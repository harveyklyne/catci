"""Propensity learners: the ``Z -> P(label | Z)`` interface.

A learner is any callable ``fit(z_train, labels_train, num_class) -> predict``,
where ``predict(z_new)`` returns an ``(m, num_class)`` row-stochastic matrix.
:func:`fit_propensities` turns a learner into the propensity matrix ``f`` / ``g``
and knows nothing about any specific ML backend.

The propensities are fitted on the full sample -- there is no cross-fitting.
The test calibrates on the fitted propensities themselves, so the sample-splitting
machinery bought nothing and cost an nfolds-fold slowdown plus a fold-assignment
RNG stream in every replicate.

The **oracle** learner returns the true propensities and lets you separate "the
test calibrates" from "the regression fit well";
:func:`xgboost_learner` and :func:`mlp_learner` are the real ones used by the experiments.
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np

__all__ = [
    "Learner",
    "Predict",
    "oracle_learner",
    "xgboost_learner",
    "mlp_learner",
    "fit_propensities",
]


class Predict(Protocol):
    def __call__(self, z_new: np.ndarray) -> np.ndarray: ...


Learner = Callable[[np.ndarray, np.ndarray, int], Predict]


def oracle_learner(true_probs: np.ndarray) -> Learner:
    """A learner that ignores the data and returns known true propensities.

    ``true_probs`` is the full ``(n, num_class)`` matrix indexed by row id; the
    returned predictor selects rows by the integer ids passed as ``z``. This
    keeps the oracle behind the same interface as real learners.
    """
    true_probs = np.asarray(true_probs, dtype=float)

    def fit(z_train: np.ndarray, labels_train: np.ndarray, num_class: int) -> Predict:
        def predict(z_new: np.ndarray) -> np.ndarray:
            return true_probs[np.asarray(z_new, dtype=int)]

        return predict

    return fit


def xgboost_learner(params: dict) -> Learner:
    """Gradient-boosted multinomial propensity learner (ports R ``fit_xgboost``).

    ``params`` matches the R tuning JSONs: ``eta``, ``max_depth`` (``"max.depth"``
    accepted), ``gamma``, ``nrounds``. Objective ``multi:softprob``,
    ``nthread=1`` so outer parallelism owns the cores (not xgboost).
    """
    import xgboost as xgb

    p = {("max_depth" if k == "max.depth" else k): v for k, v in params.items()}
    nrounds = int(p.pop("nrounds"))
    p.update(objective="multi:softprob", eval_metric="mlogloss", nthread=1)

    def fit(z_train: np.ndarray, labels_train: np.ndarray, num_class: int) -> Predict:
        dtrain = xgb.DMatrix(np.asarray(z_train, dtype=float), label=np.asarray(labels_train) - 1)
        booster = xgb.train({**p, "num_class": num_class}, dtrain, num_boost_round=nrounds)

        def predict(z_new: np.ndarray) -> np.ndarray:
            dnew = xgb.DMatrix(np.asarray(z_new, dtype=float))
            return booster.predict(dnew).reshape(-1, num_class)

        return predict

    return fit


def fit_propensities(
    z: np.ndarray,
    labels: np.ndarray,
    num_class: int,
    learner: Learner,
) -> np.ndarray:
    """Propensity matrix ``(n, num_class)``: fit ``learner`` on the full sample.

    No cross-fitting -- the learner is trained on all of ``(z, labels)`` and
    predicts back on the same ``z``.
    """
    z = np.asarray(z)
    labels = np.asarray(labels)
    predict = learner(z, labels, num_class)
    return predict(z)


def mlp_learner(params: dict) -> Learner:
    """Multi-layer perceptron multinomial propensity learner (sklearn ``MLPClassifier``).

    The alternative to :func:`xgboost_learner`. Boosted stumps approximate a
    smooth propensity surface by a staircase; a net with a smooth activation
    approximates it smoothly, which is the whole reason to try one here.

    ``params`` mirrors the tuning JSONs' ``mlp`` block: ``hidden_layer_sizes``
    (list or tuple), ``alpha`` (L2 penalty), ``learning_rate_init``,
    ``max_iter``, ``activation``, ``random_state``. Z is standardised before the
    net sees it -- an unscaled input is the usual reason an MLP silently
    underfits, and the scaler is fitted on the training rows only.

    Threads are pinned to one inside ``fit``/``predict`` (the BLAS analogue of
    xgboost's ``nthread=1``) so the runner's process pool owns the cores.
    """
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from threadpoolctl import threadpool_limits

    p = dict(params)
    hidden = p.pop("hidden_layer_sizes", (64, 64))
    kwargs = dict(
        hidden_layer_sizes=tuple(hidden),
        activation=p.pop("activation", "tanh"),
        alpha=p.pop("alpha", 1e-2),
        learning_rate_init=p.pop("learning_rate_init", 1e-3),
        max_iter=int(p.pop("max_iter", 500)),
        random_state=int(p.pop("random_state", 0)),
        solver="adam",
    )
    kwargs.update(p)  # anything else passes through to MLPClassifier

    def fit(z_train: np.ndarray, labels_train: np.ndarray, num_class: int) -> Predict:
        z_train = np.asarray(z_train, dtype=float)
        labels_train = np.asarray(labels_train)
        scaler = StandardScaler().fit(z_train)
        clf = MLPClassifier(**kwargs)
        with threadpool_limits(limits=1):
            import warnings

            with warnings.catch_warnings():
                # adam hitting max_iter is a tuning signal, not a per-fit event
                warnings.simplefilter("ignore")
                clf.fit(scaler.transform(z_train), labels_train)

        # clf.classes_ may be a strict subset of 1..num_class; scatter into the
        # full width so f/g always have num_class columns like the oracle does.
        cols = np.asarray(clf.classes_, dtype=int) - 1

        def predict(z_new: np.ndarray) -> np.ndarray:
            z_new = np.asarray(z_new, dtype=float)
            with threadpool_limits(limits=1):
                proba = clf.predict_proba(scaler.transform(z_new))
            if proba.shape[1] == num_class and np.array_equal(cols, np.arange(num_class)):
                return proba
            out = np.zeros((z_new.shape[0], num_class))
            out[:, cols] = proba
            return out

        return predict

    return fit
