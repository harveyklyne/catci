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
:func:`xgboost_learner` is the real one used by the experiments.
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np

__all__ = ["Learner", "Predict", "oracle_learner", "xgboost_learner", "fit_propensities"]


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
