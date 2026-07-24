"""Propensity learners: the ``Z -> P(label | Z)`` interface plus cross-fitting.

A learner is any callable ``fit(z_train, labels_train, num_class) -> predict``,
where ``predict(z_new)`` returns an ``(m, num_class)`` row-stochastic matrix.
Cross-fitting (below) turns a learner into out-of-fold propensities ``f`` / ``g``
and knows nothing about any specific ML backend.

Only the **oracle** learner is implemented here -- it returns the true
propensities and lets you separate "the test calibrates" from "the regression
fit well" (CODE_REVIEW.md 4.1). Gradient-boosting / multinomial learners are
left as a later step (optional ``[learners]`` extra); the interface is fixed so
they slot in without touching :func:`crossfit`.
"""

from __future__ import annotations

from typing import Callable, Protocol

import numpy as np

__all__ = ["Learner", "Predict", "oracle_learner", "xgboost_learner", "crossfit"]


class Predict(Protocol):
    def __call__(self, z_new: np.ndarray) -> np.ndarray: ...


Learner = Callable[[np.ndarray, np.ndarray, int], Predict]


def oracle_learner(true_probs: np.ndarray) -> Learner:
    """A learner that ignores the data and returns known true propensities.

    ``true_probs`` is the full ``(n, num_class)`` matrix indexed by row id; the
    returned predictor selects rows by the integer ids passed as ``z``. This
    keeps the oracle inside the same cross-fitting machinery as real learners.
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


def crossfit(
    z: np.ndarray,
    labels: np.ndarray,
    num_class: int,
    learner: Learner,
    nfolds: int = 5,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Out-of-fold propensity matrix ``(n, num_class)`` via ``nfolds`` cross-fitting."""
    if rng is None:
        rng = np.random.default_rng()
    z = np.asarray(z)
    labels = np.asarray(labels)
    n = labels.shape[0]

    fold = rng.integers(0, nfolds, size=n)
    out = np.empty((n, num_class), dtype=float)
    for k in range(nfolds):
        test = fold == k
        train = ~test
        predict = learner(z[train], labels[train], num_class)
        out[test] = predict(z[test])
    return out
