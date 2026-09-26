"""Cross-validated hyperparameter tuning for the propensity learners.

Replaces the cluster-specific R tuner (``data-raw/tuning/`` on tag
``r-frozen-oracle``), which fitted each grid point on 100 simulated
``n_tr = 800`` training sets and scored it on fresh ``n_te = 5000`` test sets.
Here the score is K-fold cross-validated multinomial log-loss, so the same code
tunes a simulated setting and a real dataset:

* **one dataset** (an application) -- ordinary K-fold CV;
* **several datasets** (independent replicates of a simulated setting) -- K-fold
  CV *within* each replicate, with the validation curves averaged over every
  (replicate, fold) pair. Each booster only ever sees one replicate's in-fold
  rows, so the training size is ``(K-1)/K * n`` -- at ``n = 1000, K = 5`` that is
  the R tuner's ``n_tr = 800``.

``nrounds`` is not a grid axis: every (replicate, fold) booster is advanced in
lockstep and the mean validation curve read off per round, with early stopping
on that mean. So a grid point costs one boosting run per fold, not one per
candidate round count.

The cost of multiclass boosting is linear in ``num_class`` (one tree per class
per round), which is what makes tuning at ``d`` in the hundreds expensive: see
:data:`FAST_GRID`.
"""

from __future__ import annotations

import itertools
import time
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

__all__ = [
    "R_GRID",
    "FAST_GRID",
    "TuneResult",
    "kfold_indices",
    "xgboost_cv_curve",
    "tune_xgboost",
]

# The grid the R tuner searched (eta fixed at 0.01, up to 1000 rounds).
R_GRID = dict(eta=[0.01], max_depth=[1, 2, 3, 4], gamma=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])

# A cheaper grid for large ``num_class``. A 10x larger step needs ~10x fewer
# rounds, and every R-tuned setting chose depth 1 with gamma in [1, 2].
FAST_GRID = dict(eta=[0.1], max_depth=[1, 2, 3], gamma=[0.0, 1.0, 2.0])


@dataclass
class TuneResult:
    """The selected hyperparameters plus the full CV table, for provenance."""

    params: dict  # eta, max_depth, gamma, nrounds -- accepted by xgboost_learner
    cv_logloss: float
    table: list = field(default_factory=list)  # one dict per grid point
    seconds: float = 0.0

    def to_json(self) -> dict:
        """The on-disk ``tune_<setting>_results.json`` layout (``"xgb"`` key as in R)."""
        p = dict(self.params)
        return {
            "xgb": {"eta": p["eta"], "max.depth": p["max_depth"], "gamma": p["gamma"], "nrounds": p["nrounds"]},
            "cv": {"logloss": self.cv_logloss, "seconds": round(self.seconds, 1), "table": self.table},
        }


def kfold_indices(n: int, n_folds: int, rng: np.random.Generator) -> list[tuple[np.ndarray, np.ndarray]]:
    """Random (not stratified) K-fold split of ``range(n)`` as ``(train, test)`` pairs.

    Not stratified on purpose: at ``d`` in the hundreds many classes have fewer
    than ``n_folds`` members, and a learner that cannot cope with a class missing
    from its training fold is a learner the test cannot use either.
    """
    if not 2 <= n_folds <= n:
        raise ValueError("Need 2 <= n_folds <= n.")
    fold = rng.permutation(n) % n_folds
    return [(np.flatnonzero(fold != k), np.flatnonzero(fold == k)) for k in range(n_folds)]


def _stack(datasets: Sequence[tuple[np.ndarray, np.ndarray]], n_folds: int, rng):
    """Stack replicates into one design and build per-replicate (train, test) folds."""
    zs, ys, folds, offset = [], [], [], 0
    for z, labels in datasets:
        z = np.asarray(z, dtype=float)
        labels = np.asarray(labels)
        if z.ndim == 1:
            z = z[:, None]
        n = z.shape[0]
        for tr, te in kfold_indices(n, n_folds, rng):
            folds.append((tr + offset, te + offset))
        zs.append(z)
        ys.append(labels)
        offset += n
    return np.vstack(zs), np.concatenate(ys), folds


def xgboost_cv_curve(
    dmatrix,
    folds,
    num_class: int,
    *,
    eta: float,
    max_depth: int,
    gamma: float,
    max_rounds: int,
    early_stopping_rounds: int | None,
    nthread: int,
    seed: int = 0,
) -> np.ndarray:
    """Mean validation mlogloss per boosting round, averaged over ``folds``.

    Truncated at the early-stopping point when ``early_stopping_rounds`` is set.
    """
    import xgboost as xgb

    params = dict(
        objective="multi:softprob", eval_metric="mlogloss", num_class=num_class,
        eta=eta, max_depth=max_depth, gamma=gamma, nthread=nthread, seed=seed,
    )
    res = xgb.cv(
        params, dmatrix, num_boost_round=max_rounds, folds=folds,
        early_stopping_rounds=early_stopping_rounds, as_pandas=False, verbose_eval=False,
    )
    return np.asarray(res["test-mlogloss-mean"], dtype=float)


def tune_xgboost(
    datasets: Sequence[tuple[np.ndarray, np.ndarray]] | tuple[np.ndarray, np.ndarray],
    num_class: int,
    *,
    grid: dict | None = None,
    n_folds: int = 5,
    max_rounds: int = 1000,
    early_stopping_rounds: int | None = 50,
    nthread: int = 1,
    rng: np.random.Generator | None = None,
    verbose: bool = False,
) -> TuneResult:
    """Choose ``(eta, max_depth, gamma, nrounds)`` by K-fold CV log-loss.

    Parameters
    ----------
    datasets : one ``(z, labels)`` pair, or a list of them (independent
        replicates of a simulated setting). ``labels`` are 1-based in
        ``{1, ..., num_class}``, matching :func:`catci.learners.fit_propensities`.
    grid : ``{"eta": [...], "max_depth": [...], "gamma": [...]}``; defaults to
        :data:`R_GRID`. ``nrounds`` is read off the CV curve, not gridded.
    early_stopping_rounds : stop a grid point once the mean validation loss has
        not improved for this many rounds (``None`` runs all ``max_rounds``, as
        the R tuner did). If the best round is within this many of
        ``max_rounds``, the curve was still falling and ``max_rounds`` is binding;
        that is recorded as ``hit_max_rounds`` in the table.
    """
    import xgboost as xgb

    if rng is None:
        rng = np.random.default_rng()
    if isinstance(datasets, tuple):
        datasets = [datasets]
    grid = dict(R_GRID if grid is None else grid)

    z, labels, folds = _stack(datasets, n_folds, rng)
    if labels.min() < 1 or labels.max() > num_class:
        raise ValueError("labels must lie in {1, ..., num_class}.")
    dmatrix = xgb.DMatrix(z, label=labels - 1)
    seed = int(rng.integers(2**31 - 1))

    t0 = time.time()
    table, best = [], None
    for eta, depth, gamma in itertools.product(grid["eta"], grid["max_depth"], grid["gamma"]):
        t = time.time()
        curve = xgboost_cv_curve(
            dmatrix, folds, num_class, eta=eta, max_depth=depth, gamma=gamma,
            max_rounds=max_rounds, early_stopping_rounds=early_stopping_rounds,
            nthread=nthread, seed=seed,
        )
        k = int(np.argmin(curve))
        row = dict(
            eta=eta, max_depth=depth, gamma=gamma, nrounds=k + 1, logloss=float(curve[k]),
            # the curve ran to max_rounds and was still (near) its minimum there
            hit_max_rounds=bool(len(curve) == max_rounds and k + 1 > max_rounds - (early_stopping_rounds or 1)),
            seconds=round(time.time() - t, 1),
        )
        table.append(row)
        if verbose:
            print(row, flush=True)
        # strict '<' keeps the first grid point on ties (grid order is the tie-break)
        if best is None or row["logloss"] < best["logloss"]:
            best = row

    params = {k: best[k] for k in ("eta", "max_depth", "gamma", "nrounds")}
    return TuneResult(params=params, cv_logloss=best["logloss"], table=table, seconds=time.time() - t0)
