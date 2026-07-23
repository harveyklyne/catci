"""Parametric double-bootstrap calibration.

Ports ``double_bootstrap_pvalue`` (``R/bootstrap_functions.R``) with two review
fixes applied:

* **2d (fall-through):** the R single-metric branch could fall through and
  return the bootstrap vector instead of a p-value. Here the scalar path always
  returns ``1 - cdf``.
* **2e (normalisation):** the R observed statistic was ranked with ``/(B+1)``
  while the inner bootstrap ranks used ``/B``, so the two stages were not exactly
  exchangeable at finite ``B``. Here both stages use the same randomised
  ``(rank - U)/(B + 1)`` rule via :func:`_randomized_cdf`.

Calibration involves RNG, so it is *not* pinned by the R fixtures (streams will
not match). It is checked instead by a calibration-under-Gaussian test.
"""

from __future__ import annotations

import numpy as np

from .bootstrap import bootstrap_T
from .criteria import ApproxChi
from .search import greedy_search

__all__ = ["double_bootstrap_pvalue", "adaptive_pvalue"]


def _randomized_cdf(x: float, x_boot: np.ndarray, rng: np.random.Generator) -> float:
    """Randomised empirical CDF of ``x`` among ``x_boot`` (R ``bootstrap_cdf``)."""
    x_boot = np.asarray(x_boot)
    B = x_boot.shape[0]
    r = np.sum(x_boot < x) + rng.uniform() * (1 + np.sum(x_boot == x))
    return float(r / (B + 1))


def double_bootstrap_pvalue(
    metrics: np.ndarray,
    metrics_boot: np.ndarray,
    rng: np.random.Generator,
) -> float:
    """Double-bootstrap p-value.

    Parameters
    ----------
    metrics : the observed criterion path, length ``M`` (``M`` search depths).
    metrics_boot : ``(M, B)`` bootstrap criterion paths.
    rng : numpy Generator for the randomised tie-breaks.
    """
    metrics = np.atleast_1d(np.asarray(metrics, dtype=float))
    metrics_boot = np.atleast_2d(np.asarray(metrics_boot, dtype=float))
    M, B = metrics_boot.shape
    if metrics.shape[0] != M:
        raise ValueError("metrics and metrics_boot have incompatible shapes.")

    if M == 1:
        # 2d: always return a p-value.
        return 1.0 - _randomized_cdf(metrics[0], metrics_boot[0], rng)

    # Inner bootstrap: turn each depth's draws into (approximately) uniform
    # q-values, then take the max over depths. 2e: same (B+1) normalisation as
    # the outer stage below.
    q_obs = np.array([_randomized_cdf(metrics[m], metrics_boot[m], rng) for m in range(M)])
    max_q_obs = float(np.max(q_obs))

    ranks = np.empty_like(metrics_boot)
    for m in range(M):
        order = np.argsort(np.argsort(metrics_boot[m], kind="stable"))  # 0-based ranks
        u = rng.uniform(size=B)
        ranks[m] = (order + 1 - u) / (B + 1)
    max_q_boot = ranks.max(axis=0)

    return 1.0 - _randomized_cdf(max_q_obs, max_q_boot, rng)


def adaptive_pvalue(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure,
    y_structure,
    n_boot: int = 100,
    criterion=None,
    rng: np.random.Generator | None = None,
) -> float:
    """Run the greedy search on the observed and bootstrap ``T`` and calibrate.

    The bootstrap draws ``T ~ N(0, Sigma)`` share the observed ``Sigma``; each is
    put through the same greedy search, and the observed criterion path is
    compared against the bootstrap paths by :func:`double_bootstrap_pvalue`.
    """
    if criterion is None:
        criterion = ApproxChi()
    if rng is None:
        rng = np.random.default_rng()

    def path(T_vec: np.ndarray) -> np.ndarray:
        return np.asarray(
            greedy_search(T_vec, Sigma, dx, dy, x_structure, y_structure, criterion).values
        )

    metrics = path(T_vector)
    boot_T = bootstrap_T(Sigma, n_boot, rng)  # (dx*dy, n_boot)
    metrics_boot = np.column_stack([path(boot_T[:, b]) for b in range(n_boot)])  # (M, n_boot)
    return double_bootstrap_pvalue(metrics, metrics_boot, rng)
