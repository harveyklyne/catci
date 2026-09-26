"""Bootstrap calibration of the search path, as a minP test.

The search returns a path of ``L`` test statistics -- one per coarsening it
visits (:mod:`catci.statistic`). Each is a valid test statistic for the *same*
null, so each has a valid bootstrap p-value; the problem is that there are ``L``
of them and we want to use whichever is most extreme. That is a plain multiple
testing problem, and this module solves it the standard way: **minP** calibrated
by resampling (Westfall & Young 1993; Romano & Wolf 2005). The first stage is
Beran (1988) prepivoting -- it replaces each statistic by its bootstrap p-value,
putting the ``L`` levels on a common scale -- and the second stage calibrates
their minimum against the bootstrap paths.

The construction pools the observed path with the ``B`` bootstrap paths into one
set of ``B + 1`` exchangeable paths and applies one identical leave-one-out rule
to every member, so ``p`` is exactly uniform on ``{1/(B+1), ..., 1}`` under the
null at finite ``B`` -- not merely close. See :func:`_loo_pvalues`.

Two consequences of the ``1/(B+1)`` grid are worth knowing:

* **Ties must be broken at random.** Each level puts exactly one path at the floor
  ``1/(B+1)``, so up to ``L`` paths tie there and the second stage is deciding
  between them. Counting ties conservatively instead makes the test unable to
  reject at all once ``L`` is comparable to ``alpha * (B + 1)``.
* **``n_boot`` bounds the power, not just the resolution.** At ``dx = dy = 8``
  (``L = 13``) and ``B = 100`` the minP test is exact but recovers only about
  half the power it reaches at ``B = 1000``. Prefer ``n_boot`` in the hundreds.

Calibration involves RNG, so it is *not* pinned by the R fixtures (streams will
not match). It is checked instead by ``tests/test_calibrate.py``, which feeds the
calibration exchangeable paths directly and asserts uniformity level by level.
"""

from __future__ import annotations

import numpy as np

from .bootstrap import bootstrap_T
from .search import MergeSearch
from .statistic import ApproxChi

__all__ = ["double_bootstrap_pvalue", "bonferroni_pvalue", "adaptive_pvalue"]


def _loo_pvalues(statistics: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Leave-one-out Monte-Carlo p-values within an exchangeable pool.

    ``statistics`` is ``(L, K)``: ``K`` exchangeable paths, ``L`` statistics each,
    larger meaning more evidence against the null. Returns ``(L, K)`` with

        ``p[l, j] = (1 + #{k != j : statistics[l, k] >= statistics[l, j]}) / K``

    i.e. every path is ranked against the *other* ``K - 1``, by one rule. Ties are
    broken at random, so each row is a random permutation of ``{1/K, ..., 1}`` and
    each ``p[l, j]`` is exactly uniform on that grid (Dwass 1957, Barnard 1963).

    Because the rule is symmetric in the pool, any statistic computed from the
    columns of the result is exchangeable whenever the input columns are -- which
    is what makes the two-stage p-value below exact.
    """
    K = statistics.shape[-1]
    # Sort ascending on the statistic, with fresh uniforms as the tie-break key.
    order = np.lexsort((rng.random(statistics.shape), statistics), axis=-1)
    rank = np.argsort(order, axis=-1) + 1  # 1 = smallest statistic in its row
    return (K + 1 - rank) / K


def _pool(statistics: np.ndarray, statistics_boot: np.ndarray) -> np.ndarray:
    """Stack the observed path as column 0 of the ``(L, B+1)`` exchangeable pool."""
    statistics = np.atleast_1d(np.asarray(statistics, dtype=float))
    statistics_boot = np.atleast_2d(np.asarray(statistics_boot, dtype=float))
    if statistics.shape[0] != statistics_boot.shape[0]:
        raise ValueError("statistics and statistics_boot have incompatible shapes.")
    return np.column_stack([statistics, statistics_boot])


def double_bootstrap_pvalue(
    statistics: np.ndarray,
    statistics_boot: np.ndarray,
    rng: np.random.Generator,
) -> float:
    """minP p-value for a path of ``L`` statistics, calibrated by double bootstrap.

    Parameters
    ----------
    statistics : the observed statistic path, length ``L`` (``L`` search depths).
    statistics_boot : ``(L, B)`` bootstrap statistic paths.
    rng : numpy Generator for the randomised tie-breaks.

    Notes
    -----
    Stage 1 prepivots every path -- observed and bootstrap alike -- into marginal
    p-values; stage 2 is the same rule applied to the negated minima, which is
    exactly ``(1 + #{b : P_min[b] <= P_min[obs]}) / (B + 1)``. At ``L = 1`` the two
    stages cancel and this reduces to the ordinary single-bootstrap p-value.
    """
    p = _loo_pvalues(_pool(statistics, statistics_boot), rng)  # stage 1: prepivot
    p_min = p.min(axis=0)  # the minP statistic, one per path
    return float(_loo_pvalues(-p_min[None, :], rng)[0, 0])  # stage 2: same rule


def bonferroni_pvalue(
    statistics: np.ndarray,
    statistics_boot: np.ndarray,
    rng: np.random.Generator,
) -> float:
    """Simple FWER control over the same path: ``min(1, L * min_l p_l)``.

    The comparator for :func:`double_bootstrap_pvalue`. It shares stage 1 and
    skips the calibration, paying for the dependence between levels with a
    Bonferroni factor instead. Its resolution floor is ``L / (B + 1)``, so unlike
    minP it cannot reject at level ``alpha`` at all unless ``B >= L / alpha``.

    :func:`double_bootstrap_pvalue` never returns a larger value than this on the
    same input -- see ``MINP.md`` -- so the gap between them measures what the
    calibration step buys.
    """
    pool = _pool(statistics, statistics_boot)
    p_obs = _loo_pvalues(pool, rng)[:, 0]  # stage 1, observed path only
    return float(min(1.0, pool.shape[0] * p_obs.min()))


def adaptive_pvalue(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure,
    y_structure,
    n_boot: int = 100,
    statistic=None,
    rng: np.random.Generator | None = None,
    search=None,
) -> float:
    """Run the label search on the observed and bootstrap ``T`` and calibrate.

    The bootstrap draws ``T ~ N(0, Sigma)`` share the observed ``Sigma``; each is
    put through the same search, and the observed statistic path is compared
    against the bootstrap paths by :func:`double_bootstrap_pvalue`.

    ``search`` picks the direction -- :class:`~catci.search.MergeSearch` (the
    default, the paper's Algorithm 1) or :class:`~catci.search.SplitSearch`. It
    is asked to ``prepare`` against the shared ``Sigma`` once, so a direction
    with per-``Sigma`` setup pays for it here rather than per draw.
    """
    if statistic is None:
        statistic = ApproxChi()
    if rng is None:
        rng = np.random.default_rng()
    if search is None:
        search = MergeSearch()

    path = search.prepare(Sigma, dx, dy, x_structure, y_structure, statistic)

    statistics = np.asarray(path(T_vector).values)
    boot_T = bootstrap_T(Sigma, n_boot, rng)  # (dx*dy, n_boot)
    statistics_boot = np.column_stack(
        [path(boot_T[:, b]).values for b in range(n_boot)]
    )  # (L, n_boot)
    return double_bootstrap_pvalue(statistics, statistics_boot, rng)
