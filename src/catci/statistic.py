"""Test statistics as ``init`` / ``update`` / ``value`` triples.

A statistic carries just enough state to be updated cheaply after a rank-one
merge, or after the split that undoes one (``ApproxChi.split``, used by
:func:`~catci.search.divisive_search`). ``ApproxChi`` is the one live statistic (Box's chi-square CDF); the
non-adaptive comparators (``euclid``, ``max``, ``mGCM``) are depth-0 value
functions -- the same code path at search depth 0, which is what the paper
claims they are. Deliberately kept small: no statistic zoo.

The search evaluates one of these on every coarsening it visits, so a length-``L``
search path is a vector of ``L`` test statistics for the same null -- which is
what :mod:`catci.calibrate` then calibrates as a minP test.

Pinned by the ``approx_chi`` and ``scalar_methods`` fixtures.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammainc

from . import merging

__all__ = ["ApproxChi", "ChiState", "euclid", "max_abs", "mgcm"]


@dataclass(frozen=True)
class ChiState:
    normsq: float
    tr: float
    tr2: float


class ApproxChi:
    """Box (1954) approximate chi-square CDF statistic; larger favours rejection."""

    def init(self, T_vector: np.ndarray, Sigma: np.ndarray) -> ChiState:
        return ChiState(
            normsq=float(np.sum(T_vector ** 2)),
            tr=float(np.trace(Sigma)),
            tr2=float(np.sum(Sigma ** 2)),
        )

    def update(
        self,
        state: ChiState,
        T_vector: np.ndarray,
        Sigma: np.ndarray,
        index1: np.ndarray,
        index2: np.ndarray,
    ) -> ChiState:
        return ChiState(
            normsq=merging.update_normsq(state.normsq, T_vector, index1, index2),
            tr=merging.update_tr(state.tr, Sigma, index1, index2),
            tr2=merging.update_tr2(state.tr2, Sigma, index1, index2),
        )

    def split(
        self,
        state: ChiState,
        Ta: np.ndarray,
        Tb: np.ndarray,
        Sa: np.ndarray,
        Sb: np.ndarray,
        cols_a: np.ndarray,
        cols_b: np.ndarray,
    ) -> ChiState:
        """State after *refining* the partition -- the counterpart of :meth:`update`.

        ``update`` merges two groups and reads the current ``(T, Sigma)``; this
        splits one group in two and reads only the new rows. See
        :mod:`catci.merging` for what the arguments mean.
        """
        return ChiState(
            normsq=merging.split_normsq(state.normsq, Ta, Tb),
            tr=merging.split_tr(state.tr, Sa, cols_b),
            tr2=merging.split_tr2(state.tr2, Sa, Sb, cols_a, cols_b),
        )

    def value(self, state: ChiState) -> float:
        return approx_chi_statistic(state.normsq, state.tr, state.tr2)


def approx_chi_statistic(normsq: float, tr: float, tr2: float) -> float:
    """``pchisq(normsq / (tr2/tr), df = tr^2 / tr2)`` -- Box (1954).

    Written as the regularised lower incomplete gamma rather than
    ``chi2.cdf``, which it equals bit-for-bit: the chi-square CDF *is*
    ``gammainc(df/2, x/2)``, but reaching it through ``scipy.stats`` costs ~80x
    more per scalar call than calling it directly. Both searches evaluate this
    once per candidate merge or split, so it is the hottest line in the package.
    """
    g = tr2 / tr
    h = tr ** 2 / tr2
    return float(gammainc(h / 2.0, normsq / g / 2.0))


# --------------------------------------------------------------------------- #
# Non-adaptive comparators (depth-0)
# --------------------------------------------------------------------------- #
def euclid(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.sqrt(np.sum(T_vector ** 2)))


def max_abs(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector)))


def mgcm(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector / np.sqrt(np.diag(Sigma)))))
