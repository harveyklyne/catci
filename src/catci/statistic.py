"""Test statistics as ``init`` / ``update`` / ``value`` triples.

A statistic carries just enough state to be updated cheaply after a rank-one
merge. ``ApproxChi`` is the one live statistic (Box's chi-square CDF); the
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
from scipy.stats import chi2

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

    def value(self, state: ChiState) -> float:
        return approx_chi_statistic(state.normsq, state.tr, state.tr2)


def approx_chi_statistic(normsq: float, tr: float, tr2: float) -> float:
    """``pchisq(normsq / (tr2/tr), df = tr^2 / tr2)`` -- Box (1954)."""
    g = tr2 / tr
    h = tr ** 2 / tr2
    return float(chi2.cdf(normsq / g, df=h))


# --------------------------------------------------------------------------- #
# Non-adaptive comparators (depth-0)
# --------------------------------------------------------------------------- #
def euclid(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.sqrt(np.sum(T_vector ** 2)))


def max_abs(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector)))


def mgcm(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector / np.sqrt(np.diag(Sigma)))))
