"""Test statistics as ``init`` / ``update`` / ``value`` triples.

A statistic carries just enough state to be updated cheaply after a rank-one
merge. ``ApproxChi`` is the one live statistic (Box's chi-square CDF); the
non-adaptive comparators (``euclid``, ``max``, ``mGCM``) are depth-0 value
functions -- the same code path at search depth 0, which is what the paper
claims they are. ``ExactChi`` is the reference ``ApproxChi`` approximates: the
same ``||T||^2`` scored against its exact null law, roughly 30x slower. It exists
to show the approximation costs nothing, not as a default. Deliberately kept
small: no statistic zoo.

The search evaluates one of these on every coarsening it visits, so a length-``L``
search path is a vector of ``L`` test statistics for the same null -- which is
what :mod:`catci.calibrate` then calibrates as a minP test.

Pinned by the ``approx_chi`` and ``scalar_methods`` fixtures.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import gammainc

from . import merging, quadform

__all__ = [
    "ApproxChi",
    "ChiState",
    "ExactChi",
    "ExactState",
    "approx_chi_array",
    "approx_chi_statistic",
    "euclid",
    "max_abs",
    "mgcm",
]


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
    """``pchisq(normsq / (tr2/tr), df = tr^2 / tr2)`` -- Box (1954).

    Evaluated as the regularised incomplete gamma rather than ``chi2.cdf``. The two
    are bit-for-bit identical -- ``chi2.cdf(x, df)`` *is* ``gammainc(df/2, x/2)`` --
    but ``chi2.cdf`` runs the whole ``rv_continuous`` argument-validation machinery
    on every scalar call, which costs an order of magnitude more than the special
    function itself. This is called once per candidate merge per bootstrap draw,
    so it was about half of the entire search.
    """
    return float(approx_chi_array(normsq, tr, tr2))


def approx_chi_array(normsq, tr, tr2):
    """Elementwise :func:`approx_chi_statistic`, for every candidate at once."""
    g = tr2 / tr
    h = tr ** 2 / tr2
    return gammainc(h / 2.0, normsq / g / 2.0)


# --------------------------------------------------------------------------- #
# Non-adaptive comparators (depth-0)
# --------------------------------------------------------------------------- #
def euclid(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.sqrt(np.sum(T_vector ** 2)))


def max_abs(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector)))


def mgcm(T_vector: np.ndarray, Sigma: np.ndarray) -> float:
    return float(np.max(np.abs(T_vector / np.sqrt(np.diag(Sigma)))))


# --------------------------------------------------------------------------- #
# Exact reference statistic
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ExactState:
    normsq: float
    lambdas: np.ndarray
    series: object = None  # a prebuilt quadform.RubenSeries, when one is cached


class ExactChi:
    """Exact CDF of ``||T||^2`` under ``T ~ N(0, Sigma)``: ``ApproxChi`` unapproximated.

    Same statistic as :class:`ApproxChi` -- ``||t||^2`` -- but scored against its
    true null law ``sum_j lambda_j Z_j^2`` rather than Box's two-moment match. No
    matrix inversion is involved, so this is *not* the unstable pseudo-inverse
    statistic the paper rules out.

    The cost is that ``(tr, tr2)`` no longer suffice: every candidate merge needs
    the *spectrum* of the merged covariance, and there is no rank-one update for a
    spectrum. Three ways to pay for it, all returning identical values:

    ``mode="dense"``
        Form the merged ``Sigma`` and eigendecompose it: ``O(p^3)`` per candidate.
    ``mode="lowrank"``
        Keep a factor ``C`` with ``C C^T = Sigma`` and eigendecompose the ``r x r``
        matrix ``C^T A^T A C = G + U^T V + V^T U``, whose non-zero spectrum is the
        same. ``O(r^3)`` per candidate, with ``r = (dx-1)(dy-1) < p = dx dy``.
    ``mode="cached"``
        ``lowrank`` plus a spectrum cache keyed on the resulting partition. The
        spectrum is a function of ``(Sigma, partition)`` alone -- never of ``T`` --
        so every bootstrap draw that revisits a partition reuses the work.

    Even cached, the CDF itself (Ruben's series or Imhof's integral) must be
    evaluated at every draw's ``||t||^2``, and that is ~90% of the cost -- which is
    why this stays ~30x slower than :class:`ApproxChi` however it is organised.

    Uses the search's opt-in context hooks (``wants_context``): ``begin_search``
    supplies the root ``Sigma`` the cache is valid for (a different one clears it),
    ``begin_level`` the partition the upcoming candidates start from, and
    ``update`` receives ``context=(dimension, i, j)``. Not reentrant; one search at
    a time per instance.
    """

    wants_context = True

    def __init__(self, mode: str = "cached", method: str = "auto"):
        if mode not in ("dense", "lowrank", "cached"):
            raise ValueError("mode must be 'dense', 'lowrank' or 'cached'.")
        self.mode = mode
        self.method = method
        self._C = None
        self._G = None
        self._root_Sigma = None
        self._partition = None
        self._cache: dict = {}
        self.stats = {"eigen": 0, "cdf": 0, "cache_hits": 0}

    # -- search context hooks ---------------------------------------------- #
    def begin_search(self, Sigma: np.ndarray) -> None:
        """A search is starting from ``Sigma``; cached spectra are only valid for that one."""
        if self._root_Sigma is None or not np.array_equal(self._root_Sigma, Sigma):
            self._root_Sigma = np.array(Sigma, dtype=float)
            self._cache.clear()

    def begin_level(self, partition: dict) -> None:
        """The partition the upcoming candidate merges start from."""
        self._partition = partition

    @staticmethod
    def _key(partition: dict, dimension: int, i: int, j: int):
        """Canonical key for the partition that merging positions ``(i, j)`` produces.

        Eigenvalues are invariant to the order of the groups, so the key sorts
        both levels -- two draws that reach the same grouping by different routes
        share a cache entry.
        """
        groups = {k: [list(g) for g in v] for k, v in partition.items()}
        side = "x" if dimension == 1 else "y"
        g = groups[side]
        g[i - 1] = g[i - 1] + g[j - 1]
        del g[j - 1]
        return (
            tuple(sorted(tuple(sorted(t)) for t in groups["x"])),
            tuple(sorted(tuple(sorted(t)) for t in groups["y"])),
        )

    # -- statistic protocol ------------------------------------------------ #
    def init(self, T_vector: np.ndarray, Sigma: np.ndarray) -> ExactState:
        if self.mode == "dense":
            self._C = self._G = None
            lambdas = quadform.positive_eigenvalues(Sigma)
        else:
            self._C = _psd_factor(Sigma)
            self._G = self._C.T @ self._C
            lambdas = np.linalg.eigvalsh(self._G)
            lambdas = lambdas[lambdas > 0]
        self.stats["eigen"] += 1
        return ExactState(normsq=float(np.sum(T_vector ** 2)), lambdas=lambdas)

    def update(
        self,
        state: ExactState,
        T_vector: np.ndarray,
        Sigma: np.ndarray,
        index1: np.ndarray,
        index2: np.ndarray,
        context=None,
    ) -> ExactState:
        normsq = merging.update_normsq(state.normsq, T_vector, index1, index2)

        cacheable = self.mode == "cached" and context is not None
        if cacheable and self._partition is not None and self._root_Sigma is not None:
            key = self._key(self._partition, *context)
            entry = self._cache.get(key)
            if entry is None:
                # Ruben's coefficients also depend on the spectrum alone, so when
                # the spectrum is tame enough for Ruben the whole series is cached
                # and each draw costs one dot product. Imhof has no q-free part, so
                # there only the eigenvalues can be reused.
                lambdas = self._spectrum(Sigma, index1, index2)
                if self._use_ruben(lambdas):
                    entry = (None, quadform.RubenSeries(lambdas))
                else:
                    entry = (lambdas, None)
                self._cache[key] = entry
            else:
                self.stats["cache_hits"] += 1
            return ExactState(normsq=normsq, lambdas=entry[0], series=entry[1])

        return ExactState(normsq=normsq, lambdas=self._spectrum(Sigma, index1, index2))

    def value(self, state: ExactState) -> float:
        self.stats["cdf"] += 1
        if state.series is not None:
            return state.series.cdf(state.normsq)
        return quadform.exact_cdf(state.normsq, state.lambdas, method=self.method)

    # -- internals --------------------------------------------------------- #
    def _use_ruben(self, lambdas: np.ndarray) -> bool:
        if self.method == "ruben":
            return True
        if self.method == "imhof":
            return False
        return quadform.ruben_terms(lambdas) <= quadform.RUBEN_TERM_BUDGET

    def _spectrum(self, Sigma: np.ndarray, index1, index2) -> np.ndarray:
        """Non-zero eigenvalues of the merged covariance."""
        self.stats["eigen"] += 1
        if self.mode == "dense":
            return quadform.positive_eigenvalues(merging.update_Sigma(Sigma, index1, index2))
        # A Sigma A^T and C^T A^T A C share their non-zero spectrum, and
        # A^T A = I + E + E^T with E[index1_t, index2_t] = 1, so C^T E C = U^T V.
        U = self._C[index1, :]
        V = self._C[index2, :]
        W = U.T @ V
        M = self._G + W + W.T
        vals = np.linalg.eigvalsh((M + M.T) / 2.0)
        top = vals[-1] if vals.size else 0.0
        return vals[vals > quadform._EIG_REL_TOL * top] if top > 0 else vals[:0]


def _psd_factor(Sigma: np.ndarray) -> np.ndarray:
    """``C`` with ``C C^T = Sigma``; zero and negative directions dropped."""
    Sigma = np.asarray(Sigma, dtype=float)
    Sigma = (Sigma + Sigma.T) / 2.0
    vals, vecs = np.linalg.eigh(Sigma)
    top = max(vals[-1], 0.0) if vals.size else 0.0
    keep = vals > quadform._EIG_REL_TOL * top if top > 0 else vals > 0
    return vecs[:, keep] * np.sqrt(vals[keep])
