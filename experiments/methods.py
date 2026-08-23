"""Method registry: name -> p-value on a fitted dataset.

Covers the adaptive catci tests (tree/ordinal/max/euclid/mGCM, calibrated by the
shared double bootstrap) and the competitors (ankan, chi_sq, multinomial). This
replaces the R split across ``formulate_statistics`` (three unconditional
competitors) and ``evaluate_sim`` (the calibrated ones): here every method is an
explicit registry entry the runner asks for, so an expensive competitor is paid
for only when requested.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2, norm

from catci.bootstrap import bootstrap_T
from catci.calibrate import double_bootstrap_pvalue
from catci.criteria import ApproxChi, euclid, max_abs, mgcm
from catci.gcm import form_t_sigma
from catci.search import greedy_search
from catci.structure import Ordinal, Saturated, Tree

ADAPTIVE = ("tree", "ordinal", "greedy", "max", "euclid", "mGCM")
COMPETITORS = ("ankan", "chi_sq", "multinomial")


@dataclass
class Fitted:
    """Everything a method might need for one replicate."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    f: np.ndarray  # fitted P(X|Z)
    g: np.ndarray  # fitted P(Y|Z)
    dx: int
    dy: int
    T_vector: np.ndarray
    Sigma: np.ndarray

    @classmethod
    def build(cls, x, y, z, f, g, dx, dy, normalise=False):
        ts = form_t_sigma(x, y, f, g, normalise=normalise)
        return cls(x, y, z, f, g, dx, dy, ts.T_vector, ts.Sigma)


# --------------------------------------------------------------------------- #
# Adaptive methods (shared bootstrap draws, matching R evaluate_sim)
# --------------------------------------------------------------------------- #
def _metric_fn(name, dx, dy, criterion):
    """A function (T_vector, Sigma) -> criterion (scalar for depth-0, vector for a search)."""
    search_structs = {
        "tree": lambda: (Tree.binary(dx), Tree.binary(dy)),
        "ordinal": lambda: (Ordinal(), Ordinal()),
        "greedy": lambda: (Saturated(), Saturated()),
    }
    if name in search_structs:
        xs, ys = search_structs[name]()

        def fn(T_vec, Sigma):
            return np.asarray(greedy_search(T_vec, Sigma, dx, dy, xs, ys, criterion).values)

        return fn
    scalar = {"max": max_abs, "euclid": euclid, "mGCM": mgcm}[name]
    return lambda T_vec, Sigma: float(scalar(T_vec, Sigma))


def adaptive_pvalues(fitted: Fitted, method_names, n_boot: int, rng: np.random.Generator) -> dict:
    """P-values for the requested adaptive methods, sharing one set of bootstrap draws."""
    criterion = ApproxChi()
    Sigma = fitted.Sigma
    boot_T = bootstrap_T(Sigma, n_boot, rng)  # (p, n_boot)

    out = {}
    for name in method_names:
        fn = _metric_fn(name, fitted.dx, fitted.dy, criterion)
        observed = np.atleast_1d(fn(fitted.T_vector, Sigma))
        metrics_boot = np.column_stack([np.atleast_1d(fn(boot_T[:, b], Sigma)) for b in range(n_boot)])
        out[name] = double_bootstrap_pvalue(observed, metrics_boot, rng)
    return out


# --------------------------------------------------------------------------- #
# Competitors
# --------------------------------------------------------------------------- #
def _li_residual(x, f):
    """Li & Shepherd (2012) residual (ports R Li_residual)."""
    x = np.asarray(x)
    n, dx = f.shape
    xcdf = np.cumsum(f, axis=1)
    xcdf = np.hstack([np.zeros((n, 1)), xcdf[:, : dx - 1], np.ones((n, 1))])
    # Pr(X < x | Z) - Pr(X > x | Z)
    return xcdf[np.arange(n), x - 1] - (1 - xcdf[np.arange(n), x])


def ankan(fitted: Fitted, rng=None) -> float:
    """Ankan & Textor (2022) residual test."""
    rx = _li_residual(fitted.x, fitted.f)
    ry = _li_residual(fitted.y, fitted.g)
    n = fitted.x.shape[0]
    prod = rx * ry
    stat = (1 / np.sqrt(n)) * prod.sum() / np.std(prod, ddof=1)
    return float(1 - chi2.cdf(stat ** 2, df=1))


def chi_sq(fitted: Fitted, rng=None) -> float:
    """Pseudo-inverse chi-square test with (dx-1)(dy-1) df (ports R chi_sq)."""
    T = fitted.T_vector
    stat = float(T @ np.linalg.pinv(fitted.Sigma) @ T)
    return float(1 - chi2.cdf(stat, df=(fitted.dx - 1) * (fitted.dy - 1)))


def multinomial(fitted: Fitted, rng=None) -> float:
    """Nested-multinomial likelihood-ratio test y ~ z vs y ~ factor(x) + z.

    sklearn stand-in for R's ``nnet::multinom`` + ``anova``; the LR statistic is
    ``2*(ll_full - ll_reduced)`` with ``(dy-1)(dx-1)`` df.
    """
    from sklearn.linear_model import LogisticRegression

    x, y, z = fitted.x, fitted.y, fitted.z
    dx, dy = fitted.dx, fitted.dy
    n = x.shape[0]
    xoh = np.zeros((n, dx))
    xoh[np.arange(n), x - 1] = 1.0
    Xreduced = z
    Xfull = np.hstack([xoh[:, 1:], z])  # drop one X dummy (reference level)

    def loglik(features):
        clf = LogisticRegression(max_iter=2000, C=1e6, tol=1e-6)
        clf.fit(features, y)
        proba = clf.predict_proba(features)
        idx = np.searchsorted(np.unique(y), y)
        return np.sum(np.log(np.clip(proba[np.arange(n), idx], 1e-12, None)))

    stat = 2 * (loglik(Xfull) - loglik(Xreduced))
    df = (dy - 1) * (dx - 1)
    return float(1 - chi2.cdf(max(stat, 0.0), df=df))


COMPETITOR_FNS = {"ankan": ankan, "chi_sq": chi_sq, "multinomial": multinomial}


def competitor_pvalues(fitted: Fitted, names, rng: np.random.Generator) -> dict:
    return {name: COMPETITOR_FNS[name](fitted, rng) for name in names}
