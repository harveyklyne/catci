"""Method registry: name -> p-value on a fitted dataset.

Covers the adaptive catci tests (tree/ordinal/max/euclid/mGCM, calibrated by the
shared minP double bootstrap), their ``_bonf`` counterparts (the same statistic
path under simple FWER control, as a comparator for the minP calibration) and the
competitors (ankan, chi_sq, multinomial). This
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
from catci.calibrate import bonferroni_pvalue, double_bootstrap_pvalue
from catci.statistic import euclid, max_abs, mgcm
from catci.gcm import form_t_sigma
from catci.search import greedy_search_paths
from catci.structure import Ordinal, Saturated, Tree

SEARCHES = ("tree", "ordinal", "greedy")
# Each search also has a `<name>_bonf` variant: same statistic path, simple FWER
# control instead of the minP calibration. Note its resolution floor is L/(n_boot+1),
# so it cannot reject at alpha unless n_boot >= L/alpha (L = dx + dy - 3).
ADAPTIVE = SEARCHES + tuple(f"{s}_bonf" for s in SEARCHES) + ("max", "euclid", "mGCM")
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
def _statistic_paths(name, dx, dy, T, Sigma, n_jobs=1):
    """Statistic paths ``(L, B)`` for every column of ``T`` (``L = 1`` for depth-0)."""
    search_structs = {
        "tree": lambda: (Tree.binary(dx), Tree.binary(dy)),
        "ordinal": lambda: (Ordinal(), Ordinal()),
        "greedy": lambda: (Saturated(), Saturated()),
    }
    name = name.removesuffix("_bonf")  # the calibration differs, the statistic path does not
    if name in search_structs:
        xs, ys = search_structs[name]()
        return greedy_search_paths(T, Sigma, dx, dy, xs, ys, n_jobs=n_jobs)
    scalar = {"max": max_abs, "euclid": euclid, "mGCM": mgcm}[name]
    return np.array([[scalar(T[:, b], Sigma) for b in range(T.shape[1])]])


def adaptive_pvalues(fitted: Fitted, method_names, n_boot: int, rng: np.random.Generator,
                     n_jobs: int = 1) -> dict:
    """P-values for the requested adaptive methods, sharing one set of bootstrap draws."""
    Sigma = fitted.Sigma
    boot_T = bootstrap_T(Sigma, n_boot, rng)  # (p, n_boot)
    T_all = np.column_stack([fitted.T_vector, boot_T])  # observed is column 0

    out, cache = {}, {}
    for name in method_names:
        base = name.removesuffix("_bonf")
        if base not in cache:  # `x` and `x_bonf` share one statistic path
            cache[base] = _statistic_paths(base, fitted.dx, fitted.dy, T_all, Sigma, n_jobs)
        paths = cache[base]
        calibrate = bonferroni_pvalue if name.endswith("_bonf") else double_bootstrap_pvalue
        out[name] = calibrate(paths[:, 0], paths[:, 1:], rng)
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
