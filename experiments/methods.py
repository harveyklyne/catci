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
from catci.statistic import ApproxChi, euclid, max_abs, mgcm
from catci.gcm import form_t_sigma
from catci.search import beam_search, evaluate_path, greedy_search, random_merges
from catci.structure import Ordinal, Saturated, Tree

SEARCHES = ("tree", "ordinal", "greedy")
# Each search also has a `<name>_bonf` variant: same statistic path, simple FWER
# control instead of the minP calibration. Note its resolution floor is L/(n_boot+1),
# so it cannot reject at alpha unless n_boot >= L/alpha (L = dx + dy - 3).
#
# TODO item 7 variants, same structures, same minP calibration:
#   `<name>_random` -- one merge path drawn at random (ignoring the data), shared by
#                      the observed and every bootstrap draw;
#   `<name>_split`  -- greedy path chosen on a random half of the rows, then
#                      evaluated and calibrated on the other half only;
#   `greedy_beam<w>` -- beam search of width w over all pairs (w = 1 is `greedy`).
# See SEARCH.md for the Gaussian-limit power study of all of these.
ALTERNATIVES = tuple(f"{s}_{v}" for s in SEARCHES for v in ("random", "split")) + ("greedy_beam5",)
ADAPTIVE = (SEARCHES + tuple(f"{s}_bonf" for s in SEARCHES) + ALTERNATIVES
            + ("max", "euclid", "mGCM"))
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
def _structures(name, dx, dy):
    return {
        "tree": lambda: (Tree.binary(dx), Tree.binary(dy)),
        "ordinal": lambda: (Ordinal(), Ordinal()),
        "greedy": lambda: (Saturated(), Saturated()),
    }[name]()


def _statistic_fn(name, dx, dy, statistic, rng):
    """A function (T_vector, Sigma) -> statistic (scalar for depth-0, vector for a search)."""
    name = name.removesuffix("_bonf")  # the calibration differs, the statistic path does not
    base, _, variant = name.partition("_")
    if base in SEARCHES:
        xs, ys = _structures(base, dx, dy)
        if variant == "random":
            merges = random_merges(dx, dy, xs, ys, rng)  # drawn once: fixed across draws

            def fn(T_vec, Sigma):
                return np.asarray(evaluate_path(T_vec, Sigma, dx, dy, merges, statistic).values)

        elif variant.startswith("beam"):
            width = int(variant.removeprefix("beam"))

            def fn(T_vec, Sigma):
                return np.asarray(beam_search(T_vec, Sigma, dx, dy, xs, ys, width, statistic).values)

        else:

            def fn(T_vec, Sigma):
                return np.asarray(greedy_search(T_vec, Sigma, dx, dy, xs, ys, statistic).values)

        return fn
    scalar = {"max": max_abs, "euclid": euclid, "mGCM": mgcm}[name]
    return lambda T_vec, Sigma: float(scalar(T_vec, Sigma))


def split_pvalue(fitted: Fitted, name: str, n_boot: int, rng: np.random.Generator,
                 normalise: bool = False) -> float:
    """Sample-split search: choose the greedy path on half A, test it on half B only.

    The propensities are the full-sample ones in ``fitted`` (no refitting per half),
    so the halves are independent only up to that shared fit -- exact with oracle
    propensities, asymptotically so otherwise.
    """
    base = name.removesuffix("_split")
    xs, ys = _structures(base, fitted.dx, fitted.dy)
    statistic = ApproxChi()
    n = fitted.x.shape[0]
    perm = rng.permutation(n)
    A, B = perm[: n // 2], perm[n // 2:]

    def half(rows):
        return form_t_sigma(fitted.x[rows], fitted.y[rows], fitted.f[rows], fitted.g[rows],
                            normalise=normalise)

    ts_A, ts_B = half(A), half(B)
    merges = greedy_search(ts_A.T_vector, ts_A.Sigma, fitted.dx, fitted.dy, xs, ys, statistic).merges

    def fn(T_vec):
        return np.asarray(evaluate_path(T_vec, ts_B.Sigma, fitted.dx, fitted.dy, merges, statistic).values)

    boot_T = bootstrap_T(ts_B.Sigma, n_boot, rng)
    statistics_boot = np.column_stack([fn(boot_T[:, b]) for b in range(n_boot)])
    return double_bootstrap_pvalue(fn(ts_B.T_vector), statistics_boot, rng)


def adaptive_pvalues(fitted: Fitted, method_names, n_boot: int, rng: np.random.Generator) -> dict:
    """P-values for the requested adaptive methods, sharing one set of bootstrap draws."""
    statistic = ApproxChi()
    Sigma = fitted.Sigma
    boot_T = bootstrap_T(Sigma, n_boot, rng)  # (p, n_boot)

    out = {}
    for name in method_names:
        if name.endswith("_split"):  # its own half-sample (T, Sigma) and draws
            out[name] = split_pvalue(fitted, name, n_boot, rng)
            continue
        fn = _statistic_fn(name, fitted.dx, fitted.dy, statistic, rng)
        calibrate = bonferroni_pvalue if name.endswith("_bonf") else double_bootstrap_pvalue
        statistics = np.atleast_1d(fn(fitted.T_vector, Sigma))
        statistics_boot = np.column_stack(
            [np.atleast_1d(fn(boot_T[:, b], Sigma)) for b in range(n_boot)]
        )
        out[name] = calibrate(statistics, statistics_boot, rng)
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
