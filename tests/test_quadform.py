"""The exact weighted-chi-square CDF, and the exact statistic built on it.

There is no R oracle for any of this -- it is new -- so the pins are the two
things that can be checked independently: the analytic chi-square special case
(all eigenvalues equal), and agreement between two algorithms whose hard cases
are complementary (Ruben slows as the spectrum spreads, Imhof as the number of
non-zero eigenvalues falls).
"""

import numpy as np
import pytest
from scipy.stats import chi2

from catci import merging, quadform
from catci.statistic import ExactChi
from catci.search import greedy_search
from catci.structure import Ordinal, Saturated, Tree

# Spectra the merge search actually visits: well-conditioned at depth 0, then
# increasingly spread and increasingly short.
SPECTRA = {
    "depth0": np.linspace(0.86, 1.81, 49),
    "one_merge": np.concatenate([np.linspace(0.88, 2.57, 42), np.geomspace(4.9e-6, 2.8e-2, 7)]),
    "mid_illcond": np.concatenate([np.linspace(0.9, 2.6, 38), np.geomspace(6.8e-7, 3e-2, 11)]),
    "deep_r12": np.geomspace(0.039, 11.0, 12),
    "deep_r6": np.geomspace(0.119, 14.3, 6),
    "deep_r4": np.geomspace(0.121, 7.4, 4),
    "deep_r2": np.array([0.3, 5.1]),
}


@pytest.mark.parametrize("r", [1, 2, 3, 4, 6, 12, 25, 49])
def test_imhof_matches_analytic_chisq(r):
    """Equal eigenvalues: sum_j lam Z_j^2 == lam * chi^2_r, exactly."""
    lam = np.full(r, 0.37)
    for q in np.linspace(0.02, 5 * r * 0.37, 15):
        assert quadform.imhof_cdf(q, lam) == pytest.approx(chi2.cdf(q / 0.37, r), abs=1e-7)


@pytest.mark.parametrize("r", [1, 2, 3, 4, 6, 12, 25, 49])
def test_ruben_matches_analytic_chisq(r):
    lam = np.full(r, 0.37)
    for q in np.linspace(0.02, 5 * r * 0.37, 15):
        assert quadform.ruben_cdf(q, lam) == pytest.approx(chi2.cdf(q / 0.37, r), abs=1e-9)


@pytest.mark.parametrize("name", ["depth0", "deep_r6", "deep_r4", "deep_r2"])
def test_ruben_and_imhof_agree(name):
    """Where both are feasible they must agree; they share no code or method."""
    lam = SPECTRA[name]
    for frac in (0.2, 0.6, 1.0, 1.8, 3.0):
        q = frac * lam.sum()
        assert quadform.ruben_cdf(q, lam) == pytest.approx(quadform.imhof_cdf(q, lam), abs=1e-7)


@pytest.mark.parametrize("name", sorted(SPECTRA))
def test_imhof_is_a_cdf(name):
    """Monotone in q, and pinned to 0 and 1 at the ends.

    Checked over the range the search actually produces, ``q <= 4 tr(Sigma)``,
    where the loose tolerances agree with a tight integration to better than
    1e-7. Far beyond that the survival function drops below the quadrature
    tolerance and the value wobbles at the 1e-5 level, but by then it has been
    1.0 to eight decimals for a long while, so nothing downstream can see it.
    """
    lam = SPECTRA[name]
    tr = lam.sum()
    qs = np.concatenate([[0.0], np.geomspace(1e-3 * tr, 4 * tr, 40)])
    vals = np.array([quadform.imhof_cdf(q, lam) for q in qs])
    tight = np.array([quadform.imhof_cdf(q, lam, limit=400, epsabs=1e-12, epsrel=1e-11) for q in qs])
    assert vals[0] == 0.0
    assert np.abs(vals - tight).max() < 1e-7
    assert np.all(np.diff(vals) >= -1e-7)


@pytest.mark.parametrize("name", sorted(SPECTRA))
def test_loose_tolerances_do_not_cost_accuracy(name):
    """The shipped tolerances agree with a much tighter integration."""
    lam = SPECTRA[name]
    for frac in (0.05, 0.5, 1.2, 2.5):
        q = frac * lam.sum()
        tight = quadform.imhof_cdf(q, lam, limit=400, epsabs=1e-12, epsrel=1e-11)
        assert quadform.imhof_cdf(q, lam) == pytest.approx(tight, abs=1e-8)


def test_auto_dispatch_matches_whichever_method_is_feasible():
    for name, lam in SPECTRA.items():
        q = 1.2 * lam.sum()
        auto = quadform.exact_cdf(q, lam, "auto")
        assert auto == pytest.approx(quadform.imhof_cdf(q, lam), abs=1e-7)


def test_ruben_series_truncation_error_is_reported():
    """Ruben's coefficients are non-negative and sum to 1, so 1 - sum is a real bound."""
    s = quadform.RubenSeries(np.linspace(1.0, 1.5, 20))
    assert s.truncation_error < 1e-10
    assert s.a.min() >= 0.0
    # An ill-conditioned spectrum is refused by the dispatcher rather than mis-answered.
    assert quadform.ruben_terms(SPECTRA["mid_illcond"]) > quadform.RUBEN_TERM_BUDGET


# --------------------------------------------------------------------------- #
# ExactChi
# --------------------------------------------------------------------------- #
def test_lowrank_spectrum_identity(shared_TS):
    """A Sigma A^T and C^T A^T A C share their non-zero spectrum."""
    T, Sigma = shared_TS
    dx = dy = 8
    for dim, i, j in [(1, 1, 3), (1, 2, 5), (2, 1, 2), (2, 4, 7)]:
        i1 = merging.get_index(dim, i, dx, dy)
        i2 = merging.get_index(dim, j, dx, dy)
        dense = quadform.positive_eigenvalues(merging.update_Sigma(Sigma, i1, i2))
        crit = ExactChi(mode="lowrank")
        crit.init(T, Sigma)
        low = crit._spectrum(Sigma, i1, i2)
        assert dense.size == low.size
        assert np.allclose(np.sort(dense), np.sort(low), atol=1e-12)


@pytest.mark.parametrize("structure", ["ordinal", "tree", "greedy"])
def test_exact_modes_agree(shared_TS, structure):
    """dense / lowrank / cached are three ways to pay for one number."""
    T, Sigma = shared_TS
    dx = dy = 8
    sts = {"ordinal": (Ordinal(), Ordinal()),
           "tree": (Tree.binary(dx), Tree.binary(dy)),
           "greedy": (Saturated(), Saturated())}[structure]
    runs = {m: greedy_search(T, Sigma, dx, dy, *sts, ExactChi(mode=m))
            for m in ("dense", "lowrank", "cached")}
    ref = runs["dense"]
    for mode, r in runs.items():
        assert np.allclose(r.values, ref.values, atol=1e-10), mode
        assert r.partitions == ref.partitions, mode


def test_cache_is_actually_used(shared_TS):
    """Bootstrap draws revisit partitions; the cached mode must notice."""
    T, Sigma = shared_TS
    dx = dy = 8
    sts = (Tree.binary(dx), Tree.binary(dy))
    crit = ExactChi(mode="cached")
    rng = np.random.default_rng(0)
    for _ in range(4):
        greedy_search(rng.standard_normal(dx * dy), Sigma, dx, dy, *sts, crit)
    assert crit.stats["cache_hits"] > 0


def test_cache_is_reset_by_a_new_sigma(shared_TS):
    """Spectra depend on Sigma too: reusing an instance on a new Sigma must not go stale."""
    T, Sigma = shared_TS
    dx = dy = 8
    sts = (Tree.binary(dx), Tree.binary(dy))
    rng = np.random.default_rng(1)
    B = rng.standard_normal((dx * dy, dx * dy))
    Sigma2 = Sigma + 0.1 * (B @ B.T) / (dx * dy)
    shared = ExactChi(mode="cached")
    greedy_search(T, Sigma, dx, dy, *sts, shared)
    got = greedy_search(T, Sigma2, dx, dy, *sts, shared)
    ref = greedy_search(T, Sigma2, dx, dy, *sts, ExactChi(mode="dense"))
    assert np.allclose(got.values, ref.values, atol=1e-10)
