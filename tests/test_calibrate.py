"""Direct tests of the minP calibration -- no DGP, no search, no ML backend.

The calibration is exact whenever the observed path and the bootstrap paths are
exchangeable, so it can be tested by handing it exchangeable paths directly. That
is the cheapest possible setting and the one that isolates the calibration from
everything upstream of it, which is why the previous suite -- which only reached
``double_bootstrap_pvalue`` through ``adaptive_pvalue`` -- missed a size defect
that grew with the number of search levels (0.049 at L=1 up to 0.358 at L=40).

``L`` is parametrised here for exactly that reason.
"""

import numpy as np
import pytest

from catci.calibrate import bonferroni_pvalue, double_bootstrap_pvalue

B = 100


def _exchangeable_paths(L, reps, rng):
    """``reps`` draws of ``(L, B+1)`` paths, correlated across levels like a real search."""
    A = rng.standard_normal((L, L))
    return (A @ rng.standard_normal((reps, L, B + 1))).transpose(1, 0, 2)


def _pvalues(L, reps, rng, fn=double_bootstrap_pvalue):
    paths = _exchangeable_paths(L, reps, rng)  # (L, reps, B+1)
    return np.array([fn(paths[:, r, 0], paths[:, r, 1:], rng) for r in range(reps)])


# --------------------------------------------------------------------------- #
# The main guarantee: uniform under exchangeability, at every L
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("L", [1, 2, 5, 10, 20, 40])
def test_uniform_under_exchangeability(L):
    reps = 10000
    pvals = _pvalues(L, reps, np.random.default_rng(1000 + L))

    # Binomial SE is 0.0022 at alpha=0.05; the observed deviations are all under
    # 0.0035, so 0.010 leaves ~3x headroom on fixed seeds while still failing the
    # pre-minP calibration at every L >= 2 (which gave 0.064, 0.087, ..., 0.336).
    for alpha in (0.05, 0.10):
        rate = float(np.mean(pvals <= alpha))
        assert abs(rate - alpha) < 0.010, f"rejection rate {rate:.4f} at alpha={alpha}, L={L}"
    assert abs(pvals.mean() - 0.5) < 0.02


def test_size_does_not_grow_with_L():
    """The signature of the old defect: size climbing monotonically in L."""
    rates = [
        float(np.mean(_pvalues(L, 4000, np.random.default_rng(7)) <= 0.05))
        for L in (1, 10, 40)
    ]
    assert max(rates) - min(rates) < 0.02, f"size varies with L: {rates}"


# --------------------------------------------------------------------------- #
# Structural properties
# --------------------------------------------------------------------------- #
def test_collapses_to_single_bootstrap_at_L1():
    """With one statistic the two stages cancel: the ordinary bootstrap p-value."""
    rng = np.random.default_rng(0)
    for _ in range(500):
        obs = rng.standard_normal()
        boot = rng.standard_normal(B)  # distinct with probability one
        expected = (1 + np.sum(boot >= obs)) / (B + 1)
        got = double_bootstrap_pvalue(np.array([obs]), boot[None, :], rng)
        assert got == pytest.approx(expected)


def test_pvalue_on_the_grid():
    rng = np.random.default_rng(2)
    grid = np.arange(1, B + 2) / (B + 1)
    for L in (1, 5, 13):
        paths = _exchangeable_paths(L, 50, rng)
        for r in range(50):
            p = double_bootstrap_pvalue(paths[:, r, 0], paths[:, r, 1:], rng)
            assert np.isclose(grid, p).any(), f"{p} is not on the 1/(B+1) grid"


def test_monotone_in_observed():
    """Raising the whole observed path cannot raise the p-value."""
    rng = np.random.default_rng(3)
    paths = _exchangeable_paths(6, 1, rng)[:, 0, :]
    obs, boot = paths[:, 0], paths[:, 1:]
    pvals = [
        double_bootstrap_pvalue(obs + shift, boot, np.random.default_rng(11))
        for shift in (0.0, 1.0, 2.0, 5.0)
    ]
    assert pvals == sorted(pvals, reverse=True), pvals


def test_rejects_a_clear_alternative():
    rng = np.random.default_rng(4)
    paths = _exchangeable_paths(5, 200, rng)
    pvals = np.array([
        double_bootstrap_pvalue(paths[:, r, 0] + 8.0, paths[:, r, 1:], rng)
        for r in range(200)
    ])
    assert float(np.mean(pvals <= 0.05)) > 0.95


def test_rejects_shape_mismatch():
    rng = np.random.default_rng(5)
    with pytest.raises(ValueError):
        double_bootstrap_pvalue(np.zeros(3), np.zeros((4, B)), rng)


# --------------------------------------------------------------------------- #
# The simple-FWER comparator
# --------------------------------------------------------------------------- #
def test_bonferroni_is_valid_but_floored():
    L, reps = 5, 4000
    rng = np.random.default_rng(6)
    pvals = _pvalues(L, reps, rng, fn=bonferroni_pvalue)
    assert pvals.min() >= L / (B + 1) - 1e-12  # resolution floor
    assert float(np.mean(pvals <= 0.05)) <= 0.05  # conservative, never inflated


def test_bonferroni_cannot_reject_when_L_exceeds_alpha_budget():
    """Floor L/(B+1) > alpha, so no data whatsoever can make it reject."""
    L = 13  # dx = dy = 8
    rng = np.random.default_rng(8)
    paths = _exchangeable_paths(L, 100, rng)
    pvals = np.array([
        bonferroni_pvalue(paths[:, r, 0] + 10.0, paths[:, r, 1:], rng) for r in range(100)
    ])
    assert L / (B + 1) > 0.05
    assert float(np.mean(pvals <= 0.05)) == 0.0


def test_bonferroni_is_never_more_significant_than_minp():
    rng = np.random.default_rng(9)
    paths = _exchangeable_paths(5, 200, rng)
    for r in range(200):
        obs, boot = paths[:, r, 0], paths[:, r, 1:]
        seed = np.random.default_rng(r)
        minp = double_bootstrap_pvalue(obs, boot, np.random.default_rng(r))
        bonf = bonferroni_pvalue(obs, boot, seed)
        assert bonf >= minp - 1e-12
