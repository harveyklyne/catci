"""DGP tests: invariants, the review fixes 2a-2c, and binary_tree == R at d=8."""

import numpy as np
import pytest

import dgp


# The hand-coded R binary-tree interaction (column-major, nrow=8) from
# R/simulation_settings.R -- the recursive construction must reproduce it.
_R_BINARY_TREE_8 = np.array([
    -1.1, -1.1, -0.9, -0.9, 0.9, 0.9, 1.1, 1.1,
    -1.1, -1.1, -0.9, -0.9, 0.9, 0.9, 1.1, 1.1,
    -0.9, -0.9, -1.1, -1.1, 1.1, 1.1, 0.9, 0.9,
    -0.9, -0.9, -1.1, -1.1, 1.1, 1.1, 0.9, 0.9,
    0.9, 0.9, 1.1, 1.1, -1.1, -1.1, -0.9, -0.9,
    0.9, 0.9, 1.1, 1.1, -1.1, -1.1, -0.9, -0.9,
    1.1, 1.1, 0.9, 0.9, -0.9, -0.9, -1.1, -1.1,
    1.1, 1.1, 0.9, 0.9, -0.9, -0.9, -1.1, -1.1,
]).reshape(8, 8, order="F")


# --------------------------------------------------------------------------- #
# marginals
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("pdf", [dgp.lin_pdf, dgp.vee_pdf, dgp.hat_pdf])
@pytest.mark.parametrize("d", [2, 3, 4, 5, 8, 16])
def test_marginal_pmfs_sum_to_one(pdf, d):
    p = pdf(d)
    assert p.shape == (d,)
    assert np.isclose(p.sum(), 1.0)
    assert np.all(p > 0)


@pytest.mark.parametrize("setting", ["sin", "sig"])
@pytest.mark.parametrize("d", [4, 8, 16])
def test_zconditional_pmfs_sum_to_one(setting, d):
    z = np.linspace(-3, 3, 50)
    P = dgp.get_pdf(setting, d, z=z)
    assert P.shape == (50, d)
    np.testing.assert_allclose(P.sum(axis=1), 1.0, atol=1e-12)
    assert np.all(P >= 0)


def test_vee_is_symmetric():
    p = dgp.vee_pdf(8)
    np.testing.assert_allclose(p, p[::-1])


# --------------------------------------------------------------------------- #
# interactions: zero margins + binary_tree matches R at d=8 (fix 2c)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("setting,dx,dy", [
    ("step", 8, 8), ("step", 6, 4), ("alt", 8, 8), ("alt", 5, 5),
    ("binary_tree", 4, 4), ("binary_tree", 8, 8), ("binary_tree", 16, 16),
])
def test_interaction_zero_margins(setting, dx, dy):
    M = dgp.get_int(setting, dx, dy)
    assert M.shape == (dx, dy)
    np.testing.assert_allclose(M.sum(axis=0), 0.0, atol=1e-12)
    np.testing.assert_allclose(M.sum(axis=1), 0.0, atol=1e-12)


def test_binary_tree_matches_r_literal_at_d8():
    np.testing.assert_allclose(dgp.binary_tree_interaction(8), _R_BINARY_TREE_8)


def test_binary_tree_rejects_non_power_of_two():
    with pytest.raises(ValueError):
        dgp.binary_tree_interaction(6)


# --------------------------------------------------------------------------- #
# fix 2b: cross-family pairs no longer crash
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("xs,ys", [
    ("lin", "sin"), ("sin", "lin"), ("hat", "sig"), ("sig", "vee"), ("sin", "sig"),
])
def test_cross_family_pairs_run(xs, ys):
    rng = np.random.default_rng(0)
    data = dgp.simulate_data(500, 8, 8, xs, ys, strength=0.5, intsetting="step", rng=rng)
    assert data["x"].min() >= 1 and data["x"].max() <= 8
    assert data["f"].shape == (500, 8) and data["g"].shape == (500, 8)
    np.testing.assert_allclose(data["f"].sum(axis=1), 1.0, atol=1e-10)


# --------------------------------------------------------------------------- #
# fix 2a: single sampling path; marginals recovered, dependence increases with strength
# --------------------------------------------------------------------------- #
def test_marginals_recovered_and_dependence_grows():
    rng = np.random.default_rng(1)
    n, d = 4000, 8

    # At strength 0 the empirical joint should factorise (mutual information ~ 0);
    # at high strength it should not. Use a coarse 2x2 MI as a monotonicity check.
    def binary_mi(x, y, dx, dy):
        # collapse to halves and measure dependence of the halves
        xb = (x > dx // 2).astype(int)
        yb = (y > dy // 2).astype(int)
        joint = np.histogram2d(xb, yb, bins=[2, 2])[0] / len(x)
        px = joint.sum(1, keepdims=True); py = joint.sum(0, keepdims=True)
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = joint * np.log(joint / (px * py))
        return np.nansum(terms)

    d0 = dgp.simulate_data(n, d, d, "lin", "lin", 0.0, "binary_tree", rng=rng)
    d1 = dgp.simulate_data(n, d, d, "lin", "lin", 1.5, "binary_tree", rng=rng)
    mi0 = binary_mi(d0["x"], d0["y"], d, d)
    mi1 = binary_mi(d1["x"], d1["y"], d, d)
    assert mi1 > mi0
    # f is the true propensity: its column means approximate the empirical label freqs
    freq = np.bincount(d0["x"], minlength=d + 1)[1:] / n
    np.testing.assert_allclose(d0["f"].mean(axis=0), freq, atol=0.03)


def test_permute_breaks_label_order():
    rng = np.random.default_rng(2)
    base = dgp.simulate_data(300, 8, 8, "lin", "lin", 0.5, "step", permute=False, rng=np.random.default_rng(2))
    perm = dgp.simulate_data(300, 8, 8, "lin", "lin", 0.5, "step", permute=True, rng=np.random.default_rng(2))
    # same rng seed up to the extra permutation draw -> f columns are a permutation
    assert perm["f"].shape == base["f"].shape
    np.testing.assert_allclose(np.sort(perm["f"], axis=1), np.sort(base["f"], axis=1), atol=1e-12)
