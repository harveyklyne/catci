"""Pins the exact identities the semi-synthetic adult construction rests on.

Every claim in the :mod:`adult_semisynth` docstring is an identity on a finite
population, so each is checked to machine precision rather than statistically:

* ``lam = 0`` is an exact null -- the kernel does not depend on ``x``;
* ``g = P(y | w)`` for every ``lam`` (C2), and the kernel is a pmf (C1);
* the departure from independence is exactly linear in ``lam``;
* ``lam = 1`` with the real direction reproduces the population's joint law;
* the planted direction is admissible, has ``max_lambda >= 1`` and is invisible
  to Ankan & Textor's Q1 stratum by stratum (C3);
* :func:`true_propensities` is consistent with the kernel for a rich ``Z``.
"""

from __future__ import annotations

import numpy as np
import pytest

import adult
import adult_semisynth as ss

TOL = 1e-12


@pytest.fixture(scope="module")
def data():
    return adult.load()


@pytest.fixture(scope="module")
def pop(data):
    """Education x Income given {Age, Sex}: the pair the four-way comparison uses."""
    return ss.build_population(data, "Education", "Income")


@pytest.fixture(scope="module")
def pop_occ(data):
    """Occupation x Income, pooled at n = 1000 (dx >= 3 on the categorical side)."""
    pooled = ss.pool_rare(data, "Occupation", min_expected=5, n=1000)
    return ss.build_population(pooled, "Occupation", "Income")


def _reachable(pop):
    return pop.px > 0


@pytest.mark.parametrize("which", ["pop", "pop_occ"])
def test_population_tables_are_consistent(which, request):
    pop = request.getfixturevalue(which)
    np.testing.assert_allclose(pop.joint.sum(axis=(1, 2)), 1.0, atol=TOL)
    np.testing.assert_allclose(pop.pw.sum(), 1.0, atol=TOL)
    np.testing.assert_allclose(pop.px, pop.joint.sum(2), atol=TOL)
    np.testing.assert_allclose(pop.py, pop.joint.sum(1), atol=TOL)


def test_null_is_exact(pop):
    """At ``lam = 0`` every reachable row of the kernel is ``P(y | w)``, whatever ``x``."""
    k0 = pop.kernel(0.0)
    target = np.broadcast_to(pop.py[:, None, :], k0.shape)
    np.testing.assert_allclose(k0, target, atol=TOL)


@pytest.mark.parametrize("lam", [0.0, 0.3, 1.0])
def test_real_kernel_is_pmf_and_holds_g_fixed(pop, lam):
    k = pop.kernel(lam)
    reach = _reachable(pop)
    np.testing.assert_allclose(k.sum(2)[reach], 1.0, atol=TOL)  # C1
    assert k[reach].min() >= -TOL
    np.testing.assert_allclose(np.einsum("cx,cxy->cy", pop.px, k), pop.py, atol=TOL)  # C2


def test_lam_one_reproduces_the_population(pop):
    np.testing.assert_allclose(pop.px[:, :, None] * pop.kernel(1.0), pop.joint, atol=TOL)


@pytest.mark.parametrize("lam", [0.25, 0.5, 0.8])
def test_departure_is_linear_in_lam(pop, lam):
    indep = pop.px[:, :, None] * pop.py[:, None, :]
    departure = pop.px[:, :, None] * pop.kernel(lam) - indep
    np.testing.assert_allclose(departure, lam * (pop.joint - indep), atol=TOL)


def test_ncp_scales_with_direction_norm(pop):
    """``ncp_per_n`` is the lam = 1 value; halving the direction quarters it."""
    delta = pop.delta_real()
    full = pop.ncp_per_n(delta)
    assert full > 0
    np.testing.assert_allclose(pop.ncp_per_n(0.5 * delta), full / 4, rtol=1e-10)


@pytest.mark.parametrize("which", ["pop", "pop_occ"])
def test_planted_direction_is_admissible_and_q1_blind(which, request):
    pop = request.getfixturevalue(which)
    delta = ss.planted_direction(pop)
    check = ss.is_valid_direction(pop, delta)
    assert check["ok"], check
    assert check["max_lam"] >= 1.0 - 1e-12
    assert abs(ss.q1_signal(pop, delta)) < 1e-14
    # C3 holds stratum by stratum, not just on average.
    rx = ss._ls_scores(pop.px)
    ry = ss._ls_scores(pop.py)
    per_cell = np.einsum("cx,cx,cy,cxy->c", pop.px, rx, ry, delta)
    np.testing.assert_allclose(per_cell, 0.0, atol=1e-14)
    # ...and it is a real alternative, not the zero direction.
    assert pop.ncp_per_n(delta) > 1e-4


def test_project_out_visible_is_admissible_and_q1_blind(pop):
    delta = ss.project_out_visible(pop)
    assert ss.is_valid_direction(pop, delta)["ok"]
    assert abs(ss.q1_signal(pop, delta)) < 1e-14
    assert ss.visible_fraction(pop, delta) < 1e-20


def test_planted_needs_three_levels(data):
    pop = ss.build_population(data, "Sex", "Income")
    with pytest.raises(ValueError, match="dx >= 3"):
        ss.planted_direction(pop)


@pytest.mark.parametrize("lam", [0.0, 0.7])
def test_true_propensities_at_w_reduce_to_population_margins(pop, lam):
    rows = np.arange(0, pop.N, 97)
    f, g = ss.true_propensities(pop, rows, lam, pop.w_names)
    np.testing.assert_allclose(f, pop.px[pop.cell[rows]], atol=TOL)
    np.testing.assert_allclose(g, pop.py[pop.cell[rows]], atol=TOL)


@pytest.mark.parametrize("delta_kind", ["real", "planted"])
@pytest.mark.parametrize("lam", [0.0, 0.6])
def test_true_propensities_rich_z_average_back_to_w(pop, lam, delta_kind):
    """With ``Z`` finer than ``W``, ``g(y | z)`` still averages to ``P(y | w)`` within ``w``."""
    delta = pop.delta_real() if delta_kind == "real" else ss.planted_direction(pop)
    z_names = pop.w_names + ("Race",)
    rows = np.arange(pop.N)
    f, g = ss.true_propensities(pop, rows, lam, z_names, delta)
    np.testing.assert_allclose(f.sum(1), 1.0, atol=TOL)
    np.testing.assert_allclose(g.sum(1), 1.0, atol=TOL)
    n_cells = len(pop.joint)
    g_by_w = np.zeros((n_cells, pop.dy))
    np.add.at(g_by_w, pop.cell, g)
    g_by_w /= np.bincount(pop.cell, minlength=n_cells)[:, None]
    np.testing.assert_allclose(g_by_w, pop.py, atol=1e-10)


def test_draw_shapes_and_label_ranges(pop):
    rng = np.random.default_rng(0)
    rep = ss.draw(pop, 500, 0.5, rng, z_names=pop.w_names + ("Race",))
    assert rep.x.shape == rep.y.shape == (500,)
    assert rep.z.shape == (500, 3)
    assert 1 <= rep.x.min() and rep.x.max() <= pop.dx
    assert 1 <= rep.y.min() and rep.y.max() <= pop.dy
    assert rep.f_true.shape == (500, pop.dx) and rep.g_true.shape == (500, pop.dy)


def test_draw_rejects_z_missing_the_generating_stratum(pop):
    with pytest.raises(ValueError, match="generating stratum"):
        ss.draw(pop, 100, 0.0, np.random.default_rng(0), z_names=("Age",))


def test_draw_rejects_lam_beyond_the_kernel(pop):
    delta = ss.planted_direction(pop)
    lam = 1.01 * ss.max_lambda(pop, delta)
    with pytest.raises(ValueError, match="negative"):
        ss.draw(pop, pop.N, lam, np.random.default_rng(0), delta=delta)


def test_pool_rare_merges_armed_forces_into_one_last_level(data):
    pooled = ss.pool_rare(data, "Occupation", min_expected=5, n=1000)
    levels = pooled.variables["Occupation"].levels
    assert levels[-1].startswith("Other (") and "Armed-Forces" in levels[-1]
    assert len(levels) < len(data.variables["Occupation"].levels)
    codes = pooled.codes["Occupation"]
    assert codes.min() == 1 and codes.max() == len(levels)
    # Every level now has an expected count of at least min_expected at n = 1000,
    # except the pooled one, which only has to be at least as big as its parts.
    share = codes.value_counts(normalize=True).sort_index().to_numpy()
    assert (share[:-1] * 1000 >= 5).all()
    # Pooling a variable leaves every other column untouched.
    others = [c for c in data.codes.columns if c != "Occupation"]
    assert pooled.codes[others].equals(data.codes[others])
