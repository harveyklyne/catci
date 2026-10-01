"""Tests for the Ankan & Textor (AAAI-23) port.

There is no reference implementation to differential-test against -- the paper
shipped no code, and the authors' later packages implement a different statistic
(see ``experiments/ankan_textor.py``). So the port is pinned by the properties
the paper itself asserts:

* the LS residual collapses to observed-minus-expected in the binary case
  (stated in the paper);
* Q1 is the squared GCM of Shah & Peters (2020) (stated in the paper);
* Q1/Q2/Q3 are one quadratic form at three shapes, with the df of
  Propositions 1-3;
* dropping the last dummy is what makes ``Sigma_d`` invertible (stated);
* and the whole thing is calibrated under a null it should be calibrated under.

The calibration test is the one that would actually catch a wrong statistic.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import chi2, kstest

from ankan_textor import (
    GLM,
    RFT,
    design_matrix,
    ci_test,
    ls_residuals,
    q_statistic,
    residual_matrix,
)


# --------------------------------------------------------------------------- #
# Residual
# --------------------------------------------------------------------------- #
def test_ls_residual_binary_is_observed_minus_expected():
    """Paper: for binary Y, the LS residual reduces to ``y - p_hat(Y = 1)``."""
    rng = np.random.default_rng(0)
    n = 200
    p2 = rng.uniform(0.05, 0.95, n)
    proba = np.column_stack([1 - p2, p2])
    labels = rng.integers(1, 3, n)

    got = ls_residuals(labels, proba)
    expected = (labels == 2).astype(float) - p2
    np.testing.assert_allclose(got, expected, atol=1e-12)


def test_ls_residual_is_between_minus_one_and_one():
    rng = np.random.default_rng(1)
    n, k = 300, 5
    proba = rng.dirichlet(np.ones(k), n)
    labels = rng.integers(1, k + 1, n)
    r = ls_residuals(labels, proba)
    assert np.all(r >= -1.0) and np.all(r <= 1.0)


def test_ls_residual_respects_level_order():
    """A higher observed level must give a weakly larger residual, same propensities."""
    proba = np.tile(np.array([0.2, 0.3, 0.5]), (3, 1))
    r = ls_residuals(np.array([1, 2, 3]), proba)
    assert r[0] < r[1] < r[2]


def test_binary_categorical_and_ordinal_residuals_agree_up_to_sign():
    """A 2-level variable is both; the two paths must give the same information."""
    rng = np.random.default_rng(2)
    n = 100
    p2 = rng.uniform(0.1, 0.9, n)
    proba = np.column_stack([1 - p2, p2])
    labels = rng.integers(1, 3, n)

    ordinal = residual_matrix(labels, proba, "ordinal")
    categorical = residual_matrix(labels, proba, "categorical")
    assert ordinal.shape == categorical.shape == (n, 1)
    np.testing.assert_allclose(ordinal, -categorical, atol=1e-12)


# --------------------------------------------------------------------------- #
# Statistic
# --------------------------------------------------------------------------- #
def test_q1_matches_squared_gcm_closed_form():
    """Paper's Q1 = (1/n)(Rx.Ry)^2 / var(Rx Ry) -- the squared GCM."""
    rng = np.random.default_rng(3)
    n = 500
    rx = rng.standard_normal((n, 1))
    ry = rng.standard_normal((n, 1))

    q = q_statistic(rx, ry)
    prod = (rx[:, 0] * ry[:, 0])
    expected = (prod.sum() ** 2) / prod.var(ddof=1) / n

    assert q.df == 1
    assert q.statistic == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("a,b", [(1, 1), (4, 1), (1, 3), (5, 2), (12, 5)])
def test_degrees_of_freedom_are_product_of_residual_widths(a, b):
    """Propositions 1-3: df is 1, k-1, and (k-1)(r-1) -- i.e. a*b throughout."""
    rng = np.random.default_rng(4)
    n = 400
    q = q_statistic(rng.standard_normal((n, a)), rng.standard_normal((n, b)))
    assert q.df == a * b
    assert q.well_conditioned


def test_statistic_is_invariant_to_column_order():
    """The quadratic form must not depend on how d and Sigma_d are indexed."""
    rng = np.random.default_rng(5)
    n = 300
    rx = rng.standard_normal((n, 3))
    ry = rng.standard_normal((n, 2))
    base = q_statistic(rx, ry)
    permuted = q_statistic(rx[:, [2, 0, 1]], ry[:, [1, 0]])
    assert permuted.statistic == pytest.approx(base.statistic, rel=1e-8)


def test_keeping_every_dummy_makes_sigma_singular():
    """Paper: one dummy is dropped 'because otherwise Sigma_d would not be full rank'."""
    rng = np.random.default_rng(6)
    n, k = 200, 4
    labels = rng.integers(1, k + 1, n)
    proba = rng.dirichlet(np.ones(k), n)
    dummies = np.zeros((n, k))
    dummies[np.arange(n), labels - 1] = 1.0
    all_k = dummies - proba  # every indicator: rows sum to zero

    ry = rng.standard_normal((n, 1))
    m_full = (all_k[:, None, :] * ry[:, :, None]).reshape(n, -1)
    assert np.linalg.matrix_rank(np.cov(m_full, rowvar=False)) < k

    dropped = residual_matrix(labels, proba, "categorical")
    m_dropped = (dropped[:, None, :] * ry[:, :, None]).reshape(n, -1)
    assert np.linalg.matrix_rank(np.cov(m_dropped, rowvar=False)) == k - 1


# --------------------------------------------------------------------------- #
# Calibration -- the test that would catch a wrong statistic
# --------------------------------------------------------------------------- #
def _null_sample(rng, n, kx, ky, kz):
    """X and Y both depend on Z, and are independent given Z."""
    z = rng.integers(1, kz + 1, n)
    bx = rng.standard_normal((kz, kx))
    by = rng.standard_normal((kz, ky))

    def draw(beta, k):
        logits = beta[z - 1]
        p = np.exp(logits - logits.max(axis=1, keepdims=True))
        p /= p.sum(axis=1, keepdims=True)
        u = rng.random((n, 1))
        return (u < np.cumsum(p, axis=1)).argmax(axis=1) + 1

    return z, draw(bx, kx), draw(by, ky)


@pytest.mark.parametrize(
    "kind_x,kind_y,kx,ky,which",
    [
        ("ordinal", "ordinal", 3, 3, "Q1"),
        ("categorical", "ordinal", 4, 3, "Q2"),
        ("categorical", "categorical", 3, 4, "Q3"),
    ],
)
def test_calibrated_under_the_null(kind_x, kind_y, kx, ky, which):
    """Type-I error near nominal, and p-values not far from uniform, under H0.

    Reps are modest to keep the suite fast, so the assertion is loose: it is a
    guard against a statistic that is wrong (wrong df, wrong normalisation,
    un-inverted Sigma), not a precision calibration study.
    """
    rng = np.random.default_rng(7)
    kz, n, reps = 4, 600, 60
    model = GLM()

    pvals = []
    for _ in range(reps):
        z, x, y = _null_sample(rng, n, kx, ky, kz)
        design = design_matrix(z[:, None], [kz])
        res = ci_test(x, y, design, kind_x, kind_y, kx, ky, model)
        assert res.which == which
        assert res.df == (kx - 1 if kind_x == "categorical" else 1) * (
            ky - 1 if kind_y == "categorical" else 1
        )
        pvals.append(res.p_value)

    pvals = np.asarray(pvals)
    assert 0.0 <= pvals.min() and pvals.max() <= 1.0
    # rejection rate at 0.05 should be near 0.05; binomial SE at 60 reps is ~0.028
    assert (pvals < 0.05).mean() < 0.20
    # and the whole distribution should be plausibly uniform
    assert kstest(pvals, "uniform").pvalue > 0.001


def test_symmetry_in_x_and_y():
    """The paper lists the symmetry property as a desideratum of the approach."""
    rng = np.random.default_rng(8)
    kz, n = 3, 500
    z, x, y = _null_sample(rng, n, 4, 3, kz)
    design = design_matrix(z[:, None], [kz])
    model = GLM()

    forward = ci_test(x, y, design, "categorical", "categorical", 4, 3, model)
    reverse = ci_test(y, x, design, "categorical", "categorical", 3, 4, model)
    assert forward.df == reverse.df
    assert forward.statistic == pytest.approx(reverse.statistic, rel=1e-8)


def test_detects_a_real_dependence():
    """Sanity: an unconditionally strong X-Y link is rejected."""
    rng = np.random.default_rng(9)
    n, kz = 800, 3
    z = rng.integers(1, kz + 1, n)
    x = rng.integers(1, 4, n)
    flip = rng.random(n) < 0.25
    y = np.where(flip, rng.integers(1, 4, n), x)  # y tracks x closely

    design = design_matrix(z[:, None], [kz])
    res = ci_test(x, y, design, "categorical", "categorical", 3, 3, GLM())
    assert res.p_value < 1e-6


def test_rft_runs_and_agrees_in_shape_with_glm():
    rng = np.random.default_rng(10)
    kz, n = 3, 400
    z, x, y = _null_sample(rng, n, 3, 3, kz)
    design = design_matrix(z[:, None], [kz])

    glm = ci_test(x, y, design, "categorical", "ordinal", 3, 3, GLM())
    rft = ci_test(x, y, design, "categorical", "ordinal", 3, 3, RFT())
    assert glm.df == rft.df == 2
    assert 0.0 <= rft.p_value <= 1.0


# --------------------------------------------------------------------------- #
# Conditioning diagnostics
# --------------------------------------------------------------------------- #
def test_wellconditioned_flags_a_degenerate_product_column():
    """A product column that is identically zero must not pass as a valid test.

    This is the adult-income failure mode in miniature: an X level and a Y level
    that never co-occur give a zero column, so Sigma_d loses rank and the
    quadratic form is undefined -- but a plain solve may still return a number.
    """
    rng = np.random.default_rng(11)
    n = 400
    rx = rng.standard_normal((n, 2))
    ry = rng.standard_normal((n, 2))
    rx[:, 1] = np.where(np.arange(n) < n // 2, rx[:, 1], 0.0)
    ry[:, 1] = np.where(np.arange(n) < n // 2, 0.0, ry[:, 1])
    # rx[:, 1] * ry[:, 1] is now identically zero
    q = q_statistic(rx, ry, method="pinv")
    assert q.rank < q.n_columns
    assert not q.well_conditioned


def test_pinv_sets_df_to_the_rank():
    rng = np.random.default_rng(12)
    n = 300
    rx = rng.standard_normal((n, 3))
    ry = rng.standard_normal((n, 2))
    q = q_statistic(rx, ry, method="pinv")
    assert q.rank == q.n_columns == 6
    assert q.df == q.rank
    # full rank: pinv and solve must agree
    assert q.statistic == pytest.approx(q_statistic(rx, ry).statistic, rel=1e-8)


def test_unknown_method_rejected():
    rng = np.random.default_rng(13)
    with pytest.raises(ValueError, match="solve.*pinv|pinv.*solve"):
        q_statistic(rng.standard_normal((50, 1)), rng.standard_normal((50, 1)), method="ridge")
