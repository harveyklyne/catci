"""Property / statistical tests that need no R oracle.

1. Rank-one update formulae must equal dense recomputation for random inputs
   (Hypothesis) -- the place a silent numerical error would be invisible in the
   figures (CODE_REVIEW.md 5.1).
2. Calibration under a known Gaussian: with T ~ N(0, Sigma), the adaptive
   double-bootstrap p-value is (approximately) uniform (CODE_REVIEW.md 5.2).
   This exercises the whole calibration path with no ML backend or DGP.
"""

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from catci import merging
from catci.calibrate import adaptive_pvalue
from catci.structure import Ordinal, Tree


# --------------------------------------------------------------------------- #
# 1. fast update == dense recompute
# --------------------------------------------------------------------------- #
@settings(max_examples=200, deadline=None)
@given(
    dx=st.integers(min_value=2, max_value=6),
    dy=st.integers(min_value=2, max_value=6),
    seed=st.integers(min_value=0, max_value=2 ** 32 - 1),
)
def test_updates_match_dense(dx, dy, seed):
    rng = np.random.default_rng(seed)
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    T = rng.standard_normal(p)

    # merge a random pair in a random dimension
    dimension = int(rng.integers(1, 3))
    d = dx if dimension == 1 else dy
    j1, j2 = sorted(rng.choice(np.arange(1, d + 1), size=2, replace=False))
    i1 = merging.get_index(dimension, int(j1), dx, dy)
    i2 = merging.get_index(dimension, int(j2), dx, dy)

    new_T = merging.update_T(T, i1, i2)
    new_Sigma = merging.update_Sigma(Sigma, i1, i2)

    fast_normsq = merging.update_normsq(float(np.sum(T ** 2)), T, i1, i2)
    fast_tr = merging.update_tr(float(np.trace(Sigma)), Sigma, i1, i2)
    fast_tr2 = merging.update_tr2(float(np.sum(Sigma ** 2)), Sigma, i1, i2)

    assert np.isclose(fast_normsq, np.sum(new_T ** 2))
    assert np.isclose(fast_tr, np.trace(new_Sigma))
    assert np.isclose(fast_tr2, np.sum(new_Sigma ** 2))


# --------------------------------------------------------------------------- #
# 2. calibration under a known Gaussian
# --------------------------------------------------------------------------- #
def test_calibration_under_gaussian():
    rng = np.random.default_rng(0)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    sqrtS = np.linalg.cholesky(Sigma)

    reps = 400
    pvals = np.empty(reps)
    for r in range(reps):
        T = sqrtS @ rng.standard_normal(p)
        pvals[r] = adaptive_pvalue(
            T, Sigma, dx, dy, Ordinal(), Ordinal(), n_boot=100, rng=rng
        )

    # Uniform(0,1): mean ~ 0.5, and the level-alpha rejection rate ~ alpha.
    assert abs(pvals.mean() - 0.5) < 0.06
    for alpha in (0.05, 0.10, 0.20):
        rate = float(np.mean(pvals < alpha))
        assert abs(rate - alpha) < 0.05, f"rejection rate {rate:.3f} at alpha={alpha}"


def test_tree_calibrates_too():
    # tree search is also a valid test -> also calibrates (a calibration test
    # alone would NOT catch finding #1; that is why we also pin the path).
    rng = np.random.default_rng(1)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    sqrtS = np.linalg.cholesky(Sigma)

    reps = 300
    pvals = np.array([
        adaptive_pvalue(
            sqrtS @ rng.standard_normal(p), Sigma, dx, dy,
            Tree.binary(dx), Tree.binary(dy), n_boot=100, rng=rng,
        )
        for _ in range(reps)
    ])
    assert abs(pvals.mean() - 0.5) < 0.07
