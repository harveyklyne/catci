"""A label absent from the sample with propensity exactly 0 must not break the searches.

Its residual column is identically zero, so ``Sigma`` has zero-variance coordinates,
and a partition such as ``{absent} vs {rest}`` has an identically-zero statistic
(the rest is minus the absent column, since residual rows sum to 0). Found on the
semi-synthetic adult data, where an MLP fitted without the class gives it
probability 0; the merge search then scored that partition ``0/0 = nan`` for some
draws only and raised "Search paths have different lengths across draws".
"""

from __future__ import annotations

import numpy as np

from catci import statistic
from catci.bootstrap import bootstrap_T
from catci.gcm import form_t_sigma
from catci.search import MergeSearch, SplitSearch
from catci.structure import Ordinal


def test_degenerate_partition_scores_zero_and_others_are_unchanged():
    normsq = np.array([3.0, -1e-16, 2.0])
    tr = np.array([2.0, -1.6e-16, 1e-13])
    tr2 = np.array([1.5, 2.2e-18, 1e-26])
    out = statistic.approx_chi_array(normsq, tr, tr2)
    assert out[1] == 0.0 and out[2] == 0.0
    from scipy.special import gammainc
    assert out[0] == gammainc((2.0 ** 2 / 1.5) / 2, 3.0 / (1.5 / 2.0) / 2)
    assert statistic.approx_chi_statistic(-1e-16, -1e-16, 1e-18) == 0.0


def _absent_level_data(n=600, dx=4, dy=6, seed=0):
    """Y never takes level 1, and its fitted propensity is exactly 0."""
    rng = np.random.default_rng(seed)
    x = rng.integers(1, dx + 1, n)
    y = rng.integers(2, dy + 1, n)
    f = np.full((n, dx), 1.0 / dx)
    g = np.zeros((n, dy))
    g[:, 1:] = 1.0 / (dy - 1)
    return x, y, f, g


def test_searches_survive_zero_variance_coordinates():
    dx, dy = 4, 6
    x, y, f, g = _absent_level_data(dx=dx, dy=dy)
    ts = form_t_sigma(x, y, f, g, normalise=False)
    assert (np.diag(ts.Sigma) == 0).sum() == dx
    T_all = np.column_stack([ts.T_vector, bootstrap_T(ts.Sigma, 3000, np.random.default_rng(1))])
    for search in (MergeSearch(), SplitSearch(), SplitSearch(2)):
        paths = search.paths(T_all, ts.Sigma, dx, dy, Ordinal(), Ordinal())
        assert np.isfinite(paths).all()
        assert ((paths >= 0) & (paths <= 1)).all()
