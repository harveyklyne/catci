"""``experiments/search_study.KronSearch`` must reproduce :mod:`catci.search`.

The study's fast path exploits ``Sigma = C_Y (x) C_X``; if it disagreed with the
package search anywhere, the study would be measuring a different method.
"""

import numpy as np
import pytest

from catci.search import beam_search, evaluate_path, random_merges
from catci.structure import Ordinal, Saturated, Tree
from search_study import KronSearch, direction


def _C(p):
    return np.diag(p) - np.outer(p, p)


@pytest.mark.parametrize("name", ["saturated", "tree", "ordinal"])
@pytest.mark.parametrize("width", [1, 4])
def test_kron_search_matches_package(name, width):
    rng = np.random.default_rng(7)
    dx, dy = 8, 4
    st = {"saturated": (Saturated(), Saturated()), "ordinal": (Ordinal(), Ordinal()),
          "tree": (Tree.binary(dx), Tree.binary(dy))}[name]
    Cx, Cy = _C(rng.dirichlet(np.ones(dx) * 3)), _C(rng.dirichlet(np.ones(dy) * 3))
    Sigma = np.kron(Cy, Cx)
    ks = KronSearch(Cx, Cy, *st)
    vals, vecs = np.linalg.eigh(Sigma)
    root = vecs[:, vals > 1e-12] * np.sqrt(vals[vals > 1e-12])
    for _ in range(5):
        T = root @ rng.standard_normal(root.shape[1]) * 1.5
        pkg = beam_search(T, Sigma, dx, dy, *st, width=width)
        v, merges = ks.beam(T, width)
        if width == 1:
            assert merges == pkg.merges
        # For width > 1 the same pair is often reachable from two beam states with
        # values equal to the last ulp, so which merge is *reported* can differ with
        # summation order; the values cannot.
        np.testing.assert_allclose(v, pkg.values, rtol=1e-9, atol=1e-12)

        path = random_merges(dx, dy, *st, rng)
        np.testing.assert_allclose(
            ks.evaluate(T, path), evaluate_path(T, Sigma, dx, dy, path).values,
            rtol=1e-9, atol=1e-12,
        )


@pytest.mark.parametrize("name", ["greedy_trap", "binary_tree", "step"])
def test_directions_in_range_of_sigma(name):
    d = 8
    p = np.full(d, 1 / d)
    C = _C(p)
    u = direction(name, d)
    Sigma = np.kron(C, C)
    assert np.allclose(Sigma @ np.linalg.pinv(Sigma) @ u, u)
