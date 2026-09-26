"""The vectorised search must reproduce the per-candidate loop it replaced.

``_greedy_search_loop`` is the original implementation, itself pinned to the R
oracle by ``test_search.py``; here it is the oracle for the vectorised kernel on
inputs the fixture does not cover -- unequal ``dx, dy``, low-rank ``Sigma`` from
real residuals, n-ary trees (the ``permitted_merges`` fallback), and many draws
through one batched call.
"""

import numpy as np
import pytest

from catci.gcm import form_t_sigma
from catci.search import _greedy_search_loop, greedy_search, greedy_search_paths
from catci.structure import Ordinal, Saturated, Tree, TreeNode, _make_leaf, _make_parent

SHAPES = [(3, 3), (4, 7), (7, 4), (6, 6), (9, 3), (8, 8)]


def _structures(name, dx, dy):
    return {
        "ordinal": lambda: (Ordinal(), Ordinal()),
        "greedy": lambda: (Saturated(), Saturated()),
        "tree": lambda: (Tree.binary(dx), Tree.binary(dy)),
        "mixed": lambda: (Tree.binary(dx), Ordinal()),
    }[name]()


def _random_ts(dx, dy, rng, n=300):
    """A realistic (T, Sigma): residual products, so Sigma is a sample covariance."""
    x = rng.integers(1, dx + 1, n)
    y = rng.integers(1, dy + 1, n)
    f = rng.dirichlet(np.ones(dx), n)
    g = rng.dirichlet(np.ones(dy), n)
    ts = form_t_sigma(x, y, f, g, normalise=bool(rng.integers(2)))
    return ts.T_vector, ts.Sigma


def _assert_same(res, ref):
    np.testing.assert_allclose(res.values, ref.values, rtol=1e-10, atol=1e-12)
    assert res.partitions == ref.partitions


@pytest.mark.parametrize("name", ["ordinal", "greedy", "tree", "mixed"])
@pytest.mark.parametrize("dx,dy", SHAPES)
def test_matches_loop(name, dx, dy):
    rng = np.random.default_rng(1000 * dx + dy)
    xs, ys = _structures(name, dx, dy)
    for _ in range(3):
        T, Sigma = _random_ts(dx, dy, rng)
        _assert_same(greedy_search(T, Sigma, dx, dy, xs, ys),
                     _greedy_search_loop(T, Sigma, dx, dy, xs, ys))


def test_low_rank_sigma():
    # n < p: Sigma is rank-deficient, the regime the carried updates are most exposed to.
    rng = np.random.default_rng(7)
    dx, dy = 8, 6
    T, Sigma = _random_ts(dx, dy, rng, n=20)
    assert np.linalg.matrix_rank(Sigma) < dx * dy
    for name in ("greedy", "tree"):
        xs, ys = _structures(name, dx, dy)
        _assert_same(greedy_search(T, Sigma, dx, dy, xs, ys),
                     _greedy_search_loop(T, Sigma, dx, dy, xs, ys))


def test_nary_tree_fallback():
    # A ternary node means groups are not always tree nodes -> generic permitted_merges path.
    leaves = [_make_leaf(i) for i in range(1, 8)]
    tree = _make_parent([_make_parent(leaves[0:3]), _make_parent(leaves[3:5]),
                         _make_parent(leaves[5:7])])
    xs = Tree(tree)
    assert xs._sibling_table is None
    rng = np.random.default_rng(3)
    for _ in range(3):
        T, Sigma = _random_ts(7, 5, rng)
        _assert_same(greedy_search(T, Sigma, 7, 5, xs, Ordinal()),
                     _greedy_search_loop(T, Sigma, 7, 5, xs, Ordinal()))


@pytest.mark.parametrize("name", ["ordinal", "greedy", "tree"])
def test_batch_matches_columns(name):
    # Draws take different paths inside one chunk; each column must equal its own search.
    rng = np.random.default_rng(11)
    dx, dy = 6, 5
    xs, ys = _structures(name, dx, dy)
    _, Sigma = _random_ts(dx, dy, rng)
    T = np.linalg.cholesky(Sigma + 1e-9 * np.eye(dx * dy)) @ rng.standard_normal((dx * dy, 40))
    paths = greedy_search_paths(T, Sigma, dx, dy, xs, ys)
    ref = np.column_stack([_greedy_search_loop(T[:, b], Sigma, dx, dy, xs, ys).values
                           for b in range(T.shape[1])])
    np.testing.assert_allclose(paths, ref, rtol=1e-10, atol=1e-12)


def test_chunking_and_threads_do_not_change_paths():
    rng = np.random.default_rng(5)
    dx, dy = 5, 6
    xs, ys = Tree.binary(dx), Tree.binary(dy)
    _, Sigma = _random_ts(dx, dy, rng)
    T = rng.standard_normal((dx * dy, 25))
    whole = greedy_search_paths(T, Sigma, dx, dy, xs, ys)
    tiny = greedy_search_paths(T, Sigma, dx, dy, xs, ys, chunk_bytes=1)  # one draw per chunk
    threaded = greedy_search_paths(T, Sigma, dx, dy, xs, ys, chunk_bytes=8 * 30 * 30 * 4, n_jobs=3)
    np.testing.assert_array_equal(whole, tiny)
    np.testing.assert_array_equal(whole, threaded)


def test_permitted_mask_matches_permitted_merges():
    # The vectorised masks must agree with permitted_merges on arbitrary reachable partitions.
    rng = np.random.default_rng(2)
    d = 9
    for structure in (Ordinal(), Saturated(), Tree.binary(d)):
        partition = [[i] for i in range(1, d + 1)]
        while True:
            pairs = structure.permitted_merges(partition)
            sizes = np.zeros((1, d), dtype=np.int64)
            gid = np.zeros((1, d), dtype=np.int64)
            for grp in partition:
                sizes[0, min(grp) - 1] = len(grp)
                gid[0, np.asarray(grp) - 1] = min(grp) - 1
            mask = structure.permitted_mask(sizes, gid)[0]
            slot = [min(g) - 1 for g in partition]
            got = sorted((slot.index(s) + 1, slot.index(t) + 1) for s, t in zip(*np.nonzero(mask)))
            assert got == pairs  # same set, and permitted_merges order is lexicographic
            if not pairs:
                break
            i, j = pairs[rng.integers(len(pairs))]
            partition[i - 1] = partition[i - 1] + partition[j - 1]
            del partition[j - 1]


@pytest.mark.parametrize("sparse", [(False, False), (True, True), (True, False), (False, True)])
@pytest.mark.parametrize("threshold", [0.0, 0.6, 1.0])  # never / default / every level
@pytest.mark.parametrize("name", ["greedy", "tree", "ordinal"])
def test_every_mode_matches_loop(name, sparse, threshold):
    # Both carrying modes and every compaction schedule must give the same paths, on
    # any structure -- the default only picks between them for speed.
    from catci.search import _run_chunk, _Shared

    rng = np.random.default_rng(17)
    dx, dy = 7, 6
    xs, ys = _structures(name, dx, dy)
    _, Sigma = _random_ts(dx, dy, rng)
    T = np.linalg.cholesky(Sigma + 1e-9 * np.eye(dx * dy)) @ rng.standard_normal((dx * dy, 12))
    paths, _ = _run_chunk(_Shared(Sigma, dx, dy), T, xs, ys,
                          compact_threshold=threshold, sparse=sparse)
    ref = np.column_stack([_greedy_search_loop(T[:, b], Sigma, dx, dy, xs, ys).values
                           for b in range(T.shape[1])])
    np.testing.assert_allclose(paths, ref, rtol=1e-10, atol=1e-12)
