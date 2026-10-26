"""Depth-truncated merge search and the random-chooser baselines.

* ``MergeSearch(max_levels=k)`` is the first ``k + 1`` entries of the full path (and
  stops early), as ``SplitSearch(max_levels=k)`` is.
* ``randomsearch`` splits a search into a generator of partition paths (never sees
  ``T``) and a scorer. The scorer is pinned by feeding it the partitions the
  *directed* searches visit, which must reproduce their values; the generator by its
  invariants (data-free, nested, fresh per column, uniform, right length) and by
  calibration under a Gaussian ``T``.
"""

import numpy as np
import pytest
from scipy.stats import chisquare

from catci.calibrate import adaptive_pvalue
from catci.randomsearch import (
    PathSet,
    RandomMergeSearch,
    RandomSplitSearch,
    _canon,
    _merge_children,
    _split_children,
    generate_merge_paths,
    generate_split_paths,
    score_paths,
)
from catci.search import (
    MergeSearch,
    SplitSearch,
    divisive_search,
    greedy_search,
    greedy_search_paths,
)
from catci.statistic import ApproxChi, ExactChi
from catci.structure import Ordinal, Saturated, Tree

from test_divisive import dense_values, pi_matrix  # noqa: F401


def gaussian_draws(dx, dy, B, seed):
    rng = np.random.default_rng(seed)
    p = dx * dy
    A = rng.standard_normal((p + 5, p))
    Sigma = A.T @ A / (p + 5)
    return np.linalg.cholesky(Sigma) @ rng.standard_normal((p, B)), Sigma


def make_structures(kind, dx, dy):
    return {"ordinal": (Ordinal(), Ordinal()),
            "saturated": (Saturated(), Saturated()),
            "tree": (Tree.binary(dx), Tree.binary(dy))}[kind]


SPLITTABLE = ["ordinal", "tree"]
ALL_KINDS = ["ordinal", "saturated", "tree"]
SIZES = [(6, 5), (7, 2), (2, 5), (2, 2), (3, 3)]


# --------------------------------------------------------------------------- #
# Truncated merge search
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kind", ALL_KINDS)
@pytest.mark.parametrize("k", [0, 1, 3])
def test_merge_truncation_is_a_prefix(kind, k):
    dx, dy = 6, 5
    T, Sigma = gaussian_draws(dx, dy, 12, seed=1)
    xs, ys = make_structures(kind, dx, dy)
    full = greedy_search_paths(T, Sigma, dx, dy, xs, ys)
    cut = MergeSearch(max_levels=k).paths(T, Sigma, dx, dy, xs, ys)
    assert cut.shape == (k + 1, 12)
    np.testing.assert_array_equal(cut, full[:k + 1])


def test_merge_truncation_beyond_the_path_is_the_full_path():
    T, Sigma = gaussian_draws(4, 4, 5, seed=2)
    full = greedy_search_paths(T, Sigma, 4, 4, Ordinal(), Ordinal())
    cut = MergeSearch(max_levels=99).paths(T, Sigma, 4, 4, Ordinal(), Ordinal())
    np.testing.assert_array_equal(cut, full)


def test_truncated_merge_refuses_what_it_cannot_do():
    T, Sigma = gaussian_draws(4, 4, 3, seed=2)
    with pytest.raises(NotImplementedError):
        MergeSearch(max_levels=2).result(T[:, 0], Sigma, 4, 4, Ordinal(), Ordinal())


# --------------------------------------------------------------------------- #
# Generator invariants
# --------------------------------------------------------------------------- #
def generate(direction, dx, dy, kind, B, seed, max_levels=None):
    xs, ys = make_structures(kind, dx, dy)
    gen = generate_merge_paths if direction == "merge" else generate_split_paths
    return gen(dx, dy, xs, ys, B, np.random.default_rng(seed), max_levels)


def visited(ps):
    """Per column, the sequence of partition keys."""
    return [[ps.keys[i] for i in ps.node_of[:, b]] for b in range(ps.node_of.shape[1])]


CASES = [("merge", k) for k in ALL_KINDS] + [("split", k) for k in SPLITTABLE]


@pytest.mark.parametrize("direction, kind", CASES)
def test_same_seed_same_paths_and_different_seed_differs(direction, kind):
    a, b = generate(direction, 6, 5, kind, 30, 3), generate(direction, 6, 5, kind, 30, 3)
    c = generate(direction, 6, 5, kind, 30, 4)
    assert visited(a) == visited(b)
    assert visited(a) != visited(c)


@pytest.mark.parametrize("direction, kind", CASES)
def test_depths_are_nested(direction, kind):
    full = generate(direction, 7, 4, kind, 40, 9)
    for k in (0, 1, 3):
        cut = generate(direction, 7, 4, kind, 40, 9, max_levels=k)
        assert cut.node_of.shape[0] == k + 1
        assert visited(cut) == [path[:k + 1] for path in visited(full)]


@pytest.mark.parametrize("direction, kind", CASES)
def test_fresh_path_per_column(direction, kind):
    ps = generate(direction, 8, 8, kind, 50, 1)
    assert len({tuple(path) for path in visited(ps)}) > 5


@pytest.mark.parametrize("direction, kind", CASES)
@pytest.mark.parametrize("dx, dy", SIZES)
def test_path_length_equals_the_directed_length(direction, kind, dx, dy):
    T, Sigma = gaussian_draws(dx, dy, 3, seed=5)
    xs, ys = make_structures(kind, dx, dy)
    directed = (MergeSearch() if direction == "merge" else SplitSearch()).paths(
        T, Sigma, dx, dy, xs, ys)
    ps = generate(direction, dx, dy, kind, 3, 0)
    assert ps.node_of.shape[0] == directed.shape[0] == dx + dy - 3


@pytest.mark.parametrize("direction, kind", CASES)
def test_every_step_is_a_permitted_move(direction, kind):
    dx, dy = 6, 5
    xs, ys = make_structures(kind, dx, dy)
    children = _merge_children if direction == "merge" else _split_children
    for path in visited(generate(direction, dx, dy, kind, 25, 2)):
        for prev, nxt in zip(path, path[1:]):
            assert nxt in children(prev, xs, ys)
    # starting points
    starts = {path[0] for path in visited(generate(direction, dx, dy, kind, 25, 2))}
    if direction == "merge":
        assert starts == {(tuple((j,) for j in range(1, dx + 1)), tuple((k,) for k in range(1, dy + 1)))}
    else:
        assert all(len(px) == 2 and len(py) == 2 for px, py in starts)


def pooled_counts(direction, dx, dy, level, seeds, B=4000):
    """Visit counts per partition key at ``level``, pooled over seeds."""
    tally = {}
    for seed in seeds:
        ps = generate(direction, dx, dy, "ordinal", B, seed, max_levels=level)
        for key_id, c in zip(*np.unique(ps.node_of[level], return_counts=True)):
            tally[ps.keys[key_id]] = tally.get(ps.keys[key_id], 0) + int(c)
    return tally


def test_split_start_is_uniform_over_coarsest_pairs():
    tally = pooled_counts("split", 5, 5, 0, range(10))  # (dx-1)(dy-1) = 16 pairs
    assert len(tally) == 16
    assert chisquare(list(tally.values())).pvalue > 1e-3


def test_first_step_is_uniform_and_pools_x_and_y():
    dx, dy = 5, 4
    tally = pooled_counts("merge", dx, dy, 1, range(10))
    assert len(tally) == (dx - 1) + (dy - 1)  # both dimensions pooled
    assert chisquare(list(tally.values())).pvalue > 1e-3


def test_observed_column_is_drawn_like_the_others():
    """Column 0 (the observed statistic) must not be treated differently."""
    dx, dy, reps = 4, 4, 3000
    tally = {0: {}, 1: {}}
    for seed in range(reps):
        ps = generate("merge", dx, dy, "ordinal", 2, seed, max_levels=1)
        for col in (0, 1):
            key = ps.keys[ps.node_of[1, col]]
            tally[col][key] = tally[col].get(key, 0) + 1
    for col in (0, 1):
        counts = np.array(list(tally[col].values()))
        assert len(counts) == 6
        assert chisquare(counts).pvalue > 1e-3


@pytest.mark.parametrize("dx, dy", [(7, 2), (2, 6), (2, 2)])
def test_degenerate_shapes(dx, dy):
    for direction in ("merge", "split"):
        ps = generate(direction, dx, dy, "ordinal", 10, 0)
        assert ps.node_of.shape[0] == dx + dy - 3
    ps = generate("merge", dx, dy, "ordinal", 10, 0, max_levels=0)
    assert ps.node_of.shape == (1, 10)


def test_saturated_split_refuses():
    with pytest.raises(NotImplementedError):
        generate("split", 5, 5, "saturated", 4, 0)


def test_result_is_not_implemented():
    with pytest.raises(NotImplementedError):
        RandomMergeSearch().result()


# --------------------------------------------------------------------------- #
# Scorer
# --------------------------------------------------------------------------- #
def pathset_from_partitions(partitions):
    keys = [(_canon(p["x"]), _canon(p["y"])) for p in partitions]
    return PathSet(keys=keys, node_of=np.arange(len(keys))[:, None])


@pytest.mark.parametrize("kind", ALL_KINDS)
@pytest.mark.parametrize("dx, dy", [(6, 5), (7, 2), (4, 4)])
def test_scorer_reproduces_the_directed_merge_search(kind, dx, dy):
    T, Sigma = gaussian_draws(dx, dy, 1, seed=11)
    xs, ys = make_structures(kind, dx, dy)
    directed = greedy_search(T[:, 0], Sigma, dx, dy, xs, ys)
    scored = score_paths(T, Sigma, pathset_from_partitions(directed.partitions), dx, dy)
    np.testing.assert_allclose(scored[:, 0], directed.values, atol=1e-10)


@pytest.mark.parametrize("kind", SPLITTABLE)
@pytest.mark.parametrize("dx, dy", [(6, 5), (7, 2), (4, 4)])
def test_scorer_reproduces_the_directed_split_search(kind, dx, dy):
    T, Sigma = gaussian_draws(dx, dy, 1, seed=12)
    xs, ys = make_structures(kind, dx, dy)
    directed = divisive_search(T[:, 0], Sigma, dx, dy, xs, ys)
    scored = score_paths(T, Sigma, pathset_from_partitions(directed.partitions), dx, dy)
    np.testing.assert_allclose(scored[:, 0], directed.values, atol=1e-10)


@pytest.mark.parametrize("direction, kind", CASES)
def test_scorer_matches_a_dense_recompute_on_random_paths(direction, kind):
    from catci.search import SearchResult

    dx, dy = 6, 5
    T, Sigma = gaussian_draws(dx, dy, 8, seed=13)
    ps = generate(direction, dx, dy, kind, 8, 3)
    values = score_paths(T, Sigma, ps, dx, dy)
    for b in range(8):
        result = SearchResult(values=list(values[:, b]))
        for node in ps.node_of[:, b]:
            px, py = ps.keys[node]
            result.partitions.append({"x": [list(g) for g in px], "y": [list(g) for g in py]})
        np.testing.assert_allclose(values[:, b], dense_values(result, T[:, b], Sigma, dx, dy),
                                   atol=1e-10)


def test_scorer_takes_any_statistic():
    dx, dy = 5, 4
    T, Sigma = gaussian_draws(dx, dy, 6, seed=14)
    ps = generate("merge", dx, dy, "ordinal", 6, 4)

    class SumSquares:
        def init(self, Tv, S):
            return float(np.sum(Tv ** 2))

        def value(self, state):
            return state

    got = score_paths(T, Sigma, ps, dx, dy, SumSquares())
    for b in range(6):
        for l, node in enumerate(ps.node_of[:, b]):
            from catci.randomsearch import _aggregation
            Pi = _aggregation(ps.keys[node], dx, dy)
            assert got[l, b] == pytest.approx(np.sum((Pi @ T[:, b]) ** 2))

    exact = score_paths(T, Sigma, ps, dx, dy, ExactChi())
    assert np.all((exact >= 0) & (exact <= 1))
    approx = score_paths(T, Sigma, ps, dx, dy, ApproxChi())
    assert np.allclose(exact, approx, atol=0.1)  # same statistic, Box vs exact null law


def test_search_paths_use_the_scorer_and_have_the_right_shape():
    dx, dy = 6, 5
    T, Sigma = gaussian_draws(dx, dy, 9, seed=15)
    xs, ys = Ordinal(), Ordinal()
    for cls in (RandomMergeSearch, RandomSplitSearch):
        full = cls(seed=3).paths(T, Sigma, dx, dy, xs, ys)
        assert full.shape == (dx + dy - 3, 9)
        cut = cls(max_levels=2, seed=3).paths(T, Sigma, dx, dy, xs, ys)
        np.testing.assert_array_equal(cut, full[:3])
        assert np.all((full >= 0) & (full <= 1))


def test_zero_variance_coordinate_is_scored_not_nan():
    """A level absent from the sample: a zero row/column of T and Sigma."""
    dx, dy = 5, 3
    T, Sigma = gaussian_draws(dx, dy, 6, seed=16)
    dead = np.arange(dx * dy) % dx == 4  # label 5 of X has no variance
    Sigma = Sigma.copy()
    Sigma[dead, :] = 0.0
    Sigma[:, dead] = 0.0
    T = T.copy()
    T[dead, :] = 0.0
    for cls in (RandomMergeSearch, RandomSplitSearch):
        out = cls(seed=1).paths(T, Sigma, dx, dy, Ordinal(), Ordinal())
        assert np.all(np.isfinite(out))


# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("make", [
    lambda s: RandomMergeSearch(seed=s),
    lambda s: RandomMergeSearch(max_levels=2, seed=s),
    lambda s: RandomSplitSearch(seed=s),
    lambda s: RandomSplitSearch(max_levels=1, seed=s),
    lambda s: MergeSearch(max_levels=2),
], ids=["rmerge", "rmerge@2", "rsplit", "rsplit@1", "merge@2"])
def test_pvalues_are_exactly_uniform_on_the_grid(make):
    """With B = 9 the minP p-value lives on {1/10, ..., 1} and is exactly uniform when
    the columns are exchangeable -- a much sharper check than a rejection rate."""
    rng = np.random.default_rng(31)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    chol = np.linalg.cholesky(Sigma)
    xs, ys = Tree.binary(dx), Tree.binary(dy)
    reps = 1500
    grid = np.array([
        round(10 * adaptive_pvalue(chol @ rng.standard_normal(p), Sigma, dx, dy, xs, ys,
                                   n_boot=9, rng=rng,
                                   search=make(int(rng.integers(2 ** 31)))))
        for _ in range(reps)
    ])
    counts = np.bincount(grid, minlength=11)[1:]
    assert counts.sum() == reps
    assert chisquare(counts).pvalue > 1e-3, counts
