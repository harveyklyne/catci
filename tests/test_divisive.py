"""Tests for the divisive (top-down) search.

There is no R oracle for this path -- it does not exist in the R package -- so
it is pinned three ways instead:

1. Against dense recomputation. ``divisive_search`` never materialises the
   aggregated ``(T, Sigma)``; it carries ``(||T||^2, tr, tr2)`` through the
   inverse update formulae and reads two rows at a time out of a block-sum
   table. Every level it reports must equal ``phi(Pi T, Pi Sigma Pi^T)`` for the
   partition it says it is at. This is the test a silent index error fails.
2. Against ``greedy_search``, structurally: same lattice, opposite direction, so
   the two must agree on the levels they share (the coarsest and the finest) and
   visit the same number in between -- while differing somewhere in the middle,
   or there would be nothing to compare.
3. By calibration, as for the merge direction: a valid search calibrates, so
   this catches a broken *bootstrap* path even though it cannot catch a broken
   *search* path (both directions are valid tests, so both calibrate).
"""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from catci.statistic import ApproxChi
from catci.calibrate import adaptive_pvalue
from catci.search import SplitSearch, divisive_search, greedy_search
from catci.structure import Ordinal, Saturated, Tree


def random_t_sigma(dx, dy, rng):
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    return rng.standard_normal(p), Sigma


def pi_matrix(x_groups, y_groups, dx, dy):
    """Pi(A, B) of the paper's (7), in the package's X-fastest ordering."""
    kx = len(x_groups)
    P = np.zeros((kx * len(y_groups), dx * dy))
    for gy, y_group in enumerate(y_groups):
        for gx, x_group in enumerate(x_groups):
            for k in y_group:
                for j in x_group:
                    P[gy * kx + gx, (k - 1) * dx + (j - 1)] = 1.0
    return P


def dense_values(result, T, Sigma, dx, dy):
    """The statistic recomputed from scratch at each partition the search visited."""
    statistic = ApproxChi()
    values = []
    for partition in result.partitions:
        P = pi_matrix(partition["x"], partition["y"], dx, dy)
        values.append(statistic.value(statistic.init(P @ T, P @ Sigma @ P.T)))
    return np.array(values)


def structures(kind, dx, dy):
    return (Tree.binary(dx), Tree.binary(dy)) if kind == "tree" else (Ordinal(), Ordinal())


# --------------------------------------------------------------------------- #
# 1. incremental state == dense recompute
# --------------------------------------------------------------------------- #
@settings(max_examples=100, deadline=None)
@given(
    dx=st.integers(min_value=2, max_value=9),
    dy=st.integers(min_value=2, max_value=9),
    kind=st.sampled_from(["tree", "ordinal"]),
    seed=st.integers(min_value=0, max_value=2 ** 32 - 1),
)
def test_divisive_values_match_dense_recompute(dx, dy, kind, seed):
    rng = np.random.default_rng(seed)
    T, Sigma = random_t_sigma(dx, dy, rng)
    result = divisive_search(T, Sigma, dx, dy, *structures(kind, dx, dy))
    assert np.allclose(result.values, dense_values(result, T, Sigma, dx, dy))


@pytest.mark.parametrize("dx, dy", [(30, 10), (16, 4)])
def test_divisive_matches_dense_at_the_motivating_size(dx, dy):
    """The paper's dX x dY, where a wrong stride would still look plausible."""
    rng = np.random.default_rng(7)
    T, Sigma = random_t_sigma(dx, dy, rng)
    for kind in ("tree", "ordinal"):
        if kind == "tree" and (dx & (dx - 1) or dy & (dy - 1)):
            structs = (Tree.binary(dx), Tree.binary(dy))  # make_binary_tree handles non-powers
        else:
            structs = structures(kind, dx, dy)
        result = divisive_search(T, Sigma, dx, dy, *structs)
        assert np.allclose(result.values, dense_values(result, T, Sigma, dx, dy)), kind


# --------------------------------------------------------------------------- #
# 2. structural agreement with, and difference from, the merge direction
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("dx, dy", [(4, 4), (8, 8), (8, 4), (5, 3), (16, 4)])
@pytest.mark.parametrize("kind", ["tree", "ordinal"])
def test_divisive_spans_the_same_lattice_as_merging(dx, dy, kind):
    rng = np.random.default_rng(dx * 100 + dy)
    T, Sigma = random_t_sigma(dx, dy, rng)
    structs = structures(kind, dx, dy)

    merge = greedy_search(T, Sigma, dx, dy, *structs)
    split = divisive_search(T, Sigma, dx, dy, *structs)

    assert len(split.values) == len(merge.values) == dx + dy - 3
    # merging starts at the singleton partition; splitting ends there
    assert split.partitions[-1] == merge.partitions[0]
    assert split.partitions[-1]["x"] == [[j] for j in range(1, dx + 1)]
    assert split.partitions[-1]["y"] == [[k] for k in range(1, dy + 1)]
    assert np.isclose(split.values[-1], merge.values[0])


def test_divisive_differs_from_merging():
    """The assertion that makes the comparison worth running at all.

    A tree has one two-group partition, so both directions must agree at that
    end too -- but if they agreed everywhere, the divisive path would be the
    merge path reversed and there would be no separate method here.
    """
    rng = np.random.default_rng(0)
    dx = dy = 8
    T, Sigma = random_t_sigma(dx, dy, rng)
    structs = (Tree.binary(dx), Tree.binary(dy))

    merge = greedy_search(T, Sigma, dx, dy, *structs)
    split = divisive_search(T, Sigma, dx, dy, *structs)

    assert np.isclose(split.values[0], merge.values[-1])  # same unique 2x2 partition
    assert split.partitions != merge.partitions[::-1]


@pytest.mark.parametrize("kind", ["tree", "ordinal"])
def test_truncation_is_a_prefix_of_the_full_path(kind):
    rng = np.random.default_rng(3)
    dx = dy = 8
    T, Sigma = random_t_sigma(dx, dy, rng)
    structs = structures(kind, dx, dy)

    full = divisive_search(T, Sigma, dx, dy, *structs)
    for levels in range(dx + dy - 3):
        short = divisive_search(T, Sigma, dx, dy, *structs, max_levels=levels)
        assert len(short.values) == levels + 1
        assert np.allclose(short.values, full.values[: levels + 1])
        assert short.partitions == full.partitions[: levels + 1]


# --------------------------------------------------------------------------- #
# 3. the structure operations themselves
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("d", [2, 3, 4, 5, 8, 16])
def test_coarsest_partitions_have_two_groups_covering_the_labels(d):
    for structure in (Ordinal(), Tree.binary(d)):
        for partition in structure.coarsest_partitions(d):
            assert len(partition) == 2
            assert sorted(j for group in partition for j in group) == list(range(1, d + 1))
    assert len(Ordinal().coarsest_partitions(d)) == d - 1
    assert len(Tree.binary(d).coarsest_partitions(d)) == 1


@pytest.mark.parametrize("d", [4, 5, 8])
def test_splitting_then_merging_is_the_identity(d):
    """A split names a merge that undoes it -- the invariant the inverse formulae assume."""
    for structure in (Ordinal(), Tree.binary(d)):
        partition = [list(g) for g in structure.coarsest_partitions(d)[0]]
        for (i, a, b) in structure.permitted_splits(partition):
            refined = partition[: i - 1] + [a, b] + partition[i:]
            assert (i, i + 1) in structure.permitted_merges(refined)
            merged = refined[: i - 1] + [refined[i - 1] + refined[i]] + refined[i + 1:]
            assert merged == partition


def test_a_groups_splits_do_not_depend_on_the_rest_of_the_partition():
    """What ``_splits_of``'s memo relies on."""
    for structure in (Ordinal(), Tree.binary(8)):
        whole = structure.permitted_splits([[1, 2, 3, 4], [5, 6], [7, 8]])
        alone = structure.permitted_splits([[1, 2, 3, 4]])
        assert [(a, b) for (i, a, b) in whole if i == 1] == [(a, b) for (_, a, b) in alone]


def test_saturated_refuses_to_split():
    with pytest.raises(NotImplementedError, match="no tractable divisive counterpart"):
        Saturated().permitted_splits([[1, 2], [3, 4]])
    with pytest.raises(NotImplementedError):
        Saturated().coarsest_partitions(4)


def test_non_contiguous_groups_are_rejected():
    """Permuted labels break the range representation, and must say so."""
    rng = np.random.default_rng(0)
    T, Sigma = random_t_sigma(4, 4, rng)

    class Scrambled(Ordinal):
        def coarsest_partitions(self, d):
            return [[[1, 3], [2, 4]]]

    with pytest.raises(ValueError, match="contiguous in label order"):
        divisive_search(T, Sigma, 4, 4, Scrambled(), Ordinal())


# --------------------------------------------------------------------------- #
# 4. calibration
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "kind, max_levels", [("ordinal", None), ("tree", None), ("tree", 2)]
)
def test_divisive_calibrates_under_gaussian(kind, max_levels):
    rng = np.random.default_rng(11)
    dx = dy = 4
    p = dx * dy
    A = rng.standard_normal((p, p))
    Sigma = A @ A.T / p + np.eye(p)
    sqrt_Sigma = np.linalg.cholesky(Sigma)
    structs = structures(kind, dx, dy)

    pvals = np.array([
        adaptive_pvalue(
            sqrt_Sigma @ rng.standard_normal(p), Sigma, dx, dy, *structs,
            n_boot=100, rng=rng, search=SplitSearch(max_levels),
        )
        for _ in range(400)
    ])

    assert abs(pvals.mean() - 0.5) < 0.06
    for alpha in (0.05, 0.10, 0.20):
        rate = float(np.mean(pvals < alpha))
        assert abs(rate - alpha) < 0.05, f"rejection rate {rate:.3f} at alpha={alpha}"
