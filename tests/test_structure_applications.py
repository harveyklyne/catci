"""Structures the applications need (TODO item 3): n-ary taxonomies and cycles.

* ``Cyclic``: groups must stay arcs in cyclic order, so that ``(i, i+1)`` plus
  ``(1, d)`` are exactly the adjacent-arc merges at every depth.
* ``Tree`` of any arity: permitted merges must match an independent oracle
  ("the union is a node, or a union of >= 2 children of one node"), the search
  must never strand, and binary trees must behave exactly as before.
* ``Tree.from_parents``: taxonomy (code -> parent) to tree, with pruning of
  unobserved branches and collapsing of single-child chains.
"""

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from catci.search import greedy_search
from catci.structure import Cyclic, Ordinal, Tree, TreeNode, make_binary_tree


def _initial_partition(d):
    return [[j] for j in range(1, d + 1)]


def _shape(node):
    return {"root": list(node.root), "children": [_shape(c) for c in node.children]}


def _merge(partition, i, j):
    out = [list(g) for g in partition]
    out[i - 1] = out[i - 1] + out[j - 1]
    del out[j - 1]
    return out


def _random_path(structure, d, choices):
    """Apply permitted merges (picked by ``choices``) down to 2 groups."""
    partition = _initial_partition(d)
    path = [partition]
    for c in choices:
        pairs = structure.permitted_merges(partition)
        if len(partition) <= 2:
            assert pairs == []
            break
        assert pairs, f"stranded at {partition}"
        partition = _merge(partition, *pairs[c % len(pairs)])
        path.append(partition)
    return path


def _random_TS(dx, dy, rng):
    A = rng.standard_normal((dx * dy, dx * dy))
    return rng.standard_normal(dx * dy) * 3, A @ A.T / (dx * dy)


# --------------------------------------------------------------------------- #
# Cyclic
# --------------------------------------------------------------------------- #
def _is_arc(group, d):
    labels = set(group)
    if len(labels) == d:
        return True
    starts = [s for s in labels if (s - 2) % d + 1 not in labels]
    if len(starts) != 1:
        return False
    return labels == {(starts[0] - 1 + k) % d + 1 for k in range(len(labels))}


def _arc_order(group, d):
    labels = set(group)
    start = next(s for s in labels if (s - 2) % d + 1 not in labels)
    return [(start - 1 + k) % d + 1 for k in range(len(labels))]


def _assert_arcs_in_cyclic_order(partition, d):
    assert all(_is_arc(g, d) for g in partition), partition
    walk = [lab for g in partition for lab in _arc_order(g, d)]
    r = walk.index(1)
    assert walk[r:] + walk[:r] == list(range(1, d + 1)), partition


def test_cyclic_pairs():
    assert Cyclic().permitted_merges(_initial_partition(4)) == [(1, 2), (2, 3), (3, 4), (1, 4)]
    assert Cyclic().permitted_merges(_initial_partition(2)) == []
    # after the wrap merge, the arc {4, 1} sits at position 1, adjacent to both 2 and 3
    assert Cyclic().permitted_merges([[1, 4], [2], [3]]) == [(1, 2), (2, 3), (1, 3)]


@settings(max_examples=300, deadline=None)
@given(d=st.integers(3, 24), choices=st.lists(st.integers(0, 10**6), min_size=24, max_size=24))
def test_cyclic_groups_stay_arcs(d, choices):
    for partition in _random_path(Cyclic(), d, choices):
        _assert_arcs_in_cyclic_order(partition, d)
        n = len(partition)
        if n > 2:
            adjacent = [(i, j) for i in range(1, n) for j in range(i + 1, n + 1)
                        if _is_arc(partition[i - 1] + partition[j - 1], d)]
            assert sorted(Cyclic().permitted_merges(partition)) == adjacent


@settings(max_examples=100, deadline=None)
@given(d=st.integers(3, 16), choices=st.lists(st.integers(0, 10**6), min_size=16, max_size=16))
def test_ordinal_groups_stay_intervals(d, choices):
    for partition in _random_path(Ordinal(), d, choices):
        assert [lab for g in partition for lab in sorted(g)] == list(range(1, d + 1))


@pytest.mark.parametrize("dx,dy", [(7, 12), (12, 24), (3, 7)])
def test_cyclic_search_reaches_two_arcs(dx, dy):
    T, Sigma = _random_TS(dx, dy, np.random.default_rng(dx * dy))
    res = greedy_search(T, Sigma, dx, dy, Cyclic(), Cyclic())
    assert len(res.values) == (dx - 2) + (dy - 2) + 1
    for p in res.partitions:
        _assert_arcs_in_cyclic_order(p["x"], dx)
        _assert_arcs_in_cyclic_order(p["y"], dy)


# --------------------------------------------------------------------------- #
# n-ary trees
# --------------------------------------------------------------------------- #
def _nodes(node):
    yield node
    for c in node.children:
        yield from _nodes(c)


def _is_valid_group(group, tree):
    """Independent oracle: a node, or a union of >= 2 children of one node."""
    g = set(group)
    for node in _nodes(tree):
        if set(node.root) == g:
            return True
        if (node.children and g <= set(node.root)
                and all(set(c.root) <= g or not set(c.root) & g for c in node.children)):
            return True
    return False


def _random_tree(labels, rng, max_arity=4):
    if len(labels) == 1:
        return TreeNode(root=(int(labels[0]),), children=())
    k = int(rng.integers(2, min(max_arity, len(labels)) + 1))
    labels = [int(x) for x in rng.permutation(labels)]
    cuts = sorted(int(c) for c in rng.choice(np.arange(1, len(labels)), size=k - 1, replace=False))
    parts = [labels[a:b] for a, b in zip([0, *cuts], [*cuts, len(labels)])]
    return TreeNode(root=tuple(sorted(labels)),
                    children=tuple(_random_tree(p, rng, max_arity) for p in parts))


def test_nary_partial_sibling_merge_regression():
    # {1,2,3} is a ternary node. After 2 and 3 merge, the block {2,3} must still
    # be able to absorb 1; the old "equals a sibling node" rule stranded it.
    leaf = lambda k: TreeNode((k,), ())
    t = Tree(TreeNode((1, 2, 3, 4), (TreeNode((1, 2, 3), (leaf(1), leaf(2), leaf(3))), leaf(4))))
    assert t.permitted_merges([[1], [2], [3], [4]]) == [(1, 2), (1, 3), (2, 3)]
    assert t.permitted_merges([[1], [2, 3], [4]]) == [(1, 2)]


@settings(max_examples=200, deadline=None)
@given(d=st.integers(3, 20), seed=st.integers(0, 2**32 - 1),
       choices=st.lists(st.integers(0, 10**6), min_size=20, max_size=20))
def test_nary_tree_merges_match_oracle(d, seed, choices):
    tree = _random_tree(list(range(1, d + 1)), np.random.default_rng(seed))
    structure = Tree(tree)
    for partition in _random_path(structure, d, choices):  # asserts never stranded
        assert all(_is_valid_group(g, tree) for g in partition)
        n = len(partition)
        if n > 2:
            expected = [(i, j) for i in range(1, n) for j in range(i + 1, n + 1)
                        if _is_valid_group(partition[i - 1] + partition[j - 1], tree)]
            assert structure.permitted_merges(partition) == expected


@settings(max_examples=50, deadline=None)
@given(d=st.integers(3, 33), choices=st.lists(st.integers(0, 10**6), min_size=33, max_size=33))
def test_binary_tree_unchanged_under_new_rule(d, choices):
    # For binary trees the context rule is exactly the old "merge sibling nodes".
    tree = make_binary_tree(d)
    sibling_pairs = [{c.root for c in node.children} for node in _nodes(tree) if node.children]
    for partition in _random_path(Tree(tree), d, choices):
        n = len(partition)
        if n > 2:
            siblings = [(i, j) for i in range(1, n) for j in range(i + 1, n + 1)
                        if {tuple(sorted(partition[i - 1])), tuple(sorted(partition[j - 1]))}
                        in sibling_pairs]
            assert Tree(tree).permitted_merges(partition) == siblings


def test_nary_search_reaches_two_by_two():
    rng = np.random.default_rng(0)
    dx, dy = 9, 10
    tx = _random_tree(list(range(1, dx + 1)), rng)
    ty = _random_tree(list(range(1, dy + 1)), rng)
    T, Sigma = _random_TS(dx, dy, rng)
    res = greedy_search(T, Sigma, dx, dy, Tree(tx), Tree(ty))
    assert len(res.values) == (dx - 2) + (dy - 2) + 1
    assert len(res.partitions[-1]["x"]) == 2 and len(res.partitions[-1]["y"]) == 2


# --------------------------------------------------------------------------- #
# Tree.from_parents
# --------------------------------------------------------------------------- #
ICD = {
    # code -> block -> chapter; keyed on (kind, id) so namespaces cannot collide
    ("code", "250.0"): ("block", "250"), ("code", "250.1"): ("block", "250"),
    ("code", "250.9"): ("block", "250"),
    ("block", "250"): ("chap", "endocrine"),
    ("code", "401.9"): ("block", "401"), ("code", "428.0"): ("block", "428"),
    ("block", "401"): ("chap", "circulatory"), ("block", "428"): ("chap", "circulatory"),
    ("code", "V58.6"): ("block", "V58"), ("block", "V58"): ("chap", "supplementary"),
    ("code", "999.9"): ("block", "999"), ("block", "999"): ("chap", "injury"),  # unobserved
}


def test_from_parents_prunes_and_collapses():
    levels = [("code", c) for c in ["250.0", "401.9", "250.1", "428.0", "V58.6"]]
    t = Tree.from_parents(ICD, levels)
    # root -> {250: {250.0, 250.1}, circulatory: {401.9, 428.0}, V58.6}. The
    # unobserved 250.9 and injury chapter are pruned; the single-child chains
    # (endocrine -> 250, 401 -> 401.9, supplementary -> V58 -> V58.6) collapse.
    leaf = lambda k: TreeNode((k,), ())
    assert _shape(t.tree) == _shape(TreeNode((1, 2, 3, 4, 5), (
        TreeNode((1, 3), (leaf(1), leaf(3))),
        TreeNode((2, 4), (leaf(2), leaf(4))),
        leaf(5),
    )))
    # 5 hangs off the root, so it waits until a chapter is whole before merging.
    assert t.permitted_merges(_initial_partition(5)) == [(1, 3), (2, 4)]
    assert t.permitted_merges([[1, 3], [2], [4], [5]]) == [(1, 4), (2, 3)]


def test_from_parents_matches_binary_tree():
    for d in (2, 3, 5, 8, 13):
        tree = make_binary_tree(d)
        parents = {c.root: node.root for node in _nodes(tree) for c in node.children}
        t = Tree.from_parents(parents, [(k,) for k in range(1, d + 1)])
        assert _shape(t.tree) == _shape(tree)


def test_from_parents_multiple_roots_and_single_level():
    t = Tree.from_parents({"a": "X", "b": "X", "c": "Y", "d": None}, ["a", "b", "c", "d"])
    # root -> {X: {a, b}, c, d}
    assert t.permitted_merges(_initial_partition(4)) == [(1, 2), (3, 4)]
    assert t.permitted_merges([[1, 2], [3], [4]]) == [(1, 2), (1, 3), (2, 3)]
    assert Tree.from_parents({"a": None}, ["a"]).tree == TreeNode((1,), ())


@pytest.mark.parametrize("parents,levels,msg", [
    ({"a": "X"}, ["a", "b"], "no entry"),
    ({"a": "X", "X": "Y", "Y": "X"}, ["a"], "Cycle"),
    ({"a": "b", "b": None}, ["a", "b"], "ancestors"),
    ({"a": None}, ["a", "a"], "distinct"),
])
def test_from_parents_rejects(parents, levels, msg):
    with pytest.raises(ValueError, match=msg):
        Tree.from_parents(parents, levels)


def test_tree_rejects_malformed():
    leaf = lambda k: TreeNode((k,), ())
    with pytest.raises(ValueError, match="single child"):
        Tree(TreeNode((1, 2), (TreeNode((1, 2), (leaf(1), leaf(2))),)))
    with pytest.raises(ValueError, match="partition"):
        Tree(TreeNode((1, 2, 3), (leaf(1), leaf(2))))
