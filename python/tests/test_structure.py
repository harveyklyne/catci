"""The structure object's permitted_merges must reproduce the R permitted-merge maps.

This pins the load-bearing redesign: ordinal/greedy/tree permitted merges on the
initial partition, and the binary-tree shape, for d in {2, 4, 8}.
"""

import numpy as np

from catci.structure import Ordinal, Saturated, Tree, make_binary_tree


def _initial_partition(d: int):
    return [[j] for j in range(1, d + 1)]


def _pairs(structure, d):
    return [list(p) for p in structure.permitted_merges(_initial_partition(d))]


def test_permitted_merges_match_oracle(oracle):
    for entry in oracle["tree_structure"]["permitted"]:
        d = entry["d"]
        assert _pairs(Ordinal(), d) == [list(p) for p in entry["ordinal"]["pairs"]]
        assert _pairs(Saturated(), d) == [list(p) for p in entry["greedy"]["pairs"]]
        assert _pairs(Tree.binary(d), d) == [list(p) for p in entry["tree"]["pairs"]]
        # num_levels is just the count -- get_num_levels does not exist here.
        assert len(_pairs(Ordinal(), d)) == entry["ordinal"]["num_levels"]
        assert len(_pairs(Tree.binary(d), d)) == entry["tree"]["num_levels"]


def _shape(node):
    return {"root": list(node.root), "children": [_shape(c) for c in node.children]}


def test_binary_tree_shape_matches_oracle(oracle):
    for entry in oracle["tree_structure"]["trees"]:
        d = entry["d"]
        assert _shape(make_binary_tree(d)) == entry["shape"]


def test_guard_at_d_equals_2():
    # every structure returns [] at 2 groups (the (d > 2) guard, by construction)
    part2 = [[1], [2]]
    assert Ordinal().permitted_merges(part2) == []
    assert Saturated().permitted_merges(part2) == []
    assert Tree.binary(2).permitted_merges(part2) == []
