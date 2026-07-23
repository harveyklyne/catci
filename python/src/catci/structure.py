"""Merge structures: the single ``permitted_merges`` operation.

This is the load-bearing redesign from CODE_REVIEW.md 4.1. The R package
represented search structure as three things threaded in lockstep -- a search
string (``"ordinal"``/``"greedy"``/``"tree"``), a ragged ``trees`` list, and a
ragged ``categories`` list -- and a separate ``get_num_levels`` that had to
agree with the merge loop by hand. That split representation *was* the bug
behind finding #1 (the string and the tree disagreeing) and the guard mismatch.

Here a single object per variable answers one question::

    structure.permitted_merges(partition) -> list of (i, j) position pairs

The row count is just ``len(permitted_merges(...))`` -- it cannot disagree with
the loop, so ``get_num_levels`` does not exist. ``partition`` is the current
list of groups (each a list of 1-based original labels); returned ``(i, j)`` are
1-based positions in that partition, ``i < j``. Every structure returns ``[]``
once the partition has ``<= 2`` groups (the ``(d > 2)`` guard, by construction).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

__all__ = [
    "Structure",
    "Ordinal",
    "Saturated",
    "Tree",
    "TreeNode",
    "make_binary_tree",
]

Partition = Sequence[Sequence[int]]
Pair = Tuple[int, int]


class Structure:
    """Base class: a merge structure for one variable."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:  # pragma: no cover
        raise NotImplementedError


class Ordinal(Structure):
    """Adjacent merges only: (1,2), (2,3), ..., (d-1, d)."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(i, i + 1) for i in range(1, d)]


class Saturated(Structure):
    """All pairs (the 'greedy' search in the R package)."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(i, j) for i in range(1, d) for j in range(i + 1, d + 1)]


# --------------------------------------------------------------------------- #
# Tree structure
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TreeNode:
    """A node in a merge tree: ``root`` is the sorted set of original labels."""

    root: Tuple[int, ...]
    children: Tuple["TreeNode", ...]


def _make_leaf(label: int) -> TreeNode:
    return TreeNode(root=(label,), children=())


def _make_parent(children: Sequence[TreeNode]) -> TreeNode:
    root = tuple(sorted({lab for c in children for lab in c.root}))
    return TreeNode(root=root, children=tuple(children))


def _make_binary_above(level: Sequence[TreeNode]) -> List[TreeNode]:
    if len(level) < 2:
        raise ValueError("Level does not contain enough nodes to construct any parents.")
    out = [_make_parent(level[2 * j:2 * j + 2]) for j in range(len(level) // 2)]
    if len(level) % 2 == 1:
        out.append(level[-1])
    return out


def make_binary_tree(d: int) -> TreeNode:
    """Balanced binary tree over in-order leaves ``1..d`` (ports ``make_binary_tree``)."""
    level: List[TreeNode] = [_make_leaf(i) for i in range(1, d + 1)]
    while len(level) > 1:
        level = _make_binary_above(level)
    return level[0]


def _get_siblings(node_labels: Sequence[int], tree: TreeNode) -> List[TreeNode]:
    """Nodes that are siblings of the subtree whose leaf-set is ``node_labels``."""
    node = set(node_labels)
    siblings: List[TreeNode] = []
    if len(set(tree.root) - node) == 0:  # node covers the whole (sub)tree root
        return siblings
    for child in tree.children:
        croot = set(child.root)
        if node <= croot and len(croot - node) > 0:  # descend into the strictly-larger child
            return _get_siblings(node_labels, child)
        if croot & node:
            if not croot <= node:
                raise ValueError("Only part of vertex in node.")
        else:
            siblings.append(child)
    return siblings


class Tree(Structure):
    """Sibling merges only, per a fixed binary tree over the original labels."""

    def __init__(self, tree: TreeNode):
        self.tree = tree

    @classmethod
    def binary(cls, d: int) -> "Tree":
        return cls(make_binary_tree(d))

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        pairs: List[Pair] = []
        for i in range(1, d):  # 1-based position of the first group
            siblings = _get_siblings(partition[i - 1], self.tree)
            sib_sets = [set(s.root) for s in siblings]
            for j in range(i + 1, d + 1):
                if any(set(partition[j - 1]) == s for s in sib_sets):
                    pairs.append((i, j))
        return pairs
