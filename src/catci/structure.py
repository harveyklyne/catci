"""Merge structures: the single ``permitted_merges`` operation.

The load-bearing redesign of the port (README.md, "Design notes"). The R
package represented search structure as three things threaded in lockstep -- a
search string (``"ordinal"``/``"greedy"``/``"tree"``), a ragged ``trees`` list,
and a ragged ``categories`` list -- and a separate ``get_num_levels`` that had
to agree with the merge loop by hand. That split representation *was* the bug
where ``"tree"`` silently ran an ordinal search, plus the guard mismatch it hid.

Here a single object per variable answers one question::

    structure.permitted_merges(partition) -> list of (i, j) position pairs

The row count is just ``len(permitted_merges(...))`` -- it cannot disagree with
the loop, so ``get_num_levels`` does not exist. ``partition`` is the current
list of groups (each a list of 1-based original labels); returned ``(i, j)`` are
1-based positions in that partition, ``i < j``. Every structure returns ``[]``
once the partition has ``<= 2`` groups (the ``(d > 2)`` guard, by construction).

:func:`~catci.search.divisive_search` walks the same lattice top-down, so a
structure answers two more questions::

    structure.coarsest_partitions(d) -> [partition, ...]
    structure.permitted_splits(partition) -> [(i, a, b), ...]

``(i, a, b)`` means "replace the group at 1-based position ``i`` by ``a``, then
``b``", with ``a + b`` a rearrangement of that group and ``b`` inserted directly
after ``a`` -- which keeps groups ordered by least label, the same invariant the
merge loop maintains. A merge has one obvious inverse, but a *start* does not:
merging bottoms out at a unique singleton partition while there are many
two-group partitions, so ``coarsest_partitions`` returns all the structure
permits and the search picks between them by statistic. A tree names its own
coarsest split, so it returns one; an ordinal variable returns ``d - 1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

__all__ = [
    "Split",
    "Structure",
    "Ordinal",
    "Saturated",
    "Tree",
    "TreeNode",
    "make_binary_tree",
]

Partition = Sequence[Sequence[int]]
Pair = Tuple[int, int]
Split = Tuple[int, List[int], List[int]]


class Structure:
    """Base class: a merge structure for one variable."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:  # pragma: no cover
        raise NotImplementedError

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        """:meth:`permitted_merges` for a batch of partitions, as a boolean mask.

        The vectorised search (:mod:`catci.search`) stores a partition of ``d``
        labels by *slot*: a group lives at the slot of its smallest label.

        Parameters
        ----------
        sizes : ``(B, d)`` int, size of the group at each slot (0 = no group there).
        gid : ``(B, d)`` int, ``gid[b, l]`` = slot of the group holding label ``l + 1``.

        Returns
        -------
        ``(B, d, d)`` bool, ``True`` at ``[b, s, t]`` (``s < t``) when merging the
        groups at slots ``s`` and ``t`` is permitted. Slot order is partition-position
        order, so a row-major scan of the mask visits candidates in the order
        ``permitted_merges`` lists them whenever that order is lexicographic -- true
        of every structure in this module.

        This default rebuilds each partition and calls :meth:`permitted_merges`
        (memoised on the partition); subclasses override it with array code.
        """
        B, d = sizes.shape
        mask = np.zeros((B, d, d), dtype=bool)
        cache = self.__dict__.setdefault("_mask_cache", {})
        if len(cache) > 100_000:  # bound the memo; it only saves recomputation
            cache.clear()
        for b in range(B):
            key = gid[b].tobytes()
            if key not in cache:
                slots = np.flatnonzero(sizes[b])
                partition = [list(np.flatnonzero(gid[b] == s) + 1) for s in slots]
                pairs = self.permitted_merges(partition)
                cache[key] = (slots[[i - 1 for i, _ in pairs]], slots[[j - 1 for _, j in pairs]])
            s, t = cache[key]
            mask[b, s, t] = True
        return mask

    def coarsest_partitions(self, d: int) -> List[Partition]:  # pragma: no cover
        raise NotImplementedError

    def permitted_splits(self, partition: Partition) -> List[Split]:  # pragma: no cover
        raise NotImplementedError


def _guard(sizes: np.ndarray) -> np.ndarray:
    """``(B, 1, 1)`` mask of partitions with more than two groups (the ``d > 2`` guard)."""
    return (np.count_nonzero(sizes, axis=1) > 2)[:, None, None]


class Ordinal(Structure):
    """Adjacent merges only: (1,2), (2,3), ..., (d-1, d)."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(i, i + 1) for i in range(1, d)]

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        active = sizes > 0
        pos = np.cumsum(active, axis=1)
        nxt = pos[:, None, :] == pos[:, :, None] + 1
        return nxt & active[:, :, None] & active[:, None, :] & _guard(sizes)

    def coarsest_partitions(self, d: int) -> List[Partition]:
        """Every contiguous two-group partition ``{1..c} | {c+1..d}``."""
        if d < 2:
            raise ValueError("Need d >= 2.")
        return [[list(range(1, c + 1)), list(range(c + 1, d + 1))] for c in range(1, d)]

    def permitted_splits(self, partition: Partition) -> List[Split]:
        """Cut each group at every interior position; groups stay intervals."""
        out: List[Split] = []
        for i, group in enumerate(partition, start=1):
            g = list(group)
            out.extend((i, g[:c], g[c:]) for c in range(1, len(g)))
        return out


class Saturated(Structure):
    """All pairs (the 'greedy' search in the R package).

    Has no divisive counterpart: splitting a group of size ``s`` into two admits
    ``2^(s-1) - 1`` bipartitions, so the top-down candidate set is exponential
    where the bottom-up one is quadratic. This asymmetry is the whole reason
    ``divisive_search`` is restricted to Ordinal and Tree.
    """

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(i, j) for i in range(1, d) for j in range(i + 1, d + 1)]

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        active = sizes > 0
        d = sizes.shape[1]
        upper = np.triu(np.ones((d, d), dtype=bool), k=1)
        return upper & active[:, :, None] & active[:, None, :] & _guard(sizes)

    def coarsest_partitions(self, d: int) -> List[Partition]:
        raise NotImplementedError(_SATURATED_HAS_NO_SPLITS)

    def permitted_splits(self, partition: Partition) -> List[Split]:
        raise NotImplementedError(_SATURATED_HAS_NO_SPLITS)


_SATURATED_HAS_NO_SPLITS = (
    "Saturated has no tractable divisive counterpart (2^(s-1) - 1 bipartitions "
    "per group of size s). Use Ordinal or Tree with divisive_search, or "
    "Saturated with greedy_search."
)


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


def _find_node(labels: Sequence[int], tree: TreeNode) -> TreeNode:
    """The node of ``tree`` whose leaf-set is exactly ``labels``."""
    target = set(labels)
    node = tree
    while set(node.root) != target:
        for child in node.children:
            if target <= set(child.root):
                node = child
                break
        else:
            raise ValueError(f"No node of the tree has leaf-set {sorted(target)}.")
    return node


class Tree(Structure):
    """Sibling merges only, per a fixed binary tree over the original labels."""

    def __init__(self, tree: TreeNode):
        self.tree = tree
        self._sibling_table = _binary_sibling_table(tree)

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

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        if self._sibling_table is None:  # n-ary tree: groups need not be nodes
            return super().permitted_mask(sizes, gid)
        s1, n1, s2, n2 = self._sibling_table
        B, d = sizes.shape
        # In a binary tree every group is a node, and a node is fixed by (min label,
        # size): the two children of an internal node are both whole groups exactly
        # when the slots of their min labels hold groups of their sizes.
        ok = (sizes[:, s1] == n1) & (sizes[:, s2] == n2)
        b, k = np.nonzero(ok)
        mask = np.zeros((B, d, d), dtype=bool)
        mask[b, s1[k], s2[k]] = True
        return mask & _guard(sizes)

    def coarsest_partitions(self, d: int) -> List[Partition]:
        """The root's children -- a tree names its own coarsest split, so there is one."""
        if not self.tree.children:
            raise ValueError("Tree has no children to split into.")
        return [[list(c.root) for c in self.tree.children]]

    def permitted_splits(self, partition: Partition) -> List[Split]:
        """Split each group into the children of the node it corresponds to."""
        out: List[Split] = []
        for i, group in enumerate(partition, start=1):
            children = _find_node(group, self.tree).children
            if len(children) == 2:
                a, b = children
                out.append((i, list(a.root), list(b.root)))
            elif children:  # pragma: no cover - make_binary_tree only makes binary nodes
                raise ValueError("Divisive tree search needs binary (or leaf) nodes.")
        return out


def _binary_sibling_table(tree: TreeNode):
    """``(s1, n1, s2, n2)`` per internal node: 0-based min label and size of each child.

    ``None`` unless every internal node has exactly two children.
    """
    rows = []
    stack = [tree]
    while stack:
        node = stack.pop()
        if not node.children:
            continue
        if len(node.children) != 2:
            return None
        c1, c2 = sorted(node.children, key=lambda c: min(c.root))
        rows.append((min(c1.root) - 1, len(c1.root), min(c2.root) - 1, len(c2.root)))
        stack.extend(node.children)
    return tuple(np.array(col, dtype=np.int64) for col in zip(*rows)) if rows else None
