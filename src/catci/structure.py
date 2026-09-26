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
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "Structure",
    "Ordinal",
    "Saturated",
    "Cyclic",
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


class Cyclic(Structure):
    """Adjacent merges on a cycle: the ordinal pairs plus the wrap pair (1, d).

    For temporal categoricals (day-of-week, month, hour). Groups stay arcs of the
    cycle, in cyclic order by position: a merge lands at the lower position, so an
    arc through the wrap point sits at position 1 and ``(1, d)`` still joins the
    two groups either side of it. Pinned by ``test_cyclic_groups_stay_arcs``.
    """

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(i, i + 1) for i in range(1, d)] + [(1, d)]


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


def _check_tree(node: TreeNode) -> None:
    """Children must partition their parent's root; no unary nodes."""
    if not node.children:
        if len(node.root) != 1:
            raise ValueError(f"Leaf {node.root} must hold exactly one label.")
        return
    if len(node.children) < 2:
        raise ValueError(f"Node {node.root} has a single child; collapse it.")
    labels = [lab for c in node.children for lab in c.root]
    if len(labels) != len(set(labels)) or tuple(sorted(labels)) != node.root:
        raise ValueError(f"Children of {node.root} do not partition it.")
    for c in node.children:
        _check_tree(c)


class Tree(Structure):
    """Merges within one node of a fixed tree over the original labels.

    The tree may have any arity. Every group in a partition reachable by the
    search is either a whole node, or a union of two or more (but not all)
    children of one node. Call that node the group's *context*: the parent of the
    node the group equals, or else the node whose children it unions. Two groups
    may merge iff they share a context. For a binary tree this is exactly
    "merge siblings" (a union of both children is the node itself). For an n-ary
    node it lets a partially merged block of children keep absorbing the rest --
    matching only whole sibling nodes would strand ``{a, b}`` next to ``c``.
    """

    def __init__(self, tree: TreeNode):
        _check_tree(tree)
        self.tree = tree
        self._parent = {}  # node root -> parent node
        self._nodes = []   # every internal node, pre-order
        stack = [tree]
        while stack:
            node = stack.pop()
            if node.children:
                self._nodes.append(node)
            for c in node.children:
                self._parent[c.root] = node
                stack.append(c)
        self._contexts = {}  # group -> context; groups recur across bootstrap draws

    @classmethod
    def binary(cls, d: int) -> "Tree":
        return cls(make_binary_tree(d))

    @classmethod
    def from_parents(
        cls,
        parents: Mapping[Hashable, Optional[Hashable]],
        levels: Sequence[Hashable],
    ) -> "Tree":
        """Build a tree from a taxonomy given as a code -> parent-code mapping.

        ``levels[k]`` is the category encoded as label ``k + 1``: each must be a
        key of ``parents``, and none may be an ancestor of another. Codes that are
        not keys, or map to ``None``, are roots; several roots are joined under
        one synthetic root. Internal codes with no level beneath them are pruned
        and single-child chains are collapsed, so a taxonomy can be passed whole
        even when only some of its leaves occur in the data. Because internal
        and leaf codes share one namespace, taxonomies whose levels reuse codes
        (e.g. numeric CCS and chapter ids) should key on ``(level, code)``.
        """
        levels = list(levels)
        if len(set(levels)) != len(levels):
            raise ValueError("levels must be distinct.")
        missing = [lev for lev in levels if lev not in parents]
        if missing:
            raise ValueError(f"levels with no entry in parents: {missing[:10]}")

        children: dict = {}
        for code in levels:  # walk each level to its root; only touched codes matter
            seen = {code}
            node = code
            while parents.get(node) is not None:
                parent = parents[node]
                if parent in seen:
                    raise ValueError(f"Cycle in parents through {parent!r}.")
                seen.add(parent)
                kids = children.setdefault(parent, [])
                if node not in kids:
                    kids.append(node)
                else:
                    break  # the rest of this path was already walked
                node = parent
            else:
                children.setdefault(None, [])
                if node not in children[None]:
                    children[None].append(node)
        label = {lev: k + 1 for k, lev in enumerate(levels)}
        clash = [lev for lev in levels if lev in children]
        if clash:
            raise ValueError(f"levels that are ancestors of other levels: {clash[:10]}")

        def build(code) -> TreeNode:
            if code in label:
                return _make_leaf(label[code])
            kids = [build(c) for c in children[code]]
            return kids[0] if len(kids) == 1 else _make_parent(kids)

        return cls(build(None))

    def _context(self, group: Sequence[int]) -> TreeNode:
        g = tuple(sorted(group))
        if g in self._parent:
            return self._parent[g]
        if g in self._contexts:
            return self._contexts[g]
        members = set(g)
        for node in self._nodes:
            if members <= set(node.root) and all(
                set(c.root) <= members or not (set(c.root) & members) for c in node.children
            ):
                self._contexts[g] = node
                return node
        raise ValueError(f"Group {list(group)} is not a union of children of any tree node.")

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        by_context: dict = {}  # bucket positions by context: O(d + pairs), not O(d^2)
        for pos, group in enumerate(partition, start=1):
            by_context.setdefault(id(self._context(group)), []).append(pos)
        return sorted(
            (i, j) for block in by_context.values()
            for a, i in enumerate(block) for j in block[a + 1:]
        )
