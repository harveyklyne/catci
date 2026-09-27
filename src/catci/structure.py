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
from typing import Hashable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "Split",
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
Split = Tuple[int, List[int], List[int]]


class Structure:
    """Base class: a merge structure for one variable."""

    def permitted_merges(self, partition: Partition) -> List[Pair]:  # pragma: no cover
        raise NotImplementedError

    def check_labels(self, d: int) -> None:
        """Raise if this structure cannot describe a variable with ``d`` labels."""

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


class Cyclic(Structure):
    """Adjacent merges on a cycle: the ordinal pairs plus the wrap pair (1, d).

    For temporal categoricals (day-of-week, month, hour). Groups stay arcs of the
    cycle, in cyclic order by position: a merge lands at the lower position, so an
    arc through the wrap point sits at position 1 and ``(1, d)`` still joins the
    two groups either side of it. Pinned by ``test_cyclic_groups_stay_arcs``.
    Pairs are listed lexicographically, so ``(1, d)`` comes second, matching the
    row-major scan of :meth:`permitted_mask` (and hence its tie-breaking).
    """

    def permitted_merges(self, partition: Partition) -> List[Pair]:
        d = len(partition)
        if d <= 2:
            return []
        return [(1, 2), (1, d)] + [(i, i + 1) for i in range(2, d)]

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        mask = Ordinal().permitted_mask(sizes, gid)
        active = sizes > 0
        first = np.argmax(active, axis=1)
        last = active.shape[1] - 1 - np.argmax(active[:, ::-1], axis=1)
        mask[np.arange(len(sizes)), first, last] = True
        return mask & _guard(sizes)


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
        self.d = len(tree.root)
        if tree.root != tuple(range(1, self.d + 1)):
            raise ValueError(f"Tree leaves must be the labels 1..{self.d}.")
        self.tree = tree
        self._sibling_table = _binary_sibling_table(tree)
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
        self._lca_table = None  # built on first n-ary permitted_mask call

    def check_labels(self, d: int) -> None:
        if d != self.d:
            raise ValueError(
                f"Tree has {self.d} labels but the variable has {d}. Every tree leaf "
                "must be a label that is coded in the data (dx = x.max()); build the "
                "tree from the observed levels only.")

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
        not keys, or map to ``None`` or NaN (a pandas missing value), are roots;
        several roots are joined under one synthetic root. Only the ancestors of
        ``levels`` are visited, so a taxonomy can be passed whole even when only
        some of its leaves occur in the data; single-child chains are collapsed.
        Pass only the levels that occur: the search needs every label in
        ``1..len(levels)`` to be coded (see :meth:`check_labels`). Because
        internal and leaf codes share one namespace, taxonomies whose levels
        reuse codes (e.g. numeric CCS and chapter ids) should key on
        ``(level, code)``.
        """
        levels = list(levels)
        if not levels:
            raise ValueError("levels must be non-empty.")
        if any(_is_root(lev) for lev in levels):
            raise ValueError("levels may not be None or NaN.")
        if len(set(levels)) != len(levels):
            raise ValueError("levels must be distinct.")
        missing = [lev for lev in levels if lev not in parents]
        if missing:
            raise ValueError(f"levels with no entry in parents: {missing[:10]}")

        children: dict = {None: []}
        walked = set()  # codes already attached to their parent (each has one parent)
        for code in levels:  # walk each level up to its root
            seen = {code}
            node = code
            while node not in walked:
                walked.add(node)
                parent = parents.get(node)
                if _is_root(parent):
                    children[None].append(node)
                    break
                if parent in seen:
                    raise ValueError(f"Cycle in parents through {parent!r}.")
                seen.add(parent)
                children.setdefault(parent, []).append(node)
                node = parent
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
                if len(self._contexts) > 100_000:  # bound the memo; it only saves work
                    self._contexts.clear()
                self._contexts[g] = node
                return node
        raise ValueError(
            f"Group {[int(x) for x in group]} is not a union of children of any tree node.")

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

    def permitted_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        if self._sibling_table is None:  # n-ary tree: groups need not be nodes
            return self._context_mask(sizes, gid)
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

    def _context_mask(self, sizes: np.ndarray, gid: np.ndarray) -> np.ndarray:
        """The context rule as array code, for any arity.

        A group's lowest enclosing node is the LCA of its leaves of smallest and
        largest DFS rank. If the group *is* that node (same size) its context is
        the node's parent; otherwise it is a union of that node's children and the
        node is its context. Two slots may merge iff their contexts coincide.
        """
        if self._lca_table is None:
            self._lca_table = _lca_table(self.tree)
        rank, by_rank, lca, node_size, node_parent = self._lca_table
        B, d = sizes.shape
        rows = np.broadcast_to(np.arange(B)[:, None], (B, d))
        lo = np.full((B, d), d - 1, dtype=np.int64)
        hi = np.zeros((B, d), dtype=np.int64)
        np.minimum.at(lo, (rows, gid), rank)
        np.maximum.at(hi, (rows, gid), rank)
        node = lca[by_rank[lo], by_rank[hi]]
        ctx = np.where(node_size[node] == sizes, node_parent[node], node)
        active = sizes > 0
        upper = np.triu(np.ones((d, d), dtype=bool), k=1)
        same = ctx[:, :, None] == ctx[:, None, :]
        return same & upper & active[:, :, None] & active[:, None, :] & _guard(sizes)

    def coarsest_partitions(self, d: int) -> List[Partition]:
        """The root's children -- a tree names its own coarsest split, so there is one."""
        if len(self.tree.children) != 2:
            raise ValueError("Divisive tree search needs a binary tree (the root has "
                             f"{len(self.tree.children)} children).")
        return [[list(c.root) for c in self.tree.children]]

    def permitted_splits(self, partition: Partition) -> List[Split]:
        """Split each group into the children of the node it corresponds to."""
        out: List[Split] = []
        for i, group in enumerate(partition, start=1):
            children = _find_node(group, self.tree).children
            if len(children) == 2:
                a, b = children
                out.append((i, list(a.root), list(b.root)))
            elif children:  # an n-ary node, e.g. from Tree.from_parents
                raise ValueError("Divisive tree search needs binary (or leaf) nodes.")
        return out


def _is_root(parent) -> bool:
    """``None`` or NaN (``NaN != NaN``) marks a missing parent."""
    return parent is None or (isinstance(parent, float) and parent != parent)


def _lca_table(tree: TreeNode):
    """Lookup tables for :meth:`Tree._context_mask`, over 0-based labels.

    ``rank[l]`` is leaf ``l``'s DFS rank and ``by_rank`` its inverse; ``lca[a, b]``
    is the index of the lowest node holding leaves ``a`` and ``b``; ``node_size`` and
    ``node_parent`` are per node index (the root's parent is ``-1``).
    """
    d = len(tree.root)
    lca = np.zeros((d, d), dtype=np.int64)
    sizes, parents, order = [], [], []
    stack = [(tree, -1)]
    while stack:  # pre-order, children in order: descendants overwrite ancestors
        node, parent = stack.pop()
        k = len(sizes)
        sizes.append(len(node.root))
        parents.append(parent)
        idx = np.asarray(node.root) - 1
        lca[np.ix_(idx, idx)] = k
        if not node.children:
            order.append(node.root[0] - 1)
        stack.extend((c, k) for c in reversed(node.children))
    by_rank = np.asarray(order, dtype=np.int64)
    rank = np.empty(d, dtype=np.int64)
    rank[by_rank] = np.arange(d)
    return rank, by_rank, lca, np.asarray(sizes), np.asarray(parents)


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
