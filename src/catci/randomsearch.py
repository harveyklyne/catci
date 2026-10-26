"""Random-chooser baselines for the merge and split searches.

The directed searches pick each next merge/split by the criterion. These baselines
spend the same depth budget on **uniformly random** steps, so comparing the two
isolates whether choosing steps by the data helps, as opposed to the family of
nested partitions and the minP correction alone.

A search here is a *generator* of partition paths plus a *scorer* of those paths:

* :func:`generate_merge_paths` / :func:`generate_split_paths` draw, for every column
  of ``T`` (the observed statistic and each bootstrap draw), a path of partitions
  uniformly among what the structures permit, x and y candidates pooled. They never
  look at ``T``, ``Sigma`` or the statistic.
* :func:`score_paths` evaluates a statistic at every partition of every path.

Because a path is drawn independently of ``T`` and afresh for each column, the map
``T -> path values`` is the same randomised map for every column, and
:mod:`catci.calibrate` applies unchanged (its guarantee needs only exchangeable
columns). One full path is drawn and truncated at ``max_levels``, so the first
``k + 1`` entries do not depend on ``max_levels``.

Extension point
---------------
The calibrator consumes an ``(S, B + 1)`` matrix of slot values and needs only
exchangeable columns. A depth-slotted path (``S`` = depth) is the special case
implemented here. A generator may return any set of candidates per column provided

1. it depends only on ``(T, Sigma_hat, structures)`` and randomness independent of
   ``T`` -- never the raw rows, which the bootstrap draws do not have;
2. it is the same map for the observed column and every bootstrap column;
3. each row ("slot") has a column-independent meaning: a fixed named family, or an
   order-statistic rule such as "the best candidate with exactly ``g`` groups". "The
   ``l``-th candidate found" in arbitrary generation order is *not* a valid slot;
4. every column returns the same ``S`` (a missing slot needs an explicit rule, e.g.
   no evidence).

To add one, write a generator that returns a :class:`PathSet`-like object (or a
scorer's output directly) under those rules and give it a ``.paths`` method with the
signature of :class:`~catci.search.MergeSearch`. This reference scorer is
statistic-agnostic; a faster scorer must agree with it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

from .statistic import ApproxChi, approx_chi_array
from .structure import Structure

__all__ = [
    "PathSet",
    "RandomMergeSearch",
    "RandomSplitSearch",
    "generate_merge_paths",
    "generate_split_paths",
    "score_paths",
]

Groups = Tuple[Tuple[int, ...], ...]  # groups of 1-based labels, each sorted, ordered by least label
Key = Tuple[Groups, Groups]  # (x groups, y groups)


@dataclass
class PathSet:
    """``keys[node_of[l, b]]`` is the partition column ``b`` visits at depth ``l``."""

    keys: List[Key]
    node_of: np.ndarray  # (L + 1, B) ints

    @property
    def depth(self) -> int:
        return self.node_of.shape[0] - 1


# --------------------------------------------------------------------------- #
# Generator
# --------------------------------------------------------------------------- #
def _canon(groups) -> Groups:
    """Sorted labels within groups, groups ordered by least label."""
    return tuple(sorted((tuple(sorted(int(x) for x in g)) for g in groups), key=lambda g: g[0]))


def _lists(groups: Groups) -> List[List[int]]:
    return [list(g) for g in groups]


def _merge_children(key: Key, x_structure: Structure, y_structure: Structure) -> List[Key]:
    """Every merge the structures permit, x then y -- the directed candidate set."""
    px, py = key
    out: List[Key] = []
    for i, j in x_structure.permitted_merges(_lists(px)):
        out.append((_canon(px[:i - 1] + (px[i - 1] + px[j - 1],) + px[i:j - 1] + px[j:]), py))
    for i, j in y_structure.permitted_merges(_lists(py)):
        out.append((px, _canon(py[:i - 1] + (py[i - 1] + py[j - 1],) + py[i:j - 1] + py[j:])))
    return out


def _split_children(key: Key, x_structure: Structure, y_structure: Structure) -> List[Key]:
    """Every split of every group of either dimension."""
    px, py = key
    out: List[Key] = []
    for i, a, b in x_structure.permitted_splits(_lists(px)):
        out.append((_canon(px[:i - 1] + (tuple(a), tuple(b)) + px[i:]), py))
    for i, a, b in y_structure.permitted_splits(_lists(py)):
        out.append((px, _canon(py[:i - 1] + (tuple(a), tuple(b)) + py[i:])))
    return out


def _generate(start_keys: Sequence[Key], children, B: int, rng: np.random.Generator,
              max_levels: int | None) -> PathSet:
    """Level-synchronous uniform random walk of ``B`` columns down the candidate graph.

    Columns sitting on the same partition share one enumeration of its children and
    one vectorised draw. Nodes are visited in id order and columns in index order, so
    the stream consumed is a function of the (data-free) walk alone: truncating at
    ``max_levels`` stops consuming it and leaves the earlier levels unchanged.
    """
    keys: List[Key] = []
    index: dict = {}

    def intern(key: Key) -> int:
        if key not in index:
            index[key] = len(keys)
            keys.append(key)
        return index[key]

    start_ids = np.array([intern(k) for k in start_keys])
    cur = start_ids[rng.integers(len(start_ids), size=B)]
    levels = [cur]
    kids_of: dict = {}
    while max_levels is None or len(levels) - 1 < max_levels:
        nxt = np.empty(B, dtype=np.int64)
        finished = 0
        for node in np.unique(cur):  # ascending ids
            cols = np.flatnonzero(cur == node)
            if node not in kids_of:
                kids_of[node] = children(keys[node])
            kids = kids_of[node]
            if not kids:
                finished += cols.size
                continue
            pick = rng.integers(len(kids), size=cols.size)
            chosen = {c: intern(kids[c]) for c in np.unique(pick).tolist()}
            nxt[cols] = [chosen[c] for c in pick.tolist()]
        if finished:
            if finished != B:
                raise ValueError("Random paths have different lengths across draws.")
            break
        cur = nxt
        levels.append(cur)
    return PathSet(keys=keys, node_of=np.vstack(levels))


def generate_merge_paths(dx: int, dy: int, x_structure: Structure, y_structure: Structure,
                         B: int, rng: np.random.Generator,
                         max_levels: int | None = None) -> PathSet:
    """``B`` random merge paths from the finest partition, ``max_levels`` merges at most.

    Untruncated, each path runs until both dimensions have two groups, the length of
    the greedy path (``dx + dy - 3`` entries).
    """
    start = (tuple((j,) for j in range(1, dx + 1)), tuple((k,) for k in range(1, dy + 1)))
    return _generate([start], lambda k: _merge_children(k, x_structure, y_structure), B, rng,
                     max_levels)


def generate_split_paths(dx: int, dy: int, x_structure: Structure, y_structure: Structure,
                         B: int, rng: np.random.Generator,
                         max_levels: int | None = None) -> PathSet:
    """``B`` random split paths from a random coarsest partition, ``max_levels`` splits at most.

    The directed search chooses its starting pair among every pair of coarsest
    partitions the structures permit; the baseline draws that pair uniformly. Needs
    ``Ordinal`` or binary ``Tree`` on both sides (``Saturated`` raises, as in
    :class:`~catci.search.SplitSearch`).
    """
    starts = [(_canon(px), _canon(py))
              for px in x_structure.coarsest_partitions(dx)
              for py in y_structure.coarsest_partitions(dy)]
    return _generate(starts, lambda k: _split_children(k, x_structure, y_structure), B, rng,
                     max_levels)


# --------------------------------------------------------------------------- #
# Scorer (reference implementation)
# --------------------------------------------------------------------------- #
def _aggregation(key: Key, dx: int, dy: int) -> np.ndarray:
    """``Pi``: sums the X-fastest ``dx*dy`` vector over the cells of a partition."""
    px, py = key
    gx = np.empty(dx, dtype=np.int64)
    for g, group in enumerate(px):
        gx[np.asarray(group) - 1] = g
    gy = np.empty(dy, dtype=np.int64)
    for g, group in enumerate(py):
        gy[np.asarray(group) - 1] = g
    cell = np.repeat(gy, dx) * len(px) + np.tile(gx, dy)
    Pi = np.zeros((len(px) * len(py), dx * dy))
    Pi[cell, np.arange(dx * dy)] = 1.0
    return Pi


def score_paths(T: np.ndarray, Sigma: np.ndarray, paths: PathSet, dx: int, dy: int,
                statistic=None) -> np.ndarray:
    """``(L + 1, B)`` statistic values along ``paths``, for every column of ``T``.

    Reference implementation: each distinct partition is scored once. For
    :class:`~catci.statistic.ApproxChi` the ``Sigma``-only traces are computed once
    per partition and ``||Pi T||^2`` for all its columns at once; any other statistic
    goes through its ``init`` / ``value`` interface one column at a time.
    """
    statistic = ApproxChi() if statistic is None else statistic
    T = np.asarray(T, dtype=float)
    if T.ndim == 1:
        T = T[:, None]
    L1, B = paths.node_of.shape
    if T.shape[1] != B:
        raise ValueError("paths were generated for a different number of columns.")
    flat = paths.node_of.ravel()
    order = np.argsort(flat, kind="stable")
    uniq, start = np.unique(flat[order], return_index=True)
    bounds = np.append(start, flat.size)
    out = np.empty(flat.size)
    fast = isinstance(statistic, ApproxChi)
    for u, node in enumerate(uniq):
        pos = order[bounds[u]:bounds[u + 1]]
        Pi = _aggregation(paths.keys[node], dx, dy)
        S = Pi @ Sigma @ Pi.T
        PT = Pi @ T[:, pos % B]
        if fast:
            out[pos] = approx_chi_array((PT ** 2).sum(axis=0), np.trace(S), np.sum(S ** 2))
        else:
            out[pos] = [statistic.value(statistic.init(PT[:, m], S)) for m in range(PT.shape[1])]
    return out.reshape(L1, B)


# --------------------------------------------------------------------------- #
# Search objects (same interface as MergeSearch / SplitSearch)
# --------------------------------------------------------------------------- #
class _RandomSearch:
    _generate = None  # set by subclasses

    def __init__(self, max_levels: int | None = None, seed=None):
        """``seed`` may be an int or a ``SeedSequence``. Do not pass a ``Generator``:
        ``default_rng(generator)`` returns that same generator, so two instances would
        continue one stream rather than replay it. Each :meth:`paths` call advances the
        stream; build a fresh instance per replicate (and per learner arm)."""
        self.max_levels = max_levels
        self.rng = np.random.default_rng(seed)

    def result(self, *args, **kwargs):
        raise NotImplementedError(
            f"{type(self).__name__} reports statistic paths only (no single observed search).")

    def paths(self, T, Sigma, dx, dy, x_structure, y_structure, statistic=None, n_jobs=1):
        """``(L + 1, B)`` values for every column of ``T`` along fresh random paths.
        ``n_jobs`` is accepted for symmetry with :class:`~catci.search.MergeSearch`."""
        T = np.asarray(T, dtype=float)
        if T.ndim == 1:
            T = T[:, None]
        if T.shape[0] != dx * dy:
            raise ValueError("T must have dx*dy rows.")
        pathset = self._generate(dx, dy, x_structure, y_structure, T.shape[1], self.rng,
                                 self.max_levels)
        return score_paths(T, Sigma, pathset, dx, dy, statistic)


class RandomMergeSearch(_RandomSearch):
    """Uniformly random merges, from the finest partition (a baseline for ``MergeSearch``)."""

    _generate = staticmethod(generate_merge_paths)


class RandomSplitSearch(_RandomSearch):
    """Uniformly random splits, from a random coarsest partition (a baseline for ``SplitSearch``)."""

    _generate = staticmethod(generate_split_paths)
