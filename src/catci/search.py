"""Greedy label-merging search (ports ``greedy_query``).

At each level, for each dimension with more than two groups, evaluate every
permitted merge (from that dimension's :class:`~catci.structure.Structure`),
take the merge that maximises the statistic, apply it, and repeat until both
dimensions have two groups. Returns the statistic value before any merge and
after each level, plus the partition sequence.

With ``colsample_bylevel = 1`` (the only mode the paper uses) this is a
deterministic function of ``(T, Sigma)`` -- so it is pinned exactly by the
``search_paths`` fixture. Candidate evaluation order and first-max tie-breaking
match the R ``merge`` loop (dimension 1 then 2; within a dimension the order the
structure returns), so the selected path is identical to R's.

Alternatives to the greedy argmax (TODO item 7), all over the same candidate set:

* :func:`beam_search` keeps the top ``width`` partition pairs per level instead of
  one; ``width = 1`` is :func:`greedy_search` exactly.
* :func:`random_merges` draws a merge path without looking at ``T`` at all, and
  :func:`evaluate_path` scores a given merge path on any ``(T, Sigma)``. Together
  they give the random search, and -- with the path taken from
  ``greedy_search(...).merges`` on one half of the data -- the sample-split search.

Each is a map ``(T, Sigma) -> `` statistic path (randomised only through draws
independent of ``T``), so :mod:`catci.calibrate` covers them as-is.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from . import merging
from .statistic import ApproxChi
from .structure import Structure


Merge = Tuple[int, int, int]  # (dimension, i, j): 1 = X / 2 = Y, 1-based positions i < j


@dataclass
class SearchResult:
    values: List[float] = field(default_factory=list)
    partitions: List[dict] = field(default_factory=list)
    merges: List[Merge] = field(default_factory=list)  # the merge applied at each level


def _copy_partition(partition: dict) -> dict:
    return {k: [list(g) for g in groups] for k, groups in partition.items()}


def greedy_search(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    statistic=None,
) -> SearchResult:
    """Run the greedy merge search; see module docstring."""
    if statistic is None:
        statistic = ApproxChi()

    T_vector = np.array(T_vector, dtype=float)
    Sigma = np.array(Sigma, dtype=float)

    structures = {1: x_structure, 2: y_structure}
    partition = {
        "x": [[j] for j in range(1, dx + 1)],
        "y": [[k] for k in range(1, dy + 1)],
    }
    key = {1: "x", 2: "y"}
    dims = {1: dx, 2: dy}

    result = SearchResult()
    result.values.append(statistic.value(statistic.init(T_vector, Sigma)))
    result.partitions.append(_copy_partition(partition))

    while dims[1] > 2 or dims[2] > 2:
        base_state = statistic.init(T_vector, Sigma)

        best = None  # (value, dimension, i, j, index1, index2)
        for dimension in (1, 2):
            groups = partition[key[dimension]]
            for (i, j) in structures[dimension].permitted_merges(groups):
                index1 = merging.get_index(dimension, i, dims[1], dims[2])
                index2 = merging.get_index(dimension, j, dims[1], dims[2])
                state = statistic.update(base_state, T_vector, Sigma, index1, index2)
                value = statistic.value(state)
                # strict '>' keeps the first candidate in loop order on ties (R which.max).
                if best is None or value > best[0]:
                    best = (value, dimension, i, j, index1, index2)

        if best is None:
            break  # no permitted merge anywhere (both dimensions guarded)

        value, dimension, i, j, index1, index2 = best
        result.values.append(value)
        result.merges.append((dimension, i, j))

        # Apply the winning merge to the partition and to (T, Sigma).
        groups = partition[key[dimension]]
        groups[i - 1] = groups[i - 1] + groups[j - 1]
        del groups[j - 1]
        result.partitions.append(_copy_partition(partition))

        dims[dimension] -= 1
        T_vector = merging.update_T(T_vector, index1, index2)
        Sigma = merging.update_Sigma(Sigma, index1, index2)

    return result


# --------------------------------------------------------------------------- #
# Alternative searches (TODO item 7)
# --------------------------------------------------------------------------- #
def _initial_partition(dx: int, dy: int) -> dict:
    return {"x": [[j] for j in range(1, dx + 1)], "y": [[k] for k in range(1, dy + 1)]}


def _candidates(partition: dict, x_structure: Structure, y_structure: Structure) -> List[Merge]:
    """Every permitted merge, in greedy's loop order (X then Y, structure order)."""
    return [(1, i, j) for (i, j) in x_structure.permitted_merges(partition["x"])] + [
        (2, i, j) for (i, j) in y_structure.permitted_merges(partition["y"])
    ]


def _apply(partition: dict, dims: dict, merge: Merge) -> Tuple[dict, dict]:
    """The partition and group counts after ``merge`` (inputs are not modified)."""
    dimension, i, j = merge
    out = _copy_partition(partition)
    groups = out["x" if dimension == 1 else "y"]
    groups[i - 1] = groups[i - 1] + groups[j - 1]
    del groups[j - 1]
    new_dims = dict(dims)
    new_dims[dimension] -= 1
    return out, new_dims


def _canonical(partition: dict) -> tuple:
    """Hashable key for a partition pair; merge order does not change (T, Sigma)."""
    return tuple(tuple(sorted(tuple(sorted(g)) for g in partition[k])) for k in ("x", "y"))


def evaluate_path(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    merges: List[Merge],
    statistic=None,
) -> SearchResult:
    """Score a *given* merge path on ``(T, Sigma)``: ``len(merges) + 1`` values.

    No search happens -- the path is fixed in advance (drawn by
    :func:`random_merges`, or chosen on independent data), so this costs one
    update per level. ``evaluate_path(T, S, dx, dy, greedy_search(T, S, ...).merges)``
    reproduces the greedy values exactly.
    """
    if statistic is None:
        statistic = ApproxChi()
    T_vector = np.array(T_vector, dtype=float)
    Sigma = np.array(Sigma, dtype=float)
    dims = {1: dx, 2: dy}
    partition = _initial_partition(dx, dy)

    result = SearchResult()
    state = statistic.init(T_vector, Sigma)
    result.values.append(statistic.value(state))
    result.partitions.append(_copy_partition(partition))
    for merge in merges:
        dimension, i, j = merge
        index1 = merging.get_index(dimension, i, dims[1], dims[2])
        index2 = merging.get_index(dimension, j, dims[1], dims[2])
        state = statistic.update(state, T_vector, Sigma, index1, index2)
        result.values.append(statistic.value(state))
        result.merges.append(merge)
        partition, dims = _apply(partition, dims, merge)
        result.partitions.append(partition)
        T_vector = merging.update_T(T_vector, index1, index2)
        Sigma = merging.update_Sigma(Sigma, index1, index2)
    return result


def random_merges(
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    rng: np.random.Generator,
) -> List[Merge]:
    """A merge path drawn uniformly at random, without looking at the data.

    At each level one merge is drawn uniformly from the candidate set the greedy
    search would score (both dimensions pooled), until both dimensions have two
    groups -- so the path has the same length as the greedy one.
    """
    dims = {1: dx, 2: dy}
    partition = _initial_partition(dx, dy)
    merges: List[Merge] = []
    while dims[1] > 2 or dims[2] > 2:
        candidates = _candidates(partition, x_structure, y_structure)
        if not candidates:
            break
        merge = candidates[rng.integers(len(candidates))]
        merges.append(merge)
        partition, dims = _apply(partition, dims, merge)
    return merges


@dataclass
class _BeamState:
    partition: dict
    dims: dict
    T_vector: np.ndarray
    Sigma: np.ndarray


def beam_search(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    width: int,
    statistic=None,
) -> SearchResult:
    """Beam search of width ``width`` over the greedy candidate set.

    Each level scores every permitted merge of every state in the beam and keeps
    the ``width`` best *distinct* partition pairs (two merge orders reaching the
    same pair give the same ``(T, Sigma)``, so duplicates are dropped). The value
    at each depth is the best found at that depth; ``partitions`` / ``merges``
    record the state attaining it, so for ``width > 1`` they need not form one
    nested chain.

    Ties are broken by beam position, then candidate order, which makes
    ``width = 1`` identical to :func:`greedy_search` (first max wins). Cost is
    about ``width`` times the greedy search.
    """
    if width < 1:
        raise ValueError("width must be >= 1.")
    if statistic is None:
        statistic = ApproxChi()

    T_vector = np.array(T_vector, dtype=float)
    Sigma = np.array(Sigma, dtype=float)
    partition = _initial_partition(dx, dy)
    beam = [_BeamState(partition, {1: dx, 2: dy}, T_vector, Sigma)]

    result = SearchResult()
    result.values.append(statistic.value(statistic.init(T_vector, Sigma)))
    result.partitions.append(_copy_partition(partition))

    while any(s.dims[1] > 2 or s.dims[2] > 2 for s in beam):
        scored = []  # (value, beam position, merge, index1, index2), in loop order
        for b, s in enumerate(beam):
            base_state = statistic.init(s.T_vector, s.Sigma)
            for merge in _candidates(s.partition, x_structure, y_structure):
                dimension, i, j = merge
                index1 = merging.get_index(dimension, i, s.dims[1], s.dims[2])
                index2 = merging.get_index(dimension, j, s.dims[1], s.dims[2])
                state = statistic.update(base_state, s.T_vector, s.Sigma, index1, index2)
                scored.append((statistic.value(state), b, merge, index1, index2))
        if not scored:
            break

        scored.sort(key=lambda r: -r[0])  # stable: ties keep loop order
        new_beam: List[_BeamState] = []
        seen = set()
        for value, b, merge, index1, index2 in scored:
            s = beam[b]
            partition, dims = _apply(s.partition, s.dims, merge)
            key = _canonical(partition)
            if key in seen:
                continue
            seen.add(key)
            if not new_beam:  # the best at this depth
                result.values.append(value)
                result.merges.append(merge)
                result.partitions.append(_copy_partition(partition))
            new_beam.append(
                _BeamState(
                    partition,
                    dims,
                    merging.update_T(s.T_vector, index1, index2),
                    merging.update_Sigma(s.Sigma, index1, index2),
                )
            )
            if len(new_beam) == width:
                break
        beam = new_beam

    return result
