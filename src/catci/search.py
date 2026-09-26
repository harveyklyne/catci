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
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np

from . import merging
from .statistic import ApproxChi
from .structure import Structure


@dataclass
class SearchResult:
    values: List[float] = field(default_factory=list)
    partitions: List[dict] = field(default_factory=list)


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

    # Statistics that cannot be updated from summary quantities alone (ExactChi)
    # need to know *which* candidate they are scoring, to key a cache on. Opt-in,
    # so the ApproxChi path is exactly what it was.
    wants_context = getattr(statistic, "wants_context", False)
    if wants_context:
        statistic.begin_search(Sigma)

    result = SearchResult()
    result.values.append(statistic.value(statistic.init(T_vector, Sigma)))
    result.partitions.append(_copy_partition(partition))

    while dims[1] > 2 or dims[2] > 2:
        base_state = statistic.init(T_vector, Sigma)
        if wants_context:
            statistic.begin_level(partition)

        best = None  # (value, dimension, i, j, index1, index2)
        for dimension in (1, 2):
            groups = partition[key[dimension]]
            for (i, j) in structures[dimension].permitted_merges(groups):
                index1 = merging.get_index(dimension, i, dims[1], dims[2])
                index2 = merging.get_index(dimension, j, dims[1], dims[2])
                if wants_context:
                    state = statistic.update(
                        base_state, T_vector, Sigma, index1, index2, context=(dimension, i, j)
                    )
                else:
                    state = statistic.update(base_state, T_vector, Sigma, index1, index2)
                value = statistic.value(state)
                # strict '>' keeps the first candidate in loop order on ties (R which.max).
                if best is None or value > best[0]:
                    best = (value, dimension, i, j, index1, index2)

        if best is None:
            break  # no permitted merge anywhere (both dimensions guarded)

        value, dimension, i, j, index1, index2 = best
        result.values.append(value)

        # Apply the winning merge to the partition and to (T, Sigma).
        groups = partition[key[dimension]]
        groups[i - 1] = groups[i - 1] + groups[j - 1]
        del groups[j - 1]
        result.partitions.append(_copy_partition(partition))

        dims[dimension] -= 1
        T_vector = merging.update_T(T_vector, index1, index2)
        Sigma = merging.update_Sigma(Sigma, index1, index2)

    return result
