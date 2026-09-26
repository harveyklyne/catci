"""Label-merging searches: agglomerative (:func:`greedy_search`, the paper's
Algorithm 1) and divisive (:func:`divisive_search`).

``greedy_search``: at each level, for each dimension with more than two groups,
evaluate every permitted merge (from that dimension's
:class:`~catci.structure.Structure`), take the merge that maximises the
statistic, apply it, and repeat until both dimensions have two groups. Returns
the statistic value before any merge and after each level, plus the partition
sequence.

With ``colsample_bylevel = 1`` (the only mode the paper uses) this is a
deterministic function of ``(T, Sigma)`` -- so it is pinned exactly by the
``search_paths`` fixture. Candidate evaluation order and first-max tie-breaking
match the R ``merge`` loop (dimension 1 then 2; within a dimension the order the
structure returns), so the selected path is identical to R's.

``divisive_search`` walks the same lattice in the opposite direction: it starts
fully collapsed, at two groups per dimension, and splits outwards. Run to
completion it visits the same ``dx + dy - 3`` levels, and it is likewise a
deterministic function of ``(T, Sigma)``, so the calibration in
:mod:`catci.calibrate` and the size guarantee behind it carry over unchanged --
the paper's Lemma 10 is stated for an arbitrary data-dependent map into
partition space, not for merging in particular.

What differs is where the work and the power sit. Merging spends its early --
most expensive -- levels near the singleton partition, where ``Sigma`` is
``dx*dy`` square, and reaches the coarse partitions only after a chain of
locally-greedy decisions taken on fine-level statistics. Splitting reaches the
coarse partitions first and picks each by its own statistic, and truncating it
with ``max_levels`` drops the fine levels, which are both the costly ones and
the ones carrying the most degrees of freedom.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from . import merging
from .blocks import SigmaBlocks, t_blocks, t_prefix
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
# Divisive search
# --------------------------------------------------------------------------- #
Range = Tuple[int, int]


def _as_range(group: List[int]) -> Range:
    """0-based half-open range of a group of 1-based labels; requires contiguity."""
    lo, hi = min(group), max(group)
    if hi - lo + 1 != len(group):
        raise ValueError(
            f"divisive_search needs groups contiguous in label order, got {sorted(group)}. "
            "Ordinal and Tree satisfy this; Saturated has no divisive counterpart."
        )
    return lo - 1, hi


def _labels(ranges: List[Range]) -> List[List[int]]:
    return [list(range(lo + 1, hi + 1)) for (lo, hi) in ranges]


def _boxes(x_ranges: List[Range], y_ranges: List[Range]) -> np.ndarray:
    """``(ylo, yhi, xlo, xhi)`` for every (X-group, Y-group) cell, X fastest."""
    xr = np.asarray(x_ranges, dtype=np.intp).reshape(-1, 2)
    yr = np.asarray(y_ranges, dtype=np.intp).reshape(-1, 2)
    out = np.empty((len(yr), len(xr), 4), dtype=np.intp)
    out[:, :, 0:2] = yr[:, None, :]
    out[:, :, 2:4] = xr[None, :, :]
    return out.reshape(-1, 4)


def _splits_of(structure: Structure, r: Range, memo: dict) -> List[Tuple[Range, Range]]:
    """The permitted splits of a single group, in range form, memoised.

    A group's splits depend only on that group (both Ordinal and Tree honour
    this -- see ``permitted_splits``), so the answer for a given range is fixed
    for the whole search and worth caching: without this, Tree re-walks the tree
    from the root for every group at every level.
    """
    if r not in memo:
        group = [list(range(r[0] + 1, r[1] + 1))]
        memo[r] = [(_as_range(a), _as_range(b)) for (_, a, b) in structure.permitted_splits(group)]
    return memo[r]


def divisive_search(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    statistic=None,
    max_levels: int | None = None,
    sigma_blocks: SigmaBlocks | None = None,
) -> SearchResult:
    """Run the greedy search top-down, from two groups per dimension outwards.

    ``max_levels`` caps the number of splits taken, so the path returned has at
    most ``max_levels + 1`` entries; ``None`` refines all the way to the
    singleton partition and returns the ``dx + dy - 3`` levels
    :func:`greedy_search` does.

    ``sigma_blocks`` lets a caller build :class:`~catci.blocks.SigmaBlocks` once
    and reuse it across bootstrap draws, which all share ``Sigma``.
    """
    if statistic is None:
        statistic = ApproxChi()
    if sigma_blocks is None:
        sigma_blocks = SigmaBlocks(Sigma, dx, dy)
    T_prefix = t_prefix(T_vector, dx, dy)

    result = SearchResult()

    # Level 0. Merging bottoms out at one partition, but splitting has to pick a
    # starting one, so score every pair the two structures permit. A tree offers
    # one; two ordinal variables offer (dx - 1)(dy - 1).
    best_start = None  # (value, state, x_ranges, y_ranges)
    for px in x_structure.coarsest_partitions(dx):
        xr = [_as_range(g) for g in px]
        for py in y_structure.coarsest_partitions(dy):
            yr = [_as_range(g) for g in py]
            boxes = _boxes(xr, yr)
            state = statistic.init(t_blocks(T_prefix, boxes), sigma_blocks.blocks(boxes, boxes))
            value = statistic.value(state)
            # strict '>' keeps the first candidate in loop order on ties, as in greedy_search.
            if best_start is None or value > best_start[0]:
                best_start = (value, state, xr, yr)

    value, state, x_ranges, y_ranges = best_start
    result.values.append(value)
    result.partitions.append({"x": _labels(x_ranges), "y": _labels(y_ranges)})

    memo: dict = {}
    level = 0
    while max_levels is None or level < max_levels:
        best = None  # (value, state, dimension, position, range_a, range_b)
        for dimension in (1, 2):
            structure = x_structure if dimension == 1 else y_structure
            ranges = x_ranges if dimension == 1 else y_ranges
            for i, r in enumerate(ranges):
                for (ra, rb) in _splits_of(structure, r, memo):
                    refined = ranges[:i] + [ra, rb] + ranges[i + 1:]
                    if dimension == 1:
                        xr, yr = refined, y_ranges
                        n_other = len(y_ranges)
                        # X is the fast axis, so a and b sit one apart within each row.
                        cols_a = np.arange(n_other) * len(refined) + i
                        cols_b = cols_a + 1
                        rows = _boxes([ra, rb], yr).reshape(n_other, 2, 4)
                        rows = np.ascontiguousarray(rows.transpose(1, 0, 2)).reshape(-1, 4)
                    else:
                        xr, yr = x_ranges, refined
                        n_other = len(x_ranges)
                        cols_a = np.arange(n_other) + i * n_other
                        cols_b = cols_a + n_other
                        rows = _boxes(xr, [ra, rb])

                    # Only the two new rows of the finer (T, Sigma), never the whole thing.
                    S = sigma_blocks.blocks(rows, _boxes(xr, yr))
                    T_ab = t_blocks(T_prefix, rows)
                    refined_state = statistic.split(
                        state,
                        T_ab[:n_other], T_ab[n_other:],
                        S[:n_other], S[n_other:],
                        cols_a, cols_b,
                    )
                    refined_value = statistic.value(refined_state)
                    if best is None or refined_value > best[0]:
                        best = (refined_value, refined_state, dimension, i, ra, rb)

        if best is None:
            break  # both dimensions are fully refined

        value, state, dimension, i, ra, rb = best
        if dimension == 1:
            x_ranges = x_ranges[:i] + [ra, rb] + x_ranges[i + 1:]
        else:
            y_ranges = y_ranges[:i] + [ra, rb] + y_ranges[i + 1:]
        result.values.append(value)
        result.partitions.append({"x": _labels(x_ranges), "y": _labels(y_ranges)})
        level += 1

    return result


# --------------------------------------------------------------------------- #
# Choosing a direction: what `calibrate` and `api` are handed
# --------------------------------------------------------------------------- #
# The bootstrap runs the same search over `n_boot + 1` vectors that all share
# one `Sigma`, so a search gets to `prepare` against that `Sigma` once and hand
# back a function of `T` alone. `MergeSearch` has nothing to set up;
# `SplitSearch` builds its block-sum table here, which is what makes truncated
# divisive search cheap -- the table is the only part that scales with dx * dy,
# and it is paid once per test rather than once per draw.


class MergeSearch:
    """Agglomerative direction: :func:`greedy_search`. The paper's Algorithm 1."""

    def prepare(self, Sigma, dx, dy, x_structure, y_structure, statistic):
        """Return ``path(T) -> SearchResult`` for this fixed ``Sigma``."""
        def path(T_vector: np.ndarray) -> SearchResult:
            return greedy_search(T_vector, Sigma, dx, dy, x_structure, y_structure, statistic)

        return path


class SplitSearch:
    """Divisive direction: :func:`divisive_search`, optionally truncated."""

    def __init__(self, max_levels: int | None = None):
        self.max_levels = max_levels

    def prepare(self, Sigma, dx, dy, x_structure, y_structure, statistic):
        """Return ``path(T) -> SearchResult``, with the block table built once."""
        blocks = SigmaBlocks(Sigma, dx, dy)

        def path(T_vector: np.ndarray) -> SearchResult:
            return divisive_search(
                T_vector, Sigma, dx, dy, x_structure, y_structure, statistic,
                max_levels=self.max_levels, sigma_blocks=blocks,
            )

        return path
