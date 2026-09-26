"""Rectangular block sums of ``(T, Sigma)`` over groups of labels.

:func:`~catci.search.greedy_search` never needs this: it starts from the
singleton partition, where ``(T, Sigma)`` *is* the input, and every later level
is a rank-one update of the level before. :func:`~catci.search.divisive_search`
starts from a two-group partition, where ``(Pi T, Pi Sigma Pi^T)`` is a sum over
large blocks of the input, and each level *refines* -- so there is no rank-one
route and the aggregation has to be paid somewhere.

Paying it per candidate costs ``O((dx*dy)^2)`` each time. Instead we cumulative-sum
``Sigma`` along all four label axes once, after which the sum over any
``(X-range x Y-range) x (X-range x Y-range)`` block is four array lookups. Both
structures that support divisive search keep groups contiguous in label order --
Ordinal by construction, Tree because ``make_binary_tree`` puts the leaves
in order -- so every group is a range and every block is a rectangle.

The split is deliberate: ``Sigma`` is shared by the observed statistic and all
``n_boot`` bootstrap draws, so :class:`SigmaBlocks` is built once per test and
reused, while the ``T`` prefix (:func:`t_prefix`) is rebuilt per draw for
``O(dx*dy)``.

Index convention matches the rest of the package: the length ``dx*dy`` vector is
ordered X-fastest, so it reshapes to ``(dy, dx)``. Ranges are 0-based half-open
``(lo, hi)``, and a *box* is the 4-vector ``(ylo, yhi, xlo, xhi)``.
"""

from __future__ import annotations

import numpy as np

__all__ = ["SigmaBlocks", "t_prefix", "t_blocks"]


class SigmaBlocks:
    """Block sums of ``Sigma`` over rectangular label ranges, in O(1) per block."""

    def __init__(self, Sigma: np.ndarray, dx: int, dy: int):
        C = np.asarray(Sigma, dtype=float).reshape(dy, dx, dy, dx)
        for axis in range(4):
            C = np.cumsum(C, axis=axis)
        # Pad a zero plane at the front of each axis so lookups can use the
        # exclusive lower end of a range without a special case at 0.
        self.cumulative = np.pad(C, [(1, 0)] * 4)
        self.dx = dx
        self.dy = dy

    def row_prefix(self, rows: np.ndarray) -> np.ndarray:
        """Collapse the row axes of each box; returns ``(R, dy + 1, dx + 1)``.

        The result is still cumulative in the *column* axes, so
        :meth:`column_blocks` finishes the job.
        """
        ylo, yhi, xlo, xhi = np.asarray(rows).T
        C = self.cumulative
        return C[yhi, xhi] - C[ylo, xhi] - C[yhi, xlo] + C[ylo, xlo]

    @staticmethod
    def column_blocks(row_prefix: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """Collapse the column axes of a :meth:`row_prefix`; returns ``(R, C)``."""
        ylo, yhi, xlo, xhi = np.asarray(cols).T
        return (
            row_prefix[:, yhi, xhi]
            - row_prefix[:, ylo, xhi]
            - row_prefix[:, yhi, xlo]
            + row_prefix[:, ylo, xlo]
        )

    def direct_blocks(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """``(R, C)`` block sums by one 16-term inclusion-exclusion.

        The alternative -- :meth:`row_prefix` then :meth:`column_blocks` -- takes
        8 gathers instead of 16, but one of its stages is ``(R, dy + 1, dx + 1)``
        regardless of how few columns are asked for. This form never touches
        anything wider than ``(R, C)``, which is what a coarse partition wants.
        """
        rylo, ryhi, rxlo, rxhi = np.asarray(rows).T
        cylo, cyhi, cxlo, cxhi = np.asarray(cols).T
        C = self.cumulative
        total = np.zeros((len(rylo), len(cylo)))
        for sr, ry, rx in ((1, ryhi, rxhi), (-1, rylo, rxhi), (-1, ryhi, rxlo), (1, rylo, rxlo)):
            for sc, cy, cx in ((1, cyhi, cxhi), (-1, cylo, cxhi), (-1, cyhi, cxlo), (1, cylo, cxlo)):
                total += (sr * sc) * C[ry[:, None], rx[:, None], cy[None, :], cx[None, :]]
        return total

    def blocks(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """``(R, C)`` matrix of block sums for row boxes ``rows``, column boxes ``cols``.

        Picks whichever of the two routes moves less memory: the two-stage one
        reads ``4R(dy + 1)(dx + 1) + 4RC`` entries, the direct one ``16RC``.
        """
        n_cols = len(cols)
        if 4 * (self.dy + 1) * (self.dx + 1) < 12 * n_cols:
            return self.column_blocks(self.row_prefix(rows), cols)
        return self.direct_blocks(rows, cols)


def t_prefix(T_vector: np.ndarray, dx: int, dy: int) -> np.ndarray:
    """Cumulative sum of ``T`` over both label axes, zero-padded at the front."""
    A = np.asarray(T_vector, dtype=float).reshape(dy, dx).cumsum(axis=0).cumsum(axis=1)
    return np.pad(A, [(1, 0), (1, 0)])


def t_blocks(prefix: np.ndarray, boxes: np.ndarray) -> np.ndarray:
    """Block sums of ``T`` for each box; returns ``(R,)``."""
    ylo, yhi, xlo, xhi = np.asarray(boxes).T
    return prefix[yhi, xhi] - prefix[ylo, xhi] - prefix[yhi, xlo] + prefix[ylo, xlo]
