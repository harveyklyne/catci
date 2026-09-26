"""Label merging and rank-one update formulae.

Ports the load-bearing update formulae (24)-(27) from the R package
(``R/merging_functions.R``). Every function here is a pure function of its
arguments and is pinned by the ``rank_one_updates`` fixture and, transitively,
by ``search_paths``.

``split_*`` invert the same three formulae for :func:`~catci.search.divisive_search`,
which walks the partition lattice the other way. Splitting a group is exactly
un-doing the merge that would put it back, so (24)-(27) rearrange to subtract
rather than add. What changes is the input they need: a merge reads the current
``(T, Sigma)``, whereas a split reads only the two *new* rows of the finer
``(T, Sigma)`` -- ``Ta``/``Tb`` and ``Sa``/``Sb`` below -- which
:mod:`catci.blocks` supplies without ever materialising the finer matrices.

Index convention (matches the R package and the fixtures): the length ``dx*dy``
vector is ordered ``(1,1),(2,1),...,(dx,1),(1,2),...,(dx,dy)`` -- X fastest.
Labels are 1-based on the public surface; positions are 0-based numpy indices.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "x_index",
    "y_index",
    "get_index",
    "update_T",
    "update_Sigma",
    "update_normsq",
    "update_tr",
    "update_tr2",
    "split_normsq",
    "split_tr",
    "split_tr2",
]


def x_index(j: int, dx: int, dy: int) -> np.ndarray:
    """Boolean mask (length ``dx*dy``) of positions whose X-label equals ``j``."""
    if not (1 <= j <= dx):
        raise ValueError("Need 1 <= j <= dx.")
    p = np.arange(dx * dy)
    return (p % dx) == (j - 1)


def y_index(k: int, dx: int, dy: int) -> np.ndarray:
    """Boolean mask (length ``dx*dy``) of positions whose Y-label equals ``k``."""
    if not (1 <= k <= dy):
        raise ValueError("Need 1 <= k <= dy.")
    p = np.arange(dx * dy)
    return (p // dx) == (k - 1)


def get_index(dimension: int, j: int, dx: int, dy: int) -> np.ndarray:
    """Mask for label ``j`` in ``dimension`` (1 = X, 2 = Y)."""
    if dimension == 1:
        return x_index(j, dx, dy)
    if dimension == 2:
        return y_index(k=j, dx=dx, dy=dy)
    raise ValueError("dimension must be 1 (X) or 2 (Y).")


def _check_balanced(index1: np.ndarray, index2: np.ndarray) -> None:
    if index1.sum() != index2.sum():
        raise ValueError("index1 and index2 must have the same number of True entries.")


def update_T(T_vector: np.ndarray, index1: np.ndarray, index2: np.ndarray) -> np.ndarray:
    """Merge the label pair: add ``index2`` entries into ``index1``, drop ``index2``."""
    _check_balanced(index1, index2)
    new_T = T_vector.copy()
    new_T[index1] = T_vector[index1] + T_vector[index2]
    return new_T[~index2]


def update_Sigma(Sigma: np.ndarray, index1: np.ndarray, index2: np.ndarray) -> np.ndarray:
    """Merge rows/cols ``index2`` into ``index1`` and drop ``index2``.

    Row update is applied before the column update, matching the R code.
    """
    _check_balanced(index1, index2)
    new_Sigma = Sigma.copy()
    new_Sigma[index1, :] = new_Sigma[index1, :] + new_Sigma[index2, :]
    new_Sigma[:, index1] = new_Sigma[:, index1] + new_Sigma[:, index2]
    keep = ~index2
    return new_Sigma[np.ix_(keep, keep)]


def update_normsq(normsq: float, T_vector: np.ndarray, index1: np.ndarray, index2: np.ndarray) -> float:
    """Update ``||T||^2`` after merging the label pair. Formula (24)."""
    _check_balanced(index1, index2)
    return normsq + 2.0 * np.sum(T_vector[index1] * T_vector[index2])


def update_tr(tr: float, Sigma: np.ndarray, index1: np.ndarray, index2: np.ndarray) -> float:
    """Update ``tr(Sigma)`` after merging the label pair. Formula (25)."""
    _check_balanced(index1, index2)
    i1 = np.flatnonzero(index1)
    i2 = np.flatnonzero(index2)
    # Sigma[i1, i2] pairs the k-th selected row with the k-th selected col -> the diagonal.
    return tr + 2.0 * np.sum(Sigma[i1, i2])


def update_tr2(tr2: float, Sigma: np.ndarray, index1: np.ndarray, index2: np.ndarray) -> float:
    """Update ``tr(Sigma^2)`` after merging the label pair. Formulae (26)-(27)."""
    _check_balanced(index1, index2)
    i1 = np.flatnonzero(index1)
    i2 = np.flatnonzero(index2)
    cross = np.sum(Sigma[:, i1] * Sigma[:, i2])
    block = np.sum(Sigma[np.ix_(i1, i1)] * Sigma[np.ix_(i2, i2)]) + np.sum(
        Sigma[np.ix_(i1, i2)] * Sigma[np.ix_(i2, i1)]
    )
    return tr2 + 4.0 * cross + 2.0 * block


# --------------------------------------------------------------------------- #
# Inverse formulae: refining a partition (see module docstring)
# --------------------------------------------------------------------------- #
# A split replaces one group by two, ``a`` then ``b``. The arguments describe
# the *finer* partition: ``Ta``, ``Tb`` are its entries for the new groups
# (length ``r``, one per group of the other dimension), and ``Sa``, ``Sb`` are
# the corresponding row blocks of its covariance, each ``(r, m)`` over all ``m``
# entries of the finer partition. ``cols_a`` and ``cols_b`` pick out the columns
# of those blocks belonging to ``a`` and to ``b``, in the same order as the rows.


def split_normsq(normsq: float, Ta: np.ndarray, Tb: np.ndarray) -> float:
    """Update ``||T||^2`` after splitting a group. Inverse of formula (24)."""
    return normsq - 2.0 * float(np.sum(Ta * Tb))


def split_tr(tr: float, Sa: np.ndarray, cols_b: np.ndarray) -> float:
    """Update ``tr(Sigma)`` after splitting a group. Inverse of formula (25)."""
    return tr - 2.0 * float(np.sum(Sa[np.arange(Sa.shape[0]), cols_b]))


def split_tr2(
    tr2: float,
    Sa: np.ndarray,
    Sb: np.ndarray,
    cols_a: np.ndarray,
    cols_b: np.ndarray,
) -> float:
    """Update ``tr(Sigma^2)`` after splitting a group. Inverse of formulae (26)-(27)."""
    cross = float(np.sum(Sa * Sb))
    Saa, Sab = Sa[:, cols_a], Sa[:, cols_b]
    Sba, Sbb = Sb[:, cols_a], Sb[:, cols_b]
    block = float(np.sum(Saa * Sbb) + np.sum(Sab * Sba))
    return tr2 - 4.0 * cross - 2.0 * block
