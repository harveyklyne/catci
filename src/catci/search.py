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

Implementation
--------------
:func:`greedy_search_paths` runs the search for many ``T`` draws at once against
one shared ``Sigma``, in lockstep over levels. It is the only implementation used
in practice; :func:`greedy_search` is its single-draw wrapper that also records
the partitions. Two ideas make it fast (see ``VECTORISATION.md``):

* **Every candidate at a level at once.** Formulae (24)-(27) are bilinear, so with
  the two merge labels left as free indices each quantity is a ``(d, d)`` matrix
  over all label pairs, and the Y branch is the same code on a transposed view.
* **Nothing is recomputed from scratch.** The expensive part of (26)-(27) --
  ``cross``, ``term1``, ``term2`` -- is carried across levels as ``(d, d)`` state
  and updated after each merge. A merge touches only two rows of the merged
  ``Sigma``, so each level costs ``O(p d)`` memory traffic plus one ``O(p^2)``-flop
  GEMM on ``O(p d)`` data, per draw, instead of ``O(p^2 d)`` reads.

Arrays are never shrunk: a merged group lives at the slot of its smallest
original label (which is also where the partition-list convention puts it), and
the absorbed label's rows are zeroed. Every draw therefore keeps the same shape
whichever path it takes, so draws batch with no padding logic, and slot order is
position order, so "first max over row-major upper triangles, X then Y" is
exactly the R tie-break for every built-in structure.

The per-candidate loop is kept as :func:`_greedy_search_loop`: it is the
reference the vectorised code is tested against, and the fallback for a
statistic other than :class:`~catci.statistic.ApproxChi`.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import List

import numpy as np

from . import merging
from .statistic import ApproxChi, approx_chi_array
from .structure import Structure

__all__ = ["SearchResult", "greedy_search", "greedy_search_paths"]

# Budget for one chunk's per-draw merged Sigmas (``8 p^2`` bytes each). Measured best
# at every d tried (8-32): larger chunks spill the per-merge temporaries out of cache.
DEFAULT_CHUNK_BYTES = 32 * 2 ** 20
# Drop absorbed labels from the padded arrays once they are at most this full.
COMPACT_THRESHOLD = 0.6


@dataclass
class SearchResult:
    values: List[float] = field(default_factory=list)
    partitions: List[dict] = field(default_factory=list)


def _copy_partition(partition: dict) -> dict:
    return {k: [list(g) for g in groups] for k, groups in partition.items()}


# --------------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------------- #
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
    if statistic is not None and not isinstance(statistic, ApproxChi):
        return _greedy_search_loop(T_vector, Sigma, dx, dy, x_structure, y_structure, statistic)

    T_vector = np.asarray(T_vector, dtype=float)
    shared = _Shared(Sigma, dx, dy)
    values, merges = _run_chunk(shared, T_vector[:, None], x_structure, y_structure)

    partition = {
        "x": [[j] for j in range(1, dx + 1)],
        "y": [[k] for k in range(1, dy + 1)],
    }
    result = SearchResult(values=[float(v) for v in values[:, 0]])
    result.partitions.append(_copy_partition(partition))
    for dim, i, j in merges:
        groups = partition["xy"[dim]]
        groups[i - 1] = groups[i - 1] + groups[j - 1]
        del groups[j - 1]
        result.partitions.append(_copy_partition(partition))
    return result


def greedy_search_paths(
    T: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    *,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    n_jobs: int = 1,
) -> np.ndarray:
    """Statistic paths of the greedy search for every column of ``T``.

    Parameters
    ----------
    T : ``(dx*dy, B)`` -- one statistic vector per column (a length-``dx*dy``
        vector is treated as ``B = 1``). All share ``Sigma``.
    chunk_bytes : memory budget for one chunk of draws; each draw holds its own
        merged ``Sigma``, ``8 (dx dy)^2`` bytes.
    n_jobs : threads over chunks (``-1`` = all cores). numpy releases the GIL in
        the heavy kernels, so this helps most once ``dx*dy`` is in the hundreds.

    Returns
    -------
    ``(L + 1, B)`` array; column ``b`` equals ``greedy_search(T[:, b], ...).values``.
    """
    T = np.asarray(T, dtype=float)
    if T.ndim == 1:
        T = T[:, None]
    p, B = T.shape
    if p != dx * dy:
        raise ValueError("T must have dx*dy rows.")
    shared = _Shared(Sigma, dx, dy)
    if n_jobs == -1:
        n_jobs = os.cpu_count() or 1
    chunk = int(max(1, min(-(-B // n_jobs), chunk_bytes // (8 * p * p))))
    slices = [slice(s, min(s + chunk, B)) for s in range(0, B, chunk)]

    def run(sl):
        return _run_chunk(shared, T[:, sl], x_structure, y_structure)[0]

    if n_jobs == 1 or len(slices) == 1:
        outs = [run(sl) for sl in slices]
    else:
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            outs = list(pool.map(run, slices))
    return np.concatenate(outs, axis=1)


# --------------------------------------------------------------------------- #
# Vectorised kernel
# --------------------------------------------------------------------------- #
#
# Notation. ``S[b, j, k, m, l] = Sigma_b[(j, k), (m, l)]`` with ``j, m`` X labels and
# ``k, l`` Y labels (``Sigma.reshape(dy, dx, dy, dx)`` transposed, since the vector
# is X-fastest); ``M[b, j, k]`` is the ``T`` table. For the X branch and a label
# pair ``(j, m)``, the deltas of formulae (24)-(27) are
#
#     d_normsq = 2 (M M^T)[j, m]
#     d_tr     = 2 sum_k S[j, k, m, k]
#     d_tr2    = 4 C[j, m] + 2 (T1[j, m] + T2[j, m])
#     C[j, m]  = sum_k (Sigma^2)[(j, k), (m, k)]        -- "cross"
#     T1[j, m] = <S[j, :, j, :], S[m, :, m, :]>         -- within-block
#     T2[j, m] = tr(S[j, :, m, :] S[j, :, m, :])        -- cross-block
#
# and the Y branch is the same with the roles of the axes swapped, i.e. on
# ``S.transpose(0, 2, 1, 4, 3)``. ``C, T1, T2`` and ``d_tr`` are the carried state.


class _Shared:
    """Draw-independent level-0 state, computed once and copied into each chunk."""

    def __init__(self, Sigma: np.ndarray, dx: int, dy: int):
        Sigma = np.asarray(Sigma, dtype=float)
        p = dx * dy
        if Sigma.shape != (p, p):
            raise ValueError("Sigma must be (dx*dy, dx*dy).")
        Sigma = (Sigma + Sigma.T) / 2.0  # the carried updates rely on exact symmetry
        self.dx, self.dy = dx, dy
        self.S = np.ascontiguousarray(Sigma.reshape(dy, dx, dy, dx).transpose(1, 0, 3, 2))
        self.tr = float(np.trace(Sigma))
        self.tr2 = float(np.sum(Sigma ** 2))
        self.C, self.T1, self.T2, self.dT, self.D = [], [], [], [], []
        for Sv in (self.S, self.S.transpose(1, 0, 3, 2)):
            d, e = Sv.shape[0], Sv.shape[1]
            flat = np.ascontiguousarray(Sv).reshape(d, e * p)
            self.C.append(flat @ flat.T)  # the one O(p^2 d) step, shared by all draws
            D = np.ascontiguousarray(np.einsum("jkjl->jkl", Sv))
            self.D.append(D)
            self.T1.append(D.reshape(d, e * e) @ D.reshape(d, e * e).T)
            self.T2.append(np.einsum("jkml,jlmk->jm", Sv, Sv))
            self.dT.append(2.0 * np.einsum("jkmk->jm", Sv))


class _State:
    """Per-chunk search state; every array has a leading draw axis of length ``n``.

    Label axes are *padded*: index ``i`` of dimension ``dim`` holds the group whose
    smallest original label is ``lab[dim][b, i]`` (0-based), and is inactive -- all
    its rows zero -- once that group has been absorbed. :meth:`compact` drops
    inactive indices when enough have accumulated; it preserves index order, which
    is partition-position order, so the tie-break is unaffected.

    ``C, T1, T2, dT`` are carried per dimension in one of two modes. **Dense**: the
    full symmetric ``(d, d)`` matrix is kept current, at one ``O(p^2)``-flop GEMM per
    merge -- right when most pairs are candidates (``Saturated``). **Sparse**: only
    the upper-triangle entries in ``valid`` are current, and ``valid`` is pruned to
    the permitted pairs every level; a permitted pair that is not valid is computed
    fresh from two rows of ``S``. Every pair a merge newly permits contains the
    merged group, so that is ``O(1)`` fresh pairs per level and each level costs
    ``O(p d)`` -- right when only ``O(d)`` pairs are ever candidates (trees, ordinal).
    """

    def __init__(self, shared: _Shared, T: np.ndarray, sparse):
        dx, dy = shared.dx, shared.dy
        n = T.shape[1]
        self.n = n
        self.sparse = tuple(sparse)
        self.S = np.broadcast_to(shared.S, (n,) + shared.S.shape).copy()
        self.M = np.ascontiguousarray(T.T.reshape(n, dy, dx).transpose(0, 2, 1))
        self.normsq = np.sum(T ** 2, axis=0)
        self.tr = np.full(n, shared.tr)
        self.tr2 = np.full(n, shared.tr2)

        def tile(mats):
            return [np.broadcast_to(c, (n,) + c.shape).copy() for c in mats]

        self.C, self.T1, self.T2, self.dT = (tile(shared.C), tile(shared.T1),
                                             tile(shared.T2), tile(shared.dT))
        self.D = tile(shared.D)  # diagonal blocks: D[dim][b, j] = S[j, :, j, :] in view(dim)
        # Level 0 is computed in full, so every entry starts valid.
        self.valid = [np.ones((n, d, d), dtype=bool) for d in (dx, dy)]
        # Partition bookkeeping, in original-label slots (what structures consume).
        self.sizes = [np.ones((n, d), dtype=np.int64) for d in (dx, dy)]
        self.gid = [np.broadcast_to(np.arange(d), (n, d)).copy() for d in (dx, dy)]
        self.lab = [np.broadcast_to(np.arange(d), (n, d)).copy() for d in (dx, dy)]
        self.compacted = [False, False]  # False while lab is still the identity

    def view(self, dim: int):
        """``(S, M)`` with the ``dim`` axis first -- the Y branch is a transposed view."""
        if dim == 0:
            return self.S, self.M
        return self.S.transpose(0, 2, 1, 4, 3), self.M.transpose(0, 2, 1)

    def active(self, dim: int) -> np.ndarray:
        return np.take_along_axis(self.sizes[dim], self.lab[dim], axis=1) > 0

    def compact(self, dim: int, threshold: float) -> None:
        """Drop inactive indices of ``dim`` if the padded size exceeds the need enough."""
        act = self.active(dim)
        d_pad = act.shape[1]
        need = int(act.sum(axis=1).max())
        if need > threshold * d_pad:
            return
        # Active indices first, in order; the tail pads with inactive (all-zero) ones.
        keep = np.argsort(~act, axis=1, kind="stable")[:, :need]
        self.compacted[dim] = True
        axes = (1, 3) if dim == 0 else (2, 4)
        self.S = _gather(_gather(self.S, keep, axes[0]), keep, axes[1])
        self.M = _gather(self.M, keep, 1 + dim)
        for mats in (self.C, self.T1, self.T2, self.dT, self.valid):
            mats[dim] = _gather(_gather(mats[dim], keep, 1), keep, 2)
        self.D[dim] = _gather(self.D[dim], keep, 1)
        self.D[1 - dim] = _gather(_gather(self.D[1 - dim], keep, 2), keep, 3)
        self.lab[dim] = _gather(self.lab[dim], keep, 1)

    def permitted(self, dim: int, structure: Structure) -> np.ndarray:
        """The structure's mask, mapped from original-label slots to padded indices."""
        mask = structure.permitted_mask(self.sizes[dim], self.gid[dim])
        if not self.compacted[dim]:
            return mask  # indices are still original slots
        lab = self.lab[dim]
        r = np.arange(self.n)[:, None, None]
        return mask[r, lab[:, :, None], lab[:, None, :]]

    def refresh(self, dim: int, mask: np.ndarray) -> None:
        """Sparse mode: make every permitted pair current, and track only those."""
        b, j, m = np.nonzero(mask & ~self.valid[dim])
        if b.size:
            Sv = self.view(dim)[0]
            Sj, Sm = Sv[b, j], Sv[b, m]  # (E, e, d, e) rows of the two groups
            self.C[dim][b, j, m] = np.einsum("ekal,ekal->e", Sj, Sm)
            Dj, Dm = self.D[dim][b, j], self.D[dim][b, m]
            self.T1[dim][b, j, m] = np.einsum("ekl,ekl->e", Dj, Dm)
            blk = Sj[np.arange(b.size), :, m]  # S[j, :, m, :]
            self.T2[dim][b, j, m] = np.einsum("ekl,elk->e", blk, blk)
            self.dT[dim][b, j, m] = 2.0 * np.einsum("ekk->e", blk)
        self.valid[dim] = mask.copy()


def _gather(arr: np.ndarray, keep: np.ndarray, axis: int) -> np.ndarray:
    """``arr`` restricted to indices ``keep[b]`` along ``axis``, per draw ``b``; contiguous.

    Batch-adjacent fancy indexing -- numpy's fast path, unlike ``take_along_axis``,
    which broadcasts an index array over every axis.
    """
    moved = np.moveaxis(arr, axis, 1)
    out = moved[np.arange(arr.shape[0])[:, None], keep]
    return np.ascontiguousarray(np.moveaxis(out, 1, axis))


def _deltas(st: _State, dim: int):
    """``(d_normsq, d_tr, d_tr2)`` for every label pair of ``dim``, each ``(n, d, d)``."""
    Mv = st.view(dim)[1]
    d_normsq = 2.0 * (Mv @ Mv.transpose(0, 2, 1))
    d_tr2 = 4.0 * st.C[dim] + 2.0 * (st.T1[dim] + st.T2[dim])
    return d_normsq, st.dT[dim], d_tr2


def _merge(st: _State, dim: int, idx: np.ndarray, u: np.ndarray, v: np.ndarray) -> None:
    """Merge label ``v`` into ``u`` along ``dim`` for draws ``idx`` (``u < v``, per draw).

    Everything is derived from ``Sigma' = P Sigma P^T`` with ``P^T P = I + E``,
    ``E = sum_k (e_uk e_vk^T + e_vk e_uk^T)``, so ``Sigma'^2 = P (Sigma^2 + Sigma E
    Sigma) P^T`` -- a rank-``2 d_other`` correction built from the two rows ``u, v``.
    Only those two rows of ``Sigma`` are read, and by symmetry every column block
    needed is a transpose of one of them.
    """
    o = 1 - dim
    Sv, Mv = st.view(dim)
    n = idx.shape[0]
    r = np.arange(n)

    # Rows u, v of the pre-merge Sigma: Su[., k, a, l] = S[u, k, a, l].
    rows = _Rows(Sv[idx, u], Sv[idx, v], u, v)
    Do = st.D[o][idx]  # (n, e, d, d): blocks S[:, h, :, h] of the other axis
    rows.Pu, rows.Pv = Do[r, :, u], Do[r, :, v]  # Pu[., j, x] = D_j[u, x]

    same = _same_axis_sparse if st.sparse[dim] else _same_axis_dense
    other = _other_axis_sparse if st.sparse[o] else _other_axis_dense
    CQ = same(st, dim, idx, rows)
    other(st, o, idx, rows)

    # ---- apply the merge ------------------------------------------------------
    # New row u (rows then columns, as merging.update_Sigma), written back as both
    # row u and, transposed, column u; row and column v are zeroed.
    new = rows.Su + rows.Sw
    new[r, :, u] += new[r, :, v]
    new[r, :, v] = 0.0
    Sv[idx, u] = new
    Sv[idx, v] = 0.0
    Sc = Sv.transpose(0, 3, 1, 2, 4)  # Sc[., m, j, k, l] = S[j, k, m, l]
    Sc[idx, u] = new.transpose(0, 2, 3, 1)
    Sc[idx, v] = 0.0
    Mv[idx, u] += Mv[idx, v]
    Mv[idx, v] = 0.0

    Do[r, :, u] += Do[r, :, v]
    Do[r, :, :, u] += Do[r, :, :, v]
    Do[r, :, v] = 0.0
    Do[r, :, :, v] = 0.0
    st.D[o][idx] = Do
    Dn = st.D[dim][idx]
    Dn[r, u] = new[r, :, u]
    Dn[r, v] = 0.0
    st.D[dim][idx] = Dn

    if CQ is not None:  # dense: C by the merge formula, T1/T2/dT row u recomputed
        CQ[r, u] += CQ[r, v]
        CQ[r, :, u] += CQ[r, :, v]
        CQ[r, v] = 0.0
        CQ[r, :, v] = 0.0
        st.C[dim][idx] = CQ
        for Tm, row in (
            (st.T1[dim], np.einsum("nkl,njkl->nj", Dn[r, u], Dn)),
            (st.T2[dim], np.einsum("nkml,nlmk->nm", new, new)),
            (st.dT[dim], 2.0 * np.einsum("nkmk->nm", new)),
        ):
            Tm[idx, u, :] = row
            Tm[idx, :, u] = row
            Tm[idx, v, :] = 0.0
            Tm[idx, :, v] = 0.0

    # Partition bookkeeping, in original-label slots.
    ou, ov = st.lab[dim][idx, u], st.lab[dim][idx, v]
    sizes, gid = st.sizes[dim], st.gid[dim]
    sizes[idx, ou] += sizes[idx, ov]
    sizes[idx, ov] = 0
    g = gid[idx]
    gid[idx] = np.where(g == ov[:, None], ou[:, None], g)


class _Rows:
    """The two pre-merge rows and the small blocks every update formula shares."""

    def __init__(self, Su, Sw, u, v):
        r = np.arange(Su.shape[0])
        self.Su, self.Sw = Su, Sw
        self.u, self.v = u, v
        self.Buu, self.Buv = Su[r, :, u], Su[r, :, v]  # (e x e) blocks S[u, :, u, :] ...
        self.Bvu, self.Bvv = Sw[r, :, u], Sw[r, :, v]


def _same_axis_dense(st, dim, idx, R):
    """C' = merge(C + W + W^T), W[a, b] = sum_{k,l} S[u,l,a,k] S[v,l,b,k]; returns C + W + W^T."""
    n, e, d, _ = R.Su.shape
    Ut = R.Su.transpose(0, 2, 1, 3).reshape(n, d, e * e)
    Vt = R.Sw.transpose(0, 2, 1, 3).reshape(n, d, e * e)
    W = Ut @ Vt.transpose(0, 2, 1)
    return st.C[dim][idx] + W + W.transpose(0, 2, 1)


def _same_axis_sparse(st, dim, idx, R):
    """W at the tracked pairs; pairs touching u or v drop out, to be refreshed if permitted."""
    r = np.arange(idx.shape[0])
    V = st.valid[dim][idx]
    V[r, R.u] = V[r, :, R.u] = V[r, R.v] = V[r, :, R.v] = False
    i, j, m = np.nonzero(V)
    if i.size:
        W2 = (np.einsum("elk,elk->e", R.Su[i, :, j], R.Sw[i, :, m])
              + np.einsum("elk,elk->e", R.Su[i, :, m], R.Sw[i, :, j]))
        st.C[dim][idx[i], j, m] += W2
    st.valid[dim][idx] = V
    return None


def _other_axis_dense(st, o, idx, R):
    """All other-axis pairs. C' = C + Z + Z^T + X + X^T, where Z sums the rank-2
    correction over the merged axis and X = Q[(u, .), (v, .)] is the new cross-block
    of Sigma'^2. Each (d x d) block is merged on both sides by P, so with (h, g) =
    (u, v): tr(B'B') = tr(BB) + 2[(BB)_gh + (BB)_hg] + B_hg^2 + 2 B_hh B_gg + B_gh^2,
    and the analogous expansion of <D_j', D_m'>."""
    Su, Sw = R.Su, R.Sw
    n, e, d, _ = Su.shape
    Z = Su.reshape(n, e * d, e).transpose(0, 2, 1) @ Sw.reshape(n, e * d, e)
    X = (Su.reshape(n, e, d * e) @ Sw.reshape(n, e, d * e).transpose(0, 2, 1)
         + R.Buu @ R.Bvv + R.Buv @ R.Buv)
    st.C[o][idx] += Z + Z.transpose(0, 2, 1) + X + X.transpose(0, 2, 1)
    F = np.einsum("nhxg,ngxh->nhg", Sw, Su)  # (BB)_vu per block
    st.T2[o][idx] += (2.0 * (F + F.transpose(0, 2, 1))
                      + R.Buv ** 2 + 2.0 * R.Buu * R.Bvv + R.Bvu ** 2)
    PP = R.Pv @ R.Pu.transpose(0, 2, 1)
    r = np.arange(n)
    a, b, c, f = R.Pu[r, :, R.v], R.Pu[r, :, R.u], R.Pv[r, :, R.v], R.Pv[r, :, R.u]
    st.T1[o][idx] += (2.0 * (PP + PP.transpose(0, 2, 1))
                      + a[:, :, None] * a[:, None, :] + b[:, :, None] * c[:, None, :]
                      + c[:, :, None] * b[:, None, :] + f[:, :, None] * f[:, None, :])
    st.dT[o][idx] += 2.0 * (R.Buv + R.Bvu)


def _other_axis_sparse(st, o, idx, R):
    """:func:`_other_axis_dense`, entry by entry at the tracked pairs only."""
    i, h, g = np.nonzero(st.valid[o][idx])
    if not i.size:
        return
    Su, Sw = R.Su, R.Sw
    E = i.size
    dot = lambda x, y: np.einsum("ek,ek->e", x.reshape(E, -1), y.reshape(E, -1))
    Z = dot(Su[i, :, :, h], Sw[i, :, :, g]) + dot(Su[i, :, :, g], Sw[i, :, :, h])
    X = (dot(Su[i, h], Sw[i, g]) + dot(Su[i, g], Sw[i, h])
         + dot(R.Buu[i, h], R.Bvv[i, g]) + dot(R.Buu[i, g], R.Bvv[i, h])
         + dot(R.Buv[i, h], R.Buv[i, :, g]) + dot(R.Buv[i, g], R.Buv[i, :, h]))
    rows = idx[i]
    st.C[o][rows, h, g] += Z + X
    F = dot(Sw[i, h, :, g], Su[i, g, :, h]) + dot(Sw[i, g, :, h], Su[i, h, :, g])
    st.T2[o][rows, h, g] += (2.0 * F + R.Buv[i, h, g] ** 2
                             + 2.0 * R.Buu[i, h, g] * R.Bvv[i, h, g] + R.Bvu[i, h, g] ** 2)
    Pu, Pv = R.Pu, R.Pv
    uu, vv = R.u[i], R.v[i]
    PP = dot(Pv[i, h], Pu[i, g]) + dot(Pv[i, g], Pu[i, h])
    st.T1[o][rows, h, g] += (2.0 * PP + Pu[i, h, vv] * Pu[i, g, vv]
                             + Pu[i, h, uu] * Pv[i, g, vv] + Pv[i, h, vv] * Pu[i, g, uu]
                             + Pv[i, h, uu] * Pv[i, g, uu])
    st.dT[o][rows, h, g] += 2.0 * (R.Buv[i, h, g] + R.Bvu[i, h, g])


def _run_chunk(shared: _Shared, T: np.ndarray, x_structure: Structure, y_structure: Structure,
               compact_threshold: float = COMPACT_THRESHOLD, sparse=None):
    """Search every column of ``T``; returns ``(values (L+1, n), merges)``.

    ``merges`` is filled only when ``n == 1`` (for :func:`greedy_search`): per
    level, ``(dim, i, j)`` with ``dim`` 0/1 and ``i < j`` the 1-based partition
    positions merged. ``sparse`` picks the carrying mode per dimension (see
    :class:`_State`); by default, sparse when the structure permits at most ``2 d``
    pairs of the finest partition.
    """
    structures = (x_structure, y_structure)
    dims0 = (shared.dx, shared.dy)
    if sparse is None:
        sparse = tuple(
            int(s.permitted_mask(np.ones((1, d), dtype=np.int64),
                                 np.arange(d)[None, :]).sum()) <= 2 * d
            for s, d in zip(structures, dims0)
        )
    st = _State(shared, T, sparse)
    n = st.n
    r = np.arange(n)

    values = [approx_chi_array(st.normsq, st.tr, st.tr2)]
    merges = []
    with np.errstate(divide="ignore", invalid="ignore"):
        while True:
            dims = st.lab[0].shape[1], st.lab[1].shape[1]  # current padded sizes
            offset = dims[0] * dims[0]
            V = np.full((n, offset + dims[1] * dims[1]), -np.inf)
            parts = []
            for dim in (0, 1):
                mask = st.permitted(dim, structures[dim])
                if st.sparse[dim]:
                    st.refresh(dim, mask)
                deltas = _deltas(st, dim)
                b, j, m = np.nonzero(mask)
                vals = approx_chi_array(st.normsq[b] + deltas[0][b, j, m],
                                        st.tr[b] + deltas[1][b, j, m],
                                        st.tr2[b] + deltas[2][b, j, m])
                vals[np.isnan(vals)] = -np.inf
                V[b, dim * offset + j * dims[dim] + m] = vals
                parts.append(deltas)
            has = np.isfinite(V).any(axis=1)
            if not has.any():
                break
            if not has.all():
                raise ValueError("Search paths have different lengths across draws.")

            best = np.argmax(V, axis=1)  # first max = R's which.max in loop order
            values.append(V[r, best])
            is_y = best >= offset
            for dim in (0, 1):
                sel = np.flatnonzero(is_y == bool(dim))
                if not sel.size:
                    continue
                flat = best[sel] - dim * offset
                u, v = flat // dims[dim], flat % dims[dim]
                d_normsq, d_tr, d_tr2 = parts[dim]
                st.normsq[sel] += d_normsq[sel, u, v]
                st.tr[sel] += d_tr[sel, u, v]
                st.tr2[sel] += d_tr2[sel, u, v]
                if n == 1:  # 1-based position = rank among active indices
                    rank = np.cumsum(st.active(dim)[0])
                    merges.append((dim, int(rank[u[0]]), int(rank[v[0]])))
                _merge(st, dim, sel, u, v)
            for dim in (0, 1):
                st.compact(dim, compact_threshold)
    return np.vstack(values), merges


# --------------------------------------------------------------------------- #
# Reference implementation
# --------------------------------------------------------------------------- #
def _greedy_search_loop(
    T_vector: np.ndarray,
    Sigma: np.ndarray,
    dx: int,
    dy: int,
    x_structure: Structure,
    y_structure: Structure,
    statistic=None,
) -> SearchResult:
    """The original per-candidate loop: one ``statistic.update`` per permitted merge."""
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
