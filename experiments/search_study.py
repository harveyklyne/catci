"""Gaussian-limit power study for the TODO item 7 searches.

Compares, at fixed size, the greedy search against

* ``beam<w>``       -- beam search of width ``w`` (7c);
* ``random_fixed``  -- one data-independent random merge path per replicate, shared
                       by the observed and every bootstrap draw (7a);
* ``random_fresh``  -- a fresh random path for every draw, observed included (7a);
* ``split``         -- greedy path chosen on half A, the whole path evaluated and
                       minP-calibrated on half B (7b);
* ``split_one``     -- as ``split``, but only the single depth that looked best on
                       half A is tested on half B (``L = 1``);

plus the depth-0 statistic and the oracle chi-square test as non-adaptive anchors.

**Why the Gaussian limit.** ``T ~ N(mu, Sigma)`` with ``Sigma`` known is exactly what
the search and the calibration see asymptotically, it isolates the *search* from
the propensity learner, and it makes sample splitting exact: two halves give
``T_A, T_B ~ N(mu / sqrt 2, Sigma)`` independently, and the full-sample statistic is
``(T_A + T_B) / sqrt 2``. Every method in a replicate sees the same ``T_A, T_B`` and
the same bootstrap draws, so method differences are paired.

**Why Kronecker ``Sigma``.** With no ``Z`` the null covariance of the residual
products is ``C_Y (x) C_X`` with ``C = diag(p) - p p^T``. Merging preserves the
product form, so ``tr`` and ``tr(Sigma^2)`` factorise over the two dimensions and every
candidate at a level is scored by a handful of ``g x g`` array operations. That is
what makes width-25 beams at ``d = 8`` and ``n_boot = 1000`` affordable here.
:class:`KronSearch` is checked path-for-path against :mod:`catci.search` by
``tests/test_search_study.py``; it is experiment code, not a package fast path (for
general ``Sigma`` that is TODO item 8).

Usage:
    python search_study.py greedy_trap --structure saturated --deltas 0,4,5 --reps 500
    python search_study.py binary_tree --structure tree --deltas 0,3,4 --reps 500

Writes ``results/search_<dgp>_<structure>.parquet`` (one row per replicate x
method) and prints rejection rates at ``alpha = 0.05``.
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import chi2

import dgp
from catci.calibrate import double_bootstrap_pvalue
from catci.structure import Ordinal, Saturated, Structure, Tree

RESULTS_DIR = Path(__file__).resolve().parent / "results"
ALPHA = 0.05

Merge = Tuple[int, int, int]


# --------------------------------------------------------------------------- #
# Search under Kronecker Sigma
# --------------------------------------------------------------------------- #
def _box(normsq, tr, tr2):
    """``ApproxChi`` value, vectorised: ``chi2.cdf(normsq / g, h)`` = ``gammainc(h/2, normsq/2g)``."""
    g = tr2 / tr
    h = tr ** 2 / tr2
    return gammainc(h / 2.0, normsq / (2.0 * g))


@dataclass
class _State:
    Tm: np.ndarray  # (gx, gy) merged T, X along rows
    A: dict  # {1: merged C_X, 2: merged C_Y}
    partition: dict
    normsq: float
    tr: dict = field(default_factory=dict)
    tr2: dict = field(default_factory=dict)

    @classmethod
    def root(cls, T_vector, Cx, Cy):
        dx, dy = Cx.shape[0], Cy.shape[0]
        Tm = np.asarray(T_vector, dtype=float).reshape(dy, dx).T  # X fastest
        A = {1: Cx, 2: Cy}
        return cls(
            Tm=Tm, A=A,
            partition={"x": [[j] for j in range(1, dx + 1)], "y": [[k] for k in range(1, dy + 1)]},
            normsq=float(np.sum(Tm ** 2)),
            tr={d: float(np.trace(A[d])) for d in (1, 2)},
            tr2={d: float(np.sum(A[d] ** 2)) for d in (1, 2)},
        )

    def value(self) -> float:
        return float(_box(self.normsq, self.tr[1] * self.tr[2], self.tr2[1] * self.tr2[2]))

    def score(self, dimension: int, pairs) -> np.ndarray:
        """Values after each merge in ``pairs`` (1-based positions) along ``dimension``."""
        return _box(*self.moments(dimension, pairs))

    def moments(self, dimension: int, pairs):
        """``(normsq, tr, tr2)`` after each merge in ``pairs`` along ``dimension``."""
        I = np.array([i for i, _ in pairs]) - 1
        J = np.array([j for _, j in pairs]) - 1
        M = self.Tm if dimension == 1 else self.Tm.T
        A = self.A[dimension]
        other = 2 if dimension == 1 else 1
        normsq = self.normsq + 2.0 * np.einsum("ik,ik->i", M[I], M[J])
        tr = (self.tr[dimension] + 2.0 * A[I, J]) * self.tr[other]
        AA = A @ A
        tr2 = (self.tr2[dimension] + 4.0 * AA[I, J]
               + 2.0 * (A[I, I] * A[J, J] + A[I, J] ** 2)) * self.tr2[other]
        return normsq, tr, tr2

    def apply(self, merge: Merge) -> "_State":
        dimension, i, j = merge
        i, j = i - 1, j - 1
        M = self.Tm if dimension == 1 else self.Tm.T
        A = self.A[dimension]
        M2 = M.copy()
        M2[i] += M[j]
        M2 = np.delete(M2, j, axis=0)
        A2 = A.copy()
        A2[i, :] += A2[j, :]
        A2[:, i] += A2[:, j]
        A2 = np.delete(np.delete(A2, j, axis=0), j, axis=1)
        key = "x" if dimension == 1 else "y"
        part = {k: [list(g) for g in v] for k, v in self.partition.items()}
        part[key][i] = part[key][i] + part[key][j]
        del part[key][j]
        newA = dict(self.A)
        newA[dimension] = A2
        tr, tr2 = dict(self.tr), dict(self.tr2)
        tr[dimension] = float(np.trace(A2))
        tr2[dimension] = float(np.sum(A2 ** 2))
        return _State(
            Tm=M2 if dimension == 1 else M2.T, A=newA, partition=part,
            normsq=float(np.sum(M2 ** 2)), tr=tr, tr2=tr2,
        )

    def done(self) -> bool:
        return len(self.partition["x"]) <= 2 and len(self.partition["y"]) <= 2

    def key(self) -> tuple:
        return tuple(tuple(sorted(tuple(sorted(g)) for g in self.partition[k])) for k in ("x", "y"))

    def child_key(self, merge: Merge) -> tuple:
        """``self.apply(merge).key()`` without building the child."""
        dimension, i, j = merge
        out = []
        for d, k in ((1, "x"), (2, "y")):
            groups = [tuple(sorted(g)) for g in self.partition[k]]
            if d == dimension:
                merged = tuple(sorted(groups[i - 1] + groups[j - 1]))
                groups = [g for q, g in enumerate(groups) if q not in (i - 1, j - 1)] + [merged]
            out.append(tuple(sorted(groups)))
        return tuple(out)


class KronSearch:
    """Greedy / beam / fixed-path search for ``Sigma = C_Y (x) C_X`` (see module docstring)."""

    def __init__(self, Cx, Cy, xs: Structure, ys: Structure):
        self.Cx, self.Cy, self.xs, self.ys = Cx, Cy, xs, ys

    def _moments(self, state: _State):
        """(normsq, tr, tr2) arrays and merges for every permitted merge, in greedy's loop order."""
        mom, merges = [], []
        for dimension, st, key in ((1, self.xs, "x"), (2, self.ys, "y")):
            pairs = st.permitted_merges(state.partition[key])
            if pairs:
                mom.append(state.moments(dimension, pairs))
                merges.extend((dimension, i, j) for i, j in pairs)
        return mom, merges

    def beam(self, T_vector, width: int = 1):
        """Beam search; ``width = 1`` is greedy. Returns (values, merges of the best state)."""
        beam = [_State.root(T_vector, self.Cx, self.Cy)]
        values, merges = [beam[0].value()], []
        while not all(s.done() for s in beam):
            mom, cand = [], []
            for b, s in enumerate(beam):
                m, mm = self._moments(s)
                mom.extend(m)
                cand.extend((b, x) for x in mm)
            if not cand:
                break
            vals = _box(*(np.concatenate([m[q] for m in mom]) for q in range(3)))
            order = np.argsort(-vals, kind="stable")  # ties keep loop order
            new, seen = [], set()
            for k in order:
                b, merge = cand[k]
                if width > 1:
                    key = beam[b].child_key(merge)
                    if key in seen:
                        continue
                    seen.add(key)
                if not new:
                    values.append(float(vals[k]))
                    merges.append(merge)
                new.append(beam[b].apply(merge))
                if len(new) == width:
                    break
            beam = new
        return np.asarray(values), merges

    def evaluate(self, T_vector, merges) -> np.ndarray:
        s = _State.root(T_vector, self.Cx, self.Cy)
        values = [s.value()]
        for merge in merges:
            s = s.apply(merge)
            values.append(s.value())
        return np.asarray(values)


# --------------------------------------------------------------------------- #
# DGPs: a direction for mu, normalised so that mu' Sigma^+ mu = delta^2
# --------------------------------------------------------------------------- #
def trap_interaction(d: int, a: float = 1.0, b: float = 3.6) -> np.ndarray:
    """Built to trap greedy (TODO 7c): a coarse 2x2 checkerboard (amplitude ``a``)
    plus a stronger fine signal (``b``) shared by row pairs ``(k, k + d/2)`` that sit
    on *opposite* sides of the coarse split, each pair on its own column pair.

    Greedy merges the fine pairs first -- they have the largest inner products --
    and in doing so cancels the coarse signal along X for good. In a narrow window
    of ``b`` the coarse 2x2 is nonetheless the better collapse; ``b = 3.6`` at
    ``d = 8`` is the worst case for greedy found by scanning ``b`` noise-free.
    """
    h = d // 2
    s = np.where(np.arange(d) < h, 1.0, -1.0)
    D = a * np.outer(s, s)
    for k in range(h):
        for r in (k, k + h):
            D[r, 2 * k] += b
            D[r, 2 * k + 1] -= b
    return D


def direction(name: str, d: int) -> np.ndarray:
    D = trap_interaction(d) if name == "greedy_trap" else dgp.get_int(name, d, d)
    D = D - D.mean(0, keepdims=True) - D.mean(1, keepdims=True) + D.mean()  # range of Sigma
    return D.reshape(-1, order="F")  # X fastest


def uniform_C(d: int) -> np.ndarray:
    p = np.full(d, 1.0 / d)
    return np.diag(p) - np.outer(p, p)


STRUCTURES = {
    "saturated": lambda d: Saturated(),
    "tree": lambda d: Tree.binary(d),
    "ordinal": lambda d: Ordinal(),
}


# --------------------------------------------------------------------------- #
# One replicate
# --------------------------------------------------------------------------- #
def replicate(dgp_name, structure, d, delta, n_boot, widths, seed_seq):
    rng = np.random.default_rng(seed_seq)
    C = uniform_C(d)
    Sigma = np.kron(C, C)
    Sp = np.linalg.pinv(Sigma)
    u = direction(dgp_name, d)
    mu = delta * u / np.sqrt(u @ Sp @ u)

    vals, vecs = np.linalg.eigh(Sigma)
    keep = vals > 1e-12
    root = vecs[:, keep] * np.sqrt(vals[keep])
    draw = lambda k: root @ rng.standard_normal((root.shape[1], k))  # noqa: E731

    T_A = mu / np.sqrt(2) + draw(1)[:, 0]
    T_B = mu / np.sqrt(2) + draw(1)[:, 0]
    T = (T_A + T_B) / np.sqrt(2)
    Z = draw(n_boot)

    st = STRUCTURES[structure](d)
    ks = KronSearch(C, C, st, st)
    out = {}

    def minp(obs, boot):
        return double_bootstrap_pvalue(np.asarray(obs), np.column_stack(boot), rng)

    # greedy and beams: the same search on the observed and every draw
    for w in widths:
        name = "greedy" if w == 1 else f"beam{w}"
        obs = ks.beam(T, w)[0]
        out[name] = minp(obs, [ks.beam(Z[:, b], w)[0] for b in range(n_boot)])
        if w == 1:
            out["depth0"] = double_bootstrap_pvalue(
                obs[:1], np.array([[ks.evaluate(Z[:, b], [])[0] for b in range(n_boot)]]), rng
            )

    # 7a random paths
    L = 2 * d - 4
    path = random_path(ks, d, rng)
    assert len(path) == L
    out["random_fixed"] = minp(ks.evaluate(T, path), [ks.evaluate(Z[:, b], path) for b in range(n_boot)])
    out["random_fresh"] = minp(
        ks.evaluate(T, random_path(ks, d, rng)),
        [ks.evaluate(Z[:, b], random_path(ks, d, rng)) for b in range(n_boot)],
    )

    # 7b sample splitting: choose on A, test on B
    vA, pathA = ks.beam(T_A, 1)
    boot_B = [ks.evaluate(Z[:, b], pathA) for b in range(n_boot)]
    out["split"] = minp(ks.evaluate(T_B, pathA), boot_B)
    l_star = int(np.argmax(vA))
    out["split_one"] = double_bootstrap_pvalue(
        ks.evaluate(T_B, pathA)[l_star:l_star + 1],
        np.array([[bb[l_star] for bb in boot_B]]), rng,
    )

    # oracle chi-square on the full sample (non-adaptive anchor)
    out["chi_sq"] = float(chi2.sf(T @ Sp @ T, df=(d - 1) ** 2))
    return out


def random_path(ks: KronSearch, d: int, rng) -> List[Merge]:
    from catci.search import random_merges

    return random_merges(d, d, ks.xs, ks.ys, rng)


def _task(args):
    dgp_name, structure, d, delta, rep, n_boot, widths, seed_seq = args
    t0 = time.time()
    pv = replicate(dgp_name, structure, d, delta, n_boot, widths, seed_seq)
    sec = time.time() - t0
    return [dict(dgp=dgp_name, structure=structure, d=d, delta=delta, rep=rep, method=m,
                 p_value=p, seconds=sec) for m, p in pv.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dgp", choices=["greedy_trap", "binary_tree", "step", "alt"])
    ap.add_argument("--structure", default="saturated", choices=list(STRUCTURES))
    ap.add_argument("--d", type=int, default=8)
    ap.add_argument("--deltas", default="0,4,5")
    ap.add_argument("--widths", default="1,5,25")
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    deltas = [float(x) for x in args.deltas.split(",")]
    widths = [int(x) for x in args.widths.split(",")]
    children = np.random.SeedSequence(args.seed).spawn(len(deltas) * args.reps)
    tasks = [(args.dgp, args.structure, args.d, delta, rep, args.n_boot, widths, children[k * args.reps + rep])
             for k, delta in enumerate(deltas) for rep in range(args.reps)]

    t0 = time.time()
    rows = []
    if args.workers == 1:
        for t in tasks:
            rows.extend(_task(t))
    else:
        import multiprocessing as mp
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("fork")) as ex:
            for res in ex.map(_task, tasks, chunksize=2):
                rows.extend(res)
    elapsed = time.time() - t0

    df = pd.DataFrame(rows)
    RESULTS_DIR.mkdir(exist_ok=True)
    stem = f"search_{args.dgp}_{args.structure}_d{args.d}{args.tag}"
    df.to_parquet(RESULTS_DIR / f"{stem}.parquet", index=False)
    (RESULTS_DIR / f"{stem}.provenance.json").write_text(json.dumps(
        dict(vars(args), runtime_seconds=round(elapsed, 1)), indent=2))
    print(f"DONE {stem}: {len(tasks)} replicates in {elapsed / 60:.1f} min")
    report(df, args.reps)


def report(df, reps):
    se = np.sqrt(0.25 / reps)
    tab = (df.assign(reject=df.p_value < ALPHA)
           .pivot_table(index="method", columns="delta", values="reject", aggfunc="mean"))
    print(f"rejection rate at alpha={ALPHA}  (reps={reps}, SE <= {se:.3f})")
    print(tab.round(3).to_string())


if __name__ == "__main__":
    main()
