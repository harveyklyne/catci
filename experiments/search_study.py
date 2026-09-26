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

``Sigma = C_Y (x) C_X`` with ``C = diag(p) - p p^T`` and uniform ``p``: the null
covariance of the residual products when there is no ``Z``. Every search is the
package's own batched code (:mod:`catci.search`), so this measures what ships.

Usage:
    python search_study.py greedy_trap --structure saturated --deltas 0,4,5 --reps 500
    python search_study.py binary_tree --structure tree --deltas 0,3,4 --reps 500

Writes ``results/search_<dgp>_<structure>_d<d><tag>.parquet`` (one row per
replicate x method) and prints rejection rates at ``alpha = 0.05``; see
``summarise_search.py`` for paired differences.
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

import dgp
from catci.bootstrap import matrix_sqrt
from catci.calibrate import double_bootstrap_pvalue
from catci.search import (beam_search_paths, evaluate_paths, greedy_search,
                          greedy_search_paths, random_merges)
from catci.structure import Ordinal, Saturated, Tree

RESULTS_DIR = Path(__file__).resolve().parent / "results"
ALPHA = 0.05


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

    root = matrix_sqrt(Sigma)
    draw = lambda k: root @ rng.standard_normal((root.shape[1], k))  # noqa: E731

    T_A = mu / np.sqrt(2) + draw(1)[:, 0]
    T_B = mu / np.sqrt(2) + draw(1)[:, 0]
    T = (T_A + T_B) / np.sqrt(2)
    Z = draw(n_boot)
    T_all = np.column_stack([T, Z])  # observed is column 0

    st = STRUCTURES[structure](d)
    args = (Sigma, d, d, st, st)
    out = {}

    def minp(paths):
        return double_bootstrap_pvalue(paths[:, 0], paths[:, 1:], rng)

    # greedy and beams: the same search on the observed and every draw
    for w in widths:
        if w == 1:
            paths = greedy_search_paths(T_all, *args)
            out["greedy"] = minp(paths)
            out["depth0"] = minp(paths[:1])
        else:
            out[f"beam{w}"] = minp(beam_search_paths(T_all, *args, width=w))

    # 7a random paths
    out["random_fixed"] = minp(evaluate_paths(T_all, Sigma, d, d, random_merges(d, d, st, st, rng)))
    out["random_fresh"] = minp(np.column_stack([
        evaluate_paths(T_all[:, b], Sigma, d, d, random_merges(d, d, st, st, rng))[:, 0]
        for b in range(n_boot + 1)
    ]))

    # 7b sample splitting: choose on A, test on B
    chosen = greedy_search(T_A, *args)
    paths_B = evaluate_paths(np.column_stack([T_B, Z]), Sigma, d, d, chosen.merges)
    out["split"] = minp(paths_B)
    l_star = int(np.argmax(chosen.values))
    out["split_one"] = minp(paths_B[l_star:l_star + 1])

    # oracle chi-square on the full sample (non-adaptive anchor)
    out["chi_sq"] = float(chi2.sf(T @ Sp @ T, df=(d - 1) ** 2))
    return out


def _task(args):
    dgp_name, structure, d, delta, rep, n_boot, widths, seed_seq = args
    t0 = time.process_time()
    pv = replicate(dgp_name, structure, d, delta, n_boot, widths, seed_seq)
    sec = time.process_time() - t0
    return [dict(dgp=dgp_name, structure=structure, d=d, delta=delta, rep=rep, method=m,
                 p_value=p, cpu_seconds=sec) for m, p in pv.items()]


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
    print(f"DONE {stem}: {len(tasks)} replicates in {elapsed / 60:.1f} min", flush=True)
    report(df, args.reps)


def report(df, reps):
    se = np.sqrt(0.25 / reps)
    tab = (df.assign(reject=df.p_value < ALPHA)
           .pivot_table(index="method", columns="delta", values="reject", aggfunc="mean"))
    print(f"rejection rate at alpha={ALPHA}  (reps={reps}, SE <= {se:.3f})")
    print(tab.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
