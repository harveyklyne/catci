"""Power and size of the divisive search against the merge search, paired.

Every replicate draws one dataset and one set of bootstrap draws, and every
search direction is calibrated on exactly those -- so differences between
methods are paired, not two independent Monte Carlo errors. Propensities are the
oracle ``f, g`` from the DGP, which isolates the search from the regression fit
(and sidesteps ``tuning/`` having hyperparameters for ``n = 1000, d = 8`` only).

Methods, for the chosen structure:

* ``merge``      -- :class:`~catci.search.MergeSearch`, the paper's Algorithm 1
* ``split``      -- :class:`~catci.search.SplitSearch` run to the singletons
* ``split@k``    -- the same, stopped after ``k`` levels

Usage:
    python power_divisive.py --d 8 --int binary_tree --strengths 0 0.4 0.8 1.2 1.6
    python power_divisive.py --d 8 --int alt --structure tree
    python power_divisive.py --d 32 --n 8000 --no-full   # truncated splits only
"""

from __future__ import annotations

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import dgp
from catci.bootstrap import bootstrap_T
from catci.calibrate import double_bootstrap_pvalue
from catci.gcm import form_t_sigma
from catci.search import MergeSearch, SplitSearch
from catci.statistic import ApproxChi
from catci.structure import Ordinal, Tree

OUT = Path(__file__).resolve().parent / "results_divisive"


def searches(truncations, full: bool) -> dict:
    out = {"merge": MergeSearch(), "split": SplitSearch()} if full else {}
    out.update({f"split@{k}": SplitSearch(max_levels=k) for k in truncations})
    return out


def structures(name: str, d: int):
    if name == "tree":
        return Tree.binary(d), Tree.binary(d)
    return Ordinal(), Ordinal()


def one_rep(task) -> list[dict]:
    args, strength, rep, seed = task
    rng = np.random.default_rng(seed)
    data = dgp.simulate_data(
        args.n, args.d, args.d, args.margin, args.margin, strength, args.int, rng=rng,
    )
    ts = form_t_sigma(data["x"], data["y"], data["f"], data["g"], normalise=False)
    xs, ys = structures(args.structure, args.d)
    statistic = ApproxChi()
    boot_T = bootstrap_T(ts.Sigma, args.n_boot, rng)

    rows = []
    for name, search in searches(args.truncations, args.full).items():
        t0 = time.perf_counter()
        path = search.prepare(ts.Sigma, args.d, args.d, xs, ys, statistic)
        observed = path(ts.T_vector)
        boot = np.column_stack([path(boot_T[:, b]).values for b in range(args.n_boot)])
        p = double_bootstrap_pvalue(np.asarray(observed.values), boot, rng)
        rows.append(dict(
            method=name, strength=strength, rep=rep, p_value=p,
            levels=len(observed.values),
            # where the observed path peaks: 0 = the search's first level
            argmax_level=int(np.argmax(observed.values)),
            seconds=time.perf_counter() - t0,
        ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, default=8)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--int", default="binary_tree", choices=["binary_tree", "alt", "step"])
    ap.add_argument("--structure", default="tree", choices=["tree", "ordinal"])
    ap.add_argument("--margin", default="lin")
    ap.add_argument("--strengths", type=float, nargs="+", default=[0.0, 0.4, 0.8, 1.2, 1.6])
    ap.add_argument("--truncations", type=int, nargs="+", default=[2, 4])
    ap.add_argument("--no-full", dest="full", action="store_false",
                    help="skip merge and full split (for large d)")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    seeds = np.random.SeedSequence(args.seed).spawn(len(args.strengths) * args.reps)
    tasks = [
        (args, s, r, seeds[i * args.reps + r])
        for i, s in enumerate(args.strengths) for r in range(args.reps)
    ]
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        rows = [row for rep_rows in ex.map(one_rep, tasks, chunksize=1) for row in rep_rows]
    df = pd.DataFrame(rows)

    OUT.mkdir(exist_ok=True)
    stem = f"{args.structure}_{args.int}_d{args.d}_n{args.n}_B{args.n_boot}{args.tag}"
    df.to_csv(OUT / f"{stem}.csv", index=False)

    se = np.sqrt(0.05 * 0.95 / args.reps)
    print(f"\n### {stem}  reps={args.reps}  ({time.time() - t0:.0f}s)  size SE~{se:.3f}")
    df["reject"] = df.p_value < 0.05
    print(df.pivot_table(index="strength", columns="method", values="reject", sort=False)
            .round(3).to_string())
    print("\nmean seconds per test:")
    print(df.groupby("method", sort=False)["seconds"].mean().round(2).to_string())


if __name__ == "__main__":
    main()
