"""The d axis: power (and size, at strength 0) as ``dx`` grows with ``n``, ``dy`` fixed.

TODO item 1: "What does the power curve do as ``d`` grows with ``n`` fixed --
does the merging advantage widen, as claimed?" Runs one :func:`config.d_grid`
block per ``dx`` through :func:`run.run`, then prints rejection rates as a
``method x dx`` table per strength. Each block writes its own parquet under
``results/`` (names carry ``_n{n}_dx{dx}_dy{dy}`` and ``--tag``).

Every ``(n, dx)`` and ``(n, dy)`` needs a tuning first (``tune.py``).

Usage:
    python run_d_axis.py lin_lin_step --n 2000 --dy 4 --dxs 8,16,32,64 \\
        --strengths 0,0.8,1.6 --reps 30 --n-boot 200 --methods ordinal,max,chi_sq \\
        --tag pilot [--workers 5] [--report-only]
"""

from __future__ import annotations

import argparse
from dataclasses import replace

import numpy as np
import pandas as pd

import methods
import run as runner
from config import d_grid

ALPHA = 0.05


def report(df: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Rejection rate at ``alpha``, one ``method x dx`` table per strength (n_reps in brackets)."""
    d = df[df.method != "__error__"].copy()
    d["reject"] = d.p_value < alpha
    tab = d.groupby(["strength", "method", "dx"])["reject"].agg(["mean", "count"]).reset_index()
    for s, block in tab.groupby("strength"):
        wide = block.pivot(index="method", columns="dx", values="mean")
        reps = int(block["count"].min())
        se = np.sqrt(max(alpha, 0.5 if s > 0 else alpha) * (1 - alpha) / reps)
        label = "size" if s == 0 else "power"
        print(f"\n### strength {s} ({label} at alpha={alpha}; >= {reps} reps/cell; SE <= {se:.3f})")
        print(wide.round(3).to_string())
    n_err = int((df.method == "__error__").sum())
    if n_err:
        print(f"\n!! {n_err} replicates errored:", df.loc[df.method == "__error__", "error"].unique()[:3])
    return tab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", help="<x>_<y>_<int>, e.g. lin_lin_step")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--dy", type=int, required=True)
    ap.add_argument("--dxs", type=str, required=True)
    ap.add_argument("--strengths", type=str, default="0,0.8,1.6")
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--methods", type=str, default=None)
    ap.add_argument("--tag", type=str, default=None)
    ap.add_argument("--learner", choices=["xgboost", "oracle"], default="xgboost")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--report-only", action="store_true", help="re-read the parquets, do not run")
    args = ap.parse_args()

    parts = args.config.split("_")
    xs, ys, ints = parts[0], parts[1], "_".join(parts[2:])
    over = dict(learner=args.learner, reps=args.reps, n_boot=args.n_boot, strengths=[float(s) for s in args.strengths.split(",")])
    if args.methods:
        ms = args.methods.split(",")
        over["adaptive"] = [m for m in ms if m in methods.ADAPTIVE]
        over["competitors"] = [m for m in ms if m in methods.COMPETITORS]
    cfgs = d_grid(xs, ys, ints, [int(v) for v in args.dxs.split(",")], dy=args.dy, n=args.n, **over)
    if args.tag:
        cfgs = [replace(c, name=f"{c.name}_{args.tag}") for c in cfgs]

    frames = []
    for cfg in cfgs:
        path = runner.RESULTS_DIR / f"{cfg.name}.parquet"
        if args.report_only:
            if path.exists():
                frames.append(pd.read_parquet(path))
            continue
        frames.append(runner.run(cfg, workers=args.workers, seed=args.seed))
    if frames:
        report(pd.concat(frames, ignore_index=True))


if __name__ == "__main__":
    main()
