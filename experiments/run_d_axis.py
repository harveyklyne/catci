"""The d axis: power (and size, at strength 0) as ``dx`` grows with ``n``, ``dy`` fixed.

TODO item 1: "What does the power curve do as ``d`` grows with ``n`` fixed --
does the merging advantage widen, as claimed?" Runs one :func:`config.d_grid`
block per ``dx`` and learner through :func:`run.run`, then prints rejection
rates as a ``method x dx`` table per learner and strength. Each block writes its
own parquet under ``results/`` (names carry ``_n{n}_dx{dx}_dy{dy}``, ``--tag`` and
the learner, e.g. ``power_lin_lin_step_n2000_dx32_dy4_pilot__oracle``).

Every ``(n, dx)`` and ``(n, dy)`` needs a tuning of each learner first
(``tune.py``); ``--learner oracle`` needs none.

Usage:
    python run_d_axis.py lin_lin_step --n 2000 --dy 4 --dxs 8,16,32,64 \\
        --strengths 0,0.8,1.6 --reps 30 --n-boot 200 --methods ordinal,max,chi_sq \\
        --tag pilot [--learner oracle] [--workers 5] [--report-only]
"""

from __future__ import annotations

import argparse
import numpy as np
import pandas as pd

import methods
import run as runner
from config import LEARNERS, d_grid, with_tag

ALPHA = 0.05


def report(df: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Rejection rate at ``alpha``, one ``method x dx`` table per learner and strength."""
    d = df[df.method != "__error__"].copy()
    if "learner" not in d:  # parquets from before results were learner-tagged
        d["learner"] = "untagged"
    d["reject"] = d.p_value < alpha
    tab = d.groupby(["learner", "strength", "method", "dx"])["reject"].agg(["mean", "count"]).reset_index()
    for (learner, s), block in tab.groupby(["learner", "strength"]):
        wide = block.pivot(index="method", columns="dx", values="mean")
        reps = int(block["count"].min())
        se = np.sqrt(max(alpha, 0.5 if s > 0 else alpha) * (1 - alpha) / reps)
        label = "size" if s == 0 else "power"
        print(f"\n### {learner}, strength {s} ({label} at alpha={alpha}; >= {reps} reps/cell; SE <= {se:.3f})")
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
    ap.add_argument("--learner", nargs="+", choices=[*LEARNERS, "oracle"], default=list(LEARNERS),
                    help="one block per learner and dx (default: mlp xgb)")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--report-only", action="store_true", help="re-read the parquets, do not run")
    args = ap.parse_args()

    parts = args.config.split("_")
    xs, ys, ints = parts[0], parts[1], "_".join(parts[2:])
    over = dict(reps=args.reps, n_boot=args.n_boot, strengths=[float(s) for s in args.strengths.split(",")])
    if args.methods:
        ms = args.methods.split(",")
        over["adaptive"] = [m for m in ms if m in methods.ADAPTIVE]
        over["competitors"] = [m for m in ms if m in methods.COMPETITORS]
    dxs = [int(v) for v in args.dxs.split(",")]
    cfgs = [with_tag(c, args.tag)
            for learner in args.learner
            for c in d_grid(xs, ys, ints, dxs, dy=args.dy, n=args.n, learner=learner, **over)]

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
