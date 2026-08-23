"""Size (null-calibration) runs.

Runs one or more cheap size settings at strength 0 and reports the rejection
rate at alpha per method.

Usage:
    python run_size.py lin_lin lin_vee --reps 1000 --workers 5
"""

from __future__ import annotations

import argparse

import numpy as np

from config import size_config
import run as runner

ALPHA = 0.05


def report(df, cfg):
    d = df[df.method != "__error__"].copy()
    d["reject"] = d.p_value < ALPHA
    tab = d.groupby("method")["reject"].agg(["mean", "count"])
    se = np.sqrt(ALPHA * (1 - ALPHA) / cfg.reps)
    print(f"\n### size {cfg.xsetting}_{cfg.ysetting}  |  reps={cfg.reps}")
    print(f"    nominal alpha={ALPHA:.3f}  binomial SE~{se:.3f}  (|rate-0.05|>~{2*se:.3f} is notable)")
    order = ["tree", "ordinal", "max", "euclid", "mGCM", "ankan", "chi_sq"]
    for m in order:
        if m in tab.index:
            rate = tab.loc[m, "mean"]
            flag = "" if abs(rate - ALPHA) <= 2 * se else ("  <-- inflated" if rate > ALPHA else "  <-- conservative")
            print(f"    {m:9s} reject@0.05 = {rate:.3f}{flag}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("settings", nargs="+", help="e.g. lin_lin lin_vee")
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    for s in args.settings:
        xs, ys = s.split("_")
        cfg = size_config(xs, ys, reps=args.reps)
        df = runner.run(cfg, workers=args.workers, seed=args.seed)
        report(df, cfg)


if __name__ == "__main__":
    main()
