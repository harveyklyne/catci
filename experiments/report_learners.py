"""Read the sweep's parquets and print the MLP-vs-XGBoost comparison tables.

Power is compared *paired*: ``run.py`` seeds from ``SeedSequence(seed)`` and the
learner does not enter the seed stream, so the two runs of a config saw the same
simulated datasets replicate for replicate. The paired difference in rejection
therefore has far less Monte Carlo noise than the two rates do separately, and
this prints its standard error so a gap can be read as real or not.

Usage:
    python report_learners.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS_DIR = Path(__file__).resolve().parent / "results"
ALPHA = 0.05
ORDER = ["tree", "ordinal", "max", "euclid", "mGCM", "ankan", "chi_sq", "multinomial"]


def _load(stem: str) -> pd.DataFrame | None:
    path = RESULTS_DIR / f"{stem}.parquet"
    if not path.exists():
        return None
    df = pd.read_parquet(path)
    return df[df.method != "__error__"]


def size_table(settings, learners):
    print("\n=== SIZE: rejection rate at alpha = 0.05 (want 0.050) ===")
    for setting in settings:
        frames = {}
        for learner in learners:
            stem = f"size_{setting}" + ("" if learner == "xgb" else f"__{learner}")
            df = _load(stem)
            if df is not None:
                frames[learner] = df
        if not frames:
            continue
        reps = max(len(f) // f.method.nunique() for f in frames.values())
        se = np.sqrt(ALPHA * (1 - ALPHA) / reps)
        print(f"\n### {setting}   reps={reps}, binomial SE~{se:.4f} "
              f"(|rate-0.05| > {2*se:.3f} is notable)")
        print(f"    {'method':10s}" + "".join(f"{k:>10s}" for k in frames))
        for m in ORDER:
            cells = []
            for df in frames.values():
                sub = df[df.method == m]
                cells.append(f"{(sub.p_value < ALPHA).mean():10.3f}" if len(sub) else f"{'-':>10s}")
            if any(c.strip() != "-" for c in cells):
                print(f"    {m:10s}" + "".join(cells))


def power_table(configs, learners=("xgb", "mlp")):
    print("\n=== POWER: rejection rate at alpha = 0.05, paired on identical data ===")
    for cfg in configs:
        a = _load(cfg)
        b = _load(f"{cfg}__{learners[1]}")
        if a is None or b is None:
            continue
        key = ["method", "strength", "rep"]
        merged = a.merge(b, on=key, suffixes=("_xgb", "_mlp"))
        merged["rej_xgb"] = merged.p_value_xgb < ALPHA
        merged["rej_mlp"] = merged.p_value_mlp < ALPHA
        merged["diff"] = merged.rej_mlp.astype(float) - merged.rej_xgb.astype(float)

        print(f"\n### {cfg}   ({merged.rep.nunique()} reps/strength)")
        print(f"    {'method':10s}{'s':>5s}{'xgb':>8s}{'mlp':>8s}{'diff':>9s}{'se(diff)':>10s}")
        for m in ORDER:
            sub = merged[merged.method == m]
            if not len(sub):
                continue
            for s, grp in sub.groupby("strength"):
                d = grp["diff"]
                se = d.std(ddof=1) / np.sqrt(len(d))
                mark = " *" if abs(d.mean()) > 2 * se else ""
                print(f"    {m:10s}{s:5.1f}{grp.rej_xgb.mean():8.3f}"
                      f"{grp.rej_mlp.mean():8.3f}{d.mean():+9.3f}{se:10.3f}{mark}")

        # one headline number per method: mean power across the strength grid
        print(f"\n    mean power over the grid:")
        for m in ORDER:
            sub = merged[merged.method == m]
            if not len(sub):
                continue
            per_s = sub.groupby("strength")[["rej_xgb", "rej_mlp"]].mean()
            print(f"      {m:10s} xgb={per_s.rej_xgb.mean():.3f}  mlp={per_s.rej_mlp.mean():.3f}"
                  f"  ({per_s.rej_mlp.mean() - per_s.rej_xgb.mean():+.3f})")


def propensity_table():
    path = RESULTS_DIR / "bench_learners.parquet"
    if not path.exists():
        return
    df = pd.read_parquet(path)
    agg = df.groupby(["setting", "learner"]).mean(numeric_only=True)
    print("\n=== PROPENSITY QUALITY: E_f = max_j E[(f_j - fhat_j)^2] ===")
    print(f"    {'setting':8s}{'E_f xgb':>10s}{'E_f mlp':>10s}{'ratio':>8s}"
          f"{'n*Ef*Eg xgb':>13s}{'n*Ef*Eg mlp':>13s}")
    for setting in agg.index.get_level_values(0).unique():
        sub = agg.loc[setting]
        x, m = sub.loc["xgb"], sub.loc["mlp"]
        print(f"    {setting:8s}{x.E_f_out:10.5f}{m.E_f_out:10.5f}"
              f"{x.E_f_out/m.E_f_out:8.2f}{1000*x.E_f_in*x.E_g_in:13.4f}"
              f"{1000*m.E_f_in*m.E_g_in:13.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", nargs="*", default=["lin_lin", "sin_sin"])
    ap.add_argument("--power", nargs="*",
                    default=["power_lin_lin_step", "power_sin_sin_binary_tree"])
    args = ap.parse_args()

    pd.set_option("display.width", 200)
    propensity_table()
    size_table(args.size, ["oracle", "xgb", "mlp"])
    power_table(args.power)


if __name__ == "__main__":
    main()
