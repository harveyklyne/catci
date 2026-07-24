"""Compare the Python pipeline against the R FIXED run on lin_lin_binary_tree.

Runs the Python grid (tree + ordinal + comparators + competitors, same xgb
hyperparameters as R) and prints power tables side-by-side with
results/power_lin_lin_binary_tree_FIXED.csv, plus per-replicate timing vs the
R timing log (1800 reps in 11.0 min on 5 workers).
"""

from __future__ import annotations

import sys
import time
from dataclasses import replace

import numpy as np
import pandas as pd

from config import power_config, REPO_ROOT
import run as runner

ALPHA = 0.05
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 100
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 5


def power_table(df, methods):
    """Rejection rate at ALPHA per (strength, method)."""
    d = df[df.method.isin(methods)].copy()
    d["reject"] = d.p_value < ALPHA
    return d.pivot_table(index="strength", columns="method", values="reject", aggfunc="mean")


def main():
    cfg = power_config(
        "lin", "lin", "binary_tree",
        reps=REPS,
        adaptive=["tree", "ordinal", "max", "euclid", "mGCM"],
    )
    print(f"Running Python: {cfg.name}  reps={cfg.reps}  strengths={cfg.strengths}  workers={WORKERS}")
    t0 = time.time()
    py = runner.run(cfg, workers=WORKERS, seed=0)
    elapsed = time.time() - t0
    n_rep = len(cfg.strengths) * cfg.reps

    methods = ["tree", "ordinal", "max", "euclid", "mGCM", "ankan", "chi_sq", "multinomial"]
    py_pow = power_table(py, methods)

    # R FIXED reference
    r_csv = REPO_ROOT / "results" / "power_lin_lin_binary_tree_FIXED.csv"
    r = pd.read_csv(r_csv).rename(columns={"chi": "chi_sq"})  # R names the pseudo-inverse test "chi"
    r_long = r.melt(id_vars=["strength"], value_vars=[m for m in methods if m in r.columns],
                    var_name="method", value_name="p_value")
    r_long["reject"] = r_long.p_value < ALPHA
    r_pow = r_long.pivot_table(index="strength", columns="method", values="reject", aggfunc="mean")

    pd.set_option("display.width", 200, "display.max_columns", 30)

    print("\n=== Python power (reject rate @0.05) ===")
    print(py_pow.round(3).to_string())
    print("\n=== R FIXED power (200 reps) ===")
    print(r_pow.round(3).to_string())

    print("\n=== tree - ordinal (the finding-#1 gap) ===")
    cmp = pd.DataFrame({
        "py_tree": py_pow["tree"], "py_ord": py_pow["ordinal"],
        "py_gap": py_pow["tree"] - py_pow["ordinal"],
        "R_tree": r_pow["tree"], "R_ord": r_pow["ordinal"],
        "R_gap": r_pow["tree"] - r_pow["ordinal"],
    })
    print(cmp.round(3).to_string())

    print("\n=== mean power across strengths ===")
    print(pd.DataFrame({"python": py_pow.mean(), "R_fixed": r_pow.mean()}).round(3).to_string())

    print("\n=== timing ===")
    print(f"python: {n_rep} reps in {elapsed/60:.1f} min ({WORKERS} workers) "
          f"= {elapsed/n_rep:.2f} s/rep wall, {elapsed*WORKERS/n_rep:.2f} s/rep-core")
    print("R FIXED equivalent (from timing log): 1800 reps in ~11.0 min (5 workers) "
          "= 0.37 s/rep wall, 1.83 s/rep-core")


if __name__ == "__main__":
    main()
