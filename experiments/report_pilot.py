"""Combine the d-axis pilot blocks (A: ordinal search, B: search-free methods) into one table.

The default pattern matches both learner-tagged names (``..._dy4_pilotA__oracle``)
and the untagged ones written before results carried the learner; ``report``
prints one table per learner.
"""
import glob
import sys

import pandas as pd

from run_d_axis import report

DEFAULT_PATTERN = "results/power_lin_lin_step_n2000_dx*_dy4_pilot*.parquet"


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATTERN
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise SystemExit(f"no parquets match {pattern}")
    report(pd.concat([pd.read_parquet(p) for p in paths], ignore_index=True))


if __name__ == "__main__":
    main()
