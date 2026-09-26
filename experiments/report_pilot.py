"""Combine the d-axis pilot blocks (A: ordinal search, B: search-free methods) into one table."""
import glob
import sys

import pandas as pd

from run_d_axis import report

pattern = sys.argv[1] if len(sys.argv) > 1 else "results/power_lin_lin_step_n2000_dx*_dy4_pilot*.parquet"
df = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(pattern))], ignore_index=True)
report(df)
