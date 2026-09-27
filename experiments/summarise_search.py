"""Summarise ``search_study.py`` output: rejection rates, and paired differences vs greedy.

Every method in a replicate sees the same data and bootstrap draws, so the
difference in rejection rate from greedy is estimated from the per-replicate
difference of the reject indicators -- much tighter than two independent SEs.

Usage:
    python summarise_search.py results/search_*.parquet
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

ALPHA = 0.05
ORDER = ["greedy", "beam5", "beam25", "random_fixed", "random_fresh", "split", "split_one",
         "depth0", "chi_sq"]


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    df = df.assign(reject=(df.p_value < ALPHA).astype(float))
    wide = df.pivot_table(index=["delta", "rep"], columns="method", values="reject")
    rows = []
    for delta, block in wide.groupby(level="delta"):
        n = len(block)
        for m in [c for c in ORDER if c in block.columns]:
            rate = block[m].mean()
            diff = block[m] - block["greedy"]
            rows.append(dict(delta=delta, method=m, reps=n, rate=rate,
                             se=np.sqrt(rate * (1 - rate) / n),
                             vs_greedy=diff.mean(), se_diff=diff.std(ddof=1) / np.sqrt(n)))
    return pd.DataFrame(rows)


def main():
    for path in sys.argv[1:]:
        df = pd.read_parquet(path)
        head = df.iloc[0]
        print(f"\n### {head.dgp} / {head.structure}  (d={head.d})   [{path}]")
        s = summarise(df)
        fmt = s.assign(rate=s.rate.map("{:.3f}".format), se=s.se.map("{:.3f}".format),
                       vs_greedy=s.vs_greedy.map("{:+.3f}".format),
                       se_diff=s.se_diff.map("{:.3f}".format))
        print(fmt.to_string(index=False))


if __name__ == "__main__":
    main()
