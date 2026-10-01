"""Summarise ``results_semisynth/``: size, power, paired differences and CPU cost.

Reads the sweep (``sweep_*.parquet`` + ``*_timing.parquet``, oracle, depth x B) and
the learner runs (``run_semisynth.py`` output) and prints markdown tables.

Usage:
    python report_semisynth.py                     # everything found
    python report_semisynth.py --pattern 'sweep_Education_Income*'
"""

from __future__ import annotations

import argparse
import re

import numpy as np
import pandas as pd

from run_semisynth import OUT

ALPHA = 0.05


def _md(df: pd.DataFrame, digits: int = 3) -> str:
    """A markdown table (no ``tabulate`` dependency)."""
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else f"{v:.{digits}f}"
        return str(v)
    names = df.index.names
    index_head = " / ".join(str(n) for n in names if n is not None) or ""
    head = [index_head] + [str(c) for c in df.columns]
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for idx, row in df.iterrows():
        idx = " / ".join(map(str, idx)) if isinstance(idx, tuple) else str(idx)
        lines.append("| " + " | ".join([idx] + [fmt(v) for v in row.to_numpy()]) + " |")
    return "\n".join(lines)


def _label(stem: str) -> str:
    return re.sub(r"_r\d+.*$", "", stem.removeprefix("sweep_"))


def load(pattern: str) -> dict[str, pd.DataFrame]:
    out = {}
    for path in sorted(OUT.glob(f"{pattern}.parquet")):
        if path.stem.endswith("_timing"):
            continue
        out[path.stem] = pd.read_parquet(path)
    return out


def sweep_report(stem: str, df: pd.DataFrame) -> None:
    df = df.assign(reject=df.p_value < ALPHA, m=df.method + "@" + df.k.astype(str))
    order = list(dict.fromkeys(df.m))
    for lam, block in df.groupby("lam"):
        reps = block.rep.nunique()
        kind = "size" if lam == 0 else "power"
        se = np.sqrt((ALPHA if lam == 0 else 0.25) * (1 - (ALPHA if lam == 0 else 0.5)) / reps)
        print(f"\n### {_label(stem)} -- lam = {lam} ({kind} at {ALPHA}; {reps} reps; SE <= {se:.3f})\n")
        tab = block.pivot_table(index="m", columns="n_boot", values="reject", sort=False)
        print(_md(tab.loc[[m for m in order if m in tab.index]]))
        if lam > 0:
            # Paired difference vs merge at the largest B: mean and SE of the
            # per-replicate difference of reject indicators.
            B = block.n_boot.max()
            wide = block[block.n_boot == B].pivot_table(index="rep", columns="m", values="reject")
            base = [c for c in wide.columns if c.startswith("merge@")][0]
            diff = wide.sub(wide[base], axis=0)
            summ = pd.DataFrame({"vs merge": diff.mean(), "se": diff.std(ddof=1) / np.sqrt(len(diff))})
            print(f"\npaired difference vs {base} at B = {B}:\n")
            print(_md(summ.loc[[m for m in order if m in summ.index]]))


def timing_report(stem: str) -> None:
    path = OUT / f"{stem}_timing.parquet"
    if not path.exists():
        return
    t = pd.read_parquet(path)
    direct = t[t.method.str.endswith("_direct")]
    if direct.empty:
        return
    direct = direct.assign(m=direct.method.str.removesuffix("_direct") + "@" + direct.k.astype(str))
    print(f"\n### {_label(stem)} -- CPU seconds per search (median over "
          f"{direct.rep.nunique()} reps), rows method@k, columns B\n")
    print(_md(direct.pivot_table(index="m", columns="n_boot", values="cpu", aggfunc="median",
                                 sort=False), 4))
    # Per-draw cost and the exponent of cost in B (log-log slope), per method.
    rows = []
    for m, g in direct.groupby("m", sort=False):
        med = g.groupby("n_boot").cpu.median()
        slope = np.polyfit(np.log(med.index), np.log(med.values), 1)[0] if len(med) > 1 else np.nan
        rows.append(dict(m=m, ms_per_draw_at_maxB=1e3 * med.iloc[-1] / med.index[-1], B_exponent=slope))
    print("\n" + _md(pd.DataFrame(rows).set_index("m"), 3))


def learner_report(stem: str, df: pd.DataFrame) -> None:
    df = df.assign(reject=df.p_value < ALPHA)
    order = list(dict.fromkeys(df.method))
    print(f"\n### {stem}\n")
    tab = df.pivot_table(index=["learner", "lam"], columns="method", values="reject", sort=False)
    print(_md(tab[[m for m in order if m in tab.columns]]))
    cost = df.groupby(["learner", "method"], sort=False)[["search_cpu", "cal_cpu", "fit_cpu"]].mean()
    print("\nmean CPU seconds per test:\n")
    print(_md(cost, 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="*")
    args = ap.parse_args()
    for stem, df in load(args.pattern).items():
        if stem.startswith("sweep_"):
            sweep_report(stem, df)
            timing_report(stem)
        elif "search_cpu" in df.columns:
            learner_report(stem, df)


if __name__ == "__main__":
    main()
