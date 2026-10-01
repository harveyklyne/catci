"""Run Ankan & Textor's CI test on every adult-income variable pair given Z.

With ``Z = {Age, Sex}`` (the conditioning set of the paper's Fig. 1b), X and Y
range over the remaining nine variables, giving 36 pairs. Each is tested with
both conditional probability models the paper uses, GLM and RFT.

Output mirrors ``run.py``: a tidy parquet plus a JSON provenance sidecar, so a
table is reproducible from files rather than from an editing session.

Usage:
    python run_adult.py                          # full data, GLM + RFT
    python run_adult.py --models glm             # one model
    python run_adult.py --n 1000 --seed 6        # the Fig. 1 subsample setting
    python run_adult.py --keep-missing           # 32,561 rows, '?' as a level

Note on interpretation: the paper reports no p-values from its own method on
this dataset (Fig. 1b's p-values come from the stratified mutual-information
baseline it argues against; Figs. 8a/8b report a skeleton and F1 scores). These
numbers therefore have no published counterpart to be checked against. What *is*
checked is the preprocessing, via the Fig. 1b df arithmetic -- see
``adult.fig1b_df_check`` and ``tests/test_adult.py``.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import adult
from ankan_textor import GLM, RFT, design_matrix, fit_residuals, test_from_residuals
from config import REPO_ROOT

RESULTS_DIR = Path(__file__).resolve().parent / "results"
DEFAULT_Z = ("Age", "Sex")

MODELS = {"glm": GLM, "rft": RFT}


def run_pairs(data: adult.AdultData, z_names=DEFAULT_Z, models=("glm", "rft"),
              rf_seed=0, method="solve") -> pd.DataFrame:
    """Every (X, Y) pair over the variables not in ``z_names``, under each model."""
    z_names = tuple(z_names)
    z_codes = data.codes[list(z_names)].to_numpy()
    design = design_matrix(z_codes, [data.n_levels(v) for v in z_names])

    targets = [v for v in data.variables if v not in z_names]
    pairs = list(combinations(targets, 2))

    rows = []
    for key in models:
        model = MODELS[key]() if key == "glm" else MODELS[key](random_state=rf_seed)

        # A residual depends only on the variable and Z, and Z is the same for
        # every pair -- so fit each of the 9 variables once rather than once per
        # pair. With 8 pairs per variable this is an 8x saving, and the GLM fits
        # are what dominate the sweep (NativeCountry alone is ~65s at full n).
        residuals = {}
        for name in targets:
            t0 = time.time()
            residuals[name] = fit_residuals(
                data.codes[name].to_numpy(), design,
                data.kind(name), data.n_levels(name), model,
            )
            print(f"  [fit] {key:4s} {name:14s} k={data.n_levels(name):3d} "
                  f"{time.time() - t0:6.1f}s", flush=True)

        for x_name, y_name in pairs:
            t0 = time.time()
            try:
                res = test_from_residuals(
                    residuals[x_name], residuals[y_name],
                    data.kind(x_name), data.kind(y_name),
                    name_x=x_name, name_y=y_name, method=method,
                )
                rows.append(dict(
                    model=key, X=x_name, Y=y_name,
                    X_code=data.variables[x_name].code, Y_code=data.variables[y_name].code,
                    kind_X=data.kind(x_name), kind_Y=data.kind(y_name),
                    k=data.n_levels(x_name), r=data.n_levels(y_name),
                    statistic=res.statistic, df=res.df, p_value=res.p_value,
                    which=res.which, rank=res.rank, cond=res.cond,
                    well_conditioned=res.well_conditioned,
                    n=data.n, seconds=round(time.time() - t0, 2),
                ))
            except Exception as e:  # keep the sweep going; record the failure
                rows.append(dict(
                    model=key, X=x_name, Y=y_name, which="__error__",
                    statistic=float("nan"), df=-1, p_value=float("nan"),
                    rank=-1, cond=float("nan"), well_conditioned=False,
                    n=data.n, error=str(e), seconds=round(time.time() - t0, 2),
                ))
            flag = "" if rows[-1]["well_conditioned"] else "   <-- ill-conditioned, undefined"
            print(f"  {key:4s} {x_name:14s} {y_name:14s} {rows[-1].get('which'):9s} "
                  f"p={rows[-1]['p_value']:.4g} df={rows[-1]['df']}{flag}", flush=True)
    return pd.DataFrame(rows)


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def _write_provenance(df, out_parquet, meta, elapsed):
    import scipy, sklearn, statsmodels
    prov = {
        **meta,
        "git_commit": _git_commit(),
        "runtime_seconds": round(elapsed, 1),
        "n_tests": len(df),
        "n_errors": int((df["which"] == "__error__").sum()),
        "python": platform.python_version(),
        "packages": {
            "numpy": np.__version__, "scipy": scipy.__version__,
            "pandas": pd.__version__, "scikit-learn": sklearn.__version__,
            "statsmodels": statsmodels.__version__,
        },
        "platform": platform.platform(),
    }
    out_parquet.with_suffix(".provenance.json").write_text(json.dumps(prov, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="glm,rft", help="comma-separated: glm,rft")
    ap.add_argument("--z", default=",".join(DEFAULT_Z), help="comma-separated conditioning variables")
    ap.add_argument("--n", type=int, default=None, help="subsample size (default: all rows)")
    ap.add_argument("--seed", type=int, default=0, help="subsample seed, and the RF seed")
    ap.add_argument("--keep-missing", action="store_true", help="keep '?' as a category")
    ap.add_argument("--method", default="solve", choices=("solve", "pinv"),
                    help="solve = the paper's form; pinv = pseudo-inverse with df = rank(Sigma_d)")
    ap.add_argument("--tag", default=None, help="output name suffix")
    args = ap.parse_args()

    models = tuple(m.strip() for m in args.models.split(","))
    z_names = tuple(v.strip() for v in args.z.split(","))
    unknown = set(models) - set(MODELS)
    if unknown:
        ap.error(f"unknown model(s): {sorted(unknown)}; choose from {sorted(MODELS)}")

    data = adult.load(drop_missing=not args.keep_missing)
    if args.n is not None:
        data = data.subsample(args.n, seed=args.seed)

    RESULTS_DIR.mkdir(exist_ok=True)
    name = args.tag or (
        f"adult_n{data.n}"
        + ("_keepmissing" if args.keep_missing else "")
        + ("" if args.method == "solve" else f"_{args.method}")
    )
    out_parquet = RESULTS_DIR / f"{name}.parquet"

    print(f"adult income: n = {data.n}, Z = {z_names}, models = {models}")
    t0 = time.time()
    df = run_pairs(data, z_names=z_names, models=models, rf_seed=args.seed, method=args.method)
    elapsed = time.time() - t0

    df.to_parquet(out_parquet, index=False)
    _write_provenance(df, out_parquet, dict(
        dataset="uci-adult-train", n=data.n, z=list(z_names), models=list(models),
        drop_missing=not args.keep_missing, subsample_n=args.n, seed=args.seed,
        method=args.method,
        levels={v: data.n_levels(v) for v in data.variables},
        kinds={v: data.kind(v) for v in data.variables},
    ), elapsed)

    ok = df["well_conditioned"].sum()
    print(f"\nDONE {name}: {len(df)} tests in {elapsed/60:.1f} min "
          f"({ok} well-conditioned, {len(df) - ok} not) -> {out_parquet}")
    return df


if __name__ == "__main__":
    main()
