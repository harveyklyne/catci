"""Runner: a Config -> one tidy long-format table + a provenance snapshot.

One row per (setting, method, strength, rep) with its p-value (CODE_REVIEW.md
4.2). Emits parquet plus a JSON sidecar recording the resolved config, git
commit, package versions, seed and runtime -- so a figure is reproducible from
files, not from an editing session.

Usage:
    python run.py <config-name> [--reps N] [--strengths a,b,c] [--workers W] [--seed S]

<config-name> is "<x>_<y>_<int>", e.g. lin_lin_binary_tree.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

import dgp
import methods
from config import Config, REPO_ROOT, power_config
from catci.learners import crossfit, xgboost_learner

RESULTS_DIR = Path(__file__).resolve().parent / "results"


def _one_replicate(cfg: Config, strength: float, rep: int, seed_seq) -> list[dict]:
    """Simulate, cross-fit propensities, and evaluate every requested method."""
    rng = np.random.default_rng(seed_seq)

    data = dgp.simulate_data(
        cfg.n, cfg.d, cfg.d, cfg.xsetting, cfg.ysetting,
        strength=strength, intsetting=cfg.intsetting, permute=False, rng=rng,
    )

    # cross-fit f (X|Z) and g (Y|Z) with the per-setting tuned learners
    f = crossfit(data["z"], data["x"], cfg.d, xgboost_learner(cfg.xgb_params(cfg.xsetting)), cfg.nfolds, rng)
    g = crossfit(data["z"], data["y"], cfg.d, xgboost_learner(cfg.xgb_params(cfg.ysetting)), cfg.nfolds, rng)

    fitted = methods.Fitted.build(data["x"], data["y"], data["z"], f, g, cfg.d, cfg.d, cfg.normalise)

    pvals = {}
    pvals.update(methods.adaptive_pvalues(fitted, cfg.adaptive, cfg.n_boot, rng))
    pvals.update(methods.competitor_pvalues(fitted, cfg.competitors, rng))

    base = dict(name=cfg.name, n=cfg.n, d=cfg.d, xsetting=cfg.xsetting,
                ysetting=cfg.ysetting, intsetting=cfg.intsetting, strength=strength, rep=rep)
    return [dict(base, method=m, p_value=p) for m, p in pvals.items()]


def _task(args):
    cfg, strength, rep, seed_seq = args
    try:
        return _one_replicate(cfg, strength, rep, seed_seq)
    except Exception as e:  # keep the grid running; record the failure as a row
        return [dict(name=cfg.name, strength=strength, rep=rep, method="__error__", p_value=float("nan"), error=str(e))]


def run(cfg: Config, workers: int = 5, seed: int = 0) -> pd.DataFrame:
    RESULTS_DIR.mkdir(exist_ok=True)
    tasks = []
    ss = np.random.SeedSequence(seed)
    # one independent child seed per (strength, rep), spawned deterministically
    children = ss.spawn(len(cfg.strengths) * cfg.reps)
    k = 0
    for strength in cfg.strengths:
        for rep in range(1, cfg.reps + 1):
            tasks.append((cfg, strength, rep, children[k]))
            k += 1

    t0 = time.time()
    rows = []
    if workers == 1:
        for t in tasks:
            rows.extend(_task(t))
    else:
        import multiprocessing as mp
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
            for res in ex.map(_task, tasks, chunksize=4):
                rows.extend(res)
    elapsed = time.time() - t0

    df = pd.DataFrame(rows)
    out_parquet = RESULTS_DIR / f"{cfg.name}.parquet"
    df.to_parquet(out_parquet, index=False)
    _write_provenance(cfg, df, elapsed, workers, seed, out_parquet)
    print(f"DONE {cfg.name}: {len(tasks)} replicates, {len(df)} rows, {elapsed/60:.1f} min "
          f"({workers} workers) -> {out_parquet}")
    return df


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def _write_provenance(cfg, df, elapsed, workers, seed, out_parquet):
    import numpy, scipy, xgboost, sklearn
    n_err = int((df.get("method") == "__error__").sum()) if "method" in df else 0
    prov = {
        "config": cfg.to_dict(),
        "git_commit": _git_commit(),
        "seed": seed,
        "workers": workers,
        "runtime_seconds": round(elapsed, 1),
        "n_replicates": len(cfg.strengths) * cfg.reps,
        "n_errors": n_err,
        "python": platform.python_version(),
        "packages": {"numpy": numpy.__version__, "scipy": scipy.__version__,
                     "xgboost": xgboost.__version__, "scikit-learn": sklearn.__version__},
        "platform": platform.platform(),
    }
    (out_parquet.with_suffix(".provenance.json")).write_text(json.dumps(prov, indent=2))


def _parse_config(name: str, reps=None, strengths=None) -> Config:
    parts = name.split("_")
    # intsetting may itself contain '_' (binary_tree)
    xsetting, ysetting, intsetting = parts[0], parts[1], "_".join(parts[2:])
    cfg = power_config(xsetting, ysetting, intsetting)
    over = {}
    if reps is not None:
        over["reps"] = reps
    if strengths is not None:
        over["strengths"] = strengths
    return replace(cfg, **over) if over else cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--strengths", type=str, default=None, help="comma-separated, e.g. 0.6,1.0,1.4")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    strengths = [float(s) for s in args.strengths.split(",")] if args.strengths else None
    cfg = _parse_config(args.config, reps=args.reps, strengths=strengths)
    run(cfg, workers=args.workers, seed=args.seed)


if __name__ == "__main__":
    main()
