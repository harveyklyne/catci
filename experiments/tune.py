"""Tune the XGBoost propensity learner for a simulated ``(n, d, setting)``.

Writes ``tuning/n{n}_numclass{d}/tune_{setting}_results.json`` -- the file
:meth:`config.Config.xgb_params` reads -- with the chosen ``xgb`` parameters
plus the whole CV table and a provenance block.

The tuning data are drawn from :func:`dgp.simulate_marginal`, the exact law of
``X | Z`` in every simulation (margins are preserved by construction), so one
tuning per ``(n, num_class, setting)`` serves every interaction, strength and
partner dimension. ``Y`` uses the same file keyed by ``ysetting`` and ``dy``.

Usage:
    python tune.py --n 1000 --d 8 --settings sin,lin [--reps 20] [--grid r|fast]
                   [--max-rounds 1000] [--workers 5] [--seed 0]

Grid points of one setting run serially; settings (and ``d`` values) are spread
over ``--workers`` processes, each with a single xgboost thread.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import dgp
from catci.tuning import FAST_GRID, R_GRID, tune_xgboost
from config import TUNING_ROOT, tuning_path

GRIDS = {"r": R_GRID, "fast": FAST_GRID}


def tune_setting(n, d, setting, reps, grid_name, max_rounds, early_stopping_rounds, seed, out_root=TUNING_ROOT):
    rng = np.random.default_rng([seed, n, d, sum(map(ord, setting))])
    datasets = []
    for _ in range(reps):
        sim = dgp.simulate_marginal(n, d, setting, rng)
        datasets.append((sim["z"], sim["x"]))
    res = tune_xgboost(
        datasets, d, grid=GRIDS[grid_name], n_folds=5, max_rounds=max_rounds,
        early_stopping_rounds=early_stopping_rounds, nthread=1, rng=rng,
    )
    out = res.to_json()
    import xgboost
    out["provenance"] = dict(
        n=n, num_class=d, setting=setting, reps=reps, n_folds=5, grid=grid_name,
        max_rounds=max_rounds, early_stopping_rounds=early_stopping_rounds, seed=seed,
        xgboost=xgboost.__version__, python=platform.python_version(),
    )
    path = tuning_path(n, d, setting, root=out_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    print(f"DONE n={n} d={d} {setting}: {res.params}  cv={res.cv_logloss:.4f}  "
          f"{res.seconds/60:.1f} min -> {path}", flush=True)
    return path


def _task(args):
    return tune_setting(*args)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--d", type=str, required=True, help="num_class; comma-separated for several")
    ap.add_argument("--settings", type=str, default="sin,sig,lin,vee,hat")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--grid", choices=sorted(GRIDS), default="r")
    ap.add_argument("--max-rounds", type=int, default=1000)
    ap.add_argument("--early-stopping", type=int, default=50, help="0 disables")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default=str(TUNING_ROOT), help="tuning root (default: experiments/tuning)")
    args = ap.parse_args()

    esr = args.early_stopping or None
    tasks = [
        (args.n, d, s, args.reps, args.grid, args.max_rounds, esr, args.seed, Path(args.out))
        for d in (int(v) for v in args.d.split(","))
        for s in args.settings.split(",")
    ]
    t0 = time.time()
    if args.workers == 1 or len(tasks) == 1:
        for t in tasks:
            _task(t)
    else:
        import multiprocessing as mp
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("fork")) as ex:
            list(ex.map(_task, tasks))
    print(f"all done in {(time.time() - t0)/60:.1f} min")


if __name__ == "__main__":
    main()
