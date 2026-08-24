"""Head-to-head propensity benchmark: MLP vs XGBoost, on the paper's own quantity.

Assumption 1 of the paper controls type I error through

    E_f := max_j E[ (f_j(Z) - fhat_j(Z))^2 | D ],   E_f, E_g = o_P(1),
    E_f * E_g = o_P(n^{-1}),

so ``E_f`` -- not log-loss, not accuracy -- is the number a propensity learner
should be judged on here. The DGP hands back the true ``f``/``g``, so it can be
computed exactly rather than estimated.

Two versions of it are reported, because the pipeline dropped cross-fitting:

* ``E_f_in``  -- on the rows the learner was fitted on. This is what the test
  actually consumes, and it is optimistic for an overfitting learner.
* ``E_f_out`` -- on fresh rows from the same DGP. This is the honest
  generalisation error.

The gap between them is the overfitting the no-cross-fitting design exposes the
test to, so it is reported as its own column rather than left implicit. The
product column ``n * E_f * E_g`` is the Assumption 1 quantity itself: the theory
wants it heading to zero, and a value far above 1 says the rate condition is not
close to being met at this ``n``.

Every cell needs tuned params for *both* learners at that ``d`` (``tune.py
--write``), so the comparison stays tuned-vs-tuned as ``d`` grows.

A note on why strength 0 is enough: the interaction matrices have zero row and
column sums, so summing the joint over Y returns ``f`` exactly whatever the
strength is. The true X|Z propensity is the same under the null and the
alternative, and E_f measured here transfers to the power runs.

Usage:
    python bench_learners.py --reps 100 --workers 8
    python bench_learners.py --d 4 8 16 32 --settings lin sig --reps 50
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import dgp
from config import TUNING_ROOT

SETTINGS = ("lin", "vee", "hat", "sin", "sig")
N, NUM_CLASS, N_OOS = 1000, 8, 5000
INTSETTING = {"lin": "step", "vee": "step", "hat": "step", "sin": "binary_tree", "sig": "binary_tree"}
LEARNERS = ("xgb", "mlp")

RESULTS_DIR = Path(__file__).resolve().parent / "results"


def _E(true: np.ndarray, est: np.ndarray) -> float:
    """max_j mean_i (true_ij - est_ij)^2 -- the Assumption 1 remainder."""
    return float(np.max(np.mean((true - est) ** 2, axis=0)))


def _log_loss(labels: np.ndarray, proba: np.ndarray) -> float:
    p = np.clip(proba[np.arange(len(labels)), np.asarray(labels) - 1], 1e-15, None)
    return float(-np.mean(np.log(p)))


def _params(setting: str, learner: str, d: int) -> dict:
    path = TUNING_ROOT / f"n{N}_numclass{d}" / f"tune_{setting}_results.json"
    blob = json.loads(path.read_text())
    if learner not in blob:
        raise KeyError(f"No {learner!r} params in {path}; run tune.py --learner {learner} "
                       f"--d {d} --write.")
    return blob[learner]


def _one_rep(args) -> list[dict]:
    setting, d, seed_seq = args
    rng = np.random.default_rng(seed_seq)

    from catci.learners import mlp_learner, xgboost_learner

    build = {"xgb": xgboost_learner, "mlp": mlp_learner}

    # strength 0: the null, where size is decided and where Assumption 1 bites
    tr = dgp.simulate_data(N, d, d, setting, setting,
                           strength=0.0, intsetting=INTSETTING[setting], rng=rng)
    te = dgp.simulate_data(N_OOS, d, d, setting, setting,
                           strength=0.0, intsetting=INTSETTING[setting], rng=rng)

    rows = []
    for learner in LEARNERS:
        maker = build[learner](_params(setting, learner, d))
        t0 = time.perf_counter()
        predict = maker(tr["z"], tr["x"], d)
        f_in = predict(tr["z"])
        fit_s = time.perf_counter() - t0
        f_out = predict(te["z"])

        # g on the same data, so the Assumption 1 product is a real paired number
        predict_g = build[learner](_params(setting, learner, d))(tr["z"], tr["y"], d)
        g_in = predict_g(tr["z"])

        rows.append(dict(
            setting=setting, d=d, learner=learner,
            E_f_in=_E(tr["f"], f_in), E_f_out=_E(te["f"], f_out),
            E_g_in=_E(tr["g"], g_in),
            mse_f_in=float(np.mean((tr["f"] - f_in) ** 2)),
            mse_f_out=float(np.mean((te["f"] - f_out) ** 2)),
            logloss_in=_log_loss(tr["x"], f_in), logloss_out=_log_loss(te["x"], f_out),
            fit_seconds=fit_s,
        ))

    rows.append(dict(
        setting=setting, d=d, learner="oracle",
        E_f_in=0.0, E_f_out=0.0, E_g_in=0.0, mse_f_in=0.0, mse_f_out=0.0,
        logloss_in=_log_loss(tr["x"], tr["f"]), logloss_out=_log_loss(te["x"], te["f"]),
        fit_seconds=0.0,
    ))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--settings", nargs="*", default=list(SETTINGS))
    ap.add_argument("--d", nargs="*", type=int, default=[NUM_CLASS],
                    help="number of levels; each needs tuned params under tuning/nN_numclassD")
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench_learners", help="parquet stem under results/")
    args = ap.parse_args()

    import multiprocessing as mp
    import pandas as pd

    ss = np.random.SeedSequence(args.seed)
    tasks = [(s, d, child) for s in args.settings for d in args.d
             for child in ss.spawn(args.reps)]

    t0 = time.time()
    if args.workers == 1:
        out = [r for t in tasks for r in _one_rep(t)]
    else:
        with ProcessPoolExecutor(max_workers=args.workers,
                                 mp_context=mp.get_context("spawn")) as ex:
            out = [r for rows in ex.map(_one_rep, tasks, chunksize=2) for r in rows]
    elapsed = time.time() - t0

    df = pd.DataFrame(out)
    RESULTS_DIR.mkdir(exist_ok=True)
    df.to_parquet(RESULTS_DIR / f"{args.out}.parquet", index=False)

    agg = df.groupby(["d", "setting", "learner"]).mean(numeric_only=True)
    pd.set_option("display.width", 200, "display.float_format", lambda v: f"{v:.5f}")

    print(f"\n=== propensity quality, n={N}, {args.reps} reps/cell, {elapsed/60:.1f} min ===")
    print("E_f = max_j E[(f_j - fhat_j)^2]; 'in' = fitted rows (what the test uses), "
          "'out' = fresh rows\n")

    for d in args.d:
        for setting in args.settings:
            sub = agg.loc[(d, setting)]
            print(f"### d={d}  {setting}")
            print(f"    {'learner':8s} {'E_f_in':>9s} {'E_f_out':>9s} {'in/out':>7s} "
                  f"{'n*Ef*Eg':>9s} {'logloss_out':>12s} {'fit_s':>7s}")
            for learner in ["oracle"] + list(LEARNERS):
                if learner not in sub.index:
                    continue
                r = sub.loc[learner]
                ratio = r.E_f_in / r.E_f_out if r.E_f_out > 0 else float("nan")
                prod = N * r.E_f_in * r.E_g_in
                print(f"    {learner:8s} {r.E_f_in:9.5f} {r.E_f_out:9.5f} {ratio:7.2f} "
                      f"{prod:9.3f} {r.logloss_out:12.5f} {r.fit_seconds:7.3f}")
            print()

    print("Winner by E_f_out (lower is better):")
    for d in args.d:
        for setting in args.settings:
            sub = agg.loc[(d, setting)].drop(index="oracle", errors="ignore")
            best = sub.E_f_out.idxmin()
            ratio = sub.E_f_out.max() / sub.E_f_out.min()
            print(f"    d={d:<3d} {setting:5s} -> {best:4s}  ({ratio:.2f}x better)")


if __name__ == "__main__":
    main()
