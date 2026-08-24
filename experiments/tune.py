"""Tune a propensity learner, on the protocol the frozen XGBoost params came from.

The frozen ``tuning/*.json`` came from the R ``tune_xgb``: simulate a fresh
train/test pair per rep, fit every grid point on the train half, score
**held-out multiclass log-loss** on the test half, average over reps, take the
argmin. This module reproduces that protocol for *both* learners, so ``mlp`` and
``xgb`` are compared tuned-vs-tuned rather than tuned-vs-default, and so a new
``(n, d)`` can be tuned at all -- the R script that produced the ``d = 8`` JSONs
was cluster-specific and was never ported, which is the README's "No tuner" gap.

Reproducing the R numbers at ``d = 8`` is the check that the port is faithful;
see ``--learner xgb --grid r``, whose grid is the R one exactly.

Usage:
    python tune.py --learner mlp --reps 20                 # all five settings
    python tune.py --learner xgb --reps 20 --d 30 --write   # a new (n, d)
    python tune.py sig --learner mlp --reps 5 --grid pilot  # one setting
"""

from __future__ import annotations

import argparse
import itertools
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import dgp
from config import TUNING_ROOT

SETTINGS = ("lin", "vee", "hat", "sin", "sig")
N_TR, N_TE, STRENGTH = 800, 5000, 0.5

# The R tuner paired each marginal setting with the interaction it is used with.
INTSETTING = {"lin": "step", "vee": "step", "hat": "step", "sin": "binary_tree", "sig": "binary_tree"}

MLP_GRIDS = {
    # A wide first pass: is the MLP even in the running, and where does it live?
    "pilot": dict(
        hidden_layer_sizes=[(8,), (32,), (8, 8), (32, 32)],
        alpha=[1e-2, 1e-1, 1.0, 10.0],
        activation=["tanh", "relu"],
        max_iter=[400],
    ),
    # The tuning grid proper. The pilot railed at the top of its alpha range on
    # every setting, so alpha runs well past it here -- these propensities carry
    # little signal about Z, and the tuner wants heavy shrinkage.
    "full": dict(
        hidden_layer_sizes=[(4,), (8,), (32,), (8, 8)],
        alpha=[1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0],
        activation=["tanh"],
        max_iter=[400],
    ),
}

# eta and the nrounds ceiling are fixed; nrounds itself is read off the boosting
# path rather than gridded, exactly as the R tuner did.
XGB_GRIDS = {
    "r": dict(max_depth=[1, 2, 3, 4], gamma=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]),
    "pilot": dict(max_depth=[1, 2], gamma=[0.0, 1.0, 2.0]),
}
XGB_ETA, XGB_MAXROUNDS = 0.01, 1000


def grid_points(grid: dict) -> list[dict]:
    keys = list(grid)
    return [dict(zip(keys, vals)) for vals in itertools.product(*(grid[k] for k in keys))]


def _log_loss(true_labels: np.ndarray, proba: np.ndarray) -> float:
    """Multiclass log-loss, matching xgboost's ``mlogloss`` (natural log, clipped)."""
    p = np.clip(proba[np.arange(len(true_labels)), np.asarray(true_labels) - 1], 1e-15, None)
    return float(-np.mean(np.log(p)))


def _simulate_pair(setting: str, d: int, rng):
    tr = dgp.simulate_data(N_TR, d, d, setting, setting,
                           strength=STRENGTH, intsetting=INTSETTING[setting], rng=rng)
    te = dgp.simulate_data(N_TE, d, d, setting, setting,
                           strength=STRENGTH, intsetting=INTSETTING[setting], rng=rng)
    return tr, te


# --------------------------------------------------------------------------- #
# MLP: one fit per grid point
# --------------------------------------------------------------------------- #
def _mlp_rep(args) -> np.ndarray:
    setting, d, points, seed_seq = args
    from catci.learners import mlp_learner

    rng = np.random.default_rng(seed_seq)
    tr, te = _simulate_pair(setting, d, rng)

    out = np.empty(len(points))
    for i, params in enumerate(points):
        try:
            predict = mlp_learner(params)(tr["z"], tr["x"], d)
            out[i] = _log_loss(te["x"], predict(te["z"]))
        except Exception:
            out[i] = np.inf
    return out


# --------------------------------------------------------------------------- #
# XGB: one fit per (depth, gamma), with nrounds read off the boosting path
# --------------------------------------------------------------------------- #
def _xgb_rep(args) -> np.ndarray:
    """Held-out mlogloss for every (depth, gamma, nrounds), flattened.

    Boosting is sequential, so the whole ``nrounds`` axis comes free from a
    single ``xgb.train`` with an eval set -- training once per (depth, gamma) and
    reading ``evals_result`` is what makes a 1000-long nrounds grid affordable.
    This is the R tuner's ``fit$evaluation_log$test_mlogloss`` trick.
    """
    setting, d, points, seed_seq = args
    import xgboost as xgb

    rng = np.random.default_rng(seed_seq)
    tr, te = _simulate_pair(setting, d, rng)
    dtrain = xgb.DMatrix(tr["z"], label=tr["x"] - 1)
    dtest = xgb.DMatrix(te["z"], label=te["x"] - 1)

    out = np.empty((len(points), XGB_MAXROUNDS))
    for i, p in enumerate(points):
        evals: dict = {}
        xgb.train(
            {"eta": XGB_ETA, "max_depth": int(p["max_depth"]), "gamma": float(p["gamma"]),
             "objective": "multi:softprob", "eval_metric": "mlogloss",
             "num_class": d, "nthread": 1},
            dtrain, num_boost_round=XGB_MAXROUNDS,
            evals=[(dtest, "test")], evals_result=evals, verbose_eval=False,
        )
        out[i] = evals["test"]["mlogloss"]
    return out.ravel()


# --------------------------------------------------------------------------- #
def _reference_rep(args):
    """Held-out mlogloss for the anchors that make a grid number readable.

    The true propensities are the floor no learner can beat; the uniform 1/d
    predictor is the score a learner shrunk to nothing achieves. A grid winner
    between them is doing something; one at or above uniform is not.
    """
    setting, d, seed_seq = args
    rng = np.random.default_rng(seed_seq)
    _, te = _simulate_pair(setting, d, rng)
    oracle = _log_loss(te["x"], te["f"])
    unif = _log_loss(te["x"], np.full((N_TE, d), 1.0 / d))
    return oracle, unif


def _map(fn, tasks, workers):
    if workers == 1:
        return [fn(t) for t in tasks]
    import multiprocessing as mp

    # "spawn", not "fork": sklearn/BLAS start threads, and a forked child that
    # inherits them aborts in the macOS Objective-C runtime on the second pool.
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        return list(ex.map(fn, tasks))


def tune(learner: str, setting: str, d: int, points, reps: int, workers: int, seed: int = 0):
    fn = {"mlp": _mlp_rep, "xgb": _xgb_rep}[learner]
    ss = np.random.SeedSequence(seed)
    tasks = [(setting, d, points, child) for child in ss.spawn(reps)]

    t0 = time.time()
    mean = np.mean(np.vstack(_map(fn, tasks, workers)), axis=0)
    elapsed = time.time() - t0

    ref_tasks = [(setting, d, child) for child in np.random.SeedSequence(seed).spawn(reps)]
    oracle, unif = np.mean(np.array(_map(_reference_rep, ref_tasks, workers)), axis=0)
    return mean, {"oracle": float(oracle), "uniform": float(unif)}, elapsed


def describe(learner: str, points, mean: np.ndarray, index: int) -> tuple[dict, str]:
    """The grid point at flat position ``index``, as params and as a printable line."""
    if learner == "mlp":
        p = points[index]
        return (dict(p, hidden_layer_sizes=list(p["hidden_layer_sizes"])),
                f"hidden={str(p['hidden_layer_sizes']):9s} alpha={p['alpha']:<7g} "
                f"act={p['activation']:5s} iter={p['max_iter']}")
    i, r = divmod(index, XGB_MAXROUNDS)
    p = points[i]
    params = {"eta": XGB_ETA, "max.depth": int(p["max_depth"]),
              "gamma": float(p["gamma"]), "nrounds": int(r) + 1}
    return params, (f"max.depth={params['max.depth']}  gamma={params['gamma']:<4g} "
                    f"nrounds={params['nrounds']}")


def _boundary_warnings(learner: str, points, best: dict) -> list[str]:
    """Flag a winner sitting on the edge of the grid.

    An argmin at a boundary means the search was truncated, not that the boundary
    is optimal -- either the range needs extending or the mlogloss curve is still
    too noisy to have a real interior minimum (too few reps). Both happened while
    tuning here, so the check is worth printing rather than remembering.
    """
    out = []
    if learner == "xgb":
        if best["nrounds"] >= XGB_MAXROUNDS:
            out.append(f"nrounds hit the ceiling ({XGB_MAXROUNDS}); raise it or add reps.")
        depths = sorted({p["max_depth"] for p in points})
        gammas = sorted({p["gamma"] for p in points})
        if best["max.depth"] == depths[-1]:
            out.append(f"max.depth hit the top of the grid ({depths[-1]}).")
        if best["gamma"] in (gammas[0], gammas[-1]):
            out.append(f"gamma hit the edge of the grid ({best['gamma']}).")
    else:
        alphas = sorted({p["alpha"] for p in points})
        if best["alpha"] in (alphas[0], alphas[-1]):
            out.append(f"alpha hit the edge of the grid ({best['alpha']}).")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("settings", nargs="*", default=None, help=f"subset of {SETTINGS}")
    ap.add_argument("--learner", choices=["mlp", "xgb"], required=True)
    ap.add_argument("--grid", default=None, help="mlp: pilot|full;  xgb: r|pilot")
    ap.add_argument("--n", type=int, default=1000, help="n the params are filed under")
    ap.add_argument("--d", type=int, default=8)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--write", action="store_true", help="merge the winner into tuning/*.json")
    ap.add_argument("--top", type=int, default=6)
    args = ap.parse_args()

    grids = MLP_GRIDS if args.learner == "mlp" else XGB_GRIDS
    grid_name = args.grid or ("full" if args.learner == "mlp" else "r")
    points = grid_points(grids[grid_name])
    n_flat = len(points) * (1 if args.learner == "mlp" else XGB_MAXROUNDS)

    settings = args.settings or list(SETTINGS)
    print(f"{args.learner} grid '{grid_name}': {n_flat} points x {args.reps} reps, "
          f"n={args.n} d={args.d}, {len(settings)} settings")

    for setting in settings:
        mean, ref, elapsed = tune(args.learner, setting, args.d, points,
                                  args.reps, args.workers, args.seed)
        best_i = int(np.argmin(mean))
        best_params, _ = describe(args.learner, points, mean, best_i)

        print(f"\n### {setting}  ({args.reps} reps, {elapsed/60:.1f} min)")
        print(f"    reference mlogloss: oracle={ref['oracle']:.5f}  uniform={ref['uniform']:.5f}")
        for i in np.argsort(mean)[: args.top]:
            _, line = describe(args.learner, points, mean, int(i))
            star = " <-- best" if i == best_i else ""
            print(f"    mlogloss={mean[i]:.5f}  {line}{star}")

        for warning in _boundary_warnings(args.learner, points, best_params):
            print(f"    WARNING: {warning}")

        if args.write:
            root = TUNING_ROOT / f"n{args.n}_numclass{args.d}"
            root.mkdir(parents=True, exist_ok=True)
            path = root / f"tune_{setting}_results.json"
            blob = json.loads(path.read_text()) if path.exists() else {}
            blob[args.learner] = best_params
            blob[f"{args.learner}_tuning"] = {
                "grid": grid_name, "reps": args.reps, "seed": args.seed,
                "n_tr": N_TR, "n_te": N_TE, "d": args.d,
                "metric": "held-out mlogloss",
                "best_mlogloss": round(float(mean.min()), 6),
                "reference_mlogloss": {k: round(v, 6) for k, v in ref.items()},
            }
            path.write_text(json.dumps(blob, indent=2) + "\n")
            print(f"    wrote {path}")


if __name__ == "__main__":
    main()
