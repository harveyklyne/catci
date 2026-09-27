"""Tune a propensity learner for a simulated ``(n, d, setting)``.

With ``--write`` the winner is merged into
``tuning/n{n}_numclass{d}/tune_{setting}_results.json`` -- the file
:meth:`config.Config.learner_params` reads -- under the learner's key (``mlp``
or ``xgb``), plus a ``<learner>_tuning`` provenance block. The other learner's
entry in the same file is left alone.

Two protocols:

``--protocol holdout`` (default; either learner)
    The one the frozen ``d = 8`` JSONs came from, the R ``tune_xgb``: simulate a
    fresh train/test pair per rep, fit every grid point on the train half, score
    **held-out multiclass log-loss** on the test half, average over reps, take
    the argmin. Running it for *both* learners means ``mlp`` and ``xgb`` are
    compared tuned-vs-tuned rather than tuned-vs-default. Reproducing the R
    numbers at ``d = 8`` is the check that the port is faithful; see
    ``--learner xgb --grid r``, whose grid is the R one exactly. ``--n`` is only
    the label the result is filed under (the train size is ``N_TR``).

``--protocol cv`` (xgb only)
    :func:`catci.tuning.tune_xgboost`: K-fold CV log-loss on ``--reps`` draws of
    size ``n`` from :func:`dgp.simulate_marginal`, the exact law of ``X | Z`` in
    every simulation (margins are preserved by construction), so one tuning per
    ``(n, num_class, setting)`` serves every interaction, strength and partner
    dimension. ``Y`` uses the same file keyed by ``ysetting`` and ``dy``. With
    ``--grid fast`` this is the cheap route to a new ``(n, d)`` on the ``d`` axis.
    ``(d, setting)`` cells are spread over ``--workers`` processes, each with a
    single xgboost thread.

Usage:
    python tune.py --learner mlp --reps 20                  # all five settings
    python tune.py --learner xgb --reps 20 --d 30 --write    # a new (n, d)
    python tune.py sig --learner mlp --reps 5 --grid pilot   # one setting
    python tune.py --learner xgb --protocol cv --n 2000 --d 16,32 --settings lin \\
        --grid fast --write [--max-rounds 1000] [--early-stopping 50]
"""

from __future__ import annotations

import argparse
import itertools
import json
import platform
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import dgp
from config import TUNING_ROOT, tuning_path

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


# --------------------------------------------------------------------------- #
# CV protocol (xgb only): K-fold CV on draws from the X | Z margin
# --------------------------------------------------------------------------- #
def _cv_grids() -> dict:
    from catci.tuning import FAST_GRID, R_GRID

    return {"r": R_GRID, "fast": FAST_GRID}


def tune_setting_cv(n, d, setting, reps, grid_name, max_rounds, early_stopping_rounds, seed):
    """K-fold-CV xgb tuning for one ``(n, d, setting)``: ``(params, provenance)``."""
    import xgboost

    from catci.tuning import tune_xgboost

    rng = np.random.default_rng([seed, n, d, sum(map(ord, setting))])
    datasets = []
    for _ in range(reps):
        sim = dgp.simulate_marginal(n, d, setting, rng)
        datasets.append((sim["z"], sim["x"]))
    res = tune_xgboost(
        datasets, d, grid=_cv_grids()[grid_name], n_folds=5, max_rounds=max_rounds,
        early_stopping_rounds=early_stopping_rounds, nthread=1, rng=rng,
    )
    out = res.to_json()
    provenance = dict(
        protocol="cv", n=n, num_class=d, setting=setting, reps=reps, n_folds=5, grid=grid_name,
        max_rounds=max_rounds, early_stopping_rounds=early_stopping_rounds, seed=seed,
        metric="K-fold CV mlogloss", xgboost=xgboost.__version__, python=platform.python_version(),
        cv=out["cv"],
    )
    print(f"DONE n={n} d={d} {setting}: {res.params}  cv={res.cv_logloss:.4f}  "
          f"{res.seconds/60:.1f} min", flush=True)
    return out["xgb"], provenance


def _cv_task(args):
    return tune_setting_cv(*args)


def write_winner(root, n: int, d: int, setting: str, learner: str, params: dict, tuning: dict) -> Path:
    """Merge one learner's winner into its JSON, leaving the other learner's entry alone."""
    path = tuning_path(n, d, setting, root=Path(root))
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = json.loads(path.read_text()) if path.exists() else {}
    blob[learner] = params
    blob[f"{learner}_tuning"] = tuning
    path.write_text(json.dumps(blob, indent=2) + "\n")
    print(f"    wrote {path}")
    return path


def main_cv(args, settings, ds):
    if args.learner != "xgb":
        raise SystemExit("--protocol cv tunes xgb only")
    grid_name = args.grid or "r"
    if grid_name not in _cv_grids():
        raise SystemExit(f"--protocol cv takes --grid {'|'.join(sorted(_cv_grids()))}")
    esr = args.early_stopping or None
    tasks = [(args.n, d, s, args.reps, grid_name, args.max_rounds, esr, args.seed)
             for d in ds for s in settings]
    print(f"xgb cv grid '{grid_name}': {len(tasks)} (d, setting) cells x {args.reps} reps, n={args.n}")
    t0 = time.time()
    results = _map(_cv_task, tasks, min(args.workers, len(tasks)))
    if args.write:
        for (n, d, s, *_), (params, tuning) in zip(tasks, results):
            write_winner(args.out, n, d, s, "xgb", params, tuning)
    print(f"all done in {(time.time() - t0)/60:.1f} min")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("settings", nargs="*", default=None, help=f"subset of {SETTINGS}")
    ap.add_argument("--settings", dest="settings_csv", type=str, default=None,
                    help="comma-separated alternative to the positional settings")
    ap.add_argument("--learner", choices=["mlp", "xgb"], required=True)
    ap.add_argument("--protocol", choices=["holdout", "cv"], default="holdout",
                    help="holdout: the R protocol, either learner;  cv: K-fold CV, xgb only")
    ap.add_argument("--grid", default=None,
                    help="holdout: mlp pilot|full, xgb r|pilot;  cv: r|fast")
    ap.add_argument("--n", type=int, default=1000,
                    help="n the params are filed under (cv: also the simulated sample size)")
    ap.add_argument("--d", type=str, default="8", help="num_class; comma-separated for several")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--write", action="store_true", help="merge the winner into tuning/*.json")
    ap.add_argument("--out", type=str, default=str(TUNING_ROOT), help="tuning root (default: experiments/tuning)")
    ap.add_argument("--top", type=int, default=6, help="holdout: grid points to print")
    ap.add_argument("--max-rounds", type=int, default=1000, help="cv: nrounds ceiling")
    ap.add_argument("--early-stopping", type=int, default=50, help="cv: 0 disables")
    args = ap.parse_args()

    settings = list(args.settings or [])
    if args.settings_csv:
        settings += args.settings_csv.split(",")
    settings = settings or list(SETTINGS)
    ds = [int(v) for v in args.d.split(",")]

    if args.protocol == "cv":
        return main_cv(args, settings, ds)

    grids = MLP_GRIDS if args.learner == "mlp" else XGB_GRIDS
    grid_name = args.grid or ("full" if args.learner == "mlp" else "r")
    points = grid_points(grids[grid_name])
    n_flat = len(points) * (1 if args.learner == "mlp" else XGB_MAXROUNDS)

    print(f"{args.learner} grid '{grid_name}': {n_flat} points x {args.reps} reps, "
          f"n={args.n} d={args.d}, {len(settings)} settings")

    for d, setting in itertools.product(ds, settings):
        mean, ref, elapsed = tune(args.learner, setting, d, points,
                                  args.reps, args.workers, args.seed)
        best_i = int(np.argmin(mean))
        best_params, _ = describe(args.learner, points, mean, best_i)

        print(f"\n### d={d} {setting}  ({args.reps} reps, {elapsed/60:.1f} min)")
        print(f"    reference mlogloss: oracle={ref['oracle']:.5f}  uniform={ref['uniform']:.5f}")
        for i in np.argsort(mean)[: args.top]:
            _, line = describe(args.learner, points, mean, int(i))
            star = " <-- best" if i == best_i else ""
            print(f"    mlogloss={mean[i]:.5f}  {line}{star}")

        for warning in _boundary_warnings(args.learner, points, best_params):
            print(f"    WARNING: {warning}")

        if args.write:
            write_winner(args.out, args.n, d, setting, args.learner, best_params, {
                "protocol": "holdout", "grid": grid_name, "reps": args.reps, "seed": args.seed,
                "n_tr": N_TR, "n_te": N_TE, "d": d,
                "metric": "held-out mlogloss",
                "best_mlogloss": round(float(mean.min()), 6),
                "reference_mlogloss": {k: round(v, 6) for k, v in ref.items()},
            })


if __name__ == "__main__":
    main()
