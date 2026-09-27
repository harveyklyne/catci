"""Score tuned XGBoost hyperparameters against the *true* propensities.

CV log-loss ranks hyperparameters by ``entropy + KL(f || f_hat)``; with simulated
data the entropy is known, so we can report the KL itself -- the quantity the
test actually consumes. Each replicate fits on ``n`` fresh draws (the full
sample, as ``run.py`` does) and reports mean ``KL(f || f_hat)``:

* ``kl_in``  -- at the training points. This is what the test sees: there is no
  cross-fitting, the propensities are predicted back on the fitted ``Z``.
* ``kl_out`` -- at fresh points.

Also reports the ``Z``-blind baseline (the empirical class frequencies), so a KL
can be read as "how much of the available signal did the learner pick up".

Usage:
    python tune_check.py --n 1000 --d 8 --settings sin,lin \\
        --tunings experiments/tuning /path/to/other/tuning [--reps 20]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

import dgp
from catci.learners import xgboost_learner
from config import tuning_path


def _kl(f, fhat):
    fhat = np.clip(fhat, 1e-12, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.nanmean(np.sum(np.where(f > 0, f * np.log(f / fhat), 0.0), axis=1)))


def check(n, d, setting, params_by_label, reps, seed=0, n_out=5000):
    rng = np.random.default_rng([seed, n, d, sum(map(ord, setting))])
    rows = {label: dict(kl_in=[], kl_out=[], seconds=[]) for label in list(params_by_label) + ["freq"]}
    for _ in range(reps):
        tr = dgp.simulate_marginal(n, d, setting, rng)
        te = dgp.simulate_marginal(n_out, d, setting, rng)
        freq = (np.bincount(tr["x"], minlength=d + 1)[1:] + 0.5) / (n + 0.5 * d)
        rows["freq"]["kl_in"].append(_kl(tr["f"], np.broadcast_to(freq, tr["f"].shape)))
        rows["freq"]["kl_out"].append(_kl(te["f"], np.broadcast_to(freq, te["f"].shape)))
        rows["freq"]["seconds"].append(0.0)
        for label, params in params_by_label.items():
            t = time.time()
            predict = xgboost_learner(params)(tr["z"], tr["x"], d)
            rows[label]["kl_in"].append(_kl(tr["f"], predict(tr["z"])))
            rows[label]["seconds"].append(time.time() - t)
            rows[label]["kl_out"].append(_kl(te["f"], predict(te["z"])))
    out = {}
    for label, r in rows.items():
        out[label] = {k: (float(np.mean(v)), float(np.std(v) / np.sqrt(len(v)))) for k, v in r.items()}
        out[label]["params"] = params_by_label.get(label)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--d", type=int, required=True)
    ap.add_argument("--settings", type=str, default="sin,sig,lin,vee,hat")
    ap.add_argument("--tunings", nargs="+", required=True, help="tuning roots to compare")
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--json", type=str, default=None, help="write the results here too")
    args = ap.parse_args()

    allres = {}
    for s in args.settings.split(","):
        params = {}
        for root in args.tunings:
            p = tuning_path(args.n, args.d, s, root=Path(root))
            params[Path(root).name] = json.loads(p.read_text())["xgb"]
        res = check(args.n, args.d, s, params, args.reps, seed=args.seed)
        allres[s] = res
        print(f"\n## n={args.n} d={args.d} {s}   (mean over {args.reps} reps, +- SE)")
        print(f"   {'tuning':24s} {'KL in-sample':>18s} {'KL fresh':>18s} {'fit s':>7s}  params")
        for label, r in res.items():
            print(f"   {label:24s} {r['kl_in'][0]:9.4f} +- {r['kl_in'][1]:.4f} "
                  f"{r['kl_out'][0]:9.4f} +- {r['kl_out'][1]:.4f} {r['seconds'][0]:7.2f}  {r['params'] or ''}")
    if args.json:
        Path(args.json).write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()
