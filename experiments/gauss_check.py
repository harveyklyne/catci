"""Is the truncated searches' small size excess finite-n non-Gaussianity, or a bug?

On each null replicate (semi-synthetic adult, oracle propensities), two observed
statistics are calibrated against the *same* ``B`` bootstrap draws from
``N(0, Sigma_hat)``:

* ``real``  -- the actual ``T`` from the data;
* ``gauss`` -- ``T`` replaced by one more independent draw from ``N(0, Sigma_hat)``.

``gauss`` is exchangeable with the bootstrap draws by construction, so the minP
p-value is exactly uniform and its rejection rate must be 0.05 up to Monte Carlo
error, for every method. If ``gauss`` is on 0.05 and ``real`` is above it, the excess
is the gap between ``T`` and its Gaussian approximation at this ``n``; if ``gauss``
is above 0.05 too, something in the search or calibration is wrong.

Usage:
    python gauss_check.py --x Education --y Income --reps 1000 --n-boot 3000
"""

from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import adult_semisynth as ss
from catci.bootstrap import bootstrap_T
from catci.calibrate import double_bootstrap_pvalue
from catci.gcm import form_t_sigma
from catci.search import MergeSearch, SplitSearch
from catci.statistic import ApproxChi
from run_semisynth import OUT, _setup, machine_info


def one_rep(task):
    args, rep_id, seed = task
    pop, delta, _, _, (xs, ys) = _setup(args)
    data_ss, boot_ss, extra_ss, cal_ss = seed.spawn(4)
    rep = ss.draw(pop, args.n, 0.0, np.random.default_rng(data_ss), z_names=tuple(args.z),
                  replace=True, delta=delta)
    ts = form_t_sigma(rep.x, rep.y, rep.f_true, rep.g_true, normalise=False)
    boot = bootstrap_T(ts.Sigma, args.n_boot, np.random.default_rng(boot_ss))
    t_gauss = bootstrap_T(ts.Sigma, 1, np.random.default_rng(extra_ss))[:, 0]
    stat = ApproxChi()
    rows = []
    for arm, t_obs in (("real", ts.T_vector), ("gauss", t_gauss)):
        T_all = np.column_stack([t_obs, boot])
        split = SplitSearch(max(args.depths)).paths(T_all, ts.Sigma, pop.dx, pop.dy, xs, ys, stat)
        merge = MergeSearch().paths(T_all, ts.Sigma, pop.dx, pop.dy, xs, ys, stat)
        cands = {f"split@{k}": split[: k + 1] for k in args.depths}
        cands.update({f"mergecoarse@{k}": merge[-(k + 1):] for k in args.depths})
        cands["merge"] = merge
        for name, paths in cands.items():
            rng = np.random.default_rng(cal_ss)  # same tie-breaks for both arms
            p = double_bootstrap_pvalue(paths[:, 0], paths[:, 1:], rng)
            rows.append(dict(arm=arm, method=name, rep=rep_id, p_value=float(p)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", required=True)
    ap.add_argument("--y", required=True)
    ap.add_argument("--w", nargs="+", default=list(ss.DEFAULT_W))
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--pool", type=float, default=5.0)
    ap.add_argument("--depths", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--n-boot", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    args.z, args.direction = list(args.w), "real"

    seeds = np.random.SeedSequence(args.seed).spawn(args.reps)
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        rows = [r for rr in ex.map(one_rep, [(args, i, s) for i, s in enumerate(seeds)], chunksize=4)
                for r in rr]
    df = pd.DataFrame(rows)
    stem = f"gausscheck_{args.x}_{args.y}_n{args.n}_B{args.n_boot}_r{args.reps}"
    OUT.mkdir(exist_ok=True)
    df.to_parquet(OUT / f"{stem}.parquet", index=False)
    (OUT / f"{stem}.provenance.json").write_text(json.dumps(
        dict(vars(args), seconds=time.time() - t0, machine=machine_info()), indent=2, default=str))
    df["reject"] = df.p_value < 0.05
    print(f"### {stem} ({time.time() - t0:.0f}s), SE ~ {np.sqrt(0.0475 / args.reps):.3f}")
    print(df.pivot_table(index="method", columns="arm", values="reject", sort=False).round(3).to_string())


if __name__ == "__main__":
    main()
