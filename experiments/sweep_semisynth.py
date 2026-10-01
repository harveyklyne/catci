"""Size and power against bootstrap size ``B`` and truncation depth ``k``, oracle only.

Semi-synthetic adult (``adult_semisynth.py``), exact propensities, so the only
moving parts are the search, its depth and the calibration. Everything is paired
and nested, which makes the whole ``(k, B)`` grid cost one search per replicate:

* **Depth.** A truncated divisive search is exactly a prefix of the full one
  (``SplitSearch(k).paths == SplitSearch().paths[:k+1]``, checked to 0.0), so
  ``split@k`` is the minP over the first ``k + 1`` rows of the full split path.
  As a comparator, ``mergecoarse@k`` is the minP over the *coarsest* ``k + 1``
  rows of the merge path -- the same multiplicity budget spent on the merge
  search's coarse end, which tests whether truncation's gain is multiplicity
  (``DIVISIVE.md``) or the divisive direction.
* **Bootstrap size.** One set of ``max(B)`` draws per replicate; each ``B`` uses
  the first ``B`` of them. The draws are i.i.d., so a prefix is a valid ``B``-draw
  bootstrap, and every ``B`` is paired with every other.

Output, in ``results_semisynth/``:

* ``<stem>.parquet`` -- one row per (replicate, lam, method, k, B): the p-value
  and the CPU seconds of its minP calibration (``cal_cpu``).
* ``<stem>_timing.parquet`` -- CPU and wall seconds of the searches. Every
  replicate times the full split, the merge and the depth-0 statistics at
  ``max(B)``. The first ``--time-reps`` replicates per ``lam`` additionally run
  ``SplitSearch(k)`` and ``MergeSearch`` *directly* at every ``B`` -- the prefix
  trick means truncated searches are otherwise never run on their own, so this
  is what measures how cost scales with ``k`` and ``B``.

Usage:
    python sweep_semisynth.py --x Education --y Income --direction planted \\
        --lams 0 0.4 0.8 1.2 --reps 500
"""

from __future__ import annotations

import argparse
import json
import subprocess
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
from catci.structure import Saturated
from run_semisynth import DEPTH0, OUT, Timer, _setup, machine_info


METHOD_IDS = {name: i for i, name in enumerate(["split", "mergecoarse", "merge", *DEPTH0])}


def one_rep(task) -> tuple[list[dict], list[dict]]:
    args, lam, rep_id, seed = task
    pop, delta, kinds, _, (xs, ys) = _setup(args)
    dx, dy = pop.dx, pop.dy
    data_rng, boot_rng, cal_seed = np.random.default_rng(seed).spawn(3)
    rep = ss.draw(pop, args.n, lam, data_rng, z_names=tuple(args.z), replace=True, delta=delta)
    if isinstance(xs, Saturated) or isinstance(ys, Saturated):
        raise ValueError("the depth sweep needs Ordinal/Tree on both sides (SplitSearch)")

    ts = form_t_sigma(rep.x, rep.y, rep.f_true, rep.g_true, normalise=False)
    B_max = max(args.n_boots)
    T_all = np.column_stack([ts.T_vector, bootstrap_T(ts.Sigma, B_max, boot_rng)])
    stat = ApproxChi()
    timing = []

    def timed(method, k, B, fn):
        with Timer() as t:
            out = fn()
        timing.append(dict(method=method, k=k, n_boot=B, lam=lam, rep=rep_id,
                           cpu=t.cpu, wall=t.wall))
        return out

    split = timed("split", None, B_max,
                  lambda: SplitSearch().paths(T_all, ts.Sigma, dx, dy, xs, ys, stat))
    merge = timed("merge", None, B_max,  # finest level first
                  lambda: MergeSearch().paths(T_all, ts.Sigma, dx, dy, xs, ys, stat))
    depth0 = {name: timed(name, 0, B_max, lambda fn=fn: np.array(
                  [[fn(T_all[:, b], ts.Sigma) for b in range(T_all.shape[1])]]))
              for name, fn in DEPTH0.items()}
    L = split.shape[0] - 1
    ks = sorted({k for k in args.depths if k < L} | {L})
    for row in timing:  # the full searches' depth, now that it is known
        if row["k"] is None:
            row["k"] = L if row["method"] == "split" else merge.shape[0] - 1

    if rep_id < args.time_reps:
        # Direct runs: the cost of split@k and merge as if each were the only method.
        for B in args.n_boots:
            T_B = T_all[:, : B + 1]
            for k in ks:
                timed("split_direct", k, B,
                      lambda: SplitSearch(k).paths(T_B, ts.Sigma, dx, dy, xs, ys, stat))
            timed("merge_direct", merge.shape[0] - 1, B,
                  lambda: MergeSearch().paths(T_B, ts.Sigma, dx, dy, xs, ys, stat))

    # The calibration's tie-break rng is reset to the same (replicate, method, k)
    # seed for every B, so B differences are not tie-break noise.
    base = int(cal_seed.integers(2**63))
    rows = []

    def calibrate(method, k, paths):
        for B in args.n_boots:
            rng = np.random.default_rng([base, METHOD_IDS[method], k])
            with Timer() as t:
                p = double_bootstrap_pvalue(paths[:, 0], paths[:, 1:B + 1], rng)
            rows.append(dict(method=method, k=k, n_boot=B, lam=lam, rep=rep_id, p_value=float(p),
                             argmax_level=int(np.argmax(paths[:, 0])), cal_cpu=t.cpu))

    for k in ks:
        calibrate("split", k, split[: k + 1])
        calibrate("mergecoarse", k, merge[-(k + 1):])
    calibrate("merge", merge.shape[0] - 1, merge)
    for name, path in depth0.items():
        calibrate(name, 0, path)
    return rows, timing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", required=True)
    ap.add_argument("--y", required=True)
    ap.add_argument("--w", nargs="+", default=list(ss.DEFAULT_W))
    ap.add_argument("--z", nargs="+", default=None)
    ap.add_argument("--direction", default="real", choices=["real", "planted"])
    ap.add_argument("--lams", type=float, nargs="+", default=[0.0])
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--pool", type=float, default=5.0)
    ap.add_argument("--depths", type=int, nargs="+", default=[1, 2, 3, 4, 6, 8])
    ap.add_argument("--n-boots", type=int, nargs="+", default=[100, 300, 1000, 3000, 10000])
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--time-reps", type=int, default=20,
                    help="replicates per lam that also time split@k / merge directly at every B")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    args.z = args.z or list(args.w)
    args.n_boots = sorted(args.n_boots)

    pop, delta, kinds, _, (xs, ys) = _setup(args)
    print(f"{args.x} ({kinds[0]}, dx={pop.dx}, {type(xs).__name__}) x "
          f"{args.y} ({kinds[1]}, dy={pop.dy}, {type(ys).__name__}) "
          f"| Z={args.z}, direction={args.direction}, "
          f"max_lambda={ss.max_lambda(pop, delta):.3f}, ncp/n={pop.ncp_per_n(delta):.4g}")

    seeds = np.random.SeedSequence(args.seed).spawn(len(args.lams) * args.reps)
    tasks = [(args, lam, r, seeds[i * args.reps + r])
             for i, lam in enumerate(args.lams) for r in range(args.reps)]
    t0 = time.time()
    rows, timing = [], []
    with ProcessPoolExecutor(args.workers) as ex:
        for r, t in ex.map(one_rep, tasks, chunksize=1):
            rows += r
            timing += t
    df, tdf = pd.DataFrame(rows), pd.DataFrame(timing)

    OUT.mkdir(exist_ok=True)
    stem = (f"sweep_{args.x}_{args.y}_{args.direction}_Z{'-'.join(args.z)}_n{args.n}"
            f"_r{args.reps}{('_' + args.tag) if args.tag else ''}")
    df.to_parquet(OUT / f"{stem}.parquet", index=False)
    tdf.to_parquet(OUT / f"{stem}_timing.parquet", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (OUT / f"{stem}.provenance.json").write_text(json.dumps(
        dict(vars(args), git_commit=commit, seconds=time.time() - t0, machine=machine_info()),
        indent=2, default=str))

    print(f"\n### {stem}  ({time.time() - t0:.0f}s)  reps={args.reps}")
    df["reject"] = df.p_value < 0.05
    for lam, block in df.groupby("lam"):
        label = "size" if lam == 0 else "power"
        print(f"\n## lam = {lam} ({label} at 0.05)  rows: method@k, columns: B")
        block = block.assign(m=block.method + "@" + block.k.astype(str))
        print(block.pivot_table(index="m", columns="n_boot", values="reject", sort=False)
                   .round(3).to_string())

    print("\n## CPU seconds per search (direct runs), rows method@k, columns B")
    direct = tdf[tdf.method.str.endswith("_direct")].assign(
        m=lambda d: d.method.str.removesuffix("_direct") + "@" + d.k.astype(str))
    print(direct.pivot_table(index="m", columns="n_boot", values="cpu", sort=False)
                .round(3).to_string())
    print("\ncalibration only:")
    cal = df.assign(m=df.method + "@" + df.k.astype(str))
    print(cal.pivot_table(index="m", columns="n_boot", values="cal_cpu", sort=False)
             .round(4).to_string())


if __name__ == "__main__":
    main()
