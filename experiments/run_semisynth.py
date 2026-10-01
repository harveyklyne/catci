"""Size and power on the semi-synthetic adult data (``adult_semisynth.py``).

Real ``(X, Z)`` rows from UCI adult, ``Y`` simulated from the mixture kernel, so
the truth of ``X indep Y | Z`` is known exactly and ``lam`` dials the effect.

Every replicate draws one dataset and, per learner, one set of bootstrap draws;
every method is calibrated on exactly those, so differences between methods are
**paired**. The learners share their bootstrap seed too (common random numbers),
so ``oracle`` vs a fitted learner isolates the cost of estimating the nuisances
from the effect of the search.

Learners:

* ``oracle`` -- the exact population propensities (:func:`adult_semisynth.true_propensities`);
* ``mlp`` / ``xgb`` -- fitted on the one-hot ``Z`` on the full sample, as
  :func:`catci.catci_test` does.

Methods (all on the same ``T, Sigma``, ``normalise=False``):

* ``merge``            -- :class:`~catci.search.MergeSearch` (Algorithm 1)
* ``split``, ``split@k`` -- :class:`~catci.search.SplitSearch`, full / truncated at
  ``k`` levels. Needs ``Ordinal`` or ``Tree`` on both sides; skipped for ``Saturated``.
* ``max``, ``euclid``, ``mGCM`` -- depth-0 comparators, minP-calibrated (``L = 1``).
  ``mGCM`` studentises internally (Shah & Peters' normalised max).
* ``chi_sq``           -- pseudo-inverse chi-square on ``(dx-1)(dy-1)`` df.
* ``at_typed``         -- Ankan & Textor with each variable typed as the data
  types it (Q1 for ordinal x ordinal), on the same propensities.
* ``at_cat``           -- Ankan & Textor with both typed categorical (Q3).
* ``strat_chi2``       -- Pearson chi-square stratified over the ``Z`` cells.

Usage:
    python run_semisynth.py --x Education --y Income --lams 0 --reps 500
    python run_semisynth.py --x Education --y Income --direction planted \\
        --lams 0 0.5 1 1.5 2 --learners oracle mlp --truncations 2 4
    python run_semisynth.py --x Occupation --y Income --pool 5 --lams 0
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

import adult
import adult_semisynth as ss
import ankan_textor as at
import methods
import taxonomies
from catci.bootstrap import bootstrap_T
from catci.calibrate import double_bootstrap_pvalue
from catci.gcm import form_t_sigma
from catci.learners import fit_propensities, mlp_learner
from catci.search import MergeSearch, SplitSearch
from catci.statistic import ApproxChi, euclid, max_abs, mgcm
from catci.structure import Ordinal, Saturated

OUT = Path(__file__).resolve().parent / "results_semisynth"


class Timer:
    """CPU and wall seconds for a block: ``with Timer() as t: ...; t.cpu, t.wall``.

    CPU is ``time.process_time()`` -- this process only, all its threads. Workers
    run single-threaded BLAS (set ``OMP_NUM_THREADS=1`` etc.), so it is the
    method's own cost, unaffected by other load on the machine; wall time is
    kept beside it to show contention.
    """

    def __enter__(self):
        self._c, self._w = time.process_time(), time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.cpu = time.process_time() - self._c
        self.wall = time.perf_counter() - self._w
        return False


def machine_info() -> dict:
    """What a CPU time was measured on, for the provenance sidecar."""
    import os
    import platform

    def sysctl(key):
        try:
            return subprocess.run(["sysctl", "-n", key], capture_output=True, text=True).stdout.strip()
        except OSError:
            return ""

    return dict(
        platform=platform.platform(), processor=platform.processor(),
        cpu_brand=sysctl("machdep.cpu.brand_string"), n_cpu=os.cpu_count(),
        load_avg_at_start=os.getloadavg(), python=platform.python_version(),
        numpy=np.__version__,
        blas_threads={k: os.environ.get(k) for k in
                      ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                       "VECLIB_MAXIMUM_THREADS")},
    )
DEPTH0 = {"max": max_abs, "euclid": euclid, "mGCM": mgcm}


def structure_for(kind: str, d: int):
    """Ordinal variables get ``Ordinal``; unordered ones ``Saturated`` (merge only).

    :func:`_setup` overrides this with a :mod:`taxonomies` tree for an unordered
    variable that has one, which is what lets the divisive search run on it.
    """
    if kind == "ordinal" or d <= 2:
        return Ordinal()
    return Saturated()


def searches(xs, ys, truncations, full_split: bool) -> dict:
    out = {"merge": MergeSearch()}
    if not (isinstance(xs, Saturated) or isinstance(ys, Saturated)):
        if full_split:
            out["split"] = SplitSearch()
        out.update({f"split@{k}": SplitSearch(max_levels=k) for k in truncations})
    return out


def stratified_chi2(x, y, zcell) -> float:
    """Pearson chi-square summed over ``Z`` cells, df summed over the observed sub-tables."""
    stat, df = 0.0, 0
    for c in np.unique(zcell):
        m = zcell == c
        xs, xi = np.unique(x[m], return_inverse=True)
        ys, yi = np.unique(y[m], return_inverse=True)
        if len(xs) < 2 or len(ys) < 2:
            continue
        tab = np.zeros((len(xs), len(ys)))
        np.add.at(tab, (xi, yi), 1.0)
        E = tab.sum(1, keepdims=True) * tab.sum(0, keepdims=True) / tab.sum()
        stat += float(((tab - E) ** 2 / E).sum())
        df += (len(xs) - 1) * (len(ys) - 1)
    return float(chi2.sf(stat, df)) if df > 0 else 1.0


_AT_NOTE: list = []


def ankan_textor_pvalue(rx, ry, kx, ky) -> float:
    """The paper's ``solve`` form; ``pinv`` (df = rank) only when that is singular.

    A rare level whose propensity is exactly 0 in every sampled stratum gives an
    all-zero residual column, so ``Sigma_d`` is singular and the paper's statistic
    is undefined. ``pinv`` is the port's documented fallback; the row is flagged in
    ``note`` so those replicates can be excluded or reported.
    """
    try:
        res = at.test_from_residuals(rx, ry, kx, ky)
        if res.well_conditioned:
            _AT_NOTE.append("")
            return res.p_value
    except np.linalg.LinAlgError:
        pass
    _AT_NOTE.append("pinv")
    return at.test_from_residuals(rx, ry, kx, ky, method="pinv").p_value


def fitted_propensities(learner: str, rep: ss.Replicate, z_levels, dx, dy):
    if learner == "oracle":
        return rep.f_true, rep.g_true
    design = at.design_matrix(rep.z, z_levels, drop_first=False)
    if learner == "mlp":
        lr = mlp_learner()
    elif learner == "xgb":
        from catci.learners import xgboost_learner
        lr = xgboost_learner(dict(eta=0.1, max_depth=3, gamma=0.0, nrounds=100))
    else:
        raise ValueError(f"unknown learner {learner!r}")
    return fit_propensities(design, rep.x, dx, lr), fit_propensities(design, rep.y, dy, lr)


def one_rep(task) -> list[dict]:
    args, lam, rep_id, seed = task
    pop, delta, kinds, z_levels, (xs, ys) = _setup(args)
    dx, dy = pop.dx, pop.dy
    data_rng, boot_seed, method_seed = np.random.default_rng(seed).spawn(3)
    rep = ss.draw(pop, args.n, lam, data_rng, z_names=tuple(args.z), replace=True, delta=delta)
    zcell = np.unique(rep.z, axis=0, return_inverse=True)[1].ravel()
    stat = ApproxChi()

    rows = []
    for learner in args.learners:
        with Timer() as tf:
            f, g = fitted_propensities(learner, rep, z_levels, dx, dy)
        ts = form_t_sigma(rep.x, rep.y, f, g, normalise=False)
        # Common random numbers: every learner sees the same N(0, I) draws.
        boot = bootstrap_T(ts.Sigma, args.n_boot, np.random.default_rng(boot_seed))
        T_all = np.column_stack([ts.T_vector, boot])
        rng = np.random.default_rng(method_seed)  # tie-breaks in the calibration

        def record(name, p, search_t, cal_t=None, paths=None):
            rows.append(dict(
                method=name, learner=learner, lam=lam, rep=rep_id, p_value=float(p),
                # search_* = statistic path for the observed T and all B draws;
                # cal_* = the minP calibration on those paths; fit_* = propensities.
                search_cpu=search_t.cpu, search_wall=search_t.wall,
                cal_cpu=None if cal_t is None else cal_t.cpu,
                fit_cpu=tf.cpu, fit_wall=tf.wall, n_boot=args.n_boot,
                note="", levels=None if paths is None else paths.shape[0],
                argmax_level=None if paths is None else int(np.argmax(paths[:, 0])),
            ))

        for name, search in searches(xs, ys, args.truncations, args.full_split).items():
            with Timer() as t_s:
                paths = search.paths(T_all, ts.Sigma, dx, dy, xs, ys, stat)
            with Timer() as t_c:
                p = double_bootstrap_pvalue(paths[:, 0], paths[:, 1:], rng)
            record(name, p, t_s, t_c, paths)
        for name, fn in DEPTH0.items():
            with Timer() as t_s:
                path = np.array([[fn(T_all[:, b], ts.Sigma) for b in range(T_all.shape[1])]])
            with Timer() as t_c:
                p = double_bootstrap_pvalue(path[:, 0], path[:, 1:], rng)
            record(name, p, t_s, t_c)

        fitted = methods.Fitted(rep.x, rep.y, rep.z, f, g, dx, dy, ts.T_vector, ts.Sigma)
        with Timer() as t_s:
            p = methods.chi_sq(fitted)
        record("chi_sq", p, t_s)
        for name, (kx, ky) in {"at_typed": kinds, "at_cat": ("categorical", "categorical")}.items():
            with Timer() as t_s:
                rx = at.residual_matrix(rep.x, f, kx)
                ry = at.residual_matrix(rep.y, g, ky)
                p = ankan_textor_pvalue(rx, ry, kx, ky)
            record(name, p, t_s)
            rows[-1]["note"] = _AT_NOTE.pop()

    # Learner-free competitor: once per replicate.
    with Timer() as t_s:
        p = stratified_chi2(rep.x, rep.y, zcell)
    rows.append(dict(method="strat_chi2", learner="none", lam=lam, rep=rep_id, p_value=p,
                     search_cpu=t_s.cpu, search_wall=t_s.wall, cal_cpu=None, fit_cpu=0.0,
                     fit_wall=0.0, n_boot=None, note="", levels=None, argmax_level=None))
    return rows


_CACHE: dict = {}


def _setup(args):
    """Population, direction, kinds, Z level counts and structures -- once per worker.

    An unordered variable with a :mod:`taxonomies` entry is recoded to the tree's
    leaf order and searched under that tree (unless ``--no-taxonomy``); other
    unordered variables get ``Saturated``, ordinal ones ``Ordinal``.
    """
    taxonomy = getattr(args, "taxonomy", True)
    key = (args.x, args.y, tuple(args.w), args.pool, args.n, args.direction, taxonomy)
    if key not in _CACHE:
        data = adult.load()
        if args.pool:
            for name in (args.x, args.y):
                if data.kind(name) == "categorical":
                    data = ss.pool_rare(data, name, args.pool, args.n)
        trees = {}
        for name in (args.x, args.y):
            if taxonomy and data.kind(name) == "categorical" and name in taxonomies.TAXONOMIES:
                data, trees[name] = taxonomies.apply_taxonomy(data, name)
        pop = ss.build_population(data, args.x, args.y, args.w)
        delta = pop.delta_real() if args.direction == "real" else ss.planted_direction(pop)
        kinds = (data.kind(args.x), data.kind(args.y))
        z_levels = [data.n_levels(name) for name in args.z]
        structures = tuple(trees.get(name) or structure_for(kind, d) for name, kind, d in
                           ((args.x, kinds[0], pop.dx), (args.y, kinds[1], pop.dy)))
        _CACHE[key] = (pop, delta, kinds, z_levels, structures)
    return _CACHE[key]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", required=True)
    ap.add_argument("--y", required=True)
    ap.add_argument("--w", nargs="+", default=list(ss.DEFAULT_W),
                    help="generating stratum (Y depends on X only within W cells)")
    ap.add_argument("--z", nargs="+", default=None,
                    help="conditioning set tested; must contain --w (default: --w)")
    ap.add_argument("--direction", default="real", choices=["real", "planted"])
    ap.add_argument("--lams", type=float, nargs="+", default=[0.0])
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--pool", type=float, default=5.0,
                    help="pool categorical levels with expected count < this (0 = off)")
    ap.add_argument("--learners", nargs="+", default=["oracle"])
    ap.add_argument("--no-taxonomy", dest="taxonomy", action="store_false",
                    help="search unordered variables as Saturated even if a tree exists")
    ap.add_argument("--truncations", type=int, nargs="+", default=[2, 4])
    ap.add_argument("--no-full-split", dest="full_split", action="store_false")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=10000,
                    help="bootstrap draws; size checks are exact at any B, power needs 10k")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    args.z = args.z or list(args.w)

    pop, delta, kinds, _, (xs, ys) = _setup(args)
    print(f"{args.x} ({kinds[0]}, dx={pop.dx}, {type(xs).__name__}) x "
          f"{args.y} ({kinds[1]}, dy={pop.dy}, {type(ys).__name__}) "
          f"| Z={args.z}, W={args.w}, direction={args.direction}, "
          f"max_lambda={ss.max_lambda(pop, delta):.3f}, ncp/n={pop.ncp_per_n(delta):.4g}")

    seeds = np.random.SeedSequence(args.seed).spawn(len(args.lams) * args.reps)
    tasks = [(args, lam, r, seeds[i * args.reps + r])
             for i, lam in enumerate(args.lams) for r in range(args.reps)]
    t0 = time.time()
    with ProcessPoolExecutor(args.workers) as ex:
        rows = [row for rr in ex.map(one_rep, tasks, chunksize=1) for row in rr]
    df = pd.DataFrame(rows)

    OUT.mkdir(exist_ok=True)
    zs = "-".join(args.z)
    stem = (f"{args.x}_{args.y}_{args.direction}_Z{zs}_n{args.n}_B{args.n_boot}"
            f"_r{args.reps}{('_' + args.tag) if args.tag else ''}")
    df.to_parquet(OUT / f"{stem}.parquet", index=False)
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (OUT / f"{stem}.provenance.json").write_text(json.dumps(
        dict(vars(args), git_commit=commit, seconds=time.time() - t0, machine=machine_info()),
        indent=2, default=str))

    print(f"\n### {stem}  ({time.time() - t0:.0f}s)  size SE~{np.sqrt(0.05 * 0.95 / args.reps):.3f}")
    df["reject"] = df.p_value < 0.05
    order = list(dict.fromkeys(df.method))
    print(df.pivot_table(index=["learner", "lam"], columns="method", values="reject", sort=False)
            [order].round(3).to_string())
    print("\nmean CPU seconds per test (search, calibration, propensity fit):")
    print(df.groupby(["learner", "method"], sort=False)[["search_cpu", "cal_cpu", "fit_cpu"]]
            .mean().round(4).to_string())


if __name__ == "__main__":
    main()
