# MLP propensity learner — worktree notes

Scratch notes for the `MLP classifier instead of xgboost` TODO item. Untracked
(`*.md` is gitignored); the code itself is committed on the branch.

- **Worktree**: `.claude/worktrees/mlp-learner`
- **Branch**: `worktree-mlp-learner`, branched from `196cc3a` (master)
- **Commit**: `4b395b0` *Add an MLP propensity learner and tune it against XGBoost*
- **Status**: code complete, committed, 68 tests passing. The size/power sweep
  was **deliberately stopped at 2/10 runs** — it is not still running and nothing
  is half-written. See [Where this was left](#where-this-was-left).

---

## TL;DR

The MLP is a **better propensity estimator** than the tuned XGBoost on 4 of 5
marginal settings (27–72% lower `E_f`), and never worse. But it probably **does
not matter for the test**, because both learners already satisfy the paper's
Assumption 1 with ~1000x margin.

The binding constraint on size at `n = 1000` is the **test's own calibration**,
and the oracle control localises it precisely: with *exact* propensities, the
depth-0 comparators (`max`, `euclid`) calibrate correctly while the adaptive
search methods (`tree`, `ordinal`) reject at 1.5–1.9x nominal. That is a bug or
a small-`n_boot` problem in the adaptive path, and no learner can touch it.

**Power was never measured** — those runs were stopped before they started — so
the "does not matter" claim rests on the `E_f`/Assumption-1 margin argument and
is a *prediction*, not a result. Treat it as such.

Provisional recommendation, to be confirmed by the power runs: keep xgboost as
the default, cite the MLP as evidence that the *learner* is not what limits the
method, and spend the effort on the adaptive calibration gap instead.

---

## Environment gotcha (read this first)

The `catci` conda env has an **editable install pointing at the main checkout**,
not at this worktree. Plain `python` in this worktree silently imports
`/Users/harveyklyne/Documents/Github/catci/src/catci` — the *old* code, with no
`mlp_learner`. Every command below therefore sets `PYTHONPATH` explicitly:

```sh
WT=/Users/harveyklyne/Documents/Github/catci/.claude/worktrees/mlp-learner
export PYTHONPATH="$WT/src:$WT/experiments"
export OMP_NUM_THREADS=1
PY=/Users/harveyklyne/miniforge3/envs/catci/bin/python
```

`conda activate catci` is equivalent for the interpreter, but does **not** fix
the `src` shadowing — the `PYTHONPATH` entry is what does. Verify with:

```sh
$PY -c "import catci; print(catci.__file__)"   # must print the worktree path
```

---

## What changed

| File | Change |
|---|---|
| `src/catci/learners.py` | new `mlp_learner(params)` |
| `experiments/tune.py` | **new** — tuner for *both* learners |
| `experiments/bench_learners.py` | **new** — propensity quality on `E_f` |
| `experiments/report_learners.py` | **new** — prints the comparison tables |
| `experiments/bench_sweep.sh` | **new** — drives the whole size/power sweep |
| `tests/test_learners.py` | **new** — 6 interface/contract tests |
| `experiments/config.py` | `Config.learner` field + `learner_params()` |
| `experiments/run.py`, `run_size.py` | `--learner xgb\|mlp\|oracle` |
| `experiments/tuning/n1000_numclass8/*.json` | added an `mlp` block per setting |
| `pyproject.toml` | declared `threadpoolctl` (imported directly now) |

Three implementation details worth remembering:

1. **Column scattering.** `MLPClassifier.predict_proba` returns one column per
   class it *saw*. If a level is absent from the training labels, the raw output
   is `(n, d-1)` and every column past the gap is shifted left — silently wrong
   propensities. `mlp_learner` scatters through `clf.classes_` into a full
   `num_class`-wide matrix. Pinned by `test_mlp_pads_absent_classes_to_num_class`.
2. **`fork` → `spawn`.** `run.py` used `mp.get_context("fork")`. sklearn starts
   BLAS threads, and a forked child that inherits them aborts in the macOS
   Objective-C runtime (`+[NSString initialize] may have been in progress...`).
   This bit on the *second* pool in a process, so it looked intermittent. All
   pools now use `spawn`. **This was a latent bug in `run.py` before this branch**
   — it just needed a threaded learner to surface.
3. **Input scaling.** The net standardises `Z` internally. Without it the MLP
   underfits toward the marginal, which is indistinguishable from "the MLP is
   bad" unless you look. `test_mlp_beats_the_constant_predictor_when_z_is_informative`
   is the floor that would catch a regression here.

---

## Why `E_f` is the metric

The paper (`experiments/Conditional_independence_testing_with_categorical_data (3).pdf`,
Assumption 1, p. 9) controls type I error through

```
E_f := max_j E[ (f_j(Z) - fhat_j(Z))^2 | D ]
E_f = o_P(1),  E_g = o_P(1),  E_f * E_g = o_P(n^-1)
```

So `E_f` — not accuracy, not log-loss — is the quantity a propensity learner
should be judged on *for this method*. The DGP returns the true `f`/`g`, so it
is computed exactly rather than estimated.

The paper's only stated justification for xgboost is empirical: *"We have found
gradient boosting multinomial regressions (xgboost package) to do well in
practice"* (p. 5, repeated p. 8). That is the claim this branch tests.

Two useful facts that fell out:

- **Strength 0 suffices.** The interaction matrices have zero row and column
  sums, so summing the joint over `Y` returns `f` exactly at any strength. The
  true `X|Z` propensity is identical under null and alternative, so `E_f`
  measured at strength 0 transfers to the power runs.
- **In-sample ≈ out-of-sample.** Cross-fitting was dropped, so propensities are
  fitted and consumed on the same rows; `bench_learners.py` reports both. The
  ratio is ≈1.00 for both learners after tuning, so the no-cross-fitting design
  is currently costing nothing. That is a property of *these tuned
  hyperparameters*, not a general guarantee — an untuned MLP overfits hard.

---

## Fairness: both learners are now tuned

The frozen `tuning/*.json` came from the R `tune_xgb` (recovered from the
`r-frozen-oracle` tag, `data-raw/tuning/tuning_functions.R`): simulate a fresh
train/test pair per rep, fit every grid point, score **held-out mlogloss**,
average over reps, argmin. Constants `n_tr=800`, `n_te=5000`, `strength=0.5`.

`experiments/tune.py` ports that protocol for **both** learners. Comparing a
tuned xgboost against a default MLP would have been meaningless — the untuned
MLP was ~10x *worse* on `E_f` before tuning.

The xgboost side reuses the R trick of reading the whole `nrounds` axis off one
`xgb.train` eval history, so a 1000-long grid is affordable. Validated by
recovering the frozen `depth=1, gamma=2` for `lin`.

This also fills part of the README's **"No tuner"** known gap: a new `(n, d)` can
now be tuned off-cluster, for either learner.

```sh
$PY $WT/experiments/tune.py --learner mlp --reps 20 --write
$PY $WT/experiments/tune.py --learner xgb --reps 20 --d 30 --write
```

`--write` merges the winner into the per-setting JSON (both learners coexist in
one file). The tuner warns when the argmin sits on a grid boundary — which
happened twice during this work and would otherwise have been missed.

---

## Results

### 1. Tuning outcome

`hidden_layer_sizes=(8,8), alpha=3.0, activation=tanh, max_iter=400` won on
**all five** settings. The first pilot railed at the top of its alpha range
everywhere — these propensities carry little signal about `Z`, and the tuner
wants heavy shrinkage.

Held-out mlogloss (20 reps), with the fraction of the oracle→uniform gap each
learner recovers:

| setting | oracle | xgb | mlp | uniform | xgb recovers | mlp recovers |
|---|---|---|---|---|---|---|
| lin | 2.0508 | 2.0691 | **2.0652** | 2.0794 | 36% | **50%** |
| vee | 2.0390 | 2.0593 | **2.0515** | 2.0794 | 50% | **69%** |
| hat | 2.0089 | 2.0349 | **2.0215** | 2.0794 | 63% | **82%** |
| sin | 2.0631 | 2.0832 | **2.0796** | 2.0794 | **−23%** | 0% |
| sig | 2.0637 | 2.0791 | **2.0753** | 2.0794 | 2% | **26%** |

**MLP wins all five.** Note `sin`: tuned xgboost scores *worse than predicting
uniform* — the stumps are actively harmful there, and the MLP merely breaks even.

### 2. Propensity quality — `E_f`, 100 reps, n=1000, d=8

| setting | `E_f` xgb | `E_f` mlp | ratio | `n·E_f·E_g` xgb | `n·E_f·E_g` mlp |
|---|---|---|---|---|---|
| hat | 0.00133 | **0.00077** | **1.72x** | 0.0017 | 0.0006 |
| vee | 0.00114 | **0.00076** | **1.51x** | 0.0013 | 0.0006 |
| lin | 0.00106 | **0.00077** | **1.37x** | 0.0011 | 0.0006 |
| sig | 0.00092 | **0.00073** | **1.27x** | 0.0008 | 0.0005 |
| sin | **0.00112** | 0.00113 | 0.99x | 0.0013 | 0.0013 |

Fit time per call is comparable (xgb 0.14–0.95 s, mlp 0.57–0.82 s); the MLP is
*faster* on `sin`, where xgboost's tuned `nrounds=861` is expensive.

**The caveat that matters:** `n·E_f·E_g ≈ 0.001` for *both* learners. Assumption
1 asks for this to go to zero; it is already ~1000x inside the requirement. A
1.7x improvement on a remainder that is this negligible has little room to
change the test's behaviour. That is the hypothesis the size/power sweep tests.

### 3. Size — the oracle control (1000 reps, both settings) ← the important one

Run with **exact** propensities, i.e. `E_f = 0`. Binomial SE ≈ 0.007, so
|rate − 0.05| > 0.014 is notable.

| method | `lin_lin` | `sin_sin` | |
|---|---|---|---|
| max | 0.046 | 0.040 | calibrated |
| euclid | 0.042 | 0.052 | calibrated |
| ankan | 0.062 | 0.046 | ~calibrated |
| **tree** | **0.077** | **0.070** | **inflated ~1.5x** |
| **ordinal** | **0.096** | **0.091** | **inflated ~1.9x** |
| mGCM | 0.169 | 0.189 | inflated (expected) |
| chi_sq | 0.161 | 0.141 | inflated (expected) |

**This is the most consequential finding on the branch, and it is not about the
MLP.** No propensity learner can improve on `E_f = 0`, so this inflation is the
*test's* calibration at `n = 1000, d = 8`, not the regression's fault.

`mGCM` and `chi_sq` are expected to fail — the paper predicts chi_sq does poorly
when `d_X d_Y = 64` is not small relative to `n`, and mGCM inherits the same
pinv problem.

**The sharp part of the diagnosis** — the two settings agree closely, so this is
not Monte Carlo noise:

- `max` and `euclid` are **correctly calibrated**. They are scalar criteria
  evaluated at **search depth 0**, and they go through the *same*
  `double_bootstrap_pvalue` and the *same* bootstrap draws.
- `tree` and `ordinal` are **inflated**, and they differ from `max`/`euclid` in
  exactly one respect: they return a **vector** of criteria over the merge path,
  and the p-value maximises over search depth.

So the mis-calibration is specific to the **adaptive search path** — not to the
GCM statistic, not to the propensities, and not to the double bootstrap as such,
since the depth-0 path exercises all three and calibrates fine. That localises
it to the vector branch of `calibrate.double_bootstrap_pvalue` /
`adaptive_pvalue` (the paper's Algorithm 2 / Theorem 2), or to `n_boot = 100`
being too small to estimate the max-over-depth null distribution.

Worth stressing: the adaptive methods are the paper's *contribution*. An
inflated `ordinal` at 0.091–0.096 with perfect propensities undercuts the size
claim for exactly the methods the paper is about, so this deserves priority over
anything on this branch.

Implication for this branch: any size difference between xgb and mlp must be
read against this floor, and the ceiling on how much a better learner can help
is small.

---

## Where this was left

The sweep was **stopped on purpose at 2 of 10 runs**, mid-way through
`size lin_lin learner=xgb`. Nothing is still running; no parquet is
half-written (each is written once, atomically, at the end of its run).

**What exists in `experiments/results/` and is trustworthy:**

| file | what it is |
|---|---|
| `bench_learners.parquet` | complete — 100 reps x 5 settings, the `E_f` table |
| `size_lin_lin__oracle.parquet` | complete — 1000 reps |
| `size_sin_sin__oracle.parquet` | complete — 1000 reps |
| `sweep.log` | console output of the 2 finished runs |
| `SMOKE-TEST-NOT-A-RESULT_size_lin_lin__mlp.40reps.parquet` | **not a result** |

That last file was a 40-rep smoke test I ran to time the pipeline. It was the
only `mlp` size file on disk and would have been easy to mistake for a real
1000-rep result, so it is renamed out of the way — the prefix also stops
`report_learners.py` from picking it up. Delete it freely.

**Never ran at all:** the four power runs. There is **no power data on this
branch**, so no claim about the MLP's effect on power is supported yet.

**To finish it**, the sweep is idempotent — just rerun the whole thing (the two
oracle runs will be recomputed, ~20 min of the ~80):

```sh
sh $WT/experiments/bench_sweep.sh > $WT/experiments/results/sweep.log 2>&1 &
$PY $WT/experiments/report_learners.py     # tables from whatever parquets exist
```

Still outstanding: **size** `{xgb, mlp}` x `{lin_lin, sin_sin}` (4 runs), and
**power** 200 reps x strengths `0.2,0.6,1.0,1.4,1.8`, `{xgb, mlp}` x
`{lin_lin_step, sin_sin_binary_tree}` (4 runs).

Power is designed as a **paired** comparison: `run.py` seeds from
`SeedSequence(seed)` and the learner does not enter the seed stream, so both
runs of a config see the same simulated datasets replicate-for-replicate, and
the same bootstrap normals. `report_learners.py` reports the paired difference
and its SE, which is far tighter than comparing two independent rejection rates.
This matters because the expected effect is small — an unpaired design at 200
reps would likely not resolve it.

---

## Next steps, in the order I'd do them

1. **Chase the adaptive size inflation.** Promoted to first on the strength of
   the oracle control above: `ordinal` at 0.091–0.096 with perfect propensities,
   while `max`/`euclid` sit at nominal, is a clean localisation and it affects
   the paper's headline methods. Candidates, cheapest first: (a) raise `n_boot`
   from 100 and see if the inflation shrinks — if it does, it is Monte Carlo
   error in the max-over-depth null, not a bug; (b) audit the vector branch of
   `double_bootstrap_pvalue` against Algorithm 2; (c) genuine `n = 1000, d = 8`
   asymptotics gap. Deserves its own worktree.
2. **Finish the sweep and read the paired power table.** The prediction is "no
   meaningful difference". If right, the MLP's value is as *evidence* that the
   learner is not the bottleneck, rather than as a replacement.
3. **Sweep `d`.** `bench_learners.py` already takes `--d 4 8 16 32`, and
   `tune.py` can tune either learner at any `d`. The motivating example is
   `dX=30, dY=10`, and the MLP shares a hidden representation across classes
   where boosting fits effectively one-vs-rest, so its advantage plausibly
   *grows* with `d`. This is the most likely place for the MLP to matter.
   Requires tuning both learners at each `d` first — xgboost tuning at `d=32` is
   the expensive part, so consider lowering `XGB_MAXROUNDS` (equally for both).
4. **Only if 2 and 3 show a real gap**, consider making the learner a documented
   choice in the paper rather than a fixed detail.

---

## Reproducing from scratch

```sh
WT=/Users/harveyklyne/Documents/Github/catci/.claude/worktrees/mlp-learner
export PYTHONPATH="$WT/src:$WT/experiments" OMP_NUM_THREADS=1
PY=/Users/harveyklyne/miniforge3/envs/catci/bin/python

$PY -m pytest $WT/tests -q                                     # 68 pass, ~5 min
$PY $WT/experiments/tune.py --learner mlp --reps 20 --write    # ~4 min
$PY $WT/experiments/bench_learners.py --reps 100 --workers 8   # ~3 min
sh $WT/experiments/bench_sweep.sh                              # ~80 min
$PY $WT/experiments/report_learners.py
```

Raw output lands in `experiments/results/` (gitignored) as parquet plus a
`.provenance.json` sidecar recording config, git commit, package versions, seed
and runtime.

---

## Update 2026-09-26 — sweep finished on minP + vectorised master

Branch rebased onto `3770b79` (minP calibration + vectorised greedy search);
188 tests pass. Full `bench_sweep.sh` rerun at `n_boot = 1000`, all 10 runs
complete. Pre-minP results archived in `experiments/results/pre-minP/`.

**The adaptive size inflation above is gone** — it was the pre-minP calibration
bug (TODO 0b). Size at alpha = 0.05, 1000 reps, SE ~0.007:

| method | lin oracle | lin xgb | lin mlp | sin oracle | sin xgb | sin mlp |
|---|---:|---:|---:|---:|---:|---:|
| tree | 0.049 | 0.048 | 0.052 | 0.041 | 0.049 | 0.044 |
| ordinal | 0.062 | 0.057 | 0.061 | 0.057 | 0.052 | 0.053 |
| max | 0.043 | 0.043 | 0.045 | 0.043 | 0.041 | 0.042 |
| euclid | 0.041 | 0.040 | 0.042 | 0.048 | 0.050 | 0.048 |
| mGCM | 0.174 | 0.190 | 0.186 | 0.189 | 0.187 | 0.185 |
| chi_sq | 0.161 | 0.165 | 0.163 | 0.141 | 0.143 | 0.142 |

**Power: the learner does not matter**, as predicted. Paired (same datasets),
200 reps x 5 strengths. Mean power over the grid, xgb -> mlp:

| setting | method | xgb | mlp | diff |
|---|---|---:|---:|---:|
| lin_lin_step | ordinal | 0.623 | 0.621 | -0.002 |
| lin_lin_step | euclid | 0.484 | 0.490 | +0.006 |
| sin_sin_binary_tree | tree | 0.821 | 0.813 | -0.008 |
| sin_sin_binary_tree | euclid | 0.665 | 0.659 | -0.006 |

No per-strength paired difference exceeds 0.03; one cell of ~70 is > 2 SE
(mGCM, s = 1.0, lin), consistent with chance.

**Conclusion:** keep xgboost as the default. The MLP's value is as evidence that
the propensity learner is not what limits the method at d = 8. Remaining open
question is step 3 above (sweep `d`), the one place the MLP could plausibly matter.
