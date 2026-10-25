# Design history and d-axis findings

Moved out of README.md. Background on decisions and measurements that the code
alone does not explain. See `README.md` for the current design summary.

## The tree method used to ignore the tree

In the R package, `query_lookup("tree")` built real tree objects and then ran
the search with `xsearch = ysearch = "ordinal"`. The trees were silently
discarded: `"tree"` was a synonym for `"ordinal"`, and every published `tree`
curve was really an ordinal curve. Turning the tree path on then exposed a
second, latent bug — the `"tree"` branch of `get_num_levels` was missing the
`(d > 2)` guard that the `greedy`/`ordinal` branches had, so `num_rows`
over-counted once a dimension reached two still-sibling groups.

Both were fixed in R before the oracle was frozen, and the fix was validated on
`lin_lin_binary_tree` (n = 1000, d = 8, 200 reps, paired on identical data):

| strength | tree | ordinal | gap |
|---:|---:|---:|---:|
| 0.2 | 0.050 | 0.045 | +0.005 |
| 0.6 | 0.280 | 0.145 | **+0.135** |
| 1.0 | 0.785 | 0.615 | **+0.170** |
| 1.4 | 0.975 | 0.920 | +0.055 |
| 1.8 | 0.995 | 0.995 | 0.000 |

Mean power 0.628 (tree) vs 0.547 (ordinal). A gap that peaks in the middle and
vanishes at both ends is the signature of a real power gain rather than Monte
Carlo noise, so the method's central claim now has direct evidence. As a sanity
check, the old buggy `"tree"` column matched the fixed run's genuine `ordinal`
column to three decimals — confirming it really had been ordinal all along.

In this codebase the bug class is unrepresentable, per the section above. It is
pinned three ways: `test_search_paths_match_oracle` (the full 13-level tree
path, values *and* partitions, against R), `test_permitted_merges_match_oracle`
(pair lists at d in {2, 4, 8}, plus the d = 2 guard), and
`test_tree_differs_from_ordinal` (the assertion the original bug would fail).
A calibration test would *not* have caught this — ordinal and tree search are
each valid tests, so both calibrate. The diagnostic has to be that on the same
`(T, Sigma)` the two searches visit different merge paths.

One cost note from the R fix does **not** carry over: genuine tree search ran
~4x slower than the buggy path in R, because the pure-R sibling bookkeeping
re-executed inside every bootstrap draw. Measured here at d = 8, tree search is
2.44 ms vs ordinal's 4.43 ms — nearly 2x *faster*, since the tree offers far
fewer candidate pairs per level. No optimisation is needed before a full grid
run.

## DGP fixes carried into the port

- A dead resampling branch in `simulate_data` recomputed `dep_pdf` and
  reassigned `x`/`y`, discarding the preceding work.
- Cross-family settings (e.g. `xsetting="sin"`, `ysetting="lin"`) crashed.
- `d = 8` was hardcoded into the `sin`/`sig`/`hat` generators. `dgp.py` is
  parametric in `d` throughout.

## Tuning, and what breaks as `d` grows

`catci.tuning.tune_xgboost` picks `(eta, max_depth, gamma, nrounds)` by K-fold
CV log-loss. `nrounds` is read off the mean validation curve rather than
gridded. Several simulated replicates are tuned in one `xgb.cv` call, with
folds kept inside each replicate, so the training size is `(K-1)/K * n` — at
`n = 1000, K = 5` that is the R tuner's `n_tr = 800`. `experiments/tune.py
--learner xgb --protocol cv` draws tuning data from `dgp.simulate_marginal`: the
interaction preserves both margins, so one tuning per `(n, num_class, setting)`
serves every interaction, strength and partner dimension. The default
`--protocol holdout` is the R tuner's own held-out-mlogloss protocol, which is
also how the MLP was tuned; both write into the same per-setting JSON, one key
per learner.

At `n = 1000, d = 8`, the cheap `FAST_GRID` (eta 0.1, 9 points) matches the
frozen R tuning (eta 0.01, 28 points) on mean `KL(f || f_hat)` to the true
propensities, and fits 10-40x faster (`tune_check.py`, 20 reps, SE ~0.0007):

| setting | R tuning | fast grid | Z-blind (class freq) |
|---|---:|---:|---:|
| sin | 0.0173 | 0.0141 | 0.0200 |
| sig | 0.0122 | 0.0122 | 0.0192 |
| lin | 0.0176 | 0.0185 | 0.0331 |
| vee | 0.0190 | 0.0183 | 0.0424 |
| hat | 0.0236 | 0.0234 | 0.0733 |

**Learners that treat labels as unrelated classes stop working at large `d`.**
In these DGPs `X | Z` carries about 0.02 nats/observation over uniform at every
`d`. That is the same order as the `(d-1)/2n` it costs merely to estimate `d`
class frequencies. At `n_train = 1600, d = 64` the held-out log-loss is 4.137
(truth), 4.159 (uniform), 4.184 (class frequencies), and 4.184 for the best
xgboost round: xgboost initialises from the noisy frequencies and fits `d`
unrelated per-class ensembles. Multinomial logistic regression, linear or
spline, is also worse than uniform at `d >= 64` for every `C`, because its
intercepts are unpenalised. The tuner correctly reports "stop after ~4
rounds". A learner that shares strength across neighbouring labels (the `lin`
pmf is a smooth ramp in the label index) would not have this problem. That is
the same structure the test exploits, and it is item 2's problem. Item 2's MLP
learner (now the default) shares a hidden layer across classes, but it has only
been tuned and compared at `d = 8`, where it matches xgboost. Until it is swept
over `d`, large-`d` simulation runs use `learner="oracle"`, which isolates the
test from the regression.

**Does the merging advantage widen with `d`? Yes (pilot).** Setup: `lin_lin_step`,
`n = 2000`, `dy = 4`, oracle propensities, `n_boot = 200`, `alpha = 0.05`, 100
reps per ordinal cell and 200 per comparator cell. The per-observation signal
is flat in `d` for `lin`. Reproduce with `run_d_axis.py --learner oracle` and
`report_pilot.py`.

| strength | method | dx=8 | 16 | 32 | 64 | 128 | 256 |
|---|---|---:|---:|---:|---:|---:|---:|
| 0 (size) | ordinal | 0.040 | 0.030 | 0.040 | 0.060 | 0.040 | 0.033 |
| | euclid | 0.025 | 0.045 | 0.020 | 0.025 | 0.030 | 0.010 |
| | max | 0.025 | 0.035 | 0.040 | 0.010 | 0.045 | 0.010 |
| | chi_sq | 0.040 | 0.090 | 0.110 | **0.405** | **0.995** | **1.000** |
| | mGCM | 0.055 | 0.080 | 0.105 | **0.320** | **0.640** | 0.020\* |
| 0.6 | ordinal | **0.640** | **0.600** | **0.570** | **0.490** | **0.460** | **0.367** |
| | euclid | 0.450 | 0.305 | 0.230 | 0.135 | 0.065 | 0.020 |
| | max | 0.215 | 0.100 | 0.070 | 0.070 | 0.035 | 0.020 |
| 1.2 | ordinal | **1.000** | **1.000** | **1.000** | **0.980** | **0.970** | **0.633** |
| | euclid | 0.995 | 0.975 | 0.870 | 0.640 | 0.410 | 0.220 |
| | max | 0.715 | 0.430 | 0.245 | 0.145 | 0.095 | 0.020 |

The `dx = 256` column has 30 reps per ordinal cell (SE up to ~0.09) and 100
per comparator cell. Ordinal search stays calibrated across the grid and loses
little power up to `dx = 128`. At `dx = 256` (`dx*dy = 1024`, about 2
observations per cell at `n = 2000`) its power at strength 1.2 drops to 0.63.
That is still ~3x Euclid's there, and ~18x at strength 0.6. Over the grid, at
strength 0.6, its power over Euclid, the best calibrated unmerged test, grows
from 1.4x to 18x. \*`mGCM` at `dx = 256` is an artefact: some cells are
empty, so `diag(Sigma) = 0` and its studentisation divides by zero. The
pseudo-inverse `chi_sq` competitor and `mGCM` lose
size control as `dx*dy` approaches `n`, even with oracle propensities (TODO
0a's `rank(Sigma)/n` effect), so their power at large `dx` is not
interpretable. `ankan` has no power against this interaction at any `d`.

**Cost.** One replicate is dominated by the search over all `n_boot + 1`
paths. With the vectorised search, a whole replicate's ordinal paths at
`dy = 4, n_boot = 200` take ~1.4 s (`dx = 32`), ~6 s (64), ~32 s (128) and
~4.5 min (256), which is roughly cubic in `dx` (measured under heavy machine
load). Before vectorisation `dx = 64` was ~4.5 min and `dx = 256` ~4.5 h.
The dense `Sigma`
is `(dx*dy)^2`, so symmetric `d` in the hundreds is out of reach whatever the
speed. The applications are asymmetric (`dx` in the hundreds, `dy` a handful),
which is why `dx` and `dy` are separate config axes.
