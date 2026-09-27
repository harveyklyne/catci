# catci

**CAT**egorical **C**onditional **I**ndependence testing — a test of
`X ⟂ Y | Z` for categorical `X`, `Y` that adaptively gains power by greedily
merging labels, calibrated with a parametric double bootstrap.

This was originally an R package. It is now Python only; the R implementation
was removed once the port was complete and differential-tested against it. The
frozen R source lives on the `r-frozen-oracle` git tag, and the fixtures it
generated are still what the test suite checks against. A thin R wrapper over
this package may follow at publication time.

## Install

```sh
conda activate catci                 # Python 3.13
pip install -e ".[test,experiments]"
pytest
```

The core install includes scikit-learn, for the default MLP propensity learner.
Extras: `test` (pytest, hypothesis), `learners` (xgboost, for the boosted
alternative), `experiments` (xgboost + pandas, pyarrow).

## Usage

```python
from catci import catci_test
from catci.structure import Ordinal, Tree

res = catci_test(
    x, y,                                  # 1-based integer label vectors
    x_structure=Tree.binary(dx),           # or Ordinal(), Cyclic(), Saturated(), Tree.from_parents(...)
    y_structure=Tree.binary(dy),
    z=z,                                   # conditioning variables
    n_boot=1000,                           # see 'Calibration' below
)
res.p_value, res.statistics, res.partitions
```

The propensities `P(X|Z)`, `P(Y|Z)` are fitted by a small neural network,
`mlp_learner()`, at hyperparameters tuned for `n = 1000, d = 8`
(`learners.DEFAULT_MLP_PARAMS`). Pass `learner=` to override them or to swap
in gradient boosting, without touching anything else:

```python
from catci.learners import mlp_learner, xgboost_learner

learner = mlp_learner({"hidden_layer_sizes": [32], "alpha": 10.0})
learner = xgboost_learner({"eta": 0.01, "max.depth": 1, "gamma": 2, "nrounds": 163})
```

The two give the same size and power at `d = 8` (paired sweep, `report_learners.py`);
the experiments run both.

Pass `f=`/`g=` instead of `z=` to supply propensities `P(X|Z)`,
`P(Y|Z)` directly — the oracle path, which separates "does the test calibrate"
from "did the regression fit well".

## Layout

```
src/catci/
  gcm.py         form_t_sigma: (x, y, f, g) -> (T, Sigma)          [pure]
  structure.py   Ordinal / Cyclic / Saturated / Tree: permitted_merges()
  merging.py     rank-one update formulae (24)-(27)                [pure]
  statistic.py   ApproxChi init/update/value + depth-0 comparators
  search.py      greedy_search: the adaptive label-merging path
  bootstrap.py   matrix_sqrt + N(0, Sigma) sampling
  calibrate.py   minP calibration of the search path + adaptive_pvalue
  learners.py    Z -> P(label|Z) interface + mlp (default) / xgboost / oracle learners
  tuning.py      K-fold CV tuner for the xgboost learner (sims and real data)
  api.py         catci_test: the public entry point

experiments/
  config.py      config-as-data: one resolved Config per figure
  dgp.py         data-generating processes, parametric in d
  methods.py     method registry: name -> p-value on a fitted dataset
  run.py         power grids -> parquet + provenance sidecar, one per learner
  run_size.py    null-calibration (size) runs, one per learner
  run_d_axis.py  power/size as dx grows with n, dy fixed (one run.py block per dx)
  tuning/        tuned hyperparameters (mlp + xgb), one JSON per (n, num_class, setting)
  tune.py        tune either learner for a simulated (n, num_class, setting):
                 held-out mlogloss (the R protocol) or K-fold CV (xgb)
  tune_check.py  score tunings by KL to the true propensities
  bench_learners.py  propensity quality: E_f, the Assumption 1 remainder
  report_learners.py print the size/power/propensity comparison tables
  bench_sweep.sh     drive the whole learner comparison

tests/
  fixtures/      the frozen R oracle (see fixtures/README.md)
  test_*.py      differential + property tests
```

## Design notes

These are the decisions that shaped the port, kept here because the code alone
does not explain why it is shaped this way.

### `permitted_merges` is the only structure operation

The R package represented search structure as three things threaded in
lockstep — a search string (`"ordinal"`/`"greedy"`/`"tree"`), a ragged `trees`
list, and a ragged `categories` list — plus a separate `get_num_levels` that had
to agree with the merge loop by hand.

Here a single object per variable answers one question:

```python
structure.permitted_merges(partition) -> [(i, j), ...]
```

The row count is just `len(permitted_merges(...))`, so it cannot disagree with
the loop, and `get_num_levels` does not exist. Every structure returns `[]` once
the partition has 2 or fewer groups, so the `(d > 2)` guard holds by
construction rather than as three parallel copies of an `if`.

This matters because that split representation *was* the bug described below.

### The tree method used to ignore the tree

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

### No cross-fitting

Propensities are fitted on the full sample (`learners.fit_propensities`). The
earlier `nfolds` cross-fitting option was removed: the test calibrates against
the fitted propensities themselves, so sample splitting bought nothing while
costing an `nfolds`-fold slowdown and an extra RNG stream per replicate.

### Statistics are init/update/value triples

So the non-adaptive comparators (`max`, `euclid`, `mGCM`) are the same code path
at search depth 0. There is one live statistic, `ApproxChi` — no statistic zoo.

### Calibration is a minP test

The search returns a path of `L = dx + dy - 3` statistics — one per coarsening it
visits — each a valid test statistic for the same null. Using whichever is most
extreme is a multiple testing problem, so `calibrate.py` solves it as **minP**
calibrated by resampling (Westfall & Young 1993; Romano & Wolf 2005), with the
first stage being Beran (1988) prepivoting:

1. Pool the observed path with the `B` bootstrap paths into `B + 1` exchangeable
   paths, and rank every one of them against the *other* `B`, by one identical
   leave-one-out rule. That turns each level into a marginal p-value.
2. Take each path's minimum over levels, and calibrate the observed minimum
   against the bootstrap minima — the same rule again.

Because every path is treated identically, the result is exactly uniform on
`{1/(B+1), ..., 1}` under the null at finite `B`. The R implementation was not:
it ranked the observed path against `B` draws *excluding itself* but each draw
against `B` *including itself*, so every bootstrap quantile was scaled by
`B/(B+1)` and the `max` over `L` levels compounded it as `(B/(B+1))^L`. Size ran
from 0.049 at `L = 1` to 0.358 at `L = 40`. `tests/test_calibrate.py` pins the
fix, parametrised by `L`.

Two practical consequences of the `1/(B+1)` grid:

- **Ties are broken at random, and that is load-bearing.** Each level puts one
  path at the floor `1/(B+1)`, so up to `L` paths tie there. Breaking those ties
  conservatively makes the test unable to reject at all once `L` approaches
  `alpha * (B + 1)`.
- **`n_boot` bounds power, not just resolution.** The test is exact at any
  `n_boot`, but at `d = 8` (`L = 13`) and `n_boot = 100` it recovers only about
  half the power it reaches at `n_boot = 1000`, which is where the curve
  flattens. `catci_test` still defaults to 100 for a cheap smoke test; the
  experiment configs use 1000.

`bonferroni_pvalue` is the comparator: the same statistic path under simple FWER
control, `min(1, L * min_l p_l)`. Its floor is `L/(B+1)`, so unlike minP it cannot
reject at level `alpha` at all unless `B >= L/alpha`.

### DGP fixes carried into the port

- A dead resampling branch in `simulate_data` recomputed `dep_pdf` and
  reassigned `x`/`y`, discarding the preceding work.
- Cross-family settings (e.g. `xsetting="sin"`, `ysetting="lin"`) crashed.
- `d = 8` was hardcoded into the `sin`/`sig`/`hat` generators. `dgp.py` is
  parametric in `d` throughout.

### Tuning, and what breaks as `d` grows

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

## Testing

`tests/` differential-tests every deterministic seam against the frozen R
oracle in `tests/fixtures/catci_fixtures.json`, plus Hypothesis property tests
(fast rank-one update == dense recompute), a calibration-under-Gaussian check,
and `test_calibrate.py`, which feeds the calibration exchangeable paths directly
and asserts uniformity at each `L`. RNG streams do not match across languages,
so the fixtures deliberately
pin only pure input to output maps, never bootstrap or randomised tie-break
paths. See `tests/fixtures/README.md` for the JSON conventions.

## Known gaps

- **Results are stale.** Every figure in `experiments/results-r-legacy/` came
  from R with cross-fitting. Nothing in the paper yet comes from this code; the
  grid needs re-running with `experiments/run.py`.
- **`dx` in the hundreds is feasible but not cheap.** At `dx = 256, dy = 4`
  a replicate costs ~4.5 min at `n_boot = 200` and ~5x that at 1000, so a full
  `d` grid needs a cluster or reduced reps at the top end.
- **No learner for large `d`.** See "Tuning, and what breaks as `d` grows":
  at `n ~ 2000, d >= 64` both xgboost and multinomial logistic regression are
  worse than uniform. The MLP has not been tried there yet.
- **`form_t_sigma` materialises the `n x dx*dy` product matrix.** At
  application scale (`n ~ 7e4, dx*dy ~ 1e3`) that is ~1 GB, doubled by
  `np.cov`. It will need row-chunking before item 4.
- **Index convention differs from the paper.** The appendix orders `dXdY`-space
  with `k` fastest; the code uses `j` (X) fastest, inherited from R's Kronecker
  layout. Self-consistent, but update formulae (24)-(27) will not line up
  index-for-index with the code. This is a paper-side edit.
- **Paper text.** The appendix says the interaction matrix rows/columns "sum to
  one"; they sum to zero.

## License

MIT — see [LICENSE](LICENSE).
