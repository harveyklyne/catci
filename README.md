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

Extras: `test` (pytest, hypothesis), `learners` (xgboost, scikit-learn),
`experiments` (learners + pandas, pyarrow).

## Usage

```python
from catci import catci_test
from catci.structure import Ordinal, Tree
from catci.learners import xgboost_learner

res = catci_test(
    x, y,                                  # 1-based integer label vectors
    x_structure=Tree.binary(dx),           # or Ordinal(), Cyclic(), Saturated(), Tree.from_parents(...)
    y_structure=Tree.binary(dy),
    z=z,                                   # conditioning variables
    learner=xgboost_learner({"eta": 0.01, "max.depth": 1, "gamma": 2, "nrounds": 163}),
    n_boot=1000,                           # see 'Calibration' below
)
res.p_value, res.statistics, res.partitions
```

Pass `f=`/`g=` instead of `learner=`/`z=` to supply propensities `P(X|Z)`,
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
  learners.py    Z -> P(label|Z) interface + oracle / xgboost learners
  api.py         catci_test: the public entry point

experiments/
  config.py      config-as-data: one resolved Config per figure
  dgp.py         data-generating processes, parametric in d
  methods.py     method registry: name -> p-value on a fitted dataset
  run.py         power grids -> parquet + provenance sidecar
  run_size.py    null-calibration (size) runs
  tuning/        frozen XGBoost hyperparameters, one JSON per marginal setting

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
- **`d` is not a real axis.** `dgp.py` is parametric in `d`, but every config is
  hardcoded to `d = 8`. The motivating example is `dX = 30, dY = 10`.
- **No tuner.** `experiments/tuning/` holds the frozen hyperparameters for
  `n = 1000, d = 8`. The R script that produced them was cluster-specific and
  was not ported, so a new `(n, d)` cannot currently be tuned.
- **Index convention differs from the paper.** The appendix orders `dXdY`-space
  with `k` fastest; the code uses `j` (X) fastest, inherited from R's Kronecker
  layout. Self-consistent, but update formulae (24)-(27) will not line up
  index-for-index with the code. This is a paper-side edit.
- **Paper text.** The appendix says the interaction matrix rows/columns "sum to
  one"; they sum to zero.

## License

MIT — see [LICENSE](LICENSE).
