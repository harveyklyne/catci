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
  search.py      greedy_search (the default); beam, random, fixed-path and
                 divisive (split-down) searches as comparators
  blocks.py      block sums of T and Sigma for the divisive search
  bootstrap.py   matrix_sqrt + N(0, Sigma) sampling
  calibrate.py   minP calibration of the search path + adaptive_pvalue
  learners.py    Z -> P(label|Z): mlp (default) / xgboost / oracle
  tuning.py      K-fold CV tuner for the xgboost learner
  api.py         catci_test: the public entry point

experiments/     simulation + real-data drivers (run.py, run_size.py,
                 run_d_axis.py, search_study.py, power_divisive.py, tune.py,
                 adult.py, ankan_textor.py, run_adult.py, adult_semisynth.py,
                 run_semisynth.py, sweep_semisynth.py, report_*.py, ...) and
                 their write-ups:
  MINP.md               minP calibration (Algorithm 2), proofs, tie/B tables
  SEARCH_RESULTS.md     greedy vs beam / random / split
  DIVISIVE.md           divisive search, its cost, vectorised speedups
  MLP_LEARNER.md        MLP vs xgboost propensity learner
  ADULT.md              Ankan & Textor reproduction on Adult
  ADULT_SEMISYNTH.md    semi-synthetic Adult study: design
  SEMISYNTH_RESULTS.md  semi-synthetic results (post-minP)
  DESIGN_NOTES.md       tree-ignored-the-tree history, DGP fixes, d-axis study
  ankan_textor_code_sources_2026-08.md   provenance of the A&T method

tests/
  fixtures/      the frozen R oracle (see fixtures/README.md)
  test_*.py      differential + property tests
```

## Design notes

- **`permitted_merges` is the only structure operation.**
  `structure.permitted_merges(partition) -> [(i, j), ...]`. The row count is
  `len(...)`, every structure returns `[]` at 2 or fewer groups, and nothing
  threads a search string / trees / categories in lockstep (the R design, which
  hid a bug where `"tree"` silently ran as `"ordinal"`: see `DESIGN_NOTES.md`).
- **No cross-fitting.** Propensities are fitted on the full sample; the test
  calibrates against the fitted propensities, so splitting bought nothing.
- **Statistics are init/update/value triples**, so the non-adaptive comparators
  (`max`, `euclid`, `mGCM`) are the same code path at search depth 0. One live
  statistic, `ApproxChi`.
- **Calibration is a minP test** over the `L = dx + dy - 3` statistics on the
  search path: pool observed + `B` bootstrap paths, rank each against the other
  `B` by one leave-one-out rule, take path minima, rank again. Exactly uniform on
  `{1/(B+1), ..., 1}` under the null. Ties must be broken at random, and
  `n_boot` bounds power (d = 8: n_boot=100 gives about half the power of 1000;
  experiment configs use 1000). `catci_test` still defaults to 100. Derivation,
  proofs and the R bug it fixed: `experiments/MINP.md`.
- **Greedy is the search; the rest are comparators.** Any path map is valid if
  every draw goes through the same one, or is randomised independently of `T`.
  None beats greedy except truncated divisive search on coarse signals; see
  `experiments/SEARCH_RESULTS.md` and `experiments/DIVISIVE.md`.
  `catci_test(..., search=SplitSearch(k))` selects the divisive search;
  `MergeSearch()` is the default.
- **Large `d`.** Learners that treat labels as unrelated classes fail at
  `d >= 64` and the MLP is only tuned at `d = 8`, so large-`d` simulations use
  `learner="oracle"`. Merging's advantage over unmerged tests grows with `d`;
  details and tuning notes in `experiments/DESIGN_NOTES.md`.

## Testing

`tests/` differential-tests every deterministic seam against the frozen R
oracle in `tests/fixtures/catci_fixtures.json`, plus Hypothesis property tests
(fast rank-one update == dense recompute), a calibration-under-Gaussian check,
and `test_calibrate.py`, which feeds the calibration exchangeable paths directly
and asserts uniformity at each `L`. RNG streams do not match across languages,
so the fixtures pin only pure input-to-output maps, never bootstrap or
randomised tie-break paths. See `tests/fixtures/README.md`.

## Known gaps

- **`form_t_sigma` materialises the `n x dx*dy` product matrix.** At application
  scale (`n ~ 7e4, dx*dy ~ 1e3`) that is ~1 GB, doubled by `np.cov`; it needs
  row-chunking.
- **`dx` in the hundreds is feasible but not cheap** (~4.5 min per replicate at
  `dx = 256, dy = 4, n_boot = 200`).
- **Index convention differs from the paper.** The appendix orders `dXdY`-space
  with `k` fastest; the code has `j` (X) fastest (R's Kronecker layout).
  Paper-side edit.
- **Paper text.** The appendix says the interaction matrix rows/columns "sum to
  one"; they sum to zero.

Real-data and semi-synthetic status lives in `TODO.md`.

## License

MIT — see [LICENSE](LICENSE).
