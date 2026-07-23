# catci (Python)

Python port of the R `catci` package — conditional independence testing for
structured categorical data. Design and rationale: `../CODE_REVIEW.md` §4 & §7.

This lives as a subdirectory of the R repo during the port so the frozen R
oracle in `../fixtures/` stays co-located for differential testing. It can be
`git subtree split` into a standalone repo once the interface stabilises.

## Layout

```
src/catci/
  gcm.py         form_t_sigma: (x, y, f, g) -> (T, Sigma)          [pure]
  structure.py   Ordinal / Saturated / Tree: permitted_merges()    [the §4.1 redesign]
  merging.py     rank-one update formulae (24)-(27)                [pure]
  criteria.py    ApproxChi init/update/value + depth-0 comparators
  search.py      greedy_search: the adaptive label-merging path
  bootstrap.py   matrix_sqrt + N(0, Sigma) sampling
  calibrate.py   double_bootstrap_pvalue (2d/2e fixed) + adaptive_pvalue
  learners.py    Z -> P(label|Z) interface + oracle learner + crossfit
  api.py         catci_test: the public entry point
```

`experiments/` (configs, runners, figures) is deferred to a later step.

## Design notes

- **`structure.permitted_merges(partition)`** is the single operation that
  replaces the R `(search string, trees, categories)` triple and its separate
  `get_num_levels`. The row count is just `len(permitted_merges(...))`, so the
  class of bug behind finding #1 (search string and structure disagreeing) is
  unrepresentable here.
- **Criteria** are `init`/`update`/`value` triples so the non-adaptive
  comparators are the same code path at search depth 0. One live criterion
  (`ApproxChi`); no metric zoo, per §7.
- **Calibration fixes:** 2d (single-metric fall-through) and 2e (mismatched
  `B` vs `B+1` normalisation) from the review are applied in `calibrate.py`.

## Testing

```sh
conda activate catci
cd python && pip install -e ".[test]" && pytest
```

`tests/` differential-tests every deterministic seam against the frozen R
oracle (`../fixtures/catci_fixtures.json`), plus Hypothesis property tests
(fast update == dense recompute) and a calibration-under-Gaussian check. See
`../fixtures/README.md` for the oracle's JSON conventions.
