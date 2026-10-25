# TODO

Open work only. Finished work is recorded where its results live (see "Done" at the
bottom); git history has the old narrative.

**Decisions (2026-09-30):**

* **Real data before simulation.** The real-data studies decide which methods are
  useful in practice; the final simulation grid waits behind them.
* **The primary method is the truncated divisive search** (`SplitSearch(max_levels=k)`).
  The merge search (Algorithm 1) is secondary.
* **`normalise=False` for all our methods.** Only the mGCM competitor studentises,
  and it does so internally (`statistic.mgcm`). `catci_test` defaults to `False`.

## Priority order

1. **Oracle-only search study on the Adult semi-synthetic data.** Which search gives
   the best power with exact propensities: merge, full split, split@k for several `k`,
   against the depth-0 comparators and the competitors. Run over several effect
   directions (real monotone, planted non-monotone, others) and pairs, so no method is
   favoured by construction. This picks `k` and tests the primary method. Queue and
   launch commands: `experiments/SEMISYNTH_RESULTS.md` ("Not run").
2. **Fitted-learner arm** (MLP, then xgb) beside the oracle on the same data and
   bootstrap draws, so the cost of estimating the nuisances separates from the effect
   of the search. Then the rich-`Z` size arm. Needs `tune.py` at the Adult
   `(n, num_class)`.
3. **Calibration checks on real Adult** (4a.1).
4. **Real Adult analysis** (4a.2-4a.3), scored against the literature's established
   dependencies (`related-literature/ci_ground_truth_2026-08.md` §2).
5. **Diabetes-130** (4b), then **ACS** (4c).
6. **Final simulation grid** (block table below) and the 0a write-up.

## Size gate and mGCM follow-ups (from the retired PLAN.md)

Gate before the Adult deep dive: on the synthetic null and the semi-synthetic Adult
null, at oracle propensities, size within MC error of nominal across
`n in {500, 1000, 2000, 5000}`, the real Adult level counts, both `normalise` settings
and every method, OR the failing region characterised and the method documented as
restricted to its complement.

* Re-run the size grid and regenerate every power number; pre-minP numbers used the
  inflated calibration.
* Re-check semisynth cells with the mGCM mechanism in hand. Unexplained: mGCM 0.060 at
  `n = 1000` -> 0.128 at `n = 5000` (Adult Occupation x Income). `pool_rare` takes `n`,
  so `dx` differs across the rows; check that before calling it a second phenomenon.
* Decide mGCM's fate: drop the studentisation, studentise per draw, or retire it as a
  comparator with the reason documented. Published mGCM power numbers are
  uninterpretable until then.
* Settle `normalise`: the decision above says `False`; check `catci_test`'s default and
  every config agree (True was worse on unpooled Occupation x Income: greedy
  0.095 -> 0.128).
* Optionally the uncentred multiplier bootstrap (robustness to support/skewness, not a
  size fix): thread `prod_mat` through `TSigma` and `bootstrap_T`. Do not centre.
* Regression tests fed from data: parametrise by `rank(Sigma)/n` and assert size (a
  `d = 12, n = 500` test would have caught mGCM); assert the leakage diagnostic; tighten
  `test_calibration_under_gaussian` (`|rate - alpha| < 0.05` allows 0.10) and
  `test_tree_calibrates_too` (only checks mean p); parametrise by number of depths.
* Report `rank(Sigma)`, `rank(Sigma)/n` and the leakage diagnostic per test in the
  appendix.

## Blocker: divisive search on unordered variables

`SplitSearch` supports `Ordinal` and `Tree` only, and only *binary* trees. `Saturated`
raises `NotImplementedError` (`2^(s-1) - 1` bipartitions per group), so any unordered
categorical needs a tree: a published taxonomy where one exists (Census occupation
groups; ICD-9 chapters for 4b), or an exhaustive top-level split for small `d`.
Occupation has a binary 1990-Census tree (`experiments/taxonomies.py`), but two groups
had to be binarised arbitrarily.

**N-ary divisive search** (design in `experiments/DIVISIVE.md`, not implemented) removes
the restriction. Tests required: all-binary n-ary tree reproduces today's paths
bit-for-bit; batched search matches the per-draw loop; constant path length across draws;
prefix property `SplitSearch(k) == SplitSearch()[:k+1]`. Then rebuild the Occupation
tree without the binarisation choices and re-time.

## Paper

* **0a: mGCM.** Resolved, needs writing up. mGCM is our `max` with `T` normalised by the
  estimated covariance (the Shah & Peters construction), and it is miscalibrated: size
  0.215 at `d = 8, n = 1000`, 0.975 at `d = 12, n = 500`, scaling with `rank(Sigma)/n`.
  That is good for us: CLT-based normalisation is unsafe for categorical data as
  `rank(Sigma)/n` grows, whereas bootstrap-on-`Sigma_hat` calibration is not. Mechanism
  in `experiments/SEMISYNTH_RESULTS.md`.
* **0b: rewrite Algorithm 2 as minP.** Draft, proof sketch, tables and citations to
  verify: `experiments/MINP.md`. Notation `M -> S`, `F -> p`, `G -> P_min`.
* **Divisive search**: `main.tex` does not describe it yet. Include the sub-linear cost
  in `B` (`experiments/DIVISIVE.md`, cost section) as a paragraph in the appendix
  "Implementation" section beside "Fast greedy merging" (`sect:fastgreedymerging`).
* **Item 9 (one line, limitation):** multi-label categoricals break the framework rather
  than extend it (`X` becomes a point in `{0,1}^d`; power-set encoding blows up,
  per-label binary encoding discards the joint).

## Final simulation grid (0b)

Post-minP numbers so far were run piecemeal with different `reps`/`n_boot`. Re-run
everything the paper cites, once, on one commit, `n_boot = 1000`.

| block | what it shows | setting | methods | reps |
|---|---|---|---|---|
| **A. Size** | every method holds level; mGCM / `chi_sq` do not | `size_config`, margins `lin_lin`, `sin_sin`, `lin_sin`; `n ∈ {500, 1000, 2000}`, `d = 8`; learners `mlp`, `xgb`, `oracle` | `tree`, `ordinal`, `greedy`, `*_bonf`, `*_exact`, `max`, `euclid`, `mGCM`, `chi_sq`, `ankan`, `multinomial` | 1000 |
| **B. Power, d = 8** | adaptive beats depth-0 and competitors | `power_config` `lin_lin_step`, `sin_sin_binary_tree`, `n = 1000`, strengths 0.2-1.8; `mlp` + `xgb` | as A, plus misspecified cross (`tree` on `step`, `ordinal` on `binary_tree`): the "not much worse when the structure is wrong" claim has no number yet | 200 |
| **C. d axis** | merging advantage widens with `d`; mGCM / `chi_sq` break | `d_grid` `lin_lin_step`, `n = 2000`, `dy = 4`, `dx ∈ {8,…,256}`, strengths `{0, 0.6, 1.2}`; `oracle` + `mlp` (needs `tune.py` at `n = 2000` per `dx`) | `ordinal`, `euclid`, `max`, `mGCM`, `chi_sq`, `ankan` | 200 |
| **D. Exact vs approx** | "same p-values in almost every case, ~1/28 the time" | B's two DGPs at strengths `{0, 0.8, 1.4}`, paired (add `*_exact` to B) | `*` vs `*_exact`, plus timings | 200 |
| **E. Search alternatives** (Gaussian limit) | greedy vs random / split / beam; divisive truncation | `search_study.py`, `power_divisive.py` cells with uniform reps and `n_boot` | as in `SEARCH_RESULTS.md`, `DIVISIVE.md` | 500 |

C is the expensive block (`p = dx·dy` up to 1024). Re-check the `ADULT_SEMISYNTH.md`
cells once 4a is done.

## Real data

Cross-cutting, all three datasets: decide how to clean and pre-process; spot check the
existing methods and identify suspicious rejections/non-rejections; apply our method with
tuned MLP and XGBoost (CV). Do we reject a superset of the methods with good size control?

Dataset order is set by whether the dataset supplies an *external, published,
outcome-blind* partition to score the learned one against
(`related-literature/ci_ground_truth_2026-08.md` §6.1).

### 4a. UCI Adult

What exists is *their* method (`experiments/ADULT.md`) plus the semi-synthetic study
(`ADULT_SEMISYNTH.md`, `SEMISYNTH_RESULTS.md`). Three deliverables of different value:

1. **Tier-1 calibration.** `education` vs `education-num` (an exact bijection: must
   abstain or flag determinism, not silently reject) and `sex` given `relationship`
   (perfect separation on 45% of rows: must not blow up). Cheap credibility; nobody does it.
2. **Head-to-head on Fig. 1b.** Six hypotheses at `n = 1000` with their discretisation and
   variable set, our verdicts beside theirs.
3. **The `fairadapt` DAG.** All six testable implications fail; three have never been
   tested. A *finding*, and the strongest single result in the repo.

### 4b. Diabetes 130-US Hospitals, 1999-2008 (UCI)

The only candidate where the reported partition is itself the scientific output. `diag_1`
has three nested published outcome-blind partitions (ICD-9 chapters ⊂ CCS-285 ⊂
phecodes), so the learned tree can be scored against something external. Also a published
stratum-specific hypothesis (Strack's diagnosis × HbA1c interaction), a hard structural
zero for calibration, and `A1Cresult`, where the right split is almost certainly {None}
vs {Norm, >7, >8}: a selection indicator masquerading as a level. Needs the n-ary divisive
search (ICD-9 chapters have ~19 children under the root). ~30% of rows are repeat
patients: restrict to first encounters, which also reproduces Strack's sample.

### 4c. ACS PUMS via `folktables`

Real value is `OCCP` × `INDP` (SOC and NAICS): a published tree on both sides at once.
Scale: `dX·dY` in the hundreds × hundreds, `n` in the millions, rich `Z` where
stratification competitors cannot run. Costs: no external CI ground truth, and five
survey-design objections. Pre-empts "Adult is 1994" (`ACSIncome` mirrors Adult's filters).
`form_t_sigma` will need row-chunking at this scale (README, Known gaps).

## Smaller open items

* `catci_test` default `n_boot = 100` -> 1000 (now affordable).
* **Weighted minP over depths** on the ordinary merge path (more weight on coarse
  levels): directly tests the multiplicity explanation of the divisive gain; possible
  paper remark.
* The divisive `alt` (merge + split union) cell was never run; prediction: tracks merge.
* Numba merge kernel, only if `p` reaches the thousands (beyond `p ~ 4096` a factored,
  `Sigma`-free form would be needed).

## Done (results live in the linked doc)

* Calibration rewritten as minP (0b), `bonferroni_pvalue` comparator: `experiments/MINP.md`.
* Tuner, `d` as a config axis, `d`-axis pilot (advantage widens with `d`):
  `experiments/DESIGN_NOTES.md`.
* MLP learner (default; xgboost matches it at `d = 8`): `experiments/MLP_LEARNER.md`.
* `Tree.from_parents`, `Cyclic`, n-ary mask fix: `structure.py`.
* Search alternatives (random, split, beam; greedy suffices): `experiments/SEARCH_RESULTS.md`.
* Divisive search, its cost, and the vectorised search (6-58x per draw, `n_jobs`):
  `experiments/DIVISIVE.md`.
* Exact chi criterion: `catci_test(statistic="exact")`, ~28x slower, p-values equal in
  17 of 18 paired reps; the CDF evaluation (91% of cost) cannot be cached across draws.
* Adult: Ankan & Textor port (`experiments/ADULT.md`), semi-synthetic design
  (`ADULT_SEMISYNTH.md`) and post-minP results (`SEMISYNTH_RESULTS.md`).
