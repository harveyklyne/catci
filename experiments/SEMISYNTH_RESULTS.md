# Semi-synthetic adult: results log (2026-09-30, minP, normalise=False)

Running log of the semi-synthetic adult runs on branch `worktree-adult-integration`.
Design: `ADULT_SEMISYNTH.md`; code: `adult_semisynth.py`, `run_semisynth.py`
(all methods, oracle and fitted propensities), `sweep_semisynth.py` (oracle, depth `k`
x bootstrap size `B`), `report_semisynth.py` (prints the tables below from
`results_semisynth/`, which is gitignored).

Conventions: alpha = 0.05; `Z = W = {Age, Sex}`; rows drawn i.i.d. from the 30,162-row
population; rare levels of unordered variables pooled at expected count < 5.
`split@k` = divisive search truncated after `k` splits (the primary method);
`merge` = Algorithm 1; `mergecoarse@k` = minP over the merge path's coarsest `k+1`
levels (**historical, withdrawn 2026-10-25**: it kept the last `k+1` rows of the full
merge path, which is not a truncation of any search; the code is gone, replaced by
`merge@k` = `MergeSearch(k)`, the first `k` merges from the finest partition, and the
random baselines `rsplit@k` / `rmerge@k`, see `RANDOM_SEARCH.md`; the mergecoarse
numbers below are kept as recorded); `max`, `euclid`, `mGCM` depth-0; `chi_sq` pseudo-inverse chi-square;
`at_typed` / `at_cat` Ankan & Textor typed as the data types it / all categorical;
`strat_chi2` Pearson chi-square stratified over Z cells. CPU = `time.process_time`
in single-threaded workers on a 6-core Apple machine (details in each provenance
sidecar).

## 1. Null re-check, B = 1000, 500 reps (SE 0.010)

Oracle and MLP propensities on the same data. Every catci method holds level; the
competitors do not.

| pair | learner | merge | split | split@2 | split@4 | max | euclid | mGCM | chi_sq | at_typed | at_cat | strat_chi2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Education x Income | oracle | 0.048 | 0.046 | 0.048 | 0.058 | 0.044 | 0.042 | 0.064 | **0.080** | 0.032 | **0.080** | 0.060 |
| | mlp | 0.048 | 0.052 | 0.048 | 0.054 | 0.048 | 0.042 | **0.072** | **0.080** | 0.044 | **0.084** | |
| Occupation x Income (Saturated merge) | oracle | 0.048 | | | | 0.060 | 0.044 | **0.080** | **0.092** | **0.092** | **0.092** | **0.078** |
| | mlp | 0.060 | | | | 0.054 | 0.054 | **0.086** | **0.084** | **0.084** | **0.084** | |

* The pre-minP over-rejection recorded in `ADULT_SEMISYNTH.md` (0.098 on the
  synthetic DGP) is gone.
* `at_cat` needed the pinv fallback (rank-deficient `Sigma_d`) in 23% of MLP
  replicates and none of the oracle ones; it over-rejects with or without it.
* HoursPerWeek x Education crashed in the MLP arm (replicate 435): an Education level
  absent from the sample got fitted propensity exactly 0, and the merge search scored
  an identically-zero partition as `nan` for some draws only. Fixed in `a9c0ba9`
  (`approx_chi_array` scores such partitions 0); that pair is re-run in the learner arm.
* These runs predate the common-random-numbers fix (`a9c0ba9`): oracle and MLP did
  not share bootstrap draws. Each arm is valid on its own; they are just not paired.

## 2. Depth x B sweep, oracle

### Education x Income, n = 1000 -- size, 1000 reps (SE 0.007)

| | 100 | 300 | 1000 | 3000 | 10000 |
|---|---|---|---|---|---|
| split@1 | 0.051 | 0.054 | 0.058 | 0.061 | 0.061 |
| split@2 | 0.052 | 0.056 | 0.057 | 0.060 | 0.059 |
| split@3 | 0.055 | 0.058 | 0.064 | 0.062 | 0.064 |
| split@4 | 0.060 | 0.057 | 0.063 | 0.066 | 0.067 |
| split@6 | 0.057 | 0.055 | 0.055 | 0.064 | 0.064 |
| split@8 | 0.055 | 0.054 | 0.057 | 0.061 | 0.062 |
| split (k = 14) | 0.059 | 0.053 | 0.053 | 0.055 | 0.057 |
| merge | 0.057 | 0.052 | 0.053 | 0.056 | 0.058 |
| mergecoarse@1..8 (historical) | 0.054--0.072 | | | | |
| max / euclid | 0.045--0.051 | | | | |
| mGCM | 0.062--0.066 | | | | |

**Mild excess for the truncated searches (~1--1.5 points at B = 10k), growing with
B.** Each cell is only 1--2 SE high, but the truncated cells are high together (they
share data, so they are correlated). minP is exact with respect to the Gaussian
bootstrap, so excess that appears as `B` resolves the calibration is the finite-`n`
gap between `T` and `N(0, Sigma_hat)`. Coarse truncated levels plausibly feel it more
because they load on rare Education levels. Check: HoursPerWeek x Education at
n = 2000 (queued).

CPU per search (median of 20 reps):

| | B = 100 | 1000 | 10000 | ms/draw at 10k | cost ~ B^ |
|---|---|---|---|---|---|
| split@1 | 0.016 | 0.023 | 0.134 | 0.013 | 0.47 |
| split@2 | 0.054 | 0.103 | 0.274 | 0.027 | 0.33 |
| split@4 | 0.148 | 0.548 | 1.381 | 0.138 | 0.48 |
| split@8 | 0.304 | 1.378 | 4.735 | 0.473 | 0.59 |
| split (k = 14) | 0.385 | 1.674 | 5.652 | 0.565 | 0.58 |
| merge | 0.032 | 0.255 | 2.656 | 0.266 | 0.97 |

The divisive search scores each partition once however many draws reach it, so its
cost grows sub-linearly in `B`; the merge search is linear. Calibration is ~3 ms.

### Education x Income, n = 1000 -- power, real (monotone) direction, 300 reps per lam

Rejection at B = 10000 (SE <= 0.020; paired differences vs merge in brackets):

| lam | split@1 | split@2 | split@4 | split (14) | merge | mergecoarse@2 (historical) | max | euclid | mGCM |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 0.317 (+.033 ± .019) | 0.310 | 0.310 | 0.280 | 0.283 | 0.293 | 0.157 | 0.173 | 0.153 |
| 0.3 | 0.613 (+.047 ± .019) | 0.610 | 0.603 | 0.570 | 0.567 | 0.600 | 0.303 | 0.400 | 0.263 |
| 0.4 | 0.903 (+.057 ± .014) | 0.897 (+.050 ± .014) | 0.890 | 0.857 | 0.847 | 0.873 (+.027 ± .010) | 0.540 | 0.693 | 0.480 |
| 0.5 | 0.987 | 0.987 | 0.987 | 0.977 | 0.980 | 0.977 | 0.753 | 0.917 | 0.713 |

* **Truncated split beats the full merge by 4--6 points at lam = 0.3--0.4** (about 3--4
  paired SEs), +3 at 0.2, nothing at the ceiling. Shallower is better: split@1 >=
  split@2 >= split@4 > split@8 > full split.
* ~~**It also beats the merge path truncated to the same depth** (split@2 vs
  mergecoarse@2: +2.4 points at lam = 0.4), so on this effect the divisive *direction*
  adds something beyond the multiplicity saving that DIVISIVE.md found at d = 8.~~
  **Withdrawn (2026-10-25):** `mergecoarse@2` is not a truncated merge search, so this
  comparison does not show what it was read as showing. Pending `merge@k` and the
  random-baseline runs (`split@k - rsplit@k`, `merge@k - rmerge@k`); the numbers are
  kept as recorded.
* **Adaptive searches beat every depth-0 statistic by 11--37 points.**
* **B barely matters for power**: 1--3 points from B = 100 to 1000, flat beyond.
* Caveat: truncated split runs ~1 point over nominal size (section above), worth
  ~1 point of this advantage. Size-adjusted, the gain is ~3--5 points.

### Education x Income, n = 1000 -- power, planted (non-monotone, Q1-blind) direction, 300 reps per lam

Rejection at B = 10000 (SE <= 0.020):

| lam | split@1 | split@2 | split@4 | split (14) | merge | mergecoarse@2 (historical) | max | euclid | mGCM |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.440 | 0.430 | 0.440 | 0.443 | 0.450 | 0.463 | 0.227 | 0.250 | 0.130 |
| 1.4 | 0.693 | 0.710 | 0.723 | 0.720 | 0.723 | **0.760** | 0.413 | 0.503 | 0.290 |
| 1.7 | 0.907 | 0.903 | 0.910 | 0.913 | 0.927 | **0.937** | 0.697 | 0.783 | 0.520 |
| 2.0 | 0.977 | 0.970 | 0.977 | 0.967 | 0.970 | 0.973 | 0.810 | 0.897 | 0.700 |

Paired vs merge: split@1 -0.010 / -0.030 / -0.020 / +0.007 (SE 0.009--0.021);
split@2 -0.020 / -0.013 / -0.023 / 0.000; mergecoarse@2 (historical, withdrawn)
+0.013 / **+0.037 ± 0.013** / +0.010 / +0.003.

* **The real-direction ranking reverses.** On this middle-vs-extremes contrast,
  truncated split is 1--3 points *below* the full merge (each ~1--1.5 SE). (The
  further claim that "the merge path truncated to the same depth is best" rested on
  `mergecoarse@2` and is withdrawn pending `merge@k`.)
* Likely mechanism: on an ordinal variable a divisive split cuts into contiguous
  pieces, so isolating the middle levels takes two splits from the top, whereas
  merging builds the middle group from below.
* Adaptive searches still beat every depth-0 statistic by 20--30 points.
* **Net so far:** truncated split wins on monotone effects (+4--6) and concedes a
  little on non-monotone ones (-1 to -3). Not uniformly best. ~~The merge path's
  coarse end (`mergecoarse@k`) is the most robust of the three across both
  directions on this pair.~~ **Withdrawn (2026-10-25):** that rested on
  `mergecoarse@k`, which is not a search; no conclusion about truncated merge or
  about the value of direction stands until the `merge@k` and random-baseline runs
  are in.

### Is the size excess a bug? Gaussian-vs-real check (`gauss_check.py`)

Education x Income, n = 1000, oracle, 1000 null reps, B = 3000 (SE 0.007). Each
replicate scores the real `T` and a fresh draw from `N(0, Sigma_hat)` against the
*same* bootstrap draws and tie-breaks.

| | real T | Gaussian T |
|---|---|---|
| split@1 / @2 / @3 / @4 | 0.060 / 0.059 / 0.061 / 0.061 | 0.043 / 0.044 / 0.047 / 0.048 |
| mergecoarse@1 / @2 / @3 / @4 (historical; `gauss_check.py` now reports `merge@k`, `rsplit@k`, `rmerge@k`) | 0.066 / 0.065 / 0.068 / 0.065 | 0.049 / 0.051 / 0.049 / 0.051 |
| merge | 0.059 | 0.047 |

* **The search and minP calibration are correct:** the Gaussian arm is nominal.
* **The ~1--1.5 point excess comes from the data:** the finite-`n` gap between `T` and
  `N(0, Sigma_hat)` at n = 1000 (rare Education levels the likely driver), shared by
  every method including the full merge. Expected to shrink with `n`; the n = 2000
  HoursPerWeek x Education size sweep that would test that was stopped.

## Not run (queues stopped at the user's request, 2026-09-30 ~22:45)

HoursPerWeek x Education real/planted power sweeps; the Occupation (Census tree)
sweep; the oracle + MLP learner arm at B = 10000; Occupation tree vs Saturated.
Launch commands are in the session's queue scripts and are reproduced by:

    python sweep_semisynth.py --x HoursPerWeek --y Education --n 2000 --direction real --lams 0.2 0.3 0.4 0.5 --reps 300 --time-reps 0
    python sweep_semisynth.py --x HoursPerWeek --y Education --n 2000 --direction planted --lams 0.55 0.7 0.85 1.0 --reps 300 --time-reps 0
    python sweep_semisynth.py --x Occupation --y Income --n 1000 --lams 0 --reps 1000 --time-reps 20
    python sweep_semisynth.py --x Occupation --y Income --n 1000 --direction real --lams 0.2 0.3 0.4 0.5 --reps 300 --time-reps 0
    python sweep_semisynth.py --x Occupation --y Income --n 1000 --direction planted --lams 0.6 0.9 1.25 1.55 --reps 300 --time-reps 0
    python run_semisynth.py --x Education --y Income --lams 0 --reps 500 --n-boot 10000 --learners oracle mlp --truncations 1 2 4

Also not run: the `merge@k` and random-baseline arms (`rsplit@k`, `rmerge@k`;
`RANDOM_SEARCH.md` Section 6). The sweep now includes them by default, so re-running
the Education x Income sweeps with the same `--lams`, `--reps`, `--seed` reproduces
the existing rows exactly and pairs the new ones with them:

    python sweep_semisynth.py --x Education --y Income --n 1000 --lams 0 --reps 1000 --n-boots 1000 --depths 1 2 4 8 --time-reps 0
    python sweep_semisynth.py --x Education --y Income --n 1000 --direction real --lams 0.2 0.3 0.4 0.5 --reps 300 --depths 1 2 4 8 --time-reps 0
    python sweep_semisynth.py --x Education --y Income --n 1000 --direction planted --lams 1.0 1.4 1.7 2.0 --reps 300 --depths 1 2 4 8 --time-reps 0
    python gauss_check.py --x Education --y Income --reps 1000 --n-boot 3000 --depths 1 2 4 8
    python run_semisynth.py --x Education --y Income --lams 0 --reps 500 --n-boot 10000 --learners oracle mlp --truncations 1 2 4 8

---

## Background findings (from the retired PLAN.md)

All measured at an exact null, oracle propensities, depth-0 criteria,
`normalise=False` unless stated. These explain the pre-minP size anomalies in
the pilot tables above.

### Finding 1 -- `rank(Sigma) = (dx-1)(dy-1)` is structural and harmless

`Sigma` is always rank-deficient (e.g. rank 3 of 8 at `dx=4, dy=2`, every
replicate). Reason: `sum_k (1{X=k} - f_k) = 0` and `sum_l (1{Y=l} - g_l) = 0`
hold in every row, so `prod_mat` has exactly `dx+dy-1` linear dependencies among
its columns. They hold for the row mean too, so the observed `T` satisfies the
same relations as the bootstrap draws. Measured leakage of `||T||^2` outside the
retained span, `1 - ||Q_r' T||^2 / ||T||^2`: **0.0000 in every cell of every
run**. `matrix_sqrt` discards exactly the structurally null directions, which `T`
does not occupy either. Leakage could only be nonzero for non-structural
degeneracy (an estimated propensity that collapsed a cell, a `Z`-stratum with no
variation); that case is unobserved. Worth reporting in the appendix per test
(one QR) alongside `rank(Sigma)` and `rank(Sigma)/n`, as a check that passes and
an alarm when it does not. A min-eigenvalue "necessity" claim would need a
construction with non-structural degeneracy; it does not exist yet.

### Finding 2 -- mGCM's studentisation over-rejects; excess scales with rank/n

Exact null, no rare cells, `dx = dy = d`, 800 reps, SE 0.008:

| d | n | rank/n | euclid | max | **mGCM** |
|---:|---:|---:|---:|---:|---:|
| 4 | 500 | 0.018 | 0.046 | 0.037 | 0.059 |
| 4 | 1000 | 0.009 | 0.066 | 0.059 | 0.059 |
| 4 | 5000 | 0.002 | 0.051 | 0.049 | 0.055 |
| 8 | 500 | 0.098 | 0.037 | 0.039 | **0.389** |
| 8 | 1000 | 0.049 | 0.054 | 0.045 | **0.215** |
| 8 | 5000 | 0.010 | 0.048 | 0.051 | 0.071 |
| 12 | 500 | 0.242 | 0.034 | 0.025 | **0.975** |
| 12 | 1000 | 0.121 | 0.040 | 0.039 | **0.790** |
| 12 | 5000 | 0.024 | 0.046 | 0.039 | 0.191 |

euclid and max are nominal everywhere. All three share `form_t_sigma`,
`matrix_sqrt`, `bootstrap_T` and the p-value code, so the fault is isolated to
the one thing mGCM does that they do not: divide by `sqrt(diag(Sigma))`.

Mechanism: `mgcm = max_j |T_j| / sqrt(Sigma_jj)`. Bootstrap draws have covariance
exactly `Sigma_hat`, so `T*_j / sqrt(Sigma_hat_jj)` is exactly standardised. The
observed `T_j` has variance `Sigma_true_jj != Sigma_hat_jj`, a ratio fluctuating
`O(1/sqrt(n))` per coordinate; the max over `p` coordinates then selects the one
where `Sigma_hat` most underestimated, a selection effect the draws cannot
reproduce. It grows with the number of coordinates per observation (`rank/n`).
Elimination argument: at depth 0 with a scalar criterion the p-value is exactly
uniform for any fixed criterion function if observed `T` and draws share a law,
so ApproxChi, Box, greedy search and the double bootstrap cannot be responsible;
only the law of `T` vs `N(0, Sigma_hat)` can. The repo's own `d=8, n=1000`
synthetic null (0.166) matches d=8, n=1000 here (0.215). A necessity claim is
available on a rate: size 0.055 -> 0.215 -> 0.975 as `(d-1)^2/n` goes
0.01 -> 0.05 -> 0.24.

Remedies for mGCM: drop the studentisation, give each draw its own variance
estimate (studentised bootstrap), or retire it as a comparator with the reason
documented. Any published power number involving mGCM is uninterpretable until
then.

Why the test suite missed it: `test_calibration_under_gaussian` feeds
`T ~ N(0, Sigma)` with `Sigma` known and fixed, so it tests the machinery given a
correct reference law and cannot see that the real `Sigma` is estimated from the
same rows. Regression tests must be fed from data.

### Finding 3 -- rare cells alone do not break size

Varying the rare level's expected count over {0.5, 1, 2, 5, 20, 100} at
`n = 1000`, `dx=4, dy=2`: euclid 0.043-0.060, max 0.037-0.061, mGCM 0.037-0.056,
leakage 0.0000, under both `normalise` settings. So the pooling effect on
Occupation x Income (0.095-0.128 unpooled -> 0.052-0.068 pooled) is likelier a
`rank/n` effect (`dx` 14 -> ~8) than a rare-cell one; pooling's real payoff is
controlling `rank/n`. Unconfirmed. Related unexplained cell: mGCM 0.060 at
n = 1000 -> 0.128 at n = 5000 on Occupation x Income contradicts "size improves
with n at fixed d"; `pool_rare` takes `n`, so `dx` may differ between the rows;
check before treating as a second phenomenon.

### The old double-bootstrap bug

Observed value ranked against `B` values excluding itself, bootstrap values
against `B` including themselves, so bootstrap quantiles were scaled by
`c = B/(B+1)` and the max over `M` depths gave `size = 1 - 0.95 c^M` (matched the
measurements to three decimals). Fixed by the minP rewrite; the arithmetic and
before/after table are already in `MINP.md`.

### Idea: uncentered multiplier bootstrap (not adopted; does not fix mGCM)

`T* = prod_mat' eps / sqrt(n)`, `eps` iid mean-zero unit-variance. Observed `T`
and every draw are the same linear functional of the same `n` rows with
different weights, so both lie in the row space of `prod_mat` by construction (verified to 1e-16 over
4000 draws), whereas Gaussian draws sit in a smaller subspace that a
deliberately degenerate `T` leaks out of by 1.5%. Do not centre: centred draws
lie in the span of centred rows, which reintroduces the mismatch. Uncentered has
conditional covariance `prod_mat' prod_mat / n`, right under the null and
slightly conservative under alternatives (`~ T T'/n` extra). Cost: thread
`prod_mat` through `TSigma`/`bootstrap_T` (`n*dx*dy` floats, 7.7 MB at
n = 30,162, dx*dy = 32). It does not fix mGCM, whose draws share the same
`Sigma_hat_jj`. Adopt only for robustness/skewness.

### Adult deep-dive design (inform the method; run once size is settled)

Why Adult (before Diabetes-130, which stays first for the paper because it has an
external outcome-blind partition, ICD-9 chapters in CCS-285 in phecodes): no data
engineering, an exact-null twin exists (`adult_semisynth.py`), and a near-oracle
nuisance arm is free (`Z = {Age, Sex}`: 14 strata, min occupancy 127 over 30,162
rows, so full-file stratum propensities proxy the truth on a subsample).

Change to the existing setup (`ADULT.md`): at n = 30,162 all 36 pairs reject
(largest p 1.1e-12), so rejection counts carry no information.
- **Subsample:** B independent draws of `n in {500, 1000, 2000}` without
  replacement; report a rejection rate over subsamples per pair.
- **Drop ill-conditioned tests, do not average over them:** 20 of 72 tests are
  ill-conditioned at n = 1000, all NativeCountry pairs. Compare on the 26
  well-conditioned pairs and report the excluded set.

Methods: `mGCM`/`euclid`/`max` (depth 0) -> `Tree.binary` -> `Saturated` (~10x
cost at dx = 13) -> `Ordinal` on ordinal pairs; against Ankan & Textor Q1/Q2/Q3
under both estimators (GLM, RFT) and both typings, plus `chi_sq` and
`multinomial`. Nuisance arm A: empirical-stratum propensities from the full file
(primary); arm B: `xgboost_learner` on the subsample (needs honest tuning; label
untuned).

Staging: A pilot, 6 Fig-1b pairs, B = 50, n = 1000, arm A, all methods (~30 min,
8 cores). B main, 26 pairs x B = 200 x n in {1000, 2000} (Saturated 7.5 s/rep at
`n_boot = 100`, Tree 0.8 s; ~22 CPU-hours, ~5.5 h on 8 cores). C learned
nuisances and partition analysis.

Record per test beyond the p-value: criterion path and selected partition; number
of groups at the criterion peak (the planted-2-group result above peaked at 7-9
groups, never 2, so if it reproduces the peak is not a usable selector and
"interpretable coarsening" needs a different stopping rule); merge-selection
frequency across subsamples (consensus partition with stability); seconds per
test; `rank` and `cond`.

Four questions: (1) does merging buy rejections over depth 0 (if mGCM and Tree
tie, the search is not justified)? (2) Tree vs Saturated: does 10x cost buy
rejections (if not, Tree.binary is default and divisive search loses urgency; if
yes, it becomes the priority since it truncates the fine end at d = 13-16)?
(3) where does catci beat A&T, and does the typing blind-spot story (ordinal Q1
df = 1 vs categorical Q2 df = 15) appear on real data (expect a loss to Q1 on
monotone pairs like Education x Income, a win on non-monotone)? (4) is the
criterion peak a usable partition selector?
