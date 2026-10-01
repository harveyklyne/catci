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
levels; `max`, `euclid`, `mGCM` depth-0; `chi_sq` pseudo-inverse chi-square;
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
| mergecoarse@1..8 | 0.054--0.072 | | | | |
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

| lam | split@1 | split@2 | split@4 | split (14) | merge | mergecoarse@2 | max | euclid | mGCM |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 0.317 (+.033 ± .019) | 0.310 | 0.310 | 0.280 | 0.283 | 0.293 | 0.157 | 0.173 | 0.153 |
| 0.3 | 0.613 (+.047 ± .019) | 0.610 | 0.603 | 0.570 | 0.567 | 0.600 | 0.303 | 0.400 | 0.263 |
| 0.4 | 0.903 (+.057 ± .014) | 0.897 (+.050 ± .014) | 0.890 | 0.857 | 0.847 | 0.873 (+.027 ± .010) | 0.540 | 0.693 | 0.480 |
| 0.5 | 0.987 | 0.987 | 0.987 | 0.977 | 0.980 | 0.977 | 0.753 | 0.917 | 0.713 |

* **Truncated split beats the full merge by 4--6 points at lam = 0.3--0.4** (about 3--4
  paired SEs), +3 at 0.2, nothing at the ceiling. Shallower is better: split@1 >=
  split@2 >= split@4 > split@8 > full split.
* **It also beats the merge path truncated to the same depth** (split@2 vs
  mergecoarse@2: +2.4 points at lam = 0.4), so on this effect the divisive *direction*
  adds something beyond the multiplicity saving that DIVISIVE.md found at d = 8.
* **Adaptive searches beat every depth-0 statistic by 11--37 points.**
* **B barely matters for power**: 1--3 points from B = 100 to 1000, flat beyond.
* Caveat: truncated split runs ~1 point over nominal size (section above), worth
  ~1 point of this advantage. Size-adjusted, the gain is ~3--5 points.

### Education x Income, n = 1000 -- power, planted (non-monotone, Q1-blind) direction, 300 reps per lam

Rejection at B = 10000 (SE <= 0.020):

| lam | split@1 | split@2 | split@4 | split (14) | merge | mergecoarse@2 | max | euclid | mGCM |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.440 | 0.430 | 0.440 | 0.443 | 0.450 | 0.463 | 0.227 | 0.250 | 0.130 |
| 1.4 | 0.693 | 0.710 | 0.723 | 0.720 | 0.723 | **0.760** | 0.413 | 0.503 | 0.290 |
| 1.7 | 0.907 | 0.903 | 0.910 | 0.913 | 0.927 | **0.937** | 0.697 | 0.783 | 0.520 |
| 2.0 | 0.977 | 0.970 | 0.977 | 0.967 | 0.970 | 0.973 | 0.810 | 0.897 | 0.700 |

Paired vs merge: split@1 -0.010 / -0.030 / -0.020 / +0.007 (SE 0.009--0.021);
split@2 -0.020 / -0.013 / -0.023 / 0.000; mergecoarse@2 +0.013 / **+0.037 ± 0.013** /
+0.010 / +0.003.

* **The real-direction ranking reverses.** On this middle-vs-extremes contrast,
  truncated split is 1--3 points *below* the full merge (each ~1--1.5 SE), and the
  merge path truncated to the same depth is best.
* Likely mechanism: on an ordinal variable a divisive split cuts into contiguous
  pieces, so isolating the middle levels takes two splits from the top, whereas
  merging builds the middle group from below.
* Adaptive searches still beat every depth-0 statistic by 20--30 points.
* **Net so far:** truncated split wins on monotone effects (+4--6) and concedes a
  little on non-monotone ones (-1 to -3). Not uniformly best; the merge path's
  coarse end (`mergecoarse@k`) is the most robust of the three across both
  directions on this pair.

### Is the size excess a bug? Gaussian-vs-real check (`gauss_check.py`)

Education x Income, n = 1000, oracle, 1000 null reps, B = 3000 (SE 0.007). Each
replicate scores the real `T` and a fresh draw from `N(0, Sigma_hat)` against the
*same* bootstrap draws and tie-breaks.

| | real T | Gaussian T |
|---|---|---|
| split@1 / @2 / @3 / @4 | 0.060 / 0.059 / 0.061 / 0.061 | 0.043 / 0.044 / 0.047 / 0.048 |
| mergecoarse@1 / @2 / @3 / @4 | 0.066 / 0.065 / 0.068 / 0.065 | 0.049 / 0.051 / 0.049 / 0.051 |
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
