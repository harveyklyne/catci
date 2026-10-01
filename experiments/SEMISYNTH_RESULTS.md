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
