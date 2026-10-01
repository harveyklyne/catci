# Paper vs code review

A review of the paper draft
(`Conditional_independence_testing_with_categorical_data/main.tex`, 1324 lines)
against master at `19d6d3f` and every worktree/branch, 2026-09-26. Refreshed
2026-09-30 against master `a2c23a2`, after TODOs 1, 2, 5 and 7 were merged. The
draft has no marked TODOs; everything below comes from reading the prose,
algorithms and proofs. Nothing here has been started.

Line numbers are `main.tex` lines. Lemma/theorem numbers use the shared counter
(`main.tex:73-74`), so "Lemma 10" is `lem:adaptivePi` (line 1030), "Lemma 5" is
`lem:greedyquerycontinuity` (876), and Theorems 1/2/3 are at 388/542/601.

- **C** items are **code TODOs**: ideas in the text not represented in the code.
- **T** items are **text TODOs**: places where the code is ahead of the text.

## Three things that matter most

1. **Cross-fitting (C1 / T7).** The estimator (lines 212-229), Assumption 1 (349,
   conditioning on `D^{(n,1)}`) and the Lemma 4 proof (749-874, "i.i.d.
   conditionally on `D^{(n,r)}`") all assume N-fold cross-fitting. The code fits
   propensities on the full sample and predicts in-sample, deliberately
   (`learners.py:8-11`, commit 970732f). Theorem 1 does not cover what the code
   does. Either restore an `nfolds` option, or replace Lemma 4 with Donsker /
   algorithmic-stability conditions.
2. **Calibration rewrite (T1).** Algorithm 2 and Theorem 2 still describe the old
   calibration. The only draft of the minP replacement is `MINP.md`, now tracked
   (it was nearly lost).
3. **Every figure and size/power claim predates the minP fix (T2).** The only
   minP-era numbers are in the study docs (SEARCH*.md, DIVISIVE.md, README,
   NOTES-mlp-learner.md); there is no size grid regenerating Fig 2.

## a) Code TODOs

### Not in TODO.md and not in any worktree

- **C1. No cross-fitting.** See above.
- **C2. The permuted-label figures (Figs 5-6, line 475) can't be regenerated.**
  `dgp.permute_labels` exists, but `experiments/run.py:42` hardcodes
  `permute=False`. There is no config field or CLI flag.
- **C3. No script draws any figure or runs the grid of settings.** The paper uses
  `size_plot.pdf` and four `power_plot_*.pdf` (lines 432-487). `run.py` runs one
  config at a time, and `experiments/results/` holds only Adult runs.
- **C4. The power DGP's λ is not the paper's λ.**
  - Paper line 443: `p_j p_k + λ d_{jk}`.
  - Code `dgp.py:253`: `indep + (strength * indep.min()/2) * interaction`.
  - The effective λ depends on the data in each replicate. The x-axis of Figs 3-6
    is really `strength ∈ {0.2,…,1.8}` (`config.py:80`). Fix one side.
- **C5. The `c^sig` and `c^sin` breakpoints differ.**
  - Paper lines 415-418: `(−2,−1,0,0.5,1,2,3)` and `(−2.3,−1.2,−0.6,0,0.6,1.2,2.3)`.
  - Code `dgp.py:88`: `center + 2Φ⁻¹(k/8)`, giving sig ≈ `(−1.80,−0.85,−0.14,0.5,1.14,1.85,2.80)`
    and sin ≈ `(−2.30,−1.35,−0.64,0,0.64,1.35,2.30)`.
  - The code's are deliberate ("fix 2c"). lin/vee/hat, `d^tree` and `d^step` match
    exactly.
- **C6. `multinomial` is in the Fig 2 caption (line 433) but not in the size
  config** (`config.py:59`, `run_size.py:29`).
- **C7. `catci_test` defaults to `normalise=True`** (`api.py:43`).
  - The paper's remark at lines 236-238 says normalising hurts, every experiment
    config uses False (`config.py:40`), and PLAN.md finds True worse (greedy
    0.095 → 0.128).
  - The theory does not cover normalising. Rescaling by `D^{-1/2}` moves the
    kernel from `𝒦` to `D^{1/2}𝒦`: Lemma 6 (917) survives, but Lemma 7's
    construction of `y` (990-1003) breaks. Nothing guards against `Σ_jj = 0`.
  - No config runs a normalised tree or ordinal test, so "see Figure size" is
    backed only by mGCM.
- **C8. Search limits are not implemented.** Line 255 mentions "fixing a maximum
  group size or a total number of merges", and Algorithm 4 takes a "search depth
  `L`" (626). The merge search always runs to 2×2. Only the divisive search can
  be truncated (`SplitSearch(max_levels)`). DIVISIVE.md finds truncation's gain is
  entirely the smaller minP multiplicity penalty, and suggests testing a weighted
  minP over the ordinary merge path instead.
- **C9. The abstract's Bonferroni baseline is not the one implemented.** Line 109
  means "testing at every possible granularity" with Bonferroni, i.e. a fixed set
  of granularities such as each tree depth. `*_bonf` (`calibrate.py:101`) applies
  Bonferroni over the *adaptive* path.
- **C10. `bonferroni_pvalue` does not share tie-breaks with minP**
  (`calibrate.py:101-119`). The docstring, and MINP.md's pathwise-dominance proof,
  say minP is never larger on the same input. The proof needs both to use the same
  stage-1 ranks, but each call draws fresh uniforms from `rng`. Either share the
  ranks or weaken the claim; check what `tests/test_calibrate.py` actually asserts.
- **C11. Weighted GCM `E[Cov(X̃_j, Ỹ_k | Z) w(Z)]`** (line 210, Scheidegger): no
  `w(Z)` option in `gcm.py`.
- **C12. No regularity diagnostics.**
  - Nothing warns when `rank(Σ̂) ≠ (dx−1)(dy−1)` (Assumption 2, line 381).
  - `bootstrap.py:25` drops eigenvalues silently.
  - A zero-trace merge scores −inf (`search.py:521`), which Lemma 6 says cannot
    happen under Assumption 2. An all-NaN draw fails with a confusing "different
    lengths" error (`search.py:525`).
  - The Adult work found p-values being returned at cond ≈ 1e21.
  - Cheap fixes: a rank warning, reporting `‖T‖²` leakage outside `range(Σ̂)` and
    `rank/n` (PLAN.md "leakage diagnostic"), and a clear error for zero-trace
    merges.
- **C13. Chi-square degrees of freedom.** Paper eq. (6), line 262: `χ²_{rank(σ)}`.
  Code `methods.py:140`: fixed `(dx−1)(dy−1)`. They agree only at generic rank.
- **C14. Resolved on master.** Oracle propensities are now `run.py --learner
  oracle`. (Was: oracle experiments existed only as an API path.)
- **C15. The multiplier bootstrap** (PLAN.md: uncentred `T* = prod_matᵀε/√n`) is
  unimplemented. It would put the observed `T` and the draws in the same row space
  by construction.
- **C16. Draft asides with nothing behind them in the code:**
  - `pcalg` structure learning on real data (line 101);
  - a formal uniformity test for the size simulation (102);
  - high-dimensional results for specific `q` (103);
  - a hierarchical-classifier learner (115). The MLP on master is a plain
    multinomial classifier, not hierarchical;
  - the Li–Shepherd parametric tests, Liu's latent-variable approach, and
    stratification for discrete `Z` (129-153). Only Ankan & Textor is implemented.

### Already tracked, with current status

| Text | TODO | Status |
|---|---|---|
All of these are now on master; the open part is mostly the write-up (T18).

| Text | TODO | Status on master |
|---|---|---|
| Sample splitting (line 138) | 7b | `<search>_split` methods; loses 8-35 points against greedy (SEARCH_COMPARISON.md) |
| Random search / sanity check | 4.1, 7a | `<search>_random` methods; greedy beats random by 7-25 points where the structure leaves the final partition free, 1-2 under a binary tree |
| Wider search | 7c | `greedy_beam5`; no gain even on the greedy trap (SEARCH.md) |
| Motivating example `d = 30 × 10` (131); "fix d = 8" (405) | 1 | CV tuner (`tuning.py`, `experiments/tune.py`), separate dx/dy, d-axis pilot up to dx = 256 (README). Learners are worse than uniform at `d ≥ 64`, so large d uses oracle propensities |
| Pre-tuning (615-617) | 1 | CV tuner on master; see T16 |
| Deep / hierarchical classifier (115) | 2 | MLP learner merged and now the **default**; every experiment runs both MLP and xgboost (NOTES-mlp-learner.md) |
| Temporal / taxonomy structures (112) | 3 | `Cyclic`, `Tree.from_parents` |
| Divisive search | 5 | `SplitSearch` with truncation; d = 8 and d = 16 answered, union of paths gives no gain (DIVISIVE.md) |
| Exact CDF criterion | 6 | `ExactChi`; paired power study never finished |
| `n_boot` default of 100 vs 1000 in the configs | 0b / 8 | Still open: `api.py:42` defaults to 100 |

## b) Text TODOs

### Calibration

- **T1. Replace Algorithm 2 (306-322), Algorithm 5 (693-713), Theorem 2 (542-548)
  and its proof (1159-1219, Lemmas 13-14 at 1221-1320) with minP.** Draft in
  `MINP.md`.
  - Code (`calibrate.py:44-98`): pool `B+1` exchangeable paths, leave-one-out ranks
    with random tie-breaks, `P_min` per path, the same rule again at stage 2.
    Exactly uniform on the `1/(B+1)` grid at finite `B`.
  - Old Algorithm 2 ranks the observed value among `B` draws but each draw among a
    pool containing itself, and can return `p = 0`. Algorithm 5 has the same
    asymmetry, plus a typo (`U` drawn as `U_l` but used as `U`). Line 304's
    continuous ECDF tie handling is superseded by random tie-breaking.
  - Theorem 2 needs `n, B → ∞` and a DKW bound. A simpler proof: exact validity for
    the Gaussian problem at fixed `B`, then joint convergence as `n → ∞` with `B`
    fixed (Theorems 3 and 4), with the p-value map a.s. continuous off ties. That
    drops Lemma 14 and the DKW bound; Lemma 13 then goes unused.
  - Also add: the Bonferroni floor `L/(B+1)` and `B ≥ L/α` (lines 193 and 301 only
    call Bonferroni "overly conservative"); that `B` bounds power; and a footnote
    on random ties (conservative ties give size 0 at `L ≥ 13, B = 100`).
  - Theorem 1 (388) says `lim_{n,B→∞}`; with minP, `B` need not grow.
  - Notation: `M → S`, `F → p`, `G → P_min`. The five citations in MINP.md are from
    memory and need checking.
- **T2. All size and power claims need regenerating.** "Correctly calibrated in
  all settings" (line 437) and "controls size in finite-sample settings" (494)
  were produced with the inflated calibration (tree 0.080, ordinal 0.072 at
  `d = 8`). Figs 2-6 and their prose must be redone.

### Estimator and theory

- **T3. √n scaling.** Line 222 defines `T̂ = (1/n)Σ…`, but lines 193 and 335 draw
  `T ~ N(0, Σ̂)`, where `Σ̂` (226) is the per-observation covariance. That is off
  by a factor of `n`, and because φ is not scale-equivariant (`g` scales, `h`
  does not) the algorithm is wrong as written. Code `gcm.py:76-77`:
  `√n·mean`, `np.cov(ddof=1)`. The code is right.
- **T4. Index order.** The paper (209, 185-186) is Y-fastest; the code
  (`gcm.py:52`, `dgp.py:20`, `merging.py:8`) is X-fastest.
- **T5. "Rows and columns of each d sum to one"** (454) should say *zero*
  (`dgp.py:141`).
- **T6. mGCM.** Replace "we were surprised that the mGCM also failed" (437) with
  the diagnosis (TODO 0a, PLAN Finding 2): oracle-propensity size 0.215 at
  `d = 8, n = 1000` and 0.975 at `d = 12, n = 500`, scaling with `rank(Σ)/n`.
  Studentising by `diag(Σ̂)` lets the max pick the coordinate whose variance was
  most underestimated, which the draws cannot reproduce. Consider a remark that
  `rank(Σ)/n`, not a minimum-eigenvalue bound, is the regularity condition that
  actually matters. Assumption 2's "fixed Σ" is already strong.
- **T7. Cross-fitting** (see C1). If full-sample fitting stays, rewrite lines
  212-229, Assumption 1 and Lemma 4; otherwise the rationale in `learners.py:8-11`
  has no theoretical backing.
- **T8. Normalised statistic** (see C7). If `normalise` stays, add it to
  Algorithm 3 (329-341) and extend Lemmas 6-7.
- **T9. `L` must not depend on `(t, σ)`.** Line 185 and Algorithm 1 leave it
  implicit. The code requires every draw to have the same path length
  (`search.py:524`), and the theorem implicitly needs it.
- **T10. Other search maps.**
  - Lemma 10 is stated for any map `t ↦ 𝒞_*`, so it covers beam and divisive
    search.
  - Lemma 5 (a.s. continuity, 876-913) is proved only for the greedy argmax, and
    only in `t` for fixed `Σ`, although Assumption 4 is joint in `(t, σ)`. Extend
    it for beam, divisive and truncated searches, and close the σ-direction gap.
  - Random search: condition on the external randomness.
  - Sample splitting: needs a conditional version of Theorem 3.

### Search, structures and implementation

- **T11. The fast greedy merging appendix (619-690, Algorithm 4) describes the old
  scheme.** It updates `T` and `Σ` explicitly and recomputes each candidate's
  trace increments per level, one draw at a time. The code (`search.py:15-36,
  196-500`):
  - scores every candidate at once as bilinear forms;
  - carries `cross`, `term1`, `term2` and `d_tr` across levels (`_other_axis_dense`,
    427-452, has no paper counterpart);
  - uses dense mode (one O(p²) GEMM per merge) for Saturated and sparse mode
    (O(pd) per level) for trees and ordinal;
  - never shrinks arrays (zeroed slots plus compaction), batches all `B+1` draws in
    lockstep, and has `n_jobs` threading.

  Also, the equation labels are misattached: `eqn:fastnormsq` sits on the trace
  line (670) and `eqn:fasttr` on the first `tr²` line (672), so "via
  (fastnormsq)" at 679 cites the trace formula.
- **T12. Bootstrap generation (Algorithm 6, 717-731)** says "chol in base R". The
  code uses an eigendecomposition with relative tolerance `p·eps·λ_max`
  (`bootstrap.py:16-27`). The law is the same.
- **T13. `ExactChi`** (`statistic.py:121-270`, `statistic="exact"`) has no
  counterpart in the paper. It is the exact weighted-χ² CDF of `‖t‖²`
  (Ruben/Imhof), not the pseudo-inverse eq. (6) the paper dismisses (263). About
  30× slower, same partitions, p-values match in 17 of 18 cases. The appendix
  claim that knowing `‖t‖², tr σ, tr σ²` suffices (622) is Box-specific. The
  proof of Lemma 5 (888-890) needs one line that the exact CDF is continuous in σ.
- **T14. Structures.**
  - Line 255 says tree merges are "between siblings". The code (`structure.py:
    205-375`) allows n-ary trees with a context rule: a partial union of children
    keeps absorbing its siblings. Pure sibling merges would strand `{a, b}` next
    to `c`.
  - Undescribed: `Tree.from_parents` (code → parent mapping, `structure.py:250`),
    `Cyclic` (ordinal plus the `(1, d)` wrap, `structure.py:111-134`), and
    `Saturated` (unrestricted "greedy", in the method registry but not in
    Section 3).
  - The theory needs no binary-tree or ordinal assumption: Lemma 5 needs only a
    finite candidate set in `𝒞_*` per stage (guard at `structure.py:90-92`).
    Say so in one sentence.
- **T15. "R package catci"** (lines 142, 402, 494) should now say Python (ff19ca5).
- **T16. Learners and tuning.** The text (217, 402) names xgboost as the learner.
  On master the **MLP is the default** in `catci_test` and every experiment runs
  both learners, paired on the same data. Hyperparameters now come from a K-fold
  CV tuner (`tuning.py`, `experiments/tune.py`), not frozen R tunings; line 615's
  "pre-tune … on 1000 datasets" should describe it. The text should also say that
  propensities beat uniform only up to moderate `d`.

### Experiments section

- **T17. Unstated settings.** "All X and Y settings under consideration" (456,
  475): the legacy runs were size hat_hat, lin_hat, lin_lin, lin_vee, sig_sig,
  sin_sig, sin_sin, vee_hat and vee_vee; power (tree) lin_lin, lin_vee, sin_sig and
  sin_sin; power (step) lin_lin and sin_sin. The strength grid 0.2-1.8 and
  `n_boot = 1000` are not stated. The power configs also compute mGCM, chi_sq and
  multinomial, which the figures don't show.
- **T18. Results on master ready to become paper sentences** (minP-era unless
  noted; the source docs have the current numbers):
  - Beam does not help, even on a setting built to trap greedy, so "greedy
    suffices because the calibration charges for search breadth" (SEARCH.md).
  - Greedy beats random search by 7-25 points where the structure leaves the
    final partition free; under a binary tree the structure, not the search,
    supplies the adaptivity. Sample splitting loses 8-35 points, often below
    plain χ² (SEARCH_COMPARISON.md, which flags its numbers as exploratory, with
    a checklist before they go in the paper).
  - Divisive search: the full split matches merging. Truncated splitting helps
    only when the signal is coarse, and that gain is entirely multiplicity:
    calibrating over 3 levels instead of 13. Union of merge and split paths gives
    no gain (DIVISIVE.md). Suggests a weighted-minP remark or experiment.
  - d-axis pilot (oracle propensities, up to dx = 256): ordinal power decays
    gently while euclid and max collapse and chi_sq/mGCM size blows up; cost is
    roughly cubic in dx (README).
  - MLP vs xgboost: learner choice doesn't affect size (NOTES-mlp-learner.md).
  - Adult and the semi-synthetic study: 5 of 6 Fig 1b pairs flip to rejection; a
    planted direction A&T cannot see (A&T 0.060 vs catci 1.000). **Pre-minP, still
    on the unmerged `adult-application` / `worktree-adult-semisynth` branches;
    needs re-running.**

### Abstract and discussion

- **T19.** Remove the draft material: "From Gemini – do not use verbatim"
  (111-112) and the CMStatistics text (117).
- **T20.** "Optimally aggregate levels" (117) should be "greedily".
- **T21.** "More powerful than … existing methods" (494) is only backed by euclid,
  max and Ankan & Textor. The discussion also says "calibrated using a single
  parametric bootstrap sample", which is still true under minP but worth
  rephrasing.
- **T22.** Add the multi-label limitation (TODO 9) as one line.

## Housekeeping

- TODO.md (untracked) still lists items 1, 2, 3, 5 and 7 without a DONE marker,
  and item 6 points to the deleted `exact-pvalue-criteria` worktree.
- Stale branches: `refactor` (fully merged); remote copies of the now-merged
  research branches (`origin/worktree-{alt-search,divisive-search,
  divisive-search-vec,mlp-learner,structures,tuner-d-axis}`).
- Still unmerged: `adult-application` and `worktree-adult-semisynth` (TODO 4a,
  pre-minP, with an uncommitted `experiments/adult_semisynth.py`).
