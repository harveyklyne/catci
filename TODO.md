# TODO

Ordered so that the things gating everything else come first. The numbered items are
unchanged so their cross-references still hold.

**Status (2026-09-30).** Items 1, 2, 3, 4.1, 5, 6, 7 and 8 are DONE and merged to
master (`51f4e09`).

**Decisions (Harvey, 2026-09-30):**

* **Real data before simulation.** The real-data studies decide which methods are
  useful in practice; the final simulation grid (0b, item 2) waits behind them.
* **The primary method is the truncated divisive search** (`SplitSearch(max_levels=k)`).
  The merge search (Algorithm 1) is secondary.
* **`normalise=False` for all our methods.** Only the mGCM competitor studentises,
  and it does so internally (`statistic.mgcm`), so it too is fed the unnormalised `T`.
  `catci_test` now defaults to `False`.

**Priority order — real data first:**

1. **Integrate the Adult work** (branch `worktree-adult-integration`): Ankan & Textor
   port merged; `adult_semisynth.py` committed with its identities pinned
   (`tests/test_adult_semisynth.py`); `run_semisynth.py` runs every method paired, for
   oracle and fitted propensities. Remaining: the minP null re-check on three pairs
   (running), then merge to master.
2. **Oracle-only search study on the Adult semi-synthetic data.** Which search gives
   the best power with the exact propensities: merge, full split, split@k for several
   `k`, against the depth-0 comparators and the competitors. Run over several
   effect directions (real monotone, planted non-monotone, and others) and pairs, so
   no method is favoured by construction. This picks `k` and tests the primary method.
3. **Add the fitted-learner arm** (MLP, then xgb) beside the oracle on the same data
   and bootstrap draws, so the cost of estimating the nuisances separates from the
   effect of the search. Then the rich-`Z` size arm. Needs `tune.py` at the Adult
   `(n, num_class)`.
4. **Calibration checks on real Adult** (4a.1): the education/education-num bijection
   and sex-given-relationship separation.
5. **Real Adult analysis** (4a.2–4a.3), scored against the literature's established
   dependencies (`related-literature/ci_ground_truth_2026-08.md` §2): Fig. 1b
   head-to-head, then the `fairadapt` DAG.
6. **Diabetes-130** (4b), then **ACS** (4c).
7. Then the simulation work: 0b's final grid, 0a's write-up.

**Blocker for the primary method on unordered variables.** `SplitSearch` supports
`Ordinal` and `Tree` only; `Saturated` raises `NotImplementedError` (`2^(s-1) - 1`
bipartitions per group). So truncated divisive search needs a tree over the levels of
any unordered categorical (Occupation, Workclass, Race, ...). Options: a published
taxonomy where one exists (Census occupation major groups; ICD-9 chapters for 4b),
or an exhaustive top-level split for small `d` — at `d = 14` the first level is 8191
candidates, which is feasible. Decide before step 2 covers categorical pairs.

Other open items:

* **0a** — resolved in principle (see below); needs writing up in the paper.
* **0b** — the paper's Algorithm 2 rewrite (`main.tex` has no minP yet; draft in
  `MINP.md`, now on master), and the final simulation grid (0b, item 2).
* **Paper: the divisive search's sub-linear cost in `B`** (found 2026-09-30, write-up
  in `experiments/DIVISIVE_COST.md`). The costly trace terms depend only on the
  partition and `Sigma`, so each partition is scored once however many bootstrap
  draws reach it, and at the coarse end the draws concentrate on few partitions.
  Truncated divisive search therefore costs a near-fixed amount plus a cheap
  per-draw update: split@2 takes 0.10 s at `B = 1000` and 0.27 s at `B = 10,000`,
  against the merge search's 0.25 s and 2.7 s. That makes the larger `B` the minP
  calibration rewards almost free. Add a paragraph to the appendix
  "Implementation" section beside "Fast greedy merging" (`sect:fastgreedymerging`),
  as part of writing up the divisive search, which `main.tex` does not describe yet.
* Small leftovers: `catci_test` default `n_boot = 100` → 1000; the divisive `alt`
  union cell; the numba merge kernel (only if `p` reaches the thousands); item 9.

## 0. Two calibration defects found 2026-09-02

Both were found while scoping the application study; see `PLAN.md` Priority 0 for
the full evidence and the scripts. They are independent of each other. 0b is fixed
(2026-09-20); 0a is resolved — it is a finding, not a bug. Results produced since
the 0b fix (the per-learner `d = 8` size/power runs, the `d`-axis pilot, the
search and divisive studies) are on minP code; anything older is not.

### 0a. Work out what `mGCM` is actually supposed to be — RESOLVED, write up

**Harvey: mGCM is the same as our 'max' (max norm with no merging) when the test statistic T is normalized using the inverse estimated covatiance matrix. This is exactly as intended - the method is from Shah and Peters. It is miscalibrated - this is good for us. We show you need to be careful with assuming CLT is in effect with categorical data.**

* **As currently written it is badly miscalibrated**, and only in the
  `normalise=False` setting where it is distinct: exact null, oracle propensities,
  size 0.215 at `d = 8, n = 1000` and **0.975** at `d = 12, n = 500`, scaling with
  `rank(Sigma)/n`. Studentising per coordinate promotes the worst-estimated
  (lowest-variance) coordinates into contention, then the max selects on their
  estimation error — something the draws cannot reproduce, since their variance is
  `Sigma_hat` by construction. `euclid` and `max` are clean everywhere.
* **More evidence from the `d`-axis pilot** (README, "Does the merging advantage
  widen with `d`?"): at `n = 2000, dy = 4`, oracle propensities, mGCM size climbs
  0.055 → 0.320 → 0.640 for `dx = 8, 64, 128`, and `chi_sq` reaches 0.405 → 1.000;
  `ordinal`, `euclid` and `max` stay at or below nominal throughout.
* **Open: the paper paragraph.** CLT-based normalisation is unsafe for categorical
  data as `rank(Sigma)/n` grows; the bootstrap-on-`Sigma_hat` calibration is not.
  The final grid (0b, item 2) supplies the table.

### 0b. Fix the bootstrap p-value — DONE (2026-09-20), rewritten as minP

Landed. `calibrate.py` is now a minP test: pool the observed path with the `B`
bootstrap paths into `B + 1` exchangeable paths, rank every one against the *other*
`B` by one identical leave-one-out rule (stage 1 = Beran prepivoting), take each
path's min over levels, and calibrate the observed min by the same rule (stage 2).
Exactly uniform on the `1/(B+1)` grid at finite `B`, for any `L`, any dependence
between levels, and any deterministic search map. Size at `alpha = 0.05`, `B = 100`:

| L | 1 | 2 | 5 | 10 | 20 | 40 |
|---|---:|---:|---:|---:|---:|---:|
| before | 0.049 | 0.064 | 0.087 | 0.129 | 0.203 | 0.336 |
| after | 0.052 | 0.047 | 0.050 | 0.049 | 0.049 | 0.049 |

Confirmed end-to-end on the PLAN.md Finding 4 setting (exact null, oracle
propensities, `d = 8`, `n = 1000`, `n_boot = 100`, 1000 reps, SE 0.0069):

| method | size@0.05 | size@0.10 | before |
|---|---:|---:|---:|
| tree | 0.051 | 0.092 | 0.080 |
| ordinal | 0.048 | 0.096 | 0.072 |
| greedy | 0.051 | 0.103 | — |
| `*_bonf` | 0.000 | 0.000 | — |

All three searches are inside 1 SE of nominal. The `_bonf` zeros are the expected
degenerate regime: floor `13/101 = 0.129 > 0.05`. At the new `n_boot = 1000`
default the floor is `13/1001 = 0.013` and the comparison becomes meaningful.

Pinned by the new `tests/test_calibrate.py`, parametrised by `L`, which reaches
`double_bootstrap_pvalue` directly with no DGP or search. `test_properties.py`
tolerances tightened from 0.05 to 0.03 and `test_tree_calibrates_too` now asserts
rejection rates, not just the mean.

Two things the fix surfaced:

* **Ties must be broken at random.** Each level puts one path at the floor
  `1/(B+1)`, so up to `L` paths tie there. Conservative tie handling gives size
  **0.000** at `L >= 13, B = 100` — it cannot reject at all.
* **`n_boot = 100` was costing real power**, ~half of what is available at
  `n_boot = 1000` for `d = 8` (`L = 13`). The experiment configs now use 1000.
  `catci_test` still defaults to 100 — worth revisiting.

`bonferroni_pvalue` added as the simple-FWER comparator, exposed as the `*_bonf`
entries in the method registry. Its floor is `L/(B+1)`, so it needs
`n_boot >= L/alpha` (~260 at `d = 8`) before the comparison is non-degenerate —
itself the cleanest argument for minP over simple FWER control.

**Still open, both for Harvey:**

1. **Rewrite Algorithm 2 in the paper.** Draft, proof sketch, tables and the
   citations to verify are in `MINP.md`. Notation: `M -> S` (test statistics),
   `F -> p`, `G -> P_min`. `MINP.md` lives on branch `worktree-track-minp`
   (with `TEXT.md`, a paper-vs-code review) — merge that first.
2. **Run the final simulation grid.** The post-minP numbers that exist were run
   piecemeal with different `reps`/`n_boot` (the `d`-axis pilot: 100 reps,
   `n_boot = 200`; `SEARCH_COMPARISON.md` flags its own cells as inconsistent).
   Re-run everything the paper will cite, once, on one commit, with
   `n_boot = 1000`, so every table shares one provenance. Proposed blocks:

   | block | what it shows | setting | methods | reps |
   |---|---|---|---|---|
   | **A. Size** | every method holds level; mGCM / `chi_sq` do not | `size_config` for margins `lin_lin`, `sin_sin` (+ one cross pair, e.g. `lin_sin`), `n ∈ {500, 1000, 2000}`, `d = 8`; learners `mlp`, `xgb`, `oracle` | `tree`, `ordinal`, `greedy`, `*_bonf`, `*_exact`, `max`, `euclid`, `mGCM`, `chi_sq`, `ankan`, `multinomial` | 1000 (SE 0.007) |
   | **B. Power, d = 8** | adaptive beats depth-0 and competitors | `power_config` `lin_lin_step` and `sin_sin_binary_tree`, `n = 1000`, strengths 0.2–1.8; `mlp` + `xgb` | as A, plus the **misspecified** cross (`tree` on `step`, `ordinal` on `binary_tree`) — the paper's "not much worse when the structure is wrong" claim currently has no number behind it | 200 |
   | **C. d axis** | merging advantage widens with `d`; mGCM / `chi_sq` break | `d_grid` `lin_lin_step`, `n = 2000`, `dy = 4`, `dx ∈ {8,…,256}`, strengths `{0, 0.6, 1.2}`; `oracle` + `mlp` (needs `tune.py` at `n = 2000` for each `dx`) | `ordinal`, `euclid`, `max`, `mGCM`, `chi_sq`, `ankan` | 200 |
   | **D. Exact vs approx** | "same answers, ~1/28 the time" | B's two DGPs at strengths `{0, 0.8, 1.4}`, paired | `*` vs `*_exact`, plus timings | 200 |
   | **E. Search alternatives** (Gaussian limit) | greedy vs random / split / beam; divisive truncation | `search_study.py` and `power_divisive.py` cells as now, with uniform reps and `n_boot` | as in `SEARCH.md`, `DIVISIVE.md` | 500, `n_boot = 1000` |

   D needs no extra data: add `*_exact` to B's method list, so the comparison is
   paired on the same replicates and draws. C is the
   expensive block (`p = dx·dy` up to 1024). Then re-check the
   `ADULT_SEMISYNTH.md` cells once 4a is merged.

## 1. Port a tuner, and make `d` a real axis — DONE (merged `51f4e09`)

`experiments/tune.py` tunes both learners by CV (`--protocol holdout|cv`, writes
need `--write`); `config.py` has separate `dx`/`dy` axes and `d_grid`;
`run_d_axis.py` runs it. **Answer: yes, the advantage widens.** Pilot (`lin_lin_step`,
`n = 2000`, `dy = 4`, oracle, README): at strength 0.6, `ordinal` falls only
0.64 → 0.37 from `dx = 8` to `256` while `euclid` falls 0.45 → 0.02 and `max`
0.22 → 0.02. Only `n = 1000, d = 8` is tuned so far; the final run is 0b item 2,
block C.

Original text:

Nothing below can run at a realistic size until this is done. `experiments/tuning/`
holds frozen hyperparameters for `n = 1000, d = 8` only; the R script that produced
them was cluster-specific and was not ported, so a new `(n, d)` cannot be tuned at
all. Every application in item 4 wants `d` in the hundreds.

Coupled to it: `dgp.py` is parametric in `d` but every config is hardcoded to
`d = 8`, while the motivating example is `dX = 30, dY = 10`. If the applications
run at `d` in the hundreds, the simulations need to cover that range too, or
there is a visible gap between what the paper demonstrates and what it applies.

  1. Port / rewrite the tuner so a new `(n, d)` can be tuned by CV.
  1. Add `d` to the config grid and re-run. What does the power curve do as `d`
     grows with `n` fixed — does the merging advantage widen, as claimed?

## 2. MLP classifier instead of xgboost — DONE (merged `51f4e09`)

The MLP is the default learner (decision 2026-09-26); every experiment runs `mlp`
and `xgb` on identical data, with `oracle` as a control. The MLP has 27–72% lower
propensity error on 4 of 5 settings, but both learners clear Assumption 1 by a wide
margin, so the learner is not what limits the test. `NOTES-mlp-learner.md` still
recommends xgboost as the default, which is stale.

## 3. Structures the applications need — DONE (merged `19d6d3f`)

`Tree.from_parents` (code→parent mapping, n-ary) and `Cyclic` are in
`structure.py`, with the n-ary mask fix.

Both are small, and both are prerequisites for item 4 rather than research.

  1. **Taxonomy → `TreeNode` constructor.** `TreeNode` is already n-ary
     (`_make_parent` takes any number of children, `_get_siblings` does not
     assume arity) — only the `Tree.binary` *constructor* is binary. Real
     taxonomies (ICD-9 chapters ⊂ CCS-285, SOC, NAICS) need a constructor from a
     code→parent mapping, not new search code.
  1. **`Cyclic` structure** for temporal categoricals — day-of-week,
     month-of-year, hour-of-day. `Ordinal` assumes a line; these wrap. Should be
     adjacent pairs plus `(1, d)`, and the invariant looks like it survives
     merging because `search.py` lands the merged group at the lower position, so
     an arc containing the wrap point stays at position 1. Verify that with a
     property test rather than assuming it. Low power ceiling (`d` = 7, 12, 24),
     so this is structure coverage, not an application.

## 4. Apply to real-world datasets. Do we actually see the method looks powerful?

**Status (2026-09-30): OPEN — the main remaining work.** The prerequisites (items 1,
2 and 3) are all done. Our method has not been applied to any dataset yet. The
Adult work so far (the Ankan & Textor port, `ADULT.md`, and the semi-synthetic
design, `ADULT_SEMISYNTH.md`) sits on the unmerged branches `adult-application` and
`worktree-adult-semisynth`, so merging those is the first step for 4a.

Reordered. The old ranking had ACS second on the grounds that it is "a more
modern, high-dimensional version of UCI Adult"; the deciding criterion is really
whether a dataset supplies an *external, published, outcome-blind* partition to
report the learned one against, and on that criterion Diabetes-130 beats ACS.
See `related-literature/ci_ground_truth_2026-08.md` §6.1.

Cross-cutting, applies to all three:

  1. Decide how to clean / pre-process the data. Spot check the existing methods
     on each dataset and identify which of their rejections / non-rejections are
     suspicious.
  1. Apply my method to the data, using appropriately tuned MLP and XGBoost
     models (CV). Do I reject a superset of the ones that will have good size
     control?

### 4a. UCI Adult — same as Ankan and Textor

`sklearn.datasets.fetch_openml("adult")`. Note that what exists so far is *their*
method, not mine: `ADULT.md` is scoped to "their method's p-values only", and the
semi-synthetic study in `ADULT_SEMISYNTH.md` has its design settled and the
separation result established but has never been run. Three deliverables of quite
different value, worth keeping separate:

  1. **Tier-1 calibration.** `education` vs `education-num` (an exact bijection —
     must abstain or flag determinism, not silently reject) and `sex` given
     `relationship` (perfect separation on 45% of rows — must not blow up).
     Cheap credibility; nobody does it.
  1. **Head-to-head on Fig. 1b.** Six hypotheses at `n = 1000` with their
     discretisation and variable set, our verdicts beside theirs.
  1. **The `fairadapt` DAG.** All six testable implications fail; three have never
     been tested. This is a *finding*, not an illustration, and it is the
     strongest single result in the repo — it should not stay buried inside
     "same as Ankan and Textor".
  1. Run the semi-synthetic study to completion.

### 4b. Diabetes 130-US Hospitals, 1999-2008 (UCI) — was 3rd, now 2nd

The only candidate where the partition the method *reports* is itself the
scientific output. `diag_1` has three nested published outcome-blind partitions
(ICD-9 chapters ⊂ CCS-285 ⊂ phecodes), all pre-dating the data, so the learned
tree has something external to be scored against. Also: a published
stratum-specific hypothesis (Strack's diagnosis × HbA1c interaction) to test, a
hard structural zero for calibration, and `A1Cresult`, where the right split is
almost certainly {None} vs {Norm, >7, >8} — a selection indicator masquerading as
a level, and something a p-value alone cannot say.

Watch: ~30% of rows are repeat patients. Restrict to first encounters, which also
reproduces Strack's analytic sample.

### 4c. ACS PUMS via `folktables` — was 2nd, now 3rd

Still worth doing, but justified differently. Its real value is not "modern
Adult" — it is `OCCP` × `INDP`, i.e. SOC and NAICS, **a published tree on both
sides at once**, which nothing else on the list gives. Scale is the other draw:
`dX·dY` in the hundreds × hundreds, `n` in the millions, and a genuinely rich `Z`
where stratification-based competitors cannot run at all.

Costs: no external CI ground truth whatsoever, and five survey-design objections
that each cost a paragraph. It does pre-empt the "Adult is 1994" objection, since
`ACSIncome`'s filters deliberately mirror Adult's.


## 4.1 Sanity checks — DONE
1) Does the adaptive search method beat a random set of partitions?

**Yes** (`SEARCH_COMPARISON.md`): greedy beats random by 3–25 points of power
wherever the structure leaves the final partition free; only 1–2 points under a
binary tree, where every path ends at the same root split.

## 5. Reverse the greedy procedure — DONE (merged `51f4e09`, see `DIVISIVE.md`)

Start with fully collapsed levels (ie just binary options), then expand outwards.
Should be more efficient in high dimensions?

**Answer.** Divisive search run to completion is not faster than merging, and it
equals merge in power (±0.015). Truncating it at `k` levels is very cheap (split@2:
0.06 s vs 15.6 s per test at `d = 16`) and calibrated (size 0.048–0.052, 1000 reps).
It gains up to +0.15 when the signal is coarse but has **no power** when the signal
is only between sibling labels. Pooling merge and split@k into one minP gives
exactly merge, so the gain is **entirely multiplicity**: calibrating over 3 levels
instead of 13.

Open follow-ups, both optional:

  1. **Weighted minP over depths** on the ordinary merge path (more weight on coarse
     levels). It tests the multiplicity explanation directly and needs no divisive
     search; a possible paper remark.
  1. The `alt` union cell was never run; the prediction is that it tracks merge.

## 6. Exact p-value criteria instead of approximate chisq — DONE (`0908335`)

Is it actually that slow to compute exact p-values rather than using the
approximation + trace updates?

`ExactChi` is on master (`catci_test(statistic="exact")`, registry entries
`*_exact`), and the `gammainc` speedup landed in `approx_chi_statistic`. Exact is
~28x slower than approx; p-values were identical in 17 of 18 paired reps (the miss
was one `1/(B+1)` grid step), although 14–32% of searches take a different path.
Phrase the paper claim as "same p-values in almost every case". The paired
power/size table is 0b item 2, block D. The original analysis follows.

Answered (`EXACT_CRITERIA_NOTES.md`, deleted with its worktree): yes, ~30x slower per test, and the cost is
**structural**. It is the weighted-chi-square CDF evaluation (91% of cost), not
the eigendecomposition (9%), and the CDF depends on `||t||²` so it cannot be
cached across bootstrap draws. Recommendation: do not switch.

## 7. Alternative search procedures — DONE (merged `dce3c8c`, see `SEARCH.md`)

**Answer.** All alternatives hold their size.

* **7a, random search:** costs 3–25 points of power, but only 1–2 under a binary
  tree. Fixing one random path across draws vs drawing a fresh path per draw makes
  no consistent difference (±3.5 points), so the conjecture below is not supported.
* **7b, sample splitting:** costs 8–35 points (`split_one` up to 47), usually
  falling below even the non-adaptive tests; minP pays for selection more cheaply
  than half the sample.
* **7c, beam search:** no gain even on the greedy-trap DGP (`w = 25`: +1.0 ± 1.0
  points). Greedy suffices.

The original designs follow.

### 7a. Random search

This is a simple default we should absolutely test in a couple of ways. We can randomly select various groups. Since these are done without conditioning on the test data, we will probably improve power by fixing the same groups on all the bootstrap samples (but we should check this).

### 7b. Sample splitting

Related to 7a, it is possible to split the data in two and use the first to choose the groups and the second to evaluate the test statistic and boostraps. There is some literature doing similar to this, I think the GCM extension by Anton Rask Lundborg.

### 7c. Beam search instead of greedy — is a wider search worth the null inflation?

Raised as "should this be Monte Carlo tree search?". It should not be. MCTS pays
for machinery this problem does not need: the criterion is available at *every*
node (no delayed reward), it is deterministic given `(t, Sigma)` (nothing to
average over rollouts), and the criterion is itself the heuristic. Worse, the
value at a node depends only on the partition pair reached, not the merge order
— merging row-sums commute — so the state space is a **DAG, not a tree**, and
naive UCT would spend most of its simulations re-scoring permutations of states
it had already visited. Add the transposition table that fixes that and MCTS
degenerates into a randomised best-first search, which a deterministic one
dominates. It would also make the statistic randomised, which the fixtures
cannot pin.

The question underneath it is real, though, and the right knob is **beam width
`w`**, with `w = 1` recovering today's greedy exactly:

  1. Keep the top `w` states per level in `search.py` instead of the argmax,
     dedup on canonicalised partitions, report the per-depth max as the
     criterion path. Deterministic, exactly `w`x cost, still pinnable.
  1. No new theory: like the divisive search, this is another deterministic
     `(t, sigma) -> ` sequence-of-partitions map, so Lemma 10 covers it.
  1. **The experiment that settles it.** Power at fixed size for
     `w in {1, 5, 25}`, on a DGP built to trap greedy — one where the best 2x2
     collapse needs a first merge that ranks poorly at depth 1. Note that
     `calibrate.py` ranks the observed criterion *within its own bootstrap null
     at each depth*, so absolute optimality of the search is irrelevant: a
     better optimiser lifts the criterion under `H0` and `H1` alike. Null
     inflation grows like `sqrt(2 log K)` in the number of partitions
     effectively explored; the alternative-side gain is `O(1)` only where greedy
     is genuinely trapped.

Either outcome is a paper sentence. If beam does not beat greedy even on the
trap, MCTS cannot either, and "greedy suffices because the calibration charges
for search breadth" is a stronger claim than an extra section. If it wins only
on traps, that characterises when richer search pays. Greedy's plausible failure
mode is narrow and specific — it will not merge cancelling labels, since that
lowers the criterion; the myopia is at the *coarse* end, where every candidate
touches 2 of `d` rows, differences are noise-dominated, and the choice binds
everything downstream. Which is also where item 5 says the cost lives.

## 8. Vectorise the search — DONE (2026-09-26)

Landed on master (`d5ff8fb`..`3770b79`). `search.greedy_search_paths(T (p, B), Sigma,
...) -> (L+1, B)` runs every bootstrap draw in lockstep and scores every candidate at a
level at once; `calibrate.adaptive_pvalue`, `experiments/methods.py` and `catci_test`
(new `n_jobs=`, `-1` = all cores) all go through it. `greedy_search` is its `B = 1`
wrapper. The old loop is kept as `_greedy_search_loop`: it is the oracle for
`tests/test_search_vectorised.py` and the path `ExactChi` still takes (no rank-one
update for a spectrum, and its cost is the per-draw CDF anyway).

Beyond the original plan: `cross`, `term1`, `term2` and `d_tr` are *carried* across
levels and updated from the two merged rows, rather than recomputed. There are two
modes, picked per dimension: dense (one `O(p^2)` GEMM per merge) for `Saturated`, and
sparse for trees and ordinal. Sparse tracks only the permitted pairs; a newly permitted
pair always contains the merged group and is computed fresh, which makes a level
`O(p d)` per draw. The search-path fixture passes with identical partitions, and every
experiment-registry p-value is identical to the old code on a fixed seed.

Measured, CPU per draw, `B = 1000`, vs the loop (the machine was heavily loaded, so the
absolute numbers are ~2x inflated; the ratios are the signal):

| structure | 8x8 | 16x16 | 32x32 |
|---|---:|---:|---:|
| greedy | 58x | 45x | — |
| tree | 22x | 8x | 6x |
| ordinal | 16x | 9x | 7x |

`n_jobs = 6` adds ~3x wall. `B = 10^4` at `p = 1024` is now minutes, not ~35 min.

Left open, only if it is ever needed:

  1. **`p` in the thousands.** At large `d` the cost is writing column `u` back into
     each draw's merged `Sigma` — for Y merges, a scatter of single elements. The
     lever is a compiled (numba) merge kernel. Memory is `8 p^2` bytes per draw in
     flight (chunked, so `B` is unconstrained); beyond `p ~ 4096` a factored,
     `Sigma`-free form would be needed (sketched, never verified).
  1. **Sharing work across draws does not pay.** Partition-only state could be shared,
     but by level 4 at `d = 16` nearly every null draw is on its own partition.
  1. `catci_test` still defaults to `n_boot = 100`; `1000` is now affordable.

**Forward-pass framing** (from the deleted planning doc, still relevant to item 7): a
search level is a batched bilinear form (all candidate values) followed by a pool over
candidates, and the criterion path is the sequence of pool outputs. The pooling
operator is the one knob that unifies the variants: `max` = greedy, `top-k` = beam
width `k` (item 7c), `softmax` at temperature `tau` = a relaxed search. Now that
candidates are scored as arrays, beam search is a change to the pool plus per-draw
state duplication, not a new search.

## 9. Out of scope, but worth one line in the paper

Multi-label categoricals (a movie tagged both "Action" and "Sci-Fi") break the
framework rather than extend it: the object here is a partition of a label set,
whereas multi-label makes `X` a point in `{0,1}^d`. Power-set encoding blows up;
per-label binary encoding discards the joint. State it as a limitation.
