# Random-search baselines and depth-truncated merge: implementation plan

Status: **plan, revised after review (Section 8); awaiting go-ahead.** Branch `worktree-random-search-baselines`.
Motivating question: does the *directed* choice of each merge/split buy power, or
does the benefit come from the multiplicity budget (the depth) and the adaptive
family alone? A baseline that spends the same depth budget on data-independent
random steps answers it.

Context: `experiments/SEMISYNTH_RESULTS.md` (Adult semi-synthetic, TODO priority 1),
`experiments/SEARCH_RESULTS.md` (earlier Gaussian-limit random-merge study),
`src/catci/search.py`, `src/catci/calibrate.py`.

## 0. Decisions already made (do not reopen without cause)

1. **Depth `k` keeps the first `k + 1` path entries, for merge and split alike**, and
   the later steps are never computed. `max_levels=None` is the full path (the
   internal sentinel stays `None`).
   * Split path: entry 0 is the starting partition, then `k` splits (coarse -> fine).
   * Merge path: entry 0 is the finest partition `dx x dy`, then `k` merges (fine -> coarse).
   * So `merge@k` and `split@k` score the same *number* of partitions but different
     ones. Cross-route comparison is not like-for-like; the fair comparisons are
     directed vs random *within* a route.
2. **`mergecoarse@k` is deleted** (it computed the full merge path and kept its last
   `k + 1` entries, which is not a truncation of the search). Its historical rows in
   `SEMISYNTH_RESULTS.md` are annotated, not removed.
3. **Random searches draw a fresh path per column** (observed and every bootstrap
   draw), from one rng, independently of `T`, `Sigma` and the statistic. No fixed
   shared path.
4. **One full random path is drawn, then truncated.** The first `k + 1` entries
   must not depend on `max_levels`, so depths are nested.
5. **Candidates are pooled over x and y** at each step, as in the directed search,
   and the random choice is uniform over what the *structure* permits (not over
   what the criterion scores finite).
6. **Path generation and path scoring are separate.** A random search never scores a
   candidate it does not take. The scorer is parameterised by the statistic, so the
   Box approximation can be swapped (exact, ...) without touching the generator.
7. **Depth-slotted only for now**, with a documented extension point for general
   candidate sets (Section 5).

## 1. Validity argument (what the design must preserve)

`double_bootstrap_pvalue` (`calibrate.py:76`) takes an `(L, B+1)` matrix and ranks
each row across columns (prepivot), takes the column minimum, and ranks that.
Its guarantee needs only: **the columns are exchangeable under the null.** The
observed `T` and the bootstrap draws are exchangeable in the Gaussian limit, so
any map `M(T)` -- possibly randomised with randomness independent of `T`, and
using `Sigma_hat` (shared by all columns) -- applied identically to every column
gives exchangeable columns. A path is not required.

Invariants every generator must satisfy (and tests must check):

* I1. Depends only on `(T, Sigma_hat, structures)` and independent randomness. Never
  raw rows (the bootstrap draws have none).
* I2. The same map for the observed column and each bootstrap column.
* I3. Row `l` means the same kind of object in every column (here: depth `l`).
* I4. Every column returns the same number of rows (the existing code raises otherwise).

Random path generation depends on none of `T`, `Sigma`, `statistic`, so I1-I2 hold
trivially for the generator; the scorer is a deterministic function of `(T, Sigma,
partition, statistic)`.

## 2. Code changes in `src/catci/`

### 2.1 `MergeSearch(max_levels=None)` (in `search.py`)

* `_run_chunk(..., max_levels=None)`: loop `while max_levels is None or level < max_levels`.
* `greedy_search_paths(..., max_levels=None)` passes it through; `MergeSearch.paths`
  forwards `self.max_levels`.
* `MergeSearch.result` with `max_levels` set: `NotImplementedError` (no caller needs it).
* Non-`ApproxChi` statistics (`_greedy_search_loop`): `NotImplementedError` when
  truncated, unless trivial to support.
* Tests: `paths(max_levels=k)` equals `paths()[:k+1]` bit-for-bit for Ordinal,
  Saturated, Tree; `k` beyond path length returns the full path; shape `(k+1, B)`.
* **No other change to `_run_chunk`, `greedy_search_paths`, `divisive_search_paths`.**
  (A first attempt put a random chooser inside these kernels; it scored every
  candidate to discard it, tied the baseline to the Box criterion, and made the
  choice depend on the scores. It has been reverted.)

### 2.2 New module `src/catci/randomsearch.py`

Public: `RandomMergeSearch(max_levels=None, seed=None)`,
`RandomSplitSearch(max_levels=None, seed=None)`; both expose
`.paths(T, Sigma, dx, dy, x_structure, y_structure, statistic=None, n_jobs=1)`
returning `(k+1, B)` -- the same contract as `MergeSearch`/`SplitSearch`.

Internals, three pieces:

**(a) Generator.** `generate(...) -> PathSet`, `PathSet = (node_of: int array (L+1, B),
partitions: list of canonical partition keys)`. A partition key is
`(tuple of x groups, tuple of y groups)`, each group a tuple of 1-based labels in a
canonical order.

* Level-synchronous. At each level, group columns by current partition key (data-free),
  enumerate that partition's permitted candidates once, and draw every column's
  choice in one vectorised `rng.integers(n_cand, size=m)`. Distinct keys are visited
  in sorted order and columns in index order, so the stream consumed is a function
  of the data-free path structure only. This gives nesting: stopping at `max_levels`
  stops consuming the stream, the first `k+1` levels are unchanged.
* *Merge:* start = finest partition (no choice). Candidates = `x_structure.permitted_merges(px)`
  tagged dim 1, plus `y_structure.permitted_merges(py)` tagged dim 2, pooled
  (same set as `random_merges`, `search.py:867`). Stops when both dimensions have 2
  groups: `dx + dy - 4` merges, `dx + dy - 3` entries -- equal to the greedy length
  (to be asserted in a test).
* *Split:* start is **not fixed** in the directed search -- level 0 picks the best of
  all pairs in `x_structure.coarsest_partitions(dx) x y_structure.coarsest_partitions(dy)`
  (`(dx-1)(dy-1)` for two ordinals, 1 for a binary tree pair). The random start is
  uniform over that product. Each later step is uniform over the splits of every
  group of both dimensions, using the same enumeration as the directed loop
  (`_splits_of`, `search.py:1019`; it requires Ordinal / binary Tree). Ends when
  every group is a singleton: `dx + dy - 3` entries untruncated.
* Candidate enumeration is memoised per partition key (`permitted_*` are pure
  functions of the partition).

**(b) Scorer.** `score(T, Sigma, pathset, statistic) -> (L+1, B)`.

* Only distinct partitions are scored (key: groups as sorted label tuples, groups
  ordered by least label, so the same partition reached by different merge orders
  shares one entry).
* **Two routes.** For `Ordinal` and binary `Tree` every partition is a set of
  contiguous label ranges (merges join adjacent groups; trees are recoded to leaf
  order), so the *range route* reuses the batched `SigmaBlocks` / prefix-sum
  machinery of `_DivisiveTable.expand` (`t_blocks`, `sigma_blocks.blocks`) for
  `tr`, `tr2`, `||Pi T||^2`. `Saturated` merge (non-contiguous groups) uses a *dense
  route*: build `Pi`, form `Pi Sigma Pi^T`. Both agree with `dense_values` (test 10).
* Cost is a design constraint, not an afterthought: on Education x Income (dy = 2) a
  random merge path can visit ~C(15, l) distinct partitions per level, ~32k per
  replicate at B = 10000. The per-partition work must be the batched range route;
  a Python loop building `Pi` is not acceptable. Profile before and after.
* `ApproxChi` route: group columns by partition and call `approx_chi_array`
  (`statistic.py:120`) vectorised.
* Generic route: any `statistic` with `init(PT, PSigmaPt)` / `value(state)` (the
  interface `dense_values` in `tests/test_divisive.py` uses). Slower, correct. Add a
  fast route for another statistic only if profiling at `B = 10000` says so.
* A degenerate partition (what directed search maps to `-inf`/0): the scorer's
  convention applies (score 0, as `approx_chi_array` after `a9c0ba9`); the generator
  does not know about it.

**(c) Seeding.** `seed` -> `np.random.default_rng(seed)` held by the instance; each
`.paths` call advances the stream. Callers pass a `SeedSequence` child (not a
`Generator`: `default_rng(generator)` returns the same generator, the trap documented
at `run_semisynth.py:186-189`) and construct a **fresh instance per learner**, so the
oracle and fitted arms see identical random paths. Chunking and threads play no role
in generation, so output does not depend on `n_jobs`.

**(d) `result` and `api`.** `api.py:97` calls `search.result(...)`. The random classes
and a truncated `MergeSearch` raise a clear `NotImplementedError` from `result`
(column 0 of the path set would give the partitions; implement only if needed).
`Saturated`/n-ary on the split side raises the same error `SplitSearch` does.

### 2.3 `__all__`/docs

Export from `catci/__init__` (or `catci.randomsearch`) only; `randomsearch.py`
imports `search` internals, so re-exporting from `search.py` would be circular.
`random_merges`/`evaluate_paths` (fixed-path, `search_study.py`) are kept untouched.

## 3. Harness changes in `experiments/`

* `sweep_semisynth.py`
  * Remove `mergecoarse`. **Keep its slot in `METHOD_IDS`** (rename the entry to a
    placeholder such as `"_retired_mergecoarse"`) so every other method's tie-break
    stream is unchanged; append `"rsplit"`, `"rmerge"` at the **end** of the list.
  * `merge@k` is `merge[:k+1]` of the full merge path the sweep already computes
    (prefix property); a separate truncated run is only for timing. `ks` always
    contains `L`, so skip `merge@L` (duplicate of the existing full `merge` row).
    `rsplit@k`/`rmerge@k`: one random path per replicate generated to the **full**
    length `L`, then prefixes.
  * Fix the timing fix-up at line 89, which assigns `k = merge.shape[0]-1` to any
    non-split method: key it on method name for the new rows.
  * Random path seed: a child of the replicate seed that does **not** shift the
    existing `spawn(3)` children (check: `Generator.spawn(4)[:3]` reproduce
    `spawn(3)`; assert it in a test, do not assume).
  * `argmax_level` and timing columns for the new methods as for the others. Direct
    timing runs (`--time-reps`) for `merge@k` (now genuinely truncated) and for the
    random searches.
* `run_semisynth.py`: add `rmerge@k`, `rsplit@k`, truncated `merge@k` as a **separate
  block after the existing methods**. The existing shared tie-break `rng`
  (`run_semisynth.py:202`) is consumed in dict order, then by `DEPTH0` (216-227);
  inserting methods before `DEPTH0` would change `max`/`euclid`/`mGCM` p-values. The
  new block uses its **own** tie-break rng from a separate spawned child (`spawn(4)`,
  child 3; children 0-2 verified unchanged). Random-path seeding as in 2.2(c).
* `gauss_check.py:55`: drop `mergecoarse`, use `merge[:k+1]` / the new classes.
* `report_semisynth.py` **needs code**: line 68 picks the paired baseline as the
  first column starting with `merge@`, alphabetically, so adding `merge@1..` silently
  makes `merge@1` the baseline for every "vs merge" column. Select the full-depth
  merge column explicitly, and add the within-route paired differences
  (`split@k - rsplit@k`, `merge@k - rmerge@k`, paired SE) that Section 6 reports.
* `SEMISYNTH_RESULTS.md`: annotate every `mergecoarse` mention as historical -- the
  bullets at lines ~96-98 and ~126-129, the "net so far" sentence, and the
  `gauss_check` table (~140) -- and withdraw the conclusions resting on it until
  `merge@k` has run. `merge@k` at small `k` sits at the fine, high-df end, so it
  says little about direction; the coarse-end comparator is gone by decision. If a
  coarse-end random reference is wanted later, a random merge path scored on its
  coarsest rows is valid (row meaning is column-independent).
* Delete the stale untracked `tests/test_random_search.py` (the reverted design;
  it fails to import) and write it afresh to Section 4.

## 4. Tests (`tests/test_random_search.py`, extend `test_search_alternatives.py`)

Generator (no data needed):
1. Same seed, different `T`/`Sigma` -> identical `node_of`/partitions (I1).
2. Nesting: `max_levels=k` equals the first `k+1` levels of the full run, same seed.
3. Fresh per column: identical columns of `T` get differing paths.
4. Seeded: same seed identical, different seed differs; independent of `n_jobs`.
5. Lengths: untruncated length equals the directed length, for Ordinal, Saturated
   (merge only), binary Tree; both dimensions.
6. Every step is a permitted merge/split of the current partition (replay validator).
7. **Uniformity:** over many seeds, the split start frequencies over
   `(dx-1)(dy-1)` pairs, and the first-step choice frequencies, pass a chi-square
   test against uniform (generic Pearson test, not a hand tolerance).
8. Pooling: with `dx > 2, dy > 2` both dimensions are chosen with probability
   proportional to their candidate counts.

Scorer:
9. **Equivalence to the directed code.** Take the partitions a *directed* search
   visits for one draw (merge: `greedy_search` partitions; split: `divisive_search`
   partitions) and feed them as a degenerate `PathSet`; scorer values equal the
   directed values (<= 1e-10). This is the test that catches an index or ordering
   error in the new scorer.
10. Scorer values equal a dense recompute (`dense_values`) at every random partition.
11. Statistic swap: a stub statistic through the generic route gives expected values.

Calibration (generic, no hand tolerances beyond 3 binomial SE, as
`_check_uniform` in `test_search_alternatives.py`):
12. Under Gaussian `T`, `p`-values uniform for `RandomMergeSearch`, `RandomSplitSearch`,
    each full and truncated, Tree and Ordinal. 400 replicates is marginal
    (SE of mean p = 0.014); use 800+ or fixed seeds chosen before looking.
13. Truncated `MergeSearch` calibrates.

Harness:
14. Existing methods' p-values in `sweep_semisynth` are bit-for-bit unchanged by
    adding the new methods (golden replicate before/after).
15. `spawn(4)` children 0-2 equal `spawn(3)` children.

Additional tests from review:
16. `dy = 2` (y has no steps) and `L = 0` (`dx = dy = 2`); `k = 0`.
17. Column 0 has the same first-step distribution as the others (catches an
    index-dependent draw).
18. Saturated / n-ary on the split side raises a clear error.
19. Add the new searches to `tests/test_degenerate.py:50` (zero-variance coordinates).
20. Test 3 (fresh per column): use large `dx`, `dy` and assert paths are not all equal.
21. Calibration (test 12) sharpened by an exact-grid check at small `B` (e.g. `B = 9`,
    p-values on a known grid, chi-square on counts), after the precedent of
    `test_random_fresh_path_calibrates`; keep the rejection-rate check as a smoke test.
22. `report_semisynth` baseline selection with `merge@1..L` present picks `merge@L`.

Regression: full existing suite passes (`pytest tests/`), with `PYTHONPATH=src`
(the editable install points at the main checkout, not this worktree).

## 5. Extension point (documented, not implemented)

Module docstring of `randomsearch.py` and the `.paths` contract in `search.py`
state the general interface:

> A search is `(generator, scorer)`. The calibrator consumes an `(S, B+1)` matrix of
> slot values and requires only exchangeable columns (Section 1, I1-I4). A
> depth-slotted path is the special case S = depth. A generator may return any set of
> candidates per column provided each row ("slot") has a column-independent meaning:
> a fixed named family, or an order-statistic rule such as "the best candidate with
> exactly g groups" or "the l-th best candidate". "The l-th candidate found" in an
> arbitrary generation order is **not** a valid slot. Missing slots need an explicit
> rule (p = 1) and every column must return the same S. The generator may use
> `Sigma_hat` but not the raw rows, so sample-splitting selectors are out of scope.

Not planned here: directed search with a different selection statistic from the
evaluation statistic (requires the directed kernels to emit partitions per column
for the scorer). Record as a follow-up.

## 6. Experiments (after code and tests pass; fold into TODO priority 1)

* Size, null, `Education x Income`, `n = 1000`, 1000 reps, `B = 1000`: random arms at
  `k in {1, 2, 4, 8, full}`. Judge against the directed methods' finite-n excess
  (~1-1.5 points on real `T`), not against 0.05 alone; the Gaussian arm of
  `gauss_check.py` is the exact check.
* Power, real and planted, same `lam` grid and **same data/bootstrap draws** as
  the existing directed rows (paired -- requires the same `--lams`, `--reps`, `--seed`,
  since the replicate seed index is `i * reps + r`), reporting paired differences with SE:
  `split@k - rsplit@k`, `merge@k - rmerge@k`.
* Then `HoursPerWeek x Education` (n = 2000) and `Occupation x Income` (after the
  n-ary divisive search).
* Report both within-route gaps and the depth profile of each (does random power rise
  or fall with `k`).
* Expectation to check, not assume: little directed advantage under a binary tree
  for *split* (every divisive path starts at the tree's one root split,
  `SEARCH_RESULTS.md`); this does not apply to merge.
* "Matched depth" matches the number of partitions scored (rows entering minP), not
  their degrees of freedom, cell counts or effective multiplicity. State this in any
  write-up: a random path may spend its rows on partitions that are nearly collinear
  or low-signal.

## 7. Risks and open points for the reviewer

* The `Generator.spawn` independence claim (3.) is asserted, not yet verified.
* Cost of the scorer at `B = 10000`, `k` large: distinct partitions grow; the
  `Pi Sigma Pi^T` per distinct partition may dominate. Mitigation: cache by key;
  profile before optimising.
* `Cyclic` structure: out of scope unless `permitted_*` / `coarsest_partitions` already
  support it; the split side requires Ordinal / binary Tree.
* Does uniform-over-permitted at the split *start* for Ordinal (uniform over
  `(dx-1)(dy-1)` first cuts) make `rsplit@0` a different object from "no information"?
  It is what the directed search chooses among, so it is the matched baseline.
* The first random draw for merge pools x and y candidates by count, so a long
  dimension gets proportionally more steps. This matches `random_merges` and the
  directed pooled argmax, but it is a modelling choice.

## 8. Review log (independent review, 2026-10-25)

Verified correct: the exchangeability argument; level-synchronous grouped draws are
iid uniform per column and nest; `Generator.spawn(4)[:3] == spawn(3)` (numpy 2.5.1,
both `SeedSequence` and `default_rng`); `METHOD_IDS` placeholder preserves other
streams; merge and split path lengths are both `dx + dy - 3`; split start is the
`coarsest_partitions` product; `_splits_of` restricts to Ordinal / binary Tree;
`random_merges` pools as claimed; truncated merge is bit-for-bit a prefix; no other
caller of `_run_chunk` / `greedy_search_paths` breaks.

Changes made in response: `report_semisynth` baseline bug and paired-difference code
(Section 3); stale test file (Section 3); seeding through `SeedSequence` and a fresh
instance per learner (2.2c); separate tie-break rng for new methods (Section 3);
`mergecoarse` annotations broadened (Section 3); range-route scorer for Ordinal/Tree
with the cost estimate (2.2b); `result`/export/path-length/timing nits (2.2d, 2.3,
Section 3); sorted-label partition keys (2.2b); extra tests 16-22 (Section 4);
size-judgement, pairing, binary-tree and matched-depth wording (Section 6).
