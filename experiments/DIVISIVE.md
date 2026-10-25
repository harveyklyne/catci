# Divisive (reverse) label search

Merged to master. Code and tests: `src/catci/search.py` (`divisive_search_paths`,
`MergeSearch`, `SplitSearch`), `blocks.py`, `structure.py`, `tests/test_divisive.py`;
their docstrings are the reference. This file keeps the conclusions.

Question addressed: *reverse the greedy procedure -- start fully collapsed (just
binary options), then expand outwards. Should be more efficient in high dimensions?*

```sh
python -m pytest tests/test_divisive.py -q
python experiments/power_divisive.py       # power/size, CSVs in experiments/results_divisive/
python experiments/bench_search.py --structure tree --dx 8 16 32 --levels 2 4 8
```

## Conclusions

* **Direction alone changes nothing.** Full split equals merge in power to within
  +-0.015 (d = 8, both DGPs) and is not faster: it visits the same `dx + dy - 3`
  levels at comparable cost. "Divisive finds better coarse partitions" is not
  supported.
* **Truncated split@k is a bet on where the signal lives.** Coarse signal: it wins
  by up to +0.15 and is far cheaper. Fine signal (`alt`, signal only between
  sibling labels): it has no power, since it never reaches that level.
* **The gain is smaller minP multiplicity.** Calibrating over 3 levels instead of
  13. Pooling merge and split@k into one minP (still exact) gives merge to the rep,
  even though split@2's level-1/2 partitions lie on the merge path only ~53% /
  ~34% of the time: the different coarse partitions carry no extra signal. (`alt`
  union not run; prediction is it tracks merge.)
* **Size is fine** everywhere (within 1 SE of 0.05; d = 16 truncated splits at 1000
  reps: 0.051 / 0.048 / 0.052 at nominal 0.05).
* **Hypothesis for follow-up:** the lever is how minP spreads its budget over
  depths. A *weighted minP over depths* (more weight on coarse levels) on the
  ordinary merge path tests this directly, needs no divisive search, and belongs
  with the calibration work (`MINP.md`).
* **Truncated divisive search is a sound option** when the analyst believes the
  effect is coarse (tree from domain knowledge), not a replacement for the full
  search.

### Power and size, d = 8

`power_divisive.py`: oracle propensities, `lin` margins, n = 1000, B = 1000, 200
reps, minP calibration, all methods on the same data and bootstrap draws.

`binary_tree` (coarse signal), rejection at 0.05:

| strength | merge | split | split@2 | split@4 | split@8 | merge+split@2 | merge+split@4 | merge+split@8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 (size) | 0.060 | 0.055 | 0.045 | 0.050 | 0.060 | 0.055 | 0.060 | 0.060 |
| 0.2 | 0.065 | 0.080 | 0.105 | 0.085 | 0.080 | | | |
| 0.4 | 0.280 | 0.270 | 0.350 | 0.365 | 0.325 | 0.280 | 0.280 | 0.280 |
| 0.6 | 0.520 | 0.520 | 0.665 | 0.645 | 0.585 | 0.520 | 0.525 | 0.510 |
| 0.8 | 0.820 | 0.820 | 0.890 | 0.865 | 0.850 | 0.820 | 0.820 | 0.815 |

Paired difference vs merge at 0.6: split@2 +0.145 +- 0.029, split@4 +0.125 +- 0.026,
split@8 +0.065 +- 0.021, full split +0.000 +- 0.012.

`alt` (fine-level signal):

| strength | merge | split | split@2 | split@4 | split@8 |
|---|---:|---:|---:|---:|---:|
| 0 (size) | 0.060 | 0.055 | 0.045 | 0.050 | 0.060 |
| 1.2 | 0.445 | 0.450 | 0.035 | 0.030 | 0.035 |
| 1.6 | 0.840 | 0.840 | 0.040 | 0.050 | 0.125 |
| 2.0 | 0.985 | 0.985 | 0.040 | 0.030 | 0.175 |

### d = 16 (binary_tree, n = 4000, B = 1000, 100 reps)

| strength | merge | split | split@2 | split@4 | split@8 |
|---|---:|---:|---:|---:|---:|
| 0.4 | 0.78 | 0.80 | 0.89 | 0.89 | 0.88 |
| 0.8 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

Paired vs merge at 0.4: split@2 / @4 +0.110 +- 0.031, split@8 +0.100 +- 0.030, full
split +0.020 +- 0.014. Seconds per test: merge 15.6, full split 58.8, split@2 0.06,
split@4 0.28, split@8 5.6. The 100-rep size excess (0.09-0.11) was noise; at 1000
reps (SE 0.007) split@2 / @4 / @8 reject 0.051 / 0.048 / 0.052 at 0.05 (0.010 /
0.010 / 0.009 at 0.01, 0.103 / 0.097 / 0.101 at 0.10).

## Why the test is valid

Theorem 1 does not need the search to be agglomerative. Lemma 10 is stated for *any*
deterministic map `Pi : t -> Pi(t)` into `C*`, and Lemma 5's continuity argument
needs only a finite candidate set per stage and an argmax selection. A divisive
search is another such map, so the double bootstrap calibrates it. Truncation is
allowed too (the paper: "`L` is assumed to be such that we stop after the `L`th
step"). Checked empirically by `test_divisive_calibrates_under_gaussian` (ordinal,
tree, truncated tree). There is no R oracle for this path; it is also pinned
against dense recomputation of `phi(Pi T, Pi Sigma Pi^T)` at every level (max abs
error 2.1e-15) and structurally against `greedy_search`.

## Limits

* **No `Saturated` divisive search.** Splitting a group of size `s` admits
  `2^(s-1) - 1` bipartitions (exponential, vs quadratic bottom-up), so
  `divisive_search` is restricted to `Ordinal` and `Tree`; `Saturated` raises.
* **Ordinal start-search cost.** Merging bottoms out at a unique singleton
  partition; splitting must choose where to start. A tree names its own coarsest
  split (one candidate); two ordinal variables offer `(dX-1)(dY-1)` (261 at
  `dX=30, dY=10`), and scoring them dominates: ordinal split@2 at `dX=30` is
  26.8 ms vs 1.25 ms for tree `d=32`. Partly an implementation artefact (a Python
  loop of ~6 numpy calls per candidate; batching into one `SigmaBlocks.blocks` call
  should collapse it), the `O(dX dY)` candidate count is intrinsic.
* **`SigmaBlocks` memory ceiling.** It allocates `(dy+1)^2 (dx+1)^2` floats: ~9.5 MB
  at `d=32`, ~143 MB at `d=64` (about 2x `Sigma` itself), capping usable `d`.
* **Contiguity.** Groups must be contiguous in label order (true for Ordinal and
  Tree, which keeps every group a range and every block a rectangle). A structure
  with non-contiguous groups would need a mask-based fallback in `blocks.py`.
* **Side change worth keeping:** `approx_chi_metric` uses
  `scipy.special.gammainc(h/2, x/2)` instead of `scipy.stats.chi2.cdf(x, df=h)`.
  Bit-identical (20,000 random pairs, max abs diff exactly 0.0) but 82.5 us vs
  1.00 us per scalar call; a free ~1.3x on the merge path too.
* Still open: d = 32 power; `alt` at d = 16 (expect the same collapse as d = 8);
  ordinal-structure power; `bench_search.py` divisive columns were not ported at
  the rebase.

## Cost: why it grows sub-linearly in `B`

**Finding (2026-09-30).** In the batched divisive search, CPU time per test grows
roughly like `B^0.4-0.6` in the number of bootstrap draws `B`, while the merge
search (`greedy_search_paths`) grows linearly (`B^0.97`). For truncated search the
cost is close to flat: split@2 takes 0.10 s at `B = 1000` and 0.27 s at
`B = 10,000`. Reason: **draws share partitions, and the expensive part of scoring a
partition does not depend on the draw.** Truncation keeps the search at the coarse
end, where partitions are shared most.

### Mechanism

The approximate chi-square criterion scores a partition pair from three numbers:
`||T||^2` (squared norm of the coarsened statistic, **depends on the draw**) and
`tr`, `tr^2` of the coarsened covariance (**depend only on the partition and on
`Sigma`**, shared by every draw). Of `(normsq, tr, tr2)` only `normsq` is
per-draw. `tr` and `tr^2` are the expensive part: every candidate split needs
block algebra on `Sigma`. The search keeps a table keyed by partition
(`_DivisiveTable`):

* **once per distinct partition**: `expand()` scores every permitted split for
  `(tr, tr^2)` in a few batched block calls and caches the result;
* **once per draw**: only the update `normsq - 2 sum(T_a T_b)` for each candidate,
  batched over all draws at that partition.

    cost ~ c_expand x (distinct partitions visited) + c_draw x B x (candidates per level)

and only the second term is linear in `B`. The per-draw loop is kept as
`_divisive_search_loop`, the reference; `test_vectorised_matches_loop` pins them
(max abs diff ~7e-15, identical partitions).

### Draws pile up on the same partitions

Distinct partitions occupied per split depth, one Education x Income null replicate
of the semi-synthetic adult data (`n = 1000`, oracle propensities; Education has 16
ordinal levels so the search chooses cut points, Income has 2):

| split depth | 0 | 1 | 2 | 3 | 4 | 5 | 6 | ... | 14 | partitions expanded |
|---|---|---|---|---|---|---|---|---|---|---|
| B = 100 | 15 | 53 | 62 | 64 | 70 | 65 | 61 | ... | 1 | 631 |
| B = 1000 | 15 | 103 | 284 | 363 | 380 | 361 | 328 | ... | 1 | 2,759 |
| B = 10,000 | 15 | 105 | 442 | 1,076 | 1,502 | 1,515 | 1,353 | ... | 1 | 8,901 |
| possible, C(15, k+1) | 15 | 105 | 455 | 1,365 | 3,003 | 5,005 | 6,435 | ... | 1 | |

At depth `k` a draw has chosen `k + 1` of 15 cut points, so at most `C(15, k+1)`
partitions are reachable. Top levels **saturate** (by `B = 1000`, depth 1 is 103 of
105), so further draws add no expansion work there; middle depths are still filling
at `B = 10,000`. Partitions expanded go 631 -> 2,759 -> 8,901, about `B^0.6` then
`B^0.5`, matching the timings.

### Measured cost

CPU seconds per search (`time.process_time`, single-threaded, median of 20
replicates; Education x Income, `n = 1000`, oracle), from `sweep_semisynth.py` on a
6-core Apple machine:

| | B = 100 | 1000 | 10,000 | ms per draw at 10k | fitted exponent of B |
|---|---|---|---|---|---|
| split@1 | 0.016 | 0.023 | 0.134 | 0.013 | 0.47 |
| split@2 | 0.054 | 0.103 | 0.274 | 0.027 | 0.33 |
| split@4 | 0.148 | 0.548 | 1.381 | 0.138 | 0.48 |
| split@8 | 0.304 | 1.378 | 4.735 | 0.473 | 0.59 |
| full split (k = 14) | 0.385 | 1.674 | 5.652 | 0.565 | 0.58 |
| full merge | 0.032 | 0.255 | 2.656 | 0.266 | 0.97 |

minP calibration adds ~3 ms. At `B = 10,000`, split@2 is ~10x cheaper than the full
merge and split@4 ~2x.

### Why truncation is the source of the saving

A truncated search stops at the coarse end, where partitions saturate: split@2 only
visits depths 0-2, at most 15 + 105 + 455 partitions however large `B`. Its
expansion cost is essentially fixed and only the cheap per-draw update grows with
`B`. The merge search keeps separate state per draw and updates it after each
merge, with no shared table; sharing was judged not to pay (TODO 8) since by level
4 at `d = 16` nearly every null draw is on its own partition. Its cost per draw is
constant (0.27 ms), total cost linear.

### Caveats

* **Partly an implementation asymmetry.** For an ordinal variable, merging from 16
  singletons faces the same `C(15, k)` counts mirrored, so a shared table could
  help its first few levels. But reaching the coarse levels means passing through
  the middle depths, where partitions are most numerous. Truncated divisive search
  never goes there, which is the real reason it is cheap.
* **0.4-0.6 is a fitted slope, not an asymptotic rate.** Near-fixed expansion plus a
  linear per-draw part, so the log-log slope rises towards 1 as `B` grows; at
  `B = 10,000` the crossover has not been reached for shallow truncations.
* **Structure-dependent.** Counts are for an ordinal variable. A binary tree has
  fewer reachable partitions (splits forced by the tree), so saturation is faster.
* One pair and one replicate for the partition counts; timings are medians of 20
  replicates on one pair.

### For the paper

A sentence or short paragraph in the appendix "Implementation" section, beside
"Fast greedy merging" (`sect:fastgreedymerging`, `alg:fastgreedymerging` in
`main.tex`). The paper does not yet describe the divisive search at all, so this
belongs with that write-up. The point: the divisive search scores each partition
once however many bootstrap draws reach it, because the costly trace terms depend
only on the partition and `Sigma`; at the coarse end, where a truncated search
stays, draws concentrate on few partitions. So a truncated divisive search costs
close to a fixed amount plus a cheap per-draw update, and increasing `B` -- which
the minP calibration rewards with power at large `L` -- is almost free.

### Reproducing

The CPU table comes from
`results_semisynth/sweep_Education_Income_real_ZAge-Sex_n1000_r1000_size_timing.parquet`:

    python report_semisynth.py --pattern 'sweep_Education_Income*size'

The partition counts (run from `experiments/` with `PYTHONPATH=../src:.`):

```python
import argparse, numpy as np
import run_semisynth as r, adult_semisynth as ss
from catci.gcm import form_t_sigma
from catci.bootstrap import bootstrap_T
from catci.search import divisive_search_paths

a = argparse.Namespace(x="Education", y="Income", w=["Age", "Sex"], z=["Age", "Sex"],
                       direction="real", n=1000, pool=5.0)
pop, delta, _, _, (xs, ys) = r._setup(a)
rep = ss.draw(pop, 1000, 0.0, np.random.default_rng(1), replace=True, delta=delta)
ts = form_t_sigma(rep.x, rep.y, rep.f_true, rep.g_true, normalise=False)
boot = bootstrap_T(ts.Sigma, 10000, np.random.default_rng(2))
for B in (100, 1000, 10000):
    T = np.column_stack([ts.T_vector, boot[:, :B]])
    _, table, visited = divisive_search_paths(T, ts.Sigma, 16, 2, xs, ys, return_ids=True)
    print(B, [len(np.unique(v)) for v in visited], len(table.nodes))
```

## Vectorised search: speedups and design

### Per-draw timings (2026-09-26 benchmark)

ms per draw, B = 1001 draws in one batch, tree, random Sigma (d = 32 and ordinal not
measured; run cut short):

| d | merge | full split | split@2 | split@4 | split@8 | old loop split@4 |
|---|---:|---:|---:|---:|---:|---:|
| 8 | 2.3 | 2.3 | 0.015 | 0.33 | 1.6 | 21 |
| 16 | 23.7 | 64.6 | 0.063 | 0.25 | 6.9 | 15 |

Full split is slower than vectorised merge at d = 16: it carries no state across
levels and fine levels have little partition sharing. Pre-vectorisation (per-draw
loop, tree, 100 draws amortised; machine contended, trust ratios to ~20%):

| structure | d | merge | split | split@2 | split@4 | split@8 |
|---|---|---:|---:|---:|---:|---:|
| tree | 8 | 2.56 | 4.34 | 0.93 | 2.02 | 3.65 |
| tree | 16 | 15.39 | 17.64 | 0.96 | 2.25 | 6.15 |
| tree | 32 | 188.80 | 199.32 | 1.25 | 2.48 | 5.38 |
| ordinal (dy=10) | 8 | 5.88 | 14.61 | 7.19 | 9.42 | 13.23 |
| ordinal (dy=10) | 16 | 14.51 | 39.02 | 14.24 | 18.55 | 26.69 |
| ordinal (dy=10) | 30 | 45.46 | 123.78 | 26.80 | 33.78 | 47.71 |

Truncated tree split is near constant in `d` (0.93 -> 0.96 -> 1.25 ms for d = 8 ->
32), ~150x cheaper than full merge at d = 32, because it never handles a partition
with more than `L + 2` groups per dimension; only the one-off block-sum table
scales with `dx*dy` (29 ms at d = 32), shared across draws.

### Merge-search vectorisation (TODO item 8, master `d5ff8fb`..`3770b79`)

`greedy_search_paths(T (p, B), Sigma, ...) -> (L+1, B)` runs every bootstrap draw in
lockstep and scores every candidate at a level at once; `adaptive_pvalue`,
`experiments/methods.py` and `catci_test` (`n_jobs=`, `-1` = all cores) use it.
`greedy_search` is its `B = 1` wrapper. The old loop is `_greedy_search_loop`: the
oracle for `tests/test_search_vectorised.py`, and the path `ExactChi` still takes
(no rank-one update for a spectrum; its cost is the per-draw CDF anyway).
`cross`, `term1`, `term2`, `d_tr` are *carried* across levels and updated from the
two merged rows, not recomputed. Two modes per dimension: dense (one `O(p^2)` GEMM
per merge) for `Saturated`; sparse for trees and ordinal, tracking only permitted
pairs (a newly permitted pair always contains the merged group, so is computed
fresh), making a level `O(p d)` per draw. Fixture passes with identical partitions;
every registry p-value identical to the old code on a fixed seed.

CPU per draw speedup vs the loop, `B = 1000` (machine loaded, absolute ~2x inflated;
ratios are the signal):

| structure | 8x8 | 16x16 | 32x32 |
|---|---:|---:|---:|
| greedy | 58x | 45x | -- |
| tree | 22x | 8x | 6x |
| ordinal | 16x | 9x | 7x |

`n_jobs = 6` adds ~3x wall. `B = 10^4` at `p = 1024` is minutes, not ~35 min. End to
end, the union run above took 237 s vs 1284 s before vectorised merge landed.

Left open, only if needed:

* **`p` in the thousands.** Cost is writing column `u` back into each draw's merged
  `Sigma` (a scatter of single elements for Y merges); lever is a compiled (numba)
  merge kernel. Memory is `8 p^2` bytes per draw in flight (chunked, so `B` is
  unconstrained); beyond `p ~ 4096` a factored, `Sigma`-free form would be needed
  (sketched, never verified).
* **Sharing work across draws does not pay for merge** (see above).
* `catci_test` still defaults to `n_boot = 100`; `1000` is now affordable.

**Forward-pass framing:** a search level is a batched bilinear form (all candidate
values) followed by a pool over candidates; the criterion path is the sequence of
pool outputs. The pooling operator unifies variants: `max` = greedy, `top-k` = beam
width `k` (TODO 7c), `softmax` at temperature `tau` = a relaxed search. Beam search
is then a change to the pool plus per-draw state duplication, not a new search.

### Planned: n-ary divisive search (TODO 5, follow-up 3; not implemented)

Today `Tree.coarsest_partitions` / `Tree.permitted_splits` raise unless every node
is binary, so taxonomies must be binarised by hand (two arbitrary choices in the
Occupation tree) and wide ones (ICD-9 chapters, ~19 under the root) cannot be used.
Theory does not need binary: minP needs only a deterministic `(T, Sigma) -> path`,
and every split adds one group, so path length is `d - 2` per variable regardless.

* *Definition:* n-ary merging run backwards. A group that is a union of `c >= 2`
  sibling subtrees (or a whole node with `c` children) may split into any two
  non-empty unions: `2^(c-1) - 1` candidates; starting partitions are the same
  bipartitions at the root. Binary nodes give exactly today's search.
* *Implementation:* keep the binary-tree speedup. The table stores groups as
  contiguous ranges with O(1) block sums from cumulative sums of `Sigma` (`_as_range`,
  `SigmaBlocks.cumulative`, `_t_prefixes`); that stays the path for `Ordinal` and
  binary `Tree`. A split like {child 1, child 3} vs {child 2} is not a range, so
  generalise a group to a sorted tuple of ranges (block sums are sums of the same
  range-sum terms; cost grows only in ranges per group). Generalise
  `Tree.permitted_splits` / `coarsest_partitions` to enumerate sibling
  bipartitions. Wide nodes: exhaustive up to arity ~10-12 (511-2047 candidates);
  beyond that a heuristic bipartition or an error naming the node. `Saturated` is
  the limiting case (one node with `d` children) and could reuse this for small `d`.
* *Tests:* all-binary n-ary tree reproduces today's paths bit-for-bit; batched
  matches per-draw loop on an n-ary tree; constant path length across draws; prefix
  property `SplitSearch(k) == SplitSearch()[:k+1]`. Then rebuild the Occupation tree
  from the 1990 Census classification without the binarisation choices and re-time.
