# Why the divisive search's cost grows sub-linearly in `B`

**Finding (2026-09-30).** In the batched divisive search
(`divisive_search_paths`, `src/catci/search.py`), CPU time per test grows roughly
like `B^0.4–0.6` in the number of bootstrap draws `B`, while the merge search
(`greedy_search_paths`) grows linearly (`B^0.97`). For a truncated divisive search
the cost is close to flat: split@2 takes 0.10 s at `B = 1000` and 0.27 s at
`B = 10,000`. The reason is that **draws share partitions, and the expensive part of
scoring a partition does not depend on the draw.** Truncation keeps the search at
the coarse end, which is exactly where partitions are shared most.

## The mechanism

The approximate chi-square criterion scores a pair of partitions from three
numbers:

* `‖T‖²`, the squared norm of the coarsened statistic -- **depends on the draw**;
* `tr` and `tr²` of the coarsened covariance -- **depend only on the partition and
  on `Sigma`**, which every bootstrap draw shares.

`tr` and `tr²` are the expensive part: for every candidate split they need block
algebra on `Sigma` (`_DivisiveTable.expand`). The search therefore keeps a table
keyed by partition:

* **once per distinct partition**: `expand()` scores every permitted split for
  `(tr, tr²)` and caches the result;
* **once per draw**: only the `‖T‖²` update for each candidate,
  `normsq − 2 Σ T_a T_b`, done in one batched array call for all draws at that
  partition.

So, roughly,

    cost ≈ c_expand × (distinct partitions visited) + c_draw × B × (candidates per level)

and only the second term is linear in `B`.

## Draws pile up on the same partitions

Distinct partitions occupied at each split depth, on one Education × Income null
replicate of the semi-synthetic adult data (`n = 1000`, oracle propensities;
Education has 16 ordinal levels, so the divisive search chooses cut points, and
Income has 2):

| split depth | 0 | 1 | 2 | 3 | 4 | 5 | 6 | … | 14 | partitions expanded |
|---|---|---|---|---|---|---|---|---|---|---|
| B = 100 | 15 | 53 | 62 | 64 | 70 | 65 | 61 | … | 1 | 631 |
| B = 1000 | 15 | 103 | 284 | 363 | 380 | 361 | 328 | … | 1 | 2,759 |
| B = 10,000 | 15 | 105 | 442 | 1,076 | 1,502 | 1,515 | 1,353 | … | 1 | 8,901 |
| possible, C(15, k+1) | 15 | 105 | 455 | 1,365 | 3,003 | 5,005 | 6,435 | … | 1 | |

At depth `k` a draw has chosen `k + 1` of the 15 cut points, so at most
`C(15, k+1)` partitions are reachable. The top levels **saturate**: by `B = 1000`,
depth 1 is full (103 of 105), and beyond that further draws add no expansion work
there. The middle depths are still filling at `B = 10,000`. Partitions expanded go
631 → 2,759 → 8,901, about `B^0.6` and then `B^0.5` -- the exponent seen in the
timings.

## Measured cost

CPU seconds per search (`time.process_time`, single-threaded, median of 20
replicates; Education × Income, `n = 1000`, oracle), from
`sweep_semisynth.py` on a 6-core Apple machine:

| | B = 100 | 1000 | 10,000 | ms per draw at 10k | fitted exponent of B |
|---|---|---|---|---|---|
| split@1 | 0.016 | 0.023 | 0.134 | 0.013 | 0.47 |
| split@2 | 0.054 | 0.103 | 0.274 | 0.027 | 0.33 |
| split@4 | 0.148 | 0.548 | 1.381 | 0.138 | 0.48 |
| split@8 | 0.304 | 1.378 | 4.735 | 0.473 | 0.59 |
| full split (k = 14) | 0.385 | 1.674 | 5.652 | 0.565 | 0.58 |
| full merge | 0.032 | 0.255 | 2.656 | 0.266 | 0.97 |

The minP calibration adds about 3 ms. At `B = 10,000`, split@2 is about 10× cheaper
than the full merge and split@4 about 2× cheaper.

## Why truncation is the source of the saving

A truncated search stops at the coarse end, where partitions saturate: split@2
only ever visits depths 0–2, at most 15 + 105 + 455 partitions however large `B`
gets. Its expansion cost is essentially fixed, and only the cheap per-draw update
grows with `B`.

The merge search keeps separate state for every draw and updates it after each
merge, with no shared table. When it was vectorised (TODO 8), sharing was judged
not to pay, since by level 4 at `d = 16` nearly every null draw is on its own
partition. Its cost per draw is therefore constant (0.27 ms) and its total cost
linear.

## Caveats

* **Partly an implementation asymmetry, not purely a property of the direction.**
  For an ordinal variable, merging from 16 singletons faces the same `C(15, k)`
  partition counts, mirrored: it starts at the fine end, where they are also
  small, so a shared table could help the merge search's first few levels. But to
  reach the coarse levels, merging must pass through the middle depths, where
  partitions are most numerous. Truncated divisive search never goes there, and
  that is the real reason it is cheap.
* **0.4–0.6 is a fitted slope, not an asymptotic rate.** The cost is a
  near-fixed expansion part plus a linear per-draw part, so the log-log slope
  rises towards 1 as `B` grows and the per-draw term dominates. At `B = 10,000`
  that crossover has not been reached for shallow truncations.
* **Structure-dependent.** The counts above are for an ordinal variable, where
  depth `k` has `C(d−1, k+1)` reachable partitions. A binary tree has fewer (each
  split is forced by the tree), so saturation is faster still. `Saturated` has no
  divisive search at all (`2^(s−1) − 1` bipartitions of a group of size `s`).
* One pair and one replicate for the partition counts; the timings are medians of
  20 replicates on one pair.

## For the paper

A sentence or short paragraph in the appendix "Implementation" section, beside
"Fast greedy merging" (`sect:fastgreedymerging`, `alg:fastgreedymerging` in
`main.tex`). The paper does not yet describe the divisive search at all, so this
belongs with that write-up. The point to make: the divisive search scores each partition once however many
bootstrap draws reach it, because the costly trace terms depend only on the
partition and `Sigma`; at the coarse end, where a truncated search stays, the draws
concentrate on a small number of partitions. So a truncated divisive search costs
close to a fixed amount plus a cheap per-draw update, and increasing `B` -- which
the minP calibration rewards with power at large `L` -- is almost free.

## Reproducing

The CPU table comes from `results_semisynth/sweep_Education_Income_real_ZAge-Sex_n1000_r1000_size_timing.parquet`:

    python report_semisynth.py --pattern 'sweep_Education_Income*size'

The partition counts:

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

(run from `experiments/` with `PYTHONPATH=../src:.`).
