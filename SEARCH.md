# Alternative search procedures (TODO item 7)

Branch `worktree-alt-search`. Three alternatives to the greedy argmax, each run
through the same minP calibration and compared with greedy at fixed size:

* **7a random search** — a merge path drawn without looking at the data, either
  *fixed* (one path per replicate, shared by the observed statistic and every
  bootstrap draw) or *fresh* (an independent path for every draw);
* **7b sample splitting** — greedy picks the path on half the data; the path
  (`split`) or just its best-looking depth (`split_one`) is tested on the other half;
* **7c beam search** — keep the top `w` partition pairs per level; `w = 1` is greedy.

## Answer

| | verdict |
|---|---|
| 7c beam | **No gain, even on a DGP built to trap greedy.** `w = 25` is +1.0 ± 1.0 points; `w = 5` is ±0.5. Greedy suffices. |
| 7a random | Valid, but **costs 3–25 points of power** (7–25 at the higher signal strength) — except under a binary tree structure, where it costs 1–2, because every tree path ends at the same root split. |
| 7a fixed vs fresh | **No consistent difference** (within ±3.5 points, sign flips between cells). The conjecture that fixing the path helps is not supported. |
| 7b split | **Costs 8–35 points everywhere** (`split_one` up to 47), usually falling below the non-adaptive tests. `split_one` is no better. Not worth pursuing: minP already pays for selection, more cheaply than half the sample. |
| 4.1 sanity check | Adaptive search beats random by 3–25 points wherever the structure leaves the endpoint free (10–25 at the higher strength). |

All methods hold size (table below). Every alternative is covered by the existing
validity argument: the minP calibration needs only that every draw goes through
the same map `(T, Sigma) -> path`, randomised if at all independently of `T`.

## Study design

`experiments/search_study.py`, Gaussian limit: `T ~ N(mu, Sigma)` with `Sigma`
known. This is what the search and calibration see asymptotically, it removes the
propensity learner from the comparison, and it makes splitting exact
(`T_A, T_B ~ N(mu/sqrt 2, Sigma)` independently, full sample `(T_A + T_B)/sqrt 2`).

* `Sigma = C (x) C`, `C = diag(p) - p p^T`, uniform `p`, `dx = dy = 8` (`L = 13`).
* `mu` points along an interaction matrix, scaled so `mu' Sigma^+ mu = delta^2`:
  `binary_tree` and `step` from `experiments/dgp.py`, and `greedy_trap` (below).
* Every method in a replicate sees the same data and bootstrap draws, so
  differences from greedy are **paired**: `vs greedy ± se` below is the mean and SE
  of the per-replicate difference of reject indicators (`summarise_search.py`).
* α = 0.05. 300 reps and `n_boot = 300` per cell unless noted.

**The trap.** `greedy_trap` is a coarse 2×2 checkerboard (amplitude `a`) plus a
stronger fine signal (`b`) on row pairs `(k, k+4)` that straddle the coarse split.
Greedy merges the straddling pairs first (largest inner products) and so cancels
the coarse signal along X for good, yet in a narrow window of `b` the coarse 2×2
is the better final collapse. Scanning `b` noise-free, `b = 3.6a` is the worst case:
a width-25 beam beats greedy's criterion by 0.35 (CDF units) at depth 6 of 13.
Outside `3.2a–3.6a` the noise-free gap is ≤ 0.06, and at `d = 6` there is
essentially no trap at all — itself evidence that greedy is rarely trapped.

## Results

Power, and difference from greedy (paired). Rows at `delta = 0` are size.

### 7c beam, greedy_trap / all pairs (δ = 5, 200 reps)

| method | rate | vs greedy |
|---|---:|---:|
| greedy | 0.595 | — |
| beam5 | 0.600 | +0.005 ± 0.011 |
| beam25 | 0.605 | +0.010 ± 0.010 |

beam5's size, from the next cell: 0.050 (greedy 0.053).

### greedy_trap / all pairs

| method | δ=0 | δ=4 | vs greedy | δ=5 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.373 | — | 0.673 | — |
| beam5 | 0.050 | 0.370 | −0.003 ± 0.007 | 0.670 | −0.003 ± 0.010 |
| random_fixed | 0.063 | 0.340 | −0.033 ± 0.023 | 0.547 | −0.127 ± 0.022 |
| random_fresh | 0.060 | 0.330 | −0.043 ± 0.022 | 0.560 | −0.113 ± 0.025 |
| split | 0.027 | 0.187 | −0.187 ± 0.027 | 0.327 | −0.347 ± 0.031 |
| split_one | 0.043 | 0.140 | −0.233 ± 0.027 | 0.200 | −0.473 ± 0.030 |
| depth0 | 0.030 | 0.400 | +0.027 ± 0.021 | 0.720 | +0.047 ± 0.019 |
| chi_sq | 0.027 | 0.413 | +0.040 ± 0.021 | 0.733 | +0.060 ± 0.019 |

The trap spreads signal over most cells, so the non-adaptive tests win here.

### binary_tree / tree structure (400 reps, `n_boot = 500`)

| method | δ=0 | δ=2 | vs greedy | δ=3 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.393 | — | 0.762 | — |
| random_fixed | 0.062 | 0.378 | −0.015 ± 0.013 | 0.752 | −0.010 ± 0.008 |
| random_fresh | 0.055 | 0.370 | −0.022 ± 0.013 | 0.743 | −0.020 ± 0.008 |
| split | 0.065 | 0.185 | −0.207 ± 0.023 | 0.435 | −0.328 ± 0.026 |
| split_one | 0.077 | 0.205 | −0.188 ± 0.025 | 0.412 | −0.350 ± 0.025 |
| chi_sq | 0.033 | 0.115 | −0.278 ± 0.025 | 0.247 | −0.515 ± 0.026 |

This cell ran on the pre-merge code (the study's own Kronecker-`Sigma` search,
tested path-for-path against `catci.search`); same methods, same design.

### step / ordinal structure (`n_boot = 500`)

| method | δ=0 | δ=3 | vs greedy | δ=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.397 | — | 0.713 | — |
| random_fixed | 0.040 | 0.263 | −0.133 ± 0.028 | 0.533 | −0.180 ± 0.028 |
| random_fresh | 0.043 | 0.297 | −0.100 ± 0.027 | 0.557 | −0.157 ± 0.030 |
| split | 0.057 | 0.207 | −0.190 ± 0.029 | 0.450 | −0.263 ± 0.029 |
| split_one | 0.050 | 0.237 | −0.160 ± 0.029 | 0.423 | −0.290 ± 0.029 |
| chi_sq | 0.027 | 0.227 | −0.170 ± 0.026 | 0.467 | −0.247 ± 0.030 |

### step / all pairs

| method | δ=0 | δ=3 | vs greedy | δ=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.257 | — | 0.533 | — |
| random_fixed | 0.060 | 0.160 | −0.097 ± 0.023 | 0.330 | −0.203 ± 0.030 |
| random_fresh | 0.033 | 0.187 | −0.070 ± 0.023 | 0.297 | −0.237 ± 0.028 |
| split | 0.023 | 0.157 | −0.100 ± 0.024 | 0.307 | −0.227 ± 0.028 |
| split_one | 0.043 | 0.127 | −0.130 ± 0.025 | 0.267 | −0.267 ± 0.029 |
| chi_sq | 0.027 | 0.227 | −0.030 ± 0.018 | 0.467 | −0.067 ± 0.025 |

### binary_tree / all pairs

| method | δ=0 | δ=3 | vs greedy | δ=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.253 | — | 0.577 | — |
| random_fixed | 0.060 | 0.167 | −0.087 ± 0.025 | 0.323 | −0.253 ± 0.028 |
| random_fresh | 0.033 | 0.183 | −0.070 ± 0.021 | 0.337 | −0.240 ± 0.029 |
| split | 0.023 | 0.170 | −0.083 ± 0.023 | 0.303 | −0.273 ± 0.030 |
| split_one | 0.043 | 0.110 | −0.143 ± 0.024 | 0.300 | −0.277 ± 0.029 |
| chi_sq | 0.027 | 0.237 | −0.017 ± 0.020 | 0.473 | −0.103 ± 0.024 |

The δ = 0 rows of the two "all pairs" cells are the same replicates (at `mu = 0`
the DGP does not enter, and both cells used seed 0) — one size check, not two.

**Size.** Across the independent null cells every method is within ~2 SE of 0.05
(largest: `split_one` 0.077 ± 0.013 in the tree cell; it is exact by construction
in this setting). `depth0` and `chi_sq` run slightly conservative, as expected of a
Box-approximated / oracle χ² at `d = 8`.

## Side observations worth a paper sentence

* **The structure, not the search, carries the tree result.** Under a binary tree
  a random path is within 1–2 points of greedy. Under ordinal or all-pairs
  structures the search is doing real work (7–25 points).
* **With all pairs permitted, greedy's edge over the non-adaptive χ² is small at
  these strengths** (2–10 points on the structured DGPs, negative on the dense
  trap). The adaptive gain lives mostly in the structured searches.

## Caveats

* Gaussian limit with oracle `Sigma`, `d = 8`, uniform marginals, one seed per
  cell. The data-level registry entries exist (`<search>_random`, `<search>_split`,
  `greedy_beam5` in `experiments/methods.py`) but were not run through the fitted
  propensity pipeline; `_split` there reuses full-sample propensities, so its halves
  are independent only up to that shared fit.
* The machine was heavily contended throughout; the plan was cut from 400–500 reps
  to 200–300 and from `n_boot = 1000` to 300–500. Paired SEs are 1–3 points, which
  is enough for every conclusion above.
* The trap is one construction. That greedy is trappable only in a narrow window,
  and that a 25-wide beam still gains nothing there, is the evidence; it is not a
  proof that no DGP rewards breadth.

## Code

* `src/catci/search.py`: `beam_search` / `beam_search_paths` (on the vectorised
  greedy kernel; `width = 1` is `greedy_search_paths` bit for bit), `random_merges`,
  `evaluate_path` / `evaluate_paths`, `SearchResult.merges`.
* `src/catci/calibrate.py`: `adaptive_pvalue(..., search=)` takes any batched
  path map.
* Tests: `tests/test_search_alternatives.py` (w=1 identity on the R fixture,
  exhaustive check of a wide beam at `d = 4`, batching, replay, calibration of
  beam / fixed / fresh random under a Gaussian null).

Reproduce (from the repo root, `PYTHONPATH=src:experiments`):

```
python experiments/search_study.py greedy_trap --structure saturated --deltas 5 --widths 1,5,25 --reps 200 --n-boot 300 --tag _beam
python experiments/search_study.py greedy_trap --structure saturated --deltas 0,4,5 --widths 1,5 --reps 300 --n-boot 300
python experiments/search_study.py step --structure ordinal --deltas 0,3,4 --widths 1 --reps 300 --n-boot 500
python experiments/search_study.py step --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py binary_tree --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py binary_tree --structure tree --deltas 0,2,3 --widths 1 --reps 400 --n-boot 500
python experiments/summarise_search.py experiments/results/search_*.parquet
```
