# Search procedures: greedy vs beam vs random vs sample-split

Merged from the former `SEARCH.md` (TODO 7c, beam) and `SEARCH_COMPARISON.md`
(TODO 7a, 7b, 4.1). Question: does the adaptive greedy search earn its keep,
compared with a wider beam, a data-independent random path, or a path chosen on an
independent half of the data?

**Status: exploratory, not final.** Cells differ in reps and bootstrap size
(300-400 reps, `B` = 300-500; the machine was heavily contended, so the plan was
cut from 400-500 reps and `B` = 1000), one cell ran on older (equivalent) code,
and the two all-pairs cells share null replicates. Qualitative conclusions are
stable at the reported paired SEs (1-3 points); do not quote the numbers in the
paper before the rerun in [Before the numbers go in the paper](#before-the-numbers-go-in-the-paper).

## Answer

1. **Every method holds its size** (within ~2 SE of alpha = 0.05; largest
   `split_one` 0.077 +- 0.013 in the tree cell). `depth0` and `chi_sq` run slightly
   conservative, as expected of a Box-approximated / oracle chi-square at `d = 8`.
2. **Beam (7c): no gain, even on a DGP built to trap greedy.** `w = 25` is
   +1.0 +- 1.0 points, `w = 5` is +-0.5. Greedy suffices.
3. **Random (7a): valid, but costs 3-25 points of power** (7-25 at the higher
   strength) wherever the structure leaves the endpoint free (ordinal, all pairs;
   3-13 on the diffuse signal). Under a binary tree it costs 1-2, because every
   tree path ends at the same root split: the structure, not the search, supplies
   the adaptivity. This is the sanity check of TODO 4.1.
4. **Fixed vs fresh random path: no consistent difference** (within +-3.5 points,
   sign flips between cells). The conjecture that fixing the path helps is not
   supported; fixed is simpler to state.
5. **Sample splitting (7b) costs 8-35 points everywhere** (`split_one` up to 47),
   usually falling below the non-adaptive tests. `split_one` is no better. minP
   already pays for selection, more cheaply than half the sample.
6. With all pairs permitted, greedy's edge over the non-adaptive chi-square is
   small (2-10 points on structured signals, negative on the diffuse trap signal).
   Under tree / ordinal searches greedy beats it by 17-52 points.

**Validity.** Every alternative is covered by the existing argument: the minP
calibration needs only that every draw (observed and bootstrap) goes through the
same map `(T, Sigma) -> path`, randomised if at all independently of `T`.

## Setup

`experiments/search_study.py`, Gaussian limit: `T ~ N(mu, Sigma)` with `Sigma`
known. This is what the search and calibration see asymptotically, it removes the
propensity learner from a comparison about the *search*, and it makes splitting
exact: `T_A, T_B ~ N(mu/sqrt 2, Sigma)` independently, full sample
`T = (T_A + T_B)/sqrt 2`.

* `Sigma = C (x) C`, `C = diag(p) - p p^T`, uniform `p`, `dx = dy = d = 8`
  (`L = 2d - 4 = 12` merges, 13 levels); `rank(Sigma) = (d-1)^2 = 49`.
* `mu = delta * u / ||u||_{Sigma^+}`, so the noncentrality `mu' Sigma^+ mu = delta^2`
  for every shape `u` (vectorised, X-fastest, double-centred interaction matrix):
  power is comparable across shapes at equal `delta`.
* Criterion: Box's approximate chi-square CDF of the merged statistic,
  `S_pi(T) = F_Box(||M_pi T||^2; tr(Sigma_pi), tr(Sigma_pi^2))`, `Sigma_pi = M_pi Sigma M_pi^T`.
* Calibration: the same minP double bootstrap for every method: `B` draws
  `Z_b ~ N(0, Sigma)`, the same path construction applied to each, stage 1
  prepivots each level against the other paths, stage 2 calibrates the minimum.
  alpha = 0.05.
* Within a replicate every method sees the same data and bootstrap draws, so
  differences from greedy are **paired**: "vs greedy" is the mean +- SE of the
  per-replicate difference of reject indicators (`summarise_search.py`). Rate SEs
  are `sqrt(r(1-r)/reps)`: <= 0.029 at 300 reps, <= 0.025 at 400.

### Methods

| name | path for the observed `T` | path for draw `Z_b` | validity |
|---|---|---|---|
| `greedy` | greedy on `T`: each level takes the permitted merge maximising `S` | greedy on `Z_b` | same deterministic map for every draw |
| `beamW` | keep the top `W` partition pairs per level (`W = 1` is greedy) | same | same deterministic map |
| `random_fixed` | `pi` drawn at random, independent of the data | the same `pi` | conditional on `pi`, a fixed map |
| `random_fresh` | `pi^(0)` drawn at random | independent `pi^(b)` | `(T, pi^(0))`, `(Z_b, pi^(b))` i.i.d. under the null |
| `split` | greedy on `T_A`; whole path scored on `T_B` | same path, scored on `Z_b` | path depends on `T_A` only, independent of `T_B`, `Z_b` |
| `split_one` | as `split`, but only `l* = argmax_l S(T_A)` tested on `T_B` (`L = 1`) | level `l*` of the same path | as `split` |
| `depth0` | finest partition only (`L = 1`) | same | non-adaptive anchor |
| `chi_sq` | `T' Sigma^+ T` vs `chi^2_{(d-1)^2}` | none | oracle non-adaptive anchor |

"Permitted" merges are set by the structure: **tree** = sibling merges in a
balanced binary tree over ordered labels; **ordinal** = adjacent groups only;
**all pairs** (`saturated`) = any two groups.

### Signals

| name | shape `u` | structure(s) |
|---|---|---|
| binary tree | hierarchical interaction on a balanced binary tree: coarse 2x2 split, finer levels weighted `0.1^(level-1)`, finest level 0 (`dgp.binary_tree_interaction`) | tree; all pairs |
| step | +-1 blocks: X split low / high halves, sign pattern flipping across the quarters of Y (`dgp.get_int("step")`) | ordinal; all pairs |
| greedy_trap ("diffuse") | coarse 2x2 checkerboard plus a stronger fine signal on row pairs straddling it (`search_study.trap_interaction`) | all pairs |

**The trap.** `greedy_trap` is a coarse 2x2 checkerboard (amplitude `a`) plus a
stronger fine signal (`b`) on row pairs `(k, k+4)` that straddle the coarse split.
Greedy merges the straddling pairs first (largest inner products), cancelling the
coarse signal along X for good, yet in a narrow window of `b` the coarse 2x2 is
the better final collapse. Scanning `b` noise-free, `b = 3.6a` is the worst case:
a width-25 beam beats greedy's criterion by 0.35 (CDF units) at depth 6 of 13.
**Narrow window caveat:** outside `3.2a-3.6a` the noise-free gap is <= 0.06, and
at `d = 6` there is essentially no trap at all, itself evidence that greedy is
rarely trapped.

| cell | reps per `delta` | `B` | code |
|---|---:|---:|---|
| binary tree / tree | 400 | 500 | pre-merge study code (Kronecker-`Sigma` search, tested path-for-path against `catci.search`) |
| step / ordinal | 300 | 500 | current |
| binary tree / all pairs | 300 | 300 | current |
| step / all pairs | 300 | 300 | current |
| greedy_trap / all pairs | 300 | 300 | current |
| greedy_trap beam cell | 200 | 300 | current |

## Results

Rejection rate at alpha = 0.05; the `delta = 0` column is size.

### Beam: greedy_trap / all pairs (delta = 5, 200 reps)

| method | rate | vs greedy |
|---|---:|---:|
| greedy | 0.595 | n/a |
| beam5 | 0.600 | +0.005 +- 0.011 |
| beam25 | 0.605 | +0.010 +- 0.010 |

beam5's size, from the next cell: 0.050 (greedy 0.053).

### greedy_trap / all pairs

| method | size | power d=4 | vs greedy | power d=5 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.373 | n/a | 0.673 | n/a |
| beam5 | 0.050 | 0.370 | -0.003 +- 0.007 | 0.670 | -0.003 +- 0.010 |
| random_fixed | 0.063 | 0.340 | -0.033 +- 0.023 | 0.547 | -0.127 +- 0.022 |
| random_fresh | 0.060 | 0.330 | -0.043 +- 0.022 | 0.560 | -0.113 +- 0.025 |
| split | 0.027 | 0.187 | -0.187 +- 0.027 | 0.327 | -0.347 +- 0.031 |
| split_one | 0.043 | 0.140 | -0.233 +- 0.027 | 0.200 | -0.473 +- 0.030 |
| depth0 | 0.030 | 0.400 | +0.027 +- 0.021 | 0.720 | +0.047 +- 0.019 |
| chi_sq | 0.027 | 0.413 | +0.040 +- 0.021 | 0.733 | +0.060 +- 0.019 |

The trap spreads signal over most cells, so the non-adaptive tests win here.

### binary_tree / tree structure

| method | size | power d=2 | vs greedy | power d=3 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.393 | n/a | 0.762 | n/a |
| random_fixed | 0.062 | 0.378 | -0.015 +- 0.013 | 0.752 | -0.010 +- 0.008 |
| random_fresh | 0.055 | 0.370 | -0.022 +- 0.013 | 0.743 | -0.020 +- 0.008 |
| split | 0.065 | 0.185 | -0.207 +- 0.023 | 0.435 | -0.328 +- 0.026 |
| split_one | 0.077 | 0.205 | -0.188 +- 0.025 | 0.412 | -0.350 +- 0.025 |
| depth0 | 0.033 | 0.117 | -0.275 +- 0.025 | 0.247 | -0.515 +- 0.026 |
| chi_sq | 0.033 | 0.115 | -0.278 +- 0.025 | 0.247 | -0.515 +- 0.026 |

### step / ordinal structure

| method | size | power d=3 | vs greedy | power d=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.397 | n/a | 0.713 | n/a |
| random_fixed | 0.040 | 0.263 | -0.133 +- 0.028 | 0.533 | -0.180 +- 0.028 |
| random_fresh | 0.043 | 0.297 | -0.100 +- 0.027 | 0.557 | -0.157 +- 0.030 |
| split | 0.057 | 0.207 | -0.190 +- 0.029 | 0.450 | -0.263 +- 0.029 |
| split_one | 0.050 | 0.237 | -0.160 +- 0.029 | 0.423 | -0.290 +- 0.029 |
| depth0 | 0.030 | 0.223 | -0.173 +- 0.025 | 0.467 | -0.247 +- 0.030 |
| chi_sq | 0.027 | 0.227 | -0.170 +- 0.026 | 0.467 | -0.247 +- 0.030 |

### binary_tree / all pairs

| method | size | power d=3 | vs greedy | power d=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.253 | n/a | 0.577 | n/a |
| random_fixed | 0.060 | 0.167 | -0.087 +- 0.025 | 0.323 | -0.253 +- 0.028 |
| random_fresh | 0.033 | 0.183 | -0.070 +- 0.021 | 0.337 | -0.240 +- 0.029 |
| split | 0.023 | 0.170 | -0.083 +- 0.023 | 0.303 | -0.273 +- 0.030 |
| split_one | 0.043 | 0.110 | -0.143 +- 0.024 | 0.300 | -0.277 +- 0.029 |
| depth0 | 0.030 | 0.223 | -0.030 +- 0.020 | 0.453 | -0.123 +- 0.024 |
| chi_sq | 0.027 | 0.237 | -0.017 +- 0.020 | 0.473 | -0.103 +- 0.024 |

### step / all pairs

| method | size | power d=3 | vs greedy | power d=4 | vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.257 | n/a | 0.533 | n/a |
| random_fixed | 0.060 | 0.160 | -0.097 +- 0.023 | 0.330 | -0.203 +- 0.030 |
| random_fresh | 0.033 | 0.187 | -0.070 +- 0.023 | 0.297 | -0.237 +- 0.028 |
| split | 0.023 | 0.157 | -0.100 +- 0.024 | 0.307 | -0.227 +- 0.028 |
| split_one | 0.043 | 0.127 | -0.130 +- 0.025 | 0.267 | -0.267 +- 0.029 |
| depth0 | 0.030 | 0.217 | -0.040 +- 0.018 | 0.460 | -0.073 +- 0.026 |
| chi_sq | 0.027 | 0.227 | -0.030 +- 0.018 | 0.467 | -0.067 +- 0.025 |

The `delta = 0` rows of the two non-trap all-pairs cells and the trap cell are the
same replicates (at `mu = 0` the DGP does not enter, and all used seed 0): one
size check, not three.

## Interpretation

**Greedy vs random.** A balanced binary tree admits many merge orders but they all
end at the root split, and the coarse levels, where this signal lives, are reached
by every path, so random matches greedy to within 1-2 points. Under ordinal or
all-pairs structures a random path usually coarsens across the signal (for the
step signal, cuts the ordering in the wrong place), and greedy's advantage is 7-13
points at the lower strength and 16-25 at the higher (11-13 on the diffuse
signal). The adaptive search does real work exactly when the structure does not do
it for it.

**Greedy vs sample splitting.** Splitting is valid here, but the test half carries
only half the noncentrality (`delta^2 / 2`). The cells are consistent with a split
test behaving roughly like a good fixed path at signal `delta / sqrt 2` (e.g. tree
structure: split at delta = 3, effective 2.1, gets 0.435, against greedy's 0.393 at
delta = 2). That is an interpretation, not something the design measures directly.
Greedy with minP uses all the data both to choose and to test, and pays for the
choice through the calibration, which is much cheaper than half the sample.
`split_one` avoids multiplicity entirely but is no better, often worse: the level
that looked best on half the data is a noisy choice.

**Beam.** Even where greedy is trappable (one narrow window of one construction),
a 25-wide beam gains nothing at the test level. This is evidence, not proof, that
no DGP rewards breadth.

**Non-adaptive anchors.** The oracle chi-square is the benchmark a practitioner
would otherwise use. Under a diffuse signal it wins by 4-6 points, as it should:
a signal spread over many cells has nothing to coarsen.

## Possible paper sentences

> Replacing the greedy search by a random merge path, drawn independently of the
> data, preserves validity but costs 7-25 percentage points of power whenever the
> merge structure leaves the final partition free; under a binary tree, where every
> path ends at the root split, the two agree to within 2 points.

> Choosing the path on one half of the sample and testing on the other is also
> valid, but loses 8-35 points relative to greedy with the minP calibration, and
> in most settings falls below the non-adaptive chi-square test: paying for
> selection through the calibration is far cheaper than paying with half the sample.

> Widening the greedy search to a beam of 5 or 25 candidate partition pairs per
> level does not improve power, even on a construction designed to defeat greedy.

## Caveats

* Gaussian limit with oracle `Sigma`, `d = 8`, uniform marginals, one seed per
  cell. The data-level registry entries (`<search>_random`, `<search>_split`,
  `greedy_beam5` in `experiments/methods.py`) were not run through the fitted
  propensity pipeline; there `_split` reuses full-sample propensities, so its halves
  are independent only up to that shared fit.
* The trap is one construction (see the narrow-window caveat above).

## Before the numbers go in the paper

Regenerate on current code at uniform settings:

* **Uniform reps and `B`** (e.g. 1000 and 1000, matching `experiments/config.py`).
* **Distinct seeds per cell** (`--seed`), so the null rows are independent.
* **Rerun binary tree / tree on current code**; the pre-merge code drew bootstrap
  samples through a different square root of `Sigma`, so the current code reproduces
  it in distribution, not draw for draw. The paper should cite one code version.
* **A `delta` grid** (0, 1, ..., 6) if this becomes a figure.
* **Optionally a data-level confirmation** via the registry entries above.

Uniform rerun (about 40k CPU-seconds: under two hours on six free cores):

```
cd <repo>; export PYTHONPATH=src:experiments
D="--deltas 0,1,2,3,4,5,6 --widths 1 --reps 1000 --n-boot 1000 --workers 6 --tag _paper"
python experiments/search_study.py binary_tree --structure tree      --seed 1 $D
python experiments/search_study.py step        --structure ordinal   --seed 2 $D
python experiments/search_study.py binary_tree --structure saturated --seed 3 $D
python experiments/search_study.py step        --structure saturated --seed 4 $D
python experiments/search_study.py greedy_trap --structure saturated --seed 5 $D
python experiments/summarise_search.py experiments/results/search_*_paper.parquet
```

Reproduce the current (exploratory) numbers, seed 0 throughout:

```
export PYTHONPATH=src:experiments
python experiments/search_study.py greedy_trap --structure saturated --deltas 5 --widths 1,5,25 --reps 200 --n-boot 300 --tag _beam
python experiments/search_study.py greedy_trap --structure saturated --deltas 0,4,5 --widths 1,5 --reps 300 --n-boot 300
python experiments/search_study.py step        --structure ordinal   --deltas 0,3,4 --widths 1 --reps 300 --n-boot 500
python experiments/search_study.py step        --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py binary_tree --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py binary_tree --structure tree      --deltas 0,2,3 --widths 1 --reps 400 --n-boot 500
python experiments/summarise_search.py experiments/results/search_*.parquet
```

## Code

See the docstrings in `src/catci/search.py` (`beam_search`, `random_merges`,
`evaluate_path`) and `adaptive_pvalue(..., search=)` in `src/catci/calibrate.py`;
tests in `tests/test_search_alternatives.py`.
