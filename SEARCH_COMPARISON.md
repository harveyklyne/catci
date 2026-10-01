# Greedy vs random vs sample-split search: setup and results

Paper-facing write-up of the comparison in TODO items 7a, 7b and 4.1: does the
adaptive (greedy) search earn its keep, compared with a data-independent random
search and with choosing the search path on an independent half of the data?
The beam-width question (7c) is in `SEARCH.md`.

**Status of these numbers:** exploratory, not final. Cells differ in reps and
bootstrap size, one cell ran on older (equivalent) code, and two cells share their
null replicates — see [Before the numbers go in the paper](#before-the-numbers-go-in-the-paper).
The qualitative conclusions are stable at the reported paired standard errors.

## Summary

1. **Every method holds its size** (all within ~2 SE of α = 0.05).
2. **Greedy beats random search by 7–25 points of power** on the structured signals
   wherever the structure leaves the final partition free (ordinal structure,
   all-pairs structure; 3–13 on the diffuse signal). Under a
   binary tree structure the gap is 1–2 points: every merge path in a binary tree
   ends at the same root split, so there the structure, not the search, supplies
   the adaptivity.
3. **Sample splitting loses 8–35 points** against greedy (up to 47 when only the
   single best-looking level is tested), and in most cells falls below even the
   non-adaptive χ² test. Choosing the partition on half the data costs half the
   noncentrality in the test half; the minP calibration pays for selection more
   cheaply.
4. **Fixing one random path across the bootstrap draws vs a fresh path per draw
   makes no consistent difference** (differences within ±3.5 points, changing sign
   between cells).
5. With all label pairs permitted, greedy's advantage over the non-adaptive χ²
   test is modest (2–10 points on structured signals) and reverses on a diffuse
   signal. The large adaptive gains are under the structured (tree, ordinal) searches.

## Setup

### The Gaussian limit

The study works directly with the limiting law of the statistic rather than with
simulated `(X, Y, Z)` data:

* `T ~ N(mu, Sigma)`, `T` of length `dx * dy`, `dx = dy = d = 8`, with `Sigma` known.
* `Sigma = C (x) C`, `C = diag(p) - p p^T`, `p` uniform on the `d` labels — the null
  covariance of the residual products `(e_X - p)(e_Y - q)` when there is no `Z`.
  `rank(Sigma) = (d - 1)^2 = 49`.
* `mu = delta * u / ||u||_{Sigma^+}`, so the noncentrality is `mu' Sigma^+ mu = delta^2`
  for every signal shape `u` — power is comparable across shapes at equal `delta`.
  `u` is the vectorised (X-fastest), double-centred interaction matrix of the shape.

This is what the search and calibration see asymptotically. It removes the
propensity learner (and its tuning) from a comparison that is about the *search*,
and it makes sample splitting exact: the two halves give independent
`T_A, T_B ~ N(mu / sqrt 2, Sigma)`, and the full-sample statistic is
`T = (T_A + T_B) / sqrt 2`.

### Criterion and calibration

For a pair of partitions `pi = (pi_X, pi_Y)`, let `M_pi` sum the cells of each
block. The criterion is Box's approximate chi-square CDF of the merged statistic,

    S_pi(T) = F_Box( ||M_pi T||^2 ; tr(Sigma_pi), tr(Sigma_pi^2) ),   Sigma_pi = M_pi Sigma M_pi^T.

A search path `pi_0, ..., pi_L` starts at the finest partition and merges one pair
of groups per level until both variables have two groups: `L = 2d - 4 = 12`
merges, so `L + 1 = 13` statistics `S_l = S_{pi_l}(T)`.

Every method is calibrated by the same minP double bootstrap: `B` draws
`Z_b ~ N(0, Sigma)`, the same path construction applied to each, stage 1
prepivots each level against the other paths, and stage 2 calibrates the minimum
over levels. α = 0.05 throughout.

### Methods

| name | path for the observed `T` | path for draw `Z_b` | validity |
|---|---|---|---|
| `greedy` | greedy on `T`: each level takes the permitted merge maximising `S` | greedy on `Z_b` | same deterministic map for every draw |
| `random_fixed` | `pi` drawn at random, independently of the data: each level merges a uniformly chosen permitted pair | the same `pi` | conditional on `pi`, a fixed map |
| `random_fresh` | `pi^(0)` drawn at random | independent `pi^(b)` | `(T, pi^(0))`, `(Z_b, pi^(b))` i.i.d. under the null |
| `split` | greedy on `T_A`; the whole path scored on `T_B` | the same path, scored on `Z_b` | path is a function of `T_A`, independent of `T_B` and `Z_b` |
| `split_one` | as `split`, but only level `l* = argmax_l S_{pi_l}(T_A)` is tested on `T_B` (`L = 1`) | level `l*` of the same path | as `split` |
| `depth0` | finest partition only (`L = 1`) | same | non-adaptive anchor |
| `chi_sq` | `T' Sigma^+ T` against `chi^2_{(d-1)^2}` | — | oracle non-adaptive anchor |

Every test is exact at finite `B` under the null (the Monte Carlo test argument:
the observed and bootstrap paths are exchangeable). "Permitted" is set by the
structure:

* **tree** — sibling merges in a balanced binary tree over the labels in order;
* **ordinal** — adjacent groups only;
* **all pairs** — any two groups.

### Signals

| name | shape `u` | structure(s) |
|---|---|---|
| binary tree | hierarchical interaction on a balanced binary tree: a coarse 2×2 split, finer levels weighted `0.1^(level-1)`, the finest level 0 (`dgp.binary_tree_interaction`) | tree; all pairs |
| step | ±1 blocks: X split into low / high halves, the sign pattern flipping across the quarters of Y (`dgp.get_int("step")`) | ordinal; all pairs |
| diffuse | a coarse 2×2 checkerboard plus a stronger fine signal on row pairs straddling it (`search_study.trap_interaction`, built for the beam question) | all pairs |

Each structure is paired with the signal it is designed for, and each signal is
also run under all pairs — the structure-agnostic search.

### Replicates and standard errors

Within a replicate every method sees the same `T_A`, `T_B` and bootstrap draws,
so differences from greedy are **paired**. `Δ vs greedy` below is the mean of the
per-replicate difference of reject indicators, `± ` its standard error — much
tighter than the difference of two independent rates. Rate SEs are `sqrt(r(1-r)/reps)`:
≤ 0.029 at 300 reps, ≤ 0.025 at 400.

| cell | reps per `delta` | `B` (bootstrap draws) | code |
|---|---:|---:|---|
| binary tree / tree | 400 | 500 | pre-merge study code (see below) |
| step / ordinal | 300 | 500 | current |
| binary tree / all pairs | 300 | 300 | current |
| step / all pairs | 300 | 300 | current |
| diffuse / all pairs | 300 | 300 | current |

## Results

Rejection rate at α = 0.05. The `δ = 0` column is size.

#### Binary-tree signal, tree structure

| method | size (δ=0) | power δ=2 | Δ vs greedy | power δ=3 | Δ vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.393 | — | 0.762 | — |
| random_fixed | 0.062 | 0.378 | −0.015 ± 0.013 | 0.752 | −0.010 ± 0.008 |
| random_fresh | 0.055 | 0.370 | −0.022 ± 0.013 | 0.743 | −0.020 ± 0.008 |
| split | 0.065 | 0.185 | −0.207 ± 0.023 | 0.435 | −0.328 ± 0.026 |
| split_one | 0.077 | 0.205 | −0.188 ± 0.025 | 0.412 | −0.350 ± 0.025 |
| depth0 | 0.033 | 0.117 | −0.275 ± 0.025 | 0.247 | −0.515 ± 0.026 |
| chi_sq | 0.033 | 0.115 | −0.278 ± 0.025 | 0.247 | −0.515 ± 0.026 |

#### Step signal, ordinal structure

| method | size (δ=0) | power δ=3 | Δ vs greedy | power δ=4 | Δ vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.050 | 0.397 | — | 0.713 | — |
| random_fixed | 0.040 | 0.263 | −0.133 ± 0.028 | 0.533 | −0.180 ± 0.028 |
| random_fresh | 0.043 | 0.297 | −0.100 ± 0.027 | 0.557 | −0.157 ± 0.030 |
| split | 0.057 | 0.207 | −0.190 ± 0.029 | 0.450 | −0.263 ± 0.029 |
| split_one | 0.050 | 0.237 | −0.160 ± 0.029 | 0.423 | −0.290 ± 0.029 |
| depth0 | 0.030 | 0.223 | −0.173 ± 0.025 | 0.467 | −0.247 ± 0.030 |
| chi_sq | 0.027 | 0.227 | −0.170 ± 0.026 | 0.467 | −0.247 ± 0.030 |

#### Binary-tree signal, all pairs

| method | size (δ=0) | power δ=3 | Δ vs greedy | power δ=4 | Δ vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.253 | — | 0.577 | — |
| random_fixed | 0.060 | 0.167 | −0.087 ± 0.025 | 0.323 | −0.253 ± 0.028 |
| random_fresh | 0.033 | 0.183 | −0.070 ± 0.021 | 0.337 | −0.240 ± 0.029 |
| split | 0.023 | 0.170 | −0.083 ± 0.023 | 0.303 | −0.273 ± 0.030 |
| split_one | 0.043 | 0.110 | −0.143 ± 0.024 | 0.300 | −0.277 ± 0.029 |
| depth0 | 0.030 | 0.223 | −0.030 ± 0.020 | 0.453 | −0.123 ± 0.024 |
| chi_sq | 0.027 | 0.237 | −0.017 ± 0.020 | 0.473 | −0.103 ± 0.024 |

#### Step signal, all pairs

| method | size (δ=0) | power δ=3 | Δ vs greedy | power δ=4 | Δ vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.257 | — | 0.533 | — |
| random_fixed | 0.060 | 0.160 | −0.097 ± 0.023 | 0.330 | −0.203 ± 0.030 |
| random_fresh | 0.033 | 0.187 | −0.070 ± 0.023 | 0.297 | −0.237 ± 0.028 |
| split | 0.023 | 0.157 | −0.100 ± 0.024 | 0.307 | −0.227 ± 0.028 |
| split_one | 0.043 | 0.127 | −0.130 ± 0.025 | 0.267 | −0.267 ± 0.029 |
| depth0 | 0.030 | 0.217 | −0.040 ± 0.018 | 0.460 | −0.073 ± 0.026 |
| chi_sq | 0.027 | 0.227 | −0.030 ± 0.018 | 0.467 | −0.067 ± 0.025 |

#### Diffuse signal, all pairs

| method | size (δ=0) | power δ=4 | Δ vs greedy | power δ=5 | Δ vs greedy |
|---|---:|---:|---:|---:|---:|
| greedy | 0.053 | 0.373 | — | 0.673 | — |
| random_fixed | 0.063 | 0.340 | −0.033 ± 0.023 | 0.547 | −0.127 ± 0.022 |
| random_fresh | 0.060 | 0.330 | −0.043 ± 0.022 | 0.560 | −0.113 ± 0.025 |
| split | 0.027 | 0.187 | −0.187 ± 0.027 | 0.327 | −0.347 ± 0.031 |
| split_one | 0.043 | 0.140 | −0.233 ± 0.027 | 0.200 | −0.473 ± 0.030 |
| depth0 | 0.030 | 0.400 | +0.027 ± 0.021 | 0.720 | +0.047 ± 0.019 |
| chi_sq | 0.027 | 0.413 | +0.040 ± 0.021 | 0.733 | +0.060 ± 0.019 |

## Interpretation

**Greedy vs random.** Random search asks whether the data-driven choice of merges
matters, or whether *any* coarsening path would do once the calibration is
exact. The answer depends on how much freedom the structure leaves. A balanced
binary tree admits many merge orders but they all end at the root split, and the
coarse levels — where this signal lives — are reached by every path; so random
matches greedy to within 1–2 points. Under ordinal or all-pairs structures a
random path usually coarsens across the signal (for the step signal, cuts the
ordering in the wrong place), and greedy's advantage is 7–13 points at the lower
strength and 16–25 at the higher (11–13 on the diffuse signal). This is the sanity check of TODO 4.1: the adaptive search is
doing real work exactly when the structure does not do it for it.

**Greedy vs sample splitting.** Splitting is the textbook way to make a
data-chosen partition valid, and it is valid here. But the test half carries only
half the noncentrality (`delta^2 / 2`). The cells are consistent with a
split test behaving roughly like a good fixed path at signal `delta / sqrt 2`
(e.g. tree structure: split at δ = 3, effective 2.1, gets 0.435, against
greedy's 0.393 at δ = 2). That is an interpretation, not something the design
measures directly. Greedy with the minP calibration uses all of the data both to
choose and to test, and pays for the choice through the calibration — which
these results say is much cheaper than giving up half the sample.
`split_one` — testing only the level that looked best on half A, which avoids any
multiplicity at all — is no better than testing the whole path, and often worse:
the level that looked best on half the data is a noisy choice.

**Fixed vs fresh random paths.** One might expect a single path shared by all
draws to calibrate more sharply than a new path per draw. In these cells the two
differ by at most 3.5 points with no consistent sign, so the choice does not
matter for power; fixed is simpler to state.

**Non-adaptive anchors.** The oracle χ² test is the benchmark a practitioner
would otherwise use. Under structured searches greedy beats it by 17–52 points;
under all pairs by 2–10 points on structured signals; on a diffuse signal it wins
by 4–6 points, as it should — a signal spread over many cells has nothing to
coarsen.

### Possible paper sentences

> Replacing the greedy search by a random merge path, drawn independently of the
> data, preserves validity but costs 7–25 percentage points of power whenever the
> merge structure leaves the final partition free; under a binary tree, where every
> path ends at the root split, the two agree to within 2 points.

> Choosing the path on one half of the sample and testing on the other is also
> valid, but loses 8–35 points relative to greedy with the minP calibration, and
> in most settings falls below the non-adaptive χ² test: paying for selection
> through the calibration is far cheaper than paying with half the sample.

## Before the numbers go in the paper

These results were produced on a heavily loaded machine and the design was cut to
fit. For the paper they should be regenerated on the current code at uniform settings:

* **Uniform reps and `B`**: e.g. 1000 reps per `delta` and `B = 1000`, matching
  `experiments/config.py`. (Current cells: 300–400 reps, `B = 300–500`.)
* **Distinct seeds per cell.** The two all-pairs cells used seed 0 with the same
  `Sigma`, so their `δ = 0` rows are the same null replicates — one size check,
  not two. Pass `--seed` per cell.
* **Rerun the binary tree / tree cell on current code.** It ran on the study's
  earlier Kronecker-`Sigma` search, which was tested path-for-path against
  `catci.search` — so the numbers should reproduce, but the paper should cite one
  code version.
* **A `delta` grid for power curves** (e.g. 0, 1, …, 6) if the comparison becomes a
  figure rather than a table.
* **Optionally, a data-level confirmation**: the registry entries
  `<search>_random` and `<search>_split` in `experiments/methods.py` run the same
  methods on simulated data with fitted propensities, but have not been run at
  scale. There, `_split` reuses the full-sample propensities, so its halves are
  independent only up to that shared fit.

A uniform rerun (≈ 40k CPU-seconds: under two hours on six free cores, much
longer on a loaded machine):

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

## Reproduce the current numbers

```
export PYTHONPATH=src:experiments
python experiments/search_study.py binary_tree --structure tree      --deltas 0,2,3 --widths 1 --reps 400 --n-boot 500
python experiments/search_study.py step        --structure ordinal   --deltas 0,3,4 --widths 1 --reps 300 --n-boot 500
python experiments/search_study.py binary_tree --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py step        --structure saturated --deltas 0,3,4 --widths 1 --reps 300 --n-boot 300
python experiments/search_study.py greedy_trap --structure saturated --deltas 0,4,5 --widths 1,5 --reps 300 --n-boot 300
python experiments/summarise_search.py experiments/results/search_*.parquet
```

(Seed 0 throughout. The binary tree / tree numbers above came from the pre-merge
code, which drew the bootstrap samples through a different square root of
`Sigma`; the current code reproduces them in distribution, not draw for draw.)
