# Divisive (reverse) label search — worktree notes

**Worktree:** `.claude/worktrees/divisive-search`
**Branch:** `worktree-divisive-search` (branched from `master` @ `196cc3a`)
**Status (2026-09-26):** rebased onto master (minP calibration, vectorised merge
search), divisive search vectorised, **power/size run at d = 8** -- see the
update section directly below. The rest of the file is the 2026-08-23 write-up;
where it says "not done", check the update first.
**Branch now:** `worktree-divisive-search-vec` (the old `worktree-divisive-search`
on origin predates the rebase).
**Date:** 2026-08-23, updated 2026-09-26

TODO item this addresses: *"Reverse the greedy procedure — start with fully
collapsed levels (ie just binary options), then expand outwards. Should be more
efficient in high dimensions?"*

```sh
conda activate catci
pip install -e ".[test,experiments]"
python -m pytest tests/test_divisive.py -q          # 31 passed (~7 min)
python experiments/bench_search.py --structure tree --dx 8 16 32 --levels 2 4 8
python experiments/bench_search.py --structure ordinal --dx 8 16 30 --dy 10
```

---

## Update 2026-09-26: vectorised, and the power question answered at d = 8

### Vectorisation

`divisive_search_paths(T, Sigma, ...)` runs every draw at once, mirroring
master's `greedy_search_paths`. Idea: of `(normsq, tr, tr2)` only `normsq`
depends on the draw; `tr, tr2` are functions of the partition, and all draws
share `Sigma`. So each partition any draw reaches is expanded once
(`_DivisiveTable.expand`: every candidate split scored in a few batched block
calls), and each draw only pays the `normsq - 2 sum(Ta Tb)` update, batched over
all draws at the same partition. The per-draw loop is kept as
`_divisive_search_loop`, the reference; `test_vectorised_matches_loop` pins
them (max abs diff ~7e-15, identical partitions). Full suite: 228 passed.

`MergeSearch` / `SplitSearch` now expose `paths(T, Sigma, ...)` (batched, what
`adaptive_pvalue(search=...)` uses) and `result(T_vector, ...)` (observed path
with partitions, what `catci_test(search=...)` uses); `prepare` is gone.

ms per draw, B = 1001 draws in one batch, tree, random Sigma:

| d | merge (master) | full split | split@2 | split@4 | split@8 | old loop split@4 |
|---|---:|---:|---:|---:|---:|---:|
| 8 | 2.3 | 2.3 | 0.015 | 0.33 | 1.6 | 21 |
| 16 | 23.7 | 64.6 | 0.063 | 0.25 | 6.9 | 15 |

(d = 32 and ordinal rows not measured -- the run was cut short.) Full split is
slower than the vectorised merge at d = 16: it does not carry state across
levels, and fine levels have little partition sharing between draws.

### Power and size, d = 8

`experiments/power_divisive.py`: oracle propensities, `lin` margins,
n = 1000, B = 1000, 200 reps, minP calibration; every method on the same data
and the same bootstrap draws. CSVs in `experiments/results_divisive/`.

**`binary_tree` (signal at the coarse levels):** rejection at 0.05

| strength | merge | split | split@2 | split@4 | split@8 |
|---|---:|---:|---:|---:|---:|
| 0 (size) | 0.060 | 0.055 | 0.045 | 0.050 | 0.060 |
| 0.2 | 0.065 | 0.080 | 0.105 | 0.085 | 0.080 |
| 0.4 | 0.280 | 0.270 | 0.350 | 0.365 | 0.325 |
| 0.6 | 0.520 | 0.520 | 0.665 | 0.645 | 0.585 |
| 0.8 | 0.820 | 0.820 | 0.890 | 0.865 | 0.850 |

Paired difference vs merge at 0.6: split@2 +0.145 +- 0.029, split@4
+0.125 +- 0.026, split@8 +0.065 +- 0.021, full split +0.000 +- 0.012.

**`alt` (signal only between sibling labels -- the finest level):**

| strength | merge | split | split@2 | split@4 | split@8 |
|---|---:|---:|---:|---:|---:|
| 0 (size) | 0.060 | 0.055 | 0.045 | 0.050 | 0.060 |
| 1.2 | 0.445 | 0.450 | 0.035 | 0.030 | 0.035 |
| 1.6 | 0.840 | 0.840 | 0.040 | 0.050 | 0.125 |
| 2.0 | 0.985 | 0.985 | 0.040 | 0.030 | 0.175 |

### Reading

* **Size is fine everywhere** (all within 1 SE of 0.05).
* **Direction alone changes nothing.** Full split equals merge to within
  +-0.015 in both DGPs. Hypothesis "divisive finds better coarse partitions" is
  *not* supported: whatever it finds, it does not show up as power.
* **Truncation is a bet on where the signal lives.** Coarse signal: it wins,
  up to +0.15, because dropping fine levels shrinks the minP multiplicity
  penalty -- and it is ~60x cheaper than merge. Fine signal: it has *no* power,
  because it never reaches the level the signal is on.
* So truncated divisive search is a sound *option* when the analyst believes
  the effect is coarse (a tree built from domain knowledge), not a replacement
  for the full search. A natural follow-up: minP over the union of a truncated
  split path and the merge path, to see whether it keeps most of both gains.

### d = 16 (binary_tree, n = 4000, B = 1000, 100 reps)

| strength | merge | split | split@2 | split@4 | split@8 |
|---|---:|---:|---:|---:|---:|
| 0.4 | 0.78 | 0.80 | 0.89 | 0.89 | 0.88 |
| 0.8 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

Paired vs merge at 0.4: split@2 / @4 +0.110 +- 0.031, split@8 +0.100 +- 0.030,
full split +0.020 +- 0.014. Same picture as d = 8. Seconds per test: merge 15.6,
full split 58.8, split@2 **0.06**, split@4 0.28, split@8 5.6.

Size: the 100-rep run showed 0.09-0.11 for the truncated splits at strength 0,
so it was re-run at 1000 reps (`..._size.csv`, SE 0.007): split@2 / @4 / @8
reject **0.051 / 0.048 / 0.052** at 0.05 (and 0.010 / 0.010 / 0.009 at 0.01,
0.103 / 0.097 / 0.101 at 0.10). Calibrated; the 100-rep excess was noise.

### Union of paths (d = 8, binary_tree) -- no gain

`merge+split@k` stacks the merge path and the split@k path into one minP (still
exact: a deterministic map of T). Same seed as the table above, so the base
methods reproduce it exactly:

| strength | merge | split@2 | merge+split@2 | merge+split@4 | merge+split@8 |
|---|---:|---:|---:|---:|---:|
| 0 (size) | 0.060 | 0.045 | 0.055 | 0.060 | 0.060 |
| 0.4 | 0.280 | 0.350 | 0.280 | 0.280 | 0.280 |
| 0.6 | 0.520 | 0.665 | 0.520 | 0.525 | 0.510 |
| 0.8 | 0.820 | 0.890 | 0.820 | 0.820 | 0.815 |

The union is merge, to the rep. It is *not* because the two searches visit the
same partitions: split@2's level-1 / level-2 partitions are also on the merge
path only ~53% / ~34% of the time (200 datasets, null and strength 0.6; level 0
is the tree root, shared by construction). The divisive search does find
different coarse partitions -- they just carry no extra signal once pooled with
merge's 13 levels. Truncation's gain is **entirely multiplicity**: calibrating
over 3 levels instead of 13. `alt` union not run; since the union contains the
whole merge path it should track merge there too (prediction, not measured).

Implication: the lever is not the search direction but *how the minP spreads
its budget over depths*. A weighted minP (more weight on coarse levels) run on
the ordinary merge path would test that directly and needs no divisive search.

End-to-end cost of this run: 237 s vs 1284 s for the same run before master's
vectorised merge landed.

### Still open

0. Weighted minP over depths on the merge path (see above) -- likely the real
   follow-up; belongs with the calibration work (TODO 0b / `MINP.md`).
1. d = 32 power; `alt` at d = 16 (expect the same collapse as d = 8).
2. Ordinal structure power (step interaction) -- not run.
3. `experiments/bench_search.py`: master's version was kept at the rebase; the
   divisive columns from the old benchmark were not ported.
4. Items 5 (memory) and 6 (non-contiguous structures) below are unchanged.

---

## Headline result

The hypothesis is **half right, and the half that works is worth having.**

Running the divisive search *to completion* is not faster than merging — it
visits the same `dx + dy - 3` levels and does comparable work at each. The win
is that the divisive direction lets you **truncate the fine end**, which
merging cannot do: merging's expensive levels come first, splitting's come last.

At `dX = dY = 32`, tree structure:

| | full merge | full split | split, 2 levels | split, 4 levels |
|---|---:|---:|---:|---:|
| ms per bootstrap draw | 188.8 | 199.3 | **1.25** | **2.48** |

A truncated divisive search is **~150x cheaper** than the full merge path and is
close to *constant in `d`* (0.93 → 0.96 → 1.25 ms as `d` goes 8 → 16 → 32),
because it never handles a partition with more than `L + 2` groups per
dimension. The only part that scales with `dx*dy` is the one-off block-sum
table, and that is shared across all `n_boot + 1` draws.

**Whether that truncation costs power is the open question — see "Not done".**

## The catch: it depends on the structure naming its own starting partition

| structure | full merge | full split | split@2 | why |
|---|---:|---:|---:|---|
| tree, `dX=dY=32` | 188.8 | 199.3 | **1.25** | root gives *one* two-group partition |
| ordinal, `dX=30, dY=10` | 45.5 | 123.8 | 26.8 | `(dX-1)(dY-1) = 261` starting partitions to score |

Merging bottoms out at a unique singleton partition; splitting has to *choose*
where to start. A tree names its own coarsest split, so there is exactly one
candidate. Two ordinal variables offer `(dX-1)(dY-1)`, and scoring them
dominates everything else — `split@2` at `dX=30` is 26.8 ms, of which almost all
is the start search.

This is **partly an implementation artefact**: the 261 starting partitions are
each a 4-element aggregation, evaluated in a Python loop with ~6 numpy round
trips apiece. Batching them into one block call should collapse this to
near-nothing. Not done — see below. The `O(dX dY)` candidate *count* is
intrinsic, but the per-candidate constant is not.

**`Saturated` has no divisive counterpart at all** and raises. Splitting a group
of size `s` admits `2^(s-1) - 1` bipartitions, so the top-down candidate set is
exponential where the bottom-up one is quadratic. This asymmetry is the reason
`divisive_search` is restricted to `Ordinal` and `Tree`.

---

## Full benchmark output

Both machines-shared caveat applies: a **peer Claude session was running the
test suite on the same machine** during earlier runs, and those numbers were
visibly contaminated (`d=32` full split read 337 ms, then 745 ms, then 199 ms
across three runs). The tables below are the quietest run obtained, `--repeats
5`, best-of. Treat ratios as reliable to maybe ±20%, absolute ms less so.
**Re-run before quoting anywhere.**

```
### structure=tree, amortised over n_boot=100 draws  (ms per draw)
  dx   dy    L      merge      split   ratio     table   split@2   split@4   split@8
   8    8   13       2.56       4.34    0.6x       0.1      0.93      2.02      3.65
  16   16   29      15.39      17.64    0.9x       1.0      0.96      2.25      6.15
  32   32   61     188.80     199.32    0.9x      29.3      1.25      2.48      5.38

### structure=ordinal, dy=10, amortised over n_boot=100 draws  (ms per draw)
  dx   dy    L      merge      split   ratio     table   split@2   split@4   split@8
   8   10   15       5.88      14.61    0.4x       0.1      7.19      9.42     13.23
  16   10   23      14.51      39.02    0.4x       0.2     14.24     18.55     26.69
  30   10   37      45.46     123.78    0.4x       0.6     26.80     33.78     47.71
```

`table` = one-off `SigmaBlocks` build (shared by all draws); the `split` columns
include it amortised over `n_boot + 1 = 101` paths.

---

## Why this is a valid test (theory)

Theorem 1 does **not** need the search to be agglomerative. Lemma 10 in the
paper is stated for *any* map `Pi : t -> Pi(t)` into `C*`, and Lemma 5's
continuity argument only needs the candidate set at each stage to be finite and
the selection to be an argmax. A divisive search is another deterministic
`(t, sigma) -> ` sequence-of-partitions map, so the double bootstrap calibrates
it for the same reason. Truncation is explicitly allowed too — the paper already
says "`L` is assumed to be such that we stop after the `L`th step, if at all".

This is checked empirically as well: `test_divisive_calibrates_under_gaussian`
covers ordinal, tree, and truncated tree.

---

## What was built

### `src/catci/structure.py` — split operations
Two new methods on `Structure`, alongside `permitted_merges`:

```python
structure.coarsest_partitions(d) -> [partition, ...]     # the two-group partitions to start from
structure.permitted_splits(partition) -> [(i, a, b), ...]  # replace group i by a then b
```

- `Ordinal`: `d - 1` starts; a group of size `s` has `s - 1` cut points.
- `Tree`: one start (root's children); a group splits into its node's children.
- `Saturated`: raises `NotImplementedError` with the reason.

A group's splits depend only on that group (not the rest of the partition) —
`_splits_of`'s memo in `search.py` relies on this, and
`test_a_groups_splits_do_not_depend_on_the_rest_of_the_partition` pins it.

### `src/catci/blocks.py` — new module
`SigmaBlocks` cumulative-sums `Sigma` along all four label axes once, after
which the sum over any `(X-range x Y-range) x (X-range x Y-range)` block is a
handful of lookups. Both supported structures keep groups contiguous in label
order (Ordinal by construction, Tree because `make_binary_tree` puts leaves
in order), so every group is a range and every block is a rectangle.

Two routes to the same answer, picked by which moves less memory:
`row_prefix` + `column_blocks` (8 gathers, but one stage is `(R, dy+1, dx+1)`
regardless of how few columns are wanted) vs `direct_blocks` (16 gathers, never
wider than `(R, C)`). Coarse partitions want the direct one.

The `Sigma` table is built once per test and reused across draws; the `T` prefix
is `O(dx*dy)` and rebuilt per draw.

### `src/catci/merging.py` — inverse formulae
`split_normsq` / `split_tr` / `split_tr2` invert (24)–(27). A merge reads the
current `(T, Sigma)`; a split reads only the two *new* rows of the finer
`(T, Sigma)`, which `blocks.py` supplies. **`divisive_search` never materialises
an aggregated `Sigma` at any level** — that is what makes truncation cheap.

### `src/catci/criteria.py` — `ApproxChi.split`
The counterpart of `ApproxChi.update`.

### `src/catci/search.py` — `divisive_search` + direction objects
`divisive_search(..., max_levels=None, sigma_blocks=None)`.

`MergeSearch` / `SplitSearch(max_levels)` are what `calibrate` and `api` are
handed. Each has `prepare(Sigma, dx, dy, xs, ys, criterion) -> path(T) -> SearchResult`,
so per-`Sigma` setup (the block table) is paid once per test rather than once
per draw.

### `calibrate.py` / `api.py`
`adaptive_pvalue(..., search=)` and `catci_test(..., search=)`, defaulting to
`MergeSearch()` — existing behaviour unchanged.

```python
from catci.search import SplitSearch
catci_test(x, y, Tree.binary(dx), Tree.binary(dy), f=f, g=g,
           search=SplitSearch(max_levels=4))
```

---

## Side change, unrelated to the search direction

`approx_chi_metric` now calls `scipy.special.gammainc(h/2, x/2)` instead of
`scipy.stats.chi2.cdf(x, df=h)`. These are **bit-identical** (verified over
20 000 random `(x, df)` pairs, max abs diff exactly 0.0, plus extremes), but the
`scipy.stats` route costs **82.5 us vs 1.00 us** per scalar call.

This was 31% of divisive search runtime before the swap and it would have masked
the comparison. It speeds up **both** directions equally, so the benchmark stays
fair — and it is a free ~1.3x on the existing merge path. All 56 pre-existing
tests pass unchanged.

**This is worth keeping regardless of what happens to the divisive work.**

---

## Test status

`tests/test_divisive.py` — **31 passed** (~7 min; the calibration tests dominate).

There is no R oracle for this path (it does not exist in the R package), so it
is pinned three other ways:

1. **Against dense recomputation.** Every level the search reports must equal
   `phi(Pi T, Pi Sigma Pi^T)` recomputed from scratch for the partition it says
   it is at. Hypothesis over `dx, dy in [2,9] x` {tree, ordinal}, plus a fixed
   check at `dX=30, dY=10` where a wrong stride would still look plausible.
   Max abs error observed: **2.1e-15**.
2. **Structurally against `greedy_search`.** Same lattice, opposite direction:
   same number of levels, agreeing at the singleton end (and, for a tree, at the
   two-group end too, since a tree has only one), while differing in between —
   `test_divisive_differs_from_merging` is the assertion that makes the whole
   comparison worth running.
3. **By calibration**, as for the merge direction.

Plus: truncation is a prefix of the full path; split-then-merge is the identity;
non-contiguous groups are rejected with a clear message; `Saturated` refuses.

Pre-existing suite: **56 passed** after the `gammainc` swap (verified before the
`search=` rewiring; re-run `pytest tests/ -q` to confirm the final state — it was
not re-run after the last edits).

---

## Not done — pick up here

1. **Re-run `pytest tests/ -q`** on the full suite. The 56 pre-existing tests
   passed after the `gammainc` change but were not re-run after `calibrate.py` /
   `api.py` gained the `search=` parameter. The public API was smoke-tested by
   hand (merge / split / split@3 all return sane p-values and level counts) but
   not under pytest.

2. **Power and size comparison — the actual research question.** Nothing here
   yet. The runtime story is only interesting if truncated divisive does not
   lose power. Three specific hypotheses worth testing, all with oracle
   propensities (`f=`, `g=`) so the search is isolated from the regression fit:
   - *Truncation should cost little power*, because power lives at the coarse
     end (`dX dY = 300` components at the finest level vs 4 at the coarsest),
     and dropping fine levels also shrinks the multiplicity penalty in the
     `max`-over-depths calibration — so it could even *gain*.
   - *Divisive should find better coarse partitions*, because it picks them by
     their own criterion rather than arriving via a chain of locally-greedy
     merges taken on noisy fine-level criteria.
   - Suggested grid: `experiments/dgp.py` with `intsetting="binary_tree"`, tree
     structure, `dX = dY in {8, 16, 32}`, strengths 0.2–1.8, ~200 reps, comparing
     `MergeSearch()` vs `SplitSearch(None)` vs `SplitSearch(2/4/8)` paired on
     identical data. Size run at strength 0 alongside.
   - Note `experiments/tuning/` only has hyperparameters for `n=1000, d=8`, which
     is why oracle propensities are the right call for a first pass.

3. **Vectorise the ordinal start search.** The `(dX-1)(dY-1)` starting
   partitions are scored in a Python loop, ~6 numpy calls each. Batching them
   into a single `SigmaBlocks.blocks` call should make ordinal `split@L`
   competitive with tree's. This is the single biggest implementation win left
   and it is what makes the ordinal row of the table look bad.

4. **`experiments/methods.py` registry entries** for the divisive variants, so
   `run.py` can drive them. Currently `_metric_fn` only knows `greedy_search`.

5. **Memory ceiling.** `SigmaBlocks` allocates `(dy+1)^2 (dx+1)^2` floats —
   ~9.5 MB at `d=32`, ~143 MB at `d=64` (vs `Sigma`'s own 134 MB, so ~2x, not
   catastrophic, but it caps the usable `d`). Not addressed.

6. **Non-contiguous structures.** Divisive search requires groups contiguous in
   label order. Fine for Ordinal and Tree as implemented, and fine under the
   label-permutation misspecification experiments (which permute the data, not
   the structure's positions) — but a structure with genuinely non-contiguous
   groups would need a mask-based fallback in `blocks.py`.

---

## Files

| file | change |
|---|---|
| `src/catci/blocks.py` | **new** — `SigmaBlocks`, `t_prefix`, `t_blocks` |
| `src/catci/structure.py` | `coarsest_partitions`, `permitted_splits` on Ordinal/Tree/Saturated |
| `src/catci/merging.py` | `split_normsq`, `split_tr`, `split_tr2` |
| `src/catci/criteria.py` | `ApproxChi.split`; `gammainc` swap |
| `src/catci/search.py` | `divisive_search`, `MergeSearch`, `SplitSearch` |
| `src/catci/calibrate.py` | `adaptive_pvalue(..., search=)` |
| `src/catci/api.py` | `catci_test(..., search=)` |
| `tests/test_divisive.py` | **new** — 31 tests |
| `experiments/bench_search.py` | **new** — the runtime benchmark |

Nothing is committed. `DIVISIVE.md` is gitignored (`*.md` with a `README.md`
negation), which is why these notes live here rather than in `README.md`.
