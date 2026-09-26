# Batched search — one contraction per level, not one per candidate

**Status:** steps 1 and 2 **landed** on branch `worktree-vectorise-search`
(2026-09-26) -- see [Outcome](#outcome-2026-09-26) at the end; the rest of this
document is the original 2026-09-02 proposal, kept for the derivation.
**Date:** 2026-09-02 (proposal), 2026-09-26 (implementation)
**TODO item:** 8.

Motivating complaint: `n_boot = 100` is too small, and the price of raising it is
linear in `B` because `calibrate.py` runs `B` sequential Python searches. The
finding below is that **the linear-in-`B` term is not the problem** — the problem
is a factor of 15–80x sitting inside each individual search, and once you take
it, you can spend it on `B`.

---

## Diagnosis: the search is overhead-bound, not FLOP-bound

`DIVISIVE.md` measured full merge paths at `dx = dy = 8, 16, 32` (tree): 2.56,
15.39, 188.80 ms per draw. That is an empirical exponent of **~3.1** in `d`.

Count the arithmetic instead. `Saturated` evaluates `~d^3/3` candidates over a
full path, and `update_tr2` ([merging.py:96](src/catci/merging.py#L96)) reads all
of `Sigma` for each one, at `O(p*d) = O(d^3)`. So the FLOP count is `O(d^6)` for
`Saturated`, `O(d^5)` for tree — far above the measured 3.1.

**The gap is the tell.** At the sizes tested, wall time is dominated by
per-candidate Python and numpy dispatch, not by arithmetic: ~1000 small numpy
calls per search, each doing very little work. That is the regime where
restructuring wins big and where extrapolating measured timings to `d = 100` is
untrustworthy in both directions.

There is a second, physical version of the same point. At `d = 32` a single
candidate's `cross` term reads `2 * p * dy * 8 = 524 kB`; with ~500 candidates at
the widest level that is **262 MB of memory traffic per level per draw**, to
produce 500 scalars. The batched form below reads `Sigma` **once** — 8.4 MB — for
the entire level.

---

## The identity: every candidate at a level, in six contractions

All four quantities the search tracks are bilinear in `Sigma` and `T`, so the
whole candidate set is a tensor contraction with the two merge labels left as
*free* indices instead of being fixed per call.

Let `p = dx*dy`, with the X-fastest ordering used throughout the package:

```python
M  = T.reshape(dy, dx).T            # (dx, dy) table, M[j, k]
S4 = Sigma.reshape(dy, dx, dy, dx)  # [y_row, x_row, y_col, x_col]
```

Then for **every** pair of X-labels `(j, m)` at once:

```python
def all_pairs(M, S4):
    """Deltas for every merge of two fast-axis labels.

    M  : (d_fast, d_slow) table
    S4 : (d_slow, d_fast, d_slow, d_fast), Sigma reshaped
    Returns (d_normsq, d_tr, d_tr2), each (d_fast, d_fast); upper triangle is
    the candidate set. Add to the current (normsq, tr, tr2) to get each
    candidate's state — i.e. exactly formulae (24)-(27), all pairs at once.
    """
    d_normsq = 2.0 * (M @ M.T)                                    # formula (24)
    d_tr     = 2.0 * np.einsum('ajam->jm', S4)                    # formula (25)
    cross    = np.einsum('abkj,abkm->jm', S4, S4)                 # (26): the dominant term
    D        = np.einsum('ajkj->jak', S4)
    term1    = np.einsum('jak,mak->jm', D, D)                     # (27), within-block
    term2    = np.einsum('ajkm,amkj->jm', S4, S4)                 # (27), cross-block
    return d_normsq, d_tr, 4.0 * cross + 2.0 * (term1 + term2)
```

**The Y branch is the same kernel on a transposed view** — no second
implementation:

```python
x_cands = all_pairs(M,   S4)
y_cands = all_pairs(M.T, np.ascontiguousarray(S4.transpose(1, 0, 3, 2)))
```

`cross` is the only term with a non-trivial cost: free `(j, m)`, contracted over
`(a, b, k)`, i.e. `O(dx^2 * p * dy) = O(d^5)` for the whole level. That is the
*same* FLOP count as the current per-candidate loop — this is not an asymptotic
improvement, it is the same arithmetic delivered as one GEMM-shaped contraction
rather than `d^2/2` dispatches. Reshaping `S4` to `(dx, p*dy)` makes `cross`
literally `A @ A.T`.

**Verified.** Every entry of all three returned matrices matches
`merging.update_{normsq,tr,tr2}` to `rtol=1e-10`, checked pairwise on both
dimensions at `(dx, dy) = (5,4), (4,7), (6,6), (3,9)`.

### Measured, numpy CPU, one full level, all candidates

| `dx=dy` | `p` | candidates | per-candidate loop | batched | speedup |
|---:|---:|---:|---:|---:|---:|
| 8 | 64 | 56 | 1.97 ms | 0.13 ms | **16x** |
| 16 | 256 | 240 | 11.60 ms | 0.68 ms | **17x** |
| 32 | 1024 | 992 | 157.22 ms | 10.32 ms | **15x** |
| 48 | 2304 | 2256 | 4931.17 ms | 60.70 ms | **81x** |

The jump at `d = 48` is the loop falling out of cache — precisely the traffic
argument above. No GPU, no draw batching, no torch: this is numpy on the same
machine.

---

## The `B` axis: honest answer

**Wall time cannot be made sub-linear in `B`.** The arithmetic is genuinely
linear in the number of draws and nothing restructures that away. Adding a
leading draw axis (`S4` becomes `(B, dy, dx, dy, dx)`) only helps while
*dispatch overhead* dominates, and that regime is narrow:

| `dx=dy` | per-draw cost at `B=1` | at `B=16` | at `B=256` | best gain |
|---:|---:|---:|---:|---:|
| 8 | 0.088 ms | 0.015 ms | 0.015 ms | **6.0x**, saturated by `B=16` |
| 16 | 0.495 ms | 0.320 ms | 0.383 ms | **1.6x**, saturated by `B=16`, then decays |

At `d = 16` the batched kernel is already arithmetic-bound at `B = 1`, so the
draw axis buys almost nothing and starts to *lose* at large `B` from memory
pressure.

**So the answer to "B is too small" is the constant factor, not sub-linearity.**
The 15–80x from candidate batching applies at every `B`. Spend it: `B = 100 ->
1000` at roughly one tenth of today's wall time. That matters for real reasons —
p-value resolution is floored at `1/(B+1)`, which is 0.0099 today and useless
under any multiplicity correction across the item 4 applications; and
`calibrate.py` reuses the same `B` draws for both bootstrap stages, so the
max-over-depth calibration is estimated coarsely at `B = 100`.

The draw axis has a different job on a GPU: it is how you get **occupancy**.
Kernel launch overhead is larger there and throughput is far higher, so the
overhead-dominated regime is much wider and batching draws is what fills the
device. It is worth building the draw axis in from the start for that reason,
not for CPU sub-linearity.

---

## torch / MPS, measured

Same kernel, `py313-torch` (torch 2.10, Apple silicon MPS), one level, X branch
only, draw axis included:

| `d` | `B` | torch CPU f64 | torch CPU f32 | MPS f32 | MPS vs CPU f64 |
|---:|---:|---:|---:|---:|---:|
| 16 | 16 | 12.3 ms | 4.9 ms | 3.9 ms | 3.2x |
| 16 | 64 | 28.0 ms | 16.7 ms | 5.8 ms | 4.8x |
| 16 | 256 | 97.6 ms | 63.8 ms | 23.0 ms | 4.2x |
| 32 | 16 | 134.8 ms | 84.5 ms | 21.1 ms | 6.4x |
| 32 | 64 | 548.4 ms | 335.0 ms | 84.0 ms | 6.5x |
| 32 | 256 | 20107.8 ms | 1300.5 ms | 376.8 ms | 53.4x |

Three findings, and they argue for doing step 1 in numpy and stopping there
until proven otherwise:

* **MPS has no float64 at all** — `Cannot convert a MPS Tensor to float64 dtype`.
  Not a tuning choice; on this machine the GPU path is fp32 or nothing.
* **MPS buys ~4–6x over CPU** in the normal regime. The 53x row is CPU fp64
  falling apart at 20 s (memory pressure at 256 x 8.4 MB), not GPU brilliance.
* **torch is not free.** torch CPU f64 at `d=16, B=16` is 12.3 ms against numpy's
  5.1 ms for the same contraction — torch's einsum picks a worse path here. Do
  not assume porting to torch is a speedup on its own; it is a speedup only once
  you are on the device, in fp32.

Net: the ~15–80x from candidate batching is the real prize and it is available in
numpy today. The GPU adds maybe another 5x, costs you fp64, and costs you the
exact fixtures.

---

## Memory budget (this is what caps the design)

`(B, p, p)` float64 = `B * p^2 * 8` bytes, and the Y branch wants a transposed
contiguous copy, so budget **2x**:

| `dx=dy` | `p` | per draw | `B=100` | `B=1000` |
|---:|---:|---:|---:|---:|
| 8 | 64 | 65 kB | 6.5 MB | 65 MB |
| 16 | 256 | 1.0 MB | 105 MB | 1.0 GB |
| 32 | 1024 | 17 MB | 1.7 GB | 17 GB |
| 64 | 4096 | 268 MB | 27 GB | — |
| 100 | 10000 | 1.6 GB | — | — |

Chunk over draws beyond `d ~ 32`; that restores linearity in `B` but with the
good constant. Beyond `d ~ 64` the per-draw `Sigma` itself is the wall, and the
fix is not batching but the factored form — `Sigma` is never needed explicitly,
since each row of `prod_mat` ([gcm.py:73](src/catci/gcm.py#L73)) is
`outer(Xres[i], Yres[i])`, rank one, giving

```
tr(Sigma_pi) = n/(n-1) * [ mean_i ( ||Ax @ Xres[i]||^2 * ||Ay @ Yres[i]||^2 )
                           - ||Ax @ Mbar @ Ay.T||_F^2 ]
```

with no `p`-by-`p` object at all (`Ax`, `Ay` = 0/1 partition membership
matrices, `Mbar` = mean of the outer products). The `tr2` analogue is an
`n`-by-`n` Gram matrix at `O(n^2 d)`, so pick whichever of `n` or `p` is
smaller. **Not yet verified** — derive it properly before relying on it. This is
the piece that decides whether item 4c (ACS, `dx`, `dy` in the hundreds) can run
at all.

---

## The forward-pass framing

Worth recording, because it is the right mental model for the implementation and
it unifies items 7 and 8.

One level of `greedy_search` is: compute all candidate values (a batched
bilinear form — the kernel above), then take the max. That is **a linear layer
followed by a max-pool**, and the criterion path `result.values` is exactly the
sequence of pool outputs. The whole search is a forward pass whose intermediate
activations *are* the statistics of interest — a multi-exit network, with the
double bootstrap's max-over-depths being a max over layer outputs.

Two caveats keep the analogy honest:

* The weights are not learned; the argmax **routes** which aggregation is
  applied next, so it is a dynamic-routing / hard-attention network, not
  `nn.Sequential`. Routing diverges per draw, which is the one real
  implementation subtlety in the batched version (per-draw `gather`, not a
  shared weight matrix).
* The `T` branch is linear but the `Sigma` branch is a congruence,
  `Sigma -> P Sigma P.T`. Both are contractions; neither is an `nn.Linear`.
* Nothing is trained and there is no loss. This is a framing that buys
  **tooling** — batching, autodiff, `torch.compile`, GPU — not a learning claim.
  Say so explicitly if it reaches the paper, or a reviewer will ask what is being
  fit.

What it buys conceptually: **the pooling operator is the single knob that
unifies every search variant under discussion.**

| pooling over the candidate axis | search |
|---|---|
| `max` | greedy — today's method |
| `top-k` | beam search, width `k` — TODO item 7 |
| `softmax` at temperature `tau` | the relaxed / "search softer" variant |

Which means one implementation covers the whole item 7 power study, and the
`tau > 0` case is differentiable, so it comes with gradients for free. Note that
max-pool is already differentiable a.e. (subgradient through the argmax), so
even the greedy path is autodiff-able if that ever becomes useful.

---

## Implementation sketch

1. **Land the kernel in numpy first, behind the existing API.** `all_pairs` +
   transposed view, replacing the inner double loop of
   [search.py:71-82](src/catci/search.py#L71-L82). Single-draw, no batch axis.
   The `search_paths` fixture must still pass **bit-exactly** — same candidate
   order, same first-max tie-breaking (`np.argmax` on the flattened upper
   triangle in dimension-then-position order reproduces R's `which.max`).
   This alone is the 15–80x.
2. **Add the draw axis**, chunked, with the chunk size derived from the memory
   table. `calibrate.adaptive_pvalue` loses its Python loop.
3. **Then** consider torch/GPU, only if 1–2 leave `d` in the hundreds out of
   reach.
4. Re-run `experiments/bench_search.py` from the divisive worktree for a
   like-for-like comparison against the numbers in `DIVISIVE.md`.

Steps 1 and 2 are pure refactors with an exact oracle to check against, and they
are worth doing regardless of what happens with items 7 and 5.

---

## Risks / what this does not change

* **fp32 breaks bit-exactness.** Any GPU path will want fp32, and on MPS fp64 is
  not available at all (confirmed above). fp32 changes tie-breaking in the argmax
  and therefore the selected path. A GPU path needs tolerance-based differential
  tests against the fp64 CPU path, not the existing exact fixtures. Keep the
  observed statistic in fp64 regardless.
* **`optimize=True` on `np.einsum` is load-bearing** and its contraction-order
  choice can change with numpy versions. Pin the path explicitly
  (`np.einsum_path`) if timings need to be reproducible across environments.
* **No power is gained.** This is wall time only. It does not touch the null
  inflation argument in item 7, and it does not make the search find better
  partitions — it makes the *existing* search cheap enough that item 7's power
  study, item 5's power/size comparison, and item 1's `d` grid become affordable.
* **`bootstrap.matrix_sqrt` is `O(p^3)` once per test** and untouched by any of
  this. At `d = 100` that eigendecomposition is itself minutes on CPU. Check
  where the balance lands before assuming the search is still the bottleneck.

---

## Outcome (2026-09-26)

Landed on branch `worktree-vectorise-search`: `catci.search.greedy_search_paths(T
(p, B), Sigma, ...) -> (L+1, B)`, used by `calibrate.adaptive_pvalue`,
`experiments/methods.py` and (via `n_jobs`) `api.catci_test`. `greedy_search` is its
`B = 1` wrapper. The old loop survives as `_greedy_search_loop`, the oracle for
`tests/test_search_vectorised.py` (every structure, both carrying modes, every
compaction schedule, unequal `dx, dy`, low-rank `Sigma`, n-ary trees, batch vs
column-wise) and the fallback for a non-`ApproxChi` statistic. The R `search_paths`
fixture passes with identical partitions; `experiments/methods.adaptive_pvalues`
gives identical p-values to the old code for a given seed.

### What changed relative to the proposal above

* **Carried, not recomputed.** The proposal recomputed `cross` each level at
  `O(p^2 d)`. Instead `C`, `T1`, `T2`, `d_tr` are carried as `(d, d)` state and
  updated after each merge from the two merged rows, using `Sigma'^2 = P(Sigma^2 +
  Sigma E Sigma)P^T` with `P^T P = I + E`. The other-axis updates of `T1`/`T2` are the
  closed-form expansions of `tr(B(I+E)B(I+E))`. All verified to `rtol=1e-10`.
* **Two carrying modes, chosen per dimension.** *Dense* (full matrices, one
  `O(p^2)`-flop GEMM per merge) for `Saturated`. *Sparse* for trees and ordinal:
  only currently permitted pairs are tracked, and every pair a merge newly permits
  contains the merged group, so it is computed fresh from two rows. Level cost
  `O(p d)` per draw.
* **Padded slots, not shrinking arrays.** A group lives at its smallest label's
  slot, so every draw keeps one shape whatever its path, and slot order is
  position order (first max over row-major upper triangles, X then Y = the R
  tie-break). Absorbed labels are compacted away at 60% occupancy. Structures
  gained a vectorised `permitted_mask(sizes, gid)` (binary trees via a
  `(min label, size)` sibling table; n-ary trees fall back to `permitted_merges`).
* **Draw axis from the start** (step 2), chunked at 32 MB of per-draw `Sigma`,
  which the chunk sweep found best at every `d` tried (8-32): larger chunks spill
  the per-merge temporaries out of cache. `n_jobs` threads over chunks.

### Measured

CPU time per draw, single-threaded BLAS, **on a heavily loaded machine** (load
average 130-370 on 6 cores throughout). Absolute numbers are inflated ~2x relative
to `DIVISIVE.md`'s (its loop at `d = 32` was 189 ms, here 216-390 ms); the ratios
are the meaningful part. `B = 1000`, `experiments/bench_search.py --cpu`.

| structure | `dx x dy` | `p` | loop ms/draw | batched ms/draw | speedup |
|---|---|---:|---:|---:|---:|
| tree | 8 x 8 | 64 | 5.9 | 0.27 | **22x** |
| tree | 16 x 16 | 256 | 25.4 | 3.05 | **8x** |
| tree | 24 x 24 | 576 | 93.6 | 14.0 | **7x** |
| tree | 32 x 32 | 1024 | 216 | 39.2 | **6x** |
| tree | 30 x 10 | 300 | 42.1 | 4.43 | **9x** |
| ordinal | 8 x 8 | 64 | 5.7 | 0.35 | **16x** |
| ordinal | 16 x 16 | 256 | 36.5 | 4.12 | **9x** |
| ordinal | 32 x 32 | 1024 | 391 | 53.4 | **7x** |
| greedy | 8 x 8 | 64 | 19.8 | 0.34 | **58x** |
| greedy | 16 x 16 | 256 | 194 | 4.30 | **45x** |

`B = 10000`, tree: 2.9 s at `8 x 8`, 32 s at `16 x 16` (CPU, one thread).
Threads: `n_jobs = 6` cut wall time 3.1x (`16 x 16`) and 3.3x (`32 x 32`) at
`B = 600`, even under that load. So at `p = 1024`, `B = 10000` is roughly 2-3 min of
wall time on this 6-core machine (idle), against ~35 min for the loop.

### Where the time goes now, and the next levers

* Per merge, each draw reads two merged rows of `Sigma` (`p d` doubles) and writes
  row `u` and column `u` back, zeroing row/column `v`. **The column writes
  dominate at large `d`.** In X-major storage a Y column is a scatter of single
  elements, ~130 us per draw at `d = 32` whether done by fancy indexing or a
  per-draw loop. Removing it needs either a lazy-column scheme (rows
  authoritative, columns aggregated on read; clean only for interval partitions,
  i.e. in-order trees and ordinal) or a compiled kernel (numba, fused and
  `prange` over draws). A compiled kernel is the obvious next step if `p` in the
  thousands is needed; it would also remove the ~15 temporaries per merge.
* **Tree and ordinal gain less than greedy** because the loop was already cheap
  for them (`O(d)` candidates per level); the batched cost there is dominated by
  the `O(p d)` row traffic that any per-draw-`Sigma` scheme must pay.
* **Sharing across draws does not pay at large `d`.** `Sigma`-side state depends
  only on the partition, but among 512 null draws at `d = 16` the number of
  distinct partitions is ~500 by level 4 (tree, ordinal and greedy alike). At
  `d = 8` tree it peaks at 116/512, so sharing would help only in the regime that
  is already cheap.
* **Memory**: `8 p^2` bytes per draw in flight (8 MB at `p = 1024`, 134 MB at
  `p = 4096`), chunked, so the draw count is unconstrained. Beyond `p ~ 4096` the
  factored form in the memory section above is still the only route; it is
  unverified.
* `bootstrap.matrix_sqrt` (`O(p^3)`, once per test) is still untouched.
