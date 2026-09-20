# Differential-test fixtures

`catci_fixtures.json` freezes the **deterministic seams** of the original R
method so the Python implementation can be tested against R as an oracle.

RNG streams will not match across languages, so we deliberately pin only the
input → output maps that are *pure functions of their arguments* — never the
bootstrap or randomised tie-break paths. Given the same inputs, the Python
implementation must reproduce every value here to numerical tolerance.

## Provenance

The R package has been removed from this repo. Its source, and the
`data-raw/fixtures/export_fixtures.R` script that produced this file, are
preserved on the **`r-frozen-oracle`** git tag:

```sh
git show r-frozen-oracle:data-raw/fixtures/export_fixtures.R
```

The exact R commit is recorded in `metadata.git_commit` inside the JSON. This
oracle is frozen: it pins behaviour the Python implementation must reproduce,
and there is no longer an R method for it to track.

## Conventions (read before parsing in Python)

- **Matrices** are objects `{nrow, ncol, data}` where `data` is **row-major**
  (`data[i]` is row `i`). Reconstruct with
  `np.array(m["data"]).reshape(m["nrow"], m["ncol"])`.
- **Every array is a JSON array**, even length 1 (so a singleton partition group
  `[3]` never collapses to a scalar `3`).
- **Label / dimension indices are 1-based** (as in R): `dimension` 1 = X, 2 = Y;
  merge pairs `[ind1, ind2]` and partition groups use 1-based original labels.
- **Index masks** (`index1`/`index2`) are 0/1 integer vectors of length `dx*dy`
  over the label order `(1,1),(2,1),…,(dx,1),(1,2),…,(dx,dy)` (X fastest).
- **Tolerance:** compare floats with `atol=1e-9, rtol=1e-7`. Verified: scipy
  `chi2.cdf` reproduces R `pchisq` on the `approx_chi` cases to < 1e-9.

## Sections

| key | seam | what it pins |
|---|---|---|
| `form_T_Sigma` | `(x, y, f, g, normalise) → (T_vector, Sigma)` | GCM construction; both `normalise` settings. Frozen `(x,y,f,g)` inputs included. |
| `search_paths` | `greedy_query` on a shared `(T, Sigma)`, `dx=dy=8` | statistic `values` **and** the partition sequence for `ordinal`, `tree`, `greedy`. Catches finding #1: `tree` ≠ `ordinal`. |
| `scalar_methods` | `query_lookup` depth-1 methods | `mGCM`, `max`, `euclid` on the shared `(T, Sigma)`. |
| `approx_chi` | Box (1954) chi-square CDF | `approx_chi_metric(normsq, tr, tr2) = pchisq(normsq/(tr2/tr), df=tr²/tr2)`. The one real cross-language numeric risk. |
| `rank_one_updates` | update formulae (24)–(27) | fast update **and** dense recomputation of `new_T`, `new_Sigma`, `normsq/tr/tr2`; they must agree. |
| `tree_structure` | the structure object | `make_binary_tree(d)` shape and the `permitted_merges` map (`get_ind2`/`get_num_levels`) for `ordinal`/`greedy`/`tree` at `d ∈ {2,4,8}`, on the initial partition. |

## Suggested Python usage

The `search_paths`, `rank_one_updates`, and `tree_structure` sections are the
load-bearing ones for the port. In particular, `tree_structure.permitted`
directly specifies the `permitted_merges(structure, partition)` operation the
redesign centres on — the Python `structure` module should reproduce those pair
lists exactly. Pair `rank_one_updates` (a fixed oracle) with a Hypothesis
property test (`fast == dense` for random `(T, Σ, i, j)`) for full coverage.
