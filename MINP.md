# Algorithm 2, rewritten as a minP test

Working note for the paper rewrite (TODO 0b). Replaces the old `M` / `F` / `G`
notation and fixes the calibration defect at the same time.

## What changed and why

The old Algorithm 2 computed, for each search level `l`, the empirical CDF value
`F_{n,l}` of the observed criterion among the bootstrap criteria, took
`G_n = max_l F_{n,l}`, and calibrated `G_n` against bootstrap copies `G_n^(b)`.

Two problems:

1. **It was not exchangeable.** `F_{n,l}` ranked the observed value against `B`
   bootstrap values *excluding itself*, while `F_{n,l}^(b')` ranked each bootstrap
   value against `B` values *including itself*. Every bootstrap quantile was
   therefore scaled by `c = B/(B+1)`, and the max over `L` levels compounded it:
   `size = 1 - (1-alpha) c^L`. Measured size at `alpha = 0.05`, `B = 100`:

   | L | 1 | 2 | 5 | 10 | 20 | 40 |
   |---|---:|---:|---:|---:|---:|---:|
   | old | 0.049 | 0.064 | 0.087 | 0.129 | 0.203 | 0.336 |
   | **new** | **0.052** | **0.047** | **0.050** | **0.049** | **0.049** | **0.049** |

2. **The framing hid what the second stage is doing.** After the first
   quantilisation the `F_{n,l}` are already (one minus) p-values. Maximising them
   is just a multiple test across `L` dependent hypotheses, which is the setting
   of **minP** (Westfall & Young 1993; Romano & Wolf 2005). Saying so lets us cite
   the literature instead of re-deriving it, and makes the comparison to simple
   FWER control obvious.

Notation, accordingly: the per-level values are **test statistics** `S_l` (they
are statistics for the same null, computed on successively coarser tables), their
bootstrap p-values are `p_l`, and the minP statistic is `P_min`. The old `M`
("criteria") is retired; the code calls them `statistics` throughout.

## The algorithm

```
Input:  Observed statistics S-hat in R^L; bootstrap statistics S^(1),...,S^(B) in R^L.
Output: P-value p in [0, 1].

Set S^(0) = S-hat.                                   // B+1 exchangeable paths
for j = 0,...,B do
    for l = 1,...,L do
        p^(j)_l = (1 + #{k != j : S^(k)_l >= S^(j)_l}) / (B+1)     // leave-one-out
    end
    P^(j)_min = min_{l=1,...,L} p^(j)_l
end
Set p = (1 + #{b = 1,...,B : P^(b)_min <= P^(0)_min}) / (B+1).
```

Caption:

> The first stage turns each of the `L` candidate test statistics into a marginal
> bootstrap p-value, each individually valid. Because the statistics are computed
> on nested coarsenings of the same table they are heavily dependent, so their
> minimum is not itself a p-value; we calibrate it by the same construction
> applied to the bootstrap paths. Every one of the `B+1` paths is ranked against
> the other `B` by an identical rule, so the `P_min` are exchangeable under the
> null and `p` is exactly uniform on the `1/(B+1)` grid at finite `B`.

## Why it is exact (finite `B`, no asymptotics)

Under the null the `B+1` vectors `S^(j)`, `j = 0,...,B`, are i.i.d., hence
exchangeable. The map `S^(j) -> P^(j)_min` is of the form
`Psi(S^(j); {S^(k)}_{k != j})` with `Psi` symmetric in its second argument, so
`(P^(0)_min, ..., P^(B)_min)` is exchangeable. The rank of `P^(0)_min` among
`B+1` exchangeable values is therefore uniform on `{1,...,B+1}` (ties broken at
random), and `P(p <= alpha) <= alpha` for every `alpha`, with equality on the
grid. This is the standard Monte-Carlo test argument (Dwass 1957; Barnard 1963;
Hope 1968) applied to the second stage.

Note this needs nothing about `L`, about the dependence between levels, or about
the search being agglomerative — only that the same deterministic map takes each
`(T, Sigma)` draw to a path. The same argument therefore covers the divisive
search (see `DIVISIVE.md`) and any truncation of either.

At `L = 1` the two stages cancel and the construction reduces exactly to the
ordinary bootstrap p-value `(1 + #{b : S^(b) >= S-hat})/(B+1)`, which is why the
depth-0 comparators were never affected by the old defect.

## Two consequences of the `1/(B+1)` grid

**Ties must be broken at random.** Each level places exactly one path at the floor
`1/(B+1)`, so up to `L` paths tie there and the second stage is choosing among
them. Counting ties conservatively instead gives, at `B = 100`:

| L | 2 | 5 | 13 | 20 |
|---|---:|---:|---:|---:|
| randomised | 0.052 | 0.048 | 0.053 | 0.051 |
| conservative | 0.046 | 0.045 | **0.000** | **0.000** |

i.e. the conservative version cannot reject at all once `L` is comparable to
`alpha (B+1)`. This is worth a footnote in the paper.

**`B` bounds power, not just resolution.** The test is exact at every `B`, but the
tie at the floor is resolved by a coin flip, so small `B` throws away information.
At `L = 13` (`dx = dy = 8`), fixed alternative, `alpha = 0.05`:

| B | 100 | 200 | 400 | 1000 |
|---|---:|---:|---:|---:|
| power | 0.522 | 0.878 | 0.930 | 0.952 |

The experiment configs now use `n_boot = 1000`, where the curve has flattened. No
tie rule recovers this; only more draws do.

## Comparison: simple FWER control

The obvious alternative to calibrating `P_min` is to Bonferroni-correct it:
`p_bonf = min(1, L * min_l p_l)`. It is valid and needs no second stage, but its
resolution floor is `L/(B+1)` rather than `1/(B+1)` — so at `dx = dy = 8` and
`B = 100` it cannot reject at `alpha = 0.05` under *any* data, and it needs
`B >= L/alpha ~ 260` before the comparison is even non-degenerate. It is also
conservative for the usual reason: the `L` levels are nested coarsenings of one
table and hence strongly dependent, which is exactly the dependence minP's
resampling step accounts for and Bonferroni discards.

**minP dominates it pathwise, not just on average.** Write `P_min = m/(B+1)` for the
observed minimum. For each level `l`, exactly `m` of the `B+1` paths have
`p_l <= m/(B+1)` (stage 1 makes each row a permutation of the grid), so across `L`
levels at most `Lm` paths have `P_min <= m/(B+1)` — one of which is the observed
path. Hence

    p_minP = (1 + #{b : P_min[b] <= P_min[0]}) / (B+1) <= Lm/(B+1) = p_bonf

for every realisation, with no null assumption and no appeal to dependence. So minP
is never less significant than Bonferroni on the same data, and the gap is the value
of the resampling step. (`tests/test_calibrate.py` asserts this.)

Implemented as `catci.calibrate.bonferroni_pvalue` and as the `*_bonf` entries in
the experiment registry, so the paper's comparison has numbers behind it.

## References to cite

- Westfall, P.H. & Young, S.S. (1993). *Resampling-Based Multiple Testing*. Wiley.
  — minP.
- Romano, J.P. & Wolf, M. (2005). Exact and approximate stepdown methods for
  multiple hypothesis testing. *JASA* 100(469), 94–108. — and the Econometrica
  companion, "Stepwise multiple testing as formalized data snooping".
- Ge, Y., Dudoit, S. & Speed, T.P. (2003). Resampling-based multiple testing for
  microarray data analysis. *TEST* 12, 1–77. — the pooled `B+1` ranking used here.
- Beran, R. (1988). Prepivoting test statistics: a bootstrap view of asymptotic
  refinements. *JASA* 83, 687–697. — stage 1 is prepivoting; the bridge from
  "double bootstrap" to "minP".
- Dwass, M. (1957); Barnard, G.A. (1963); Hope, A.C.A. (1968). — exactness of the
  `(1 + #)/(B+1)` Monte-Carlo p-value at finite `B`.

All five should be checked against the actual volumes before they go in the
bibliography — they are from memory.
