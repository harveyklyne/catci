# A semi-synthetic UCI adult study, after the 401k design

Status: **design settled and verified; the separation result established; the
full study not yet run.** Worktree `.claude/worktrees/adult-semisynth`, branch
`worktree-adult-semisynth`, based on `adult-application`. All new code is
`experiments/adult_semisynth.py`; nothing else in the worktree is modified.

## The design being ported

In Klyne & Shah (2025, *AoS* 53(1)) §5.1.2, the 401k data supplies the real
predictors `(X, Z)` and only the response `Y` is simulated, from a known `f_P`.
That buys a known target parameter `θ_P` without giving up a realistic predictor
distribution — the point being to show coverage where the truth is checkable,
*then* apply the method to the real response in §5.3.

The catci analogue has to differ in one respect: the target is not a scalar
parameter but the truth of the null `X ⟂ Y | Z`. So what must be known exactly is
**whether conditional independence holds, and under the alternative by how much**.
"Coverage" becomes size at λ = 0 and power at λ > 0.

### The desirable outcome

The 401k study is not only a calibration check. Its regression functions are
*chosen* (§5.1.3: `f_plm`, `f_add`, `f_int`), and `f_int` exists specifically so
that a competitor breaks on it — the PLR method "risks completely losing coverage
when the response is nonlinear in `X`". So the target here is the same shape: **a
semi-synthetic setting, real predictors, exactly known truth, in which the
competing method has essentially no power and catci does.**

That target is met, and the construction is analytic rather than empirical —
Ankan & Textor's Q1 population signal is driven to exactly zero. At `n = 2000`
against a true alternative, Q1 rejects 6.0% of the time (its size is 5.2%) while
catci rejects 100%. See [Separating catci from Ankan & Textor](#separating-catci-from-ankan--textor).

## The construction

Population = the full 30,162-row empirical distribution of the adult data, the
same posture the 401k study took toward its 9,915 rows. Every "true" quantity is
then a finite sum, known exactly rather than estimated.

Fix a **generating stratum** `W` — discrete conditioning variables, default
`{Age, Sex}`. One replicate:

1. Draw `n` rows. `X` and every `Z` variable are the **real, untouched** values.
2. Draw `Y_i` independently from the mixture kernel

   ```
   h_λ(y | x, w)  =  (1 − λ) · P(y | w)  +  λ · P(y | x, w)
   ```

   both tables being the population's own.

λ = 0 draws `Y` from `Z` alone using fresh randomness, so `Y ⟂ (X, everything) | W`.
λ = 1 reproduces the population's own joint law of `(X, Y) | W` — the **real**
dependence, not an invented one. λ interpolates.

### Answering the permutation question directly

Yes — this *is* the permutation idea, in its cleaner with-replacement form. The
three variants, and why I ended up here:

| variant | null exact? | i.i.d.? | true `f`, `g` |
|---|---|---|---|
| permute `Y` within `W`-cells | yes | **no** — sampling without replacement makes `Y_i, Y_j` in a cell negatively dependent | real, nonparametric |
| **resample `Y` within `W`-cells (used)** | yes | yes | real, nonparametric |
| draw `Y ~ ĝ(Z)` from a fitted model | yes | yes | a model, not the data |

Straight permutation is the natural first thought and it does enforce
conditional independence — but it breaks the i.i.d. assumption the test's
calibration is derived under, and does so worst exactly where cells are small.
Resampling within cells fixes that at no cost, since the within-cell empirical
distribution is the same object either way. The third variant is what you would
need if `Z` were continuous; on adult it is unnecessary and strictly worse,
because it replaces a real conditional with a fitted one.

### Two properties that made me choose the mixture over a tilt

Taking `Z = W`:

* `f(x | w) = P(x | w)` is untouched — `X` is never modified.
* `g(y | w) = Σ_x P(x | w) h_λ(y | x, w) = P(y | w)` for **every** λ, since the
  second component averages back.

So λ is a *pure conditional-dependence knob*: the two nuisance functions a
GCM-type test must estimate do not move as it turns. An exponential tilt moves
`g` while it moves the dependence, which would confound "the test lost
calibration" with "the test gained power".

Second, the departure from independence is exactly **linear** in λ,

```
P_λ(x,y|w) − P(x|w)P(y|w)  =  λ · [P(x,y|w) − P(x|w)P(y|w)]
```

so any chi-square non-centrality is `n · λ² · ncp_per_n`, and a power grid can be
placed without a pilot run.

### Making the nuisances hard without breaking the null

`W` need not be all of `Z`. Since λ = 0 gives `Y ⟂ (X, Z) | W`, weak union gives
`Y ⟂ X | Z` for **any** `Z ⊇ W`. So generating under `W = {Age, Sex}` while
*testing* with a rich `Z` keeps the null exact and makes the estimation genuinely
hard: `f(x | Z)` is the real, unknown, high-dimensional conditional, and the
learner must discover for itself that `g` depends on only two `Z` variables.

Cost: for λ > 0 the identity `g = P(y | w)` no longer holds exactly, because the
mixture's second component then averages over `P(x | Z)` rather than `P(x | w)`.
Both remain closed-form on a finite population; `true_propensities` returns
whichever applies. This suggests running **size with rich `Z`, power with `Z = W`**.

## Verification

All four checks pass (`scratchpad/verify.py`, Occupation × Income, `W = {Age, Sex}`):

| claim | result |
|---|---|
| `g_true` invariant in λ at `Z = W` | max abs difference `2.2e−16` |
| drawn `Y` matches `P(y|w)` at every λ | max abs error 0.0069 / 0.0093 / 0.0060 on 200k draws at λ = 0 / 0.5 / 1 |
| departure exactly linear in λ | ratio to `λ·δ₁` in `[0.100000, 0.100000]` at λ = 0.1, likewise 0.3, 0.7 |
| λ = 0 really is conditionally independent | stratified χ² = 169.4 on df = 182, **p = 0.74**, on a 300k draw |

The λ² law is accurate to five significant figures: at λ = 1, n = 300k, predicted
non-centrality 30399.9 against observed χ² of 30399.4.

## The data

`n = 30,162` after listwise deletion (the `adult.py` default, itself pinned by the
Fig. 1b df arithmetic — see `ADULT.md`).

`Z = {Age, Sex}` is the only conditioning set with no sparsity problem at all —
14 strata, minimum occupancy **127**. Richer sets degrade fast:

| Z | possible | occupied | min | median | rows in cells < 20 |
|---|---:|---:|---:|---:|---:|
| Age+Sex | 14 | 14 | 127 | 1518 | 0.0% |
| Age+Sex+Race | 70 | 68 | 1 | 55 | 0.5% |
| Age+Sex+MaritalStatus | 98 | 89 | 1 | 51 | 0.6% |
| Age+Sex+Education | 224 | 215 | 1 | 33 | 2.2% |
| Age+Sex+Race+MaritalStatus | 490 | 307 | 1 | 7 | 3.7% |
| Age+Sex+Education+MaritalStatus | 1568 | 888 | 1 | 5 | 10.3% |

Real conditional dependence given `{Age, Sex}`, as non-centrality per unit `n`
at λ = 1 (so `ncp = n·λ²·`this) — the top and the interesting tail:

| X | Y | dx | dy | CMI | ncp_per_n |
|---|---|---:|---:|---:|---:|
| MaritalStatus | Relationship | 7 | 6 | 0.508 | 1.030 |
| Education | Occupation | 16 | 14 | 0.258 | 0.580 |
| Workclass | Occupation | 7 | 14 | 0.124 | 0.316 |
| Education | Income | 16 | 2 | 0.056 | 0.113 |
| **Occupation** | **Income** | **14** | **2** | **0.053** | **0.101** |
| HoursPerWeek | Income | 4 | 2 | 0.018 | 0.035 |
| **Income** | **Workclass** | **2** | **7** | **0.008** | **0.017** |

Suggested headline pairs:

* **Occupation × Income** — the interpretable one. At `n = 1000`, λ ∈ [0, 0.5]
  traces a full power curve. The scientific question ("which occupations group
  together as high-earning, given age and sex?") is exactly what the merging
  search is supposed to answer, so the recovered partition is a second, richer
  readout than the p-value.
* **Income × Workclass** — the honest exception from the `ADULT.md` sweep, the
  one pair that did not reject at `n = 1000` under either Ankan–Textor estimator.
  Weakest real signal in the table, so it is where adaptive merging should earn
  its keep against a stratified test.
* **Education × Income** — the separation pair. Both are typed ordinal, which is
  what sends Ankan & Textor down the one-degree-of-freedom Q1 branch and makes
  the blinding of the next section possible. `dx = 16` gives ample room to plant
  a direction orthogonal to the Li–Shepherd score.

### Rare levels have to be pooled, and pooling is free

At `n = 1000`, Occupation's `Armed-Forces` has an expected count of **0.3** and
`Priv-house-serv` **4.7**; Workclass's `Without-pay` **0.5**. A level that is
usually absent contributes a near-zero column to `Σ`, and `ApproxChi` divides by
it — numerically degenerate, not merely uninformative.

`pool_rare` merges them **once, before anything is simulated**, so it is a
coarsening of the population and every claim above survives verbatim. The
alternative — dropping unobserved levels per replicate — would make the level set
random and the "true" propensities replicate-dependent. Pooling costs essentially
nothing: `ncp_per_n` for Occupation × Income moves 0.1013 → 0.1012.

## Separating catci from Ankan & Textor

### Why a blind spot has to exist

`ankan_textor.residual_matrix` returns **one** column for a variable typed
ordinal — the Li–Shepherd residual — and `k−1` columns for one typed categorical.
So for a pair of ordinal variables, A&T's Q1 is one number on one degree of
freedom, and its population signal is the single scalar

```
S = E[ r_X(X, W) · r_Y(Y, W) ],    r_X(x, w) = P(X < x | w) − P(X > x | w).
```

A test whose entire signal is one linear functional is blind to every direction
orthogonal to it. Since `Population.kernel` now takes an arbitrary direction
`delta`, we can plant a departure from independence with `S = 0` — at every λ and
every `n` — while leaving intact the rest of the dependence, which catci's
full-table criterion does see.

### The construction, and a false start worth recording

Three constraints must hold at once. C1 (`h` is a pmf) and C2 (`g` stays fixed)
are what make λ a pure dependence knob; C3 is the new one:

```
C0  delta = 0 wherever P(y | w) = 0          -- support
C1  Σ_y  delta(x, y, w) = 0                  -- h is a pmf
C2  Σ_x  P(x|w) delta(x, y, w) = 0           -- g fixed at P(y|w)
C3  ⟨G, delta⟩ = 0,  G = r_X ⊗ r_Y           -- Q1 blind
```

The obvious route — take the *real* dependence and project out its `G`
component (`project_out_visible`) — is exactly orthogonal and **practically
useless**. The leftover direction puts its mass on cells where `P(y|w)` is nearly
zero, so `max_lambda` is 0.004 for Education × Income and the attainable effect
size is negligible. That is why `planted_direction` *constructs* the direction
instead, carrying a `P(y|w)` factor:

```
delta(x, y, w) = c · P(y | w) · u_w(x) · v_w(y)
```

C0 is then automatic, non-negativity is a bound on `c·u·v` alone, and scaling
makes `max_lambda` exactly 1. C2 and C3 become two conditions on `u_w`, imposed
by Gram–Schmidting a template against `{1, r_X}` in the `P(x|w)` inner product,
**separately in each stratum**, so they hold exactly rather than on average.

Two design points:

* **Blinding goes through `X`, not `Y`.** With `dy = 2` the space of
  `P(y|w)`-centred functions of `y` is one-dimensional and `r_Y` spans it, so no
  admissible `v` is orthogonal to it. Blinding through `X` needs only `dx ≥ 3`,
  works for any `dy`, and is the interpretable half — it is a statement about
  which levels of `X` group together.
* **`v_w = r_Y`.** The entire `Y`-side signal is put in the direction A&T would
  have been *most* sensitive to, so the blindness cannot be dismissed as an
  unlucky choice on that side.

The default template is the second harmonic `cos(2π(x−1)/(dx−1))` — the canonical
non-monotone contrast, separating the middle levels from both extremes. On
Education × Income given `{Age, Sex}`, in the largest stratum, it produces:

| Education | P(>50K \| x, w) at λ=1 |
|---|---:|
| Preschool | 0.572 |
| 9th | 0.452 |
| **HS-grad** | **0.320** |
| **Some-college** | **0.292** |
| Bachelors | 0.360 |
| Doctorate | 0.401 |

against a stratum baseline of 0.338 — a U-shape in education. **Be honest about
what is realistic here**: `f`, `g`, and the predictor distribution are exactly
the real ones, but the *sign pattern* of the interaction is designed, not real.
That is precisely the status of `f_int` in the 401k study.

### It works, exactly

| quantity | real direction | planted direction |
|---|---:|---:|
| `q1_signal` | 0.0694 | **4.9e−17** |
| `visible_fraction` (squared cosine with `G*`) | 0.717 | **5.3e−33** |
| C1 / C2 residuals | — | 2.8e−16 / 1.3e−16 |
| `max_lambda` | 1 | 2.05 |

Empirically, on a **300,000**-row draw at λ = 1: mean `r_X·r_Y` = +0.00050
(SE 0.00042), giving **Q1 = 1.41 on df = 1**. The same draw from the real
direction gives Q1 = 28,904. Q1 is at its nominal level against an alternative
this large, at any sample size.

### The four-way comparison

Education × Income given `{Age, Sex}`, **every method given the exact true
propensities** so this compares statistics and not nuisance fits. 40 reps at
`n = 2000` (SE 0.079) for the catci columns; 400 reps (SE 0.011) for the A&T
columns from a separate, cheaper run.

| direction | λ | A&T Q1 | stratified χ² (df≈210) | catci ordinal | catci greedy |
|---|---:|---:|---:|---:|---:|
| planted | 0 | 0.052 | 0.100 | 0.025 | 0.125 |
| **planted** | **2.0** | **0.060** | 0.975 | **1.000** | **1.000** |
| real | 0.3 | 1.000 | 0.250 | 0.900 | 0.850 |

**The two directions reverse the ranking completely.** On the real Education →
Income effect, which is strongly monotone, A&T's df = 1 is perfectly aimed and it
beats catci. On the planted direction it has no power at all. That contrast is a
better result than a uniform win would have been: it identifies the mechanism
rather than just the outcome.

### The caveat that makes the claim sharper, not weaker

The blindness is entirely a consequence of typing Education **ordinal**. Typed
categorical, A&T uses Q2 with df = 15, which spans all directions in `X` and is
not blind. 400 reps at `n = 2000`, oracle propensities, SE 0.011:

| direction | λ | Q1 (ordinal, df=1) | Q2 (categorical, df=15) |
|---|---:|---:|---:|
| planted | 0 | 0.052 | 0.058 |
| **planted** | **2.0** | **0.060** | **1.000** |
| real | 0 | 0.052 | 0.058 |
| **real** | **0.1** | **0.275** | **0.140** |

So A&T has a knob — the ordinal/categorical typing — that **the paper never
specifies** (already flagged in `ADULT.md` as a choice left open), and *each
setting has a blind spot*: ordinal buys power against monotone alternatives
(0.275 vs 0.140) at the price of zero power against non-monotone ones; categorical
buys omnidirectional power at the price of df = 15. catci requires no such choice
and is competitive in both regimes.

That is the claim worth making in the paper — not "our method is better", but
"their method requires an unstated choice, either setting of which is exploitable,
and ours does not."

### Honest negative: the partition is not a clean two-group readout

I expected the greedy search to collapse Education into `{middle}` vs
`{extremes}`, recovering the planted contrast. It does not. Across 40 replicates
at λ = 2 the criterion peaks at 7–9 groups, and never at 2:

```
n_groups   5  6  7  8  9 10 11 12 13 14 15
count      1  2  6  9  6  3  4  1  3  2  3
```

Power is unaffected — the p-value is fine — but the "interpretable coarsening"
selling point does not come for free from `ApproxChi`'s peak. Worth understanding
before leaning on partition recovery as a headline.

## What the pilot runs found

### 1. The construction is clean; the mild over-rejection is not ours

`lam = 0`, **oracle** propensities (so nothing here is about nuisance
estimation), `normalise=False`, 500 reps, binomial SE 0.010, α = 0.05:

| population | n | tree | mGCM |
|---|---:|---:|---:|
| adult Income × Race (dx 2, dy 5) | 1000 | 0.052 | 0.072 |
| adult Occupation × Income (pooled) | 1000 | 0.068 | 0.060 |
| adult Occupation × Income (pooled) | 5000 | 0.058 | 0.128 |
| **the repo's own synthetic DGP, strength 0** | 1000 | **0.098** | **0.166** |

I started this pilot because a 20-rep look suggested the adult null was
over-rejecting, and an exactly-known null that over-rejects would have meant a
broken construction. It is the other way round: **the repo's own synthetic null
control over-rejects roughly twice as much as the adult semi-synthetic one.**
Whatever the inflation is, it is a property of the statistic/calibration path,
present on the package's own DGP, and not introduced by this design.

That is worth chasing separately — it bears on the stale-results item in the
README's "Known gaps" — but it is not a blocker here, and the semi-synthetic
study will in fact *measure* it against a null that is exact by construction
rather than by a DGP argument.

An earlier, unpooled 400-rep sweep (Occupation with all 14 levels including
`Armed-Forces`) gave 0.095–0.128 for the saturated and tree searches, versus
0.052–0.068 pooled — consistent with the degenerate-column story above, and the
reason `pool_rare` exists.

### 2. `normalise` disagrees between the API and the experiments

`catci_test` defaults to `normalise=True`; every config in `experiments/config.py`
uses `normalise=False`. In the 400-rep unpooled sweep the True setting was
consistently worse calibrated (e.g. Occupation × Income at n = 1000: greedy 0.095
→ 0.128). Not necessarily a bug, but the two defaults should not disagree
silently.

### 3. Saturated search is ~10× the cost of tree search

At `dx = 13`, `dy = 2`, measured on one replicate:

| search | per search | per replicate at `n_boot = 100` |
|---|---:|---:|
| `Saturated()` (greedy) | 73.8 ms | **7.5 s** |
| `Tree.binary(13)` | 7.8 ms | 0.8 s |

This is what made the first pilot take 43 minutes for 3200 fits. It is direct
evidence for TODO item 2 (the reverse/divisive greedy procedure): the saturated
search's `C(d,2)` candidate count per level is the whole cost, and a real-data
`X` with 13–16 levels is precisely where it bites. Budget accordingly — a 6-point
λ grid × 500 reps × saturated search is ~6 CPU-hours per pair.

## What is not done

* **The full λ grid.** A 200-rep, 9-point grid at `n = 1000` (λ ∈ {0, 0.5, 1,
  1.5, 2} planted and {0, 0.05, 0.1, 0.2} real) was launched and **killed before
  it finished**, so the four-way table above rests on 40 reps (SE 0.079) for the
  catci columns. The separation is far larger than that error bar, but the
  numbers should be re-run before they go in a paper. `scratchpad/gap.py` takes
  `reps` and `n` as arguments and is ready to relaunch; budget ~25 min on 8 cores.
* **A smaller-effect grid point.** At λ = 2, `n = 2000` the stratified χ² also
  rejects (0.975), so that row separates catci from A&T but *not* from the
  stratified baseline. A lower λ or `n` should separate all three — worth finding
  the point where it does, since "beats both competitors at once" is the stronger
  figure.
* **The learned-nuisance arm.** Everything above is oracle-propensity. A 20-rep
  look with `xgboost_learner` on the unpooled pair gave 0.25 at α = 0.05, but at
  20 reps that is ±0.10 and it was the unpooled setting; it needs redoing at
  proper rep count on the pooled population before it means anything.
* **Hyperparameters.** `experiments/tuning/` only has `n1000_numclass8`, and the
  README notes the tuner was never ported. Occupation is `numclass=13`, Education
  16. Either reuse the `d = 8` parameters and say so, or port a tuner.
* **The rich-`Z` size arm**, which is the part that actually stresses nuisance
  estimation, and the whole reason the weak-union trick is in the module.
* **Other pairs for the blinding result.** Only Education × Income was run.
  `HoursPerWeek × Education` is the other ordinal pair with room (`dx = 4`,
  `visible_fraction` of the *real* dependence only 0.068 — A&T is already nearly
  blind there without any planting, which is a result in itself).
* **Nothing is pytest-pinned.** `pool_rare`, `planted_direction` and the C1/C2/C3
  identities are all verified by hand in scratch scripts. The three identities in
  particular are exact and cheap, so they should be unit tests.

## Files

| file | what |
|---|---|
| `experiments/adult_semisynth.py` | the whole construction — see below |
| `ADULT.md` | the earlier Ankan–Textor reproduction on the same data |
| `experiments/adult.py` | loading and preprocessing (unchanged) |
| `experiments/ankan_textor.py` | the competitor, unchanged |

`adult_semisynth.py` in two halves:

| | |
|---|---|
| the population | `Population`, `build_population`, `pool_rare`, `delta_real` |
| sampling | `draw`, `true_propensities`, `Population.kernel`, `ncp_per_n` |
| blinding | `planted_direction`, `visible_direction`, `project_out_visible`, `q1_signal`, `visible_fraction`, `is_valid_direction`, `max_lambda` |

Scratch scripts (session scratchpad, not in the worktree): `verify.py` (the four
invariant checks), `size2.py` (400-rep unpooled sweep), `size3.py` (the
three-population control), `blind.py` (visible fractions and the failed
projection route), `blind2.py` (planted-direction validation and the 300k draw),
`gap.py` (the four-way comparison — killed mid-grid), `kindcheck.py` (Q1 vs Q2
typing). These are worth promoting into `experiments/` if the result is kept.
