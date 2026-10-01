# Ankan & Textor on UCI adult income — what we reproduced

Reproduction of the conditional independence tests of *A Simple Unified Approach to
Testing High-Dimensional Conditional Independences for Categorical and Ordinal Data*
(Ankan & Textor, AAAI-23, pp. 12180–12188) on the adult-income dataset.

Scope: **their method's p-values only**. The learned skeleton (Fig. 8a), the F1
comparison (Fig. 8b) and the competing baselines are out of scope.

## The target had to be reframed

The paper prints **no p-values from its own method** — on this dataset or any
other. The six p-values in Fig. 1b come from the *baseline it argues against*:

> **(a)** Skeleton estimated by the stable PC algorithm … from 1000 samples of the
> adult income data **using the default CI test in the R package 'bnlearn', a
> stratified mutual information test**. … **(b)** A closer inspection of test
> results reveals high degrees of freedom that sometimes exceed the sample size.

Their method's only appearances on this data are a skeleton drawing and an F1
curve. So there is no published number for our p-values to be checked against,
and this is a **port validated by properties and by one numeric anchor**, not a
value-for-value replication. Nothing here should be cited as "matching the paper".

There is also no replication package, and the descendant implementations the
authors ship are a different statistic — pgmpy's `pillai` and dagitty's
`cis.pillai` both swap the chi-square quadratic form for a canonical-correlation
trace under an F-approximation, and dagitty's drops the Li–Shepherd residual
entirely. See `ankan_textor_code_sources_2026-08.md`.

## The one thing that *is* checkable: Fig. 1b's df

Degrees of freedom are a deterministic function of the coding,
`df = (k−1)(r−1) × n_strata`, so the six df in Fig. 1b pin the preprocessing.
All six reproduce, which fixes `n_strata = 14` (Age in 7 bins × Sex in 2) and the
level counts:

| pair | k, r | (k−1)(r−1) | ×14 | paper |
|---|---|---|---|---|
| Edct–Wrkc | 16, 6 | 75 | 1050 | 1050 ✓ |
| Occp–Wrkc | 13, 6 | 60 | 840 | 840 ✓ |
| Rltn–HrPW | 6, 4 | 15 | 210 | 210 ✓ |
| Incm–Occp | 2, 13 | 12 | 168 | 168 ✓ |
| Incm–Wrkc | 2, 6 | 5 | 70 | 70 ✓ |
| Incm–HrPW | 2, 4 | 3 | 42 | 42 ✓ |

This settled a choice the paper never states. Workclass must show **6** levels and
Occupation **13**, below their full-data counts. Keeping `?` as a category can
only *add* a level, so it cannot produce those numbers; **listwise deletion can,
and does** — 30,162 rows, and a 1000-row subsample missing the rarest levels
(`Without-pay`, `Never-worked`, `Armed-Forces`) lands exactly on 6 and 13. Two of
ten seeds tried reproduce all six df; seed 6 is the one used below.

Hence `adult.load()` defaults to `drop_missing=True`. Pinned by
`tests/test_adult.py`, including the negative case that retaining `?` cannot match.

## Headline result

The paper's introduction argues that stratification-based tests fail on this data
because df explodes, not because the variables are independent. On the **same 1000
rows**, that argument is now quantified — Fig. 1b's baseline beside our port:

| X | Y | stat | df (MI) | p (MI) | df (A&T) | p (GLM) | p (RFT) |
|---|---|---|---:|---:|---:|---:|---:|
| Education | Workclass | Q2 | 1050 | 1.00 | 5 | 1.4e−06 | 4.9e−06 |
| Occupation | Workclass | Q3 | 840 | 1.00 | 60 | 2.1e−28 | 3.9e−26 |
| Relationship | HoursPerWeek | Q2 | 210 | 0.99 | 5 | 3.1e−05 | 5.9e−05 |
| Income | Occupation | Q2 | 168 | 0.08 | 12 | 6.9e−16 | 2.3e−15 |
| Income | Workclass | Q2 | 70 | 0.50 | 5 | **0.168** | **0.170** |
| Income | HoursPerWeek | Q1 | 42 | 0.003 | 1 | 1.7e−10 | 2.0e−10 |

Five of six flip from "no evidence" to overwhelming, with df collapsing by one to
two orders of magnitude. Income–Workclass does not reject under either test — the
honest exception, and worth keeping in view.

## Full sweep

`Z = {Age, Sex}`, X and Y ranging over the other nine variables: **36 pairs**,
each under both of the paper's estimators.

**At n = 30,162 (all rows): every one of the 36 pairs rejects at α = 0.05, under
both GLM and RFT.** The largest p-value across the whole sweep is 1.1e−12
(HoursPerWeek–NativeCountry); eight underflow to exactly 0 in float64. GLM and RFT
agree on every pair, and disagree on nothing.

This is what the paper asserts should happen — "any reasonable structure should be
dense" — but note it is *unfalsifiable in this form*: a test that rejected
everything unconditionally would look identical. The n = 1000 run is the more
informative one, where 20 of 26 well-conditioned pairs reject and six do not.

Reproduce with:

```sh
python run_adult.py                       # n = 30,162, GLM + RFT, ~1.4 min
python run_adult.py --n 1000 --seed 6     # the Fig. 1 setting, ~8 s
```

Results land in `results/*.parquet` with a provenance sidecar. `results/` is
gitignored, so the numbers above are the tracked record.

## A failure mode the paper does not discuss

Dropping one dummy makes `Sigma_d` full rank in principle, but not in practice.
When a level of X and a level of Y never co-occur within a stratum of Z, the
product column `R_I(x=l) · R_I(y=m)` is identically zero and `Sigma_d` loses rank.
With `Z = {Age, Sex}` giving only 14 strata, this bites on every pair involving
NativeCountry (41 levels).

It matters because **the failure is silent**. `np.linalg.solve` raises only on
exact singularity, so a condition number of 1e21 returns a plausible-looking
number rather than an error. At n = 2000 we measured GLM returning p-values with
`cond ≈ 1.3e21`, no warning of any kind. RFT at least errors, because random-forest
probabilities hit exact zeros.

`q_statistic` therefore always reports `rank` and `cond`, and results carry a
`well_conditioned` flag. At full n only 4 of 72 tests are ill-conditioned (the
NativeCountry pairs, each rank-deficient by exactly one); at n = 1000 it is 20 of
72. Treat those as undefined, not as evidence. `--method pinv` substitutes the
pseudo-inverse with `df = rank(Sigma_d)`, which keeps them computable but departs
from the paper's stated null.

## Port fidelity

Implemented from the paper text (pp. 12182–12183); see `ankan_textor.py`.

* **LS residual** `p̂(Y<y_i|z) − p̂(Y>y_i|z)`, collapsing to `y_i − p̂(Y=1|z)` for
  binary variables, and applied per dummy for categorical ones.
* **Q1/Q2/Q3 are one statistic at three shapes.** With residual matrices of width
  `a` and `b`, `Q = (1/n) d Σ_d⁻¹ dᵀ` on the `a·b` elementwise-product columns has
  `df = a·b`, which is 1, `k−1` and `(k−1)(r−1)` exactly as Propositions 1–3 state.
  One implementation, three propositions.
* **Estimators**: GLM = multinomial logistic (`nnet::multinom`) / proportional-odds
  (`VGAM`); RFT = 50-tree probability forest (`ranger`). Python equivalents are
  scikit-learn `LogisticRegression`, statsmodels `OrderedModel`, and
  `RandomForestClassifier(n_estimators=50)`.

No reference implementation exists to differential-test against, so
`tests/test_ankan_textor.py` pins the properties the paper itself asserts: the
binary residual collapse, Q1 equalling the squared GCM of Shah & Peters (2020),
the df of Propositions 1–3, symmetry in X and Y, the singularity that dropping a
dummy avoids — and calibration under a null, which is the test that would actually
catch a wrong statistic.

## Choices the paper leaves open

Marked `NOTE:` in the source.

* **Z encoding** — unstated. We one-hot encode and drop the first level per
  variable. R would use polynomial contrasts for ordered factors.
* **Ordinal vs categorical typing** — unstated. We take Age, HoursPerWeek,
  Education and Income as ordinal, the rest categorical.
* **HoursPerWeek bin 3** — the paper's own list overlaps at 30 ("21-30, 30-40").
  Read as 31-40.
* **Penalisation** — `nnet::multinom` is unpenalised; scikit-learn always
  penalises, so `C = 1e6` approximates the MLE. Some 41-class fits still hit the
  iteration cap.
* **Seeds** — none given, so RFT cannot match the paper's RNG.

## Files

| file | what |
|---|---|
| `adult.py` | data loading, the paper's preprocessing, `fig1b_df_check` |
| `ankan_textor.py` | LS residuals, GLM/RFT estimators, Q1/Q2/Q3 |
| `run_adult.py` | the 36-pair sweep → parquet + provenance |
| `data/adult-train.csv` | vendored UCI train split (32,561 rows) |
| `ankan_textor_code_sources_2026-08.md` | why there is no code to port from |
| `../tests/test_adult.py`, `../tests/test_ankan_textor.py` | the pins |
