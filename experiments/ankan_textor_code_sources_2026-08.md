# Code sources for Ankan & Textor (AAAI-23) — reproducing the Adult-income results

Research note, 23 Aug 2026. Target: *A Simple Unified Approach to Testing High-Dimensional
Conditional Independences for Categorical and Ordinal Data*, Ankan & Textor, AAAI-23
(pp. 12180–12188; arXiv:2206.04356v3). Goal is reproducing **their own method** on Adult
(Fig. 1a/1b, Fig. 8a/8b) — baselines out of scope.

## Bottom line

**There is no replication package.** No code-availability statement in the AAAI camera-ready or in
arXiv v3; no supplementary link; nothing under either author's GitHub matching the naming pattern
Ankan uses for paper repos (`2020-dagitty-manual`-style — probed ~30 candidate names, control repo
confirmed the probe works). The experiments were run in R against `bnlearn`/`pcalg`/`nnet`/`VGAM`/
`ranger` and appear never to have been released.

**But both authors ship a descendant of the method in their maintained packages**, and neither is
the paper's statistic. That distinction is the main thing to get right before building:

| Source | Language | Residuals | Statistic | Null |
|---|---|---|---|---|
| **Paper** | R (unreleased) | Li–Shepherd LS residuals; multinomial logit (`nnet`) or proportional-odds (`VGAM`), or `ranger` probability forest (50 trees) | Hotelling T²-type quadratic form Q1/Q2/Q3 on residual dot-products | χ²(1), χ²(k−1), χ²((k−1)(r−1)) |
| **pgmpy** `pillai` | Python | dummy − predicted class probability, drop last col; sklearn RF | Pillai's trace of canonical correlations between residual *matrices* | F-approximation (Muller & Peterson 1984) |
| **dagitty** `cis.pillai` | R | **linear** `lm` residuals on dummy/integer codes | Pillai's trace via `cancor` | `CCP::p.asym` asymptotic |

So: pgmpy is the right scaffolding and the right *family*, but a faithful Fig. 8 needs the paper's
Q-statistics written from scratch. Good news is that pgmpy's PC accepts a custom callable CI test,
so the port only has to be the statistic — see below.

---

## 1. pgmpy — Ankan's own Python implementation (closest runnable code)

Repo: <https://github.com/pgmpy/pgmpy> (HEAD 1.1.2, Aug 2026). Relevant files:

- `pgmpy/ci_tests/pillai_trace.py` — registry name `"pillai"`, class `PillaiTrace`. Docstring cites
  `ankan_textor_2023` and `li_shepherd_2010`.
- Siblings sharing the same residual machinery: `hotelling_lawley.py` (`HotellingLawley` — the
  Hotelling–Lawley *trace*, still an F-approximation, not the paper's T²), `wilks_lambda.py`,
  `roys_largest_root.py`, `gcm.py`, `generalized_cov.py`.
- `pgmpy/ci_tests/_base.py` → `_ResidualMixin.get_residuals()` is where residualization lives:
  Z one-hot encoded with an appended `_intercept_Z` column of ones; `RandomForestClassifier(random_state=0)`
  for categorical/ordinal targets, `RandomForestRegressor` for numeric; residual matrix =
  `get_dummies(X) − predict_proba(Z)` with the **last column dropped**.
- `docs/references.bib` line 81: the paper is `ankan_textor_2023`, keyed `metrics_and_independence_tests`.

**Provenance.** Introduced as `ci_pillai` in PR #1779, commit `724e90eb`, 17 Jun 2024, by Ankur Ankan
("Residualization based CI test") — that first version used **XGBoost** (`XGBClassifier`/`XGBRegressor`),
not random forests. Renamed `ci_pillai` → `pillai_trace` in `bee9971b` (Sep 2024); the estimator later
became sklearn RF, which is *closer* to the paper's RFT than the XGBoost original.

**Where it departs from the paper** (all four matter for reproducing Fig. 8):

1. **Statistic and null.** Paper: `Q(x,y) = (1/n)·d·Σ̂_d⁻¹·dᵀ` where `d` is the vector of dot products
   between residual indicator columns, asymptotically χ² with (k−1)(r−1) df. pgmpy: sum of squared
   canonical correlations with an F-approximation, effect size = partial η² = V/s. Different df,
   different p-values, different effect size — a direct numeric comparison against the paper's
   Fig. 1b table will *not* match.
2. **Ordinality is dropped.** The paper's LS residual `r_i = p̂(Y<y_i|z) − p̂(Y>y_i|z)` is what makes
   the ordinal case work; pgmpy's dummy-minus-probability residual discards the order (and pgmpy's
   `Adult` dataset is tagged `is_ordinal: False`).
3. **Estimator defaults.** Paper's RFT = `ranger` probability forest with **50 trees**, otherwise
   defaults. pgmpy = sklearn RF with 100 trees. Both accept an injected estimator
   (`PillaiTrace(data, estimator=...)`), so this one is cheap to align.
4. **`gcm` is the one honest overlap.** For two *ordinal* variables the paper's Q1 is exactly the
   squared generalized covariance measure (Shah & Peters 2020) — so pgmpy's `gcm` is a legitimate
   cross-check for the ordinal–ordinal case only.

**The highest-leverage fact for the build:** `pgmpy/causal_discovery/PC.py` takes
`ci_test: str | Callable`, with `variant="stable"` (order-independent, Colombo & Maathuis — exactly
what the paper used via `pcalg::pc(stable.fast)`) or `"parallel"` (same result, multicore). A faithful
Q1/Q2/Q3 implementation can be dropped in as a callable and the skeleton search, v-structures,
Meek rules and CPDAG come for free. Note the API moved: old `pgmpy.estimators.CITests` functions →
`pgmpy.ci_tests` classes.

**pgmpy also ships Adult**: `pgmpy/datasets/adult.py` → `load_dataset("adult")`, 32,561 rows, mixed
types, ordinal `educ` (16 levels, order spelled out in the file) and `Income`, categorical
`workclass, mar-stat, occup, relat, race, sex, nat-count`, plus an expert-knowledge file. It pulls
from the HF hub repo `pgmpy/example_datasets`, so it needs network at build time.

## 2. dagitty — Textor's own R implementation (second descendant, further from the paper)

Repo: <https://github.com/jtextor/dagitty>, `r/R/ci-tests.R` → `.ci.test.pillai`, exposed as
`localTests(..., type = "cis.pillai")`. Added Oct 2023 by Textor. Categorical variables are dummy-coded
via `fastDummies` (dropping the *most frequent* category), ordinals are integer-coded, X and Y are
residualized by **`lm`**, then `cancor` + `CCP::p.asym(tstat="Pillai")`; effect size is the RMS canonical
correlation. Documented in `r/man/localTests.Rd`.

Useful as a second reading of the canonical-correlation branch and as a sanity oracle on toy data —
not usable for Fig. 8, since linear residualization of dummies is not the paper's GLM/RF LS residual.

## 3. Building blocks for a faithful port

The paper is self-contained enough to port; the pieces and their Python equivalents:

- **LS residual** (Li & Shepherd 2012): `p̂(Y<y_i|z) − p̂(Y>y_i|z)`; for binary Y it collapses to
  `y_i − p̂(Y=1|z)`. Reference implementation to unit-test against: R package **`PResiduals`**
  (Li & Shepherd's own — `presid()`, COBOT), and `VGAM`.
- **GLM estimator**: multinomial logit ≈ sklearn `LogisticRegression(multinomial)`; the ordinal
  proportional-odds model (`VGAM::vglm(cumulative(parallel=TRUE))`) ≈
  `statsmodels.miscmodels.ordinal_model.OrderedModel(distr="logit")` or `mord`.
- **RFT estimator**: `ranger` probability forest ≈ `RandomForestClassifier(n_estimators=50).predict_proba`.
- **Σ̂_d**: sample covariance of the per-observation product vectors (element-wise products of the
  residual indicator columns); statistic `(1/n)·d·Σ̂_d⁻¹·dᵀ`. Watch the dropped dummy — Σ̂_d is
  rank-deficient without it.
- **PC-stable**: pgmpy `PC(variant="stable")`, or `causal-learn` if a second opinion is wanted.

## 4. The Adult recipe as the paper states it

- 11 variables in Fig. 1a / 8a: Income (binarized at $50K), Workclass, Education, Marital Status,
  Occupation, Relationship, Race, Sex, HoursPerWeek, NativeCountry, Age.
- Age binned `<21, 21-30, …, 61-70, >70`; HoursPerWeek binned `≤20, 21-30, 30-40, >40`
  (the paper's own bin list is what it is — `21-30` and `30-40` overlap at 30; pick a convention and
  note it).
- **Fig. 8a**: PC-stable + RFT on n = 1000 → densely connected skeleton (contrast with Fig. 1a, same
  n, stratified MI → near-empty).
- **Fig. 8b**: 10 subsamples per point, x-axis ≈ 200–800. F1 compares d-connected pairs in the learned
  CPDAG against pairs called "dependent" in the data by a chi-square effect size
  `RMSEA = sqrt((χ² − df)/(n·df)) > 0.05`.
- **Fig. 1b** is a table of p and df for six pairs given Z = {Age, Sex} — the df-inflation illustration.

**Underspecified, so exact numbers won't reproduce:** PC significance level (presumably 0.05), which
variables were declared ordinal vs categorical, the exact subsample sizes, and all RNG seeds. Aim for
qualitative reproduction of the skeleton density and the F1 ordering.

## 5. Environment notes (this cloud container)

- Network: `git clone` from GitHub works. **PyPI, UCI, and Hugging Face are all unreachable** — so
  `pip install pgmpy` fails here; clone the repo instead, or build on a machine with PyPI access.
  numpy 2.4.4 / pandas 3.0.2 / scikit-learn 1.8.0 / scipy 1.17.1 are preinstalled.
- **No R** in the container and none on the Mac VM, and no CRAN access — the exact R stack
  (`bnlearn`, `pcalg`, `VGAM`, `nnet`, `ranger`) is not a viable route. Python port is the path.
- Raw Adult data is mirrorable from GitHub: `jbrownlee/Datasets` → `adult-train.csv` (32,561 rows,
  matches the UCI train split used by pgmpy) and `adult-all.csv` (48,842). Headerless; column names
  in `adult.names`.

## Sources

- Paper: <https://ojs.aaai.org/index.php/AAAI/article/view/26436> · arXiv <https://arxiv.org/abs/2206.04356>
- pgmpy: <https://github.com/pgmpy/pgmpy> (`pgmpy/ci_tests/`, `pgmpy/causal_discovery/PC.py`, `pgmpy/datasets/adult.py`)
- dagitty: <https://github.com/jtextor/dagitty> (`r/R/ci-tests.R`)
- Adult mirror: <https://github.com/jbrownlee/Datasets>

---

## 6. Addendum — do the descendant methods have papers? (asked 24 Aug 2026)

**Short answer: no. Both descendants exist only as code.** The published record still describes the
AAAI-23 χ² statistic; the canonical-correlation variants were never written up.

### 6.1 The pgmpy trace family — unpublished

`pillai`, `hotelling_lawley`, `wilks_lambda`, `roys_largest_root` are the four classical MANOVA
statistics (sum of squared canonical correlations; sum of ρ²/(1−ρ²); product of (1−ρ²); largest ρ²),
each with the F-approximation from **Muller & Peterson (1984)**, *Practical Methods for Computing
Power in Testing the Multivariate General Linear Hypothesis*, CSDA 2(2):143–158 — a power-calculation
paper for the multivariate GLM, not a CI test. In pgmpy's own bibliography:

- `pillai_trace.py` cites `ankan_textor_2023` + `li_shepherd_2010` + `muller_peterson_1984`
- `hotelling_lawley.py`, `wilks_lambda.py`, `roys_largest_root.py` cite **only** `muller_peterson_1984`
- `generalized_cov.py` (residual cross-covariance determinant, permutation p-values) cites
  `ankan_textor_2023` + `phipson_smyth_2010` (the permutation-p-value correction)
- `gcm.py` cites `shah_peters_2020` — the one genuinely published statistic in the family

So three of them have no CI-test citation at all, and the two that do point back at the AAAI-23 paper
whose statistic they do not implement.

**Ankan says as much himself.** In his post *Causal Discovery with Mixed Data using pgmpy* (Medium),
the only citation given is the AAAI-23 paper, described as: pgmpy "implements a variant of this method
by utilizing the residualizing idea with a measure of association based on canonical correlations."
No preprint, no arXiv link.

**Their own later papers still describe the original test.** In Ankan & Textor, *Expert-In-The-Loop
Causal Discovery* (UAI 2025, PMLR v286:3465–3479 — not on arXiv), the test is introduced as
"a residualization-based CI test [Ankan and Textor, 2023] **that returns a chi-square distributed
test statistic**" — the AAAI-23 Q-statistic, not the F-approximated trace. The pgmpy software paper
(Ankan & Textor, JMLR 25(265):1–8, 2024) predates the variant and lists only "residualization test
(Ankan and Textor, 2022)" among PC's CI tests.

Nothing in the recent literature documents it either: the Aug 2026 TMLR survey *Conditional
Independence Tests for Constraint-Based Causal Discovery* (arXiv:2608.11156) does not discuss the test,
and its library table doesn't even credit pgmpy with a regression-based CI test.

*(Erratum worth knowing if you cite via pgmpy: its bib entry for `ankan_textor_2023` carries
`note = {arXiv:2306.01638}`, but this paper's preprint is **arXiv:2206.04356**.)*

### 6.2 dagitty's `cis.pillai` — unpublished, and not LS residuals

Documented only in `?localTests` (`r/man/localTests.Rd`): ordinals integer-coded, categoricals
dummy-coded dropping the most frequent level, residualized by **`lm`**, then `cancor` +
`CCP::p.asym(tstat="Pillai")`, effect size = RMS canonical correlation. No accompanying paper; it does
not use Li–Shepherd residuals at all.

### 6.3 What *is* published in the Li–Shepherd line (the useful citations)

The residual itself has a well-developed published lineage out of Vanderbilt — this, not the pgmpy
variants, is the literature to position against:

- **Li & Shepherd (2010)**, *Test of Association Between Two Ordinal Variables While Adjusting for
  Covariates*, JASA 105(490):612–620 — effectively the ancestor CI test for ordinal data; the AAAI-23
  proofs are explicitly adapted from it.
- **Li & Shepherd (2012)**, *A New Residual for Ordinal Outcomes*, Biometrika 99(2):473–480 — the LS
  residual `p̂(Y<y|z) − p̂(Y>y|z)`.
- **Shepherd, Li & Liu (2016)**, *Probability-Scale Residuals for Continuous, Discrete, and Censored
  Data*, Canadian J. Statistics 44:463–476 — the general PSR framework.
- **Liu, Shepherd, Wanga & Li (2018)**, *Covariate-Adjusted Spearman's Rank Correlation with
  Probability-Scale Residuals*, Biometrics 74:595–605 — fit models of X and Y on Z, take PSRs, correlate
  them. Structurally the same move as the paper's Q1, published eight years earlier, in the
  rank-correlation framing.
- R implementation of all of the above: **`PResiduals`** (Dupont, Horner, Li, Liu, Shepherd) —
  `presid()`, `partial_Spearman()`, `conditional_Spearman()`. Best oracle for unit-testing a port of
  the residual.

Adjacent published CI tests the paper itself positions against: **Shah & Peters (2020)**, AoS
48(3):1514–1538 (GCM — the paper's Q1 is its square, for the ordinal–ordinal case), and
**Petersen & Hansen (2021)**, JMLR 22(70):1–47 (partial copula — the paper shows the LS residual is a
discrete limit of it).
