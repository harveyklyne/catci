# MLP propensity learner

Notes on the `MLP classifier instead of xgboost` item. Code is merged
(`mlp_learner` in `src/catci/learners.py`, tuner in `experiments/tune.py`).

## TL;DR

The MLP is a **better propensity estimator** than the tuned XGBoost on 4 of 5
marginal settings (27-72% lower `E_f`), and never worse. But it **does not matter
for the test**: both learners already satisfy the paper's Assumption 1 with ~1000x
margin (`n*E_f*E_g ~ 0.001`), size is identical across oracle / xgb / mlp, and
paired power differences are within noise (tables below).

Decision: **keep xgboost as the default.** The MLP's value is as evidence that
the propensity learner is not what limits the method at `d = 8`. The one open
place it might matter is **sweeping `d`** (see the end).

## Gotchas

1. **Tune both learners.** `experiments/tune.py` ports the R `tune_xgb` protocol
   (fresh train/test pair per rep, held-out mlogloss, argmin) for both. Comparing
   a tuned xgboost against a default MLP is meaningless: the untuned MLP was ~10x
   *worse* on `E_f`, and overfits hard in-sample. The tuner warns when the argmin
   sits on a grid boundary. `--write` merges the winner into
   `experiments/tuning/n1000_numclass8/*.json` (both learners coexist per file).
2. **`predict_proba` drops absent classes.** `MLPClassifier.predict_proba` returns
   one column per class it *saw*; if a level is absent from the training labels
   the output is `(n, d-1)` and later columns shift left, silently giving wrong
   propensities. `mlp_learner` scatters through `clf.classes_` into a full
   `num_class`-wide matrix. Pinned by `test_mlp_pads_absent_classes_to_num_class`.
3. **`fork` -> `spawn` on macOS.** `run.py` used `mp.get_context("fork")`; sklearn
   starts BLAS threads and a forked child that inherits them aborts in the macOS
   Objective-C runtime (`+[NSString initialize] may have been in progress...`),
   intermittently (second pool in a process). All pools now use `spawn`. This was
   a latent bug in `run.py` that a threaded learner surfaced.
4. **Input scaling.** The net standardises `Z` internally; without it the MLP
   underfits toward the marginal, indistinguishable from "the MLP is bad".
   `test_mlp_beats_the_constant_predictor_when_z_is_informative` is the floor.
5. **`E_f` is the right metric.** Assumption 1 (paper p. 9) needs
   `E_f := max_j E[(f_j(Z) - fhat_j(Z))^2 | D] = o_P(1)`, likewise `E_g`, and
   `E_f*E_g = o_P(n^-1)`. Strength 0 suffices to measure it (the interaction
   matrices have zero row/column sums, so the true `X|Z` propensity is identical
   under null and alternative). In-sample vs out-of-sample `E_f` ratio is ~1.00
   for both tuned learners, so skipping cross-fitting costs nothing here (a
   property of these hyperparameters, not a general guarantee).

## Results

### Tuning outcome

`hidden_layer_sizes=(8,8), alpha=3.0, activation=tanh, max_iter=400` won on all
five settings (heavy shrinkage: these propensities carry little signal about
`Z`). Held-out mlogloss, 20 reps, with the fraction of the oracle->uniform gap
recovered:

| setting | oracle | xgb | mlp | uniform | xgb recovers | mlp recovers |
|---|---|---|---|---|---|---|
| lin | 2.0508 | 2.0691 | **2.0652** | 2.0794 | 36% | **50%** |
| vee | 2.0390 | 2.0593 | **2.0515** | 2.0794 | 50% | **69%** |
| hat | 2.0089 | 2.0349 | **2.0215** | 2.0794 | 63% | **82%** |
| sin | 2.0631 | 2.0832 | **2.0796** | 2.0794 | **-23%** | 0% |
| sig | 2.0637 | 2.0791 | **2.0753** | 2.0794 | 2% | **26%** |

On `sin`, tuned xgboost scores worse than predicting uniform; the MLP merely
breaks even.

### Propensity quality, `E_f`, 100 reps, n=1000, d=8

| setting | `E_f` xgb | `E_f` mlp | ratio | `n*E_f*E_g` xgb | `n*E_f*E_g` mlp |
|---|---|---|---|---|---|
| hat | 0.00133 | **0.00077** | **1.72x** | 0.0017 | 0.0006 |
| vee | 0.00114 | **0.00076** | **1.51x** | 0.0013 | 0.0006 |
| lin | 0.00106 | **0.00077** | **1.37x** | 0.0011 | 0.0006 |
| sig | 0.00092 | **0.00073** | **1.27x** | 0.0008 | 0.0005 |
| sin | **0.00112** | 0.00113 | 0.99x | 0.0013 | 0.0013 |

Fit time per call is comparable (xgb 0.14-0.95 s, mlp 0.57-0.82 s). Assumption 1
asks `n*E_f*E_g -> 0`; it is already ~1000x inside the requirement for both.

### Size, minP master (`3770b79`), `n_boot = 1000`, 1000 reps, SE ~0.007

| method | lin oracle | lin xgb | lin mlp | sin oracle | sin xgb | sin mlp |
|---|---:|---:|---:|---:|---:|---:|
| tree | 0.049 | 0.048 | 0.052 | 0.041 | 0.049 | 0.044 |
| ordinal | 0.062 | 0.057 | 0.061 | 0.057 | 0.052 | 0.053 |
| max | 0.043 | 0.043 | 0.045 | 0.043 | 0.041 | 0.042 |
| euclid | 0.041 | 0.040 | 0.042 | 0.048 | 0.050 | 0.048 |
| mGCM | 0.174 | 0.190 | 0.186 | 0.189 | 0.187 | 0.185 |
| chi_sq | 0.161 | 0.165 | 0.163 | 0.141 | 0.143 | 0.142 |

mGCM and chi_sq fail by design at `d_X d_Y = 64`, `n = 1000` (see
`SEMISYNTH_RESULTS.md`, background findings, for the mGCM mechanism).

### Power: the learner does not matter

Paired (same simulated datasets and bootstrap normals; the learner is not in the
seed stream), 200 reps x 5 strengths. Mean power over the grid, xgb -> mlp:

| setting | method | xgb | mlp | diff |
|---|---|---:|---:|---:|
| lin_lin_step | ordinal | 0.623 | 0.621 | -0.002 |
| lin_lin_step | euclid | 0.484 | 0.490 | +0.006 |
| sin_sin_binary_tree | tree | 0.821 | 0.813 | -0.008 |
| sin_sin_binary_tree | euclid | 0.665 | 0.659 | -0.006 |

No per-strength paired difference exceeds 0.03; one cell of ~70 is > 2 SE
(mGCM, s = 1.0, lin), consistent with chance.

## Only open question

**Sweep `d`.** `bench_learners.py` takes `--d 4 8 16 32` and `tune.py` can tune
either learner at any `d`. The MLP shares a hidden representation across classes
where boosting fits effectively one-vs-rest, so its advantage plausibly grows
with `d` (motivating case `dX=30, dY=10`). This is the only place the MLP might
matter. Both learners must be tuned at each `d` first; xgboost tuning at `d=32`
is the expensive part, so consider lowering `XGB_MAXROUNDS` equally for both.
