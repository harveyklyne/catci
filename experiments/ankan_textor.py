"""Ankan & Textor (AAAI-23): the residualization CI test, ported from the paper.

There is no replication package for the paper, and the two descendant
implementations shipped by the authors are *not* this statistic -- pgmpy's
``pillai`` and dagitty's ``cis.pillai`` both replace the chi-square quadratic
form with a canonical-correlation trace under an F-approximation, and dagitty's
drops the Li-Shepherd residual entirely. So this is a port from the paper text
(pp. 12182-12183), not a wrapper.

The method has two independent parts.

**The residual** (Li & Shepherd 2012). For an ordinal variable,

    R_i = p_hat(Y < y_i | z_i) - p_hat(Y > y_i | z_i)

For a binary variable this collapses to ``y_i - p_hat(Y = 1 | z_i)``, which the
paper notes is the ordinary observed-minus-expected residual. A categorical
variable with k levels is handled by applying that binary form to each of its
first k-1 dummy indicators, giving ``I(x_i = j) - p_hat(X = j | z_i)``.

**The statistic.** The paper gives three, but they are one statistic at three
shapes. Let ``Rx`` be an ``(n, a)`` residual matrix and ``Ry`` an ``(n, b)`` one,
where a variable contributes ``a = 1`` column if ordinal and ``a = k - 1``
columns if categorical. Form the ``(n, a*b)`` matrix ``M`` of elementwise
products, and let ``d = M.sum(0)`` be the vector of dot products. Then

    Q = (1 / n) * d @ inv(Sigma_d) @ d.T,    Sigma_d = cov(M),    df = a * b

recovers all three propositions exactly:

===========  =============  =============  ==========  ==================
statistic    X              Y              df          paper
===========  =============  =============  ==========  ==================
Q1           ordinal        ordinal        1           chi2(1)
Q2           categorical    ordinal        k-1         chi2(k-1)
Q3           categorical    categorical    (k-1)(r-1)  chi2((k-1)(r-1))
===========  =============  =============  ==========  ==================

Q1 is also the squared generalized covariance measure of Shah & Peters (2020),
which :func:`q_statistic` reproduces as the ``a = b = 1`` case -- one of the
port's unit tests.

Dropping the last dummy is not cosmetic: with all k indicators the columns of
``M`` are linearly dependent and ``Sigma_d`` is singular.

Conditioning
------------
Dropping the dummy is necessary but not sufficient. When a variable has many
levels relative to the sample, ``Sigma_d`` goes rank-deficient anyway: if level
``l`` of X and level ``m`` of Y never co-occur within any stratum of Z, the
product column ``R_I(x=l) * R_I(y=m)`` is identically zero (or numerically so),
and the quadratic form is undefined. On adult income with ``Z = {Age, Sex}``
(14 strata) this bites on every pair involving ``NativeCountry`` (41 levels).

The paper does not discuss this. It matters because the failure is *silent* under
a plain solve: ``np.linalg.solve`` raises only on exact singularity, so a
condition number of 1e21 returns a number rather than an error. :func:`q_statistic`
therefore always reports ``rank`` and ``cond`` alongside the statistic, and
callers should treat ``well_conditioned=False`` results as undefined rather than
as evidence. Passing ``method="pinv"`` substitutes the Moore-Penrose inverse with
``df = rank(Sigma_d)``, which is a deviation from the paper (it changes the null)
but is what makes the high-cardinality pairs computable at all.

Conditional probability models
------------------------------
The paper uses ``nnet::multinom`` for categorical targets, ``VGAM`` proportional
odds for ordinal targets, and a 50-tree ``ranger`` probability forest for RFT.
The Python equivalents here are :class:`GLM` (scikit-learn multinomial logistic /
statsmodels ``OrderedModel``) and :class:`RFT`
(``RandomForestClassifier(n_estimators=50)``).

Choices the paper leaves open are marked ``NOTE:`` below.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np
from scipy.stats import chi2

Kind = Literal["ordinal", "categorical"]

__all__ = [
    "GLM", "RFT", "ls_residuals", "residual_matrix", "q_statistic",
    "fit_residuals", "test_from_residuals", "ci_test", "TestResult",
    "QResult", "COND_LIMIT",
]


# --------------------------------------------------------------------------- #
# Design matrix for the conditioning set
# --------------------------------------------------------------------------- #
def design_matrix(z: np.ndarray, n_levels: list[int], drop_first: bool = True) -> np.ndarray:
    """One-hot encode categorical/ordinal conditioning variables.

    NOTE: the paper does not state how Z is encoded. R would treatment-code
    unordered factors and use polynomial contrasts for ordered ones; pgmpy
    one-hot encodes everything. We one-hot encode everything and drop the first
    level of each variable for identifiability, which matches R's treatment
    coding for unordered factors and is the more conservative reading for
    ordered ones (it imposes no shape on the Z effect).
    """
    z = np.asarray(z)
    if z.ndim == 1:
        z = z[:, None]
    blocks = []
    for j, k in enumerate(n_levels):
        oh = np.zeros((z.shape[0], k))
        oh[np.arange(z.shape[0]), z[:, j] - 1] = 1.0
        blocks.append(oh[:, 1:] if drop_first else oh)
    return np.hstack(blocks) if blocks else np.zeros((z.shape[0], 0))


# --------------------------------------------------------------------------- #
# Conditional probability models: (design, labels, n_levels) -> (n, n_levels)
# --------------------------------------------------------------------------- #
class ProbabilityModel(Protocol):
    def fit_predict(self, design: np.ndarray, labels: np.ndarray, n_levels: int, kind: Kind) -> np.ndarray: ...


@dataclass(frozen=True)
class GLM:
    """The paper's GLM estimator.

    Categorical targets: multinomial logistic regression (``nnet::multinom``).
    Ordinal targets: proportional-odds logistic regression (``VGAM``), which
    takes the level order into account.

    NOTE: ``nnet::multinom`` is unpenalised by default (``decay = 0``).
    scikit-learn always penalises, so ``C`` is set very large to approximate the
    unpenalised MLE while keeping the fit numerically stable on rare levels.
    """

    C: float = 1e6
    max_iter: int = 5000

    def fit_predict(self, design, labels, n_levels, kind):
        if n_levels < 2:
            raise ValueError("a variable needs at least 2 observed levels")
        if kind == "ordinal" and n_levels > 2:
            return self._proportional_odds(design, labels, n_levels)
        return self._multinomial(design, labels, n_levels)

    def _multinomial(self, design, labels, n_levels):
        from sklearn.linear_model import LogisticRegression

        clf = LogisticRegression(C=self.C, max_iter=self.max_iter, tol=1e-8)
        clf.fit(design, labels)
        return _expand_to_levels(clf.predict_proba(design), clf.classes_, n_levels)

    def _proportional_odds(self, design, labels, n_levels):
        from statsmodels.miscmodels.ordinal_model import OrderedModel

        # OrderedModel supplies its own thresholds, so the design must not carry
        # an intercept column; drop_first in design_matrix already ensures that.
        model = OrderedModel(labels, design, distr="logit")
        res = model.fit(method="bfgs", disp=False, maxiter=self.max_iter)
        proba = np.asarray(res.model.predict(res.params, exog=design))
        return _expand_to_levels(proba, np.unique(labels), n_levels)


@dataclass(frozen=True)
class RFT:
    """The paper's random-forest estimator: ``ranger`` probability forest, 50 trees.

    NOTE: order is not passed to the forest -- ``ranger`` treats the target as a
    factor for probability estimation regardless, so ordinal and categorical
    targets are fitted identically. The order still enters through the LS
    residual, which is where it matters.

    NOTE: the paper gives no seed. ``random_state`` is exposed so a run is
    reproducible; results will not match the paper's RNG.
    """

    n_estimators: int = 50
    random_state: int = 0
    n_jobs: int = 1

    def fit_predict(self, design, labels, n_levels, kind):
        from sklearn.ensemble import RandomForestClassifier

        clf = RandomForestClassifier(
            n_estimators=self.n_estimators,
            random_state=self.random_state,
            n_jobs=self.n_jobs,
        )
        clf.fit(design, labels)
        return _expand_to_levels(clf.predict_proba(design), clf.classes_, n_levels)


def _expand_to_levels(proba: np.ndarray, classes: np.ndarray, n_levels: int) -> np.ndarray:
    """Place a model's per-class columns into a full ``(n, n_levels)`` matrix.

    A level absent from the fitted labels gets a probability of exactly 0 rather
    than shifting the remaining columns left, so column j always means level j+1.
    """
    full = np.zeros((proba.shape[0], n_levels))
    full[:, np.asarray(classes, dtype=int) - 1] = proba
    return full


# --------------------------------------------------------------------------- #
# Residuals
# --------------------------------------------------------------------------- #
def ls_residuals(labels: np.ndarray, proba: np.ndarray) -> np.ndarray:
    """Li-Shepherd residual ``p(Y < y_i) - p(Y > y_i)`` for 1-based ``labels``."""
    labels = np.asarray(labels)
    n, k = proba.shape
    cum = np.cumsum(proba, axis=1)
    rows = np.arange(n)
    below = np.where(labels > 1, cum[rows, np.clip(labels - 2, 0, k - 1)], 0.0)
    above = 1.0 - cum[rows, labels - 1]
    return below - above


def residual_matrix(labels: np.ndarray, proba: np.ndarray, kind: Kind) -> np.ndarray:
    """The residual columns a variable contributes to the statistic.

    Ordinal -> one LS-residual column. Categorical with k levels -> k-1 columns
    of ``I(x = j) - p_hat(X = j | z)``, the binary LS residual per dummy. The
    last dummy is dropped to keep ``Sigma_d`` full rank.
    """
    labels = np.asarray(labels)
    if kind == "ordinal":
        return ls_residuals(labels, proba)[:, None]
    k = proba.shape[1]
    dummies = np.zeros((labels.shape[0], k))
    dummies[np.arange(labels.shape[0]), labels - 1] = 1.0
    return (dummies - proba)[:, : k - 1]


# --------------------------------------------------------------------------- #
# The statistic
# --------------------------------------------------------------------------- #
# A condition number above this means the solve is not meaningfully invertible:
# float64 carries ~16 digits, so 1e12 already loses most of them.
COND_LIMIT = 1e12


@dataclass(frozen=True)
class QResult:
    """The statistic plus the diagnostics needed to know whether to believe it."""

    statistic: float
    df: int
    rank: int
    cond: float
    n_columns: int

    @property
    def well_conditioned(self) -> bool:
        return self.rank == self.n_columns and self.cond < COND_LIMIT


def q_statistic(rx: np.ndarray, ry: np.ndarray, method: str = "solve") -> QResult:
    """``Q = (1/n) d Sigma_d^-1 d^T`` for residual matrices ``rx`` (n, a), ``ry`` (n, b).

    ``d`` is the vector of dot products between residual columns and ``Sigma_d``
    the covariance of their elementwise products. This is Q1/Q2/Q3 of the paper
    at ``(a, b) = (1, 1) / (k-1, 1) / (k-1, r-1)``.

    ``method="solve"`` is the paper's form, with ``df = a*b``. ``method="pinv"``
    uses the pseudo-inverse and sets ``df = rank(Sigma_d)``, which keeps
    rank-deficient pairs computable at the cost of departing from the stated null.
    """
    n = rx.shape[0]
    # column (l, m) of M is Rx[:, l] * Ry[:, m]; ordering is Ry-major, matching
    # the paper's d = (R_I(x=1).R_I(y=1), ..., R_I(x=k-1).R_I(y=1), ...)
    m = (rx[:, None, :] * ry[:, :, None]).reshape(n, -1)
    d = m.sum(axis=0)
    sigma = np.atleast_2d(np.cov(m, rowvar=False, ddof=1))

    rank = int(np.linalg.matrix_rank(sigma))
    cond = float(np.linalg.cond(sigma))
    n_columns = m.shape[1]

    if method == "pinv":
        stat = float(d @ np.linalg.pinv(sigma) @ d / n)
        df = rank
    elif method == "solve":
        stat = float(d @ np.linalg.solve(sigma, d) / n)
        df = n_columns
    else:
        raise ValueError(f"method must be 'solve' or 'pinv', got {method!r}")

    return QResult(statistic=stat, df=df, rank=rank, cond=cond, n_columns=n_columns)


@dataclass(frozen=True)
class TestResult:
    x: str
    y: str
    statistic: float
    df: int
    p_value: float
    which: str  # "Q1" | "Q2" | "Q3"
    rank: int
    cond: float
    well_conditioned: bool


def _which_q(kind_x: Kind, kind_y: Kind) -> str:
    if kind_x == "ordinal" and kind_y == "ordinal":
        return "Q1"
    if kind_x == "categorical" and kind_y == "categorical":
        return "Q3"
    return "Q2"


def fit_residuals(
    labels: np.ndarray,
    z_design: np.ndarray,
    kind: Kind,
    n_levels: int,
    model: ProbabilityModel,
) -> np.ndarray:
    """Fit ``p_hat(V | Z)`` and turn it into V's residual columns.

    Split out from :func:`ci_test` because a residual depends only on the
    variable and Z, never on the variable it is tested against. A sweep over many
    pairs that share one Z should fit each variable once and reuse the result --
    see ``run_adult.py``, where doing so is an 8x saving.
    """
    proba = model.fit_predict(z_design, labels, n_levels, kind)
    return residual_matrix(labels, proba, kind)


def test_from_residuals(
    rx: np.ndarray,
    ry: np.ndarray,
    kind_x: Kind,
    kind_y: Kind,
    name_x: str = "X",
    name_y: str = "Y",
    method: str = "solve",
) -> TestResult:
    """The statistic and its p-value, given residuals already fitted."""
    q = q_statistic(rx, ry, method=method)
    return TestResult(
        x=name_x, y=name_y, statistic=q.statistic, df=q.df,
        p_value=float(chi2.sf(q.statistic, q.df)),
        which=_which_q(kind_x, kind_y),
        rank=q.rank, cond=q.cond, well_conditioned=q.well_conditioned,
    )


def ci_test(
    x: np.ndarray,
    y: np.ndarray,
    z_design: np.ndarray,
    kind_x: Kind,
    kind_y: Kind,
    n_levels_x: int,
    n_levels_y: int,
    model: ProbabilityModel,
    name_x: str = "X",
    name_y: str = "Y",
    method: str = "solve",
) -> TestResult:
    """Test ``X _||_ Y | Z`` by the paper's method.

    ``z_design`` is the already-encoded conditioning matrix (see
    :func:`design_matrix`); it is passed in rather than built here so that a
    sweep over many pairs encodes Z once.
    """
    rx = fit_residuals(x, z_design, kind_x, n_levels_x, model)
    ry = fit_residuals(y, z_design, kind_y, n_levels_y, model)
    return test_from_residuals(rx, ry, kind_x, kind_y, name_x, name_y, method)
