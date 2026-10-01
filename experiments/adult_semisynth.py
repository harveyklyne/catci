"""Semi-synthetic UCI adult income: real predictors, a simulated response.

This is the categorical analogue of the 401k study in Klyne & Shah (2025, AoS
53(1)), where the real ``(X, Z)`` predictors were kept and only ``Y`` was
simulated, so that the target parameter was known exactly while the predictor
distribution stayed realistic. Here the target is not a scalar parameter but the
*truth of the null* ``X indep Y | Z``, so what has to be known exactly is whether
conditional independence holds -- and, under the alternative, by how much.

The population
--------------
The population is the **full 30,162-row empirical distribution** of the adult
data (all eleven variables), exactly as the 401k study treated its 9,915 rows as
the predictor distribution. Everything below is a statement about that finite
population, so every "true" quantity is a finite sum and is known exactly rather
than estimated.

The construction
----------------
Fix a *generating stratum* ``W`` -- a tuple of discrete conditioning variables,
by default ``{Age, Sex}``. Write ``w(i)`` for row ``i``'s cell. A replicate:

1. Draw ``n`` rows from the population. ``X`` and every ``Z`` variable are the
   **real, untouched** values on those rows.
2. Draw ``Y_i`` independently from the mixture kernel

       h_lam(y | x, w) = (1 - lam) * P(y | w)  +  lam * P(y | x, w)

   where both ``P(y | w)`` and ``P(y | x, w)`` are the population's own
   conditional tables.

``lam = 0`` gives ``Y_i ~ P(y | w(i))`` using randomness independent of ``X``, so

    Y  indep  (X, all other variables)  |  W,

and hence, by weak union, ``Y indep X | Z`` for **any** ``Z`` containing ``W``.
The null is exact by construction, not by approximation. ``lam = 1`` reproduces
the population's own joint law of ``(X, Y) | W`` -- the real dependence, not an
invented one. ``lam`` interpolates.

Why the mixture, and not a tilt
-------------------------------
Because it moves the dependence while holding both nuisances fixed. Taking
``Z = W``:

* ``f(x | w) = P(x | w)`` is untouched -- ``X`` is never modified.
* ``g(y | w) = sum_x P(x | w) h_lam(y | x, w) = P(y | w)`` for **every** ``lam``,
  since the mixture's second component averages back to ``P(y | w)``.

So ``lam`` is a pure conditional-dependence knob: the marginals a GCM-type test
has to estimate do not move as it turns. Moreover the departure from the null is
*exactly linear* in ``lam``,

    P_lam(x, y | w) - P(x | w) P(y | w) = lam * [P(x, y | w) - P(x | w) P(y | w)],

so any chi-square-type non-centrality scales as ``n * lam**2``, which makes the
power grid cheap to place -- see :meth:`Population.ncp_per_n`.

An exponential tilt has neither property: it perturbs ``g`` as it perturbs the
dependence, confounding "the test lost calibration" with "the test gained power".

Generating with ``W`` coarser than ``Z``
----------------------------------------
The weak-union argument above means ``W`` need not be all of ``Z``. Generating
under ``W = {Age, Sex}`` while *testing* with a rich ``Z`` keeps the null exact
and makes the nuisance estimation genuinely hard: ``f(x | Z)`` is the real,
unknown, high-dimensional conditional, and the learner must discover for itself
that ``g`` depends on only two of the ``Z`` variables. The cost is that for
``lam > 0`` the identity ``g = P(y | w)`` no longer holds exactly, because the
mixture's second component then averages over ``P(x | Z)`` rather than
``P(x | w)``. Both are still computable in closed form on a finite population;
:func:`true_propensities` returns whichever applies.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import adult

DEFAULT_W = ("Age", "Sex")


def _cell_index(codes: pd.DataFrame, names) -> tuple[np.ndarray, int]:
    """Map rows to a contiguous 0-based cell id over the product of ``names``."""
    idx = np.zeros(len(codes), dtype=np.int64)
    size = 1
    for name in names:
        k = int(codes[name].max())
        idx += (codes[name].to_numpy() - 1) * size
        size *= k
    _, inv = np.unique(idx, return_inverse=True)
    return inv, int(inv.max()) + 1


@dataclass(frozen=True)
class Population:
    """The finite population every "true" quantity below is defined against.

    ``joint[c]`` is the ``dx x dy`` table ``P(x, y | W = c)``; ``px[c]`` and
    ``py[c]`` are its margins. ``cell`` gives each of the ``N`` population rows
    its ``W`` cell. Because the population is finite and fully enumerated, these
    tables *are* the truth -- there is no estimation error in them by definition.
    """

    data: adult.AdultData
    x_name: str
    y_name: str
    w_names: tuple[str, ...]
    cell: np.ndarray  # (N,) cell id per population row
    joint: np.ndarray  # (n_cells, dx, dy)
    px: np.ndarray  # (n_cells, dx)
    py: np.ndarray  # (n_cells, dy)

    @property
    def dx(self) -> int:
        return self.joint.shape[1]

    @property
    def dy(self) -> int:
        return self.joint.shape[2]

    @property
    def N(self) -> int:
        return len(self.cell)

    @property
    def pw(self) -> np.ndarray:
        """``P(W = c)`` over the population."""
        return np.bincount(self.cell, minlength=len(self.joint)) / self.N

    def delta_real(self) -> np.ndarray:
        """``P(y | x, w) - P(y | w)``: the population's own departure from independence.

        Rows for an ``(x, w)`` combination that never occurs are undefined; they
        are set to zero, which is unreachable during sampling and keeps the
        result a valid direction (see :meth:`kernel`).
        """
        with np.errstate(invalid="ignore", divide="ignore"):
            cond = self.joint / self.px[:, :, None]
        cond[self.px == 0] = 0.0
        out = cond - self.py[:, None, :]
        out[self.px == 0] = 0.0
        return out

    def kernel(self, lam: float, delta: np.ndarray | None = None) -> np.ndarray:
        """``h(y | x, w) = P(y | w) + lam * delta``, an ``(n_cells, dx, dy)`` kernel.

        ``delta`` defaults to :meth:`delta_real`, which recovers the mixture of
        the module docstring exactly: ``P(y|w) + lam * [P(y|x,w) - P(y|w)]``.

        Any ``delta`` satisfying the two constraints of :func:`is_valid_direction`
        gives a kernel with the same guarantees -- ``f`` untouched, ``g`` fixed at
        ``P(y | w)`` for every ``lam``, and departure from independence exactly
        linear in ``lam``. That freedom is what :func:`project_out_visible` uses.
        """
        delta = self.delta_real() if delta is None else delta
        return self.py[:, None, :] + lam * delta

    def ncp_per_n(self, delta: np.ndarray | None = None) -> float:
        """Stratified chi-square non-centrality per unit ``n`` at ``lam = 1``.

        The non-centrality at general ``lam`` and sample size ``n`` is
        ``n * lam**2 * ncp_per_n()``, by the linearity noted in the module
        docstring. Use it to place a power grid without a pilot run.

        This is the non-centrality of the *stratified* test -- the one that pays
        ``(dx-1)(dy-1) * n_cells`` degrees of freedom. It is a fair summary of
        "how much dependence is there", not of "how much power will any given
        test have"; see :func:`visible_fraction`.
        """
        delta = self.delta_real() if delta is None else delta
        joint = self.px[:, :, None] * self.kernel(1.0, delta)
        pw = self.pw
        total = 0.0
        for c in range(len(self.joint)):
            E = np.outer(self.px[c], joint[c].sum(0))
            ok = E > 0
            total += pw[c] * np.sum((joint[c][ok] - E[ok]) ** 2 / E[ok])
        return float(total)


def pool_rare(data: adult.AdultData, name: str, min_expected: float, n: int) -> adult.AdultData:
    """Merge levels of ``name`` expected to appear fewer than ``min_expected`` times.

    At ``n = 1000`` the adult data's rarest Occupation level (``Armed-Forces``,
    9 rows in 30,162) has an expected count of 0.3, and Workclass's
    ``Without-pay`` 0.5. A level that is usually absent contributes a near-zero
    column to ``Sigma``, and ``ApproxChi`` divides by it -- so those coordinates
    are numerically degenerate rather than merely uninformative.

    Pooling them is a **coarsening of the population**, applied once, before
    anything is simulated. Everything the module docstring establishes -- the
    exactness of the ``lam = 0`` null, the invariance of ``g``, the linearity in
    ``lam`` -- is a statement about whatever population is handed to
    :func:`build_population`, so it survives unchanged. The alternative, dropping
    unobserved levels per replicate, would not: it makes the level set random and
    the "true" propensities replicate-dependent.

    The pooled levels become a single ``Other (...)`` level appended last.
    """
    counts = data.codes[name].value_counts().sort_index()
    rare = counts.index[(counts / len(data.codes)) * n < min_expected].to_numpy()
    if len(rare) == 0:
        return data

    var = data.variables[name]
    keep = [i for i in counts.index if i not in set(rare.tolist())]
    remap = {old: new for new, old in enumerate(keep, start=1)}
    for old in rare.tolist():
        remap[old] = len(keep) + 1

    codes = data.codes.copy()
    codes[name] = codes[name].map(remap).astype(np.int64)
    pooled_label = "Other (" + ", ".join(var.levels[i - 1] for i in rare.tolist()) + ")"
    levels = tuple(var.levels[i - 1] for i in keep) + (pooled_label,)
    variables = dict(data.variables)
    # Pooling destroys any ordering among the merged levels, so an ordinal
    # variable that loses interior levels can no longer claim to be ordinal.
    variables[name] = adult.Variable(var.name, var.code, var.kind, levels)
    return adult.AdultData(codes=codes, variables=variables)


def build_population(data: adult.AdultData, x_name: str, y_name: str,
                     w_names=DEFAULT_W) -> Population:
    """Enumerate the population's ``(X, Y) | W`` tables."""
    w_names = tuple(w_names)
    cell, n_cells = _cell_index(data.codes, w_names)
    x = data.codes[x_name].to_numpy() - 1
    y = data.codes[y_name].to_numpy() - 1
    dx, dy = data.n_levels(x_name), data.n_levels(y_name)

    joint = np.zeros((n_cells, dx, dy))
    np.add.at(joint, (cell, x, y), 1.0)
    joint /= joint.sum(axis=(1, 2), keepdims=True)
    return Population(data=data, x_name=x_name, y_name=y_name, w_names=w_names,
                      cell=cell, joint=joint, px=joint.sum(2), py=joint.sum(1))


# --------------------------------------------------------------------------- #
# Planting a direction that Ankan & Textor's Q1 cannot see
# --------------------------------------------------------------------------- #
# For a pair of variables both typed *ordinal*, Ankan & Textor reduce each to a
# single Li-Shepherd residual column (``ankan_textor.residual_matrix``), so their
# Q1 statistic is one number with one degree of freedom, and its population value
# is the single scalar
#
#     S = E[ r_X(X, W) * r_Y(Y, W) ],   r_X(x, w) = P(X < x | w) - P(X > x | w).
#
# A test whose entire signal is one linear functional is blind to every direction
# orthogonal to it. Since :meth:`Population.kernel` accepts an arbitrary ``delta``,
# we can plant a departure from independence that is exactly orthogonal to that
# functional -- ``S = 0`` identically, at every ``lam`` and every ``n`` -- while
# leaving the rest of the dependence, which catci's full-table criterion does see,
# untouched. This is the categorical analogue of the interaction setting (22) in
# Klyne & Shah (2025), which exists precisely to break the partially-linear method.
#
# Three constraints have to hold simultaneously, and the projection below is the
# orthogonal projection that enforces all three:
#
#   (C1)  sum_y delta(x, y, w) = 0              -- h is a probability mass function
#   (C2)  sum_x P(x|w) delta(x, y, w) = 0       -- g stays fixed at P(y|w)
#   (C3)  <G, delta> = 0                        -- Q1's population signal vanishes
#
# C1 and C2 are what make ``lam`` a pure dependence knob (module docstring); C3 is
# the new one. All three are linear, so they compose.


def _ls_scores(probs: np.ndarray) -> np.ndarray:
    """Li-Shepherd score ``P(< l | w) - P(> l | w)`` for each level, per cell.

    ``probs`` is ``(n_cells, d)``. This is the population version of
    ``ankan_textor.ls_residuals``: that function evaluates the score at each
    observation's own label, this one tabulates it at every level.
    """
    cum = np.cumsum(probs, axis=1)
    below = np.concatenate([np.zeros((len(probs), 1)), cum[:, :-1]], axis=1)
    above = 1.0 - cum
    return below - above


def visible_direction(pop: Population) -> np.ndarray:
    """``G*``: the only direction Ankan & Textor's Q1 can see, projected into C1/C2.

    ``G(x, y, w) = r_X(x, w) * r_Y(y, w)`` represents Q1's population signal as
    ``S = <G, delta>`` under the inner product

        <A, B> = sum_w P(w) sum_{x,y} P(x|w) A(x,y,w) B(x,y,w).

    Projecting ``G`` into the C1/C2 subspace does not change ``<G, delta>`` for any
    admissible ``delta`` -- the projection is orthogonal and ``delta`` already lies
    in that subspace -- but it makes ``G*`` itself admissible, so it can be
    subtracted from a direction without breaking C1 or C2.
    """
    rx = _ls_scores(pop.px)  # (n_cells, dx)
    ry = _ls_scores(pop.py)  # (n_cells, dy)
    return _project_c1_c2(pop, rx[:, :, None] * ry[:, None, :])


def _project_c1_c2(pop: Population, arr: np.ndarray) -> np.ndarray:
    """Orthogonal projection onto ``{C0} ∩ {C1} ∩ {C2}``, i.e. double centering on the support.

    Unweighted over ``y`` (C1 is a plain sum) and ``P(x|w)``-weighted over ``x``
    (C2 is a weighted sum). Each is the orthogonal projection for its constraint
    under the inner product of :func:`visible_direction`, and they act on
    different indices, so they commute and their product is the projection onto
    the intersection.

    C0 is the support constraint ``delta(x, y, w) = 0 wherever P(y | w) = 0``. It
    is needed because ``h = P(y|w) + lam * delta`` must stay non-negative: a level
    of ``Y`` unobserved in stratum ``w`` has no probability mass to absorb a
    perturbation, so any non-zero ``delta`` there drives ``max_lambda`` to zero.
    Several ``(Age, Sex)`` strata of the adult data do have empty Education and
    HoursPerWeek levels, so this is not hypothetical. Centering therefore runs
    over the *supported* ``y`` in each stratum rather than over all of them.
    """
    support = pop.py > 0  # (n_cells, dy)
    n_sup = support.sum(axis=1, keepdims=True)  # (n_cells, 1)

    out = np.where(support[:, None, :], arr, 0.0)
    out = out - (out.sum(axis=2, keepdims=True) / n_sup[:, None, :]) * support[:, None, :]
    out = out - np.einsum("cx,cxy->cy", pop.px, out)[:, None, :] * support[:, None, :]
    out[pop.px == 0] = 0.0  # unreachable rows: zero satisfies C1 and is inert in C2
    return out


def _inner(pop: Population, a: np.ndarray, b: np.ndarray) -> float:
    return float(np.einsum("c,cx,cxy,cxy->", pop.pw, pop.px, a, b))


def q1_signal(pop: Population, delta: np.ndarray) -> float:
    """Q1's population signal ``S`` at ``lam = 1`` -- the number to drive to zero."""
    rx = _ls_scores(pop.px)
    ry = _ls_scores(pop.py)
    return _inner(pop, rx[:, :, None] * ry[:, None, :], delta)


def visible_fraction(pop: Population, delta: np.ndarray) -> float:
    """Fraction of ``delta``'s energy lying in the one direction Q1 can see.

    The squared cosine between ``delta`` and ``G*``. Small values mean Q1 is
    already nearly blind to the *real* dependence, and that projecting the
    remainder out costs almost nothing.
    """
    g = visible_direction(pop)
    gg = _inner(pop, g, g)
    if gg <= 0:
        return 0.0
    return float(_inner(pop, g, delta) ** 2 / (gg * _inner(pop, delta, delta)))


def project_out_visible(pop: Population, delta: np.ndarray | None = None) -> np.ndarray:
    """Remove from ``delta`` the single component Ankan & Textor's Q1 responds to.

    The result still satisfies C1 and C2, so it is a valid direction for
    :meth:`Population.kernel` with every guarantee intact, and additionally has
    ``q1_signal == 0`` to machine precision. Defaults to blinding the population's
    own dependence, which keeps the alternative as close to the real data as an
    Ankan-&-Textor-invisible alternative can be.
    """
    delta = pop.delta_real() if delta is None else delta
    g = visible_direction(pop)
    gg = _inner(pop, g, g)
    if gg <= 0:
        return delta
    return delta - (_inner(pop, g, delta) / gg) * g


def planted_direction(pop: Population, shape: np.ndarray | str = "u_shape") -> np.ndarray:
    """A well-conditioned direction with ``q1_signal == 0`` exactly.

    :func:`project_out_visible` gets the orthogonality right but is useless in
    practice: subtracting ``G*`` from the real dependence leaves a direction whose
    mass sits on cells where ``P(y | w)`` is nearly zero, so ``max_lambda`` is
    ~0.004 and the attainable effect size is negligible. Constructing the
    direction instead fixes that, by carrying a ``P(y | w)`` factor:

        delta(x, y, w) = c * P(y | w) * u_w(x) * v_w(y).

    Then ``h = P(y|w) * (1 + lam * c * u_w(x) * v_w(y))``, so non-negativity is a
    bound on ``c * u * v`` alone and the scaling below guarantees ``max_lambda >= 1``
    (only the negative products bind, so it is usually larger -- 1.72 for
    Education x Income given {Age, Sex}). The constraints become conditions on
    the factors:

    * C0 (support) -- automatic, from the ``P(y | w)`` factor.
    * C1 (``h`` sums to 1) -- ``sum_y P(y|w) v_w(y) = 0``.
    * C2 (``g`` fixed) -- ``sum_x P(x|w) u_w(x) = 0``.
    * C3 (Q1 blind) -- ``sum_x P(x|w) u_w(x) r_X(x, w) = 0``, per stratum.

    C3 is imposed on the ``X`` side because the ``Y`` side may have no room: with
    ``dy = 2`` the space of ``P(y|w)``-centred functions of ``y`` is one
    dimensional and ``r_Y`` spans it, so no admissible ``v`` is orthogonal to it.
    Blinding through ``X`` works for any ``dy``, needs only ``dx >= 3``, and is the
    more interpretable half anyway -- it is a statement about which levels of ``X``
    group together.

    ``u_w`` is ``shape`` Gram-Schmidted against ``{1, r_X}`` in the ``P(x|w)`` inner
    product, separately in each stratum, so C2 and C3 hold exactly rather than on
    average. The default ``"u_shape"`` template is the second harmonic
    ``cos(2*pi*(x-1)/(dx-1))``, the canonical non-monotone contrast: it separates
    the middle levels of ``X`` from both extremes. That is a shape no single
    monotone score can represent, and a two-group split catci's merging search
    can. ``v_w`` is taken to be ``r_Y`` itself, which puts the whole ``Y``-side
    signal in the direction A&T would have been most sensitive to -- so the
    resulting blindness cannot be dismissed as an unlucky choice on that side.
    """
    dx, dy = pop.dx, pop.dy
    if dx < 3:
        raise ValueError(f"blinding through X needs dx >= 3, got dx={dx}")
    if isinstance(shape, str):
        if shape != "u_shape":
            raise ValueError(f"unknown shape {shape!r}")
        s = np.cos(2.0 * np.pi * np.arange(dx) / (dx - 1))
    else:
        s = np.asarray(shape, dtype=float)
        if s.shape != (dx,):
            raise ValueError(f"shape must have length dx={dx}, got {s.shape}")

    rx = _ls_scores(pop.px)  # (n_cells, dx)
    ry = _ls_scores(pop.py)  # (n_cells, dy)

    u = np.zeros((len(pop.joint), dx))
    for c in range(len(pop.joint)):
        w = pop.px[c]
        vec = s.copy()
        # Gram-Schmidt against 1 then r_X, in the <a, b> = sum_x P(x|w) a b product.
        for basis in (np.ones(dx), rx[c]):
            nrm = float(w @ (basis * basis))
            if nrm > 1e-14:
                vec = vec - (float(w @ (vec * basis)) / nrm) * basis
        u[c] = vec

    delta = pop.py[:, None, :] * (u[:, :, None] * ry[:, None, :])
    delta[pop.px == 0] = 0.0

    # Scale so that |lam * u * v| <= 1 for lam <= 1, hence max_lambda >= 1.
    reachable = pop.px > 0
    prod = np.abs(u[:, :, None] * ry[:, None, :])
    prod = np.where(reachable[:, :, None] & (pop.py > 0)[:, None, :], prod, 0.0)
    peak = prod.max()
    if peak > 0:
        delta = delta / peak
    return delta


def is_valid_direction(pop: Population, delta: np.ndarray, tol: float = 1e-10) -> dict:
    """Check C1, C2, and the non-negativity of ``P(y|w) + lam * delta`` at ``lam = 1``."""
    reachable = pop.px > 0
    c1 = np.abs(delta.sum(axis=2)[reachable]).max()
    c2 = np.abs(np.einsum("cx,cxy->cy", pop.px, delta)).max()
    return dict(
        c1=float(c1), c2=float(c2), ok=bool(c1 < tol and c2 < tol),
        max_lam=float(max_lambda(pop, delta)),
    )


def max_lambda(pop: Population, delta: np.ndarray) -> float:
    """Largest ``lam`` keeping ``P(y|w) + lam * delta`` a non-negative kernel.

    Projecting out the visible component can push a near-zero cell negative, so
    a planted direction generally admits a smaller ``lam`` range than the real
    one (which is valid up to ``lam = 1`` by construction). The power grid has to
    respect this bound.
    """
    reachable = pop.px > 0
    py = np.broadcast_to(pop.py[:, None, :], delta.shape)
    neg = (delta < 0) & reachable[:, :, None]
    if not neg.any():
        return np.inf
    return float(np.min(py[neg] / -delta[neg]))


@dataclass(frozen=True)
class Replicate:
    """One simulated dataset, plus the propensities that are true for it."""

    rows: np.ndarray  # indices into the population
    x: np.ndarray  # 1-based, real
    y: np.ndarray  # 1-based, simulated
    z: np.ndarray  # (n, len(z_names)), real
    z_names: tuple[str, ...]
    f_true: np.ndarray  # (n, dx) true P(X | Z)
    g_true: np.ndarray  # (n, dy) true P(Y | Z)
    lam: float


def draw(pop: Population, n: int, lam: float, rng: np.random.Generator,
         z_names=None, replace: bool = False,
         delta: np.ndarray | None = None) -> Replicate:
    """One replicate: real ``X, Z`` on ``n`` sampled rows, ``Y`` from ``h_lam``.

    ``z_names`` defaults to ``pop.w_names``. Any ``z_names`` **containing**
    ``pop.w_names`` leaves the ``lam = 0`` null exact (weak union); anything else
    does not, and is rejected.

    ``delta`` defaults to the population's own dependence. Pass
    :func:`project_out_visible` to plant one Ankan & Textor's Q1 cannot see.
    """
    z_names = tuple(pop.w_names if z_names is None else z_names)
    missing = set(pop.w_names) - set(z_names)
    if missing:
        raise ValueError(
            f"z_names must contain the generating stratum {pop.w_names}; missing "
            f"{sorted(missing)}. Otherwise conditional independence at lam=0 does "
            "not hold for the Z being tested."
        )
    delta = pop.delta_real() if delta is None else delta

    rows = rng.choice(pop.N, size=n, replace=replace)
    if not replace:
        rows = np.sort(rows)
    cell = pop.cell[rows]
    x = pop.data.codes[pop.x_name].to_numpy()[rows]

    kernel = pop.kernel(lam, delta)  # (n_cells, dx, dy)
    probs = kernel[cell, x - 1, :]  # (n, dy)
    if probs.min() < -1e-12:
        raise ValueError(
            f"lam={lam} makes the kernel negative; max_lambda for this delta is "
            f"{max_lambda(pop, delta):.4f}."
        )
    u = rng.random(n)[:, None]
    y = (u > np.clip(probs, 0.0, None).cumsum(1)).sum(1) + 1

    f_true, g_true = true_propensities(pop, rows, lam, z_names, delta)
    z = pop.data.codes[list(z_names)].to_numpy()[rows]
    return Replicate(rows=rows, x=x, y=y, z=z, z_names=z_names,
                     f_true=f_true, g_true=g_true, lam=lam)


def true_propensities(pop: Population, rows: np.ndarray, lam: float,
                      z_names: tuple[str, ...],
                      delta: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """``P(X | Z)`` and ``P(Y | Z)`` on the given rows, exactly.

    Both are population quantities: ``f`` is the empirical conditional of the
    real ``X`` given the real ``Z``, and ``g`` is obtained by averaging the
    generating kernel over that same ``f``. When ``z_names == pop.w_names`` this
    reduces to ``g = P(y | w)`` for every ``lam`` (see the module docstring); the
    general form is computed here so that a richer ``Z`` is handled too.
    """
    zcell, n_z = _cell_index(pop.data.codes, z_names)
    x_all = pop.data.codes[pop.x_name].to_numpy() - 1

    # f(x | z): population counts over Z cells.
    fx = np.zeros((n_z, pop.dx))
    np.add.at(fx, (zcell, x_all), 1.0)
    fx /= fx.sum(1, keepdims=True)

    # g(y | z) = sum_x f(x | z) * h_lam(y | x, w). W is a coarsening of Z -- draw()
    # enforces that -- so every Z cell sits inside exactly one W cell.
    kernel = pop.kernel(lam, pop.delta_real() if delta is None else delta)
    w_of_z = np.zeros(n_z, dtype=np.int64)
    w_of_z[zcell] = pop.cell
    gy = np.einsum("zx,zxy->zy", fx, kernel[w_of_z])

    return fx[zcell[rows]], gy[zcell[rows]]
