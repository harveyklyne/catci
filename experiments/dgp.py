"""Data-generating processes for the catci simulations.

Rewritten port of the R ``simulate_data`` and its helpers, with the review's
DGP fixes applied:

* **2a** -- the R lin/vee/hat branch drew ``x, y`` via a ``u/v`` mixture and then
  *unconditionally overwrote* them from a second block with a different
  dependence scaling, so the first draws (and one of the two scalings) were dead
  code. Here there is a single sampling path: compute the marginals ``f, g``,
  form the joint ``dep = indep + scale * interaction`` once, and sample from it.
* **2b** -- each variable's marginal now depends only on its own setting, so
  cross-family pairs (e.g. ``x="lin", y="sin"``) no longer crash.
* **2c** -- nothing is hardcoded to ``d = 8``: sin/sig breakpoints come from
  ``norm.ppf`` on an equally-spaced grid, and the binary-tree interaction is
  built from tree contrasts for any power-of-two ``dx``, ``dy`` (not
  necessarily equal).
* The ``step`` interaction no longer vanishes for ``dy in {2, 3}``.

This module lives in ``experiments/`` (not the ``catci`` method library): the
method knows nothing about how data is simulated.

Index convention matches the method library: the length ``dx*dy`` joint vector
is ordered ``(1,1),(2,1),...,(dx,1),(1,2),...`` -- X fastest.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

__all__ = [
    "simulate_data",
    "simulate_marginal",
    "z_correlated_normal",
    "lin_pdf",
    "vee_pdf",
    "hat_pdf",
    "sin_pdf",
    "sig_pdf",
    "get_pdf",
    "marginal_matrix",
    "get_int",
    "binary_tree_interaction",
    "permute_labels",
]

_MARGINAL_MIX = ("lin", "vee", "hat")
_MARGINAL_Z = ("sin", "sig")


# --------------------------------------------------------------------------- #
# Z
# --------------------------------------------------------------------------- #
def z_correlated_normal(n: int, p: int, corr: float, rng: np.random.Generator) -> np.ndarray:
    """``Z ~ N_p(0, Sigma)`` with unit variances and all off-diagonals ``corr``."""
    Sigma = np.full((p, p), corr)
    np.fill_diagonal(Sigma, 1.0)
    L = np.linalg.cholesky(Sigma)
    return rng.standard_normal((n, p)) @ L.T


# --------------------------------------------------------------------------- #
# Marginal pmfs (parametric in d)
# --------------------------------------------------------------------------- #
def lin_pdf(d: int) -> np.ndarray:
    p = np.linspace(1.0, 3.0, d)
    return p / p.sum()


def vee_pdf(d: int) -> np.ndarray:
    half = np.linspace(1.0, 3.0, int(np.ceil(d / 2)))
    p = np.concatenate([half, half[: d // 2][::-1]])
    return p / p.sum()


def hat_pdf(d: int) -> np.ndarray:
    p = np.ones(d)
    p[d // 3 : int(np.ceil(2 * d / 3))] = 3.0  # middle third raised
    return p / p.sum()


def _pdf_from_breakpoints(f0: np.ndarray, d: int, center: float, std: float) -> np.ndarray:
    """Bin a latent ``f0`` into ``d`` ordinal classes using ``norm.ppf`` breakpoints.

    Interior breakpoints ``center + std * ppf(k/d)`` make the classes equiprobable
    when ``f0 == center`` and generalise to any ``d`` (fix 2c; supersedes the
    hand-tuned length-9 vectors the R code hardcoded for ``d = 8``).
    """
    f0 = np.asarray(f0, dtype=float)
    n = f0.shape[0]
    interior = center + std * norm.ppf(np.arange(1, d) / d)  # length d-1
    z = (interior[None, :] - f0[:, None]) / std
    cdf = norm.cdf(z)
    cdf = np.hstack([np.zeros((n, 1)), cdf, np.ones((n, 1))])
    return cdf[:, 1:] - cdf[:, :-1]


def sin_pdf(z: np.ndarray, d: int, std: float = 2.0) -> np.ndarray:
    f0 = np.exp(-np.asarray(z) ** 2 / 2) * np.sin(np.asarray(z))  # sinusoidal, a=1
    return _pdf_from_breakpoints(f0, d, center=0.0, std=std)


def sig_pdf(z: np.ndarray, d: int, std: float = 2.0) -> np.ndarray:
    f0 = 1.0 / (1.0 + np.exp(-3.0 * np.asarray(z)))  # sigmoidal, s=3, in (0,1)
    return _pdf_from_breakpoints(f0, d, center=0.5, std=std)


def get_pdf(setting: str, d: int, z=None, modify: bool = False) -> np.ndarray:
    if setting == "lin":
        out = lin_pdf(d)
    elif setting == "vee":
        out = vee_pdf(d)
    elif setting == "hat":
        out = hat_pdf(d)
    elif setting == "sin":
        out = sin_pdf(z, d)
    elif setting == "sig":
        out = sig_pdf(z, d)
    else:
        raise ValueError(f"Setting not recognised: {setting!r}")
    if modify:
        out = 2.0 / d - out
    return out


def marginal_matrix(setting: str, d: int, z_col: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """The true conditional propensity matrix ``(n, d)`` for one variable (fix 2b).

    Depends only on this variable's own ``setting`` and ``z_col``.
    """
    z_col = np.asarray(z_col, dtype=float)
    if setting in _MARGINAL_MIX:
        base = get_pdf(setting, d, modify=False)
        mod = get_pdf(setting, d, modify=True)
        q = 1.0 / (1.0 + np.exp(-3.0 * z_col))  # latent mixing weight from z
        return (1.0 - q)[:, None] * base[None, :] + q[:, None] * mod[None, :]
    if setting in _MARGINAL_Z:
        return get_pdf(setting, d, z=z_col, modify=False)
    raise ValueError(f"Setting not recognised: {setting!r}")


# --------------------------------------------------------------------------- #
# Interaction matrices (parametric in d), all with zero row/column sums
# --------------------------------------------------------------------------- #
def _get_int_step(dx: int, dy: int) -> np.ndarray:
    """X: a step (low half vs high half). Y: a vee over quarters, ``+ - (0) - +``.

    The vee needs ``dy >= 4``; below that the quarter blocks are empty and the R
    construction returned the zero matrix -- no signal at any strength. For
    ``dy in {2, 3}`` Y falls back to a single step, ``+ (0) -``, so the pattern
    is still a 2x2 block contrast.
    """
    if dx < 2 or dy < 2:
        raise ValueError("step interaction needs dx, dy >= 2.")
    outcol = [-1] * (dx // 2) + [0] * (dx % 2) + [1] * (dx // 2)
    rc = outcol[::-1]
    if dy < 4:
        outall = outcol + [0] * (dx * (dy % 2)) + rc
        return np.array(outall, dtype=float).reshape(dx, dy, order="F")
    outall = (
        outcol * (dy // 4)
        + rc * (dy // 4)
        + [0] * (dx * (dy % 4))
        + rc * (dy // 4)
        + outcol * (dy // 4)
    )
    return np.array(outall, dtype=float).reshape(dx, dy, order="F")


def _get_int_alt(dx: int, dy: int) -> np.ndarray:
    outcol = [-1, 1] * (dx // 2)
    if dx % 2 == 1:
        outcol = outcol[: dx // 2] + [0] + outcol[dx // 2 :]
    neg = [-v for v in outcol]
    outall = (outcol + neg) * (dy // 2)
    if dy % 2 == 1:
        cut = dx * (dy // 2)
        outall = outall[:cut] + [0] * dx + outall[cut:]
    return np.array(outall, dtype=float).reshape(dx, dy, order="F")


def _haar(d: int, ell: int) -> np.ndarray:
    """Level-``ell`` tree contrast on ``d = 2^L`` leaves: ``1^{(ell-1)} (x) h (x) 1^{(L-ell)}``, ``h = (1, -1)``.

    ``+1`` under the left child and ``-1`` under the right child of every
    depth-``(ell-1)`` node, so level 1 is low half vs high half.
    """
    L = int(round(np.log2(d)))
    return np.kron(np.kron(np.ones(2 ** (ell - 1)), [1.0, -1.0]), np.ones(2 ** (L - ell)))


def binary_tree_interaction(dx: int, dy: int | None = None) -> np.ndarray:
    """Hierarchical interaction over balanced binary trees on X and Y.

    ``T = - sum_ell w_ell a_ell^X (a_ell^Y)^T`` over tree levels
    ``ell = 1, ..., min(Lx, Ly)``, where ``a_ell`` is the level-``ell`` contrast
    (:func:`_haar`). The coarsest level has weight 1, finer levels decay
    geometrically (``0.1^{ell-1}``), and the finest shared level has weight 0 --
    on the smaller side those are sibling leaves, which should be exchangeable so
    that merging them loses nothing. Rows and columns sum to zero.

    For ``dx == dy`` this is the recursive Kronecker construction of fix 2c
    (``kron(T, J) - w_ell kron(J^{(ell-1)}, H)`` with ``H = h h^T``,
    ``J = 1 1^T``) written as a sum of outer products, and at ``d = 8`` it
    reproduces the hand-coded R matrix exactly. ``dx != dy`` needs both to be
    powers of two, ``>= 2``, with ``max(dx, dy) >= 4``.
    """
    dy = dx if dy is None else dy
    pow2 = lambda d: d >= 2 and (d & (d - 1)) == 0
    if not (pow2(dx) and pow2(dy)) or max(dx, dy) < 4:
        raise ValueError("binary_tree interaction needs dx, dy powers of two, >= 2, max(dx, dy) >= 4.")
    L = int(round(np.log2(min(dx, dy))))
    T = np.zeros((dx, dy))
    for ell in range(1, L + 1):
        w = 1.0 if ell == 1 else (0.0 if ell == L else 0.1 ** (ell - 1))
        T -= w * np.outer(_haar(dx, ell), _haar(dy, ell))
    return T


def get_int(setting: str, dx: int, dy: int) -> np.ndarray:
    if setting == "step":
        return _get_int_step(dx, dy)
    if setting == "alt":
        return _get_int_alt(dx, dy)
    if setting == "binary_tree":
        return binary_tree_interaction(dx, dy)
    raise ValueError(f"Setting not recognised: {setting!r}")


# --------------------------------------------------------------------------- #
# Label permutation (breaks ordinal structure)
# --------------------------------------------------------------------------- #
def permute_labels(x: np.ndarray, f: np.ndarray, rng: np.random.Generator):
    """Randomly permute the labels of ``x`` and the columns of ``f`` consistently."""
    d = f.shape[1]
    perm = rng.permutation(d)  # 0-based
    order = np.argsort(perm)
    x_new = perm[np.asarray(x) - 1] + 1
    return x_new, f[:, order]


# --------------------------------------------------------------------------- #
# Joint simulation (single sampling path -- fix 2a)
# --------------------------------------------------------------------------- #
def _sample_categorical(probs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Row-wise categorical draw; returns 0-based class indices."""
    cum = np.cumsum(probs, axis=1)
    u = rng.random(probs.shape[0])
    return (u[:, None] < cum).argmax(axis=1)


def simulate_marginal(n: int, d: int, setting: str, rng: np.random.Generator, column: int = 0) -> dict:
    """Draw ``(Z, X)`` from the law of one variable given ``Z`` in :func:`simulate_data`.

    Every interaction has zero row and column sums, so the joint in
    :func:`simulate_data` preserves both margins exactly: ``X | Z`` has pmf
    ``marginal_matrix(xsetting, dx, Z[:, 0])`` whatever ``strength``,
    ``intsetting`` or ``dy`` are (and ``Y | Z`` likewise with ``Z[:, 1]``).
    So a propensity learner can be tuned on this, without ever forming the
    ``dx * dy`` joint -- which is what makes tuning at large ``d`` cheap.
    """
    z = z_correlated_normal(n, 5, 0.5, rng)
    f = marginal_matrix(setting, d, z[:, column], rng)
    x = _sample_categorical(f, rng) + 1
    return {"z": z, "x": x, "f": f}


def simulate_data(
    n: int,
    dx: int,
    dy: int,
    xsetting: str,
    ysetting: str,
    strength: float,
    intsetting: str,
    permute: bool = False,
    rng: np.random.Generator | None = None,
) -> dict:
    """Generate an i.i.d. ``(X, Y, Z)`` dataset with the true propensities ``f, g``.

    ``strength`` scales the interaction; at ``strength = 0`` the data are
    conditionally independent. Returns a dict with ``x, y`` (1-based int arrays),
    ``z`` ``(n, 5)``, and ``f``/``g`` (the true ``P(.|Z)`` matrices).
    """
    if rng is None:
        rng = np.random.default_rng()

    z = z_correlated_normal(n, 5, 0.5, rng)
    f = marginal_matrix(xsetting, dx, z[:, 0], rng)  # (n, dx)
    g = marginal_matrix(ysetting, dy, z[:, 1], rng)  # (n, dy)

    indep = (g[:, :, None] * f[:, None, :]).reshape(n, dx * dy)  # X fastest
    interaction = get_int(intsetting, dx, dy).reshape(-1, order="F")  # X fastest
    dep = indep + (strength * indep.min() / 2.0) * interaction[None, :]

    if np.any(dep < -1e-12):
        raise ValueError("Negative probabilities: strength too large for this interaction.")
    dep = np.clip(dep, 0.0, None)

    xy = _sample_categorical(dep, rng)
    x = xy % dx + 1
    y = xy // dx + 1

    if permute:
        x, f = permute_labels(x, f, rng)
        y, g = permute_labels(y, g, rng)

    if np.any(np.isnan(f)) or np.any(np.isnan(g)):
        raise ValueError("Data not generated correctly (NaN in f or g).")

    return {"x": x, "y": y, "z": z, "f": f, "g": g}
