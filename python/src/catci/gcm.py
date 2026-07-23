"""Generalised covariance measure: build ``(T, Sigma)`` from labels + propensities.

Pure port of ``form_T_Sigma`` (``R/simulation_settings.R``); no model fitting
happens here (that is the ``learners`` layer). Pinned by the ``form_T_Sigma``
fixture for both ``normalise`` settings.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["TSigma", "form_t_sigma"]


@dataclass(frozen=True)
class TSigma:
    """The generalised-covariance statistic and its covariance estimate."""

    T_vector: np.ndarray  # length dx*dy
    Sigma: np.ndarray  # (dx*dy, dx*dy)

    @property
    def dx_dy(self) -> int:
        return self.T_vector.shape[0]


def _one_hot(labels: np.ndarray, d: int) -> np.ndarray:
    """n x d indicator matrix for 1-based integer ``labels`` in ``{1,...,d}``."""
    labels = np.asarray(labels)
    out = np.zeros((labels.shape[0], d), dtype=float)
    out[np.arange(labels.shape[0]), labels - 1] = 1.0
    return out


def form_t_sigma(
    x: np.ndarray,
    y: np.ndarray,
    f: np.ndarray,
    g: np.ndarray,
    normalise: bool,
) -> TSigma:
    """Form ``(T_vector, Sigma)`` from labels ``x, y`` and propensity matrices ``f, g``.

    Parameters
    ----------
    x, y : integer arrays of length ``n``, values in ``{1,...,dx}`` / ``{1,...,dy}``.
    f, g : arrays ``(n, dx)`` / ``(n, dy)``; row ``i`` estimates ``P(X=.|Z_i)`` / ``P(Y=.|Z_i)``.
    normalise : if True, use the normalised generalised covariance (unit ``diag(Sigma)``).

    The per-row product is ``kron(Y_row, X_row)`` -> length ``dx*dy`` with X fastest,
    matching the ``(1,1),(2,1),...`` ordering used throughout.
    """
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f, dtype=float)
    g = np.asarray(g, dtype=float)

    n = x.shape[0]
    dx = f.shape[1]
    dy = g.shape[1]
    if y.shape[0] != n:
        raise ValueError("x and y must have the same length.")
    if f.shape[0] != n or g.shape[0] != n:
        raise ValueError("f and g must have the same number of rows as x and y.")
    if x.min() < 1 or y.min() < 1 or x.max() > dx or y.max() > dy:
        raise ValueError("Need 1 <= x <= dx and 1 <= y <= dy.")

    X = _one_hot(x, dx) - f  # residuals
    Y = _one_hot(y, dy) - g

    # prod[i] = kron(Y[i], X[i]) -> X varies fastest, matching R's `Y[row,] %x% X[row,]`.
    prod_mat = (Y[:, :, None] * X[:, None, :]).reshape(n, dx * dy)

    T_vector = np.sqrt(n) * prod_mat.mean(axis=0)
    Sigma = np.cov(prod_mat, rowvar=False, ddof=1)  # matches R stats::var (N-1 denominator)

    if normalise:
        scale = np.sqrt(np.diag(Sigma))
        T_vector = T_vector / scale
        Sigma = Sigma / np.outer(scale, scale)

    return TSigma(T_vector=T_vector, Sigma=Sigma)
