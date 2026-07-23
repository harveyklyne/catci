"""Parametric bootstrap sampling ``T ~ N(0, Sigma)``.

Ports ``matrix_sqrt`` / ``bootstrap_T`` (``R/bootstrap_functions.R``). ``Sigma``
is typically low-rank (it is a covariance of merged residual products), so we
use an eigendecomposition-based PSD square root rather than a pivoted Cholesky;
the sampling law ``N(0, Sigma)`` is identical.
"""

from __future__ import annotations

import numpy as np

__all__ = ["matrix_sqrt", "bootstrap_T"]


def matrix_sqrt(Sigma: np.ndarray) -> np.ndarray:
    """A matrix ``A`` with ``A A^T = Sigma`` for symmetric PSD ``Sigma``.

    Returns an ``(p, r)`` factor where ``r`` is the numerical rank; columns
    beyond the rank (zero/negative eigenvalues) are dropped.
    """
    Sigma = np.asarray(Sigma, dtype=float)
    Sigma = (Sigma + Sigma.T) / 2.0
    vals, vecs = np.linalg.eigh(Sigma)
    tol = max(Sigma.shape) * np.finfo(float).eps * max(vals.max(), 0.0)
    keep = vals > tol
    return vecs[:, keep] * np.sqrt(vals[keep])


def bootstrap_T(Sigma: np.ndarray, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Draw ``n_boot`` samples ``T ~ N(0, Sigma)`` as columns of a ``(p, n_boot)`` array."""
    A = matrix_sqrt(Sigma)  # (p, r)
    r = A.shape[1]
    return A @ rng.standard_normal(size=(r, n_boot))
