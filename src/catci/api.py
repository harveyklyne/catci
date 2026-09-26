"""Public entry point: ``catci_test``.

This is what Algorithm 3 describes and what the README promised but the R
package never exposed -- a single function from data to a p-value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .calibrate import adaptive_pvalue
from .statistic import ApproxChi
from .gcm import form_t_sigma
from .learners import Learner, fit_propensities
from .search import greedy_search
from .structure import Structure

__all__ = ["CatciResult", "catci_test"]


@dataclass(frozen=True)
class CatciResult:
    p_value: float
    statistics: np.ndarray  # observed statistic path
    partitions: list  # partition sequence visited by the observed search


def catci_test(
    x: np.ndarray,
    y: np.ndarray,
    x_structure: Structure,
    y_structure: Structure,
    z: Optional[np.ndarray] = None,
    f: Optional[np.ndarray] = None,
    g: Optional[np.ndarray] = None,
    learner: Optional[Learner] = None,
    n_boot: int = 100,
    normalise: bool = True,
    rng: Optional[np.random.Generator] = None,
    n_jobs: int = 1,
) -> CatciResult:
    """Test conditional independence ``X _||_ Y | Z`` for categorical ``X, Y``.

    Provide propensities directly (``f``, ``g`` -- the oracle path) or a
    ``learner`` plus ``z`` to fit them on the full sample. Returns the p-value
    together with the observed statistic path and the partitions the search
    visited. ``n_jobs`` threads share the bootstrap searches (``-1`` = all cores).
    """
    x = np.asarray(x)
    y = np.asarray(y)
    dx = int(x.max())
    dy = int(y.max())
    if rng is None:
        rng = np.random.default_rng()

    if f is None or g is None:
        if learner is None or z is None:
            raise ValueError("Provide either (f, g) or (learner, z).")
        f = fit_propensities(z, x, dx, learner)
        g = fit_propensities(z, y, dy, learner)

    ts = form_t_sigma(x, y, f, g, normalise=normalise)
    statistic = ApproxChi()

    observed = greedy_search(ts.T_vector, ts.Sigma, dx, dy, x_structure, y_structure, statistic)
    p = adaptive_pvalue(
        ts.T_vector, ts.Sigma, dx, dy, x_structure, y_structure,
        n_boot=n_boot, statistic=statistic, rng=rng, n_jobs=n_jobs,
    )
    return CatciResult(p_value=p, statistics=np.asarray(observed.values), partitions=observed.partitions)
