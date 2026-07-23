"""Public entry point: ``catci_test``.

This is what Algorithm 3 describes and what the README promised but the R
package never exposed -- a single function from data to a p-value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .calibrate import adaptive_pvalue
from .criteria import ApproxChi
from .gcm import form_t_sigma
from .learners import Learner, crossfit
from .search import greedy_search
from .structure import Structure

__all__ = ["CatciResult", "catci_test"]


@dataclass(frozen=True)
class CatciResult:
    p_value: float
    criteria: np.ndarray  # observed criterion path
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
    nfolds: int = 5,
    rng: Optional[np.random.Generator] = None,
) -> CatciResult:
    """Test conditional independence ``X _||_ Y | Z`` for categorical ``X, Y``.

    Provide propensities directly (``f``, ``g`` -- the oracle path) or a
    ``learner`` plus ``z`` to cross-fit them. Returns the p-value together with
    the observed criterion path and the partitions the search visited.
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
        f = crossfit(z, x, dx, learner, nfolds=nfolds, rng=rng)
        g = crossfit(z, y, dy, learner, nfolds=nfolds, rng=rng)

    ts = form_t_sigma(x, y, f, g, normalise=normalise)
    criterion = ApproxChi()

    observed = greedy_search(ts.T_vector, ts.Sigma, dx, dy, x_structure, y_structure, criterion)
    p = adaptive_pvalue(
        ts.T_vector, ts.Sigma, dx, dy, x_structure, y_structure,
        n_boot=n_boot, criterion=criterion, rng=rng,
    )
    return CatciResult(p_value=p, criteria=np.asarray(observed.values), partitions=observed.partitions)
