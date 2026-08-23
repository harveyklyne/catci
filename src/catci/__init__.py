"""catci: conditional independence testing for structured categorical data.

See README.md for the design rationale and ``tests/fixtures/`` for the frozen
R oracle this implementation is differential-tested against.
"""

from __future__ import annotations

from .api import CatciResult, catci_test
from .criteria import ApproxChi, approx_chi_metric, euclid, max_abs, mgcm
from .gcm import TSigma, form_t_sigma
from .search import greedy_search
from .structure import Ordinal, Saturated, Structure, Tree, make_binary_tree

__all__ = [
    "catci_test",
    "CatciResult",
    "form_t_sigma",
    "TSigma",
    "greedy_search",
    "Structure",
    "Ordinal",
    "Saturated",
    "Tree",
    "make_binary_tree",
    "ApproxChi",
    "approx_chi_metric",
    "euclid",
    "max_abs",
    "mgcm",
]
