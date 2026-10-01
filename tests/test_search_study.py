"""``experiments/search_study.py``: every DGP's mean must lie in the range of ``Sigma``.

A component outside it would be invisible to the bootstrap draws (which live in
the range) and make the "power" partly an artefact.
"""

import numpy as np
import pytest

from search_study import direction, uniform_C


@pytest.mark.parametrize("name", ["greedy_trap", "binary_tree", "step"])
def test_directions_in_range_of_sigma(name):
    d = 8
    C = uniform_C(d)
    Sigma = np.kron(C, C)
    u = direction(name, d)
    assert np.allclose(Sigma @ np.linalg.pinv(Sigma) @ u, u)
