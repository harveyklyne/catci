"""form_t_sigma must reproduce the R form_T_Sigma oracle for both normalise settings."""

import numpy as np
from conftest import mat, vec

from catci.gcm import form_t_sigma


def test_form_t_sigma_matches_oracle(oracle):
    inp = oracle["form_T_Sigma"]["inputs"]
    x = np.asarray(inp["x"], dtype=int)
    y = np.asarray(inp["y"], dtype=int)
    f = mat(inp["f"])
    g = mat(inp["g"])

    for case in oracle["form_T_Sigma"]["cases"]:
        ts = form_t_sigma(x, y, f, g, normalise=case["normalise"])
        np.testing.assert_allclose(ts.T_vector, vec(case["T_vector"]), atol=1e-9, rtol=1e-7)
        np.testing.assert_allclose(ts.Sigma, mat(case["Sigma"]), atol=1e-9, rtol=1e-7)
