"""Contract tests for the propensity learners.

The learner interface is ``fit(z, labels, num_class) -> predict``, with
``predict(z_new)`` an ``(m, num_class)`` row-stochastic matrix. These pin that
contract for :func:`mlp_learner` -- in particular the two ways an sklearn
classifier breaks it: a column count set by the labels it happened to see rather
than by ``num_class``, and a column *order* keyed to ``classes_``.
"""

import numpy as np
import pytest

from catci.learners import fit_propensities, mlp_learner, xgboost_learner

FAST_MLP = {"hidden_layer_sizes": (4,), "max_iter": 60, "alpha": 3.0}


def _toy(rng, n=200, num_class=4, p=3):
    """Labels that genuinely depend on z, so a fitted learner has something to find."""
    z = rng.standard_normal((n, p))
    logits = z[:, [0]] * np.linspace(-2, 2, num_class)[None, :]
    probs = np.exp(logits)
    probs /= probs.sum(1, keepdims=True)
    labels = np.array([rng.choice(num_class, p=probs[i]) + 1 for i in range(n)])
    return z, labels, probs


def test_mlp_returns_row_stochastic_matrix_of_the_right_shape():
    rng = np.random.default_rng(0)
    z, labels, _ = _toy(rng)
    out = fit_propensities(z, labels, 4, mlp_learner(FAST_MLP))
    assert out.shape == (200, 4)
    assert np.all(out >= 0.0)
    np.testing.assert_allclose(out.sum(axis=1), 1.0, atol=1e-8)


def test_mlp_pads_absent_classes_to_num_class():
    """num_class, not the observed label set, fixes the width -- and the padding
    lands in the absent column rather than shifting the present ones left."""
    rng = np.random.default_rng(1)
    z = rng.standard_normal((120, 3))
    labels = rng.choice([1, 2, 4], size=120)  # class 3 never observed

    out = fit_propensities(z, labels, 4, mlp_learner(FAST_MLP))
    assert out.shape == (120, 4)
    assert np.all(out[:, 2] == 0.0)  # the absent class, not a shifted neighbour
    np.testing.assert_allclose(out.sum(axis=1), 1.0, atol=1e-8)


def test_mlp_is_deterministic_given_random_state():
    rng = np.random.default_rng(2)
    z, labels, _ = _toy(rng)
    a = fit_propensities(z, labels, 4, mlp_learner({**FAST_MLP, "random_state": 7}))
    b = fit_propensities(z, labels, 4, mlp_learner({**FAST_MLP, "random_state": 7}))
    np.testing.assert_allclose(a, b)


def test_mlp_beats_the_constant_predictor_when_z_is_informative():
    """A sanity floor: if the net were mis-wired (e.g. unscaled inputs collapsing
    it to the marginal) it would score no better than predicting 1/num_class."""
    rng = np.random.default_rng(3)
    z, labels, probs = _toy(rng, n=800)
    out = fit_propensities(z, labels, 4, mlp_learner({"hidden_layer_sizes": (8, 8),
                                                     "alpha": 3.0, "max_iter": 300}))
    fitted_err = np.mean((probs - out) ** 2)
    constant_err = np.mean((probs - 0.25) ** 2)
    assert fitted_err < constant_err


@pytest.mark.parametrize(
    "learner",
    [mlp_learner(FAST_MLP),
     xgboost_learner({"eta": 0.3, "max.depth": 2, "gamma": 0, "nrounds": 10})],
)
def test_learners_share_one_interface(learner):
    """Whatever the backend, the runner may swap one for the other blind."""
    rng = np.random.default_rng(4)
    z, labels, _ = _toy(rng)
    predict = learner(z, labels, 4)
    fresh = rng.standard_normal((17, 3))
    out = predict(fresh)
    assert out.shape == (17, 4)
    np.testing.assert_allclose(out.sum(axis=1), 1.0, atol=1e-6)
