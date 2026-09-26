"""Exact CDF of a weighted sum of chi-squares, ``Q = sum_j lambda_j Z_j^2``.

This is the quantity Box (1954) *approximates* in
:func:`catci.statistic.approx_chi_statistic`. Under the Gaussian limit
``T ~ N(0, Sigma)``, ``||T||^2`` has exactly this law with ``lambda = eig(Sigma)``,
so ``P(Q <= ||t||^2)`` is the exact version of that statistic -- same ``||T||^2``,
no matrix inversion, and therefore none of the instability the paper attributes
to the pseudo-inverse statistic. Used by :class:`catci.statistic.ExactChi`.

Two standard algorithms, both taking the eigenvalues as input:

``imhof``
    Numerical inversion of the characteristic function, Imhof (1961). What
    ``CompQuadForm::imhof`` computes. Accurate, but the integrand oscillates
    with frequency ``q/2`` and the truncation point blows up as the number of
    non-zero eigenvalues falls, so it is slowest exactly where the search ends up.

``ruben``
    Ruben (1962) / Farebrother AS 204: a central-chi-square mixture whose
    coefficients satisfy a recursion. No oscillation, geometric convergence at
    rate set by the spread of the eigenvalues.

Both are pure functions of ``(q, lambdas)``.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import gammainc, gammaln
from scipy.stats import chi2

__all__ = ["imhof_cdf", "ruben_cdf", "RubenSeries", "exact_cdf", "positive_eigenvalues"]

_EIG_REL_TOL = 1e-10


def positive_eigenvalues(Sigma: np.ndarray, rel_tol: float = _EIG_REL_TOL) -> np.ndarray:
    """Eigenvalues of symmetric ``Sigma`` above ``rel_tol * max``; the rest are numerical zeros."""
    Sigma = np.asarray(Sigma, dtype=float)
    vals = np.linalg.eigvalsh((Sigma + Sigma.T) / 2.0)
    if vals.size == 0:
        return vals
    top = vals[-1]
    if top <= 0:
        return vals[:0]
    return vals[vals > rel_tol * top]


# --------------------------------------------------------------------------- #
# Imhof (1961)
# --------------------------------------------------------------------------- #
def imhof_cdf(
    q: float,
    lambdas: np.ndarray,
    u0: float = 1.0,
    limit: int = 100,
    epsabs: float = 1e-7,
    epsrel: float = 1e-5,
) -> float:
    """``P(sum_j lambda_j Z_j^2 <= q)`` by Imhof's characteristic-function inversion.

    ``P(Q > q) = 1/2 + (1/pi) int_0^inf sin(theta(u)) / (u rho(u)) du`` with
    ``theta(u) = A(u) - q u / 2``, ``A(u) = 0.5 sum_j arctan(lambda_j u)`` and
    ``rho(u) = prod_j (1 + lambda_j^2 u^2)^{1/4}``.

    Integrating that directly is what ``CompQuadForm::imhof`` does, and it fails
    whenever few eigenvalues are non-zero: the integrand decays only as
    ``u^{-1-r/2}``, so the truncation point needed for a given accuracy explodes as
    ``r`` falls, while the integrand oscillates at frequency ``q/2`` throughout.
    At ``r = 3`` the adaptive quadrature simply gives up. That is precisely the
    regime the merge search ends in, so the oscillation is peeled out instead:

        sin(A - wu) = sin(A) cos(wu) - cos(A) sin(wu),   w = q/2

    Both ``sin(A(u))/(u rho(u))`` and ``cos(A(u))/(u rho(u))`` are smooth and
    decaying, so ``[u0, inf)`` becomes a pair of Fourier integrals that QUADPACK's
    QAWF evaluates without resolving every oscillation. Only ``[0, u0]``, which
    carries at most a few periods, goes through ordinary adaptive quadrature.
    The split is needed because ``cos(A)/(u rho)`` alone diverges at the origin.
    """
    lam = np.asarray(lambdas, dtype=float)
    lam = lam[lam > 0]
    if lam.size == 0:
        return 1.0
    if q <= 0:
        return 0.0
    if lam.size == 1:
        return float(chi2.cdf(q / lam[0], df=1))

    w = 0.5 * q
    lam2 = lam * lam
    lam_sum = float(lam.sum())

    # Every integrand evaluation is a Python callback over `lam`, and the whole
    # cost of this function is the number of those callbacks -- hence the hoists.
    def head(u: float) -> float:
        if u == 0.0:
            return 0.5 * (lam_sum - q)
        A = 0.5 * np.arctan(lam * u).sum()
        return np.sin(A - w * u) * np.exp(-0.25 * np.log1p(lam2 * u * u).sum()) / u

    def amp_cos(u: float) -> float:
        A = 0.5 * np.arctan(lam * u).sum()
        return np.sin(A) * np.exp(-0.25 * np.log1p(lam2 * u * u).sum()) / u

    def amp_sin(u: float) -> float:
        A = 0.5 * np.arctan(lam * u).sum()
        return np.cos(A) * np.exp(-0.25 * np.log1p(lam2 * u * u).sum()) / u

    kw = dict(limit=limit, epsabs=epsabs, epsrel=epsrel)
    val = quad(head, 0.0, u0, **kw)[0]
    val += quad(amp_cos, u0, np.inf, weight="cos", wvar=w, **kw)[0]
    val -= quad(amp_sin, u0, np.inf, weight="sin", wvar=w, **kw)[0]

    sf = 0.5 + val / np.pi
    return float(min(max(1.0 - sf, 0.0), 1.0))


# --------------------------------------------------------------------------- #
# Ruben (1962) / Farebrother AS 204
# --------------------------------------------------------------------------- #
class RubenSeries:
    """Ruben's central-chi-square mixture for a fixed spectrum.

    ``F(q) = sum_k a_k P(chi^2_{r + 2k} <= q / beta)`` with

    * ``beta = min_j lambda_j``, which makes every ``a_k >= 0`` and ``sum_k a_k = 1``,
      so the truncated series has the rigorous remainder bound :attr:`truncation_error`;
    * ``a_0 = prod_j sqrt(beta / lambda_j)``, ``g_i = 0.5 sum_j (1 - beta/lambda_j)^i``,
      and ``k a_k = sum_{i=1}^{k} g_i a_{k-i}``.

    The coefficients depend on the spectrum alone, never on ``q``. That split is the
    point of this class: in the merge search the spectrum is a function of the
    partition, so a whole series can be built once per partition and then evaluated
    at each bootstrap draw's ``q`` for the cost of one dot product.
    """

    __slots__ = ("r", "beta", "a", "half_r", "truncation_error")

    def __init__(self, lambdas: np.ndarray, eps: float = 1e-12, max_terms: int = 100_000):
        lam = np.asarray(lambdas, dtype=float)
        lam = lam[lam > 0]
        self.r = lam.size
        if self.r == 0:
            self.beta = 1.0
            self.a = np.ones(1)
            self.half_r = 0.5
            self.truncation_error = 0.0
            return

        self.beta = float(lam.min())
        self.half_r = self.r / 2.0
        ratio = 1.0 - self.beta / lam  # in [0, c], c = 1 - lambda_min / lambda_max

        c = float(ratio.max())
        if c <= 0.0:  # all eigenvalues equal: the series is one exact term
            self.a = np.ones(1)
            self.truncation_error = 0.0
            return

        # c^k is only a rough guide to how fast the mass accumulates, so start from
        # it and extend until the *measured* remaining mass is below eps.
        n_terms = min(max_terms, int(np.ceil(np.log(eps) / np.log(c))) + 20)
        log_a0 = 0.5 * float(np.sum(np.log(self.beta / lam)))

        while True:
            powers = ratio[:, None] ** np.arange(1, n_terms + 1)[None, :]
            g = 0.5 * powers.sum(axis=0)  # g_i = 0.5 sum_j ratio_j^i

            a = np.empty(n_terms + 1)
            a[0] = np.exp(log_a0)
            for k in range(1, n_terms + 1):
                a[k] = float(g[:k] @ a[k - 1::-1]) / k

            # Coefficients are non-negative and sum to 1, so the mass we dropped is
            # itself the worst-case error of the truncated CDF.
            err = float(max(0.0, 1.0 - a.sum()))
            if err <= eps or n_terms >= max_terms:
                break
            n_terms = min(max_terms, 2 * n_terms)

        self.a = a
        self.truncation_error = err

    def cdf(self, q: float) -> float:
        """``P(sum_j lambda_j Z_j^2 <= q)``."""
        if self.r == 0:
            return 1.0
        if q <= 0:
            return 0.0

        z = 0.5 * q / self.beta
        K = self.a.size - 1
        # G_0 = P(chi^2_r <= q/beta); G_{k+1} = G_k - t_k with
        # t_k = z^{r/2+k} e^{-z} / Gamma(r/2+k+1), so t_{k+1} = t_k * z / (r/2+k+1).
        G0 = float(gammainc(self.half_r, z))
        if K == 0:
            return float(min(max(G0, 0.0), 1.0))

        log_z = np.log(z)
        log_t0 = self.half_r * log_z - z - gammaln(self.half_r + 1.0)
        # log t_k = log t_0 + sum_{m=1}^{k} (log z - log(r/2 + m)), for k = 0..K-1.
        incr = log_z - np.log(self.half_r + np.arange(1, K))
        log_t = log_t0 + np.concatenate(([0.0], np.cumsum(incr)))
        t = np.where(log_t < -745.0, 0.0, np.exp(np.minimum(log_t, 709.0)))

        G = np.empty(K + 1)
        G[0] = G0
        G[1:] = G0 - np.cumsum(t)
        np.clip(G, 0.0, 1.0, out=G)

        return float(min(max(float(self.a @ G), 0.0), 1.0))


def ruben_cdf(q: float, lambdas: np.ndarray, eps: float = 1e-12) -> float:
    """``P(sum_j lambda_j Z_j^2 <= q)`` by Ruben's series; builds the series each call."""
    return RubenSeries(lambdas, eps=eps).cdf(q)


def ruben_terms(lambdas: np.ndarray, eps: float = 1e-12) -> float:
    """How many series terms Ruben would need -- the cost of the cheap method.

    The coefficient mass accumulates at rate ``c = 1 - lambda_min / lambda_max``, so
    the term count grows with the condition number. Merging label pairs drives the
    condition number up by orders of magnitude (``A^T A`` is singular), which is why
    this has to be checked rather than assumed.
    """
    lam = np.asarray(lambdas, dtype=float)
    lam = lam[lam > 0]
    if lam.size <= 1:
        return 1.0
    c = 1.0 - lam.min() / lam.max()
    if c <= 0.0:
        return 1.0
    return float(np.ceil(np.log(eps) / np.log(c)))


#: Above this many Ruben terms, fall back to Imhof.
RUBEN_TERM_BUDGET = 400


def exact_cdf(q: float, lambdas: np.ndarray, method: str = "auto") -> float:
    """``P(sum_j lambda_j Z_j^2 <= q)``, exactly.

    ``method="auto"`` picks Ruben when its series is short and Imhof otherwise. The
    two methods' hard cases are complementary: Ruben slows down as the spectrum
    spreads, Imhof as the number of non-zero eigenvalues falls.
    """
    if method == "auto":
        method = "ruben" if ruben_terms(lambdas) <= RUBEN_TERM_BUDGET else "imhof"
    if method == "ruben":
        return ruben_cdf(q, lambdas)
    if method == "imhof":
        return imhof_cdf(q, lambdas)
    raise ValueError("method must be 'auto', 'ruben' or 'imhof'.")
