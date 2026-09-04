"""
theory_glauber.py — exact finite-N correlators of the 1D Glauber-Ising chain.

Conventions are those of docs/section1.tex (Section 1 of the paper):

    gamma = (N^2 - nu^2) / (N^2 + nu^2)      beta = log(N/nu) / 2
    eta   = (N - nu)     / (N + nu)          p    = nu / (N + nu)

All correlators below are functions of *Glauber time* t (hop rate 1/2 per
direction).  The Monte-Carlo engine works in macroscopic time tau = micro/N^2
with hop rate 1 per direction, so its generator is twice the Glauber one and
the formulas must be evaluated at

    t = t_formula(tau, N) = 2 * N^2 * tau .

Numerics.  Modified Bessel functions are evaluated through the exponentially
scaled  ive(k, x) = exp(-|x|) I_k(x).  For x = g t >= 0 this gives, e.g.,

    exp(-2t) I_a(g t) I_b(g t) = exp(-2(1-g) t) ive(a, g t) ive(b, g t),

which never forms the astronomically large raw I_a(g t) at large argument.
Every series in the lattice index is truncated at |index| <= K, with K doubled
(default start 400) until the result stops changing by more than `tol`.
"""

import numpy as np
from scipy.special import ive

__all__ = [
    "gamma_of", "eta_of", "p_of", "beta_of", "t_formula",
    "C_eq", "G_two_time", "C_single_disordered", "C_single_thesis",
    "H_interface", "C_eq_ring", "G_two_time_ring",
]


# ---------------------------------------------------------------------------
# Derived parameters
# ---------------------------------------------------------------------------

def gamma_of(N, nu=1.0):
    """gamma = (N^2 - nu^2)/(N^2 + nu^2) = tanh(2 beta)."""
    N = float(N); nu = float(nu)
    return (N * N - nu * nu) / (N * N + nu * nu)


def eta_of(N, nu=1.0):
    """eta = (N - nu)/(N + nu) = tanh(beta)."""
    N = float(N); nu = float(nu)
    return (N - nu) / (N + nu)


def p_of(N, nu=1.0):
    """Equilibrium interface density p = nu/(N + nu) = (1 - eta)/2."""
    N = float(N); nu = float(nu)
    return nu / (N + nu)


def beta_of(N, nu=1.0):
    """beta = log(N/nu)/2."""
    return 0.5 * np.log(float(N) / float(nu))


def t_formula(tau_macro, N):
    """Glauber time corresponding to engine macroscopic time tau: 2 N^2 tau."""
    return 2.0 * float(N) ** 2 * np.asarray(tau_macro, dtype=float)


def _gamma_from_eta(eta):
    """gamma = 2 eta / (1 + eta^2)  (inverse of eta = (1-sqrt(1-g^2))/g)."""
    return 2.0 * eta / (1.0 + eta * eta)


def _resolve(eta, gamma):
    eta = float(eta)
    g = _gamma_from_eta(eta) if gamma is None else float(gamma)
    return eta, g


def _converge(compute, K0=400, tol=1e-12, Kmax=1 << 22):
    """Evaluate ``compute(K)`` with K doubled from K0 until the result settles."""
    K = int(K0)
    prev = np.asarray(compute(K), dtype=float)
    while K < Kmax:
        K *= 2
        cur = np.asarray(compute(K), dtype=float)
        if np.max(np.abs(cur - prev)) < tol:
            return cur
        prev = cur
    return prev


def _as_array(n):
    a = np.atleast_1d(np.asarray(n))
    return a, np.isscalar(n) or (np.ndim(n) == 0)


# ---------------------------------------------------------------------------
# Equilibrium pair correlator (Glauber 63 ; thesis 39 ; Henkel 2.10)
# ---------------------------------------------------------------------------

def C_eq(n, eta):
    """<sigma_0 sigma_n>_eq = eta^|n| on the infinite line."""
    a, scalar = _as_array(n)
    out = float(eta) ** np.abs(a).astype(float)
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Two-time correlator at equilibrium (thesis 43 ; solution of Henkel 2.4)
#   G_n(t) = exp(-t) sum_k eta^|k| I_{n-k}(g t)
# ---------------------------------------------------------------------------

def G_two_time(n, t, eta, gamma=None, K0=400, tol=1e-12):
    eta, g = _resolve(eta, gamma)
    a, scalar = _as_array(n)
    a = a.astype(int)
    gt = g * float(t)
    pref = np.exp(-(1.0 - g) * float(t))

    def compute(K):
        k = np.arange(-K, K + 1)
        w = eta ** np.abs(k).astype(float)              # (2K+1,)
        order = a[:, None] - k[None, :]                 # (Nn, 2K+1)
        return pref * (ive(order, gt) * w[None, :]).sum(axis=1)

    out = _converge(compute, K0, tol)
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Single-time correlator after a disordered (T = inf) quench
# Henkel (2.9):  C_n(t) = eta^n - sum_{m>=1} eta^m exp(-2t)[I_{n-m}(2gt)-I_{n+m}(2gt)]
# ---------------------------------------------------------------------------

def C_single_disordered(n, t, eta, gamma=None, K0=400, tol=1e-12):
    eta, g = _resolve(eta, gamma)
    a, scalar = _as_array(n)
    a = np.abs(a.astype(int))
    g2t = 2.0 * g * float(t)
    pref = np.exp(-2.0 * (1.0 - g) * float(t))
    eq = eta ** a.astype(float)

    def compute(K):
        m = np.arange(1, K + 1)
        w = eta ** m.astype(float)                       # (K,)
        term = ive(a[:, None] - m[None, :], g2t) - ive(a[:, None] + m[None, :], g2t)
        return eq - pref * (w[None, :] * term).sum(axis=1)

    out = _converge(compute, K0, tol)
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Single-time correlator, thesis (41) double-Bessel form, disordered data.
# Equals C_single_disordered by the Bessel addition theorem (A&S 9.6.33):
#   r_{0,n}(t) = eta^n - exp(-2t) sum_{j>=1} eta^j sum_{l in Z}
#                [ I_{-l}(gt) I_{n-l-j}(gt) - I_{-l-j}(gt) I_{n-l}(gt) ]
# ---------------------------------------------------------------------------

def C_single_thesis(n, t, eta, gamma=None, K0=400, tol=1e-12):
    eta, g = _resolve(eta, gamma)
    a, scalar = _as_array(n)
    a = a.astype(int)
    gt = g * float(t)
    pref = np.exp(-2.0 * (1.0 - g) * float(t))

    def compute(K):
        l = np.arange(-K, K + 1)
        j = np.arange(1, K + 1)
        wj = eta ** j.astype(float)                      # (J,)
        I_ml = ive(-l, gt)                               # (L,)
        res = np.empty(len(a), dtype=float)
        for idx, nn in enumerate(a):
            # T1 = sum_j eta^j sum_l ive(-l) ive(n-l-j)
            M1 = ive(nn - l[None, :] - j[:, None], gt)   # (J, L)
            inner1 = (I_ml[None, :] * M1).sum(axis=1)    # (J,)
            T1 = (wj * inner1).sum()
            # T2 = sum_j eta^j sum_l ive(-l-j) ive(n-l)
            I_nl = ive(nn - l, gt)                       # (L,)
            M2 = ive(-l[None, :] - j[:, None], gt)       # (J, L)
            inner2 = (M2 * I_nl[None, :]).sum(axis=1)    # (J,)
            T2 = (wj * inner2).sum()
            res[idx] = eta ** float(nn) - pref * (T1 - T2)
        return res

    out = _converge(compute, K0, tol)
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Time-delayed interface correlator at equilibrium, thesis (45) regrouped:
#   H_j(t) = ((1-eta)/2)^2
#          + (1/4)(1/eta - eta) exp(-2t)
#            sum_{l<=0} sum_{m>=1} eta^{m-l}
#               [ I_{j-l} I_{j+1-m} - I_{j-m} I_{j+1-l} ](g t)
# The double sum factorizes:  (w_l . A)(w_m . B) - (w_m . C)(w_l . D).
# ---------------------------------------------------------------------------

def H_interface(j, t, eta, gamma=None, K0=400, tol=1e-12):
    eta, g = _resolve(eta, gamma)
    a, scalar = _as_array(j)
    a = a.astype(int)
    gt = g * float(t)
    pref = np.exp(-2.0 * (1.0 - g) * float(t))
    plateau = ((1.0 - eta) / 2.0) ** 2
    coef = 0.25 * (1.0 / eta - eta)

    def compute(K):
        li = np.arange(0, K + 1)            # |l|, l = 0,-1,...,-K
        m = np.arange(1, K + 1)
        wl = eta ** li.astype(float)        # w_l = eta^{-l} = eta^{|l|}
        wm = eta ** m.astype(float)         # w_m = eta^{m}
        res = np.empty(len(a), dtype=float)
        for idx, jj in enumerate(a):
            A = ive(jj + li, gt)            # I_{j-l},   l=-li
            D = ive(jj + 1 + li, gt)        # I_{j+1-l}
            B = ive(jj + 1 - m, gt)         # I_{j+1-m}
            C = ive(jj - m, gt)             # I_{j-m}
            dbl = (wl @ A) * (wm @ B) - (wm @ C) * (wl @ D)
            res[idx] = plateau + coef * pref * dbl
        return res

    out = _converge(compute, K0, tol)
    return float(out[0]) if scalar else out


# ---------------------------------------------------------------------------
# Exact finite-ring objects (no infinite-line truncation)
# ---------------------------------------------------------------------------

def C_eq_ring(l, eta):
    """Ring equilibrium correlator C_eq(m) = (eta^m + eta^{l-m})/(1 + eta^l)."""
    m = np.arange(l, dtype=float)
    eta = float(eta)
    return (eta ** m + eta ** (l - m)) / (1.0 + eta ** l)


def G_two_time_ring(n, t, l, eta, gamma=None):
    """Exact finite-ring two-time correlator (thesis/section1 eq. 23):

        G_n(t) = (1/l) sum_q exp(-(1 - g cos(2 pi q/l)) t) cos(2 pi q n/l) Chat(q),

    with Chat the DFT of the ring equilibrium correlator C_eq_ring.  Returns the
    value(s) at the requested separation(s) n.
    """
    eta, g = _resolve(eta, gamma)
    Ceq = C_eq_ring(l, eta)
    q = 2.0 * np.pi * np.arange(l) / l
    decay = np.exp(-(1.0 - g * np.cos(q)) * float(t))
    G_full = np.real(np.fft.ifft(np.fft.fft(Ceq) * decay))
    a, scalar = _as_array(n)
    out = G_full[np.mod(a.astype(int), l)]
    return float(out[0]) if scalar else out
