"""
glauber_exact_ring.py — exact correlators of the 1D Glauber-Ising model
on a periodic ring of l sites, in the rate convention of the C engine:

  spin-flip rates (alpha = 2):
    both neighbours disagree -> 1 + gamma   (annihilation)
    one  neighbour disagrees -> 1           (movement)
    no   neighbour disagrees -> 1 - gamma   (creation)

  micro time t = N^2 * T_macro  (engine rescaling, interface_state.h)

Conventions match Henkel, arXiv:2501.05912 up to a factor 2 in time
(Henkel normalises the movement flip rate to 1/2): t_Henkel = 2 * t_engine.

Implemented (all exact on the discrete finite ring):
  C_eq_ring        equilibrium correlator (eta^n + eta^(l-n)) / (1 + eta^l)
  two_time_stationary   G_n(tau) in the stationary state (FFT spectral)
  single_time_quench    C_n(t) after a quench, DST solver with C_0 = 1 pinned
  two_time_from_profile C_n(tau, s) = ring kernel applied to C(s) (FFT)
  theta_single_time / theta_two_time   Henkel's continuum finite-N formulas
                                       (3.7) and (3.10), as cross-checks
"""

import numpy as np
from numpy.fft import fft, ifft
from scipy.fft import dst, idst


# ─────────────────────────────────────────────────────────────────────────────
# Equilibrium / invariant correlator on the ring
# ─────────────────────────────────────────────────────────────────────────────

def gamma_from_eta(eta: float) -> float:
    """gamma = tanh(2 beta J) given eta = tanh(beta J)."""
    return 2.0 * eta / (1.0 + eta * eta)


def C_eq_ring(l: int, eta: float) -> np.ndarray:
    """Ising-ring pair correlator <sigma_0 sigma_n>, n = 0..l-1.
    Equals the conditioned-even Bernoulli(p) bond measure with p=(1-eta)/2."""
    n = np.arange(l, dtype=float)
    return (eta**n + eta**(l - n)) / (1.0 + eta**l)


# ─────────────────────────────────────────────────────────────────────────────
# Stationary two-time correlator  G_n(tau) = <sigma_0(0) sigma_n(tau)>
# ─────────────────────────────────────────────────────────────────────────────

def two_time_stationary(l: int, eta: float, t_micro: float) -> np.ndarray:
    """Exact stationary G_n(tau) on the ring, engine micro time t_micro.

    Conditional magnetisation obeys ds_n/dt = -2 s_n + gamma (s_{n-1}+s_{n+1}),
    so G(tau) = K_tau * C_eq with symbol exp(-2 t (1 - gamma cos q)),
    q_k = 2 pi k / l.  Returns full array n = 0..l-1.
    """
    gamma = gamma_from_eta(eta)
    q = 2.0 * np.pi * np.arange(l) / l
    decay = np.exp(-2.0 * t_micro * (1.0 - gamma * np.cos(q)))
    return np.real(ifft(fft(C_eq_ring(l, eta)) * decay))


# ─────────────────────────────────────────────────────────────────────────────
# Quench: single-time correlator with the C_0 = 1 constraint (discrete, exact)
# ─────────────────────────────────────────────────────────────────────────────

def single_time_quench(l: int, gamma: float, t_micro_arr,
                       C0: np.ndarray = None) -> np.ndarray:
    """Exact C_n(t) on the ring after a quench, engine time units.

    EOM (n = 1..l-1):  dC_n/dt = -4 C_n + 2 gamma (C_{n-1} + C_{n+1}),
    with C_0 = C_l = 1 pinned (Dirichlet).  This is Henkel eq. (2.3) times
    the engine factor 2.  Solved by DST-I eigenmodes sin(pi k n / l) with
    eigenvalues lambda_k = -4 (1 - gamma cos(pi k / l)).

    For gamma = 1 the equilibrium profile is C_eq = 1 identically; for
    gamma < 1 we relax towards the finite-l equilibrium of the SAME pinned
    problem, which for the quench observable is handled by evolving the
    deviation from the particular solution.  Here we only need gamma = 1
    (T = 0 quench), where the stationary solution of the pinned problem is
    exactly 1, so u = C - 1 obeys the homogeneous equation.

    C0: initial correlator C_n(0), length l (default: disordered delta).
    Returns array (len(t_arr), l//2 + 1).
    """
    if C0 is None:
        C0 = np.zeros(l)
        C0[0] = 1.0

    if abs(gamma - 1.0) > 1e-12:
        raise NotImplementedError("quench solver implemented for gamma=1")

    u0 = C0[1:l].copy() - 1.0          # deviation from pinned value 1
    coef = dst(u0, type=1)             # DST-I, modes k = 1..l-1
    k = np.arange(1, l, dtype=float)
    lam = -4.0 * (1.0 - np.cos(np.pi * k / l))

    t_arr = np.atleast_1d(np.asarray(t_micro_arr, dtype=float))
    out = np.empty((len(t_arr), l // 2 + 1))
    for i, t in enumerate(t_arr):
        u_t = idst(coef * np.exp(lam * t), type=1)
        full = np.empty(l)
        full[0] = 1.0
        full[1:] = 1.0 + u_t
        out[i] = full[: l // 2 + 1]
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Quench: two-time correlator from the single-time profile (discrete, exact)
# ─────────────────────────────────────────────────────────────────────────────

def two_time_from_profile(C_s_full: np.ndarray, gamma: float,
                          tau_micro: float) -> np.ndarray:
    """C_n(tau, s) = <sigma_n(s+tau) sigma_0(s)> on the ring.

    In tau the conditional magnetisation evolves freely (no pinning):
    dC_n/dtau = -2 C_n + gamma (C_{n-1} + C_{n+1})   (engine units;
    this is Henkel eq. (2.4) times the engine factor 2).
    C_s_full: single-time correlator at waiting time s, full length-l array.
    """
    l = len(C_s_full)
    q = 2.0 * np.pi * np.arange(l) / l
    decay = np.exp(-2.0 * tau_micro * (1.0 - gamma * np.cos(q)))
    return np.real(ifft(fft(C_s_full) * decay))


def symmetrize(C_half: np.ndarray, l: int) -> np.ndarray:
    """Extend C_n given on n = 0..l/2 to the full ring by C_n = C_{l-n}."""
    full = np.empty(l)
    half = l // 2
    full[: half + 1] = C_half
    full[half + 1:] = C_half[1:half][::-1]
    return full


# ─────────────────────────────────────────────────────────────────────────────
# Henkel's continuum finite-N formulas (theta functions) — cross-checks
# ─────────────────────────────────────────────────────────────────────────────

def _theta3(z, q, kmax=60):
    """Jacobi theta_3(z, q) = 1 + 2 sum_{k>=1} q^{k^2} cos(2 k z)."""
    z = np.asarray(z, dtype=float)
    out = np.ones_like(z)
    for k in range(1, kmax + 1):
        qk = q ** (k * k)
        if qk < 1e-300:
            break
        out = out + 2.0 * qk * np.cos(2.0 * k * z)
    return out


def theta_single_time(x_arr, t_henkel: float, Nring: int, n_quad=2000):
    """Henkel (3.7): C(t; x) = int_0^1 du theta3(pi/2 (u + |x|/N), e^{-pi^2 t/N^2}).
    x in lattice units, t in Henkel time units (= 2 * engine micro time)."""
    qq = np.exp(-np.pi**2 * t_henkel / Nring**2)
    u = (np.arange(n_quad) + 0.5) / n_quad
    x_arr = np.atleast_1d(np.asarray(x_arr, dtype=float))
    out = np.empty(len(x_arr))
    for i, x in enumerate(x_arr):
        out[i] = np.mean(_theta3(np.pi / 2 * (u + abs(x) / Nring), qq))
    return out


def theta_two_time(x_arr, tau_henkel: float, s_henkel: float,
                   Nring: int, n_quad=400):
    """Henkel (3.10): two-time correlator on the finite ring (continuum).
    Note the tau/2 in the second theta: the two-time kernel diffuses with
    D = 1/2 in Henkel units."""
    qs = np.exp(-np.pi**2 * s_henkel / Nring**2)
    qt = np.exp(-np.pi**2 * (tau_henkel / 2.0) / Nring**2)
    u = (np.arange(n_quad) + 0.5) / n_quad
    v = (np.arange(n_quad) + 0.5) / n_quad
    U, V = np.meshgrid(u, v, indexing="ij")
    TS = _theta3(np.pi / 2 * (U + V), qs)        # waiting-time factor
    x_arr = np.atleast_1d(np.asarray(x_arr, dtype=float))
    out = np.empty(len(x_arr))
    for i, x in enumerate(x_arr):
        xi = abs(x) / Nring
        T1 = _theta3(np.pi / 2 * (xi - V), qt)
        T2 = _theta3(np.pi / 2 * (xi + V), qt)
        out[i] = 0.5 * np.mean(TS * (T1 + T2))
    return out
