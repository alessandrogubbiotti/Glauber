"""
exactdiag_check.py — gold-standard validation of the closed formulas
against the FULL generator of the spin dynamics on a small ring.

Builds the 2^l x 2^l generator L with the engine's rates
    both neighbours disagree -> 1 + gamma
    one  neighbour disagrees -> 1
    no   neighbour disagrees -> 1 - gamma
and verifies, at machine precision:

  1. REVERSIBILITY: detailed balance mu_x L_xy = mu_y L_yx w.r.t. the Ising
     ring measure mu (= Bernoulli(p) bonds conditioned even, p = (1-eta)/2),
     i.e. the Dirichlet form E(f,g) = -<f, Lg>_mu is symmetric.
  2. STATIONARY two-time correlator G_n(tau) vs the FFT spectral formula.
  3. QUENCH (gamma=1) single-time C_n(t) vs the DST solver with C_0=1 pinned.
  4. QUENCH two-time C_n(tau,s) vs ring-kernel propagation of C(s).

Run:  python3 python_scripts/theory/exactdiag_check.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from scipy.linalg import expm

from glauber_exact_ring import (C_eq_ring, gamma_from_eta, two_time_stationary,
                                single_time_quench, two_time_from_profile,
                                symmetrize)


def build_generator(l: int, gamma: float):
    """Full generator on {-1,1}^l, engine rate convention."""
    nstates = 1 << l
    spins = np.empty((nstates, l), dtype=np.int8)
    for x in range(nstates):
        for i in range(l):
            spins[x, i] = 1 if (x >> i) & 1 else -1

    L = np.zeros((nstates, nstates))
    for x in range(nstates):
        s = spins[x]
        for i in range(l):
            d = int(s[i] != s[(i - 1) % l]) + int(s[i] != s[(i + 1) % l])
            rate = (1.0 + gamma, 1.0, 1.0 - gamma)[2 - d]
            y = x ^ (1 << i)
            L[x, y] += rate          # column convention: dP/dt = P L? see below
    # We use row convention: (Lf)(x) = sum_y L[x,y] (f(y)-f(x))
    np.fill_diagonal(L, -L.sum(axis=1))
    return L, spins


def ising_measure(spins: np.ndarray, eta: float):
    """mu(sigma) propto exp(beta J sum sigma_i sigma_{i+1}); eta = tanh(beta J)."""
    K = np.arctanh(eta)
    energy = np.sum(spins * np.roll(spins, -1, axis=1), axis=1)
    w = np.exp(K * energy.astype(float))
    return w / w.sum()


def main(l=10, eta=0.6):
    gamma = gamma_from_eta(eta)
    print(f"l={l}, eta={eta}, gamma={gamma:.6f}")
    L, spins = build_generator(l, gamma)
    mu = ising_measure(spins, eta)

    # 1 ── detailed balance / Dirichlet-form symmetry
    DB = mu[:, None] * L - (mu[:, None] * L).T
    print(f"1. detailed balance  max |mu_x L_xy - mu_y L_yx| = "
          f"{np.abs(DB).max():.3e}")

    # stationarity of mu
    print(f"   stationarity      max |mu L|                 = "
          f"{np.abs(mu @ L).max():.3e}")

    # 2 ── stationary two-time correlator
    tau = 0.7
    P = expm(tau * L)
    sig = spins.astype(float)
    G_ed = np.array([np.sum(mu * sig[:, 0] * (P @ sig[:, n]))
                     for n in range(l)])
    G_th = two_time_stationary(l, eta, tau)
    print(f"2. stationary G_n(tau={tau})   max |exactdiag - FFT formula| = "
          f"{np.abs(G_ed - G_th).max():.3e}")
    # equilibrium correlator while we are at it
    C_ed = np.array([np.sum(mu * sig[:, 0] * sig[:, n]) for n in range(l)])
    print(f"   equilibrium C_eq          max diff = "
          f"{np.abs(C_ed - C_eq_ring(l, eta)).max():.3e}")

    # 3 ── quench at gamma=1 from uniform (disordered) measure
    Lq, _ = build_generator(l, 1.0)
    pi0 = np.full(1 << l, 1.0 / (1 << l))
    for t in (0.3, 1.0):
        Pt = expm(t * Lq)
        pit = pi0 @ Pt
        C_ed = np.array([np.sum(pit * sig[:, 0] * sig[:, n])
                         for n in range(l // 2 + 1)])
        C_th = single_time_quench(l, 1.0, [t])[0]
        print(f"3. quench C_n(t={t})        max |exactdiag - DST solver| = "
              f"{np.abs(C_ed - C_th).max():.3e}")

    # 4 ── quench two-time correlator
    s_, tau_ = 0.3, 0.5
    Ps   = expm(s_ * Lq)
    Ptau = expm(tau_ * Lq)
    pis  = pi0 @ Ps
    C2_ed = np.array([np.sum(pis * sig[:, 0] * (Ptau @ sig[:, n]))
                      for n in range(l)])
    C_s_full = symmetrize(single_time_quench(l, 1.0, [s_])[0], l)
    C2_th = two_time_from_profile(C_s_full, 1.0, tau_)
    print(f"4. quench C_n(tau={tau_},s={s_})  max |exactdiag - kernel formula| = "
          f"{np.abs(C2_ed - C2_th).max():.3e}")


if __name__ == "__main__":
    main()
