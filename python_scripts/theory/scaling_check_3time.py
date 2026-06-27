"""
scaling_check_3time.py
======================

!! SUPERSEDED by cont3_final.py (see ../../VERIFICATION.md #1).  The continuum
!! kernel used below, K = [[H,-H'],[H',-H'']], is the one PRINTED in the paper
!! (e:Kblock) and is WRONG for the connected k>=3 correlator: it gives an
!! identically-zero connected 3-point.  The validated complex kernel
!! K = [[P, i(2H-P-H')],[i(P-2H-H'), P]] lives in cont3_final.py.  Keep this file
!! only as the discrete-side exact-generator reference; do not trust its
!! continuum comparison for k>=3.

Does the THREE-time interface correlator obey the scaling limit predicted by
Theorem 1.11(ii) of glauber_paper.tex?  Two things to check:

  (a) SCALING RATE.  The theorem says the connected k-point interface correlator
      rescales with one power of N per point:  N^k <eta...eta>_c -> finite.  For
      k=3 we verify N^3 gives a finite nonzero limit while N^2 -> 0 and N^4 -> inf.

  (b) LIMIT VALUE.  That limit is the continuum "two oriented triangles": the
      connected (single-loop) part of the Pfaffian of the 2x2 block kernel
          K(z,t) = [[ H, -H' ],[ H', -H'' ]](z,t),   H'' = 4(H - P),
      assembled block-antisymmetrically.  We anchor the continuum formula by
      checking the k=2 case reproduces f = HP - H^2 + 1/4 (H')^2, then compare
      the N^3-rescaled exact correlator (calibrated ring, eta_N=(N-1)/(N+1)) to
      the continuum two-triangles, and extrapolate in 1/N.

Discrete side: exact generator (free of any closed-form assumption), sped up by
the reversible eigendecomposition L~ = D^{1/2} L D^{-1/2} (symmetric).

Run:  python scaling_check_3time.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))                       # theory/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))      # python_scripts/

import numpy as np
from math import factorial
from itertools import combinations
from scipy.linalg import eigh

from glauber_exact_ring import gamma_from_eta
from multitime_interface_dpp import build_generator, ising_measure, pfaffian
import theory_continuum as tc


# ─────────────────────────────────────────────────────────────────────────────
# Fast exact model (reversible eigendecomposition)
# ─────────────────────────────────────────────────────────────────────────────
class FastModel:
    def __init__(self, l, eta):
        self.l = l
        self.eta = eta
        self.gamma = gamma_from_eta(eta)
        L, spins = build_generator(l, self.gamma)
        mu = ising_measure(spins, eta)
        sq = np.sqrt(mu)
        Lsym = (sq[:, None] * L) / sq[None, :]      # D^{1/2} L D^{-1/2}, symmetric
        Lsym = 0.5 * (Lsym + Lsym.T)                # kill 1e-15 asymmetry
        self.w, self.V = eigh(Lsym)
        self.sq = sq
        self.mu = mu
        self.spins = spins
        self.wall = np.empty((l, 1 << l))
        for x in range(l):
            self.wall[x] = 0.5 * (1.0 - spins[:, x].astype(float) *
                                  spins[:, (x + 1) % l].astype(float))
        self.sig = spins.astype(float)

    def P_apply(self, dt, vec):
        """Return expm(dt L) @ vec via the symmetric eigendecomposition."""
        if dt <= 1e-14:
            return vec
        x = self.sq * vec
        x = self.V @ (np.exp(dt * self.w) * (self.V.T @ x))
        return x / self.sq

    def moment_general(self, obs_points):
        pts = sorted(obs_points, key=lambda p: p[1])
        g = pts[-1][0].copy()
        for k in range(len(pts) - 2, -1, -1):
            dt = pts[k + 1][1] - pts[k][1]
            g = pts[k][0] * self.P_apply(dt, g)
        return float(self.mu @ g)

    def boxvec(self, center, w):
        """Box wall-density observable: mean of eta over 2w+1 bonds about center."""
        idx = [(center + i) % self.l for i in range(-w, w + 1)]
        return self.wall[idx].mean(axis=0)

    def cumulant(self, obs_points):
        """Joint cumulant of diagonal observables (vec, t) at distinct times."""
        n = len(obs_points)
        total = 0.0
        for part in _set_partitions(list(range(n))):
            sz = len(part)
            coef = (-1) ** (sz - 1) * factorial(sz - 1)
            prod = 1.0
            for block in part:
                prod *= self.moment_general([obs_points[i] for i in block])
            total += coef * prod
        return total

    def box_cumulant(self, boxes):
        """Connected correlator of box observables.  boxes = [(center,w,t),...]."""
        return self.cumulant([(self.boxvec(c, w), t) for (c, w, t) in boxes])

    def spin_two_point(self, n, t):
        return self.moment_general([(self.sig[:, 0], 0.0), (self.sig[:, n], t)])


def _set_partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        yield [collection]; return
    first = collection[0]
    for smaller in _set_partitions(collection[1:]):
        for i, subset in enumerate(smaller):
            yield smaller[:i] + [[first] + subset] + smaller[i + 1:]
        yield [[first]] + smaller


# ─────────────────────────────────────────────────────────────────────────────
# Continuum two-triangles:  connected k-point of the interface field
# ─────────────────────────────────────────────────────────────────────────────
def _img(func, z, tau, L, nimg):
    z = np.atleast_1d(np.asarray(z, float))
    acc = np.zeros_like(z)
    for m in range(-nimg, nimg + 1):
        acc = acc + func(z + m * L, tau)
    return acc


def Kblock(z, tau, L=None, nimg=6):
    """2x2 kernel [[H, -H'],[H', -H'']](z,tau).  H''=4(H-P).
    If L is given, use the length-L torus image sums (with the tanh(L) parity
    factor on H, H', H'')."""
    if L is None:
        H = tc.H(z, tau)[0]; H1 = tc.dxH(z, tau)[0]; Pp = tc.P(z, tau)[0]
        H2 = 4.0 * (H - Pp)
    else:
        tnh = np.tanh(L)
        H = tnh * _img(lambda zz, tt: tc.H(zz, tt), z, tau, L, nimg)[0]
        H1 = tnh * _img(lambda zz, tt: tc.dxH(zz, tt), z, tau, L, nimg)[0]
        Hl = tnh * _img(lambda zz, tt: tc.H(zz, tt), z, tau, L, nimg)[0]
        Pl = tnh * _img(lambda zz, tt: tc.P(zz, tt), z, tau, L, nimg)[0]
        H2 = 4.0 * (Hl - Pl)
    return np.array([[H, -H1], [H1, -H2]])


def cont_interface_cumulant(points, L=None):
    """Connected k-point interface correlator (k=2 or 3) of the continuum field:
        kappa_k = -(-1/2)^k * Pf(A),
    A the 2k x 2k block-antisymmetric matrix with off-diagonal blocks
    K(x_i - x_j, |tau_i - tau_j|) (i<j) and zero diagonal blocks.
    Anchored so that k=2 gives f (validated in main)."""
    k = len(points)
    A = np.zeros((2 * k, 2 * k))
    for i in range(k):
        for j in range(i + 1, k):
            xi, ti = points[i]; xj, tj = points[j]
            B = Kblock(xi - xj, abs(ti - tj), L=L)
            A[2 * i:2 * i + 2, 2 * j:2 * j + 2] = B
            A[2 * j:2 * j + 2, 2 * i:2 * i + 2] = -B.T
    pf = float(np.real(pfaffian(A)))
    return -((-0.5) ** k) * pf


def cont_box_cumulant(centers, taus, h, L=None, nq=7):
    """Box-averaged continuum two-triangles: average f_k over a width-h box in
    EACH spatial argument (uniform per position), matching the box observable
    eta_box = (1/h) int_box eta.  centers/taus length k (=2 or 3)."""
    k = len(centers)
    s = (np.arange(nq) + 0.5) / nq * h - h / 2.0      # uniform box samples
    acc = 0.0; cnt = 0
    if k == 2:
        for s0 in s:
            for s1 in s:
                acc += cont_interface_cumulant(
                    [(centers[0] + s0, taus[0]), (centers[1] + s1, taus[1])], L=L)
                cnt += 1
    else:
        for s0 in s:
            for s1 in s:
                for s2 in s:
                    acc += cont_interface_cumulant(
                        [(centers[0] + s0, taus[0]),
                         (centers[1] + s1, taus[1]),
                         (centers[2] + s2, taus[2])], L=L)
                    cnt += 1
    return acc / cnt


# ─────────────────────────────────────────────────────────────────────────────
def main():
    Ns = [8, 10, 12]
    L_RING = 1                      # ring of l = N sites (unit macroscopic torus)
    H_BOX = 0.25                    # macroscopic box width (so eta is smeared)

    # macroscopic 3-box configurations (distinct times). centers as fractions of N.
    configs = [
        dict(name="A", fracs=[0.0, 0.25, 0.5], taus=[0.04, 0.09, 0.16]),
        dict(name="B", fracs=[0.0, 0.30, 0.55], taus=[0.06, 0.12, 0.20]),
    ]

    # ---- 0. convention sanity: continuum k=2 (pointwise) == f_torus ----
    print("== continuum k=2 anchor:  cont_interface_cumulant vs f_torus ==")
    for (x, tau) in [(0.3, 0.1), (0.5, 0.2), (0.2, 0.05)]:
        c2 = cont_interface_cumulant([(0.0, 0.0), (x, tau)], L=L_RING)
        fT = tc.f_torus(x, tau, L_RING)[0]
        print(f"   x={x}, tau={tau}:  cont k=2 = {c2:+.6e}   f_torus = {fT:+.6e}"
              f"   ratio = {c2/fT:.4f}")
    print()

    # ---- 1. spin two-point sanity: discrete -> H_torus ----
    print("== discrete spin two-point -> H_torus (calibration sanity) ==")
    for N in Ns:
        m = FastModel(N, (N - 1) / (N + 1))
        tau = 0.1; t = N * N * tau; nb = N // 4
        d = m.spin_two_point(nb, t)
        h = tc.H_torus(nb / N, tau, L_RING)[0]
        print(f"   N={N}: <s0 s_{nb}>(tau=0.1) = {d:+.5f}   H_torus = {h:+.5f}")
    print()

    def boxes_for(cfg, N):
        w = max(int(round(H_BOX * N / 2.0)), 1)          # box half-width in sites
        centers = [int(round(fr * N)) for fr in cfg['fracs']]
        return w, centers, (2 * w + 1) / N               # actual macroscopic width

    # ---- 2. scaling RATE on the BOX correlator: N^a, a=2,3,4 ----
    print(f"== scaling rate of the BOX 3-correlator (box width h~{H_BOX}):  "
          f"N^a * <eta_box eta_box eta_box>_c ==")
    print("   (a=3 should converge; a=2 -> 0, a=4 -> diverge)")
    for cfg in configs:
        print(f"  config {cfg['name']}  centers(frac)={cfg['fracs']}  tau={cfg['taus']}")
        print(f"   {'N':>3} {'h_N':>6} {'conn':>13} {'N^2':>11} {'N^3':>11} {'N^4':>11}")
        for N in Ns:
            m = FastModel(N, (N - 1) / (N + 1))
            w, centers, hN = boxes_for(cfg, N)
            bx = [(c, w, N * N * tau) for c, tau in zip(centers, cfg['taus'])]
            c = m.box_cumulant(bx)
            print(f"   {N:>3} {hN:>6.3f} {c:>13.3e} {N**2*c:>11.3e} "
                  f"{N**3*c:>11.3e} {N**4*c:>11.3e}")
        print()

    # ---- 3. LIMIT: N^3 * discrete BOX  vs  box-averaged continuum two-triangles ----
    print("== N^3 * box 3-correlator  vs  box-averaged continuum two-triangles ==")
    for cfg in configs:
        print(f"  config {cfg['name']}")
        print(f"   {'N':>3} {'h_N':>6} {'N^3*disc':>13} {'cont box':>13} {'ratio':>8}")
        vals = []
        for N in Ns:
            m = FastModel(N, (N - 1) / (N + 1))
            w, centers, hN = boxes_for(cfg, N)
            xs = [c / N for c in centers]
            bx = [(c, w, N * N * tau) for c, tau in zip(centers, cfg['taus'])]
            disc = N ** 3 * m.box_cumulant(bx)
            cont = cont_box_cumulant(xs, cfg['taus'], hN, L=L_RING)
            vals.append((1.0 / N, disc, cont))
            print(f"   {N:>3} {hN:>6.3f} {disc:>13.4e} {cont:>13.4e} {disc/cont:>8.3f}")
        inv = np.array([v[0] for v in vals]); dd = np.array([v[1] for v in vals])
        cc = np.array([v[2] for v in vals])
        A = np.vstack([np.ones_like(inv), inv]).T
        d0 = np.linalg.lstsq(A, dd, rcond=None)[0][0]
        c0 = np.linalg.lstsq(A, cc, rcond=None)[0][0]
        print(f"   Richardson 1/N->0:  discrete = {d0:+.4e}   continuum = {c0:+.4e}"
              f"   ratio = {d0/c0:.3f}")
        print()


if __name__ == "__main__":
    main()
