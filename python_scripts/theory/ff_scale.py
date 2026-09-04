"""
ff_scale.py — diffusive scaling limit of the connected k-point INTERFACE
correlator via the exact free-fermion engine on a long OPEN chain.

Calibration (engine units, glauber_exact_ring):
    eta_N = (N-1)/(N+1),   x = site/N,   t_micro = N^2 * tau.
Pointwise limits (no box-averaging needed; the engine is exact):
    N^2 <eta_0(0) eta_{Nx}(N^2 tau)>_c        -> f(x,tau)              (k=2)
    N^3 <eta(0) eta(N^2 t2) eta(N^2 t3)>_c     -> f_3 continuum density (k=3)
The continuum reference is the LINE Pfaffian of the VALIDATED complex kernel,
cont3_final.fk(., Kcont), anchored by f at k=2 (ratio 1.0 verified).

NOTE (fixed 2026-08-12): this script previously took its k=3 continuum reference
from scaling_check_3time.cont_interface_cumulant, which is built from the real
block [[H,-H'],[H',-H'']] and is WRONG for k>=3 -- at the configs below it is
wrong-signed and off by a config-dependent factor (ratio to the validated
reference: -1.63 at config A, -0.92 at config B).  The discrete side was never
affected; only the comparison target was.  See VERIFICATION.md.

Run:  python ff_scale.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))   # python_scripts/
import numpy as np
import ff_engine as fe
import theory_continuum as tc
from cont3_final import Kcont, fk


def richardson(invN, vals):
    """Linear fit vals ~ a + b/N, return intercept a (1/N -> 0)."""
    A = np.vstack([np.ones_like(invN), invN]).T
    return np.linalg.lstsq(A, vals, rcond=None)[0][0]


def k2_anchor(L=6):
    print(f"== k=2 anchor: N^2 <eta_0 eta_Nx>_c  vs  f(x,tau)   (open chain, L={L}) ==")
    x, tau = 0.3, 0.10
    print(f"   x={x}, tau={tau};   f(x,tau) = {tc.f(x,tau)[0]:+.6e}")
    print(f"   {'N':>4} {'eta':>8} {'M':>5} {'N^2*conn':>13} {'ratio->f':>10}")
    rows = []
    for N in (10, 20, 40, 80):
        eta = (N - 1) / (N + 1)
        M = N * L
        ff = fe.FF(M, eta)
        c0 = M // 2
        n1 = c0 + int(round(N * x))
        disc = N**2 * ff.iface_cumulant([(c0, 0.0), (n1, N*N*tau)]).real
        rows.append((1.0/N, disc))
        print(f"   {N:>4} {eta:>8.5f} {M:>5} {disc:>13.6e} {disc/tc.f(x,tau)[0]:>10.4f}")
    inv = np.array([r[0] for r in rows]); dd = np.array([r[1] for r in rows])
    print(f"   Richardson 1/N->0: {richardson(inv,dd):+.6e}   "
          f"(f = {tc.f(x,tau)[0]:+.6e}, ratio {richardson(inv,dd)/tc.f(x,tau)[0]:.4f})\n")


def k3_limit(L=6):
    print(f"== k=3: N^3 <eta eta eta>_c  vs  continuum two-triangles   (open chain, L={L}) ==")
    configs = [
        dict(name="A", xs=[-0.2, 0.0, 0.25], taus=[0.05, 0.10, 0.18]),
        dict(name="B", xs=[0.0, 0.15, 0.35], taus=[0.04, 0.09, 0.15]),
    ]
    for cfg in configs:
        xs, taus = cfg['xs'], cfg['taus']
        cont = fk(list(zip(xs, taus)), Kcont).real
        print(f"  config {cfg['name']}  xs={xs}  taus={taus}")
        print(f"   continuum f_3 (LINE) = {cont:+.6e}")
        print(f"   {'N':>4} {'M':>5} {'N^3*conn':>14} {'ratio->cont':>12}")
        rows = []
        for N in (10, 20, 40, 80):
            eta = (N - 1) / (N + 1)
            M = N * L
            ff = fe.FF(M, eta)
            c0 = M // 2
            pts = [(c0 + int(round(N*x)), N*N*t) for x, t in zip(xs, taus)]
            disc = N**3 * ff.iface_cumulant(pts).real
            rows.append((1.0/N, disc))
            print(f"   {N:>4} {M:>5} {disc:>14.6e} {disc/cont:>12.4f}")
        inv = np.array([r[0] for r in rows]); dd = np.array([r[1] for r in rows])
        ext = richardson(inv, dd)
        print(f"   Richardson 1/N->0: {ext:+.6e}   ratio {ext/cont:.4f}\n")


if __name__ == "__main__":
    k2_anchor()
    k3_limit()
