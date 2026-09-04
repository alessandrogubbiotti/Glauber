"""
cont3_final.py — the continuum connected k-point INTERFACE density as the exact
scaling limit of the validated discrete free-fermion Pfaffian.

Discrete CENTERED bond moment = i^k Pf[zero-diagonal-block 2k x 2k Majorana
contraction] (paper Cor. c:centered); it coincides with the connected cumulant
for k <= 3 only (k4_zerodiag_check.py), so "connected" below means k <= 3.
Scaling each off-diagonal 2x2 block by N gives the continuum block (extracted and
identified):

    K(z,dt) = [[ P,            i(2H - P - H') ],
               [ i(P - 2H - H'), P            ]]   (z, dt),   H' = d_z H,

with H,H',P the line functions of theory_continuum.  Then

    f_k(points) = (-i/2)^k * Pf[ off-diagonal, time-ordered K-blocks ].

For k=2 this is exactly f = HP - H^2 + 1/4 H'^2 (W=2H-P forced by the anchor).
This module tests f_3 against the discrete N^3 limit, and, as a NEGATIVE
control, the real block [[H,-H'],[H',-H'']] -- the kernel the paper printed in
the theorem before 2026-06-27, which now survives only as the k=2 covariance
"shadow" e:Kblock (Cor. c:k2).  It is kept under its historical name `Kpaper`
(alias `Kblock`).  Note that fk(real kernel).real == 0 at odd k is structural
(imaginary prefactor times a real Pfaffian, Remark r:Kstruct), so the
meaningful negative control is the historical REAL-prefactor evaluation in
scaling_check_3time, which is nonzero and wrong-signed (gate [5] of
run_checks.py).
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import numpy as np
import theory_continuum as tc
from multitime_interface_dpp import pfaffian


def Kcont(z, dt):
    H = tc.H(z, dt)[0]; H1 = tc.dxH(z, dt)[0]; P = tc.P(z, dt)[0]
    return np.array([[P, 1j * (2*H - P - H1)],
                     [1j * (P - 2*H - H1), P]], dtype=complex)


def Kpaper(z, dt):
    """Real block e:Kblock = [[H, -H'], [H', -H'']] (H'' = 4(H-P)).  Historical
    name: it was the block printed in thm:pfaffian before 2026-06-27; it is the
    k=2 covariance only.  Negative control, never a reference for k >= 3."""
    H = tc.H(z, dt)[0]; H1 = tc.dxH(z, dt)[0]; P = tc.P(z, dt)[0]
    return np.array([[H, -H1], [H1, -4.0 * (H - P)]], dtype=complex)


Kblock = Kpaper                           # the paper's current name for it


def fk(points, kernel):
    """(-i/2)^k Pf[off-diagonal time-ordered k-blocks]: the continuum CENTERED
    k-point density (e:pfhier); equals the connected one for k <= 3."""
    pts = sorted(points, key=lambda p: p[1])         # tau ascending
    k = len(pts)
    A = np.zeros((2 * k, 2 * k), dtype=complex)
    for p in range(k):
        for q in range(p + 1, k):
            xp, tp = pts[p]; xq, tq = pts[q]
            B = kernel(xp - xq, tq - tp)
            A[2*p:2*p+2, 2*q:2*q+2] = B
            A[2*q:2*q+2, 2*p:2*p+2] = -B.T
    return ((-0.5j) ** k) * pfaffian(A)


def box_fk(centers, taus, h, kernel=Kcont, nq=7):
    """Box-averaged continuum density: average f_k over a width-h box in EACH
    spatial argument, matching the box observable eta_box = (1/h) int_box eta.

    This is the validated replacement for scaling_check_3time.cont_box_cumulant,
    which built the same average from the REAL paper block and is wrong (indeed
    wrong-signed) for k >= 3.  Any k (centered density, as fk); midpoint rule
    with nq nodes per axis.  LINE geometry only: the torus image-sum option
    (L=) of the old function is deliberately not carried over -- a torus-
    geometry caller must not migrate here.  Gated in run_checks.py [6].
    """
    from itertools import product
    centers = list(centers); taus = list(taus)
    if len(centers) != len(taus):
        raise ValueError("centers and taus must have the same length")
    s = (np.arange(nq) + 0.5) / nq * h - h / 2.0      # uniform box samples
    acc = 0.0
    for shifts in product(s, repeat=len(centers)):
        pts = [(c + d, t) for c, d, t in zip(centers, shifts, taus)]
        acc += fk(pts, kernel).real
    return acc / (nq ** len(centers))


if __name__ == "__main__":
    import ff_engine as fe

    print("k=2 anchor: fk(Kcont) and fk(Kpaper) vs f")
    for (x, tau) in [(0.3, 0.1), (0.5, 0.2), (0.2, 0.05)]:
        a = fk([(0.0, 0.0), (x, tau)], Kcont)
        b = fk([(0.0, 0.0), (x, tau)], Kpaper)
        print(f"  x={x},tau={tau}: Kcont={a.real:+.6e}(Im{a.imag:+.0e})  "
              f"Kpaper={b.real:+.6e}  f={tc.f(x,tau)[0]:+.6e}")

    print("\nk=3: discrete N^3 limit (Richardson)  vs  fk(Kcont)  vs  fk(Kpaper)")
    configs = [
        dict(name="A", xs=[-0.2, 0.0, 0.25], taus=[0.05, 0.10, 0.18]),
        dict(name="B", xs=[0.0, 0.15, 0.35], taus=[0.04, 0.09, 0.15]),
        dict(name="C", xs=[0.1, -0.15, 0.3], taus=[0.07, 0.12, 0.05]),
        dict(name="D", xs=[-0.1, 0.05, 0.22], taus=[0.03, 0.11, 0.20]),
    ]
    for cfg in configs:
        xs, taus = cfg['xs'], cfg['taus']
        rows = []
        for N in (20, 40, 80):
            eta = (N - 1) / (N + 1); M = N * 8
            ff = fe.FF(M, eta); c0 = M // 2
            pts = [(c0 + int(round(N*x)), N*N*t) for x, t in zip(xs, taus)]
            rows.append((1.0/N, N**3 * ff.iface_cumulant(pts).real))
        inv = np.array([r[0] for r in rows]); dd = np.array([r[1] for r in rows])
        disc = np.linalg.lstsq(np.vstack([np.ones(3), inv]).T, dd, rcond=None)[0][0]
        cc = fk(list(zip(xs, taus)), Kcont)
        cp = fk(list(zip(xs, taus)), Kpaper)
        rp = f"{disc/cp.real:.4f}" if cp.real != 0 else "n/a (exactly 0)"
        print(f"  {cfg['name']}: discrete={disc:+.6e}  "
              f"fk(Kcont)={cc.real:+.6e} (r={disc/cc.real:.4f})  "
              f"fk(Kblock)={cp.real:+.6e} (r={rp})")
