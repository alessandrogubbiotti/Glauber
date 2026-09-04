"""
paper_kernel_check.py — verify that the EXACT convention written in Theorem
thm:pfaffian of glauber_paper.tex (eqs e:Kpf, e:pfhier) reproduces the validated
continuum interface density.

The theorem states: order the points so tau_1 >= ... >= tau_k (DESCENDING); form
the 2k x 2k antisymmetric matrix whose upper 2x2 block (i<j) is
        K(x_i - x_j, tau_i - tau_j),
lower blocks the negative transposes, diagonal blocks zero, with the complex kernel

        K(z, dt) = [[ P,              i(2H - P - H') ],
                    [ i(P - 2H - H'),  P             ]]   (z, dt),  H' = d_z H,

and the connected k-point function is  (-i/2)^k Pf[...].  This module checks that
construction against fk_code (the cont3_final.py convention, already matched to the
discrete N^3 free-fermion limit) for k=2 and k=3, and confirms it is independent of
the z-sign orientation.

Companion to cont3_final.py / ff_engine.py.  See the %NOTE at e:pfhier in
Interfacce/glauber_paper.tex.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))                       # theory/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))      # python_scripts/
import numpy as np
import theory_continuum as tc
from multitime_interface_dpp import pfaffian


def K(z, dt):
    """The complex kernel e:Kpf (P on the diagonal)."""
    H = tc.H(z, dt)[0]; H1 = tc.dxH(z, dt)[0]; P = tc.P(z, dt)[0]
    return np.array([[P, 1j * (2*H - P - H1)],
                     [1j * (P - 2*H - H1), P]], dtype=complex)


def fk_code(points):
    """Validated cont3_final.py convention (ascending sort), the reference."""
    pts = sorted(points, key=lambda p: p[1])           # tau ascending
    k = len(pts); A = np.zeros((2*k, 2*k), dtype=complex)
    for p in range(k):
        for q in range(p + 1, k):
            xp, tp = pts[p]; xq, tq = pts[q]
            B = K(xp - xq, tq - tp)                     # p the earlier time
            A[2*p:2*p+2, 2*q:2*q+2] = B
            A[2*q:2*q+2, 2*p:2*p+2] = -B.T
    return ((-0.5j) ** k) * pfaffian(A)


def fk_paper(points, zsign=+1):
    """The literal Theorem thm:pfaffian convention (e:pfhier): DESCENDING time
    order, upper block (i<j) = K(x_i - x_j, tau_i - tau_j)."""
    pts = sorted(points, key=lambda p: p[1], reverse=True)   # tau DESCENDING
    k = len(pts); A = np.zeros((2*k, 2*k), dtype=complex)
    for i in range(k):
        for j in range(i + 1, k):
            xi, ti = pts[i]; xj, tj = pts[j]           # ti >= tj
            B = K(zsign * (xi - xj), ti - tj)
            A[2*i:2*i+2, 2*j:2*j+2] = B
            A[2*j:2*j+2, 2*i:2*i+2] = -B.T
    return ((-0.5j) ** k) * pfaffian(A)


def _close(a, b, tol=1e-9):
    return abs(a - b) < tol


if __name__ == "__main__":
    ok = True

    print("== k=2 ==  paper convention vs fk_code vs tc.f")
    for (x, tau) in [(0.3, 0.1), (0.5, 0.2), (0.2, 0.05)]:
        pts = [(0.0, 0.0), (x, tau)]
        c = fk_code(pts); a = fk_paper(pts, +1); b = fk_paper(pts, -1)
        ref = tc.f(x, tau)[0]
        ok &= _close(a.real, c.real) and _close(b.real, c.real) and _close(a.real, ref)
        ok &= abs(a.imag) < 1e-12 and abs(b.imag) < 1e-12
        print(f"  x={x},tau={tau}: paper(+)={a.real:+.6e}(Im{a.imag:+.0e})  "
              f"paper(-)={b.real:+.6e}  code={c.real:+.6e}  f={ref:+.6e}")

    print("\n== k=3 ==  paper convention vs fk_code (matched to discrete N^3)")
    configs = [
        dict(name="A", xs=[-0.2, 0.0, 0.25], taus=[0.05, 0.10, 0.18]),
        dict(name="B", xs=[0.0, 0.15, 0.35], taus=[0.04, 0.09, 0.15]),
        dict(name="C", xs=[0.1, -0.15, 0.3], taus=[0.07, 0.12, 0.05]),
        dict(name="D", xs=[-0.1, 0.05, 0.22], taus=[0.03, 0.11, 0.20]),
    ]
    for cfg in configs:
        pts = list(zip(cfg['xs'], cfg['taus']))
        c = fk_code(pts); a = fk_paper(pts, +1); b = fk_paper(pts, -1)
        ok &= _close(a.real, c.real) and _close(b.real, c.real)
        ok &= abs(a.imag) < 1e-12 and abs(b.imag) < 1e-12
        print(f"  {cfg['name']}: paper(+)={a.real:+.6e}(Im{a.imag:+.0e})  "
              f"paper(-)={b.real:+.6e}  code={c.real:+.6e}")

    print("\nALL CHECKS PASS" if ok else "\n*** MISMATCH ***")
    sys.exit(0 if ok else 1)
