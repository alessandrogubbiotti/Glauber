"""
ff_extract_check.py — CORRECTNESS GATE for the free-fermion 3-time engine.

The Glauber chain is a free-fermion (quadratic Majorana) system, so the two-time
Majorana contraction

    Phi_{ab}(s) = <psi0| a_a e^{s Lt} a_b |psi0>          (s = t_b - t_a >= 0)

must factorise through a single-particle (2M x 2M) propagator:

    a_b(s) = e^{s Lt} a_b e^{-s Lt} = sum_g R(s)_{bg} a_g,   R(s) = e^{s Gsp},
    Phi(s)_{ab} = <psi0| a_a a_b(s) |psi0> = (C @ R(s).T)_{ab},

with  C_{ab} = <psi0| a_a a_b |psi0>  (static Majorana covariance) and
Gsp_{bg} from  [Lt, a_b] = sum_g Gsp_{bg} a_g.

This script extracts C and Gsp from the BRUTE-FORCE operators (pfaffian_verify),
and checks:
  (1) [Lt, a_b] really is linear in the Majoranas (=> Lt is quadratic);
  (2) Phi(s) = C @ e^{s Gsp}.T reproduces the brute-force contraction for several s;
  (3) i^k Pf[assembled Phi] reproduces the exact bond correlator at k=2,3.
  (4) report the band structure of C and Gsp (for the scalable builder).

Run:  python ff_extract_check.py [M]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.linalg import expm
import pfaffian_verify as pv
from multitime_interface_dpp import pfaffian


def extract(M, eta):
    m = pv.FermiModel(M, eta)
    n = 2 * M
    dim = 1 << M
    a = m.a
    Lt = m.Lt
    psi0 = m.psi0.astype(complex)
    apsi = [a[i] @ psi0 for i in range(n)]

    # static covariance C_{ab} = <psi0| a_a a_b |psi0> = (a_a psi0)^H (a_b psi0)
    C = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            C[i, j] = np.vdot(apsi[i], apsi[j])

    # single-particle generator Gsp_{bg} from [Lt, a_b] = sum_g Gsp_{bg} a_g
    # Gsp_{bg} = Tr[[Lt,a_b] a_g] / Tr[a_g a_g],  Tr[a_g a_d] = dim * delta
    Gsp = np.zeros((n, n), dtype=complex)
    lin_resid = 0.0
    for b in range(n):
        comm = Lt @ a[b] - a[b] @ Lt
        for g in range(n):
            Gsp[b, g] = np.trace(comm @ a[g]) / dim
        # residual: how much of comm is NOT a single Majorana (quadratic test)
        recon = sum(Gsp[b, g] * a[g] for g in range(n))
        lin_resid = max(lin_resid, np.abs(comm - recon).max())
    return m, C, Gsp, lin_resid


def phi_brute(m, a_idx, b_idx, s):
    """<psi0| a_{a_idx} e^{s Lt} a_{b_idx} |psi0>, s>=0."""
    psi0 = m.psi0.astype(complex)
    v = m.a[b_idx] @ psi0
    if s > 1e-14:
        v = expm(s * m.Lt) @ v
    v = m.a[a_idx] @ v
    return complex(psi0 @ v)              # psi0 real => <psi0|.> = psi0 . v


def assemble_phi(C, Gsp, legs):
    """legs = [(leg_index, time), ...] time-sorted ascending. Build the 2k x 2k
    antisymmetric contraction matrix Phi using the single-particle factorisation."""
    nlegs = len(legs)
    Phi = np.zeros((nlegs, nlegs), dtype=complex)
    # precompute R(s) for the needed time gaps
    for p in range(nlegs):
        for q in range(p + 1, nlegs):
            lp, tp = legs[p]
            lq, tq = legs[q]
            s = tq - tp                    # >= 0 since sorted
            R = expm(s * Gsp) if s > 1e-14 else np.eye(Gsp.shape[0])
            Phi[p, q] = (C @ R.T)[lp, lq]
            Phi[q, p] = -Phi[p, q]
    return Phi


def bond_corr_ff(C, Gsp, points):
    """i^k Pf[Phi] for the bond correlator, points=[(x,t),...]."""
    pts = sorted(points, key=lambda p: p[1])
    legs = []
    for (x, t) in pts:
        legs.append((2 * x + 1, t))
        legs.append((2 * x + 2, t))
    Phi = assemble_phi(C, Gsp, legs)
    k = len(pts)
    return (1j ** k) * pfaffian(Phi)


def main(M=6):
    eta = 1.0 / 3.0
    m, C, Gsp, lin_resid = extract(M, eta)
    print(f"OPEN chain M={M}, eta={eta:.4f}")
    print(f"(1) [Lt,a_b] linear-in-Majorana residual (quadratic test): {lin_resid:.2e}")
    print(f"    C antisym-part check |C+C.T - 2I|? C_ab+C_ba should be 2 delta: "
          f"{np.abs(C + C.T - 2*np.eye(2*M)).max():.2e}")

    # (2) Phi(s) factorisation vs brute force
    worst = 0.0
    for s in (0.0, 0.2, 0.5, 1.0):
        R = expm(s * Gsp)
        Phi_fac = C @ R.T
        for (aa, bb) in [(1, 4), (2, 5), (1, 2), (3, 6), (0, 7)]:
            pb = phi_brute(m, aa, bb, s)
            worst = max(worst, abs(Phi_fac[aa, bb] - pb))
    print(f"(2) max |Phi_factorised - Phi_brute|  over s,legs: {worst:.2e}")

    # (3) i^k Pf vs exact bond correlator
    tests = {2: [(0, 0.0), (2, 0.6)],
             3: [(0, 0.0), (2, 0.45), (3, 0.95)]}
    print("(3) bond correlator: free-fermion i^k Pf  vs  exact operator product")
    for k, pts in tests.items():
        ff = bond_corr_ff(C, Gsp, pts)
        ex = m.bond_corr(pts)
        print(f"    k={k}: ff={ff.real:+.10f}  exact={ex.real:+.10f}  "
              f"|diff|={abs(ff-ex):.2e}")

    # (4) band structure
    np.set_printoptions(precision=3, suppress=True, linewidth=200)
    print("\n(4) Gsp (single-particle generator), real part:")
    print(Gsp.real)
    print("    Gsp imag part max:", np.abs(Gsp.imag).max())
    print("\n    C (static covariance) imag part (the Gamma covariance):")
    print(C.imag)
    print("    C real part max off-diag:", np.abs(C.real - np.eye(2*M)).max())


if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    main(M)
