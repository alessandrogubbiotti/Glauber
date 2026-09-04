"""
ff_engine.py — scalable free-fermion engine for multi-time Glauber-Ising BOND
(interface) correlators on a long OPEN chain.

The Glauber chain is quadratic in Majoranas, so the two-time contraction
factorises (validated in ff_extract_check.py):

    Phi(s) = C @ R(s).T ,   C = I + i*Gamma ,   R(s) = expm(s*Gsp) ,
    <mu_{x1}(t1)...mu_{xk}(tk)> = i^k Pf[ Phi assembled over the 2k legs ].

C (static covariance) and Gsp (single-particle generator) are LOCAL/structured,
with closed forms read off and validated against the exact 2^M operators:

  Gsp = i*G,  G real antisymmetric, banded:
    G[2x,2x+1]   = -2/(1+eta^2)        (within-site; boundary ends: -2*sqrt(1-eta^2))
    G[2x+1,2x+2] = -2*gamma            (inter-site;  boundary ends: extracted)
    G[2x+1,2x+4] = +gamma*eta          (skip coupling)
  Gamma (= Im C) antisymmetric:
    Gamma[2x,2y+1] = -cL(2x)*cR(2y+1)*eta^(y-x)   for y>=x
                     cL(a)=1 if a==0 else sqrt(1-eta^2);  cR(b)=1 if b==2M-1 else sqrt(1-eta^2)
    Gamma[2x,2x-1] = +eta              (x>=1)

This reaches N in the hundreds (Pf of 6x6, expm of banded 2M x 2M).

Run:  python ff_engine.py        # validates builders vs brute force, then 3-pt
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.linalg import expm
from glauber_exact_ring import gamma_from_eta
from multitime_interface_dpp import pfaffian


def build_Gsp(M, eta, boundary=None):
    g = gamma_from_eta(eta)
    n = 2 * M
    G = np.zeros((n, n))
    a = 2.0 / (1.0 + eta * eta)          # within-site bulk
    b = 2.0 * g                          # inter-site bulk
    c = g * eta                          # skip coupling
    for x in range(M):
        if 2 * x + 1 < n:
            G[2 * x, 2 * x + 1] = -a
        if 2 * x + 2 < n:
            G[2 * x + 1, 2 * x + 2] = -b
        if 2 * x + 4 < n:
            G[2 * x + 1, 2 * x + 4] = +c
    # boundary within-site (both ends): -2 sqrt(1-eta^2)
    s = np.sqrt(max(1.0 - eta * eta, 0.0))
    G[0, 1] = -2.0 * s
    G[n - 2, n - 1] = -2.0 * s
    # boundary inter-site entries G[1,2] and G[n-3,n-2] (extracted, M-independent)
    if boundary is not None:
        G[1, 2] = boundary
        G[n - 3, n - 2] = boundary
    # Gsp = i * G, with G real antisymmetric (entries set in the upper triangle)
    return 1j * _to_antisym_upper(G)


def _to_antisym_upper(G):
    """Given a matrix whose intended nonzeros are all in the upper triangle (i<j),
    return the antisymmetric matrix A with A[i,j]=G[i,j], A[j,i]=-G[i,j]."""
    U = np.triu(G, 1)
    return U - U.T


def build_C(M, eta):
    n = 2 * M
    Gam = np.zeros((n, n))
    s = np.sqrt(max(1.0 - eta * eta, 0.0))
    for x in range(M):
        # backward coupling Gamma[2x, 2x-1] = +eta
        if x >= 1:
            Gam[2 * x, 2 * x - 1] = +eta
        # long tail Gamma[2x, 2y+1] = -cL*cR*eta^(y-x), y>=x
        cL = 1.0 if (2 * x == 0) else s
        for y in range(x, M):
            beta = 2 * y + 1
            cR = 1.0 if (beta == n - 1) else s
            Gam[2 * x, beta] = -cL * cR * eta ** (y - x)
    # antisymmetrise from the entries we set: some are upper (2x<2y+1) and the
    # backward ones have 2x>2x-1 (lower). Build full antisym carefully.
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if Gam[i, j] != 0.0:
                A[i, j] = Gam[i, j]
                A[j, i] = -Gam[i, j]
    return np.eye(n) + 1j * A, A


def extract_boundary_b(eta):
    """Get the exact boundary inter-site Gsp entry G[1,2] from a tiny brute force."""
    import ff_extract_check as ff
    _, _, Gsp, _ = ff.extract(5, eta)
    return Gsp.imag[1, 2]


# ---- correlators -----------------------------------------------------------
# Gsp = i*(real antisym) is HERMITIAN, eigenvalues come in +/- mu pairs.  The
# contraction Phi(s)=C e^{s Gsp^T} only involves the decaying (lambda<=0) sector;
# the growing (+mu) modes carry zero physical weight but make a naive expm(s*Gsp)
# overflow at s=N^2 tau.  We diagonalise once and drop the growing modes.
class FF:
    def __init__(self, M, eta, tol=1e-9):
        bnd = extract_boundary_b(eta)
        self.M = M
        self.C, _ = build_C(M, eta)
        self.Gsp = build_Gsp(M, eta, boundary=bnd)
        lam, V = np.linalg.eigh(self.Gsp)          # Hermitian -> real lam, unitary V
        self.lam = lam
        self.V = V
        self.tol = tol
        self._Rt = {}

    def _RT(self, s):
        """e^{s Gsp}^T with the growing (lambda>tol) modes dropped."""
        key = round(float(s), 12)
        if key not in self._Rt:
            e = np.where(self.lam > self.tol, 0.0, np.exp(np.minimum(s * self.lam, 0.0)))
            # R(s) = V diag(e) V^H ;  R(s)^T = conj(V) diag(e) V^T
            Vc = self.V.conj()
            self._Rt[key] = (Vc * e) @ self.V.T
        return self._Rt[key]

    def phi(self, legs):
        nlegs = len(legs)
        Phi = np.zeros((nlegs, nlegs), dtype=complex)
        for p in range(nlegs):
            for q in range(p + 1, nlegs):
                lp, tp = legs[p]; lq, tq = legs[q]
                s = tq - tp
                RT = self._RT(s) if s > 1e-14 else np.eye(2 * self.M)
                Phi[p, q] = (self.C @ RT)[lp, lq]
                Phi[q, p] = -Phi[p, q]
        return Phi

    def bond_corr(self, points):
        pts = sorted(points, key=lambda p: p[1])
        legs = []
        for (x, t) in pts:
            legs.append((2 * x + 1, t)); legs.append((2 * x + 2, t))
        k = len(pts)
        return (1j ** k) * pfaffian(self.phi(legs))

    def iface_cumulant(self, points):
        """Connected interface correlator <eta...eta>_c, eta=(1-mu)/2,
        kappa_eta = (-1/2)^k kappa_mu, via moment-cumulant on bond moments."""
        from math import factorial
        n = len(points)

        def parts(coll):
            coll = list(coll)
            if len(coll) == 1:
                yield [coll]; return
            first = coll[0]
            for sm in parts(coll[1:]):
                for i, ss in enumerate(sm):
                    yield sm[:i] + [[first] + ss] + sm[i + 1:]
                yield [[first]] + sm
        kap_mu = 0.0 + 0j
        for part in parts(range(n)):
            sz = len(part); coef = (-1) ** (sz - 1) * factorial(sz - 1)
            pr = 1.0 + 0j
            for blk in part:
                pr *= self.bond_corr([points[i] for i in blk])
            kap_mu += coef * pr
        return ((-0.5) ** n) * kap_mu


# thin wrappers kept for the validation harness (small s, direct expm is fine)
def assemble_phi(C, Gsp, legs):
    nlegs = len(legs)
    Phi = np.zeros((nlegs, nlegs), dtype=complex)
    Rcache = {}
    for p in range(nlegs):
        for q in range(p + 1, nlegs):
            lp, tp = legs[p]; lq, tq = legs[q]
            s = tq - tp; key = round(s, 12)
            if key not in Rcache:
                Rcache[key] = expm(s * Gsp) if s > 1e-14 else np.eye(Gsp.shape[0])
            Phi[p, q] = (C @ Rcache[key].T)[lp, lq]
            Phi[q, p] = -Phi[p, q]
    return Phi


def bond_corr(C, Gsp, points):
    pts = sorted(points, key=lambda p: p[1])
    legs = []
    for (x, t) in pts:
        legs.append((2 * x + 1, t)); legs.append((2 * x + 2, t))
    return (1j ** len(pts)) * pfaffian(assemble_phi(C, Gsp, legs))


# ---- validation ------------------------------------------------------------
def validate(models=((6, 1 / 3.), (8, 0.6), (8, 0.8), (10, 0.9)), verbose=True):
    """Closed-form builders (Gsp, C) vs brute-force 2^M extraction, and a 3-point
    bond correlator vs the exact operator route.  Returns the worst deviation
    over all models; run_checks.py gates on it."""
    import ff_extract_check as ff
    if verbose:
        print("== builder validation vs brute-force operators ==")
    worst = 0.0
    for (M, eta) in models:
        m, Cb, Gb, _ = ff.extract(M, eta)
        bnd = Gb.imag[1, 2]
        Gsp = build_Gsp(M, eta, boundary=bnd)
        C, _ = build_C(M, eta)
        eG = np.abs(Gsp - Gb).max()
        eC = np.abs(C - Cb).max()
        # and the correlators themselves
        pts3 = [(1, 0.0), (3, 0.3), (4, 0.7)]
        ex = m.bond_corr(pts3).real
        bf = bond_corr(C, Gsp, pts3).real
        worst = max(worst, eG, eC, abs(bf - ex))
        if verbose:
            print(f"  M={M:2d} eta={eta:.3f}: max|Gsp-brute|={eG:.1e}  max|C-brute|={eC:.1e}  "
                  f"3pt bond ff={bf:+.8f} exact={ex:+.8f} d={abs(bf-ex):.1e}")
    return worst


if __name__ == "__main__":
    worst = validate()
    ok = worst < 1e-12
    print(f"  [{'PASS' if ok else 'FAIL'}] worst deviation {worst:.2e} (< 1e-12)")
    sys.exit(0 if ok else 1)
