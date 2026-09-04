"""
pfaffian_verify.py
==================

Machine-precision check of the central finite-N theorem of glauber_paper.tex
(Cor. cor:multitime, eq. e:multipf):

    <mu_{x1}(t1) ... mu_{xk}(tk)>  =  i^k * Pf[ Phi ]_{2k x 2k}

where mu_x = sigma_x sigma_{x+1} = i a_{2x+1} a_{2x+2} is the bond (interface)
operator, a_alpha are the Jordan-Wigner Majorana operators, and Phi is the
2k x 2k antisymmetric matrix of the elementary TWO-TIME Majorana contractions
in the Gaussian (free-fermion) equilibrium state.

We build everything explicitly on a small OPEN chain (where Jordan-Wigner has no
boundary/winding term, so the identity is exact):

  * Pauli operators, Majoranas a_alpha  (verify Clifford algebra + mu_x = i a a),
  * Glauber generator L (reversible open-chain rates), symmetrised L~ = pi^{1/2} L pi^{-1/2},
    null vector |psi0> = sqrt(pi)  (verify L~ symmetric, L~|psi0>=0),
  * the contraction  Phi_{pq} = <psi0| a_{lp} e^{(tq-tp)L~} a_{lq} |psi0>  (p before q),
  * compare i^k Pf[Phi] with the exact operator-product correlator, at k = 2, 3, 4.

The 2-time case is the determinant; the 3-time case is the genuinely new content
(the two oriented triangles).  We also print the connected (cumulant) parts.

This module is BOTH a report (run it and read the numbers) and a gate: the
`*_deviations` / `fuzz_deviations` functions return the measured errors and are
consumed by run_checks.py, so a regression fails CI instead of being printed to a
terminal nobody is watching.

Run:  python pfaffian_verify.py [M]        (exit code 0 iff all checks pass)
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from scipy.linalg import expm

from glauber_exact_ring import gamma_from_eta
from multitime_interface_dpp import pfaffian


# Standard multi-time configurations used by every check below and by the gates
# in run_checks.py.  Bonds (0,1) and (3,4) are adjacent (they share a site), so
# the overlapping-bond corner is covered here and not only in the fuzz.
STD_TESTS = {
    2: [(0, 0.0), (2, 0.6)],
    3: [(0, 0.0), (2, 0.45), (3, 0.95)],
    4: [(0, 0.0), (1, 0.3), (3, 0.7), (4, 1.15)],
}

DEFAULT_ETA = 1.0 / 3.0


# ── Pauli operators in the sigma^z (config) basis ──────────────────────────────
I2 = np.eye(2)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)


def op_at(P, j, M):
    """Single-site Pauli P at site j (0=leftmost) on M sites, via Kronecker."""
    mats = [I2] * M
    mats[j] = P
    out = mats[0]
    for k in range(1, M):
        out = np.kron(out, mats[k])
    return out


def build_paulis(M):
    sx = [op_at(SX, j, M) for j in range(M)]
    sy = [op_at(SY, j, M) for j in range(M)]
    sz = [op_at(SZ, j, M) for j in range(M)]
    return sx, sy, sz


def build_majoranas(M, sx, sy, sz):
    """a_{2j} = (prod_{k<j} sx_k) sz_j,  a_{2j+1} = (prod_{k<j} sx_k) sy_j."""
    dim = 1 << M
    a = [None] * (2 * M)
    string = np.eye(dim, dtype=complex)
    for j in range(M):
        a[2 * j] = string @ sz[j]
        a[2 * j + 1] = string @ sy[j]
        string = string @ sx[j]          # extend string by sx_j for next sites
    return a


# ── Glauber generator on the OPEN chain (functions / row convention) ───────────
def build_generator_open(M, eta):
    """Reversible Glauber rates  w_j = 1 - sigma_j * tanh(K * h_j),  K=arctanh(eta),
    h_j = sum of EXISTING neighbours.  Interior: tanh(K*(+-2))=+-gamma (engine rule
    1+-gamma, 1); boundary: tanh(K*(+-1))=+-zeta (one neighbour).  This is the rate
    that satisfies detailed balance w.r.t. the open Ising measure exp(K sum s_j s_{j+1}).
    Returns L with (Lf)(x) = sum_y L[x,y](f(y)-f(x))."""
    K = np.arctanh(eta)
    dim = 1 << M
    spins = np.empty((dim, M), dtype=np.int8)
    for x in range(dim):
        for i in range(M):
            spins[x, i] = 1 if (x >> i) & 1 else -1
    L = np.zeros((dim, dim))
    for x in range(dim):
        s = spins[x]
        for j in range(M):
            h = 0
            if j > 0:
                h += s[j - 1]
            if j < M - 1:
                h += s[j + 1]
            rate = 1.0 - s[j] * np.tanh(K * h)
            y = x ^ (1 << j)
            L[x, y] += rate
    np.fill_diagonal(L, -L.sum(axis=1))
    return L, spins


def gibbs_open(spins, eta):
    """Open Ising chain measure pi ∝ exp(K sum_{j} s_j s_{j+1}), K=arctanh(eta)."""
    K = np.arctanh(eta)
    M = spins.shape[1]
    energy = np.zeros(spins.shape[0])
    for j in range(M - 1):
        energy += spins[:, j] * spins[:, j + 1]
    w = np.exp(K * energy)
    return w / w.sum()


# ── correlators ────────────────────────────────────────────────────────────────
class FermiModel:
    def __init__(self, M, eta):
        self.M = M
        self.eta = eta
        self.gamma = gamma_from_eta(eta)
        sx, sy, sz = build_paulis(M)
        self.sx, self.sy, self.sz = sx, sy, sz
        self.a = build_majoranas(M, sx, sy, sz)
        self.L, self.spins = build_generator_open(M, eta)
        self.pi = gibbs_open(self.spins, eta)
        sq = np.sqrt(self.pi)
        # symmetrised generator L~ = pi^{1/2} L pi^{-1/2}
        self.Lt = (sq[:, None] * self.L) / sq[None, :]
        self.psi0 = sq.copy()                     # null vector of L~
        self._Pcache = {}
        # bond operator mu_x = sigma^z_x sigma^z_{x+1} as a diagonal matrix
        self.mu = [self.sz[x] @ self.sz[x + 1] for x in range(M - 1)]

    def U(self, dt):
        key = round(float(dt), 10)
        if key not in self._Pcache:
            self._Pcache[key] = expm(key * self.Lt)
        return self._Pcache[key]

    # exact multi-time bond correlator via the operator product
    def bond_corr(self, points):
        """<mu_{x1}(t1) ... mu_{xk}(tk)>, points=[(x,t),...] time-sorted."""
        pts = sorted(points, key=lambda p: p[1])
        v = self.psi0.copy()                      # ket, build right-to-left
        prev_t = pts[-1][1]
        for (x, t) in reversed(pts):
            dt = prev_t - t
            if dt > 1e-14:
                v = self.U(dt) @ v
            v = self.mu[x] @ v
            prev_t = t
        # remaining evolution to the left end is trivial on <psi0|
        return complex(self.psi0 @ v)

    # elementary two-time Majorana contraction Phi between two legs
    def contraction(self, leg_p, t_p, leg_q, t_q):
        """<psi0| a_{leg_p} e^{(t_q-t_p)L~} a_{leg_q} |psi0>, assuming t_p <= t_q."""
        v = self.a[leg_q] @ self.psi0
        dt = t_q - t_p
        if dt > 1e-14:
            v = self.U(dt) @ v
        v = self.a[leg_p] @ v
        return complex(self.psi0 @ v)

    # Pfaffian prediction for the multi-time bond correlator
    def bond_corr_pfaffian(self, points):
        pts = sorted(points, key=lambda p: p[1])
        k = len(pts)
        legs = []                                  # (leg_index, time) in order
        for (x, t) in pts:
            legs.append((2 * x + 1, t))
            legs.append((2 * x + 2, t))
        n = 2 * k
        Phi = np.zeros((n, n), dtype=complex)
        for p in range(n):
            for q in range(p + 1, n):
                lp, tp = legs[p]
                lq, tq = legs[q]
                # p is earlier-or-equal in time order (legs already time-sorted)
                Phi[p, q] = self.contraction(lp, tp, lq, tq)
                Phi[q, p] = -Phi[p, q]
        return (1j ** k) * pfaffian(Phi), Phi


def connected_from_moments(moment_fn, k):
    """joint cumulant of k single-occurrence variables via moment-cumulant."""
    from math import factorial

    def parts(coll):
        coll = list(coll)
        if len(coll) == 1:
            yield [coll]; return
        first = coll[0]
        for sm in parts(coll[1:]):
            for i, s in enumerate(sm):
                yield sm[:i] + [[first] + s] + sm[i + 1:]
            yield [[first]] + sm
    tot = 0.0 + 0j
    for part in parts(range(k)):
        sz = len(part)
        coef = (-1) ** (sz - 1) * factorial(sz - 1)
        pr = 1.0 + 0j
        for b in part:
            pr *= moment_fn(tuple(b))
        tot += coef * pr
    return tot


# ── gate-facing measurements (return numbers; run_checks.py consumes these) ────
def pfaffian_deviations(m, tests=None):
    """|<mu..mu> - i^k Pf[Phi]| per k, on the exact 2^M operator route.
    This is the finite-N (pre-scaling) identity e:multipf / cor:multitime."""
    tests = tests or STD_TESTS
    return {k: abs(m.bond_corr(pts) - m.bond_corr_pfaffian(pts)[0])
            for k, pts in tests.items()}


def algebra_deviations(m):
    """Structural prerequisites of the identity: Clifford algebra, mu = i a a,
    reversibility of L~, and stationarity of |psi0>."""
    dim = 1 << m.M
    cliff = 0.0
    for al in range(2 * m.M):
        for be in range(2 * m.M):
            anti = m.a[al] @ m.a[be] + m.a[be] @ m.a[al]
            cliff = max(cliff, np.abs(anti - 2 * (al == be) * np.eye(dim)).max())
    return {
        "clifford": cliff,
        "mu_is_i_a_a": max(np.abs(1j * m.a[2 * x + 1] @ m.a[2 * x + 2] - m.mu[x]).max()
                           for x in range(m.M - 1)),
        "L~_symmetric": np.abs(m.Lt - m.Lt.T).max(),
        "L~_psi0_zero": np.abs(m.Lt @ m.psi0).max(),
    }


def fuzz_deviations(seed=20260812, trials_per_k=25, ks=(2, 3, 4),
                    models=((6, 1.0 / 3.0), (7, 0.6), (8, 0.9), (8, 0.45))):
    """Randomised sweep of the finite-N identity over the corners the fixed
    configurations above do not reach: coinciding times, the SAME bond observed
    at two times, adjacent/overlapping bonds, identical (bond, time) pairs,
    several (M, eta).

    Checks the operator route AND the closed-form free-fermion engine against the
    same exact ground truth.  Returns (worst_deviation, coverage_counts).
    """
    import ff_engine as fe
    rng = np.random.default_rng(seed)
    worst = 0.0
    cov = {"configs": 0, "equal_times": 0, "repeated_bond": 0, "both": 0,
           "identical_pair": 0}
    for M, eta in models:
        m = FermiModel(M, eta)
        ff = fe.FF(M, eta)
        for k in ks:
            for trial in range(trials_per_k):
                xs = rng.integers(0, M - 1, size=k)
                mode = trial % 3          # 0 generic, 1 some ties, 2 all equal
                if mode == 0:
                    ts = np.round(rng.uniform(0, 1.2, size=k), 3)
                elif mode == 1:
                    base = np.round(rng.uniform(0, 1.2, size=max(1, k - 1)), 3)
                    ts = np.array([base[i % len(base)] for i in range(k)])
                else:
                    ts = np.full(k, round(float(rng.uniform(0, 1.2)), 3))
                pts = [(int(x), float(t)) for x, t in zip(xs, ts)]
                # identical (bond, time) pairs are KEPT: mu^2 = 1 on the operator
                # side and the 2k-leg Pfaffian degenerates consistently
                # (<a_l a_l> = 1), so the identity holds there too.
                exact = m.bond_corr(pts)
                worst = max(worst,
                            abs(exact - m.bond_corr_pfaffian(pts)[0]),
                            abs(exact - ff.bond_corr(pts)))
                eqt = len({t for _, t in pts}) < k
                rep = len({x for x, _ in pts}) < k
                cov["configs"] += 1
                cov["equal_times"] += int(eqt)
                cov["repeated_bond"] += int(rep)
                cov["both"] += int(eqt and rep)
                cov["identical_pair"] += int(len(set(pts)) < k)
    return worst, cov


def kernel_report(m):
    """Print and interpret the elementary kernel the Pfaffian is built from."""
    print("=" * 70)
    print("THE KERNEL  Phi_{ab}(t,t') = <psi0| a_a(t) a_b(t') |psi0>")
    print("=" * 70)
    M = m.M
    # spin two-point D_ij(t) = <sigma_i(t) sigma_j(0)>  (the dynamical building block)
    def D(i, j, t):
        v = m.sz[j] @ m.psi0
        if t > 1e-14:
            v = m.U(t) @ v
        v = m.sz[i] @ v
        return float((m.psi0 @ v).real)
    t = 0.6
    print(f"\n(1) Static (equal-time) input: <sigma_0 sigma_j> should be zeta^j = "
          f"{m.eta:.4f}^j")
    print("    j :   exact        zeta^j")
    for j in range(min(M, 4)):
        print(f"    {j} : {D(0,j,0.0):+.6f}    {m.eta**j:+.6f}")

    # the 4x4 kernel block for two bonds (x at t=0, y at time t)
    x, y = 0, 2
    legs = [(2*x+1, 0.0), (2*x+2, 0.0), (2*y+1, t), (2*y+2, t)]
    Phi = np.zeros((4, 4), dtype=complex)
    for p in range(4):
        for q in range(p+1, 4):
            lp, tp = legs[p]; lq, tq = legs[q]
            Phi[p, q] = m.contraction(lp, tp, lq, tq)
            Phi[q, p] = -Phi[p, q]
    print(f"\n(2) The 4x4 antisymmetric kernel  Phi  for bond {x}@t=0 and bond {y}@t={t}")
    print("    legs = [a_{2x+1}@0, a_{2x+2}@0 | a_{2y+1}@t, a_{2y+2}@t]")
    np.set_printoptions(precision=4, suppress=True)
    print("    Re Phi =\n", Phi.real)
    print("    Im Phi =\n", Phi.imag)
    print(f"    Pf[Phi] = {pfaffian(Phi):+.6f}   "
          f"i^2 Pf = {(1j**2*pfaffian(Phi)).real:+.6f}  "
          f"= <mu_{x}(0) mu_{y}({t})> = {m.bond_corr([(x,0.0),(y,t)]).real:+.6f}")

    # diagonal block = equal-time: Phi[0,1] = <a_{2x+1} a_{2x+2}> = mu_x / i = -i*zeta
    print(f"\n    equal-time block Phi[0,1] = {Phi[0,1]:+.4f}  (= -i*<mu_x> = -i*zeta = "
          f"{-1j*m.eta:+.4f})")

    # connected two-time wall correlator: SIGN is the diagnostic
    mu_c = (m.bond_corr([(x,0.0),(y,t)]) - m.bond_corr([(x,0.0)])*m.bond_corr([(y,t)]))
    print(f"\n(3) Connected two-time wall correlator <mu_{x}(0)mu_{y}({t})>_c = {mu_c.real:+.5e}")
    print(f"    sign = {'POSITIVE' if mu_c.real>0 else 'negative'}: a pure number-conserving")
    print("    (determinantal) wall process forces ANTI-correlation; a positive value")
    print("    is the fingerprint of the ANOMALOUS (pair create/annihilate) contractions.")
    print("=" * 70 + "\n")


def verify_2time_determinant(m):
    """Check eq.(2time):  <mu_x(t) mu_y(0)>_c = det[[D_xy,D_x,y+1],[D_x+1,y,D_x+1,y+1]],
    D_ij(t) = <sigma_i(t) sigma_j(0)>."""
    print("=" * 70)
    print("CHECK eq.(2time): connected bond corr = det of spin two-point D")
    print("=" * 70)

    def D(i, j, t):
        v = m.sz[j] @ m.psi0
        if t > 1e-14:
            v = m.U(t) @ v
        v = m.sz[i] @ v
        return float((m.psi0 @ v).real)

    for (x, y, t) in [(0, 2, 0.6), (1, 3, 0.4), (0, 1, 0.8)]:
        det = (D(x, y, t) * D(x+1, y+1, t) - D(x, y+1, t) * D(x+1, y, t))
        full = m.bond_corr([(x, 0.0), (y, t)]).real
        conn = full - m.bond_corr([(x, 0.0)]).real * m.bond_corr([(y, t)]).real
        print(f"  x={x},y={y},t={t}:  det(D-block)={det:+.8f}   "
              f"<mu mu>_c={conn:+.8f}   diff={abs(det-conn):.1e}")
    print()


def D_universality_gaps(m, tests=None):
    """NEGATIVE CONTROL.  Is the SCALAR spin two-point D the full kernel, or is the
    matrix-valued (anomalous) Majorana contraction Phi genuinely needed?  Tests
        <mu_{x1}(t1)...mu_{xk}(tk)>  =?  Pf[ D-matrix of the 2k spin legs ].
    A MISMATCH at every k is the expected, wanted outcome: it is what rules out a
    number-conserving (determinantal) description of the wall process.

    Returns k -> (exact, Pf[D], |signed gap|, |gap between magnitudes|).  Both
    gaps are reported so that a mismatch cannot be dismissed as a Pfaffian sign
    convention: the magnitudes differ too.
    """
    tests = tests or STD_TESTS

    def D(i, j, t):
        v = m.sz[j] @ m.psi0
        if t > 1e-14:
            v = m.U(t) @ v
        v = m.sz[i] @ v
        return float((m.psi0 @ v).real)

    out = {}
    for k, pts in tests.items():
        legs = []
        for (x, t) in pts:
            legs.append((x, t)); legs.append((x + 1, t))
        n = 2 * k
        Mt = np.zeros((n, n))
        for p in range(n):
            for q in range(p + 1, n):
                ip, tp = legs[p]; iq, tq = legs[q]
                Mt[p, q] = D(ip, iq, abs(tp - tq))
                Mt[q, p] = -Mt[p, q]
        pf = pfaffian(Mt).real
        exact = m.bond_corr(pts).real
        out[k] = (exact, pf, abs(pf - exact), abs(abs(pf) - abs(exact)))
    return out


def D_universality_test(m):
    print("=" * 70)
    print("IS THE SCALAR SPIN TWO-POINT D THE WHOLE KERNEL?  (Pf[D] vs exact)")
    print("  [negative control: MISMATCH at every k is the wanted outcome]")
    print("=" * 70)
    for k, (exact, pf, gap, maggap) in D_universality_gaps(m).items():
        match = "MATCH" if gap < 1e-9 else "*** MISMATCH (expected) ***"
        print(f"  k={k}: exact<mu..mu>={exact:+.8f}   Pf[D]={pf:+.8f}   "
              f"|diff|={gap:.2e}  (magnitudes differ by {maggap:.2e})   {match}")
    print()


def spin_wick_gaps(m):
    """Does the multi-time SPIN correlator equal Pf of the spin two-point D?
       <sigma_{i1}(t1)...sigma_{i2m}(t2m)>  =?  Pf[ D_{pq} ].

    Outcome (this is why the theory is stated for BONDS and not for spins): the
    EQUAL-TIME spin correlator is the Pfaffian of D exactly -- the static Ising
    Wick theorem -- but the MULTI-TIME one is not.  sigma_x is a Jordan-Wigner
    string operator, not a Majorana bilinear, so it is not quasi-free in time;
    mu_x = sigma_x sigma_{x+1} is.  Returns label -> (exact, Pf[D], |gap|).
    """
    def Dt(s_i, t_i, s_j, t_j):
        if t_i < t_j:
            s_i, t_i, s_j, t_j = s_j, t_j, s_i, t_i
        v = m.sz[s_j] @ m.psi0
        if t_i - t_j > 1e-14:
            v = m.U(t_i - t_j) @ v
        v = m.sz[s_i] @ v
        return float((m.psi0 @ v).real)

    def spin_corr(legs):
        pts = sorted(legs, key=lambda p: p[1])
        v = m.psi0.copy(); prev = pts[-1][1]
        for (s, t) in reversed(pts):
            dt = prev - t
            if dt > 1e-14:
                v = m.U(dt) @ v
            v = m.sz[s] @ v; prev = t
        return float((m.psi0 @ v).real)

    cases = [
        ("4 spins, equal time", [(0, 0.0), (1, 0.0), (3, 0.0), (4, 0.0)]),
        ("4 spins, multi-time", [(0, 0.0), (1, 0.3), (3, 0.6), (4, 0.9)]),
        ("6 spins, multi-time", [(0, 0.0), (1, 0.2), (2, 0.5), (4, 0.7), (5, 0.9), (3, 1.1)]),
    ]
    out = {}
    for label, legs in cases:
        n = len(legs)
        order = sorted(range(n), key=lambda k: legs[k][1])
        ol = [legs[k] for k in order]
        Mt = np.zeros((n, n))
        for p in range(n):
            for q in range(p + 1, n):
                (sp, tp), (sq, tq) = ol[p], ol[q]
                Mt[p, q] = Dt(sp, tp, sq, tq); Mt[q, p] = -Mt[p, q]
        pf = pfaffian(Mt).real
        exact = spin_corr(legs)
        out[label] = (exact, pf, abs(pf - exact))
    return out


def spin_wick_test(m):
    print("=" * 70)
    print("MULTI-TIME SPIN WICK:  <sigma...sigma>  vs  Pf[D]")
    print("  [expected: equal time MATCHES (static Ising Wick); multi-time FAILS")
    print("   -- sigma is a JW string, not a bilinear; this is why the theory")
    print("   is stated for BONDS mu = sigma sigma]")
    print("=" * 70)
    for label, (exact, pf, gap) in spin_wick_gaps(m).items():
        match = "MATCH" if gap < 1e-9 else "*** MISMATCH (expected) ***"
        print(f"  {label:22s}: exact={exact:+.8f}  Pf[D]={pf:+.8f}  "
              f"|diff|={gap:.2e}  {match}")
    print()


def main(M=6, eta=DEFAULT_ETA):
    m = FermiModel(M, eta)
    print(f"OPEN chain M={M}, eta=zeta={eta:.4f}, gamma={m.gamma:.6f}\n")

    # ── sanity: Clifford algebra, mu = i a a, reversibility ──
    alg = algebra_deviations(m)
    print(f"  Clifford {{a,a}}=2d   : max err {alg['clifford']:.2e}")
    print(f"  mu_x = i a a        : max err {alg['mu_is_i_a_a']:.2e}")
    print(f"  L~ symmetric (revers): max err {alg['L~_symmetric']:.2e}")
    print(f"  L~|psi0>=0 (station) : max err {alg['L~_psi0_zero']:.2e}\n")

    # ── k = 2, 3, 4 : exact operator product vs i^k Pf[Phi] ──
    for k, pts in STD_TESTS.items():
        exact = m.bond_corr(pts)
        pf, _ = m.bond_corr_pfaffian(pts)
        # connected interface (eta) cumulant: eta=(1-mu)/2 -> kappa_eta = (-1/2)^k kappa_mu
        kap_mu = connected_from_moments(
            lambda S: m.bond_corr([pts[i] for i in S]), k)
        kap_eta = ((-0.5) ** k) * kap_mu
        bonds = ",".join(str(p[0]) for p in pts)
        times = ",".join(f"{p[1]:g}" for p in pts)
        print(f"  k={k}  bonds=({bonds}) times=({times})")
        print(f"     exact <mu..mu>   = {exact.real:+.10f}  (Im {exact.imag:+.1e})")
        print(f"     i^k Pf[Phi]      = {pf.real:+.10f}  (Im {pf.imag:+.1e})")
        print(f"     |exact - Pf|     = {abs(exact - pf):.3e}")
        print(f"     connected <eta..eta>_c = {kap_eta.real:+.6e}\n")
    return m


if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    eta = float(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_ETA
    m = main(M, eta)                      # one model, reused by every report
    verify_2time_determinant(m)
    D_universality_test(m)
    spin_wick_test(m)
    kernel_report(m)

    # ── gate summary (machine-checkable; run_checks.py uses the same functions)
    alg = algebra_deviations(m)
    dev = pfaffian_deviations(m)
    Dgaps = D_universality_gaps(m)
    wick = spin_wick_gaps(m)
    gates = [
        ("algebra prerequisites < 1e-12", max(alg.values()) < 1e-12),
        ("<mu..mu> = i^k Pf[Phi], k=2,3,4 < 1e-12", max(dev.values()) < 1e-12),
        ("negative control: Pf[D] FAILS at every k (signed AND magnitude gap)",
         min(g[2] for g in Dgaps.values()) > 1e-8 and min(g[3] for g in Dgaps.values()) > 1e-8),
        ("equal-time spin Wick holds < 1e-12", wick["4 spins, equal time"][2] < 1e-12),
        ("negative control: multi-time spin Wick FAILS",
         min(g[2] for lbl, g in wick.items() if "multi-time" in lbl) > 1e-8),
    ]
    print("=" * 70)
    n_ok = 0
    for name, ok in gates:
        n_ok += bool(ok)
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"  {n_ok}/{len(gates)} gates passed.")
    sys.exit(0 if n_ok == len(gates) else 1)
