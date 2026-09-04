"""
multitime_interface_dpp.py
==========================

Question (from interface_dpp note): the TWO-time connected interface (domain-wall)
correlator of the Glauber-Ising chain is a single 2x2 determinant of one-wall
propagators G, i.e. an *extended determinantal* structure.  Does a 3-time / 4-time
interface correlator exist, and does it STILL carry a determinantal structure?

This script answers it on the exact finite-ring generator (gold standard, the same
build_generator used in exactdiag_check.py).  Strategy:

  * Compute genuine *connected* multi-time wall correlators (joint cumulants of the
    0/1 wall indicators eta_x(t) = (1 - sigma_x sigma_{x+1})/2) by exact
    diagonalisation -- no ansatz, machine precision.

  * Test the rigid structure a DETERMINANTAL point process must obey.  For a DPP
    with kernel K the truncated (connected) correlation functions are sums over
    single cycles of kernel entries:
        rho_2^T(1,2)   = -K12 K21
        rho_3^T(1,2,3) = +(K12 K23 K31 + K13 K32 K21)
        rho_4^T        = -(sum over 6 directed 4-cycles)
    These imply *gauge-invariant* constraints that do not depend on knowing K:
      (3pt, real kernel) discriminant  D3 = (rho3T)^2 - 4 * C12 C23 C31  >= 0,
          with C_ij = -rho2T(ij)   [must hold for ANY real kernel].
      (3pt, complex kernel) the two 3-cycles are roots of
          z^2 - rho3T z + C12 C23 C31 = 0  (always solvable -> 3pt alone is weak).
      (4pt) over-determined consistency of the cycle data -> decisive even for
          complex kernels.
    The Glauber generator is reversible (detailed balance) so its natural extended
    kernel is real; a negative D3 therefore rules out the determinantal structure.

  * Cross-check against the PFAFFIAN alternative (annihilating walks are Pfaffian
    point processes, Tribe-Zaboronski): build the 2x2-matrix extended kernel from
    the measured one- and two-point data and verify it reproduces rho_3^T, rho_4^T.

Run:  python multitime_interface_dpp.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from math import factorial
from itertools import combinations
from scipy.linalg import expm

from glauber_exact_ring import gamma_from_eta, two_time_stationary


# ────────────────────────────────────────────────────────────────────────────
# Exact generator (engine rate convention), identical to exactdiag_check.py
# ────────────────────────────────────────────────────────────────────────────
def build_generator(l, gamma):
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
            L[x, y] += rate
    np.fill_diagonal(L, -L.sum(axis=1))               # (Lf)(x)=sum_y L[x,y](f(y)-f(x))
    return L, spins


def ising_measure(spins, eta):
    K = np.arctanh(eta)
    energy = np.sum(spins * np.roll(spins, -1, axis=1), axis=1)
    w = np.exp(K * energy.astype(float))
    return w / w.sum()


# ────────────────────────────────────────────────────────────────────────────
# Multi-time machinery
# ────────────────────────────────────────────────────────────────────────────
class Model:
    def __init__(self, l, eta):
        self.l = l
        self.eta = eta
        self.gamma = gamma_from_eta(eta)
        self.L, self.spins = build_generator(l, self.gamma)
        self.mu = ising_measure(self.spins, eta)
        # wall (bond) indicators eta_x = (1 - s_x s_{x+1})/2,  x = 0..l-1
        self.wall = np.empty((l, 1 << l))
        for x in range(l):
            self.wall[x] = 0.5 * (1.0 - self.spins[:, x].astype(float) *
                                  self.spins[:, (x + 1) % l].astype(float))
        # spin observable s_x for self-test
        self.sig = self.spins.astype(float)
        self._Pcache = {}

    def P(self, dt):
        key = round(float(dt), 10)
        if key not in self._Pcache:
            self._Pcache[key] = expm(key * self.L)
        return self._Pcache[key]

    def moment_general(self, obs_points):
        """E[ prod_i f_{a_i}(t_i) ] for general diagonal observables.
        obs_points: list of (vec, t); vec is a length-nstates observable."""
        pts = sorted(obs_points, key=lambda p: p[1])
        g = pts[-1][0].copy()
        for k in range(len(pts) - 2, -1, -1):
            dt = pts[k + 1][1] - pts[k][1]
            if dt > 1e-14:
                g = pts[k][0] * (self.P(dt) @ g)
            else:
                g = pts[k][0] * g
        return float(self.mu @ g)

    def wall_moment(self, points):
        """E[ prod eta_{b_i}(t_i) ],  points = [(bond, t), ...]."""
        return self.moment_general([(self.wall[b], t) for (b, t) in points])

    def wall_cumulant(self, points):
        """Joint cumulant (= truncated correlation = connected correlator) of the
        0/1 wall indicators eta_{b_i}(t_i)."""
        n = len(points)
        idx = list(range(n))
        total = 0.0
        for part in _set_partitions(idx):
            size = len(part)
            coef = (-1) ** (size - 1) * factorial(size - 1)
            prod = 1.0
            for block in part:
                prod *= self.wall_moment([points[i] for i in block])
            total += coef * prod
        return total


def _set_partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        yield [collection]
        return
    first = collection[0]
    for smaller in _set_partitions(collection[1:]):
        for i, subset in enumerate(smaller):
            yield smaller[:i] + [[first] + subset] + smaller[i + 1:]
        yield [[first]] + smaller


def cumulant_from_moments(moment_fn, n):
    """Joint cumulant of n variables given a moment_fn(subset_indices)->moment."""
    total = 0.0
    for part in _set_partitions(list(range(n))):
        size = len(part)
        coef = (-1) ** (size - 1) * factorial(size - 1)
        prod = 1.0
        for block in part:
            prod *= moment_fn(tuple(block))
        total += coef * prod
    return total


# ────────────────────────────────────────────────────────────────────────────
# SINGLE-TIME (possibly non-equilibrium) wall point process
# ────────────────────────────────────────────────────────────────────────────
class SingleTimeState:
    """Wall (domain-wall) point process in a fixed single-time distribution `pi`.
    Used for the T=0 quench (pure annihilating walks): the cleanest setting in
    which the determinantal-vs-Pfaffian distinction is sharp, since at a single
    time a real reversible dynamics has a symmetric correlation kernel."""
    def __init__(self, l, pi, spins):
        self.l = l
        self.pi = pi
        self.wall = np.empty((l, 1 << l))
        for x in range(l):
            self.wall[x] = 0.5 * (1.0 - spins[:, x].astype(float) *
                                  spins[:, (x + 1) % l].astype(float))

    def moment(self, bonds):
        v = np.ones(len(self.pi))
        for b in bonds:
            v = v * self.wall[b]
        return float(self.pi @ v)

    def cumulant(self, bonds):
        bonds = list(bonds)
        return cumulant_from_moments(
            lambda idx: self.moment([bonds[i] for i in idx]), len(bonds))


def quench_state(l, t):
    """pi_t after a T=0 (gamma=1, pure annihilation) quench from full disorder."""
    Lq, spins = build_generator(l, 1.0)
    pi0 = np.full(1 << l, 1.0 / (1 << l))
    pit = pi0 @ expm(t * Lq)
    return SingleTimeState(l, pit, spins), spins


def single_time_determinantal_test(l=12, t=0.25):
    """Sharp test.  For a single-time DETERMINANTAL process with a REAL SYMMETRIC
    kernel (the natural kernel of a real reversible dynamics):

        rho2T(i,j) = -K_ij^2  <= 0                         (anti-correlation)
        rho3T(i,j,k) = 2 K_ij K_jk K_ki
            =>  (rho3T)^2 = 4 (-rho2T_ij)(-rho2T_jk)(-rho2T_ki)   EXACTLY.

    Pure annihilating random walks (Tribe-Zaboronski) are Pfaffian, NOT
    determinantal, so this equality must FAIL.  We report the relative defect
        R = (rho3T)^2 / (4 prod(-rho2T))  -  1 ;   R==0 <=> determinantal."""
    st, spins = quench_state(l, t)
    print(f"── SINGLE-TIME determinantal test  (T=0 quench, l={l}, t={t}) ──")
    # equal-time wall density and pair anti-correlation
    p = st.moment([0])
    print(f"   wall density p(t) = {p:.5f}")
    # check rho2T <= 0 for adjacent / near bonds
    triples = []
    for b1 in range(l):
        for b2 in range(b1 + 1, l):
            for b3 in range(b2 + 1, l):
                triples.append((b1, b2, b3))
    worst_R = None
    n_bad = 0
    shown = 0
    pos_pair = 0
    for (b1, b2, b3) in triples:
        k12 = st.cumulant([b1, b2]); k23 = st.cumulant([b2, b3]); k13 = st.cumulant([b1, b3])
        if max(k12, k23, k13) > 1e-12:
            pos_pair += 1
        k3 = st.cumulant([b1, b2, b3])
        denom = 4.0 * (-k12) * (-k23) * (-k13)
        if abs(denom) < 1e-30:
            continue
        R = (k3 * k3) / denom - 1.0
        if abs(R) > 1e-6:
            n_bad += 1
        if worst_R is None or abs(R) > abs(worst_R[0]):
            worst_R = (R, (b1, b2, b3), k3, (k12, k23, k13))
        if abs(R) > 1e-3 and shown < 8:
            print(f"   bonds {(b1,b2,b3)}: (rho3T)^2/(4∏(-rho2T)) - 1 = {R:+.4e}   "
                  f"[rho3T={k3:+.3e}, rho2T={k12:+.2e},{k23:+.2e},{k13:+.2e}]")
            shown += 1
    print(f"   pairs with rho2T>0 (would already break a Hermitian DPP): "
          f"{pos_pair}")
    print(f"   triples failing the determinantal equality (|R|>1e-6): "
          f"{n_bad} / {len(triples)}")
    Rw, bw, k3w, k2w = worst_R
    print(f"   largest relative defect  R = {Rw:+.4e}  at bonds {bw}")
    print(f"   => equality holds (determinantal)?  {'YES' if abs(Rw)<1e-6 else 'NO'}\n")
    return worst_R


# ────────────────────────────────────────────────────────────────────────────
# Self-test: 2-time spin correlator vs the validated FFT spectral formula
# ────────────────────────────────────────────────────────────────────────────
def self_test(m: Model):
    print("── self-test: <sigma_0(0) sigma_n(t)> exact-diag vs FFT spectral formula ──")
    for t in (0.3, 0.9):
        G_ed = np.array([m.moment_general([(m.sig[:, 0], 0.0), (m.sig[:, n], t)])
                         for n in range(m.l)])
        G_th = two_time_stationary(m.l, m.eta, t)
        err = np.abs(G_ed - G_th).max()
        print(f"   t={t}:  max|exactdiag - formula| = {err:.2e}")
    print()


# ────────────────────────────────────────────────────────────────────────────
# DETERMINANTAL test (kernel-agnostic)
# ────────────────────────────────────────────────────────────────────────────
def determinantal_3pt_test(m: Model, times, verbose=True):
    """For every triple of distinct (bond,time) points with three DISTINCT times,
    compute the real-kernel determinantal discriminant
        D3 = (rho3T)^2 - 4 * C12 C23 C31,   C_ij = -rho2T(ij).
    A real determinantal point process REQUIRES D3 >= 0.  Report the minimum."""
    t1, t2, t3 = times
    bonds = range(m.l)
    worst = None
    n_neg = 0
    n_tot = 0
    examples = []
    for b1 in bonds:
        for b2 in bonds:
            for b3 in bonds:
                pts = [(b1, t1), (b2, t2), (b3, t3)]
                k2_12 = m.wall_cumulant([pts[0], pts[1]])
                k2_23 = m.wall_cumulant([pts[1], pts[2]])
                k2_13 = m.wall_cumulant([pts[0], pts[2]])
                C12, C23, C13 = -k2_12, -k2_23, -k2_13
                k3 = m.wall_cumulant(pts)
                Q = C12 * C23 * C13
                D3 = k3 * k3 - 4.0 * Q
                n_tot += 1
                if D3 < -1e-12:
                    n_neg += 1
                if worst is None or D3 < worst[0]:
                    worst = (D3, pts, k3, (k2_12, k2_23, k2_13), Q)
                if D3 < -1e-9 and len(examples) < 6:
                    examples.append((D3, pts, k3, (k2_12, k2_23, k2_13), Q))
    if verbose:
        print(f"── DETERMINANTAL 3-point test  (times {t1},{t2},{t3}; "
              f"{n_tot} triples) ──")
        print(f"   triples with D3 < 0 (real-DPP forbidden): {n_neg} / {n_tot}")
        D3w, pw, k3w, k2w, Qw = worst
        print(f"   most negative discriminant D3 = {D3w:+.6e}")
        print(f"       at bonds/times {pw}")
        print(f"       rho3T = {k3w:+.6e},  rho2T(12,23,13) = "
              f"({k2w[0]:+.3e},{k2w[1]:+.3e},{k2w[2]:+.3e}),  4*C12C23C13 = {4*Qw:+.3e}")
        if examples:
            print("   sample DPP-violating configurations:")
            for D3e, pe, k3e, k2e, Qe in examples:
                print(f"       {pe}: D3={D3e:+.3e}  rho3T={k3e:+.3e}  "
                      f"4Q={4*Qe:+.3e}")
        print()
    return worst, n_neg, n_tot


# ────────────────────────────────────────────────────────────────────────────
# Pfaffian (Parlett-Reid, O(n^3)) and a self-test
# ────────────────────────────────────────────────────────────────────────────
def pfaffian(A):
    """Pfaffian of a (real or complex) antisymmetric matrix via Parlett-Reid."""
    A = np.array(A, dtype=complex)
    N = A.shape[0]
    if N % 2 == 1:
        return 0.0 + 0.0j
    if N == 0:
        return 1.0 + 0.0j
    A = A.copy()
    pf = 1.0 + 0.0j
    for k in range(0, N - 1, 2):
        # pivot: largest |A[i,k]| for i>k
        col = np.abs(A[k + 1:, k])
        kp = k + 1 + int(np.argmax(col))
        if kp != k + 1:
            A[[k + 1, kp], :] = A[[kp, k + 1], :]
            A[:, [k + 1, kp]] = A[:, [kp, k + 1]]
            pf = -pf
        piv = A[k, k + 1]
        if piv == 0:
            return 0.0 + 0.0j
        pf *= piv
        if k + 2 < N:
            tau = A[k, k + 2:] / piv          # row-k tail
            col2 = A[k + 2:, k + 1].copy()    # col-(k+1) tail
            A[k + 2:, k + 2:] += np.outer(tau, col2) - np.outer(col2, tau)
    return pf


def _pfaffian_self_test():
    rng = np.random.default_rng(0)
    worst = 0.0
    for N in (2, 4, 6, 8):
        for _ in range(20):
            M = rng.standard_normal((N, N))
            A = M - M.T
            pf = pfaffian(A)
            worst = max(worst, abs(pf * pf - np.linalg.det(A)))
    print(f"── pfaffian self-test: max |Pf(A)^2 - det(A)| = {worst:.2e} ──\n")


# ────────────────────────────────────────────────────────────────────────────
# Over-determined structure fits: does a DETERMINANTAL / PFAFFIAN kernel exist?
# ────────────────────────────────────────────────────────────────────────────
from scipy.optimize import least_squares


def _all_subset_moments(moment_fn, n):
    """Return list of (subset_tuple, moment) for all nonempty subsets of range(n)."""
    out = []
    for r in range(1, n + 1):
        for S in combinations(range(n), r):
            out.append((S, moment_fn(S)))
    return out


def fit_determinantal(moment_fn, n, complexK=True, restarts=6, seed=0):
    """Try to find a kernel K (n x n) with det K[S,S] = moment(S) for ALL subsets S.
    Returns the minimal RMS residual found.  ~0 => determinantal; bounded away => not."""
    subs = _all_subset_moments(moment_fn, n)
    M = len(subs)
    rng = np.random.default_rng(seed)

    def unpack(v):
        K = np.zeros((n, n), dtype=complex)
        idx = 0
        for i in range(n):
            K[i, i] = v[idx]; idx += 1               # real diagonal
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if complexK:
                    K[i, j] = v[idx] + 1j * v[idx + 1]; idx += 2
                else:
                    K[i, j] = v[idx]; idx += 1
        return K

    nvar = n + (2 if complexK else 1) * n * (n - 1)

    def resid(v):
        K = unpack(v)
        r = np.empty(M, dtype=complex)
        for q, (S, mom) in enumerate(subs):
            sub = K[np.ix_(S, S)]
            r[q] = np.linalg.det(sub) - mom
        return np.concatenate([r.real, r.imag])

    best = np.inf
    for s in range(restarts):
        v0 = 0.1 * rng.standard_normal(nvar)
        # seed diagonal with the measured densities
        for i in range(n):
            v0[i] = subs[i][1]
        sol = least_squares(resid, v0, method="lm" if M >= nvar else "trf",
                            max_nfev=4000)
        rms = np.sqrt(np.mean(sol.fun ** 2))
        best = min(best, rms)
    return best, M, nvar


def fit_pfaffian(moment_fn, n, complexK=False, restarts=6, seed=0):
    """Try to find a 2x2-matrix antisymmetric kernel with Pf = moment(S) for ALL S."""
    subs = _all_subset_moments(moment_fn, n)
    M = len(subs)
    rng = np.random.default_rng(seed)
    pairs = list(combinations(range(n), 2))

    def unpack(v):
        # diagonal antisym 2x2 blocks: param d_i  -> [[0,d],[-d,0]]
        d = v[:n]
        blocks = {}
        idx = n
        for (i, j) in pairs:
            if complexK:
                B = (v[idx:idx + 4] + 1j * v[idx + 4:idx + 8]).reshape(2, 2)
                idx += 8
            else:
                B = v[idx:idx + 4].reshape(2, 2); idx += 4
            blocks[(i, j)] = B
        return d, blocks

    nvar = n + (8 if complexK else 4) * len(pairs)

    def assemble(S, d, blocks):
        m = len(S)
        A = np.zeros((2 * m, 2 * m), dtype=complex)
        for a in range(m):
            di = d[S[a]]
            A[2 * a, 2 * a + 1] = di
            A[2 * a + 1, 2 * a] = -di
        for a in range(m):
            for b in range(a + 1, m):
                i, j = S[a], S[b]
                B = blocks[(i, j)] if (i, j) in blocks else -blocks[(j, i)].T
                A[2 * a:2 * a + 2, 2 * b:2 * b + 2] = B
                A[2 * b:2 * b + 2, 2 * a:2 * a + 2] = -B.T
        return A

    def resid(v):
        d, blocks = unpack(v)
        r = np.empty(M, dtype=complex)
        for q, (S, mom) in enumerate(subs):
            r[q] = pfaffian(assemble(S, d, blocks)) - mom
        return np.concatenate([r.real, r.imag])

    best = np.inf
    for s in range(restarts):
        v0 = 0.1 * rng.standard_normal(nvar)
        for i in range(n):
            v0[i] = subs[i][1]      # seed diagonal block with density
        sol = least_squares(resid, v0, method="lm" if M >= nvar else "trf",
                            max_nfev=6000)
        rms = np.sqrt(np.mean(sol.fun ** 2))
        best = min(best, rms)
    return best, M, nvar


def structure_fit_report(label, moment_fn, n):
    print(f"── STRUCTURE FIT ({label}), n={n} points, "
          f"{2**n - 1} subset-moments ──")
    det_r_real, Mr, vr = fit_determinantal(moment_fn, n, complexK=False)
    det_r_cplx, _, vc = fit_determinantal(moment_fn, n, complexK=True)
    pf_r, _, vp = fit_pfaffian(moment_fn, n, complexK=False)
    print(f"   DETERMINANTAL (real kernel,    {vr:3d} params): "
          f"min RMS residual = {det_r_real:.3e}")
    print(f"   DETERMINANTAL (complex kernel, {vc:3d} params): "
          f"min RMS residual = {det_r_cplx:.3e}")
    print(f"   PFAFFIAN      (real 2x2 kernel,{vp:4d} params): "
          f"min RMS residual = {pf_r:.3e}")
    verdict = "PFAFFIAN (determinant fails)" if (pf_r < 1e-7 < det_r_cplx) \
        else ("determinantal" if det_r_cplx < 1e-7 else "neither fit converged")
    print(f"   => {verdict}\n")
    return det_r_real, det_r_cplx, pf_r


# ────────────────────────────────────────────────────────────────────────────
# Display: genuine connected 3-time and 4-time EQUILIBRIUM interface correlators
# ────────────────────────────────────────────────────────────────────────────
def show_multitime_correlators(m: Model, times3, times4):
    print("── connected 3-time interface correlators  "
          f"<eta(t1)eta(t2)eta(t3)>_c  (times {times3}) ──")
    t1, t2, t3 = times3
    for (b1, b2, b3) in [(0, 1, 2), (0, 2, 4), (0, 1, 3), (1, 3, 5)]:
        k3 = m.wall_cumulant([(b1, t1), (b2, t2), (b3, t3)])
        print(f"   bonds ({b1},{b2},{b3}): kappa_3 = {k3:+.6e}")
    print("── connected 4-time interface correlators  "
          f"<eta eta eta eta>_c  (times {times4}) ──")
    t1, t2, t3, t4 = times4
    for (b1, b2, b3, b4) in [(0, 1, 2, 3), (0, 2, 4, 6), (0, 1, 2, 4)]:
        k4 = m.wall_cumulant([(b1, t1), (b2, t2), (b3, t3), (b4, t4)])
        print(f"   bonds ({b1},{b2},{b3},{b4}): kappa_4 = {k4:+.6e}")
    print()


def run_all(L_RING, ETA, stage="all"):
    print(f"Glauber-Ising ring l={L_RING}, eta(=zeta)={ETA}, "
          f"gamma={gamma_from_eta(ETA):.6f}, p_wall={(1-ETA)/2:.4f}\n")
    if stage in ("all", "pf"):
        _pfaffian_self_test()
    m = Model(L_RING, ETA)
    if stage in ("all", "basic"):
        self_test(m)
        determinantal_3pt_test(m, times=(0.0, 0.5, 1.1))
        show_multitime_correlators(m, times3=(0.0, 0.4, 0.9),
                                   times4=(0.0, 0.3, 0.7, 1.2))
    if stage in ("all", "single"):
        single_time_determinantal_test(l=min(L_RING + 2, 14), t=0.20)
        single_time_determinantal_test(l=min(L_RING + 2, 14), t=0.45)
    if stage in ("all", "fit_single"):
        # decisive over-determined fit on the single-time quench walls
        st, _ = quench_state(min(L_RING + 2, 14), 0.30)
        n = 7
        bonds = list(range(n))
        structure_fit_report("single-time quench walls",
                             lambda S: st.moment([bonds[i] for i in S]), n)
    if stage in ("all", "fit_multi"):
        # decisive over-determined fit on the EQUILIBRIUM multi-time walls
        # 7 points, 7 distinct times, walking bonds -> genuine 7-time correlator
        n = 7
        ts = [0.0, 0.25, 0.55, 0.9, 1.3, 1.8, 2.4][:n]
        bds = [0, 1, 2, 3, 4, 5, 6][:n]
        pts = list(zip(bds, ts))
        structure_fit_report("equilibrium 7-time walls",
                             lambda S: m.wall_moment([pts[i] for i in S]), n)


if __name__ == "__main__":
    L_RING = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    ETA = float(sys.argv[2]) if len(sys.argv) > 2 else 0.45
    STAGE = sys.argv[3] if len(sys.argv) > 3 else "all"
    run_all(L_RING, ETA, STAGE)
