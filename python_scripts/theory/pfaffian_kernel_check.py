"""pfaffian_kernel_check.py — the checks of pfaffian_kernel.py.

  [1] the Pfaffian: Pf(A)^2 = det(A) and Pf(B^T A B) = det(B) Pf(A), a zero pivot;
  [2] the limit: k solves its equation with its data, the stationary case is the
      Poisson process, the value at contact of rho_2 is nu^2;
  [3] the stationary case against the kernels that were already in the
      repository: the correlation functions at two times against
      theory_continuum.f and at three times against cont3_final.fk (the complex
      kernel of the paper); the blocks, in the limit and on the lattice, against
      the ones of certify_note_kernel.py;
  [4] the lattice: the equation of the two-point function, the start from
      independent spins against theory_glauber.C_single_disordered, two times in
      the stationary case against theory_glauber.H_interface;
  [5] the lattice against the master equation: the generator of the Glauber
      dynamics on a ring of 14 sites is written explicitly, and the correlation
      functions of the interfaces after the quench, at one, two and three times,
      are compared with the Pfaffians;
  [6] the Glauber scaling: the rescaled blocks of the lattice and N^3 times the
      three-point functions converge to the ones of the limit, with an error of
      order 1/N.

Run:  python3 python_scripts/theory/pfaffian_kernel_check.py     (exit code 0 iff all pass)
"""
import ast
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))                    # theory/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))   # python_scripts/

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

import pfaffian_kernel as pk
import theory_continuum as tc
import theory_glauber as tg
from cont3_final import Kcont, fk
from multitime_interface_dpp import pfaffian as pfaffian_dpp

_results = []


def check(name, value, bound, detail=""):
    ok = bool(abs(value) <= bound)
    _results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value:.1e} (bound {bound:.0e})"
          + (f"\n         {detail}" if detail else ""))


def accepted(call):
    """0 if the call raises ValueError, as it has to, and 1 if it is accepted."""
    try:
        call()
    except ValueError:
        return 0.0
    return 1.0


# ---------------------------------------------------------------------------
# The kernel of certify_note_kernel.py
# ---------------------------------------------------------------------------

def note_kernel():
    """The definitions of certify_note_kernel.py, as a dictionary.

    That script runs its certificates when it is imported.  Here only its imports,
    its functions and classes and the two constants kap and chat are executed.
    Its kappa is 2 nu and its times are the Glauber ones, twice the ones of the engine.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "certify_note_kernel.py")
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read())
    keep = [node for node in tree.body
            if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.ClassDef))
            or (isinstance(node, ast.Assign)
                and any(isinstance(t, ast.Name) and t.id in ("kap", "chat") for t in node.targets))]
    namespace = {}
    exec(compile(ast.Module(body=keep, type_ignores=[]), path, "exec"), namespace)
    return namespace


def note_kernel_deviations():
    """The largest differences between the blocks of the note and the ones of pfaffian_kernel.py."""
    old = note_kernel()
    nu = old["kap"] / 2.0
    limit = 0.0
    for v in (0.05, 0.5, 0.7, 2.0):
        for w in (-1.0, -0.3, 0.0, 0.3, 1.0, 2.5):
            limit = max(limit, np.max(np.abs(old["Kcont"](v, w) - pk.extended_block(w, 0.3, v / 2.0, nu))))
    for w in (0.01, 0.3, 1.0, 2.5):
        limit = max(limit, np.max(np.abs(old["Kcont0"](0.0, w) - pk.kernel_block(w, 0.4, nu))))
    lattice = pfaffians = 0.0
    for N in (20, 50):
        D = old["Discrete"]((2 * N - old["kap"]) / (2 * N + old["kap"]))
        D.ZMAX = 12 * N
        rescale = np.array([[1.0, N], [N, N * N]])
        for theta in (0.02, 0.25):
            for r in (-N, -3, -1, 0, 1, 2, N // 3, N):
                lattice = max(lattice, np.max(np.abs(D.K_ext(2.0 * N * N * theta, 0, r) * rescale
                                                     - pk.lattice_extended_block(r, 0.3, theta, N, nu))))
        for r in (1, 2, 7, N):
            lattice = max(lattice, np.max(np.abs(D.K_static(0, r) * rescale - pk.lattice_block(r, 0.3, N, nu))))
        for points in ([(0.0, 0), (0.7, 3), (1.2, 2)], [(0.0, 0), (1.0, 1), (2.5, -1)],
                       [(0.0, 0), (0.0, 2), (2.0, 1), (2.0, 4)]):
            a = D.pf_intensity(points)                 # the times of the note are the Glauber ones
            b = pk.lattice_correlation([(0.1 + u / (2.0 * N * N), x) for u, x in points], N, nu)
            pfaffians = max(pfaffians, abs(b / N ** len(points) / a - 1.0))
    return limit, lattice, pfaffians


# ---------------------------------------------------------------------------
# The master equation on a ring
# ---------------------------------------------------------------------------

class Ring:
    """The Glauber dynamics on a ring of l sites, all the 2^l configurations.

    The spin at x flips at rate 1 - (gamma/2) sigma_x (sigma_{x-1} + sigma_{x+1}).  The
    start is the measure with independent interfaces of probability lam on each bond.
    """

    def __init__(self, l, gamma, lam):
        states = np.arange(1 << l)
        self.spin = [2 * ((states >> x) & 1) - 1 for x in range(l)]
        self.eta = [(1 - self.spin[x] * self.spin[(x + 1) % l]) // 2 for x in range(l)]
        rows, cols, vals = [], [], []
        total = np.zeros(1 << l)
        for x in range(l):
            w = 1.0 - 0.5 * gamma * self.spin[x] * (self.spin[x - 1] + self.spin[(x + 1) % l])
            rows.append(states ^ (1 << x))           # the configuration after the flip
            cols.append(states)
            vals.append(w)
            total += w
        rows.append(states)
        cols.append(states)
        vals.append(-total)
        # the transpose of the generator: it acts on the measures
        self.Lt = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(1 << l, 1 << l)).tocsr()
        walls = sum(self.eta)
        mu = lam ** walls * (1.0 - lam) ** (l - walls)
        self.mu0 = mu / mu.sum()

    def correlation(self, points):
        """E[eta_{x_1}(t_1) ... eta_{x_n}(t_n)], points = [(t, x), ...], t the time of the lattice."""
        v, now = self.mu0, 0.0
        for t, x in sorted(points):
            if t > now:
                v = expm_multiply(self.Lt * (t - now), v)
                now = t
            v = v * self.eta[x]
        return float(v.sum())


def master_equation_deviation(N, nu, nu_prime, l=14):
    """The largest difference between the ring and the Pfaffians, over some sets of points."""
    gamma = (N * N - nu * nu) / (N * N + nu * nu)
    ring = Ring(l, gamma, nu_prime / (N + nu_prime))
    t1, t2, t3 = 0.4, 0.9, 1.5                       # times of the lattice
    cases = [
        [(t2, 3)],
        [(t1, 3), (t1, 4)], [(t2, 3), (t2, 5)], [(t3, 2), (t3, 5)],
        [(t2, 3), (t2, 4), (t2, 5)], [(t3, 2), (t3, 4), (t3, 5)], [(t1, 2), (t1, 3), (t1, 5), (t1, 6)],
        [(t1, 3), (t2, 3)], [(t1, 3), (t2, 4)], [(t1, 4), (t3, 2)], [(0.0, 3), (t2, 5)],
        [(t1, 3), (t2, 5), (t3, 4)], [(t1, 4), (t2, 4), (t3, 4)], [(t1, 5), (t1, 3), (t3, 4)],
        [(t1, 2), (t2, 5), (t2, 3), (t3, 4)],
    ]
    worst = 0.0
    for points in cases:
        exact = ring.correlation(points)
        scaled = pk.lattice_correlation([(t / (N * N), x) for t, x in points], N, nu, nu_prime)
        worst = max(worst, abs(scaled / N ** len(points) - exact))
    return worst, len(cases)


# ---------------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(0)
    print("=" * 74)
    print(" The checks of pfaffian_kernel.py")
    print("=" * 74)

    print("\n[1] the Pfaffian")
    worst = 0.0
    for n in (2, 4, 6, 8, 10):
        M = rng.normal(size=(n, n))
        M = M - M.T
        worst = max(worst, abs(pk.pfaffian(M) ** 2 - np.linalg.det(M)) / abs(np.linalg.det(M)))
    check("Pf(A)^2 = det(A), sizes 2, ..., 10", worst, 1e-10)
    B = rng.normal(size=(6, 6))
    M = rng.normal(size=(6, 6))
    M = M - M.T
    check("Pf(B^T A B) = det(B) Pf(A)",
          (pk.pfaffian(B.T @ M @ B) - np.linalg.det(B) * pk.pfaffian(M)) / pk.pfaffian(M), 1e-10)
    worst = 0.0
    for n in (2, 4, 6, 8, 10):
        M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        M = M - M.T
        worst = max(worst, abs(pk.pfaffian(M) / pfaffian_dpp(M) - 1.0))
    check("the same as multitime_interface_dpp.pfaffian, complex matrices", worst, 1e-10)
    M = np.array([[0.0, 0.0, 1.0, 2.0], [0.0, 0.0, 3.0, 4.0], [-1.0, -3.0, 0.0, 5.0], [-2.0, -4.0, -5.0, 0.0]])
    check("a zero in the place of the first pivot: Pf = a12 a34 - a13 a24 + a14 a23 = 2", pk.pfaffian(M) - 2.0, 1e-12)
    check("an odd size gives 0; refused: a matrix that is not antisymmetric, also a small one, or not square",
          abs(pk.pfaffian(np.zeros((3, 3)))) + accepted(lambda: pk.pfaffian(np.ones((4, 4))))
          + accepted(lambda: pk.pfaffian(4e-9 * np.ones((4, 4)))) + accepted(lambda: pk.pfaffian(np.zeros((2, 3)))), 0.5)
    M = rng.normal(size=(6, 6))
    M = 1e12 * (M - M.T)
    check("a large antisymmetric matrix is accepted: Pf(A)^2 = det(A)",
          pk.pfaffian(M) ** 2 / np.linalg.det(M) - 1.0, 1e-10)

    print("\n[2] the limit")
    pde = der = data = 0.0
    h = 1e-5
    d = np.array([0.1, 0.4, 1.0, 2.5])
    for nu, nup, tau in ((1.0, 3.0, 0.05), (1.0, 0.4, 0.3), (2.0, 0.7, 0.1), (0.5, 2.5, 1.0), (1.0, 100.0, 0.02)):
        k, k1, k2 = pk.k_limit(d, tau, nu, nup)
        dtau = (pk.k_limit(d, tau + h, nu, nup)[0] - pk.k_limit(d, tau - h, nu, nup)[0]) / (2 * h)
        pde = max(pde, np.max(np.abs(dtau - (2 * k2 - 8 * nu * nu * k))))
        up, down = pk.k_limit(d + h, tau, nu, nup), pk.k_limit(d - h, tau, nu, nup)
        der = max(der, np.max(np.abs((up[0] - down[0]) / (2 * h) - k1)),
                  np.max(np.abs((up[1] - down[1]) / (2 * h) - k2)))
        data = max(data, abs(pk.k_limit(0.0, tau, nu, nup)[0][0] - 1.0),
                   np.max(np.abs(pk.k_limit(d, 0.0, nu, nup)[0] - np.exp(-2 * nup * d))))
    check("d/dtau k = 2 k'' - 8 nu^2 k", pde, 1e-5)
    check("k' and k'' are the derivatives of k", der, 1e-5)
    check("k(0) = 1 and k = exp(-2 nu' d) at tau = 0", data, 1e-12)
    k, k1, k2 = pk.k_limit(d, 0.0, 1.0, 3.0)
    check("at tau = 0 also k' = -2 nu' k and k'' = 4 nu'^2 k",
          max(np.max(np.abs(k1 + 6.0 * k)), np.max(np.abs(k2 - 36.0 * k))), 1e-12)
    check("refused: a block at distance 0, two times that coincide, two equal points, a negative distance",
          accepted(lambda: pk.kernel_block(0.0, 0.1)) + accepted(lambda: pk.extended_block(0.3, 0.1, 0.0))
          + accepted(lambda: pk.kernel_matrix([(0.1, 0.2), (0.1, 0.2)])) + accepted(lambda: pk.k_limit(-0.3, 0.1)), 0.5)
    check("refused: a negative time and a time that is not a number, in every function of the limit",
          accepted(lambda: pk.k_limit(0.3, -0.1)) + accepted(lambda: pk.k_limit(0.3, np.nan))
          + accepted(lambda: pk.density(-0.1)) + accepted(lambda: pk.density(np.nan))
          + accepted(lambda: pk.kernel_block(0.3, -0.1)) + accepted(lambda: pk.kernel_block(np.nan, 0.1))
          + accepted(lambda: pk.extended_block(0.3, -0.1, 0.2)) + accepted(lambda: pk.extended_block(0.3, 0.1, np.nan)), 0.5)
    B = pk.kernel_block(0.3, 0.1, 1.0, 3.0)
    B[0, 0] = 99.0
    check("the block that is returned can be modified without changing the next ones",
          pk.kernel_block(0.3, 0.1, 1.0, 3.0)[0, 0] + 0.5 * pk.k_limit(0.3, 0.1, 1.0, 3.0)[0][0], 1e-15)
    points = [(0.3, 0.1), (0.3, 0.45), (0.3, 1.2), (0.3, 1.3)]
    check("nu' = nu, one time: rho_4 = nu^4", pk.correlation(points, 1.7) - 1.7 ** 4, 1e-10)
    points = [(0.3, 0.1), (0.5, 0.45), (0.9, -0.2)]
    check("nu' = nu: rho_1 = nu at every time", pk.kernel_matrix(points, 1.7)[4, 5] - 1.7, 1e-12)
    check("the density is nu' at tau = 0", pk.density(0.0, 1.0, 3.0) - 3.0, 1e-12)
    check("the density tends to nu", pk.density(50.0, 1.0, 3.0) - 1.0, 1e-12)
    check("at contact rho_2(0+) = nu^2", pk.pair_density(1e-9, 0.2, 1.3, 4.0)[0] - 1.3 ** 2, 1e-6)
    worst = 0.0
    for nu, nup in ((1.0, 3.0), (1.5, 0.5)):
        worst = max(worst, abs(pk.correlation([(0.2, 0.0), (0.2 + 1e-9, 0.6)], nu, nup)
                               - pk.correlation([(0.2, 0.0), (0.2, 0.6)], nu, nup)))
    check("two times: as the times get close, the correlation at one time", worst, 1e-7)
    a = pk.correlation([(0.05, 0.0), (0.2, 0.3), (0.2, 0.7), (0.4, 0.1)], 1.0, 3.0)
    b = pk.correlation([(0.4, 0.1), (0.2, 0.7), (0.05, 0.0), (0.2, 0.3)], 1.0, 3.0)
    check("the correlation does not depend on the order of the points", a - b, 1e-12)

    print("\n[3] the stationary case against the kernels that were already in the repository (nu = 1)")
    worst = 0.0
    for x, theta in ((0.3, 0.1), (-0.5, 0.2), (0.2, 0.05), (1.5, 0.7)):
        worst = max(worst, abs(pk.correlation([(0.4, 0.0), (0.4 + theta, x)], 1.0) - 1.0 - tc.f(x, theta)[0]))
    check("two times: rho - nu^2 = theory_continuum.f", worst, 1e-10)
    worst = 0.0
    for xs, taus in (([-0.2, 0.0, 0.25], [0.05, 0.10, 0.18]), ([0.0, 0.15, 0.35], [0.04, 0.09, 0.15]),
                     ([0.1, -0.15, 0.3], [0.07, 0.12, 0.05]), ([-0.1, 0.05, 0.22], [0.03, 0.11, 0.20])):
        P = list(zip(taus, xs))
        pairs = sum(pk.correlation([P[i], P[j]], 1.0) for i, j in ((0, 1), (0, 2), (1, 2)))
        centered = pk.correlation(P, 1.0) - pairs + 2.0          # E prod (eta_i - nu), nu = 1
        worst = max(worst, abs(centered - fk(list(zip(xs, taus)), Kcont).real))
    check("three times: the centered function = cont3_final.fk", worst, 1e-10)
    limit, lattice, pfaffians = note_kernel_deviations()
    check("certify_note_kernel.py: the blocks of the limit, at one and at two times", limit, 1e-12)
    check("certify_note_kernel.py: the blocks of the lattice, at one and at two times", lattice, 1e-10)
    check("certify_note_kernel.py: its Pfaffians on the lattice, three and four points (relative)", pfaffians, 1e-10)

    print("\n[4] the lattice")
    N, nu, nup, tau = 40, 1.0, 3.0, 0.07
    gamma = (N * N - nu * nu) / (N * N + nu * nu)
    r = np.arange(1, 60)
    h = 1e-6
    k = pk.k_lattice(np.arange(0, 62), tau, N, nu, nup)
    dt = (pk.k_lattice(r, tau + h, N, nu, nup) - pk.k_lattice(r, tau - h, N, nu, nup)) / (2 * h) / (N * N)
    rhs = 2.0 * gamma * (k[r + 1] + k[r - 1] - 2 * k[r]) - 4.0 * (1 - gamma) * k[r]
    check("dk/dt = 2 gamma (k(r+1) + k(r-1) - 2 k(r)) - 4 (1 - gamma) k", np.max(np.abs(dt - rhs)), 1e-8)
    check("k(0) = 1 and the datum at t = 0",
          max(abs(k[0] - 1.0), np.max(np.abs(pk.k_lattice(r, 0.0, N, nu, nup) - ((N - nup) / (N + nup)) ** r))),
          1e-13)
    worst = 0.0
    for N, nu, tau in ((20, 1.0, 0.05), (50, 2.0, 0.3), (30, 0.5, 1.0)):
        zeta, gamma = (N - nu) / (N + nu), (N * N - nu * nu) / (N * N + nu * nu)
        r = np.arange(0, 80)
        ref = tg.C_single_disordered(r, 2 * N * N * tau, zeta, gamma)   # the Glauber time is 2 N^2 tau
        worst = max(worst, np.max(np.abs(pk.k_lattice(r, tau, N, nu, N) - ref)))
    check("nu' = N, independent spins: theory_glauber.C_single_disordered", worst, 1e-11)
    worst = 0.0
    for N, nu, theta in ((20, 1.0, 0.05), (50, 2.0, 0.3), (30, 0.5, 1.0)):
        zeta, gamma = (N - nu) / (N + nu), (N * N - nu * nu) / (N * N + nu * nu)
        for j in (-7, -1, 0, 1, 2, 15):
            ref = tg.H_interface(j, 2 * N * N * theta, zeta, gamma)
            mine = pk.lattice_correlation([(0.3, 5), (0.3 + theta, 5 + j)], N, nu) / (N * N)
            worst = max(worst, abs(mine - ref) / ref)
    check("nu' = nu, two times: theory_glauber.H_interface (relative)", worst, 1e-9)
    N, nu, nup = 30, 1.0, 3.0
    P = [(0.05, 0), (0.2, 9), (0.2, 21), (0.4, 3)]
    a = pk.lattice_correlation(P, N, nu, nup)
    worst = max(abs(pk.lattice_correlation([P[i] for i in order], N, nu, nup) / a - 1.0)
                for order in ((3, 2, 0, 1), (1, 0, 3, 2), (2, 3, 1, 0), (3, 0, 2, 1)))
    check("the correlation on the lattice does not depend on the order of the points", worst, 1e-12)
    k = pk.k_lattice(np.arange(0, 3000), 0.07, N, nu, nup)             # one table for all the distances
    G, p = pk._propagated_table(0.07, 0.05, N, nu, nup, 4096)
    rescale = np.array([[1.0, N], [N, N * N]])
    worst = 0.0
    for r in (254, 255, 256, 257, 510, 511, 512, 1022, 1023, 1024):    # around the sizes of the tables
        B = -0.5 * np.array([[k[r], k[r + 1] - k[r]], [k[r - 1] - k[r], 2 * k[r] - k[r - 1] - k[r + 1]]]) * rescale
        worst = max(worst, np.max(np.abs(pk.lattice_block(r, 0.07, N, nu, nup) - B)))
        for sign in (1, -1):
            g = [sign * G[r + m] for m in (-1, 0, 1)] if sign > 0 else [-G[r - m] for m in (-1, 0, 1)]
            q = [p[r + sign * m] for m in (-1, 0, 1)]                  # G and p at sign r - 1, sign r, sign r + 1
            B = -0.5 * np.array([[g[1] - q[1], g[2] - g[1] - q[2] + q[1]],
                                 [g[0] - g[1] + q[0] + q[1], 2 * g[1] - g[0] - g[2] + q[2] - q[0]]]) * rescale
            worst = max(worst, np.max(np.abs(pk.lattice_extended_block(sign * r, 0.07, 0.05, N, nu, nup) - B)))
    check("the blocks at the distances where a table ends are the ones of a longer table", worst, 1e-15)
    check("refused: a block of the lattice at distance 0, two equal points, a negative distance",
          accepted(lambda: pk.lattice_block(0, 0.1, N)) + accepted(lambda: pk.lattice_kernel_matrix([(0.1, 2), (0.1, 2)], N))
          + accepted(lambda: pk.k_lattice(-1, 0.1, N)), 0.5)
    check("refused: a negative time on the lattice, also in the stationary case, and two times that coincide",
          accepted(lambda: pk.k_lattice(1, -0.1, N)) + accepted(lambda: pk.lattice_block(3, -0.1, N, 1.0, 3.0))
          + accepted(lambda: pk.lattice_block(3, -0.1, N, 1.0)) + accepted(lambda: pk.lattice_density(np.nan, N))
          + accepted(lambda: pk.lattice_extended_block(3, -0.1, 0.2, N, 1.0))
          + accepted(lambda: pk.lattice_extended_block(1, 0.1, 0.0, N)), 0.5)

    print("\n[5] the lattice against the master equation on a ring of 14 sites")
    for N, nu, nup in ((4, 3.0, 6.0), (5, 3.5, 4.4)):
        worst, n = master_equation_deviation(N, nu, nup)
        check(f"N = {N}, nu = {nu}, nu' = {nup}: {n} correlation functions, up to four points and three times",
              worst, 1e-6)

    print("\n[6] the Glauber scaling")
    nu, nup = 1.0, 3.0
    Ns = (100, 200, 400)
    one, two, three, three_times = [], [], [], []
    for N in Ns:
        one.append(np.max(np.abs(pk.lattice_block(round(0.3 * N), 0.07, N, nu, nup)
                                 - pk.kernel_block(0.3, 0.07, nu, nup))))
        two.append(max(np.max(np.abs(pk.lattice_extended_block(round(d * N), 0.07, 0.05, N, nu, nup)
                                     - pk.extended_block(d, 0.07, 0.05, nu, nup))) for d in (-0.4, 0.0, 0.3)))
        P = [(0.07, 0.0), (0.07, 0.3), (0.07, 0.8)]
        three.append(abs(pk.lattice_correlation([(t, round(N * u)) for t, u in P], N, nu, nup)
                         / pk.correlation(P, nu, nup) - 1.0))
        P = [(0.03, 0.0), (0.07, 0.3), (0.12, -0.2)]
        three_times.append(abs(pk.lattice_correlation([(t, round(N * u)) for t, u in P], N, nu, nup)
                               / pk.correlation(P, nu, nup) - 1.0))
    for name, err in (("the block at one time", one), ("the block at two times", two),
                      ("N^3 rho_3^N / rho_3 - 1, one time", three),
                      ("N^3 rho_3^N / rho_3 - 1, three times", three_times)):
        detail = "N = 100, 200, 400: " + ", ".join(f"{e:.1e}" for e in err)
        halves = all(1.6 < err[i] / err[i + 1] < 2.4 for i in range(2))
        check(f"{name}: the error at N = 400", err[-1] if halves else np.inf, 2e-2, detail)

    print("\n" + "=" * 74)
    print(f" {sum(_results)}/{len(_results)} checks passed.")
    print("=" * 74)
    return 0 if all(_results) else 1


if __name__ == "__main__":
    sys.exit(main())
