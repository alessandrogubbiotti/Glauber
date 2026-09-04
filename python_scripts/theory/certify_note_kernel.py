"""Numerical certificates for note_collaborators.tex

A. Theorem 2.5 (Garrod Thm 7, eq:G7) two-time 2-point intensity vs exact
   stationary Glauber chain on a ring (exact diagonalisation / expm).
B. Same at three times (6x6 Pfaffian).
C. Pfaffian 'collapse' claims: static discrete table -> theta^n, continuum
   equal-time kernel -> (kappa/2)^n, for n = 3,4,5.
D. Lemma 5.6 (extended entries): gauged discrete entries at large N vs the
   erfc closed forms of Psi^{+-}; error should be O(1/N).
"""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import ive, erfc
from math import atanh, tanh, sqrt, pi, exp

np.set_printoptions(precision=10, linewidth=140)

# ---------------------------------------------------------------- utilities
def pfaffian(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    if n == 0:
        return 1.0
    if n % 2:
        return 0.0
    tot = 0.0
    for j in range(1, n):
        if A[0, j] == 0.0:
            continue
        idx = [k for k in range(n) if k not in (0, j)]
        sub = A[np.ix_(idx, idx)]
        tot += ((-1) ** (j + 1)) * A[0, j] * pfaffian(sub)
    return tot

def assemble(points, block, diag12):
    """points: list of labels; block(i,j) -> 2x2 K(p_i,p_j) for i<j;
    diag12(i) -> K_12(p_i,p_i). Returns 2n x 2n antisymmetric matrix."""
    n = len(points)
    M = np.zeros((2 * n, 2 * n))
    for i in range(n):
        d = diag12(points[i])
        M[2 * i, 2 * i + 1] = d
        M[2 * i + 1, 2 * i] = -d
        for j in range(i + 1, n):
            B = block(points[i], points[j])
            M[2 * i:2 * i + 2, 2 * j:2 * j + 2] = B
            M[2 * j:2 * j + 2, 2 * i:2 * i + 2] = -B.T
    return M

# ---------------------------------------------------------------- discrete model
class Discrete:
    def __init__(self, zeta):
        self.zeta = zeta
        self.beta = atanh(zeta)
        self.gamma = tanh(2 * self.beta)
        self.theta = (1 - zeta) / 2
        self.ZMAX = 80

    def Rt(self, x, z):  # odd extension of zeta^{|z-x|}, \tilde R(x,x)=0
        d = z - x
        return np.sign(d) * self.zeta ** abs(d)

    def G(self, u, n):  # Green's function of L^A at the Glauber rates
        return ive(n, self.gamma * u) * np.exp(-(1 - self.gamma) * u)  # = e^{-u} I_n(gamma u)

    def conv(self, u, x, y):  # (p_u * \tilde R)(x,y) = sum_z G_u(y-z) R~(x,z)
        zs = np.arange(x - self.ZMAX, x + self.ZMAX + 1)
        return np.sum(self.G(u, y - zs) * self.Rt(x, zs))

    def K_static(self, x, y):  # eq:G2kernel with R = zeta^{y-x}, x<y
        n = y - x
        z = self.zeta
        return np.array([[-0.5 * z ** n, 0.5 * (1 - z) * z ** n],
                         [-0.5 * (1 - z) * z ** (n - 1), 0.5 * (1 - z) ** 2 * z ** (n - 1)]])

    def K_ext(self, u, x, y):  # eq:G7, s<t, u=t-s, any x,y
        c = lambda xx, yy: self.conv(u, xx, yy)
        p = lambda xx, yy: self.G(u, yy - xx)
        C = c(x, y)
        D2C = c(x, y + 1) - C
        D1C = c(x + 1, y) - C
        D12C = c(x + 1, y + 1) - c(x + 1, y) - c(x, y + 1) + C
        K11 = -0.5 * (C - p(x, y))
        K12 = -0.5 * (D2C - p(x, y + 1) + p(x, y))
        K21 = -0.5 * (D1C + p(x, y - 1) + p(x, y))
        K22 = -0.5 * (D12C + p(x, y + 1) - p(x, y - 1))
        return np.array([[K11, K12], [K21, K22]])

    def K_block(self, s, x, t, y):
        """general block for space-time points, using antisymmetry for s>t"""
        if s < t:
            return self.K_ext(t - s, x, y)
        if s > t:
            return -self.K_ext(s - t, y, x).T
        # equal times
        if x < y:
            return self.K_static(x, y)
        if x > y:
            return -self.K_static(y, x).T
        return np.array([[0.0, self.theta], [-self.theta, 0.0]])

    def pf_intensity(self, pts):  # pts = list of (t,x)
        M = assemble(pts, lambda a, b: self.K_block(a[0], a[1], b[0], b[1]),
                     lambda a: self.theta)
        return pfaffian(M)

# ---------------------------------------------------------------- exact ring
class Ring:
    def __init__(self, L, beta):
        self.L, self.beta = L, beta
        gamma = tanh(2 * beta)
        S = 1 << L
        states = np.arange(S, dtype=np.int64)
        bits = ((states[:, None] >> np.arange(L)) & 1).astype(np.int8)
        sig = 1 - 2 * bits  # spin +-1
        self.sig = sig
        # Gibbs weights
        E = np.sum(sig * np.roll(sig, -1, axis=1), axis=1)
        w = np.exp(beta * E)
        self.pi = w / w.sum()
        rows, cols, vals = [], [], []
        diag = np.zeros(S)
        for x in range(L):
            nb = sig[:, (x - 1) % L] + sig[:, (x + 1) % L]
            rate = 0.5 * (1 - 0.5 * gamma * sig[:, x] * nb)
            tgt = states ^ (1 << x)
            rows.append(states); cols.append(tgt); vals.append(rate)
            diag -= rate
        rows.append(states); cols.append(states); vals.append(diag)
        self.Q = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(S, S))

    def eta(self, x):  # wall on bond (x,x+1)
        return 0.5 * (1 - self.sig[:, x % self.L] * self.sig[:, (x + 1) % self.L])

    def evolve(self, u, h):
        return expm_multiply(self.Q * u, h)

    def two_time(self, a, u, b):
        return float(np.dot(self.pi * self.eta(a), self.evolve(u, self.eta(b))))

    def three_time(self, a, t1, b, t2, c):
        inner = self.eta(b) * self.evolve(t2 - t1, self.eta(c))
        return float(np.dot(self.pi * self.eta(a), self.evolve(t1, inner)))

# ================================================================ A + B
zeta = 0.35
D = Discrete(zeta)
print(f"zeta={zeta}, beta={D.beta:.6f}, gamma={D.gamma:.6f}, theta={D.theta}")
L = 16
R = Ring(L, D.beta)
print(f"ring L={L}; stationary check |pi Q|_inf = {np.abs(R.pi @ R.Q).max():.2e}")
print(f"E[eta] ring = {float(np.dot(R.pi, R.eta(0))):.10f}  vs theta = {D.theta:.10f}  (ring correction ~ zeta^L = {zeta**L:.1e})")

print("\n=== A. two-time intensities E[eta_0(0) eta_n(u)]: exact ring vs Pf(eq:G7) ===")
print(f"{'u':>5} {'n':>3} {'exact':>16} {'Pf(G7)':>16} {'diff':>11} {'theta^2':>10}")
for u in (0.5, 2.0, 5.0):
    for n in (-3, -1, 0, 1, 2, 4):
        ex = R.two_time(0, u, n)
        pf = D.pf_intensity([(0.0, 0), (u, n)])
        print(f"{u:5.1f} {n:3d} {ex:16.10f} {pf:16.10f} {ex-pf:11.2e} {D.theta**2:10.6f}")

print("\n=== B. three-time intensities E[eta_a(0) eta_b(t1) eta_c(t2)]: exact ring vs 6x6 Pfaffian ===")
print(f"{'(a,t1,b,t2,c)':>22} {'exact':>16} {'Pf':>16} {'diff':>11}")
for (a, t1, b, t2, c) in [(0, 1.0, 1, 2.5, -1), (0, 0.7, 3, 1.2, 2), (0, 2.0, 0, 4.0, 0),
                          (0, 1.0, -2, 3.0, 3), (0, 0.4, 1, 0.9, 1)]:
    ex = R.three_time(a, t1, b, t2, c)
    pf = D.pf_intensity([(0.0, a), (t1, b), (t2, c)])
    print(f"{str((a,t1,b,t2,c)):>22} {ex:16.10f} {pf:16.10f} {ex-pf:11.2e}")

# ================================================================ C
print("\n=== C. Pfaffian collapse: static discrete table, random points -> theta^n ? ===")
rng = np.random.default_rng(1)
for n in (3, 4, 5):
    xs = sorted(rng.choice(np.arange(-8, 9), size=n, replace=False))
    pf = D.pf_intensity([(0.0, int(x)) for x in xs])
    print(f"n={n} points={xs}: Pf={pf:.12f}, theta^n={D.theta**n:.12f}, diff={pf-D.theta**n:.1e}")

print("\n--- continuum equal-time kernel (eq:Kcmulti at v=0) -> (kappa/2)^n ? ---")
kap = 2.0
def Kcont0(x, y):  # x<y
    w = y - x
    e = np.exp(-kap * w)
    return np.array([[-0.5 * e, 0.5 * kap * e], [-0.5 * kap * e, 0.5 * kap ** 2 * e]])
for n in (3, 4, 5):
    xs = sorted(rng.uniform(-2, 2, size=n))
    M = assemble(xs, lambda a, b: Kcont0(a, b), lambda a: kap / 2)
    pf = pfaffian(M)
    print(f"n={n}: Pf={pf:.12f}, (kappa/2)^n={(kap/2)**n:.12f}, diff={pf-(kap/2)**n:.1e}")

# ================================================================ D
print("\n=== D. Lemma 5.6: gauged discrete extended entries vs continuum closed forms ===")
chat = kap ** 2 / 2
def Psi_plus(v, w):
    s = sqrt(2 * v)
    return 0.5 * (np.exp(-kap * w) * erfc((kap * v - w) / s) + np.exp(kap * w) * erfc((kap * v + w) / s))
def Psi_minus(v, w):
    s = sqrt(2 * v)
    return 0.5 * (np.exp(-kap * w) * erfc((kap * v - w) / s) - np.exp(kap * w) * erfc((kap * v + w) / s))
def gch(v, w):
    return exp(-chat * v) * exp(-w * w / (2 * v)) / sqrt(2 * pi * v)
def Kcont(v, w):
    return np.array([[-0.5 * Psi_minus(v, w), 0.5 * kap * Psi_plus(v, w) - gch(v, w)],
                     [-0.5 * kap * Psi_plus(v, w), 0.5 * kap ** 2 * Psi_minus(v, w)]])

# Wronskian / consistency checks of the closed forms
v0, w0, h = 0.7, 0.4, 1e-5
dPp = (Psi_plus(v0, w0 + h) - Psi_plus(v0, w0 - h)) / (2 * h)
dPm = (Psi_minus(v0, w0 + h) - Psi_minus(v0, w0 - h)) / (2 * h)
print(f"Wronskian: dPsi+ + kap Psi- = {dPp + kap*Psi_minus(v0,w0):.2e};  dPsi- - (2g - kap Psi+) = {dPm - (2*gch(v0,w0) - kap*Psi_plus(v0,w0)):.2e}")
print(f"Psi+_v(0) = {Psi_plus(v0,0.0):.10f} vs erfc(kap sqrt(v/2)) = {erfc(kap*sqrt(v0/2)):.10f}")
# direct numerical convolution check of Psi+
from scipy.integrate import quad
num = quad(lambda z: gch(v0, w0 - z) * np.exp(-kap * abs(z)), -30, 30, limit=400)[0]
print(f"Psi+ closed form {Psi_plus(v0,w0):.10f} vs numerical convolution {num:.10f}")

v = 0.5
for N in (50, 200, 800):
    zN = (2 * N - kap) / (2 * N + kap)
    DN = Discrete(zN)
    DN.ZMAX = int(12 * N)
    u = N ** 2 * v
    print(f"\nN={N}: zeta_N={zN:.6f}, gamma_N={DN.gamma:.8f}, N^2(1-gamma_N)={N**2*(1-DN.gamma):.6f} (chat={chat})")
    print(f"{'w':>6} {'entry':>6} {'gauged discrete':>18} {'continuum':>16} {'diff':>10}")
    for w in (-1.0, -0.3, 0.0, 0.3, 1.0):
        n = int(round(N * w))
        Kd = DN.K_ext(u, 0, n)
        Kg = np.array([[Kd[0, 0], N * Kd[0, 1]], [N * Kd[1, 0], N ** 2 * Kd[1, 1]]])
        Kc = Kcont(v, w)
        for (i, j, name) in ((0, 0, "11"), (0, 1, "12"), (1, 0, "21"), (1, 1, "22")):
            print(f"{w:6.2f} {name:>6} {Kg[i,j]:18.10f} {Kc[i,j]:16.10f} {Kg[i,j]-Kc[i,j]:10.2e}")

# also: N^2 * discrete 2-point intensity vs continuum Pfaffian (the actual claim of Thm 5.1)
print("\n--- N^2 rho^(2)_N(0,0; v, w) vs continuum Pfaffian ---")
for N in (50, 200, 800):
    zN = (2 * N - kap) / (2 * N + kap)
    DN = Discrete(zN); DN.ZMAX = int(12 * N)
    for w in (0.0, 0.3, 1.0):
        n = int(round(N * w))
        pfN = N ** 2 * DN.pf_intensity([(0.0, 0), (N ** 2 * v, n)])
        Kc = Kcont(v, w)
        pfc = (kap / 2) ** 2 - Kc[0, 0] * Kc[1, 1] + Kc[0, 1] * Kc[1, 0]
        print(f"N={N:4d} w={w:4.1f}: N^2 rho_N = {pfN:.8f}, continuum = {pfc:.8f}, diff = {pfN-pfc:.2e}")
