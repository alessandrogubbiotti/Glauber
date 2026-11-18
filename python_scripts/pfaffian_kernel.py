"""pfaffian_kernel.py — the Pfaffian kernel of the interfaces, on the lattice and in the limit.

The interfaces (walls) of the Glauber dynamics are at every time a Pfaffian point
process (Garrod): their correlation functions are Pfaffians of a 2x2 block kernel
built from the spin two-point function K(x,y) = <sigma_x sigma_y>.  Here the
interfaces start independent with intensity nu' and the dynamics has the rates of
intensity nu (hop 1 per direction, annihilation 1+gamma, creation 1-gamma,
gamma = (N^2-nu^2)/(N^2+nu^2)).  The case nu' = nu is the stationary one, and
nu' = N is the start from independent spins.

In the limit (mesh 1/N, time tau = t/N^2) the two-point function is K(u,v) = k(v-u),
where k is the bounded solution of

    d/dtau k = 2 k'' - 8 nu^2 k   for d > 0,    k(0) = 1,    k = exp(-2 nu' d) at tau = 0,

and the interfaces are the Pfaffian point process with the kernel, for u < v,

    K(u,v) = -1/2 [[ k,   k'  ],
                   [-k', -k'' ]] (v-u),       K_12(u,u) = -k'(0+)/2 = rho_1 .

For two times sigma < tau the kernel is the one at the time sigma, with the second
argument propagated by the heat kernel killed at rate 4 nu^2 (extended kernel).
The correlation functions are

    rho(tau_1,u_1; ...; tau_n,u_n) = Pf [ K(tau_i,u_i; tau_j,u_j) ]_{i,j<=n} ,

and on the lattice N^n times the correlation function converges to it.

Functions (nu_prime is the initial intensity nu'; without it the start is stationary):

    k_limit(d, tau, nu, nu_prime)             k, k', k'' in the limit
    density(tau, nu, nu_prime)                rho_1
    pair_density(d, tau, nu, nu_prime)        rho_2 at one time
    kernel_block(d, tau, ...)                 the 2x2 block at one time
    extended_block(d, sigma, theta, ...)      the 2x2 block between the times sigma, sigma+theta
    kernel_matrix(points, nu, nu_prime)       the antisymmetric matrix of the blocks
    correlation(points, nu, nu_prime)         its Pfaffian; points = [(tau, u), ...]

    k_lattice(r, tau, N, nu, nu_prime)        <sigma_0 sigma_r> on the lattice, exact
    lattice_density(tau, N, ...)              N E[eta_x]
    lattice_block, lattice_extended_block     the blocks on the lattice, rescaled
    lattice_kernel_matrix(points, N, ...)     the rescaled matrix; points = [(tau, x), ...]
    lattice_correlation(points, N, ...)       N^n E[eta_x1(t_1) ... eta_xn(t_n)], exact

    pfaffian(A)                               the Pfaffian of an antisymmetric matrix

The blocks on the lattice are rescaled as [[K_11, N K_12], [N K_21, N^2 K_22]], so that
they converge entry by entry to the ones of the limit.  The time is always the
macroscopic one, tau = t/N^2, which is the time of the programs, and x is the
integer index of the bond: the interface eta_x is between the sites x and x+1.

Usage:
    python3 python_scripts/pfaffian_kernel.py 0.1,0 0.1,0.3 0.25,0.1 --nu 1 --nu-prime 3
    python3 python_scripts/pfaffian_kernel.py 0.1,0 0.1,0.3 --nu 1 --nu-prime 3 --N 100

Each point is tau,u.  The first prints the matrix of the kernel at the points and the
correlation function of the limit.  With --N it prints also the rescaled ones of the
lattice of mesh 1/N, at the bonds floor(N u).  The checks of this module are in
theory/pfaffian_kernel_check.py.
"""
import argparse
from functools import lru_cache

import numpy as np
from scipy.special import erfcx, ive

__all__ = ["k_limit", "density", "pair_density", "kernel_block", "extended_block",
           "kernel_matrix", "correlation", "k_lattice", "lattice_density",
           "lattice_block", "lattice_extended_block", "lattice_kernel_matrix",
           "lattice_correlation", "pfaffian"]


# ---------------------------------------------------------------------------
# Three functions of the Gaussian.  For Z normal with variance 2 s,
#
#     odd_part(a, d, s)  = E[ sgn(d+Z) exp(-a |d+Z|) ]     (Phi_a(d,s) of the notes)
#     even_part(a, d, s) = E[ exp(-a |d+Z|) ]
#     heat(d, s)         = the density of Z at d
#
# and  d/dd odd = -a even + 2 heat,   d/dd even = -a odd.  All of them are for
# s > 0 and are multiplied by exp(-c), with c given: a^2 s can be large while the
# product is not.
# ---------------------------------------------------------------------------

def _exp_erfc(a, d, s, c, sign):
    """exp(-c + a^2 s + sign a d) erfc((2 a s + sign d) / (2 sqrt(s))), without overflow."""
    z = (2.0 * a * s + sign * d) / (2.0 * np.sqrt(s))
    gauss = np.exp(-c - d * d / (4.0 * s))            # exp(-c + a^2 s + sign a d - z^2)
    out = gauss * erfcx(np.abs(z))
    neg = z < 0                                       # erfc(z) = 2 - erfc(-z)
    out[neg] = 2.0 * np.exp(-c + a * a * s + sign * a * d[neg]) - out[neg]
    return out


def odd_part(a, d, s, c=0.0):
    """exp(-c) E[sgn(d+Z) exp(-a|d+Z|)], Z normal with variance 2s, s > 0."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    return 0.5 * (_exp_erfc(a, d, s, c, -1.0) - _exp_erfc(a, d, s, c, +1.0))


def even_part(a, d, s, c=0.0):
    """exp(-c) E[exp(-a|d+Z|)], Z normal with variance 2s, s > 0."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    return 0.5 * (_exp_erfc(a, d, s, c, -1.0) + _exp_erfc(a, d, s, c, +1.0))


def heat(d, s, c=0.0):
    """exp(-c) times the density at d of the normal with variance 2s, s > 0."""
    d = np.atleast_1d(np.asarray(d, dtype=float))
    return np.exp(-c - d * d / (4.0 * s)) / np.sqrt(4.0 * np.pi * s)


def _key(tau, nu, nu_prime):
    """The arguments as they enter the tables: floats, and nu' = nu when it is not given."""
    if not tau >= 0:                                  # written so that a NaN is refused too
        raise ValueError("the time must be >= 0")
    return float(tau), float(nu), float(nu if nu_prime is None else nu_prime)


# ---------------------------------------------------------------------------
# The limit, one time
# ---------------------------------------------------------------------------

def k_limit(d, tau, nu=1.0, nu_prime=None):
    """The two-point function k(d) of the limit and its first two derivatives, for d >= 0.

    k(d) = exp(-2 nu d) + exp(-8 nu^2 tau) [Phi_{2nu'}(d, 2tau) - Phi_{2nu}(d, 2tau)].
    The three arrays k, k', k'' are returned; at d = 0 the derivatives are the ones
    from the right.
    """
    nu_prime = nu if nu_prime is None else nu_prime
    d = np.atleast_1d(np.asarray(d, dtype=float))
    if not np.all(d >= 0):
        raise ValueError("k_limit is defined for d >= 0")
    if not tau >= 0:
        raise ValueError("the time must be >= 0")
    a, b, s, c = 2.0 * nu_prime, 2.0 * nu, 2.0 * tau, 8.0 * nu * nu * tau
    if tau == 0:
        k = np.exp(-a * d)
        return k, -a * k, a * a * k
    # The stationary profile and the two terms of Phi_{2nu} are taken together,
    #     exp(-b d) - exp(-c) Phi_b = (um + up)/2,
    # with um = exp(-b d) erfc((d - 2bs)/(2 sqrt s)) and up = exp(b d) erfc((d + 2bs)/(2 sqrt s)):
    # a sum of positive terms, so that k keeps its relative precision where it is small.
    um = _exp_erfc(-b, d, s, c, +1.0)
    up = _exp_erfc(b, d, s, c, +1.0)
    k = 0.5 * (um + up) + odd_part(a, d, s, c)
    k1 = -0.5 * b * (um - up) - a * even_part(a, d, s, c)
    k2 = 0.5 * b * b * (um + up) + a * a * odd_part(a, d, s, c)
    return k, k1, k2


@lru_cache(maxsize=4096)
def _density(tau, nu, nu_prime):
    return float(-0.5 * k_limit(0.0, tau, nu, nu_prime)[1][0])


def density(tau, nu=1.0, nu_prime=None):
    """The density of the interfaces at the time tau, rho_1 = -k'(0+)/2.

    rho_1 = nu erf(2 nu sqrt(2 tau)) + nu' exp(-8 nu^2 tau) erfcx(2 nu' sqrt(2 tau)).
    """
    return _density(*_key(tau, nu, nu_prime))


def pair_density(d, tau, nu=1.0, nu_prime=None):
    """The two-point function of the interfaces at one time, at distance d > 0.

    rho_2 = rho_1^2 + (k k'' - k'^2)/4, the Pfaffian of the 4x4 matrix.
    """
    k, k1, k2 = k_limit(d, tau, nu, nu_prime)
    rho1 = density(tau, nu, nu_prime)
    return rho1 * rho1 + 0.25 * (k * k2 - k1 * k1)


@lru_cache(maxsize=4096)
def _kernel_block(d, tau, nu, nu_prime):
    k, k1, k2 = (x[0] for x in k_limit(d, tau, nu, nu_prime))
    return -0.5 * np.array([[k, k1], [-k1, -k2]])


def kernel_block(d, tau, nu=1.0, nu_prime=None):
    """The 2x2 block K(u,v) at one time, for d = v - u > 0."""
    if not d > 0:
        raise ValueError("kernel_block needs d > 0")
    return _kernel_block(float(d), *_key(tau, nu, nu_prime)).copy()


# ---------------------------------------------------------------------------
# The limit, two times
# ---------------------------------------------------------------------------

def extended_block(d, sigma, theta, nu=1.0, nu_prime=None):
    """The 2x2 block K(sigma,u; sigma+theta,v) of the extended kernel, d = v - u real, theta > 0.

    Let G(d) be the odd extension of the two-point function at the time sigma, propagated
    for a time theta by the heat kernel killed at rate 4 nu^2, and p this heat kernel:

        G = exp(-4 nu^2 theta) Phi_{2nu}(d, theta)
            + exp(-4 nu^2 (theta + 2 sigma)) [Phi_{2nu'} - Phi_{2nu}](d, theta + 2 sigma).

    Then G' = g1 + 2p and G'' = g2 + 2p', with g1 and g2 without the heat kernel, and

        K = [[ -G/2,   -g1/2 - p ],
             [  g1/2,   g2/2     ]].
    """
    nu_prime = nu if nu_prime is None else nu_prime
    if not theta > 0:
        raise ValueError("extended_block needs theta > 0; at equal times use kernel_block")
    if not sigma >= 0:
        raise ValueError("the time must be >= 0")
    a, b = 2.0 * nu_prime, 2.0 * nu
    c1 = 4.0 * nu * nu * theta                       # the killing between the two times
    c2 = c1 + 8.0 * nu * nu * sigma
    s2 = theta + 2.0 * sigma
    G = odd_part(b, d, theta, c1) + odd_part(a, d, s2, c2) - odd_part(b, d, s2, c2)
    g1 = -b * even_part(b, d, theta, c1) - a * even_part(a, d, s2, c2) + b * even_part(b, d, s2, c2)
    g2 = b * b * odd_part(b, d, theta, c1) + a * a * odd_part(a, d, s2, c2) - b * b * odd_part(b, d, s2, c2)
    p = heat(d, theta, c1)
    return np.array([[-0.5 * G[0], -0.5 * g1[0] - p[0]], [0.5 * g1[0], 0.5 * g2[0]]])


def _assemble(points, diagonal, block):
    """The antisymmetric matrix with the given 2x2 blocks; the rows 2i, 2i+1 are the point i."""
    n = len(points)
    A = np.zeros((2 * n, 2 * n))
    for i in range(n):
        rho = diagonal(points[i])
        A[2 * i, 2 * i + 1], A[2 * i + 1, 2 * i] = rho, -rho
        for j in range(i + 1, n):
            B = block(points[i], points[j])
            A[2 * i:2 * i + 2, 2 * j:2 * j + 2] = B
            A[2 * j:2 * j + 2, 2 * i:2 * i + 2] = -B.T
    return A


def kernel_matrix(points, nu=1.0, nu_prime=None):
    """The 2n x 2n antisymmetric matrix of the blocks at the space-time points (tau_i, u_i).

    The points must be pairwise distinct.  The block of two points is the one of the
    theorem if the first comes before, in time or, at the same time, in space, and
    minus the transpose of the block of the exchanged points otherwise.
    """
    nu_prime = nu if nu_prime is None else nu_prime

    def block(P, Q):
        (tp, up), (tq, uq) = P, Q
        if tp == tq:
            if up == uq:
                raise ValueError("the points must be pairwise distinct")
            B = kernel_block(abs(uq - up), tp, nu, nu_prime)
            return B if up < uq else -B.T
        if tp < tq:
            return extended_block(uq - up, tp, tq - tp, nu, nu_prime)
        return -extended_block(up - uq, tq, tp - tq, nu, nu_prime).T

    pts = [(float(t), float(u)) for t, u in points]
    return _assemble(pts, lambda P: density(P[0], nu, nu_prime), block)


def correlation(points, nu=1.0, nu_prime=None):
    """rho(tau_1,u_1; ...; tau_n,u_n) of the limit: the Pfaffian of kernel_matrix(points)."""
    return pfaffian(kernel_matrix(points, nu, nu_prime))


# ---------------------------------------------------------------------------
# The lattice, exact.  x is the index of the bond and tau the macroscopic time.
# ---------------------------------------------------------------------------

def _lattice_parameters(N, nu, nu_prime):
    zeta = (N - nu) / (N + nu)                       # the spin correlation across a bond
    zeta0 = (N - nu_prime) / (N + nu_prime)          # the same in the initial measure
    gamma = (N * N - nu * nu) / (N * N + nu * nu)
    return zeta, zeta0, gamma


def _size(n):
    """The size of a table that contains the indices 0, ..., n: a power of 2, at least 256."""
    return max(256, 1 << int(np.ceil(np.log2(n + 2))))


@lru_cache(maxsize=256)
def _k_table(tau, N, nu, nu_prime, size):
    """k(r) = <sigma_0 sigma_r> at the time N^2 tau, for r = 0, ..., size-1.

    k(r) = zeta^r + exp(-4t) sum_j I_{r-j}(4 gamma t) sgn(j) (zeta'^|j| - zeta^|j|),  t = N^2 tau:
    the stationary profile plus the initial discrepancy, carried by the difference of
    two killed walks.
    """
    zeta, zeta0, gamma = _lattice_parameters(N, nu, nu_prime)
    r = np.arange(size)
    if tau == 0:
        out = zeta0 ** r
    elif nu_prime == nu:
        out = zeta ** r
    else:
        t = N * N * tau
        x = 4.0 * gamma * t
        W = int(10.0 * np.sqrt(x)) + 50               # beyond it the Bessel weights are negligible
        q = np.exp(-4.0 * t + x) * ive(np.arange(-W, W + 1), x)        # exp(-4t) I_n(4 gamma t)
        j = np.arange(-W, size + W)
        datum = np.sign(j) * (zeta0 ** np.abs(j) - zeta ** np.abs(j))  # odd in j
        out = zeta ** r + np.convolve(datum, q, mode="valid")
    out.setflags(write=False)
    return out


@lru_cache(maxsize=256)
def _propagated_table(sigma, theta, N, nu, nu_prime, size):
    """G(r) = sum_j p(r-j) sgn(j) k_sigma(|j|) and p(r) = exp(-2t) I_r(2 gamma t), t = N^2 theta.

    Both for r = 0, ..., size-1; G is odd and p is even.
    """
    _, _, gamma = _lattice_parameters(N, nu, nu_prime)
    t = N * N * theta
    x = 2.0 * gamma * t
    W = int(10.0 * np.sqrt(x)) + 50
    weight = np.exp(-2.0 * t + x)
    p = weight * ive(np.arange(-W, W + 1), x)
    j = np.arange(-W, size + W)
    k = _k_table(sigma, N, nu, nu_prime, _size(size + W))
    G = np.convolve(np.sign(j) * k[np.abs(j)], p, mode="valid")
    p = weight * ive(np.arange(size), x)
    G.setflags(write=False)
    p.setflags(write=False)
    return G, p


def k_lattice(r, tau, N, nu=1.0, nu_prime=None):
    """<sigma_0 sigma_r> on the lattice of mesh 1/N at the time N^2 tau, r >= 0 integer.

    The error is of the order of 1e-16 times zeta^r, the stationary value.  Far in
    the tail, where the function is much smaller than zeta^r (it happens for nu' > nu),
    its relative precision is lost; the correlation functions do not feel it.
    """
    r = np.atleast_1d(np.asarray(r, dtype=int))
    if np.any(r < 0):
        raise ValueError("k_lattice is defined for r >= 0")
    tau, nu, nu_prime = _key(tau, nu, nu_prime)
    return _k_table(tau, int(N), nu, nu_prime, _size(int(r.max())))[r]


def lattice_density(tau, N, nu=1.0, nu_prime=None):
    """N times the probability of an interface on a bond at the time N^2 tau: N (1 - k(1))/2."""
    return float(N * 0.5 * (1.0 - k_lattice(1, tau, N, nu, nu_prime)[0]))


def lattice_block(r, tau, N, nu=1.0, nu_prime=None):
    """The block K(x, x+r) of Garrod at one time, r >= 1, rescaled.

    K = -1/2 [[k, D2 k], [D1 k, D1 D2 k]], where D2 moves the second bond to the right
    (r+1) and D1 the first one (r-1).
    """
    if r < 1:
        raise ValueError("lattice_block needs r >= 1")
    tau, nu, nu_prime = _key(tau, nu, nu_prime)
    k = _k_table(tau, int(N), nu, nu_prime, _size(r + 1))
    d2 = k[r + 1] - k[r]
    d1 = k[r - 1] - k[r]
    d12 = 2.0 * k[r] - k[r - 1] - k[r + 1]
    return -0.5 * np.array([[k[r], N * d2], [N * d1, N * N * d12]])


def lattice_extended_block(r, sigma, theta, N, nu=1.0, nu_prime=None):
    """The block K(sigma, x; sigma+theta, x+r) of the extended kernel of Garrod, rescaled.

    r is any integer and theta > 0.  With p(r) = exp(-2t) I_r(2 gamma t), t = N^2 theta, the
    Green function of gamma Delta - 2 (1 - gamma), and G = p * (odd extension of k_sigma),

        K_11 = -1/2 [ G - p ](r)
        K_12 = -1/2 [ D2 G - p(r+1) + p(r) ]
        K_21 = -1/2 [ D1 G + p(r-1) + p(r) ]
        K_22 = -1/2 [ D1 D2 G + p(r+1) - p(r-1) ].
    """
    if not theta > 0:
        raise ValueError("lattice_extended_block needs theta > 0; at equal times use lattice_block")
    sigma, nu, nu_prime = _key(sigma, nu, nu_prime)
    Gt, pt = _propagated_table(sigma, float(theta), int(N), nu, nu_prime, _size(abs(r) + 1))

    def G(m):
        return np.sign(m) * Gt[abs(m)]

    def p(m):
        return pt[abs(m)]

    K11 = -0.5 * (G(r) - p(r))
    K12 = -0.5 * (G(r + 1) - G(r) - p(r + 1) + p(r))
    K21 = -0.5 * (G(r - 1) - G(r) + p(r - 1) + p(r))
    K22 = -0.5 * (2.0 * G(r) - G(r - 1) - G(r + 1) + p(r + 1) - p(r - 1))
    return np.array([[K11, N * K12], [N * K21, N * N * K22]])


def lattice_kernel_matrix(points, N, nu=1.0, nu_prime=None):
    """The rescaled antisymmetric matrix of the blocks at the points (tau_i, x_i) of the lattice.

    x_i is the index of a bond.  The Pfaffian of this matrix is N^n times the
    correlation function, see lattice_correlation.
    """
    def block(P, Q):
        (tp, xp), (tq, xq) = P, Q
        if tp == tq:
            if xp == xq:
                raise ValueError("the points must be pairwise distinct")
            B = lattice_block(abs(xq - xp), tp, N, nu, nu_prime)
            return B if xp < xq else -B.T
        if tp < tq:
            return lattice_extended_block(xq - xp, tp, tq - tp, N, nu, nu_prime)
        return -lattice_extended_block(xp - xq, tq, tp - tq, N, nu, nu_prime).T

    pts = [(float(t), int(x)) for t, x in points]
    return _assemble(pts, lambda P: lattice_density(P[0], N, nu, nu_prime), block)


def lattice_correlation(points, N, nu=1.0, nu_prime=None):
    """N^n E[eta_{x_1}(N^2 tau_1) ... eta_{x_n}(N^2 tau_n)] on the lattice: a Pfaffian.

    points = [(tau_1, x_1), ...], with x_i the index of a bond.  As N grows, at the
    bonds x_i = floor(N u_i), it converges to correlation([(tau_i, u_i), ...]).
    """
    return pfaffian(lattice_kernel_matrix(points, N, nu, nu_prime))


# ---------------------------------------------------------------------------
# The Pfaffian
# ---------------------------------------------------------------------------

def pfaffian(A):
    """The Pfaffian of an antisymmetric matrix, real or complex.

    The matrix is reduced to the tridiagonal form by congruences, with pivoting
    (Parlett and Reid); the Pfaffian is then the product of the entries (2i, 2i+1).
    """
    A = np.array(A, dtype=complex if np.iscomplexobj(A) else float)
    n = A.shape[0]
    if A.shape != (n, n):
        raise ValueError("the matrix must be square")
    if n == 0:
        return 1.0
    if np.max(np.abs(A + A.T)) > 1e-10 * np.max(np.abs(A)):
        raise ValueError("the matrix must be antisymmetric")
    if n % 2 == 1:
        return 0.0
    pf = 1.0
    for k in range(0, n - 1, 2):
        # the largest entry of the row k, to the right of the diagonal, goes to the column k+1
        kp = k + 1 + int(np.argmax(np.abs(A[k, k + 1:])))
        if kp != k + 1:
            A[[k + 1, kp], :] = A[[kp, k + 1], :]
            A[:, [k + 1, kp]] = A[:, [kp, k + 1]]
            pf = -pf
        if A[k, k + 1] == 0:
            return 0.0
        pf = pf * A[k, k + 1]
        if k + 2 < n:
            row = A[k, k + 2:] / A[k, k + 1]
            A[k + 2:, k + 2:] += np.outer(row, A[k + 2:, k + 1]) - np.outer(A[k + 2:, k + 1], row)
    return pf


# ---------------------------------------------------------------------------
# From the command line
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="The Pfaffian kernel of the interfaces at the given space-time points.")
    parser.add_argument("points", nargs="+", metavar="tau,u", help="a space-time point")
    parser.add_argument("--nu", type=float, default=1.0, help="the intensity of the dynamics")
    parser.add_argument("--nu-prime", type=float, default=None,
                        help="the initial intensity (default: nu, the stationary case)")
    parser.add_argument("--N", type=int, default=None,
                        help="print also the rescaled kernel of the lattice of mesh 1/N")
    args = parser.parse_args()

    points = [tuple(float(x) for x in point.split(",")) for point in args.points]
    np.set_printoptions(precision=6, suppress=True, linewidth=150)
    print("points (tau, u):", points)
    print("the kernel of the limit:")
    print(kernel_matrix(points, args.nu, args.nu_prime))
    print(f"rho = {correlation(points, args.nu, args.nu_prime):.10g}")
    if args.N:
        N = args.N
        sites = [(tau, int(np.floor(N * u + 1e-9))) for tau, u in points]      # 100 * 0.29 is 28.999...
        print(f"\nthe rescaled kernel of the lattice, N = {N}, at the bonds {[x for _, x in sites]}:")
        print(lattice_kernel_matrix(sites, N, args.nu, args.nu_prime))
        print(f"N^{len(points)} rho^N = {lattice_correlation(sites, N, args.nu, args.nu_prime):.10g}")


if __name__ == "__main__":
    main()
