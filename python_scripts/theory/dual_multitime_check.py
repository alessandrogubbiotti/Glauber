"""
dual_multitime_check.py  (2026-07-05)

Route-B validation of the multi-time bond Pfaffian (paper Cor. [multi-time
Pfaffian], eq. e:multipf) via the ANNIHILATING-RANDOM-WALK DUAL -- fully
independent of the Jordan-Wigner / Majorana / Pfaffian pipeline.

Duality (paper sec. s:dual), derived from the exact operator identity
    L sigma_A = sum_{x in A} [ -2 sigma_A
                 + gamma ( sigma_{A xor {x-1,x}} + sigma_{A xor {x,x+1}} ) ]:
    e^{tL} sigma_A = E_A[ exp(-2c int_0^t |A_s| ds) sigma_{A_t} ],
where A_s = annihilating random walks: each walker hops left/right at rate
gamma, a walker jumping onto an occupied site annihilates together with the
occupant, and c = 1 - gamma is the per-walker mass rate.

Multi-time nesting (ascending times t_1 <= ... <= t_k, Delta_j = t_j-t_{j-1}):
    G_1 = phi ,   G_{j+1} = T_{Delta_{j+1}} X_{A_j} G_j ,
    E[ sigma_{A_1}(t_1) ... sigma_{A_k}(t_k) ] = G_k(A_k),
with phi(Z) = <sigma_Z>_pi = zeta^{(z_2-z_1)+(z_4-z_3)+...} (0 for odd |Z|),
(X_A g)(S) = g(A xor S)  (sigma_A sigma_S = sigma_{A xor S}), and
T_t = exp(t*Agen), Agen the weighted dual generator above, acting on functions
of the subsets (size <= 2k, fixed parity) of a window of length L.  The window
is a TRUNCATION OF Z, not a physical boundary: the margin is chosen so that a
walker reaches the edge with probability ~< 1e-9 within the total time span.

Naming: the chain/window length is L (no scaling is performed here; N is
reserved in the paper for the scaling parameter).

Validation ladder (each rung independent of the Pfaffian pipeline except [4]):
  [1] k=1 :  <mu> = zeta                                (statics)
  [2] sigma-sigma two-time vs the Z Bessel formula D(n,t)  (fixes conventions)
  [3] k=2 :  <mu mu> = zeta^2 + det(D-block)            (Glauber determinant)
  [4] k=2,3,4 : vs ff_engine i^k Pf[Phi] on the bulk of a long open chain
"""

import sys, time
import numpy as np
from itertools import combinations
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import ive

# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------

def eta_of_gamma(gamma):
    """zeta = tanh(beta) from gamma = tanh(2 beta)."""
    return (1.0 - np.sqrt(1.0 - gamma * gamma)) / gamma


# --------------------------------------------------------------------------
# subset space on a window of length L: all subsets with |S| <= smax and
# |S| = parity (mod 2), encoded as bitmasks (L <= 32).
# --------------------------------------------------------------------------

def build_space(L, smax):
    """All subsets of {0..L-1} with |S| <= smax (both parities: the xor step
    flips parity when |A_j| is odd, e.g. in the sigma-sigma rung)."""
    masks = []
    for m in range(0, smax + 1):
        for comb in combinations(range(L), m):
            mk = 0
            for u in comb:
                mk |= (1 << u)
            masks.append(mk)
    masks = np.array(sorted(masks), dtype=np.int64)
    return masks


def popcount(arr, L):
    pc = np.zeros_like(arr)
    for i in range(L):
        pc += (arr >> i) & 1
    return pc


def build_generator(masks, L, gamma):
    """Sparse weighted dual generator Agen on the subset space:
    hop gamma per direction, annihilation on contact (mask XOR does both:
    target = mask ^ bit_u ^ bit_v covers move AND mutual annihilation),
    mass -2c|S| on the diagonal.  Window edges: jumps out of the window are
    simply absent (truncation of Z; margin makes this irrelevant)."""
    c = 1.0 - gamma
    n = len(masks)
    pc = popcount(masks, L)
    rows, cols, vals = [], [], []
    nmoves = np.zeros(n, dtype=np.int64)
    for u in range(L):
        bu = np.int64(1 << u)
        has_u = (masks & bu) != 0
        for dv in (-1, +1):
            v = u + dv
            if v < 0 or v >= L:
                continue
            bv = np.int64(1 << v)
            src = np.nonzero(has_u)[0]
            tgt_masks = masks[src] ^ bu ^ bv
            tgt = np.searchsorted(masks, tgt_masks)
            assert np.all(masks[tgt] == tgt_masks), "target outside space"
            rows.append(src); cols.append(tgt)
            vals.append(np.full(len(src), gamma))
            nmoves[src] += 1
    rows.append(np.arange(n)); cols.append(np.arange(n))
    vals.append(-gamma * nmoves - 2.0 * c * pc)
    A = csr_matrix((np.concatenate(vals),
                    (np.concatenate(rows), np.concatenate(cols))),
                   shape=(n, n))
    return A


def build_phi(masks, L, zeta):
    """phi(Z) = <sigma_Z>_pi = zeta^{alternating gap sum} (0 for odd |Z|)."""
    out = np.zeros(len(masks))
    for i, mk in enumerate(masks):
        pos = [u for u in range(L) if (mk >> u) & 1]
        if len(pos) % 2 == 1:
            continue
        s = 0
        for j in range(0, len(pos), 2):
            s += pos[j + 1] - pos[j]
        out[i] = zeta ** s
    return out


def xor_op(masks, w, maskA, smax, L):
    """(X_A w)(S) = w(A xor S); entries whose image leaves the space (size
    > smax) are set to 0 -- they are never read downstream (the dual
    generator is block lower-triangular in |S|, and the final evaluation
    only queries small sets)."""
    t = masks ^ np.int64(maskA)
    pc = popcount(t, L)
    ok = pc <= smax
    out = np.zeros_like(w)
    idx = np.searchsorted(masks, t[ok])
    assert np.all(masks[idx] == t[ok])
    out[ok] = w[idx]
    return out


def dual_multitime(sets, times, gamma, L, smax=None):
    """E[ sigma_{sets[0]}(times[0]) ... sigma_{sets[-1]}(times[-1]) ]
    in equilibrium on (a window of) Z, by the nested annihilating-walk dual.
    sets: list of site tuples (window coordinates); times ascending."""
    zeta = eta_of_gamma(gamma)
    k = len(sets)
    tot = sum(len(A) for A in sets)
    if smax is None:
        smax = tot
    masks = build_space(L, smax)
    A = build_generator(masks, L, gamma)
    w = build_phi(masks, L, zeta)
    for j in range(k - 1):
        maskA = 0
        for u in sets[j]:
            maskA |= (1 << u)
        w = xor_op(masks, w, maskA, smax, L)
        dt = times[j + 1] - times[j]
        if dt > 1e-14:
            w = expm_multiply(A * dt, w)
    maskF = 0
    for u in sets[-1]:
        maskF |= (1 << u)
    i = np.searchsorted(masks, np.int64(maskF))
    return w[i]


# --------------------------------------------------------------------------
# independent Z references
# --------------------------------------------------------------------------

def D_bessel(n, t, gamma, K=80):
    """Two-time spin two-point on Z (paper e:Dline):
    D(n,t) = sum_m zeta^{|m|} e^{-2t} I_{n-m}(2 gamma t),
    via scaled Bessel: e^{-2t} I_m(2 g t) = ive(m, 2gt) e^{-2(1-g)t}."""
    zeta = eta_of_gamma(gamma)
    c = 1.0 - gamma
    m = np.arange(-K, K + 1)
    return float(np.sum(zeta ** np.abs(m)
                        * ive(n - m, 2.0 * gamma * t)) * np.exp(-2.0 * c * t))


def mu2_bessel(d, t, gamma):
    """<mu_x(0) mu_{x+d}(t)> = zeta^2 + det [[D(d),D(d+1)],[D(d-1),D(d)]]
    (paper e:pf2time; Glauber's two-time determinant)."""
    zeta = eta_of_gamma(gamma)
    Dd = D_bessel(d, t, gamma)
    return zeta ** 2 + (Dd * Dd - D_bessel(d + 1, t, gamma)
                        * D_bessel(d - 1, t, gamma))


# --------------------------------------------------------------------------
# the ladder
# --------------------------------------------------------------------------

def run(gamma, do_k4=False):
    zeta = eta_of_gamma(gamma)
    print(f"\n===== gamma = {gamma}   (zeta = {zeta:.6f}) =====")

    # [1] k=1 statics --------------------------------------------------------
    L, ctr = 15, 7
    v = dual_multitime([(ctr, ctr + 1)], [0.0], gamma, L)
    print(f"[1] k=1  <mu> dual = {v:+.12f}   zeta = {zeta:+.12f}"
          f"   |diff| = {abs(v - zeta):.2e}")

    # [2] sigma-sigma two-time vs Bessel D ----------------------------------
    print("[2] sigma-sigma two-time vs Z Bessel formula (t = 0.5):")
    t = 0.5
    margin = 13
    for d in (0, 1, 3):
        L = d + 1 + 2 * margin
        ctr = margin
        v = dual_multitime([(ctr,), (ctr + d,)], [0.0, t], gamma, L)
        ref = D_bessel(d, t, gamma)
        print(f"      d={d}: dual = {v:+.10f}   Bessel = {ref:+.10f}"
              f"   |diff| = {abs(v - ref):.2e}")

    # [3] k=2 bond moment vs Glauber determinant ----------------------------
    print("[3] k=2 bond moment vs zeta^2 + det(D-block)  (Bessel):")
    for (d, t) in ((2, 0.6), (1, 0.8), (3, 0.4)):
        margin = 13
        L = d + 2 + 2 * margin
        b1 = margin
        v = dual_multitime([(b1, b1 + 1), (b1 + d, b1 + d + 1)],
                           [0.0, t], gamma, L)
        ref = mu2_bessel(d, t, gamma)
        print(f"      d={d}, t={t}: dual = {v:+.10f}   Bessel = {ref:+.10f}"
              f"   |diff| = {abs(v - ref):.2e}")

    # [4] k=2,3(,4) vs the free-fermion Pfaffian (bulk, long open chain) ----
    from ff_engine import FF
    Mff = 201
    ff = FF(Mff, zeta)
    x0 = Mff // 2

    cases = [("k=2", (0, 2), (0.0, 0.6)),
             ("k=3", (0, 2, 3), (0.0, 0.45, 0.95))]
    if do_k4:
        cases.append(("k=4", (0, 1, 3, 4), (0.0, 0.2, 0.45, 0.7)))

    for name, rel, ts in cases:
        margin = 13 if len(rel) <= 3 else 9
        span = max(rel) + 2
        L = span + 2 * margin
        b0 = margin
        sets = [(b0 + r, b0 + r + 1) for r in rel]
        t0 = time.time()
        v = dual_multitime(sets, list(ts), gamma, L)
        pts = [(x0 + r, t) for r, t in zip(rel, ts)]
        ref = ff.bond_corr(pts)
        print(f"[4] {name}: bonds {rel} times {ts}  (window L={L})")
        print(f"      dual (annihilating walks) = {v:+.10f}")
        print(f"      i^k Pf[Phi]  (ff_engine)   = {ref.real:+.10f}"
              f"  (Im {ref.imag:+.1e})")
        print(f"      |diff| = {abs(v - ref.real):.2e}"
              f"      [{time.time() - t0:.1f}s]")


if __name__ == "__main__":
    do_k4 = "--k4" in sys.argv
    for gamma in (0.6, 0.9):
        run(gamma, do_k4=do_k4 and gamma == 0.6)
