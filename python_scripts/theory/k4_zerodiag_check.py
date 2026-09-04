"""
k4_zerodiag_check.py — what does the ZERO-DIAGONAL multi-time Pfaffian equal:
the joint CUMULANT or the CENTERED MOMENT of the bond/interface indicators?

Background.  The moment identity (e:multipf)

    <mu_{x1}(t1) ... mu_{xk}(tk)> = i^k Pf[Phi]        (full 2k x 2k contraction)

is validated to 1e-15 (operator route) / 4e-16 (dual route).  An EARLIER DRAFT
of the paper's scaling theorem (thm:pfaffian(ii), e:pfhier) stated that the
CONNECTED (cumulant) k-point density is the Pfaffian with the k within-bond 2x2
diagonal blocks set to zero, "the within-bond self-contractions are exactly what
the cumulant subtraction removes, at every k".  That is false from k = 4 on;
since 2026-07-05 the paper states e:pfhier as the CENTERED moment (c:centered),
which is what this script established.  It is kept as the regression guard.

Combinatorial issue: perfect matchings of the 2k legs with no within-bond pair
form a single 2k-cycle only for k = 2, 3.  From k = 4 on, disconnected
matchings survive (e.g. bonds 1-2 paired among themselves and 3-4 among
themselves), so the zero-diagonal Pfaffian should equal the CENTERED moment

    E[ prod_i (mu_i - <mu_i>) ]  =  kappa_k  +  (disconnected cumulant products),

which coincides with the cumulant kappa_k only for k <= 3.  The cumulant itself
should be the single-2k-cycle part of the matching sum.

This script settles it EXACTLY at finite open chain (no scaling, no MC):
  ground truth  : FermiModel.bond_corr (2^M operator route, no fermionic input)
                  -> all subset moments -> cumulants (Moebius) and centered moments;
  Pfaffian side : ff_engine.FF.phi contraction (validated vs brute force),
                  full Pf, zero-diagonal Pf, and a full signed-matching
                  enumeration classified by the induced partition of the bonds
                  (self-pairs / disconnected / single cycle).

Checks, k = 3, 4, 5:
  [A] i^k Pf[Phi_full]      == moment                    (e:multipf sanity)
  [B] i^k Pf[Phi_zerodiag]  == centered moment           (the corrected reading)
  [C] i^k (single-cycle sum)== cumulant                  (the connected object)
  [D] zerodiag - cumulant   == sum over pair-partitions of cumulant products
      (identically 0 at k=3; NONZERO at k=4,5 -> the old draft's "at every k" fails)

Run:  python k4_zerodiag_check.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from itertools import combinations
from math import factorial
import numpy as np

import pfaffian_verify as pv
import ff_engine as fe
from multitime_interface_dpp import pfaffian, _set_partitions

TOL_ID = 5e-12          # identities expected at machine precision
FAILURES = []


def check(label, val, tol=TOL_ID):
    ok = abs(val) < tol
    if not ok:
        FAILURES.append(label)
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}: {val:.3e}")
    return ok


# ---- moment-side objects (operator route, exact) ---------------------------
def subset_moments(m, points):
    """All subset moments E[prod_{i in S} mu_i(t_i)] from the 2^M operator route."""
    n = len(points)
    mom = {(): 1.0 + 0j}
    for r in range(1, n + 1):
        for S in combinations(range(n), r):
            mom[S] = m.bond_corr([points[i] for i in S])
    return mom


def cumulant(mom, S):
    tot = 0.0 + 0j
    for part in _set_partitions(list(S)):
        coef = (-1) ** (len(part) - 1) * factorial(len(part) - 1)
        pr = 1.0 + 0j
        for blk in part:
            pr *= mom[tuple(sorted(blk))]
        tot += coef * pr
    return tot


def centered_moment(mom, S):
    """E[prod_{i in S}(mu_i - <mu_i>)] by inclusion-exclusion over subsets."""
    S = tuple(S)
    tot = 0.0 + 0j
    for r in range(len(S) + 1):
        for T in combinations(S, r):
            rest = tuple(x for x in S if x not in T)
            pr = (-1) ** len(T)
            for i in T:
                pr *= mom[(i,)]
            tot += pr * mom[rest]
    return tot


# ---- Pfaffian-side objects --------------------------------------------------
def phi_matrix(ff, points):
    """Time-sorted 2k x 2k contraction, legs (2x+1, 2x+2) per bond (as FF.bond_corr)."""
    pts = sorted(points, key=lambda p: p[1])
    legs = []
    for (x, t) in pts:
        legs.append((2 * x + 1, t)); legs.append((2 * x + 2, t))
    return ff.phi(legs)


def zero_diag(Phi):
    A = Phi.copy()
    k = A.shape[0] // 2
    for p in range(k):
        A[2 * p:2 * p + 2, 2 * p:2 * p + 2] = 0.0
    return A


def matchings(n):
    """All perfect matchings of range(n) with their permutation signs."""
    def rec(rem):
        if not rem:
            yield [], 1
            return
        i = rem[0]
        for jx, j in enumerate(rem[1:]):
            rest = rem[1:jx + 1] + rem[jx + 2:]
            # sign: moving j next to i across jx intervening elements
            for sub, s in rec(rest):
                yield [(i, j)] + sub, s * ((-1) ** jx)
    yield from rec(list(range(n)))


def classified_matching_sums(Phi):
    """Split Pf[Phi] = sum over signed matchings into classes by the induced
    multigraph on the k bonds: 'self' (has a within-bond pair), 'cycle'
    (single cycle through all k bonds), 'disc' (no self-pair, disconnected)."""
    n = Phi.shape[0]
    sums = {'self': 0j, 'cycle': 0j, 'disc': 0j}
    for pairs, sgn in matchings(n):
        term = sgn * np.prod([Phi[i, j] for (i, j) in pairs])
        bonds = [(i // 2, j // 2) for (i, j) in pairs]
        if any(a == b for (a, b) in bonds):
            cls = 'self'
        else:
            # union-find on k bonds
            k = n // 2
            parent = list(range(k))
            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]; x = parent[x]
                return x
            for (a, b) in bonds:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
            cls = 'cycle' if len({find(x) for x in range(k)}) == 1 else 'disc'
        sums[cls] += term
    return sums


# ---- the test ---------------------------------------------------------------
def run(M, eta, points, label):
    k = len(points)
    print(f"\n== k={k}  {label}  (open chain M={M}, eta={eta}) ==")
    m = pv.FermiModel(M, eta)
    ff = fe.FF(M, eta)

    mom = subset_moments(m, points)
    S = tuple(range(k))
    kap = cumulant(mom, S)
    cen = centered_moment(mom, S)

    Phi = phi_matrix(ff, points)
    pf_full = (1j ** k) * pfaffian(Phi)
    pf_zd = (1j ** k) * pfaffian(zero_diag(Phi))
    cls = classified_matching_sums(Phi)
    cyc = (1j ** k) * cls['cycle']
    dsc = (1j ** k) * cls['disc']

    # disconnected pair-partition products of cumulants (blocks of size >= 2)
    disc_prod = 0j
    for part in _set_partitions(list(S)):
        if len(part) == 1 or any(len(b) < 2 for b in part):
            continue
        pr = 1.0 + 0j
        for blk in part:
            pr *= cumulant(mom, tuple(sorted(blk)))
        disc_prod += pr

    print(f"  moment          = {mom[S].real:+.12e}")
    print(f"  cumulant        = {kap.real:+.12e}")
    print(f"  centered moment = {cen.real:+.12e}")
    print(f"  i^k Pf[full]    = {pf_full.real:+.12e}")
    print(f"  i^k Pf[zerodiag]= {pf_zd.real:+.12e}")
    print(f"  i^k cycle sum   = {cyc.real:+.12e}")
    print(f"  i^k 2+2/3+2 sum = {dsc.real:+.12e}   (Sum prod kappa = {disc_prod.real:+.12e})")

    check("[A] Pf[full] - moment       ", abs(pf_full - mom[S]))
    check("[B] Pf[zerodiag] - centered ", abs(pf_zd - cen))
    check("[C] cycle sum - cumulant    ", abs(cyc - kap))
    check("[D'] disc class - Sum kappa-products", abs(dsc - disc_prod))
    gap = abs(pf_zd - kap)
    rel = gap / max(abs(kap), abs(pf_zd), 1e-300)
    verdict = "COINCIDE (k<=3 accident)" if gap < TOL_ID else \
        f"DIFFER: zerodiag Pf is NOT the cumulant (rel. offset {rel:.2f})"
    print(f"    [D] |Pf[zerodiag] - cumulant| = {gap:.3e}  -> {verdict}")
    return gap


if __name__ == "__main__":
    # two adjacent close pairs, staggered times, well separated along the chain
    A, B = (1, 0.00), (2, 0.35)
    C, D = (4, 0.05), (5, 0.40)
    E = (3, 0.20)

    for eta in (0.6, 0.9):
        g3 = run(8, eta, [A, B, C],       "three points (control: no 2+2 exists)")
        g4 = run(8, eta, [A, B, C, D],    "two close pairs (2+2 matchings exist)")
        g5 = run(8, eta, [A, B, C, D, E], "five points (2+3 matchings exist)")

    print("\n================================================================")
    if FAILURES:
        print("SOME IDENTITY CHECKS FAILED:", FAILURES)
        sys.exit(1)
    print("ALL IDENTITY CHECKS PASS:")
    print("  zero-diagonal Pf = CENTERED moment (every k);")
    print("  cumulant = single-cycle sum;  they coincide only for k <= 3.")
    print("  -> thm:pfaffian(ii)/e:pfhier is (correctly) STATED AS CENTERED in the")
    print("     paper since 2026-07-05; this check guards against reverting it to")
    print("     'cumulant at every k' (the old draft error).")
