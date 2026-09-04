"""
run_checks.py — explicit, runnable verification gates for the scaling-limit
claims of glauber_paper.tex: the finite-N Pfaffian identity (`e:multipf`,
`c:centered`), the boxed density `e:fmain`, and the hierarchy kernel `e:Kpf`
of Theorem `thm:pfaffian`.

Each gate prints PASS/FAIL with the measured number and names the paper claim it
backs.  These are the FAST analytic / free-fermion gates (~10 s on an unloaded
machine); the heavy box Monte-Carlo collapse (Fig `f:collapse`) and the 3/8
persistence run are separate long campaigns documented in ../../VERIFICATION.md,
and the annihilating-walk dual cross-check is `dual_multitime_check.py` (~2 min,
add --k4 for the k=4 rung).

Layout:
  [1] finite-N (pre-scaling): algebra prerequisites; <mu..mu> = i^k Pf[Phi] on
      the exact 2^M operator route AND the closed-form free-fermion engine;
      a randomised 300-configuration fuzz over coinciding times, repeated and
      adjacent bonds, identical (bond,time) pairs, four (M, eta).
  [2] negative controls: the scalar spin kernel D must FAIL, in sign AND in
      magnitude (the anomalous matrix contraction is genuinely needed);
      multi-time spin Wick must FAIL while the equal-time (static Ising) case
      holds.
  [3] centered vs connected (`c:centered`, `e:pfhier` as stated): zero-diagonal
      Pf = CENTERED moment; cumulant = single-2k-cycle sum; they must DIFFER at
      k=4 (an early draft claimed they agree at every k -- regression guard).
  [4] continuum k=2 anchor: fk(Kcont) reproduces the boxed f (`e:fmain`).  This
      is an ALGEBRAIC anchor (Kcont and tc.f share the H, P primitives); the
      closed forms themselves are checked against MC / finite-N exact forms in
      VERIFICATION.md (rows e:Hlimit, e:fmain, s:numerics).
  [5] continuum k=3: (a) the complex kernel `e:Kpf` (printed in the theorem
      since 2026-06-27) matches the discrete N^3 Richardson limit at three
      configurations; (b) NEGATIVE CONTROL -- the superseded real-kernel
      reference of scaling_check_3time (real prefactor, the pre-2026-06-27
      reading that contaminated ff_scale / freefermion_largeN, VERIFICATION.md
      R2) has the WRONG SIGN relative to the discrete limit; (c) its k>=3 guard
      trips.  (fk(real kernel).real == 0 at odd k is structural -- imaginary
      prefactor times real Pfaffian, Remark r:Kstruct -- and is therefore NOT a
      gate: any real kernel passes it.)
  [6] box_fk (the validated box-averaged continuum density, replacement for the
      superseded cont_box_cumulant): equals the old function at k=2, where both
      are valid, and converges to the pointwise fk as h -> 0 at k=3.

Run:  python run_checks.py        (exit code 0 iff all gates pass)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))                    # theory/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))   # python_scripts/

import numpy as np
import theory_continuum as tc
from cont3_final import Kcont, Kblock, fk, box_fk
import ff_engine as fe
import pfaffian_verify as pv
import scaling_check_3time as sc

_results = []


def check(name, ok, detail):
    _results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n         {detail}")


def discrete_k3_limit(xs, taus, Ns=(20, 40)):
    """N^3 x connected 3-point from the free-fermion engine, Richardson in 1/N."""
    rows = []
    for N in Ns:
        eta = (N - 1) / (N + 1)
        ff = fe.FF(N * 8, eta)
        c0 = (N * 8) // 2
        pts = [(c0 + int(round(N * x)), N * N * t) for x, t in zip(xs, taus)]
        rows.append((1.0 / N, N ** 3 * ff.iface_cumulant(pts).real))
    inv = np.array([r[0] for r in rows])
    dd = np.array([r[1] for r in rows])
    return np.linalg.lstsq(np.vstack([np.ones(len(Ns)), inv]).T, dd, rcond=None)[0][0]


def main():
    print("=" * 74)
    print(" Glauber scaling-limit verification gates (fast analytic / free-fermion)")
    print("=" * 74)

    # -- Gate 1: finite-N identity e:multipf (the theorem BEFORE the limit) -----
    print("\n[1] e:multipf   <mu_x1(t1)..mu_xk(tk)> = i^k Pf[Phi]   (finite N)")
    m = pv.FermiModel(6, pv.DEFAULT_ETA)
    alg = pv.algebra_deviations(m)
    check("algebra prerequisites (Clifford, mu=iaa, reversibility, stationarity)",
          max(alg.values()) < 1e-12,
          "  ".join(f"{k}={v:.1e}" for k, v in alg.items()))
    dev = pv.pfaffian_deviations(m)
    check("operator route: |<mu..mu> - i^k Pf[Phi]| at k=2,3,4",
          max(dev.values()) < 1e-12,
          "  ".join(f"k={k}: {v:.1e}" for k, v in dev.items()))
    ffworst = fe.validate(models=((6, 1 / 3.), (8, 0.9)), verbose=False)
    check("free-fermion engine builders (closed-form C, Gsp) vs brute force",
          ffworst < 1e-12, f"worst over (6,1/3),(8,0.9): {ffworst:.1e}")
    worst, cov = pv.fuzz_deviations()
    check(f"randomised fuzz ({cov['configs']} configs: {cov['equal_times']} w/ "
          f"coinciding times, {cov['repeated_bond']} w/ repeated bond, "
          f"{cov['identical_pair']} w/ identical (bond,time) pair)",
          worst < 1e-12 and cov["equal_times"] > 0 and cov["repeated_bond"] > 0
          and cov["identical_pair"] > 0,
          f"worst |exact - i^k Pf| = {worst:.1e} over operator AND engine routes, "
          f"4 (M,eta), k=2,3,4")

    # -- Gate 2: negative controls ----------------------------------------------
    print("\n[2] negative controls   (a wrong kernel MUST fail; if these 'pass'")
    print("    accidentally, the positive checks above lose their meaning)")
    Dg = pv.D_universality_gaps(m)
    check("scalar spin kernel: Pf[D] FAILS at every k, in sign-insensitive "
          "magnitude too (anomalous part needed, not a sign convention)",
          min(g[2] for g in Dg.values()) > 1e-8 and min(g[3] for g in Dg.values()) > 1e-8,
          "  ".join(f"k={k}: gap {g[2]:.1e} (|mag| {g[3]:.1e})" for k, g in Dg.items()))
    wick = pv.spin_wick_gaps(m)
    eq_gap = wick["4 spins, equal time"][2]
    mt_gap = min(g[2] for lbl, g in wick.items() if "multi-time" in lbl)
    check("spin Wick: equal-time HOLDS (static Ising), multi-time FAILS "
          "(sigma is a JW string; the bond mu is the quasi-free object)",
          eq_gap < 1e-12 and mt_gap > 1e-8,
          f"equal-time gap {eq_gap:.1e}; smallest multi-time gap {mt_gap:.1e}")

    # -- Gate 3: centered vs connected (c:centered / e:pfhier reading) ----------
    print("\n[3] c:centered   zero-diagonal Pf = CENTERED moment; cumulant = cycle sum")
    import k4_zerodiag_check as kz
    from multitime_interface_dpp import pfaffian
    m8 = pv.FermiModel(8, 0.6)
    ff8 = fe.FF(8, 0.6)
    pts4 = [(1, 0.00), (2, 0.35), (4, 0.05), (5, 0.40)]
    mom = kz.subset_moments(m8, pts4)
    S = (0, 1, 2, 3)
    kap = kz.cumulant(mom, S)
    cen = kz.centered_moment(mom, S)
    Phi = kz.phi_matrix(ff8, pts4)
    pf_zd = (1j ** 4) * pfaffian(kz.zero_diag(Phi))
    cyc = (1j ** 4) * kz.classified_matching_sums(Phi)["cycle"]
    check("k=4: i^4 Pf[zerodiag] = centered moment  AND  cycle sum = cumulant",
          abs(pf_zd - cen) < 5e-12 and abs(cyc - kap) < 5e-12,
          f"|Pf_zd - centered| = {abs(pf_zd - cen):.1e}   "
          f"|cycle - cumulant| = {abs(cyc - kap):.1e}")
    check("k=4: centered and cumulant DIFFER (early 'at every k' claim was wrong;"
          " e:pfhier must stay stated as centered)",
          abs(pf_zd - kap) > 1e-5,
          f"|Pf_zd - cumulant| = {abs(pf_zd - kap):.1e} "
          f"(the 2+2 pair products; k4_zerodiag_check.py for k=3,5 and eta=0.9)")

    # -- Gate 4: two-point continuum density f (e:fmain), quartered Wronskian ----
    print("\n[4] e:fmain   f = 1/4[(d_xH)^2 - H d_x^2 H]  (quartered Wronskian)")
    ok = True
    for (x, tau) in [(0.3, 0.1), (0.5, 0.2), (0.2, 0.05)]:
        a = fk([(0.0, 0.0), (x, tau)], Kcont).real
        f = tc.f(x, tau)[0]
        ok &= abs(a - f) < 1e-9 * max(1.0, abs(f))
    check("k=2 Pfaffian anchor reproduces the closed-form f (algebraic anchor: "
          "same H,P primitives; closed forms vs MC are in VERIFICATION.md)",
          ok, "fk(Kcont) == tc.f to < 1e-9 at three (x,tau)")

    # -- Gate 5: k=3 -- complex kernel right, superseded reference wrong-signed --
    print("\n[5] e:Kpf   k=3 connected interface 3-point  (thm:pfaffian)")
    configs = [dict(name="A", xs=[-0.2, 0.0, 0.25], taus=[0.05, 0.10, 0.18]),
               dict(name="B", xs=[0.0, 0.15, 0.35], taus=[0.04, 0.09, 0.15]),
               dict(name="C", xs=[0.1, -0.15, 0.3], taus=[0.07, 0.12, 0.05])]
    cont_matches, rows = True, []
    disc_of, bad_of = {}, {}
    for cfg in configs:
        xs, taus = cfg["xs"], cfg["taus"]
        disc = discrete_k3_limit(xs, taus)
        cc = fk(list(zip(xs, taus)), Kcont).real
        bad = sc.cont_interface_cumulant(list(zip(xs, taus)), L=None, allow_wrong_k3=True)
        disc_of[cfg["name"]], bad_of[cfg["name"]] = disc, bad
        cont_matches &= abs(disc - cc) < 0.08 * abs(cc)
        rows.append(f"cfg {cfg['name']}: discrete={disc:+.3e}  Kpf={cc:+.3e} "
                    f"(r={disc/cc:.3f})  superseded={bad:+.3e} (r={bad/disc:+.3f})")
    detail = "\n         ".join(rows)
    check("complex kernel e:Kpf (printed in the theorem) matches the discrete "
          "N^3 limit at A, B, C", cont_matches,
          detail + "\n         Richardson (N=20,40) within 8% at every config")
    check("NEGATIVE CONTROL: superseded real-kernel reference (scaling_check_3time, "
          "real prefactor) has the WRONG SIGN vs the discrete limit at A and B "
          "(the R2 contamination)",
          bad_of["A"] * disc_of["A"] < 0 and bad_of["B"] * disc_of["B"] < 0,
          f"sign(superseded*discrete): A {np.sign(bad_of['A']*disc_of['A']):+.0f}, "
          f"B {np.sign(bad_of['B']*disc_of['B']):+.0f}  (C is only 10% off -- "
          f"the error is config-dependent, which is what made it hard to spot)")
    pts3 = [(0.0, 0.0), (0.1, 0.05), (0.2, 0.1)]
    tripped = tripped_box = False
    try:
        sc.cont_interface_cumulant(pts3)
    except RuntimeError:
        tripped = True
    try:
        sc.cont_box_cumulant([0.0, 0.1, 0.2], [0.0, 0.05, 0.1], 0.1)
    except RuntimeError:
        tripped_box = True
    k2_ok = abs(sc.cont_interface_cumulant([(0.0, 0.0), (0.3, 0.1)], L=None)
                - tc.f(0.3, 0.1)[0]) < 1e-9
    check("guard: scaling_check_3time k>=3 continuum functions RAISE without "
          "allow_wrong_k3; k=2 still equals f",
          tripped and tripped_box and k2_ok,
          f"cont_interface_cumulant raised: {tripped}; cont_box_cumulant raised: "
          f"{tripped_box}; k=2 == f: {k2_ok}")

    # -- Gate 6: box_fk -----------------------------------------------------------
    print("\n[6] box_fk   box-averaged continuum density (replacement for the "
          "superseded cont_box_cumulant)")
    dmax = 0.0
    for c, t, h in [((0.0, 0.3), (0.0, 0.1), 0.25), ((0.0, 0.5), (0.0, 0.2), 0.4)]:
        a = box_fk(c, t, h, nq=5)
        b = sc.cont_box_cumulant(list(c), list(t), h, L=None, nq=5)   # k=2: valid
        dmax = max(dmax, abs(a - b))
    check("k=2: box_fk == old cont_box_cumulant (same nodes, same average) where "
          "both are valid", dmax < 1e-12, f"max |diff| = {dmax:.1e} at two (c,t,h)")
    xs, taus = configs[0]["xs"], configs[0]["taus"]
    pw = fk(list(zip(xs, taus)), Kcont).real
    rel = [(h, (box_fk(xs, taus, h) - pw) / pw) for h in (0.1, 0.02)]
    check("k=3: box_fk -> pointwise fk as h -> 0 (monotone, < 1% at h=0.02)",
          abs(rel[1][1]) < 1e-2 and abs(rel[1][1]) < abs(rel[0][1]),
          "  ".join(f"h={h}: rel {r:+.1e}" for h, r in rel) + f"   (pointwise {pw:+.4e})")

    # -- summary ----------------------------------------------------------------
    print("\n" + "=" * 74)
    n_pass, n = sum(_results), len(_results)
    print(f" {n_pass}/{n} gates passed.")
    print(" Not gated here (separate campaigns / independent routes) -- see")
    print(" ../../VERIFICATION.md:")
    print("   * dual (annihilating-walk) cross-check: python dual_multitime_check.py")
    print("   * k>=4 continuum limit; regime-family f & tagged persistence: not run")
    print("   * box-MC collapse (Fig f:collapse) & 3/8 persistence: long campaigns")
    print("=" * 74)
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
