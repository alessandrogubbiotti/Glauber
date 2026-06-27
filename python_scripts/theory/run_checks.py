"""
run_checks.py — explicit, runnable verification gates for the scaling-limit
claims of glauber_paper.tex (Theorem `thm:pfaffian`, boxed density `e:fmain`,
kernel `e:Kblock`).

Each gate prints PASS/FAIL with the measured number and names the paper claim it
backs. These are the FAST analytic / free-fermion gates only (a few seconds);
the heavy box Monte-Carlo collapse (Fig `f:collapse`) and the 3/8 persistence
run are separate long campaigns documented in ../../VERIFICATION.md.

The point of this script is to make the open k>=3 correctness issue impossible to
miss: gate [2] reproduces, on every run, that the kernel currently PRINTED in the
theorem gives an identically-zero connected 3-point, while the corrected complex
kernel matches the discrete limit.

Run:  python run_checks.py        (exit code 0 iff all gates pass)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))                    # theory/
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))   # python_scripts/

import numpy as np
import theory_continuum as tc
from cont3_final import Kcont, Kpaper, fk
import ff_engine as fe

_results = []


def check(name, ok, detail):
    _results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n         {detail}")


def main():
    print("=" * 74)
    print(" Glauber scaling-limit verification gates (fast analytic / free-fermion)")
    print("=" * 74)

    # -- Gate 1: two-point continuum density f (e:fmain), quartered Wronskian ----
    print("\n[1] e:fmain   f = 1/4[(d_xH)^2 - H d_x^2 H]  (quartered Wronskian)")
    ok = True
    for (x, tau) in [(0.3, 0.1), (0.5, 0.2), (0.2, 0.05)]:
        a = fk([(0.0, 0.0), (x, tau)], Kcont).real
        f = tc.f(x, tau)[0]
        ok &= abs(a - f) < 1e-9 * max(1.0, abs(f))
    check("k=2 Pfaffian anchor reproduces the closed-form f",
          ok, "fk(Kcont) == tc.f to < 1e-9 at three (x,tau)")

    # -- Gate 2: k=3 — the printed kernel is identically zero (the open bug) -----
    print("\n[2] e:Kblock  k=3 connected interface 3-point  (Theorem thm:pfaffian)")
    configs = [dict(name="A", xs=[-0.2, 0.0, 0.25], taus=[0.05, 0.10, 0.18]),
               dict(name="C", xs=[0.1, -0.15, 0.3], taus=[0.07, 0.12, 0.05])]
    paper_zero, cont_matches = True, True
    detail_rows = []
    for cfg in configs:
        xs, taus = cfg["xs"], cfg["taus"]
        rows = []
        for N in (20, 40):
            eta = (N - 1) / (N + 1)
            ff = fe.FF(N * 8, eta)
            c0 = (N * 8) // 2
            pts = [(c0 + int(round(N * x)), N * N * t) for x, t in zip(xs, taus)]
            rows.append((1.0 / N, N ** 3 * ff.iface_cumulant(pts).real))
        inv = np.array([r[0] for r in rows])
        dd = np.array([r[1] for r in rows])
        disc = np.linalg.lstsq(np.vstack([np.ones(2), inv]).T, dd, rcond=None)[0][0]
        cc = fk(list(zip(xs, taus)), Kcont).real
        cp = fk(list(zip(xs, taus)), Kpaper).real
        paper_zero &= abs(cp) < 1e-12
        cont_matches &= abs(disc - cc) < 0.08 * abs(cc)
        detail_rows.append(f"cfg {cfg['name']}: discrete={disc:+.3e}  "
                           f"Kcont={cc:+.3e}  Kpaper={cp:+.3e}")
    detail = "\n         ".join(detail_rows)
    check("PRINTED kernel [[H,-H'],[H',-H'']] gives connected-3pt == 0 "
          "(reproduces the open bug)", paper_zero,
          detail + "\n         => Theorem thm:pfaffian is false as printed for k>=3; "
                   "see VERIFICATION.md #1")
    check("CORRECTED complex kernel matches the discrete N^3 limit",
          cont_matches, "Richardson (N=20,40) within 8% at both configs")

    # -- summary ----------------------------------------------------------------
    print("\n" + "=" * 74)
    n_pass, n = sum(_results), len(_results)
    print(f" {n_pass}/{n} gates passed.")
    print(" Not gated here (separate campaigns / open) -- see ../../VERIFICATION.md:")
    print("   * k>=4 continuum; regime-family f & tagged-colour persistence: not run")
    print("   * box-MC collapse (Fig f:collapse) & 3/8 persistence: long campaigns")
    print("=" * 74)
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main())
