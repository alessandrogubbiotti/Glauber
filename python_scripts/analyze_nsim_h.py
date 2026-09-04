"""Analyze a set of BoxCorr run dirs: chi^2/dof vs EXACT finite-N theory
(ring G for spin, H_interface for iface) and how the error bar scales with the
box width h at fixed nsims => the nsim<->h law.  Usage: analyze_nsim_h.py DIR...
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_glauber as tg


def spin_theory(cfg, centers):
    N, l, eta, gamma = cfg["N"], cfg["l"], cfg["eta"], cfg["gamma"]
    w = (cfg["box_sites"] - 1) // 2; c0 = cfg["reference_site"]
    offs = np.arange(-w, w + 1); sep = offs[:, None] - offs[None, :]
    out = []
    for T in cfg["T_values"]:
        g = tg.G_two_time_ring(np.arange(l), 2.0 * N * N * T, l, eta, gamma)
        out.append(np.array([g[(cj - c0 + sep) % l].mean() for cj in centers]))
    return np.array(out)


def iface_theory(cfg, centers):
    N, l, eta, gamma, p = cfg["N"], cfg["l"], cfg["eta"], cfg["gamma"], cfg["p_equilibrium"]
    w = (cfg["box_sites"] - 1) // 2; c0 = cfg["reference_site"]
    off = np.arange(-w, w + 1); sep2 = off[:, None] - off[None, :]
    out = []
    for T in cfg["T_values"]:
        gall = tg.H_interface(np.arange(l // 2 + 1), 2.0 * N * N * T, eta, gamma) - p * p
        col = []
        for cj in centers:
            dd = np.abs(((cj - c0 + sep2 + l // 2) % l) - l // 2)
            col.append(N * N * gall[dd].mean())
        out.append(np.array(col))
    return np.array(out)


rows = []
for d in sys.argv[1:]:
    cfg = json.load(open(os.path.join(d, "configuration.json")))
    raw = np.atleast_2d(np.loadtxt(os.path.join(d, "box_corr.txt")))
    K = len(cfg["T_values"]); N = cfg["N"]; l = cfg["l"]
    x = raw[:, 0]
    centers = np.round(x * N).astype(int) % l
    mc  = np.stack([raw[:, 1 + 2 * k] for k in range(K)])
    err = np.stack([raw[:, 2 + 2 * k] for k in range(K)])
    obs = cfg["observable"]; h = cfg["box_width_h"]
    th = iface_theory(cfg, centers) if obs == "iface" else spin_theory(cfg, centers)
    for k, T in enumerate(cfg["T_values"]):
        m = err[k] > 0
        chi2 = np.sum(((mc[k][m] - th[k][m]) / err[k][m]) ** 2) / m.sum()
        rows.append((obs, h, T, chi2, err[k][m].max(), cfg["box_sites"], cfg.get("nsims")))

print(f"{'obs':>6} {'h':>5} {'|B|':>5} {'T':>6} {'chi2/dof':>9} {'max_stderr':>12}")
for r in sorted(rows):
    print(f"{r[0]:>6} {r[1]:>5g} {r[5]:>5} {r[2]:>6g} {r[3]:>9.2f} {r[4]:>12.3e}")

print("\nerror-bar scaling at fixed nsims (log-log fit of max_stderr vs h):")
for obs in ("spin", "iface"):
    for T in sorted({r[2] for r in rows if r[0] == obs}):
        pts = sorted((r[1], r[4]) for r in rows if r[0] == obs and r[2] == T)
        hs = np.array([p[0] for p in pts]); es = np.array([p[1] for p in pts])
        if len(hs) >= 2 and (es > 0).all():
            s = np.polyfit(np.log(hs), np.log(es), 1)[0]
            print(f"  {obs:>6} T={T:<5g}  stderr ~ h^{s:+.2f}"
                  f"   => nsim ~ h^{2*s:+.2f} to hold fixed precision")
