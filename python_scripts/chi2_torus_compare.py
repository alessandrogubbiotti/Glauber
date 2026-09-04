"""Full-curve chi2/dof of the box iface correlator vs the bare line f and vs
the image-summed ring f_torus, per h and tau.  Quantifies how much of the
'h=2 small-tau deviation' is really the finite-ring wrap-around."""
import os, sys, glob, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_continuum as tc

combined = sys.argv[1] if len(sys.argv) > 1 else "grid_T05/combined"
for d in sorted(glob.glob(os.path.join(combined, "h*"))):
    bp = os.path.join(d, "box_corr.txt")
    if not os.path.exists(bp):
        continue
    cfg = json.load(open(os.path.join(d, "configuration.json")))
    raw = np.atleast_2d(np.loadtxt(bp)); x = raw[:, 0]
    L = cfg["L"]; h = cfg["box_width_h"]; K = len(cfg["T_values"])
    sep = np.where((x - L / 2) % L > L / 2, (x - L / 2) % L - L, (x - L / 2) % L)
    print(f"\n{os.path.basename(d)}  h={h:g}  (chi2/dof: line -> torus)")
    for k, T in enumerate(cfg["T_values"]):
        mc = raw[:, 1 + 2 * k]; se = raw[:, 2 + 2 * k]
        cl = tc.box_average(tc.f, sep, T, h)
        ct = tc.box_average(lambda xx, tt: tc.f_torus(xx, tt, L), sep, T, h)
        m = se > 0; dof = max(m.sum(), 1)
        chil = np.sum(((mc[m] - cl[m]) / se[m]) ** 2) / dof
        chit = np.sum(((mc[m] - ct[m]) / se[m]) ** 2) / dof
        print(f"  tau={T:>4.3g}:  {chil:7.2f} -> {chit:7.2f}")
