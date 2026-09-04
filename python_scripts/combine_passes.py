"""Combine multiple multi-h BoxCorr passes (same params, different seeds) into a
high-statistics estimate, PER box width.  Each PASS_DIR holds h<value>/ subdirs
(box_corr.txt + configuration.json).  Output: OUT_DIR/h<value>/{box_corr.txt,
configuration.json}, usable by the existing validators/plotters.
Usage: combine_passes.py OUT_DIR PASS_DIR...
"""
import json, os, sys, glob, subprocess
import numpy as np

out = sys.argv[1]
passdirs = sys.argv[2:]

# discover the set of h-subdirs from any pass that has completed ones
hnames = []
for d in passdirs:
    for s in sorted(glob.glob(os.path.join(d, "h*"))):
        hn = os.path.basename(s)
        if os.path.exists(os.path.join(s, "box_corr.txt")) and hn not in hnames:
            hnames.append(hn)
if not hnames:
    print("no completed passes yet"); sys.exit(0)

for hn in hnames:
    arrs, cfg = [], None
    for d in passdirs:
        bp = os.path.join(d, hn, "box_corr.txt")
        if os.path.exists(bp):
            arrs.append(np.atleast_2d(np.loadtxt(bp)))
            if cfg is None:
                cfg = json.load(open(os.path.join(d, hn, "configuration.json")))
    if not arrs:
        continue
    K = len(cfg["T_values"]); npass = len(arrs); M = arrs[0].shape[0]
    comb = np.zeros((M, 1 + 2 * K)); comb[:, 0] = arrs[0][:, 0]
    for k in range(K):
        mcs = np.stack([a[:, 1 + 2 * k] for a in arrs])
        ses = np.stack([a[:, 2 + 2 * k] for a in arrs])
        comb[:, 1 + 2 * k] = mcs.mean(axis=0)
        comb[:, 2 + 2 * k] = np.sqrt((ses ** 2).sum(axis=0)) / npass
    od = os.path.join(out, hn); os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "box_corr.txt"), "w") as f:
        f.write("# x")
        for T in cfg["T_values"]:
            f.write(f"\tMC_T{T:g}\tstderr_T{T:g}")
        f.write("\n")
        for row in comb:
            f.write(f"{row[0]:.6f}")
            for k in range(K):
                f.write(f"\t{row[1 + 2 * k]:.8f}\t{row[2 + 2 * k]:.8f}")
            f.write("\n")
    cfg2 = dict(cfg); cfg2["nsims"] = int(cfg["nsims"]) * npass
    json.dump(cfg2, open(os.path.join(od, "configuration.json"), "w"), indent=2)
    print(f"{hn}: combined {npass} passes -> {od}  (total nsims/pt = {cfg2['nsims']})")

# Refresh the per-h figures from the freshly combined data (MC vs continuum f).
hdirs = sorted(d for d in glob.glob(os.path.join(out, "h*")) if os.path.isdir(d))
if hdirs:
    plotdir = os.path.join(os.path.dirname(out) or ".", "plots")
    os.makedirs(plotdir, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    r = subprocess.run([sys.executable, os.path.join(here, "plot_box_vary_h.py"),
                        *hdirs, "--out", os.path.join(plotdir, "box.png")],
                       capture_output=True, text=True)
    print(f"plots -> {plotdir}/box_h*.png"
          + ("" if r.returncode == 0 else f"  (plot warning: {r.stderr.strip()[:200]})"))
