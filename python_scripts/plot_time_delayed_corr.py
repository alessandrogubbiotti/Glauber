"""
plot_time_delayed_corr.py  —  G_n(T) = <sigma_0(0) sigma_n(T)>
Called as:  python3 python_scripts/plot_time_delayed_corr.py <outdir>
Reads <outdir>/corr2pt_multi.txt and saves <outdir>/time_delayed_corr.png.
"""

import sys, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: plot_time_delayed_corr.py <outdir>")
    sys.exit(1)

outdir = sys.argv[1]
data_path = os.path.join(outdir, "corr2pt_multi.txt")

# parse header to get T values
with open(data_path) as f:
    header = f.readline().strip()

tokens = header.lstrip("#").split("\t")
T_vals = []
for tok in tokens:
    if tok.startswith("G_T"):
        T_vals.append(float(tok[3:]))

raw = np.loadtxt(data_path)
n_over_l = raw[:, 1]   # n/l

# regime label from directory name
ann_match = re.search(r"Ann([\d.]+)", outdir)
ann = float(ann_match.group(1)) if ann_match else 0.0
if ann > 0.05:
    regime = r"Supercritical ($a = -0.5$)"
elif ann > 0.003:
    regime = r"Critical ($a = -1$)"
else:
    regime = r"Subcritical ($a = -1.5$)"

colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(T_vals)))

fig, ax = plt.subplots(figsize=(7, 5))

for i, T in enumerate(T_vals):
    col_idx = 2 + 2 * i
    err_idx = 2 + 2 * i + 1
    G   = raw[:, col_idx]
    err = raw[:, err_idx]
    ax.plot(n_over_l, G,
            color=colors[i], lw=1.8, label=f"$T={T}$")
    ax.fill_between(n_over_l, G - err, G + err,
                    color=colors[i], alpha=0.20)

ax.axhline(0, color="k", lw=0.6, ls="--")
ax.set_xlabel(r"$n / l$", fontsize=13)
ax.set_ylabel(r"$G_n(T) = \langle\sigma_0(0)\,\sigma_n(T)\rangle$", fontsize=12)
ax.set_title(
    f"Time-delayed spin correlation — {regime}\n"
    r"$N=150$, Poissonian IC ($p=1/N$), Gillespie interface model",
    fontsize=10)
ax.legend(fontsize=10, loc="upper right")
ax.set_xlim(0, 0.5)

out_png = os.path.join(outdir, "time_delayed_corr.png")
fig.tight_layout()
fig.savefig(out_png, dpi=160)
print(f"Saved {out_png}")
