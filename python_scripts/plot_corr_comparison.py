"""
plot_corr_comparison.py  —  3-panel comparison of G_n(T) for sub/crit/super.
Usage:  python3 python_scripts/plot_corr_comparison.py
Saves:  corr_comparison.png in the Glauber root.
"""

import os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DIRS = [
    ("results_N150_T4.0_L4_Ann0.08_create3.63e-06_corr2pt", r"Supercritical  ($a = -0.5$)"),
    ("results_N150_T4.0_L4_Ann0.01_create2.96e-07_corr2pt", r"Critical  ($a = -1$)"),
    ("results_N150_T4.0_L4_Ann0.00_create2.42e-08_corr2pt", r"Subcritical  ($a = -1.5$)"),
]

def load(outdir):
    path = os.path.join(outdir, "corr2pt_multi.txt")
    with open(path) as f:
        header = f.readline().strip()
    tokens = header.lstrip("#").split("\t")
    T_vals = [float(tok[3:]) for tok in tokens if tok.startswith("G_T")]
    raw = np.loadtxt(path)
    return raw[:, 1], T_vals, raw  # n_over_l, T_vals, full_array

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False,
                         gridspec_kw={"wspace": 0.28})
fig.suptitle(
    r"Time-delayed spin correlation  $G_n(T)=\langle\sigma_0(0)\,\sigma_n(T)\rangle$"
    "\nPoissonian IC, $N=150$, Gillespie interface model",
    fontsize=12, y=1.02)

for ax, (ddir, title) in zip(axes, DIRS):
    n_over_l, T_vals, raw = load(ddir)
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(T_vals)))

    for i, T in enumerate(T_vals):
        G   = raw[:, 2 + 2 * i]
        err = raw[:, 2 + 2 * i + 1]
        ax.plot(n_over_l, G, color=colors[i], lw=1.8, label=f"$T={T}$")
        ax.fill_between(n_over_l, G - err, G + err,
                        color=colors[i], alpha=0.20)

    ax.axhline(0, color="k", lw=0.6, ls="--")
    ax.set_title(title, fontsize=12, pad=6)
    ax.set_xlabel(r"$n / l$", fontsize=12)
    ax.set_xlim(0, 0.5)
    ax.legend(fontsize=8.5, loc="upper right", framealpha=0.85)

axes[0].set_ylabel(r"$G_n(T)$", fontsize=12)

out = "corr_comparison.png"
fig.savefig(out, dpi=160, bbox_inches="tight")
print(f"Saved {out}")
