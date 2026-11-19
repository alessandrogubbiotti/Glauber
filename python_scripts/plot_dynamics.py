"""
plot_dynamics.py  —  space-time kymograph for subcritical / critical / supercritical.

Usage:
    python3 python_scripts/plot_dynamics.py [SUPER CRIT SUB [OUT.png]]

Reads trajectory/spins_output.bin from three results directories of Gillespie,
one per regime (a = -0.5, a = -1, a = -1.5, in this order), and saves the three
pictures side by side.  Without arguments the directories are the three
N=150, T=4 ones below and the figure is dynamics_comparison.png in the Glauber
root.
"""

import os, json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

DIRS = [
    ("results_N150_T4.0_L4_Ann0.08_create3.63e-06_gillespie", "Supercritical  ($a = -0.5$)"),
    ("results_N150_T4.0_L4_Ann0.01_create2.96e-07_gillespie", "Critical  ($a = -1$)"),
    ("results_N150_T4.0_L4_Ann0.00_create2.42e-08_gillespie", "Subcritical  ($a = -1.5$)"),
]
out = "dynamics_comparison.png"

# The three directories, and the name of the figure, can be given on the command line.
if len(sys.argv) >= 4:
    DIRS = [(d, title) for d, (_, title) in zip(sys.argv[1:4], DIRS)]
if len(sys.argv) >= 5:
    out = sys.argv[4]

MARKER = -9999

def load_traj(traj_dir):
    cfg_path   = os.path.join(traj_dir, "configuration.json")
    spins_path = os.path.join(traj_dir, "spins_output.bin")
    with open(cfg_path) as f:
        cfg = json.load(f)
    T, N, L, res = cfg["T"], cfg["N"], cfg["L"], cfg["resolution"]
    l = N * L
    # The frames are the ones of the times j/res not after T; the factor is for the
    # rounding of T, as in gillespie_main.c: 0.29 * 100 is 28.999999999999996.
    Mac_T = int(np.floor(T * res * (1.0 + 1e-12)))
    raw = np.fromfile(spins_path, dtype=np.int32)
    expected = Mac_T * (l + 1)
    if len(raw) != expected:
        raise ValueError(f"{traj_dir}: expected {expected} int32s, got {len(raw)}")
    # reshape: (Mac_T, l+1) then drop the marker column
    frames = raw.reshape(Mac_T, l + 1)[:, :l]   # shape (Mac_T, l)
    return frames, T, L, N, res


cmap = ListedColormap(["#2166ac", "#d73027"])   # blue = −1, red = +1

# The gap leaves room for the labels 4.0 and 0.0 of two neighbouring panels.
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True,
                         gridspec_kw={"wspace": 0.1})

for ax, (results_dir, title) in zip(axes, DIRS):
    traj_dir = os.path.join(results_dir, "trajectory")
    frames, T, L, N, res = load_traj(traj_dir)     # (Mac_T, l)
    fig.suptitle(f"Interface dynamics — Glauber model, $N={N}$", fontsize=14, y=1.01)

    im = ax.imshow(
        frames,
        cmap=cmap, vmin=-1, vmax=1,
        aspect="auto",
        extent=[0, L, 0, T],
        origin="lower",
        interpolation="nearest",
    )
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Space $x$", fontsize=11)

axes[0].set_ylabel("Macroscopic time $t$", fontsize=11)

# The colorbar takes its room from the three panels, so that they have the same width.
cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                    fraction=0.015, pad=0.013, ticks=[-1, 1], shrink=0.85)
cbar.ax.set_yticklabels([r"$\sigma = -1$", r"$\sigma = +1$"], fontsize=10)

fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved {out}")
