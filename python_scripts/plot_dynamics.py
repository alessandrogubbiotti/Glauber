"""
plot_dynamics.py  —  space-time kymograph for subcritical / critical / supercritical.

Usage:
    python3 python_scripts/plot_dynamics.py

Reads trajectory/spins_output.bin from the three N=150 regime directories
and saves dynamics_comparison.png in the Glauber root.
"""

import os, json
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

MARKER = -9999

def load_traj(traj_dir):
    cfg_path   = os.path.join(traj_dir, "configuration.json")
    spins_path = os.path.join(traj_dir, "spins_output.bin")
    with open(cfg_path) as f:
        cfg = json.load(f)
    T, N, L, res = cfg["T"], cfg["N"], cfg["L"], cfg["resolution"]
    l = N * L
    Mac_T = int(T * res)
    raw = np.fromfile(spins_path, dtype=np.int32)
    expected = Mac_T * (l + 1)
    if len(raw) != expected:
        raise ValueError(f"{traj_dir}: expected {expected} int32s, got {len(raw)}")
    # reshape: (Mac_T, l+1) then drop the marker column
    frames = raw.reshape(Mac_T, l + 1)[:, :l]   # shape (Mac_T, l)
    return frames, T, L, N, res


cmap = ListedColormap(["#2166ac", "#d73027"])   # blue = −1, red = +1

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True,
                         gridspec_kw={"wspace": 0.04})
fig.suptitle("Interface dynamics — Glauber model, $N=150$", fontsize=14, y=1.01)

for ax, (results_dir, title) in zip(axes, DIRS):
    traj_dir = os.path.join(results_dir, "trajectory")
    frames, T, L, N, res = load_traj(traj_dir)     # (Mac_T, l)

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

# colorbar to the right of last axis only
cbar = fig.colorbar(im, ax=axes[-1], orientation="vertical",
                    fraction=0.046, pad=0.04, ticks=[-1, 1], shrink=0.85)
cbar.ax.set_yticklabels(["$\sigma = -1$", "$\sigma = +1$"], fontsize=10)

out = "dynamics_comparison.png"
fig.savefig(out, dpi=180, bbox_inches="tight")
print(f"Saved {out}")
