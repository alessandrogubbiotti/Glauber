"""plot_voter_comparison.py — the noisy voter model next to the Glauber dynamics.

Usage:  python3 python_scripts/plot_voter_comparison.py <output.txt> <results> [OUT.png]

<output.txt> is the file written by nvm_sim, the program of the submodule
Continuum-Noisy-Voter-Model, and <results> is a results directory of Gillespie,
whose trajectory/ is read.  The two space-time pictures are drawn side by side
with the same colours.

The time of the noisy voter model is divided by 2: its genealogies have variance
1 in a unit of time and the ones of the Glauber dynamics have variance 2.  With
lambda = 2 nu^2 the two dynamics have then the same limit.  The default name of
the picture is voter_comparison.png.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

MARKER = -9999


def load_voter(path):
    """The colours of nvm_sim as spins, with the length of the ring, the final time and lambda."""
    with open(path) as f:
        header = dict(item.split("=") for item in f.readline().lstrip("# ").split())
    spins = 2 * np.loadtxt(path, ndmin=2) - 1
    last_time = (spins.shape[0] - 1) / float(header["resolution"])
    return spins, float(header["L"]), last_time, float(header["lambda"])


def load_glauber(results_dir):
    """The trajectory saved by Gillespie, with the length of the ring, the final time and nu."""
    trajectory = os.path.join(results_dir, "trajectory")
    with open(os.path.join(trajectory, "configuration.json")) as f:
        cfg = json.load(f)
    l, frames = cfg["N"] * cfg["L"], int(cfg["T"] * cfg["resolution"])
    raw = np.fromfile(os.path.join(trajectory, "spins_output.bin"), dtype=np.int32)
    if len(raw) != frames * (l + 1) or np.any(raw.reshape(frames, l + 1)[:, l] != MARKER):
        raise ValueError(f"{trajectory}: spins_output.bin is not the one of configuration.json")
    return raw.reshape(frames, l + 1)[:, :l], cfg["L"], cfg["T"], cfg.get("nu")


def main(voter_path, results_dir, out="voter_comparison.png"):
    voter, L_voter, T_voter, lam = load_voter(voter_path)
    glauber, L_glauber, T_glauber, nu = load_glauber(results_dir)

    cmap = ListedColormap(["#2166ac", "#d73027"])   # blue = -1, red = +1
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, gridspec_kw={"wspace": 0.1})
    panels = ((voter, L_voter, T_voter / 2.0, rf"Noisy voter model, $\lambda = {lam:g}$"),
              (glauber, L_glauber, T_glauber, "Glauber dynamics" + (rf", $\nu = {nu:g}$" if nu else "")))
    for ax, (spins, length, last_time, title) in zip(axes, panels):
        im = ax.imshow(spins, cmap=cmap, vmin=-1, vmax=1, aspect="auto", origin="lower",
                       interpolation="nearest", extent=[0, length, 0, last_time])
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Space $x$", fontsize=11)
    axes[0].set_ylabel("Macroscopic time of the Glauber dynamics", fontsize=11)

    cbar = fig.colorbar(im, ax=axes[-1], orientation="vertical", fraction=0.046, pad=0.04,
                        ticks=[-1, 1], shrink=0.85)
    cbar.ax.set_yticklabels([r"$\sigma = -1$", r"$\sigma = +1$"], fontsize=10)

    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"Saved {out}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
