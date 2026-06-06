"""
plot_corr2pt.py  —  Plot two-point time-delayed correlation C(r, t1).

Usage (called automatically from corr2pt_main):
    python3 python_scripts/plot_corr2pt.py <outdir> <T_burnin> <t1> <x0>

Reads:  <outdir>/corr2pt.txt
Writes: <outdir>/corr2pt.png  (linear scale with error bars)
        <outdir>/corr2pt_log.png (log scale with exponential fit)
"""

import sys, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def exp_decay(r, A, xi):
    return A * np.exp(-r / xi)

def main():
    if len(sys.argv) < 5:
        print(__doc__); sys.exit(1)

    outdir   = sys.argv[1]
    T_burnin = float(sys.argv[2])
    t1       = float(sys.argv[3])
    x0       = int(sys.argv[4])

    data_path = os.path.join(outdir, "corr2pt.txt")
    data = np.loadtxt(data_path, comments="#")
    r       = data[:, 0]
    r_over_l = data[:, 1]
    C       = data[:, 2]
    sterr   = data[:, 3]

    # Attempt exponential fit on positive part
    pos = C > sterr * 0.5
    xi_fit = None
    A_fit  = None
    if pos.sum() >= 4:
        try:
            popt, pcov = curve_fit(exp_decay, r[pos], C[pos],
                                   p0=[C[0], r[-1] / 4],
                                   bounds=([0, 0.1], [2, r[-1]]))
            A_fit, xi_fit = popt
            perr = np.sqrt(np.diag(pcov))
            print(f"  Exp fit: A={A_fit:.4f}±{perr[0]:.4f}, "
                  f"xi={xi_fit:.2f}±{perr[1]:.2f} lattice units")
        except Exception as e:
            print(f"  Fit failed: {e}")

    # Try to load N and l from configuration.json
    cfg_path = os.path.join(outdir, "configuration.json")
    N, l = None, None
    if os.path.exists(cfg_path):
        with open(cfg_path) as f:
            cfg = json.load(f)
        N = cfg.get("N", "?")
        l = int(N) * int(cfg.get("L", 4)) if isinstance(N, int) else None

    title = f"C(r, t₁={t1:.2f}) after burn-in T={T_burnin:.2f},  x₀={x0}"
    if N: title += f",  N={N}"

    # ---- Linear scale ----
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(r_over_l, C, yerr=sterr, fmt="o-", ms=5, lw=1.4,
                capsize=3, label="Empirical  $C(r, t_1)$")
    if xi_fit and l:
        rf = np.linspace(0, r_over_l[-1], 400)
        ax.plot(rf, exp_decay(rf * l, A_fit, xi_fit), "--", lw=1.4,
                label=f"Exp. fit  $\\xi={xi_fit:.1f}$,  A={A_fit:.3f}")
    ax.axhline(0, color="gray", lw=0.8, ls=":")
    ax.set_xlabel("$r / l$")
    ax.set_ylabel("$C(r, t_1)$")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = os.path.join(outdir, "corr2pt.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved {out}")

    # ---- Log scale ----
    fig, ax = plt.subplots(figsize=(8, 5))
    mask = C > sterr   # only plot points well above noise
    if mask.any():
        ax.errorbar(r_over_l[mask], C[mask], yerr=sterr[mask],
                    fmt="o-", ms=5, lw=1.4, capsize=3,
                    label="Empirical  $|C(r, t_1)|$")
        if xi_fit and l:
            rf = np.linspace(r_over_l[mask][0], r_over_l[mask][-1], 400)
            ax.semilogy(rf, exp_decay(rf * l, A_fit, xi_fit), "--", lw=1.4,
                        label=f"Exp. fit  $\\xi={xi_fit:.1f}$")
    ax.set_yscale("log")
    ax.set_xlabel("$r / l$")
    ax.set_ylabel("$C(r, t_1)$  [log scale]")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = os.path.join(outdir, "corr2pt_log.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved {out}")

if __name__ == "__main__":
    main()
