"""
plot_theory_vs_mc.py — overlay the finite-N Monte-Carlo correlators produced by
./Corr2pt with the exact theory of theory_glauber.py (docs/section1.tex).

Usage:  python python_scripts/plot_theory_vs_mc.py <results_dir>

The results directory must contain configuration.json plus, depending on the
run's corr_start:
  equilibrium :  corr2pt_multi.txt (G_n)  and  iface_corr2pt_multi.txt (H_n)
  disordered  :  corr_eqtime_multi.txt (C_n, single-time after the quench)

For each macroscopic time T the script reads N, nu, L, ann, create and
time_scale_formula (= 2 N^2) from configuration.json, evaluates the theory at
the Glauber time t = time_scale_formula * T, overlays it on the MC with its
error band, and prints a chi^2/dof per curve.
"""

import sys
import os
import json
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_glauber as th


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------

def load_config(results_dir):
    with open(os.path.join(results_dir, "configuration.json")) as f:
        return json.load(f)


def read_columns(path, sym):
    """Read a tidy correlator file.

    Returns (n, {T: (mean, stderr)}) for every column header  f"{sym}_T<value>".
    """
    with open(path) as f:
        header = f.readline().lstrip("#").strip()
    tokens = header.split("\t")
    raw = np.loadtxt(path)
    n = raw[:, 0].astype(int)
    pat = re.compile(r"^%s_T([0-9.]+)$" % re.escape(sym))
    out = {}
    for j, tok in enumerate(tokens):
        m = pat.match(tok.strip())
        if m:
            T = float(m.group(1))
            out[T] = (raw[:, j], raw[:, j + 1])
    return n, out


def chi2_per_dof(mc, err, theory):
    """Diagonal chi^2/dof over points with a positive error bar."""
    mask = err > 0
    if not np.any(mask):
        return float("nan")
    r = (mc[mask] - theory[mask]) / err[mask]
    return float(np.sum(r * r) / mask.sum())


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------

def panel_two_curve(ax, n, cols, N, nu, eta, g, l, time_scale, kind):
    """G-style panel: MC vs infinite-line and exact-ring theory.

    kind = 'G' (two-time spin) -> overlay G_two_time and G_two_time_ring.
    """
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, max(len(cols), 1)))
    for i, T in enumerate(sorted(cols)):
        mean, err = cols[T]
        t = time_scale * T
        line = th.G_two_time(n, t, eta, g)
        ring = th.G_two_time_ring(n, t, l, eta, g)
        c = colors[i]
        ax.errorbar(n, mean, yerr=err, fmt="o", ms=3, color=c, alpha=0.6,
                    capsize=0, lw=0.8, label=f"MC  T={T:g}")
        ax.plot(n, ring, "-", color=c, lw=1.6)
        ax.plot(n, line, "--", color=c, lw=1.0)
        chi_line = chi2_per_dof(mean, err, line)
        chi_ring = chi2_per_dof(mean, err, ring)
        print(f"  {kind}  T={T:<5g}  chi2/dof(ring)={chi_ring:8.3f}  "
              f"chi2/dof(line)={chi_line:8.3f}")
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$G_n(T)=\langle\sigma_0(0)\,\sigma_n(T)\rangle$")
    ax.set_title("Spin two-time correlator\n(— exact ring, -- infinite line)")
    ax.legend(fontsize=8, loc="upper right")


def panel_interface(ax, n, cols, N, nu, eta, g, l, p, time_scale):
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, max(len(cols), 1)))
    for i, T in enumerate(sorted(cols)):
        mean, err = cols[T]
        t = time_scale * T
        theory = th.H_interface(n, t, eta, g)
        c = colors[i]
        ax.errorbar(n, mean, yerr=err, fmt="o", ms=3, color=c, alpha=0.6,
                    capsize=0, lw=0.8, label=f"MC  T={T:g}")
        ax.plot(n, theory, "-", color=c, lw=1.6)
        chi = chi2_per_dof(mean, err, theory)
        print(f"  H  T={T:<5g}  chi2/dof={chi:8.3f}")
    ax.axhline(p * p, color="k", lw=0.7, ls=":", label=r"plateau $p^2$")
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$H_n(T)=\langle\eta_0(0)\,\eta_n(T)\rangle$")
    ax.set_title("Interface two-time correlator")
    ax.legend(fontsize=8, loc="upper right")


def panel_single(ax, n, cols, N, nu, eta, g, l, time_scale):
    colors = plt.cm.cividis(np.linspace(0.1, 0.85, max(len(cols), 1)))
    for i, T in enumerate(sorted(cols)):
        mean, err = cols[T]
        t = time_scale * T
        theory = th.C_single_disordered(n, t, eta, g)
        c = colors[i]
        ax.errorbar(n, mean, yerr=err, fmt="o", ms=3, color=c, alpha=0.6,
                    capsize=0, lw=0.8, label=f"MC  T={T:g}")
        ax.plot(n, theory, "-", color=c, lw=1.6)
        chi = chi2_per_dof(mean, err, theory)
        print(f"  C  T={T:<5g}  chi2/dof={chi:8.3f}")
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$C_n(T)=\langle\sigma_0(T)\,\sigma_n(T)\rangle$")
    ax.set_title("Single-time correlator (disordered quench)")
    ax.legend(fontsize=8, loc="upper right")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(results_dir):
    cfg = load_config(results_dir)
    N  = int(cfg["N"])
    nu = float(cfg.get("nu", 1.0))
    L  = int(cfg["L"])
    l  = N * L
    eta = float(cfg.get("eta", th.eta_of(N, nu)))
    g   = float(cfg.get("gamma", th.gamma_of(N, nu)))
    p   = float(cfg.get("p_equilibrium", th.p_of(N, nu)))
    time_scale = float(cfg.get("time_scale_formula", 2.0 * N * N))
    start = cfg.get("corr_start", "equilibrium")

    print(f"results: {results_dir}")
    print(f"N={N} nu={nu} L={L} l={l}  gamma={g:.6f} eta={eta:.6f} p={p:.6f}")
    print(f"corr_start={start}  time_scale=2N^2={time_scale:g}")

    panels = []
    g_path = os.path.join(results_dir, "corr2pt_multi.txt")
    h_path = os.path.join(results_dir, "iface_corr2pt_multi.txt")
    c_path = os.path.join(results_dir, "corr_eqtime_multi.txt")

    if start != "disordered":
        if os.path.exists(g_path):
            n, cols = read_columns(g_path, "G")
            if cols:
                panels.append(("G", n, cols))
        if os.path.exists(h_path):
            n, cols = read_columns(h_path, "H")
            if cols:
                panels.append(("H", n, cols))
    else:
        if os.path.exists(c_path):
            n, cols = read_columns(c_path, "C")
            if cols:
                panels.append(("C", n, cols))

    if not panels:
        print("No matching observable files found; nothing to plot.")
        return

    fig, axes = plt.subplots(1, len(panels), figsize=(6.2 * len(panels), 4.8),
                             squeeze=False)
    suptitle = (rf"$N={N}$, $\nu={nu:g}$, $\gamma={g:.4f}$, "
                rf"$\eta={eta:.4f}$, $p={p:.4f}$")
    fig.suptitle(suptitle, fontsize=12)

    for ax, (kind, n, cols) in zip(axes[0], panels):
        print(f"[{kind}]")
        if kind == "G":
            panel_two_curve(ax, n, cols, N, nu, eta, g, l, time_scale, "G")
        elif kind == "H":
            panel_interface(ax, n, cols, N, nu, eta, g, l, p, time_scale)
        elif kind == "C":
            panel_single(ax, n, cols, N, nu, eta, g, l, time_scale)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(results_dir, "theory_vs_mc.png")
    fig.savefig(out, dpi=160)
    print(f"Saved {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: plot_theory_vs_mc.py <results_dir>")
        sys.exit(1)
    main(sys.argv[1])
