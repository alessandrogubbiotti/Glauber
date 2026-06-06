"""
plot_gillespie.py  —  visualise output from the Gillespie interface simulation.

Usage:
    python3 python_scripts/plot_gillespie.py <results_dir> [<results_dir2> ...]

Per directory, produces:
    mean_spin.png          <sigma(x0,t)> vs macroscopic time
    corr_pair.png          C(x1,x2,t) = <sigma(x1,t)sigma(x2,t)> vs t
    spin_corr_static.png   time-averaged C(r) on log-y axis with exponential fit
    spin_corr_dynamics.png C(r,t) at multiple selected times (stacked curves)
    iface_corr.png         interface pair correlation g(r)
    timing.png             histogram of per-simulation wall times
    overview.png           4-panel summary

Multi-directory: comparison.png overlaying all dirs.
"""

import sys, os, json, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_json(path):
    with open(path) as f: return json.load(f)

def load_txt(path, skip_header=True):
    """Load a tab-separated file, skip comment/header lines starting with #."""
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1: data = data.reshape(1, -1)
    return data

def load_sel_t(path):
    """Load spin_corr_sel_t.txt: first col r/l, remaining cols are C(r,t_k)."""
    with open(path) as f:
        header = f.readline().strip()          # "# r/l  C_t=...  C_t=..."
    t_vals = [float(s.split("=")[1]) for s in header.split("\t")[1:]]
    data   = np.loadtxt(path, comments="#")
    r_over_l = data[:, 0]
    C_rt     = data[:, 1:]                    # shape (n_r, n_t)
    return r_over_l, C_rt, np.array(t_vals)

def load_timing(path):
    d = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            key, val = line.split("\t")
            try:    d[key] = float(val)
            except: d[key] = val
    return d

def exp_decay(r, xi):
    return np.exp(-r / xi)

def fit_xi(r_over_l, C, l):
    r = r_over_l * l
    mask = (r >= 1) & (C > 1e-10)
    if mask.sum() < 3: return None
    try:
        popt, _ = curve_fit(exp_decay, r[mask], C[mask], p0=[5.], bounds=(0.1, l))
        return popt[0]
    except: return None


# ---------------------------------------------------------------------------
# Per-directory analysis
# ---------------------------------------------------------------------------

def analyse(results_dir):
    cfg_path  = os.path.join(results_dir, "configuration.json")
    if not os.path.exists(cfg_path):
        print(f"  Skipping {results_dir}: no configuration.json"); return None

    cfg = load_json(cfg_path)
    N, L = int(cfg.get("N", 0)), int(cfg.get("L", 0))
    l = N * L
    ann = cfg.get("annihilation", "?")
    cre = cfg.get("creation", "?")
    label = f"N={N}, ann={ann:.3f}, cre={cre:.2e}"

    def p(name): return os.path.join(results_dir, name)
    def exists(name): return os.path.exists(p(name))

    result = {"label": label, "l": l, "cfg": cfg, "dir": results_dir}

    # --- mean spin evolution ---
    if exists("mean_spin_t.txt"):
        d = load_txt(p("mean_spin_t.txt"))
        result["mean_spin_t"] = d[:, 0]
        result["mean_spin"]   = d[:, 1]
        result["mean_spin_std"] = d[:, 2]

    # --- two-site correlation ---
    if exists("corr_pair_t.txt"):
        d = load_txt(p("corr_pair_t.txt"))
        result["corr_pair_t"]   = d[:, 0]
        result["corr_pair"]     = d[:, 1]
        result["corr_pair_std"] = d[:, 2]

    # --- time-averaged spin correlation ---
    if exists("spin_corr.txt"):
        d = load_txt(p("spin_corr.txt"))
        result["r_sc"] = d[:, 0]
        result["C_sc"] = d[:, 1]
        result["xi"]   = fit_xi(d[:, 0], d[:, 1], l)
        if result["xi"]:
            print(f"  [{os.path.basename(results_dir)}] "
                  f"xi = {result['xi']:.2f} lattice units")

    # --- time-dependent spin correlation C(r,t) ---
    if exists("spin_corr_sel_t.txt"):
        r_over_l, C_rt, t_vals = load_sel_t(p("spin_corr_sel_t.txt"))
        result["sel_r"] = r_over_l
        result["C_rt"]  = C_rt
        result["t_vals"] = t_vals

    # --- time autocorrelation ---
    if exists("time_autocorr.txt"):
        with open(p("time_autocorr.txt")) as f:
            header = f.readline()
        grand_mean = float(re.search(r"<s>=([-\d.e+]+)", header).group(1)) \
                     if re.search(r"<s>=([-\d.e+]+)", header) else 0.0
        d = load_txt(p("time_autocorr.txt"))
        result["autocorr_tau"] = d[:, 0]
        result["autocorr_A"]   = d[:, 1]   # raw
        result["autocorr_G"]   = d[:, 2]   # connected
        result["grand_mean"]   = grand_mean

    # --- interface pair correlation ---
    if exists("iface_corr.txt"):
        d = load_txt(p("iface_corr.txt"))
        result["r_ic"] = d[:, 0]
        result["g_ic"] = d[:, 1]

    # --- timing ---
    if exists("timing_stats.txt"):
        result["timing"] = load_timing(p("timing_stats.txt"))

    return result


# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------

CMAP = plt.get_cmap("plasma")

def plot_mean_spin(res, outdir):
    if "mean_spin" not in res: return
    t, m, s = res["mean_spin_t"], res["mean_spin"], res["mean_spin_std"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, m, lw=1.5, label=r"$\langle\sigma(x_0,t)\rangle$")
    ax.fill_between(t, m - s, m + s, alpha=0.25, label="±1 std")
    ax.axhline(0, color="gray", lw=0.7, ls="--")
    ax.set_xlabel("Macroscopic time $t$")
    ax.set_ylabel(r"$\langle\sigma(x_0,t)\rangle$")
    ax.set_title(f"Mean spin at $x_0$ — {res['label']}")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "mean_spin.png"), dpi=150); plt.close(fig)


def plot_corr_pair(res, outdir):
    if "corr_pair" not in res: return
    t, c, s = res["corr_pair_t"], res["corr_pair"], res["corr_pair_std"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, c, lw=1.5, color="tab:orange",
            label=r"$\langle\sigma(x_1,t)\sigma(x_2,t)\rangle$")
    ax.fill_between(t, c - s, c + s, alpha=0.25, color="tab:orange")
    ax.axhline(0, color="gray", lw=0.7, ls="--")
    ax.set_xlabel("Macroscopic time $t$")
    ax.set_ylabel(r"$C(x_1,x_2,t)$")
    ax.set_title(f"Two-site correlation — {res['label']}")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "corr_pair.png"), dpi=150); plt.close(fig)


def plot_spin_corr_static(res, outdir):
    if "r_sc" not in res: return
    r, C = res["r_sc"], res["C_sc"]
    xi = res.get("xi")
    fig, ax = plt.subplots(figsize=(8, 4))
    mask = C > 0
    ax.semilogy(r[mask], C[mask], "o-", ms=3, lw=1.2,
                label="Empirical $C(r)$ (time avg.)")
    if xi and res["l"]:
        l = res["l"]
        rf = np.linspace(r[mask][0], r[mask][-1], 300)
        ax.semilogy(rf, np.exp(-rf * l / xi), "--", lw=1.2,
                    label=f"Exp. fit  $\\xi={xi:.1f}$")
    ax.set_xlabel("$r/l$")
    ax.set_ylabel("$C(r)$ [log scale]")
    ax.set_title(f"Time-averaged spin-spin correlation — {res['label']}")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "spin_corr_static.png"), dpi=150)
    plt.close(fig)


def plot_spin_corr_dynamics(res, outdir):
    if "C_rt" not in res: return
    r, C_rt, t_vals = res["sel_r"], res["C_rt"], res["t_vals"]
    n_t = len(t_vals)
    colors = [CMAP(i / max(n_t - 1, 1)) for i in range(n_t)]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"Spin correlation $C(r,t)$ — {res['label']}", fontsize=11)

    ax = axes[0]
    for i, (ti, col) in enumerate(zip(t_vals, colors)):
        ax.plot(r, C_rt[:, i], color=col, lw=1.2, alpha=0.85,
                label=f"$t={ti:.2f}$" if i % max(1, n_t // 6) == 0 else "")
    ax.set_xlabel("$r/l$"); ax.set_ylabel("$C(r,t)$")
    ax.set_title("Linear scale"); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    ax = axes[1]
    for i, (ti, col) in enumerate(zip(t_vals, colors)):
        mask = C_rt[:, i] > 1e-10
        if mask.any():
            ax.semilogy(r[mask], C_rt[mask, i], color=col, lw=1.2, alpha=0.85)
    sm = plt.cm.ScalarMappable(cmap=CMAP,
                                norm=plt.Normalize(t_vals[0], t_vals[-1]))
    fig.colorbar(sm, ax=ax, label="Macroscopic time $t$")
    ax.set_xlabel("$r/l$"); ax.set_ylabel("$C(r,t)$ [log]")
    ax.set_title("Log scale"); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "spin_corr_dynamics.png"), dpi=150)
    plt.close(fig)


def plot_autocorr(res, outdir):
    if "autocorr_A" not in res: return
    tau, A, G = res["autocorr_tau"], res["autocorr_A"], res["autocorr_G"]
    gm = res.get("grand_mean", 0.0)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"Time autocorrelation $A(\\tau)$ — {res['label']}", fontsize=11)

    ax = axes[0]
    ax.plot(tau, A, lw=1.5, label=r"$A(\tau)=\langle\sigma(t)\sigma(t+\tau)\rangle$")
    ax.axhline(gm**2, color="gray", ls="--", lw=1,
               label=f"$\\langle\\sigma\\rangle^2={gm**2:.4f}$")
    ax.set_xlabel(r"$\tau$ (macroscopic time)")
    ax.set_ylabel(r"$A(\tau)$")
    ax.set_title("Raw autocorrelation")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    ax = axes[1]
    mask = np.abs(G) > 1e-12
    if mask.any():
        ax.semilogy(tau[mask], np.abs(G[mask]), lw=1.5, color="tab:orange",
                    label=r"$|G(\tau)|=|A(\tau)-\langle\sigma\rangle^2|$")
        # Fit exponential decay to connected part
        pos = G[mask] > 0
        if pos.sum() > 5:
            try:
                popt, _ = curve_fit(exp_decay, tau[mask][pos], G[mask][pos],
                                    p0=[1.0], bounds=(0.01, 100))
                tf = np.linspace(0, tau[mask][-1], 300)
                ax.semilogy(tf, exp_decay(tf, popt[0]), "--", lw=1.2,
                            label=f"Exp fit  $\\tau_c={popt[0]:.3f}$")
            except: pass
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$|G(\tau)|$ [log]")
    ax.set_title("Connected autocorrelation (log)")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "time_autocorr.png"), dpi=150); plt.close(fig)
    print(f"  Saved {os.path.join(outdir, 'time_autocorr.png')}")


def plot_iface_corr(res, outdir):
    if "r_ic" not in res: return
    r, g = res["r_ic"], res["g_ic"]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(r, g, "s-", ms=3, lw=1.2, color="tab:green", label="$g(r)$ empirical")
    far = np.median(g[max(0, len(g) - len(g) // 4):])
    ax.axhline(far, color="gray", ls="--", lw=1,
               label=f"Far-field ≈ {far:.3e} (Poisson baseline)")
    ax.set_xlabel("$r/l$"); ax.set_ylabel("Pair density")
    ax.set_title(f"Interface pair correlation — {res['label']}")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "iface_corr.png"), dpi=150); plt.close(fig)


def plot_timing(res, outdir):
    if "timing" not in res: return
    tim = res["timing"]
    fig, ax = plt.subplots(figsize=(6, 3))
    keys = ["mean_wall_time_s", "std_wall_time_s", "mean_n_events", "std_n_events"]
    labels = ["Mean wall time (s)", "Std wall time (s)",
              "Mean events", "Std events"]
    vals = [tim.get(k, 0) for k in keys]
    bars = ax.bar(labels, vals, color=["tab:blue","tab:cyan","tab:orange","tab:red"])
    ax.bar_label(bars, fmt="%.2f", fontsize=8)
    ax.set_title(f"Timing stats — {res['label']}")
    ax.set_ylabel("Value"); ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "timing.png"), dpi=150); plt.close(fig)


def plot_overview(res, outdir):
    """4-panel summary."""
    has_ms  = "mean_spin"  in res
    has_cp  = "corr_pair"  in res
    has_sc  = "r_sc"       in res
    has_dyn = "C_rt"       in res
    if not any([has_ms, has_cp, has_sc, has_dyn]): return

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(f"Overview — {res['label']}", fontsize=12)

    # Panel 1: mean spin
    ax = axes[0, 0]
    if has_ms:
        t, m, s = res["mean_spin_t"], res["mean_spin"], res["mean_spin_std"]
        ax.plot(t, m, lw=1.5); ax.fill_between(t, m-s, m+s, alpha=0.25)
        ax.axhline(0, color="gray", lw=0.7, ls="--")
        ax.set_xlabel("$t$"); ax.set_ylabel(r"$\langle\sigma(x_0,t)\rangle$")
        ax.set_title("Mean spin at $x_0$"); ax.grid(True, alpha=0.3)

    # Panel 2: pair correlation
    ax = axes[0, 1]
    if has_cp:
        t, c, s = res["corr_pair_t"], res["corr_pair"], res["corr_pair_std"]
        ax.plot(t, c, lw=1.5, color="tab:orange")
        ax.fill_between(t, c-s, c+s, alpha=0.25, color="tab:orange")
        ax.axhline(0, color="gray", lw=0.7, ls="--")
        ax.set_xlabel("$t$"); ax.set_ylabel(r"$C(x_1,x_2,t)$")
        ax.set_title("Two-site correlation"); ax.grid(True, alpha=0.3)

    # Panel 3: static spin correlation (log scale)
    ax = axes[1, 0]
    if has_sc:
        r, C = res["r_sc"], res["C_sc"]; xi = res.get("xi")
        mask = C > 0
        ax.semilogy(r[mask], C[mask], "o-", ms=2, lw=1)
        if xi and res["l"]:
            l = res["l"]; rf = np.linspace(r[mask][0], r[mask][-1], 300)
            ax.semilogy(rf, np.exp(-rf*l/xi), "--", lw=1, label=f"$\\xi={xi:.1f}$")
            ax.legend(fontsize=8)
        ax.set_xlabel("$r/l$"); ax.set_ylabel("$C(r)$  [log]")
        ax.set_title("Time-avg spin correlation"); ax.grid(True, alpha=0.3)

    # Panel 4: C(r,t) dynamics
    ax = axes[1, 1]
    if has_dyn:
        r, C_rt, t_vals = res["sel_r"], res["C_rt"], res["t_vals"]
        n_t = len(t_vals)
        for i, ti in enumerate(t_vals):
            col = CMAP(i / max(n_t-1, 1))
            ax.plot(r, C_rt[:, i], color=col, lw=1, alpha=0.8)
        sm = plt.cm.ScalarMappable(cmap=CMAP,
                                    norm=plt.Normalize(t_vals[0], t_vals[-1]))
        fig.colorbar(sm, ax=ax, label="$t$", shrink=0.8)
        ax.set_xlabel("$r/l$"); ax.set_ylabel("$C(r,t)$")
        ax.set_title("Spin correlation vs time"); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "overview.png"), dpi=150); plt.close(fig)
    print(f"  Saved {os.path.join(outdir, 'overview.png')}")


def process_dir(results_dir):
    print(f"Processing {results_dir} ...")
    res = analyse(results_dir)
    if res is None: return None
    out = results_dir
    plot_mean_spin(res, out)
    plot_corr_pair(res, out)
    plot_spin_corr_static(res, out)
    plot_spin_corr_dynamics(res, out)
    plot_autocorr(res, out)
    plot_iface_corr(res, out)
    plot_timing(res, out)
    plot_overview(res, out)
    return res


# ---------------------------------------------------------------------------
# Multi-directory comparison
# ---------------------------------------------------------------------------

def plot_comparison(results):
    results = [r for r in results if r]
    if len(results) < 2: return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Regime comparison — Gillespie", fontsize=13)

    for res in results:
        lbl = res["label"]
        if "r_sc" in res:
            axes[0,0].semilogy(res["r_sc"], np.clip(res["C_sc"], 1e-12, None),
                               "o-", ms=2, lw=1, label=lbl)
        if "mean_spin" in res:
            axes[0,1].plot(res["mean_spin_t"], res["mean_spin"], lw=1, label=lbl)
        if "corr_pair" in res:
            axes[1,0].plot(res["corr_pair_t"], res["corr_pair"], lw=1, label=lbl)
        if "r_ic" in res:
            axes[1,1].plot(res["r_ic"], res["g_ic"], "-", lw=1, label=lbl)

    axes[0,0].set_title("C(r) time-avg [log]")
    axes[0,0].set_xlabel("$r/l$"); axes[0,0].legend(fontsize=7)
    axes[0,0].grid(True, alpha=0.3)

    axes[0,1].set_title(r"$\langle\sigma(x_0,t)\rangle$")
    axes[0,1].set_xlabel("$t$"); axes[0,1].legend(fontsize=7)
    axes[0,1].grid(True, alpha=0.3)

    axes[1,0].set_title(r"$C(x_1,x_2,t)$")
    axes[1,0].set_xlabel("$t$"); axes[1,0].legend(fontsize=7)
    axes[1,0].grid(True, alpha=0.3)

    axes[1,1].set_title("Interface pair correlation $g(r)$")
    axes[1,1].set_xlabel("$r/l$"); axes[1,1].legend(fontsize=7)
    axes[1,1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig("comparison.png", dpi=150); plt.close(fig)
    print("Saved comparison.png")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    dirs = sys.argv[1:]
    results = [process_dir(d) for d in dirs]
    if len(results) > 1:
        plot_comparison(results)
