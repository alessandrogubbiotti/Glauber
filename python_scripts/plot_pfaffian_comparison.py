"""plot_pfaffian_comparison.py — the estimate of ./Pfaffian against the Pfaffian formulas.

Usage:  python3 python_scripts/plot_pfaffian_comparison.py <dir> [OUT.png]

<dir> is the output directory of ./Pfaffian.  For the four entries of the
rescaled kernel, for the density rho_1, for the two-point and the three-point
functions at one time and for the two-point function at two times, the script
compares three things:

  the estimate    the average over the simulations, with its standard error;
  the lattice     the exact value at the same N, the Pfaffian of the kernel of
                  Garrod, averaged over the same bins of the estimate;
  the limit       the Pfaffian of the limit kernel.

The last two are computed by pfaffian_kernel.py.  The script prints, for each
observable, the chi^2 per point of the estimate against the lattice and the
largest distance between the lattice and the limit, and draws the picture
(default <dir>/pfaffian_comparison.png).  The table has all the times; the
picture has at most six of them, and for the two-time function the pairs that
start from the first time.  The chi^2 is the diagonal one: the
bins of one observable are estimated on the same simulations, so they are not
independent.  rho_2, rho_3, K_22 and the two-time function are counts of pairs
and of triples: with few of them in a bin (some tens) the standard error of the
bin is not reliable and the chi^2 neither.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pfaffian_kernel as pk

# One colour for each time, always in this order.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
ENTRIES = ["K_11", "K_12", "K_21", "K_22"]


def chi2(mc, err, exact):
    """The sum of ((mc - exact)/err)^2 and the number of points with a positive error."""
    mask = err > 0
    r = (mc[mask] - exact[mask]) / err[mask]
    return float(np.sum(r * r)), int(mask.sum())


def kernel_on_bins(cfg, tau, limit):
    """The entries of the kernel and rho_2 averaged over the bins: a dictionary of arrays."""
    N, nu, nup, b, nb = cfg["N"], cfg["nu"], cfg["nu_prime"], cfg["bin"], cfg["n_bins"]
    out = {name: np.zeros(nb) for name in ENTRIES + ["rho_2"]}
    rho1 = pk.density(tau, nu, nup) if limit else pk.lattice_density(tau, N, nu, nup)
    for j in range(nb):
        for r in range(j * b + 1, (j + 1) * b + 1):
            K = pk.kernel_block(r / N, tau, nu, nup) if limit else pk.lattice_block(r, tau, N, nu, nup)
            for name, value in zip(ENTRIES, K.ravel()):
                out[name][j] += value / b
            out["rho_2"][j] += (rho1 * rho1 - K[0, 0] * K[1, 1] + K[0, 1] * K[1, 0]) / b
    return out


def rho3_on_bins(cfg, tau, limit):
    """rho_3 averaged over the pairs of bins of the two gaps: an array nb3 x nb3."""
    N, nu, nup, b, nb3 = cfg["N"], cfg["nu"], cfg["nu_prime"], cfg["bin"], cfg["n_bins_rho3"]
    out = np.zeros((nb3, nb3))
    for i in range(nb3):
        for j in range(nb3):
            for r1 in range(i * b + 1, (i + 1) * b + 1):
                for r2 in range(j * b + 1, (j + 1) * b + 1):
                    if limit:
                        value = pk.correlation([(tau, 0.0), (tau, r1 / N), (tau, (r1 + r2) / N)], nu, nup)
                    else:
                        value = pk.lattice_correlation([(tau, 0), (tau, r1), (tau, r1 + r2)], N, nu, nup)
                    out[i, j] += value / (b * b)
    return out


def two_times_on_bins(cfg, sigma, tau, limit):
    """rho(sigma, 0; tau, r): the value at r = 0 and then the averages over the bins of |r|."""
    N, nu, nup, b, nb = cfg["N"], cfg["nu"], cfg["nu_prime"], cfg["bin"], cfg["n_bins"]

    def rho(r):
        if limit:
            return pk.correlation([(sigma, 0.0), (tau, r / N)], nu, nup)
        return pk.lattice_correlation([(sigma, 0), (tau, r)], N, nu, nup)

    out = np.zeros(nb + 1)
    out[0] = rho(0)
    for j in range(nb):
        for r in range(j * b + 1, (j + 1) * b + 1):
            out[j + 1] += 0.5 * (rho(r) + rho(-r)) / b
    return out


def style(ax, title, xlabel):
    ax.set_title(title, fontsize=11, color=INK, loc="left")
    ax.set_xlabel(xlabel, fontsize=10, color=MUTED)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def draw(ax, x, mc, err, lattice, x_fine, limit, color):
    """Points for the estimate, a dashed line for the lattice, a full line for the limit."""
    ax.plot(x_fine, limit, "-", color=color, lw=1.8, zorder=2)
    ax.plot(x, lattice, "--", color=INK, lw=0.9, zorder=3)
    ax.errorbar(x, mc, yerr=err, fmt="o", ms=3.8, color=color, mec="white", mew=0.5,
                elinewidth=0.9, capsize=0, zorder=4)


def main(results_dir, out=None):
    with open(os.path.join(results_dir, "configuration.json")) as f:
        cfg = json.load(f)
    N, nu, nup, times = cfg["N"], cfg["nu"], cfg["nu_prime"], cfg["T_values"]
    out = out or os.path.join(results_dir, "pfaffian_comparison.png")
    print(f"N = {N}, nu = {nu:g}, nu' = {nup:g}, {cfg['nsims']} simulations, bins of {cfg['bin']} bonds")
    shown = times
    if len(times) > 6:                       # six times, spread from the first to the last
        shown = [times[i] for i in sorted(set(np.round(np.linspace(0, len(times) - 1, 6)).astype(int)))]
        print("the picture shows the times", ", ".join(f"{tau:g}" for tau in shown))
    colour = {tau: COLORS[i] for i, tau in enumerate(shown)}
    print(f"{'observable':<28}{'chi^2/points (estimate vs lattice)':>36}{'max |lattice - limit|':>26}")

    def report(name, mc, err, lattice, limit):
        c, n = chi2(mc.ravel(), err.ravel(), lattice.ravel())
        print(f"{name:<28}{c / max(n, 1):>24.2f}  ({n:4d} points){np.max(np.abs(lattice - limit)):>22.4f}")

    fig, axes = plt.subplots(2, 4, figsize=(17, 8), gridspec_kw={"wspace": 0.26, "hspace": 0.36})
    fig.suptitle(rf"The Pfaffian kernel and the correlation functions of the interfaces:  $N={N}$,  "
                 rf"$\nu={nu:g}$,  $\nu'={nup:g}$,  {cfg['nsims']} simulations", fontsize=13, color=INK, y=0.97)

    # ---- the four entries of the kernel and rho_2, at each time ----
    titles = {"K_11": r"$K_{11}$", "K_12": r"$K_{12}$", "K_21": r"$K_{21}$", "K_22": r"$K_{22}$"}
    for tau in times:
        data = np.loadtxt(os.path.join(results_dir, f"kernel_T{tau:g}.txt"), ndmin=2)
        d = data[:, 2]
        lattice, limit = kernel_on_bins(cfg, tau, False), kernel_on_bins(cfg, tau, True)
        d_fine = np.linspace(d[0] / 4, d[-1], 300)
        fine = {name: [] for name in ENTRIES}
        for x in d_fine if tau in colour else []:
            for name, value in zip(ENTRIES, pk.kernel_block(x, tau, nu, nup).ravel()):
                fine[name].append(value)
        fine["rho_2"] = pk.pair_density(d_fine, tau, nu, nup)
        for q, name in enumerate(ENTRIES + ["rho_2"]):
            mc, err = data[:, 3 + 2 * q], data[:, 4 + 2 * q]
            report(f"{name}  tau={tau:g}", mc, err, lattice[name], limit[name])
            if tau in colour:
                draw(axes[0, q] if q < 4 else axes[1, 1], d, mc, err, lattice[name], d_fine, fine[name], colour[tau])
    for q, name in enumerate(ENTRIES):
        style(axes[0, q], titles[name] + r"$(0,d)$, rescaled", r"$d$")
    style(axes[1, 1], r"$\rho_2(0,d)$ at one time", r"$d$")

    # ---- the density ----
    ax = axes[1, 0]
    data = np.loadtxt(os.path.join(results_dir, "density.txt"), ndmin=2)
    tau_grid = data[:, 0]
    lattice = np.array([pk.lattice_density(x, N, nu, nup) for x in tau_grid])
    limit = np.array([pk.density(x, nu, nup) for x in tau_grid])
    report("rho_1 on the grid of times", data[:, 1], data[:, 2], lattice, limit)
    tau_fine = np.linspace(0.0, tau_grid[-1], 400)
    ax.plot(tau_fine, [pk.density(x, nu, nup) for x in tau_fine], "-", color=MUTED, lw=1.8, zorder=2)
    ax.plot(tau_fine, [pk.lattice_density(x, N, nu, nup) for x in tau_fine], "--", color=INK, lw=0.9, zorder=3)
    ax.errorbar(tau_grid, data[:, 1], yerr=data[:, 2], fmt="o", ms=3.8, color=MUTED, mec="white",
                mew=0.5, elinewidth=0.9, zorder=4)
    ax.axhline(nu, color=MUTED, lw=0.8, ls=":")
    ax.annotate(r"$\nu$", (0.0, nu), textcoords="offset points", xytext=(6, 4), fontsize=10, color=MUTED)
    style(ax, r"The density $\rho_1(\tau)$", r"$\tau$")

    # ---- rho_3 at one time: the picture shows the equal gaps ----
    ax = axes[1, 2]
    nb3 = cfg["n_bins_rho3"]
    for tau in times if nb3 else []:
        data = np.loadtxt(os.path.join(results_dir, f"rho3_T{tau:g}.txt"), ndmin=2)
        mc, err = data[:, 4].reshape(nb3, nb3), data[:, 5].reshape(nb3, nb3)
        lattice, limit = rho3_on_bins(cfg, tau, False), rho3_on_bins(cfg, tau, True)
        report(f"rho_3  tau={tau:g}", mc, err, lattice, limit)
        if tau in colour:
            d = data[:, 2].reshape(nb3, nb3)[:, 0]
            d_fine = np.linspace(d[0] / 4, d[-1], 120)
            fine = [pk.correlation([(tau, 0.0), (tau, x), (tau, 2 * x)], nu, nup) for x in d_fine]
            draw(ax, d, np.diag(mc), np.diag(err), np.diag(lattice), d_fine, fine, colour[tau])
    if nb3:
        style(ax, r"$\rho_3(0,d,2d)$ at one time", r"$d$")
    else:                                    # one bin only: no room for the two gaps
        ax.set_axis_off()

    # ---- two times: all the pairs are in the table, the picture shows the ones from the first time ----
    ax = axes[1, 3]
    for a in range(len(times)):
        for b in range(a + 1, len(times)):
            sigma, tau = times[a], times[b]
            data = np.loadtxt(os.path.join(results_dir, f"twotime_S{sigma:g}_T{tau:g}.txt"), ndmin=2)
            d, mc, err = data[:, 2], data[:, 3], data[:, 4]
            lattice, limit = two_times_on_bins(cfg, sigma, tau, False), two_times_on_bins(cfg, sigma, tau, True)
            report(f"two times  {sigma:g}, {tau:g}", mc, err, lattice, limit)
            if a == 0 and tau in colour:
                d_fine = np.linspace(0.0, d[-1], 300)
                fine = [pk.correlation([(sigma, 0.0), (tau, x)], nu, nup) for x in d_fine]
                draw(ax, d, mc, err, lattice, d_fine, fine, colour[tau])
    if len(times) > 1:
        style(ax, rf"$\rho(\sigma,0;\tau,d)$ at two times, $\sigma={times[0]:g}$", r"$d$")
    else:
        ax.set_axis_off()

    # ---- the legends, once for the whole picture: the three marks and the colours of the times ----
    marks = [plt.Line2D([], [], marker="o", ls="", ms=5, color=MUTED, mec="white", label="estimate from the simulation"),
             plt.Line2D([], [], ls="--", lw=0.9, color=INK, label=f"lattice, exact at $N={N}$"),
             plt.Line2D([], [], ls="-", lw=1.8, color=MUTED, label="limit")]
    fig.legend(handles=marks, loc="lower center", ncol=3, fontsize=10, frameon=False, labelcolor=INK,
               bbox_to_anchor=(0.27, 0.0))
    colours = [plt.Line2D([], [], marker="o", ls="-", lw=1.8, ms=5, color=colour[tau], mec="white",
                          label=rf"$\tau={tau:g}$") for tau in shown]
    fig.legend(handles=colours, loc="lower center", ncol=len(shown), fontsize=10, frameon=False,
               labelcolor=INK, bbox_to_anchor=(0.72, 0.0))

    fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
    print(f"Saved {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:3])
