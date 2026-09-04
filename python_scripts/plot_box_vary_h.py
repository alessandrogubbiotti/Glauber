"""
plot_box_vary_h.py — fixed N; ONE FIGURE PER box size h (not overlaid).

For each result directory (one box size h) it produces a separate figure with
K panels (one per macroscopic time tau), each showing the box two-time
correlator (points: MC, line: continuum) for that single h.

Spin:      continuum = box-average of H_torus(.,tau,L)              (e:Herfc/e:Htorus)
Interface: continuum = box-average of f(.,tau)  [corrected f]       (e:f)

The reference box is at x=L/2; N is fixed (in the title). Output files are
named <out_base>_h<h>.{png,pdf}, one per directory.

Usage: python3 plot_box_vary_h.py DIR1 [DIR2 ...] [--out base.png]
"""
import argparse, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_continuum as tc


def load(d):
    cfg = json.load(open(os.path.join(d, "configuration.json")))
    raw = np.atleast_2d(np.loadtxt(os.path.join(d, "box_corr.txt")))
    K = len(cfg["T_values"]); per = (raw.shape[1] - 1) // K
    x = raw[:, 0]
    mc = np.stack([raw[:, 1 + per * k] for k in range(K)])
    err = np.stack([raw[:, 2 + per * k] for k in range(K)])
    return cfg, x, mc, err


def continuum(obs, sep, T, h, L):
    if obs == "iface":
        # ring (length-L torus) continuum: image-summed f, NOT the bare line f.
        # The wrap-around roughly doubles the antipode value (see tc.f_torus).
        return tc.box_average(lambda xx, tt: tc.f_torus(xx, tt, L), sep, T, h)
    return tc.box_average(lambda xx, tt: tc.H_torus(xx, tt, L), sep, T, h)


def plot_one(cfg, x, mc, err, out_base):
    L, N = cfg["L"], cfg["N"]
    obs = cfg.get("observable", "spin")
    h = cfg["box_width_h"]
    T_values = cfg["T_values"]; K = len(T_values)
    nsims = cfg.get("nsims", "?")
    ylab = (r"$N^2\,\langle(\eta_{\rm box}{-}p)(\eta_{\rm box}{-}p)\rangle \to \bar f$"
            if obs == "iface" else
            r"$C(x,\tau)=\langle m_{\rm box}(L/2,0)\,m_{\rm box}(x,\tau)\rangle$")
    name = "interface memory $f$" if obs == "iface" else "spin correlator"

    sep = np.where((x - L / 2) % L > L / 2, (x - L / 2) % L - L, (x - L / 2) % L)
    xf = np.linspace(0, L, 500, endpoint=False)
    sepf = np.where((xf - L / 2) % L > L / 2, (xf - L / 2) % L - L, (xf - L / 2) % L)

    fig, axes = plt.subplots(1, K, figsize=(4.3 * K, 4.1), sharey=True)
    if K == 1:
        axes = [axes]
    chis = []
    for k, T in enumerate(T_values):
        cont = continuum(obs, sep, T, h, L)
        m = err[k] > 0
        chi = np.sum(((mc[k][m] - cont[m]) / err[k][m]) ** 2) / max(m.sum(), 1)
        chis.append(chi)
        ax = axes[k]
        ax.plot(xf, continuum(obs, sepf, T, h, L), "-", color="k", lw=1.6, zorder=5)
        ax.errorbar(x, mc[k], yerr=err[k], fmt="o", ms=3.6, color="#1f77b4",
                    capsize=1.6, lw=0.9, zorder=6)
        ax.axvline(L / 2, color="0.6", lw=0.7, ls="--")
        ax.set_title(rf"$\tau={T:g}$   ($\chi^2/\mathrm{{dof}}={chi:.2f}$)", fontsize=10)
        ax.set_xlabel(r"macroscopic position $x$")
    axes[0].set_ylabel(ylab)
    fig.suptitle(rf"Continuum limit of the box {name} at $N={N}$, $L={L}$, "
                 rf"$h={h:g}$ (equilibrium, {nsims} sims/pt); points: MC, line: continuum",
                 y=1.02)
    fig.tight_layout()
    base = f"{out_base}_h{h:g}"
    fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(base + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return h, obs, N, chis, base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", default="box")
    args = ap.parse_args()
    base = os.path.splitext(args.out)[0]
    for d in args.dirs:
        cfg, x, mc, err = load(d)
        h, obs, N, chis, outbase = plot_one(cfg, x, mc, err, base)
        print(f"{obs} N={N} h={h:g}: chi2/dof = "
              + " ".join(f"T={t:g}:{c:.2f}" for t, c in zip(cfg["T_values"], chis))
              + f"  -> {outbase}.png/.pdf")


if __name__ == "__main__":
    main()
