"""Validate the interface (eta) box engine: MC vs exact finite-N (H_interface)
vs continuum box-averaged f and f/4.  Decides the prefactor of the interface
two-time correlator (the '1/4' flagged in the paper). Usage: validate_iface.py DIR"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_glauber as tg, theory_continuum as tc

d = sys.argv[1]
cfg = json.load(open(os.path.join(d, "configuration.json")))
N, L, l = cfg["N"], cfg["L"], cfg["l"]
eta, gamma, p = cfg["eta"], cfg["gamma"], cfg["p_equilibrium"]
h = cfg["box_width_h"]; w = (cfg["box_sites"] - 1) // 2; c0 = cfg["reference_site"]
T_values = cfg["T_values"]; K = len(T_values)
raw = np.atleast_2d(np.loadtxt(os.path.join(d, "box_corr.txt")))
x = raw[:, 0]
mc = np.stack([raw[:, 1 + 2 * k] for k in range(K)])     # N^2 * connected eta-box corr
err = np.stack([raw[:, 2 + 2 * k] for k in range(K)])
centers = np.round(x * N).astype(int) % l
sep = np.where((x - L / 2) % L > L / 2, (x - L / 2) % L - L, (x - L / 2) % L)


def exact_box(T):                       # N^2 * box-averaged connected H_interface (ring)
    t = 2.0 * N * N * T
    gall = tg.H_interface(np.arange(l // 2 + 1), t, eta, gamma) - p * p
    off = np.arange(-w, w + 1); sep2 = off[:, None] - off[None, :]
    out = np.empty(len(centers))
    for i, cj in enumerate(centers):
        dd = np.abs(((cj - c0 + sep2 + l // 2) % l) - l // 2)
        out[i] = N * N * gall[dd].mean()
    return out


xf = np.linspace(0, L, 400, endpoint=False)
sepf = np.where((xf - L / 2) % L > L / 2, (xf - L / 2) % L - L, (xf - L / 2) % L)
print(f"iface validation: N={N} l={l} h={h} |B|={2*w+1} nsims={cfg.get('nsims')}")
print(f"{'T':>6} {'chi2/dof vs exact':>18} {'vs f/4':>10} {'vs f':>10}")
fig, axes = plt.subplots(1, K, figsize=(4.3 * K, 4.0), sharey=True)
if K == 1:
    axes = [axes]
for k, T in enumerate(T_values):
    ex = exact_box(T)
    cf = tc.box_average(tc.f, sep, T, h)
    m = err[k] > 0
    chi_ex = np.sum(((mc[k][m] - ex[m]) / err[k][m]) ** 2) / m.sum()
    chi_f4 = np.sum(((mc[k][m] - cf[m] / 4) / err[k][m]) ** 2) / m.sum()
    chi_f = np.sum(((mc[k][m] - cf[m]) / err[k][m]) ** 2) / m.sum()
    print(f"{T:>6} {chi_ex:>18.2f} {chi_f4:>10.2f} {chi_f:>10.2f}")
    ax = axes[k]
    ax.errorbar(x, mc[k], yerr=err[k], fmt="o", ms=3, color="C0", capsize=1.5, label="MC")
    ax.plot(xf, tc.box_average(tc.f, sepf, T, h) / 4, "k-", lw=1.6, label="cont. f/4")
    ax.plot(xf, tc.box_average(tc.f, sepf, T, h), "r--", lw=1.0, alpha=0.6, label="cont. f (paper e:f)")
    ax.set_title(f"T={T:g}"); ax.set_xlabel("x")
axes[0].set_ylabel(r"$N^2\,\langle(\eta_{\rm box}-p)(\eta_{\rm box}-p)\rangle$")
axes[-1].legend(fontsize=8)
fig.suptitle(f"Interface memory: MC vs continuum f/4 (N={N}, h={h})", y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(d, "iface_validate.png"), dpi=160, bbox_inches="tight")
print("saved", os.path.join(d, "iface_validate.png"))
