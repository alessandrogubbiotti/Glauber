"""
mc_box_correlation.py — Monte-Carlo vs exact box-magnetization correlation.

Fix a number of macroscopic time points T_1..T_K, a number of macroscopic space
points x_1..x_M, and a small box width h.  The local (box) magnetization at a
macroscopic point x is the average spin over a box of width h centred at x,

    m_box(x, t) = (1/|B|) sum_{i in B(x)} sigma_i(t),   |B| = 2w+1,  w = round(hN/2).

We estimate, by Monte Carlo from the reversible (equilibrium) measure, the
two-time correlation of the box magnetization at the centre L/2 and time 0 with
the box magnetization at every x_j and time T_k,

    C(x_j, T_k) = < m_box(L/2, 0) m_box(x_j, T_k) > ,

and overlay the exact solution

    C_th(x_j, T_k) = (1/|B|^2) sum_{i in B(L/2)} sum_{p in B(x_j)} G_ring(p-i, 2 N^2 T_k),

where G_ring(n, t) = <sigma_0(0) sigma_n(t)>_eq is the exact finite-ring two-time
spin correlator of theory_glauber.py (Glauber time t = 2 N^2 tau).

Sampling modes (--mode):
  independent (default): a SEPARATE independent ensemble of `--nsims`
      realizations for EVERY (x_j, T_k) grid point — fresh reference draw, fresh
      trajectories, independent RNG stream.  The estimates at different grid
      points are then statistically independent, so the per-point error bars are
      independent and chi^2/dof per time concentrates around 1.  Cost scales as
      M * K * nsims trajectories.
  joint: ONE ensemble of `--nsims` realizations measured at every (x_j, T_k)
      jointly (M*K times cheaper).  Still unbiased, but the points at a fixed T
      share the same reference and trajectories, so they are correlated (the
      whole curve shifts together) and the diagonal chi^2/dof is not meaningful.

Equilibrium IC (matches src/spin_init.c): number of interfaces k ~ Binomial(l,p)
conditioned even, placed on a uniform k-subset of bonds (exclusion), integrated
to spins.  Dynamics: spin-flip BKL at the closing "glauber" rates in code units
(hop rate 1/dir); macroscopic time tau = micro/N^2.

Usage:
    python python_scripts/mc_box_correlation.py [--N 12] [--L 4] [--nu 1.0]
        [--h 0.5] [--M 13] [--T 0.05,0.1,0.2,0.4] [--nsims 800]
        [--mode independent|joint] [--seed 0] [--out box_corr]
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theory_glauber as th


# ---------------------------------------------------------------------------
# Equilibrium initial condition (particle sampling; matches src/spin_init.c)
# ---------------------------------------------------------------------------

def equilibrium_ic(l, p, rng):
    """Reversible (Ising-ring) measure, sampled the "particle" way (Prop 1.2):

      1. number of interfaces k ~ Binomial(l, p) CONDITIONED on k even
         (the exact finite-l law; Poisson is only the N->inf limit);
      2. place the k interfaces on a uniformly random k-subset of the l bonds
         (distinct bonds = exclusion rule);
      3. integrate bonds -> spins with sigma_0 = +/-1 uniform.

    Exact conditioned-even Bernoulli(p) field (pair correlator C_eq_ring)."""
    while True:
        k = int(rng.binomial(l, p))
        if (k & 1) == 0:
            break
    bond = np.zeros(l, dtype=np.int64)
    if k:
        bond[rng.choice(l, size=k, replace=False)] = 1     # uniform k-subset
    prefix = np.concatenate(([0], np.cumsum(bond[:-1]))) & 1
    s0 = 1.0 if rng.random() < 0.5 else -1.0
    return s0 * (1.0 - 2.0 * prefix)                        # sigma_i in {+1,-1}


# ---------------------------------------------------------------------------
# Spin-flip BKL: advance sigma in place to micro time u_target
# ---------------------------------------------------------------------------

def evolve_to(sigma, u_target, la, lc, lm, Lidx, Ridx, rng):
    """Continuous-time spin-flip dynamics at the glauber rates (code units):
       flip rate of site i = 1 - (gamma/2) sigma_i(sigma_{i-1}+sigma_{i+1})
         = 1+gamma (lone dissenter), 1 (wall), 1-gamma (aligned bulk)."""
    l = sigma.size
    u = 0.0
    while u < u_target:
        s = sigma * (sigma[Lidx] + sigma[Ridx])            # in {-2, 0, +2}
        rates = np.where(s < 0, la, np.where(s > 0, lc, lm))
        cum = np.cumsum(rates)
        R = cum[-1]
        dt = rng.exponential(1.0 / R)
        if u + dt >= u_target:
            break
        u += dt
        j = int(np.searchsorted(cum, rng.random() * R))
        if j >= l:
            j = l - 1
        sigma[j] = -sigma[j]


# ---------------------------------------------------------------------------
# Independent sampling: one fresh ensemble per (x_j, T_k) grid point
# ---------------------------------------------------------------------------

def _grid_point_worker(pk):
    """One independent ensemble for a single (k, j) grid point.  Top-level and
    self-contained (no module globals) so it pickles for multiprocessing."""
    k, j, l, p, gamma, u_target, ref_box, box, nsims, seedseq = pk
    la, lc, lm = 1.0 + gamma, 1.0 - gamma, 1.0
    Lidx = (np.arange(l) - 1) % l
    Ridx = (np.arange(l) + 1) % l
    sim_seeds = seedseq.spawn(nsims)
    acc = acc2 = 0.0
    for s in range(nsims):
        rng = np.random.default_rng(sim_seeds[s])
        sigma = equilibrium_ic(l, p, rng)
        m0 = sigma[ref_box].mean()                         # fresh reference
        evolve_to(sigma, u_target, la, lc, lm, Lidx, Ridx, rng)
        pr = m0 * sigma[box].mean()
        acc += pr
        acc2 += pr * pr
    mean = acc / nsims
    return k, j, mean, np.sqrt(max(acc2 / nsims - mean * mean, 0.0) / nsims)


def run_independent(l, p, gamma, ref_box, boxes, t_targets, nsims, seed, jobs=1):
    K, M = len(t_targets), boxes.shape[0]
    Cmc = np.zeros((K, M))
    Cerr = np.zeros((K, M))
    # Independent RNG sub-stream for every (k, j) grid point.
    point_seeds = np.random.SeedSequence(seed).spawn(K * M)
    tasks = [(k, j, l, p, gamma, t_targets[k], ref_box, boxes[j], nsims,
              point_seeds[k * M + j])
             for k in range(K) for j in range(M)]
    # Longest-processing-time-first: dispatch the largest-T (slowest) ensembles
    # first so the parallel makespan isn't dominated by a late long task.
    tasks.sort(key=lambda tk: -tk[5])

    done = 0
    if jobs > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            for k, j, mean, se in pool.imap_unordered(_grid_point_worker, tasks):
                Cmc[k, j] = mean
                Cerr[k, j] = se
                done += 1
                if done % M == 0:
                    print(f"  {done}/{K * M} grid points done", flush=True)
    else:
        for t in tasks:
            k, j, mean, se = _grid_point_worker(t)
            Cmc[k, j] = mean
            Cerr[k, j] = se
            done += 1
            if done % M == 0:
                print(f"  {done}/{K * M} grid points done", flush=True)
    return Cmc, Cerr


# ---------------------------------------------------------------------------
# Joint sampling: one ensemble measured at all (x_j, T_k) — cheaper, correlated
# ---------------------------------------------------------------------------

def run_joint(l, p, gamma, ref_box, boxes, t_targets, nsims, seed):
    K, M = len(t_targets), boxes.shape[0]
    la, lc, lm = 1.0 + gamma, 1.0 - gamma, 1.0
    Lidx = (np.arange(l) - 1) % l
    Ridx = (np.arange(l) + 1) % l

    Csum = np.zeros((K, M))
    Csq = np.zeros((K, M))
    sim_seeds = np.random.SeedSequence(seed).spawn(nsims)
    for s in range(nsims):
        rng = np.random.default_rng(sim_seeds[s])
        sigma = equilibrium_ic(l, p, rng)
        m0 = sigma[ref_box].mean()
        u = 0.0
        for k in range(K):
            # advance the SAME trajectory to each successive target
            evolve_to_from = t_targets[k] - u
            evolve_to(sigma, evolve_to_from, la, lc, lm, Lidx, Ridx, rng)
            u = t_targets[k]
            pr = m0 * sigma[boxes].mean(axis=1)
            Csum[k] += pr
            Csq[k] += pr * pr
    Cmc = Csum / nsims
    Cerr = np.sqrt(np.maximum(Csq / nsims - Cmc * Cmc, 0.0) / nsims)
    return Cmc, Cerr


# ---------------------------------------------------------------------------
# Exact theory: box-averaged two-time spin correlator on the ring
# ---------------------------------------------------------------------------

def theory_curve(centers, c0, w, l, eta, gamma, T_values, N):
    offs = np.arange(-w, w + 1)
    sep = offs[:, None] - offs[None, :]
    C = np.empty((len(T_values), len(centers)))
    for k, T in enumerate(T_values):
        t = 2.0 * N * N * T
        gfull = th.G_two_time_ring(np.arange(l), t, l, eta, gamma)
        for j, cj in enumerate(centers):
            idx = (cj - c0 + sep) % l
            C[k, j] = gfull[idx].mean()
    return C


# ---------------------------------------------------------------------------
# Plot-only mode: overlay the exact theory on an existing results directory
# (e.g. produced by the C program ./BoxCorr, whose txt has no theory columns)
# ---------------------------------------------------------------------------

def plot_from_dir(outdir):
    with open(os.path.join(outdir, "configuration.json")) as f:
        cfg = json.load(f)
    N, L, nu = cfg["N"], cfg["L"], cfg["nu"]
    l = cfg["l"]
    eta, gamma = cfg["eta"], cfg["gamma"]
    w = (cfg["box_sites"] - 1) // 2
    c0 = cfg["reference_site"]
    T_values = cfg["T_values"]
    K = len(T_values)

    raw = np.loadtxt(os.path.join(outdir, "box_corr.txt"))
    raw = np.atleast_2d(raw)
    xpos = raw[:, 0]
    per = (raw.shape[1] - 1) // K        # 2 cols/T (C output) or 3 (Python)
    Cmc = np.stack([raw[:, 1 + per * k] for k in range(K)])
    Cerr = np.stack([raw[:, 2 + per * k] for k in range(K)])
    centers = np.round(xpos * N).astype(int) % l

    Cth = theory_curve(centers, c0, w, l, eta, gamma, T_values, N)
    centers_fine = np.arange(l)
    Cth_fine = theory_curve(centers_fine, c0, w, l, eta, gamma, T_values, N)

    print(f"plot-only: {outdir}  N={N} l={l} nu={nu} |B|={cfg['box_sites']} "
          f"nsims={cfg.get('nsims', '?')} engine={cfg.get('engine', 'python')}")
    print("chi^2/dof (MC vs exact):")
    for k, T in enumerate(T_values):
        mask = Cerr[k] > 0
        chi2 = np.sum(((Cmc[k][mask] - Cth[k][mask]) / Cerr[k][mask]) ** 2) / mask.sum()
        print(f"  t={T:<6g}  chi2/dof={chi2:7.3f}")

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, K))
    for k, T in enumerate(T_values):
        c = colors[k]
        ax.plot(centers_fine / N, Cth_fine[k], "-", color=c, lw=1.7,
                label=f"exact  t={T:g}")
        ax.errorbar(xpos, Cmc[k], yerr=Cerr[k], fmt="o", ms=3.4, color=c,
                    alpha=0.85, lw=0.9, capsize=1.5, label=f"MC  t={T:g}")
    ax.axvline(L / 2.0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$C(x,t)=\langle m_{\rm box}(L/2,0)\,m_{\rm box}(x,t)\rangle$")
    ax.set_title(rf"Box-magnetization correlation (— exact, $\circ$ MC)"
                 "\n"
                 rf"$N={N},\ \nu={nu:g},\ \gamma={gamma:.4f},\ \eta={eta:.4f},"
                 rf"\ |B|={cfg['box_sites']},\ {cfg.get('nsims', '?')}$ sims/pt")
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    fig.tight_layout()
    png = os.path.join(outdir, "box_corr.png")
    fig.savefig(png, dpi=160)
    print(f"Saved {png}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=12)
    ap.add_argument("--L", type=int, default=4)
    ap.add_argument("--nu", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.5, help="macroscopic box width")
    ap.add_argument("--M", type=int, default=13,
                    help="number of macroscopic space points (0 => every site)")
    ap.add_argument("--T", type=str, default="0.05,0.1,0.2,0.4")
    ap.add_argument("--nsims", type=int, default=800,
                    help="realizations PER grid point (independent) or total (joint)")
    ap.add_argument("--mode", choices=["independent", "joint"], default="independent")
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallel worker processes (independent mode)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="box_corr")
    ap.add_argument("--plot-only", type=str, default=None, metavar="DIR",
                    help="no MC: overlay exact theory on an existing results "
                         "directory (e.g. from the C program ./BoxCorr)")
    args = ap.parse_args()

    if args.plot_only:
        plot_from_dir(args.plot_only)
        return

    N, L, nu = args.N, args.L, args.nu
    l = N * L
    eta, gamma, p = th.eta_of(N, nu), th.gamma_of(N, nu), th.p_of(N, nu)
    T_values = [float(x) for x in args.T.split(",")]
    K = len(T_values)

    w = max(int(round(args.h * N / 2.0)), 0)
    boxw = 2 * w + 1
    c0 = l // 2
    M = args.M if args.M > 0 else l
    centers = np.unique(np.round(np.linspace(0, l, M, endpoint=False)).astype(int) % l)
    M = centers.size

    offs = np.arange(-w, w + 1)
    boxes = (centers[:, None] + offs[None, :]) % l
    ref_box = (c0 + offs) % l
    t_targets = [T * N * N for T in T_values]

    print(f"N={N} L={L} l={l} nu={nu}  gamma={gamma:.6f} eta={eta:.6f} p={p:.6f}")
    print(f"box h={args.h} -> w={w} ({boxw} sites)  reference site {c0} (x=L/2)")
    print(f"mode={args.mode}  M={M} space pts  K={K} times {T_values}  nsims={args.nsims}")
    if args.mode == "independent":
        print(f"  => {M*K} independent ensembles of {args.nsims} "
              f"({M*K*args.nsims} trajectories total)")

    t0 = time.time()
    if args.mode == "independent":
        Cmc, Cerr = run_independent(l, p, gamma, ref_box, boxes, t_targets,
                                    args.nsims, args.seed, jobs=args.jobs)
    else:
        Cmc, Cerr = run_joint(l, p, gamma, ref_box, boxes, t_targets,
                              args.nsims, args.seed)
    print(f"MC done in {time.time() - t0:.1f}s")

    Cth = theory_curve(centers, c0, w, l, eta, gamma, T_values, N)
    # Exact solution on a FINE grid (every integer box centre) for a smooth
    # curve, evaluated independently of the coarse MC sampling grid.
    centers_fine = np.arange(l)
    xfine = centers_fine / N
    Cth_fine = theory_curve(centers_fine, c0, w, l, eta, gamma, T_values, N)

    print("\nchi^2/dof (MC vs exact)" +
          ("  [independent: should concentrate near 1]" if args.mode == "independent"
           else "  [joint: correlated points, not meaningful]") + ":")
    for k, T in enumerate(T_values):
        mask = Cerr[k] > 0
        chi2 = np.sum(((Cmc[k][mask] - Cth[k][mask]) / Cerr[k][mask]) ** 2) / mask.sum()
        print(f"  T={T:<6g}  chi2/dof={chi2:7.3f}")

    # ---- outputs ----
    xpos = centers / N
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    cfg = dict(N=N, L=L, nu=nu, l=l, gamma=gamma, eta=eta, p_equilibrium=p,
               box_width_h=args.h, box_sites=boxw, reference_site=int(c0),
               reference_x=L / 2.0, M=int(M), T_values=T_values,
               nsims=args.nsims, mode=args.mode, time_scale_formula=2.0 * N * N,
               seed=args.seed)
    with open(os.path.join(outdir, "configuration.json"), "w") as f:
        json.dump(cfg, f, indent=2)

    with open(os.path.join(outdir, "box_corr.txt"), "w") as f:
        f.write("# x")
        for T in T_values:
            f.write(f"\tMC_T{T:g}\tstderr_T{T:g}\tTH_T{T:g}")
        f.write("\n")
        for j in range(M):
            f.write(f"{xpos[j]:.6f}")
            for k in range(K):
                f.write(f"\t{Cmc[k, j]:.8f}\t{Cerr[k, j]:.8f}\t{Cth[k, j]:.8f}")
            f.write("\n")
    print(f"Wrote {os.path.join(outdir, 'box_corr.txt')}")

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, K))
    for k, T in enumerate(T_values):
        c = colors[k]
        ax.plot(xfine, Cth_fine[k], "-", color=c, lw=1.7,
                label=f"exact  T={T:g}")
        ax.errorbar(xpos, Cmc[k], yerr=Cerr[k], fmt="o", ms=3.4, color=c,
                    alpha=0.85, lw=0.9, capsize=1.5, label=f"MC  T={T:g}")
    ax.axvline(L / 2.0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlabel(r"macroscopic position $x$")
    ax.set_ylabel(r"$C(x,T)=\langle m_{\rm box}(L/2,0)\,m_{\rm box}(x,T)\rangle$")
    ax.set_title(rf"Box-magnetization correlation (— exact, $\circ$ MC, {args.mode})"
                 "\n"
                 rf"$N={N},\ \nu={nu:g},\ \gamma={gamma:.4f},\ \eta={eta:.4f},"
                 rf"\ h={args.h:g}\,(|B|={boxw}),\ {args.nsims}$ sims/pt")
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    png = os.path.join(outdir, "box_corr.png")
    fig.savefig(png, dpi=160)
    print(f"Saved {png}")


if __name__ == "__main__":
    main()
