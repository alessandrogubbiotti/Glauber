"""
finiteN_validation.py — finite-N validation suite for the Glauber dynamics
on the discrete torus T_N (paper draft, Sections 1-4), at N = 20.

Model and units are those of the paper:
    calibration (1.7):  e^{2 beta_N} = N,  zeta = tanh(beta_N) = (N-1)/(N+1)
    Glauber rates (2.5)-(2.7):  gamma = tanh(2 beta_N) = (N^2-1)/(N^2+1),
        flip rates  w in {1-gamma, 1, 1+gamma}  (hop rate 1 per direction)
    micro (lattice) time t throughout; macroscopic time tau = t / N^2.

Every observable is computed by THREE independent routes:
  (M) matrix diagonalisation — numerical eigh of the closed-hierarchy
      N x N matrix  A = -2 I + gamma (S + S^T)  (eq. (4.2)); statics by
      numerical diagonalisation of the 2 x 2 transfer matrix (eq. (1.2)).
  (B) exact formulas — wrapped Bessel kernel (4.1)/(4.3) via scipy ive,
      static correlator (4.5), boxed two-time formula (4.9).
  (MC) continuous-time Monte Carlo of the single-spin dynamics on the full
      torus (no box), exact uniformised jump chain.  Per the validation
      protocol, EVERY data point (every separation d, site x, or time t at
      which a quantity is plotted) uses its OWN independent ensemble of
      simulations, with initial conditions freshly drawn from the same
      distribution (exact conditioned-Bernoulli sampling of mu_{N,beta},
      Lemma 1.1 + Lemma 2.1, respecting the even-parity conditioning of
      Remark 4.2).

Observables:
  1. stationary space-time trajectory (single Gillespie run, figure only);
  2. equal-time correlator <sigma_0 sigma_d> at stationarity, plus the same
     correlator measured after a long run from the all-plus state
     (convergence to the invariant measure);
  3. magnetisation: profile m_x(t) from a two-domain initial condition and
     uniform-mode decay m(t) = e^{-2ct} from the all-plus state;
  4. stationary two-time correlator <sigma_y(0) sigma_x(t)>.

Outputs (in <repo>/validation_finiteN/):
  fig1_trajectory.{pdf,png}     fig2_static.{pdf,png}
  fig3_magnetization.{pdf,png}  fig4_twotime.{pdf,png}
  appendix_info.txt   appendix_validation.tex   data_finiteN.npz

Run:   python python_scripts/theory/finiteN_validation.py [--quick]
       (--quick: small ensembles, pipeline smoke test only)
"""

import os
import sys
import time

import numpy as np
from scipy.special import ive
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

N = 20
ZETA = (N - 1.0) / (N + 1.0)                # tanh(beta_N), calibration (1.7)
BETA = np.arctanh(ZETA)                     # = (1/2) log N
GAMMA = 2.0 * ZETA / (1.0 + ZETA ** 2)      # tanh(2 beta_N) = (N^2-1)/(N^2+1)
A_RATE = 1.0 + GAMMA                        # annihilation
C_RATE = 1.0 - GAMMA                        # creation
P_WALL = (1.0 - ZETA) / 2.0                 # = 1/(N+1)

QUICK = "--quick" in sys.argv

# ensembles (one INDEPENDENT ensemble of this size per plotted point)
M_STATIC = 2000 if QUICK else 20000
M_PROFILE = 2000 if QUICK else 20000
M_DECAY = 5000 if QUICK else 50000
M_TWOTIME = 2000 if QUICK else 20000

# evaluation grids (micro time)
D_STATIC = np.arange(N // 2 + 1)            # separations 0..10
T_EV_STATIC = 25.0                          # evolution of stationary start
T_EQ_FROM_PLUS = 100.0 if QUICK else 600.0  # burn-in from the all-plus state
T_PROFILE = [2.0, 10.0, 50.0]
T_DECAY = [1, 2, 5, 10, 20, 50, 100, 200, 400]
T_TWOTIME = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]
D_TWOTIME = [0, 1, 3, 6]
T_TRAJ = 2 * N * N                          # trajectory window: tau in [0,2]

SEED = 20260612
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "..", "..", "validation_finiteN")


# ─────────────────────────────────────────────────────────────────────────────
# (M) matrix-diagonalisation route
# ─────────────────────────────────────────────────────────────────────────────

def transfer_static_corr(d_arr):
    """<sigma_0 sigma_d> by numerical diagonalisation of the 2x2 transfer
    matrix (1.2): C_d = Tr(S T^d S T^{N-d}) / Tr(T^N), S = diag(1,-1)."""
    T = np.array([[np.exp(BETA), np.exp(-BETA)],
                  [np.exp(-BETA), np.exp(BETA)]])
    lam, V = np.linalg.eigh(T)
    S = V.T @ np.diag([1.0, -1.0]) @ V       # spin operator in the eigenbasis
    Z = np.sum(lam ** N)
    out = np.empty(len(d_arr))
    for i, d in enumerate(d_arr):
        out[i] = np.einsum("i,ij,j,ji->", lam ** float(d), S,
                           lam ** float(N - d), S) / Z
    return out


def hierarchy_propagator():
    """Eigendecomposition of A = -2 I + gamma (S + S^T) of eq. (4.2)."""
    A = -2.0 * np.eye(N)
    idx = np.arange(N)
    A[idx, (idx + 1) % N] += GAMMA
    A[idx, (idx - 1) % N] += GAMMA
    evals, evecs = np.linalg.eigh(A)
    return evals, evecs


_EVALS, _EVECS = hierarchy_propagator()


def propagate_matrix(v, t):
    """e^{tA} v  through the numerical eigendecomposition."""
    return _EVECS @ (np.exp(_EVALS * t) * (_EVECS.T @ v))


# ─────────────────────────────────────────────────────────────────────────────
# (B) exact-formula route (Bessel / closed sums)
# ─────────────────────────────────────────────────────────────────────────────

def static_corr_formula(d_arr):
    """Eq. (4.5):  C_d = (zeta^d + zeta^{N-d}) / (1 + zeta^N)."""
    d = np.asarray(d_arr, dtype=float)
    return (ZETA ** d + ZETA ** (N - d)) / (1.0 + ZETA ** N)


def kernel_bessel(t, m_images=60):
    """Wrapped heat kernel (4.1)/(4.3):
    K_z(t) = e^{-2t} sum_m I_{z+mN}(2 gamma t),  z = 0..N-1,
    computed with the scaled Bessel ive to avoid overflow:
    e^{-2t} I_n(2 gamma t) = e^{-2(1-gamma) t} ive(n, 2 gamma t)."""
    z = np.arange(N)
    m = np.arange(-m_images, m_images + 1)
    orders = z[:, None] + N * m[None, :]
    return np.exp(-2.0 * C_RATE * t) * ive(orders, 2.0 * GAMMA * t).sum(axis=1)


def magnetization_bessel(m0, t):
    """m_x(t) = sum_y K_{x-y}(t) m_y(0), direct circular sum (eq. (4.3))."""
    K = kernel_bessel(t)
    x = np.arange(N)
    return np.array([np.sum(K[(xi - x) % N] * m0) for xi in range(N)])


def two_time_formula(t, d):
    """Boxed eq. (4.9), stationary <sigma_y(0) sigma_x(t)> at distance d."""
    k = 2.0 * np.pi * np.arange(N) / N
    pref = (np.sqrt(1.0 - GAMMA ** 2) * (1.0 - ZETA ** N)
            / (N * (1.0 + ZETA ** N)))
    rates = 1.0 - GAMMA * np.cos(k)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = pref * np.sum(np.cos(k * d)[None, :] / rates[None, :]
                        * np.exp(-2.0 * t[:, None] * rates[None, :]), axis=1)
    return out if out.size > 1 else out[0]


def two_time_bessel(t, d):
    """Cross-check: convolution (4.8) of the wrapped kernel with (4.5)."""
    K = kernel_bessel(t)
    C = static_corr_formula(np.arange(N))
    z = np.arange(N)
    return np.sum(K[(d - z) % N] * C)


# ─────────────────────────────────────────────────────────────────────────────
# Monte Carlo: exact sampling of mu_{N,beta} and uniformised dynamics
# ─────────────────────────────────────────────────────────────────────────────

def sample_equilibrium(M, rng):
    """M independent exact samples of mu_{N,beta} (Lemma 1.1 + Lemma 2.1):
    i.i.d. Bernoulli(p) walls conditioned on even parity (rejection), then a
    uniform global colour.  Returns (spins (M,N) int8, acceptance rate)."""
    rows = []
    n_drawn = n_kept = 0
    while n_kept < M:
        batch = max(M - n_kept, 1000)
        walls = (rng.random((batch, N)) < P_WALL)
        even = walls.sum(axis=1) % 2 == 0
        n_drawn += batch
        n_kept += int(even.sum())
        rows.append(walls[even])
    walls = np.concatenate(rows, axis=0)[:M]
    flips = np.zeros((M, N), dtype=np.int64)
    flips[:, 1:] = np.cumsum(walls[:, :-1], axis=1)
    sign0 = (1 - 2 * rng.integers(0, 2, size=M)).astype(np.int8)
    spins = sign0[:, None] * np.where(flips % 2 == 0, 1, -1).astype(np.int8)
    return spins, n_kept / n_drawn


def mc_evolve(spins, t, rng):
    """Evolve M chains (rows of spins, in place) for micro time t with the
    Glauber rates w_x = 1 - (gamma/2) sigma_x (sigma_{x-1} + sigma_{x+1}).

    Exact uniformisation: dominating rate B = 2 >= 1+gamma per site; events
    arrive at rate B*N, hit a uniform site, flip with probability w/B.
    No time-discretisation error."""
    M = spins.shape[0]
    B = 2.0
    n_ev = rng.poisson(B * N * t, size=M)
    kmax = int(n_ev.max()) if M else 0
    rows = np.arange(M)
    for k in range(kmax):
        x = rng.integers(0, N, size=M)
        u = rng.random(M)
        s = spins[rows, x]
        nb = spins[rows, (x + 1) % N] + spins[rows, (x - 1) % N]
        w = 1.0 - 0.5 * GAMMA * s * nb
        flip = np.nonzero((n_ev > k) & (u * B < w))[0]
        spins[flip, x[flip]] = -s[flip]
    return spins


def gillespie_trajectory(t_max, rng):
    """Single-chain exact Gillespie run from an equilibrium sample.
    Returns (event times, spin configuration after each event incl. t=0)."""
    s, _ = sample_equilibrium(1, rng)
    s = s[0].astype(np.int8)
    t, times, configs = 0.0, [0.0], [s.copy()]
    while True:
        nb = np.roll(s, 1) + np.roll(s, -1)
        w = 1.0 - 0.5 * GAMMA * s * nb
        R = w.sum()
        t += rng.exponential(1.0 / R)
        if t > t_max:
            break
        x = rng.choice(N, p=w / R)
        s[x] = -s[x]
        times.append(t)
        configs.append(s.copy())
    return np.array(times), np.array(configs)


def mean_stderr(g):
    """Ensemble mean and standard error of per-chain estimates g (len M)."""
    return g.mean(), g.std(ddof=1) / np.sqrt(len(g))


# ─────────────────────────────────────────────────────────────────────────────
# Observable runners — one independent ensemble (and fresh seed) per point
# ─────────────────────────────────────────────────────────────────────────────

def run_static(seeds):
    """<sigma_0 sigma_d>: (a) stationary start evolved for T_EV_STATIC,
    (b) all-plus start evolved for T_EQ_FROM_PLUS.  One ensemble per d."""
    res_a = np.empty((len(D_STATIC), 2))
    res_b = np.empty((len(D_STATIC), 2))
    acc = []
    for i, d in enumerate(D_STATIC):
        rng = np.random.default_rng(seeds.spawn(1)[0])
        spins, a = sample_equilibrium(M_STATIC, rng)
        acc.append(a)
        mc_evolve(spins, T_EV_STATIC, rng)
        g = np.mean(spins * np.roll(spins, -d, axis=1), axis=1)
        res_a[i] = mean_stderr(g)

        rng = np.random.default_rng(seeds.spawn(1)[0])
        spins = np.ones((M_STATIC, N), dtype=np.int8)
        mc_evolve(spins, T_EQ_FROM_PLUS, rng)
        g = np.mean(spins * np.roll(spins, -d, axis=1), axis=1)
        res_b[i] = mean_stderr(g)
        print(f"  static d={d:2d}  stat={res_a[i,0]:+.4f}({res_a[i,1]:.4f})"
              f"  from-plus={res_b[i,0]:+.4f}({res_b[i,1]:.4f})", flush=True)
    return res_a, res_b, float(np.mean(acc))


def domain_ic():
    """Two-domain initial condition: +1 on [0, N/2), -1 on [N/2, N)."""
    m0 = np.ones(N)
    m0[N // 2:] = -1.0
    return m0


def run_profile(seeds):
    """m_x(t) from the two-domain IC; one ensemble per (x, t) point."""
    m0 = domain_ic().astype(np.int8)
    out = {}
    for t in T_PROFILE:
        res = np.empty((N, 2))
        for x in range(N):
            rng = np.random.default_rng(seeds.spawn(1)[0])
            spins = np.tile(m0, (M_PROFILE, 1))
            mc_evolve(spins, t, rng)
            col = spins[:, x].astype(float)
            res[x] = mean_stderr(col)
        out[t] = res
        print(f"  profile t={t:5.1f} done", flush=True)
    return out


def run_decay(seeds):
    """Uniform magnetisation from the all-plus state; one ensemble per t."""
    res = np.empty((len(T_DECAY), 2))
    for i, t in enumerate(T_DECAY):
        rng = np.random.default_rng(seeds.spawn(1)[0])
        spins = np.ones((M_DECAY, N), dtype=np.int8)
        mc_evolve(spins, float(t), rng)
        g = spins.mean(axis=1)
        res[i] = mean_stderr(g)
        print(f"  decay t={t:4d}  m={res[i,0]:+.5f}({res[i,1]:.5f})",
              flush=True)
    return res


def run_twotime(seeds):
    """Stationary <sigma_y(0) sigma_x(t)> at distance d; one ensemble and
    fresh equilibrium initial conditions per (d, t) point."""
    out = {}
    for d in D_TWOTIME:
        res = np.empty((len(T_TWOTIME), 2))
        for i, t in enumerate(T_TWOTIME):
            rng = np.random.default_rng(seeds.spawn(1)[0])
            s0, _ = sample_equilibrium(M_TWOTIME, rng)
            s1 = s0.copy()
            mc_evolve(s1, float(t), rng)
            g = np.mean(s0 * np.roll(s1, -d, axis=1), axis=1)
            res[i] = mean_stderr(g)
        out[d] = res
        print(f"  twotime d={d} done", flush=True)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Cross-checks between the deterministic routes
# ─────────────────────────────────────────────────────────────────────────────

def crosschecks():
    checks = {}
    d_all = np.arange(N)
    # statics: transfer matrix vs formula (4.5)
    checks["static: transfer vs (4.5)"] = np.abs(
        transfer_static_corr(d_all) - static_corr_formula(d_all)).max()
    # (4.9) at t=0 returns (4.5)
    checks["(4.9) at t=0 vs (4.5)"] = max(
        abs(two_time_formula(0.0, d) - static_corr_formula([d])[0])
        for d in d_all)
    # two-time: formula (4.9) vs matrix propagation vs Bessel convolution
    C_eq = transfer_static_corr(d_all)
    e_mat, e_bes = 0.0, 0.0
    for t in T_TWOTIME:
        Gm = propagate_matrix(C_eq, float(t))
        for d in D_TWOTIME:
            f = two_time_formula(float(t), d)
            e_mat = max(e_mat, abs(f - Gm[d]))
            e_bes = max(e_bes, abs(f - two_time_bessel(float(t), d)))
    checks["two-time: (4.9) vs matrix e^{tA}"] = e_mat
    checks["two-time: (4.9) vs Bessel conv (4.8)"] = e_bes
    # magnetisation: matrix vs Bessel kernel, and uniform mode vs e^{-2ct}
    m0 = domain_ic()
    e_pm = max(np.abs(propagate_matrix(m0, t) - magnetization_bessel(m0, t)).max()
               for t in T_PROFILE)
    checks["profile: matrix vs Bessel kernel"] = e_pm
    ones = np.ones(N)
    e_dm = max(abs(propagate_matrix(ones, float(t))[0]
                   - np.exp(-2.0 * C_RATE * t)) for t in T_DECAY)
    checks["decay: matrix vs e^{-2ct}"] = e_dm
    return checks


def chi2(mc, exact):
    """chi^2/dof of MC points (mean, stderr) against exact values;
    deterministic points (zero stderr, e.g. d=0) are skipped."""
    ok = mc[:, 1] > 0
    dev = (mc[ok, 0] - np.asarray(exact)[ok]) / mc[ok, 1]
    return float(np.mean(dev ** 2)), float(np.abs(dev).max())


# ─────────────────────────────────────────────────────────────────────────────
# Figures (article titles; run details go to the appendix files)
# ─────────────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "legend.fontsize": 8.5, "figure.dpi": 110, "savefig.dpi": 220,
    "axes.grid": True, "grid.alpha": 0.3,
})

CLR_MC = "#1f77b4"


def savefig(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUTDIR, f"{name}.{ext}"),
                    bbox_inches="tight")
    plt.close(fig)


def fig_trajectory(times, configs):
    fig, ax = plt.subplots(figsize=(7.0, 2.9))
    edges_t = np.concatenate([times, [T_TRAJ]])
    edges_x = np.arange(N + 1)
    cmap = matplotlib.colors.ListedColormap(["#2166ac", "#b2182b"])
    norm = matplotlib.colors.BoundaryNorm([-2, 0, 2], 2)
    ax.pcolormesh(edges_t, edges_x, configs.T, cmap=cmap, norm=norm,
                  shading="flat", rasterized=True)
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"site $x$")
    ax.set_title(r"A stationary space–time trajectory of the Glauber dynamics"
                 f" on $\\mathbb{{T}}_{{{N}}}$")
    ax.grid(False)
    savefig(fig, "fig1_trajectory")


def fig_static(exact, transfer, mc_a, mc_b):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    d = D_STATIC
    ax.plot(d, exact, "k-", lw=1.2, label=r"exact, eq. (4.5)")
    ax.plot(d, transfer, "s", ms=6, mfc="none", mec="crimson",
            label="transfer matrix (numerical)")
    ax.errorbar(d, mc_a[:, 0], yerr=mc_a[:, 1], fmt="o", ms=3.5,
                color=CLR_MC, capsize=2, lw=1,
                label="MC, stationary start")
    ax.errorbar(d, mc_b[:, 0], yerr=mc_b[:, 1], fmt="^", ms=3.5,
                color="darkorange", capsize=2, lw=1,
                label=f"MC, all-plus start, $t={T_EQ_FROM_PLUS:g}$")
    ax.set_xlabel(r"separation $d$")
    ax.set_ylabel(r"$\langle \sigma_0\, \sigma_d \rangle$")
    ax.set_title(f"Equal-time correlator at stationarity, $N={N}$")
    ax.legend()
    savefig(fig, "fig2_static")


def fig_magnetization(prof_mc, decay_mc):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.4))
    m0 = domain_ic()
    x = np.arange(N)
    colors = plt.cm.viridis(np.linspace(0.0, 0.75, len(T_PROFILE)))
    ax1.step(np.arange(N + 1) - 0.5, np.append(m0, m0[-1]), where="post",
             color="0.6", lw=1, label=r"$t=0$")
    for t, col in zip(T_PROFILE, colors):
        ax1.plot(x, propagate_matrix(m0, t), "-", color=col, lw=1.2,
                 label=rf"$t={t:g}$")
        ax1.plot(x, magnetization_bessel(m0, t), "--", color="k", lw=0.7)
        mc = prof_mc[t]
        ax1.errorbar(x, mc[:, 0], yerr=mc[:, 1], fmt="o", ms=2.8,
                     color=col, capsize=1.5, lw=0.8)
    ax1.plot([], [], "--", color="k", lw=0.7, label="Bessel, eq. (4.3)")
    ax1.set_xlabel(r"site $x$")
    ax1.set_ylabel(r"$m_x(t)$")
    ax1.set_title("Magnetisation profile, two-domain initial state")
    ax1.legend(ncols=2)

    t_dense = np.logspace(0, np.log10(T_DECAY[-1]), 200)
    ax2.semilogy(t_dense, np.exp(-2.0 * C_RATE * t_dense), "k-", lw=1.2,
                 label=r"$e^{-2ct}$ (slowest mode)")
    t_mark = np.array(T_DECAY, dtype=float)
    ax2.semilogy(t_mark,
                 [propagate_matrix(np.ones(N), t)[0] for t in t_mark],
                 "s", ms=6, mfc="none", mec="crimson",
                 label="matrix diagonalisation")
    ax2.errorbar(t_mark, decay_mc[:, 0], yerr=decay_mc[:, 1], fmt="o",
                 ms=3.5, color=CLR_MC, capsize=2, lw=1, label="MC")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"time $t$")
    ax2.set_ylabel(r"$m(t)$")
    ax2.set_title("Uniform-mode decay from the all-plus state")
    ax2.legend()
    fig.suptitle(f"Magnetisation sector of the Glauber dynamics, $N={N}$",
                 y=1.02)
    fig.tight_layout()
    savefig(fig, "fig3_magnetization")


def fig_twotime(tt_mc):
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    t_dense = np.logspace(np.log10(0.1), np.log10(300.0), 250)
    colors = plt.cm.tab10(np.arange(len(D_TWOTIME)))
    C_eq = transfer_static_corr(np.arange(N))
    t_mark = np.array(T_TWOTIME, dtype=float)
    for d, col in zip(D_TWOTIME, colors):
        ax.plot(t_dense, two_time_formula(t_dense, d), "-", color=col,
                lw=1.2, label=rf"$d={d}$")
        ax.plot(t_mark, [propagate_matrix(C_eq, t)[d] for t in t_mark],
                "s", ms=5, mfc="none", mec="k", mew=0.7)
        mc = tt_mc[d]
        ax.errorbar(t_mark, mc[:, 0], yerr=mc[:, 1], fmt="o", ms=3.2,
                    color=col, capsize=1.8, lw=0.9)
    ax.plot([], [], "s", ms=5, mfc="none", mec="k", mew=0.7,
            label="matrix diagonalisation")
    ax.errorbar([], [], fmt="o", ms=3.2, color="0.4", label="MC")
    ax.set_xscale("log")
    ax.set_xlabel(r"time lag $t$")
    ax.set_ylabel(r"$\langle \sigma_y(0)\, \sigma_{y+d}(t) \rangle$")
    ax.set_title(f"Stationary two-time correlator, $N={N}$"
                 "  (lines: eq. (4.9))")
    ax.legend(ncols=2)
    savefig(fig, "fig4_twotime")


# ─────────────────────────────────────────────────────────────────────────────
# Appendix output
# ─────────────────────────────────────────────────────────────────────────────

def write_appendix(checks, accept, chi2_table, elapsed):
    txt = os.path.join(OUTDIR, "appendix_info.txt")
    with open(txt, "w") as f:
        f.write("Finite-N validation of the Glauber dynamics on T_N "
                "(appendix data)\n")
        f.write("=" * 70 + "\n\n")
        f.write("Model parameters (calibration (1.7), Glauber rates "
                "(2.5)-(2.7)):\n")
        f.write(f"  N        = {N}\n")
        f.write(f"  beta_N   = (1/2) log N = {BETA:.10f}\n")
        f.write(f"  zeta     = tanh beta   = (N-1)/(N+1) = {ZETA:.10f}\n")
        f.write(f"  gamma    = tanh 2beta  = (N^2-1)/(N^2+1) = {GAMMA:.10f}\n")
        f.write(f"  a = 1+gamma = {A_RATE:.10f}   (annihilation)\n")
        f.write(f"  c = 1-gamma = {C_RATE:.10f}   (creation)\n")
        f.write(f"  p = (1-zeta)/2 = 1/(N+1) = {P_WALL:.10f}\n")
        f.write(f"  slowest relaxation rate 2c = {2*C_RATE:.10f}\n\n")
        f.write("Monte Carlo protocol:\n")
        f.write("  exact uniformised jump chain (dominating rate B=2 per "
                "site);\n  every plotted point uses its own independent "
                "ensemble with fresh\n  initial conditions; equilibrium "
                "starts use exact conditioned-Bernoulli\n  sampling of "
                "mu_{N,beta} (even-parity rejection + uniform colour).\n")
        f.write(f"  ensemble sizes: static {M_STATIC}, profile {M_PROFILE},"
                f" decay {M_DECAY}, two-time {M_TWOTIME} (per point)\n")
        f.write(f"  static stationary start evolved for t = {T_EV_STATIC}\n")
        f.write(f"  static all-plus start evolved for t = {T_EQ_FROM_PLUS}"
                f"  (2c t = {2*C_RATE*T_EQ_FROM_PLUS:.2f})\n")
        f.write(f"  sampler acceptance: measured {accept:.4f},"
                f" exact (1+zeta^N)/2 = {(1+ZETA**N)/2:.4f}\n")
        f.write(f"  master seed {SEED}; per-point seeds spawned by "
                "numpy SeedSequence\n\n")
        f.write("Cross-checks between deterministic routes "
                "(max abs deviation):\n")
        for k, v in checks.items():
            f.write(f"  {k:42s} {v:.3e}\n")
        f.write("\nMC vs exact (chi^2/dof, max |dev|/stderr per curve):\n")
        for k, (c2, mx) in chi2_table.items():
            f.write(f"  {k:42s} chi2/dof = {c2:6.3f}   max = {mx:5.2f} sigma\n")
        f.write(f"\nTotal run time: {elapsed:.1f} s"
                f"{'  (QUICK mode)' if QUICK else ''}\n")

    tex = os.path.join(OUTDIR, "appendix_validation.tex")
    with open(tex, "w") as f:
        f.write("% appendix_validation.tex — generated by "
                "finiteN_validation.py\n")
        f.write("\\begin{table}[ht]\\centering\n"
                "\\begin{tabular}{ll}\n\\hline\n")
        f.write(f"$N$ & ${N}$ \\\\\n")
        f.write(f"$\\zeta=\\tanh\\beta_N$ & $(N-1)/(N+1)={ZETA:.6f}$ \\\\\n")
        f.write(f"$\\gamma=\\tanh 2\\beta_N$ & "
                f"$(N^2-1)/(N^2+1)={GAMMA:.6f}$ \\\\\n")
        f.write(f"$(a,c)$ & $({A_RATE:.6f},\\,{C_RATE:.6f})$ \\\\\n")
        f.write(f"ensembles per point & {M_STATIC} (correlators), "
                f"{M_DECAY} (magnetisation decay) \\\\\n\\hline\n")
        f.write("\\multicolumn{2}{l}{max deviations between the "
                "deterministic routes:}\\\\\n")
        for k, v in checks.items():
            f.write(f"\\multicolumn{{2}}{{l}}{{\\quad {k}: "
                    f"$\\;{v:.1e}$}}\\\\\n")
        f.write("\\hline\n\\end{tabular}\n")
        f.write("\\caption{Validation of the finite-$N$ Monte Carlo against "
                "matrix diagonalisation and the exact Bessel formulas; every "
                "data point uses an independent ensemble with freshly "
                "sampled initial conditions.}\n\\end{table}\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)
    print(f"N={N}  zeta={ZETA:.6f}  gamma={GAMMA:.6f}  "
          f"a={A_RATE:.6f}  c={C_RATE:.6f}"
          f"{'  [QUICK]' if QUICK else ''}", flush=True)

    print("cross-checks between deterministic routes:", flush=True)
    checks = crosschecks()
    for k, v in checks.items():
        print(f"  {k:42s} {v:.3e}", flush=True)

    seeds = np.random.SeedSequence(SEED)

    print("trajectory ...", flush=True)
    rng = np.random.default_rng(seeds.spawn(1)[0])
    times, configs = gillespie_trajectory(T_TRAJ, rng)
    print(f"  {len(times)-1} events in t in [0, {T_TRAJ}]", flush=True)

    print("static correlator ...", flush=True)
    mc_a, mc_b, accept = run_static(seeds)
    print("magnetisation profile ...", flush=True)
    prof_mc = run_profile(seeds)
    print("magnetisation decay ...", flush=True)
    decay_mc = run_decay(seeds)
    print("two-time correlator ...", flush=True)
    tt_mc = run_twotime(seeds)

    # chi^2 of MC against the exact route
    exact_static = static_corr_formula(D_STATIC)
    chi2_table = {
        "static, stationary start": chi2(mc_a, exact_static),
        "static, all-plus start": chi2(mc_b, exact_static),
        "magnetisation decay": chi2(
            decay_mc, np.exp(-2.0 * C_RATE * np.array(T_DECAY, float))),
    }
    m0 = domain_ic()
    for t in T_PROFILE:
        chi2_table[f"profile t={t:g}"] = chi2(prof_mc[t],
                                              propagate_matrix(m0, t))
    for d in D_TWOTIME:
        chi2_table[f"two-time d={d}"] = chi2(
            tt_mc[d], two_time_formula(np.array(T_TWOTIME, float), d))

    print("MC vs exact:", flush=True)
    for k, (c2, mx) in chi2_table.items():
        print(f"  {k:42s} chi2/dof = {c2:6.3f}   max = {mx:5.2f} sigma",
              flush=True)

    print("figures ...", flush=True)
    fig_trajectory(times, configs)
    fig_static(exact_static, transfer_static_corr(D_STATIC), mc_a, mc_b)
    fig_magnetization(prof_mc, decay_mc)
    fig_twotime(tt_mc)

    np.savez(os.path.join(OUTDIR, "data_finiteN.npz"),
             N=N, zeta=ZETA, gamma=GAMMA, seed=SEED,
             d_static=D_STATIC, static_exact=exact_static,
             static_mc_stationary=mc_a, static_mc_fromplus=mc_b,
             t_profile=np.array(T_PROFILE),
             profile_mc=np.array([prof_mc[t] for t in T_PROFILE]),
             profile_exact=np.array([propagate_matrix(m0, t)
                                     for t in T_PROFILE]),
             t_decay=np.array(T_DECAY, float), decay_mc=decay_mc,
             t_twotime=np.array(T_TWOTIME, float),
             d_twotime=np.array(D_TWOTIME),
             twotime_mc=np.array([tt_mc[d] for d in D_TWOTIME]),
             twotime_exact=np.array(
                 [two_time_formula(np.array(T_TWOTIME, float), d)
                  for d in D_TWOTIME]),
             traj_times=times, traj_configs=configs)

    elapsed = time.time() - t0
    write_appendix(checks, accept, chi2_table, elapsed)
    print(f"done in {elapsed:.1f} s — output in {os.path.abspath(OUTDIR)}",
          flush=True)


if __name__ == "__main__":
    main()
