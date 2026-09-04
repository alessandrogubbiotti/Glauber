# encoding: utf-8
"""
Interface particle survival function.

For each interface particle track birth_time and death_time (macroscopic).
Censored = still alive at T_max.
Plot: Kaplan-Meier survival function S(tau) = P(lifetime > tau).

Speed: two-tier sampling separates the ~2*n_ifaces boundary sites (rate lm=1)
from the l-2*n_ifaces bulk sites (rate lc << 1).  Most steps O(n_ifaces) ~ O(3).
"""

import math
import bisect
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Fast BKL with particle tracking
# ---------------------------------------------------------------------------

def bkl_survival(sigma0, T_max, N, a, rng):
    """
    Returns list of (birth_macro, lifetime_macro, censored:bool).
    Two-tier sampling: boundary sites (rate lm or la, ~6 total) vs bulk (rate lc).
    """
    la = float(N) ** a
    lc = float(N) ** (a - 2)
    lm = 1.0

    l = len(sigma0)
    sig = list(sigma0.astype(int))   # pure Python list

    def s_at(k):
        return sig[k] * (sig[(k - 1) % l] + sig[(k + 1) % l])

    def rate_of_s(s):
        if s < 0: return la
        if s > 0: return lc
        return lm

    def rate_of(k):
        return rate_of_s(s_at(k))

    # Initialise rate structures
    # fast = boundary (lm) or isolated (la); slow = bulk (lc)
    rates    = [rate_of(k) for k in range(l)]
    fst      = sorted(k for k in range(l) if rates[k] != lc)
    fst_set  = set(fst)
    fst_r    = [rates[k] for k in fst]
    R_fst    = sum(fst_r)
    n_bulk   = l - len(fst)
    R_bulk   = lc * n_bulk
    R        = R_fst + R_bulk

    def _update(k):
        nonlocal R_fst, R_bulk, R, n_bulk
        old_r = rates[k]
        new_r = rate_of(k)
        if old_r == new_r:
            return
        rates[k] = new_r
        was_fst = old_r != lc
        now_fst = new_r != lc
        if was_fst and now_fst:               # fast -> fast (rate changed)
            idx = fst.index(k)
            fst_r[idx] = new_r
            R_fst += new_r - old_r
        elif was_fst and not now_fst:         # fast -> bulk
            idx = fst.index(k)
            fst.pop(idx); fst_set.discard(k); fst_r.pop(idx)
            R_fst -= old_r; R_bulk += lc; n_bulk += 1
        elif not was_fst and now_fst:         # bulk -> fast
            idx = bisect.bisect_left(fst, k)
            fst.insert(idx, k); fst_set.add(k); fst_r.insert(idx, new_r)
            R_fst += new_r; R_bulk -= lc; n_bulk -= 1
        # else: bulk -> bulk with same lc (impossible since lc==lc)
        R = R_fst + R_bulk

    # Particle tracking
    next_id = [0]
    p_birth  = {}
    p_pos    = {}
    pos2pid  = {}

    def spawn(pos, t):
        pid = next_id[0]; next_id[0] += 1
        p_birth[pid] = t; p_pos[pid] = pos; pos2pid[pos] = pid

    def kill(pos, t, rec):
        pid = pos2pid.pop(pos, None)
        if pid is None: return
        birth = p_birth.pop(pid); del p_pos[pid]
        rec.append((birth, t - birth, False))

    def move(old, new):
        pid = pos2pid.pop(old, None)
        if pid is None: return
        pos2pid[new] = pid; p_pos[pid] = new

    for j in range(l):
        if sig[j] != sig[(j + 1) % l]:
            spawn(j, 0.0)

    rec      = []
    t_micro  = 0.0
    T_end_mi = T_max * N * N

    while t_micro < T_end_mi:
        if R <= 0.0: break
        dt = rng.exponential(1.0 / R)
        if t_micro + dt > T_end_mi: break
        t_micro += dt
        t_mac    = t_micro / (N * N)

        # Two-tier sampling
        u = rng.random() * R
        if u < R_fst:
            # Sample from fast list (typically ~6 entries)
            v = rng.random() * R_fst
            cum = 0.0; k = fst[0]
            for site, rate in zip(fst, fst_r):
                cum += rate
                if cum >= v: k = site; break
        else:
            # Sample uniformly from bulk (acceptance ≈ n_bulk/l ~ 97%)
            while True:
                k = int(rng.integers(l))
                if k not in fst_set: break

        s_k = s_at(k)

        if s_k < 0:                          # annihilation
            kill((k - 1) % l, t_mac, rec)
            kill(k,            t_mac, rec)
        elif s_k > 0:                        # creation
            spawn((k - 1) % l, t_mac)
            spawn(k,            t_mac)
        else:                                # movement
            if sig[(k - 1) % l] != sig[k]:  # interface on left -> slides right
                move((k - 1) % l, k)
            else:                            # interface on right -> slides left
                move(k, (k - 1) % l)

        sig[k] = -sig[k]
        for j in ((k - 1) % l, k, (k + 1) % l):
            _update(j)

    # Censored particles still alive
    T_mac_end = T_end_mi / (N * N)
    for pid in list(p_pos.keys()):
        birth = p_birth[pid]
        rec.append((birth, T_mac_end - birth, True))

    return rec


# ---------------------------------------------------------------------------
# Kaplan-Meier estimator
# ---------------------------------------------------------------------------

def kaplan_meier(records):
    """
    Proper KM estimator — handles censored observations correctly.
    Returns (tau, S) arrays.
    """
    obs = sorted((lt, cens) for _, lt, cens in records)
    n   = len(obs)
    S   = 1.0
    tau_arr = [0.0]; S_arr = [1.0]
    n_at_risk = n
    i = 0
    while i < n:
        t_i = obs[i][0]
        j   = i; d = 0
        while j < n and obs[j][0] == t_i:
            if not obs[j][1]: d += 1   # death (not censored)
            j += 1
        if d > 0:
            S *= 1.0 - d / n_at_risk
            tau_arr.append(t_i); S_arr.append(S)
        n_at_risk -= (j - i)
        i = j
    return np.array(tau_arr), np.array(S_arr)


# ---------------------------------------------------------------------------
# Initial condition: Poisson with given number of interfaces
# ---------------------------------------------------------------------------

def fixed_ic(l, n_ifaces, rng):
    """Place n_ifaces interfaces at random positions (evenly drawn without replacement)."""
    n = n_ifaces
    if n % 2 != 0: n += 1
    pos   = sorted(rng.choice(l, size=n, replace=False))
    sigma = np.ones(l, dtype=int); color = 1; prev = 0
    for p in pos:
        sigma[prev:p] = color; color = -color; prev = p
    sigma[prev:] = color
    return sigma


def poisson_ic(l, rho_0, rng):
    n = max(int(round(rng.poisson(l * rho_0))), 2)
    if n % 2 != 0: n += 1
    return fixed_ic(l, n, rng)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run():
    L      = 4
    N      = 50        # N^2 = 2500 micro steps per macro time unit
    T_max  = 6.0       # macroscopic time
    n_runs = 60        # trajectories per a-value
    a_list = [-2, -1, 0]

    rho_0 = math.tanh(1.0) / N
    l     = L * N

    print(f"N={N}  L={L}  l={l}  rho_0={rho_0:.5f}  T_max={T_max}  n_runs={n_runs}")
    for a in a_list:
        print(f"  a={a:+d}:  la={N**a:.3g}  lc={N**(a-2):.3g}  lm=1")

    all_rec = {}
    for a in a_list:
        recs = []
        for run_i in range(n_runs):
            rng   = np.random.default_rng((a + 10) * 10000 + run_i)
            sigma = poisson_ic(l, rho_0, rng)
            recs += bkl_survival(sigma, T_max, N, a, rng)
        all_rec[a] = recs

        n_tot  = len(recs)
        n_dead = sum(1 for _, _, c in recs if not c)
        lts    = [lt for _, lt, c in recs if not c]
        med    = float(np.median(lts)) if lts else float('nan')
        print(f"  a={a:+d}: {n_tot} particles | {n_dead} deaths "
              f"| median lifetime = {med:.4f}")

    # ---- Kaplan-Meier plot ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {-2: 'C0', -1: 'C1', 0: 'C2'}
    labels = {
        -2: r'$a=-2$  subcritical  [form E]',
        -1: r'$a=-1$  CRITICAL  [form $D^{c_1}$]',
         0: r'$a=0$   Glauber  [form D]',
    }

    for a in a_list:
        tau, S = kaplan_meier(all_rec[a])
        c = colors[a]; lab = labels[a]
        axes[0].step(tau, S, where='post', lw=2, color=c, label=lab)
        m = (tau > 1e-5) & (S > 1e-3)
        if m.sum() > 2:
            axes[1].loglog(tau[m], S[m], lw=2, color=c, label=lab)

    # Power-law references
    tau_r = np.geomspace(1e-3, T_max, 400)
    axes[1].loglog(tau_r, 0.3  / np.sqrt(tau_r), 'k--', lw=1.4,
                   label=r'$\propto\tau^{-1/2}$')
    axes[1].loglog(tau_r, 0.05 / tau_r,           'k:',  lw=1.4,
                   label=r'$\propto\tau^{-1}$')

    for ax in axes:
        ax.set_xlabel(r'lifetime $\tau$  (macroscopic)', fontsize=12)
        ax.set_ylabel(r'$S(\tau)=P(\mathrm{lifetime}>\tau)$', fontsize=12)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

    axes[0].set_xlim(0, T_max)
    axes[0].set_title('Kaplan-Meier survival function — linear scale', fontsize=11)
    axes[1].set_xlim(1e-3, T_max)
    axes[1].set_title('Kaplan-Meier survival function — log-log scale', fontsize=11)

    fig.suptitle(
        f'Interface particle lifetimes  (N={N}, L={L}, {n_runs} runs, '
        f'$T_{{\\rm max}}={T_max}$)\n'
        r'Birth: $T=0$ or bulk creation.   Death: annihilation with adjacent particle.',
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig('interface_survival.png', dpi=150, bbox_inches='tight')
    print('Saved  interface_survival.png')


if __name__ == '__main__':
    run()
