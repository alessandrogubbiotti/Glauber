"""
droplet_lifetime_check.py -- numerical gates for the DROPLET-LIFETIME INTENSITY
claim of note_collaborators.tex, Sec. "The limit object: the CNVM of the marked
Brownian web" (items (a)-(d) and the paragraph "Continuum: marks sit on the
backward web ..."):

    births of droplets of (web-)lifetime > a have intensity  c^/(2 sqrt(pi a))
    per unit macroscopic space-time                        [note's units]

where the note's conventions are: the NOISY VOTER model of FINR06 (rate 1 per
site; with prob 1-p adopt a random neighbour, with prob p resample uniformly;
flip rates 1/2 next to one wall, p/2 next to none, 1-p/2 between two; gamma =
1-p; p = delta^2 lambda, lambda = c^), backward web of DIFFUSIVITY 1, effective
flip marks at rate c^/2 per unit dual-path time, c^ = lim N^2 (1/2)(1-gamma_N)
* 2 = 2.

The GLAUBER PAPER (../Interfacce/glauber_paper.tex, d:rates and e:calib) is the
SAME process at DOUBLE speed: a site flips at rate 1 next to one wall (hop rate
1 per direction per wall, diffusivity 2 in t = N^2 tau), at rate a_N =
2N^2/(N^2+1) between two walls (pair annihilation), at rate c_N = 2/(N^2+1)
next to none (pair creation); e^{2 beta_N} = N, p_N = 1/(N+1).  Hence

    tau_voter = 2 tau_glauber,   a_voter = 2 a_glauber,
    intensity per unit voter macro time = 1/2 x intensity per Glauber macro time.

What the claim reduces to (pair-coalescence survival).  A nucleated pair is two
adjacent walls; their relative coordinate d is the birth-death chain with rate
2r per direction for d >= 2 (r = hop rate per direction = 1 in Glauber units)
and, from d = 1, rate 2r to d = 2 and absorption (annihilation) at rate A = a_N.
Its survival is asymptotically P_1(T > t) ~ h(1) sqrt(2/(pi 4 r t)), h(1) =
2r/A (the harmonic function of the defect chain), so the finite-N intensity

    I_N(a) = N^3 c_N P_1(T > N^2 a)   [Glauber macro units, births/(dx dt)]

tends to (2 c^/a^)/sqrt(2 pi a) = sqrt(2/(pi a)) with c^ = lim N^2 c_N = 2,
a^ = lim a_N = 2; converted to voter units this is 1/sqrt(pi a_voter) =
c^/(2 sqrt(pi a_voter)), i.e. the note's formula, PROVIDED the lifetime of a
mark is the coalescence time of its two released forward-web paths ("web-age",
FINR Rem. 2.2) -- which is what the note asserts.

Gates:
  [A] EXACT pair-coalescence survival P_1(T > t): matrix exponential of the
      master equation on a truncated lattice (spectral decomposition of the
      symmetric tridiagonal generator; truncation checked by enlarging the box;
      cross-checked against scipy expm_multiply and against an independent
      Laplace-transform route, Talbot inversion of the renewal formula
      phi(s) = A/(s + 2r + A - 2r psi(s)) in mpmath), versus
        * the closed form e^{-4rt}(I_0(4rt) + I_1(4rt)) when A = 2r  (1e-10),
        * the asymptotic h(1) sqrt(2/(pi 4rt)), h(1) = 2r/A, for A != 2r.
  [B] FINITE-N INTENSITY I_N(a) for N in (25,50,100,200,400), a in (0.01, 0.05,
      0.2, 1.0); ratio to the Glauber-unit constant 2/sqrt(2 pi a) and, in the
      note's units, to c^/(2 sqrt(pi a_voter)) = 1/sqrt(pi a_voter); gate:
      within 2% at N = 400 for a >= 0.05.
  [C] MONTE CARLO of the actual lattice model (Glauber units) on a ring of
      L*N bonds with background walls in equilibrium, event-driven (exact
      Gillespie), FAITHFUL forward-web ghost paths: walls and ghosts are the
      same kind of object (coalescing rate-1-per-direction walkers); a wall-wall
      meeting annihilates the colour (both merge into one ghost path); each
      creation releases two labelled paths and its WEB-AGE is the time the two
      labels meet on one path; its VISIBLE lifetime is the time until either of
      its two walls is annihilated (by its partner or by any other wall).
      Output: empirical intensity of births with lifetime > N^2 a for both
      definitions vs the exact I_N(a) of [B] (times the exact bulk-site
      fraction (N/(N+1))^2, because creation only happens at sites with no
      adjacent wall), with Poisson error bars.

Run:  python droplet_lifetime_check.py [--fast]      (exit code 0 iff all pass)
      --fast : smaller Monte-Carlo budget (N = 50 only, 64 runs, ~10 s)
      full   : N = 50, 100, 200 with ~1e6 tracked creations each (~5 min on 8 cores)
"""
import sys, os, time
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from scipy.special import ive

_results = []
INF = float("inf")


def check(name, ok, detail):
    _results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n         {detail}")


# ============================================================================
# [A]  exact pair-coalescence survival
# ============================================================================
def survival_spectral(r, A, ts, D):
    """P_1(T > t) for the relative-coordinate chain (rate 2r per direction for
    d >= 2; from d = 1: to 2 at rate 2r, absorbed at rate A), by the exact
    exponential of the generator truncated to d in {1..D} (absorbing beyond D;
    the chain is reversible w.r.t. counting measure, so the generator is a
    symmetric tridiagonal matrix and e^{tQ} 1 is a spectral sum)."""
    diag = np.full(D, -4.0 * r)
    diag[0] = -(2.0 * r + A)
    off = np.full(D - 1, 2.0 * r)
    w, V = eigh_tridiagonal(diag, off)
    amp = V[0, :] * V.sum(axis=0)          # <e_1, v_k> <v_k, 1>
    ts = np.atleast_1d(np.asarray(ts, dtype=float))
    return np.array([float(np.dot(amp, np.exp(w * t))) for t in ts])


def box_size(r, t_max, k=8.0):
    """truncation D = k standard deviations of the relative coordinate."""
    return int(k * np.sqrt(4.0 * r * t_max)) + 100


def survival_expm(r, A, t, D):
    """same quantity via scipy expm_multiply (Al-Mohy--Higham), independent code."""
    diag = np.full(D, -4.0 * r)
    diag[0] = -(2.0 * r + A)
    off = np.full(D - 1, 2.0 * r)
    Q = diags([off, diag, off], [-1, 0, 1], format="csr")
    e1 = np.zeros(D)
    e1[0] = 1.0
    return float(expm_multiply(Q.T * t, e1).sum())       # row 1 of e^{tQ}, summed


def survival_laplace(r, A, t, dps=30):
    """Independent route, no lattice truncation: renewal formula in Laplace
    space.  From d = 1 the chain leaves at rate 2r + A; with prob A/(2r+A) it is
    absorbed, otherwise it jumps to 2 and returns to 1 after a free first-passage
    time (transform psi); hence phi(s) = E_1 e^{-sT} = A / (s + 2r + A - 2r psi(s)),
    psi(s) = [(s+4r) - sqrt((s+4r)^2 - 16 r^2)]/(4r), and P_1(T > t) is the Talbot
    inverse of (1 - phi)/s."""
    import mpmath as mp
    mp.mp.dps = dps
    r_, A_ = mp.mpf(r), mp.mpf(A)

    def F(s):
        psi = ((s + 4 * r_) - mp.sqrt((s + 4 * r_) ** 2 - 16 * r_ ** 2)) / (4 * r_)
        phi = A_ / (s + 2 * r_ + A_ - 2 * r_ * psi)
        return (1 - phi) / s
    return float(mp.invertlaplace(F, mp.mpf(t), method="talbot"))


def survival_bessel(r, t):
    """closed form for A = 2r (free walk absorbed at 0): e^{-4rt}(I_0 + I_1)(4rt)."""
    x = 4.0 * r * t
    return ive(0, x) + ive(1, x)


def asymptotic(r, A, t):
    return (2.0 * r / A) * np.sqrt(2.0 / (np.pi * 4.0 * r * t))


def gate_A():
    print("\n[A] exact pair-coalescence survival P_1(T > t)")
    r = 1.0
    ts = np.array([0.125, 1.25, 12.5, 125.0, 1250.0, 12500.0, 25000.0])   # 4rt up to 1e5
    D = box_size(r, ts.max())
    # A = 2r: Bessel closed form
    P = survival_spectral(r, 2 * r, ts, D)
    B = survival_bessel(r, ts)
    err = np.abs(P - B)
    check("A = 2r: spectral master equation == e^{-4rt}(I_0+I_1)(4rt) to 1e-10 "
          f"(7 times, 4rt in [0.5, 1e5], D = {D})", err.max() < 1e-10,
          "max |diff| = %.1e   (P_1 ranges %.3e .. %.3e)" % (err.max(), P.min(), P.max()))
    # truncation: enlarge the box
    P2 = survival_spectral(r, 2 * r, ts[-1:], int(1.5 * D))
    trunc = abs(P2[0] - P[-1])
    check("truncation: D -> 1.5 D changes P_1(T > t_max) by < 1e-12",
          trunc < 1e-12, f"|diff| = {trunc:.1e} at 4rt = {4*r*ts[-1]:.0f}")
    # expm_multiply cross-check (independent code path), moderate t
    t0 = 125.0
    D0 = box_size(r, t0)
    for A in (2.0, 1.5):
        Pe = survival_expm(r, A, t0, D0)
        Ps = survival_spectral(r, A, [t0], D0)[0]
        check(f"A = {A}: expm_multiply == spectral at 4rt = {4*r*t0:.0f} (1e-10)",
              abs(Pe - Ps) < 1e-10, f"expm {Pe:.12f}  spectral {Ps:.12f}  diff {abs(Pe-Ps):.1e}")
    # Laplace/Talbot cross-check, A != 2r, no truncation at all
    okL, rows = True, []
    for A, t in ((1.5, 125.0), (2.4, 1250.0), (2 * 2500 / 2501.0, 2500.0 * 0.2)):
        Pl = survival_laplace(r, A, t)
        Ps = survival_spectral(r, A, [t], box_size(r, t))[0]
        okL &= abs(Pl - Ps) < 1e-9 * max(1.0, Ps)
        rows.append(f"A={A:.6f} t={t:g}: Talbot {Pl:.12e} spectral {Ps:.12e} diff {abs(Pl-Ps):.1e}")
    check("A != 2r: Talbot inversion of the renewal formula == spectral (1e-9)", okL,
          "\n         ".join(rows))
    # asymptotic h(1) sqrt(2/(pi 4rt)), h(1) = 2r/A
    okas, rows = True, []
    for A in (1.0, 2.0, 3.0):
        P = survival_spectral(r, A, ts, D)
        ratio = P / asymptotic(r, A, ts)
        dev = np.abs(ratio - 1)
        okas &= dev[-1] < 1e-3 and np.all(np.diff(dev[2:]) < 0)
        rows.append(f"A={A}: P/asym at 4rt=50,500,5e3,5e4,1e5 = " +
                    " ".join(f"{x:.6f}" for x in ratio[2:]))
    check("P_1(T>t) / [h(1) sqrt(2/(pi 4rt))] -> 1 for A = 1, 2, 3 (h(1) = 2r/A): "
          "within 1e-3 at 4rt = 1e5, deviation decreasing", okas, "\n         ".join(rows))


# ============================================================================
# [B]  finite-N intensity
# ============================================================================
def rates_N(N):
    return 2.0 * N * N / (N * N + 1.0), 2.0 / (N * N + 1.0)        # a_N, c_N


def I_N_exact(N, a_list, k=8.0):
    """I_N(a) = N^3 c_N P_1(T > N^2 a), Glauber macro units, exact (spectral)."""
    aN, cN = rates_N(N)
    a_arr = np.atleast_1d(np.asarray(a_list, dtype=float))
    ts = N * N * a_arr
    D = box_size(1.0, ts.max(), k)
    return N ** 3 * cN * survival_spectral(1.0, aN, ts, D)


def gate_B():
    print("\n[B] finite-N intensity I_N(a) = N^3 c_N P_1(T > N^2 a)  (Glauber units)")
    Ns = (25, 50, 100, 200, 400)
    a_list = (0.01, 0.05, 0.2, 1.0)
    c_hat, a_hat = 2.0, 2.0
    print("  ratio_G = I_N(a) sqrt(2 pi a) / (2 c^/a^)      [-> 1; Glauber units, c^=a^=2]")
    print("  ratio_V = I_V(a_v) / [c^/(2 sqrt(pi a_v))]  [-> 1; note's units: a_v = 2a,"
          " I_V = I_N/2]")
    print("  (the two ratios coincide identically: the unit conversion is exact)")
    hdr = "   N   " + "".join(f"{'a=%g' % a:>14}" for a in a_list)
    print("  I_N(a):"); print(hdr)
    table = {}
    for N in Ns:
        I = I_N_exact(N, a_list)
        table[N] = I
        print(f"  {N:4d}  " + "".join(f"{x:14.6f}" for x in I))
    print("  ratio_G (= ratio_V):"); print(hdr)
    ratios = {}
    for N in Ns:
        I = table[N]
        rG = I * np.sqrt(2 * np.pi * np.asarray(a_list)) / (2 * c_hat / a_hat)
        a_v = 2 * np.asarray(a_list)
        rV = (I / 2.0) / (c_hat / (2 * np.sqrt(np.pi * a_v)))
        assert np.allclose(rG, rV, rtol=1e-12)
        ratios[N] = rG
        print(f"  {N:4d}  " + "".join(f"{x:14.8f}" for x in rG))
    dev = {N: np.abs(ratios[N] - 1) for N in Ns}
    print("  |ratio - 1| * N^2 a   (the finite-size correction is O(1/(N^2 a)), not O(1/N)):")
    print(hdr)
    for N in Ns:
        print(f"  {N:4d}  " + "".join(f"{dev[N][j]*N*N*a:14.5f}" for j, a in enumerate(a_list)))
    ok400 = all(dev[400][j] < 0.02 for j, a in enumerate(a_list) if a >= 0.05)
    check("N = 400: ratio within 2% for a >= 0.05", ok400,
          "max |ratio-1| = %.2e" % max(dev[400][j] for j, a in enumerate(a_list) if a >= 0.05))
    trend = all(np.all(np.diff([dev[N][j] for N in Ns]) < 0) for j in range(len(a_list)))
    check("deviation decreases monotonically in N at every a (trend ~ 1/(N^2 a))", trend,
          "limit constant: sqrt(2/(pi a)) [Glauber]  <->  1/sqrt(pi a_v) = c^/(2 sqrt(pi a_v)) [note]")
    return table


# ============================================================================
# [C]  Monte Carlo of the lattice model with faithful forward-web ghosts
# ============================================================================
def mc_run(args):
    """One ring.  Returns (birth, vis, web) arrays (micro time) and diagnostics.

    Walls live on bonds of a ring of l = L N bonds.  Sites: no adjacent wall ->
    creation at rate c_N; one adjacent wall -> flip at rate 1 (= wall hop at
    rate 1 per direction); both adjacent bonds walls -> flip at rate a_N (pair
    annihilation).  Forward web: every path (wall or ghost) hops at rate 1 per
    direction; landing on an occupied bond coalesces the two paths.  Wall+wall
    -> (accept w.p. a_N/2, so that the pair rate is a_N) both colours vanish
    and ONE ghost path continues at the target bond; wall+ghost -> wall;
    ghost+ghost -> ghost.  Every path carries a dict {mark: sides}; when both
    sides of a mark sit on one path its web-age is recorded and the label is
    dropped; label-free ghosts are deleted.  A creation on a bond carrying a
    ghost turns that ghost into the new wall (same web path)."""
    N, L, T_burn, T_end, seed, a_max = args
    rng = np.random.default_rng(seed)
    l = L * N
    aN, cN = rates_N(N)
    acc = aN / 2.0
    pN = 1.0 / (N + 1.0)
    while True:
        w0 = np.nonzero(rng.random(l) < pN)[0]
        if len(w0) % 2 == 0:
            break
    occ = [-1] * l
    pos, isw, lab, wm, slot = [], [], [], [], []
    active, free = [], []
    birth, vis, web = [], [], []

    def new_path(b, w, labels, mark):
        if free:
            i = free.pop()
            pos[i] = b; isw[i] = w; lab[i] = labels; wm[i] = mark
        else:
            i = len(pos)
            pos.append(b); isw.append(w); lab.append(labels); wm.append(mark); slot.append(-1)
        slot[i] = len(active); active.append(i); occ[b] = i
        return i

    def kill(i):
        k = slot[i]; last = active.pop()
        if last != i:
            active[k] = last; slot[last] = k
        slot[i] = -1; free.append(i); lab[i] = None

    def merge(dst, src, t):
        for m, bits in src.items():
            if m in dst:
                if dst[m] | bits == 3:
                    web[m] = t
                    del dst[m]
                else:                       # cannot happen (each label exists once)
                    dst[m] |= bits
            else:
                dst[m] = bits

    for b in w0:
        new_path(int(b), True, {}, -1)

    t = 0.0
    t_track = N * N * a_max
    sweep_dt = 0.05 * N * N
    next_sweep = sweep_dt
    cr_rate = cN * l
    BATCH = 1 << 16
    ub = rng.random(BATCH); eb = rng.exponential(size=BATCH); kb = 0
    n_ev = 0; n_cre = 0; n_cre_att = 0
    path_samples = []; wall_samples = []
    while True:
        if kb == BATCH:
            ub = rng.random(BATCH); eb = rng.exponential(size=BATCH); kb = 0
        if t >= next_sweep:                  # prune BEFORE the rate is built
            next_sweep += sweep_dt
            path_samples.append(len(active))
            wall_samples.append(sum(1 for i in active if isw[i]))
            for i in list(active):
                d = lab[i]
                if d:
                    stale = [m for m in d if t - birth[m] > t_track]
                    for m in stale:
                        del d[m]
                if not isw[i] and not d:
                    occ[pos[i]] = -1; kill(i)
        n = len(active)
        R = 2.0 * n + cr_rate
        t += eb[kb] / R
        if t >= T_end:
            break
        u0 = ub[kb] * R
        kb += 1
        n_ev += 1
        if u0 < 2.0 * n:
            k = int(u0)
            i = active[k >> 1]
            b = pos[i]
            b2 = (b + 1) if (k & 1) else (b - 1)
            if b2 == l:
                b2 = 0
            elif b2 < 0:
                b2 = l - 1
            j = occ[b2]
            if j == -1:
                occ[b] = -1; occ[b2] = i; pos[i] = b2
                continue
            wi, wj = isw[i], isw[j]
            if wi and wj:
                if rng.random() >= acc:
                    continue                               # null event: pair rate a_N
                for q in (i, j):
                    m = wm[q]
                    if m >= 0 and vis[m] == INF:
                        vis[m] = t
                occ[b] = -1
                di, dj = lab[i], lab[j]
                kill(i)
                isw[j] = False; wm[j] = -1
                merge(dj, di, t)
                if not dj:
                    occ[b2] = -1; kill(j)
            elif wi:                                       # wall onto ghost
                merge(lab[i], lab[j], t)
                kill(j)
                occ[b] = -1; occ[b2] = i; pos[i] = b2
            elif wj:                                       # ghost onto wall
                merge(lab[j], lab[i], t)
                occ[b] = -1; kill(i)
            else:                                          # ghost onto ghost
                dj = lab[j]
                merge(dj, lab[i], t)
                occ[b] = -1; kill(i)
                if not dj:
                    occ[b2] = -1; kill(j)
        else:
            x = int((u0 - 2.0 * n) / cN)
            if x >= l:
                x = l - 1
            n_cre_att += 1
            bl = x - 1 if x > 0 else l - 1
            br = x
            jl, jr = occ[bl], occ[br]
            if (jl != -1 and isw[jl]) or (jr != -1 and isw[jr]):
                continue                                   # not a bulk site
            n_cre += 1
            if t >= T_burn:
                m = len(birth)
                birth.append(t); vis.append(INF); web.append(INF)
                labs = ({m: 1}, {m: 2})
            else:
                m = -1
                labs = ({}, {})
            for bond, jj, lb in ((bl, jl, labs[0]), (br, jr, labs[1])):
                if jj == -1:
                    new_path(bond, True, lb, m)
                else:
                    isw[jj] = True; wm[jj] = m
                    merge(lab[jj], lb, t)
    diag = dict(n_ev=n_ev, n_cre=n_cre, n_cre_att=n_cre_att,
                mean_paths=float(np.mean(path_samples)) if path_samples else 0.0,
                max_paths=max(path_samples) if path_samples else 0,
                mean_walls=float(np.mean(wall_samples)) if wall_samples else 0.0)
    return (np.array(birth), np.array(vis), np.array(web), diag)


def mc_campaign(N, L, n_runs, T_burn_macro, T_obs_macro, a_list, seed0, workers):
    T_burn = T_burn_macro * N * N
    T_end = (T_burn_macro + T_obs_macro) * N * N
    a_max = max(a_list)
    jobs = [(N, L, T_burn, T_end, seed0 + k, a_max) for k in range(n_runs)]
    t0 = time.time()
    if workers > 1:
        import multiprocessing as mp
        with mp.Pool(workers) as pool:
            outs = pool.map(mc_run, jobs)
    else:
        outs = [mc_run(j) for j in jobs]
    wall = time.time() - t0
    birth = np.concatenate([o[0] for o in outs])
    vis = np.concatenate([o[1] for o in outs])
    web = np.concatenate([o[2] for o in outs])
    diag = dict(n_ev=sum(o[3]["n_ev"] for o in outs),
                n_cre=sum(o[3]["n_cre"] for o in outs),
                n_cre_att=sum(o[3]["n_cre_att"] for o in outs),
                mean_paths=np.mean([o[3]["mean_paths"] for o in outs]),
                max_paths=max(o[3]["max_paths"] for o in outs),
                mean_walls=np.mean([o[3]["mean_walls"] for o in outs]),
                wall=wall, n_marks=len(birth))
    rows = []
    for a in a_list:
        ta = N * N * a
        elig = birth <= T_end - ta
        n_el = int(elig.sum())
        n_vis = int((elig & (vis - birth > ta)).sum())
        n_web = int((elig & (web - birth > ta)).sum())
        vol = n_runs * L * (T_obs_macro - a)            # macro space-time volume
        rows.append(dict(a=a, n_el=n_el, n_vis=n_vis, n_web=n_web, vol=vol,
                         I_vis=n_vis / vol, e_vis=np.sqrt(max(n_vis, 1)) / vol,
                         I_web=n_web / vol, e_web=np.sqrt(max(n_web, 1)) / vol,
                         I_birth=n_el / vol, e_birth=np.sqrt(max(n_el, 1)) / vol))
    return rows, diag


def gate_C(N, L, n_runs, a_list, workers, seed0=12345):
    print(f"\n[C] Monte Carlo, N = {N}, ring of {L*N} bonds (L = {L}), {n_runs} runs, "
          f"burn-in 1, observation window 3 (macro units), {workers} worker(s)")
    aN, cN = rates_N(N)
    f_bulk = (N / (N + 1.0)) ** 2                 # P(both adjacent bonds empty) in equilibrium
    I_exact = I_N_exact(N, a_list)
    rows, diag = mc_campaign(N, L, n_runs, 1.0, 3.0, a_list, seed0, workers)
    print(f"  wall time {diag['wall']:.1f} s, {diag['n_ev']:.3e} events, "
          f"{diag['n_marks']} tracked creations (+{diag['n_cre']-diag['n_marks']} in burn-in), "
          f"mean paths {diag['mean_paths']:.1f} (walls {diag['mean_walls']:.1f}, "
          f"max {diag['max_paths']}); events/s = {diag['n_ev']/diag['wall']/max(1,workers):.2e} per worker")
    Ib_th = N ** 3 * cN * f_bulk
    print(f"  birth intensity: MC {rows[0]['I_birth']:.3f} +- {rows[0]['e_birth']:.3f}  vs  "
          f"N^3 c_N (N/(N+1))^2 = {Ib_th:.3f}   (N^3 c_N = {N**3*cN:.3f})")
    print("  intensity of births with lifetime > N^2 a, per unit macro space-time (Glauber units):")
    print(f"  {'a':>6} {'n_elig':>7} {'I_web (MC)':>19} {'I_vis (MC)':>19} {'I_N(a)':>9} "
          f"{'I_N fb':>9} {'web/(I_N fb)':>13} {'vis/web':>11} {'-ln(r)/sqrt(a)':>14}")
    ok_web, ok_birth = True, True
    ratios_vw, err_vw = [], []
    for k, (row, Ie) in enumerate(zip(rows, I_exact)):
        ref = Ie * f_bulk
        rw = row["I_web"] / ref
        ew = row["e_web"] / ref
        rv = row["I_vis"] / row["I_web"] if row["I_web"] > 0 else float("nan")
        ev = rv * np.sqrt(1.0 / max(row["n_vis"], 1) + 1.0 / max(row["n_web"], 1))
        ratios_vw.append(rv); err_vw.append(ev)
        ok_web &= abs(rw - 1) < 3 * ew + 0.01
        kap = -np.log(rv) / np.sqrt(row['a']) if rv > 0 else float('inf')
        print(f"  {row['a']:6.3f} {row['n_el']:7d} {row['I_web']:9.4f} +- {row['e_web']:7.4f} "
              f"{row['I_vis']:9.4f} +- {row['e_vis']:7.4f} {Ie:9.4f} {ref:9.4f} "
              f"{rw:7.4f}+-{ew:5.4f} {rv:5.3f}+-{ev:5.3f} {kap:14.2f}")
    zb = abs(rows[0]["I_birth"] - Ib_th) / rows[0]["e_birth"]
    check(f"N={N}: birth intensity == N^3 c_N (N/(N+1))^2 within 3 sigma", zb < 3,
          f"{zb:.2f} sigma")
    check(f"N={N}: WEB-AGE intensity == exact I_N(a) x bulk fraction at every a "
          "(|ratio-1| < 3 sigma + 1%)", ok_web,
          "the web-age of a mark is the pair-coalescence time of its two released "
          "paths, independent of the background: the note's identification")
    k0 = list(a_list).index(0.005) if 0.005 in a_list else 0
    r0, rN = ratios_vw[k0], ratios_vw[-1]
    e0, eN = err_vw[k0], err_vw[-1]
    mono = all(ratios_vw[k + 1] <= ratios_vw[k] + 2 * np.hypot(err_vw[k], err_vw[k + 1])
               for k in range(len(ratios_vw) - 1))
    check(f"N={N}: VISIBLE/WEB ratio tends to 1 as a -> 0 (> 0.6 at a = {a_list[k0]}, "
          f"increasing towards 1 below), falls below 0.1 at a = {a_list[-1]}, "
          "non-increasing in a (within 2 sigma)",
          r0 > 0.6 and rN < 0.1 and (r0 - rN) > 5 * np.hypot(e0, eN) and mono,
          f"ratio at a={a_list[0]}: {ratios_vw[0]:.3f}+-{err_vw[0]:.3f};  a={a_list[k0]}: "
          f"{r0:.3f}+-{e0:.3f};  a={a_list[-1]}: {rN:.3f}+-{eN:.3f}   "
          "(visible lifetime < web-age whenever a droplet wall is annihilated by a "
          "background wall, unit density: O(sqrt a) loss at small a, total at a = O(1))")
    return rows


def main(argv):
    fast = "--fast" in argv
    print("=" * 78)
    print(" Droplet-lifetime intensity  c^/(2 sqrt(pi a))  (note_collaborators, CNVM section)")
    print("=" * 78)
    t0 = time.time()
    gate_A()
    gate_B()
    a_list = (0.001, 0.002, 0.005, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
    L = 8
    try:
        workers = max(1, min(8, os.cpu_count() or 1))
    except Exception:
        workers = 1
    # budgets: ~1e6 tracked creations per N.  Measured (8 logical cores, 2026-08):
    # N=50 ~40 s, N=100 ~80 s, N=200 ~165 s (a single core does ~1.4e6 events/s,
    # the pool only ~1.5x that in total); a 2x budget gave identical conclusions.
    gate_C(50, L, 64 if fast else 400, a_list, workers)
    if not fast:
        gate_C(100, L, 200, a_list, workers, seed0=777)
        gate_C(200, L, 100, a_list, workers, seed0=31415)
    print("\n" + "=" * 78)
    n_pass, n = sum(_results), len(_results)
    print(f" {n_pass}/{n} gates passed.   total wall time {time.time()-t0:.0f} s")
    print("=" * 78)
    return 0 if n_pass == n else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
