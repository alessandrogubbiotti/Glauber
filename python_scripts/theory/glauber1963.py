"""
Numerical implementation of Glauber (1963) stochastic Ising dynamics.
All correlators for a periodic ring of N spins.

Key parameter relations:
  gamma = tanh(2*beta*J)   (rate parameter)
  eta   = tanh(beta*J)     (correlation length parameter)
  gamma = 2*eta / (1 + eta^2)

Methods implemented
-------------------
1. DST  -- exact spectral solver for equal-time r_n(t)
2. FFT  -- exact spectral solver for two-time G_n(tau)
3. Bessel -- direct evaluation via exponentially-scaled ive()
4. KMC  -- kinetic Monte Carlo validation
"""

import numpy as np
from numpy.fft import fft, ifft
from scipy.fft import dst, idst
from scipy.special import ive          # ive(n,z) = exp(-|z|)*I_n(z)
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# Physical parameter helpers
# ─────────────────────────────────────────────────────────────────────────────

def beta_to_gamma_eta(beta_J: float):
    """Return (gamma, eta) from dimensionless beta*J."""
    eta   = np.tanh(beta_J)
    gamma = 2.0 * eta / (1.0 + eta**2)   # = tanh(2*beta*J)
    return gamma, eta


# ─────────────────────────────────────────────────────────────────────────────
# Static (equilibrium) correlator on a ring of N sites
# ─────────────────────────────────────────────────────────────────────────────

def r_eq_ring(N: int, eta: float) -> np.ndarray:
    """
    r_n^eq = (eta^n + eta^{N-n}) / (1 + eta^N),  n = 0, …, N-1.
    Returns array of length N.
    """
    n     = np.arange(N, dtype=float)
    etaN  = eta**N
    # use log-space to avoid overflow for large N and eta->1
    # but eta<1 always so eta^N < 1 => no overflow
    return (eta**n + eta**(N - n)) / (1.0 + etaN)


# ─────────────────────────────────────────────────────────────────────────────
# Method 1: DST spectral solver — equal-time pair correlator r_n(t)
# ─────────────────────────────────────────────────────────────────────────────

def equal_time_dst(r0_arr: np.ndarray, t_arr: np.ndarray,
                   gamma: float, eta: float, alpha: float = 1.0) -> np.ndarray:
    """
    Evolve the equal-time correlator r_n(t) using DST-I eigenmodes.

    r0_arr : initial condition r_n(0), length N  (r_0 = 1 enforced internally)
    t_arr  : times at which to evaluate, shape (T,)
    Returns : array shape (T, N) with r_n(t) for n = 0, …, N-1.

    The equation for u_n = r_n - r_n^eq  (n=1..N-1) is
        du/dt = M u,   M tridiagonal, u_0 = u_N = 0  (Dirichlet BC)
    Eigenmodes: sin(n k pi/N), eigenvalue lambda_k = -2*alpha*(1 - gamma*cos(k*pi/N))
    """
    N      = len(r0_arr)
    r_inf  = r_eq_ring(N, eta)
    u0     = r0_arr.copy()
    u0[0]  = 0.0           # deviation at n=0 is always 0 (r_0 = 1)

    # DST-I of u0[1..N-1]  (length N-1 transform)
    u_inner = u0[1:]                                   # shape (N-1,)
    u_inf_inner = r_inf[1:]
    dev_inner   = u_inner - u_inf_inner

    # scipy dst type=1: X_k = 2 sum_{n=1}^{N-1} x_n sin(n k pi / N), k=1..N-1
    # inverse: x_n = (1/N) sum_{k=1}^{N-1} X_k sin(n k pi / N)
    C = dst(dev_inner, type=1)                         # shape (N-1,)

    k      = np.arange(1, N, dtype=float)
    lam_k  = -2.0 * alpha * (1.0 - gamma * np.cos(k * np.pi / N))  # eigenvalues

    T     = len(t_arr)
    result = np.empty((T, N))
    result[:, 0] = 1.0                                 # r_0 = 1 always

    for i, t in enumerate(t_arr):
        decay = np.exp(lam_k * t)                      # shape (N-1,)
        Ct    = C * decay
        dev_t = idst(Ct, type=1)                       # idst is already the true inverse of dst
        result[i, 1:] = r_inf[1:] + dev_t

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Method 2: FFT spectral solver — two-time correlator G_n(tau)
# ─────────────────────────────────────────────────────────────────────────────

def two_time_fft(n_val: int, tau_arr: np.ndarray,
                 N: int, gamma: float, eta: float,
                 alpha: float = 1.0) -> np.ndarray:
    """
    G_n(tau) = (1/N) sum_m  S_N(q_m) exp(-alpha*(1-gamma*cos(q_m))*tau) * exp(i q_m n)

    Returns the real part (G_n is real by time-reversal symmetry).
    """
    r_inf = r_eq_ring(N, eta)            # equilibrium correlator = G_n(0)
    S_q   = np.real(fft(r_inf))         # S_N(q_m), real since r_n = r_{N-n}

    m     = np.arange(N, dtype=float)
    q_m   = 2.0 * np.pi * m / N
    rate  = alpha * (1.0 - gamma * np.cos(q_m))   # decay rate for each mode

    result = np.empty(len(tau_arr))
    for i, tau in enumerate(tau_arr):
        decay  = np.exp(-rate * tau)
        # G_n(tau) = IFFT[ S_q * decay ]  evaluated at index n
        G_full = np.real(ifft(S_q * decay))
        result[i] = G_full[n_val % N]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Method 3: Direct Bessel series (infinite chain or ring with winding sum)
# ─────────────────────────────────────────────────────────────────────────────

def _bessel_propagator(n_vals: np.ndarray, t: float,
                       gamma: float, alpha: float,
                       tol: float = 1e-12) -> np.ndarray:
    """
    K_n(t) = exp(-alpha*t) * I_n(alpha*gamma*t)

    Uses ive(n, z) = exp(-z)*I_n(z), so
        exp(-alpha*t)*I_n(alpha*gamma*t) = exp(-alpha*(1-gamma)*t) * ive(n, alpha*gamma*t)

    Returns K_n for the given integer n_vals array.
    """
    z      = alpha * gamma * t
    prefac = np.exp(-alpha * (1.0 - gamma) * t)
    return prefac * ive(n_vals, z)


def two_time_bessel_ring(n_val: int, tau_arr: np.ndarray,
                         N: int, gamma: float, eta: float,
                         alpha: float = 1.0,
                         n_windings: int = 3) -> np.ndarray:
    """
    G_n(tau) = sum_l K_{n-l}(tau) * eta^|l|      (infinite chain sum)
    For the ring: replace K with wrapped kernel K^N_n = sum_p K_{n+pN}(tau).

    The sum over l is truncated when ive drops below tol.
    """
    result = np.empty(len(tau_arr))

    for i, tau in enumerate(tau_arr):
        if tau == 0.0:
            r_inf = r_eq_ring(N, eta)
            result[i] = r_inf[n_val % N]
            continue

        z      = alpha * gamma * tau
        prefac = np.exp(-alpha * (1.0 - gamma) * tau)

        # Determine truncation: I_n(z) is negligible for |n| >> z + few*sqrt(z)
        l_max = int(z + 10.0 * np.sqrt(z) + 50)
        l_arr = np.arange(-l_max, l_max + 1, dtype=int)

        # Equilibrium weights eta^|l| (l from -l_max to l_max)
        weights = eta**np.abs(l_arr)

        # Wrapped kernel: sum over windings p in {-n_windings..n_windings}
        G = 0.0
        for p in range(-n_windings, n_windings + 1):
            # K_{n-l+pN}(tau)
            indices = (n_val - l_arr + p * N)
            K_vals  = prefac * ive(indices, z)
            G      += np.dot(K_vals, weights)

        result[i] = G

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Method 4: Kinetic Monte Carlo
# ─────────────────────────────────────────────────────────────────────────────

def kmc_glauber(N: int, beta_J: float, t_max: float,
                n_samples: int = 200,
                t_measure_eq: np.ndarray = None,
                t_measure_2t: np.ndarray = None,
                alpha: float = 1.0,
                rng_seed: int = 42) -> dict:
    """
    Continuous-time (BKL) KMC simulation of the Glauber chain.
    Flip rate: w_j = (alpha/2) * [1 - (gamma/2)*sigma_j*(sigma_{j-1}+sigma_{j+1})]

    Returns dict with:
        r_n_eq  : equal-time correlator (time-averaged after equilibration)
        G_0_tau : autocorrelation G_0(tau) estimated from equilibrated configs
    """
    gamma, eta = beta_to_gamma_eta(beta_J)
    rng        = np.random.default_rng(rng_seed)

    # Storage
    r_accum = np.zeros(N)
    r_count = 0
    G0_tau  = None
    if t_measure_2t is not None:
        G0_tau = np.zeros(len(t_measure_2t))
        G0_cnt = np.zeros(len(t_measure_2t), dtype=int)

    t_equil = t_max * 0.3    # first 30% for equilibration

    for _ in range(n_samples):
        # Random initial configuration
        sigma = rng.choice([-1, 1], size=N)
        t     = 0.0

        # Track a reference snapshot for two-time measurements
        sigma_ref     = None
        t_ref         = None

        while t < t_max:
            # Compute all flip rates
            left  = np.roll(sigma, +1)
            right = np.roll(sigma, -1)
            rates = (alpha / 2.0) * (1.0 - (gamma / 2.0) * sigma * (left + right))
            rates = np.maximum(rates, 0.0)    # numerical safety
            R_tot = rates.sum()

            if R_tot == 0:
                break

            # Time to next event (exponential)
            dt = rng.exponential(1.0 / R_tot)
            t_new = t + dt

            # Equal-time measurement (after equilibration)
            if t >= t_equil and t_new >= t_equil:
                if t_measure_eq is not None:
                    pass    # handled below
                else:
                    # accumulate continuously
                    for dn in range(N):
                        r_accum[dn] += np.mean(sigma * np.roll(sigma, -dn)) * dt
                    r_count += dt

            # Two-time measurement: snapshot at t_equil, then track G_0(tau)
            if sigma_ref is None and t_new >= t_equil:
                sigma_ref = sigma.copy()
                t_ref     = t_equil
            if sigma_ref is not None and t_measure_2t is not None:
                for k, tau in enumerate(t_measure_2t):
                    t_target = t_ref + tau
                    if t <= t_target < t_new:
                        G0_tau[k] += np.mean(sigma * sigma_ref)
                        G0_cnt[k] += 1

            # Choose spin to flip
            j = rng.choice(N, p=rates / R_tot)
            sigma[j] *= -1
            t = t_new

    # Normalize
    if r_count > 0:
        r_eq_kmc = r_accum / r_count
    else:
        r_eq_kmc = np.full(N, np.nan)

    G0_kmc = None
    if t_measure_2t is not None:
        G0_kmc = np.where(G0_cnt > 0, G0_tau / G0_cnt, np.nan)

    return {"r_eq": r_eq_kmc, "G0_tau": G0_kmc, "n_samples": n_samples}


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: dynamic structure factor S(q, omega)
# ─────────────────────────────────────────────────────────────────────────────

def dynamic_structure_factor(q_arr: np.ndarray, omega_arr: np.ndarray,
                              N: int, gamma: float, eta: float,
                              alpha: float = 1.0) -> np.ndarray:
    """
    S(q, omega) = integral e^{i omega tau} G~(q, tau) d tau
                = S_N(q) / [ alpha^2*(1-gamma cos q)^2 + omega^2 ]
                           * alpha*(1-gamma cos q)

    Returns shape (len(q_arr), len(omega_arr)).
    """
    r_inf = r_eq_ring(N, eta)
    n_arr = np.arange(N)

    out = np.empty((len(q_arr), len(omega_arr)))
    for i, q in enumerate(q_arr):
        Sq   = np.real(np.sum(r_inf * np.exp(-1j * q * n_arr)))
        rate = alpha * (1.0 - gamma * np.cos(q))
        # Lorentzian: S(q,omega) = S(q) * rate / (rate^2 + omega^2)
        out[i] = Sq * rate / (rate**2 + omega_arr**2)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Plotting / demo
# ─────────────────────────────────────────────────────────────────────────────

def run_demo(beta_J: float = 1.0, N: int = 64, alpha: float = 1.0):
    gamma, eta = beta_to_gamma_eta(beta_J)
    print(f"beta*J={beta_J:.2f}  gamma={gamma:.4f}  eta={eta:.4f}  xi={-1/np.log(eta):.3f}")

    r_inf = r_eq_ring(N, eta)

    # ── Equal-time evolution (DST) ────────────────────────────────────────────
    t_arr = np.array([0.0, 1.0, 5.0, 20.0, 100.0])
    # start from uncorrelated (all zero for n>=1)
    r0 = np.zeros(N); r0[0] = 1.0
    r_dst = equal_time_dst(r0, t_arr, gamma, eta, alpha)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    ax = axes[0]
    n_plot = np.arange(N // 2 + 1)
    for i, t in enumerate(t_arr):
        ax.plot(n_plot, r_dst[i, n_plot], label=f"t={t:.0f}")
    ax.plot(n_plot, r_inf[n_plot], 'k--', lw=2, label="eq")
    ax.set_xlabel("n"); ax.set_ylabel(r"$r_n(t)$")
    ax.set_title(f"Equal-time correlator (DST), N={N}, βJ={beta_J}")
    ax.legend(fontsize=7)

    # ── Two-time correlator (FFT & Bessel) ────────────────────────────────────
    tau_arr = np.linspace(0, 20, 200)
    G0_fft    = two_time_fft(0, tau_arr, N, gamma, eta, alpha)
    G0_bessel = two_time_bessel_ring(0, tau_arr, N, gamma, eta, alpha)

    ax = axes[1]
    ax.plot(tau_arr, G0_fft,    label="FFT", lw=2)
    ax.plot(tau_arr, G0_bessel, '--', label="Bessel", lw=2)
    ax.set_xlabel(r"$\tau$"); ax.set_ylabel(r"$G_0(\tau)$")
    ax.set_title(f"Autocorrelation (FFT vs Bessel), βJ={beta_J}")
    ax.legend()

    # ── Dynamic structure factor ──────────────────────────────────────────────
    q_vals   = np.array([0.0, np.pi / 4, np.pi / 2, np.pi])
    omega_arr = np.linspace(-5, 5, 500)
    Sqw = dynamic_structure_factor(q_vals, omega_arr, N, gamma, eta, alpha)

    ax = axes[2]
    labels = ["q=0", "q=π/4", "q=π/2", "q=π"]
    for i, lbl in enumerate(labels):
        ax.plot(omega_arr, Sqw[i], label=lbl)
    ax.set_xlabel(r"$\omega$"); ax.set_ylabel(r"$S(q,\omega)$")
    ax.set_title(f"Dynamic structure factor, βJ={beta_J}")
    ax.legend()

    plt.tight_layout()
    plt.savefig("glauber_demo.png", dpi=150)
    print("Saved glauber_demo.png")

    # ── Verify DST vs Bessel for several n ───────────────────────────────────
    print("\nVerification: max |DST - Bessel| for G_n(tau), tau=0..20")
    for n_test in [0, 1, 5]:
        G_fft    = two_time_fft(n_test, tau_arr, N, gamma, eta, alpha)
        G_bessel = two_time_bessel_ring(n_test, tau_arr, N, gamma, eta, alpha)
        print(f"  n={n_test}: max diff = {np.max(np.abs(G_fft - G_bessel)):.2e}")

    # ── Static correlator self-consistency ────────────────────────────────────
    print("\nSelf-consistency: G_n(tau=0) == r_n^eq")
    G_tau0 = np.array([two_time_fft(n, np.array([0.0]), N, gamma, eta, alpha)[0]
                       for n in range(N // 2)])
    print(f"  max |G_n(0) - r_n^eq| = {np.max(np.abs(G_tau0 - r_inf[:N//2])):.2e}")

    # ── Long-time DST: should recover equilibrium ─────────────────────────────
    r_longt = equal_time_dst(r0, np.array([1000.0]), gamma, eta, alpha)[0]
    print(f"\nDST long-time vs eq: max diff = {np.max(np.abs(r_longt - r_inf)):.2e}")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# KMC comparison (separate, slower)
# ─────────────────────────────────────────────────────────────────────────────

def run_kmc_check(beta_J: float = 0.8, N: int = 32, alpha: float = 1.0,
                  n_samples: int = 500):
    gamma, eta = beta_to_gamma_eta(beta_J)
    r_inf      = r_eq_ring(N, eta)
    tau_arr    = np.linspace(0, 10, 30)

    print(f"Running KMC: N={N}, beta*J={beta_J}, {n_samples} samples …")
    kmc = kmc_glauber(N, beta_J, t_max=60.0,
                      n_samples=n_samples,
                      t_measure_2t=tau_arr,
                      alpha=alpha)

    G0_fft = two_time_fft(0, tau_arr, N, gamma, eta, alpha)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    n_plot = np.arange(N // 2 + 1)
    ax.plot(n_plot, r_inf[n_plot], 'b-', lw=2, label="exact")
    ax.plot(n_plot, kmc["r_eq"][n_plot], 'ro', ms=5, label="KMC")
    ax.set_xlabel("n"); ax.set_ylabel(r"$r_n^{\rm eq}$")
    ax.set_title(f"Equilibrium correlator: exact vs KMC")
    ax.legend()

    ax = axes[1]
    ax.plot(tau_arr, G0_fft, 'b-', lw=2, label="FFT (exact)")
    ax.plot(tau_arr, kmc["G0_tau"], 'ro', ms=5, label="KMC")
    ax.set_xlabel(r"$\tau$"); ax.set_ylabel(r"$G_0(\tau)$")
    ax.set_title(f"Autocorrelation: exact vs KMC (n_samples={n_samples})")
    ax.legend()

    plt.tight_layout()
    plt.savefig("glauber_kmc.png", dpi=150)
    print("Saved glauber_kmc.png")

    r_err = np.max(np.abs(kmc["r_eq"][1:N//2] - r_inf[1:N//2]))
    print(f"Max |KMC r_eq - exact|  (n=1..N/2) = {r_err:.4f}")
    valid_G = ~np.isnan(kmc["G0_tau"])
    if valid_G.any():
        G_err = np.max(np.abs(kmc["G0_tau"][valid_G] - G0_fft[valid_G]))
        print(f"Max |KMC G0(tau) - exact|           = {G_err:.4f}")

    return fig


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    run_demo(beta_J=1.0, N=64)

    if "--kmc" in sys.argv:
        run_kmc_check(beta_J=0.8, N=32, n_samples=300)
    else:
        print("\n(pass --kmc to also run the KMC validation, ~30 s)")
