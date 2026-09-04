"""
theory_continuum.py — N -> infinity scaling limits of the equilibrium Glauber
correlation functions (the closed forms of glauber_paper, Sec. "Scaling limits
of the correlation functions", ss:macform).

Macroscopic variables: space x = site/N, time tau = micro/N^2.  All formulas are
on the line R (the L-torus corrections are O(e^{-2L}) and tanh(L)~1 for L>=4;
use H_torus for the exact periodic version).

    H(x,tau)  = <sigma_0(0) sigma_{Nx}(N^2 tau)>  (e:Herfc)
              = e^{-4tau} E[e^{-2|x+W_2tau|}]
              = 1/2 [ e^{-2|x|} erfc((4tau-|x|)/2sqrt(tau))
                    + e^{ 2|x|} erfc((4tau+|x|)/2sqrt(tau)) ]
    H(x,0)    = e^{-2|x|},   H(0,tau) = erfc(2 sqrt(tau))           (e:auto)
    P(x,tau)  = e^{-4tau} e^{-x^2/4tau} / sqrt(4 pi tau)  (killed heat kernel)
    d_x H     = sgn(x) [ e^{2|x|} erfc(a+) - e^{-2|x|} erfc(a-) ],
                a± = (4tau±|x|)/(2 sqrt(tau))   (heat-kernel terms cancel)
    S         = -1/4 d_x H
    f(x,tau)  = 4 H P - 4 H^2 + (d_x H)^2 = 4[HP - H^2 + 4 S^2]      (e:f)
                interface memory; f(x,0)=0, f != 0 for tau>0.

box_average() convolves any of these with the tent kernel of a width-h box
(the law of the difference of two uniforms on [-h/2,h/2]), matching the box
observable m_box of BoxCorr / mc_box_correlation.py.
"""
import numpy as np
from scipy.special import erfc

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))   # numpy 2.x rename

__all__ = ["H", "P", "dxH", "S", "f", "f_torus", "autocorr", "box_average",
           "H_torus", "static_line"]


def _ax(x):
    return np.atleast_1d(np.asarray(x, dtype=float))


def H(x, tau):
    """Two-time spin correlator H(x,tau) on the line (e:Herfc)."""
    x = _ax(x); ax = np.abs(x)
    if tau <= 0:
        return np.exp(-2.0 * ax)
    s = 2.0 * np.sqrt(tau)
    am = (4.0 * tau - ax) / s
    ap = (4.0 * tau + ax) / s
    return 0.5 * (np.exp(-2.0 * ax) * erfc(am) + np.exp(2.0 * ax) * erfc(ap))


def P(x, tau):
    """Heat kernel of d_x^2 killed at rate 4 (e:HSP)."""
    x = _ax(x)
    if tau <= 0:
        return np.zeros_like(x)
    return np.exp(-4.0 * tau) * np.exp(-x * x / (4.0 * tau)) / np.sqrt(4.0 * np.pi * tau)


def dxH(x, tau):
    """d_x H(x,tau) (closed form; the 1/sqrt(pi tau) terms cancel)."""
    x = _ax(x); ax = np.abs(x)
    if tau <= 0:
        return -2.0 * np.sign(x) * np.exp(-2.0 * ax)
    s = 2.0 * np.sqrt(tau)
    am = (4.0 * tau - ax) / s
    ap = (4.0 * tau + ax) / s
    return np.sign(x) * (np.exp(2.0 * ax) * erfc(ap) - np.exp(-2.0 * ax) * erfc(am))


def S(x, tau):
    return -0.25 * dxH(x, tau)


def f(x, tau):
    """Interface memory function: the VALIDATED continuum limit of
    N^2 ( <eta_0(0) eta_{Nx}(N^2 tau)> - p^2 ),

        f(x,tau) = H P - H^2 + (1/4) (d_xH)^2 .

    This is 1/4 of the paper's boxed e:f (= 4HP - 4H^2 + (d_xH)^2 = f_paper).
    Both the box Monte Carlo and the exact finite-N thesis-45 correlator
    (theory_glauber.H_interface) converge to THIS, not to f_paper: the ratio
    f_paper/(N^2 box-conn) -> 4.00 across N and tau, and the MC chi^2/dof ~ 1
    against f/4 (986/147/8.7 against f). This resolves the '1/4' prefactor
    flagged at eq (4.14) of the paper."""
    h = H(x, tau); pp = P(x, tau); d = dxH(x, tau)
    return h * pp - h * h + 0.25 * d * d


def f_paper(x, tau):
    """The paper's boxed e:f = 4HP - 4H^2 + (d_xH)^2 = 4*f(x,tau) (4x too large;
    kept only to display the discrepancy)."""
    return 4.0 * f(x, tau)


def autocorr(tau):
    """Autocorrelation H(0,tau) = erfc(2 sqrt(tau)) (e:auto)."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    return erfc(2.0 * np.sqrt(np.clip(tau, 0, None)))


def static_line(x):
    """Equal-time macroscopic static correlator H(x,0) = e^{-2|x|}."""
    return np.exp(-2.0 * np.abs(_ax(x)))


def H_torus(xi, tau, Lper, n_img=6):
    """Periodic (length-Lper torus) version: tanh(Lper) sum_m H(xi+m Lper, tau).
    For Lper>=4 this equals H(xi,tau) to O(e^{-2 Lper}); tanh(Lper)~1."""
    xi = _ax(xi)
    acc = np.zeros_like(xi)
    for m in range(-n_img, n_img + 1):
        acc = acc + H(xi + m * Lper, tau)
    return np.tanh(Lper) * acc


def f_torus(xi, tau, Lper, n_img=8):
    """Periodic (length-Lper torus) interface memory: sum_m f(xi + m Lper, tau).

    The line f decays like the spin correlator, so on a ring of circumference
    Lper the two points at the antipode (sep = Lper/2) are correlated via BOTH
    arcs: the m=+/-1 images are ~equal to the m=0 term there and roughly DOUBLE
    f.  This is the interface analogue of H_torus for the spin correlator; near
    sep=0 the images are O(e^{-2 Lper}) and f_torus ~ f.  The ring-normalisation
    prefactor (analogue of H_torus's tanh(Lper)) is 1 - O(e^{-2 Lper}) ~ 1 for
    Lper>=4 and is omitted.  Validated: at h=2 (box = Lper/2) the box-average of
    f_torus tracks the ring MC at the antipode, whereas the bare line f is ~2x
    too low (up to 12 sigma)."""
    xi = _ax(xi)
    acc = np.zeros_like(xi)
    for m in range(-n_img, n_img + 1):
        acc = acc + f(xi + m * Lper, tau)
    return acc


def box_average(func, x, tau, h, nquad=401):
    """Box-average of func(., tau) over a width-h box: convolve with the tent
    kernel tri(s) = (h-|s|)/h^2 on [-h,h] (difference of two Unif[-h/2,h/2]).
    func is called as func(x_shifted, tau) and must accept an array."""
    x = _ax(x)
    s = np.linspace(-h, h, nquad)
    tent = np.maximum(h - np.abs(s), 0.0) / (h * h)
    vals = np.array([func(x + ss, tau) for ss in s])     # (nquad, len(x))
    return _trapz(tent[:, None] * vals, s, axis=0)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import theory_glauber as th

    print("== self-tests ==")
    taus = np.array([0.05, 0.1, 0.2, 0.4])

    # 1. autocorrelation matches erfc(2 sqrt(tau)) and H(0,tau)
    a1 = autocorr(taus); a2 = np.array([H(0.0, t)[0] for t in taus])
    print(f"1. autocorr vs H(0,tau): max|diff| = {np.max(np.abs(a1-a2)):.2e}")
    print(f"   H(0,tau) = {a2}")

    # 2. PDE residual d_tau H = d_xx H - 4 H  (central differences)
    x0, dt, dx = 0.7, 1e-4, 1e-4
    dtauH = (H(x0, 0.2 + dt)[0] - H(x0, 0.2 - dt)[0]) / (2 * dt)
    dxxH = (H(x0 + dx, 0.2)[0] - 2 * H(x0, 0.2)[0] + H(x0 - dx, 0.2)[0]) / dx**2
    print(f"2. PDE residual at (x=0.7,tau=0.2): "
          f"{abs(dtauH - (dxxH - 4*H(x0,0.2)[0])):.2e}")

    # 3. d_xH closed form vs finite difference
    fd = (H(x0 + dx, 0.2)[0] - H(x0 - dx, 0.2)[0]) / (2 * dx)
    print(f"3. dxH closed vs FD at (0.7,0.2): {abs(dxH(x0,0.2)[0]-fd):.2e}")

    # 4. f(x,0+) -> 0
    print(f"4. f(0.5, 1e-6) = {f(0.5, 1e-6)[0]:.2e}  (should ~0)")

    # 5. CONVENTION: box-avg line H vs box-avg exact finite ring at large N (L=4)
    L, h = 4, 0.5
    for N in (100, 200, 400):
        l = N * L; eta = th.eta_of(N); gamma = th.gamma_of(N)
        w = int(round(h * N / 2)); c0 = l // 2
        centers = np.round(np.arange(15) * l / 15).astype(int) % l
        offs = np.arange(-w, w + 1); sep = offs[:, None] - offs[None, :]
        # exact ring, box-averaged, at tau=0.2
        t = 2.0 * N * N * 0.2
        gfull = th.G_two_time_ring(np.arange(l), t, l, eta, gamma)
        ring = np.array([gfull[(cj - c0 + sep) % l].mean() for cj in centers])
        # continuum L-torus H, box-averaged at the same macroscopic separations
        xsep = ((centers - c0 + l // 2) % l - l // 2) / N
        cont_line = box_average(H, xsep, 0.2, h)
        cont_tor = box_average(lambda xx, tt: H_torus(xx, tt, L), xsep, 0.2, h)
        print(f"5. N={N:3d}: max|box-ring - box-H_torus| = {np.max(np.abs(ring-cont_tor)):.2e}"
              f"   (line-H would be off by {np.max(np.abs(ring-cont_line)):.2e})")
    print("done.")
