"""sqrt(v) law from the continuum kernel:
E[(Y_{t+v}-Y_t)^2] ~ ||f||^2 [kappa - 2 * int C_v(w) dw]  (f smooth),
C_v = -K11 K22 + K12 K21 = (k^2/4)[(Psi-)^2 - (Psi+)^2] + (k/2) g Psi+ .
"""
import numpy as np
from scipy.special import erfc
from scipy.integrate import quad
from math import sqrt, pi, exp

kap = 2.0
chat = kap ** 2 / 2
def Pp(v, w):
    s = sqrt(2 * v)
    return 0.5 * (np.exp(-kap * w) * erfc((kap * v - w) / s) + np.exp(kap * w) * erfc((kap * v + w) / s))
def Pm(v, w):
    s = sqrt(2 * v)
    return 0.5 * (np.exp(-kap * w) * erfc((kap * v - w) / s) - np.exp(kap * w) * erfc((kap * v + w) / s))
def g(v, w):
    return exp(-chat * v) * exp(-w * w / (2 * v)) / sqrt(2 * pi * v)
def C(v, w):
    return (kap ** 2 / 4) * (Pm(v, w) ** 2 - Pp(v, w) ** 2) + (kap / 2) * g(v, w) * Pp(v, w)

print(f"{'v':>8} {'kappa-2*intC':>14} {'/sqrt(v)':>10} {'/v':>10}")
for v in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
    I = quad(lambda w: C(v, w), -12, 12, points=[-5 * sqrt(v), 0, 5 * sqrt(v)], limit=800)[0]
    val = kap - 2 * I
    print(f"{v:8.0e} {val:14.8f} {val/sqrt(v):10.5f} {val/v:10.2f}")
print("prediction if only the same-particle term mattered: 2*(k/2)*(1-erfc(k sqrt(v/2)))/sqrt(v) -> k*k*sqrt(2/pi) =", kap * kap * sqrt(2 / pi))
