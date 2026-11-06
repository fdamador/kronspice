# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador

"""Skyrme knots: B = 1 hedgehog (ANW 1983) and B = 4 cube (rational map, Houghton-Manton-Sutcliffe 1998).

Units: lengths in 1/(e F_pi), energies in F_pi/e (ANW). With e = 5.45, F_pi = 129.4 MeV:
length unit hbar c/(e F_pi) = 0.2798 fm, energy unit 23.74 MeV.
Rational-map energy (HMS eq. 2.9, same units after rescaling):
  E = 4pi * Int [ (1/8)(r^2 F'^2 + 2B sin^2F) + (1/2)(2B sin^2F F'^2 + I sin^4F / r^2) ] dr  (ANW form for B = I = 1)
"""
import json
import numpy as np
from scipy.optimize import minimize

HBARC = 197.3269804
E_SK, FPI = 5.45, 129.4
LU = HBARC / (E_SK * FPI)      # fm
EU = FPI / E_SK                # MeV


def grid(n=400, s=1.5):
    x = (np.arange(n) + 0.5) / n          # (0,1)
    r = s * np.tan(np.pi * x / 2)          # (0, inf)
    w = s * np.pi / 2 / np.cos(np.pi * x / 2) ** 2 / n
    return r, w


def energy(F, r, w, B, I):
    # interior nodes F, boundary F(0) = pi, F(inf) = 0 via extended arrays
    Fe = np.concatenate([[np.pi], F, [0.0]])
    re = np.concatenate([[0.0], r, [r[-1] * 50]])
    dF = np.diff(Fe) / np.diff(re)
    rm = 0.5 * (re[1:] + re[:-1])
    Fm = 0.5 * (Fe[1:] + Fe[:-1])
    s2 = np.sin(Fm) ** 2
    dr = np.diff(re)
    dens = (1 / 8) * (rm ** 2 * dF ** 2 + 2 * B * s2) + 0.5 * (2 * B * s2 * dF ** 2 + I * s2 ** 2 / np.maximum(rm, 1e-12) ** 2)
    return 4 * np.pi * np.sum(dens * dr)


def solve(B, I, n=400):
    r, w = grid(n)
    a = np.sqrt(B)
    F0 = 2 * np.arctan((a / r) ** 2) if B > 1 else 2 * np.arctan(1 / r ** 2)
    res = minimize(energy, F0, args=(r, w, B, I), method="L-BFGS-B", options={"maxiter": 20000, "maxfun": 10 ** 7, "ftol": 1e-15, "gtol": 1e-10})
    return r, res.x, res.fun


def baryon_number(r, F):
    Fe = np.concatenate([[np.pi], F, [0.0]])
    return (Fe[0] - Fe[-1]) / np.pi  # hedgehog: B = [F(0) - F(inf)]/pi


if __name__ == "__main__":
    out = {}
    r, F1, E1 = solve(1, 1)
    # moment of inertia integral Lambda = Int r^2 sin^2F [1 + 4(F'^2 + sin^2F/r^2)] dr
    dF = np.gradient(F1, r)
    s2 = np.sin(F1) ** 2
    Lam = np.trapezoid(r ** 2 * s2 * (1 + 4 * (dF ** 2 + s2 / r ** 2)), r)
    # baryon density b(r) = -(1/2pi^2) sin^2F F'/r^2 ; isoscalar radius
    b = -s2 * dF / (2 * np.pi ** 2 * r ** 2)
    Bnum = np.trapezoid(4 * np.pi * r ** 2 * b, r)
    r2 = np.trapezoid(4 * np.pi * r ** 4 * b, r)
    # fix e, F_pi from M_N = 939, M_Delta = 1232 (ANW eqs.): M = E*Fpi/e + 3/(8 lam)*... use ANW closed form
    lam_unit = 2 * np.pi / 3 * Lam   # lambda = lam_unit / (e^3 F_pi)
    # M_N = E Fpi/e + (3/8)/lambda, M_D = E Fpi/e + (15/8)/lambda -> 1/lambda = (M_D - M_N)/1.5
    inv_lam = (1232 - 939) / 1.5
    Mcl = 939 - 3 / 8 * inv_lam
    # Mcl = E1 * Fpi / e and 1/lambda = e^3 Fpi / lam_unit  (in MeV, hbar = c = 1)
    # => Fpi/e = Mcl/E1 ; e^3 Fpi = inv_lam * lam_unit -> e^4 = inv_lam*lam_unit/(Mcl/E1)
    e4 = inv_lam * lam_unit / (Mcl / E1)
    e = e4 ** 0.25
    Fpi = Mcl / E1 * e
    LUf = HBARC / (e * Fpi)
    rad = np.sqrt(r2) * LUf
    MN = 938.9
    lam_fm = lam_unit / (e ** 3 * Fpi)      # MeV^-1
    # ANW moments (in nuclear magnetons, M_N = 939): mu_p,n = M_N <r2>_{I=0}/(9 lam)... ANW eq 34-35:
    # mu_{I=0} = <r^2>_{I=0} M_N / (9 lambda) ... with <r^2> in MeV^-2 ; mu_{I=1} = M_N lambda / 3
    r2_mev = r2 / (e * Fpi) ** 2
    mu_s = r2_mev * MN / (6 * lam_fm)
    mu_v = MN * lam_fm / 3
    mup, mun = mu_s + mu_v, mu_s - mu_v
    out["B1"] = dict(E=E1, Lambda=Lam, Bnum=Bnum, e=e, Fpi=Fpi, radius_fm=rad, mu_p=mup, mu_n=mun, ratio=abs(mup / mun), length_fm=LUf)
    print("B=1: E %.3f (ANW 36.5)  Lambda %.2f (50.9)  B %.5f  e %.3f  Fpi %.1f MeV  r %.3f fm  mu_p %.3f mu_n %.3f |ratio| %.3f" %
          (E1, Lam, Bnum, e, Fpi, rad, mup, mun, abs(mup / mun)))
    # B = 4 cube, rational map I = 20.65
    r4, F4, E4 = solve(4, 20.65)
    hms1, hms4 = E1 / (12 * np.pi ** 2 / 4 * 1) * 1, E4 / (12 * np.pi ** 2 / 4 * 4)
    print("per-winding energy (12 pi^2 units): B=1 %.4f (1.232)  B=4 %.4f (1.137)" % (E1 / (3 * np.pi ** 2), E4 / (4 * 3 * np.pi ** 2)))
    M1, M4 = E1 * Fpi / e, E4 * Fpi / e
    bind = 4 * M1 - M4
    print("classical masses: B=1 %.0f MeV, B=4 %.0f MeV; binding %.0f MeV (%.1f %%) vs measured 28.3 MeV (0.75 %%)" % (M1, M4, bind, 100 * bind / (4 * M1)))
    out["B4"] = dict(E=E4, M1=M1, M4=M4, bind=bind, frac=bind / (4 * M1))
    out["profiles"] = dict(r1=r.tolist(), F1=F1.tolist(), r4=r4.tolist(), F4=F4.tolist(), length_fm=LUf)
    json.dump(out, open("skyrme.json", "w"))
