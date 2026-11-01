# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent references for the M4 diode golden files.

The SPICE3 level-1 diode equations are written out again here, in Python,
from the published model (SPICE3f3 User's Manual §4.2; Massobrio & Antognetti,
Semiconductor Device Modeling with SPICE, ch. 1), without sharing code with
the Go implementation. DC points are solved by bisection to machine
precision; transients are integrated with classical fourth-order Runge-Kutta
at a step far below the Go simulator's.

Run from the repository root:  python3 internal/golden/ref/diode_ref.py
It needs only the Python standard library.
"""
import math
import os

K, Q, T = 1.380649e-23, 1.602176634e-19, 300.15
GMIN = 1e-12


class Diode:
    def __init__(self, IS=1e-14, N=1.0, RS=0.0, CJO=0.0, VJ=1.0, M=0.5, TT=0.0, FC=0.5, BV=None, IBV=1e-3):
        self.IS, self.N, self.RS, self.CJO, self.VJ, self.M, self.TT, self.FC, self.BV = IS, N, RS, CJO, VJ, M, TT, FC, BV
        self.vte = N * K * T / Q
        if BV is not None:
            # Breakdown voltage adjusted so the current at BV is IBV
            # (SPICE3f5 diode temperature setup), with vt = kT/q.
            vt = K * T / Q
            if IBV < IS * BV / vt:
                self.xbv = BV
            else:
                g = lambda x: -(IS * (math.exp(min((BV - x) / vt, 700)) - 1 + x / vt) - IBV)
                self.xbv = bisect(g, -100.0, BV)

    def id(self, v):
        """Junction current without gmin, and its derivative."""
        vte, IS = self.vte, self.IS
        if v >= -3 * vte:
            e = math.exp(min(v / vte, 700))
            return IS * (e - 1), IS * e / vte
        if self.BV is None or v >= -self.xbv:
            a = (3 * vte / (v * math.e)) ** 3
            return -IS * (1 + a), 3 * IS * a / v
        e = math.exp(min(-(self.xbv + v) / vte, 700))
        return -IS * e, IS * e / vte

    def cap(self, v):
        """Total small-signal capacitance dq/dv (diffusion + depletion)."""
        _, gd = self.id(v)
        c = self.TT * gd
        if self.CJO:
            fcv = self.FC * self.VJ
            if v < fcv:
                c += self.CJO * (1 - v / self.VJ) ** (-self.M)
            else:
                f2 = (1 - self.FC) ** (1 + self.M)
                f3 = 1 - self.FC * (1 + self.M)
                c += self.CJO / f2 * (f3 + self.M * v / self.VJ)
        return c

    def terminal_current(self, v):
        """Current for terminal voltage v, including RS and gmin."""
        if self.RS == 0:
            return self.id(v)[0] + GMIN * v
        # i = id(vj) + gmin vj with vj = v - RS i; f is monotone in vj.
        f = lambda vj: (v - vj) / self.RS - self.id(vj)[0] - GMIN * vj
        vj = bisect(f, -1e3, 1e3)
        return (v - vj) / self.RS


def bisect(f, lo, hi):
    """Root of a decreasing function f on [lo, hi], to machine precision."""
    flo = f(lo)
    for _ in range(400):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)


def frange(a, b, s):
    n = int(math.floor((b - a) / s + 1e-9)) + 1
    return [a + k * s for k in range(n)]


def write(name, header, rows, note):
    path = os.path.join("testdata", "golden", name)
    with open(path, "w") as fh:
        fh.write("# reference: " + note + " (internal/golden/ref/diode_ref.py)\n")
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(repr(float(x)) for x in r) + "\n")


def diode_fwd():
    d = Diode(IS=1e-14, N=1.05, RS=10)
    V, R = 5.0, 1e3
    va = bisect(lambda va: (V - va) / R - d.terminal_current(va), -10, 10)
    i = (V - va) / R
    path = os.path.join("testdata", "golden", "diode_fwd.csv")
    with open(path, "w") as fh:
        fh.write("# reference: bisection on the SPICE3 diode equations (internal/golden/ref/diode_ref.py)\n")
        fh.write("name,value\n")
        for k, v in [("i(d1)", i), ("i(r1)", i), ("i(v1)", -i), ("v(a)", va), ("v(in)", V)]:
            fh.write("%s,%r\n" % (k, v))


def clipper():
    d = Diode(IS=2.52e-9, N=1.752, RS=0.568)
    R = 1e3
    rows = []
    for vin in frange(-5, 5, 0.25):
        f = lambda vo: (vin - vo) / R - d.terminal_current(vo) + d.terminal_current(-vo)
        vo = bisect(f, -10, 10)
        rows.append((vin, vin, vo, d.terminal_current(vo)))
    write("diode_clipper_dc.csv", ["sweep", "v(in)", "v(out)", "i(d1)"], rows,
          "bisection on the SPICE3 diode equations")


def zener():
    d = Diode(IS=1e-14, N=1, RS=2, BV=5.1, IBV=1e-3)
    R, RL = 470.0, 1e4
    rows = []
    for vin in frange(-3, 12, 0.5):
        f = lambda vo: (vin - vo) / R + d.terminal_current(-vo) - vo / RL
        vo = bisect(f, -20, 20)
        rows.append((vin, vo, d.terminal_current(-vo)))
    write("zener_dc.csv", ["sweep", "v(out)", "i(d1)"], rows, "bisection on the SPICE3 diode equations")


def rk4(f, y, t0, t1, h, record, every):
    n = int(round((t1 - t0) / h))
    out = [(t0, y)]
    for k in range(1, n + 1):
        t = t0 + (k - 1) * h
        k1 = f(t, y)
        k2 = f(t + h / 2, y + h / 2 * k1)
        k3 = f(t + h / 2, y + h / 2 * k2)
        k4 = f(t + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        if k % every == 0:
            out.append((t0 + k * h, y))
    return out


def rectifier():
    d = Diode(IS=1e-14, N=1)
    R, C = 1e3, 10e-6
    vin = lambda t: 10 * math.sin(2 * math.pi * 1e3 * t)
    f = lambda t, vo: (d.terminal_current(vin(t) - vo) - vo / R) / C
    pts = rk4(f, 0.0, 0.0, 3e-3, 1e-8, None, 1000)  # record every 10 us
    write("rectifier_rc.csv", ["time", "v(out)"], [(t, v) for t, v in pts],
          "RK4 integration of the circuit ODE, h = 10 ns")


def charge():
    d = Diode(IS=1e-14, N=1, CJO=10e-12, VJ=0.7, M=0.4, FC=0.5, TT=20e-9)
    R = 1e3
    vin = lambda t: 2 * math.sin(2 * math.pi * 1e6 * t)
    f = lambda t, va: ((vin(t) - va) / R - d.id(va)[0] - GMIN * va) / d.cap(va)
    pts = rk4(f, 0.0, 0.0, 3e-6, 1e-12, None, 2000)  # record every 2 ns
    write("diode_charge_tran.csv", ["time", "v(a)"], [(t, v) for t, v in pts],
          "RK4 integration of the circuit ODE, h = 1 ps")


if __name__ == "__main__":
    diode_fwd()
    clipper()
    zener()
    rectifier()
    charge()
