# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Circuit schematics for examples/magnetic_gas, drawn with matplotlib.

    python3 examples/magnetic_gas/viz/schematic.py

Writes schematic_marx.png and schematic_pump.png beside this script. Part
names and node names are those of the netlists in export/ (circuits.go).
"""
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyBboxPatch, Polygon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, SURF, ACC = "#0b0b0b", "#52514e", "#fcfcfb", "#2a78d6"
LW = 1.6


class Sheet:
    def __init__(self, w, h):
        self.fig, self.ax = plt.subplots(figsize=(w, h))
        self.fig.set_facecolor(SURF)
        self.ax.set_facecolor(SURF)
        self.ax.set_aspect("equal")
        self.ax.axis("off")

    def line(self, pts, color=INK, lw=LW, ls="-"):
        xs, ys = zip(*pts)
        self.ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="round")

    def dot(self, p):
        self.ax.add_patch(Circle(p, 0.07, color=INK, zorder=5))

    def text(self, p, s, ha="center", va="center", size=9, color=INK, **kw):
        self.ax.text(p[0], p[1], s, ha=ha, va=va, fontsize=size, color=color, **kw)

    def ground(self, p):
        x, y = p
        self.line([(x, y), (x, y - 0.25)])
        for i, w in enumerate((0.3, 0.2, 0.1)):
            yy = y - 0.25 - 0.1 * i
            self.line([(x - w, yy), (x + w, yy)])

    # Two-terminal parts drawn from p to q; the symbol (length 1) sits in
    # the middle and leads fill the rest.
    def part(self, p, q, kind, label="", value="", side=1, color=INK):
        (x0, y0), (x1, y1) = p, q
        L = math.hypot(x1 - x0, y1 - y0)
        ux, uy = (x1 - x0) / L, (y1 - y0) / L
        vx, vy = -uy, ux
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2

        def P(u, v):
            return (mx + u * ux + v * vx, my + u * uy + v * vy)

        half = 0.5
        self.line([p, P(-half, 0)], color=color)
        self.line([P(half, 0), q], color=color)
        if kind == "R":
            n = 6
            pts = [P(-half, 0)]
            for i in range(n):
                pts.append(P(-half + (i + 0.5) / n, 0.15 * (1 if i % 2 == 0 else -1)))
            pts.append(P(half, 0))
            self.line(pts, color=color)
        elif kind in ("C", "CV"):
            self.line([P(-half, 0), P(-0.08, 0)], color=color)
            self.line([P(0.08, 0), P(half, 0)], color=color)
            for u in (-0.08, 0.08):
                self.line([P(u, -0.3), P(u, 0.3)], color=color, lw=2.2)
            if kind == "CV":  # variable: arrow through the plates
                a, b = P(-0.32, -0.38), P(0.32, 0.38)
                self.ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))
        elif kind == "L":
            n, r = 4, 0.125
            for i in range(n):
                cu = -half + r + 2 * r * i
                th = [math.pi * k / 20 for k in range(21)]
                self.line([P(cu - r * math.cos(t), r * math.sin(t)) for t in th], color=color)
        elif kind == "GAP":  # spark gap: two pointed electrodes
            self.line([P(-half, 0), P(-0.28, 0)], color=color)
            self.line([P(0.28, 0), P(half, 0)], color=color)
            self.ax.add_patch(Polygon([P(-0.28, 0.16), P(-0.28, -0.16), P(-0.07, 0)], closed=True, color=color))
            self.ax.add_patch(Polygon([P(0.28, 0.16), P(0.28, -0.16), P(0.07, 0)], closed=True, color=color))
        elif kind == "D":  # anode at p, cathode at q
            self.line([P(-half, 0), P(-0.2, 0)], color=color)
            self.line([P(0.2, 0), P(half, 0)], color=color)
            self.ax.add_patch(Polygon([P(-0.2, 0.22), P(-0.2, -0.22), P(0.2, 0)], closed=True, color=color))
            self.line([P(0.2, -0.22), P(0.2, 0.22)], color=color, lw=2.2)
        elif kind == "V":  # + at q
            self.line([P(-half, 0), P(-0.3, 0)], color=color)
            self.line([P(0.3, 0), P(half, 0)], color=color)
            self.ax.add_patch(Circle((mx, my), 0.3, fill=False, ec=color, lw=LW))
            self.text(P(0.14, 0), "+", size=9)
            self.text(P(-0.14, 0), "−", size=9)
        if label:
            off = 0.42 * side
            lp = P(0, off)
            ha = "center"
            if abs(vx) > 0.7:  # vertical part: label to the side
                ha = "left" if off * vx > 0 else "right"
            txt = label if not value else f"{label}\n{value}"
            self.text(lp, txt, ha=ha, size=8.5, linespacing=1.15)

    def save(self, name, title, sub):
        self.ax.set_title(title, loc="left", fontsize=13, color=INK, pad=24)
        self.ax.text(0, 1.015, sub, transform=self.ax.transAxes, fontsize=9, color=INK2)
        self.fig.tight_layout()
        self.fig.savefig(os.path.join(HERE, name), dpi=150, facecolor=SURF, bbox_inches="tight", pad_inches=0.25)


def marx():
    s = Sheet(12, 6.2)
    N, dx, x0 = 4, 2.6, 2.6
    yH, yL = 3.0, 0.0
    X = [x0 + k * dx for k in range(N)]
    # supply into h1
    s.part((-0.6, yL), (-0.6, yH), "V", "VS", "9.63 kV\n(0.9 × breakdown)", side=-1)
    s.ground((-0.6, yL))
    s.part((-0.6, yH), (X[0], yH), "R", "RS", "1 MΩ")
    for k, x in enumerate(X):
        n = k + 1
        s.part((x, yL), (x, yH), "C", f"C{n}", "10 nF", side=-1)
        s.dot((x, yH))
        s.text((x, yH + 0.3), f"h{n}", size=8, color=INK2)
        if n > 1:
            s.dot((x, yL))
            s.text((x, yL - 0.32), f"l{n}", size=8, color=INK2)
        if n < N:
            xn = X[k + 1]
            s.part((x, yH), (xn, yH), "R", f"RH{n}", "100 kΩ")
            s.part((x, yL), (xn, yL), "R", f"RL{n}", "100 kΩ", side=-1)
            s.part((x, yH), (xn, yL), "GAP", f"RG{n}", "", side=1, color=ACC)
    s.ground((X[0], yL))
    # output gap, loop inductance and load
    xo = X[-1] + 2.4
    s.part((X[-1], yH), (xo - 1.2, yH), "GAP", "RG4", "", color=ACC)
    s.dot((xo - 1.2, yH))
    s.text((xo - 1.2, yH + 0.3), "x", size=8, color=INK2)
    s.part((xo - 1.2, yH), (xo + 0.6, yH), "L", "LS", "1 µH (loop)")
    s.dot((xo + 0.6, yH))
    s.text((xo + 0.6, yH + 0.3), "out", size=8, color=INK2)
    s.part((xo + 0.6, yH), (xo + 0.6, yL), "R", "RLOAD", "50 Ω", side=-1)
    s.ground((xo + 0.6, yL))
    # legend for the gaps
    s.text((X[0] - 0.2, -1.25), "blue: spark gaps RG1–RG4, 0.2 mm, filled with Ar, H₂, O₂ or HHO at 500 psi "
           "(open while charging; 0.5 Ω arcs once fired)", ha="left", size=8.5, color=INK2)
    s.text((X[0] - 0.2, -1.65), "charging: C1–C4 fill in parallel through RS and the 100 kΩ ladders · firing: the gaps "
           "erect C1–C4 in series into LS and RLOAD", ha="left", size=8.5, color=INK2)
    s.ax.set_xlim(-2.6, xo + 2.4)
    s.ax.set_ylim(-1.9, 3.6)
    s.save("schematic_marx.png", "Four-stage Marx generator with gas spark gaps",
           "Node and part names as in export/marx_argon_charge.cir and marx_argon_fire.cir; values for argon")


def pump():
    s = Sheet(10, 5.4)
    yT, yB = 3.0, 0.0
    s.part((0, yB), (0, yT), "V", "VIN", "100 V", side=-1)
    s.ground((0, yB))
    s.text((0, yT + 0.3), "in", size=8, color=INK2)
    s.dot((0, yT))
    s.part((0, yT), (2.4, yT), "D", "D1", "")
    s.dot((2.4, yT))
    s.text((2.4, yT + 0.3), "v", size=8, color=INK2)
    s.part((2.4, yB), (2.4, yT), "CV", "CV  gas gap", "10 nF; ×3 with\nplasma lit", side=1, color=ACC)
    s.ground((2.4, yB))
    s.part((2.4, yT), (5.0, yT), "D", "D2", "")
    s.dot((5.0, yT))
    s.text((5.0, yT + 0.3), "out", size=8, color=INK2)
    s.part((5.0, yB), (5.0, yT), "C", "COUT", "100 nF", side=1)
    s.line([(5.0, yT), (6.8, yT)])
    s.part((6.8, yB), (6.8, yT), "R", "RLOAD", "1 MΩ", side=-1)
    s.line([(5.0, yB), (6.8, yB)])
    s.ground((5.0, yB))
    s.dot((5.0, yB))
    # RF supply that lights the plasma
    box = FancyBboxPatch((-0.9, -1.75), 4.9, 0.85, boxstyle="round,pad=0.05,rounding_size=0.12",
                         fc=SURF, ec=INK2, lw=1.0, ls=(0, (4, 3)))
    s.ax.add_patch(box)
    s.text((1.55, -1.33), "RF supply that lights the plasma (not modelled):\npays the work that changes CV (EXT=electric-field)",
           size=8, color=INK2)
    s.ax.annotate("", xy=(2.7, 1.1), xytext=(2.9, -0.88),
                  arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0, ls=(0, (4, 3))))
    s.text((4.2, -1.33), "D1, D2: Si diodes, IS = 1e-14 A, RS = 0.1 Ω", ha="left", size=8, color=INK2)
    s.ax.set_xlim(-2.6, 9.2)
    s.ax.set_ylim(-2.0, 3.8)
    s.save("schematic_pump.png", "Plasma-varactor charge pump",
           "Fill CV at 100 V while the plasma is lit (C high), let it go out (C ×⅓, V ×3), D2 hands the charge to COUT")


if __name__ == "__main__":
    marx()
    pump()
