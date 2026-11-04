# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Schematic of the stacked dual-motor machine with one RLC tank between the
two machines' output coils (dual.go, Output = OutSeries; output_tank.go).

One current runs shaft → disk 1 → R_int → coil 1 → tank → coil 2 → R_int →
disk 2 → shaft. The two shafts are the common (ground); shaft 2 counter-
rotates and coil 2 is wound to aid coil 1, so the two disks' EMFs add round
the loop. The tank is R_L ∥ C ∥ L_p, L_p the primary of one transformer whose
secondary has its common end tied to the shafts' common and is rectified by a
bridge (DC side floating) into a capacitor bank for DC loads. Values: the
first-pass design of output_tank.go, to be optimized for a commercial device.

Usage: python viz/plot_schematic.py [out_dir]   (schemdraw, matplotlib)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import schemdraw
import schemdraw.elements as elm

here = os.path.dirname(os.path.abspath(__file__))
out = sys.argv[1] if len(sys.argv) > 1 else here
TEAL, ORANGE, QUIET = "#127b8a", "#d9772b", "#5b5b66"
FS = 9
Y0, YG2, YB = 0.0, -5.4, -11.0  # tank top (g1), tank bottom (g2), common (the shafts)
X1, X2 = 2.0, 6.6               # the machines' vertical branches
XR, XC, XP, XS = 10.0, 12.0, 14.6, 16.0  # R_L, C, primary, secondary


def disk_branch(d, n, x, ytop, col, motor_left, vfd):
    """Shaft (common) → disk → R_int, up to ytop."""
    d.add(elm.Line().at((x, YB)).to((x, YB + 0.6)).linewidth(3.5).color(QUIET))
    d.add(elm.Label().at((x + (0.15 if motor_left else -0.15), YB + 0.3)).label(
        f"shaft {n}", fontsize=8, color=QUIET, halign="left" if motor_left else "right"))
    d.add(elm.SourceV().at((x, YB + 0.6)).to((x, YB + 3.1)).color(col))
    d.add(elm.Resistor().at((x, YB + 3.1)).to((x, YB + 5.3)).label(
        "R_int\n0.27 mΩ", loc="bot" if motor_left else "top", fontsize=FS))
    d.add(elm.Line().at((x, YB + 5.3)).to((x, ytop)))
    side, ha = (-0.35, "right") if motor_left else (0.35, "left")
    d.add(elm.Label().at((x - side * 1.7, YB + 1.85)).label(  # inner side, away from the motor
        f"disk {n}", fontsize=FS, color=col, halign="left" if motor_left else "right"))
    d.add(elm.Label().at((x - side * 0.45, YB + 3.35)).label(
        f"rim {n}", fontsize=8, color=QUIET, halign="left" if motor_left else "right"))
    xm = x - 3.0 if motor_left else x + 3.0
    d.add(elm.Motor().at((xm, YB + 0.6)).up(2.5).label(f"motor {n}", loc="top" if motor_left else "bot", fontsize=FS))
    xe, xs = (xm + 0.55, x - 0.6) if motor_left else (xm - 0.55, x + 0.6)
    d.add(elm.Line().at((xe, YB + 1.85)).to((xs, YB + 1.85)).linestyle("--").color(QUIET))
    d.add(elm.Label().at((xm, YB - 0.55)).label(vfd, fontsize=8))


with schemdraw.Drawing(show=False) as d:
    d.config(fontsize=10, unit=2.5)
    d.add(elm.Label().at((-1.5, 3.0)).label(
        "Stacked dual-motor homopolar machine: ONE parallel RLC tank between the two machines' output coils\n"
        "values: the first-pass design, to be optimized for a commercial device",
        fontsize=12, halign="left"))
    # machine 1: up to the top rail, coil 1 along it into the tank's top (g1)
    disk_branch(d, 1, X1, Y0, TEAL, True, "VFD 1: +ω0·(1 ± sin 2π·7.8t)")
    d.add(elm.Inductor2(loops=4).at((X1, Y0)).to((X1 + 3.0, Y0)).color(TEAL)
          .label("coil 1, 43 µH, around disk 2", fontsize=FS, color=TEAL))
    d.add(elm.Line().at((X1 + 3.0, Y0)).to((XP, Y0)))
    # machine 2: up to the tank's bottom rail (g2), coil 2 along it
    disk_branch(d, 2, X2, YG2, ORANGE, False, "VFD 2: −ω0·(1 ± sin 2π·7.8t), counter-rotating")
    d.add(elm.Inductor2(loops=4).at((X2, YG2)).to((X2 + 3.0, YG2)).color(ORANGE)
          .label("coil 2, 43 µH, around disk 1\n(wound to aid coil 1)", fontsize=FS, color=ORANGE))
    d.add(elm.Line().at((X2 + 3.0, YG2)).to((XP, YG2)))
    # the common: both shafts, grounded
    d.add(elm.Line().at((X1, YB)).to((X2, YB)).linewidth(3.5).color(QUIET))
    d.add(elm.Dot().at((X1, YB)))
    d.add(elm.Dot().at((X2, YB)))
    d.add(elm.Ground().at((4.3, YB)))
    d.add(elm.Label().at((4.3, YB - 1.25)).label("common: shaft 1 = shaft 2", fontsize=9, color=QUIET))
    # the one tank: R_L ∥ C ∥ L_p between coil 1's end (g1) and coil 2's end (g2)
    for x in (XR, XC, XP):
        d.add(elm.Dot().at((x, Y0)))
        d.add(elm.Dot().at((x, YG2)))
    d.add(elm.Resistor().at((XR, Y0)).to((XR, YG2)).label("R_L\n127 mΩ", loc="bot", fontsize=FS))
    d.add(elm.Capacitor().at((XC, Y0)).to((XC, YG2)).label("C\n4.4 F", loc="bot", fontsize=FS))
    d.add(elm.Inductor2(loops=5).at((XP, Y0)).to((XP, YG2)))
    d.add(elm.Label().at((XP - 0.4, (Y0 + YG2) / 2)).label("L_p\n433 µH", fontsize=FS, halign="right"))
    d.add(elm.Label().at(((XR + XP) / 2, Y0 + 0.45)).label("tank", fontsize=9, color=QUIET))
    for dx in (0.6, 0.8):
        d.add(elm.Line().at((XP + dx, Y0 - 0.2)).to((XP + dx, YG2 + 0.2)))
    d.add(elm.Inductor2(loops=5).at((XS, YG2)).to((XS, Y0)))
    d.add(elm.Label().at(((XP + XS) / 2, Y0 + 0.9)).label("1 : 12.5, k = 0.98", fontsize=FS))
    # secondary: its common end on the shafts' common (ground), the other end to the bridge
    xsn, xa, xb, xbank, xdc = XS + 1.2, XS + 4.0, XS + 5.6, XS + 7.4, XS + 9.2
    yhi, ylo = Y0 - 1.0, YG2 + 1.0
    ya, yb2 = Y0 - 2.4, YG2 + 2.4
    d.add(elm.Line().at((XS, Y0)).to((xa - 0.7, Y0)))
    d.add(elm.Line().at((xa - 0.7, Y0)).to((xa - 0.7, ya)))
    d.add(elm.Line().at((xa - 0.7, ya)).to((xa, ya)))
    d.add(elm.Line().at((XS, YG2)).to((xb + 0.7, YG2)))
    d.add(elm.Line().at((xb + 0.7, YG2)).to((xb + 0.7, yb2)))
    d.add(elm.Line().at((xb + 0.7, yb2)).to((xb, yb2)))
    d.add(elm.Dot().at((xsn, Y0)))
    d.add(elm.Dot().at((xsn, YG2)))
    d.add(elm.Dot().at((XS + 2.6, YG2)))
    d.add(elm.Ground().at((XS + 2.6, YG2)))
    d.add(elm.Label().at((XS + 2.9, YG2 - 0.95)).label("secondary common = shafts' common", fontsize=8, color=QUIET, halign="left"))
    d.add(elm.Resistor().at((xsn, Y0)).to((xsn, YG2)).color(QUIET).label("snubber\n20 Ω + 50 µF", loc="bot", fontsize=8))
    for x, ym in ((xa, ya), (xb, yb2)):
        d.add(elm.Diode().at((x, ylo)).to((x, ym)))
        d.add(elm.Diode().at((x, ym)).to((x, yhi)))
        d.add(elm.Dot().at((x, ym)))
    d.add(elm.Line().at((xa, yhi)).to((xdc, yhi)))
    d.add(elm.Line().at((xa, ylo)).to((xdc, ylo)))
    d.add(elm.Label().at(((xa + xb) / 2, yhi + 0.4)).label("bridge", fontsize=FS))
    d.add(elm.Label().at((xdc + 0.25, yhi)).label("+", fontsize=11, halign="left"))
    d.add(elm.Label().at((xdc + 0.25, ylo)).label("−", fontsize=11, halign="left"))
    for x in (xa, xb, xbank, xdc):
        d.add(elm.Dot().at((x, yhi)))
        d.add(elm.Dot().at((x, ylo)))
    d.add(elm.Capacitor(polar=True).at((xbank, yhi)).to((xbank, ylo)).label("bank\n25 mF", loc="bot", fontsize=FS))
    d.add(elm.Resistor().at((xdc, yhi)).to((xdc, ylo)).label("DC load\n40 Ω", loc="bot", fontsize=FS))
    d.add(elm.Label().at((-1.5, YB - 2.3)).label(
        "one current I: shaft 1 → disk 1 → R_int → coil 1 → tank → coil 2 → R_int → disk 2 → shaft 2 → common → shaft 1;"
        " each coil's field drives the other disk (cross-excitation, M) and reaches its own disk (K_S = 0.52)\n"
        "the bridge's DC side floats (only the secondary's common end is tied to the shafts); ground symbols are one node",
        fontsize=8, color=QUIET, halign="left"))
    d.save(os.path.join(out, "dual_motor_schematic.png"), dpi=200)
print("wrote dual_motor_schematic.png")
