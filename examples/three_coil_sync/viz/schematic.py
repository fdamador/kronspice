# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Schematic and coil-construction drawings for the three-coil example.

    python3 examples/three_coil_sync/viz/schematic.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz

Writes schematic.png (the circuit kronspice simulates, with the evolved
values from data/design.csv) and coils.png (side section of the three coils
to scale, with the winding data).
"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

INK, MUTED, HVC, ROY, OUT = "#0b0b0b", "#52514e", "#2a78d6", "#eb6834", "#1baf7a"
data, outdir = sys.argv[1], sys.argv[2]


def load_design():
    with open(os.path.join(data, "design.csv")) as f:
        return {r["name"]: float(r["value"]) for r in csv.DictReader(f)}


D = load_design()


def eng(x, unit):
    if x == 0:
        return "0 " + unit
    for p, s in ((1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n"), (1e-12, "p")):
        if abs(x) >= p * 0.9999:
            return f"{x / p:.3g} {s}{unit}"
    return f"{x:.3g} {unit}"


# ---------------------------------------------------------------- primitives
def wire(ax, pts, c=INK, lw=1.4):
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=c, lw=lw, solid_capstyle="round")


def dot(ax, x, y, c=INK):
    ax.add_patch(Circle((x, y), 0.06, color=c, zorder=5))


def resistor(ax, x, y, vertical=True, length=1.0, c=INK, label=None, lpos="r"):
    n, w = 6, 0.15
    pts = []
    for i in range(n + 1):
        t = 0.2 + 0.6 * i / n
        off = 0 if i in (0, n) else (w if i % 2 else -w)
        pts.append((x + off, y - t * length) if vertical else (x + t * length, y + off))
    head = [(x, y), pts[0]] if vertical else [(x, y), pts[0]]
    tail = [pts[-1], (x, y - length)] if vertical else [pts[-1], (x + length, y)]
    wire(ax, head, c), wire(ax, pts, c), wire(ax, tail, c)
    if label:
        tx, ty = (x + 0.3, y - length / 2) if vertical else (x + length / 2, y + 0.3)
        if lpos == "l" and vertical:
            tx = x - 0.3
        ax.text(tx, ty, label, fontsize=7.5, color=c, ha="left" if lpos == "r" else "right", va="center")


def inductor(ax, x, y, vertical=True, length=1.0, n=4, c=INK, label=None, lpos="r", lw=1.4):
    a = 0.8 * length / n / 2
    s = 0.1 * length
    if vertical:
        wire(ax, [(x, y), (x, y - s)], c, lw)
        th = [math.pi / 2 - math.pi * t / 30 for t in range(31)]
        for i in range(n):
            cy = y - s - a - 2 * a * i
            ax.plot([x + a * math.cos(t) for t in th], [cy + a * math.sin(t) for t in th], color=c, lw=lw)
        wire(ax, [(x, y - length + s), (x, y - length)], c, lw)
        if label:
            ax.text(x + (0.35 if lpos == "r" else -0.35), y - length / 2, label, fontsize=7.5, color=c,
                    ha="left" if lpos == "r" else "right", va="center")
    else:
        wire(ax, [(x, y), (x + s, y)], c, lw)
        th = [math.pi - math.pi * t / 30 for t in range(31)]
        for i in range(n):
            cx = x + s + a + 2 * a * i
            ax.plot([cx + a * math.cos(t) for t in th], [y + a * math.sin(t) for t in th], color=c, lw=lw)
        wire(ax, [(x + length - s, y), (x + length, y)], c, lw)
        if label:
            ax.text(x + length / 2, y + 0.35, label, fontsize=7.5, color=c, ha="center", va="bottom")


def capacitor(ax, x, y, vertical=True, length=1.0, c=INK, label=None, lpos="r"):
    g, w = 0.08, 0.25
    if vertical:
        m = y - length / 2
        wire(ax, [(x, y), (x, m + g)], c), wire(ax, [(x, m - g), (x, y - length)], c)
        wire(ax, [(x - w, m + g), (x + w, m + g)], c, 2), wire(ax, [(x - w, m - g), (x + w, m - g)], c, 2)
        if label:
            ax.text(x + (0.35 if lpos == "r" else -0.35), m, label, fontsize=7.5, color=c,
                    ha="left" if lpos == "r" else "right", va="center")
    else:
        m = x + length / 2
        wire(ax, [(x, y), (m - g, y)], c), wire(ax, [(m + g, y), (x + length, y)], c)
        wire(ax, [(m - g, y - w), (m - g, y + w)], c, 2), wire(ax, [(m + g, y - w), (m + g, y + w)], c, 2)
        if label:
            ax.text(m, y + 0.35, label, fontsize=7.5, color=c, ha="center", va="bottom")


def ground(ax, x, y, c=INK):
    wire(ax, [(x, y), (x, y - 0.15)], c)
    for i, w in enumerate((0.25, 0.16, 0.07)):
        wire(ax, [(x - w, y - 0.15 - 0.08 * i), (x + w, y - 0.15 - 0.08 * i)], c)


def source(ax, x, y, kind, c=INK, label=None):
    """Circle source between (x, y) and (x, y − 1), + at the top."""
    wire(ax, [(x, y), (x, y - 0.25)], c), wire(ax, [(x, y - 0.75), (x, y - 1)], c)
    ax.add_patch(Circle((x, y - 0.5), 0.25, fill=False, color=c, lw=1.4))
    if kind == "sin":
        t = [i / 20 for i in range(21)]
        ax.plot([x - 0.15 + 0.3 * u for u in t], [y - 0.5 + 0.09 * math.sin(2 * math.pi * u) for u in t], color=c, lw=1)
    else:
        ax.text(x, y - 0.4, "+", fontsize=8, ha="center", va="center", color=c)
        ax.text(x, y - 0.62, "−", fontsize=8, ha="center", va="center", color=c)
    if label:
        ax.text(x - 0.35, y - 0.5, label, fontsize=7.5, color=c, ha="right", va="center")


def diode(ax, x0, y0, x1, y1, c=INK, label=None):
    """Diode from (x0, y0) anode to (x1, y1) cathode."""
    mx, my = (x0 + x1) / 2, (y0 + y1) / 2
    L = math.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / L, (y1 - y0) / L
    px, py = -uy, ux
    s = 0.17
    wire(ax, [(x0, y0), (mx - ux * s, my - uy * s)], c), wire(ax, [(mx + ux * s, my + uy * s), (x1, y1)], c)
    tri = [(mx - ux * s + px * s, my - uy * s + py * s), (mx - ux * s - px * s, my - uy * s - py * s),
           (mx + ux * s, my + uy * s)]
    ax.add_patch(plt.Polygon(tri, closed=True, fill=True, color=c))
    wire(ax, [(mx + ux * s + px * s, my + uy * s + py * s), (mx + ux * s - px * s, my + uy * s - py * s)], c, 2)
    if label:
        ax.text(mx + px * 0.35, my + py * 0.35, label, fontsize=7, color=c, ha="center", va="center")


def mosfet(ax, x, y, c=INK, label=None, mirror=False):
    """N-MOSFET: drain at (x, y), source at (x, y − 1.2), gate lead to the outside."""
    wire(ax, [(x, y), (x, y - 0.35), (x - 0.25 if not mirror else x + 0.25, y - 0.35)], c)
    wire(ax, [(x, y - 1.2), (x, y - 0.85), (x - 0.25 if not mirror else x + 0.25, y - 0.85)], c)
    bx = x - 0.25 if not mirror else x + 0.25
    wire(ax, [(bx, y - 0.25), (bx, y - 0.95)], c, 2)
    gx = bx - 0.12 if not mirror else bx + 0.12
    wire(ax, [(gx, y - 0.3), (gx, y - 0.9)], c, 1.6)
    g_end = gx - 0.45 if not mirror else gx + 0.45
    wire(ax, [(gx, y - 0.6), (g_end, y - 0.6)], c)
    if label:
        ax.text(x + (0.15 if not mirror else -0.15), y - 0.6, label, fontsize=7.5, color=c,
                ha="left" if not mirror else "right", va="center")
    return (g_end, y - 0.6)


def coupling(ax, x, y0, y1, c=MUTED, label=None):
    for dx in (-0.05, 0.05):
        wire(ax, [(x + dx, y0), (x + dx, y1)], c, 1.0)
    if label:
        ax.text(x, y0 + 0.15, label, fontsize=7, color=c, ha="center", va="bottom")


def box(ax, x0, y0, x1, y1, c, title):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0.05,rounding_size=0.2",
                                fill=False, ec=c, lw=1.0, ls=(0, (4, 3))))
    ax.text(x0 + 0.15, y1 - 0.15, title, fontsize=8.5, color=c, ha="left", va="top", weight="bold")


# ---------------------------------------------------------------- schematic
def schematic(path):
    fig, ax = plt.subplots(figsize=(13, 6.4))
    ax.set_aspect("equal")
    ax.axis("off")

    # ----- coil 1: HV generator and air-core step-up transformer
    box(ax, -0.4, -1.6, 5.2, 5.0, HVC, "Coil 1 — HV step-up (air core)")
    if D.get("Sync", 0) > 0:
        source(ax, 0.6, 3.3, "sin", HVC, f"HV amplifier\n{D['VgPk']:.0f} V pk, {eng(D['FRoyer'], 'Hz')}\n"
                                         f"clocked by coil 2\nphase {D['SyncAngle']:.0f}°, gain {D['SyncGain']:.2f}")
        ax.annotate("", xy=(0.25, 3.05), xytext=(0.25, 5.5), arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.7))
        wire(ax, [(0.25, 5.5), (9.6, 5.5)], MUTED, 0.7)
        ax.text(1.0, 5.58, "sync: pickup loop on coil 2 → phase shifter → HV amplifier input", fontsize=7, color=MUTED)
    else:
        source(ax, 0.6, 3.3, "sin", HVC, f"Vg 48 V pk\n{eng(D['FHV'], 'Hz')}, {D['PhiHV']:.0f}°")
    wire(ax, [(0.6, 3.3), (0.6, 3.8), (1.0, 3.8)], HVC)
    resistor(ax, 1.0, 3.8, vertical=False, length=0.9, c=HVC, label="Rg 1 Ω")
    wire(ax, [(1.9, 3.8), (2.1, 3.8), (2.1, 3.4)], HVC)
    inductor(ax, 2.1, 3.4, length=2.2, n=4, c=HVC, label=f"L1p\n{int(D['N1p'])} t\n{eng(D['L1p'], 'H')}", lpos="l")
    wire(ax, [(2.1, 1.2), (2.1, 0.8)], HVC)
    wire(ax, [(0.6, 2.3), (0.6, 0.8), (2.1, 0.8)], HVC)
    ground(ax, 1.35, 0.8, HVC)
    coupling(ax, 2.85, 3.2, 1.4, MUTED, f"k {D['K1']:.2f}")
    inductor(ax, 3.6, 3.4, length=2.2, n=6, c=HVC, label=f"L1s\n{int(D['N1s'])} t\n{eng(D['L1s'], 'H')}")
    wire(ax, [(3.6, 1.2), (3.6, 0.8)], HVC)
    ground(ax, 3.6, 0.8, HVC)
    wire(ax, [(3.6, 3.4), (3.6, 4.3), (6.6, 4.3)], HVC)
    dot(ax, 3.6, 4.3, HVC)
    ax.text(4.3, 4.45, "HV terminal", fontsize=7.5, color=HVC)
    ax.text(0.0, -0.4, f"self-capacitance C1self {eng(D['C1self'], 'F')} (HV terminal to ground)\n"
                       f"winding resistances R1p {eng(D['R1p'], 'Ω')}, R1s {eng(D['R1s'], 'Ω')}",
            fontsize=7, color=MUTED, va="top")

    # ----- coil 3: resonant output loop (series chain to ground)
    box(ax, 5.6, -1.6, 9.0, 5.0, OUT, "Coil 3 — resonant output")
    inductor(ax, 6.6, 4.3, length=2.0, n=5, c=OUT, label=f"L3  {int(D['N3'])} t\n{eng(D['L3'], 'H')}", lpos="l")
    resistor(ax, 6.6, 2.3, length=0.7, c=OUT, label=f"R3 {eng(D['R3'], 'Ω')}", lpos="l")
    capacitor(ax, 6.6, 1.6, length=0.8, c=OUT, label=f"C3 {eng(D['C3'], 'F')}", lpos="l")
    dot(ax, 6.6, 0.8, OUT)
    ax.text(6.8, 0.8, "out", fontsize=7.5, color=OUT, va="center")
    resistor(ax, 6.6, 0.8, length=1.2, c=OUT, label=f"RL {eng(D['RL'], 'Ω')}", lpos="l")
    ground(ax, 6.6, -0.4, OUT)
    coupling(ax, 7.35, 4.1, 2.5, MUTED, f"k {D['K23']:.2f}")
    wire(ax, [(7.5, 3.3), (10.6, 3.3), (10.6, 3.9)], MUTED, 0.8)
    ax.text(8.2, 3.05, "coil 2 wound\nover coil 3", fontsize=7, color=MUTED, va="top")

    # ----- coil 2: ZVS Royer (Mazzilli), bifilar centre-tapped tank
    box(ax, 9.6, -1.6, 15.6, 6.0, ROY, "Coil 2 — ZVS Royer, bifilar tank")
    yc, xa, xt, xb = 4.0, 11.0, 12.6, 14.2
    inductor(ax, xa, yc, vertical=False, length=xt - xa, n=3, c=ROY)
    inductor(ax, xt, yc, vertical=False, length=xb - xt, n=3, c=ROY)
    ax.text((xa + xt) / 2, yc - 0.45, f"L2a {int(D['N2'])} t", fontsize=7.5, color=ROY, ha="center")
    ax.text((xt + xb) / 2, yc - 0.45, f"L2b {int(D['N2'])} t", fontsize=7.5, color=ROY, ha="center")
    ax.text(9.8, 5.25, f"strands {eng(D['L2strand'], 'H')} each,\nk(strands) {D['K2ab']:.2f}", fontsize=7, color=ROY, va="top")
    dot(ax, xt, yc, ROY)
    # choke from the centre tap up to +24 V
    inductor(ax, xt, yc + 1.3, length=1.3, n=3, c=ROY, label="Lchoke\n100 µH", lpos="r")
    ax.text(xt, yc + 1.45, "+24 V (Vcc, 20 A)", fontsize=7.5, color=ROY, ha="center", va="bottom")
    # drains: c1 and c2 down to the MOSFETs, C2 across them
    wire(ax, [(xa, yc), (xa, 1.6)], ROY)
    wire(ax, [(xb, yc), (xb, 1.6)], ROY)
    capacitor(ax, xa, 3.0, vertical=False, length=xb - xa, c=ROY)
    ax.text(xt, 3.25, f"C2 {eng(D['C2'], 'F')}", fontsize=7.5, color=ROY, ha="center", va="bottom")
    dot(ax, xa, 3.0, ROY), dot(ax, xb, 3.0, ROY)
    ax.text(xa - 0.12, 3.6, "c1", fontsize=7, color=ROY, ha="right"), ax.text(xb + 0.12, 3.6, "c2", fontsize=7, color=ROY)
    g1 = mosfet(ax, xa, 1.6, ROY, "M1", mirror=True)
    g2 = mosfet(ax, xb, 1.6, ROY, "M2", mirror=False)
    ground(ax, xa, 0.4, ROY), ground(ax, xb, 0.4, ROY)
    ax.text(xa - 0.2, 0.55, "IRFP260N", fontsize=6.5, color=ROY, ha="right")
    ax.text(xb + 0.2, 0.55, "IRFP260N", fontsize=6.5, color=ROY, ha="left")
    # crossed diodes: D1 g1 → c2, D2 g2 → c1
    dot(ax, *g1, ROY), dot(ax, *g2, ROY)
    dot(ax, xa, 2.4, ROY), dot(ax, xb, 2.4, ROY)
    diode(ax, g1[0], g1[1], xb, 2.4, ROY)
    diode(ax, g2[0], g2[1], xa, 2.4, ROY)
    ax.text(xt, 2.25, "D1, D2 UF4007", fontsize=7, color=ROY, ha="center")
    # gate pull-ups from the 12 V rail
    resistor(ax, g1[0], g1[1], length=1.4, c=ROY)
    resistor(ax, g2[0], g2[1], length=1.4, c=ROY)
    wire(ax, [(g1[0], g1[1] - 1.4), (g2[0], g2[1] - 1.4)], ROY)
    ax.text(xt, g1[1] - 1.3, "2 × 470 Ω", fontsize=7, color=ROY, ha="center", va="bottom")
    ax.text(xt, g1[1] - 1.55, "+12 V gate rail", fontsize=7.5, color=ROY, ha="center", va="top")
    ax.text(9.8, -1.1, "gate C 3 nF, C_oss 1 nF per MOSFET (in the model)", fontsize=6.8, color=MUTED)

    ax.set_xlim(-0.6, 15.9)
    ax.set_ylim(-1.9, 6.2)
    fig.suptitle("Circuit simulated by kronspice — evolved values", fontsize=11, x=0.02, ha="left")
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- coils
def coils(path):
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.25]})
    # coils 2 and 3: side section
    ax = axs[0]
    ax.set_aspect("equal")
    r3 = 25.0 + 0.5
    n3, n2 = int(D["N3"]), int(D["N2"])
    l3 = n3 * 1.1
    gap = D["gap23"] * 1e3
    r2 = 25.0 + 1.0 + gap + 1.0
    l2 = 2 * n2 * 2.2
    z2 = 0.5 * (l3 - l2)
    ax.add_patch(Rectangle((-25, -5), 50, l3 + 10, fill=True, fc="#eeeeec", ec=MUTED, lw=0.8))
    ax.text(0, l3 + 7, "Ø50 mm former (PVC / PTFE)", ha="center", fontsize=8, color=MUTED)
    for s in (-1, 1):
        for i in range(n3):
            ax.add_patch(Circle((s * r3, i * 1.1 + 0.55), 0.5, color=OUT))
        if gap > 0.5:
            ax.add_patch(Rectangle((s * (26.0 if s > 0 else 26.0 + gap) if s > 0 else -(26.0 + gap), z2 - 3), gap, l2 + 6,
                                   fc="#f6f1e7", ec=MUTED, lw=0.5))
        for i in range(2 * n2):
            ax.add_patch(Circle((s * r2, z2 + i * 2.2 + 1.1), 1.0, color=ROY if i % 2 == 0 else "#b8461d"))
    ax.set_xlim(-r2 - 22, r2 + 22)
    ax.set_ylim(-12, l3 + 14)
    ax.axis("off")
    ax.set_title("Coils 3 and 2 (section, to scale, mm)", fontsize=10, loc="left")
    ax.annotate("", xy=(r2 + 8, z2), xytext=(r2 + 8, z2 + l2), arrowprops=dict(arrowstyle="<->", color=ROY, lw=0.8))
    ax.text(r2 + 10, z2 + l2 / 2, f"coil 2: {n2}+{n2} t bifilar\nØ2.0 mm, Ø{2 * r2:.0f} mm\nlength {l2:.0f} mm",
            fontsize=7.5, color=ROY, va="center")
    ax.annotate("", xy=(-r3 - 8 - gap, 0), xytext=(-r3 - 8 - gap, l3), arrowprops=dict(arrowstyle="<->", color=OUT, lw=0.8))
    ax.text(-r3 - 10 - gap, l3 / 2, f"coil 3: {n3} t\nØ1.0 mm, Ø{2 * r3:.0f} mm\nlength {l3:.0f} mm",
            fontsize=7.5, color=OUT, va="center", ha="right")
    ax.text(0, -9, f"spacer {gap:.1f} mm between the windings · k(2–3) {D['K23']:.2f}", ha="center", fontsize=7.5,
            color=MUTED)

    # coil 1: HV transformer
    ax = axs[1]
    ax.set_aspect("equal")
    n1s, n1p = int(D["N1s"]), int(D["N1p"])
    l1 = n1s * 0.45
    r1 = 37.5
    rp = r1 + 10
    lp = n1p * 3.0
    ax.add_patch(Rectangle((-r1 + 0.3, -5), 2 * r1 - 0.6, l1 + 10, fc="#eeeeec", ec=MUTED, lw=0.8))
    for s in (-1, 1):
        ax.add_patch(Rectangle((s * r1 - 0.2 if s > 0 else -r1 - 0.2, 0), 0.4, l1, color=HVC))
        for i in range(n1p):
            ax.add_patch(Circle((s * rp, i * 3.0 + 1.5), 1.0, color="#1d5ea8"))
    ax.set_xlim(-rp - 60, rp + 60)
    ax.set_ylim(-14, max(l1, lp) + 14)
    ax.axis("off")
    ax.set_title("Coil 1 — HV transformer (section, to scale, mm)", fontsize=10, loc="left")
    ax.text(r1 + 3, l1 * 0.6, f"secondary {n1s} t, Ø0.4 mm\nclose-wound on Ø75 mm\nlength {l1:.0f} mm\n"
                              f"{eng(D['L1s'], 'H')}; top = HV terminal → coil 3", fontsize=7.5, color=HVC, va="center")
    ax.text(-rp - 3, lp / 2 + 2, f"primary {n1p} t, Ø2.0 mm\non Ø95 mm, 3 mm pitch\nat the grounded end\n"
                                 f"{eng(D['L1p'], 'H')}, k {D['K1']:.2f}", fontsize=7.5, color="#1d5ea8", ha="right",
            va="center")
    ax.text(0, -10, "keep coil 1 ≥ 3 diameters away from coils 2–3 (no stray coupling is modelled)",
            ha="center", fontsize=7.5, color=MUTED)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)


schematic(os.path.join(outdir, "schematic.png"))
coils(os.path.join(outdir, "coils.png"))
