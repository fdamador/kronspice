# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figures for the harmonic_cavitation report, from the CSVs that

    go run ./examples/harmonic_cavitation -viz examples/harmonic_cavitation/viz

writes. python examples/harmonic_cavitation/viz/plot_all.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Wedge, FancyArrow  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
STEEL, STEELDK, WATER, CAV, STEAM, INK, MUTED = "#9aa7ad", "#6f7c82", "#2f7fb8", "#d0612a", "#c43b2c", "#1d2a30", "#5d6d74"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})


def read(name):
    with open(os.path.join(HERE, name), newline="") as f:
        r = csv.reader(f)
        head = next(r)
        return head, [[float(x) for x in row] for row in r]


def col(rows, i):
    return [r[i] for r in rows]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, name), dpi=200)
    plt.close(fig)
    print("wrote", name)


def machine():
    """Section of the heater design: drive in the middle, a dimple drum at each outlet (R 84 mm, 2 rows of
    Ø 6 mm), fixed housing, outlet at each end cap. The Tesla turbine's disks reach the drums' radius, so the
    gap over the drums runs on over the turbine's rim: water entering at each drum's outer top feeds the
    turbine at its periphery. The gap is drawn ×4."""
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    R, dl, gap, hw, drv, shaft = 84.3, 19.2, 3.2, 8, 34, 6
    ax.plot([-(drv + dl + 70), drv + dl + 70], [0, 0], ls=(0, (6, 3, 1, 3)), color=MUTED, lw=.6)
    # drive (motor, or Tesla disk stack) in the middle
    # the turbine's disks reach the drums' radius R: one rim, one gap
    ax.add_patch(Rectangle((-drv, -R), 2 * drv, 2 * R, fc=STEELDK, ec=INK, lw=.6))
    for i in range(7):
        x = -drv + 5 + i * (2 * drv - 10) / 6
        ax.add_patch(Rectangle((x - .8, -R + 4), 1.6, 2 * R - 8, fc=STEEL, ec="none"))
    ax.add_patch(Rectangle((-(drv + dl + 30), -shaft), 2 * (drv + dl + 30), 2 * shaft, fc=INK, ec="none"))
    # the gap continues over the turbine's rim: water from both drums reaches its periphery
    xp = drv + 4
    for y0 in (R, -(R + gap)):
        ax.add_patch(Rectangle((-xp, y0), 2 * xp, gap, fc=WATER, ec="none"))
    # housing shell over the plenum, joining the two drum shells
    ax.add_patch(Rectangle((-xp + 2, R + gap), 2 * xp - 4, hw, fc=STEELDK, ec=INK, lw=.6))
    ax.add_patch(Rectangle((-xp + 2, -R - gap - hw), 2 * xp - 4, hw, fc=STEELDK, ec=INK, lw=.6))
    for s in (-1, 1):
        x0 = drv + 4 if s > 0 else -drv - 4 - dl  # drum's left edge
        # housing: shell over the drum and an end cap with the outlet
        cap = x0 + dl + gap if s > 0 else x0 - gap - hw
        ax.add_patch(Rectangle((x0 - 2, R + gap), dl + 4, hw, fc=STEELDK, ec=INK, lw=.6))
        ax.add_patch(Rectangle((x0 - 2, -R - gap - hw), dl + 4, hw, fc=STEELDK, ec=INK, lw=.6))
        ax.add_patch(Rectangle((cap, shaft + 6), hw, R + gap + hw - shaft - 6, fc=STEELDK, ec=INK, lw=.6))
        ax.add_patch(Rectangle((cap, -(R + gap + hw)), hw, R + gap + hw - shaft - 6, fc=STEELDK, ec=INK, lw=.6))
        # water: the gap over the drum, joined at the drum's inner corner to the plenum round the turbine
        ax.add_patch(Rectangle((x0, R), dl, gap, fc=WATER, ec="none"))
        ax.add_patch(Rectangle((x0, -R - gap), dl, gap, fc=WATER, ec="none"))
        # drum
        ax.add_patch(Rectangle((x0, -R), dl, 2 * R, fc=STEEL, ec=INK, lw=.6))
        for k in range(2):
            x = x0 + dl * (k + .5) / 2
            for y in (R, -R):
                sg = 1 if y > 0 else -1
                ax.add_patch(Wedge((x, y), 3.2, 180 if sg > 0 else 0, 360 if sg > 0 else 180, fc=WATER, ec="none"))
                ax.add_patch(Wedge((x + 3, y + sg * gap), 3.2, 0 if sg > 0 else 180, 180 if sg > 0 else 360, fc=WATER, ec="none"))
                ax.plot(x, y - sg * 1.4, "o", ms=2.6, color=CAV)
        # exhaust steam out along the shaft through the end cap, to the condenser
        xs = cap + (hw if s > 0 else 0)
        ax.add_patch(FancyArrow(xs, shaft + 3, s * 34, 0, width=1.6, head_width=7, head_length=7, color=STEAM))
        ax.text(xs + s * 40, shaft + 12, "STEAM OUT", color=STEAM, ha="left" if s > 0 else "right", va="center", fontsize=8, weight="bold")
        ax.text(xs + s * 40, shaft - 4, "→ condenser (home heat)\n→ pump back to WATER IN", color=MUTED,
                ha="left" if s > 0 else "right", va="top", fontsize=6.5)
        # water in at the drum's outer top corner
        xo = cap + (hw if s > 0 else 0)
        ax.add_patch(FancyArrow(xo + s * 34, R + gap / 2, -s * 30, 0, width=1.4, head_width=7, head_length=6, color=WATER))
        ax.text(xo + s * 38, R + gap / 2, "WATER IN", ha="left" if s > 0 else "right", va="center", fontsize=8, weight="bold")
        # from the plenum into the turbine at its outer rim, toward the axis
        for y in (1, -1):
            ax.add_patch(FancyArrow(s * drv * .5, y * (R + gap / 2), 0, -y * 15, width=1.2, head_width=5,
                                    head_length=4, color=STEAM))
    ax.text(0, R * .3, "TESLA\nTURBINE\n(steam)", ha="center", va="center", fontsize=8, weight="bold", color="white",
            bbox=dict(fc=STEELDK, ec="none", pad=1.5))
    ax.text(0, -R - gap - hw - 4, "one shaft; a motor supplies the rest", ha="center", va="top", fontsize=7, color=MUTED)
    ax.annotate("dimple drum (R 84 mm)", xy=(drv + 4 + dl / 2, -R), xytext=(drv + dl + 50, -R + 6),
                fontsize=7, color=MUTED, arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6))
    ax.annotate("fixed housing", xy=(-drv - 4 - dl / 2, R + gap + hw), xytext=(-drv - dl - 70, R + 22),
                fontsize=7, color=MUTED, arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6))
    ax.set_xlim(-(drv + dl + 150), drv + dl + 150)
    ax.set_ylim(-R - 32, R + 30)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "machine.png")


def calibration():
    _, cb = read("calib.csv")
    _, gc = read("gcurve.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    rpm = col(cb, 0)
    rate = col(cb, 2)
    mx = max(rate)
    a.plot(rpm, [r / mx for r in rate], color=CAV, label="model collapse rate (normalized)")
    for x, lab in ((3000, "no bubbles"), (3400, "bubbles"), (3800, "best")):
        a.axvline(x, color=MUTED, ls="--", lw=.8)
        a.text(x + 30, 0.05, lab, fontsize=7, color=MUTED, rotation=90, va="bottom")
    a.set_xlabel("Xue's rotor speed [rpm]")
    a.set_ylabel("collapse rate, relative")
    a.set_title("calibration on Xue et al. 2024", fontsize=9)
    b.plot(col(gc, 0), col(gc, 1), color=WATER)
    b.axvline(1.105, color=MUTED, ls="--", lw=.8)
    b.text(1.115, 0.9, "inception", fontsize=7, color=MUTED)
    b.axvline(0.775, color=CAV, ls="--", lw=.8)
    b.text(0.785, 0.05, "Xue best", fontsize=7, color=CAV)
    b.set_xlabel("cavitation number σ")
    b.set_ylabel("vapour law g(σ)")
    b.set_title("vapour per cavity, relative", fontsize=9)
    save(fig, "calibration.png")


def front():
    _, fr = read("front.csv")
    _, se = read("sensitivity.csv")
    labels = ["Cd 0.10", "Cd 0.40", "pulse 0.15", "pulse 0.50", "κv 0.10", "St 0.446 only"][:len(se)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    sc = a.scatter(col(fr, 1), col(fr, 0), c=col(fr, 3), cmap="viridis", s=18)
    fig.colorbar(sc, ax=a, label="σ")
    a.set_xlabel("collapse power [W]")
    a.set_ylabel("η (collapse / shaft)")
    a.set_title("Pareto front, shaft ≤ 3 kW", fontsize=9)
    b.bar(range(len(se)) , col(se, 1), color=WATER)
    b.axhline(fr[0][0], color=CAV, ls="--", lw=1)
    b.text(len(se) - .5, fr[0][0] + .01, f"base {fr[0][0]:.3f}", color=CAV, ha="right", fontsize=7)
    b.set_xticks(range(len(se)))
    b.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
    b.set_ylabel("best η")
    b.set_title("declared constants changed", fontsize=9)
    save(fig, "front.png")


def ring():
    _, rg = read("ring.csv")
    fig, a = plt.subplots(figsize=(7.2, 2.4))
    a.semilogy([r[0] / 1e3 for r in rg], [r[1] for r in rg], color=WATER, lw=.8)
    a.set_xlabel("frequency [kHz]")
    a.set_ylabel("|p| at the driven node (scaled)")
    a.set_title("gap ring in kronspice AC, clear water: azimuthal modes m = 1, 2, 3, 4", fontsize=9)
    save(fig, "ring.png")


def collapse():
    _, tr = read("collapse.csv")
    t = [r[0] * 1e6 for r in tr]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.6))
    a.plot(t, [r[1] * 1e3 for r in tr], color=WATER)
    a.set_xlabel("time [µs]")
    a.set_ylabel("cavity radius [mm]")
    a.set_title("Rayleigh–Plesset collapse", fontsize=9)
    b.semilogy(t, [max(abs(r[3]), 1) for r in tr], color=CAV)
    b.set_xlabel("time [µs]")
    b.set_ylabel("|radiated p| at the axis [Pa]")
    b.set_title("one collapse, seen 36 mm away", fontsize=9)
    save(fig, "collapse.png")


def drives():
    _, dr = read("drives.csv")
    fig, a = plt.subplots(figsize=(7.2, 2.5))
    labels = ["motor on the shaft", "pump + Tesla turbine"]
    for i, d in enumerate(dr):
        shaft, water, room, coll = d[2], d[3], d[4], d[5]
        a.barh(i, coll / 1e3, color=CAV, label="collapse energy" if i == 0 else None)
        a.barh(i, (shaft - coll) / 1e3, left=coll / 1e3, color=WATER, label="other shaft work → heat" if i == 0 else None)
        a.barh(i, (water - shaft) / 1e3, left=shaft / 1e3, color="#8fb8d6", label="pump/turbine losses → water" if i == 0 else None)
        a.barh(i, room / 1e3, left=water / 1e3, color=STEEL, label="motor losses → room" if i == 0 else None)
        a.text(d[1] / 1e3 + .2, i, f"{d[1] / 1e3:.1f} kW electric, {d[6]:.3f} collapse/J", va="center", fontsize=7)
    a.set_yticks([0, 1])
    a.set_yticklabels(labels)
    a.set_xlabel("power [kW]")
    a.set_xlim(0, 19)
    a.legend(fontsize=7, frameon=False, loc="lower right")
    save(fig, "drives.png")


def bulb():
    _, pn = read("pern.csv")
    _, best = read("best.csv")
    fig = plt.figure(figsize=(7.2, 2.9))
    for k, n in enumerate((3, 8, 12)):
        _, pr = read(f"bulb_n{n}.csv")
        ax = fig.add_subplot(1, 4, k + 1, projection="polar")
        ax.fill([r[0] for r in pr], [r[1] for r in pr], color=STEEL, ec=INK, lw=.5)
        ax.set_title(f"n = {n}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(0, 1.45)
    a = fig.add_subplot(1, 4, 4)
    a.bar(col(pn, 0), col(pn, 1), color=WATER)
    a.axhline(best[0][1], color=CAV, ls="--", lw=1)
    a.text(16.4, best[0][1] + .005, "dimples", color=CAV, ha="right", fontsize=7)
    a.set_xlabel("n")
    a.set_ylabel("best η")
    a.set_ylim(0, .22)
    save(fig, "bulb.png")


def loop():
    _, le = read("loop_eta.csv")
    _, lp = read("loop.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for pd, c in ((1.18e5, CAV), (1.5e5, MUTED), (3e5, WATER)):
        pts = [(r[1], r[2]) for r in le if r[0] == pd and r[2] == r[2]]
        a.plot([x for x, _ in pts], [y for _, y in pts], color=c, label=f"drum {pd / 1e5:.2f} bar")
    a.set_xlabel("water temperature [°C]")
    a.set_ylabel("η (collapse / shaft)")
    a.set_title("cavitation fades as the water nears boiling", fontsize=9)
    a.legend(fontsize=7, frameon=False)
    labels = ["A: boil in drums", "B: subcooled + flash"]
    for i, r in enumerate(lp):
        a2 = b
        a2.barh(i, r[7] / 1e3, color=CAV, label="collapse" if i == 0 else None)
        a2.barh(i, (r[3] + r[2] - r[7]) / 1e3, left=r[7] / 1e3, color=WATER, label="other drum heat" if i == 0 else None)
        a2.barh(i + 0.38, r[2] / 1e3, height=0.3, color="#8fb8d6", label="turbine returns" if i == 0 else None)
        a2.text(3.15, i, f"{r[4] / 1e3:.2f} kW electric", va="center", fontsize=7)
    b.set_yticks([0.15, 1.15])
    b.set_yticklabels(labels)
    b.set_xlabel("power [kW] (3 kW of shaft work)")
    b.set_xlim(0, 4.3)
    b.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.45, 1.2), ncol=3)
    save(fig, "loop.png")


if __name__ == "__main__":
    machine()
    calibration()
    front()
    ring()
    collapse()
    drives()
    bulb()
    loop()
