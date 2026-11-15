# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Pictures for examples/void_charge, drawn from the Go run (not schematic).

    go run ./examples/void_charge -viz examples/void_charge/viz
    python3 examples/void_charge/viz/plot_charge.py

Reads viz/*.csv and viz/frames/*.bin (frames are not kept in git) and writes
void_charge_graph.png (graph), void_charge_collision.png (image) and
void_charge.mp4 / void_charge.gif (e+e− and e−e− meeting, in real time).
Cuts through the axis are mirrored (−ρ … ρ); the axis is horizontal.
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

HERE = os.path.dirname(os.path.abspath(__file__))
FR = os.path.join(HERE, "frames")
H, LR, LZ = 0.4, 24.0, 24.0


def table(name):
    rows = list(csv.DictReader(open(os.path.join(HERE, name + ".csv"))))
    out = {}
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except ValueError:
            out[k] = [r[k] for r in rows]
    return out


def read_frames(name):
    raw = open(os.path.join(FR, name), "rb").read()
    nr, nz = np.frombuffer(raw[:8], dtype="<i4")
    per = 1 + 2 * nr * nz
    a = np.frombuffer(raw[8:], dtype="<f4").reshape(-1, per)
    return a[:, 0], a[:, 1:1 + nr * nz].reshape(-1, nr, nz), a[:, 1 + nr * nz:].reshape(-1, nr, nz)


def mirror(f):
    return np.vstack([f[::-1], f])


def draw(ax, dens, light, lim=16):
    ax.imshow(mirror(dens), origin="lower", extent=[-LZ, LZ, -LR, LR], cmap="Blues_r", vmin=0, vmax=1.3)
    if light.max() > 0:
        zz = np.arange(dens.shape[1]) * H + H / 2 - LZ
        rr = np.concatenate([-(np.arange(dens.shape[0])[::-1] * H + H / 2), np.arange(dens.shape[0]) * H + H / 2])
        m = light.max()
        ax.contourf(zz, rr, mirror(light), levels=[0.08 * m, 0.35 * m, 0.7 * m, 1.01 * m], cmap="autumn_r", alpha=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim * 0.6, lim * 0.6)
    ax.set_xticks([])
    ax.set_yticks([])


def graph():
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.9))
    e = table("enclosed")
    a = ax[0]
    a.plot(e["r_over_lbar"], e["deficit"], label="charge where ether is missing")
    a.plot(e["r_over_lbar"], e["wall"], label="charge on the wall (Tewari's interface)")
    a.plot(e["r_over_lbar"], e["light"], label="charge with the light")
    a.plot(e["r_over_lbar"], e["qed"], "k--", label="QED (Uehling): rises above e inside ƛ")
    a.axvline(0.866, color="grey", lw=0.8, ls=":")
    a.text(0.9, 0.04, "Darwin rms 0.866 ƛ", color="grey", fontsize=8)
    a.axvline(1.043, color="purple", lw=0.8, ls=":")
    a.text(1.08, 0.12, "Tewari r_e", color="purple", fontsize=8)
    a.set_xscale("log")
    a.set_xlabel("radius r / ƛ")
    a.set_ylabel("charge seen inside r (units of e)")
    a.set_title("Where the charge sits: void bubble vs QED")
    a.legend(fontsize=8, loc="upper left")
    a.set_ylim(0, 1.1)

    p = table("pair")
    a = ax[1]
    a.plot(p["D"], p["e_pm"], "o-", label="e+e− (shared light mode)")
    a.plot(p["D"], p["e_mm"], "s--", label="e−e− (one photon per mode)")
    a.plot(p["D"], -p["coulomb"], "k:", label="Coulomb, opposite charges −αε/D")
    a.plot(p["D"], p["coulomb"], "r:", label="Coulomb, like charges +αε/D")
    a.axhline(0, color="grey", lw=0.6)
    a.set_xlabel("centre separation D (ξ)")
    a.set_ylabel("interaction energy (units of μ)")
    a.set_title("Two voids: short-range pull, no sign, no 1/r")
    a.legend(fontsize=8)
    a.set_yscale("symlog", linthresh=1e-3)

    b = table("breathe")
    q = table("quad")
    a = ax[2]
    a.plot(b["t"], b["local_deficit"], label="missing ether within 7 ξ, rest = 1 (3 % breathing start)")
    a2 = a.twinx()
    a2.plot(q["t"], q["quadrupole"], color="C3", lw=0.8, label="elongation (stretched start)")
    a2.set_ylabel("elongation ∫(1−|ψ|²)(2z²−ρ²)", color="C3")
    a.set_xlabel("time (ħ/μ)")
    a.set_ylabel("charge if carried by missing ether (rest = 1)")
    a.set_title("A charge made of missing ether flickers")
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "void_charge_graph.png"), dpi=130)
    plt.close(fig)


def collision_image():
    t1, d1, x1 = read_frames("rt_pm.bin")
    t2, d2, x2 = read_frames("rt_mm.bin")
    picks = [0, 6, 10, 20, 80]
    fig, ax = plt.subplots(2, len(picks), figsize=(3.4 * len(picks), 4.6))
    for c, k in enumerate(picks):
        k = min(k, len(t1) - 1)
        draw(ax[0, c], d1[k], x1[k])
        draw(ax[1, c], d2[k], x2[k])
        ax[0, c].set_title(f"t = {t1[k]:.0f}")
    ax[0, 0].set_ylabel("e+e−")
    ax[1, 0].set_ylabel("e−e−")
    fig.suptitle("Electron and positron voids meeting: e+e− gather both photons in one void and the other collapses; e−e− coalesce; the light never escapes (blue: ether, orange: light)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "void_charge_collision.png"), dpi=120)
    plt.close(fig)


def movie():
    t1, d1, x1 = read_frames("rt_pm.bin")
    t2, d2, x2 = read_frames("rt_mm.bin")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))

    def frame(k):
        for a in ax:
            a.clear()
        draw(ax[0], d1[k], x1[k])
        draw(ax[1], d2[k], x2[k])
        ax[0].set_title(f"e+e−  t = {t1[k]:.1f}")
        ax[1].set_title(f"e−e−  t = {t2[k]:.1f}")
        return []

    anim = animation.FuncAnimation(fig, frame, frames=min(len(t1), len(t2)), interval=80)
    try:
        anim.save(os.path.join(HERE, "void_charge.mp4"), writer=animation.FFMpegWriter(fps=12))
    except Exception as err:  # ffmpeg missing
        print("mp4 skipped:", err)
    anim.save(os.path.join(HERE, "void_charge.gif"), writer=animation.PillowWriter(fps=12))
    plt.close(fig)


if __name__ == "__main__":
    graph()
    collision_image()
    movie()
