# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Pictures for examples/void_shells, drawn from the Go run (not schematic).

    go run ./examples/void_shells -viz examples/void_shells/viz
    python3 examples/void_shells/viz/plot_shells.py

Reads viz/states.csv and viz/frames/*.bin (the frames are not kept in git) and
writes shells_g.png (graph), shells_gallery.png (image) and shells.mp4 /
shells.gif (the prepared shell and the round bubble in real time). Cuts
through the axis are mirrored (−ρ … ρ); the axis is horizontal.
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


def states():
    return list(csv.DictReader(open(os.path.join(HERE, "states.csv"))))


def read_frames(name):
    raw = open(os.path.join(FR, name), "rb").read()
    nr, nz = np.frombuffer(raw[:8], dtype="<i4")
    per = 1 + 2 * nr * nz
    a = np.frombuffer(raw[8:], dtype="<f4").reshape(-1, per)
    return a[:, 0], a[:, 1:1 + nr * nz].reshape(-1, nr, nz), a[:, 1 + nr * nz:].reshape(-1, nr, nz)


def mirror(f):
    return np.vstack([f[::-1], f])


def draw(ax, dens, light, lim=14):
    ax.imshow(mirror(dens), origin="lower", extent=[-LZ, LZ, -LR, LR], cmap="Blues_r", vmin=0, vmax=1.25)
    if light.max() > 0:
        zz = np.arange(dens.shape[1]) * H + H / 2 - LZ
        rr = np.concatenate([-(np.arange(dens.shape[0])[::-1] * H + H / 2), np.arange(dens.shape[0]) * H + H / 2])
        m = light.max()
        ax.contourf(zz, rr, mirror(light), levels=[0.08 * m, 0.35 * m, 0.7 * m, 1.01 * m], cmap="autumn_r", alpha=0.8)
    ax.axhline(0, color="w", lw=0.5, ls=":")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xticks([])
    ax.set_yticks([])


def graph(S):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    groups = [("1", "round bubble vs ε", "eps"), ("2", "one photon, mode k", "mode"), ("4", "N photons, one per mode", "N")]
    for a, (part, title, kind) in zip(ax, groups):
        rows = [r for r in S if r["part"] == part]
        xs = []
        for r in rows:
            n = r["name"]
            xs.append(float(n.split("=")[1]) if "=" in n else float(n.split()[-1]))
        g = [float(r["g"]) for r in rows]
        a.plot(xs, g, "o-", color="#1f77b4", lw=2, ms=6, label="computed g = M/E_L")
        a.axhline(2, color="k", ls="--", lw=1, label="g = 2 (Dirac)")
        a.axhline(4 / 3, color="#d62728", ls=":", lw=1.2, label="4/3: held by pressure")
        a.axhline(1.5, color="#ff7f0e", ls=":", lw=1.2, label="3/2: held by a wall")
        a.set_ylim(1.15, 2.1)
        if kind == "eps":
            a.set_xscale("log")
            a.set_xlabel("light strength ε (plays ħc)")
        else:
            a.set_xlabel("mode index k" if kind == "mode" else "photons N")
        voids = [int(r["voids"]) for r in rows]
        for x, y, v in zip(xs, g, voids):
            a.annotate(f"{v} void", (x, y), (0, 8), textcoords="offset points", ha="center", fontsize=7, color="#555")
        a.set_ylabel("g")
        a.set_title(title)
        a.legend(fontsize=7, loc="upper right")
    fig.suptitle("void_shells (H50): every state is one void, and g stays at 1.31–1.37, never near 2",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "shells_g.png")
    fig.savefig(out, dpi=125)
    plt.close(fig)
    print("wrote", out)


def gallery():
    names = [("state_ball.bin", "1 photon, lowest mode:\nthe round bubble"),
             ("state_mode3.bin", "1 photon in mode 3:\nthe void stretches, no shell"),
             ("state_mode6.bin", "1 photon in mode 6"),
             ("state_N4.bin", "4 photons, one per mode:\none shared void"),
             ("state_shell_start.bin", "prepared shell (start):\nether core in a void shell"),
             ("state_onion_start.bin", "prepared onion (start):\nvoid, ether layer, void")]
    fig, ax = plt.subplots(1, len(names), figsize=(19, 4.3))
    for a, (f, lab) in zip(ax, names):
        _, D, C = read_frames(f)
        draw(a, D[0], C[0], lim=15)
        a.set_title(lab, fontsize=9)
    fig.text(0.5, 0.02, "Cuts through the axis (axis horizontal, mirrored). Blue: ether density (white = empty). Orange–red: the "
             "light. The two prepared starts on the right do not survive: relaxed, both become the round bubble.",
             ha="center", fontsize=9)
    fig.suptitle("void_shells: the states the light and the ether settle into", fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    out = os.path.join(HERE, "shells_gallery.png")
    fig.savefig(out, dpi=115)
    plt.close(fig)
    print("wrote", out)


def movie():
    A = read_frames("rt_shell.bin")
    B = read_frames("rt_ball.bin")
    n = min(len(A[0]), len(B[0]))
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.6))
    title = fig.suptitle("")

    def frame(k):
        for a, (t, D, C), lab in zip(ax, (A, B), ("prepared shell: ether core inside a void shell", "the round bubble")):
            a.cla()
            draw(a, D[k], C[k], lim=16)
            a.set_title(lab, fontsize=10)
        title.set_text(f"void_shells, real time, no friction: t = {A[0][k]:.1f}   (blue: ether, white: empty, orange: light)")
        return []

    anim = animation.FuncAnimation(fig, frame, frames=n, interval=80)
    out = os.path.join(HERE, "shells.mp4")
    anim.save(out, writer=animation.FFMpegWriter(fps=12, bitrate=2000))
    print("wrote", out)
    out = os.path.join(HERE, "shells.gif")
    anim.save(out, writer=animation.PillowWriter(fps=10), dpi=60)
    print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    graph(states())
    gallery()
    movie()
