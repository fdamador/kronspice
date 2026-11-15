# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Part 13 pictures: the full ether simulation (Go run, not schematic).

    go run ./examples/ether_horizon -part13 -viz examples/ether_horizon/viz
    python3 examples/ether_horizon/viz/plot_ether13.py

Reads viz/part13_runs.csv and viz/part13_frames/*.bin (written by the Go run;
the frames are not kept in git) and writes ether13_graph.png,
ether13_snapshots.png, ether13.mp4 and ether13.gif. Cuts through the axis are
mirrored (−ρ … ρ) so the donut shows both sides.
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

HERE = os.path.dirname(os.path.abspath(__file__))
FR = os.path.join(HERE, "part13_frames")
H, LR, LZ = 0.4, 20.0, 26.0
AE = 0.00115965218059
COL = {"thin, ε = 100": "#9ecae1", "thin, ε = 300": "#1f77b4", "thin, ε = 1000": "#d62728",
       "half, ε = 300": "#ff7f0e", "horn, ε = 300": "#2ca02c"}


def runs():
    out = {}
    for r in csv.DictReader(open(os.path.join(HERE, "part13_runs.csv"))):
        d = out.setdefault(r["case"], {k: [] for k in r if k != "case"})
        for k, v in r.items():
            if k == "case":
                continue
            d[k].append(1.0 if v == "true" else 0.0 if v == "false" else float(v))
    return {k: {a: np.array(b) for a, b in v.items()} for k, v in out.items()}


def frames(name):
    raw = open(os.path.join(FR, name), "rb").read()
    nr, nz = np.frombuffer(raw[:8], dtype="<i4")
    per = 1 + 2 * nr * nz
    a = np.frombuffer(raw[8:], dtype="<f4").reshape(-1, per)
    return a[:, 0], a[:, 1:1 + nr * nz].reshape(-1, nr, nz), a[:, 1 + nr * nz:].reshape(-1, nr, nz)


def graph(R):
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    for k, d in R.items():
        c = COL.get(k, "k")
        t = d["t"]
        ax[0, 0].plot(t, d["inner"], color=c, lw=1.8, label=k)
        ax[0, 1].plot(t, d["z_c"], color=c, lw=1.8, label=k)
        ax[1, 0].plot(t, d["winding"].clip(-2.2, 2.2) + 0.03 * list(R).index(k), color=c, lw=1.4, label=k)
        ax[1, 1].plot(t, d["anomaly"] / AE, color=c, lw=1.6, label=k)
    ax[0, 0].axhline(0.2, color="k", lw=0.6, ls=":")
    ax[0, 0].text(1, 0.6, "axis (hole closed when the inner edge reaches it)", fontsize=8)
    ax[0, 0].set(xlabel="time [ξ/c_s units]", ylabel="inner edge of the light's void, ρ [ξ]",
                 title="(a) the hole: thin donuts keep it open; strong light, the half-fat\ndonut and the horn close it (inner edge reaches the axis)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].set(xlabel="time", ylabel="z of the void's centre",
                 title="(b) every donut that keeps its vortex travels (~⅓ of sound speed)")
    ax[1, 0].set(xlabel="time", ylabel="phase winding round the void",
                 title="(c) winding round the light's void: thin donuts keep 1; once the void\nsits on the axis the count jumps (the vortex has left the light)")
    ax[1, 0].set_ylim(-2.5, 2.5)
    ax[1, 1].axhline(1, color="m", ls=":", lw=1.2, label="measured")
    ax[1, 1].axhline(0.959, color="gray", ls="--", lw=1, label="Part 12 (thin-ring formulas, horn torus)")
    ax[1, 1].set(xlabel="time", ylabel="(g − 2)/2 from the wall charge / measured", ylim=(0, 1.1),
                 title="(d) the anomaly from the computed walls\n(means 0.28–0.67 of measured; Part 12 said 0.959)")
    ax[1, 1].legend(fontsize=7, loc="lower right")
    fig.suptitle("H49 Part 13: full ether simulation (Gross–Pitaevskii ether + the m = ½ light), real time, no friction",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "ether13_graph.png")
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print("wrote", out)


def mirror(f):
    return np.vstack([f[::-1], f])  # ρ from −LR to LR


CASES = [("frames_thin_e300.bin", "thin donut, ε = 300"), ("frames_thin_e1000.bin", "thin donut, strong light ε = 1000"),
         ("frames_horn.bin", "Part 12's horn torus, ε = 300")]


def draw(axh, dens, chi2, zc=None):
    ext = [-LZ, LZ, -LR, LR]
    axh.imshow(mirror(dens), origin="lower", extent=ext, cmap="Blues_r", vmin=0, vmax=1.25, aspect="equal")
    lvl = chi2.max() * np.array([0.08, 0.35, 0.7])
    zz = np.arange(dens.shape[1]) * H + H / 2 - LZ
    rr = np.concatenate([-(np.arange(dens.shape[0])[::-1] * H + H / 2), np.arange(dens.shape[0]) * H + H / 2])
    axh.contourf(zz, rr, mirror(chi2), levels=list(lvl) + [chi2.max() * 1.01], cmap="autumn_r", alpha=0.75)
    axh.axhline(0, color="w", lw=0.5, ls=":")
    axh.set_xlim(-LZ, LZ)
    axh.set_ylim(-LR, LR)
    axh.set_xticks([])
    axh.set_yticks([])


def snapshots():
    F = [frames(f) for f, _ in CASES]
    times = [0, 10, 20, 40, 60]
    fig, ax = plt.subplots(len(CASES), len(times), figsize=(17, 7.6))
    for r, ((t, D, C), (_, lab)) in enumerate(zip(F, CASES)):
        for c, tt in enumerate(times):
            k = int(np.argmin(abs(t - tt)))
            draw(ax[r, c], D[k], C[k])
            ax[r, c].set_title(f"{lab}\nt = {t[k]:.0f}" if c == 0 else f"t = {t[k]:.0f}", fontsize=9)
    fig.text(0.5, 0.01, "Cuts through the axis (axis horizontal, white dotted; z left–right, ρ up–down, mirrored). Blue: ether "
             "density (white = empty). Orange–red: the light. Sound ripples fill the box after the start.", ha="center",
             fontsize=9)
    fig.suptitle("H49 Part 13: what the donut void does in the full ether simulation", fontsize=12)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    out = os.path.join(HERE, "ether13_snapshots.png")
    fig.savefig(out, dpi=110)
    plt.close(fig)
    print("wrote", out)


def movie():
    F = [frames(f) for f, _ in CASES]
    n = min(len(f[0]) for f in F)
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    title = fig.suptitle("", fontsize=11)

    def frame(k):
        for a, (t, D, C), (_, lab) in zip(ax, F, CASES):
            a.cla()
            draw(a, D[k], C[k])
            a.set_title(lab, fontsize=10)
        title.set_text(f"H49 Part 13 — full ether simulation, real time, no friction:  t = {F[0][0][k]:.0f}   "
                       "(blue: ether, white: empty, orange: light; axis horizontal)")
        return []

    anim = animation.FuncAnimation(fig, frame, frames=n, interval=120)
    out = os.path.join(HERE, "ether13.mp4")
    anim.save(out, writer=animation.FFMpegWriter(fps=8, bitrate=2400))
    print("wrote", out)
    out = os.path.join(HERE, "ether13.gif")
    anim.save(out, writer=animation.PillowWriter(fps=8), dpi=55)
    print("wrote", out)
    plt.close(fig)


if __name__ == "__main__":
    graph(runs())
    snapshots()
    movie()
