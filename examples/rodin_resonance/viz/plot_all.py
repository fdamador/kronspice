# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Figures for examples/rodin_resonance (H43), from the CSV files the Go
program writes to viz/:

    python3 examples/rodin_resonance/viz/plot_all.py

  geometry.png  the {12/5} Rodin pair (A, B) on the former, top view and 3-D,
                beside the (2, 1) winding of photon_electron's class, and the
                evolved output groups
  fields_xy.png, fields_xz.png
                field maps in the mid-plane (z = 0) and the vertical plane
                y = 0 for four drives: |B| peak over a cycle, rms, the
                time-averaged energy density and the ellipticity (how much the
                field vector turns), with the direction of the major axis
  energy.png    magnetic energy over the last cycles of the 1 kHz drives;
                B(t) at points; the 100 kHz build-ups
  map.png       open-circuit EMF and load power of a reference coil against
                (frequency, angle round the toroid)
  modes.png     self-resonances: open, series, plain and crossed rings, for the
                Rodin (5, 12) and the (2, 1) windings, coloured by trip parity
  bank.png      the winner's bank voltage and energy against time, and the GA
                history
  toroidal.png  Part 3t: the 12-section toroidal output winding on the pair
                (top view and cross-section), the EMF of one turn round the
                tube against angle and against (frequency, angle), and the
                winners' bank charging and GA histories

Sequential maps use one hue (light to dark); categorical lines use the
reference palette's fixed order (dataviz skill).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SEQ = LinearSegmentedColormap.from_list("seq", ["#fcfcfb", "#cfe0f6", "#7fb0ea", "#2a78d6", "#123f78"])
plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": False, "figure.dpi": 150, "axes.titlesize": 8.5})

R0 = (0.0995 / 2 + 0.029 / 2) / 2
A0 = (0.0995 / 2 - 0.029 / 2) / 2


def rows(name):
    with open(os.path.join(HERE, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def have(name):
    return os.path.exists(os.path.join(HERE, name))


def tor_sections():
    """The Part 3t winding (the in-phase winner if present, else the sin/cos one)."""
    for tag in ("_phase", ""):
        if have(f"tor_sections{tag}.csv"):
            return rows(f"tor_sections{tag}.csv")
    return []


def draw_tor3d(ax, P, sec, every=1):
    """The pair and the toroidal output winding in 3-D: every turn a circle round
    the tube in its meridian plane, through the hole and over the outside."""
    ax.plot(*(P["A"].T * 100), color=SER[0], lw=0.8)
    ax.plot(*(P["B"].T * 100), color=SER[1], lw=0.8)
    u = np.linspace(0, 2 * np.pi, 48)
    for x in sec:
        c, arc, a = np.radians(float(x["centre_deg"])), np.radians(float(x["arc_deg"])), float(x["radius"])
        n = int(x["turns"])
        col = SER[2] if float(x["sign"]) > 0 else SER[4]
        for j in range(0, n, every):
            ph = c - arc / 2 + arc * (j + 0.5) / n
            r = R0 + a * np.cos(u)
            ax.plot(r * np.cos(ph) * 100, r * np.sin(ph) * 100, a * np.sin(u) * 100, color=col, lw=0.5, alpha=0.8)
    ax.view_init(elev=32, azim=-62)
    try:
        ax.set_box_aspect((1, 1, 0.45), zoom=1.35)
    except TypeError:  # matplotlib < 3.6
        ax.set_box_aspect((1, 1, 0.45))
    ax.set_axis_off()


def geometry():
    r = rows("paths.csv")
    P = {}
    for x in r:
        P.setdefault(x["winding"], []).append((float(x["x"]), float(x["y"]), float(x["z"])))
    P = {k: np.array(v) for k, v in P.items()}
    sec = tor_sections()
    ncol = 4 if sec else 3
    fig = plt.figure(figsize=(3 * ncol, 3.4))
    ax = fig.add_subplot(1, ncol, 1)
    th = np.linspace(0, 2 * np.pi, 200)
    for rr in (R0 - A0, R0 + A0):
        ax.plot(rr * np.cos(th) * 100, rr * np.sin(th) * 100, color=GRID, lw=1)
    ax.plot(P["A"][:, 0] * 100, P["A"][:, 1] * 100, color=SER[0], lw=1, label="A (1→6→11…)")
    ax.plot(P["B"][:, 0] * 100, P["B"][:, 1] * 100, color=SER[1], lw=1, label="B (mirror, 1→8→3…)")
    for k in range(12):
        a = 2 * np.pi * k / 12
        ax.text(5.6 * np.cos(a), 5.6 * np.sin(a), str(k + 1), ha="center", va="center", color=INK2, fontsize=7)
    for tag, col, lab in (("", SER[2], "free search"), ("_dist", SER[6], "distributed groups")):
        if not have(f"groups{tag}.csv"):
            continue
        for g in rows(f"groups{tag}.csv"):
            t = np.radians(float(g["theta_deg"]))
            rho, rad, beta = float(g["rho"]), float(g["radius"]), np.radians(float(g["beta_deg"]))
            # the coil's outline seen from above: an ellipse of half-axes rad (along ρ) and rad·sin β
            a = np.linspace(0, 2 * np.pi, 64)
            u = rad * np.cos(a)
            w = rad * np.sin(a) * np.sin(beta)
            x = rho * np.cos(t) + u * np.cos(t) - w * np.sin(t)
            y = rho * np.sin(t) + u * np.sin(t) + w * np.cos(t)
            ax.plot(x * 100, y * 100, color=col, lw=1.2)
        ax.plot([], [], color=col, label=f"Part 3 output coil(s), {lab}")
    ax.set_aspect("equal")
    ax.set_title("Rodin pair, top view: the {12/5} star")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, frameon=False)
    ax = fig.add_subplot(1, ncol, 2, projection="3d")
    ax.plot(*(P["A"].T * 100), color=SER[0], lw=0.8)
    ax.plot(*(P["B"].T * 100), color=SER[1], lw=0.8)
    ax.view_init(elev=35, azim=-60)
    ax.set_title("3-D: 5 times round, 12 times through")
    ax.set_box_aspect((1, 1, 0.45))
    ax.set_axis_off()
    ax = fig.add_subplot(1, ncol, 3, projection="3d")
    ax.plot(*(P["K21A"].T * 100), color=SER[6], lw=1)
    ax.view_init(elev=35, azim=-60)
    ax.set_title("(2, 1) winding (photon_electron's class)")
    ax.set_box_aspect((1, 1, 0.45))
    ax.set_axis_off()
    if sec:
        ax = fig.add_subplot(1, ncol, 4, projection="3d")
        draw_tor3d(ax, P, sec, every=4)
        ax.set_title(f"Part 3t: toroidal output winding\n12 sections × {sec[0]['turns']} turns (every 4th drawn)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "geometry.png"))
    plt.close(fig)


def load_planes():
    r = rows("field_planes.csv")
    out = {}
    for pl in ("xy", "xz"):
        sub = [x for x in r if x["plane"] == pl]
        n = int(round(np.sqrt(len(sub))))
        pts = np.array([[float(x["x"]), float(x["y"]), float(x["z"])] for x in sub])
        bA = np.array([[float(x["bAx"]), float(x["bAy"]), float(x["bAz"])] for x in sub])
        bB = np.array([[float(x["bBx"]), float(x["bBy"]), float(x["bBz"])] for x in sub])
        near = np.array([x["near"] == "1" for x in sub])
        out[pl] = (n, pts, bA, bB, near)
    return out


def drives():
    d = {}
    for x in rows("drives.csv"):
        d[x["case"]] = (complex(float(x["IA_re"]), float(x["IA_im"])), complex(float(x["IB_re"]), float(x["IB_im"])))
    return d


def ellipse(ia, ib, bA, bB):
    c = ia * bA + ib * bB
    mag2 = np.sum(np.abs(c) ** 2, axis=1)
    cc = np.abs(np.sum(c * c, axis=1))
    a2 = 0.5 * (mag2 + cc)
    b2 = np.maximum(0, 0.5 * (mag2 - cc))
    # direction of the major axis: Re(c·e^{-jψ}) with ψ = ½·arg(c·c)
    psi = 0.5 * np.angle(np.sum(c * c, axis=1))
    major = np.real(c * np.exp(-1j * psi)[:, None])
    return np.sqrt(a2), np.sqrt(mag2 / 2), np.sqrt(b2 / np.maximum(a2, 1e-300)), major


CASES = ["A only", "A + B in phase, toroidal sense", "A + B quadrature (sin/cos), toroidal sense",
         "A + B quadrature (sin/cos), dipole sense"]
SHORT = ["A only", "in phase, toroidal", "sin/cos, toroidal", "sin/cos, dipole"]


def fields():
    P = load_planes()
    D = drives()
    mu0 = 4e-7 * np.pi
    for pl, (u, v) in (("xy", (0, 1)), ("xz", (0, 2))):
        n, pts, bA, bB, near = P[pl]
        ext = [pts[:, u].min() * 100, pts[:, u].max() * 100, pts[:, v].min() * 100, pts[:, v].max() * 100]
        fig, axs = plt.subplots(4, 4, figsize=(9, 9.3))
        for i, case in enumerate(CASES):
            ia, ib = D[case]
            pk, rms, ell, maj = ellipse(ia, ib, bA, bB)
            for arr in (pk, rms, ell):
                arr[near] = np.nan
            en = rms ** 2 / (2 * mu0)
            vmax = np.nanpercentile(pk, 99) * 1e3
            imgs = [(pk * 1e3, "|B| peak [mT]", 0, vmax), (rms * 1e3, "|B| rms [mT]", 0, vmax / np.sqrt(2)),
                    (en, "⟨B²⟩/2μ0 [J/m³]", 0, np.nanpercentile(en, 99)), (ell, "ellipticity (0 swings, 1 turns)", 0, 1)]
            for j, (img, lab, lo, hi) in enumerate(imgs):
                ax = axs[i, j]
                im = ax.imshow(img.reshape(n, n).T, origin="lower", extent=ext, cmap=SEQ, vmin=lo, vmax=hi, interpolation="nearest")
                if j == 0:  # direction of the major axis (an axis: no arrowheads)
                    s = slice(None, None, max(1, n // 14))
                    U = maj[:, u].reshape(n, n).T
                    V = maj[:, v].reshape(n, n).T
                    L = np.hypot(U, V) + 1e-30
                    X = pts[:, u].reshape(n, n).T * 100
                    Y = pts[:, v].reshape(n, n).T * 100
                    ax.quiver(X[s, s], Y[s, s], (U / L)[s, s], (V / L)[s, s], color=INK, headwidth=0, headlength=0,
                              headaxislength=0, pivot="middle", scale=22, width=0.004)
                cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
                cb.ax.tick_params(labelsize=6)
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(lab)
                if j == 0:
                    ax.set_ylabel(SHORT[i])
        fig.suptitle(f"Field of the Rodin pair in the {'mid-plane z = 0' if pl == 'xy' else 'vertical plane y = 0'} (±7 cm; ±3 mm round the wires blank); "
                     "1 kHz, 1.25 A in A alone", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, f"fields_{pl}.png"))
        plt.close(fig)


def energy():
    fig, axs = plt.subplots(1, 3, figsize=(9, 2.9))
    r = rows("wm_lowf.csv")
    ax = axs[0]
    names = ["A only", "A + B in phase, toroidal sense", "A + B quadrature (sin/cos), toroidal sense", "A + B quadrature (sin/cos), dipole sense"]
    for i, nm in enumerate(names):
        t = np.array([float(x["t"]) for x in r if x["case"] == nm])
        w = np.array([float(x["Wm"]) for x in r if x["case"] == nm])
        m = t > t.max() - 2e-3
        ax.plot((t[m] - t[m][0]) * 1e3, w[m] * 1e6, color=SER[i], lw=1.4, label=SHORT[i])
    ax.set_xlabel("t [ms] (last two cycles at 1 kHz)")
    ax.set_ylabel("magnetic energy [µJ]")
    ax.set_title("Stored magnetic energy")
    ax.legend(fontsize=6, frameon=False)
    # B(t) at the tube point for the toroidal quadrature case
    D = drives()
    pts = rows("points.csv")
    ax = axs[1]
    wt = np.linspace(0, 2 * np.pi, 200)
    ia, ib = D["A + B quadrature (sin/cos), toroidal sense"]
    p = pts[1]
    bA = np.array([float(p["bAx"]), float(p["bAy"]), float(p["bAz"])])
    bB = np.array([float(p["bBx"]), float(p["bBy"]), float(p["bBz"])])
    b = np.real(np.outer(np.exp(1j * wt), ia * bA + ib * bB))
    ax.plot(wt / (2 * np.pi), b[:, 1] * 1e3, color=SER[0], lw=1.4, label="B_φ (round the axis)")
    ax.plot(wt / (2 * np.pi), b[:, 2] * 1e3, color=SER[1], lw=1.4, label="B_z")
    ax.plot(wt / (2 * np.pi), b[:, 0] * 1e3, color=SER[2], lw=1.4, label="B_ρ")
    ax.set_xlabel("t / T")
    ax.set_ylabel("B [mT]")
    ax.set_title("sin/cos, toroidal: B(t) in the tube at pole 1")
    ax.legend(fontsize=6, frameon=False)
    ax = axs[2]
    r = rows("ringup.csv")
    for i, nm in enumerate(["A only, tuned", "in phase, toroidal, tuned to L + M", "quadrature, toroidal, tuned to L + M",
                            "quadrature, toroidal sense, tuned", "A only, not tuned"]):
        t = np.array([float(x["t"]) for x in r if x["case"] == nm])
        w = np.array([float(x["Wm"]) for x in r if x["case"] == nm])
        if len(t):
            ax.plot(t * 1e3, w * 1e6, color=SER[i], lw=0.8, label=nm)
    ax.set_xlabel("t [ms] at 100 kHz")
    ax.set_ylabel("magnetic energy, cycle maximum [µJ]")
    ax.set_title("Build-up at resonance, 100 kHz")
    ax.legend(fontsize=5.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "energy.png"))
    plt.close(fig)


def readmap(name):
    with open(os.path.join(HERE, name), newline="", encoding="utf-8") as f:
        r = list(csv.reader(f))
    th = np.array([float(x) for x in r[0][1:]])
    f = np.array([float(x[0]) for x in r[1:]])
    z = np.array([[float(v) for v in x[1:]] for x in r[1:]])
    return f, th, z


def maps():
    fig, axs = plt.subplots(2, 4, figsize=(9.5, 5.4))
    k = 0
    for shape in ("tube", "top"):
        for sense in ("toroidal", "dipole"):
            for row, kind, lab in ((0, "emf", "open-circuit EMF [V]"), (1, "power", "load power [W]")):
                f, th, z = readmap(f"map_{kind}_{shape}_{sense}.csv")
                ax = axs[row, k]
                im = ax.pcolormesh(th, f, np.log10(np.maximum(z, 1e-12)), cmap=SEQ, shading="nearest")
                ax.set_yscale("log")
                if k == 0:
                    ax.set_ylabel("drive frequency [Hz]")
                ax.set_xlabel("angle round the toroid [°]")
                ax.set_xticks([0, 90, 180, 270, 360])
                ax.set_title(f"{shape} coil, {sense}: log10 {lab}", fontsize=7)
                cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
                cb.ax.tick_params(labelsize=6)
            k += 1
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "map.png"))
    plt.close(fig)


def modes():
    r = rows("modes.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3), sharey=False)
    tops = ["winding A open", "series, toroidal sense", "series, dipole sense", "plain ring", "crossed ring"]
    for ax, kind in zip(axs, ["Rodin (5, 12)", "(2, 1) torus knot"]):
        for j, top in enumerate(tops):
            sub = [x for x in r if x["kind"] == kind and x["topology"] == top]
            for x in sub:
                fr = float(x["f"]) / 1e6
                par = float(x["parity"])
                col = INK2
                if not np.isnan(par):
                    col = SER[7] if par < -0.5 else (SER[0] if par > 0.5 else SER[3])
                ax.plot([fr, fr], [j - 0.35, j + 0.35], color=col, lw=2)
        ax.set_yticks(range(len(tops)))
        ax.set_yticklabels(tops)
        ax.set_xscale("log")
        ax.set_xlabel("frequency [MHz]")
        ax.set_title(f"{kind}: first six modes")
    axs[1].plot([], [], color=SER[7], lw=2, label="inverted after one trip (720°)")
    axs[1].plot([], [], color=SER[0], lw=2, label="same way up after one trip")
    axs[1].plot([], [], color=SER[3], lw=2, label="mixed")
    axs[1].plot([], [], color=INK2, lw=2, label="open line (no ring)")
    axs[1].legend(fontsize=6, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "modes.png"))
    plt.close(fig)


def bank():
    if not have("bank.csv"):
        return
    fig, axs = plt.subplots(1, 3, figsize=(9, 2.8))
    for i, (tag, lab) in enumerate((("", "free search"), ("_dist", "distributed groups"))):
        if not have(f"bank{tag}.csv"):
            continue
        r = rows(f"bank{tag}.csv")
        t = np.array([float(x["t"]) for x in r])
        v = np.array([float(x["vbank"]) for x in r])
        e = np.array([float(x["ebank"]) for x in r])
        axs[0].plot(t / t[-1], v, color=SER[i], lw=1.4, label=f"{lab} ({t[-1]*1e3:.3g} ms)")
        axs[1].plot(t / t[-1], e * 1e6, color=SER[i], lw=1.4, label=lab)
        h = rows(f"ga_history{tag}.csv")
        axs[2].plot([int(x["generation"]) for x in h], [10 ** (6 - float(x["best"])) for x in h], color=SER[i], lw=1.4, label=lab)
    axs[0].set_xlabel("t / run length")
    axs[0].set_ylabel("bank voltage [V]")
    axs[0].set_title("Bank voltage (kronspice transient)")
    axs[0].legend(fontsize=6, frameon=False)
    axs[1].set_xlabel("t / run length")
    axs[1].set_ylabel("bank energy [µJ]")
    axs[1].set_title("Bank energy")
    axs[2].set_yscale("log")
    axs[2].set_xlabel("generation")
    axs[2].set_ylabel("P_load × bridge factor [W]")
    axs[2].set_title("Genetic searches, best so far")
    axs[2].legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "bank.png"))
    plt.close(fig)


def toroidal():
    if not have("tor_turn_emf.csv"):
        return
    fig = plt.figure(figsize=(9.5, 6.2))
    # (a) 3-D: the pair and the 12 sections, every turn through the hole
    ax = fig.add_subplot(2, 3, 1, projection="3d")
    P = {}
    for x in rows("paths.csv"):
        P.setdefault(x["winding"], []).append((float(x["x"]), float(x["y"]), float(x["z"])))
    P = {k: np.array(v) for k, v in P.items()}
    sec = rows("tor_sections.csv") if have("tor_sections.csv") else []
    draw_tor3d(ax, P, sec, every=4)
    ax.set_title("(a) 12 toroidal sections (green) wound\nround the tube over the pair (every 4th turn drawn)", fontsize=7)
    # (b) cross-section through one section: tube, bundles shell, turn
    ax = fig.add_subplot(2, 3, 2)
    for x in sec[:1]:
        a = float(x["radius"])
        u = np.linspace(0, 2 * np.pi, 200)
        for rr, col, lw, lab in ((A0, GRID, 1.5, "former"), (a, SER[2], 1.5, "output turn (round the tube)")):
            ax.plot((R0 + rr * np.cos(u)) * 100, rr * np.sin(u) * 100, color=col, lw=lw, label=lab)
            ax.plot((-R0 + rr * np.cos(u)) * 100, rr * np.sin(u) * 100, color=col, lw=lw)
        ax.axvline(0, color=INK2, lw=0.5, ls=":")
    for k, col in (("A", SER[0]), ("B", SER[1])):
        a3 = np.array(P[k])
        rho = np.hypot(a3[:, 0], a3[:, 1])
        sel = np.abs(np.arctan2(a3[:, 1], a3[:, 0])) < np.radians(4)
        ax.scatter(rho[sel] * 100, a3[sel, 2] * 100, s=4, color=col, label=f"winding {k} near 0°")
    ax.set_aspect("equal")
    ax.set_xlabel("ρ [cm] (meridian plane)")
    ax.set_ylabel("z [cm]")
    ax.set_ylim(-4.2, 2.8)
    ax.set_title("(b) cross-section: each turn goes through the hole", fontsize=7)
    ax.legend(fontsize=5.5, frameon=False, loc="lower center", ncol=2)
    # (c) EMF of one turn vs angle, 1 kHz
    ax = fig.add_subplot(2, 3, 3)
    r = rows("tor_turn_emf.csv")
    phi = np.array([float(x["phi_deg"]) for x in r])
    for key, col, lab in (("tor_sincos_amp", SER[0], "toroidal sense, sin/cos"), ("tor_inphase_amp", SER[0], "toroidal sense, in phase"),
                          ("dip_sincos_amp", SER[1], "dipole sense, sin/cos"), ("dip_inphase_amp", SER[1], "dipole sense, in phase"),
                          ("tor_A_amp", SER[3], "A alone"), ("tor_B_amp", SER[2], "B alone")):
        ls = "--" if "inphase" in key else "-"
        ax.plot(phi, np.array([float(x[key]) for x in r]) * 1e3, color=col, ls=ls, lw=1.1, label=lab)
    for k in range(13):
        ax.axvline(30 * k, color=GRID, lw=0.5)
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 60))
    ax.set_xlabel("azimuth of the turn [°] (poles at 0°, 30° …)")
    ax.set_ylabel("|EMF| per turn [mV]")
    ax.set_ylim(0, 23)
    ax.set_title("(c) one turn round the tube, 1 kHz", fontsize=7)
    ax.legend(fontsize=5.2, frameon=False, ncol=2, loc="upper center")
    # (d) map
    ax = fig.add_subplot(2, 3, 4)
    f, thm, z = readmap("tor_map_toroidal.csv")
    im = ax.pcolormesh(thm, f, np.log10(np.maximum(z * 1e3, 1e-6)), cmap=SEQ, shading="nearest")
    ax.set_yscale("log")
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel("azimuth of the turn [°]")
    ax.set_ylabel("drive frequency [Hz]")
    ax.set_title("(d) log10 |EMF| per turn [mV], toroidal sense, sin/cos", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    # (e) bank, (f) GA
    axb = fig.add_subplot(2, 3, 5)
    axg = fig.add_subplot(2, 3, 6)
    for i, (tag, lab) in enumerate((("", "sin/cos drive"), ("_phase", "drive phase evolved"))):
        if have(f"tor_bank{tag}.csv"):
            b = rows(f"tor_bank{tag}.csv")
            t = np.array([float(x["t"]) for x in b])
            axb.plot(t * 1e3, [float(x["vbank"]) for x in b], color=SER[i], lw=1.3, label=lab)
        if have(f"tor_ga_history{tag}.csv"):
            h = rows(f"tor_ga_history{tag}.csv")
            axg.plot([int(x["generation"]) for x in h], [10 ** (6 - float(x["best"])) for x in h], color=SER[i], lw=1.3, label=lab)
    axb.set_xlabel("t [ms]")
    axb.set_ylabel("bank voltage [V]")
    axb.set_title("(e) bank charging (kronspice transient)", fontsize=7)
    axb.legend(fontsize=6, frameon=False)
    axg.set_yscale("log")
    axg.set_xlabel("generation")
    axg.set_ylabel("P_load × bridge factor [W]")
    axg.set_title("(f) genetic searches, best so far", fontsize=7)
    axg.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "toroidal.png"))
    plt.close(fig)


if __name__ == "__main__":
    geometry()
    fields()
    energy()
    if have("map_emf_tube_toroidal.csv"):
        maps()
    modes()
    bank()
    toroidal()
    print("figures written to", HERE)
