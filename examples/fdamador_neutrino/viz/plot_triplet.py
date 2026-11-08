# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""The three spin states of the ring (H41, Part 8): m = +1, 0, −1.

    go run ./examples/fdamador_neutrino -only 8 -field -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/plot_triplet.py

Writes viz/structure_m+1.png, viz/structure_m0.png, viz/structure_m-1.png in
the layout of viz/structure.png. Every field drawn here is the field
solver's (em/fdtd.go), not a sketch: the filtered phasors of Part 8d on the
finest grid run (r/3 by default, dx = 0.833 mm), the E-ring (TE01δ) for
m = 0 and the turning patterns (cos φ ± i sin φ)/√2 of the lowest odd
lattice mode for m = ±1, sampled by trilinear interpolation (FieldAt) at
t = T/8 (θ = 45°) for panels (a) and (b):

  (a) 3-D view of the bench ring (a = 10 mm, r = 2.5 mm, εr = 100): E on the
      ring's centre circle (orange) and B beside the tube and on the axis
      (blue), arrows scaled per colour; the turning sense (from the
      computed rotation of the dipole, which agrees with the sign of J_z)
      as a curved arrow above the ring;
  (b) the z = 0 plane: E in the plane (orange arrows) and B_z (blue
      contours, solid > 0, dashed < 0) at the same instant;
  (c) two periods at the state's own frequency, from the phasors: electric
      and magnetic energy / Ū, J_z·ω/Ū (Minkowski), and the dipole that
      the state carries (axial magnetic for m = 0, the in-plane electric
      dipole for m = ±1, whose magnetic dipole is zero to round-off).
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
A_MM, R_MM = 10.0, 2.5
E_COL, B_COL = "#d9480f", "#1c64b4"
TAGS = {"p1": ("structure_m+1.png", "+1"), "0": ("structure_m0.png", "0"), "m1": ("structure_m-1.png", "−1")}


def read(name):
    with open(os.path.join(HERE, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def field_meta():
    with open(os.path.join(HERE, "p8_field_meta.csv")) as f:
        rows = list(csv.DictReader(f))
    cells = max(float(r["cells"]) for r in rows)
    return {r["state"]: r for r in rows if float(r["cells"]) == cells}, cells


def turning_sense(cy, tag):
    """+1 counter-clockwise seen from +z, −1 clockwise, 0 none: the rotation
    of the dipole the state carries (p in the plane for ±1)."""
    if tag == "0":
        return 0
    px, py = cy["px"], cy["py"]
    cross = px * np.roll(py, -1) - py * np.roll(px, -1)
    return 1 if cross.mean() > 0 else -1


def torus(ax):
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, 2 * np.pi, 30))
    X = (A_MM + R_MM * np.cos(v)) * np.cos(u)
    Y = (A_MM + R_MM * np.cos(v)) * np.sin(u)
    Z = R_MM * np.sin(v)
    ax.plot_surface(X, Y, Z, color="#f2c6a8", alpha=0.3, linewidth=0)


def panel_a(ax, pts, sense, label):
    torus(ax)
    x, y, z = pts["x"] * 1e3, pts["y"] * 1e3, pts["z"] * 1e3
    E = np.vstack([pts["ex"], pts["ey"], pts["ez"]]).T
    B = np.vstack([pts["bx"], pts["by"], pts["bz"]]).T
    ring = np.abs(np.hypot(x, y) - A_MM) < 1e-6
    ring &= np.abs(z) < 1e-6
    em = np.abs(E[ring]).max() if ring.any() else 1
    bm = np.linalg.norm(B, axis=1).max()
    for k in np.where(ring)[0]:
        e = E[k] / em * 5
        ax.quiver(x[k], y[k], z[k], e[0], e[1], e[2], color=E_COL, lw=1.8, arrow_length_ratio=0.35)
    for k in np.where(~ring)[0]:
        b = B[k] / bm * 5
        if np.linalg.norm(b) < 0.25:
            continue
        ax.quiver(x[k], y[k], z[k], b[0], b[1], b[2], color=B_COL, lw=1.3, arrow_length_ratio=0.35)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(A_MM * np.cos(t), A_MM * np.sin(t), 0, color=E_COL, lw=0.8, alpha=0.6)
    if sense != 0:
        th = np.linspace(0.15 * np.pi, 1.55 * np.pi, 80)
        if sense < 0:
            th = th[::-1]
        rr = A_MM + 6
        ax.plot(rr * np.cos(th), rr * np.sin(th), 9 * np.ones_like(th), color="#333", lw=2)
        d = np.array([-np.sin(th[-1]), np.cos(th[-1]), 0]) * sense
        ax.quiver(rr * np.cos(th[-1]), rr * np.sin(th[-1]), 9, d[0] * 4, d[1] * 4, 0, color="#333", lw=2, arrow_length_ratio=0.6)
        ax.text(0, 0, 19, "pattern turns " + ("counter-clockwise" if sense > 0 else "clockwise") + "\n(seen from +z)",
                ha="center", fontsize=8.5, color="#333")
    else:
        ax.text(0, 0, 19, "no turning: E and B swap\nenergy in place", ha="center", fontsize=8.5, color="#333")
    lim = 22
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-16, 22)
    ax.set_box_aspect((1, 1, 0.85))
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
    ax.view_init(elev=28, azim=-60)
    ax.set_title(f"(a) m = {label}: E on the ring (orange), B beside the\ntube and on the axis (blue), at t = T/8", fontsize=10)


def panel_b(ax, pl, label):
    x, y = pl["x"] * 1e3, pl["y"] * 1e3
    n = int(round(np.sqrt(len(x))))
    X = x.reshape(n, n)
    Y = y.reshape(n, n)
    Bz = pl["bz"].reshape(n, n)
    Ex, Ey = pl["ex"].reshape(n, n), pl["ey"].reshape(n, n)
    bm = np.abs(Bz).max()
    lev = np.linspace(-0.9, 0.9, 10) * bm
    lev = lev[np.abs(lev) > 0.05 * bm]
    cs = ax.contour(X, Y, Bz, levels=lev[lev > 0], colors=B_COL, linewidths=1.0)
    ax.contour(X, Y, Bz, levels=lev[lev < 0], colors=B_COL, linewidths=1.0, linestyles="--")
    del cs
    s = 3
    em = np.hypot(Ex, Ey).max()
    ax.quiver(X[::s, ::s], Y[::s, ::s], Ex[::s, ::s] / em, Ey[::s, ::s] / em, color=E_COL, scale=14, width=0.004)
    t = np.linspace(0, 2 * np.pi, 200)
    for r in (A_MM - R_MM, A_MM + R_MM):
        ax.plot(r * np.cos(t), r * np.sin(t), color="#b08060", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(x.min(), x.max()); ax.set_ylim(y.min(), y.max())
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title(f"(b) z = 0 plane at t = T/8: E in the plane (arrows),\nB_z (blue contours, dashed < 0); m = {label}", fontsize=10)


def panel_c(ax, cy, meta, tag, label):
    th = np.concatenate([cy["theta"], cy["theta"] + 2 * np.pi]) / (2 * np.pi)
    rep = lambda v: np.concatenate([v, v])  # noqa: E731
    u = (cy["ue"] + cy["ub"]).mean()
    ax.plot(th, rep(cy["ue"] / u), color=E_COL, lw=1.6, label="electric energy / Ū")
    ax.plot(th, rep(cy["ub"] / u), color=B_COL, lw=1.6, label="magnetic energy / Ū")
    ax.plot(th, rep(cy["jzM_w"] / u), color="#222", lw=1.4, label="J_z·ω/Ū (Minkowski)")
    if tag == "0":
        mz = cy["mz"] / np.abs(cy["mz"]).max()
        ax.plot(th, rep(mz), color=B_COL, lw=1.0, ls="--", label="magnetic moment m_z (axial)")
    else:
        pm = np.hypot(cy["px"], cy["py"]).max()
        ax.plot(th, rep(cy["px"] / pm), color=E_COL, lw=1.0, ls="--", label="electric dipole p_x")
        ax.plot(th, rep(cy["py"] / pm), color=E_COL, lw=1.0, ls=":", label="electric dipole p_y")
    f = float(meta["f"]) / 1e9
    ax.axhline(0, color="#999", lw=0.6)
    ax.set_xlabel("time [periods]")
    ax.set_ylim(-1.15, 1.9)
    ax.set_xlim(0, 2)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9, ncol=2)
    if tag == "0":
        tt = "(c) two periods at f = %.4g GHz: B passes through zero\ntwice a period; J_z = 0; m along the axis, mean 0" % f
    else:
        tt = ("(c) two periods at f = %.4g GHz: U_B stays Ū/2, J_z·ω/Ū = %+.2f\n(lattice m = 1 mode); no magnetic dipole, p turns"
              % (f, float(meta["jzM"])))
    ax.set_title(tt, fontsize=10)


def main():
    meta, cells = field_meta()
    for tag, (out, label) in TAGS.items():
        pts = read(f"p8_pts_{tag}.csv")
        pl = read(f"p8_plane_{tag}.csv")
        cy = read(f"p8_cycle_{tag}.csv")
        sense = turning_sense(cy, tag)
        fig = plt.figure(figsize=(15, 5.6))
        head = {"p1": "turning one way", "0": "not turning", "m1": "turning the other way"}[tag]
        fig.suptitle(f"Part 8 (H41): the ring state m = {label}, {head} — field solver, grid r/{cells:g}", fontsize=13)
        ax = fig.add_subplot(1, 3, 1, projection="3d")
        panel_a(ax, pts, sense, label)
        ax = fig.add_subplot(1, 3, 2)
        panel_b(ax, pl, label)
        ax = fig.add_subplot(1, 3, 3)
        panel_c(ax, cy, meta[tag], tag, label)
        foot = ("Dielectric torus εr = 100 (a = 10 mm, r = 2.5 mm) in em/fdtd.go; fields are the filtered phasors of Part 8d. "
                "J_z labels the state; its sign is set by how the field is started, nothing in Maxwell's equations selects it.")
        fig.text(0.5, 0.01, foot, ha="center", fontsize=8.5, color="#444")
        fig.tight_layout(rect=(0, 0.04, 1, 0.94))
        path = os.path.join(HERE, out)
        fig.savefig(path, dpi=140)
        plt.close(fig)
        print("wrote", path)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        HERE = sys.argv[1]
    main()
