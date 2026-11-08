# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""What the fdamador neutrino looks like: the E ring and the B round it (H28).

    python3 examples/fdamador_neutrino/viz/plot_structure.py [OUT.png]

Writes viz/structure.png (default):
  (a) 3-D view of the README's bench ring (a = 10 mm, r = 2.5 mm): E runs
      round the ring (toroidal), B loops through the hole and round the tube
      (poloidal), drawn in eight meridional planes;
  (b) the same B lines in the (ρ, z) half-plane, with E going into and out
      of the page through the tube's cross-section;
  (c) one ring-down at the README's f0 = 2.0606 GHz and Q = 40.52: the
      ring's E (displacement current round the ring) and its magnetic
      moment (B through the hole) a quarter period apart, their energy
      sloshing between electric and magnetic, and m(t) averaging to zero.

B is the field of a circular current loop of radius a (the ring's
displacement current I = C dV/dt treated as a thin ring): field lines are
the contours of ψ = ρ·A_φ, with A_φ from Smythe's elliptic-integral form,
K and E computed by the arithmetic-geometric mean. Amplitudes are
normalised; only the geometry and f0, Q are physical.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
A_MM, R_MM = 10.0, 2.5           # ring radius and tube radius (README bench ring)
F0, Q = 2.0606e9, 40.52          # README §1, εr = 100
E_COL, B_COL = "#d9480f", "#1c64b4"


def ellipke(m):
    """Complete elliptic integrals K(m), E(m) by the arithmetic-geometric mean."""
    a = np.ones_like(m)
    b = np.sqrt(1.0 - m)
    c = np.sqrt(m)
    s = 0.5 * c * c
    p = 1.0
    for _ in range(30):
        a, b, c = 0.5 * (a + b), np.sqrt(a * b), 0.5 * (a - b)
        p *= 2.0
        s += 0.5 * p * c * c
    K = np.pi / (2.0 * a)
    return K, K * (1.0 - s)


def psi(rho, z, a=A_MM):
    """ρ·A_φ of a unit current loop of radius a (μ0 = 1 units); B lines are its contours."""
    m = 4.0 * a * rho / ((a + rho) ** 2 + z ** 2)
    m = np.clip(m, 1e-12, 1 - 1e-12)
    K, E = ellipke(m)
    k = np.sqrt(m)
    return rho * np.sqrt(a / rho) / (np.pi * k) * ((1 - m / 2) * K - E)


def b_lines():
    rho = np.linspace(0.05, 32, 400)
    z = np.linspace(-24, 24, 400)
    RR, ZZ = np.meshgrid(rho, z)
    P = psi(RR, ZZ)
    levels = psi(np.array([A_MM + 1.2, A_MM + 3.0, A_MM + 6.0, A_MM + 11.0, A_MM + 18.0]), np.zeros(5))
    return rho, z, P, np.sort(levels)


def main(out):
    rho, z, P, levels = b_lines()
    fig = plt.figure(figsize=(15, 5.6))
    fig.suptitle("The fdamador neutrino (H28): the toroidal field inverted — E is the ring, B circles round it",
                 fontsize=13)

    # (a) 3-D
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, 2 * np.pi, 30))
    X = (A_MM + R_MM * np.cos(v)) * np.cos(u)
    Y = (A_MM + R_MM * np.cos(v)) * np.sin(u)
    Z = R_MM * np.sin(v)
    ax.plot_surface(X, Y, Z, color="#f2c6a8", alpha=0.35, linewidth=0)
    cs = plt.figure().add_subplot().contour(rho, z, P, levels=levels)
    plt.close()
    for phi in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        for seg_list in cs.allsegs:
            for seg in seg_list:
                r_, z_ = seg[:, 0], seg[:, 1]
                ax.plot(r_ * np.cos(phi), r_ * np.sin(phi), z_, color=B_COL, lw=0.9, alpha=0.75)
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(A_MM * np.cos(t), A_MM * np.sin(t), 0, color=E_COL, lw=3)
    for th in np.linspace(0, 2 * np.pi, 10, endpoint=False):
        x0, y0 = A_MM * np.cos(th), A_MM * np.sin(th)
        ax.quiver(x0, y0, 0, -np.sin(th) * 4, np.cos(th) * 4, 0, color=E_COL, lw=2, arrow_length_ratio=0.45)
    ax.quiver(0, 0, -14, 0, 0, 28, color=B_COL, lw=2, arrow_length_ratio=0.12)
    ax.set_box_aspect((1, 1, 0.75))
    lim = 22
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-16, 16)
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.set_zlabel("z [mm]")
    ax.view_init(elev=24, azim=-60)
    ax.set_title("(a) E round the ring (orange), B through the hole\nand round the tube (blue); a = 10 mm, r = 2.5 mm",
                 fontsize=10)

    # (b) meridional cross-section: B = curl A from ψ, traced as streamlines on the full plane
    ax = fig.add_subplot(1, 3, 2)
    xs = np.linspace(-28, 28, 281)
    zs = np.linspace(-20, 20, 201)
    XX, ZZ = np.meshgrid(xs, zs)
    h = 1e-3
    rr = np.maximum(np.abs(XX), 2e-3)
    Bz = (psi(rr + h, ZZ) - psi(rr - h, ZZ)) / (2 * h) / rr
    Brho = -(psi(rr, ZZ + h) - psi(rr, ZZ - h)) / (2 * h) / rr
    Bx = np.sign(XX) * Brho
    inside = np.hypot(np.abs(XX) - A_MM, ZZ) < R_MM * 0.35
    Bx[inside] = np.nan
    Bz[inside] = np.nan
    seeds = [(sx * (A_MM + d), 0.0) for sx in (1, -1) for d in (3.0, 4.2, 5.8, 8.0, 11.0, 15.0)]
    ax.streamplot(xs, zs, Bx, Bz, start_points=seeds, color=B_COL, linewidth=1.0,
                  density=10, arrowsize=0.9, integration_direction="both", broken_streamlines=False)
    for sx, mark in ((1, "x"), (-1, "o")):
        ax.add_patch(plt.Circle((sx * A_MM, 0), R_MM, color="#f2c6a8", ec=E_COL, lw=1.5, zorder=3))
        ax.plot(sx * A_MM, 0, marker=mark, color=E_COL, ms=11, mew=2.5,
                mfc="none" if mark == "o" else E_COL, zorder=4)
        if mark == "o":
            ax.plot(sx * A_MM, 0, ".", color=E_COL, ms=6, zorder=4)
    ax.annotate("", xy=(0, 9), xytext=(0, -9), arrowprops=dict(arrowstyle="->", color=B_COL, lw=2))
    ax.text(0.8, 9.5, "B through the hole", color=B_COL, fontsize=9)
    ax.text(A_MM, -15.5, "E into the page", color=E_COL, fontsize=8.5, ha="center")
    ax.text(-A_MM, -15.5, "E out of the page", color=E_COL, fontsize=8.5, ha="center")
    ax.axvline(0, color="#999", lw=0.6, ls=":")
    ax.set_aspect("equal")
    ax.set_xlim(-28, 28); ax.set_ylim(-20, 20)
    ax.set_xlabel("ρ [mm] (cut through the ring's axis)"); ax.set_ylabel("z [mm]")
    ax.set_title("(b) cross-section: B lines circle the tube and thread\nthe hole, as for a ring current (TE01δ pattern)",
                 fontsize=10)

    # (c) ring-down
    ax = fig.add_subplot(1, 3, 3)
    w0 = 2 * np.pi * F0
    tt = np.linspace(0, 12e-9, 4000)
    env = np.exp(-w0 * tt / (2 * Q))
    e = env * np.cos(w0 * tt)
    m = env * np.sin(w0 * tt)
    ax.plot(tt * 1e9, e, color=E_COL, lw=1.0, label="E round the ring (capacitor voltage)")
    ax.plot(tt * 1e9, m, color=B_COL, lw=1.0, label="magnetic moment m(t) (B through the hole)")
    ax.plot(tt * 1e9, env, color="#555", lw=0.8, ls="--", label=f"envelope, Q = {Q}")
    ax.plot(tt * 1e9, -env, color="#555", lw=0.8, ls="--")
    ax.set_xlabel("time [ns]"); ax.set_ylabel("normalised amplitude")
    ax.set_title(f"(c) it exists only as an oscillation (f0 = {F0 / 1e9:.4f} GHz):\n"
                 "E and B a quarter period apart; m(t) averages to 0", fontsize=10)
    ax.legend(fontsize=7.5, loc="lower right", framealpha=0.9)
    ax.set_xlim(0, 12)
    ax.grid(alpha=0.3)

    fig.text(0.5, 0.01, "Bench ring of the README (§1): dielectric ring εr = 100. B is the field of the ring's "
             "displacement current as a thin ring (field lines = contours of ρ·A_φ). In vacuum the same ring "
             "has Q ≈ 0.04 and is gone in under one period (§2).", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "structure.png"))
