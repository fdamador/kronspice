# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""What the photon electron looks like: light twice round a torus (H23).

    python3 examples/photon_electron/viz/plot_structure.py [OUT.png]

Writes viz/structure.png (default), the companion of
examples/fdamador_neutrino/viz/structure.png (the same field inverted, H28):
  (a) the paper's (2, 1) torus-knot path, twice round a ring of radius
      R = λ/4π = 193 fm while once round a tube of radius R/2, with one
      Compton wavelength of light on it: E (orange) and B (blue) across the
      path, E×B along it, drawn exactly as field.go launches them
      (E ∝ −sin(ks)·n2, B ∝ sin(ks)·n1/c, n2 = ẑ taken off the tangent,
      n1 = n2 × t);
  (b) the 720° property: E·n2 along the path against the angle round the
      ring; one wavelength fills two trips, so after one trip the field is
      inverted and only after two does it match itself. Markers: the
      pulse measured at the start point after each trip on the plain and
      crossed two-wire ring lines (README §2, kronspice transient);
  (c) at one point on the path E and c·B rise and fall together (a
      travelling wave, unlike the neutrino ring's quarter-period swap), at
      f = m_e c²/h = 1.24×10²⁰ Hz; markers: the energy still within 2R in
      the field-solver run (README §3), which is gone within half a period.
Amplitudes are normalised; geometry, frequency and the plotted data are the
example's own.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LAMBDA_FM = 2426.3102           # Compton wavelength h/(m_e c) in fm (README §1: 2.426×10⁻¹² m)
R = LAMBDA_FM / (4 * np.pi)      # 193 fm
A = R / 2
F0 = 1.2356e20                   # m_e c²/h, Hz (README §4: 1.24×10²⁰ Hz)
E_COL, B_COL, P_COL = "#d9480f", "#1c64b4", "#7b2d8b"

# README §2: pulse back at the start after each trip, fraction of launch
PLAIN = [0.90, 0.69, 0.64, 0.57]
CROSSED = [-0.91, 0.70, -0.59, 0.52]
# README §3: energy within 2R (%) at time in periods
KNOT_T = [0, 0.125, 0.25, 0.5, 1, 2]
KNOT_E = [98.9, 65.1, 5.8, 1.6e-6, 1e-4, 6e-7]


def knot(u):
    """(2, 1) torus knot as in field.go: φ = u, θ = u/2; point and unit tangent."""
    th = u / 2
    rr = R + A * np.cos(th)
    p = np.stack([rr * np.cos(u), rr * np.sin(u), A * np.sin(th)], -1)
    d = np.stack([-A / 2 * np.sin(th) * np.cos(u) - rr * np.sin(u),
                  -A / 2 * np.sin(th) * np.sin(u) + rr * np.cos(u),
                  A / 2 * np.cos(th) * np.ones_like(u)], -1)
    return p, d / np.linalg.norm(d, axis=-1, keepdims=True)


def frame(u):
    p, t = knot(u)
    n2 = np.stack([-t[..., 2] * t[..., 0], -t[..., 2] * t[..., 1], 1 - t[..., 2] ** 2], -1)
    n2 /= np.linalg.norm(n2, axis=-1, keepdims=True)
    n1 = np.cross(n2, t)
    return p, t, n1, n2


def arc(u):
    """Arc length along the knot from u = 0 and the total length L."""
    uu = np.linspace(0, 4 * np.pi, 20001)
    p, _ = knot(uu)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))])
    return np.interp(u, uu, s), s[-1]


def main(out):
    fig = plt.figure(figsize=(15, 6.2))
    fig.suptitle("The photon electron (H23, Williamson & van der Mark 1997): one wavelength of light "
                 "twice round a torus", fontsize=13)

    # (a) 3-D
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, 2 * np.pi, 30))
    ax.plot_surface((R + A * np.cos(vv)) * np.cos(uu), (R + A * np.cos(vv)) * np.sin(uu), A * np.sin(vv),
                    color="#d9d2e9", alpha=0.25, linewidth=0)
    u = np.linspace(0, 4 * np.pi, 800)
    p, _ = knot(u)
    ax.plot(p[:, 0], p[:, 1], p[:, 2], color=P_COL, lw=2.5)
    us = np.linspace(0, 4 * np.pi, 28, endpoint=False)
    ps, ts, n1, n2 = frame(us)
    s, L = arc(us)
    amp = np.sin(2 * np.pi * s / L)
    sc = 70.0
    e = (-amp)[:, None] * n2 * sc
    b = amp[:, None] * n1 * sc
    ax.quiver(ps[:, 0], ps[:, 1], ps[:, 2], e[:, 0], e[:, 1], e[:, 2], color=E_COL, lw=1.6, arrow_length_ratio=0.3)
    ax.quiver(ps[:, 0], ps[:, 1], ps[:, 2], b[:, 0], b[:, 1], b[:, 2], color=B_COL, lw=1.6, arrow_length_ratio=0.3)
    k = 3
    ax.quiver(ps[k, 0], ps[k, 1], ps[k, 2], *(ts[k] * 90), color=P_COL, lw=2.5, arrow_length_ratio=0.35)
    lim = 330
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-200, 200)
    ax.set_box_aspect((1, 1, 0.6))
    ax.set_xlabel("x [fm]"); ax.set_ylabel("y [fm]"); ax.set_zlabel("z [fm]")
    ax.view_init(elev=28, azim=-58)
    ax.set_title("(a) the (2, 1) path (purple): twice round the ring, once round\n"
                 f"the tube; E (orange), B (blue) across it; R = λ/4π = {R:.0f} fm", fontsize=10)

    # (b) 720°
    ax = fig.add_subplot(1, 3, 2)
    u = np.linspace(0, 8 * np.pi, 4000)
    s_, L = arc(np.mod(u, 4 * np.pi))
    s_ = s_ + L * np.floor(u / (4 * np.pi))
    for k0 in (0, 720):
        ax.axvspan(k0, k0 + 360, color="#f3f0f8")
    ax.plot(np.degrees(u), -np.sin(2 * np.pi * s_ / L), color=E_COL, lw=1.6, label="E·n2 along the path (t = 0)")
    ax.axhline(0, color="#999", lw=0.6)
    for k0, lab in ((180, "trip 1"), (540, "trip 2"), (900, "trip 3"), (1260, "trip 4")):
        ax.text(k0, 1.14, lab, ha="center", fontsize=8.5, color="#555")
    ax.plot([360, 720, 1080, 1440], PLAIN, "s", color="#555", ms=7, label="plain ring line, pulse after each trip")
    ax.plot([360, 720, 1080, 1440], CROSSED, "o", color=P_COL, ms=8, label="crossed ring line (kronspice, README §2)")
    ax.set_xlim(0, 1460)
    ax.set_ylim(-1.25, 1.3)
    ax.set_xticks([0, 360, 720, 1080, 1440])
    ax.set_xlabel("angle round the ring along the path [°]")
    ax.set_ylabel("field / returned pulse (normalised)")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, framealpha=0.9)
    ax.set_title("(b) spin ½ as the paper reads it: after 360° the field is inverted,\n"
                 "after 720° it matches itself (crossed line: −, +, −, +)", fontsize=10)

    # (c) time at one point, and the field solver's verdict
    ax = fig.add_subplot(1, 3, 3)
    T = 1 / F0
    tt = np.linspace(0, 2, 1000)
    ax.plot(tt, np.sin(2 * np.pi * tt), color=E_COL, lw=1.6, label="E at one point on the path")
    ax.plot(tt, np.sin(2 * np.pi * tt), color=B_COL, lw=1.2, ls="--", label="c·B at the same point (in phase)")
    ax.set_xlabel(f"time [periods of the light; one period = {T * 1e21:.2f}×10⁻²¹ s]")
    ax.set_ylabel("normalised field")
    ax.set_ylim(-1.3, 1.3)
    ax3 = ax.twinx()
    ax3.plot(KNOT_T, KNOT_E, "o-", color="#2b8a3e", lw=1.4, ms=6, label="energy still within 2R, field solver (README §3)")
    ax3.set_ylabel("energy within 2R [%]", color="#2b8a3e")
    ax3.set_ylim(-5, 105)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, framealpha=0.9)
    ax.set_xlim(0, 2)
    ax.grid(alpha=0.3)
    ax.set_title("(c) a travelling wave at f = m_e c²/h = 1.24×10²⁰ Hz;\n"
                 "under linear Maxwell it leaves within half a period", fontsize=10)

    fig.text(0.5, 0.01, "Path and field frame as launched in field.go (λ = Compton wavelength, R = λ/4π, tube R/2). "
             "Inverted (E as the ring, B round it) this becomes the fdamador neutrino, H28 "
             "(examples/fdamador_neutrino/viz/structure.png).", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 0.94))
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "structure.png"))
