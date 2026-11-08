# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""One cycle of the photon electron (H23), animated, beside its circuit mirror.

    python3 examples/photon_electron/viz/plot_cycle.py [OUT.gif]

Writes viz/cycle.gif: one period of the light (f = m_e c²/h = 1.24×10²⁰ Hz,
T = 8.09×10⁻²¹ s), slowed so the loop lasts 3 s (about 4×10²⁰×). Same path
and field frame as viz/structure.png (field.go):

  (a) 3-D: the (2, 1) path twice round the ring, the wave cos(ks − ωt)
      travelling along it, E (orange) and B (blue) across the path;
  (b) the field along the path over its two trips (0–720°) at this moment:
      the pattern moves on by one wavelength (both trips) per period, and
      at every moment the value after 360° is minus the value at the start;
  (c) the circuit mirror (viz/circuit.png): the parallel tank, its voltage
      (E through the hole) and inductor current (B round the ring) a
      quarter period apart, and the energy passing between C and L.

Panels (a, b) are the model's own travelling wave, in which E and B rise and
fall together; panel (c) is the lumped mirror of the fdamador neutrino, in
which they alternate. The two agree on frequency per mode, not on the
phase between E and B. The escape of the light (README §3) is left out.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_structure import A, B_COL, E_COL, F0, P_COL, R, arc, frame, knot  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
NF = 48


def main(out):
    us = np.linspace(0, 4 * np.pi, 32, endpoint=False)
    ps, ts, n1, n2 = frame(us)
    s_us, Lp = arc(us)
    upath = np.linspace(0, 4 * np.pi, 600)
    ppath, _ = knot(upath)
    uu, vv = np.meshgrid(np.linspace(0, 2 * np.pi, 48), np.linspace(0, 2 * np.pi, 20))
    torus = ((R + A * np.cos(vv)) * np.cos(uu), (R + A * np.cos(vv)) * np.sin(uu), A * np.sin(vv))
    ang = np.linspace(0, 4 * np.pi, 800)
    s_ang, _ = arc(ang)
    tt = np.linspace(0, 1, 400)

    fig = plt.figure(figsize=(15, 5.4))
    fig.suptitle("The photon electron (H23): one cycle of light twice round the torus, "
                 "beside its circuit mirror (parallel tank)", fontsize=12.5)
    ax3 = fig.add_subplot(1, 3, 1, projection="3d")
    axb = fig.add_subplot(1, 3, 2)
    axc = fig.add_subplot(1, 3, 3)

    def draw(k):
        ph = 2 * np.pi * k / NF
        amp = np.cos(2 * np.pi * s_us / Lp - ph)
        ax3.cla()
        ax3.plot_surface(*torus, color="#d9d2e9", alpha=0.22, linewidth=0)
        ax3.plot(ppath[:, 0], ppath[:, 1], ppath[:, 2], color=P_COL, lw=1.6, alpha=0.7)
        e = (-amp)[:, None] * n2 * 75
        b = amp[:, None] * n1 * 75
        ax3.quiver(ps[:, 0], ps[:, 1], ps[:, 2], e[:, 0], e[:, 1], e[:, 2], color=E_COL, lw=1.5, arrow_length_ratio=0.3)
        ax3.quiver(ps[:, 0], ps[:, 1], ps[:, 2], b[:, 0], b[:, 1], b[:, 2], color=B_COL, lw=1.5, arrow_length_ratio=0.3)
        j = int(np.argmax(amp))
        ax3.scatter(*ps[j], color="k", s=30)
        ax3.set_xlim(-330, 330); ax3.set_ylim(-330, 330); ax3.set_zlim(-200, 200)
        ax3.set_box_aspect((1, 1, 0.6))
        ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
        ax3.view_init(elev=28, azim=-58 + 90 * k / NF)
        ax3.set_title(f"(a) the wave travels along the (2, 1) path; E, B across it\n"
                      f"R = λ/4π = {R:.0f} fm (black dot: a crest)", fontsize=9.5)

        axb.cla()
        f = -np.cos(2 * np.pi * s_ang / Lp - ph)
        axb.axvspan(0, 360, color="#f3f0f8")
        axb.plot(np.degrees(ang), f, color=E_COL, lw=1.8)
        axb.axhline(0, color="#999", lw=0.6)
        axb.axvline(360, color="#555", lw=0.8, ls=":")
        i0 = 0
        i1 = np.argmin(np.abs(np.degrees(ang) - 360))
        axb.plot([0, 360], [f[i0], f[i1]], "o", color=P_COL, ms=8)
        axb.text(10, 1.12, "trip 1", fontsize=8.5, color="#555"); axb.text(370, 1.12, "trip 2", fontsize=8.5, color="#555")
        axb.set_xlim(0, 720); axb.set_ylim(-1.3, 1.3); axb.set_xticks([0, 180, 360, 540, 720])
        axb.set_xlabel("angle round the ring along the path [°]"); axb.set_ylabel("E across the path (normalised)")
        axb.set_title("(b) the field along both trips now: the value at 360°\n"
                      f"is always minus the value at 0° (purple); t = {k / NF:.2f} T", fontsize=9.5)

        axc.cla()
        v, il = np.cos(ph), np.sin(ph)
        axc.plot(tt, np.cos(2 * np.pi * tt), color=E_COL, lw=1.6, label="tank voltage V (E through the hole)")
        axc.plot(tt, np.sin(2 * np.pi * tt), color=B_COL, lw=1.6, label="inductor current I_L (B round the ring)")
        axc.axvline(k / NF, color="#333", lw=1)
        axc.plot(k / NF, v, "o", color=E_COL); axc.plot(k / NF, il, "o", color=B_COL)
        axc.bar([1.12, 1.26], [v * v, il * il], width=0.11, color=[E_COL, B_COL], alpha=0.8)
        axc.text(1.12, -0.13, "U_C", ha="center", fontsize=8.5); axc.text(1.26, -0.13, "U_L", ha="center", fontsize=8.5)
        axc.axvline(1.04, color="#bbb", lw=0.6)
        axc.set_xlim(0, 1.34); axc.set_ylim(-1.2, 1.25); axc.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        axc.set_xlabel("time [periods]"); axc.grid(alpha=0.3)
        axc.legend(fontsize=7.5, loc="lower left")
        axc.set_title("(c) the circuit mirror (parallel tank): energy C → L → C;\n"
                      "the model's own E and B (a, b) rise and fall together", fontsize=9.5)
        return []

    fig.text(0.5, 0.015, "f = m_e c²/h = 1.24×10²⁰ Hz, T = 8.09×10⁻²¹ s, slowed about 4×10²⁰× to 3 s. "
             "Escape left out: under linear Maxwell the light leaves within ¼–½ period (README §3). "
             "Mirror of examples/fdamador_neutrino/viz/cycle.gif.", ha="center", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    FuncAnimation(fig, draw, frames=NF).save(out, writer=PillowWriter(fps=16), dpi=80)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cycle.gif"))
