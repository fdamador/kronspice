# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes

d = json.load(open("skyrme.json"))
LU = d["profiles"]["length_fm"]
INK, INK2, GRID, BG = "#1d1d1b", "#5c5b57", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREEN, GREY, RED = "#2a78d6", "#d0622a", "#2f8f5b", "#9a988f", "#b8322a"
plt.rcParams.update({"font.size": 13, "figure.facecolor": BG, "savefig.facecolor": BG, "axes.facecolor": BG,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID})

r1, F1 = np.array(d["profiles"]["r1"]), np.array(d["profiles"]["F1"])
r4, F4 = np.array(d["profiles"]["r4"]), np.array(d["profiles"]["F4"])


def bdens(r, F):
    dF = np.gradient(F, r)
    return -np.sin(F) ** 2 * dF / (2 * np.pi ** 2 * r ** 2)


# 1. profiles and baryon density
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for r, F, c, lab in ((r1, F1, BLUE, "nucleon knot (B = 1)"), (r4, F4, ORANGE, "helium-4 knot (B = 4)")):
    ax[0].plot(r * LU, F / np.pi, color=c, lw=2.5, label=lab)
    b = bdens(r, F)
    B = 4 if F is F4 else 1
    ax[1].plot(r * LU, 4 * np.pi * r ** 2 * b / B / LU, color=c, lw=2.5, label=lab)
ax[0].set_xlim(0, 3); ax[0].set_xlabel("radius (fm)"); ax[0].set_ylabel("twist F(r) / π")
ax[0].set_title("How far the field has turned: from π at the centre to 0 outside", fontsize=12, loc="left")
ax[1].set_xlim(0, 3); ax[1].set_xlabel("radius (fm)"); ax[1].set_ylabel("knot density per unit winding, 4πr²b(r)/B (1/fm)")
ax[1].set_title("Where the knot's winding sits", fontsize=12, loc="left")
for a in ax:
    a.legend(frameon=False)
fig.suptitle("Skyrme knots: the nucleon and helium-4 (e = 5.44, Fπ = 129 MeV, fitted to N and Δ masses)", fontsize=15)
fig.tight_layout(); fig.savefig("fig1_profiles.png", dpi=120); plt.close(fig)

# 2. proton vs neutron charge density
dF = np.gradient(F1, r1); s2 = np.sin(F1) ** 2
rhoS = bdens(r1, F1)                                     # normalized to 1
rhoV = s2 * (1 + 4 * (dF ** 2 + s2 / r1 ** 2))
rhoV = rhoV / np.trapezoid(4 * np.pi * r1 ** 2 * rhoV, r1)
p = 0.5 * rhoS + 0.5 * rhoV
n = 0.5 * rhoS - 0.5 * rhoV
fig, ax = plt.subplots(figsize=(10, 5.6))
ax.plot(r1 * LU, 4 * np.pi * r1 ** 2 * p / LU, color=RED, lw=2.6, label="proton (charge +1)")
ax.plot(r1 * LU, 4 * np.pi * r1 ** 2 * n / LU, color=BLUE, lw=2.6, label="neutron (charge 0)")
ax.axhline(0, color=INK2, lw=0.8)
ax.fill_between(r1 * LU, 0, 4 * np.pi * r1 ** 2 * n / LU, where=n > 0, color=BLUE, alpha=0.12)
ax.fill_between(r1 * LU, 0, 4 * np.pi * r1 ** 2 * n / LU, where=n < 0, color=BLUE, alpha=0.25)
ax.set_xlim(0, 2.5); ax.set_xlabel("radius (fm)"); ax.set_ylabel("charge per fm of radius, 4πr²ρ(r) (e/fm)")
ax.annotate("neutron: positive core", xy=(0.25, 0.25), xytext=(0.55, 0.55), color=BLUE, fontsize=12, arrowprops=dict(arrowstyle="->", color=BLUE))
ax.annotate("… and a negative skin", xy=(1.1, -0.11), xytext=(1.45, 0.3), color=BLUE, fontsize=12, arrowprops=dict(arrowstyle="->", color=BLUE))
ax.set_title("One knot, two spin states: where the charge sits in the proton and the neutron", fontsize=13, loc="left")
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig("fig2_charge.png", dpi=120); plt.close(fig)

# 3. comparison of models
models = ["one light ring", "quark light rings", "Skyrme knot", "measured"]
cols = [GREY, GREEN, ORANGE, INK]
mup = [1.000, 2.793, d["B1"]["mu_p"], 2.793]
mun = [None, -1.862, d["B1"]["mu_n"], -1.913]
rad = [0.210, None, d["B1"]["radius_fm"], 0.72]
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
for k, (vals, title, unit) in enumerate([(mup, "proton magnetic moment", "μN"), (mun, "neutron magnetic moment", "μN"), (rad, "size (isoscalar / charge radius)", "fm")]):
    for i, v in enumerate(vals):
        if v is None:
            ax[k].text(i, 0.02 if k != 1 else -0.05, "n/a" if k != 1 else "", ha="center", va="bottom", color=INK2, fontsize=10)
            continue
        ax[k].bar(i, v, color=cols[i], width=0.65)
        ax[k].text(i, v, f"{v:+.2f}" if unit == "μN" else f"{v:.2f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=11, color=INK)
    ax[k].set_xlim(-0.5, 3.5); ax[k].set_xticks(range(4)); ax[k].set_xticklabels(models, rotation=20, ha="right", fontsize=10)
    ax[k].set_title(f"{title} ({unit})", fontsize=12, loc="left"); ax[k].axhline(0, color=INK2, lw=0.8)
ax[1].text(0, -1.0, "p + e + ν̄\nrings:\n−1835.5\n(off scale)", ha="center", va="center", fontsize=10, color=INK2)
fig.suptitle("The nucleon in each picture we have built, against measurement", fontsize=15)
fig.tight_layout(); fig.savefig("fig3_compare.png", dpi=120); plt.close(fig)

# 4. helium binding
M1, M4 = d["B4"]["M1"], d["B4"]["M4"]
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].bar([0, 1], [M1, M4 / 4], color=[BLUE, ORANGE], width=0.6)
for i, v in enumerate([M1, M4 / 4]):
    ax[0].text(i, v, f"{v:.0f} MeV", ha="center", va="bottom")
ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["one knot (B = 1)", "helium-4 knot, per winding"])
ax[0].set_ylim(700, 900); ax[0].set_ylabel("classical knot mass (MeV)")
ax[0].set_title("The cube knot is lighter per winding: helium holds together", fontsize=12, loc="left")
ax[1].bar([0, 1], [100 * d["B4"]["frac"], 0.75], color=[ORANGE, INK], width=0.6)
ax[1].text(0, 100 * d["B4"]["frac"], f"{d['B4']['bind']:.0f} MeV\n{100 * d['B4']['frac']:.1f} %", ha="center", va="bottom")
ax[1].text(1, 0.75, "28.3 MeV\n0.75 %", ha="center", va="bottom")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["Skyrme knot", "measured"]); ax[1].set_ylim(0, 10)
ax[1].set_ylabel("binding energy (% of mass)")
ax[1].set_title("…but about ten times too tightly", fontsize=12, loc="left")
fig.suptitle("Helium-4 as one Skyrme knot (rational-map approximation)", fontsize=15)
fig.tight_layout(); fig.savefig("fig4_helium.png", dpi=120); plt.close(fig)

# 5. 3-D knots: baryon density isosurfaces
def density3d(r, F, B, n=130, L=2.6):
    xs = np.linspace(-L, L, n)
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    R = np.sqrt(X ** 2 + Y ** 2 + Z ** 2) + 1e-9
    Fr = np.interp(R, r, F)
    dFr = np.interp(R, r, np.gradient(F, r))
    radial = -np.sin(Fr) ** 2 * dFr / (2 * np.pi ** 2 * R ** 2)
    if B == 1:
        ang = np.ones_like(R)
    else:
        th = np.arccos(np.clip(Z / R, -1, 1)); ph = np.arctan2(Y, X)
        z = np.tan(th / 2) * np.exp(1j * ph)
        a = 2 * np.sqrt(3) * 1j
        num, den = z ** 4 + a * z ** 2 + 1, z ** 4 - a * z ** 2 + 1
        Rm = num / den
        dR = ((4 * z ** 3 + 2 * a * z) * den - num * (4 * z ** 3 - 2 * a * z)) / den ** 2
        ang = ((1 + np.abs(z) ** 2) * np.abs(dR) / (1 + np.abs(Rm) ** 2)) ** 2
    return xs, radial * ang


fig = plt.figure(figsize=(14, 7))
for k, (r, F, B, c, title) in enumerate(((r1, F1, 1, BLUE, "nucleon (proton or neutron): B = 1, a ball"),
                                         (r4, F4, 4, ORANGE, "helium-4: B = 4, a cube"))):
    xs, b = density3d(r, F, B, L=3.0 if B == 1 else 6.0)
    h = xs[1] - xs[0]
    tot = b.sum() * h ** 3
    lev = (0.35 if B == 1 else 0.5) * b.max()
    verts, faces, _, _ = marching_cubes(b, lev, spacing=(h, h, h))
    verts = (verts + xs[0]) * LU
    ax = fig.add_subplot(1, 2, k + 1, projection="3d")
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], color=c, alpha=0.9, lw=0, shade=True)
    lim = 4.6 * LU if B == 4 else 1.6 * LU
    for s in (ax.set_xlim, ax.set_ylim, ax.set_zlim):
        s(-lim, lim)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=22, azim=35)
    ax.set_xlabel("x (fm)", fontsize=10); ax.set_ylabel("y (fm)", fontsize=10); ax.set_zlabel("z (fm)", fontsize=10)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor("#f1f0ec"); a.pane.set_edgecolor(GRID)
    ax.set_title(f"{title}\n(surface of high winding density; total winding {tot:.2f})", fontsize=12)
fig.suptitle("Skyrme knots in 3-D: where the field's winding is concentrated", fontsize=15)
fig.tight_layout(); fig.savefig("fig5_knots3d.png", dpi=120); plt.close(fig)
print("ok")
