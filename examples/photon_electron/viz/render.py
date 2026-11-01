"""Picture for examples/photon_electron (run `go run . -viz viz` first).

    python3 viz/render.py

Writes
  viz/round_trip.png   the voltage between the two wires at the start point,
                       plain ring (top) and crossed ring (bottom), each trip marked
  viz/knot_light.png   the double-loop path, and how much of the light's energy
                       stays near it under linear Maxwell equations (field solver)
  viz/pumps.png        the three electron pumps of Part 4
  viz/magnetron.png    Part 5: the evolution of the magnetron-driven pump, and the
                       evolved design's electron bursts
Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 17, "axes.titlesize": 20, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})


def load(name):
    with open(os.path.join(HERE, f"{name}_trace.csv")) as f:
        rows = list(csv.reader(f))[1:]
    return [float(r[0]) for r in rows], [float(r[1]) for r in rows]


TRIP = 9.17  # µs, one trip round the ring
fig, axs = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
for ax, (name, title, colr) in zip(axs, [
        ("plain", "Plain ring: the pulse comes back the same way up every trip (360°)", "#1f5fbf"),
        ("crossed", "Crossed ring: upside down after one trip, upright after two (720°)", "#c0392b")]):
    t, v = load(name)
    ax.plot(t, v, color=colr, lw=2.5)
    ax.axhline(0, color="#555", lw=1)
    for n in range(1, 5):
        ax.axvline(n * TRIP, color="#888", ls=":", lw=1.5)
        ax.text(n * TRIP, 3.3, f"trip {n}", ha="center", fontsize=15)
    ax.set_ylim(-3.6, 3.8)
    ax.set_ylabel("volts between\nthe wires")
    ax.set_title(title, loc="left")
    ax.grid(alpha=0.3)
axs[1].set_xlabel("time (µs)")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "round_trip.png"), dpi=100)
print("wrote", os.path.join(HERE, "round_trip.png"))


# --- Part 3: light on the double loop ---------------------------------------
import numpy as np  # noqa: E402

rows = list(csv.reader(open(os.path.join(HERE, "knot_energy.csv"))))[1:]
T = [float(r[0]) for r in rows]
inside = [100 * float(r[1]) for r in rows]
gone = [100 * float(r[3]) for r in rows]
fig = plt.figure(figsize=(16, 7.5))
ax = fig.add_subplot(1, 2, 1, projection="3d")
u = np.linspace(0, 4 * np.pi, 800)
R, a = 1.0, 0.5
rr = R + a * np.cos(u / 2)
ax.plot(rr * np.cos(u), rr * np.sin(u), a * np.sin(u / 2), color="#c0392b", lw=4)
th = np.linspace(0, 2 * np.pi, 100)
ax.plot(R * np.cos(th), R * np.sin(th), 0 * th, color="#888", ls=":", lw=1.5)
ax.set_box_aspect((1, 1, 0.45))
ax.set_xticks([]), ax.set_yticks([]), ax.set_zticks([])
ax.view_init(elev=28, azim=30)
ax.set_title("The paper's path: twice round the ring\n(one wavelength of light)", fontsize=18)
ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(T, inside, "o-", color="#c0392b", lw=3, ms=9, label="energy still near the loop (within 2R)")
ax2.plot(T, gone, "s-", color="#1f5fbf", lw=3, ms=8, label="energy that has left the box (±6R)")
ax2.axhline(100, color="#27864a", ls="--", lw=2.5)
ax2.text(1.0, 104, "what the model needs: 100 % stays", color="#27864a", fontsize=16, ha="center")
ax2.set_ylim(-5, 115)
ax2.set_xlabel("time (in periods of the light)")
ax2.set_ylabel("percent of the starting energy")
ax2.set_title("Ordinary Maxwell equations: the light leaves", loc="left")
ax2.grid(alpha=0.3)
ax2.legend(loc="center right", fontsize=15)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "knot_light.png"), dpi=100)
print("wrote", os.path.join(HERE, "knot_light.png"))


# --- Part 4: the three pumps --------------------------------------------------
def rd(name):
    rows = list(csv.reader(open(os.path.join(HERE, name))))[1:]
    return [[float(x) for x in r] for r in rows]

fig, axs = plt.subplots(1, 3, figsize=(20, 7))
se = rd("pump_single.csv")
ax = axs[0]
ax.plot([r[0] / 1e6 for r in se], [r[1] for r in se], "o-", color="#1f5fbf", lw=3, ms=10)
ax.axhline(1, color="#27864a", ls="--", lw=2)
ax.set_xscale("log")
ax.set_ylim(0.6, 1.05)
ax.set_xlabel("drive frequency (MHz)")
ax.set_ylabel("electrons moved per cycle")
ax.set_title("1. Single-electron pump:\nexactly one per cycle until the drive is too fast", fontsize=17)
ax.grid(alpha=0.3, which="both")

mp = rd("pump_multipactor.csv")
ax = axs[1]
ax.semilogy([r[0] for r in mp], [max(r[1], 1e-3) for r in mp], color="#c0392b", lw=2.5, label="electrons in the gap")
ax.set_ylim(1e-3, 1e10)
ax.set_xlabel("time (ns)")
ax.set_ylabel("electrons in the gap", color="#c0392b")
ax2 = ax.twinx()
ax2.plot([r[0] for r in mp], [r[2] for r in mp], color="#1f5fbf", lw=2.5)
ax2.set_ylabel("gap voltage (V)", color="#1f5fbf")
ax2.set_ylim(0, 160)
ax.set_title("2. Multipactor pump: bursts —\nelectrons multiply, drain the resonator, restart", fontsize=17)
ax.grid(alpha=0.3)

gp = rd("pump_ground.csv")
ax = axs[2]
ax.plot([r[0] for r in gp], [r[1] for r in gp], color="#8e44ad", lw=2.5)
ax.axhline(-1199.5, color="#27864a", ls="--", lw=2)
ax.text(3, -1150, "ladder's limit −2·n·peak", color="#27864a", fontsize=15, ha="center")
ax.set_xlabel("time (ms)")
ax.set_ylabel("bank voltage (V)")
ax.set_title("3. Ground charge pump:\n10 V → 200 V by resonance → −1,195 V bank", fontsize=17)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "pumps.png"), dpi=100)
print("wrote", os.path.join(HERE, "pumps.png"))


# --- Part 5: magnetron-driven pump, evolved -----------------------------------
h = rd("magnetron_history.csv")
mb = rd("magnetron_best.csv")
fig, axs = plt.subplots(1, 2, figsize=(17, 7))
ax = axs[0]
ax.semilogy([r[0] for r in h], [r[1] for r in h], "o-", color="#1f5fbf", lw=3, ms=8)
ax.set_xlabel("generation")
ax.set_ylabel("mains watts per collected microamp (lower is better)")
ax.set_title("Genetic algorithm: 40 designs per generation", loc="left")
ax.grid(alpha=0.3, which="both")
ax = axs[1]
ax.semilogy([r[0] for r in mb], [max(r[1], 1e-3) for r in mb], color="#c0392b", lw=2)
ax.set_ylim(1e-3, 1e11)
ax.set_xlabel("time (ns)")
ax.set_ylabel("electrons in the gap", color="#c0392b")
ax2 = ax.twinx()
ax2.plot([r[0] for r in mb], [r[2] for r in mb], color="#1f5fbf", lw=2)
ax2.set_ylabel("gap voltage (V)", color="#1f5fbf")
ax2.set_ylim(0, max(r[2] for r in mb) * 1.15)
ax.set_title("Evolved design: bursts that pump about 0.85 mA", loc="left")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "magnetron.png"), dpi=100)
print("wrote", os.path.join(HERE, "magnetron.png"))
