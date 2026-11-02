"""Pictures for examples/kaluza_kron (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/kaluza_kron.png: orbits from 5-D geodesics, charge-to-mass
against fifth-dimension momentum, and the Kron machine solved three ways.
Light background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

fig, axs = plt.subplots(1, 3, figsize=(24, 8))
orb = {}
for r in csv.DictReader(open(os.path.join(HERE, "orbits.csv"))):
    orb.setdefault(r["pw"], []).append((float(r["x"]), float(r["y"])))
ax = axs[0]
for (pw, pts), col in zip(orb.items(), ["#7f8c8d", "#e67e22", "#1f5fbf", "#16a085", "#c0392b"]):
    if float(pw) < 0.2:
        continue
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, lw=4, label="p_w = %s" % pw)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("5-D free fall seen in 4-D: circles in a magnetic field", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13)

rows = list(csv.DictReader(open(os.path.join(HERE, "qm.csv"))))
ax = axs[1]
pw = [float(r["pw"]) for r in rows]
ax.plot(pw, [float(r["exact"]) for r in rows], color="#1f5fbf", lw=6, label="exact geodesic: 2p/√(1+p²)")
ax.plot(pw, [float(r["kaluza_linear"]) for r in rows], color="#c0392b", lw=3, ls="--", label="Kaluza eq. 8a: 2p")
ax.axhline(2, color="black", lw=2, ls=":", label="ceiling q/(m√G) = 2")
ax.set_ylim(0, 6)
ax.set_xlabel("momentum along the fifth dimension p_w")
ax.set_ylabel("charge ÷ mass, q/(m√G)")
ax.set_title("Charge is fifth-dimension motion — capped at 2\n(the electron needs 2×10²¹)", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13, loc="upper left")

rows = list(csv.DictReader(open(os.path.join(HERE, "machine.csv"))))
ax = axs[2]
t = [float(r["t"]) for r in rows]
ax.plot(t, [float(r["geodesic"]) for r in rows], color="#1f5fbf", lw=8, label="Kron geodesic (q1, q2, θ)")
ax.plot(t, [float(r["reduced"]) for r in rows], color="#f1c40f", lw=4, ls="--", label="rotor alone, fluxes conserved")
ax.plot(t, [float(r["kronspice"]) for r in rows], color="#c0392b", lw=2, ls=":", label="kronspice machine element")
ax.set_xlabel("time (s)")
ax.set_ylabel("rotor speed (rad/s)")
ax.set_title("Kron machine: winding charges as the fifth dimension", loc="left")
ax.grid(alpha=0.3)
ax.legend(fontsize=13, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "kaluza_kron.png"), dpi=90)
print("wrote", os.path.join(HERE, "kaluza_kron.png"))
