"""Pictures for examples/forward_casimir (run `go run . -viz viz` first).

    python3 viz/render.py

Writes viz/forward_casimir.png: Casimir pressure for ideal mirrors and real
aluminium, and the servo voltage and charge over a stroke. Light
background and large type.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

fig, axs = plt.subplots(1, 2, figsize=(18, 8))
rows = list(csv.DictReader(open(os.path.join(HERE, "casimir.csv"))))
a = [float(r["gap_nm"]) for r in rows]
ax = axs[0]
ax.plot(a, [float(r["ideal"]) for r in rows], color="#7f8c8d", lw=6, label="perfect mirrors (Forward eq. 2)")
ax.plot(a, [float(r["al"]) for r in rows], color="#1f5fbf", lw=4, label="aluminium (Lifshitz, λp = 107 nm)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("gap (nm)")
ax.set_ylabel("attractive pressure (Pa)")
ax.set_title("Casimir pull between two plates", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13)

st = {}
for r in csv.DictReader(open(os.path.join(HERE, "stroke.csv"))):
    st.setdefault(r["case"], []).append((float(r["gap_nm"]), float(r["volt"])))
ax = axs[1]
for (name, pts), col, ls, lw in zip(st.items(), ["#7f8c8d", "#1f5fbf", "#e67e22", "#c0392b"], ["-", "-", "--", "--"], [6, 6, 3, 3]):
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=col, ls=ls, lw=lw, label=name)
ax.set_xscale("log")
ax.set_xlabel("gap (nm)")
ax.set_ylabel("servo voltage on the leaves (V)")
ax.set_title("Forward's servo: voltage needed to balance the pull", loc="left")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "forward_casimir.png"), dpi=100)
print("wrote", os.path.join(HERE, "forward_casimir.png"))
