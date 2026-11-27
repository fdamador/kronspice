"""Pictures for examples/pulsed_bifilar (run the Go example first).

    python3 viz/render.py [viz-dir] [run_output.txt]

Writes into viz/:
  coil_voltage.png         the two solenoids (ordinary, bifilar), every turn
                           coloured by its voltage at the moment of largest
                           voltage, with the voltage along the wire and
                           between neighbouring turns
  coil_voltage_shapes.png  the same for the conical and flat-spiral bifilar coils
  coil_voltage_counter.png Tesla (same-direction) and opposite-current bifilar solenoids
  coil_voltage_parallel.png the two parallel bifilar solenoids (same direction, opposite)
  coil_current.png         the series solenoids (ordinary, Tesla bifilar, opposite-current
                           bifilar), every turn coloured by the current circulating in it
                           (red one way round, blue the other) at the moment of largest
                           current, with each section's current along the wire
  coil_current_parallel.png the same for the two parallel bifilar solenoids
  coil_amps.png            wire and primary current building up pulse by pulse, all seven
  coil_vip.png             volts, amps and power at f1, all seven, side by side
  vi_frequency.png         V and I against frequency (AC, per volt of drive) for the
                           five solenoid windings, and V/I
  pulse_sweep.png          coil-end voltage against pulse rate (rate / f1), all seven
  ring_up.png              the coil-end voltage building up pulse by pulse, all seven
  ring_up.gif              the two solenoids glowing as the pulses build them up
  ring_up_shapes.gif       the same for the cone and the flat spiral
Light background and large type throughout.
"""
import csv
import os
import re
import sys
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import cm, colors  # noqa: E402
from PIL import Image  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
VIZ = sys.argv[1] if len(sys.argv) > 1 else HERE
RUN = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "..", "run_output.txt")
ONLY = sys.argv[3:]  # optional figure-function names, e.g. vi_frequency

plt.rcParams.update({"font.size": 15, "axes.titlesize": 17, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})
CMAP = matplotlib.colormaps["coolwarm"]
GLOW = matplotlib.colormaps["inferno"]

# tag, label, colour, the heading the Go program prints for it
COILS = [
    ("ordinary", "Ordinary solenoid", "#1f5fbf", "ordinary cylindrical solenoid"),
    ("bifilar", "Bifilar solenoid", "#c0392b", "Tesla bifilar cylindrical solenoid"),
    ("cone_bifilar", "Bifilar cone", "#27864a", "Tesla bifilar conical"),
    ("flat_bifilar", "Bifilar flat spiral", "#8e44ad", "Tesla bifilar flat spiral"),
    ("counter", "Series bifilar, opposite currents", "#d68910", "opposite-current bifilar cylindrical solenoid"),
    ("parallel", "Parallel bifilar, same direction", "#16a085", "parallel bifilar (same direction)"),
    ("parallel_counter", "Parallel bifilar, opposite currents", "#7f8c8d", "parallel bifilar (opposite currents)"),
]


def label_of(tag):
    return next(c[1] for c in COILS if c[0] == tag)


def load(tag, kind):
    """kind: envelope (time, coil-end volts), glow (time, |V| per position per
    frame), peak (time, V per position at the largest coil-end voltage) or
    geom (pos, radius, height)."""
    with open(os.path.join(VIZ, f"{tag}_{kind}.csv")) as f:
        rows = list(csv.reader(f))
    data = np.array([[float(x) for x in r] for r in rows[1:]])
    return data[:, 0], data[:, 1:]


def geometry(tag):
    _, g = load(tag, "geom")
    return g[:, 0], g[:, 1]  # radius, height per position


def wire_order(tag):
    """Physical position of each turn in wire order (from geom's turn column)."""
    _, g = load(tag, "geom")
    return np.argsort(g[:, 2]).astype(int)


def positions(n, every):
    """Neighbouring pairs (p, p+1) every `every` turns, so both bifilar strands show."""
    return sorted(set(range(0, n, every)) | set(range(1, n, every)))


def draw_coil(ax, tag, values, cmap, norm, every, lw=1.6, ticks=True, azim=35):
    r, z = geometry(tag)
    th = np.linspace(0, 2 * np.pi, 72)
    flat = np.ptp(z) < 1e-9
    for p in positions(len(r), every):
        ax.plot(r[p] * np.cos(th) * 100, r[p] * np.sin(th) * 100, np.full_like(th, z[p] * 100),
                color=cmap(norm(values[p])), lw=lw)
    rmax = r.max() * 100
    height = max(np.ptp(z) * 100, 1)
    ax.set_xlim(-rmax, rmax), ax.set_ylim(-rmax, rmax)
    if flat:
        ax.set_zlim(-rmax / 2, rmax / 2)
        ax.set_box_aspect((1, 1, 0.5))
        ax.view_init(elev=35, azim=azim)
    else:
        ax.set_box_aspect((1, 1, min(3, height / (2 * rmax))))
        ax.view_init(elev=14, azim=azim)
    ax.set_xticks([]), ax.set_yticks([])
    if not ticks or flat:
        ax.set_zticks([])
    else:
        ax.set_zlabel("height (cm)", labelpad=8)


def coil_voltage(tags, out, title):
    fig = plt.figure(figsize=(16, 11))
    for col, tag in enumerate(tags):
        _, v = load(tag, "peak")
        snap = v[0]
        vmax = np.max(np.abs(snap))
        norm = colors.Normalize(-vmax, vmax)
        ax = fig.add_subplot(2, 2, 1 + col, projection="3d")
        draw_coil(ax, tag, snap, CMAP, norm, every=2)
        ax.set_title(f"{label_of(tag)}\nturns coloured by voltage (±{vmax:,.0f} V)", pad=2)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=CMAP), ax=ax, shrink=0.55, pad=0.08, label="volts")
        ax2 = fig.add_subplot(2, 2, 3 + col)
        n = len(snap)
        along = snap[wire_order(tag)]
        ax2.plot(np.arange(n), along, color="#1f5fbf", lw=2.5, label="voltage along the wire")
        adj = np.abs(np.diff(snap))  # between physically adjacent turns
        ax2.plot(np.arange(n - 1), adj, color="#c0392b", lw=2, label="between neighbouring turns")
        ax2.set_xlabel("turn (along the wire for blue, inside/bottom outward for red)")
        ax2.set_ylabel("volts")
        ax2.grid(alpha=0.3)
        ax2.legend(loc="best", fontsize=13)
        ax2.set_title(f"largest neighbour-to-neighbour voltage {adj.max():,.0f} V")
    fig.suptitle(title, fontsize=19)
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, out), dpi=100)
    plt.close(fig)


def parse_run():
    with open(RUN) as f:
        txt = f.read()
    out = {}
    for block in re.split(r"\n== ", txt)[1:]:
        head = block.split(" winding ==")[0]
        tag = next(c[0] for c in COILS if head.startswith(c[3]))
        f1 = float(re.search(r"f1 = ([0-9.e+]+) Hz", block).group(1))
        num = r"\s+([0-9.e+-]+)"
        rows = re.findall(r"^" + num * 11 + r"\s*$", block, re.M)
        # rate/f1, coil ends V, adjacent turns V, wire A, primary A, in W, heat W
        out[tag] = (f1, [tuple(float(r[i]) for i in (0, 2, 3, 4, 5, 6, 7)) for r in rows])
    return out


def pulse_sweep():
    data = parse_run()
    fig, ax = plt.subplots(figsize=(12, 7.5))
    for tag, label, colr, _ in COILS:
        f1, rows = data[tag]
        ax.plot([r[0] for r in rows], [r[1] for r in rows], "o-", color=colr, lw=2.5, ms=8,
                label=f"{label} (f1 = {f1 / 1e3:,.1f} kHz)")
    for r, t in ((1 / 3, "f1/3"), (0.5, "f1/2"), (1, "f1")):
        ax.axvline(r, color="#888", ls=":", lw=1.5)
        ax.text(r, 1.2e4, t, ha="center")
    ax.set_yscale("log")
    ax.set_ylim(0.5, 2e4)
    ax.set_xlabel("pulse rate ÷ the coil's resonance f1")
    ax.set_ylabel("voltage between the coil's ends (V)")
    ax.set_title("One-sign 100 V pulses: build-up only in step with resonance")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=13, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, "pulse_sweep.png"), dpi=100)
    plt.close(fig)


def ring_up():
    fig, axs = plt.subplots(2, 4, figsize=(26, 10))
    axs.flat[-1].axis("off")
    for ax, (tag, label, colr, _) in zip(axs.flat, COILS):
        m, e = load(tag, "envelope")
        e = e[:, 0]
        ax.fill_between(m * 1e6, -e, e, color=colr, alpha=0.35)
        ax.plot(m * 1e6, e, color=colr, lw=2.5)
        ax.set_title(textwrap.fill(f"{label}: coil-end voltage, pulses at f1", 34), fontsize=15)
        ax.set_xlabel("time (µs)")
        ax.set_ylabel("volts")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, "ring_up.png"), dpi=100)
    plt.close(fig)


def ring_up_gif(tags, out):
    sets = [(tag, label_of(tag)) + load(tag, "glow") for tag in tags]
    images = []
    for f in range(len(sets[0][2])):
        fig = plt.figure(figsize=(12, 8))
        for i, (tag, label, t, env) in enumerate(sets):
            ax = fig.add_subplot(1, 2, 1 + i, projection="3d")
            top = env.max()
            norm = colors.Normalize(0, top)
            draw_coil(ax, tag, 0.15 * top + 0.85 * env[f], GLOW, norm, every=6, lw=2, ticks=False, azim=35 + 2 * f)
            ax.set_title(f"{label}\n{t[f] * 1e6:,.0f} µs: up to {env[f].max():,.0f} V", fontsize=16)
        fig.suptitle("Pulses in step with resonance build the coil up (brighter = higher voltage)", fontsize=16)
        fig.tight_layout()
        fig.canvas.draw()
        img = Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).convert("RGB")
        images.append(img.resize((img.width * 3 // 4, img.height * 3 // 4)))
        plt.close(fig)
    images += [images[-1]] * 12
    images[0].save(os.path.join(VIZ, out), save_all=True, append_images=images[1:], duration=120, loop=0)


def coil_current(tags, out, title):
    """Turns coloured by the current circulating in them; below, each section's current."""
    fig = plt.figure(figsize=(7 * len(tags), 11))
    for col, tag in enumerate(tags):
        _, cur = load(tag, "current")  # wire_amps, circulating_amps, peak_amps per section
        circ = cur[:, 1]
        r, _ = geometry(tag)
        n = len(r)
        turn = load(tag, "geom")[1][:, 2].astype(int)
        vals = circ[turn * len(circ) // n]  # each position takes its section's current
        imax = max(np.max(np.abs(circ)), 1e-12)
        norm = colors.Normalize(-imax, imax)
        ax = fig.add_subplot(2, len(tags), 1 + col, projection="3d")
        draw_coil(ax, tag, vals, CMAP, norm, every=2)
        ax.set_title(f"{label_of(tag)}\nturns coloured by circulating current (±{imax:.3g} A)", pad=2)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=CMAP), ax=ax, shrink=0.55, pad=0.08,
                     label="amps (red one way round, blue the other)")
        ax2 = fig.add_subplot(2, len(tags), 1 + len(tags) + col)
        s = np.arange(len(circ))
        ax2.bar(s - 0.2, cur[:, 0], width=0.4, color="#1f5fbf", label="along the wire")
        ax2.bar(s + 0.2, circ, width=0.4, color="#c0392b", label="circulating (sense of the turns)")
        ax2.plot(s, cur[:, 2], "k_", ms=14, mew=2, label="peak |I|, last 10 %")
        ax2.axhline(0, color="#555", lw=1)
        ax2.set_xlabel("section along the wire (20 sections)")
        ax2.set_ylabel("amps")
        ax2.grid(alpha=0.3)
        ax2.legend(loc="best", fontsize=11)
        ax2.set_title(f"mean circulating current {np.mean(circ):+.3g} A")
    fig.suptitle(title, fontsize=19)
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, out), dpi=100)
    plt.close(fig)


def coil_amps():
    fig, axs = plt.subplots(2, 4, figsize=(26, 10))
    axs.flat[-1].axis("off")
    for ax, (tag, label, colr, _) in zip(axs.flat, COILS):
        m, e = load(tag, "ienv")
        ax.plot(m * 1e6, e[:, 0], color=colr, lw=2.5, label="coil wire (largest section)")
        ax.plot(m * 1e6, e[:, 1], color="#555", lw=2, ls="--", label="primary")
        ax.set_title(textwrap.fill(f"{label}: current, pulses at f1", 34), fontsize=15)
        ax.set_xlabel("time (µs)")
        ax.set_ylabel("amps (peak per period)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=11, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, "coil_amps.png"), dpi=100)
    plt.close(fig)


def coil_vip():
    data = parse_run()
    labels = [textwrap.fill(c[1], 13) for c in COILS]
    cols = [c[2] for c in COILS]
    at_f1 = [next(r for r in data[c[0]][1] if abs(r[0] - 1) < 1e-9) for c in COILS]
    panels = [("coil-end voltage (V)", [r[1] for r in at_f1], True),
              ("neighbouring-turn voltage (V)", [r[2] for r in at_f1], True),
              ("peak wire current (A)", [r[3] for r in at_f1], False),
              ("peak primary current (A)", [r[4] for r in at_f1], False),
              ("mean power from the source (W)", [r[5] for r in at_f1], False),
              ("mean power into heat (W)", [r[6] for r in at_f1], False)]
    fig, axs = plt.subplots(2, 3, figsize=(24, 12))
    for ax, (name, vals, log) in zip(axs.flat, panels):
        ax.bar(range(len(vals)), vals, color=cols)
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.3g}", ha="center", va="bottom", fontsize=12)
        if log:
            ax.set_yscale("log")
        ax.set_xticks(range(len(vals)), labels, fontsize=10)
        ax.set_title(name)
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle("Pulses at each coil's own f1 (100 V, 50 Ω, through the 5-turn primary): volts, amps and power",
                 fontsize=18)
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, "coil_vip.png"), dpi=100)
    plt.close(fig)


def vi_frequency():
    with open(os.path.join(VIZ, "vi_sweep.csv")) as f:
        rows = list(csv.reader(f))
    head = rows[0]
    d = np.array([[float(x) for x in r] for r in rows[1:]])
    col = {h: d[:, i] for i, h in enumerate(head)}
    f = col["freq"]
    sets = [(c[0], c[1], c[2]) for c in COILS if c[0] + "_v_ends" in col]
    fig, axs = plt.subplots(3, 1, figsize=(14, 17), sharex=True)
    for k, label, c in sets:
        v, i = col[k + "_v_ends"], col[k + "_i_wire"]
        pk = int(np.argmax(v))
        axs[0].loglog(f, v, color=c, lw=2.5, label=f"{label}: peak {v[pk]:.3g} V at {f[pk] / 1e3:,.1f} kHz")
        axs[1].loglog(f, i * 1e3, color=c, lw=2.5, label=f"{label}: {i[pk] * 1e3:.3g} mA at that peak")
        axs[2].loglog(f, v / i, color=c, lw=2.5, label=f"{label}: V/I {v[pk] / i[pk]:,.3g} Ω at the peak")
        for ax in axs:
            ax.axvline(f[pk], color=c, ls=":", lw=1.3)
    axs[0].set_ylabel("voltage between the coil's ends (V)")
    axs[1].set_ylabel("largest wire current (mA)")
    axs[2].set_ylabel("V / I (Ω)")
    axs[2].set_xlabel("frequency (Hz)")
    axs[0].set_title("Per volt of drive on the 5-turn primary (AC): same L, so at resonance the current is alike\n"
                     "and the bifilar's voltage is lower by √(C_bifilar / C_ordinary)")
    for ax in axs:
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=11, loc="best")
    fig.tight_layout()
    fig.savefig(os.path.join(VIZ, "vi_frequency.png"), dpi=100)
    plt.close(fig)


if __name__ == "__main__" and ONLY:
    for name in ONLY:
        globals()[name]()
    print("wrote", ", ".join(ONLY), "to", VIZ)
elif __name__ == "__main__":
    coil_voltage(["ordinary", "bifilar"], "coil_voltage.png",
                 "Pulses in step with resonance (rate = f1): the voltage pattern on each solenoid")
    coil_voltage(["cone_bifilar", "flat_bifilar"], "coil_voltage_shapes.png",
                 "Conical and flat-spiral bifilar coils at f1: the voltage pattern")
    coil_voltage(["bifilar", "counter"], "coil_voltage_counter.png",
                 "Bifilar solenoid, strands in series: currents the same way (Tesla) and opposite")
    coil_voltage(["parallel", "parallel_counter"], "coil_voltage_parallel.png",
                 "Bifilar solenoid, strands in parallel: currents the same way and opposite")
    coil_current(["ordinary", "bifilar", "counter"], "coil_current.png",
                 "Current in the turns at f1, strands in series: ordinary, Tesla bifilar (same way), opposite currents")
    coil_current(["parallel", "parallel_counter"], "coil_current_parallel.png",
                 "Current in the turns at f1, strands in parallel: same way and opposite")
    coil_amps()
    coil_vip()
    vi_frequency()
    pulse_sweep()
    ring_up()
    ring_up_gif(["ordinary", "bifilar"], "ring_up.gif")
    ring_up_gif(["cone_bifilar", "flat_bifilar"], "ring_up_shapes.gif")
    print("wrote the pictures and animations to", VIZ)
