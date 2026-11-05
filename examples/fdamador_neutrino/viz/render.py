"""Chart and sound for examples/fdamador_neutrino.

    go run ./examples/fdamador_neutrino -viz examples/fdamador_neutrino/viz
    python3 examples/fdamador_neutrino/viz/render.py

Writes
  viz/resonance.png            the resonance chart: E-ring resonance curves by
                               permittivity, Q against permittivity (formula and
                               kronspice), the inverted ring lines' spectra, and the
                               ring held by the zero-point field
  viz/fdamador_neutrino.wav    the whole sound, five parts with gaps (44.1 kHz, 16-bit mono)
  viz/sound_*.wav              each part on its own
  viz/harmonics.wav            the ring lines' modes added one at a time (plain, then crossed)

Sound: every real frequency here is far above hearing (the ceramic ring 2.06 GHz,
the ring lines 20-160 MHz, the 0.1 eV neutrino ring 24 THz). Each sound is the
kronspice waveform slowed down by one factor, so its pitch lands in the audible
range and every frequency ratio, decay rate and Q is kept:
  1. ceramic ring, εr 100 (Q 40.5): kronspice ring-down, f0 → 440 Hz (×2.135e-7)
  2. ceramic ring, εr 1000 (Q 1140): kronspice ring-down, f0 → 440 Hz (×6.752e-7)
  3. the ring held by the zero-point field (Q 40.7): kronspice noise run, ω0 → 440 Hz
  4. plain inverted ring line: its four highest modes (kronspice AC peaks), top → 880 Hz
  5. crossed inverted ring line: the same, top → 880 Hz
Parts 4 and 5 are synthesized from the measured mode frequencies (equal loudness,
1.2 s decay); parts 1-3 are the simulated waveforms themselves.
"""
import csv
import os
import wave

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
SEQ = {30: "#86b6ef", 100: "#3987e5", 300: "#1c5cab", 1000: "#0d366b"}  # blue ramp 250/400/550/700
plt.rcParams.update({"font.size": 15, "axes.titlesize": 17, "axes.titleweight": "normal",
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb", "axes.spines.top": False, "axes.spines.right": False})


def load(name):
    with open(os.path.join(HERE, name)) as f:
        rows = list(csv.reader(f))[1:]
    return np.array([[float(x) for x in r] for r in rows])


def grid(ax):
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------- chart
fig, axs = plt.subplots(5, 1, figsize=(11, 25), gridspec_kw={"height_ratios": [1.1, 1, 0.8, 0.8, 1]})
fig.subplots_adjust(hspace=0.55, left=0.12, right=0.95, top=0.945, bottom=0.04)
fig.suptitle("fdamador neutrino: resonances of the inverted ring (kronspice)", fontsize=19, color=INK, y=0.985)

ax = axs[0]
qs = {r[0]: r[2] for r in load("ering_q.csv")}
for er, col in SEQ.items():
    d = load(f"ering_er{er:g}.csv")
    ax.plot(d[:, 0] * 100, d[:, 1], color=col, lw=2, label=f"εr {er:g}  (Q {qs[er]:.3g})")
ax.legend(frameon=False, fontsize=13, loc="upper right")
ax.set_xlim(-12, 12)
ax.set_xlabel("frequency offset from resonance (% of f₀)")
ax.set_ylabel("ring current / peak")
ax.set_title("1. Ceramic E-ring: resonance narrows as permittivity rises")
grid(ax)

ax = axs[1]
q = load("ering_q.csv")
er = np.logspace(np.log10(20), np.log10(1500), 100)
e = 8.8541878128e-12 * er * 2.5e-3 ** 2 / (2 * 10e-3)
lam = np.log(8 * 4) - 1.75
L = 1.25663706212e-6 * 10e-3 * lam
w0 = 1 / np.sqrt(L * e)
ka = w0 * 10e-3 / 299792458.0
qrad = 6 * lam / (np.pi * ka ** 3)
qd = 1 / 1e-4
qtot = 1 / (1 / qrad + 1 / qd)
ax.loglog(er, qtot, color=INK2, lw=2, label="formula: Q_rad ∝ εr^1.5, with tan δ = 1e-4")
ax.loglog(q[:, 0], q[:, 2], "o", ms=10, color=BLUE, markeredgecolor="#fcfcfb", mew=2, label="kronspice AC bandwidth")
for x, y in zip(q[:, 0], q[:, 2]):
    ax.annotate(f"{y:.3g}", (x, y), xytext=(8, -16), textcoords="offset points", color=INK2, fontsize=13)
ax.set_xlabel("relative permittivity εr")
ax.set_ylabel("Q")
ax.set_title("2. The permittivity is the confinement")
ax.legend(frameon=False, fontsize=13, loc="upper left")
grid(ax)

for ax, name, col, title, frac in [
        (axs[2], "plain", BLUE, "3. Inverted (CL) ring line, plain: whole wavelengths", ["1", "2", "3", "4"]),
        (axs[3], "crossed", ORANGE, "4. Inverted ring line, crossed: half wavelengths (720°, spin ½)", ["½", "1½", "2½", "3½"])]:
    d = load(f"dual_{name}.csv")
    m = load(f"dual_{name}_modes.csv")[:, 0]
    ax.semilogy(d[:, 0] / 1e6, d[:, 1], color=col, lw=1.6)
    for j, (f, lab) in enumerate(zip(m, frac)):
        k = np.argmin(np.abs(d[:, 0] - f))
        text = f"{lab} λ\n{f / 1e6:.1f} MHz" if j < 2 else f"{lab} λ"
        ax.annotate(text, (f / 1e6, d[k, 1]), xytext=(0, 8 if j != 3 else 34), textcoords="offset points",
                    ha="center", va="bottom", color=INK2, fontsize=12)
    ax.set_xlim(10, 200)
    ax.set_ylim(top=d[:, 1].max() * 200)
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("|Z| (Ω)")
    ax.set_title(title)
    grid(ax)

ax = axs[4]
z = load("zpf_trace.csv")
t, en = z[:, 0], z[:, 2]
period = 2 * np.pi
sel = t > 10 * 40.68
tp = t[sel] / period
ax.plot(tp, en[sel], color=BLUE, lw=0.6, alpha=0.8)
win = max(1, int(np.sum(sel) / 60))
ax.plot(tp, np.convolve(en[sel], np.ones(win) / win, mode="same"), color=INK, lw=2.2)
ax.axhline(0.5, color=ORANGE, lw=2, ls="--")
ax.text(tp[-1], 0.5, " ½ħω₀", color=INK, va="center", fontsize=14)
ax.text(tp[0], ax.get_ylim()[1] * 0.92, "blue: energy each instant   black: running mean", color=INK2, fontsize=12)
ax.set_xlabel("time (periods of the ring)")
ax.set_ylabel("stored energy (ħω₀)")
ax.set_title("5. Zero-point inflow: the ring holds ½ħω₀ on average, as noise")
grid(ax)

out = os.path.join(HERE, "resonance.png")
fig.savefig(out, dpi=110)
print("wrote", out)

# ---------------------------------------------------------------- sound
RATE = 44100


def resample(t_sec, y, dur):
    ts = np.arange(int(dur * RATE)) / RATE
    return np.interp(ts, t_sec, y)


def fade(x, ms=15):
    n = int(RATE * ms / 1000)
    r = np.linspace(0, 1, n)
    x[:n] *= r
    x[-n:] *= r[::-1]
    return x


def norm(x, peak=0.8):
    return x / np.max(np.abs(x)) * peak


parts = {}
for er, name in [(100, "ring_eps100"), (1000, "ring_eps1000")]:
    d = load(f"ringdown_er{er}.csv")
    t_sec = d[:, 0] / 440.0  # one period of the ring = one period at 440 Hz
    parts[name] = fade(norm(resample(t_sec, d[:, 1], t_sec[-1])))

z = load("zpf_trace.csv")
sel = z[:, 0] > 10 * 40.68
t_sec = (z[sel, 0] - z[sel, 0][0]) / (2 * np.pi) / 440.0
parts["zero_point_ring"] = fade(norm(resample(t_sec, z[sel, 1], t_sec[-1])), ms=60)

for name in ["plain", "crossed"]:
    m = load(f"dual_{name}_modes.csv")[:, 0]
    f_audio = m / m.max() * 880.0
    ts = np.arange(int(2.6 * RATE)) / RATE
    x = sum(np.sin(2 * np.pi * f * ts) for f in f_audio) * np.exp(-ts / 1.2)
    parts[f"ring_line_{name}"] = fade(norm(x))
    print(f"  {name} ring line chord: " + ", ".join(f"{f:.1f} Hz" for f in f_audio))


def write(path, x):
    pcm = (np.clip(x, -1, 1) * 32767).astype("<i2")
    with wave.open(path, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(RATE)
        w.writeframes(pcm.tobytes())


gap = np.zeros(int(0.7 * RATE))
order = ["ring_eps100", "ring_eps1000", "zero_point_ring", "ring_line_plain", "ring_line_crossed"]
for k, name in enumerate(order, 1):
    write(os.path.join(HERE, f"sound_{k}_{name}.wav"), parts[name])
allx = np.concatenate([np.concatenate([parts[n], gap]) for n in order])
write(os.path.join(HERE, "fdamador_neutrino.wav"), allx)
print("wrote", os.path.join(HERE, "fdamador_neutrino.wav"), f"({len(allx) / RATE:.1f} s)")
for n in order:
    print(f"  {n}: {len(parts[n]) / RATE:.2f} s")

# ---------------------------------------------------------------- harmonics
# The ring lines' modes one at a time: each step adds the next mode and
# holds the growing chord for 1.6 s, first the plain ring (undertones
# 1, ½, ⅓, ¼ of the top mode), then the crossed ring (1, ⅓, ⅕, ⅐: every
# even member gone). The top mode is 880 Hz in both.
step = 1.6
hx = []
for name in ["plain", "crossed"]:
    m = np.sort(load(f"dual_{name}_modes.csv")[:, 0])[::-1]
    f_audio = m / m.max() * 880.0
    ts = np.arange(int(step * RATE)) / RATE
    for n in range(1, len(f_audio) + 1):
        x = sum(np.sin(2 * np.pi * f * ts) for f in f_audio[:n]) / np.sqrt(n)  # about equal loudness per step
        hx.append(fade(x, ms=25))
    hx.append(np.zeros(int(0.8 * RATE)))
harm = np.concatenate(hx)
harm = harm / np.max(np.abs(harm)) * 0.8
write(os.path.join(HERE, "harmonics.wav"), harm)
print("wrote", os.path.join(HERE, "harmonics.wav"), f"({len(harm) / RATE:.1f} s)")
