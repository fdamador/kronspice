#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Charts of the four RNA bases against measurement, and the single-molecule
sound files (modes/ and summary/ from the Go run, summary/measured.csv from
prep/measured.py). Run from examples/rna_vibration.

Sound: every frequency divided by 10 billion (1 cm^-1 -> 3 Hz), loudness by
the square root of the infrared intensity: each vibration as a short tone
from lowest to highest, then all of them together.
"""
import csv
import wave

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SR, K = 44100, 3.0
NAMES = {"A": "adenine", "U": "uracil", "G": "guanine", "C": "cytosine"}
LABEL = {"A": "A (9-methyladenine)", "U": "U (1-methyluracil)", "G": "G (9-methylguanine)", "C": "C (1-methylcytosine)"}
COL = {"A": "#1f5fbf", "U": "#c0392b", "G": "#1e8449", "C": "#8e44ad"}


def modes(key):
    return [(float(r["wavenumber_cm-1"]), float(r["ir_intensity_km_mol"])) for r in csv.DictReader(open("modes/%s.csv" % key.lower(), encoding="utf-8"))]


def measured():
    rows = []
    for r in csv.DictReader(open("summary/measured.csv", encoding="utf-8")):
        if r["model_cm-1"]:
            rows.append(r)
    return rows


def tones(md, step=0.08, chord=5.0):
    imax = max(i for _, i in md)
    seq = []
    for w, i in sorted(md):
        t = np.arange(int(SR * step * 3)) / SR
        seq.append((0.25 + 0.75 * np.sqrt(i / imax)) * np.sin(2 * np.pi * K * w * t) * np.exp(-t / 0.08))
    n = int(SR * step)
    arp = np.zeros(n * len(seq) + len(seq[-1]))
    for k, s in enumerate(seq):
        arp[k * n:k * n + len(s)] += s
    arp /= np.max(np.abs(arp))
    t = np.arange(int(SR * chord)) / SR
    rng = np.random.default_rng(1)
    ch = sum(np.sqrt(i) * np.sin(2 * np.pi * K * w * t + rng.uniform(0, 2 * np.pi)) for w, i in md)
    ch *= np.minimum(1, np.minimum(t / 0.15, (chord - t) / 0.6))
    ch /= np.max(np.abs(ch))
    return 0.8 * np.concatenate([arp, np.zeros(int(SR * 0.5)), ch])


def write(fn, y):
    with wave.open(fn, "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(SR)
        f.writeframes((y * 32767).astype("<i2").tobytes())


def base_chart():
    meas = measured()
    key = {"Adenine": "A", "Uracil": "U", "Guanine": "G", "Cytosine": "C"}
    fig, axs = plt.subplots(4, 1, figsize=(13, 10), sharex=True, facecolor="white")
    for ax, b in zip(axs, "AUGC"):
        md = modes(b)
        imax = max(i for _, i in md)
        for w, i in md:
            ax.plot([w, w], [0, i / imax], color=COL[b], lw=2)
        for r in meas:
            if key[r["molecule"]] == b:
                ax.plot(float(r["measured_cm-1"]), 1.08, marker="v", color="black", ms=8)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([])
        ax.set_title(LABEL[b], loc="left", fontsize=12, color=COL[b], fontweight="bold")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    axs[-1].set_xlim(0, 3700)
    axs[-1].set_xlabel("wavenumber (cm⁻¹)")
    fig.suptitle("RNA bases as circuits: model vibrations (bar height = IR intensity) and laboratory band positions (triangles)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig("summary/base_vibrations.png", dpi=110, facecolor="white")


def mvm_chart():
    meas = measured()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 6), facecolor="white", gridspec_kw={"width_ratios": [1.3, 1]})
    labels = []
    for k, r in enumerate(meas):
        m, x = float(r["measured_cm-1"]), float(r["model_cm-1"])
        ax.plot(m, k, marker="v", color="black", ms=8)
        ax.plot(x, k, marker="o", color="#c0392b", ms=7)
        ax.plot([m, x], [k, k], color="#999999", lw=1)
        labels.append(r["molecule"] + " " + r["band"])
        bx.barh(k, float(r["difference_pct"]), color="#c0392b" if r["molecule"] == "Uracil" else "#1f5fbf")
    for a in (ax, bx):
        a.set_yticks(range(len(meas)))
        a.invert_yaxis()
    ax.set_yticklabels(labels, fontsize=9)
    bx.set_yticklabels([])
    ax.set_xlabel("wavenumber (cm⁻¹): triangle measured, dot model")
    bx.set_xlabel("model − measured (%)")
    bx.axvline(0, color="black", lw=0.8)
    fig.suptitle("Model (harmonic AMBER OL3, gas phase) against measurement; uracil rows in red", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("summary/model_vs_measured.png", dpi=110, facecolor="white")


if __name__ == "__main__":
    base_chart()
    mvm_chart()
    allb = []
    for b in "AUGC":
        y = tones(modes(b))
        write("sound/%s.wav" % NAMES[b], y)
        allb += [y, np.zeros(int(SR * 0.8))]
    write("sound/four_bases.wav", np.concatenate(allb))
    for k, n in (("AU", "AU_pair"), ("GC", "GC_pair"), ("GU", "GU_wobble_pair")):
        write("sound/%s.wav" % n, tones(modes(k)))
    print("ok")
