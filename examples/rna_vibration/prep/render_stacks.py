#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Charts and sound files for the stacked base pairs (stacks/ from the Go run).

Sound: every frequency divided by 10 billion (1 cm^-1 -> 3 Hz), loudness by
the square root of the infrared intensity, as in the single-base files.
No measured low-frequency band of an A-form RNA duplex was found to mark
on the chart. Run from examples/rna_vibration.
"""
import csv
import wave

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SR, K = 44100, 3.0
plt.rcParams.update({"font.size": 15, "font.family": "DejaVu Sans"})


def modes(key):
    return list(csv.DictReader(open("stacks/modes_%s.csv" % key.lower())))


def tone_file(rows, fn, step=0.06, chord=6.0):
    md = [(float(r["wavenumber_cm-1"]), float(r["ir_intensity_km_mol"])) for r in rows]
    imax = max(i for _, i in md)
    seq = []
    for w, i in sorted(md):
        t = np.arange(int(SR * step * 3)) / SR
        seq.append((0.25 + 0.75 * np.sqrt(i / imax)) * np.sin(2 * np.pi * K * w * t) * np.exp(-t / 0.06))
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
    y = 0.8 * np.concatenate([arp, np.zeros(int(SR * 0.5)), ch])
    with wave.open(fn, "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(SR)
        f.writeframes((y * 32767).astype("<i2").tobytes())


def low_file(rows, fn, dur=6.0):
    """Only the vibrations below 200 cm^-1 (the stacking and pairing motions)."""
    md = [(float(r["wavenumber_cm-1"]), float(r["ir_intensity_km_mol"])) for r in rows if float(r["wavenumber_cm-1"]) < 200]
    t = np.arange(int(SR * dur)) / SR
    rng = np.random.default_rng(2)
    y = sum((0.3 + np.sqrt(i)) * np.sin(2 * np.pi * K * w * t + rng.uniform(0, 2 * np.pi)) for w, i in md)
    y *= np.minimum(1, np.minimum(t / 0.2, (dur - t) / 0.8))
    y = 0.8 * y / np.max(np.abs(y))
    with wave.open(fn, "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(SR)
        f.writeframes((y * 32767).astype("<i2").tobytes())


def chart():
    fig, axs = plt.subplots(2, 1, figsize=(14, 11), facecolor="white")
    cols = {1: "#7a7a7a", 2: "#1f5fbf", 3: "#1e8449", 4: "#c0392b"}
    for ax, fam in zip(axs, ("AU", "GC")):
        for n, key in enumerate([fam, fam + "2", fam + "3", fam + "4"], start=1):
            r = modes(key)
            w = [float(x["wavenumber_cm-1"]) for x in r if float(x["wavenumber_cm-1"]) < 260]
            I = [float(x["ir_intensity_km_mol"]) for x in r if float(x["wavenumber_cm-1"]) < 260]
            rigid = [x["stack_motion"] != "" for x in r if float(x["wavenumber_cm-1"]) < 260]
            y = 5 - n
            for wi, ri in zip(w, rigid):
                ax.plot([wi, wi], [y - 0.35, y + 0.35], color=cols[n], lw=3.2 if ri else 1.6, alpha=1 if ri else 0.55)
            ax.text(-4, y, "%d pair%s" % (n, "" if n == 1 else "s"), ha="right", va="center", fontsize=15, color=cols[n], fontweight="bold")
        ax.set_xlim(0, 260)
        ax.set_ylim(0.4, 5.1)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="x", color="#cccccc")
        ax.set_title(("A·U" if fam == "AU" else "G·C") + " pairs stacked (A-form step of PDB 1EKA): low-frequency vibrations", fontsize=18, fontweight="bold", loc="left")
    axs[1].set_xlabel("Wavenumber, cm⁻¹   (in the sound files × 3 = hertz)", fontsize=16)
    fig.text(0.01, 0.005, "Thick lines: the pairs move as solid blocks against each other (stacking motions).  Thin lines: motion inside the pairs.\n"
             "No measured low-frequency band of an A-form RNA duplex was found to compare with.",
             fontsize=12)
    fig.tight_layout(rect=(0.06, 0.06, 1, 1))
    fig.savefig("summary/stack_low_modes.png", dpi=110, facecolor="white")

    # Band splitting at high frequency: the N–H and C=O regions.
    fig, axs = plt.subplots(2, 1, figsize=(14, 9), facecolor="white", sharex=False)
    for ax, (lo, hi, title) in zip(axs, ((1650, 1900, "Ring and C=O stretch region"), (3200, 3400, "N–H stretch region"))):
        for n, key in enumerate(["AU", "AU2", "AU3", "AU4"], start=1):
            r = [x for x in modes(key) if lo <= float(x["wavenumber_cm-1"]) <= hi]
            imax = max(float(x["ir_intensity_km_mol"]) for x in modes(key))
            for x in r:
                w, I = float(x["wavenumber_cm-1"]), float(x["ir_intensity_km_mol"])
                ax.plot([w, w], [5 - n - 0.4, 5 - n - 0.4 + 0.8 * np.sqrt(I / imax)], color=cols[n], lw=2.6)
            ax.text(lo - 3, 5 - n, "%d A·U" % n, ha="right", va="center", fontsize=14, color=cols[n], fontweight="bold")
        ax.set_xlim(lo, hi)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.set_title(title + ": each line splits into a cluster as pairs are added", fontsize=16, fontweight="bold", loc="left")
        ax.grid(axis="x", color="#dddddd")
    axs[1].set_xlabel("Wavenumber, cm⁻¹", fontsize=15)
    fig.tight_layout(rect=(0.06, 0, 1, 1))
    fig.savefig("summary/stack_band_splitting.png", dpi=110, facecolor="white")


def chord_file(rows, fn, dur=8.0):
    """All vibrations at once, as one chord (no tone-by-tone opening)."""
    md = [(float(r["wavenumber_cm-1"]), float(r["ir_intensity_km_mol"])) for r in rows]
    t = np.arange(int(SR * dur)) / SR
    rng = np.random.default_rng(1)
    y = sum(np.sqrt(i) * np.sin(2 * np.pi * K * w * t + rng.uniform(0, 2 * np.pi)) for w, i in md)
    y *= np.minimum(1, np.minimum(t / 0.15, (dur - t) / 0.8))
    y = 0.8 * y / np.max(np.abs(y))
    with wave.open(fn, "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(SR)
        f.writeframes((y * 32767).astype("<i2").tobytes())


if __name__ == "__main__":
    chart()
    single = {"a": "adenine", "u": "uracil", "g": "guanine", "c": "cytosine", "au": "AU_pair", "gc": "GC_pair", "gu": "GU_wobble_pair"}
    for k, n in single.items():
        chord_file(list(csv.DictReader(open("modes/%s.csv" % k))), "sound/%s_chord.wav" % n)
    for fam in ("AU", "GC"):
        for n in (2, 3, 4):
            chord_file(modes("%s%d" % (fam, n)), "sound/%s_stack_%d_chord.wav" % (fam, n))
        low_file(modes(fam), "sound/%s_low_1pair.wav" % fam)
        low_file(modes(fam + "4"), "sound/%s_low_4pairs.wav" % fam)
    # Runs and mixed stacks, named by sequence as in examples/dna_vibration.
    for key, name in (("polyA2", "common_AA-UU_2pairs"), ("polyA3", "common_AAA-UUU_3pairs"), ("polyA4", "common_AAAA-UUUU_4pairs"),
                      ("ApU", "mixed_AU_2pairs"), ("AU3", "mixed_AAU-AUU_3pairs"), ("AU4", "mixed_AAUU_4pairs"),
                      ("GU2", "tandem_UG-GU_wobble_2pairs")):
        chord_file(modes(key), "sound/%s_chord.wav" % name)
