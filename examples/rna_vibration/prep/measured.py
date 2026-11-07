#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Write summary/measured.csv: published band positions against the model
(summary/bands.csv from the Go run). Run from examples/rna_vibration.

Adenine, guanine and cytosine are the same molecules as in
examples/dna_vibration (9-methyladenine, 9-methylguanine and
1-methylcytosine are capped identically), so their measured values and
sources are reused from that example's summary/measured.csv; only the model
column is this run's. Uracil values are new. No gas-phase or droplet IR
spectrum of an A·U or G·U pair, and no measured low-frequency mode of an
A-form RNA duplex, was found, so none is listed.
"""
import csv

bands = {(r["molecule"], r["band"]): float(r["model_cm-1"]) for r in csv.DictReader(open("summary/bands.csv", encoding="utf-8"))}
share = {(r["molecule"], r["band"]): float(r["energy_share"]) for r in csv.DictReader(open("summary/bands.csv", encoding="utf-8"))}

IVANOV = "Ivanov, Krasnokutski, Sheina, Blagoi, Biophysical Bulletin (Visnyk Kharkiv Univ.) 2(11), 14 (2002), Table 1 and text p. 16"
COLARUSSO = "Colarusso, Zhang, Guo, Bernath, Chem. Phys. Lett. 269, 39 (1997), values read from Colarusso's PhD thesis (Univ. Waterloo 1997) Table 4.1"
PLUTZER01 = "Plützer, Nir, de Vries, Kleinermanns, PCCP 3, 5466 (2001)"
NOWAK = "Nowak, Reva, Lopes Jesus, Lapinski, Fausto, PCCP 21, 22857 (2019), ESI Table S2"
MIN = "Min, Lee, Choi, Miller, Bull. Korean Chem. Soc. 30, 3039 (2009)"
NIR = "Nir, Janzen, Imhof, Kleinermanns, de Vries, PCCP 4, 732 (2002)"
CB = "Cavaluzzi & Borer, Nucleic Acids Res. 32, e13 (2004)"
ONI = "Onidas et al., J. Phys. Chem. B 106, 11367 (2002), via Markovitsi et al. 2012"
REUSE = "; value and source reused from examples/dna_vibration"

rows = [
    # molecule, band label, model key (molecule, band), measured, conditions, source
    ("Adenine", "NH₂ antisymmetric stretch", ("A", "A NH₂ antisymmetric stretch"), 3569, "9H-adenine, gas, jet-cooled IR–UV; model is 9-methyladenine" + REUSE, PLUTZER01),
    ("Adenine", "NH₂ symmetric stretch", ("A", "A NH₂ symmetric stretch"), 3452, "9H-adenine, gas, jet-cooled IR–UV" + REUSE, PLUTZER01),
    ("Uracil", "N3–H stretch", ("U", "U N3–H stretch"), 3428.6, "1-methyluracil, Ar matrix 12 K, main band (site component 3433.4)", IVANOV),
    ("Uracil", "N3–H stretch (uracil, gas)", ("U", "U N3–H stretch"), 3436, "uracil (not methylated), hot gas 325 °C, FT absorption/emission; model is 1-methyluracil", COLARUSSO),
    ("Uracil", "C2=O stretch", ("U", "U C2=O stretch"), 1737.5, "1-methyluracil, Ar matrix 12 K (quoted in the text)", IVANOV),
    ("Uracil", "C4=O stretch", ("U", "U C4=O stretch"), 1700, "1-methyluracil, Ar matrix 12 K (quoted in the text)", IVANOV),
    ("Guanine", "NH₂ antisymmetric stretch", ("G", "G NH₂ antisymmetric stretch"), 3543, "9-methylguanine amino-oxo, Ar matrix 14 K (site 3543/3535)" + REUSE, NOWAK),
    ("Guanine", "N1–H stretch", ("G", "G N1–H stretch"), 3440, "9-methylguanine amino-oxo, Ar matrix (3440/3435)" + REUSE, NOWAK),
    ("Guanine", "NH₂ symmetric stretch", ("G", "G NH₂ symmetric stretch"), 3430, "9-methylguanine amino-oxo, Ar matrix (3430/3427)" + REUSE, NOWAK),
    ("Guanine", "C6=O stretch", ("G", "G C6=O stretch"), 1752, "9-methylguanine amino-oxo, Ar matrix (1752/1741)" + REUSE, NOWAK),
    ("Cytosine", "NH₂ antisymmetric stretch", ("C", "C NH₂ antisymmetric stretch"), 3572.7, "cytosine keto form, helium nanodroplets 0.37 K; model is 1-methylcytosine" + REUSE, MIN),
    ("Cytosine", "NH₂ symmetric stretch", ("C", "C NH₂ symmetric stretch"), 3451.7, "cytosine keto form, helium nanodroplets" + REUSE, MIN),
    ("Cytosine", "C2=O stretch", ("C", "C C2=O stretch"), 1740, "cytosine keto, gas GC-FTIR 300 °C (as quoted)" + REUSE, NIR),
]

with open("summary/measured.csv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["molecule", "band", "measured_cm-1", "model_cm-1", "difference_pct", "conditions", "source"])
    for mol, band, key, meas, cond, src in rows:
        m = bands[key]
        cond = cond + "; model mode holds %.0f%% of the bond's energy" % (100 * share[key])
        w.writerow([mol, band, meas, "%.1f" % m, "%.1f" % (100 * (m - meas) / meas), cond, src])
    for mol, band, meas, cond, src in [
        ("Adenine", "UV absorption maximum (259 nm)", 38610, "dAMP in water pH 7; electronic, not in the vibration circuit" + REUSE, CB),
        ("Guanine", "UV absorption maximum (252 nm)", 39683, "dGMP in water pH 7; electronic, not in the vibration circuit" + REUSE, CB),
        ("Uracil", "UV absorption maximum (262 nm)", 38168, "pU (UMP) in water or 0.10 M phosphate pH 6.9–7.0, Table 1; electronic, not in the vibration circuit", CB),
        ("Cytosine", "UV absorption maximum (271 nm)", 36900, "dCMP in water pH 7; electronic, not in the vibration circuit" + REUSE, CB),
        ("Adenine", "Fluorescence maximum (306 nm)", 32680, "dAMP in water, quantum yield ~1e-4; electronic, not in the vibration circuit" + REUSE, ONI),
        ("Guanine", "Fluorescence maximum (334 nm)", 29940, "dGMP in water, quantum yield ~1e-4; electronic, not in the vibration circuit" + REUSE, ONI),
        ("Cytosine", "Fluorescence maximum (328 nm)", 30488, "dCMP in water, quantum yield ~1e-4; electronic, not in the vibration circuit" + REUSE, ONI),
    ]:
        w.writerow([mol, band, meas, "", "", cond, src])
print("ok")
