# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""power_filter report: harmonic filters for a six-pulse rectifier, scored by geometric-algebra power (M9f).

python examples/power_filter/viz/make_report.py examples/power_filter/power_filter_report.pdf

Numbers come from README.md (default run, 10 x 24, seed 1) and doc/M9f.md
(the checks of the GA power code and the longer 30 x 20 run). The example has
no figures and no run_output.txt. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/power_filter", "Harmonic filters for a six-pulse rectifier",
               "Can a GA scored by geometric-algebra power raise a rectifier's power factor? — milestone M9f, "
               "CLAUDE.md §3b item 3",
               "CLAUDE.md §3b, M9f", "10 generations × 24, seed 1")

    r.h1("Summary")
    r.p("A three-phase diode bridge draws distorted current from its supply. This example evolves two series L–C "
        "filter branches per phase at the point of common coupling (PCC) so that the power factor rises. The "
        "fitness is a power-quality measure computed with geometric algebra (GA): the geometric product of the "
        "harmonic voltage and current multivectors splits into active, reactive and distortion parts, and names "
        "the pairs of harmonics that carry the distortion. The result is checked against IEEE Std 1459 and "
        "against the quaternion instantaneous power of the three phases.")
    r.bullets([
        "<b>Power factor 0.9033 → 0.9559</b> with the default run (10 generations × 24, 222 evaluations); total active "
        "power 7389 W → 7673 W.",
        "<b>Evolved branches:</b> L5 = 3.9 mH with C5 = 82 µF (281 Hz) and L7 = 0.15 mH with C7 = 10 µF (4109 Hz) — not at "
        "300 Hz and 420 Hz.",
        "<b>Current distortion falls:</b> D<sub>I</sub> 1018 → 784 VA and THD<sub>I</sub> 40.4 % → 30.6 % on phase a, while "
        "THD<sub>V</sub> rises from 9.9 % to 12.8 %.",
        "<b>The GA names the culprits:</b> the largest distortion pair 1–5 falls from 961 VA to 578 VA; fundamental "
        "reactive power swings from +515 var to −73 var.",
        "<b>Twins agree:</b> GA and IEEE 1459 per phase to 0.3 % in S; the quaternion mean power equals ΣP (7389 W and 7673 W).",
        "<b>A longer run</b> (30 × 20, 562 evaluations, <i>doc/M9f.md</i>) reached PF 0.957 with branches at 284 Hz and 673 Hz.",
    ])

    r.h1("Background")
    r.p("CLAUDE.md §3b item 3 describes geometric-algebra power theory: for non-sinusoidal waveforms each harmonic's "
        "cosine and sine components are basis vectors of a geometric algebra, voltage V and current I become "
        "multivectors, and the geometric product M = V·I splits into a scalar (active power) and bivector parts "
        "(reactive and distortion power), each tied to a specific pair of harmonics. These measures are meant to "
        "feed GA fitness functions for power-quality filters. Milestone M9f builds this in <i>analysis/power.go</i> and "
        "checks it against IEEE Std 1459-2010; this example is its demonstration. No hypothesis is involved.")

    r.h1("The model")
    r.bullets([
        "<b>Supply:</b> 60 Hz, 120 V (phase), three-phase, with 0.05 Ω and 0.5 mH per line.",
        "<b>Load:</b> a diode bridge into a 2200 µF, 10 Ω DC load.",
        "<b>Filter:</b> at the PCC each phase carries two series L–C branches to a grounded star point, nominally for "
        "the 5th and 7th harmonics, each with a fixed branch resistance of 0.05 Ω.",
        "<b>Genome:</b> L5, L7 (0.1–10 mH) and C5, C7 (10–1000 µF), all on the E12 series.",
        "<b>Harmonic analysis:</b> for each phase the last two of 14 simulated cycles are resolved into harmonics up to "
        "the 25th (<i>analysis.HarmonicsOf</i>).",
        "<b>GA power:</b> the geometric product M = V·I in Cl(2H), with H the number of harmonics "
        "(<i>analysis.GeometricPower</i>), gives the scalar part P, the same-harmonic bivectors Q and the "
        "cross-harmonic bivectors D, grouped by pair of harmonics.",
        "<b>Fitness:</b> 1 − PF, with PF = ΣP/ΣS over the three phases (S the apparent power).",
    ])

    r.h1("Results")
    r.p("Default run: 10 generations × 24, 222 evaluations. Per-phase quantities are for phase a. "
        "Q<sub>1</sub> is the fundamental reactive power, D<sub>I</sub> the current-distortion power of IEEE 1459 "
        "(|V<sub>1</sub>I<sub>H</sub>|), THD total harmonic distortion.")
    r.table([
        ["", "No filter", "Evolved"],
        ["PF", "0.9033", "<b>0.9559</b>"],
        ["P (three phases)", "7389 W", "7673 W"],
        ["Phase a: S / P / Q<sub>1</sub>", "2731 VA / 2463 W / +515 var", "2701 VA / 2559 W / −73 var"],
        ["D<sub>I</sub>", "1018 VA", "784 VA"],
        ["THD<sub>I</sub> / THD<sub>V</sub>", "40.4 % / 9.9 %", "30.6 % / 12.8 %"],
        ["Largest distortion pairs", "1–5 961 VA, 1–7 333 VA, 1–11 192 VA", "1–5 578 VA, 1–7 461 VA, 1–11 156 VA"],
    ], [0.3, 0.35, 0.35])
    r.p("Evolved branches: L5 3.9 mH with C5 82 µF (281 Hz) and L7 0.15 mH with C7 10 µF (4109 Hz).")
    r.h2("Why the branches are not at 300 Hz and 420 Hz")
    r.bullets([
        "Maximizing PF trades harmonic absorption against the leading fundamental reactive power the capacitors add.",
        "The 5th-harmonic branch is tuned a little low, which removes most of the 1–5 distortion and swings Q<sub>1</sub> "
        "from +515 var to −73 var.",
        "The second branch in this run is mostly a capacitor.",
    ])
    r.h2("Longer run (doc/M9f.md)")
    r.table([
        ["Run", "PF", "THD<sub>I</sub>", "Largest distortion pairs (phase a)"],
        ["No filter", "0.903", "40 %", "1–5 961 VA, 1–7 333 VA, 1–11 192 VA"],
        ["Evolved, 30 × 20, 562 evaluations", "0.957", "28 %", "1–5 611 VA, 1–7 402 VA, 1–13 151 VA"],
    ], [0.32, 0.1, 0.1, 0.48])
    r.p("Its branches sit at 284 Hz and 673 Hz, Q<sub>1</sub> goes from +515 var to −226 var per phase, and the "
        "quaternion mean power equals ΣP (7389 W and 7736 W). The one-line summary in examples/README.md "
        "(\"PF 0.903 to 0.957\") matches this longer run.")

    n0 = len(r.s)
    r.h1("Validation")
    r.h2("In this example")
    r.bullets([
        "The GA and IEEE 1459 columns agree per phase, S to 0.3 % (the commutation notches above the 25th harmonic are "
        "outside the GA's harmonic set).",
        "The quaternion instantaneous power of the three phases (scalar part of v·i*) averages to ΣP: 7389 W and 7673 W.",
        "<i>main_test.go</i>: GA against IEEE 1459, largest pair 1–5, quaternion mean = ΣP, and a short GA run that must "
        "improve PF by at least 0.01. The program exits non-zero if the evolved filter does not improve the power factor.",
    ])
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.h2("Checks of the GA power code (doc/M9f.md)")
    r.table([
        ["Check", "Result"],
        ["GA vs IEEE 1459 on random band-limited waveforms (9 harmonics, DC offset): P, S, S<sub>1</sub>, D<sub>I</sub>, "
         "D<sub>V</sub>, S<sub>H</sub>, P<sub>1</sub>, Q<sub>1</sub>; S² = P² + Q² + D²; S² = S<sub>1</sub>² + D<sub>I</sub>² + "
         "D<sub>V</sub>² + S<sub>H</sub>²; Σ pair² = D²", "50 cases, 10<super>−9</super>"],
        ["R–L from the engine's transient: P = I²R, Q<sub>1</sub> = I²ωL, no distortion", "10<super>−4</super>"],
        ["Single-phase bridge rectifier from the engine: GA and IEEE 1459 agree (S to 1 %: harmonics above the 32nd); "
         "largest distortion pair V<sub>1</sub> × I<sub>3</sub>", "PF 0.56, THD<sub>I</sub> 140 %"],
        ["Quaternion power s = v·i*: scalar Σv·i, vector −v×i; p invariant and vector co-rotating under Clarke and Park",
         "10<super>−13</super>, 10<super>−10</super>"],
        ["Clarke and Park rotors = matrices", "500 random angles and vectors, 10<super>−14</super>"],
    ], [0.72, 0.28])
    n1 = len(r.s)
    r.h2("Exported circuits")
    r.p("<i>export_test.go</i> writes the unfiltered rectifier and one with 1 mH, 270 µF / 150 µF branches (300 Hz, "
        "420 Hz) as netlists, KiCad schematics and results sheets in <i>export/</i>. Kron's torn solution is shown there "
        "for linear circuits only; this circuit contains diodes, so the sheet reports the MNA operating point alone "
        "(the torn Newton solve is checked in kron's tests).")
    r.s[n1:] = [KeepTogether(r.s[n1:])]

    r.h1("Limits of the model")
    r.bullets([
        "Harmonics are resolved only up to the 25th, so commutation notches above it are left out of the GA power "
        "(the 0.3 % gap in S).",
        "IEEE 1459 three-phase effective quantities (S<sub>e</sub>, positive-sequence P<sub>1</sub><super>+</super>, "
        "unbalance power) and GA power for three-phase sets are not done (doc/M9f.md, \"Not done here\").",
        "The supply, rectifier and branches are ideal apart from the stated resistances; the fitness is PF only.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/power_filter    # about 2 min; -gens 10 -pop 24 -seed 1 are the defaults\n"
           "go test ./examples/power_filter   # GA and IEEE 1459 agree; a short GA improves PF (-short skips the GA)\n"
           "KRONSPICE_EXPORT=1 go test ./examples/power_filter -run TestExport")
    r.small("Sources: IEEE Std 1459-2010; M. Castro-Núñez, R. Castro-Puche, \"The IEEE Standard 1459, the CPC power theory, "
            "and geometric algebra in circuits with nonsinusoidal sources and linear loads,\" <i>IEEE Trans. Circuits Syst. I</i>, "
            "2012; doc/M9f.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "power_filter_report.pdf"))
