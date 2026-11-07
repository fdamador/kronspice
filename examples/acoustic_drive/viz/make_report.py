# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/acoustic_drive (piezo ring driving a liquid sphere).

python examples/acoustic_drive/viz/make_report.py examples/acoustic_drive/acoustic_drive_report.pdf

Numbers come from examples/acoustic_drive/README.md (results table) and the
results sheets in examples/acoustic_drive/export/*_results.md (torn solve).
The example has no run_output.txt and no figures.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/acoustic_drive", "Acoustic drive of a sonoluminescence flask",
               "Electrical → mechanical → acoustic in one circuit, in H<sub>2</sub>O and D<sub>2</sub>O "
               "— CLAUDE.md §3.1, §3h", "CLAUDE.md §3.1, §3h",
               "AC sweep, 2001 points over ±3 % of f<sub>0</sub>")

    r.h1("Summary")
    r.p("One netlist joins three energy domains of CLAUDE.md §3.1: a 50 Ω generator, an optional series inductor "
        "that tunes out the transducer's capacitance, a piezoelectric ring, its coupling into the breathing mode of a "
        "3 cm liquid sphere, and that mode itself, in light water (H<sub>2</sub>O) and heavy water (D<sub>2</sub>O). "
        "The question is how much drive and power are needed for a 1.3 atm pressure amplitude at the centre of the "
        "flask, the range where single-bubble sonoluminescence is observed, and where every watt goes.")
    r.bullets([
        "<b>Power into the liquid is fixed by the mode:</b> 6.96 mW in H<sub>2</sub>O and 6.72 mW in D<sub>2</sub>O at "
        "1.3 atm, equal to P²/(2·R<sub>par</sub>), with or without matching.",
        "<b>Matching raises the pressure per volt about 21 times:</b> from 3.22 to 69.4 kPa/V in H<sub>2</sub>O, so the "
        "generator EMF falls from 40.9 V to 1.90 V.",
        "<b>Matching moves the loss, it does not remove it:</b> without the inductor the transducer's loss takes about a "
        "third of the input (3.50 of 10.9 mW); with it, most of the input (19.0 of 26.2 mW) is burned in the "
        "generator's own 50 Ω.",
        "<b>D<sub>2</sub>O resonates 6.5 % lower</b> than H<sub>2</sub>O (23.40 against 25.03 kHz; sound speed 1404 "
        "against 1502 m/s) and takes slightly less power.",
        "<b>The ledger closes in every run:</b> source power in equals heat out, residual 0 of the scale; the tests "
        "require below 10<super>−9</super>.",
    ])

    r.h1("What is being tested")
    r.p("This is a standard-physics construction, not a register hypothesis. It is the first circuit that joins the "
        "electrical, mechanical and acoustic channels through Kron's network analogies (§3a.1 items 3 and 7): each "
        "non-electrical channel has an across and a through variable whose product is power, and different domains "
        "are joined by connection tensors (Kron 1939, Ch. V p. 133). The liquids come from the fluid/ package (§3h). "
        "The bubble itself is out of scope here; it enters as a cited equivalent circuit in examples/bubble_collapse, "
        "on this chain's pressure node p.")

    r.small("Sources: Brenner, Hilgenfeldt &amp; Lohse, Rev. Mod. Phys. 74, 425 (2002); Kinsler et al., Fundamentals of "
            "Acoustics, 4th ed. (2000), §7.5 and §9.3; Firestone, J. Acoust. Soc. Am. 4, 249 (1933); Wagner &amp; Pruß, "
            "J. Phys. Chem. Ref. Data 31, 387 (2002); Herrig et al., J. Phys. Chem. Ref. Data 47, 043102 (2018); "
            "Kron, Tensor Analysis of Networks (1939), Ch. V p. 133.")
    r.h1("The model")
    r.p("Every node belongs to one domain, and ground is the common reference of all of them (zero velocity, ambient "
        "pressure). The mechanical domain uses the mobility analogy (velocity across, force through), so a mass is a "
        "capacitor to ground and Kirchhoff's current law is the force balance at a point. The acoustic domain is "
        "pressure across, volume velocity through.")
    r.table([
        ["stage", "elements", "domain"],
        ["generator", "VGEN 1 V AC, RGEN 50 Ω", "electrical"],
        ["matching", "LMATCH = 1/(ω<sub>0</sub>²·C<sub>0</sub>) with 0.5 Ω, or RWIRE 0.5 Ω alone", "electrical"],
        ["piezo ring, electrical side", "C<sub>0</sub> 4 nF clamped capacitance", "electrical"],
        ["transducer", "gyrator YPZT, force factor φ = 0.6 N/V", "electrical ↔ mechanical"],
        ["piezo ring, mechanical side", "CMASS 4 g moving mass; LSPRING 1/k with k = ω<sub>0</sub>²·m (resonance placed at "
         "the mode); RLOSS 1/b, b = 3 N·s/m", "mechanical"],
        ["coupling", "gyrator YCOUPLE, effective area S<sub>eff</sub> = 2 cm²", "mechanical ↔ acoustic"],
        ["flask breathing mode", "CMODE C<sub>a</sub>, LMODE M<sub>a</sub>, RMODE R<sub>par</sub> from fluid.SphericalResonator", "acoustic"],
    ], [0.24, 0.56, 0.20])
    r.p("A transducer is an ideal gyrator, i<sub>1</sub> = g·V<sub>2</sub>, i<sub>2</sub> = −g·V<sub>1</sub>, which passes "
        "power from one port to the other without loss: a piezo element of force factor φ has g = φ, a piston of area S "
        "has g = S. In Kron's primitive network it is two axes coupled by an antisymmetric mutual admittance.")
    r.p("The breathing mode of a liquid sphere of radius R with a pressure-release wall sits at f<sub>0</sub> = c/(2R). It "
        "is a parallel acoustic resonator whose across variable is the centre pressure: C<sub>a</sub> = 2R³/(πρc²) from "
        "the mode's potential energy, M<sub>a</sub> = 1/(ω<sub>0</sub>²·C<sub>a</sub>), R<sub>par</sub> = "
        "Q/(ω<sub>0</sub>·C<sub>a</sub>), with an assumed Q of 1500. The program sweeps 2001 points over ±3 % of "
        "f<sub>0</sub> with analysis.AC, finds the pressure peak, scales the drive to 1.3 atm at the centre and reads the ledger.")

    r.h1("Results")
    r.p("All four runs at a 1.3 atm centre amplitude (the sonoluminescence range is about 1.2–1.5 atm at 20–40 kHz, "
        "Brenner, Hilgenfeldt &amp; Lohse 2002). f<sub>0</sub> is the mode frequency; the pressure peak of the chain is in brackets.")
    r.table([
        ["run", "f<sub>0</sub> (peak)", "pressure per volt", "generator EMF", "source power", "into the liquid",
         "transducer loss", "generator's own 50 Ω"],
        ["H<sub>2</sub>O, no matching", "25.03 kHz (24.27)", "3.22 kPa/V", "40.9 V", "10.9 mW", "6.96 mW", "3.50 mW", "0.42 mW"],
        ["H<sub>2</sub>O, matched", "25.03 kHz (25.02)", "69.4 kPa/V", "1.90 V", "26.2 mW", "6.96 mW", "0.42 µW", "19.0 mW"],
        ["D<sub>2</sub>O, no matching", "23.40 kHz (22.69)", "3.20 kPa/V", "41.2 V", "10.5 mW", "6.72 mW", "3.26 mW", "0.51 mW"],
        ["D<sub>2</sub>O, matched", "23.40 kHz (23.40)", "72.2 kPa/V", "1.82 V", "23.5 mW", "6.72 mW", "0.39 µW", "16.6 mW"],
    ], [0.16, 0.14, 0.11, 0.11, 0.11, 0.11, 0.12, 0.14])
    r.bullets([
        "Without the inductor the chain's peak sits below f<sub>0</sub> (24.27 kHz in H<sub>2</sub>O) and the transducer "
        "loss takes a third of the input. With it the chain is a low impedance at resonance, most of the input is "
        "dissipated in the generator's 50 Ω, and the peak returns to f<sub>0</sub>.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "Ledger: source power in equals heat out in every run, residual 0 of the scale.",
        "main_test.go checks, for both liquids: liquid heat equals P²/(2·R<sub>par</sub>) to 10<super>−9</super>; ledger "
        "residual below 10<super>−9</super>; matching puts the peak back at f<sub>0</sub> within 0.1 % and raises the "
        "pressure per volt by more than ten times; D<sub>2</sub>O's f<sub>0</sub> is below H<sub>2</sub>O's.",
        "Kron twin (export sheets): the matched chains are torn automatically at LMATCH into two pieces; torn and MNA "
        "node voltages agree to 4.6×10<super>−12</super> (H<sub>2</sub>O) and 2.6×10<super>−12</super> (D<sub>2</sub>O) "
        "of the largest node voltage, and the power V·I* is the same torn and untorn (Roth 1959 §5). The unmatched "
        "H<sub>2</sub>O chain is torn at RGEN and agrees to 4.4×10<super>−15</super>.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The flask Q (1500), the coupling area S<sub>eff</sub> and the transducer values (C<sub>0</sub>, φ, m, b) are "
        "declared assumptions in main.go, not measurements. A hydrophone calibration of the real flask sets them; after "
        "that the chain is a design tool, and an evolvable one (matching network, transducer and drive frequency for "
        "the least power at a target pressure).",
        "The bubble dynamics (Rayleigh–Plesset or Keller–Miksis, collapse, light) stay out of the engine (CLAUDE.md §1, "
        "non-goals); they enter as a cited equivalent circuit in examples/bubble_collapse.",
        "These acoustic modes, about 25 kHz, sit ten orders of magnitude below molecular vibration frequencies "
        "(doc/future-work.md).",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/acoustic_drive\n"
           "go test ./examples/acoustic_drive ./fluid\n"
           "python examples/acoustic_drive/viz/make_report.py examples/acoustic_drive/acoustic_drive_report.pdf")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "acoustic_drive_report.pdf"))
