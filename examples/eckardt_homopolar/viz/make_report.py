# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H15 report: Eckardt's self-excited twin-disk homopolar generator, as a PDF.

python examples/eckardt_homopolar/viz/make_report.py examples/eckardt_homopolar/eckardt_homopolar_report.pdf

Numbers come from examples/eckardt_homopolar/README.md (the example has no
run_output.txt and no figures). Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1  # keep headings with the text that follows
kronreport.H2.keepWithNext = 1


def build(out):
    r = Report("examples/eckardt_homopolar", "Eckardt's self-excited homopolar generator",
               "Does the twin-disk machine's “giant resonance” draw energy from spacetime, or from the shaft? — "
               "register entry H15, CLAUDE.md §3f", "register entry H15")

    r.h1("Summary")
    r.p("H. Eckardt (A.I.A.S./UPITEC, 2012) describes a twin-disk homopolar machine that excites its own field and "
        "reports an energy surplus at a “giant resonance” where the effective resistance vanishes. This example runs "
        "the paper's equations as printed, then the same machine as a kronspice circuit built on the N machine "
        "element, and compares their energy books.")
    r.bullets([
        "<b>The resonance is reproduced:</b> eq. (9) gives ω<sub>res</sub> = 29.84 rad/s (paper: about 30). It is the "
        "self-excitation threshold of a series-wound dynamo.",
        "<b>The paper's equations run away:</b> the current passes 10<super>6</super> A at 0.463 s, and the energy "
        "into R and the coil is 1.42×10<super>9</super> J against 10.7 J of drive work.",
        "<b>The surplus is the omitted back torque:</b> the missing 1.423×10<super>9</super> J equals, to "
        "10<super>−3</super>, the work the back torque ½r²(B<sub>m</sub> + B<sub>I</sub>)·I would have taken from "
        "the shaft, which eq. (12) leaves out.",
        "<b>With the back torque the current stays bounded:</b> in kronspice the speed settles at ω<sub>res</sub> with "
        "I ≈ 9.949 A and P = α·ω<sub>res</sub> = 29.84 W from the drive; the ledger closes to 3×10<super>−5</super> J.",
        "<b>The paper's flywheel pulse is the flywheel's energy:</b> what R receives is what the shaft gave (books "
        "close to 10<super>−8</super> J).",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H15).</b> ECE theory (Evans; A.I.A.S. UFT paper 107) holds that homopolar resonances draw "
        "energy from spacetime. Eckardt's 2012 paper, “Simulation of resonances in the homopolar generator”, models a "
        "self-excited twin-disk machine and reports a “giant resonance” at R<sub>eff</sub> = 0 with an energy surplus.")
    r.p("<b>Standard prediction.</b> R<sub>eff</sub> = 0 is the self-excitation threshold of a series-wound dynamo "
        "(Siemens 1867; Bullard 1955; Moffatt 1978, ch. 6). Above it the current grows by taking power from the "
        "drive. The Lorentz force that gives the disk EMF ½ωr²B also gives a back torque ½r²B·I per disk, so "
        "ω·(torque) = EMF·I exactly and the shaft pays for every joule.")

    r.h1("The machine and the paper's equations")
    r.p("Two disks of radius r on one shaft are in series with a resistance R and a coil of N turns, length l and "
        "core permeability μ<sub>r</sub> (paper Fig. 1). The “static” disk sits in a permanent-magnet field "
        "B<sub>m</sub>; the “dynamic” disk sits in the coil's own field B<sub>I</sub> = μ<sub>0</sub>μ<sub>r</sub>NI/l, "
        "so the loop current I excites it. ω is the shaft speed, J the inertia and α the drive torque.")
    r.code("L·dI/dt + R·I = ½·ω·r²·(Bm + B_I)          (6)\n"
           "R_eff = R − ½·ω·r²·μ0·μr·N/l               (8)\n"
           "ω_res = 2·R·l/(r²·μ0·μr·N)                 (9)\n"
           "J·dω/dt = α   (constant drive torque)       (12)  as printed\n"
           "J·dω/dt = α − ½·r²·(Bm + B_I)·I                   with the back torque")
    r.table([
        ["parameter (paper p. 4)", "value", "parameter", "value"],
        ["R", "0.3 Ω", "r", "10 cm"],
        ["L", "1 mH", "μ<sub>r</sub>", "100"],
        ["N", "800", "B<sub>m</sub>", "0.1 T"],
        ["l", "5 cm", "J", "0.01 kg·m²"],
        ["α", "1 N·m", "ω<sub>res</sub> from eq. (9)", "29.84 rad/s (paper: about 30)"],
    ], [0.25, 0.2, 0.25, 0.3])
    r.p("In kronspice the machine is one N element with four windings (static disk path, dynamic disk path, coil, "
        "and a magnet equivalent at 1 A). Its speed-matrix entries are K1 = ½r²B<sub>m</sub> = 0.0005 V·s and "
        "K2 = ½r²μ<sub>0</sub>μ<sub>r</sub>N/l = 0.01005 V·s/A; the shaft is a mechanical node with inertia J and "
        "the drive torque α. The back torque is part of the element (Kron's speed-voltage tensor), not an added term.")

    r.h1("Results")
    r.h2("1. The paper's equations as printed (RK4, no back torque)")
    r.p("The speed rises linearly, and once ω passes ω<sub>res</sub> the current runs away: 0.56 A at 0.30 s, 211 A "
        "at 0.40 s, 1.2×10<super>5</super> A at 0.45 s, over 10<super>6</super> A at 0.463 s. The energy into R and "
        "the coil is 1.42×10<super>9</super> J against 10.7 J of drive work, all of it spin-up. The missing "
        "1.423×10<super>9</super> J equals, to 10<super>−3</super>, the work the omitted back torque "
        "½r²(B<sub>m</sub> + B<sub>I</sub>)·I would have taken from the shaft. The “energy from spacetime” in these "
        "equations is the shaft work they leave out.")
    r.h2("2. The same machine in kronspice (back torque inside the element), L = 1 mH, 3 s")
    r.table([
        ["t", "ω", "I", "B<sub>I</sub>"],
        ["0.30 s", "30.08 rad/s", "0.579 A", "1.16 T"],
        ["0.60 s", "25.77 rad/s", "5.05 A", "10.2 T"],
        ["1.20 s", "27.37 rad/s", "12.69 A", "25.5 T"],
        ["3.00 s", "30.33 rad/s", "9.68 A", "19.5 T"],
    ], [0.25, 0.25, 0.25, 0.25])
    r.p("The current stays bounded and the speed settles at ω<sub>res</sub>. Energy ledger: drive work 84.75 J = "
        "R 80.1 J + coil 0.047 J + spin-up 4.6 J + reference resistors 0.0002 J; it closes to 3×10<super>−5</super> J. "
        "Steady state: α = K1·I + K2·I² gives I = 9.949 A and P = α·ω<sub>res</sub> = 29.84 W = R·I² (29.69 W).")
    r.p("That steady state needs B<sub>I</sub> = 20 T in the iron, which saturates near 1.5–2 T, so a real coil limits "
        "the excitation well before. With the coil inductance the stated geometry implies (50.5 H if B<sub>I</sub> "
        "fills the disk area) the run is much slower and the shaft reverses as the field hunts; the books still close.")
    r.h2("3. The paper's second run: disk spun to 40 rad/s, then connected (Figs. 6–9)")
    r.table([
        ["drive", "peak current", "speed after 0.5 s", "energy into R"],
        ["α = 0", "32.34 A at 0.057 s", "19.53 rad/s", "6.092 J = flywheel 6.092 J"],
        ["α = 1 N·m", "46.32 A at 0.049 s", "25.50 rad/s",
         "19.197 J = flywheel 4.749 J + drive 14.466 J − coil 0.017 J"],
    ], [0.13, 0.2, 0.17, 0.5])
    r.p("The pulse is the flywheel's energy discharged through the self-excited dynamo: what R receives is what the "
        "shaft gave (books close to 10<super>−8</super> J).")

    r.h1("Ledger and validation")
    r.bullets([
        "ω<sub>res</sub> from eq. (9) against the paper's value.",
        "The energy missing from the paper's books equals the work of the omitted back torque, to 10<super>−3</super>.",
        "With the back torque the current is bounded and the ledger closes (3×10<super>−5</super> J in run 2; "
        "10<super>−8</super> J in the flywheel runs).",
        "<i>main_test.go</i> covers all of these; <i>ECKARDT_DUMP=1</i> prints the comparison series. "
        "<i>export_test.go</i> writes the driven and flywheel netlists, KiCad schematics and results sheets to "
        "<i>export/</i>.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The coil field is linear (B<sub>I</sub> = μ<sub>0</sub>μ<sub>r</sub>NI/l); the 20 T steady state is far "
        "above iron saturation (1.5–2 T), so a real machine would limit the excitation much earlier.",
        "The main runs use the paper's L = 1 mH, not the 50.5 H its geometry implies.",
        "The two ECE-only resonances of UFT 107 are not modelled, because no equations are supplied for them.",
    ])

    r.h1("Status")
    r.p("Register status: <b>untested</b>. The classical resonance (R<sub>eff</sub> = 0) is the series-dynamo "
        "self-excitation threshold and is reproduced. The paper's energy surplus equals, to 10<super>−3</super>, the "
        "work of the back torque its eq. (12) omits; with it the current stays bounded (≈ 9.9 A at ω<sub>res</sub>, "
        "30 W from the drive) and the ledger closes. The two ECE-only resonances are not modelled.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/eckardt_homopolar\n"
           "go test ./examples/eckardt_homopolar\n"
           "python examples/eckardt_homopolar/viz/make_report.py examples/eckardt_homopolar/eckardt_homopolar_report.pdf")
    r.small("Sources: H. Eckardt, “Simulation of resonances in the homopolar generator,” A.I.A.S./UPITEC, 13 February "
            "2012 (hypothesis source); E. C. Bullard, Proc. Camb. Phil. Soc. 51, 744 (1955); H. K. Moffatt, Magnetic "
            "Field Generation in Electrically Conducting Fluids, Cambridge 1978, ch. 6; G. Kron, Equivalent Circuits "
            "of Electric Machinery, Wiley 1951, Epilogue.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "eckardt_homopolar_report.pdf"))
