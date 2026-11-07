# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H16 report: the Amador & Eckardt 2012 parametric LCR circuit with the pump on the ledger, as a PDF.

    python examples/parametric_lcr/viz/make_report.py examples/parametric_lcr/parametric_lcr_report.pdf

Numbers come from examples/parametric_lcr/README.md (the example has no
run_output.txt and no figures) and export/parametric_lcr_fig4_constant_c_results.md
(MNA against Kron's torn solve for the constant-C case). Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/parametric_lcr", "Parametric LCR circuit with C(t)",
               "Does a pumped capacitor draw energy \"from the time domain\", or from the pump? — register entry H16",
               "register entry H16")

    r.h1("Summary")
    r.p("F. D. Amador and H. Eckardt (\"Simulation of a Parametric Resonance Circuit\", A.I.A.S./UPITEC, 30 March 2012) "
        "drive a series R, L, C(t) loop and report a coefficient of performance (COP) of about 200 against the source's "
        "input, read as energy \"from the time domain\" (spin connection resonance). This example reproduces every case "
        "in the paper with kronspice's time-varying capacitor, and puts on the ledger the one agent the paper's energy "
        "sum leaves out: the pump that varies C. The paper's reading is hypothesis H16; the standard prediction is "
        "parametric excitation (Mandelstam &amp; Papalexi 1932; Kouřil &amp; Vrba §10.2).")
    r.bullets([
        "<b>Every case is reproduced</b>, and in every one the source plus the pump pays for the heat in R and the "
        "stored energy; the books close to 0.001–0.4 % where resolved.",
        "<b>Tuned pumping at 2ω<sub>0</sub> grows without limit</b> (87.3 A at 30 ms, σ = 245 1/s, as Kouřil &amp; Vrba "
        "predict), paid by the pump: 256 J in against −0.085 J from the source.",
        "<b>Detuned to 0.97 f<sub>0</sub> at φ = −90°</b>, the case the paper highlights: bounded 16.1 A and 2.43 kW; "
        "the source gives 6.75 J, the pump 1139 J, R takes 1143 J.",
        "<b>The paper's COP divides heat in R by the source energy alone:</b> 169 for that case. Divided by source plus "
        "pump it is 0.998 in the bounded cases and below 1 everywhere.",
        "<b>Pumping at the drive frequency</b> (§3.2) grows only for a<sub>0</sub> of about 0.7 or more, as the second Mathieu tongue "
        "requires.",
        "<b>Bench test:</b> measure the power that drives C(t) — the motor shaft or the varactor bias supply.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis H16</b> (Amador &amp; Eckardt 2012): a parametric LCR circuit with C(t) draws energy \"from the "
        "time domain\", with a COP of about 200 by source input. <b>Standard prediction beside it:</b> a capacitor whose "
        "value is varied is driven by an outside agent (the motor that moves the plates, or a varactor's bias supply), "
        "and that agent does work on the circuit, P<sub>pump</sub> = −ṡ·½·C·v². Parametric excitation (Mandelstam and "
        "Papalexi 1932, the paper's ref. [13]; Kouřil &amp; Vrba §10.2.2) says the circuit self-excites when the "
        "elastance modulation depth m exceeds a threshold, and that all its energy comes from the source and the pump.")

    r.h1("The circuit")
    r.p("Paper Fig. 2, eqs. 19–27: a sine source U = U<sub>0</sub>·sin(ωt), inductance L, resistance R and a capacitance "
        "C<sub>eff</sub>(t) = C·(1 + a<sub>0</sub>·sin(ω<sub>1</sub>·t + φ)), with q = C<sub>eff</sub>·u<sub>C</sub> "
        "(eq. 25). kronspice's VaryingCapacitor, q = s(t)·C·v, is eq. (25) exactly. Every case is checked as source + "
        "pump = heat in R + change of stored energy.")
    r.table([
        ["parameter (paper eqs. 28–36)", "value"],
        ["U<sub>0</sub>", "5 V"],
        ["C, L, R", "230 nF, 20 mH, 20 Ω"],
        ["resonance f<sub>0</sub>", "2346.6 Hz"],
        ["quality factor Q", "14.74"],
        ["damping α = R/2L", "500 1/s"],
        ["pump", "ω<sub>1</sub> = 2ω<sub>0</sub>, a<sub>0</sub> = 0.2, φ = 0"],
    ], [0.45, 0.55])
    r.h2("Standard theory beside the paper")
    r.bullets([
        "a<sub>0</sub> = 0.2 gives an elastance depth m = 0.2020 against the tuned self-excitation threshold "
        "m<sub>c</sub> = 2/Q = 0.1356 (K&amp;V 10.101), so pumping at 2ω<sub>0</sub> self-excites.",
        "C(1 + a sin) raises the mean elastance by 1/√(1 − a²), so the circuit's own frequency is 1.0103 f<sub>0</sub>.",
        "Pumped at 2 × 0.97 f<sub>0</sub> the detuning is 4.2 % and the threshold rises to 0.2144 &gt; m: below "
        "self-excitation, the circuit is a parametric amplifier of the drive (K&amp;V 10.103).",
    ])

    r.h1("Results")
    r.p("\"Books\" is the closure of source + pump against heat in R + stored energy.")
    r.table([
        ["paper case", "kronspice", "books"],
        ["Fig. 4, constant C at f<sub>0</sub>", "peak 0.25 A = U<sub>0</sub>/R; source 0.0175 J → R 0.0169 J + stored "
         "0.0006 J", "0.001 %"],
        ["§3.1 Fig. 5, C(1 + 0.2 sin 2ω<sub>0</sub>t) at f<sub>0</sub>, 30 ms", "87.3 A and growing without limit, "
         "σ = 245 1/s as K&amp;V predict; source −0.085 J, pump 256 J → R 165 J + stored 91 J", "0.07 %"],
        ["§3.1 Figs. 6–7, drive and pump at 0.97 f<sub>0</sub>, φ = 0", "bounded 7.34 A, 504 W; pump 241.5 J, source "
         "−3.8 J", "0.04 %"],
        ["same, φ = −90°", "bounded 16.1 A, 2.43 kW; source 6.75 J, pump 1139 J → R 1143 J", "0.04 %"],
        ["same, φ = +90°", "bounded 5.72 A, 307 W; pump 138 J", "0.04 %"],
        ["§3.2 Fig. 8, C(1 + 0.2 sin<super>7</super> ω<sub>1</sub>t), ω<sub>1</sub> = drive = 1.07 ω<sub>0</sub>",
         "0.123 A, no growth", "0.4 %"],
        ["§3.2 Fig. 8, a<sub>0</sub> = 0.8", "207 A; pump 607 J → R 262 J + stored 345 J", "0.3 %"],
        ["§3.2 Fig. 9, rectangular C 0.8C ↔ 1.2C at 0.95 f<sub>0</sub>, drive off after 2 periods", "decays, 1.1 mA",
         "0.4 %"],
        ["§3.2 Fig. 9, rectangular C 0.2C ↔ 1.8C", "7.6 kA from a 2-period drive, or 443 A from a 1 mA start with no "
         "drive; pump 2360 J → R 614 J + stored 1854 J", "few % (integration error: 0.97 % → 0.32 % from 2000 to 8000 "
         "steps per period)"],
    ], [0.33, 0.47, 0.2])
    r.h2("The COP, with and without the pump")
    r.bullets([
        "The paper's COP divides the heat in R by the source energy alone: 0.964 for Fig. 4, 169 for the φ = −90° case "
        "the paper highlights, and negative where the source absorbs.",
        "Divided by source plus pump it is 0.998 in the bounded cases and below 1 everywhere.",
        "In every case the heat in R and the stored energy are paid by the source plus the pump that varies C. The "
        "pump's share is what the paper's \"energy from spacetime\" counts.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "The pump's power P<sub>pump</sub> = −ṡ·½·C·v² is on the ledger in every run; closure 0.001 % to 0.4 %, and a "
        "few % only for the 0.2C ↔ 1.8C rectangular case, where it falls with the step as integration error does.",
        "Tests (main_test.go): classical amplitude (TestClassical), tuned growth paid by the pump (TestTunedGrowth), "
        "detuned bounded case (TestDetunedBounded), second tongue a<sub>0</sub> = 0.2 against 0.8 (TestSecondTongue).",
        "Export sheet, Fig. 4 constant-C case at 2346.61 Hz: 0.25 A through R (5 V across), MNA and Kron's torn solve "
        "agree to 1.3×10<super>−14</super> of the largest node voltage, and the power taken by the R, L, C branches is "
        "the same torn and untorn (Roth 1959 §5). The pumped-C circuit is nonlinear, so its sheet gives MNA only.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The capacitor variation s(t) is imposed; the pump is ideal, and only its power on the circuit is booked, not "
        "the losses of a real motor or bias supply.",
        "Lumped series R, L, C; no losses beyond R.",
    ])

    r.h1("Status")
    r.p("<b>H16: untested</b> (reproduced; the COP omits the pump). All of the paper's cases are reproduced: "
        "U<sub>0</sub>/R = 0.25 A classical case; unbounded growth at 2ω<sub>0</sub> pumping (σ = 245 1/s); bounded "
        "≈ 16 A, 2.4 kW at 0.97 f<sub>0</sub>, φ = −90° (source-only ratio 169); pumping at ω<sub>0</sub> grows only "
        "for a<sub>0</sub> of about 0.7 or more. In every case the source plus the pump pays for R and the stored energy. The bench "
        "test is to measure the power that drives C(t), at the motor shaft or the varactor bias supply.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/parametric_lcr\n"
           "go test ./examples/parametric_lcr\n"
           "KRONSPICE_EXPORT=1 go test ./examples/parametric_lcr -run TestExport   # export/ sheets")
    r.small("Related: the same parametric pumping in Kron's machine frame, examples/kron_machine and "
            "examples/parametric_pump.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "parametric_lcr_report.pdf"))
