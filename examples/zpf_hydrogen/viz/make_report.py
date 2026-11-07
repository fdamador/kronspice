# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/zpf_hydrogen: Puthoff 1987, hydrogen's ground state and the zero-point field (H24).

Run from the repository root:

    python examples/zpf_hydrogen/viz/make_report.py examples/zpf_hydrogen/zpf_hydrogen_report.pdf

Numbers come from examples/zpf_hydrogen/run_output.txt (the program's output;
where the README rounds or differs, run_output.txt is used) and README.md
(construction, declared choices, sources); the figure is viz/zpf_hydrogen.png.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/zpf_hydrogen", "Hydrogen held up by the zero-point field",
               "Does the orbiting electron absorb from the zero-point field exactly what it radiates, and does that fix "
               "hydrogen's ground state? (Puthoff 1987) — register entry H24",
               "register entry H24", "CODATA 2018; 16 random-phase runs")

    r.h1("Summary")
    r.p("H. E. Puthoff (Phys. Rev. D 35, 3266, 1987) works in <b>stochastic electrodynamics</b> (SED), which treats empty "
        "space as filled with a real, random classical electromagnetic field carrying ½ħω of energy per mode, the "
        "zero-point field. A classical electron circling a proton radiates and should spiral in within about 16 ps. "
        "Puthoff shows that on the orbit with m·ω·r² = ħ (m the electron mass, ω the orbit's angular frequency, r its "
        "radius) the electron absorbs from the zero-point field exactly what it radiates, so it stays. This example "
        "recomputes the paper's balance with CODATA values, checks whether the balance is restoring, simulates the "
        "paper's oscillator in a random zero-point field, and builds the same oscillator as a kronspice circuit driven "
        "by its resistor's zero-point noise.")
    r.bullets([
        "<b>The balance gives Bohr's numbers:</b> r = 5.291772×10<super>−11</super> m (the Bohr radius a<sub>0</sub>) and "
        "−13.6057 eV, with 4.669288×10<super>−8</super> W absorbed and the same radiated. This is expected: the paper's "
        "eq. 20 is Bohr's condition.",
        "<b>The averaged balance is restoring:</b> absorbed/radiated is 1.414 at 0.5 a<sub>0</sub> and 0.707 at "
        "2 a<sub>0</sub>, pulling the orbit back to a<sub>0</sub> with a time constant of 9.34×10<super>−11</super> s "
        "(≈ 614,349 orbits).",
        "<b>The oscillator in a simulated zero-point field holds ½ħω<sub>0</sub>:</b> mean energy 0.5037 ± 0.0009 "
        "ħω<sub>0</sub> (expected 0.5036); absorbed 0.005076 and radiated 0.005077 agree to 0.03 %.",
        "<b>The kronspice tank in zero-point resistor noise</b> stores 0.4916 ± 0.0035 ħω<sub>0</sub> (expected 0.4935); "
        "its ledger closes (noise source 645.0 J = resistor 639.3 J + stored 4.59 J, mismatch 0.002).",
        "<b>Zero net energy exchange:</b> what the zero-point field puts in comes back out. Nothing here is an energy source.",
        "<b>Status:</b> reproduced at Bohr-theory level, as the paper claims; beyond it untested here.",
    ])

    r.h1("The hypothesis and the standard prediction")
    r.p("<b>H24 (hypothesis):</b> hydrogen's ground state is a zero-point-fluctuation-determined state: in SED the "
        "orbiting electron absorbs from the random zero-point field exactly what it radiates when m·ω·r² = ħ.")
    r.p("<b>Standard prediction:</b> the quantum-mechanical ground state has zero orbital angular momentum and a spread "
        "of radii (most probable a<sub>0</sub>, mean 1.5 a<sub>0</sub>); a single circular orbit with angular momentum ħ "
        "is not it. Puthoff's result is at the level of Bohr theory, as he says; extending SED to the Schrödinger ground "
        "state was \"yet to be achieved\" in 1987 (his footnote 17).")

    r.h1("Results")
    r.h2("1. The paper's numbers with CODATA 2018 values")
    r.table([
        ["Quantity", "Result", "Compare"],
        ["Orbit radius from m·ω·r² = ħ and Coulomb (eq. 20)", "5.291772×10<super>−11</super> m",
         "Bohr radius a<sub>0</sub> = 5.291772×10<super>−11</super> m"],
        ["Energy", "−13.6057 eV", "Bohr and Schrödinger, infinitely heavy proton: −13.6057 eV (measured ionization "
         "energy 13.598 eV includes the proton's recoil)"],
        ["Power absorbed from the zero-point field (eq. 18)", "4.669288×10<super>−8</super> W",
         "radiated (eq. 19): 4.669288×10<super>−8</super> W"],
        ["Orbit frequency", "6.5797×10<super>15</super> Hz", ""],
        ["Collapse time with no zero-point field", "1.56×10<super>−11</super> s", "the classical problem Puthoff answers"],
    ], [0.36, 0.22, 0.42])

    r.h2("2. Does the balance hold the orbit?")
    r.p("Net power P<sub>abs</sub> − P<sub>rad</sub> on Coulomb orbits of other radii, averaged as in the paper:")
    r.table([
        ["Radius", "Absorbed ÷ radiated", "Net power"],
        ["0.5 a<sub>0</sub>", "1.414", "+3.095×10<super>−7</super> W (gains energy, grows)"],
        ["0.7 a<sub>0</sub>", "1.195", "+3.797×10<super>−8</super> W"],
        ["0.9 a<sub>0</sub>", "1.054", "+3.850×10<super>−9</super> W"],
        ["1.0 a<sub>0</sub>", "1.000", "−2.847×10<super>−17</super> W"],
        ["1.1 a<sub>0</sub>", "0.953", "−1.484×10<super>−9</super> W"],
        ["1.5 a<sub>0</sub>", "0.816", "−1.693×10<super>−9</super> W"],
        ["2.0 a<sub>0</sub>", "0.707", "−8.548×10<super>−10</super> W (loses energy, shrinks)"],
    ], [0.2, 0.3, 0.5], bold_rows=(4,))
    r.p("The ratio absorbed/radiated is ħ/(mωr²), which on Coulomb orbits falls as 1/√r (0.707 at 2 a<sub>0</sub>, "
        "1.414 at 0.5 a<sub>0</sub>). The averaged balance pulls the orbit back to a<sub>0</sub> with a time constant of "
        "9.34×10<super>−11</super> s (≈ 614,349 orbits).")

    r.h2("3. The paper's oscillator in a simulated zero-point field")
    r.p("Puthoff builds the orbit from two one-dimensional oscillators (eqs. 4–17). Here that oscillator is simulated "
        "directly: a random field with the zero-point spectrum is made from random phases (eq. 2) and a Fourier "
        "transform, and the radiating oscillator is integrated in it (own Runge–Kutta), 16 independent runs.")
    r.table([
        ["", "Simulated", "Expected for this field", "Formula"],
        ["Mean energy", "0.5037 ± 0.0009 ħω<sub>0</sub>", "0.5036", "½ħω<sub>0</sub> = 0.5000"],
        ["Power absorbed", "0.005076", "0.005076", "eq. 16: ħΓω<sub>0</sub>³/2 = 0.005000"],
        ["Power radiated", "0.005077", "", ""],
    ], [0.2, 0.3, 0.22, 0.28])
    r.p("In steady state the oscillator radiates what it absorbs (difference 0.03 %) and holds ½ħω<sub>0</sub>, the known "
        "SED result for a harmonic oscillator (Boyer 1975). Declared choices:")
    r.bullets([
        "<b>Units:</b> scaled, m = ω<sub>0</sub> = ħ = 1.",
        "<b>Damping:</b> Γω<sub>0</sub> = 0.01 instead of the electron's 2.6×10<super>−7</super> at the hydrogen orbit "
        "frequency, so a run spans 524 relaxation times. Every result is compared with the formula for that damping.",
        "<b>Field cut:</b> the field stops at 2ω<sub>0</sub>. The full zero-point spectrum makes a classical particle's "
        "kinetic energy diverge.",
        "<b>Radiation reaction:</b> used in its reduced (Landau–Lifshitz) form, which is equal at resonance.",
    ])
    r.figure(os.path.join(HERE, "zpf_hydrogen.png"),
             "Figure 1. Left: power absorbed from the zero-point field (eq. 18, blue) and radiated (Larmor, eq. 19, red "
             "dashed) on Coulomb orbits from 0.3 to 3 a<sub>0</sub>; they cross at the Bohr radius, with absorption "
             "winning inside and radiation outside. Right: the oscillator's energy in one simulated zero-point-field "
             "run (scaled time 1/ω<sub>0</sub>); its time average, 0.503 ħω<sub>0</sub>, lies on ½ħω<sub>0</sub>.",
             max_height_mm=90)

    r.h2("4. The same oscillator as a kronspice circuit")
    r.p("A series tank of L = 1 H, C = 1 F and R = 0.01 Ω is driven by the resistor's own zero-point noise. By Callen "
        "and Welton's fluctuation–dissipation theorem (1951, cited in Puthoff's footnote 13), a resistor at absolute zero "
        "has a one-sided voltage noise S<sub>V</sub> = ħωR/π per unit angular frequency. ħ is scaled to 1 J·s; the "
        "circuit is linear, so energies scale with ħ.")
    r.table([
        ["", "kronspice", "Expected"],
        ["Energy stored in the tank", "0.4916 ± 0.0035 ħω<sub>0</sub>", "0.4935"],
        ["Ledger", "noise source 645.0 J = resistor 639.3 J + stored 4.59 J", "mismatch 0.002"],
    ], [0.3, 0.45, 0.25])
    r.bullets([
        "The expected value allows for kronspice joining the noise waveform's points by straight lines.",
        "A kronspice lesson from this run: after every corner of a stepwise (PWL) waveform kronspice takes a short "
        "backward-Euler step (standard SPICE practice). With 65,536 corners that step damped the tank about as much as "
        "its resistor did, and the energy came out 7 % low. A maximum step of 0.02 s makes the extra damping negligible.",
        "The noise source puts in what the resistor takes out: net, the zero-point field supplies nothing, as in "
        "Puthoff's balance.",
    ])

    r.h1("Where this stands")
    r.bullets([
        "<b>Bohr-theory level.</b> Eq. 20 is Bohr's condition, so energy and size agree with Bohr and Schrödinger. The "
        "test of SED is everything beyond that.",
        "<b>Later full SED simulations:</b> Cole &amp; Zou (Phys. Lett. A, 2003) report the Schrödinger ground-state "
        "probability density; Nieuwenhuizen &amp; Liska (Phys. Scr. T165, 014006, 2015), with longer runs, find that the "
        "atom eventually ionizes itself.",
        "<b>Energy.</b> The zero-point exchange is zero net: absorbed equals radiated.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "Paper's balance: absorbed and radiated power equal at a<sub>0</sub> (net −2.847×10<super>−17</super> W).",
        "Oscillator: absorbed and radiated agree to 0.03 %; mean energy matches the expectation for the cut field.",
        "kronspice tank: ledger noise source = resistor + stored to a mismatch of 0.002.",
        "Tests (<i>main_test.go</i>): Bohr numbers, stability ratio, oscillator energy and balance, tank energy and ledger, FFT.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "A single circular orbit, averaged as in the paper; no orbital spread and no full three-dimensional SED atom.",
        "The oscillator uses a damping 0.01 instead of 2.6×10<super>−7</super> and a field cut at 2ω<sub>0</sub> "
        "(declared).",
    ])

    r.h1("Status")
    r.p("H24 is <b>reproduced at Bohr-theory level, as the paper claims; beyond it untested here.</b> The quantum ground "
        "state has L = 0 and a spread of radii; later full SED simulations report the Schrödinger density (Cole &amp; Zou) "
        "and, with longer runs, self-ionization (Nieuwenhuizen &amp; Liska). The zero-point exchange is zero net, not an energy source.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/zpf_hydrogen -viz examples/zpf_hydrogen/viz   # about 15 seconds\n"
           "python3 examples/zpf_hydrogen/viz/render.py\n"
           "go test ./examples/zpf_hydrogen")
    r.small("Sources: H. E. Puthoff, Phys. Rev. D 35, 3266 (1987); T. H. Boyer, Phys. Rev. D 11, 790 (1975); "
            "H. B. Callen, T. A. Welton, Phys. Rev. 83, 34 (1951); D. C. Cole, Y. Zou, Phys. Lett. A (2003), "
            "doi:10.1016/j.physleta.2003.08.022; T. M. Nieuwenhuizen, M. T. P. Liska, Phys. Scr. T165, 014006 (2015); "
            "CODATA 2018.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "zpf_hydrogen_report.pdf"))
