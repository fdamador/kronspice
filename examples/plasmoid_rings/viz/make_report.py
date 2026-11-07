# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H26 report: Nachamkin's force-free plasmoids as toroidal current rings.

python examples/plasmoid_rings/viz/make_report.py examples/plasmoid_rings/plasmoid_rings_report.pdf

Numbers come from run_output.txt (the program's report) and README.md; the
figure is viz/plasmoid_rings.png made by viz/render.py. Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
from reportlab.platypus import KeepTogether  # noqa: E402


def build(out):
    r = Report("examples/plasmoid_rings", "Plasmoids as toroidal current rings",
               "Does Nachamkin's force-free plasmoid hold together and stop radiating? — register entry H26",
               "register entry H26", "about 8 seconds; deterministic")

    r.h1("Summary")
    r.p("J. Nachamkin's 1992 Phillips Laboratory report proposes a plasmoid made of a sphere of plasma holding a "
        "standing electromagnetic wave whose electric and magnetic fields are parallel everywhere (a force-free, or "
        "Beltrami, field), a swirling plasma flow that balances the field's pressure, and a \"resonant size\" at which "
        "it should stop radiating. He offers it as a model of ball lightning, with a worked example of 1 MJ in a "
        "10 cm sphere. This example checks his field and numbers, computes the sphere's radiation lifetime, and "
        "models the two real rings his field corresponds to with kronspice: a spheromak decaying through its "
        "resistance, and an electron ring compressed by a capacitor bank.")
    r.bullets([
        "<b>The field mathematics is sound.</b> curl F = αF to 2×10<super>−10</super>, div F = 0 to 10<super>−11</super>; field lines "
        "stay on nested tori around a current ring at 0.611 of the radius (61.1 mm for R = 10 cm).",
        "<b>His density condition is ordinary plasma physics:</b> it is exactly the cold-plasma dispersion relation "
        "ω² = ω<sub>p</sub>² + α²c² (to −4×10<super>−17</super>). His Appendix B numbers are reproduced exactly.",
        "<b>His worked example breaks his own model:</b> the electrons would oscillate at 7.92×10<super>8</super> m/s = 2.64 c, "
        "in a field about 1,000 times air's breakdown strength.",
        "<b>The resonant size does not stop radiation:</b> the plasma sphere's energy lifetime is 0.6 µs (TE) and 1.2 µs (TM), "
        "and the two halves of the force-free field need different sizes (αR = 4.4933 and 2.7437).",
        "<b>A resistive spheromak ring</b> of 1 MJ in 10 cm lasts 0.598 µs, 18.9 µs and 0.598 ms at 1, 10 and 100 eV; "
        "lasting one second needs 14082 eV (163 million K).",
        "<b>A collapsing electron ring</b> driven by a 50 kJ bank goes from 16.4 cm to 2.60 cm and 2.0 to 15.0 MeV, as "
        "R proportional to B<super>−½</super> predicts; it takes 21 J of the bank's energy and the ledger closes to 10<super>−9</super>.",
    ])

    r.h1("What is being tested")
    r.h2("The claim (hypothesis H26)")
    r.p("Long-lived plasmoids (ball lightning) are force-free, time-harmonic electromagnetic fields in a plasma, "
        "stabilized by vortical flow and non-radiating at a resonant size (Nachamkin, <i>Force-Free Time-Harmonic "
        "Plasmoids</i>, PL-TR-92-3044, 1992).")
    r.h2("The standard prediction beside it")
    r.p("A standing wave in a cold plasma obeys ω² = ω<sub>p</sub>² + α²c², and a finite plasma sphere with no outside "
        "feed has only leaky (radiating) modes. At zero frequency the same l = 1 field is the spheromak, a laboratory "
        "force-free current ring that decays through its resistance (µs to ms in the lab). A ring of electrons in a "
        "rising magnetic field is the electron-ring compressor of the electron-ring accelerators, whose radius follows "
        "R proportional to B<super>−½</super> from conservation of canonical angular momentum. The example sits beside the other "
        "electron-ring models: <i>kron_electron</i> (H21), <i>space_vortex</i> (H13) and <i>photon_electron</i> (H23).")

    r.page_break()
    r.h1("Results")
    r.figure(os.path.join(HERE, "plasmoid_rings.png"),
             "Figure 1. Left: field lines of Nachamkin's l = 1 force-free field, lying on nested tori around the current "
             "ring (section 1). Middle: the 1 MJ, 10 cm spheromak ring in kronspice, energy left against time on a log "
             "axis at electron temperatures 1, 10 and 100 eV (section 4). Right: the collapsing electron ring, ring "
             "radius (cm), field × 10 (T) and electron energy (γ) against time as the bank's field rises (section 5).")

    r.h2("1. The force-free field")
    r.p("His l = 1 field (eq. 8) is identical to the spheromak field. It obeys curl F = αF to 2×10<super>−10</super> and "
        "div F = 0 to 10<super>−11</super>. A traced field line keeps its poloidal flux ψ to 10<super>−14</super>, so it stays on "
        "its torus, and the tori nest around a current ring at 0.611 of the sphere's radius (61.1 mm for R = 10 cm).")

    r.h2("2. His critical frequency is the cold-plasma cutoff")
    r.p("Appendix B recomputed (R = 10 cm, wavenumber k = 2686 /m):")
    r.table([
        ["Quantity", "Here", "Nachamkin"],
        ["Free charge density ρ = (m/e)ε<sub>0</sub>c²(k² − α²)", "−32.63 C/m³", "−32.63 C/m³"],
        ["Free electrons", "8.532×10<super>17</super>", "8.532×10<super>17</super>"],
        ["Field E<sub>0</sub> = C·R·ω²/|e/m| for his C = 0.01 at 126 GHz", "3.566×10<super>9</super> V/m", "3.566×10<super>9</super> V/m"],
    ], [0.56, 0.22, 0.22])
    r.bullets([
        "That density's plasma frequency ω<sub>p</sub> satisfies ω² = ω<sub>p</sub>² + α²c² to −4×10<super>−17</super>: the plasmoid "
        "is a standing wave of a cold plasma with wavenumber α, and his critical frequency ω<sub>0</sub> = αc is that wave's cutoff.",
        "At his chosen field the electrons oscillate at ω·C·R = 7.92×10<super>8</super> m/s = 2.64 c, above the speed of light, "
        "so his non-relativistic electron model does not hold for his own example. The field is also about 1,000 times "
        "air's breakdown strength.",
    ])

    r.h2("3. Does the resonant size stop radiation?")
    r.p("Exact l = 1 resonances of a plasma sphere with his radius and density, with outgoing waves only (nothing "
        "feeding it from outside):")
    r.table([
        ["Part of his field", "Frequency", "Size it needs (αR)", "Q", "Energy lifetime"],
        ["TE (his M part)", "128.1583 GHz", "4.4933 (his choice)", "4.8×10<super>5</super>", "6×10<super>−7</super> s"],
        ["TM (his N part)", "128.1471 GHz", "2.7437", "9.4×10<super>5</super>", "1.2×10<super>−6</super> s"],
    ], [0.22, 0.18, 0.22, 0.16, 0.22])
    r.bullets([
        "The plasmoid radiates: its energy leaks away in about a microsecond. At 1 MJ the initial radiated power is "
        "about 2×10<super>12</super> W.",
        "The two halves need different sizes. An equal mix of them, which is what makes E parallel to H, is not one "
        "free mode of the sphere; holding it needs outside driving or surface currents, the weak point Nachamkin names himself.",
    ])

    r.h2("4. The spheromak ring in kronspice")
    r.p("At zero frequency the field is a force-free current ring: there is no Lorentz force on the current, so the ring "
        "neither pinches nor bursts, and electrical resistance ends it. The field decays as one mode with lifetime "
        "μ<sub>0</sub>/(ηα²), with η the Spitzer resistivity. As a circuit the ring is an inductor L = 2W/I² carrying the "
        "ring current I, in series with a resistor R = P/I² (W the field energy, P the ohmic power). For 1 MJ in a 10 cm "
        "sphere: field scale B<sub>0</sub> = 80 T, ring current 7.43 MA, magnetic pressure 2.5×10<super>4</super> atm, held only "
        "by a wall or an outside field.")
    r.table([
        ["Electron temperature", "Ring L", "Ring R", "Energy lifetime (theory)", "kronspice", "Heat"],
        ["1 eV", "3.62×10<super>−8</super> H", "0.0303 Ω", "5.98×10<super>−7</super> s", "5.98×10<super>−7</super> s", "9.997×10<super>5</super> J of 10<super>6</super> J"],
        ["10 eV", "3.62×10<super>−8</super> H", "0.000957 Ω", "1.89×10<super>−5</super> s", "1.89×10<super>−5</super> s", "9.997×10<super>5</super> J of 10<super>6</super> J"],
        ["100 eV", "3.62×10<super>−8</super> H", "3.03×10<super>−5</super> Ω", "0.000598 s", "0.000598 s", "9.997×10<super>5</super> J of 10<super>6</super> J"],
    ], [0.15, 0.14, 0.14, 0.19, 0.16, 0.22])
    r.p("To last one second the ring would need T = 14082 eV (163 million K), fusion-plasma temperatures. Ball "
        "lightning reported to last seconds cannot be a resistive spheromak in air.")

    r.h2("5. A collapsing electron ring, driven by kronspice")
    r.p("Declared setup: a 1000 µF bank at 10 kV is discharged into a 10 µH solenoid with 16 turns per metre (kronspice "
        "solves this circuit); 2 MeV electrons circle the axis in a 0.05 T bias field. As the coil field rises, the induced "
        "electric field speeds the electrons up and the stronger field bends them tighter, so the ring collapses.")
    r.table([
        ["Quantity", "Start", "End", "Law"],
        ["Field", "0.050 T", "1.984 T (in a quarter period)", ""],
        ["Ring radius", "16.4 cm", "2.60 cm", "R proportional to B<super>−½</super>: 2.60 cm"],
        ["Momentum", "1", "×6.3007", "×6.2989"],
        ["Electron energy", "2.0 MeV", "15.0 MeV", ""],
        ["Ring current (10<super>13</super> electrons)", "456 A", "2933 A", ""],
    ], [0.3, 0.17, 0.28, 0.25])
    r.bullets([
        "The electrons' angle round the axis is a cyclic coordinate (as in <i>kaluza_kron</i>), so their canonical angular "
        "momentum is conserved, here to 6×10<super>−4</super>.",
        "The work of the induced field, 2.0808×10<super>−12</super> J per electron, equals the energy gain, 2.0808×10<super>−12</super> J.",
        "The bank (5×10<super>4</super> J) drives the coil and the ring takes 21 J of it. The \"void\" inside the ring is field, "
        "not a source: every joule the ring gains comes from the bank.",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "Spheromak decay: the kronspice energy lifetime equals the theoretical μ<sub>0</sub>/(ηα²) at all three temperatures, "
        "and the heat (9.997×10<super>5</super> J) accounts for the 10<super>6</super> J stored.",
        "Electron ring: kronspice ledger −10<super>−9</super>; induced-field work equals the electron energy gain.",
        "Program conclusion: no energy appears from the plasmoid or the void inside the ring; radiation, resistance and "
        "the bank account for all of it.",
        "<i>main_test.go</i> checks the field, Appendix B, the modes, the spheromak decay and the ring compression.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "The plasma sphere's resonances are computed for a uniform cold plasma with only outgoing waves; outside "
        "driving, surface currents and the swirling flow are not modelled.",
        "The spheromak is treated as one resistive decay mode with Spitzer resistivity; the wall or outside field that "
        "would hold its 2.5×10<super>4</super> atm pressure is not modelled.",
        "Nachamkin's own suggested extra ingredients (quantum effects, outside fields) have no equations in the paper "
        "and are not modelled.",
    ])

    r.h1("Status")
    r.p("Register entry H26. The paper's field mathematics is sound and its electron density follows from standard "
        "plasma physics. Two parts fail: the worked example needs electron quiver speeds of 2.6 c, so the "
        "non-relativistic model breaks; and the resonant size does not stop radiation (energy lifetime about 1 µs), "
        "with the TE and TM halves needing different sizes. A 1 MJ, 10 cm spheromak lives 19 µs at 10 eV and needs "
        "14 keV to last 1 s — not ball lightning in air. Long-lived ball-lightning plasmoids, if electromagnetic, need "
        "something the paper lacks. The real, measured force-free rings are spheromaks and compressed electron rings, "
        "both modelled here with kronspice.")

    n0 = len(r.s)
    r.h1("Reproducing the results")
    r.code("go run ./examples/plasmoid_rings -viz examples/plasmoid_rings/viz   # about 8 seconds\n"
           "python3 examples/plasmoid_rings/viz/render.py\n"
           "go test ./examples/plasmoid_rings")
    r.small("Sources: J. Nachamkin, PL-TR-92-3044 (1992), DTIC AD-A257 765; S. Chandrasekhar, P. C. Kendall, <i>Astrophys. J.</i> "
            "126, 457 (1957); M. N. Rosenbluth, M. N. Bussac, <i>Nucl. Fusion</i> 19, 489 (1979); F. F. Chen, <i>Introduction to "
            "Plasma Physics and Controlled Fusion</i>; J. A. Stratton, <i>Electromagnetic Theory</i> (1941), §9.22; NRL Plasma "
            "Formulary; V. I. Veksler, Proc. CERN Symposium on High Energy Accelerators (1956).")
    r.s[n0:] = [KeepTogether(r.s[n0:])]
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "plasmoid_rings_report.pdf"))
