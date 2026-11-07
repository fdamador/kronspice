# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Report: electrodynamic tether in low Earth orbit, as a PDF.

python examples/tether/viz/make_report.py examples/tether/tether_report.pdf

Numbers come from examples/tether/README.md (the example has no
run_output.txt and no figures) and the results sheets in export/.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/tether", "Electrodynamic tether in low Earth orbit",
               "Pushing and pulling against a planetary magnetic field, with energy and momentum ledgers — "
               "CLAUDE.md §3f", "CLAUDE.md §3f", "milestone M9n; declared v, B and contactor values")

    r.h1("Summary")
    r.p("A vertical aluminium wire moves east at orbital speed through the northward geomagnetic field. The Lorentz "
        "force on its charges, computed by grade selection (<i>em.MotionalEMF</i>), gives an EMF (v×B)·L along it, and "
        "the circuit closes through the ionosphere by plasma contactors at the two ends. The same loop runs as a "
        "generator (the orbit pays for the load) or as a thruster (a supply drives current against the EMF). The "
        "force's reaction is booked in the momentum ledger, as every push–pull design must (CLAUDE.md §3f). NASA flew "
        "this circuit three times (PMG 1993, TSS-1 1992, TSS-1R 1996) and designed it twice more (ProSEDS, the ISS "
        "reboost study); their numbers are set beside the example's.")
    r.bullets([
        "<b>EMF</b> from the Lorentz force: 1155.0 V, equal to vBL.",
        "<b>Generator</b> (500 Ω load): 1.330 A, 884 W to the load, drag 0.199 N; F·v = −1536 W, so the orbit pays.",
        "<b>Thruster</b> (3000 V supply): −5.004 A, 15,012 W from the supply, thrust 0.751 N, 5780 W into the orbit: "
        "0.05 N/kW.",
        "<b>Ledgers:</b> every energy and momentum identity holds to 10<super>−9</super>; the energy ledger closes to "
        "10<super>−12</super>; the reaction force lands on the field and ionosphere.",
        "<b>Against NASA:</b> the example's 231 V/km is an upper end (v ⟂ B at 30 µT; TSS-1R measured 177 V/km), and "
        "its 0.05 N/kW is half the ISS study's 0.1 N/kW, because 9.2 kW of the 15 kW is I²R heat in the linear "
        "contactor stand-ins.".replace("⟂", "perpendicular to"),
    ])

    r.h1("The construction")
    r.p("Declared, not looked up:")
    r.table([
        ["input", "value"],
        ["orbital speed v", "7.7 km/s"],
        ["magnetic field B", "30 µT (order of magnitude in LEO; a real orbit takes B from a cited field model such as "
         "IGRF)"],
        ["tether", "5 km of 1 mm aluminium wire (169 Ω)"],
        ["contactors", "100 Ω each, a linear stand-in for their nonlinear current–voltage characteristics"],
    ], [0.3, 0.7])
    r.p("Two modes of the same loop, each solved as a kronspice operating point with the ledger on:")
    r.code("EMF - tether resistance - contactors - load          (generator)\n"
           "EMF - tether resistance - contactors - supply > EMF  (thruster)")
    r.p("The current feels the force ∫I·dl×B. As a generator the force is drag: the power the load receives comes "
        "from the orbit's kinetic energy, F·v = −EMF·I. With a supply driving the current against the EMF the force "
        "turns into thrust. <i>em.MomentumLedger</i> books the force on the tether with its reaction on the "
        "geomagnetic field and ionosphere, the exterior source.")
    r.p("This is the circuit NASA's <i>Tethers in Space Handbook</i> writes for the electrodynamic tether (Cosmo &amp; "
        "Lorenzini 1997, §4): EMF V = ∫(v×B)·dl, which for a straight tether perpendicular to the field is "
        "V = L·v·B·sin θ; force F = ∫I·dl×B; and loop equation V<sub>ind</sub> = I·R + ΔV<sub>low</sub> + "
        "ΔV<sub>up</sub> + ΔV<sub>ion</sub> + ΔV<sub>load</sub>, where the end terms are the drops across the "
        "contactors' space-charge sheaths or plasma clouds and ΔV<sub>ion</sub> is the ionosphere's drop, of order "
        "1 to 20 Ω. The handbook notes that the sheaths are the largest impedances in the loop; the two 100 Ω "
        "resistors stand in for them, with the ionospheric drop folded in.")

    r.h1("Results")
    r.table([
        ["", "generator, 500 Ω load", "thruster, 3000 V supply"],
        ["current up the tether", "1.330 A", "−5.004 A"],
        ["electrical power", "884 W to the load", "15,012 W from the supply"],
        ["force along x", "−0.199 N (drag)", "+0.751 N (thrust)"],
        ["F·v", "−1536 W: the orbit pays", "+5780 W into the orbit"],
        ["thrust per power", "—", "0.05 N/kW"],
        ["momentum ledger", "reaction +0.199 N on the field and ionosphere",
         "reaction −0.751 N on the field and ionosphere"],
    ], [0.26, 0.37, 0.37])
    r.p("Checks in the test (all to 10<super>−9</super>): F·v = −EMF·I in both modes; generator load power plus I²R "
        "losses equals EMF·I; thruster supply power equals losses plus F·v; the momentum ledger's net force equals the "
        "force on the tether; the energy ledger closes to 10<super>−12</super>.")
    r.p("The export sheet for the generator (<i>export/tether_generator_results.md</i>) gives the MNA node voltages "
        "−132.957 V (bot), −357.26 V (m), 797.74 V (top), 664.783 V (x) and 0 (y). Automatic tearing finds no cut of "
        "this single loop whose pieces each keep a reference (Aitchison 1983, Lemma 2), so only the MNA column is "
        "filled.")

    r.h1("Against NASA's flights and studies")
    r.table([
        ["", "this example", "PMG (1993)", "TSS-1R (1996)", "ProSEDS (designed)", "ISS reboost study (1998)"],
        ["tether", "5 km, 1 mm Al, 169 Ω", "500 m, #18 AWG Teflon-insulated Cu",
         "19.7 km deployed of 20.7 km; 10 strands of 34 AWG Cu on a Nomex core, 0.305 mm FEP insulation, 2.54 mm "
         "overall", "5 km of 7-strand 28 AWG 1350-O Al, mostly bare, plus 10 km Dyneema",
         "10 km, bare lower section, &lt; 200 kg"],
        ["orbit", "v = 7.7 km/s, B = 30 µT, v perpendicular to B", "193 × 869 km, 25.7°",
         "297 km at the break, near the equator", "400 km circular, 36°", "ISS orbit"],
        ["EMF", "1155 V (231 V/km)", "bias ±65 to −130 V applied in series with the I×B EMF", "3482 V (177 V/km)",
         "—", "about 100 V/km"],
        ["current", "1.33 A generator, 5.0 A thruster", "100–300 mA by day, 10–50 mA at night; 0.3 A at 130 V bias",
         "0.97 A just before the break; 1.1 A for 75 s after it", "—", "—"],
        ["force", "0.199 N drag, 0.751 N thrust", "—", "—", "about 0.4 N drag; 14 km/day descent",
         "0.5 N from 5 kW, 0.8 N from 10 kW"],
        ["power", "884 W load; 15 kW supply", "—", "—", "about 100 W to the spacecraft", "5–10 kW from the arrays"],
        ["ends", "two 100 Ω linear contactors",
         "hollow-cathode xenon plasma contactors at both ends, current fully reversible",
         "1.6 m conducting sphere as anode, electron guns and the Orbiter ground as cathode",
         "bare wire collects, hollow cathode emits", "bare wire collects, plasma contactor on the ISS emits"],
    ], [0.09, 0.15, 0.18, 0.21, 0.18, 0.19])
    r.p("What the comparison says about the example's declared inputs:")
    r.bullets([
        "<b>EMF per kilometre.</b> The example's 231 V/km is the upper end: it takes v perpendicular to B at 30 µT. "
        "TSS-1R measured 177 V/km near the equator at 297 km, and the ISS study quotes about 100 V/km for the "
        "station's 51.6° orbit. The handbook's sin θ is the difference, which is why a real orbit must take B from a "
        "field model.",
        "<b>Thrust per watt.</b> The example's thruster gives 0.05 N/kW against the study's 0.1 N/kW at 5 kW and "
        "0.08 N/kW at 10 kW. Of the 15 kW drawn, 9.2 kW is I²R heat in the 369 Ω of tether and contactor resistance "
        "and 5.8 kW goes into the orbit. NASA's design cuts both: a bare tether collects electrons along its length "
        "with no end sheath to cross (Sanmartín, Martínez-Sánchez &amp; Ahedo 1993), and a hollow cathode emits them "
        "at a drop of tens of volts.",
        "<b>The contactors</b> are the weakest part of the model. TSS-1R's sphere collected currents up to three times "
        "what the Parker–Murphy model predicted, at lower satellite potentials than predicted, and Stone &amp; Raitt "
        "conclude that the system's motion through the plasma affects collection. A nonlinear contactor "
        "characteristic fitted to these flight curves is the next step for this example.",
        "<b>Reversibility.</b> PMG demonstrated in flight what the two columns of the results table assume: the same "
        "tether and contactors run as a generator and as a motor.",
    ])
    r.h2("The TSS-1R break")
    r.p("On 25 February 1996 (057/01:29:16 GMT), during STS-75, with 19.7 km of the planned 20.7 km deployed and the "
        "EMF at 3482 V, the tether parted at the Orbiter end. The failure board (Szalai 1996) traced it to an arc "
        "through a flaw in the insulation that burned the tether through. Gilchrist et al. (1998) report that the "
        "current, which had reached 0.97 A, rose to 1.1 A as the bare broken end entered the plasma, and stayed there "
        "for about 75 s: a few short strands of copper emitted more current than the Orbiter's conducting surfaces, by "
        "a gas-enhanced discharge later reproduced in a plasma chamber. The lesson for this example is the one §3e "
        "applies to coils: every design is checked against insulation breakdown at its peak voltage, here the full "
        "1155 V between the tether and the surrounding plasma when the far end is open.")

    r.h1("Limits of the model")
    r.bullets([
        "v and B are declared (v perpendicular to B at 30 µT), not taken from an orbit and a field model such as IGRF.",
        "The contactors and the ionosphere are two linear 100 Ω resistors; real sheaths are nonlinear and depend on "
        "the motion through the plasma.",
        "Each mode is a single kronspice DC operating point of one loop at fixed v and B."
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/tether\n"
           "go test ./examples/tether\n"
           "KRONSPICE_EXPORT=1 go test ./examples/tether -run TestExport    # netlists, KiCad schematics, results sheets in export/")
    r.small("Sources: M. L. Cosmo, E. C. Lorenzini (eds.), <i>Tethers in Space Handbook</i>, 3rd ed., 1997, NTRS "
            "19980018321; L. Johnson, M. Herrmann, NASA/TM-1998-208538, NTRS 19980223479; N. H. Stone, W. J. Raitt, "
            "Adv. Space Res. (1998), NTRS 19990084046; B. E. Gilchrist et al., Geophys. Res. Lett. 25(4), 437 (1998); "
            "K. Szalai, TSS-1R Mission Failure Investigation Board, 1996; L. Johnson et al., AIAA 1998, NTRS "
            "19990107381; J. Ballance, L. Johnson, NTRS 20010020385; L. Curtis et al., STAIF 2002, NTRS 20020048608; "
            "I. Katz et al., J. Geophys. Res. 100 (1995); J. R. Sanmartín et al., J. Propulsion and Power 9(3), 353 "
            "(1993); L. W. Parker, B. L. Murphy, J. Geophys. Res. 72, 1631 (1967); P. Joot, Theorem 3.16 (reference "
            "only); P. W. Aitchison, Quart. Appl. Math. 41, 265 (1983); doc/em-push-pull.md; doc/M9n.md.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "tether_report.pdf"))
