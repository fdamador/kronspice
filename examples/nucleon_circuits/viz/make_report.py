# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H31 report: the proton and the neutron as circuits of light, as a PDF.

    python examples/nucleon_circuits/viz/make_report.py examples/nucleon_circuits/nucleon_circuits_report.pdf

Numbers come from examples/nucleon_circuits/run_output.txt (the Go program),
skyrme/run_output.txt (the Skyrme knot solver) and the two README.md files
(measured values, tests, method); figures are skyrme/fig1_profiles.png to
skyrme/fig5_knots3d.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401

SKY = os.path.join(HERE, "..", "skyrme")


def build(out):
    r = Report("examples/nucleon_circuits", "Nucleon circuits",
               "Are the proton and the neutron circuits of light? — register entry H31",
               "register entry H31")

    r.h1("Summary")
    r.p("The light-ring picture used in this repository for the electron (H23) and the neutrino (H28) is carried to "
        "the nucleons. A light ring is light of energy mc² going once round a ring of radius ħ/(mc); a charge q carried "
        "round it at c gives the Dirac moment qħ/2m exactly. Each ring is built as a kronspice inductor with Maxwell "
        "mutual coupling to its neighbours. The proton is tried as one ring, the neutron as the rings of its decay "
        "products (p + e + ν̄), and the standard accounts — three quark rings, and Skyrme's knot, computed here — are "
        "set beside them. <b>H31 is a hypothesis</b>; the measured values and the standard models are the yardstick.")
    r.bullets([
        "<b>One light ring is not the proton:</b> radius 0.2103 fm against the measured 0.84075 fm (4.0×), moment "
        "1.0000 μN against 2.7928 μN (2.79×).",
        "<b>The neutron as p + e + ν̄ rings</b> gets charge 0, spin ½ and instability (+0.7823 MeV) right, but with "
        "each ring at its own size its moment is −1835.5 μN, 959× too large; the electron ring that fits −1.9130 μN "
        "(0.9897 fm) carries 199.4 MeV against 1.293 MeV available (154×).",
        "<b>Three quark rings with one number</b> (ring energy 335.96 MeV, fixed by μ<sub>p</sub>) predict "
        "μ<sub>n</sub> = −1.8619 μN against the measured −1.9130 (−2.7 %), and electromagnetism makes the neutron "
        "1.053 MeV lighter (lattice QED part −1.00 MeV).",
        "<b>Skyrme's knot, computed:</b> μ<sub>p</sub> +1.880, μ<sub>n</sub> −1.324, |μ<sub>p</sub>/μ<sub>n</sub>| "
        "1.419 (measured 1.460), radius 0.595 fm (0.72); helium-4 as one B = 4 cube binds 267 MeV (7.7 %) against "
        "28.3 MeV (0.75 %).",
        "<b>Lifetime:</b> the light rings carry no weak coupling, so they give none; the standard model gives 880.5 s "
        "against the measured 878.4 s.",
        "<b>Status: bounded.</b> The decay-products picture fails on moment and energy by factors of 150 to 1,000.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis H31</b> (user hypothesis, 2026-11-05): the proton and the neutron are circuits of light, as the "
        "electron and the neutrino are pictured in H23 and H28, the neutron being the rings of its decay products. "
        "<b>Standard prediction beside it:</b> the constituent-quark picture (three spin-weighted quark rings) and "
        "Skyrme's knot of the pion field, with the measured nucleon properties any model has to reproduce:")
    r.table([
        ["", "proton", "neutron"],
        ["charge", "+e", "0"],
        ["spin", "½", "½"],
        ["mass", "938.272 MeV", "939.565 MeV (1.293 MeV heavier)"],
        ["magnetic moment", "+2.7928 μN", "−1.9130 μN"],
        ["size", "charge radius 0.84075 fm", "&lt;r²&gt; = −0.1155 fm²"],
        ["lifetime", "stable", "878.4 s (n → p e ν̄, 0.782 MeV released)"],
    ], [0.25, 0.35, 0.4])
    r.small("Measured values: CODATA 2022 (masses, moments, r<sub>p</sub>); PDG 2024 (&lt;r²&gt;<sub>n</sub>, τ<sub>n</sub>, "
            "V<sub>ud</sub>, λ). μN is the nuclear magneton.")

    r.h1("The model: a light ring as a circuit")
    r.p("Light of energy mc² going once round a ring of circumference h/(mc) sits at radius R = ħ/(mc). A charge q "
        "carried round at c is a current qc/(2πR) on an area πR², so its moment is qcR/2 = qħ/(2m): the Dirac moment, "
        "exactly (TestLightRingDirac). In kronspice each ring is an inductor (a thin ring of core radius R/2, "
        "L = μ<sub>0</sub>R(ln 16 − 7/4)) carrying the ring's current from a current source, coupled to its neighbours "
        "by Maxwell's coaxial-circle mutual inductance. The ledger books every ½LI² and every M·I<sub>1</sub>·I<sub>2</sub>, "
        "matching the analytic sum to 1e-6 (TestCircuitTwin). The four circuits (one-ring proton, p + e rings, uud "
        "and udd) are exported as netlists, KiCad schematics and results sheets in export/.")

    r.h1("Results")
    r.h2("1. The proton as one light ring")
    r.table([
        ["quantity", "one light ring", "measured"],
        ["radius ħ/(m<sub>p</sub>c)", "0.2103 fm", "0.84075 fm (4.0× larger)"],
        ["current", "3.635×10<super>4</super> A", ""],
        ["moment (kronspice currents)", "1.0000 μN", "2.7928 μN (2.79× larger)"],
        ["magnetic self-energy ½LI²", "1.1143 MeV (ledger 1.1143 MeV, residual 0.0e+00)", ""],
        ["ring carrying +e at c with the measured moment", "radius 0.5874 fm, a light ring of 335.96 MeV", "= m<sub>p</sub>c²/2.7928"],
    ], [0.36, 0.4, 0.24])
    r.p("One ring reproduces the moment of a point particle; the measured moment is 2.79 times that. The ring that "
        "would carry the measured moment holds about a third of the proton's energy: <b>one ring cannot be the whole "
        "proton.</b>")

    r.h2("2. The neutron as the light rings of its decay products")
    r.p("The neutron decays to a proton, an electron and an antineutrino; in the light model these are three rings — "
        "the proton ring, the electron's photon ring (H23) and the inverted neutrino ring (H28).")
    r.table([
        ["requirement", "result", "verdict"],
        ["charge 0", "+e − e + 0 = 0", "passes"],
        ["spin ½", "three spin-½ rings can combine to ½ (Rutherford's 1920 bare p + e could not)", "passes"],
        ["unstable", "m<sub>n</sub> − (m<sub>p</sub> + m<sub>e</sub>) = +0.7823 MeV: more energy than its parts", "passes"],
        ["moment −1.9130 μN", "each ring at its own size (electron ring 386.2 fm): −1835.5 μN, 959× too large", "<b>fails</b>"],
        ["moment, refitted", "electron ring of radius 0.9897 fm gives −1.9130 μN (kronspice: currents +1.302×10<super>4</super> A "
         "and −7724 A; magnetic energy 0.1300 MeV, analytic 0.1300, residual 0.0e+00)", "fits by construction"],
        ["energy of that ring", "light on a 0.9897 fm ring carries 199.4 MeV against 1.293 MeV available: 154× too much", "<b>fails</b>"],
        ["binding", "Coulomb attraction of the two rings −1.617 MeV (own Coulomb energies +1.724 MeV)", "far too weak"],
    ], [0.18, 0.62, 0.2])
    r.p("The antineutrino ring is neutral (an E ring: no current, no moment, energy &lt; 0.45 eV), so it changes none "
        "of this. The decay-products picture gets the charge, the spin and the instability right, and fails on the "
        "moment and the energy by two to three orders of magnitude — the same reason Rutherford's 1920 neutron was "
        "abandoned after Chadwick's 1932 discovery.")

    r.h2("3. Standard counterpart: three quark light rings")
    r.p("Each quark is a light ring (u: +2/3 e, d: −1/3 e), and the nucleon's spin-½ state is a weighted mix of spin "
        "arrangements (for the proton u↑u↑d↓ with weight 2/3, u↑u↓d↑ and u↓u↑d↑ with 1/6 each: the SU(6) wavefunction). "
        "Each arrangement is a three-ring kronspice circuit, and its moment is the I·πR² of the solved currents.")
    r.table([
        ["", "quark rings", "measured"],
        ["quark ring energy (the one free number, from μ<sub>p</sub>)", "335.96 MeV, radius 0.5874 fm", "—"],
        ["proton moment", "+2.7928 μN (fitted)", "+2.7928"],
        ["<b>neutron moment</b>", "<b>−1.8619 μN (predicted)</b>", "−1.9130 (−2.7 %)"],
        ["μ<sub>n</sub>/μ<sub>p</sub>", "−0.6667 (= −2/3)", "−0.6850"],
        ["three rings' light energy", "1007.9 MeV", "938.3 (p) / 939.6 (n): 7.4 % high; binding not in the model"],
        ["electromagnetic energy (rings 0.587 fm apart)", "p 1.7075 MeV, n 0.6547 MeV", ""],
        ["m<sub>n</sub> − m<sub>p</sub> from electromagnetism", "−1.053 MeV; at 0.5R −1.285, 2.0R −0.844, 4.0R −0.688 MeV",
         "lattice QED part −1.00 MeV"],
    ], [0.38, 0.32, 0.30])
    r.p("The ring that fits the proton's moment in Part 1 is exactly one quark ring; with that one number the neutron's "
        "moment follows to 2.7 %. The ring spacing is a declared choice, not fitted. Electromagnetism alone makes the "
        "neutron about 1 MeV <i>lighter</i>, in line with the lattice's QED part; the measured +1.293 MeV needs the d "
        "ring heavier than the u ring (lattice QCD part +2.52 MeV; Borsanyi et al., <i>Science</i> 347, 1452, 2015). "
        "Ledger residuals 0.0e+00, 0.0e+00.")

    r.h2("4. The knot picture: Skyrme's model, computed")
    r.p("In Skyrme's model the proton and the neutron are one knot of the pion field (winding number B = 1) in two "
        "states of internal rotation; helium-4 is one knot of winding 4, with cubic symmetry. skyrme/skyrme.py solves "
        "the hedgehog profile F(r) (F(0) = π, F(∞) = 0) on a stretched grid r = 1.5·tan(πx/2) of 400 points by "
        "L-BFGS-B, in units 1/(eF<sub>π</sub>) for length and F<sub>π</sub>/e for energy; the B = 4 cube uses the "
        "rational-map ansatz of Houghton, Manton &amp; Sutcliffe (1998), angular integral I = 20.65. The constants e and "
        "F<sub>π</sub> are fitted to the nucleon and Δ masses, as Adkins, Nappi &amp; Witten (ANW 1983) did.")
    r.table([
        ["", "this computation", "ANW 1983", "measured"],
        ["B = 1 energy (F<sub>π</sub>/e)", "36.462", "36.5", "—"],
        ["e, F<sub>π</sub>", "5.437, 129.1 MeV", "5.45, 129 MeV", "F<sub>π</sub> = 186 MeV (this convention)"],
        ["isoscalar radius", "0.595 fm", "0.59 fm", "0.72 fm"],
        ["μ<sub>p</sub>", "+1.880 μN", "+1.87", "+2.793"],
        ["μ<sub>n</sub>", "−1.324 μN", "−1.31", "−1.913"],
        ["|μ<sub>p</sub>/μ<sub>n</sub>|", "1.419", "1.43", "1.460"],
        ["g<sub>A</sub>", "not computed", "0.61", "1.23"],
        ["B = 4 energy per winding (12π² units)", "1.1366", "HMS 1.137", "—"],
        ["helium-4 binding", "267 MeV (7.7 %)", "—", "28.3 MeV (0.75 %)"],
    ], [0.32, 0.22, 0.18, 0.28])
    r.p("The solver also prints Λ = 50.70 (ANW 50.9), winding B = 1.00005, per-winding energy B = 1 1.2315 (1.232), and "
        "classical masses B = 1 866 MeV, B = 4 3196 MeV. The B = 1 knot reproduces ANW to 0.1–1 %. It gives the neutron "
        "a negative moment and the moment ratio to 3 % with no charged parts at all, but absolute moments and sizes "
        "come out about 30 % low. The cube is lighter per winding than four separate knots, so helium-4 holds together, "
        "but about ten times too tightly: the classical knot over-binds, a known shortcoming of the original model.")

    r.figure(os.path.join(SKY, "fig3_compare.png"),
             "Figure 1. The nucleon in each picture against measurement: proton moment, neutron moment and size for one "
             "light ring (grey), quark light rings (green), the Skyrme knot (orange) and measurement (black). The "
             "p + e + ν̄ rings' neutron moment, −1835.5 μN, is off the scale. No size is given for the quark rings.")
    r.figure(os.path.join(SKY, "fig1_profiles.png"),
             "Figure 2. Skyrme knots with e = 5.44, F<sub>π</sub> = 129 MeV. Left: the profile F(r)/π, how far the field "
             "has turned, from π at the centre to 0 outside, for the nucleon (B = 1) and helium-4 (B = 4). Right: the "
             "radial winding density per unit winding, 4πr²b(r)/B.")
    r.figure(os.path.join(SKY, "fig2_charge.png"),
             "Figure 3. Charge per fm of radius, 4πr²ρ(r), in the proton and the neutron, from ρ = ½ρ<sub>S</sub> ± "
             "½ρ<sub>V</sub> (ρ<sub>S</sub> the winding density, ρ<sub>V</sub> from the moment of inertia density, each "
             "normalized). The neutron has a positive core and a negative skin, the sign pattern electron scattering "
             "measures.", width_frac=0.85)
    r.figure(os.path.join(SKY, "fig4_helium.png"),
             "Figure 4. Helium-4 as one Skyrme knot (rational-map approximation). Left: classical knot mass, 866 MeV for "
             "one B = 1 knot against 799 MeV per winding in the B = 4 cube, so helium holds together. Right: binding "
             "energy as a fraction of mass, 267 MeV (7.7 %) for the knot against the measured 28.3 MeV (0.75 %).")
    r.figure(os.path.join(SKY, "fig5_knots3d.png"),
             "Figure 5. Surfaces of high winding density in 3-D: a ball for the nucleon (B = 1) and a cube for helium-4 "
             "(B = 4). The integrated winding inside the plotted box is 0.95 and 3.99.", width_frac=0.9)

    r.h2("5. Neutron decay")
    r.bullets([
        "<b>Energy budget:</b> 939.5654 − 938.2721 − 0.5110 = 0.7823 MeV, shared by the electron and the antineutrino.",
        "<b>Lifetime:</b> the light rings carry no weak coupling, so the light model gives no lifetime. The standard "
        "model gives τ = 4908.6 s/(V<sub>ud</sub>²(1 + 3λ²)) = 880.5 s with V<sub>ud</sub> = 0.97373 and λ = 1.2754 "
        "(Czarnecki, Marciano &amp; Sirlin, <i>PRL</i> 120, 202002, 2018), against the measured 878.4 s.",
    ])

    r.h1("What it means")
    r.bullets([
        "The light ring captures one thing exactly: a ring of light carrying a charge has the Dirac moment qħ/2m.",
        "The proton is not one ring, and the neutron is not its decay products held together: that picture fails on "
        "moment and energy by factors of 150 to 1,000.",
        "Three rings work: the quark-ring circuit predicts the neutron's moment to 2.7 % from one number, and gets the "
        "electromagnetic part of the mass difference in sign and size.",
        "The knot (Skyrme) gets the same pattern with no charged parts.",
        "Neither the light model nor the circuits supply the binding that holds three rings together (in the standard "
        "account, gluons) or the weak interaction that turns d into u in the decay.",
    ])

    r.h1("Validation")
    r.table([
        ["test", "checks", "held to"],
        ["TestLightRingDirac", "ring of radius ħ/(m<sub>p</sub>c) with +e at c: moment 1 μN from the circuit's current; "
         "radius 0.2103 fm", "1e-9; 1e-6"],
        ["TestCircuitTwin", "three coupled rings, both current signs: ledger stored energy = ½ΣLI² + ΣMI<sub>1</sub>I<sub>2</sub>; "
         "currents as set; ledger", "1e-6; 1e-9"],
        ["TestQuarkRings", "μ<sub>p</sub> fitted, μ<sub>n</sub>/μ<sub>p</sub> = −2/3; electromagnetic m<sub>n</sub> − "
         "m<sub>p</sub> near −1 MeV", "1e-9; ±0.5 MeV"],
        ["TestNeutronFromDecayRings", "electron ring for the measured moment 0.990 fm; its energy &gt; 100× the 1.293 MeV "
         "available", "1e-3"],
        ["TestLifetime", "standard-model τ within 0.5 % of 878.4 s", "0.5 %"],
        ["TestExport", "the four ring circuits as netlists, KiCad schematics and results sheets", "—"],
    ], [0.26, 0.58, 0.16])

    r.h1("Limits of the model")
    r.bullets([
        "Rings are thin coaxial rings with a declared core radius R/2 and declared spacings; the spacing of the quark "
        "rings is chosen, not fitted.",
        "No binding force and no weak interaction in the ring circuits.",
        "The Skyrme knot is classical (no quantized vibrations), and the B = 4 cube is the rational-map approximation.",
    ])

    r.h1("Status")
    r.p("<b>H31: bounded.</b> One ring gives the proton 1 μN (measured 2.793) at 0.210 fm (0.841); the neutron as "
        "p + e + ν̄ rings gets charge 0, spin ½ and instability right but its moment 959× too large, and the electron "
        "ring that fits the moment carries 154× the available energy. The standard counterparts — quark rings "
        "(μ<sub>n</sub> −1.862 from one number) and Skyrme's knot — are computed beside it, and the lifetime comes "
        "only from the standard weak interaction.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/nucleon_circuits      # under a second\n"
           "go test ./examples/nucleon_circuits\n"
           "KRONSPICE_EXPORT=1 go test ./examples/nucleon_circuits   # export/: netlists, KiCad, results sheets\n"
           "cd examples/nucleon_circuits/skyrme\n"
           "python3 skyrme.py      # solves both knots, writes skyrme.json (a few seconds)\n"
           "python3 render.py      # fig1-fig5 from skyrme.json")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "nucleon_circuits_report.pdf"))
