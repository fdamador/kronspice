# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/toroidal_coherence: Konstapel 2025 toroidal coherence propulsion (H35).

Run from the repository root:

    python examples/toroidal_coherence/viz/make_report.py examples/toroidal_coherence/toroidal_coherence_report.pdf

Numbers come from examples/toroidal_coherence/README.md and the H35 row of the
register in CLAUDE.md. The example has no run_output.txt and no figures.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/toroidal_coherence", "Toroidal coherence propulsion",
               "What does standard physics predict for the Tier 1 bench of Konstapel's toroidal-coherence paper, "
               "and for Hall's reflection reading of lensed galaxies? — register entry H35",
               "register entry H35")

    r.h1("Summary")
    r.p("J. Konstapel, \"Toroidal Coherence Engineering: From Hall's Photon Theory to Practical Propulsion Systems\" "
        "(Leiden, 1 Nov 2025, 17 pp.) proposes that heated fiber-optic toroidal coils, phase-locked to the Earth's "
        "electromagnetic field, invert a scalar \"Sn field\" (identified with the scalar part of Maxwell's quaternion "
        "equations, after V. Robinson and C. Hall) and reduce the inertial mass of a craft. The paper gives no equation "
        "with a value for its coupling σ, so the hypothesis cannot be put into the simulator as a hypothetical model. "
        "This example computes the <b>standard prediction</b> for each of the paper's three Tier 1 measurements (its "
        "§6.1), the ordinary effects a bench must remove, the paper's own figures, and Hall's explanation of lensed "
        "galaxy images.")
    r.bullets([
        "<b>The (3n ± 1) harmonic signature is ordinary symmetry.</b> Three coils at 120° with an ordinary diode "
        "nonlinearity, in a kronspice transient, show the multiples of 3 only at the 10<super>−10</super> level; a single "
        "coil shows every harmonic.",
        "<b>Discrete lines are ordinary loop modes:</b> a closed 10 m silica loop has longitudinal modes every "
        "c/(n<sub>g</sub>·ℓ) = 20.53 MHz.",
        "<b>On the balance, heat makes the coil heavier by E/c² (+0.39 µg),</b> while air buoyancy (16.7 g) and feed-lead "
        "forces (1.3 g per metre) alone reach the claimed 500 mg – 50 g.",
        "<b>The paper's figures do not hold:</b> 1 kN on a 1000 t craft gives 0.001 m/s², not 1 m/s².",
        "<b>Hall's reflection at hydrogen clouds turns light by at most 1.3×10<super>−3</super> arcsec</b> (dense core), "
        "against 1–30 arcsec observed in strong lensing, and does so chromatically (2.1 % more for blue than for 1 µm).",
        "<b>Status: untested.</b> Tests proposed: a vacuum balance with a dummy coil at equal power and a written "
        "coherence protocol; multiband image positions of a lensed quasar.",
    ])

    r.h1("The claim and the standard prediction")
    r.table([
        ["Tier 1 measurement or claim", "Paper (hypothesis)", "Standard prediction (this example)"],
        ["Harmonics of the coil array", "a (3n ± 1) \"topological signature\"",
         "any balanced three-coil array with an ordinary nonlinearity radiates only 3n ± 1; one coil radiates every harmonic"],
        ["Infrared emission", "discrete lines from a coherent coil", "longitudinal loop modes every c/(n<sub>g</sub>·ℓ); "
         "thermal light filtered, not phase-locked"],
        ["Weight on a balance", "0.001 %–0.1 % lighter (500 mg – 50 g)", "+0.39 µg (E/c²); ordinary forces of mg to g size"],
        ["Galaxy lensing", "frustrated total internal reflection at hydrogen clouds (Hall)",
         "gravitational lensing: achromatic, 1–30 arcsec separations"],
    ], [0.24, 0.3, 0.46])
    r.p("CLAUDE.md §3c states that the quaternion form of Maxwell's equations is equivalent to the vector form; its "
        "scalar part adds no physics. Related entries: H23 (the toroidal electron the paper builds on), H18 "
        "(inertial-mass reduction of a craft), H8 and H10 (scalar-field electrodynamics and Bearden's claims, untested).")

    r.h1("Part 1 — the (3n ± 1) \"topological signature\"")
    r.p("Three coils at 120°, each 1 mH + 2 Ω, closed through 20 Ω in parallel with a diode (an ordinary rectifying "
        "contact), driven at 10 kHz (inside the paper's 1–100 kHz), solved by kronspice transient. The array field is "
        "the space vector (2/3)(i<sub>a</sub> + a·i<sub>b</sub> + a²·i<sub>c</sub>), with a = e<super>j2π/3</super>, the "
        "field on the axis of three coils at 120°. Harmonic amplitudes relative to the fundamental:")
    r.table([
        ["n", "One coil", "Array, forward", "Array, backward", "Linear control"],
        ["1", "1", "1", "1.7×10<super>−9</super>", "1"],
        ["2", "3.25×10<super>−2</super>", "1.5×10<super>−10</super>", "3.25×10<super>−2</super>", "2.4×10<super>−15</super>"],
        ["3", "1.06×10<super>−3</super>", "5.5×10<super>−10</super>", "5.3×10<super>−10</super>", "1.6×10<super>−15</super>"],
        ["4", "3.12×10<super>−3</super>", "3.12×10<super>−3</super>", "1.0×10<super>−10</super>", "1.0×10<super>−15</super>"],
        ["5", "2.58×10<super>−4</super>", "3.0×10<super>−10</super>", "2.58×10<super>−4</super>", "1.3×10<super>−15</super>"],
        ["6", "8.50×10<super>−4</super>", "1.0×10<super>−10</super>", "1.2×10<super>−10</super>", "1.7×10<super>−16</super>"],
        ["7", "1.01×10<super>−4</super>", "1.01×10<super>−4</super>", "2.1×10<super>−10</super>", "8.0×10<super>−16</super>"],
    ], [0.08, 0.23, 0.23, 0.23, 0.23])
    r.bullets([
        "Harmonics 3n + 1 (4, 7) rotate forward, 3n − 1 (2, 5) backward; multiples of 3 cancel by symmetry. The linear "
        "control shows no harmonics above round-off.",
        "Ledger: 11.98 mJ in − 8.65 mJ returned = 10.48 mJ dissipated − 7.16 mJ released from storage (start-up); "
        "integration mismatch 1.6 nJ, worst point residual 0.",
        "The paper's Tier 2 array would show 3n ± 1 with no new physics. For the single Tier 1 coil the only unexpected "
        "result would be the <i>absence</i> of multiples of 3, which the paper does not state.",
    ])

    r.h1("Part 2 — discrete infrared lines")
    r.p("A closed 10 m silica loop (group index n<sub>g</sub> 1.46, declared) has longitudinal modes every "
        "c/(n<sub>g</sub>·ℓ) = 20.53 MHz; every closed loop shows discrete lines (compare the ring lines of "
        "<i>examples/photon_electron</i>, AC peaks at n·v/ℓ). Thermal light at 1000 K peaks at 2.90 µm with a coherence "
        "length of order that wavelength, 3×10<super>−7</super> of the loop: the loop filters it into modes but cannot "
        "phase-lock it. A test needs a stated line spacing or linewidth that differs from c/(n<sub>g</sub>·ℓ).")

    r.h1("Part 3 — the balance")
    r.p("Rig (declared where the paper is silent): 50 kg on the balance, 100 kW heater (paper), coil at 1000 K (paper), "
        "10 m (paper) × 50 mm bundle, emissivity 0.8, silica c<sub>p</sub> ≈ 1000 J/(kg·K), 250 A DC feed, 50 µT.")
    r.table([
        ["Force on the balance", "Equivalent mass"],
        ["Claim: 0.001 %–0.1 % lighter", "500 mg – 50 g"],
        ["Standard prediction: heat adds E/c² (heavier)", "+0.39 µg"],
        ["Radiation recoil, all emission one way (P/c)", "≤ 34 mg"],
        ["Feed lead in the geomagnetic field", "1.3 g per metre"],
        ["Air buoyancy lost round the hot coil (in air)", "16.7 g"],
    ], [0.6, 0.4], bold_rows=(2,))
    r.bullets([
        "The coil radiates 71 kW at 1000 K, so the heater holds the temperature.",
        "Ordinary effects reach the claimed band. A bench must run in vacuum, cancel lead forces (twisted, reversed "
        "feed), compare against a dummy coil without fiber at the same power, cycle the drive on and off, and fix the "
        "\"coherence protocol\" in writing before the run (the paper attributes null results to poor phase-lock, §5.3, "
        "§9.1). A well-shielded balance can then bound the effect below 10<super>−6</super> of the rig mass.",
        "Existing bound: MICROSCOPE finds no composition-dependent free fall at the 10<super>−15</super> level, "
        "η(Ti, Pt) = [−1.5 ± 2.3 (stat) ± 1.5 (syst)]×10<super>−15</super> (Touboul et al., PRL 129, 121102, 2022), which "
        "leaves little room for gravity as an electromagnetic coupling.",
    ])

    r.h1("Part 4 — the paper's figures")
    r.bullets([
        "1 kN on a 1000 t craft gives 0.001 m/s², not 1 m/s² (§4.3); the 5 and 100 m/s² rows inherit the ×1000 error.",
        "6000 g on 1000 t with 1 kN needs a residual mass fraction 1.7×10<super>−8</super>, not the \"&gt;99.99 %\" "
        "(10<super>−4</super>) the paper states.",
        "Coil length: 3–10 km per coil (§4.1), 3,000–5,000 km per coil (§4.3), ~2,000 km in total (§6.3).",
        "Aluminium cladding (melts at 933 K) is run at 1,500–2,000 K; fused silica softens near 1,940 K (paper: melts at "
        "2,400 K).",
        "Lock to the 7.8 Hz Schumann mode within 0.1 %: that mode has Q ≈ 5 (Sentman 1995; <i>examples/earth_machine</i>), "
        "a linewidth of about 20 %.",
        "E<sub>eff</sub> = (1 − f)·γ·m<sub>0</sub>·c² for v &gt; c (§8.1): γ is not real above c, whatever the mass.",
    ])

    r.h1("Part 5 — Hall's \"FTIR lensing\" of galaxy images")
    r.p("Hall's photon theory explains lensed galaxy images as total (and frustrated total) internal reflection at thin "
        "hydrogen clouds with an empty gap between them, not as gravitational lensing. Total reflection from hydrogen to "
        "vacuum happens only for rays within the critical grazing angle √(2(n − 1)) of a sharp boundary, and turns them "
        "by at most twice that. With H<sub>2</sub>'s measured refractivity (Peck &amp; Huang, JOSA 67, 1550, 1977; "
        "n − 1 = 1.388×10<super>−4</super> at STP, 589 nm), scaled to cloud density:")
    r.table([
        ["H<sub>2</sub> density", "n − 1", "Largest deflection"],
        ["1 cm<super>−3</super> (diffuse)", "5.2×10<super>−24</super>", "1.3×10<super>−6</super> arcsec"],
        ["10<super>4</super> cm<super>−3</super> (molecular cloud)", "5.2×10<super>−20</super>", "1.3×10<super>−4</super> arcsec"],
        ["10<super>6</super> cm<super>−3</super> (dense core)", "5.2×10<super>−18</super>", "1.3×10<super>−3</super> arcsec"],
    ], [0.4, 0.3, 0.3])
    r.bullets([
        "Observed strong-lensing images are separated by about 1–30 arcsec, 10³–10<super>7</super> times more.",
        "Real cloud edges are graded over parsecs, not the sub-wavelength step that total reflection needs, and "
        "frustrated total reflection needs a gap no wider than about a wavelength.",
        "The reflection is chromatic: blue light is turned 2.1 % more than 1 µm light, while gravitational lensing turns "
        "every wavelength equally. <b>Test:</b> compare the image positions of a multiply imaged quasar across bands; "
        "positions that agree across bands rule out the reflection reading.",
    ])

    r.h1("Limits")
    r.bullets([
        "No hypothetical model is simulated: the paper gives no equation for its coupling σ, so only the standard "
        "prediction and the ordinary bench effects are computed.",
        "The balance rig's mass, bundle size, emissivity, feed current and field are declared where the paper is silent.",
        "Only Part 1 is a circuit simulation (kronspice transient); Parts 2–5 are closed-form estimates.",
    ])

    r.h1("Status")
    r.p("Register entry H35 is <b>untested</b>: no equation gives the coupling, so there is no hypothetical model to "
        "compare. Every claimed Tier 1 signature has an ordinary explanation of at least the claimed size. What would move "
        "the entry: a vacuum balance run with a dummy coil at equal power and a coherence protocol fixed in writing "
        "beforehand, or lensed-quasar image positions that shift with wavelength.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/toroidal_coherence\n"
           "go test ./examples/toroidal_coherence")
    r.small("Sources: J. Konstapel, \"Toroidal Coherence Engineering\" (Leiden, 1 Nov 2025); Touboul et al., PRL 129, "
            "121102 (2022); Peck &amp; Huang, JOSA 67, 1550 (1977); Sentman (1995).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "toroidal_coherence_report.pdf"))
