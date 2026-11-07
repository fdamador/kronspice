# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Light hydrogen (register entry H32) as a PDF report.

    python examples/light_hydrogen/viz/make_report.py [out.pdf]

Numbers come from examples/light_hydrogen/run_output.txt (the recorded run),
README.md and REPORT.md (left unchanged; this PDF covers it). The figures are
the SVGs in viz/ (fig1–fig9, written by `go run ./examples/light_hydrogen -svg`);
they are rasterized to PNG in a temporary directory with PyMuPDF at build time,
so nothing is written to the repository except the PDF.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import pymupdf  # noqa: E402


def svg_png(name, tmp, dpi=200):
    """Rasterize viz/<name>.svg to a PNG in tmp; return its path."""
    doc = pymupdf.open(os.path.join(HERE, name + ".svg"))
    out = os.path.join(tmp, name + ".png")
    doc[0].get_pixmap(dpi=dpi).save(out)
    return out


def main(out):
    tmpdir = tempfile.TemporaryDirectory()
    tmp = tmpdir.name

    def fig(name, caption, **kw):
        r.figure(svg_png(name, tmp), caption, **kw)

    r = Report("examples/light_hydrogen", "Light hydrogen",
               "The hydrogen atom built from the repository's light models — toroidal-photon electron, "
               "light-ring, quark-ring or Skyrme proton, ether-knot confinement — against measured hydrogen "
               "(register entry H32)",
               "register entry H32", "recorded run: run_output.txt (about 40 s)")

    r.h1("Summary")
    r.p("Hypothesis H32 (user hypothesis, 2026-11-05) builds hydrogen from three earlier light models: the "
        "electron as the toroidal photon of Williamson &amp; van der Mark (H23), the proton as a circuit of light "
        "(H31: one light ring, three quark rings, or Skyrme's knot), and the ether knot of H29 as what would hold "
        "the electron's light together. The binding is solved with Kron's 1945 radial Schrödinger network, both "
        "as a matrix and as a kronspice circuit whose frequency response carries the spectrum; the 21 cm line is "
        "a separate kronspice circuit of current loops. Every prediction is set beside measured hydrogen and the "
        "standard account (the Dirac electron with QED and the measured proton).")
    r.bullets([
        "<b>The rings are points at the atom's scale.</b> Electron ring 193.1 fm, proton light ring 0.2103 fm, "
        "Bohr radius 52918 fm (ratios 0.0036 and 4.0×10<super>−6</super>), so the atom is the Coulomb problem "
        "with Z = q<sub>e</sub>q<sub>p</sub>/e².",
        "<b>The frequency response contains the spectrum.</b> Each level is a pole of the network's "
        "driving-point impedance; poles equal the network matrix to all ten printed digits and reach "
        "−Z²/2n² to 3.1×10<super>−6</super>; with q = e the Lyman and Balmer lines fall within 0.01 nm of "
        "measured, and 1S–2S is 13.4 ppm low.",
        "<b>With the toroidal-photon charge (0.910 e) the atom fails.</b> It carries +0.0897 e (measured limit "
        "|q| &lt; 10<super>−21</super> e), binds 1S at 11.2686 eV and puts Lyman α at 146.701 nm "
        "(measured 121.5668 nm): every line is 1/Z² = 1.207 times too long.",
        "<b>With the charge set to e the levels are Bohr's</b>, and the proton is tested: the light-ring proton "
        "gives the 21 cm line at 508.623 MHz (0.35808 of measured) and shifts 1S–2S by −0.907 MHz; the quark-ring "
        "proton (moment fitted) gives 1420.506 MHz but shifts 1S–2S by −0.338 MHz, against a measurement good "
        "to 10 Hz.",
        "<b>The 21 cm circuit equals Fermi's contact energy</b> to 1.3×10<super>−5</super>; its ledger stored "
        "energy equals ½ΣLI² + ΣMII to 4.4×10<super>−16</super>.",
        "<b>The confining knot predicts an unseen particle.</b> Sized to the electron (ℓ = 160.9 fm, energy "
        "m<sub>e</sub>c²) its squashing tone is ħω = 1.737 MeV, an excited electron at 2.248 MeV; none has been "
        "seen.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H32).</b> The hydrogen atom is the toroidal-photon electron (H23) bound to a proton of "
        "light (H31: light ring, quark rings or Skyrme knot), with the electron's light held by the ether knot "
        "(H29). <b>Standard prediction.</b> The Dirac electron with QED, bound to the measured proton, gives every "
        "measured column, including the 23.0 GHz of 1S–2S and the 0.0069 % of the 21 cm line that the Bohr and "
        "Fermi formulas leave out.")
    r.table([
        ["Part", "Model", "Numbers used"],
        ["Electron (H23)", "one Compton wavelength of light on a double loop",
         "ring radius λ<sub>C</sub>/4π = 193.1 fm; charge q = √(3ε<sub>0</sub>ħc)/2π = 0.910 e (eq. 5); "
         "g = 2(1 + α′/2π), α′ = (q/e)²α (eqs. 8, 12)"],
        ["Electron, charge set to e", "the same, with the paper's eq. 13 comparison", "q = e, g = 2(1 + α/2π)"],
        ["Proton, light ring (H31)", "+e carried at c round radius ħ/(m<sub>p</sub>c)", "0.2103 fm, 1 μ<sub>N</sub>"],
        ["Proton, quark rings (H31)", "uud rings fitted to the proton's moment",
         "2.7928 μ<sub>N</sub>, rms charge radius 0.678 fm"],
        ["Proton, Skyrme knot (H31)", "B = 1 hedgehog (nucleon_circuits/skyrme)",
         "1.880 μ<sub>N</sub>; no charge radius computed"],
        ["Confinement (H29)", "the ether knot",
         "energy 547.5√(κ<sub>2</sub>κ<sub>4</sub>), core radius 1.2ℓ, squashing tone ω = 1.416 c/ℓ"],
    ], [0.2, 0.33, 0.47])
    fig("fig6_scales",
        "Figure 1. Sizes on one logarithmic axis: the electron size bound from g − 2 (0.001 fm), the proton light "
        "ring (0.21 fm), quark rings rms (0.678 fm), the measured proton charge radius (0.841 fm), the knot scale "
        "ℓ (161 fm), the electron ring (193 fm) and the Bohr radius (5.29×10<super>4</super> fm). Both rings "
        "are far smaller than the atom.")

    r.h1("The model")
    r.h2("Kron's radial Schrödinger network")
    r.p("Kron's 1945 network (<i>Phys. Rev.</i> 67, 39) for s states is a radial ladder: series inductors carry "
        "the kinetic energy, shunt inductors the potential −Z/r, unit capacitors the energy, and the wave "
        "function u = rψ is the node voltage, with short circuits at r = 0 and at the edge. Eigenvalues come from "
        "a Sturm count with Richardson extrapolation to zero spacing (h → 0). The same network is built in "
        "kronspice in the dual form (capacitors and inductors swapped, because AC analysis needs a DC operating "
        "point); the network admittance vanishes at ω² = 1/(E + V0), so each state is a pole of the "
        "driving-point impedance. The light models change only Z = q<sub>e</sub>q<sub>p</sub>/e².")
    r.table([
        ["Z = q<sub>e</sub>q<sub>p</sub>/e²", "E<sub>1s</sub> (network, h → 0)", "Exact −Z²/2",
         "E<sub>2s</sub> (h → 0)", "Exact −Z²/8"],
        ["1.000000", "−0.499999997 Ha", "−0.500000000", "−0.125000000", "−0.125000000"],
        ["0.910318 (H23)", "−0.414339490 Ha", "−0.414339492", "−0.103584873", "−0.103584873"],
    ], [0.2, 0.22, 0.18, 0.2, 0.2])
    r.p("kronspice twin (h = 0.1 a0, 299 nodes, AC resonance by bisection on the sign of the reactance): "
        "E<sub>1s</sub> = −0.4987562112 Ha against the network matrix −0.4987562112 (difference "
        "7.9×10<super>−15</super>).")
    r.h2("The 21 cm line as a circuit")
    r.p("The electron's spin, spread over the 1s cloud, is a magnetization M(r) = μ<sub>e</sub>|ψ|² ẑ. A "
        "magnetization is a current, J = curl M = (2/a)M sinθ φ̂, so the cloud is 32 coaxial current loops round the "
        "proton ring (8 radial steps × 4 Gauss–Legendre angles, exact quadrature weights). Each loop and the ring "
        "is a kronspice inductor driven to its current, coupled to the ring by Maxwell's coaxial-circle mutual "
        "inductance. The circuit gives the aligned and anti-aligned (zz) energy ±U; spin ½ (S·S rather than "
        "S<sub>z</sub>S<sub>z</sub>) makes the splitting 4U.")
    fig("fig5_loops",
        "Figure 2. The 21 cm circuit in a half-plane through the axis: shading is the 1s density |ψ|², dots are "
        "the 32 current loops (area proportional to current), and the proton ring sits at the origin (drawn "
        "enlarged; true size 4×10<super>−6</super> a0).", max_height_mm=110)

    r.h1("Results")
    r.h2("The frequency response of the whole network")
    r.p("The whole network was swept at 8000 frequencies (1299 nodes for Z = 1, 1427 for Z = 0.910) with a 1 A "
        "source at r = 1 a0. Each pole was pinned by bisection, and two networks (spacing h = 0.1 and 0.05 a0) "
        "were extrapolated to h → 0. Energies in hartree:")
    r.table([
        ["n", "Pole, h = 0.1", "Matrix, h = 0.1", "Pole, h → 0", "Exact −Z²/2n²", "Rel. diff"],
        ["<b>Z = 1</b>", "", "", "", "", ""],
        ["1", "−0.4987562112", "−0.4987562112", "−0.4999984496", "−0.5000000000", "−3.1×10<super>−6</super>"],
        ["2", "−0.1249219725", "−0.1249219725", "−0.1249999756", "−0.1250000000", "−1.9×10<super>−7</super>"],
        ["3", "−0.0555401320", "−0.0555401320", "−0.0555555534", "−0.0555555556", "−3.9×10<super>−8</super>"],
        ["4", "−0.0312451187", "−0.0312451187", "−0.0312499996", "−0.0312500000", "−1.2×10<super>−8</super>"],
        ["5", "−0.0199980004", "−0.0199980004", "−0.0199999999", "−0.0200000000", "−5.2×10<super>−9</super>"],
        ["<b>Z = 0.910318</b>", "", "", "", "", ""],
        ["1", "−0.4134846441", "−0.4134846441", "−0.4143386084", "−0.4143394918", "−2.1×10<super>−6</super>"],
        ["2", "−0.1035312793", "−0.1035312793", "−0.1035848591", "−0.1035848730", "−1.3×10<super>−7</super>"],
        ["3", "−0.0460271288", "−0.0460271288", "−0.0460377201", "−0.0460377213", "−2.6×10<super>−8</super>"],
        ["4", "−0.0258928660", "−0.0258928660", "−0.0258962180", "−0.0258962182", "−8.4×10<super>−9</super>"],
        ["5", "−0.0165722065", "−0.0165722065", "−0.0165735796", "−0.0165735797", "−3.6×10<super>−9</super>"],
    ], [0.14, 0.18, 0.18, 0.18, 0.18, 0.14])
    fig("fig1_response",
        "Figure 3. Frequency response of Kron's network: |Z<sub>in</sub>| from the kronspice AC sweep, plotted "
        "against x = 1/√(−2E), so a level n of charge product Z sits at x = n/Z. Blue: hydrogen network "
        "(q = e, Z = 1); orange: light atom (q = 0.910 e, Z = 0.910). Vertical lines are the poles pinned by "
        "bisection.")
    r.p("The lines are differences of the extrapolated levels (non-relativistic Coulomb levels depend on n alone, "
        "so the s levels also fix the Lyman lines). Measured values are NIST vacuum wavelengths:")
    r.table([
        ["Line", "Measured (nm)", "q = e, from the poles (nm)", "q = e − measured", "q = 0.910 e (nm)"],
        ["Lyman α", "121.57", "121.5689", "−0.0011 nm", "146.7019"],
        ["Lyman β", "102.57", "102.5737", "+0.0037 nm", "123.7797"],
        ["Lyman γ", "97.25", "97.2551", "+0.0051 nm", "117.3615"],
        ["Lyman δ", "94.97", "94.9757", "+0.0057 nm", "114.6108"],
        ["Balmer α (Hα)", "656.46", "656.4698", "+0.0098 nm", "792.1882"],
        ["Balmer β (Hβ)", "486.27", "486.2739", "+0.0039 nm", "586.8061"],
        ["Balmer γ (Hγ)", "434.17", "434.1731", "+0.0031 nm", "523.9340"],
    ], [0.2, 0.17, 0.25, 0.18, 0.2])
    r.p("1S–2S from the poles: 2466.028388 THz with q = e and 2043.548496 THz with q = 0.910 e, against the "
        "measured 2466.061413 THz (Parthey et al. 2011, ±10 Hz). With q = e the model is 13.4 ppm low, of which "
        "about 9 ppm is fine structure and Lamb shift (the Dirac and QED part) and the rest the extrapolation.")
    fig("fig2_spectrum",
        "Figure 4. Lyman and Balmer lines read from the poles (h → 0), against measured hydrogen (NIST, black). "
        "Blue (q = e) lies on black to the 0.01 nm shown; orange (q = 0.910 e) is the same pattern stretched by "
        "1/Z² = 1.207.")
    fig("fig3_frequencies",
        "Figure 5. Where the frequencies are, on one logarithmic axis: the network poles in circuit units "
        "(0.0491 Hz), the measured 21 cm, Hα and 1S–2S lines, the electron ring m<sub>e</sub>c²/h "
        "(1.24×10<super>20</super> Hz), the knot squashing tone (4.2×10<super>20</super> Hz) and the proton "
        "light ring m<sub>p</sub>c²/h (2.27×10<super>23</super> Hz). The atom's lines come from the Coulomb "
        "binding, not from the parts' internal oscillation; the poles are an analog carrying the levels through "
        "E = 1/ω² − V0.")

    r.h2("Spectrum, neutrality and proton size, model by model")
    r.p("The levels are Bohr's with reduced mass, scaled by Z². The last column is the proton-size shift of "
        "1S–2S, (2/3n³)(Zα)<super>4</super>m<sub>r</sub>c²(r/ƛ<sub>r</sub>)², relative to the measured radius.")
    r.table([
        ["Electron + proton", "Net charge", "1S binding", "Lyman α", "1S–2S size shift"],
        ["toroidal photon (H23) + light ring (H31)", "+0.0897 e", "11.2686 eV", "146.701 nm", "−0.926 MHz"],
        ["toroidal photon (H23) + quark rings (H31)", "+0.0897 e", "11.2686 eV", "146.701 nm", "−0.535 MHz"],
        ["q set to e + light ring (H31)", "0", "13.5983 eV", "121.568 nm", "−0.907 MHz"],
        ["q set to e + quark rings (H31)", "0", "13.5983 eV", "121.568 nm", "−0.338 MHz"],
        ["measured (Dirac + QED) + measured", "0", "13.5983 eV", "121.568 nm", "0"],
        ["<b>Measured hydrogen</b>", "|q| &lt; 10<super>−21</super> e", "13.5984 eV", "121.5668 nm",
         "1S–2S known to 10 Hz"],
    ], [0.36, 0.15, 0.14, 0.15, 0.20], bold_rows=(6,))
    r.small("Measured: neutrality, Bressi et al., <i>Phys. Rev. A</i> 83, 052101 (2011); ionization and Lyman α, "
            "NIST ASD; 1S–2S, Parthey et al., <i>Phys. Rev. Lett.</i> 107, 203001 (2011). The 23.0 GHz between "
            "the reduced-mass Bohr 1S–2S (2466.038424 THz) and the measured value is the Dirac and QED (Lamb) "
            "part, which no row models.")

    r.h2("The 21 cm line")
    r.table([
        ["Electron + proton", "μ<sub>e</sub> (μ<sub>B</sub>)", "μ<sub>p</sub> (μ<sub>N</sub>)",
         "kronspice zz energy", "ν (MHz)", "/ measured"],
        ["toroidal photon + light ring", "0.911194", "1.0000", "5.785×10<super>−26</super> J", "349.218", "0.24586"],
        ["toroidal photon + quark rings", "0.911194", "2.7928", "1.616×10<super>−25</super> J", "975.311", "0.68664"],
        ["toroidal photon + Skyrme knot", "0.911194", "1.8800", "1.088×10<super>−25</super> J", "656.529", "0.46221"],
        ["q set to e + light ring", "1.001161", "1.0000", "8.425×10<super>−26</super> J", "508.623", "0.35808"],
        ["q set to e + quark rings", "1.001161", "2.7928", "2.353×10<super>−25</super> J", "1420.506", "1.00007"],
        ["q set to e + Skyrme knot", "1.001161", "1.8800", "1.584×10<super>−25</super> J", "956.211", "0.67320"],
        ["measured (Dirac + QED) + measured", "1.001160", "2.7928", "2.353×10<super>−25</super> J", "1420.503",
         "1.00007"],
        ["<b>Measured</b> (Hellwig et al. 1970)", "", "", "", "1420.405752", "1"],
    ], [0.3, 0.12, 0.12, 0.18, 0.13, 0.15], bold_rows=(8,))
    r.p("The 0.0069 % left on the standard row is the Breit, recoil and proton-structure terms beyond Fermi's "
        "formula. The quark-ring proton matches only because its moment was fitted to the measured one.")
    fig("fig4_hyperfine",
        "Figure 6. The 21 cm line from the current-loop circuit, model by model, as ν / ν<sub>measured</sub> "
        "(measured 1420.406 MHz; kronspice loop interaction × 4 / h). Orange: toroidal-photon charge 0.910 e; "
        "blue: charge set to e; black: measured particles.")
    fig("fig7_atom3d",
        "Figure 7. Light hydrogen in 3-D. Left, at the atom's scale: the 1s cloud (sampled from |ψ|²) with the "
        "21 cm spin-current loops round the proton (scale bar 1 a0 = 52 918 fm). Right, at the particles' scale: "
        "the electron's double loop of light (ring 193 fm) and the proton ring (0.21 fm, a dot); the whole right "
        "panel is 0.006 a0 wide, a point in the left panel.")

    r.h2("What would hold the electron's light")
    r.bullets([
        "<b>Scale:</b> the knot's core radius 1.2ℓ set to the ring's 193.1 fm gives ℓ = 160.9 fm.",
        "<b>Energy:</b> 547.5√(κ<sub>2</sub>κ<sub>4</sub>) set to m<sub>e</sub>c² gives κ<sub>2</sub> = "
        "0.0009294 J/m and κ<sub>4</sub> = 2.406×10<super>−29</super> J·m.",
        "<b>Prediction:</b> the squashing tone ω = 1.416 c/ℓ gives ħω = 1.737 MeV = 3.40 m<sub>e</sub>c², an "
        "excited electron at 2.248 MeV.",
        "<b>Against measurement:</b> collider searches exclude excited electrons to the TeV scale (PDG 2024); "
        "g − 2 limits the electron's size to about 10<super>−18</super> m against the knot's "
        "10<super>−13</super> m; the knot carries no charge and no spin ½, which still come from H23's assumption.",
    ])

    r.h2("Hydrogen under pressure and the 2p<sub>z</sub> state")
    r.p("The network's outer short circuit is a hard wall; moving it in compresses the atom (Michels, de Boer "
        "&amp; Bijl, <i>Physica</i> 4, 981, 1937). This is Walter Russell's \"high pressure, short waves\" read as "
        "the confined atom (doc/russell-octave-report.md).")
    r.table([
        ["Wall (a0)", "E<sub>1s</sub> (Ha)", "E<sub>2s</sub> (Ha)", "E<sub>2p</sub> (Ha)", "Lyman α (nm)",
         "Pressure"],
        ["60.000", "−0.500000", "−0.125000", "−0.125000", "121.57", "≈ 0"],
        ["8.000", "−0.499975", "−0.084739", "−0.104450", "115.26", "0.00157 GPa"],
        ["5.000", "−0.496417", "+0.141254", "+0.007594", "90.45", "0.525 GPa"],
        ["3.000", "−0.423967", "+1.111685", "+0.481250", "50.36", "30.1 GPa"],
        ["2.500", "−0.334910", "+1.865480", "+0.851978", "38.41", "97.7 GPa"],
        ["2.000", "−0.125000", "+3.327509", "+1.576019", "26.80", "375 GPa"],
        ["1.835", "+0.000219", "+4.113959", "+1.967265", "23.18", "618 GPa"],
        ["1.500", "+0.437018", "+6.644122", "+3.231051", "16.32", "1.93×10<super>3</super> GPa"],
    ], [0.14, 0.15, 0.15, 0.15, 0.17, 0.24])
    r.small("Selected rows of the 16 in run_output.txt (part 7).")
    r.bullets([
        "Compression raises every level and shortens every line; the free atom's 2s = 2p splits (2p stays lower).",
        "At a 2 a0 wall the ground state is exactly −1/8 Ha and Lyman α is 26.80 nm (375 GPa; the pressure at "
        "the centre of the Earth is about 364 GPa, PREM). At 1.835 a0 the ground state reaches zero: the atom "
        "is ionized by pressure.",
        "With the centrifugal term l(l+1)/2r² in the shunts the network gives E<sub>2p</sub> = −0.125000000 Ha "
        "(h → 0), exact −1/8; its kronspice twin (h = 0.1 a0, 399 nodes) gives −0.125026081484 Ha against the "
        "matrix −0.125026081484 (difference 1.8×10<super>−14</super>). The density r²e<super>−r</super>cos²θ "
        "has two lobes on the axis and is zero on the central plane.",
    ])
    fig("fig8_pressure",
        "Figure 8. Hydrogen compressed by a wall in Kron's network. Top: the 1s, 2s and 2p levels against the "
        "wall radius r<sub>c</sub> (E<sub>1s</sub> = −1/8 exactly at 2.0 a0; ionized at 1.835 a0). Bottom: Lyman "
        "α against the pressure on the wall, P = −dE<sub>1s</sub>/dV, from 121.57 nm for the free atom; the "
        "vertical line is the Earth's centre, 364 GPa.", max_height_mm=150)
    fig("fig9_orbitals",
        "Figure 9. Electron density in a plane through the axis, each panel scaled to its own peak: 1s "
        "(E = −1/2 Ha, one sphere) and 2p<sub>z</sub> (E = −1/8 Ha, two lobes, nodal plane z = 0). Kron's "
        "network gives both levels.")

    r.h1("Validation")
    r.table([
        ["Test", "Checks", "Tolerance"],
        ["TestNetworkCoulomb", "network 1s and 2s against −Z²/2 and −Z²/8 for Z = 1 and 0.910",
         "10<super>−7</super>"],
        ["TestNetworkTwin", "kronspice resonance against the network matrix", "10<super>−9</super> Ha"],
        ["TestToroidalElectron", "charge 0.910318 e, g from eq. 12, ring radius 193.08 fm", "10<super>−6</super>"],
        ["TestStandardHydrogen", "ionization, Lyman α and the 21 cm line of the standard atom against measurement",
         "10<super>−4</super>"],
        ["TestHyperfineCircuit", "loop interaction = Fermi's contact energy; ledger = ½ΣLI² + ΣMII; 4U/h = "
         "Fermi's formula", "10<super>−4</super>, 10<super>−12</super>"],
        ["TestLightAtom", "net charge 0.09 e, far above 10<super>−21</super> e; levels scale as Z²", "10<super>−12</super>"],
        ["TestKnotExcitation", "knot energy = m<sub>e</sub>c²; squashing tone 1.737 MeV", "0.002 MeV"],
        ["TestSpectrum", "poles = matrix (10<super>−9</super>); levels → −Z²/2n² (5×10<super>−6</super>); lines "
         "within 0.01 nm of measured; light-atom lines 1/Z² longer (10<super>−5</super>)", "as stated"],
        ["TestConfined", "E<sub>1s</sub> = −1/8 at a 2 a0 wall; zero crossing in (1.83, 1.84) a0; 2p below 2s "
         "under compression; free 2s = 2p = −1/8", "10<super>−6</super>"],
        ["TestTwin2p", "the l = 1 network in kronspice = its matrix", "10<super>−9</super>"],
        ["TestExport", "the radial network and the hyperfine loops export, and their results sheets solve", "—"],
    ], [0.24, 0.6, 0.16])
    r.p("Energy ledger of the 21 cm circuit: stored energy against ½ΣLI² + ΣMII, 4.4×10<super>−16</super>; "
        "ledger residual 0.0.")

    r.h1("Limits of the model")
    r.bullets([
        "The network's poles are an analog: its frequencies stand for energies through E = 1/ω² − V0; they are "
        "not the light the atom emits.",
        "The levels are non-relativistic Bohr levels with reduced mass: fine structure and the Lamb shift (about "
        "9 ppm of 1S–2S, the 23.0 GHz) are not modeled; the 21 cm circuit is Fermi's contact term, missing the "
        "Breit, recoil and proton-structure terms (0.0069 %).",
        "The 21 cm line is not a resonance of the network; it is a magnetic interaction energy from a separate "
        "circuit. The rings' own frequencies (1.24×10<super>20</super> and 2.27×10<super>23</super> Hz) are not "
        "in the spectrum.",
        "The knot's constants κ<sub>2</sub>, κ<sub>4</sub> are fixed here by the electron's mass and ring size, "
        "not by an independent measurement; the Skyrme proton has no computed charge radius.",
    ])

    r.h1("Status")
    r.p("H32 is <i>bounded</i> in the register. With the toroidal-photon charge the atom is charged (0.090 e, "
        "excluded at 10<super>−21</super> e), every line moves by Z² = 0.829, and the 21 cm line falls to 0.69 of "
        "the proton factor. With the charge set to e the levels are Bohr's; the light-ring proton misses 21 cm by "
        "its moment (1 μ<sub>N</sub> against 2.793) and 1S–2S by its size; the quark-ring proton, fitted to the "
        "moment, gets 21 cm and misses 1S–2S by 0.34 MHz because it is smaller (0.678 fm against 0.841). The knot "
        "that would hold the electron's light predicts an excited electron at 2.25 MeV, not seen. The standard "
        "prediction (Dirac electron with QED and the measured proton) gives every measured column.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/light_hydrogen       # about 40 seconds\n"
           "go run ./examples/light_hydrogen -svg examples/light_hydrogen/viz   # also the figures and the 3-D scene\n"
           "go test ./examples/light_hydrogen\n"
           "KRONSPICE_EXPORT=1 go test ./examples/light_hydrogen -run TestExport\n"
           "python examples/light_hydrogen/viz/make_report.py examples/light_hydrogen/light_hydrogen_report.pdf")
    r.small("The interactive 3-D view is viz/atom3d.html. Sources: Williamson &amp; van der Mark, <i>Ann. Fond. "
            "Louis de Broglie</i> 22, 133 (1997); G. Kron, <i>Phys. Rev.</i> 67, 39 (1945); E. Fermi, <i>Z. "
            "Phys.</i> 60, 320 (1930); Eides, Grotch &amp; Shelyuto, <i>Phys. Rep.</i> 342, 63 (2001); Adkins, "
            "Nappi &amp; Witten, <i>Nucl. Phys. B</i> 228, 552 (1983); Faddeev &amp; Niemi, <i>Nature</i> 387, 58 "
            "(1997); CODATA 2022; PDG 2024; NIST ASD.")
    r.build(out)
    tmpdir.cleanup()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "light_hydrogen_report.pdf"))
