# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Hydrogen poles (register entry H34) as one PDF report.

    python examples/hydrogen_poles/viz/make_report.py [out.pdf]

Covers the example's three HTML reports, one section each (left unchanged):
hydrogen_atom_report.html (Parts 1–2), h2_report.html (Part 3) and
hydrogen_spectrum_report.html (Part 4). Numbers come from run_output.txt (the
recorded `go run ./examples/hydrogen_poles`), README.md and those reports. The
figures are the SVGs in viz/ (rasterized to PNG in a temporary directory with
PyMuPDF at build time) and the first frame of viz/atom_spheres.gif.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import pymupdf  # noqa: E402

J32 = "3/2"   # j = 3/2 written plainly (the fraction-slash glyph is not in every font)


def svg_png(name, tmp, dpi=200):
    """Rasterize viz/<name>.svg to a PNG in tmp; return its path."""
    doc = pymupdf.open(os.path.join(HERE, name + ".svg"))
    out = os.path.join(tmp, name + ".png")
    doc[0].get_pixmap(dpi=dpi).save(out)
    return out


def main(out):
    tmpdir = tempfile.TemporaryDirectory()
    tmp = tmpdir.name

    r = Report("examples/hydrogen_poles", "Hydrogen poles",
               "Walter Russell's \"four pairs of poles — not one pair\" read as hydrogen: one atom's magnetic poles "
               "level by level, the four magnets of H<sub>2</sub>, and hydrogen's lines as kronspice resonances "
               "(register entry H34)",
               "register entry H34", "recorded run: run_output.txt (about a minute)")

    def fig(name, caption, **kw):
        r.figure(svg_png(name, tmp), caption, **kw)

    r.h1("Summary")
    r.p("Hypothesis H34 (user hypothesis, 2026-11-05) reads Walter Russell's painting \"Analysis of the Octave Bar "
        "Magnet\", beneath which he wrote <i>\"a bar magnet has four pairs of poles — not one pair\"</i>, as "
        "hydrogen, and takes the painting's two touching spheres as two hydrogen atoms (the bar reading is H33, "
        "<i>octave_magnet</i>). The example computes one atom's magnetic field level by level from the electron's "
        "orbital current and spin magnetization, counts its N and S zones as a bench would, builds H<sub>2</sub> "
        "with the Heitler–London wavefunction and its four spin magnets, and finds hydrogen's spectral lines as "
        "resonances of circuits solved by kronspice. This PDF covers the example's three reports: one atom (Parts "
        "1–2), two atoms (Part 3) and the spectrum (Part 4).")
    r.bullets([
        "<b>Ground state: one pair.</b> 1s½ has a dipole of −1.001160 μ<sub>B</sub> (equal to the Landé value) "
        "and no octupole; so does every j = ½ level.",
        "<b>Four pairs first appear in 2p" + J32 + ".</b> The m = ½ state has octupole −54.062621 "
        "μ<sub>B</sub>·a0² and shows four N S N S zones round the atom from 0.5 to 15.4 a0 "
        "(r<sub>c</sub> = 15.61 a0); the octupole comes entirely from the electron's spin magnetization.",
        "<b>Selection rules hold numerically:</b> Wigner–Eckart ratio −3.000000; even-l moments and j = ½ "
        "octupoles vanish to round-off; series field = direct sum of 1.28 million current loops to "
        "2.4×10<super>−4</super>.",
        "<b>H<sub>2</sub>: the bond cancels the electron magnets.</b> The singlet binds (R<sub>e</sub> 1.6425 a0, "
        "D<sub>e</sub> 3.1557 eV; measured 1.4011 a0, 4.748 eV); the triplet repels at every R. The bound molecule "
        "has zero spin density, leaving only the proton dipoles (1.52×10<super>−3</super> μ<sub>B</sub> each).",
        "<b>Spectrum as resonances:</b> Lyman, Balmer, Paschen and Brackett α within 16 ppm, the 2p fine structure "
        "10968.90 MHz against 10969.10 MHz, 21 cm +56 ppm and ESR +18 ppm from the model; H<sub>2</sub>/HD/"
        "D<sub>2</sub> vibration 13.5 % low with isotope ratios to 0.07 %.",
        "<b>Status: bounded.</b> No measurement of an excited hydrogen atom's magnetic octupole, or of "
        "H<sub>2</sub>'s nuclear four-zone field, is known.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H34).</b> Russell's \"four pairs of poles — not one pair\" describes hydrogen, and his two "
        "touching spheres are two hydrogen atoms. <b>Standard prediction.</b> Dirac/QED hydrogen with these "
        "multipoles: ground-state hydrogen has one pair of poles (a dipole); four pairs (an octupole) first appear "
        "in the 2p" + J32 + " level. In H<sub>2</sub> the bond pairs the electron spins, so only the protons' tiny "
        "dipoles remain.")
    r.h2("What \"four pairs\" means for an atom")
    r.bullets([
        "<b>Poles come in powers of two.</b> A dipole has one pair; an octupole (2³ poles) has four. \"Four pairs, "
        "not one pair\" names the octupole against the dipole.",
        "<b>No two-pair state.</b> A state of definite parity has no magnetic quadrupole, so an atom's magnetic "
        "poles come as one pair or four pairs, never two.",
        "<b>Four pairs need higher states.</b> A magnetic 2<super>k</super>-pole needs 2j ≥ k. The ground state "
        "(1s, j = ½) can have only one pair; four pairs first become possible in 2p" + J32 + ".",
        "<b>On a sphere round the atom</b> an axial octupole shows four alternating zones of B<sub>r</sub> (N S N S "
        "from pole to pole); a dipole shows two.",
    ])

    # ---------------- Report 1: one atom ----------------
    r.h1("Report 1 — One hydrogen atom (Parts 1–2)")
    r.h2("The model")
    r.bullets([
        "<b>The electron:</b> Pauli spinors |n l j m&gt; built with Clebsch–Gordan coefficients from the hydrogen "
        "radial functions. Its field has two sources: the orbital current J = −2μ<sub>B</sub> Σ "
        "m<sub>l</sub>|ψ|²/ρ, and the spin magnetization M = −(g<sub>s</sub>/2)μ<sub>B</sub> ψ<super>†</super>σψ "
        "with g<sub>s</sub> = 2.002319.",
        "<b>The field:</b> both enter as one azimuthal current, J<sub>orb</sub> + curl M, whose field is computed "
        "exactly from the Legendre-series Green's function (Jackson §5.5).",
        "<b>The proton:</b> its spin dipole (+1.5210×10<super>−3</super> μ<sub>B</sub>, m<sub>I</sub> = +½) is "
        "added at the centre.",
    ])
    r.h2("The atom's magnetism, level by level")
    r.table([
        ["Level", "∫|ψ|²", "Dipole (μ<sub>B</sub>)", "Landé −g<sub>j</sub>m<sub>j</sub>",
         "Octupole (μ<sub>B</sub>·a0²)", "Pole pattern round the atom"],
        ["1s½, m = ½", "1.000000", "−1.001160", "−1.001160", "0", "2 zones (one pair) at every radius"],
        ["2p½, m = ½", "1.000000", "−0.332947", "−0.332947", "0", "2 zones at every radius"],
        ["2p" + J32 + ", m = " + J32, "1.000000", "−2.001160", "−2.001160", "+18.020874",
         "2 zones at every radius (r<sub>c</sub> = 4.25 a0 falls inside the cloud)"],
        ["2p" + J32 + ", m = ½", "1.000000", "−0.667053", "−0.667053", "−54.062621",
         "<b>4 zones 0.5–15.4 a0</b>, then 2 (r<sub>c</sub> = 15.61 a0)"],
    ], [0.15, 0.1, 0.13, 0.13, 0.15, 0.34])
    r.p("Split by source (dipole, octupole):")
    r.table([
        ["Level", "Orbital current", "Spin magnetization"],
        ["1s½, m = +½", "−0.000000, −0.000000", "−1.001160, +0.000000"],
        ["2p½, m = +½", "−0.666667, −0.000000", "+0.333720, +0.000000"],
        ["2p" + J32 + ", m = +" + J32, "−1.000000, −0.000000", "−1.001160, +18.020874"],
        ["2p" + J32 + ", m = +½", "−0.333333, −0.000000", "−0.333720, −54.062621"],
    ], [0.3, 0.35, 0.35])
    r.p("The octupole comes entirely from the electron's spin magnetization; the l = 1 orbital current carries "
        "none. The j = ½ levels have none in either part, as 2j ≥ 3 requires.")
    fig("figH1_atom_fieldlines",
        "Figure 1 (H1). Field lines (contours of ρA<sub>φ</sub>) of each level in a plane through the spin axis, "
        "over the electron density; the ring at 8 a0 is coloured by the sign of B<sub>r</sub> (orange N, blue S). "
        "Panels: 1s½ m = +½, 2p½ m = +½, 2p" + J32 + " m = +" + J32 + ", 2p" + J32 + " m = +½. Only the last "
        "shows four zones; its dashed circle marks 15.6 a0, where they end.", width_frac=0.85, max_height_mm=170)

    r.h2("Counting the poles")
    r.p("The pole count is read on spheres of growing radius; each sign change of B<sub>r</sub> from pole to pole "
        "separates an N zone from an S zone. Outside the cloud the octupole leads inside r<sub>c</sub> = "
        "√(3m<sub>3</sub>/m<sub>1</sub>), or √(−2m<sub>3</sub>/m<sub>1</sub>) when the ratio is negative; beyond "
        "it the dipole leads, because it falls as 1/r³ and the octupole as 1/r<super>5</super>. For m = " + J32 +
        ", r<sub>c</sub> = 4.25 a0 lies inside the cloud, where the field is not the exterior expansion, and no "
        "four-zone pattern appears at any radius.")
    fig("figH2_sphere_scans",
        "Figure 2 (H2). B<sub>r</sub> from the north pole (0°) to the south pole (180°) on spheres of 3 to 24 a0, "
        "each curve scaled to its own peak. Left, 2p" + J32 + " m = +½: three sign changes (four zones) out to "
        "15.6 a0, then one. Right, m = +" + J32 + ": one sign change at every radius.")
    fig("figH3_zones",
        "Figure 3 (H3). Zones against sphere radius (log scale). Two zones (one pair) for the ground state and "
        "every j = ½ level; four for 2p" + J32 + " m = ½ out to r<sub>c</sub> = 15.6 a0.")
    r.figure(os.path.join(HERE, "atom_spheres.gif"),
             "Figure 4. First frame (radius 1 a0) of the animation viz/atom_spheres.gif: spheres round 2p½ (left), "
             "2p" + J32 + " m = " + J32 + " (middle) and 2p" + J32 + " m = ½ (right), coloured by B<sub>r</sub> "
             "(orange N, blue S). The right-hand sphere carries four bands; in the animation they stay until the "
             "radius passes 15.6 a0 (the orange tick), then become two.", width_frac=0.7)
    r.p("The octave bar's own bench scan (<i>octave_magnet</i>) on the atom: sign changes of B<sub>z</sub> on a "
        "line beside the spin axis, z from −60 to 60 a0. A dipole gives 2; more mean more pole sites along the line.")
    r.table([
        ["Level", "ρ = 2 a0", "ρ = 5 a0", "ρ = 10 a0", "ρ = 20 a0"],
        ["1s½, m = +½", "2", "2", "2", "2"],
        ["2p½, m = +½", "0", "2", "2", "2"],
        ["2p" + J32 + ", m = +" + J32, "0", "0", "2", "2"],
        ["2p" + J32 + ", m = +½", "0", "<b>4</b>", "2", "2"],
    ], [0.28, 0.18, 0.18, 0.18, 0.18])
    fig("figH4_line_scan",
        "Figure 5 (H4). B<sub>z</sub> along a line 5 a0 beside the spin axis, each curve scaled to its own peak. "
        "Sign changes: 1s½ 2, 2p½ 2, 2p" + J32 + " m = +" + J32 + " 0, 2p" + J32 + " m = +½ 4 (two N and two S "
        "sites along the line).")
    r.h2("Checks")
    r.table([
        ["Check", "Result"],
        ["Normalization ∫|ψ|², every level", "1 to 10<super>−8</super>"],
        ["Dipole against the Landé value −g<sub>j</sub>m<sub>j</sub>μ<sub>B</sub>", "equal to 10<super>−8</super>"],
        ["Even-l moments (parity)", "≤ 4.2×10<super>−10</super>"],
        ["Octupole of the j = ½ levels (2j ≥ k)", "≤ 10<super>−11</super>"],
        ["Wigner–Eckart: octupole(m = ½)/octupole(m = " + J32 + "), theory −3 "
         "(proportional to m(5m² − 3j(j+1) + 1))", "−3.000000"],
        ["Series field against a direct sum of 1 280 000 current loops, 4 points outside the cloud",
         "worst 2.4×10<super>−4</super> of |B|"],
    ], [0.7, 0.3])

    # ---------------- Report 2: H2 ----------------
    r.h1("Report 2 — Two hydrogen atoms, four magnets (Part 3)")
    r.p("Russell's two touching spheres read as two hydrogen atoms bonding into H<sub>2</sub>. The molecule holds "
        "four small magnets: two electron spins and two proton spins. The model is the Heitler–London wavefunction "
        "(Heitler &amp; London 1927) with the exact one- and two-centre integrals (Sugiura 1927 for the exchange "
        "integral).")
    r.bullets([
        "<b>Singlet (spins opposite) binds:</b> R<sub>e</sub> 1.6425 a0 (0.8692 Å), D<sub>e</sub> 3.1557 eV, the "
        "classic Heitler–London result; measured 1.4011 a0 and 4.748 eV. The simple wavefunction gets two thirds of "
        "the binding.",
        "<b>Triplet (spins parallel) repels at every R:</b> E − 2E<sub>H</sub> = +10.099 eV at 1.4 a0, +4.189 eV "
        "at 2 a0, +0.2392 eV at 4 a0.",
    ])
    fig("figM1_energies",
        "Figure 6 (M1). Heitler–London energies of two hydrogen atoms against their separation, relative to the "
        "separated atoms. Only the singlet (blue) has a well, at 1.64 a0 and 3.16 eV; the measured bond (1.40 a0, "
        "4.75 eV) is marked. The triplet (orange) repels at every R.")
    fig("figM2_spin_density",
        "Figure 7 (M2). Along the bond axis at R = 1.40 a0: electron density of the singlet (grey), spin density of "
        "the triplet (orange, cusps at the protons) and of the singlet (blue, zero everywhere, because its spin part "
        "gives &lt;σ<sub>z</sub>&gt; = 0 for each electron).")
    r.p("The field of the four magnets, at the measured bond length 1.4011 a0 with spins quantized along the bond. "
        "The electron field comes from the triplet spin magnetization M<sub>z</sub> = −(g<sub>s</sub>/2)μ<sub>B</sub> "
        "s<sub>z</sub>, summed as current loops and evaluated inside the cloud on cell corners, where the "
        "neighbouring loops' near fields cancel (grid h = 0.05 against 0.025: worst 3.4×10<super>−4</super> of |B|; "
        "h = 0.1 against 0.05: 1.4×10<super>−3</super>). The proton spins are point dipoles at the nuclei.")
    r.table([
        ["Spin arrangement", "∫ spin density", "Moment (μ<sub>B</sub>)", "Bound?", "Zones round the molecule",
         "B<sub>z</sub> sign changes at ρ = 2 a0"],
        ["singlet electrons + para protons", "0.0000", "+0.00000", "yes", "none — no magnetic field at all", "—"],
        ["singlet electrons + ortho protons", "0.0000", "+0.00304", "yes",
         "2 zones 0.20–0.29 a0; 4 zones 0.31–1.56 a0; 2 zones 1.62–20.00 a0", "2"],
        ["triplet electrons + para protons", "2.0001", "−2.00232", "no", "2 zones 0.20–20.00 a0", "2"],
        ["triplet electrons + ortho protons (all four aligned)", "2.0001", "−1.99928", "no",
         "2 zones 0.20–0.59 a0; 4 zones 0.61–0.77 a0; 2 zones 0.81–20.00 a0", "2"],
    ], [0.25, 0.1, 0.1, 0.1, 0.3, 0.15])
    fig("figM3_h2_fieldlines",
        "Figure 8 (M3). Field lines (contours of ρA<sub>φ</sub>) in a plane through the bond (bond vertical), over "
        "the spin density; ring at 1 a0 coloured by the sign of B<sub>r</sub>. Panels: singlet + ortho protons "
        "(bound; electrons cancel, two proton dipoles of 0.0015 μ<sub>B</sub>), triplet + para protons (unbound; "
        "−2.0 μ<sub>B</sub> in all), triplet + ortho (all four aligned, unbound). Long arrows are the electron "
        "moments, short arrows the protons' (drawn far larger than their 1/658 share).")
    fig("figM4_h2_zones",
        "Figure 9 (M4). Zones on spheres centred between the protons (each 0.70 a0 from the centre), against "
        "radius. Ortho protons alone make four zones near the molecule (two dipoles end to end); the triplet "
        "electrons make one pair at almost every radius. Singlet + para has no field and no zones and is left off.")
    r.p("The two-atom picture does hold four small magnets, and with all four aligned the molecule would carry about "
        "2 μ<sub>B</sub>. But that needs the electron spins parallel, and the triplet does not bind. The bond forms "
        "only with the spins opposite, and then the electron magnets cancel exactly, at every point. What is left "
        "are the proton magnets, 658 times weaker. Ortho-H<sub>2</sub> shows a four-zone pattern near the protons; "
        "para-H<sub>2</sub>, the commonest form at low temperature, has no magnetic field at all. The rotating "
        "molecule's own small moment (ortho, J = 1) is not modelled.")

    # ---------------- Report 3: spectrum ----------------
    r.h1("Report 3 — Hydrogen's spectral lines as circuit resonances (Part 4)")
    r.bullets([
        "<b>Atomic levels</b> are resonances of Kron's 1945 radial network (<i>Phys. Rev.</i> 67, 39), one network "
        "per orbital (1s, 2p, 3d, 4f, 5g), solved by kronspice's AC analysis in its dual form as in "
        "<i>light_hydrogen</i>: each level is where the driving-point reactance changes sign. Grids of h = 0.1 and "
        "0.05 a0 are Richardson-extrapolated, and energies are scaled by the reduced mass μ/m = 0.999456.",
        "<b>Fine structure:</b> the p networks carry the spin–orbit term (g<sub>s</sub> − 1)(α²/2)&lt;L·S&gt;/r³ "
        "(Bethe &amp; Salpeter §12), so the j = ½ and j = " + J32 + " networks resonate apart.",
        "<b>H<sub>2</sub> vibration:</b> an LC circuit, the masses as capacitors and the Heitler–London spring "
        "(k = 0.2759 E<sub>h</sub>/a0² at 1.6425 a0) as an inductor 1/k; its single transfer pole is the vibration.",
        "<b>Two lines come from the model, not a circuit:</b> the 21 cm line uses the computed ψ(0) of Part 1 "
        "(Fermi contact), and the ESR frequency the computed dipole.",
    ])
    r.table([
        ["Level", "kronspice h", "kronspice h/2", "Richardson", "Exact", "AC solves"],
        ["1s", "−0.4987562112", "−0.4996878900", "−0.4999984496", "−0.500000", "112"],
        ["2p (no spin–orbit)", "−0.1250260815", "−0.1250065129", "−0.1249999900", "−0.125000", "108"],
        ["2p½", "−0.1250271946", "−0.1250076252", "−0.1250011021", "−0.125000", "108"],
        ["2p" + J32, "−0.1250255250", "−0.1250059568", "−0.1249994341", "−0.125000", "108"],
        ["3d", "−0.0555565843", "−0.0555558128", "−0.0555555556", "−0.055556", "104"],
        ["4f", "−0.0312501395", "−0.0312500349", "−0.0312500000", "−0.031250", "104"],
        ["5g", "−0.0200000310", "−0.0200000072", "−0.0199999993", "−0.020000", "102"],
    ], [0.2, 0.17, 0.17, 0.17, 0.15, 0.14])
    r.p("Energies in hartree. Twin: kronspice resonance against the network matrix (Sturm), h/2: worst "
        "4.7×10<super>−14</super> Ha. The fine structure 2p" + J32 + " − 2p½ from the two networks' resonances is "
        "10979.6 MHz at h, 10971.6 MHz at h/2 and 10968.9 MHz by Richardson; with the spin–orbit term off the two "
        "are one network and one resonance.")
    r.table([
        ["Line", "Found", "Measured", "Difference", "How"],
        ["Lyman α (2p → 1s)", "121.5689 nm", "121.5670 nm", "−0.0016 %", "kronspice resonance"],
        ["Balmer α, Hα (3d → 2p)", "656.4697 nm", "656.4610 nm", "−0.0013 %", "kronspice resonance"],
        ["Paschen α (4f → 3d)", "1875.6274 nm", "1875.6100 nm", "−0.0009 %", "kronspice resonance"],
        ["Brackett α (5g → 4f)", "4052.2813 nm", "4052.2700 nm", "−0.0003 %", "kronspice resonance"],
        ["fine structure 2p" + J32 + " − 2p½", "10968.90 MHz", "10969.10 MHz", "−0.0018 %", "kronspice resonance"],
        ["hyperfine 1s, 21 cm (F = 1 → 0)", "1420.4855 MHz", "1420.4058 MHz", "+0.0056 %",
         "model (Fermi contact, ψ(0) of Part 1)"],
        ["electron spin resonance of H, per tesla", "28024.95 MHz", "28024.45 MHz", "+0.0018 %",
         "model (dipole of Part 1)"],
    ], [0.3, 0.15, 0.15, 0.13, 0.27])
    r.p("The remainders are physics these models leave out: for the optical lines the Dirac level shifts and the "
        "Lamb shift; for 21 cm the Breit, recoil and proton-structure terms; for ESR the bound g-factor correction "
        "−α²/3. Measured values: NIST ASD (vacuum wavelengths), Hellwig et al. (1970), Tiedemann &amp; Robinson "
        "(1977), and the established ≈ 10 969 MHz fine-structure interval (e.g. Hagley &amp; Pipkin 1994).")
    fig("figS2_lines",
        "Figure 10 (S2). Each line's difference from its measured value in parts per million, placed by frequency "
        "from the radio to the ultraviolet. Blue: kronspice circuit resonances; orange: this example's model "
        "(Part 1). The ESR point is the frequency per tesla, plotted at 1 T. All lie within 60 ppm.")
    fig("figS1_fine_structure",
        "Figure 11 (S1). The fine structure as two kronspice resonances: Im(1/V) at the drive of Kron's 2p network "
        "(h = 0.05 a0) against the frequency offset from the 2p level without spin–orbit. The j = ½ and j = " + J32 +
        " networks cross zero 10971.6 MHz apart at this grid; without spin–orbit the crossing lies between them.")
    r.table([
        ["Molecule", "kronspice ω<sub>e</sub> (cm<super>−1</super>)", "Measured ω<sub>e</sub>", "Difference"],
        ["H<sub>2</sub>", "3805.0", "4401.21", "−13.5 %"],
        ["HD", "3295.5", "3813.15", "−13.6 %"],
        ["D<sub>2</sub>", "2691.2", "3115.50", "−13.6 %"],
    ], [0.25, 0.3, 0.25, 0.2])
    r.bullets([
        "The Heitler–London spring is too soft: every vibration comes out 13.5 % low, the same error that leaves its "
        "bond energy at two thirds of the measured value.",
        "The isotope ratios, set by the masses alone, agree to 0.07 %: HD/H<sub>2</sub> 0.8661 against 0.8664, "
        "D<sub>2</sub>/H<sub>2</sub> 0.7073 against 0.7079; the small remainders are anharmonicity.",
        "Rotation B = 1/(2μR²): 44.30 cm<super>−1</super> at the Heitler–London R<sub>e</sub>, 60.89 at the "
        "measured R<sub>e</sub>, against measured B<sub>e</sub> 60.853. The bond length, not the formula, is off.",
    ])
    fig("figS3_vibration",
        "Figure 12 (S3). The vibration circuits of H<sub>2</sub>, HD and D<sub>2</sub> (damped to Q = 200 for "
        "display): |V| at the second mass against wavenumber. Each peak sits 13.5 % below its measured line "
        "(vertical lines).")
    r.p("None of these lines involves the 2p" + J32 + " octupole of Report 1: its magnetic interactions with the "
        "proton are far below these intervals, and the proton has no octupole to couple to.")

    r.h1("Validation")
    r.table([
        ["Test", "Checks", "Tolerance"],
        ["TestNormAndDipole", "∫|ψ|² = 1; dipole = −g<sub>j</sub>m<sub>j</sub>μ<sub>B</sub>, all four levels",
         "10<super>−6</super>"],
        ["TestOctupoleSelection", "even-l moments zero; j = ½ octupole zero; j = " + J32 + " ratio −3",
         "10<super>−8</super>; 10<super>−6</super>"],
        ["TestSeriesTwin", "series field = direct loop sum", "10<super>−3</super>"],
        ["TestZones", "1s 2 zones at 1–40 a0; 2p" + J32 + " m = ½ 4 zones to 14 a0, 2 beyond 18 a0; "
         "r<sub>c</sub> = 15.61 a0", "exact; 0.05 a0"],
        ["TestE1", "E1(x) against Abramowitz &amp; Stegun Table 5.1", "10<super>−12</super>"],
        ["TestHeitlerLondon", "R<sub>e</sub> = 1.6425 a0, D<sub>e</sub> = 0.1160 E<sub>h</sub>; triplet unbound at "
         "every R", "0.002 a0, 5×10<super>−4</super> E<sub>h</sub>"],
        ["TestSpinDensity", "∫ spin density 0 (singlet), 2 (triplet)", "exact; 10<super>−3</super>"],
        ["TestSpectrum", "kronspice resonance = matrix; Lyman α to 30 ppm; fine structure, 21 cm to 100 ppm; ESR "
         "to 30 ppm; LC pole = √(k/μ); D<sub>2</sub>/H<sub>2</sub> ratio",
         "10<super>−12</super> E<sub>h</sub>; as listed; 10<super>−9</super>; 2×10<super>−3</super>"],
        ["TestMoleculeField", "triplet far field → dipole of −g<sub>s</sub>μ<sub>B</sub>; ortho protons 4 zones at "
         "1 a0; para moment 0", "10/r²; exact"],
    ], [0.22, 0.56, 0.22])

    r.h1("Limits of the model")
    r.bullets([
        "Non-relativistic Pauli spinors and Bohr-type levels with spin–orbit added; Dirac level shifts, the Lamb "
        "shift, Breit, recoil and proton-structure terms and the bound g correction are left out.",
        "H<sub>2</sub> is the Heitler–London wavefunction: its bond is too long and too weak (1.6425 a0, 3.1557 eV "
        "against 1.4011 a0, 4.748 eV) and its spring too soft; the field is evaluated at the measured bond length.",
        "The rotating molecule's own small moment (ortho, J = 1) is not modelled.",
    ])

    r.h1("Status")
    r.p("H34 is <i>bounded</i>. One atom: ground-state hydrogen has one pair, and so does every j = ½ level; four "
        "pairs (an octupole) first appear in 2p" + J32 + ", which with m = ±½ shows four N S N S zones round the "
        "atom out to 15.6 a0, from the electron's spin magnetization. Two atoms: the bond forms only with the "
        "electron spins opposite, so the electron magnets cancel exactly; the proton magnets, 658 times weaker, "
        "remain, with a four-zone pattern for ortho-H<sub>2</sub> and no field for para-H<sub>2</sub>. What would "
        "move it: an observation of an excited atom's magnetic octupole, or of H<sub>2</sub>'s nuclear four-zone "
        "field; neither has been measured, to our knowledge.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/hydrogen_poles                # about a minute\n"
           "go run ./examples/hydrogen_poles -svg viz       # also the figures and animations\n"
           "go test ./examples/hydrogen_poles\n"
           "python examples/hydrogen_poles/viz/make_report.py examples/hydrogen_poles/hydrogen_poles_report.pdf")
    r.small("Sources: W. Russell, \"Analysis of the Octave Bar Magnet\" (painting; doc/russell-octave-report.md); "
            "J. D. Jackson, <i>Classical Electrodynamics</i>, 3rd ed., §5.5–5.6; J. J. Sakurai, <i>Modern Quantum "
            "Mechanics</i>, §3.8, §3.11; Heitler &amp; London, <i>Z. Phys.</i> 44, 455 (1927); Sugiura, <i>Z. "
            "Phys.</i> 45, 484 (1927); Slater (1963), App. 6; Kołos &amp; Wolniewicz, <i>J. Chem. Phys.</i> 49, 404 "
            "(1968); G. Kron, <i>Phys. Rev.</i> 67, 39 (1945); Bethe &amp; Salpeter (1957); NIST ASD; Hellwig et "
            "al. (1970); Tiedemann &amp; Robinson (1977); Hagley &amp; Pipkin (1994); Huber &amp; Herzberg (1979); "
            "Abramowitz &amp; Stegun, Table 5.1; CODATA 2018.")
    r.build(out)
    tmpdir.cleanup()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "hydrogen_poles_report.pdf"))
