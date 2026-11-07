# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Octave magnet (register entry H33) as one PDF report.

    python examples/octave_magnet/viz/make_report.py [out.pdf]

Covers the example's two HTML reports, one section each (left unchanged):
report.html (the bar, Parts 1–3) and ring_report.html (the ring of twelve,
Part 5). Numbers come from run_output.txt (the recorded
`go run ./examples/octave_magnet`), README.md, those reports and
doc/russell-octave-report.md. Figures are the SVGs in viz/ (rasterized to PNG in
a temporary directory with PyMuPDF at build time; CSS class styles are inlined
first, since MuPDF ignores them) and the first frames of the two ring GIFs.
"""
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import pymupdf  # noqa: E402


def svg_png(name, tmp, dpi=200):
    """Rasterize viz/<name>.svg to a PNG in tmp; return its path.

    Simple class rules (.name{...}) in a <style> block are copied onto the
    elements as style attributes, because MuPDF does not apply CSS classes."""
    with open(os.path.join(HERE, name + ".svg"), encoding="utf-8") as f:
        src = f.read()
    for cls, body in re.findall(r"\.([A-Za-z_][\w-]*)\{([^}]*)\}", "".join(re.findall(r"<style>(.*?)</style>",
                                                                                         src, re.S))):
        src = src.replace(f'class="{cls}"', f'style="{body}"')
    doc = pymupdf.open(stream=src.encode("utf-8"), filetype="svg")
    out = os.path.join(tmp, name + ".png")
    doc[0].get_pixmap(dpi=dpi).save(out)
    return out


def main(out):
    tmpdir = tempfile.TemporaryDirectory()
    tmp = tmpdir.name

    r = Report("examples/octave_magnet", "Octave magnet",
               "Does a bar magnet have four pairs of poles, as Walter Russell's painting says? The exact field of a "
               "bar, the bench test that decides, and the pattern closed into a ring of twelve sections "
               "(register entry H33)",
               "register entry H33", "recorded run: run_output.txt (a few seconds)")

    def fig(name, caption, **kw):
        r.figure(svg_png(name, tmp), caption, **kw)

    r.h1("Summary")
    r.p("Hypothesis H33 (user hypothesis, 2026-11-05) comes from Walter Russell's painting \"Analysis of the Octave "
        "Bar Magnet\", which says a bar magnet has <i>\"four pairs of poles — not one pair\"</i>, drawn as a chain of "
        "lens-shaped rings strung along the magnet's axis. The example builds the standard field of a bar magnet "
        "exactly, shows what each claim looks like as a field, and fixes the bench measurement that decides between "
        "them. A second part bends the pattern into a ring magnet of twelve 30° sections. This PDF covers the "
        "example's two reports: the bar (Parts 1–3) and the ring (Part 5).")
    r.bullets([
        "<b>The model is exact to round-off:</b> the current-sheet field equals the solenoid formula on the axis to "
        "1.3×10<super>−6</super> of μ<sub>0</sub>M<sub>0</sub>; the loop field equals Biot–Savart off the axis to "
        "2.160×10<super>−6</super>; a kronspice search coil reads the same flux as ΣM·I and 2πψ to six digits.",
        "<b>One pair against four:</b> a uniform bar has 2 magnetic-charge sites and gives 2 sign changes of "
        "B<sub>z</sub> 3 mm beside it; four reversed segments (S N | N S | S N | N S) have 5 sites and give 5.",
        "<b>The lenses are real field lines</b> of a chain of poles: the four-segment bar draws the chain of lenses "
        "in Russell's painting; the uniform bar's outside lines run end to end in one family.",
        "<b>The decisive test is a Hall-probe scan</b> of an ordinary bar: two sign changes support the standard "
        "prediction, five or more on an unmodified bar support Russell. The model cannot decide which a real bar is.",
        "<b>The ring of twelve works when built that way:</b> both alternating rings (radial, or along the ring with "
        "like poles meeting at the joints) show 12 sign changes round the outside (peaks 291 and 315 mT); the "
        "all-outward control shows 0.",
        "<b>Status: untested</b> for the bar, until the measurement is made.",
    ])

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H33).</b> A bar magnet has four pairs of poles, not one; its field is a chain of lens-shaped "
        "rings on the axis. <b>Standard prediction.</b> A uniformly magnetized bar has one pair of poles (magnetic "
        "charge on the end faces only) and its outside field lines run end to end in one family. Four pairs, as a "
        "field, need three reversals of M inside the bar (five charge sites end to end, eight if four magnets are "
        "separated) and give at least five sign changes of B<sub>z</sub> beside it.")
    r.small("The painting itself is described in doc/russell-octave-report.md; it is not reproduced here. That note "
            "also records where the painting disagrees with measurement: Russell's structures repeat in octaves "
            "(frequency ratio 2), while hydrogen's lines follow 1/n² (Lyman α to the Lyman limit is 4/3, Hα to Lyman "
            "α is 5.4).")

    # ---------------- Report 1: the bar ----------------
    r.h1("Report 1 — The bar magnet (Parts 1–3)")
    r.h2("The model")
    r.bullets([
        "<b>The bar as currents.</b> A magnetization M(z)ẑ in a rod of radius a is exactly a sheet of azimuthal "
        "current M(z) on its side (J = curl M, K = M×n, n the outward unit normal; Jackson §5.8–5.9), built as 2,000 coaxial loops each "
        "carrying M·dz.",
        "<b>Poles</b> are the magnetic charge −div M: the end faces, plus every place inside where M changes.",
        "<b>Field lines</b> are the contours of the flux function ψ = ρA<sub>φ</sub>, summed from the exact "
        "elliptic-integral potential of a loop (Smythe §7.10).",
        "<b>The bar:</b> an alnico rod, 12.7 mm × 76.2 mm, B<sub>r</sub> 1.25 T (M<sub>0</sub> = "
        "B<sub>r</sub>/μ<sub>0</sub> = 9.947×10<super>5</super> A/m). Demagnetization is ignored, declared.",
    ])
    r.h2("Checks")
    r.table([
        ["Check", "Result"],
        ["Sheet field on the axis (2000 loops) against the exact solenoid formula",
         "worst 1.3×10<super>−6</super> of μ<sub>0</sub>M<sub>0</sub>"],
        ["One loop off the axis, elliptic formula against Biot–Savart (em.Loop, 2000 segments)",
         "B<sub>ρ</sub> 2.160×10<super>−6</super>, B<sub>z</sub> 3.097×10<super>−7</super>"],
        ["kronspice search coil: flux ∫V dt against ΣM·I and 2πψ", "equal to six digits"],
    ], [0.62, 0.38])
    r.h2("Poles as magnetic charge")
    r.table([
        ["Bar", "Charge sites (position, charge)"],
        ["uniform (one pair)", "2: −38.1 mm (−126 A·m), +38.1 mm (+126 A·m)"],
        ["four reversed segments (four pairs)", "5: −38.1 mm (−126 A·m), −19.1 mm (+252 A·m), +0.0 mm (−252 A·m), "
         "+19.1 mm (+252 A·m), +38.1 mm (−126 A·m)"],
    ], [0.32, 0.68])
    r.p("The uniform bar carries charge only on its two end faces: one pair. Four pairs need M to reverse three "
        "times inside; then S N | N S | S N | N S gives the two ends plus three doubled \"consequent\" poles, four "
        "dipoles end to end.")
    fig("figM1_fieldlines",
        "Figure 1 (M1). Field lines (contours of ψ = ρA<sub>φ</sub>) in a plane through the axis (axis horizontal, "
        "as in Russell's drawing); the bar is the grey outline and dots mark the magnetic charge (orange N, black S). "
        "Top: uniform bar, charge on the two end faces, outside lines end to end. Bottom: four segments S N | N S | "
        "S N | N S, the lines bow from pole to pole along the axis, a chain of lenses.",
        width_frac=0.85, max_height_mm=170)
    r.p("A field line that leaves the axis at one pole and returns at the next is a lens strung on the axis, so a "
        "chain of poles draws a chain of lenses, as in Russell's painting.")
    r.h2("The bench test")
    r.table([
        ["Bar", "Charge sites", "B<sub>z</sub> sign changes on the axis", "B<sub>z</sub> sign changes 3 mm beside it"],
        ["uniform (an ordinary bar): one pair", "2", "0", "2"],
        ["four segments S N | N S | S N | N S: four pairs end to end", "5", "3", "5"],
    ], [0.4, 0.15, 0.22, 0.23])
    fig("figM2_scan",
        "Figure 2 (M2). B<sub>z</sub> along a line 3 mm beside the bar (a Hall probe run parallel to the axis); "
        "the shaded band is the bar. Uniform bar: 2 sign changes. Four segments: 5. One sign change per pole site.")
    r.p("Search coil in kronspice (a ring round the bar, radius a + 3 mm), flux from ∫V dt against ΣM·I and 2πψ:")
    r.table([
        ["z", "circuit ∫V dt", "mutual inductances ΣM·I", "2πψ"],
        ["0.0 mm", "1.53867×10<super>−4</super> Wb", "1.53867×10<super>−4</super> Wb",
         "1.53867×10<super>−4</super> Wb"],
        ["19.1 mm", "1.49659×10<super>−4</super> Wb", "1.49659×10<super>−4</super> Wb",
         "1.49659×10<super>−4</super> Wb"],
        ["38.1 mm", "7.85866×10<super>−5</super> Wb", "7.85866×10<super>−5</super> Wb",
         "7.85866×10<super>−5</super> Wb"],
        ["57.2 mm", "7.28009×10<super>−6</super> Wb", "7.28009×10<super>−6</super> Wb",
         "7.28009×10<super>−6</super> Wb"],
    ], [0.16, 0.28, 0.28, 0.28])
    r.p("<b>Bench protocol for the bar</b> (from report.html and doc/russell-octave-report.md):")
    r.bullets([
        "Take an ordinary bar magnet and declare its measured B<sub>r</sub>.",
        "Run a Hall probe (or the search coil modelled here) parallel to the axis at a fixed gap, from beyond one "
        "end to beyond the other, in about 1 mm steps.",
        "Count the sign changes of B<sub>z</sub>. Two supports the standard model; more than two on an unmodified "
        "bar supports Russell (five for four pairs end to end, eight if separated).",
        "Compare the trace point by point with the example's scan for the declared B<sub>r</sub>.",
    ])

    # ---------------- Report 2: the ring ----------------
    r.h1("Report 2 — The ring of twelve sections (Part 5)")
    r.p("The octave bar bent into a ring: a ring magnet cut into twelve 30° sections, each magnetized opposite to "
        "its neighbours (user request, 2026-11-05). Two ways of making neighbours opposite are built, with a control:")
    r.bullets([
        "<b>Radial.</b> Each section is magnetized along its own centre line, outward and inward in turn, so the "
        "outer face reads N S N S … — the standard multipole ring of motors, magnetic couplings and encoders.",
        "<b>Along the ring.</b> Each section is magnetized along the ring, clockwise and anticlockwise in turn, so "
        "like poles meet at every joint (N|N, S|S): the four-segment bar closed into a ring, with consequent poles "
        "at the joints.",
        "<b>Control.</b> All twelve sections magnetized outward.",
    ])
    r.p("The ring is 50 mm across, 30 mm bore, 10 mm tall, B<sub>r</sub> 1.25 T (declared); each section is "
        "uniformly (\"parallel\") magnetized, as segmented magnets are made. A uniformly magnetized body carries all "
        "its magnetic charge on its surface, σ = M·n (n the outward unit normal; Jackson §5.9): on the outer arc, the inner arc and the two cut "
        "faces; the flat top and bottom carry none. Outside, B = μ<sub>0</sub>H is the Coulomb field of that charge. "
        "Twin: the same volume as point dipoles m = M·dV (Jackson eq. 5.56). Field lines are traced by RK4 along "
        "the unit vector of B from the outer N faces, at equal steps of flux, until they reach the magnet again.")
    r.table([
        ["Check", "Radial", "Along the ring"],
        ["Surface charge against the point-dipole twin, 6 points outside", "7.8×10<super>−5</super> of |B|",
         "7.9×10<super>−5</super>"],
        ["Net magnetic charge, of Σ|Q|", "4.8×10<super>−17</super>", "2.2×10<super>−17</super>"],
        ["div B and curl B outside (central differences, of the field gradient)", "1×10<super>−7</super>, 3×10<super>−7</super>",
         "3×10<super>−7</super>, 1×10<super>−7</super>"],
    ], [0.5, 0.25, 0.25])
    fig("figR1_ring3d",
        "Figure 3 (R1). The ring in 3-D, orange N and blue S, with field lines traced from the outer N faces (four "
        "per pole at three heights) until they return into an S face. Left: radial, poles at the section centres. "
        "Right: along the ring, like poles meeting at every joint. Each line arches over to the neighbouring S "
        "pole: a crown of twelve arches, flat at mid-height and bulging above and below it.")
    r.figure(os.path.join(HERE, "ring3d_radial.gif"),
             "Figure 4. First frame of the turntable animation viz/ring3d_radial.gif: the radial ring, poles at the "
             "middle of each section.", width_frac=0.5)
    r.figure(os.path.join(HERE, "ring3d_tangential.gif"),
             "Figure 5. First frame of viz/ring3d_tangential.gif: the ring magnetized along the ring; each section "
             "is N at one end and S at the other, so the poles sit at the joints.", width_frac=0.5)
    fig("figR2_midplane",
        "Figure 6 (R2). The field in the mid-plane, seen from above; by symmetry B lies in that plane, so these "
        "lines are exact. Arrows in the sections show M. Left: radial, 12 poles outside. Middle: along the ring, "
        "12 poles at the joints. Right: control, all outward; its outer face is all N and its bore all S, so its "
        "flux leaves radially and returns over the top and bottom into the bore, and no lines close in the plane.")
    r.h2("The field round the outside")
    r.p("B<sub>r</sub> 3 mm outside the outer face, in the mid-plane: the scan a Hall probe makes when carried round "
        "the ring.")
    r.table([
        ["Ring", "Sign changes", "Peak |B<sub>r</sub>|", "N poles at"],
        ["radial, alternating (N S N S … outside)", "12", "291 mT",
         "15°, 75°, 135°, 195°, 255°, 315° (section centres)"],
        ["along the ring, alternating (N|N, S|S joints)", "12", "315 mT",
         "30°, 90°, 150°, 210°, 270°, 330° (joints)"],
        ["radial, all outward (control)", "0", "223 mT", "none (B<sub>r</sub> &gt; 0 all round)"],
    ], [0.34, 0.12, 0.13, 0.41])
    fig("figR3_scan",
        "Figure 7 (R3). B<sub>r</sub> round the ring, 3 mm outside, in the mid-plane. Twelve sign changes for both "
        "alternating rings, their peaks shifted by half a section (15°); none for the all-outward control.")
    r.p("Both alternating rings have twelve poles round the outside, six N and six S; the field falls off quickly "
        "with distance because neighbouring poles cancel. They look alike from outside; their poles sit 15° apart.")
    r.h2("A pickup coil beside the turning ring, in kronspice")
    r.p("A 200-turn, 6 mm pickup coil faces the outer side, 1.5 mm out (R 8.26 Ω from the wire, L 246 µH from "
        "Wheeler's multilayer formula, both declared), driving an 8.2 Ω load while the radial ring turns at "
        "3000 rpm.")
    r.table([
        ["Quantity", "Result"],
        ["Electrical cycles per revolution", "6 (300 Hz), one per pole pair"],
        ["Open-circuit EMF", "6.23 V peak"],
        ["Across the load", "1.754 V rms, 375 mW"],
        ["Ledger over two revolutions", "in 3.01119×10<super>−2</super> J (shaft, through the EMF) = heat "
         "3.01119×10<super>−2</super> J + Δ stored −1.5×10<super>−7</super> J"],
        ["Ledger residual", "5.8×10<super>−6</super>"],
    ], [0.4, 0.6])
    r.p("Every joule the load receives is work done turning the ring; the coil's EMF source stands in for the shaft.")

    r.h1("Validation")
    r.table([
        ["Test", "Checks", "Tolerance"],
        ["TestSheetAxis", "sheet field = exact solenoid formula on the axis, both bars",
         "10<super>−5</super> of μ<sub>0</sub>M<sub>0</sub>"],
        ["TestLoopFieldTwin", "elliptic loop field = Biot–Savart, four points", "10<super>−5</super>"],
        ["TestPoles", "2 charge sites (uniform), 5 (four segments), net charge zero", "exact"],
        ["TestScan", "sign changes 2 / 0 (uniform), 5 / 3 (four segments)", "exact"],
        ["TestSearchCoil", "kronspice fluxmeter = ΣM·I = 2πψ", "10<super>−6</super>"],
        ["TestRingTwin", "surface charge = point-dipole twin, 6 points outside, both rings", "2×10<super>−4</super>"],
        ["TestRingChargeDivCurl", "net charge zero; div B, curl B zero outside",
         "10<super>−12</super>; 10<super>−5</super>"],
        ["TestRingScan", "12 / 12 / 0 sign changes; N poles at 15° + 60°k (radial), 30° + 60°k (along the ring)",
         "exact; 0.6°"],
        ["TestRingLinesClose", "every line from an outer N face ends on an S section", "exact"],
        ["TestRingPickup", "kronspice pickup coil, 6 cycles per revolution, ledger closes",
         "exact; 10<super>−4</super>"],
    ], [0.24, 0.56, 0.2])

    r.h1("Limits of the model")
    r.bullets([
        "Demagnetization is ignored and B<sub>r</sub> is declared; a real bar's measured B<sub>r</sub> must be "
        "declared for the bench comparison.",
        "The model computes either case exactly; it cannot say which one a real bar is. That is a measurement.",
        "The ring's sections are ideal uniform magnets; the pickup coil's R and L are declared, not measured.",
        "Lens shapes made from two centres, colour as state and the octave (ratio 2) in the painting have no "
        "counterpart in the model (doc/russell-octave-report.md).",
    ])

    r.h1("Status")
    r.p("H33 is <i>untested</i>: the model and the bench protocol are built; the claim needs the Hall-probe scan of an "
        "ordinary bar. Russell's like-poles-meet pattern closed into a ring gives twelve real poles, as many as the "
        "standard radial multipole ring, but the pattern has to be built into the magnetization: a ring magnetized "
        "one way all round has no poles round its side. The bar's status is unchanged by the ring. The same words "
        "read as hydrogen are register entry H34 (<i>hydrogen_poles</i>).")

    r.h1("Reproducing the results")
    r.code("go run ./examples/octave_magnet                 # a few seconds\n"
           "go run ./examples/octave_magnet -svg viz        # also the figures and the ring GIFs\n"
           "go test ./examples/octave_magnet\n"
           "python examples/octave_magnet/viz/make_report.py examples/octave_magnet/octave_magnet_report.pdf")
    r.small("Sources: W. Russell, \"Analysis of the Octave Bar Magnet\" (painting); <i>The Universal One</i> (1926); "
            "<i>The Secret of Light</i> (1947). J. D. Jackson, <i>Classical Electrodynamics</i>, 3rd ed., §5.5, §5.6 "
            "(eq. 5.56), §5.8–5.9; W. R. Smythe, <i>Static and Dynamic Electricity</i>, §7.10; H. A. Wheeler, "
            "<i>Proc. IRE</i> 16, 1398 (1928).")
    r.build(out)
    tmpdir.cleanup()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "octave_magnet_report.pdf"))
