# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H18 report: the Pais annular cavity (US 10,144,532 B2) as a Kron ring line, as a PDF.

    python examples/pais_cavity/viz/make_report.py examples/pais_cavity/pais_cavity_report.pdf

Numbers come from examples/pais_cavity/run_output.txt (the program's output),
README.md, doc/pais-us10144532-comparison.md (the patent's claims and figures)
and export/pais_ring_4emitters_results.md (torn against full solve). The
example has no figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/pais_cavity", "The Pais annular cavity",
               "US 10,144,532 B2 as a Kron ring line, and the patent's figures against physical limits — "
               "register entry H18",
               "register entry H18", "1 MW drive scaling, copper walls")

    r.h1("Summary")
    r.p("S. C. Pais's patent \"Craft using an inertial mass reduction device\" (U.S. Navy, filed 2016, issued 2018) "
        "claims that an annular resonant cavity between two hull walls, filled with xenon, driven by microwave emitters, "
        "with the outer wall charged, vibrating and rotated, creates \"a local polarized vacuum\" and reduces the craft's "
        "inertial mass. This example solves that cavity as a closed transmission-line ring with Kron's diakoptics "
        "(Happ's two-level method), checks it against MNA and the energy ledger, and sets the patent's figures against "
        "physical ceilings. The claim is a hypothesis (H18); the standard prediction is computed beside it.")
    r.bullets([
        "<b>Kron and MNA agree:</b> 3.3×10<super>−10</super> off resonance; on the Q ≈ 1.5×10<super>4</super> "
        "resonance, 5.8×10<super>−8</super> with backward errors 8.4×10<super>−14</super> (Kron) and "
        "2.4×10<super>−14</super> (MNA). The interconnection system has order 16 against 1,280 MNA unknowns.",
        "<b>Ledger closes:</b> input power equals wall loss to 2.5×10<super>−10</super>.",
        "<b>1 MW into four in-phase emitters</b> stores 6.167 J at Q 1.477×10<super>4</super>, with 148.9 kV and "
        "2.978 MV/m at most.",
        "<b>Sectional steering does not work on a high-Q ring:</b> energy fed into one section fills the ring evenly "
        "(0.250 per quadrant) on or off resonance; only a lossy ring (Q 7.422) keeps 0.817 of it near the emitter.",
        "<b>The patent's figures exceed physical limits:</b> σ = 1 C/m² needs 1.13×10<super>11</super> V/m (wall "
        "pressure 5.65×10<super>10</super> Pa against copper's 2.1×10<super>8</super> Pa); 10<super>33</super> W/m² "
        "needs a wall speed of 3×10<super>13</super> c; 30,000 rpm at 2 m bursts a carbon-fibre hoop (margin 3.34).",
        "<b>Standard prediction:</b> no net force on a closed cavity, mass unchanged; beaming all 1 MW out gives at "
        "most P/c = 3.34 mN. <b>Status: untested.</b>",
    ])

    r.h1("The claim and the standard prediction")
    r.p("From the patent (as summarized in doc/pais-us10144532-comparison.md):")
    r.bullets([
        "An inner and an outer resonant-cavity wall form an annular cavity. Microwave emitters drive it, which makes the "
        "outer wall vibrate \"in an accelerated mode\" and \"create a local polarized vacuum\" outside it.",
        "The cavity is filled with a noble gas (xenon); the outer wall is electrically charged; the cavity is rotated "
        "\"in an accelerated mode\".",
        "Figures: emitters at 300 MHz–300 GHz; cavity Q of 10<super>4</super>–10<super>9</super>; wall vibration at "
        "10<super>9</super>–10<super>18</super> Hz; σ = 50,000 C/m² on a 2 m disc at 30,000 rpm; energy flux "
        "S = f<sub>G</sub>·(σ²/ε<sub>0</sub>)·[R<sub>r</sub>·ω + R<sub>v</sub>·ν + v·R], about 10<super>24</super> W/m² "
        "for rotation and 10<super>33</super> W/m² at σ = 1 C/m² and 10<super>11</super> Hz; steering by driving "
        "sections of the shroud separately.",
        "The patent gives no wall gap, height, voltage or drive power.",
    ])
    r.p("<b>Standard prediction:</b> a closed cavity's fields exert no net time-averaged force on the craft (momentum "
        "conservation), the mass is unchanged, and the only way to push is to emit momentum, at most P/c. The gyroscope "
        "weight reduction the patent cites (Hayasaka &amp; Takeuchi 1989) was not reproduced by three groups (Faller et "
        "al. 1990; Quinn &amp; Picard 1990; Nitschke &amp; Wilmarth 1990).")

    r.h1("The construction")
    r.p("The annulus is a closed parallel-plate TEM line around the craft (Paul §4.2): L′ = μ<sub>0</sub>d/h, "
        "C′ = ε<sub>0</sub>h/d, R′ = 2R<sub>s</sub>/h, with d the gap, h the wall height and R<sub>s</sub> the copper "
        "surface resistance at the drive frequency. Because the patent gives no dimensions, the geometry is declared. "
        "Emitters are current injections. It is the same circuit as examples/ring_towers and examples/double_hull: the "
        "closed ring transmission line at craft scale and microwave frequency.")
    r.table([
        ["parameter", "value"],
        ["radius R", "2 m (the patent's example radius)"],
        ["gap d, wall height h", "5 cm, 0.5 m"],
        ["sections", "640 (40 per wavelength)"],
        ["ring mode", "n = 16 at 381.315 MHz (λ = 0.786 m), inside the patent's 300 MHz–300 GHz band"],
        ["characteristic impedance Z<sub>0</sub>", "37.7 Ω"],
        ["walls", "copper (case E: niobium, R<sub>s</sub> = 20 nΩ, declared)"],
        ["total drive", "1 MW"],
    ], [0.35, 0.65])
    r.p("<b>Kron analysis.</b> Happ's two-level diakoptics with a single common bus (kron.TwoLevel; Happ 1971 Ch. XII "
        "§4, Fig. 12.12): the ring is torn into 16 sectors, each hinged at the inner wall, joined by 16 link branches, "
        "so the interconnection system has order 16 against 1,280 MNA unknowns. The twin is the same netlist through "
        "kronspice's MNA AC analysis.")

    r.h1("Results")
    r.p("Five drive cases (1 MW scaling; \"quadrant\" is the energy fraction in each 90° quarter of the ring, the "
        "emitter's quadrant first):")
    r.table([
        ["case", "f [MHz]", "Q", "power in", "stored energy"],
        ["A. 4 emitters at 0/90/180/270°, in phase", "381.3153", "1.477×10<super>4</super>", "1×10<super>6</super> W",
         "6.167 J"],
        ["B. top emitter only, same current, on resonance", "381.3153", "1.477×10<super>4</super>",
         "6.25×10<super>4</super> W", "0.3854 J"],
        ["C. top emitter only, between modes n and n + 1 (steering attempt)", "393.2046", "1.5×10<super>4</super>",
         "0.1845 W", "1.12×10<super>−6</super> J"],
        ["D. top emitter only, walls 2000× copper's R<sub>s</sub>", "381.3153", "7.422", "107 W",
         "3.315×10<super>−7</super> J"],
    ], [0.46, 0.12, 0.13, 0.14, 0.15])
    r.table([
        ["case", "V<sub>max</sub>", "E<sub>max</sub>", "RF σ", "energy by quadrant", "ledger residual"],
        ["A", "148.9 kV", "2.978 MV/m", "26.4 µC/m²", "0.250 0.250 0.250 0.250", "2.5×10<super>−10</super>"],
        ["B", "37.22 kV", "0.7445 MV/m", "6.59 µC/m²", "0.250 0.250 0.250 0.250", "2.7×10<super>−10</super>"],
        ["C", "0.06346 kV", "0.001269 MV/m", "0.0112 µC/m²", "0.250 0.249 0.250 0.251", "1.6×10<super>−14</super>"],
        ["D", "0.06388 kV", "0.001278 MV/m", "0.0113 µC/m²", "0.817 0.089 0.006 0.088", "1.0×10<super>−14</super>"],
    ], [0.08, 0.15, 0.17, 0.15, 0.27, 0.18])
    r.p(""
        "<b>Case E</b>, niobium walls with R<sub>s</sub> = 20 nΩ (declared; superconducting RF cavities at 2 K): "
        "Q 3.76×10<super>9</super>, and the same 6.17 J stored needs 3.93 W.")
    r.bullets([
        "<b>Sectional steering does not work on a high-Q ring.</b> One wave decays over 1/α = 2Z<sub>0</sub>/R′ ≈ 3.7 km "
        "and the ring is 12.6 m around, so energy fed into one section fills the ring evenly on or off resonance "
        "(cases B and C). Only a ring with single-digit Q keeps the energy near its emitter (case D), as the ring-tower "
        "runs found: synchronized drives move the antinodes, not which sector holds the energy.",
        "Driving between modes (case C) does not steer the energy either; it only lowers the power the ring accepts, "
        "to 0.1845 W at the same current.",
    ])

    r.h2("The patent's figures against physical ceilings (copper case, 1 MW scaling)")
    r.table([
        ["figure", "value", "limit"],
        ["σ = 1 C/m²", "E = σ/ε<sub>0</sub> = 1.13×10<super>11</super> V/m",
         "Kilpatrick 1.91×10<super>7</super> V/m at 381 MHz; Fowler–Nordheim emission 2.45×10<super>15</super> A/m² "
         "drains the charge in 4.1×10<super>−16</super> s; wall pressure 5.65×10<super>10</super> Pa against copper's "
         "2.1×10<super>8</super> Pa (UTS); RF drive to reach it in this cavity 1.44×10<super>15</super> W"],
        ["σ = 50,000 C/m²", "E = 5.65×10<super>15</super> V/m",
         "4×10<super>−3</super> of the Schwinger field 1.32×10<super>18</super> V/m; pressure 1.41×10<super>20</super> Pa; "
         "emission 1.09×10<super>25</super> A/m² drains it in 4.6×10<super>−21</super> s; drive 3.6×10<super>24</super> W"],
        ["10<super>24</super> W/m² (rotation)", "reproduced: σ²/ε<sub>0</sub>·Rω = 1.77×10<super>24</super> W/m²",
         "the magnitude of E×H beside a moving charged surface; it circulates (div S = 0), so no energy leaves"],
        ["10<super>33</super> W/m² (vibration)", "needs R<sub>v</sub>·ν = 8.85×10<super>21</super> m/s",
         "3×10<super>13</super> c"],
        ["30,000 rpm at R = 2 m", "rim 6,283 m/s, (v/c)² = 4.4×10<super>−10</super>",
         "burst margin 3.34 against a T1000G carbon-fibre hoop (spacetime.AnalyzeRotor; &gt; 1: flies apart)"],
        ["wall vibration 10<super>9</super>–10<super>18</super> Hz", "—",
         "lattice vibrations end at the Debye frequency: Cu 7.1×10<super>12</super> Hz, Al 8.9×10<super>12</super> Hz"],
        ["thrust", "—", "zero net time-averaged force on a closed cavity; beaming all 1 MW out gives P/c = 3.34 mN"],
    ], [0.2, 0.28, 0.52])

    r.h1("Twin tests and ledger")
    r.table([
        ["check", "result"],
        ["Kron vs MNA at 393.2046 MHz (off resonance)", "3.3×10<super>−10</super> relative; backward error Kron "
         "3.9×10<super>−16</super>, MNA 1.3×10<super>−13</super>"],
        ["Kron vs MNA at 381.3153 MHz (on resonance)", "5.8×10<super>−8</super> relative; backward error Kron "
         "8.4×10<super>−14</super>, MNA 2.4×10<super>−14</super>"],
        ["ledger (input − wall loss)/input", "1.0×10<super>−14</super> to 2.7×10<super>−10</super> over cases A–D"],
        ["Q against ωL′/R′", "to 2×10<super>−3</super>"],
        ["resonance against the lumped-ladder dispersion", "follows sin x/x exactly"],
        ["export sheet, case A (torn into 16 pieces of 119–120 unknowns)",
         "largest node-voltage difference 1.1×10<super>−8</super> of the largest node voltage; power taken by the "
         "R, L, C branches equal torn and untorn (Roth 1959 §5)"],
    ], [0.4, 0.6])
    r.p("On the resonance the system is ill-conditioned by about Q·n, so both solutions are checked by their backward "
        "error rather than by their difference alone.")

    r.h1("Limits of the model")
    r.bullets([
        "The geometry (gap, height, drive power) is declared; the patent gives none of it.",
        "The cavity is a 1-D TEM ring line with copper surface loss; the xenon fill, the static charge, the rotation and "
        "the wall vibration are assessed against limits, not simulated.",
        "Emitters are ideal current injections.",
    ])

    r.h1("Status and bench test")
    r.p("<b>H18: untested.</b> Bench test: weigh a sealed, microwave-driven, charged cavity on a balance against a dummy "
        "load dissipating the same power (to separate heating and convection). The standard prediction is no weight "
        "change beyond P/c of any escaping radiation, with every watt booked in through the emitters and out as heat.")

    r.h1("Reproducing the results")
    r.code("go run ./examples/pais_cavity\n"
           "go test ./examples/pais_cavity\n"
           "KRONSPICE_EXPORT=1 go test ./examples/pais_cavity -run TestExport   # export/ sheet for case A")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "pais_cavity_report.pdf"))
