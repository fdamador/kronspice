# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H46 report: an efficient cavitation heater, two dimple drums on one shaft
driven by a motor or a Tesla turbine, and the Mandelbulb challenger, as a PDF.

    go run ./examples/harmonic_cavitation -viz examples/harmonic_cavitation/viz
    python examples/harmonic_cavitation/viz/plot_all.py
    python examples/harmonic_cavitation/viz/make_report.py

Numbers come from run_output.txt and the CSVs in viz/; layout from
internal/pdfreport/kronreport.py.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1
kronreport.H2.keepWithNext = 1


def fig(name):
    return os.path.join(HERE, name)


def rows(name):
    with open(fig(name), newline="") as f:
        r = csv.reader(f)
        next(r)
        return [[float(x) for x in row] for row in r]


def build(out):
    best = rows("best.csv")
    d, b = best[0], best[1]
    dr = rows("drives.csv")
    pn = rows("pern.csv")
    front = rows("front.csv")
    se = rows("sensitivity.csv")
    etas = [x[0] for x in front]
    lp = rows("loop.csv")

    r = Report("examples/harmonic_cavitation", "Harmonic cavitation",
               "An efficient cavitation heater: a dimpled metal drum at each outlet, spinning on one shaft in a fixed "
               "housing, sized from published rotational cavitation reactors, with the Mandelbulb cavity evolved as a "
               "challenger — register entry H46", "register entry H46", "NSGA-II seed 1, 80 × 80")

    r.h1("Summary")
    r.p("Two dimpled metal drums, one at each outlet, turn on one shaft driven by a motor or a Tesla turbine. Each "
        "drum carries rows of hemispherical dimples facing a matching ring of dimples in a fixed housing across a small "
        "gap. Water is cavitated in the gap and leaves through the end caps as hot water or steam, to heat homes. The "
        "measure is <b>cavitation per joule</b>: η = collapse energy, (p − p<sub>v</sub>)·V<sub>vapour</sub> per cavity "
        "cycle, divided by shaft work. Heat per joule is no measure, because all shaft work ends as heat whatever the "
        "cavity.")
    r.bullets([
        f"<b>Heater design (3 kW shaft):</b> η = {d[1]:.3f}, {d[2]:.0f} W of collapse energy from {d[3]:.0f} W. Drums of "
        "84 mm radius with 2 rows of Ø 6 mm hemispheres, 26 per row, against a housing dimple ring; gap 0.8 mm, "
        "1,830 rpm (16 m/s), outlet 1.18 bar, σ = 0.881. With a plain housing the same drums give 0.117.",
        f"<b>Size is free:</b> the evolved front is flat in η ({min(etas):.4f}–{max(etas):.4f}) from 1.1 kW to "
        "3.0 kW of shaft power, so the drums can be sized to the heating load.",
        "<b>One speed for every dimple:</b> on a drum every dimple runs at the same speed and cavitation number, so one "
        "shaft speed puts them all at their best.",
        "<b>The design rules hold</b> in every rerun with the declared constants changed: full hemispheres, dimples on "
        f"the housing, the smallest gap, σ ≈ 0.88. The absolute η is less certain: {min(x[1] for x in se):.2f}–"
        f"{max(x[1] for x in se):.2f}, mostly with the declared drag coefficient.",
        f"<b>Mandelbulb (H46):</b> evolved over its power n and the machine, the best is η = {b[1]:.3f} (n = 14). "
        f"Each n separately gives {min(x[1] for x in pn):.3f}–{max(x[1] for x in pn):.3f}. With its fine detail "
        "removed the same bulb gives 0.191: a smooth lobed ring does as well as dimples, and the fractal detail lowers it.",
        f"<b>Energy:</b> heat is never more than the electricity. A motor on the shaft gives {dr[0][6]:.3f} collapse "
        f"per electrical joule; a pump driving a Tesla turbine gives {dr[1][6]:.3f}, a quarter. A heat pump on the "
        f"motor's {dr[0][1] / 1e3:.1f} kW would deliver {dr[0][8] / 1e3:.1f} kW of heat.",
        f"<b>Closed steam loop (the user's layout, Part 10):</b> water enters at the drums' outer top, the drums make "
        f"steam, the steam drives the Tesla turbine on the same shaft, and it is condensed (the home's heat) and pumped "
        f"back. The turbine returns only {lp[0][2]:.0f} W of 3,000 W (Carnot allows 13 % between 111 °C and 60 °C, and "
        f"the turbine is declared at 30 % of ideal), so a motor supplies the rest and the loop cannot run itself. "
        f"Boiling in the drums drops η to {lp[0][8]:.3f}: steam that leaves has not collapsed. Pressurized drums (3 bar, "
        f"no boiling) with a flash tank keep η = {lp[1][8]:.3f}.",
    ])
    r.figure(fig("machine.png"), "Figure 1. The design in section, in the user's loop. Water enters at each drum's outer "
             "top; the drums (orange: cavitation in the dimples, facing the housing's dimple ring across the gap, blue: "
             "water, gap drawn ×4) heat it; the Tesla turbine's disks reach the drums' radius, so the same gap "
             "runs on over the turbine's rim and the water enters the turbine at its periphery; the steam drives the turbine on the same shaft and leaves "
             "along the shaft to the condenser (the home's heat), from which the water is pumped back.", max_height_mm=80)

    r.h1("What is being tested")
    r.p("<b>Hypothesis (H46, the user's).</b> Cavities shaped as a Mandelbulb, the 3-D power-n fractal of White and "
        "Nylander, give the most cavitation on a rotor turned by a Tesla turbine. Their harmonics and implosions gather "
        "in the middle. The heat, as steam, heats homes. The power n is to be evolved. The user then fixed the layout: "
        "a dimple drum at each outlet.")
    r.p("<b>Standard prediction.</b> Rotational cavitation reactors are well studied: cavities on a rotor passing a "
        "fixed surface or fixed cavities across a gap of a few millimetres. Their measured heaters return 80–90 % of "
        "the electricity as heat (Kwon &amp; Yoon 2012; Yoon group). All shaft work becomes heat, sound and steam. "
        "Focused collapses raise the peak pressure, not the energy.")

    r.h1("The model")
    r.p("A reduced model calibrated on published reactors, not CFD (CLAUDE.md §1 keeps CFD and bubble dynamics out "
        "of the engine; the bubble enters as Rayleigh–Plesset). Every constant is either fitted to a published number "
        "or declared, and the declared ones are run at the ends of their ranges.")
    r.table([
        ["quantity", "value", "origin"],
        ["inception σ<sub>i</sub>", "1.105 (bracket 0.968–1.243)", "Xue et al. 2024: none at 3000 rpm, bubbles at 3400 (Ø 80 mm); 1 atm declared"],
        ["best σ", "0.775", "Xue 2024: strongest pulsation at 3800 rpm"],
        ["pulsation St = fW/U", "0.446; also 0.28", "Xue 2024 (1122 Hz at 3600 rpm); cloud shedding (Stutz &amp; Legoupil)"],
        ["vapour law", "g = x·e<super>1−x</super>, x = (σ<sub>i</sub> − σ)/Δ, Δ = 0.272", "Δ fitted so Xue's collapse rate peaks at 3800 rpm"],
        ["vapour at the peak κ<sub>v</sub>", "0.2 of the cavity", "declared, order of Sun et al. 2021's 5.44×10<super>−8</super> m³ beside Ø 10 mm units"],
        ["cavity drag C<sub>d</sub>", "0.25 on the projected area", "declared, order of Sun 2021's 1.63 N·m (one 1/32 sector assumed)"],
        ["housing pulse", "0.3·½ρU²·(1.5 mm/gap)<super>0.60</super>", "exponent fitted to Petkovšek 2013 (0.7 bar at 3.5 mm, 1.7 at 0.8 mm); 0.3 declared"],
        ["wall friction", "max(μU/gap, turbulent flat plate)", "Schlichting c<sub>f</sub> = 0.0592 Re<super>−0.2</super>; drum surface and both end faces (hub 10 mm)"],
        ["gap pressure", "the outlet pressure", "set at the drum's end (a gene)"],
    ], [0.2, 0.3, 0.5])
    r.p("Per cavity, the collapse energy per event is (p<sub>local</sub> − p<sub>v</sub> + housing pulse)·V<sub>vapour</sub>. "
        "The event rate is the cloud rate St·U/W on a plain housing, and the other ring's passing rate, capped at the "
        "cloud rate, with housing dimples. It is never faster than the Rayleigh collapse allows. Each design is "
        "scored at the worse St. Collapse energy is a part of the cavity drag, and a design claiming more is rejected.")
    r.p("<b>The shaft as a kronspice circuit.</b> In the mobility analogy (angular speed across, torque through) the "
        "shaft is a speed source ω driving one conductance per loss: friction, cavity drag less collapse, and "
        "collapse. Each absorbs ω·torque. The GA ranks every design by the power the operating point books into the "
        "collapse branch over the power the source gives. The ledger closes to round-off (residual 0).")

    r.h1("Part 1–2. Calibration and drum sizing")
    r.p("The resonance rule sets each ring's dimple count so the passing frequency equals the cloud frequency: "
        "N = 2πr·St/W, a pitch of about 2.2 W. On Xue's own rotor it gives 18.7 grooves against the 20 built.")
    r.p("A dimple works only within a window of σ = (p − p<sub>v</sub>)/(½ρU²). At 1 atm cavitation starts near "
        "U = 13.3 m/s and is strongest near 15.9 m/s, so a drum of 50 mm radius cavitates from 2,550 rpm and is "
        "strongest at 3,040 rpm; one of 84 mm, from 1,520 and at 1,810 rpm. On a drum every dimple runs at the same "
        "speed, so one shaft speed suits them all.")
    r.figure(fig("calibration.png"), "Figure 2. Left: the model's collapse rate on Xue's rotor against speed, with the "
             "measured no-bubble, first-bubble and best speeds. Right: the fitted vapour law g(σ), zero above inception "
             "and falling again as the cavity chokes.", max_height_mm=75)

    r.h1("Part 3–4. The search")
    r.p(f"NSGA-II (population 80, 80 generations) over drum radius, rows, dimple diameter and depth, rotor and housing "
        f"pitch, housing dimples on or off, gap, speed and outlet pressure. Two objectives: η and collapse power. "
        f"Constraints: shaft ≤ 3 kW, land ≥ 1 mm, energy consistency. {len(front)} feasible designs on the front, all "
        f"within η {min(etas):.4f}–{max(etas):.4f}; the heater design is the one with the most collapse power within "
        f"0.5 % of the best η.")
    r.table([
        ["", "heater design (two drums)"],
        ["drum radius / rows", "84.3 mm / 2"],
        ["dimples", "Ø 6.0 mm, depth W/2 (hemisphere), 26 per row, pitch 3.3 W"],
        ["housing", "dimple ring, pitch 1.2 W (73 per ring)"],
        ["gap / speed / outlet", "0.80 mm / 1,828 rpm (16.1 m/s) / 1.18 bar"],
        ["σ", f"{d[6]:.3f}"],
        ["shaft = friction + drag", f"{d[3]:.0f} W = {d[4]:.0f} + {d[5]:.0f} W"],
        ["collapse energy / η", f"{d[2]:.0f} W / {d[1]:.4f}"],
    ], [0.3, 0.7])
    r.figure(fig("front.png"), "Figure 3. Left: the Pareto front (η against collapse power, colour σ): flat in η, so "
             "the size follows the load. Right: the best η when the search is rerun with one declared constant changed. "
             "The drag coefficient dominates.", max_height_mm=75)
    r.p("The most collapse per joule comes at σ ≈ 0.88, a slower speed than the most collapse per second (σ = 0.775). "
        "Every search drives the gap to the smallest value in the data (0.8 mm), so the true optimum may be smaller "
        "still, outside the calibrated range.")

    r.h1("Part 5. The gap as an acoustic ring")
    r.p("The annular gap is a kronspice ladder closed on itself: series inertance ρΔx/A, shunt compliance AΔx/ρc², "
        "64 sections, driven by a 1 m³/s volume-velocity source, Q 2000. Its peaks equal the lattice's exact "
        "frequencies to 10<super>−4</super>, and these tend to m·c/2πR within the lattice's sin x/x. The 26-dimple "
        "pattern passes at 792 Hz, while its own mode 26 is at 73.7 kHz in clear water. The two meet only if the "
        "bubbly mixture's sound speed falls to the rotor speed, 16.1 m/s. Wood's minimum at 118 kPa is 21.7 m/s, so "
        "this resonance is out of reach at this design.")
    r.figure(fig("ring.png"), "Figure 4. |p| at the driven node of the gap ring (kronspice AC), clear water: the first "
             "four azimuthal modes.", max_height_mm=60)

    r.h1("Part 6. One collapse, and collapses on the axis")
    r.p("The cavity's vapour (11.1 mm³, a sphere of R<sub>0</sub> = 1.38 mm) collapses under 171 kPa: the local "
        "114 kPa plus the housing pulse of 57 kPa. Empty, it collapses in 96.66 µs against Rayleigh's 96.69 µs. With "
        "10<super>−3</super> gas, surface tension and viscosity, its energy books close to 2.1×10<super>−6</super>.")
    r.figure(fig("collapse.png"), "Figure 5. The heater design's cavity: radius, and the pressure it radiates to the "
             "axis 84 mm away (incompressible Rayleigh–Plesset; the final peaks are overstated).", max_height_mm=70)
    r.p("<b>Gathering in the middle.</b> Rotor and housing rings line up gcd(N<sub>r</sub>, N<sub>h</sub>) cavities at "
        "a time. With 26 and 73 that is one; with equal counts all 26 collapse together and their pulses reach the "
        "axis together, raising the peak there by up to 26×. The energy is still only the 26 collapses' share of the "
        "shaft work. These figures are upper bounds (linear superposition of incompressible near fields). In a solid "
        "drum the axis is metal, so the focused pulse goes into the rotor; a hollow drum with a water core would put "
        "it in the water.")

    r.h1("Part 7. Energy books")
    r.table([
        ["drive", "electric", "shaft", "heat to water", "heat to room", "collapse / electric J", "steam"],
        ["motor on the shaft", f"{dr[0][1]:.0f} W", f"{dr[0][2]:.0f} W", f"{dr[0][3]:.0f} W", f"{dr[0][4]:.0f} W", f"{dr[0][6]:.3f}", f"{dr[0][7]:.2f} kg/h"],
        ["pump + Tesla turbine", f"{dr[1][1]:.0f} W", f"{dr[1][2]:.0f} W", f"{dr[1][3]:.0f} W", f"{dr[1][4]:.0f} W", f"{dr[1][6]:.3f}", f"{dr[1][7]:.2f} kg/h"],
    ], [0.22, 0.11, 0.1, 0.15, 0.14, 0.16, 0.12])
    r.figure(fig("drives.png"), "Figure 6. Where the electricity goes. Declared: motor 0.90, pump 0.70, Tesla turbine "
             "0.30 (to be measured).", max_height_mm=60)
    r.p("The turbine's and pump's losses also warm the water, so heat per joule is similar, but cavitation per joule is "
        "four times lower. A heat pump (COP 3.5) on the same electricity gives 3.5× the heat of either.")

    r.h1("Part 8. The Mandelbulb challenger (H46)")
    r.p("The bulb: v → v<super>n</super> + c in spherical form, polar convention (θ = arccos z/r), 10 iterations, "
        "detail below 0.3 mm (metal printing) removed. Tested: both conventions have (n − 1)-fold symmetry about the "
        "axis; in the latitude convention (θ = arcsin z/r) the equator is exactly the degree-n Multibrot set; n = 8 "
        "shows 7 lobes. The bulb's latitude rings, scaled to the drum radius, form the drum's surface, and each valley "
        "between lobes is a cavity with its own width, depth and volume.")
    r.figure(fig("bulb.png"), "Figure 7. Left: the equator of the n = 3, 8 and 12 bulbs (outer circle: the drum "
             "radius). Right: the best η for each n, each from its own NSGA-II search (50 × 50), against the dimple drums.",
             max_height_mm=70)
    r.bullets([
        f"NSGA-II over n and the machine (60 × 60): best η = {b[1]:.3f} at n = 14, {b[8]:.0f} cavities, on 31 mm drums at "
        "6,740 rpm.",
        "The n = 8 bulb put on the heater's 84 mm drums without re-optimizing gives 0.099 at 4.7 kW: the bulb needs its "
        "own size and speed.",
        "Fine detail removed (3 iterations against 10): η = 0.191, as good as the dimples (0.189). The small structure "
        "adds drag faster than it adds collapse.",
        "Status: <b>bounded</b> in this model. The Mandelbulb works as a ring of lobes; a smooth lobed ring matches the "
        "dimples, and the fractal detail does not help.",
    ])

    r.h1("Part 10. The closed steam loop")
    r.p("The user's layout: water in at each drum's outer top; the drums' cavitation and shaft work raise steam; the "
        "steam drives the Tesla turbine on the drums' own shaft; the exhaust is condensed at the outlet, where its heat "
        "goes to the home, and the condensate is pumped back. A motor on the same shaft supplies what the turbine does "
        "not. Saturation pressures from IAPWS-IF97 (verified to 10<super>−8</super> at 300 K); steam states from the "
        "IAPWS-95 tables; declared: turbine isentropic efficiency 0.30, condenser 0.2 bar (60 °C, warm enough for "
        "radiators), motor 0.90, pump 0.70.")
    r.p("<b>Cavitation needs subcooled water.</b> Collapse energy is (p − p<sub>v</sub>)·V. As the water warms, its "
        "vapour pressure p<sub>v</sub> rises toward the drum pressure and the collapse fades; at boiling it is gone, "
        "and steam that boils off and leaves the drum has not collapsed. So two ways to raise the steam are compared. "
        "<b>A:</b> the drums are the boiler at 1.5 bar. <b>B:</b> the drums run at 3 bar, where water at 100–120 °C is "
        "still well below boiling, and the hot water is throttled into a flash tank at 1 atm, whose steam drives the "
        "turbine while the remaining hot water returns to the drums.")
    lp = rows("loop.csv")
    r.table([
        ["cycle (3 kW shaft)", "steam", "turbine", "motor shaft", "electric", "home heat", "collapse", "η", "collapse / elec J"],
        ["A: boil in the drums, 1.5 bar", f"{lp[0][1]:.2f} kg/h", f"{lp[0][2]:.0f} W", f"{lp[0][3]:.0f} W", f"{lp[0][4]:.0f} W",
         f"{lp[0][5]:.0f} W", f"{lp[0][7]:.0f} W", f"{lp[0][8]:.3f}", f"{lp[0][9]:.3f}"],
        ["B: subcooled drums 3 bar + flash", f"{lp[1][1]:.2f} kg/h", f"{lp[1][2]:.0f} W", f"{lp[1][3]:.0f} W", f"{lp[1][4]:.0f} W",
         f"{lp[1][5]:.0f} W", f"{lp[1][7]:.0f} W", f"{lp[1][8]:.3f}", f"{lp[1][9]:.3f}"],
    ], [0.24, 0.09, 0.08, 0.1, 0.09, 0.1, 0.09, 0.07, 0.14])
    r.figure(fig("loop.png"), "Figure 8. Left: η against water temperature at three drum pressures (heater drums at the "
             "speed best for 60 °C water); it falls to zero as the water nears boiling, later at higher pressure. Right: "
             "the two loops at 3 kW of shaft work: the collapse share, the rest of the drum heat, and the little the "
             "turbine returns.", max_height_mm=75)
    r.bullets([
        "<b>The loop cannot run itself.</b> Heat becomes work at most at the Carnot efficiency, 13 % between 111 °C and "
        "60 °C, and a Tesla turbine at 30 % of ideal returns 3–4 % of the shaft work. The motor supplies about 96 %.",
        "<b>The turbine does not add heat.</b> The home receives the condenser's heat, about the electricity less the "
        "motor's losses, as with no turbine at all. The books close (electricity = home + room to 10<super>−12</super> W), "
        "and kronspice's shaft circuit, with the motor as a speed source and the turbine as a torque source, has the "
        "motor supply P<sub>shaft</sub> − P<sub>turbine</sub> with a closed ledger.",
        "<b>Boiling wastes the cavitation:</b> 91 % of the heat in A goes into boiling, where nothing collapses, so η "
        f"falls to {lp[0][8]:.3f}. Subcooled drums with a flash tank keep η = {lp[1][8]:.3f} (the hot water at 3 bar has "
        "about the same p − p<sub>v</sub> as cold water at 1 atm).",
        "Collapse in hot water is further reduced by heat transfer into the bubble (the thermal effect, Brennen 1995 "
        "§2.3), not modelled here, so the hot-water values are upper bounds.",
    ])

    r.h1("Limits of this study")
    r.bullets([
        "Reduced model, not CFD. Calibrated on one reactor's grooves (Xue 2024) and applied to hemispheres; the drag, "
        "vapour and pulse constants are declared from the order of published values.",
        "Gap at the edge of the calibrated range in every search (0.8 mm).",
        "Rayleigh–Plesset is incompressible; peak pressures and focusing are upper bounds.",
        "Turbine, pump and motor efficiencies are declared.",
        "Next: CFD of the heater drums (OpenFOAM, Schnerr–Sauer and SST k-ω, as Sun et al. 2020), then a bench with KI "
        "dosimetry, hydrophones, a power analyser and a steam meter.",
    ])

    r.h1("Run it")
    r.code("go run ./examples/harmonic_cavitation\n"
           "go run ./examples/harmonic_cavitation -viz examples/harmonic_cavitation/viz\n"
           "python examples/harmonic_cavitation/viz/plot_all.py\n"
           "python examples/harmonic_cavitation/viz/make_report.py\n"
           "go test ./examples/harmonic_cavitation")

    r.h1("Sources")
    r.bullets([
        "L. Xue, Z. Hao, W. Ren et al., <i>Ultrason. Sonochem.</i> 109, 106999 (2024).",
        "X. Sun, X. Xuan, Y. Song et al., <i>Ultrason. Sonochem.</i> 70, 105311 (2020); X. Sun et al., "
        "<i>Ultrason. Sonochem.</i> (2021), doi:10.1016/j.ultsonch.2021.105771.",
        "M. Petkovšek, M. Zupanc, M. Dular et al., <i>Sep. Purif. Technol.</i> 118, 415 (2013).",
        "S. Kwon, J. Y. Yoon, <i>Proc. IMechE E</i> 227, 67 (2012), doi:10.1177/0954408912451535.",
        "Y.-H. Zhang et al., <i>Ultrason. Sonochem.</i> 124, 107727 (2025), review of rotational cavitators.",
        "Lord Rayleigh, <i>Phil. Mag.</i> 34, 94 (1917); C. E. Brennen, <i>Cavitation and Bubble Dynamics</i>, OUP 1995.",
        "A. B. Wood, <i>A Textbook of Sound</i>, 1930; H. Schlichting, <i>Boundary-Layer Theory</i>.",
        "D. White, P. Nylander, the Mandelbulb (2007–2009); D. Hvidtfeldt, distance-estimated 3-D fractals (2011).",
        "IAPWS-95 and IAPWS-IF97 (water); IAPWS R12-08 (viscosity).",
    ])
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "harmonic_cavitation_report.pdf"))
