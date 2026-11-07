# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Tesla's homopolar dynamo (CLAUDE.md §3f, milestone M9n) as a PDF report.

    python examples/homopolar/viz/make_report.py [out.pdf]

Numbers come from examples/homopolar/README.md and examples/homopolar/REPORT.md
(both written from `go run ./examples/homopolar` and `-dual`); figures are
viz/dual_motor_layouts.png, viz/dual_motor_resonance.png and
viz/dual_motor_schematic.png. REPORT.md is left unchanged; this PDF covers it.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def main(out):
    r = Report("examples/homopolar", "Tesla's homopolar dynamo",
               "The disk EMF from the Lorentz force, Tesla's two-disk machine, and a cross-excited "
               "dual-motor machine with resonance, an output tank and rim bearings — CLAUDE.md §3f, M9n",
               "CLAUDE.md §3f, M9n", "GA: 20 generations of 30, seed 1")

    r.h1("Summary")
    r.p("Tesla's U.S. Patent 406,968 (1889) puts two conducting disks on one shaft in fields of opposite "
        "polarity, joins their rims with a flexible metal belt and takes the current from the two shafts, so "
        "no brush bears on a fast rim. This example derives the disk EMF from the Lorentz force by grade "
        "selection, builds the machine from kronspice's N machine element, compares Tesla's two-disk "
        "arrangement with a single disk and a rim brush, and evolves the two-disk machine for efficiency. "
        "A second part (the dual motor, written up in REPORT.md) removes the magnets: two disks on their "
        "own motors excite each other through each other's field coils, and the work studies resonance, an "
        "RLC output tank with a transformer and rectifier, and ball bearings as rim collectors.")
    r.bullets([
        "<b>Lorentz EMF = formula.</b> The EMF per disk from q(v×B) is 2.8149 V, equal to "
        "½Bω(r<sub>o</sub>² − r<sub>i</sub>²) (tested to 10<super>−12</super>).",
        "<b>Two disks beat one disk with a rim brush.</b> Into the same 10 mΩ load the two-disk machine gives "
        "3044.6 W at efficiency 0.9801; the single disk with a rim brush gives 711.7 W at 0.9477. Ledger "
        "residual 0 in both.",
        "<b>Evolved machine.</b> r<sub>o</sub> 0.361 m, 6.8 mm thick, 1.22 T, 1 Ω load: 24.93 V per disk, "
        "2485 W into the load from 2486 W at the shaft, efficiency 0.9998.",
        "<b>Stacking lowers the self-excitation threshold.</b> Stacked 8 cm apart the dual motor self-excites "
        "at 794 rad/s (rim 119 m/s), below the 157 m/s copper hoop limit; side by side it needs 1206 rad/s "
        "(181 m/s).",
        "<b>The machine resonates without a capacitor</b> (shaft inertia as the second store): "
        "f<sub>em</sub> = 0.587 Hz, Q = 41. At the 7.8 Hz Schumann frequency the response is forced, not "
        "resonant: ±346 A on 6.5 kA.",
        "<b>Output tank.</b> One R<sub>L</sub> || C || L<sub>p</sub> tank between the coils, with a 1 : 12.5 "
        "transformer, gives 5.34 kW to the AC load at h = 1 (scaled to 10 kA); rectified, the bank holds "
        "449 V and feeds 5.1 kW (scaled to 10 kA).",
        "<b>Rim ball bearings</b> work only with superconducting field coils: a PI-regulated design holds "
        "19.81 A at 1.60 T, 102 W to the load and 29 W in the bearings from 131.5 W at the drives. Every watt "
        "is paid by the motors; the ledgers close.",
    ])

    r.h1("Background")
    r.p("A disk spinning at angular speed ω in an axial field B carries its charges at v = ω×r. The force "
        "q(v×B), the grade-1 part of &lt;F·Q<super>†</super>&gt; (<i>em.LorentzGA</i>; Joot, Theorem 3.16), pushes "
        "them radially, so the disk is a source of EMF between axle and rim. <i>em.MotionalEMF</i> integrates "
        "v×B along a radius and must equal")
    r.code("EMF = ½ B ω (r_o² − r_i²)")
    r.p("This is standard physics (CLAUDE.md §3f, the electric pair of the Lorentz force). The mechanical "
        "torque needed to keep the disk turning equals the electrical power drawn, so the energy ledger shows "
        "inflow from the shaft and outflow to the load and losses.")

    r.h1("The model")
    r.p("Each disk is one N machine element with two windings, the disk's radial path and the field winding "
        "(driven at 1 A, the field strength in the coupling), and one speed-matrix entry "
        "K = ½·B·(r<sub>o</sub>² − r<sub>i</sub>²). The element gives the speed voltage ω·K·i<sub>f</sub> in "
        "the disk and the torque K·i<sub>disk</sub>·i<sub>f</sub> on the shaft (Kron 1951, Epilogue eq. 9). "
        "Each shaft is a mechanical node driven at ω, so the ledger books the shaft power as inflow and shows "
        "it leaving as heat in the load, the disks, the belt and the contacts.")
    r.table([
        ["Reference design", "Value"],
        ["Copper disks", "15 cm, 5 mm thick, 1 cm axle"],
        ["Field, speed", "0.8 T, 3000 r/min"],
        ["Load", "10 mΩ"],
        ["Contacts and belt (declared, illustrative, not measured)",
         "50 µΩ shaft, 500 µΩ rim brush, 100 µΩ belt"],
    ], [0.55, 0.45])

    r.h1("Results: Tesla's two-disk machine")
    r.table([
        ["", "Two disks, belt (Tesla)", "Single disk, rim brush"],
        ["EMF per disk, Lorentz force", "2.8149 V", "2.8149 V"],
        ["½Bω(r<sub>o</sub>² − r<sub>i</sub>²)", "2.8149 V", "2.8149 V"],
        ["Disk spreading resistance", "1.45 µΩ", "1.45 µΩ"],
        ["Load", "5.5178 V, 3044.6 W", "2.6678 V, 711.7 W"],
        ["Shaft (ledger inflow)", "3106.4 W, torques 4.94 + 4.94 N·m", "750.9 W, torque 2.39 N·m"],
        ["Efficiency", "0.9801", "0.9477"],
        ["Ledger residual", "0", "0"],
    ], [0.36, 0.34, 0.30])
    r.bullets([
        "The two disks' EMFs add (open circuit: 2 × EMF at the terminals, to 10<super>−6</super>).",
        "The belt costs far less than a brush on a rim moving at 47 m/s, so the two-disk machine puts more "
        "than four times the power into the same load.",
        "<b>Evolved</b> (genes: outer radius, thickness, field up to 1.5 T, load) for at least 2 kW with the "
        "least loss, 20 generations of 30, seed 1 (562 evaluations): r<sub>o</sub> 0.361 m, 6.8 mm thick, "
        "1.22 T, 1 Ω load, 24.93 V per disk, 2485 W into the load from 2486 W at the shaft, efficiency 0.9998.",
    ])

    r.h1("Results: the cross-excited dual motor")
    r.p("The same copper disks without magnets: each disk on its own motor, its current running through the "
        "field coil of the other unit (Rikitake's cross-excitation, as in <i>earth_machine</i>), optionally with "
        "a feedback coil carrying its own current. The coupling M (ω-EMF per ampere) is the coil's flux through "
        "the disk annulus over 2π (<i>em.MutualCircles</i> averaged over the turns), L is Nagaoka's, R the "
        "copper of a 4-turn coil in a 6 cm × 2 cm window. Liquid-metal shaft contacts (5 µΩ), the motor rotors "
        "and the coupling of two coils on one former (k = 0.95) are declared.")
    r.table([
        ["", "Side by side", "Stacked on one axis, 8 cm apart"],
        ["Coupling M, own former", "2.196×10<super>−7</super> V·s/A", "2.196×10<super>−7</super> V·s/A"],
        ["Coil's stray coupling to its own disk", "— (a feedback coil buys it with copper)", "0.519 M, free"],
        ["Coil–coil mutual", "—", "3.77 µH (k = 0.39)"],
        ["Self-excitation threshold", "1206 rad/s, rim 181 m/s (K = 0); 1180 rad/s with K = 1",
         "794 rad/s, rim 119 m/s"],
        ["Against the copper hoop limit (157 m/s, 220 MPa declared)", "above", "below"],
    ], [0.36, 0.34, 0.30])
    r.figure(os.path.join(HERE, "dual_motor_layouts.png"),
             "Figure 1. The two dual-motor layouts with the coils' field lines (blue, darker = stronger), copper "
             "disks and coil turns. (a) Side by side: each coil reaches only its partner's disk; threshold "
             "1206 rad/s, rim 181 m/s, above the copper limit. (b) Stacked 8 cm apart: coil 1's stray field also "
             "threads disk 1 (built-in feedback K = 0.52) and the coils couple with k = 0.39; threshold 794 rad/s, "
             "rim 119 m/s.")
    r.bullets([
        "A feedback coil side by side adds coupling and its own resistance in proportion, so the threshold "
        "hardly moves (K = −0.5 raises it to 3566 rad/s).",
        "With both shafts held at one speed the field grows at exactly "
        "(−R + (1+K)·M·ω)/(l<sub>11</sub> + l<sub>12</sub>) in both layouts (engine against exact, "
        "10<super>−3</super>).",
        "Driven by torque (20 N·m each), Rikitake's μ = R·√(J/(G·M·L)) is 9.5: the machine excites itself and "
        "settles to a steady 8 kA field; stacked it runs at 135 m/s on the rim, side by side at 206 m/s (past "
        "the limit).",
        "Reversals need μ ≈ 1: 1800 N·m and 4.7 MW at the shafts, 16 reversals in 11.5 s (mean 0.63 s) side by "
        "side, at a rim past the limit; with a feedback coil or stacked (both K &gt; 0) the field locks to one "
        "polarity.",
        "Start-up overshoots the running speed, so a real machine starts on its speed drives. The ledger closes "
        "to 10<super>−9</super>–10<super>−6</super> of the drive work.",
    ])

    r.h2("Resonance and the load")
    r.p("The stacked machine has a 50 µΩ load in each loop and bearing and air friction of "
        "5×10<super>−3</super> N·m·s (both declared). At 20 N·m it settles at ω* = R/M<sub>s</sub> = 944 rad/s "
        "(rim 142 m/s), 6.77 kA per loop, 2.29 kW into each load from 18.9 kW per motor. A homopolar loop alone "
        "(R and L) cannot resonate; the shafts' inertia is a capacitor on the mechanical channel, so the "
        "self-excited machine is an L–J oscillator. At the steady state the loop's net resistance is zero, so "
        "to first order only friction damps it:")
    r.code("f_em = √(2·M_s·G/(L_s·J)) / 2π = 0.587 Hz,   Q = J·2πf_em/b = 41")
    r.table([
        ["Oscillating speed by", "Result"],
        ["torque ±0.5 N·m (±2.5 %), swept", "resonance at 0.575 Hz (softened from 0.587 by the drive), δω 50–70 "
         "rad/s against 1.6 at half the frequency; circuit = RK4 to 0.2 %, ledger 2×10<super>−8</super>"],
        ["torque ±0.5 / ±2 / ±5 N·m at f<sub>em</sub>", "harmonics of the loop current: 2nd 15 / 24 / 32 %, "
         "3rd 2.6 / 4.3 / 7.5 %, 4th ≤ 2.5 %; no reversal. At ±5 N·m the speed peaks at 10 200 rpm, right at "
         "the copper limit"],
        ["speed drives ±30 %, no capacitor", "growth rate (−2.3 1/s at 0.9 ω*) independent of the oscillation "
         "frequency: no resonance, the loop is first order"],
        ["speed drives ±30 %, capacitor across the load", "slow oscillations (f<sub>p</sub> ≤ about f0) lift the "
         "field because the speed spends part of each cycle above the AC onset; faster ones average out"],
    ], [0.32, 0.68])
    r.figure(os.path.join(HERE, "dual_motor_resonance.png"),
             "Figure 2. (a) No capacitor: speed oscillation δω and current oscillation δI against the frequency "
             "of a ±0.5 N·m torque oscillation on 20 N·m; peak 0.575 Hz (linear f<sub>em</sub> 0.587 Hz), Q = 41. "
             "(b) Driven ±5 N·m at f<sub>em</sub>: loop current and shaft speed against time; the current is "
             "non-sinusoidal (2nd harmonic 32 %, 3rd 7.5 %). (c) Capacitor across the load, torque-driven at "
             "20 N·m: 300 F stays DC (2.3 kW in the load), 1000 F and 12 000 F lock in as AC (2.0 kW and 1.0 kW).")
    r.p("With the shafts held at speed the tank locks in once C &gt; L/R<sub>L</sub>² (5401 F); with free, "
        "torque-driven shafts the inertia joins the tank and the DC state gives way at C ≈ 454 F (bisection, "
        "from 1 % off the DC state).")
    r.table([
        ["C across each load", "tank f0", "Q", "outcome", "loop current", "load"],
        ["30 F", "7.9 Hz", "0.07", "DC", "6.77 kA", "2.29 kW"],
        ["300 F", "2.5 Hz", "0.24", "DC", "6.77 kA", "2.29 kW"],
        ["1000 F", "1.37 Hz", "0.43", "AC", "±7.6 kA at 0.22 Hz, 3rd harmonic 66 %", "2.03 kW"],
        ["3000 F", "0.79 Hz", "0.75", "AC", "±8.0 kA at 0.23 Hz, 3rd 59 %", "1.83 kW"],
        ["12 000 F", "0.40 Hz", "1.49", "AC", "±8.5 kA at 0.27 Hz, 3rd 37 %", "1.00 kW"],
        ["24 000 F", "0.28 Hz", "2.11", "AC", "±7.3 kA at 0.14 Hz, 3rd 78 %", "0.67 kW"],
    ], [0.15, 0.11, 0.08, 0.10, 0.40, 0.16])
    r.bullets([
        "Locked in, the machine is an AC generator below both the tank and f<sub>em</sub>; it gives the load "
        "less, because much of the AC circulates between the coils and the capacitor. Circuit = RK4 to 0.6 %, "
        "ledger 10<super>−8</super>.",
        "A faster speed swing does not let the capacitor shrink: Q = R<sub>L</sub>·√(C/L) falls as √C. With "
        "the rim-safe swing h = 0.16 only 3000 F grows (+0.17 1/s at f0/4).",
        "Ceilings for this copper machine: tank f0 ≤ 0.99 Hz with the speed held (Q ≥ 1), ≤ 3.4 Hz with free "
        "shafts (Q ≥ 0.29); with no capacitor, f<sub>em</sub> is 0.59 / 1.9 / 5.9 Hz at 20 / 200 / 2000 N·m.",
        "<b>Shafts oscillated at 7.8 Hz</b> (first Schumann frequency): ±101 rad/s about ω* takes ±274 N·m "
        "per shaft and gives a forced ±346 A at 7.8 Hz on the 6.5 kA DC, exactly "
        "M<sub>s</sub>·I·δω/(2πf·L<sub>s</sub>); harmonics 1.4 % and 0.5 %; the load gets only 3 W at 7.8 Hz. "
        "Resonance at 7.8 Hz would need a 31 F tank (Q = 0.07) or 3.5 kN·m, 3.3 MW per shaft.",
        "As a Schumann source the coils are a 282 A·m² dipole at 7.8 Hz, above the natural Schumann field "
        "(about 1 pT) only within 300 m, radiating 6×10<super>−22</super> W: a local 7.8 Hz magnet, not a "
        "Schumann source.",
    ])

    r.h2("rpm = 1768·(1 ± h·sin 2π·7.8t)")
    r.p("The speed drives impose the profile. For a self-excited loop at a held speed")
    r.code("I(t) ~ exp((M_s·ω0 − R)/L_s · t) · exp(−β·cos 2πft),   β = M_s·ω0·h / (L_s·2πf)    (~ : proportional to)")
    r.p("so the mean speed decides growth and the sine only ripples the field, with harmonics "
        "2·I<sub>n</sub>(β)/I<sub>0</sub>(β) (modified Bessel functions).")
    r.bullets([
        "<b>Air coils:</b> the machine self-excites at 9014 rpm, so at 1768 rpm the field decays at "
        "−18.7 1/s for h = 0, 0.5 and 1 alike.",
        "<b>Iron magnetic circuit</b> (gain k = 5.10, declared, chosen so 1768 rpm is the onset): "
        "L<sub>s</sub> = 69 µH, mean field 1.6 T at 10.5 kA; ripple 0.0093, 0.047 and 0.093 of the DC at "
        "h = 0.1, 0.5 and 1; RK4 and circuit equal the exact value to 3×10<super>−4</super>, ledger "
        "3×10<super>−9</super>.",
        "<b>Antiphase</b> shafts: the field decays slightly (−0.04 1/s at h = 1). A 6 F or 24 F capacitor "
        "(Q = 0.015–0.03) changes nothing.",
        "Drives: ±50 N·m at h = 0.1, ±503 N·m at h = 1; up to 186 kW is swapped with the rotor every cycle, so "
        "the drives must be four-quadrant (regenerative).",
    ])
    r.table([
        ["h", "DC", "AC at 7.8 Hz", "at 15.6 Hz", "load DC", "load AC", "AC share"],
        ["0.1", "10 kA", "93 A peak (66 A rms)", "0.2 A", "5.0 kW", "0.2 W", "0.004 %"],
        ["0.5", "10 kA", "466 A (330 A rms)", "5.4 A", "5.0 kW", "5.4 W", "0.11 %"],
        ["1", "10 kA", "932 A (659 A rms)", "22 A", "5.0 kW", "22 W", "0.43 %"],
    ], [0.07, 0.10, 0.25, 0.13, 0.13, 0.14, 0.18])

    r.h2("One RLC tank between the two machines' coils, transformer and rectifier")
    r.p("One current runs shaft 1 → disk 1 → R<sub>int</sub> → coil 1 → R<sub>L</sub> || C || L<sub>p</sub> → "
        "coil 2 → R<sub>int</sub> → disk 2 → shaft 2 → common. Shaft 2 counter-rotates and coil 2 is wound to "
        "aid coil 1; co-rotating, the EMFs cancel and the current decays as a plain R–L loop (test: −1.0 1/s, "
        "matching the formula to 5 %). This is a series self-excited dynamo; the two-loop Rikitake mode, "
        "reversals included, is gone. The iron gain is 4.45, set so the DC onset is 1768 rpm. At 7.8 Hz each "
        "machine is a source with R<sub>s</sub> = 0.275 mΩ behind X<sub>s</sub> = 2.945 mΩ "
        "(L<sub>s</sub> = 60 µH), and the output is conjugate-matched to it with tank Q = 3.")
    r.table([
        ["Element", "Per machine", "The one tank between the coils"],
        ["R<sub>L</sub>", "63.7 mΩ", "127 mΩ"],
        ["C", "8.79 F", "4.40 F"],
        ["L<sub>p</sub> (primary)", "216.5 µH", "433 µH"],
        ["Transformer", "1 : 12.5, k = 0.98", "1 : 12.5, k = 0.98"],
        ["AC load R<sub>ac</sub>", "10 Ω", "20 Ω"],
    ], [0.3, 0.3, 0.4])
    r.figure(os.path.join(HERE, "dual_motor_schematic.png"),
             "Figure 3. Schematic of the stacked dual-motor machine with one parallel RLC tank between the two "
             "machines' output coils: motors and speed drives (counter-rotating), disks, R<sub>int</sub> 0.27 mΩ, "
             "coils 43 µH each, the tank R<sub>L</sub> 127 mΩ || C 4.4 F || L<sub>p</sub> 433 µH, a 1 : 12.5 "
             "transformer (k = 0.98), snubber, bridge, 25 mF bank and 40 Ω DC load. The secondary's common end "
             "is tied to the shafts' common. Values are a first-pass design.", max_height_mm=120)
    r.p("Parametric pumping: with the matched output the 7.8 Hz speed sine pumps the field, so it holds below "
        "the DC onset. Scaled to a 10 kA loop current:")
    r.table([
        ["h", "steady mean rpm", "loop AC", "primary AC", "R<sub>L</sub>", "secondary", "R<sub>ac</sub>"],
        ["0.1", "1758", "1043 A", "295 A", "150 W", "76 V / 3.8 A", "144 W"],
        ["0.5", "1592", "4306 A", "1218 A", "2.55 kW", "316 V / 15.8 A", "2.45 kW"],
        ["1.0", "1337", "6359 A", "1799 A", "5.56 kW", "470 V / 23.5 A", "5.34 kW"],
    ], [0.07, 0.15, 0.13, 0.14, 0.13, 0.22, 0.16])
    r.p("The circulating current rises with the tank Q, but the power that crosses the transformer does not:")
    r.table([
        ["Tank Q", "L<sub>p</sub>", "C", "primary AC", "secondary", "AC power out"],
        ["1", "1299 µH", "3.75 F", "659 A", "23.0 A", "10.8 kW"],
        ["3", "433 µH", "4.40 F", "1799 A", "23.1 A", "10.9 kW"],
        ["10", "130 µH", "6.64 F", "5833 A", "22.8 A", "10.6 kW"],
        ["30", "43.3 µH", "13.1 F", "15 176 A", "19.8 A", "8.0 kW"],
    ], [0.12, 0.16, 0.14, 0.18, 0.18, 0.22])
    r.bullets([
        "The circulating part is the primary's magnetizing current and stays in the tank; only the in-phase "
        "part reaches the secondary, at most EMF²/(8·R<sub>s</sub>) per machine.",
        "<b>Rectifier into a bank</b> (Q = 3 design, one bridge, 25 mF bank pre-charged to 400 V, 40 Ω DC load, "
        "x = 0.756, h = 1): at 8.8 kA loop current the bank holds 396 V with 33 V p-p ripple at 15.6 Hz and the "
        "DC load takes 3.9 kW; scaled to 10 kA, 449 V and 5.1 kW. Ledger closes to 4×10<super>−8</super>.",
        "Circuit = RK4 (21 666 A against 21 669 A), ledger 10<super>−6</super> of the drive work; the primary "
        "carries the loop's 10 kA DC (21.7 kJ at 433 µH), so its core needs an air gap or a cancelled DC flux.",
    ])

    r.h2("Ball bearings at the rims and an adjustable feedback")
    r.p("Each disk's rim current goes through a ball-bearing bank (P element): measured contact drop 0.73 V "
        "per bank (Safdarzadeh et al., <i>Lubricants</i> 12, 148, 2024); damage limit 0.3 A/mm² of Hertzian "
        "contact area scaled to a 300 mm rim (a 7.5× extrapolation, labelled); Huber κψ never measured, so "
        "swept. The declared bank (four 300 mm races at 5 kN) carries 24.8 A at the damage limit. Because the "
        "bearing drop is a fixed voltage, at fixed feedback the loop has one unstable equilibrium current "
        "I* = 2·V<sub>b</sub>/(2·M<sub>s</sub>·ω − R<sub>loop</sub>).")
    r.bullets([
        "<b>Copper machine:</b> no equilibrium at 1.0 or 1.5 × onset, I* = 29.5 kA at 2.0 and 4.5 kA at 2.5. At "
        "10 kA the bearings would carry 404× their damage limit (15 kW of heat). 1.6 T costs 37 kW per copper "
        "coil against 82 W generated per disk at the bearing limit.",
        "<b>Bearing-matched design:</b> superconducting field coils of 2438 turns (loop inductance 45 H), set "
        "point 19.8 A, load 0.261 Ω; each disk makes 3.32 V, each bank drops 0.73 V. At fixed feedback 5 % off "
        "the equilibrium the current runs away to 17.19 A or 22.53 A after 600 s (RK4: 17.14 A and 22.48 A).",
        "<b>PI regulator</b> on the adjustable feedback K<sub>v</sub>: from a 7.4 A start pulse it reaches 95 % "
        "in 283 s and holds 19.81 A at K<sub>v</sub> = 0, 1.60 T; load 102 W, bearings 29 W, drives 131.5 W. "
        "Books: drive 393 kJ = load 300 kJ + bearings 86 kJ + stored 7.5 kJ (residual 2×10<super>−8</super>).",
        "<b>Huber torque:</b> in EM mode the drives pay the same 0.355 N·m per shaft at κψ = 0, "
        "10<super>−4</super> and 10<super>−3</super> (at 10<super>−2</super> the current collapses); thermal "
        "mode returns at most η = 10 % of the bearing heat (≤ 0.008 N·m).",
    ])

    r.h1("Ledger and validation")
    r.bullets([
        "Tests: Lorentz EMF = formula to 10<super>−12</super>; ledger shaft power = Σ torque·ω to "
        "10<super>−9</super>; torque per disk = K·I with I = 2·EMF/R<sub>total</sub>; ledger residual below "
        "10<super>−12</super>.",
        "Dual motor: threshold growth = exact in both layouts; circuit = RK4 with feedback coils and stacked; "
        "the one tank equals the per-machine tanks (10<super>−9</super>); I<sub>2</sub> = −I<sub>1</sub> at "
        "every step.",
        "Every watt at the outputs comes from the two motors and their speed drives; the neutrino and other "
        "unmodeled channels are zero. No step adds energy: not the stacking, the resonance, the tank Q nor the "
        "transformer.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Contact, belt and brush resistances, shaft contacts, rotor inertia, windage, hoop strength, iron gain, "
        "transformer, snubber, diode, bank and load values are declared, illustrative inputs, not measured. "
        "The component values are a first pass, to be optimized later for a commercial device.",
        "The rim-bearing damage limit is a Hertz-scaled extrapolation; the cryocooler power of the "
        "superconducting coils is outside the model; the iron is linear (no saturation to stop a runaway).",
        "The bridge is nonlinear, so only the engine models it; its check is the ledger.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/homopolar [-gens 20] [-pop 30] [-p 2000]\n"
           "go run ./examples/homopolar -dual      # all dual-motor parts, A–I\n"
           "go run ./examples/homopolar -rim       # rim bearings only\n"
           "go test ./examples/homopolar\n"
           "python examples/homopolar/viz/make_report.py examples/homopolar/homopolar_report.pdf")
    r.small("Sources: N. Tesla, U.S. Patent 406,968 (1889); P. Joot, <i>Geometric Algebra for Electrical "
            "Engineers</i>, Theorem 3.16; G. Kron, <i>Equivalent Circuits of Electric Machinery</i> (1951), "
            "Epilogue eq. 9; doc/M9n.md. The full dual-motor write-up is REPORT.md in this folder.")
    r.build(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "homopolar_report.pdf"))
