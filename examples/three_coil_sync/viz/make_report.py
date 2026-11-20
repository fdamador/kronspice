# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Three coils, two drives: the PDF report.

    go run ./examples/three_coil_sync -data examples/three_coil_sync/viz/data_free > examples/three_coil_sync/run_free.txt
    go run ./examples/three_coil_sync -sync -data examples/three_coil_sync/viz/data > examples/three_coil_sync/run_output.txt
    go run ./examples/three_coil_sync -waves -data examples/three_coil_sync/viz/data -data-free examples/three_coil_sync/viz/data_free
    python3 examples/three_coil_sync/viz/plot_report.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz
    python3 examples/three_coil_sync/viz/plot_report.py examples/three_coil_sync/viz/data_free examples/three_coil_sync/viz/free
    python3 examples/three_coil_sync/viz/schematic.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz
    python3 examples/three_coil_sync/viz/schematic.py examples/three_coil_sync/viz/data_free examples/three_coil_sync/viz/free
    python3 examples/three_coil_sync/viz/make_report.py examples/three_coil_sync/three_coil_sync_report.pdf

Every number is read from the runs' CSV files (viz/data: synchronized design;
viz/data_free: free-running design) and the run logs.
"""
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report  # noqa: E402


def rows(d, name):
    with open(os.path.join(HERE, d, name)) as f:
        return list(csv.DictReader(f))


def design(d):
    return {r["name"]: float(r["value"]) for r in rows(d, "design.csv")}


def eng(x, unit, d=3):
    if x == 0:
        return f"0 {unit}"
    for p, s in ((1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n"), (1e-12, "p")):
        if abs(x) >= p * 0.9999:
            return f"{x / p:.{d}g} {s}{unit}"
    return f"{x:.{d}g} {unit}"


def ga_stats(log):
    run = open(os.path.join(EX, log), encoding="utf-8").read()
    pop, gens, seed = re.search(r"population (\d+), (\d+) generations, seed (\d+)", run).groups()
    nevals, wall = re.search(r"(\d+) evaluations in (\S+)", run).groups()
    return pop, gens, seed, nevals, wall


def build(out):
    A, B = design("data"), design("data_free")  # A synchronized, B free-running
    cA = {r["case"]: r for r in rows("data", "compare.csv")}
    cB = {r["case"]: r for r in rows("data_free", "compare.csv")}
    chA, chB = rows("data", "check.csv")[0], rows("data_free", "check.csv")[0]
    tA, tB = rows("data", "trends.csv"), rows("data_free", "trends.csv")
    gA, gB = ga_stats("run_output.txt"), ga_stats("run_free.txt")

    def p(c, k, key="p_out"):
        return float(c[k][key])

    royA, hvA = p(cA, "evolved Royer only"), p(cA, "evolved HV only")
    royB, hvB = p(cB, "evolved Royer only"), p(cB, "evolved HV only")
    phA = [float(q["p_out"]) for q in tA if q["sweep"] == "phase"]
    detB = [q for q in tB if q["sweep"] == "detune"]
    lockB = [q for q in detB if q["locked"] == "true"]
    best_lock = max(lockB, key=lambda q: float(q["p_out"])) if lockB else None

    r = Report("examples/three_coil_sync",
               "Three coils, two drives: high voltage and high current combined in one resonant coil",
               "A high-voltage coil wired to a resonant output coil and a Royer oscillator's coil wound over it, "
               "evolved by kronspice's genetic algorithm for the most power in the load",
               "CLAUDE.md §3e, §3j, §4")

    # ------------------------------------------------------------ summary
    r.h1("Summary")
    r.p("Coil 1 is an air-core step-up transformer; its high-voltage terminal is wired to coil 3. Coil 2 is the "
        "bifilar, centre-tapped tank coil of a zero-voltage-switching (resonant) Royer oscillator on 24 V, wound over "
        "coil 3. Coil 3 is a series-resonant RLC loop feeding the load. Coil 3 is therefore driven twice: by voltage "
        "through the wire from coil 1 and by current through the magnetic coupling from coil 2. Every inductance and "
        "coupling is computed from the winding geometry, and kronspice ran every candidate as a full transient with "
        "the MOSFETs and diodes in it and the energy ledger on. Two searches were run, each 48 designs × 25 "
        "generations:")
    r.bullets([
        f"<b>A — synchronized (recommended).</b> The HV amplifier is clocked from a pickup on coil 2 through a phase "
        f"shifter, so both drives share one frequency by construction; the GA chose the phase and gain. Load "
        f"<b>{A['POut']:.0f} W</b> from {A['PInRoy']:.0f} W (Royer) + {A['PInHV']:.0f} W (HV amplifier), "
        f"<b>η = {A['Eta']:.3f}</b>, steady (output ripple {A['Beat'] * 100:.1f} %), at "
        f"{eng(A['FRoyer'], 'Hz', 4)}. Coil 3's current is set mostly by the Royer; the HV drive, phased at "
        f"{A['SyncAngle']:.0f}°, adds {A['POut'] - royA:.0f} W over the Royer alone ({royA:.0f} W).",
        f"<b>B — free-running.</b> An independent HV generator; the GA chose its frequency and phase. Load "
        f"<b>{B['POut']:.0f} W</b> from {B['PInHV']:.0f} W (HV) + {B['PInRoy']:.0f} W (Royer), η = {B['Eta']:.3f}. The "
        f"HV generator ({eng(B['FHV'], 'Hz', 4)}) and the Royer ({eng(B['FRoyer'], 'Hz', 4)}) are not locked: the "
        f"output beats at their difference and swings ±{B['Beat'] * 50:.0f} % within each beat. "
        f"{B['loss_Rg']:.0f} W of the HV input is lost in the generator's own 1 Ω source resistance.",
        ("<b>Locking B.</b> The Royer injection-locks to the HV generator for HV frequencies "
         f"{float(lockB[0]['x']):.2f}–{float(lockB[-1]['x']):.2f}× B's value; locked at {float(best_lock['x']):.2f}× it gives "
         f"{float(best_lock['p_out']):.0f} W, steady, but draws {float(best_lock['p_hv']):.0f} W from the HV generator, "
         "past its 400 W rating, which is why the GA stopped at the edge of the lock range." if best_lock else
         "<b>Locking B.</b> No lock was found in ±15 % of the HV frequency."),
        f"<b>The GA did not keep coil 1 as a high-voltage coil.</b> Both searches chose a low step-up "
        f"({int(A['N1s'])} and {int(B['N1s'])} secondary turns), so the HV terminal reaches only {A['VHVPk']:.0f} V (A) and "
        f"{B['VHVPk']:.0f} V (B) peak. A high-voltage, low-current drive cannot push much current through coil 3's "
        "loop, whose impedance at resonance is a few tens of ohms; the load wants current, so evolution traded voltage "
        "for current.",
        f"<b>Phase decides the result.</b> On design A the phase between the two drives alone moves the output from "
        f"{min(phA):.0f} to {max(phA):.0f} W. In the hand-picked starting design the drives were in the wrong phase: "
        f"both on gave {p(cA, 'baseline both'):.2f} W against {p(cA, 'baseline HV only'):.2f} W from the HV coil alone, "
        "and the Royer's supply took power back.",
        f"<b>No surplus anywhere.</b> The load never gets more than the two supplies deliver. The transient ledger "
        f"closes to {abs(float(chA['mismatch'])):.0e} (A) and {abs(float(chB['mismatch'])):.0e} (B) of the input "
        f"energy; a 4× finer time step moves the load power by {abs(float(chA['rel'])) * 100:.2f} % (A) and "
        f"{abs(float(chB['rel'])) * 100:.1f} % (B, whose beating is sensitive to the step).",
    ])

    # ------------------------------------------------------------ circuit
    r.h1("The circuit")
    r.figure(os.path.join(HERE, "schematic.png"),
             "Figure 1. Design A as kronspice simulates it, with the evolved values. In design B the HV amplifier is a "
             "free-running 48 V generator (Figure 9).", width_frac=1.0)
    r.bullets([
        "<b>Coil 1 (HV):</b> generator or amplifier with 1 Ω source resistance into an air-core transformer, the primary "
        "over the grounded end of a long single-layer secondary; the top of the secondary is the HV terminal, wired to "
        "the top of coil 3. Its self-capacitance (Medhurst) is included.",
        "<b>Coil 3 (output):</b> L3, its copper resistance R3, the resonating capacitor C3 and the load RL in series "
        "from the HV terminal to ground.",
        "<b>Coil 2 (Royer):</b> the Mazzilli ZVS form of the resonant Royer: two MOSFETs, a 100 µH feed choke to the "
        "centre tap, C2 across the tank, gates pulled up through 470 Ω from 12 V and pulled down by fast diodes to the "
        "opposite drains. Its frequency is set by the tank, 1/(2π√(L<sub>tank</sub>·C2)), pulled by the load on "
        "coil 3. The tank coil is two strands wound side by side and joined in the middle, so its halves couple "
        "tightly (k ≈ 0.96).",
        "<b>Synchronization (A):</b> the tank voltage is sensed (in hardware, a one-turn pickup loop on coil 2), "
        "shifted by a first-order all-pass (0–180°) and an optional inversion (the other 180°), and amplified with "
        "gain G into coil 1's primary. The amplifier's power comes from its own supply and is counted as input.",
        "<b>Models:</b> level-1 MOSFETs (V<sub>T</sub> 4 V, R<sub>DS</sub> ≈ 30 mΩ) with 3 nF gate and 1 nF output "
        "capacitance; diodes with 0.1 Ω series resistance; copper with skin effect, no proximity effect.",
        "<b>Start-up:</b> a symmetric ZVS can sit with both MOSFETs on; one gate resistor is 2 % high (a part "
        "tolerance), so the oscillator starts the same way every run.",
    ])

    # ------------------------------------------------------------ designs
    r.h1("The two evolved designs")

    def two(label, ka, kb=None, f=lambda x: f"{x:.3g}"):
        kb = kb or ka
        return [label, f(A[ka]), f(B[kb])]

    i = lambda x: f"{int(x)}"
    r.table([
        ["", "A — synchronized", "B — free-running"],
        two("coil 3 turns N3", "N3", f=i),
        two("coil 2 turns per strand N2", "N2", f=i),
        two("spacer coil 3 → coil 2", "gap23", f=lambda x: f"{x * 1e3:.1f} mm"),
        two("coil 1 secondary / primary turns", "N1s", f=i)[:1] + [f"{int(A['N1s'])} / {int(A['N1p'])}",
                                                                   f"{int(B['N1s'])} / {int(B['N1p'])}"],
        two("Royer tank capacitor C2", "C2", f=lambda x: eng(x, "F")),
        two("output capacitor C3", "C3", f=lambda x: eng(x, "F")),
        two("load RL", "RL", f=lambda x: eng(x, "Ω")),
        ["HV drive", f"slaved, {A['SyncAngle']:.0f}°, gain {A['SyncGain']:.3f} ({A['VgPk']:.0f} V pk)",
         f"{eng(B['FHV'], 'Hz', 4)}, 48 V pk"],
        two("L3", "L3", f=lambda x: eng(x, "H")),
        two("L2 per strand / tank", "L2strand")[:1] + [f"{eng(A['L2strand'], 'H')} / {eng(A['Ltank'], 'H')}",
                                                       f"{eng(B['L2strand'], 'H')} / {eng(B['Ltank'], 'H')}"],
        two("L1s, L1p, k1", "L1s")[:1] + [f"{eng(A['L1s'], 'H')}, {eng(A['L1p'], 'H')}, {A['K1']:.2f}",
                                         f"{eng(B['L1s'], 'H')}, {eng(B['L1p'], 'H')}, {B['K1']:.2f}"],
        two("coupling k(2–3)", "K23", f=lambda x: f"{x:.3f}"),
        two("Royer frequency, measured", "FRoyer", f=lambda x: eng(x, "Hz", 4)),
    ], [0.34, 0.33, 0.33])
    r.table([
        ["steady state", "A", "B", "rating"],
        ["load power", f"{A['POut']:.1f} W", f"{B['POut']:.1f} W", ""],
        ["HV input / Royer input", f"{A['PInHV']:.1f} / {A['PInRoy']:.1f} W", f"{B['PInHV']:.1f} / {B['PInRoy']:.1f} W",
         "≤ 400 W / 24 V × 20 A"],
        ["efficiency η", f"{A['Eta']:.3f}", f"{B['Eta']:.3f}", ""],
        ["output ripple (per cycle, max − min / mean)", f"{A['Beat']:.3f}", f"{B['Beat']:.3f}", ""],
        ["MOSFET peak drain voltage", f"{A['VdsPk']:.0f} V", f"{B['VdsPk']:.0f} V", "≤ 150 V"],
        ["supply current (mean)", f"{A['ISupply']:.1f} A", f"{B['ISupply']:.1f} A", "≤ 20 A"],
        ["coil 2 strand current", f"{A['I2rms']:.1f} A rms", f"{B['I2rms']:.1f} A rms", "≤ 30 A rms"],
        ["coil 3 current", f"{A['I3rms']:.2f} A rms", f"{B['I3rms']:.2f} A rms", "≤ 8 A rms"],
        ["HV terminal peak", f"{A['VHVPk']:.0f} V", f"{B['VHVPk']:.0f} V", "≤ 20 kV"],
        ["C3 peak", f"{A['VC3Pk']:.0f} V", f"{B['VC3Pk']:.0f} V", "≤ 2 kV"],
        ["load peak voltage", f"{A['VOutPk']:.0f} V", f"{B['VOutPk']:.0f} V", ""],
    ], [0.36, 0.2, 0.2, 0.24])
    r.small("Both designs sit at their binding limits: A at the 20 A supply and 30 A strand current, B at the supply "
            "and close to the HV generator's 400 W. A larger supply or thicker coil 2 raises both.")

    # ------------------------------------------------------------ build
    r.h1("How to build the coils (design A)")
    r.figure(os.path.join(HERE, "coils.png"),
             "Figure 2. Design A to scale. Left: coil 3 on its former, the spacer, and coil 2's two strands "
             "(alternating shades) over it. Right: coil 1, secondary close-wound with the primary over its grounded "
             "end.", width_frac=1.0)
    l3 = A["coil_3_length"] * 1e3
    l2 = 2 * A["N2"] * 2.2
    r.bullets([
        f"<b>Coil 3 (output).</b> {int(A['N3'])} turns of Ø1.0 mm enamelled copper (AWG 18), close-wound at 1.1 mm "
        f"pitch on a Ø50 mm PVC or PTFE tube: winding length {l3:.0f} mm, about {A['coil_3_wire_length']:.0f} m of "
        f"wire. Anchor both ends through holes in the former. Target {eng(A['L3'], 'H')}; check with an LCR meter "
        "before winding coil 2.",
        f"<b>Spacer.</b> Build the outside of coil 3 up by {A['gap23'] * 1e3:.1f} mm in radius with PTFE or polyimide "
        f"tape, or slide on a tube. It also insulates coil 2 from coil 3, which reaches {A['VHVPk']:.0f} V peak at its "
        "top.",
        f"<b>Coil 2 (Royer tank).</b> Two Ø2.0 mm enamelled wires (AWG 12) side by side, {int(A['N2'])} turns of the "
        f"pair over the spacer, centred along coil 3 (winding length {l2:.0f} mm). Join the end of strand A to the "
        "start of strand B: that joint is the centre tap, to the choke. The free start of A goes to M1's drain (c1), "
        f"the free end of B to M2's drain (c2). The tank carries {A['I2rms']:.0f} A rms per strand: litz wire or "
        f"Ø6 mm copper tube runs cooler and changes the inductance little. Target {eng(A['L2strand'], 'H')} per "
        f"strand, {eng(A['Ltank'], 'H')} end to end.",
        f"<b>Coil 1 (HV).</b> Close-wind {int(A['N1s'])} turns on a Ø75 mm tube "
        f"(modelled with Ø0.4 mm wire, winding length {A['coil_1_secondary_length'] * 1e3:.0f} mm; Ø1.0 mm wire makes it "
        f"{A['N1s'] * 1.1:.0f} mm long), ground its bottom end, and wind the primary, "
        f"{int(A['N1p'])} turns of Ø2.0 mm wire at 3 mm pitch, on a Ø95 mm ring over the grounded end. In design A "
        f"coil 1 is a low-ratio transformer ({eng(A['L1s'], 'H')} : {eng(A['L1p'], 'H')}); in design B it has "
        f"{int(B['N1s'])} : {int(B['N1p'])} turns. <b>Use Ø1.0 mm wire for the secondary, not the Ø0.4 mm the model "
        f"assumed:</b> it carries coil 3's full current ({A['I3rms']:.1f} A rms), which is "
        f"{A['I3rms'] / (math.pi * 0.2 ** 2):.0f} A/mm² in Ø0.4 mm wire and would overheat it. Thicker wire only "
        f"lowers the loss (the model's {A['loss_R1s']:.0f} W in it is pessimistic) and changes the inductance little.",
        "<b>Pickup loop (sync).</b> One turn of insulated wire round coil 2, into a comparator or PLL and an "
        f"adjustable phase shifter that clocks the HV amplifier; set the phase to {A['SyncAngle']:.0f}° of the tank "
        "voltage, and trim it on the bench for the largest load power (Figure 6).",
        "<b>Layout.</b> Keep coil 1 at least three diameters from coils 2–3, axes at right angles: the model has no "
        "coupling between them. Keep the Royer's high-current loop (C2, MOSFETs, coil 2's leads) short and wide.",
    ])

    # ------------------------------------------------------------ BOM
    r.h1("Bill of materials (design A)")
    r.table([
        ["ref", "part", "value / rating", "qty", "notes"],
        ["coil 3", "enamelled copper wire Ø1.0 mm (AWG 18)", f"{A['coil_3_wire_length'] * 1.2:.0f} m", "1", "20 % spare"],
        ["", "former, PVC/PTFE tube", "Ø50 mm OD", "1", f"≥ {l3 + 40:.0f} mm long"],
        ["coil 2", "enamelled copper wire Ø2.0 mm (AWG 12)", f"2 × {A['coil_2A_wire_length'] * 1.3:.1f} m", "2",
         "or litz / Ø6 mm copper tube"],
        ["", "spacer: PTFE / polyimide tape or tube", f"{A['gap23'] * 1e3:.1f} mm radial", "1", "≥ 1 kV"],
        ["coil 1", "enamelled copper wire Ø1.0 mm (AWG 18)", f"{A['coil_1_secondary_wire_length'] * 1.2:.0f} m", "1",
         "secondary (carries coil 3's current)"],
        ["", "former, PVC/acrylic tube", "Ø75 mm OD", "1", f"≥ {A['coil_1_secondary_length'] * 1e3 + 40:.0f} mm long"],
        ["", "wire Ø2.0 mm, PTFE or enamelled", f"{A['coil_1_primary_wire_length'] * 1.3:.1f} m", "1", "primary"],
        ["C2", "polypropylene film (MKP), high-current", f"{eng(A['C2'], 'F')}, ≥ 250 V", "1 set",
         f"e.g. 10 × 470 nF in parallel to share ~{A['I2rms'] * 2:.0f} A rms"],
        ["C3", "polypropylene film (MKP / pulse)", f"{eng(A['C3'], 'F')}, ≥ {max(250, math.ceil(A['VC3Pk'] * 1.5 / 50) * 50)} V",
         "1", f"rated for {A['I3rms']:.0f} A rms"],
        ["RL", "load (power resistor for testing)", f"{eng(A['RL'], 'Ω')}, ≥ {math.ceil(A['POut'] * 1.5 / 50) * 50} W", "1",
         "non-inductive"],
        ["M1, M2", "N-MOSFET IRFP260N", "200 V, 50 A, TO-247", "2", "on a heatsink"],
        ["D1, D2", "fast diode UF4007", "1 kV, 1 A", "2", ""],
        ["RG1, RG2", "resistor", "470 Ω, 1 W", "2", ""],
        ["ZD1, ZD2", "zener (gate protection)", "15 V, 1 W", "2", "not in the model"],
        ["Lchoke", "toroidal choke", "100 µH, ≥ 25 A", "1", ""],
        ["PSU1", "DC supply", "24 V, 20 A", "1", "Royer"],
        ["PSU2", "regulator", "12 V, 0.5 A (7812)", "1", "gate rail"],
        ["AMP1", "class-D bridge or audio power amplifier", f"≥ {max(20, math.ceil(A['VgPk'] * 2))} V pk, ≥ 50 W, "
         f"10–30 kHz", "1", "design B: 48 V pk, 400 W"],
        ["SYNC", "pickup loop, comparator (LM311), RC phase shifter or PLL (CD4046)", "", "1",
         "clocks AMP1 from coil 2"],
        ["", "fuses 2 × 25 A slow, HV wire, standoffs", "", "", ""],
    ], [0.09, 0.34, 0.24, 0.07, 0.26])
    r.small("Wire lengths come from the evolved geometry plus spare. Prices and suppliers are left out on purpose.")

    # ------------------------------------------------------------ results A
    r.h1("Results — design A (synchronized)")
    r.figure(os.path.join(HERE, "compare.png"),
             "Figure 3. Mean powers in steady state for design A and for the hand-picked baseline: what each source "
             "delivers and what reaches the load. \"HV only\" drives coil 1 with a sine at the same frequency and "
             "amplitude as when slaved, Royer off.", width_frac=1.0)
    r.p(f"In design A the Royer does most of the work: alone it puts {royA:.0f} W in the load. Driving coil 1 alone "
        f"gives only {hvA:.1f} W, because coil 1's path is lossy at this frequency. Together, in phase, the load gets "
        f"{A['POut']:.0f} W: more than the sum of the two alone ({royA + hvA:.0f} W) and close to the coherent sum "
        f"(√P₁ + √P₂)² = {(math.sqrt(royA) + math.sqrt(hvA)) ** 2:.0f} W, because coil 3 carries one current and its "
        "power goes as the square of it. The extra power is paid by both sources: the Royer delivers more when the HV "
        "drive pushes the same way.")
    r.figure(os.path.join(HERE, "wave.png"),
             "Figure 4. Design A, five periods of the steady state: voltages at the HV terminal, the Royer tank and the "
             "load; currents in coil 3, coil 2 (circulating) and the HV amplifier; the voltage across C3.",
             width_frac=1.0)
    r.figure(os.path.join(HERE, "startup.png"),
             "Figure 5. Design A from switch-on: peak tank and load voltage per 10 µs. Averages are taken over the "
             "last 0.5 ms.", width_frac=0.85)
    r.figure(os.path.join(HERE, "trends.png"),
             "Figure 6. Design A, one parameter at a time: HV amplifier gain, the phase between the drives, load, C3, "
             "and the spacer between coils 3 and 2 with the coupling it gives.", width_frac=1.0)
    gain = [q for q in tA if q["sweep"] == "gain"]
    gp = [q for q in tA if q["sweep"] == "gap"]
    gpp = max(gp, key=lambda q: float(q["p_out"]))
    r.bullets([
        f"<b>Gain.</b> More HV drive gives more output ({float(gain[-1]['p_out']):.0f} W at 1.5× the evolved gain), "
        f"but the Royer's supply current rises with it; at the evolved gain it is already {A['ISupply']:.1f} A of the "
        "20 A limit. At low gain the HV amplifier absorbs power (negative HV input): the phase that suits the Royer "
        "makes coil 1 a load until the HV drive is strong enough.",
        f"<b>Phase.</b> The output runs from {min(phA):.0f} to {max(phA):.0f} W round the circle; the best phase is "
        "where the HV drive's contribution to coil 3's current lines up with the Royer's. This is the bench "
        "adjustment that matters most.",
        "<b>Load and C3</b> peak near the evolved values: C3 tunes coil 3's loop to the Royer frequency, and the load "
        "is matched to that loop.",
        f"<b>Spacer.</b> Coupling falls from {float(gp[0]['k23']):.2f} at 1 mm to {float(gp[-1]['k23']):.2f} at 20 mm; "
        f"the output peaks at {float(gpp['x']):.1f} mm. Tighter is not better: a tighter coupling loads the Royer and "
        "pulls its frequency.",
    ])
    r.figure(os.path.join(HERE, "losses.png"), "Figure 7. Design A: mean power into each part.", width_frac=0.85)

    # ------------------------------------------------------------ results B
    r.h1("Results — design B (free-running)")
    r.figure(os.path.join(HERE, "free", "compare.png"),
             "Figure 8. Design B: each drive alone and both together.", width_frac=1.0)
    r.p(f"Alone, the HV coil gives {hvB:.0f} W and the Royer {royB:.0f} W; together {B['POut']:.0f} W, against "
        f"{hvB + royB:.0f} W for the plain sum and {(math.sqrt(hvB) + math.sqrt(royB)) ** 2:.0f} W for the coherent "
        "sum. The two are near enough in frequency to pull on each other but not locked, so the output beats "
        f"(ripple {B['Beat']:.2f}). The HV side is lossy: of {B['PInHV']:.0f} W from the generator, "
        f"{B['loss_Rg']:.0f} W heat its 1 Ω source resistance and {B['loss_R1s'] + B['loss_R1p']:.0f} W coil 1.")
    r.figure(os.path.join(HERE, "free", "schematic.png"), "Figure 9. Design B with its evolved values.",
             width_frac=1.0)
    r.figure(os.path.join(HERE, "free", "trends.png"),
             "Figure 10. Design B, one parameter at a time (shaded where the Royer locks to the HV generator).",
             width_frac=1.0)
    r.figure(os.path.join(HERE, "free", "losses.png"), "Figure 11. Design B: mean power into each part.",
             width_frac=0.85)

    # ------------------------------------------------------------ GA
    r.h1("The genetic search")
    r.figure(os.path.join(HERE, "convergence.png"),
             f"Figure 12. Design A's search: best score per generation (population {gA[0]}, {gA[1]} generations, "
             f"seed {gA[2]}; {gA[3]} transient simulations in {gA[4]}). Design B's search: {gB[3]} simulations in "
             f"{gB[4]}, seed {gB[2]}.", width_frac=0.8)
    r.bullets([
        "<b>Genes (10):</b> N3, N2, the spacer, N1s, N1p, C2 and C3 (E12), RL (E24); then for A the drive phase "
        "(0–360°) and amplifier gain, for B the HV frequency (as a ratio of the Royer's unloaded tank frequency) and "
        "its starting phase.",
        "<b>Fitness:</b> 100/(100 + P<sub>load</sub>) plus a penalty proportional to any rating exceeded (drain "
        "150 V, supply 20 A, strand 30 A rms, coil 3 8 A rms, C3 2 kV, HV terminal 20 kV, HV input 400 W) and a "
        "penalty when the Royer does not oscillate. The GA cannot buy power by overrating a part.",
        "<b>Each candidate</b> is a 2.5 ms transient (switch-on, settling, 0.5 ms measured) at ≤ 0.2 µs steps with "
        "the ledger on. Seeded, so each run repeats bit for bit.",
    ])

    # ------------------------------------------------------------ physics
    r.h1("What the physics allows")
    r.bullets([
        "<b>Energy:</b> every watt in the load comes from the 24 V supply or the HV amplifier, and the ledger closes. "
        "High voltage and high current are not two energy sources inside coil 3; they are two ways of pushing one "
        "current round one loop.",
        "<b>Fields:</b> coil 1 drives coil 3 through the wire (its electric field), coil 2 through the shared flux "
        "(its magnetic field). Round the loop both are EMFs in series: what adds is EMF, and so current, which is why "
        "phase matters and why the load power can exceed the sum of the two alone while the inputs rise with it.",
        "<b>Synchronization</b> is what turns two drives into one: the same frequency and a phase that adds. A "
        "free-running Royer follows the HV drive only inside its lock range; clocking the amplifier from the Royer "
        "removes that limit.",
    ])
    r.h1("Next steps")
    r.bullets([
        "<b>Bench check:</b> measure L3, L2 (per strand and end to end), k(2–3), L1s and L1p with an LCR meter and "
        "rerun the example with the measured values; measure input and output power with a power analyser.",
        "<b>Rerun with coil 1's wire as a gene</b> (it is fixed at Ø0.4 mm here) and a rating on its current.",
        "<b>Raise the limits:</b> a 24 V, 40 A supply and litz coil 2 lift design A's binding constraints; rerun the "
        "search with the new ratings.",
        "<b>Model limits:</b> no proximity effect in the copper, air-core coils only, no stray coupling between coil 1 "
        "and coils 2–3, level-1 MOSFETs, an ideal phase shifter and amplifier apart from the 1 Ω source resistance, "
        "no radiation (all coils far below λ/10).",
    ])
    r.h1("Reproducing the results")
    r.code("go run ./examples/three_coil_sync -data examples/three_coil_sync/viz/data_free > examples/three_coil_sync/run_free.txt\n"
           "go run ./examples/three_coil_sync -sync -data examples/three_coil_sync/viz/data > examples/three_coil_sync/run_output.txt\n"
           "go run ./examples/three_coil_sync -waves -data examples/three_coil_sync/viz/data \\\n"
           "    -data-free examples/three_coil_sync/viz/data_free > examples/three_coil_sync/check_output.txt\n"
           "python3 examples/three_coil_sync/viz/plot_report.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz\n"
           "python3 examples/three_coil_sync/viz/plot_report.py examples/three_coil_sync/viz/data_free examples/three_coil_sync/viz/free\n"
           "python3 examples/three_coil_sync/viz/schematic.py examples/three_coil_sync/viz/data examples/three_coil_sync/viz\n"
           "python3 examples/three_coil_sync/viz/schematic.py examples/three_coil_sync/viz/data_free examples/three_coil_sync/viz/free\n"
           "python3 examples/three_coil_sync/viz/make_report.py examples/three_coil_sync/three_coil_sync_report.pdf\n"
           "go test ./examples/three_coil_sync\n"
           "go run ./cmd/kronspice examples/three_coil_sync/viz/data/best.cir     # design A as a netlist")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(EX, "three_coil_sync_report.pdf"))
