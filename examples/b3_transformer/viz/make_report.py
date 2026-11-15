# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H52 report: the B(3) transformer, a spring as the transformer core, as a PDF.

    go run ./examples/b3_transformer -data examples/b3_transformer/viz/data > examples/b3_transformer/run_output.txt
    python3 examples/b3_transformer/viz/plot_report.py examples/b3_transformer/viz/data examples/b3_transformer/viz
    python3 examples/b3_transformer/viz/make_report.py examples/b3_transformer/b3_transformer_report.pdf

Numbers come from the run's CSV files (viz/data) and run_output.txt; the
per-core constants in the general summary are recomputed here from the same
declared geometry and materials as core.go.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401

DATA = os.path.join(HERE, "data")
MU0 = 4e-7 * math.pi


def load(name):
    with open(os.path.join(DATA, name)) as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def sci(x, d=2):
    """x as reportlab markup a×10^b."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    if 2 <= e <= 4:
        return f"{x:,.0f}"
    if -2 <= e <= 1:
        return f"{x:.{d + 1}g}"
    return f"{x / 10 ** e:.{d}g}×10<super>{'−' if e < 0 else ''}{abs(e)}</super>"


# Design geometry (core.go DefaultSpring) and declared materials.
N, R, A, P, REACH, RET = 10, 10e-3, 2e-3, 12e-3, 10e-3, 16
AREA = math.pi * A * A
LT = math.hypot(2 * math.pi * R, P)
MATS = [("nanocrystalline tape", 1e5, 1.2, "nano"), ("silicon-steel wire", 4000, 1.6, "steel"),
        ("Galfenol wire", 87, 1.5, "galfenol")]


def AL(mur):
    """Inductance per effective turn², 1/(N·R_turn + R_return) [H]."""
    rt = LT / (MU0 * mur * AREA)
    rr = (N * P + 2 * REACH) / (MU0 * mur * AREA * RET)
    return 1 / (N * rt + rr)


def band(f, g, frac=10 ** (-3 / 20)):
    """Frequencies where |ratio| stays above frac × its maximum."""
    gm = max(g)
    inside = [fi for fi, gi in zip(f, g) if gi >= frac * gm]
    return min(inside), max(inside), gm


def build(out):
    fr = load("freq.csv")
    L = load("ledger.csv")
    r = Report("examples/b3_transformer",
               "B(3) transformer: a spring as the transformer core",
               "How the turns ratio and the rate of change differ from a normal transformer when the core is an "
               "N-turn spring — register entry H52, honouring Myron Evans",
               "register entry H52")

    # ------------------------------------------------------------ summary
    r.h1("Summary")
    r.p("The core of this transformer is itself a spring: a helix of N turns closed by a return leg, so the magnetic "
        "flux runs along the core wire and goes round the spring's axis N times — the helical \"B(3)\" field of the "
        "design. The windings are wrapped round the spring's wall, and each winding turn encloses some of the core's "
        "passes. The design case has a primary of 100 turns round all 10 passes and a secondary of 100 turns round 2.")
    r.bullets([
        "<b>The turns ratio changes from N<sub>s</sub>/N<sub>p</sub> to (N<sub>s</sub>·n<sub>s</sub>)/"
        "(N<sub>p</sub>·n<sub>p</sub>).</b> Each winding turn links the flux once for every core pass it encloses. "
        "With 100:100 windings a normal transformer gives 1; the spring core gives (100·2)/(100·10) = 0.2 "
        "(magnetic network 0.200001, kronspice 0.198).",
        "<b>The rate of change is the same in every core pass, but not per winding turn.</b> In a normal transformer "
        "every turn of both windings sees the same dΦ/dt, so the volts per turn are equal. On the spring core a "
        "primary turn sees n<sub>p</sub>·dΦ/dt = 10·dΦ/dt and a secondary turn n<sub>s</sub>·dΦ/dt = 2·dΦ/dt: the "
        "volts per turn differ by n<sub>p</sub>/n<sub>s</sub> = 5 (0.200 V against 0.040 V at 20 V peak).",
        "<b>The core's turns count like winding turns:</b> with the primary round all N passes and the secondary "
        "round 2 the ratio is 2/N (N = 2 to 40, within 1×10<super>−4</super>), and per primary volt dΦ/dt falls as "
        "1/(N<sub>p</sub>·N).",
        "<b>No energy is gained.</b> Efficiency is 0.998–0.999 for every N, and the energy ledger closes to "
        "1.5×10<super>−7</super> or better: the core passes change the ratio, never the energy balance.",
        "<b>Flux must stay on the wire.</b> Wider gaps help only logarithmically; a thick return leg, a high μ<sub>r</sub> "
        "and a core wire thick compared with the spring are what keep it there. The design (nanocrystalline tape, "
        "Ø4 mm core on a Ø20 mm spring, 8 mm gaps, 16× return leg) holds every pass equal to 6×10<super>−5</super>.",
        "<b>As a magnetic amplifier</b> the same rule multiplies a control winding's ampere-turns: current gain "
        "(N<sub>c</sub>·n<sub>c</sub>)/(N<sub>g</sub>·n<sub>g</sub>), (n<sub>c</sub>/n<sub>g</sub>)× that of the same "
        "windings on an ordinary core.",
    ])

    r.h1("The general rule, for any N<sub>p</sub>, N<sub>s</sub> and core")
    r.p("Let a winding have N<sub>w</sub> turns, each enclosing n<sub>c</sub> passes of a spring core of N turns, "
        "and let Φ be the flux in every pass (no gap leakage). Then:")
    r.table([
        ["quantity", "normal transformer", "spring core"],
        ["flux linkage of a winding", "λ = N<sub>w</sub>·Φ", "λ = N<sub>w</sub>·n<sub>c</sub>·Φ"],
        ["voltage ratio v<sub>s</sub>/v<sub>p</sub>", "N<sub>s</sub>/N<sub>p</sub>",
         "(N<sub>s</sub>·n<sub>s</sub>)/(N<sub>p</sub>·n<sub>p</sub>)"],
        ["current ratio i<sub>s</sub>/i<sub>p</sub> (ideal)", "N<sub>p</sub>/N<sub>s</sub>",
         "(N<sub>p</sub>·n<sub>p</sub>)/(N<sub>s</sub>·n<sub>s</sub>)"],
        ["impedance ratio Z<sub>p</sub>/Z<sub>s</sub>", "(N<sub>p</sub>/N<sub>s</sub>)²",
         "((N<sub>p</sub>·n<sub>p</sub>)/(N<sub>s</sub>·n<sub>s</sub>))²"],
        ["rate of change in the core, per primary volt", "dΦ/dt = v<sub>p</sub>/N<sub>p</sub>",
         "dΦ/dt = v<sub>p</sub>/(N<sub>p</sub>·n<sub>p</sub>), the same in every pass"],
        ["volts per winding turn", "dΦ/dt, equal in both windings",
         "n<sub>p</sub>·dΦ/dt (primary), n<sub>s</sub>·dΦ/dt (secondary)"],
        ["peak flux density", "V<sub>pk</sub>/(ω·N<sub>p</sub>·A)", "V<sub>pk</sub>/(ω·N<sub>p</sub>·n<sub>p</sub>·A)"],
        ["magnetizing inductance", "N<sub>p</sub>²·A<sub>L</sub>", "(N<sub>p</sub>·n<sub>p</sub>)²·A<sub>L</sub>"],
        ["one wire up the bore", "—", "N·dΦ/dt (links every pass)"],
    ], [0.3, 0.28, 0.42])
    r.p("The spring transformer is therefore exactly a normal transformer with effective turns "
        "N<sub>p</sub>·n<sub>p</sub> and N<sub>s</sub>·n<sub>s</sub>: the 100:100 spring design behaves as a "
        "1000:200 normal one. The core constants for the three declared cores in the design geometry "
        "(A<sub>L</sub> = 1/(N·R<sub>turn</sub> + R<sub>return</sub>), the inductance per effective turn², from the "
        "same geometry as core.go):")
    rows = [["core", "μ<sub>r</sub>", "A<sub>L</sub> [H]", "L<sub>p</sub> for N<sub>p</sub>·n<sub>p</sub> = 1000",
             "B<sub>pk</sub> per V<sub>pk</sub> at 1 kHz", "−3 dB band of v<sub>s</sub>/v<sub>p</sub> (loaded)",
             "ratio at 1 kHz"]]
    for name, mur, bs, tag in MATS:
        al = AL(mur)
        lo, hi, gm = band(fr["f"], fr["ratio_" + tag])
        i1k = min(range(len(fr["f"])), key=lambda i: abs(fr["f"][i] - 1e3))
        rows.append([name, sci(mur, 2), sci(al, 3), f"{1e6 * al:.3g} H",
                     f"{1 / (2 * math.pi * 1e3 * 1000 * AREA) * 1e3:.3g} mT/V",
                     ("below 1 Hz" if lo <= fr["f"][0] else f"{sci(lo)} Hz") + f" – {sci(hi)} Hz",
                     f"{fr['ratio_' + tag][i1k]:.4f}"])
    r.table(rows, [0.18, 0.08, 0.12, 0.15, 0.15, 0.2, 0.12])
    r.small("Band: where |v<sub>s</sub>/v<sub>p</sub>| stays within 3 dB of its maximum with the design load of each "
            "core (its secondary's reactance at 1 kHz). The flux density per volt depends only on the effective primary "
            "turns and the core section, not on μ<sub>r</sub>.")

    # ------------------------------------------------------------ core
    r.h1("The core and the windings")
    r.figure(os.path.join(HERE, "core.png"),
             "Figure 1. The spring core: 10 turns, Ø20 mm, Ø4 mm core wire, 12 mm pitch, closed by a return leg of 16× "
             "the core section. Left: the primary (blue) wraps the spring's wall round all 10 core passes; the secondary "
             "(red) round 2. Right: section through the axis — every winding turn encloses n<sub>c</sub> core "
             "cross-sections, each carrying the same flux Φ.", width_frac=0.95)
    r.bullets([
        "<b>Magnetic network</b> (core.go): the core wire in 24 segments per turn, each a reluctance "
        "ℓ/(μ<sub>0</sub>μ<sub>r</sub>A); neighbouring turns joined at every azimuth by the air permeance of two "
        "parallel cylinders, μ<sub>0</sub>π/acosh(p/2a) per unit length; the return leg one reluctance; a winding an "
        "MMF N<sub>w</sub>·I in each pass it encloses. Solved by band Cholesky; inductances from the linkages.",
        "<b>Circuit</b> (circuit.go, kronspice): 20 V peak at 1 kHz from 0.5 Ω, the network inductances with 1 % "
        "declared winding air leakage and Ø0.5 mm copper, the secondary into a load equal to its reactance at 1 kHz "
        "(632 Ω for the design), and the bore wire into 1 MΩ.",
        "<b>Materials</b> (declared): nanocrystalline tape μ<sub>r</sub> 10<super>5</super>, B<sub>s</sub> 1.2 T; "
        "3 % Si steel μ<sub>r</sub> 4000, 1.6 T; Galfenol μ<sub>r</sub> 87, 1.5 T. Linear, no hysteresis, no eddy "
        "currents (valid for the 18 µm tape at 1 kHz, skin depth 54 µm; not for the solid steel or Galfenol wire).",
    ])

    # ------------------------------------------------------------ turns ratio
    r.h1("The turns ratio against a normal transformer")
    r.figure(os.path.join(HERE, "ratio.png"),
             "Figure 2. Left: v<sub>s</sub>/v<sub>p</sub> against the passes the secondary encloses, primary round all 10, "
             "100-turn windings: a normal 100:100 transformer stays at 1; the spring core follows n<sub>s</sub>/10. "
             "Right: against the spring's core turns N (primary round all N, secondary round 2): 2/N.")
    rr = load("ratio.csv")
    t = [["secondary passes n<sub>s</sub>", "normal N<sub>s</sub>/N<sub>p</sub>",
          "spring (N<sub>s</sub>·n<sub>s</sub>)/(N<sub>p</sub>·n<sub>p</sub>)", "network"]]
    for s, x, ra, idl, nm in zip(rr["sweep"], rr["x"], rr["ratio"], rr["ideal"], rr["normal"]):
        if s == 0 and int(x) in (1, 2, 5, 10):
            t.append([str(int(x)), f"{nm:.3f}", f"{idl:.3f}", f"{ra:.6f}"])
    r.table(t, [0.25, 0.25, 0.3, 0.2])
    t = [["core turns N", "normal", "2/N", "network", "dΦ/dt per primary volt"]]
    for s, x, ra, idl, nm in zip(rr["sweep"], rr["x"], rr["ratio"], rr["ideal"], rr["normal"]):
        if s == 1 and int(x) in (2, 4, 10, 20, 40):
            t.append([str(int(x)), f"{nm:.3f}", f"{idl:.4f}", f"{ra:.6f}", f"{1 / (100 * x) * 1e3:.3g} mWb/s"])
    r.table(t, [0.18, 0.15, 0.17, 0.2, 0.3])

    # ------------------------------------------------------------ input / output
    r.h1("Input and output, and the rate of change")
    r.figure(os.path.join(HERE, "wave.png"),
             "Figure 3. Two cycles at 1 kHz in steady state (kronspice transient). Left, spring core: primary 20 V peak, "
             "secondary 4.0 V, bore wire 0.2 V; volts per winding turn 0.200 V (primary) against 0.040 V (secondary); "
             "dΦ/dt in passes 1, 6 and 10 on top of one another (0.0198 Wb/s peak). Right, a normal 100:100 transformer "
             "with the same primary: output equals input, and the volts per turn are equal.", width_frac=1.0)
    r.p("The difference from a normal transformer is in one place: <b>how many times each winding turn links the "
        "flux</b>. The flux itself changes at the same rate everywhere in the core. A primary turn wraps 10 passes, so it "
        "sees ten times that rate; a secondary turn wraps 2 and sees twice it. A normal transformer's turns each wrap "
        "the core once, so every turn sees the rate once and the ratio can only come from the turn counts.")
    r.figure(os.path.join(HERE, "passes.png"),
             "Figure 4. dΦ/dt along the spring with the primary round all 10 passes, relative to the largest. The "
             "design holds every pass equal to 6×10<super>−5</super>; a plain spring (thin core, narrow gaps, thin return) "
             "or a lower-μ<sub>r</sub> core lets flux leak across the gaps and the passes differ.", width_frac=0.7)

    # ------------------------------------------------------------ frequency response
    r.h1("Frequency response")
    r.figure(os.path.join(HERE, "freq.png"),
             "Figure 5. Loaded AC response 1 Hz – 1 MHz (kronspice): |v<sub>s</sub>/v<sub>p</sub>|, its phase, the input "
             "impedance and the efficiency, for the three cores in the design geometry. Dashed line: 1 kHz.", width_frac=1.0)
    t = [["frequency", "|v<sub>s</sub>/v<sub>p</sub>|", "phase", "|Z<sub>in</sub>|", "η"]]
    for fq in (1, 10, 100, 1e3, 1e4, 1e5):
        i = min(range(len(fr["f"])), key=lambda k: abs(math.log(fr["f"][k] / fq)))
        t.append([f"{fr['f'][i]:.3g} Hz", f"{fr['ratio_nano'][i]:.4f}", f"{fr['phase_nano'][i]:.1f}°",
                  f"{fr['zin_nano'][i]:.3g} Ω", f"{fr['eta_nano'][i]:.3f}"])
    r.p("Design core (nanocrystalline):")
    r.table(t, [0.2, 0.2, 0.2, 0.2, 0.2])
    r.bullets([
        "<b>Mid band:</b> the ratio holds at 0.198 (the ideal 0.2 less the declared 1 % air leakage and the copper) "
        "from a few hertz to a few kilohertz, with η ≈ 0.999.",
        "<b>Low end:</b> the magnetizing reactance ω·(N<sub>p</sub>·n<sub>p</sub>)²·A<sub>L</sub> falls toward the "
        "copper and source resistance; the core's effective turns raise it by n<sub>p</sub>², so the spring core "
        "reaches lower frequencies than a 100-turn normal winding on the same core.",
        "<b>High end:</b> the declared winding leakage inductance against the load sets the roll-off (about 10 kHz for "
        "632 Ω). Winding capacitance is not modelled, so there is no self-resonance in the plot.",
        "<b>Steel and Galfenol:</b> lower A<sub>L</sub> and stronger gap leakage narrow the band; Galfenol never "
        "reaches 0.2 and its efficiency stays below 0.5.",
    ])

    # ------------------------------------------------------------ ledger
    r.h1("Energy ledger against core turns")
    r.p("The nanocrystalline spring with N core turns, primary round all N passes, secondary round 2, 100 turns each, "
        "20 V peak at 1 kHz into one fixed 632 Ω load, 20 cycles from rest. <b>E in</b> is the energy the source "
        "delivers; <b>E back</b> what the magnetizing (reactive) current returns to it in part of each cycle; "
        "<b>net in</b> = E in − E back must equal load + copper + stored, to the residual.")
    r.figure(os.path.join(HERE, "ledger.png"),
             "Figure 6. Left: the ledger's terms over the run against N. Right: efficiency, the ledger residual "
             "(×10<super>6</super>) and the peak dΦ/dt per pass.", width_frac=1.0)
    t = [["N", "ratio", "B<sub>pk</sub> [T]", "E in [mJ]", "E back", "net in", "E load", "E copper", "residual",
          "η"]]
    for i in range(len(L["N"])):
        t.append([str(int(L["N"][i])), f"{L['ratio'][i]:.4f}", f"{L['b_pk'][i]:.3f}" + (" *" if L["b_pk"][i] > 1.2 else ""),
                  f"{L['e_in'][i] * 1e3:.4f}", f"{L['e_back'][i] * 1e3:.4f}", f"{L['e_net'][i] * 1e3:.4f}",
                  f"{L['e_load'][i] * 1e3:.4f}", f"{L['e_copper'][i] * 1e3:.4f}",
                  sci(L["residual"][i] / L["e_in"][i], 1), f"{L['p_out'][i] / L['p_in'][i]:.3f}"])
    r.table(t, [0.06, 0.09, 0.1, 0.11, 0.1, 0.1, 0.1, 0.11, 0.13, 0.08])
    r.small("* above the tape's 1.2 T saturation: outside the linear model; a real core would saturate at N = 2 with "
            "this drive.")
    r.bullets([
        "<b>Every joule is accounted for:</b> net in = load + copper + stored to 1.5×10<super>−7</super> of the input "
        "or better, for every N.",
        "<b>Output power falls as 1/N²</b> at a fixed drive and load (v<sub>s</sub> = 2·v<sub>p</sub>/N); the input "
        "falls with it, and the efficiency stays 0.998–0.999. More core turns step the voltage down further; they do "
        "not change what the transformer costs or returns.",
        "<b>The reactive energy</b> (E back) peaks near N = 10, where the primary's magnetizing current is largest "
        "relative to the load current; it goes back to the source each cycle and is not lost.",
    ])

    # ------------------------------------------------------------ leakage
    r.h1("Keeping the flux off the gaps")
    r.p("The air between neighbouring turns is a second path. Two leaks are measured: the <b>spread</b> of dΦ/dt over "
        "the passes with the 10-pass primary driven, and the <b>edge leakage</b> of the 2-pass winding (its flux "
        "returning round its edges through the unwound neighbours), seen as the coupling k and the swap ratio (ideal 5). "
        "The air path is weak when the helix number H = (2πR)²/(μ<sub>r</sub>·a²·acosh(p/2a)) ≪ 1.")
    r.table([
        ["core", "spring", "H", "spread", "k", "swap"],
        ["nanocrystalline", "plain: a 1 mm, gap 1 mm, return 1×", "0.041", "2.4×10<super>−3</super>", "0.925", "4.28"],
        ["nanocrystalline", "gap 1 → 8 mm", "0.017", "2.5×10<super>−3</super>", "0.962", "4.63"],
        ["nanocrystalline", "gap 1 → 18 mm", "0.014", "3.4×10<super>−3</super>", "0.966", "4.67"],
        ["nanocrystalline", "return 1× → 16×", "0.041", "1.5×10<super>−4</super>", "0.930", "4.33"],
        ["nanocrystalline", "core wire a 1 → 2 mm", "0.014", "1.2×10<super>−3</super>", "0.971", "4.71"],
        ["<b>nanocrystalline</b>", "<b>design: a 2 mm, gap 8 mm, return 16×</b>", "<b>0.0057</b>",
         "<b>5.8×10<super>−5</super></b>", "<b>0.989</b>", "<b>4.89</b>"],
        ["silicon steel", "design", "0.14", "1.4×10<super>−3</super>", "0.814", "3.31"],
        ["Galfenol", "design", "6.6", "0.034", "0.383", "0.73"],
    ], [0.17, 0.35, 0.1, 0.16, 0.1, 0.12])
    r.p("Wider gaps help least (the gap permeance falls only as 1/acosh(p/2a)) and lengthen the return leg; the thick "
        "return leg removes the spread; the edge leakage falls with μ<sub>r</sub> and with (a/R)².")

    # ------------------------------------------------------------ mag amp
    r.h1("As a magnetic amplifier for older circuits")
    r.p("A magnetic amplifier (Mali 1960, in reference/) is a saturable reactor: a <b>gate</b> winding in series with "
        "the AC load blocks the supply while the core is unsaturated, and a small DC <b>control</b> current sets how much "
        "of each half cycle passes before the core saturates. It was the solid, tube-free amplifier of servo drives, "
        "lamp dimmers, regulated supplies, ship and aircraft controls from the 1940s to the 1960s, and it suits the "
        "staged Royer step-down of CLAUDE.md §3j. In a self-saturating amplifier the ampere-turns balance, "
        "I<sub>load</sub>·N<sub>g</sub> ≈ I<sub>c</sub>·N<sub>c</sub> (Mali pp. 34–35, 50–55).")
    r.p("<b>On the spring core each winding's ampere-turns are multiplied by the passes it encloses</b>, exactly as its "
        "volts are, so the balance becomes I<sub>load</sub>·N<sub>g</sub>·n<sub>g</sub> ≈ "
        "I<sub>c</sub>·N<sub>c</sub>·n<sub>c</sub>:")
    r.code("current gain   K_i = (N_c · n_c) / (N_g · n_g)\n"
           "power gain     K_p = K_i² · R_load / R_control\n"
           "control current to saturate   I_c = B_s · A · (N·R_turn + R_return) / (N_c · n_c)\n"
           "gate blocks up to   V_rms ≈ 4.44 · f · N_g · n_g · B_s · A")
    r.table([
        ["windings (nanocrystalline design, 60 Hz, 100 Ω load)", "I<sub>c</sub> to saturate", "gate blocks",
         "K<sub>i</sub>", "K<sub>i</sub> ordinary core", "R<sub>c</sub>", "K<sub>p</sub>"],
        ["control 100 t × 10 passes, gate 100 t × 2", "6.19 mA", "0.8 V", "5.00", "1.00", "2.07 Ω", "1.2×10<super>3</super>"],
        ["control 100 t × 10 passes, gate 100 t × 10", "6.19 mA", "4.0 V", "1.00", "1.00", "2.07 Ω", "48"],
        ["control 100 t × 2 passes, gate 100 t × 10", "31.0 mA", "4.0 V", "0.20", "1.00", "0.39 Ω", "10"],
        ["control 500 t × 10 passes, gate 50 t × 2", "1.24 mA", "0.4 V", "50.0", "10.0", "10.3 Ω", "2.4×10<super>4</super>"],
    ], [0.34, 0.12, 0.1, 0.08, 0.12, 0.1, 0.14])
    r.p("How to use it with an older circuit:")
    r.bullets([
        "<b>Control round all the passes, gate round a few.</b> A weak signal — a tube's plate current, a thermocouple, a "
        "photocell, the output of a small transistor stage — drives the control winding; its ampere-turns count "
        "n<sub>c</sub> times, so a few milliamps saturate the core and switch a load current n<sub>c</sub>/n<sub>g</sub> "
        "times larger than the same windings would on an ordinary core.",
        "<b>Two gate windings on opposite walls</b>, each with its own rectifier, make the full-wave self-saturating "
        "amplifier; the control and bias windings share the passes, so one spring carries all of them.",
        "<b>Size the gate for the supply.</b> This small design core blocks only 0.8 V at 60 Hz with a 100 × 2 gate. A "
        "120 V, 60 Hz load needs N<sub>g</sub>·n<sub>g</sub> ≈ 30,000 effective turns at 1.2 T — for example "
        "1,500 turns round 20 passes — or a core of larger section (the blocking voltage scales with A).",
        "<b>Core material:</b> use a square-loop nanocrystalline or permalloy tape for a sharp switching knee; reset the "
        "core each cycle (bias winding) when it is driven by one-sign pulses (§3j).",
        "<b>It adds no energy.</b> The load power comes from the AC supply through the gate; the control winding only "
        "steers it. The spring gives the gain of an ordinary core wound N<sub>c</sub>·n<sub>c</sub> and "
        "N<sub>g</sub>·n<sub>g</sub> turns, set with one physical winding each.",
        "<b>Not yet simulated:</b> this example's core is linear. The saturating amplifier needs the hysteretic core "
        "model of §3j (device/core.go, device/magamp.go, milestone M9r), which is the next step.",
    ])

    # ------------------------------------------------------------ B(3)
    r.h1("The helical field and B(3) (H52)")
    r.p("Myron Evans' B(3) field (<i>Physica B</i> 182, 227, 1992; AIAS Essays 30, 34–39 in reference/aias-b3/) is the "
        "longitudinal magnetic field that, in his theory, circularly polarized radiation carries, tied to the transverse "
        "fields by the B Cyclic theorem B(1)×B(2) = i·B(0)·B(3)*. The user's design honours that work: in the spring core "
        "the flux follows a helix whose transverse part turns once per pitch, a magnetic circular polarization laid out "
        "in space.")
    r.table([
        ["at 0.25 T in the wire (design drive)", "value"],
        ["transverse part B(0) = B·cos α (pitch angle α = 10.8°)", "0.249 T"],
        ["|B(3)| by the B Cyclic theorem (= B(0))", "0.249 T"],
        ["axial field in the wire, B·sin α", "0.0475 T"],
        ["|B(1)×B(2)| = B(0)², standard reading (an axial vector, not a field)", "0.062 T²"],
        ["net axial flux through any cut of the spring", "Φ (returned by the leg)"],
    ], [0.7, 0.3])
    r.bullets([
        "<b>Measurement that separates the two parts:</b> a wire up the bore reads the transverse circulation, "
        "N·dΦ/dt; a coil round the outside of the spring reads the axial flux, once.",
        "<b>Standard prediction (H52):</b> the electric field does not cancel — curl E = −∂B/∂t is the windings' EMF. At "
        "1 kHz the field is quasi-static (wavelength 300 km against a 120 mm spring): a spatial helix oscillating in "
        "place, linearly polarized in time at each point, so B(1)×B(2) = 0 there.",
        "<b>Next step toward the user's emulation:</b> two interleaved helices (a bifilar spring core) driven in "
        "quadrature make the transverse field rotate in time; a fluxgate in the bore then looks for a static axial field, "
        "which the standard prediction puts at zero in vacuum.",
    ])
    r.p("<b>Status: untested.</b> The standard prediction is computed; no bench data yet.")

    r.h1("Limits of the model")
    r.bullets([
        "Linear core: no saturation, no hysteresis, no eddy currents (the N = 2 ledger row exceeds the tape's 1.2 T).",
        "Air leakage only between neighbouring turns at the same azimuth (two-cylinder formula, an estimate for "
        "p ≈ R); none from the spring's outside or between distant turns; winding air leakage declared at 1 %.",
        "No winding capacitance, so no self-resonance above the band.",
        "Materials declared from data sheets, not fitted to a bench core.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/b3_transformer -data examples/b3_transformer/viz/data > examples/b3_transformer/run_output.txt\n"
           "python3 examples/b3_transformer/viz/plot_report.py examples/b3_transformer/viz/data examples/b3_transformer/viz\n"
           "python3 examples/b3_transformer/viz/core.py examples/b3_transformer/viz/core.png\n"
           "python3 examples/b3_transformer/viz/make_report.py examples/b3_transformer/b3_transformer_report.pdf\n"
           "go test ./examples/b3_transformer")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "b3_transformer_report.pdf"))
