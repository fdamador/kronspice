# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""M9 report: the NMOS common-source amplifier's gain / bandwidth / power Pareto front, as a PDF.

    python examples/mos_amplifier_pareto/viz/make_report.py \
        examples/mos_amplifier_pareto/mos_amplifier_pareto_report.pdf

Numbers come from examples/mos_amplifier_pareto/README.md (seed 1 results and
the seed 1-6 ranges), config.json (genes, objectives, GA settings), circuit.cir
and doc/M9.md (NSGA-II benchmarks). The example has no run_output.txt and no
figures. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/mos_amplifier_pareto", "NMOS amplifier Pareto front",
               "Trading gain, bandwidth and supply power with NSGA-II on a Level 3 NMOS stage — milestone M9 "
               "(CLAUDE.md §4, §9)",
               "CLAUDE.md §4, M9", "seed 1, NSGA-II 100 × 100")

    r.h1("Summary")
    r.p("This is the M9 acceptance example: the MOSFET levels 1–3 and the NSGA-II multi-objective search used "
        "together. A common-source NMOS stage on 5 V, driven from a 1 kΩ source into a 100 kΩ || 5 pF load, is "
        "evolved for three objectives at once — gain, bandwidth and supply power — under the constraint that the gain "
        "is at least 20 dB. No weights are traded against each other; the run returns the set of designs none of which "
        "is beaten on every objective by another (the Pareto front).")
    r.bullets([
        "<b>The front holds 99 designs</b> after 100 generations (seed 1); in generation 0 only one random design of "
        "100 met the 20 dB constraint.",
        "<b>Highest gain:</b> 34.45 dB at 1.10 MHz bandwidth and 404 µW.",
        "<b>Widest band:</b> 33.8 MHz at 20.17 dB and 19.4 mW.",
        "<b>Lowest power:</b> 22.2 µW at 24.10 dB and 344 kHz.",
        "<b>Bandwidth costs power,</b> and the gain tops out near 35 dB on every seed: the 100 kΩ load caps the drain "
        "resistance and gm/I<sub>D</sub> cannot exceed its weak-inversion value. Seeds 1–6 reach 34.5–35.2 dB, "
        "34–36 MHz at 17–21 mW, and 20–34 µW.",
        "<b>Every front netlist is re-simulated from scratch</b> by TestMOSParetoExample (constraint, recorded errors, "
        "energy ledger, mutual non-dominance).",
    ])

    r.h1("The circuit")
    r.p("circuit.cir: VDD 5 V; VIN through RSIG 1 kΩ and CIN 1 µF into the gate; the gate bias divider R1/R2; RD to "
        "the drain; RSS bypassed by CSS 100 µF at the source; M1 (drain d, gate g, source src, bulk 0, model nproc); "
        "COUT 1 µF into RL 100 kΩ || CL 5 pF; AC sweep .ac dec 10 100 1g. The starting values are R1 1 MΩ, R2 500 kΩ, "
        "RD 10 kΩ, RSS 1 kΩ and W/L = 20 µm / 2 µm.")
    r.p("The process model is SPICE <b>Level 3</b> (semi-empirical) with weak inversion (NFS), so low-power designs "
        "biased near or below threshold follow the exponential subthreshold law instead of Level 1's square law. Its "
        "parameters are derived from TOX, NSUB and UO as in SPICE3: TOX 20 nm, NSUB 3e16, UO 550, VMAX 1.5e5, "
        "THETA 0.06, ETA 0.05, KAPPA 0.4, NFS 4e11, plus overlap and junction capacitances — the same n-channel process "
        "as the mos3_* golden circuits. The Level 3 current is written in forward-mode dual numbers, so g<sub>m</sub>, "
        "g<sub>ds</sub> and g<sub>mbs</sub> are the exact derivatives of the stamped current; the gate capacitances "
        "are Meyer's.")

    r.h1("Genes, objectives and search")
    r.table([
        ["gene", "range", "series"],
        ["R1, R2 (gate bias divider)", "100 kΩ – 10 MΩ", "E24"],
        ["RD (drain load)", "100 Ω – 1 MΩ", "E24"],
        ["RSS (bypassed source resistor)", "10 Ω – 100 kΩ", "E24"],
        ["M1 W (instance parameter M1.w)", "2 – 500 µm", "continuous"],
        ["M1 L (instance parameter M1.l)", "1 – 10 µm", "continuous"],
    ], [0.45, 0.35, 0.2])
    r.table([
        ["objective", "kind", "meaning"],
        ["gain at 10 kHz, V(o)/V(in)", "gain_max, ceiling 40 dB", "more is better (error = shortfall/20)"],
        ["upper −3 dB edge", "bandwidth_max, ceiling 1 GHz", "more is better (error = decades short)"],
        ["power from VDD at the operating point", "supply_power, weight 1000 (mW)", "less is better"],
        ["gain at 10 kHz ≥ 20 dB", "gain_max with \"constraint\": true",
         "a constraint, not an objective: a design below 20 dB is dominated by every design above it"],
    ], [0.3, 0.3, 0.4])
    r.p("<b>GA settings</b> (config.json): algorithm nsga2, population 100, 100 generations, mutation 0.2, σ 0.15, "
        "seed 1, timeout 2000 ms per evaluation.")
    r.p("<b>NSGA-II</b> (Deb et al. 2002, evo/nsga2.go):")
    r.bullets([
        "<b>Ranking:</b> fast non-dominated sort (§III-A) and crowding distance (§III-B); ties are broken by index, so "
        "runs are bit-identical across worker counts and resume exactly from a checkpoint.",
        "<b>Constrained dominance (§VI):</b> feasible designs dominate violators, violators are ordered by total "
        "weighted violation, and failed simulations (singular matrix, no convergence, timeout) come last, ordered by "
        "their M5 penalty class. A failure never crashes the run.",
        "<b>Loop (§III-C):</b> binary crowded tournament, BLX-α crossover and Gaussian mutation as in the weighted GA, "
        "then the best population of parents plus offspring by rank and crowding; elitist by construction.",
        "<b>Benchmarks</b> (TestNSGA2ZDT, 10 variables, population 60, 150 generations): on ZDT1 the 60 front members "
        "lie within 0.0065 of the true front; on ZDT2 (non-convex, where a weighted sum can only find the two ends) "
        "within 0.0040 with 31 interior points.",
    ])

    r.h1("Results (seed 1)")
    r.p("Five of the 99 front designs:")
    r.table([
        ["design", "gain", "bandwidth", "power", "RD", "RSS", "W / L"],
        ["highest gain", "34.45 dB", "1.10 MHz", "404 µW", "62 kΩ", "180 Ω", "500 µm / 1.84 µm"],
        ["widest band", "20.17 dB", "33.8 MHz", "19.4 mW", "1 kΩ", "130 Ω", "499 µm / 1.22 µm"],
        ["lowest power", "24.10 dB", "344 kHz", "22.2 µW", "1 MΩ", "27 kΩ", "301 µm / 6.47 µm"],
        ["≥ 30 dB at least power", "31.41 dB", "392 kHz", "56.7 µW", "430 kΩ", "910 Ω", "427 µm / 3.43 µm"],
        ["≥ 10 MHz at least power", "23.49 dB", "11.7 MHz", "2.88 mW", "3.3 kΩ", "750 Ω", "500 µm / 1.32 µm"],
    ], [0.24, 0.11, 0.12, 0.11, 0.1, 0.1, 0.22])
    r.bullets([
        "<b>Bandwidth costs power:</b> a small RD against the 5 pF load needs current to keep the gain.",
        "<b>The gain ceiling</b> is near 35 dB because the gain is g<sub>m</sub>·(RD || RL || r<sub>o</sub>): the "
        "100 kΩ load caps the drain resistance, and g<sub>m</sub>/I<sub>D</sub> is at most its weak-inversion value "
        "1/(n·V<sub>T</sub>). The highest-gain designs use the widest device the genes allow, in moderate inversion.",
        "<b>Across seeds 1–6:</b> highest gain 34.5–35.2 dB, widest band 34–36 MHz at 17–21 mW, lowest power "
        "20–34 µW, all at gain ≥ 20 dB.",
    ])
    r.p("Besides the usual GA outputs the run writes front.csv (every front design with its three errors, the "
        "constraint and the gene values, sorted by the first objective) and front/NNN.cir, one netlist per design, each "
        "of which opens in ngspice, LTspice and QSpice as it is. No front plot is stored in the example folder.")

    r.h1("Validation")
    r.bullets([
        "<b>TestMOSParetoExample</b> (evo/mos_pareto_test.go) re-simulates every front netlist from scratch and checks "
        "the constraint, the recorded errors, the energy ledger and mutual non-dominance, and that the front spans the "
        "trade-off (a design above 34 dB, one above 30 MHz, one below 50 µW).",
        "<b>TestNSGA2ZDT</b> checks the search itself on the ZDT1 and ZDT2 benchmark fronts (numbers above).",
        "The export/ folder holds the starting circuit as .cir, .kicad_sch and a results sheet (MNA only: Kron's torn "
        "solve on that sheet is linear-only, and M1 is nonlinear).",
    ])

    r.h1("Limits")
    r.bullets([
        "One generic 1 µm-class Level 3 process; no process corners or mismatch.",
        "Gain is scored at one frequency (10 kHz) and power at the DC operating point only.",
        "Small-signal AC analysis: no distortion or output swing objective.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/mos_amplifier_pareto/config.json   # about 6 s\n"
           "go test ./evo -run 'TestMOSParetoExample|TestNSGA2'                       # acceptance and ZDT benchmarks")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "mos_amplifier_pareto_report.pdf"))
