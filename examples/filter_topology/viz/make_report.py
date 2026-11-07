# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""M8 report: a low-pass filter grown from an empty embryo, as a PDF.

python examples/filter_topology/viz/make_report.py examples/filter_topology/filter_topology_report.pdf

Numbers come from examples/filter_topology/README.md, config.json, embryo.cir
and doc/M8.md (the example has no run_output.txt and no figures). Layout from
internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401
import kronreport  # noqa: E402

kronreport.H1.keepWithNext = 1  # keep headings with the text that follows
kronreport.H2.keepWithNext = 1


def build(out):
    r = Report("examples/filter_topology", "A low-pass filter grown from an empty embryo",
               "Can the topology genome invent a filter that meets a mask, starting from nothing between input "
               "and output? — CLAUDE.md §4, milestone M8", "CLAUDE.md §4, M8",
               "population 300, 200 generations, seed 1")

    r.h1("Summary")
    r.p("This is the M8 milestone example: topology evolution. The embryo holds only the test fixture (after Koza "
        "et al., <i>Genetic Programming III</i>, 1999, ch. 25): a 1 V source with 1 kΩ source resistance driving node "
        "<i>in</i>, and a 1 kΩ load on node <i>out</i>. Nothing connects <i>in</i> to <i>out</i>. Evolution decides "
        "which parts exist (inductors, capacitors, resistors with E24 values), how they are wired among ground, "
        "<i>in</i>, <i>out</i> and eight internal nodes, and their values.")
    r.bullets([
        "<b>Spec met from nothing:</b> seed 1 grows a four-part third-order elliptic (Cauer) low-pass after 200 "
        "generations (60,000 evaluations).",
        "<b>Checked on a denser sweep:</b> at 200 points per decade (ten times the density evolution used) the "
        "passband stays between −6.02 and −7.12 dB (mask −5.5 to −7.5 dB) and the stopband reaches at most "
        "−46.77 dB (mask −46 dB).",
        "<b>Every seed succeeds:</b> seeds 1–8 all meet the mask; six of them with four parts in one of the two "
        "textbook forms of the third-order elliptic low-pass.",
        "<b>Parsimony needs time:</b> stopped at “spec met”, seed 1 had 15 parts at generation 50; run to 200 "
        "generations, the 0.002-per-part penalty pared it to 4.",
        "<b>Penalties, never a crash:</b> 15–25 % of each generation fails validation or measurement and is "
        "penalized; none panics. A run checkpointed and resumed equals the uninterrupted run.",
    ])

    r.h1("The construction")
    r.h2("Embryo (embryo.cir)")
    r.code("VIN vs 0 DC 0 AC 1\nRS vs in 1k\nRL out 0 1k\n.ac dec 20 10 100k")
    r.p("With matched 1 kΩ source and load, an ideal lossless passband sits at −6.02 dB. The export sheet in "
        "<i>export/</i> shows the bare fixture: <i>out</i> at 0 V, and no automatic tearing possible because the "
        "circuit has no path to tear.")
    r.h2("Genome (evo/encode_topology.go)")
    r.p("A variable-length list of parts, each with a kind, two nodes and a log<sub>10</sub> value, grown into a copy "
        "of the embryo as L_g1, C_g2, …; the embryo is never modified. This is a direct netlist encoding, not Koza's "
        "developmental program. At decode, self-loops and parts with an end on an internal node nothing else touches "
        "are pruned, repeatedly: they carry no current, so every node voltage is unchanged, and without the pruning "
        "debris from mutation would turn into floating-node penalties. What remains invalid (capacitor-only islands, "
        "inductor loops) is penalized.")
    r.table([
        ["setting", "value"],
        ["terminals", "in, out (plus ground)"],
        ["internal nodes", "8"],
        ["parts", "3 at start, at most 20"],
        ["L", "10 µH – 10 H, E24"],
        ["C", "100 pF – 100 µF, E24"],
        ["R", "1 Ω – 1 MΩ, E24"],
        ["GA", "population 300, 200 generations, seed 1, 2 s timeout per evaluation"],
    ], [0.3, 0.7])
    r.p("<b>Operators.</b> Change value (Gaussian, in decades), add, delete, rewire one end, change kind (keeping the "
        "relative position in the range), split a part in series through a fresh node, add a part in parallel; one "
        "operator per child plus a geometric number more (mean 1.43). Crossover is one-point on the part lists. New "
        "parts attach to nodes already in use plus one fresh internal node, so they mostly join the circuit.")
    r.h2("Objectives")
    r.table([
        ["mask band", "requirement"],
        ["passband 10 Hz – 1 kHz", "−7.5 dB ≤ gain ≤ −5.5 dB (matched: −6.02 dB)"],
        ["stopband 3 kHz – 100 kHz", "gain ≤ −46 dB (40 dB below the passband)"],
    ], [0.35, 0.65])
    r.p("Fitness is the mask error (<i>evo.FilterMask</i>) plus 0.002 per part (<i>evo.PartCount</i>). The mask error "
        "is 0 when the response stays inside the mask at every sweep point and band edge; otherwise it is 0.1 + the "
        "fraction of points outside + the mean excess in dB/20, after Koza's filter scoring. The 0.1 floor means any "
        "design that leaves the mask scores worse than any design that meets it, so “fitness &lt; 0.1” reads "
        "“spec met”, and once the mask is met evolution keeps removing parts.")

    r.h1("Results")
    r.h2("Seed 1: the evolved filter")
    r.code("L_g1 in  out 0.2       series inductor\n"
           "C_g4 in  out 11n       across it: a notch at 1/(2π√(LC)) = 3.39 kHz\n"
           "C_g2 0   in  180n      shunt capacitors at both ends\n"
           "C_g3 out 0   220n")
    r.p("That is a third-order elliptic (Cauer) low-pass: a π-section with the trap across the series arm. Re-simulated "
        "from the written <i>best.cir</i> on a sweep of 200 points per decade, the passband stays between −6.02 and "
        "−7.12 dB and the stopband reaches at most −46.77 dB, inside the mask everywhere; this re-simulation is the "
        "acceptance test. <i>best.cir</i> opens in ngspice, LTspice and QSpice as it is.")
    r.h2("Seeds 1–8 (200 generations, population 300)")
    r.table([
        ["seed", "parts", "form"],
        ["1", "4", "series L with C across it (notch 3.4 kHz), shunt C at both ends"],
        ["2", "6", "same idea, with two small resistors"],
        ["3", "4", "series L with C across it, shunt C both ends"],
        ["4", "4", "series L, shunt C, series-LC trap to ground at the input (notch 3.9 kHz)"],
        ["5", "6", "two series L, trap across, shunt C both ends, one small resistor"],
        ["6", "4", "series L with C across it, shunt C both ends"],
        ["7", "4", "series L, shunt C, series-LC trap to ground"],
        ["8", "4", "series L, shunt C, series-LC trap to ground"],
    ], [0.1, 0.1, 0.8])
    r.p("All eight seeds meet the mask. The four-part designs are the two textbook forms of the third-order elliptic "
        "low-pass: the trap across the series inductor, or a series-resonant trap in a shunt arm. Third order is the "
        "lowest that can meet about 1 dB ripple and 40 dB attenuation from 3× the passband edge (a Butterworth needs "
        "five reactive parts, a Chebyshev four); here it is reached from an empty fixture with standard values.")
    r.h2("How the search behaves")
    r.bullets([
        "<b>Parsimony needs time.</b> With an early stop at “spec met”, seed 1 stopped at generation 50 with 15 parts; "
        "run to 200 generations, the 0.002-per-part penalty pared it to 4.",
        "<b>Penalties, never a crash.</b> Every decoded circuit passes <i>circuit.Validate</i> before simulation; "
        "15–25 % of each generation fails validation or measurement and takes a penalty from the M5 classes "
        "(no-measurement &lt; no-convergence &lt; singular &lt; invalid &lt; ledger &lt; timeout &lt; other &lt; "
        "panic). None panics.",
        "<b>Deterministic and resumable.</b> Operators draw from the per-generation PCG stream; a topology run "
        "checkpointed at generation 2 and resumed equals the uninterrupted run (<i>TestTopologyResume</i>).",
    ])

    r.h1("Validation")
    r.table([
        ["M8 criterion", "status (doc/M8.md)"],
        ["filter from an empty embryo meets spec",
         "<i>TestFilterTopologyExample</i>: written best netlist re-simulated at 200 points/decade stays inside the "
         "mask (passband −6.02…−7.12 dB, stopband ≤ −46.77 dB)"],
        ["every decoded circuit validated before simulation", "<i>circuit.Validate</i> before any analysis"],
        ["penalties, never a crash", "M5 penalty classes; none panics"],
        ["deterministic, resumable", "<i>TestTopologyResume</i>"],
        ["go test -race ./...", "clean"],
    ], [0.38, 0.62])

    r.h1("Limits")
    r.bullets([
        "Koza's developmental (circuit-constructing) encoding is not implemented; the direct encoding was enough for "
        "this milestone.",
        "Only R, C and L are evolvable parts; transistors and sources are not.",
        "No multi-part local search.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./cmd/evolve -config examples/filter_topology/config.json   # about 40 s\n"
           "go test ./evo -run 'TestFilterTopologyExample|TestTopologyResume'   # the M8 acceptance tests\n"
           "python examples/filter_topology/viz/make_report.py examples/filter_topology/filter_topology_report.pdf")
    r.small("Source: J. R. Koza, F. H. Bennett III, D. Andre, M. A. Keane, Genetic Programming III: Darwinian Invention "
            "and Problem Solving, Morgan Kaufmann (1999), ch. 25 — the embryo fixture, the filter mask scoring and the "
            "evolved low-pass filters this example repeats with a direct encoding.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "filter_topology_report.pdf"))
