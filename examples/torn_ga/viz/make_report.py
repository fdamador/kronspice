# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""PDF report for examples/torn_ga: Kron diakoptics inside a running genetic algorithm (milestone M9b).

Run from the repository root:

    python examples/torn_ga/viz/make_report.py examples/torn_ga/torn_ga_report.pdf

Numbers come from examples/torn_ga/README.md, doc/M9b.md (the two-CPU timings and
cache counts) and examples/torn_ga/export/torn_ga_line40_results.md (the torn
solve against MNA). The example has no figures and no run_output.txt.
Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401


def build(out):
    r = Report("examples/torn_ga", "Torn GA: diakoptics in a running evolution",
               "Does tearing the circuit (Kron's diakoptics) with cached pieces make a value-genome GA faster, "
               "with the same result? — CLAUDE.md §3a, milestone M9b",
               "CLAUDE.md §3a, milestone M9b", "population 40, 30 generations, seed 1")

    r.h1("Summary")
    r.p("A two-pole RC pre-filter (R1, C1, R2, C2, from the E24 and E12 series) is evolved in front of a fixed plant, "
        "joined to it by a resistor RT, for a 5 kHz cutoff at the plant's output. The same run (same seed) is made "
        "three times: AC analysis by modified nodal analysis (MNA), by Kron's node-pair form untorn, and by Kron's "
        "diakoptics torn at RT with the plant's factored pieces cached and shared by all workers. A mutation changes "
        "only the pre-filter, so only its piece is re-factored and the plant's factors are reused from one candidate "
        "to the next. This is the acceptance example of milestone M9b.")
    r.bullets([
        "<b>Same best design from all three solvers</b> on both plants; the program exits non-zero if they differ. "
        "On a 16-CPU run of the line: R1 13 Ω, C1 47 nF, R2 3 kΩ, C2 390 nF, fitness 5.6×10<super>−6</super>.",
        "<b>Torn is an order of magnitude faster than untorn Kron:</b> 11.1× on the 40-section LC line and 5.6× on the "
        "12×12 RC grid (two CPUs); 8.2× on the line with 16 CPUs.",
        "<b>Torn Kron is still slower than sparse MNA:</b> 2.0× on the line and 1.5× on the grid (two CPUs), 2.3× on the "
        "line with 16 CPUs. MNA stays the default; tearing is opt-in.",
        "<b>The cache works as intended:</b> on the line every candidate factors its own 4-unknown piece and reuses the "
        "121-unknown line (233,046 pieces taken from the cache, 226,038 factored, 16-CPU run).",
        "<b>Torn equals MNA:</b> node voltages agree to 4.5×10<super>−15</super> of the largest node voltage at 10 kHz, "
        "and the power taken by the R, L, C branches is identical in both (Roth's invariance).",
    ])

    r.h1("What is being tested")
    r.p("CLAUDE.md §3a gives the reason <i>kron</i> exists in this project: a mutation usually changes one part of a "
        "circuit, so with tearing only the torn piece containing that change is re-solved and the other pieces' factored "
        "solutions are reused from the parent. M9b is done when the value-genome GA runs faster with tearing than "
        "without, with the same results.")

    r.h1("The model")
    r.h2("Circuit")
    r.p("V1 in 0 AC 1, the pre-filter R1 in m 1k, C1 m 0 10n, R2 m f 1k, C2 f 0 10n, the coupling RT f p0 10, and one "
        "of two plants from p0 to out:")
    r.bullets([
        "<b>line:</b> n sections of 1 Ω, 1 µH, 1 nF (√(L/C) ≈ 31.6 Ω), matched by RL out 0 31.6; 40 sections by default, "
        "125 unknowns;",
        "<b>grid:</b> an n×n RC mesh, 10 Ω between neighbours and 1 nF to ground at every node, driven at one corner and "
        "loaded by 100 Ω at the other; a plant whose LU factorization fills in; 12×12 by default, 144 nodes.",
    ])
    r.h2("Problem")
    r.table([
        ["Item", "Setting"],
        ["Genes", "R1, R2: 10 Ω – 100 kΩ (E24); C1, C2: 100 pF – 10 µF (E12)"],
        ["Objective", "<i>evo.Cutoff</i> at node out, target 5 kHz"],
        ["Analysis", ".ac dec 50 100 1meg (201 frequencies per candidate)"],
        ["Run", "population 40, 30 generations (1,240 evaluations), seed 1"],
    ], [0.2, 0.8])
    r.h2("Three solvers")
    r.bullets([
        "<b>MNA:</b> the ordinary <i>analysis.AC</i> path (sparse Gilbert–Peierls LU with a reused pattern, M6).",
        "<b>Kron, untorn:</b> <i>evo.Problem.Tear</i> with no torn branches, so every candidate goes through the node-pair "
        "form Y′ = A<super>T</super>·Y·A (M9a).",
        "<b>Kron, torn at RT:</b> <i>kron.TornAC</i> cuts RT out, so the node-pair matrix is block-diagonal (pre-filter "
        "in one block, plant in the other) plus RT's rank-one term.",
    ])
    r.p("The torn solve, with B the block-diagonal matrix, U the torn branch's incidence, D its admittance and c the "
        "source vector:")
    r.code("z = B^-1·c                       each piece on its own\n"
           "(I + D·U^T·B^-1·U)·w = D·U^T·z     interconnection system, order 1 here\n"
           "x = z − B^-1·U·w")
    r.p("Here w is the current through RT, Kron's intersection current through the open path that tearing created. "
        "The torn solve refines against the full matrix and checks its componentwise backward error, reporting a nearly "
        "floating piece instead of returning a wrong answer. A <i>kron.TornCache</i> keeps factored pieces by their exact "
        "content (every entry's bit pattern, confirmed entry by entry on a hit), bounded by count and by 128 MiB; a "
        "cached and an uncached solve are bit-identical.")

    r.h1("Results")
    r.h2("Milestone run, two CPUs")
    r.table([
        ["Plant", "MNA", "Kron node form, untorn", "Kron diakoptics, torn at RT", "Same best design"],
        ["40-section LC line (125 unknowns)", "6.4 s", "138.8 s", "12.5 s (11.1× faster than untorn)", "yes"],
        ["12×12 RC grid (144 nodes)", "21.3 s", "177.0 s", "31.5 s (5.6× faster than untorn)", "yes"],
    ], [0.28, 0.1, 0.18, 0.3, 0.14])
    r.p("On the line the plant's pieces were taken from the cache 235,394 times and the pre-filter's piece factored "
        "223,690 times (doc/M9b.md).")
    r.h2("Line run on a 16-CPU machine (2026-11-01)")
    r.table([
        ["Solver", "Time", "Pieces cached / factored"],
        ["MNA", "1.18 s", ""],
        ["Kron node form, untorn", "21.80 s", "9 / 229,533"],
        ["Kron diakoptics, torn at RT", "2.67 s", "233,046 / 226,038"],
    ], [0.4, 0.2, 0.4], bold_rows=(3,))
    r.p("Best design in all three: R1 13 Ω, C1 47 nF, R2 3 kΩ, C2 390 nF, fitness 5.6×10<super>−6</super>. Torn is "
        "8.2× faster than untorn Kron here.")
    r.h2("Against MNA")
    r.p("The torn Kron solve is still slower than the MNA engine: 2.0× on the line and 1.5× on the grid in the two-CPU "
        "run, and 2.3× on the line in the 16-CPU run. Once the fixed piece is cached, what remains per frequency is "
        "building the primitive network, the branch relations and the node-pair rows, which costs about as much as MNA's "
        "whole sparse factorization of a chain. On the grid, where the LU fills in, the gap narrows. Reusing the "
        "per-network structure across the frequencies of one sweep, and building only the changed piece of a mutated "
        "candidate, are the steps listed in doc/M9b.md that would let the torn path overtake sparse MNA.")
    r.h2("Torn solution against MNA (export sheet)")
    r.p("The 40-section line problem, analysed at 10 kHz, is torn at RT into 2 pieces of 4 and 121 unknowns joined by "
        "1 link; the interconnection system has order 1, against 125 unknowns for the whole circuit.")
    r.table([
        ["Quantity", "Torn (Kron)", "Untorn (MNA)"],
        ["Link current through RT", "0.000454594 A, phase −20.73°", ""],
        ["Voltage across RT", "0.00454594 V, phase −20.73°", ""],
        ["Power V·I* in the R, L, C branches", "0.000553677, phase −16.21°", "0.000553677, phase −16.21°"],
        ["Node out", "0.0143161, phase −28.08°", "0.0143161, phase −28.08°"],
        ["Node f", "0.0369572, phase −23.80°", "0.0369572, phase −23.80°"],
        ["Largest node-voltage difference", "4.5×10<super>−15</super> of the largest node voltage", ""],
    ], [0.38, 0.32, 0.3])

    r.h1("Validation")
    r.bullets([
        "<i>TestSameBestFasterTorn</i> runs a 20-section line (population 16, 6 generations, seed 3) and requires the same "
        "best values from all three solvers and a torn time below the untorn one.",
        "<i>BenchmarkGeneration</i> times one generation of 8 candidates with a fresh seed each iteration (the plant's "
        "pieces stay cached, the pre-filter's do not) for line/40 and grid/12, by MNA and torn.",
        "Milestone checks in doc/M9b.md: the torn solve equals the untorn solve and MNA on the linear golden circuits and "
        "60 random circuits, 2,289 torn solves, worst 5.1×10<super>−11</super> relative to MNA; power invariant under "
        "tearing; a cut that leaves a singular block is reported, not solved.",
    ])

    r.h1("Limits")
    r.bullets([
        "Cuts are named by hand; there is no automatic partitioner in this example.",
        "Only AC analysis of linear candidates is torn; nonlinear circuits and circuits with idealized active elements "
        "are solved by MNA.",
        "Timings depend on the machine; the README reports a two-CPU and a 16-CPU run.",
    ])

    r.h1("Reproducing the results")
    r.code("go run ./examples/torn_ga                      # 40-section LC line\n"
           "go run ./examples/torn_ga -plant grid -n 12    # 12×12 RC mesh\n"
           "go run ./examples/torn_ga -gens 30 -pop 40 -seed 1\n"
           "go test ./examples/torn_ga                     # same best design, torn faster than untorn\n"
           "go test ./examples/torn_ga -bench Generation   # one generation, fresh seeds")
    r.small("Sources: G. Kron, Diakoptics (1963); G. Kron, Tensor Analysis of Networks (1939), p. 133; H. H. Happ, "
            "Diakoptics and Networks (1971), Ch. XII; P. W. Aitchison, Quart. Appl. Math. 41, 265–272 (1983); J. P. Roth, "
            "Quart. Appl. Math. 17 (1959), Sec. 5; J. R. Gilbert, T. Peierls, SIAM J. Sci. Stat. Comput. 9, 862 (1988).")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "torn_ga_report.pdf"))
