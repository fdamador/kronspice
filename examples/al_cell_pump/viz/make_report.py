# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""H38 report: the aluminium-cell ground pump as a PDF.

python3 examples/al_cell_pump/viz/make_report.py examples/al_cell_pump/al_cell_pump_report.pdf

Numbers are those of run_output.txt (seed 1, 40 x 25 NSGA-II); figures are
viz/pump.png and viz/pipes.png. Layout from internal/pdfreport/kronreport.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
import kronreport as K  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph,  # noqa: E402
                                Spacer, Table, TableStyle)

BODY, SMALL, CAP, H1, H2 = K.BODY, K.SMALL, K.CAP, K.H1, K.H2
CELL, CELLB, BULLET = K.CELL, K.CELLB, K.BULLET


def P(t, s=BODY):
    return Paragraph(t, s)


def bullets(items):
    return [Paragraph(t, BULLET, bulletText="•") for t in items]


def table(rows, widths, bold_rows=()):
    data = [[Paragraph(str(c), CELLB if (r == 0 or r in bold_rows) else CELL) for c in row]
            for r, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    st = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#dfe7f1")),
        ("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.HexColor("#1f3a5f")),
        ("LINEBELOW", (0, -1), (-1, -1), 0.4, colors.HexColor("#999999")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]
    for r in range(1, len(rows)):
        if r % 2 == 0:
            st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#f5f7fa")))
    for r in bold_rows:
        st.append(("BACKGROUND", (0, r), (-1, r), colors.HexColor("#fff4d6")))
    t.setStyle(TableStyle(st))
    return t


def figure(name, width, caption):
    path = os.path.join(HERE, name)
    img = Image(path)
    img.drawHeight = width * img.imageHeight / img.imageWidth
    img.drawWidth = width
    return KeepTogether([img, P(caption, CAP)])


def build(out):
    W = A4[0] - 36 * mm
    s = []

    r = K.Report("examples/al_cell_pump", "Aluminium-cell ground pump",
                 "Evolved for catalytic efficiency — register entry H38", "register entry H38",
                 "seed 1, NSGA-II 40 × 25")
    W = r.W
    s = r.s

    s.append(P("Summary", H1))
    s.append(P(
        "The negative Cockcroft–Walton ground pump of <i>photon_electron</i> was rebuilt with its rectifiers "
        "replaced by aluminium–water valve cells (the Pollak/Nodon electrolytic rectifier of <i>sparkgap_diode</i>). "
        "A genetic search (NSGA-II) then designed the cells and the drive so that electrons are pumped from ground "
        "onto a negative bank with as little energy as possible. Two layouts were evolved: plate stacks, with "
        "each cell in its own tank, and concentric pipes, with one tank per rectifier."))
    s += bullets([
        "<b>Best result:</b> η = 0.857 at −321 V and 7.8 mA from 2.9 W, using one cell of six 13 mm² plates formed at "
        "545 V and driven at 471 Hz. It delivers 1.66×10<super>16</super> electrons per joule, against the "
        "e·|V| ideal of 1.94×10<super>16</super>. A single aluminium rod in an iron cup does slightly better: "
        "η = 0.869 at −197 V.",
        "<b>The baseline is lost in the tank, not the rectifiers.</b> At 10 V and 100 kHz, 98–99 % of the drive "
        "heats the tank coil, so η is 0.0006 with silicon diodes and 0.00002 with aluminium cells. Evolution moves to "
        "small tanks at hundreds of hertz.",
        "<b>The oxide film's capacitance is the stack's main barrier.</b> Every efficient design uses the least "
        "aluminium area and the highest forming voltage it can, with cells in series.",
        "<b>η falls with voltage:</b> 0.5 at −3.4 kV and 0.15 at −10 kV, because the coil's share of the loss grows again.",
        "<b>Silicon diodes are still better,</b> by 0.01–0.09 in η on every evolved design. The 2 % varactor effect of "
        "the oxide lowers η by up to 0.013, and it is far too small for parametric (2f) pumping.",
        "<b>No design goes under the e·|V| floor.</b> The ledger enforces this as a test, so the result is a measure "
        "of how close a real circuit gets to the floor, not a way around it.",
    ])

    s.append(P("What “catalytic efficiency” means", H1))
    s.append(P(
        "Lifting one electron from ground to a bank at −|V| costs at least e·|V|, whatever the circuit. That is a "
        "floor, and no resonance or barrier lowering goes under it. What a circuit <i>can</i> change is how much "
        "more than e·|V| it spends. The fitness is therefore"))
    s += bullets([
        "<b>η = P<sub>load</sub> / P<sub>in</sub></b> in steady state. The bank voltage must drift less than 1 % over "
        "the measuring window, so energy still going into capacitors is not counted as output.",
        "Two further objectives: the bank voltage and the load current.",
        "Two constraints: a stack never sees more than 90 % of its holding voltage in reverse, and the bank must be "
        "in steady state.",
    ])
    s.append(P(
        "η = 1 would mean the drive pays only e·|V| per electron. Each design's input is split into the load, coil "
        "resistance, conducting rectifiers (electrolysis onset plus electrolyte), leaking rectifiers and the rest. "
        "The parts sum to the input within 1 %."))

    s.append(P("The aluminium valve stack", H1))
    s.append(P("A stack is <i>Ns</i> cells in series × <i>Np</i> plates in parallel, with plate area <i>A</i>, a water "
               "gap and a forming voltage <i>V</i><sub>form</sub>:"))
    s.append(table([
        ["property", "value", "from"],
        ["holds (reverse)", "Ns · Vform", "oxide formed at Vform, ≤ 700 V"],
        ["forward onset", "Ns · 2 V", "electrolysis per cell"],
        ["forward resistance", "Ns · ρ · gap / (A · Np), ρ = 0.1 Ω·m", "water electrolyte"],
        ["reverse leakage", "A·Np · j0 · exp(β·E), 0.1 A/m² at Vform", "Güntherschulze–Betz, β = 3×10<super>−8</super> m/V, 1.4 nm/V"],
        ["above Vform", "sparks (2 Ω per cell)", "breakdown"],
        ["film capacitance", "εr 9.5 · ε0 · A·Np / (1.4 nm/V · Vform) / Ns", "formed oxide"],
        ["varactor", "C(v) = C0 · (1 + κ·|v|/hold), κ = 2 %", "fixed-thickness barrier oxide"],
    ], [0.2 * W, 0.42 * W, 0.38 * W]))
    s.append(Spacer(1, 3 * mm))
    s.append(P(
        "Cells in series raise the voltage the stack holds and divide its film capacitance. Plates in parallel raise "
        "the current and multiply the film capacitance. The 11 genes are drive frequency (50 Hz–100 kHz) and "
        "amplitude (1–50 V), tank capacitor, stages (1–8), ladder capacitors, load, cells in series (1–10), plates in "
        "parallel (1–20), plate area (10 mm²–100 cm²), gap (1–20 mm) and forming voltage (50–700 V). "
        "The film is evolved as a linear capacitor (κ = 0), which runs 16× faster. Every Pareto design is then "
        "re-measured with κ = 2 % and again with silicon diodes."))

    s.append(PageBreak())
    s.append(P("Results: plate stacks", H1))
    s.append(figure("pump.png", W,
                    "Figure 1. (a) The Pareto front, with the baselines at the bottom. (b) Where the input goes. "
                    "(c) The most efficient design charging its bank from rest. (d) A multi-plate stack."))
    s.append(table([
        ["design", "bank", "load", "input", "η", "coil", "forward"],
        ["baseline, silicon diodes (10 V, 100 kHz, 1 µF tank)", "−2953 V", "0.30 mA", "1.55 kW", "0.0006", "98.9 %", "0 %"],
        ["baseline, aluminium stacks", "−194 V", "0.02 mA", "210 W", "0.00002", "98.0 %", "0 %"],
        ["evolved, most efficient (1 cell × 6 plates)", "−321 V", "7.8 mA", "2.9 W", "0.857", "12.8 %", "1.4 %"],
        ["evolved, 5 cells × 3 plates", "−357 V", "36 mA", "15.6 W", "0.816", "", ""],
        ["evolved, 8 cells × 5 plates, 2 stages", "−1113 V", "12 mA", "17.9 W", "0.751", "", ""],
        ["evolved, 7 cells × 10 plates, 3 stages", "−3372 V", "69 mA", "459 W", "0.505", "", ""],
        ["evolved, 7 cells × 8 plates, 3 stages", "−10.2 kV", "140 mA", "9.5 kW", "0.150", "", ""],
    ], [0.38 * W, 0.1 * W, 0.1 * W, 0.1 * W, 0.1 * W, 0.11 * W, 0.11 * W], bold_rows=(3,)))
    s.append(Spacer(1, 3 * mm))
    s.append(P(
        "<b>Most efficient design.</b> One cell of six 13.2 mm² plates with a 7.9 mm gap, formed at 545 V, in one "
        "stage. It is driven at 471 Hz and 50 V into a 0.25 µF tank. The stack holds 545 V; its forward path is "
        "2.0 V plus 9.96 Ω, and its film is 0.00875 µF. The tank swings 179 V, and the largest reverse voltage is 60 % "
        "of the oxide's rating. The electrolysis onset alone would allow η ≤ |V|/(|V| + 4 V) = 0.988; most of what "
        "is lost goes to the coil."))

    s.append(P("What the evolution found", H2))
    s += bullets([
        "<b>The baseline's loss is the tank.</b> A 1 µF tank at 100 kHz circulates hundreds of amps to hold its swing. "
        "Evolution moves to 0.2–1 µF tanks at hundreds of hertz to a few kilohertz. There the coil takes 13–20 % "
        "on the designs above η = 0.75.",
        "<b>The film capacitance is the main barrier.</b> At 100 kHz the baseline stack's 0.05 µF film all but "
        "shorts it, the high-frequency failure the 1911 <i>Britannica</i> describes. Efficient designs push "
        "plate area to 13–60 mm² and forming voltage to 545–700 V, and stack cells in series.",
        "<b>Conduction costs little.</b> The 2 V per-cell onset and the water's resistance together take 0.4–9 % of "
        "the input.",
        "<b>Leakage is negligible</b> once reverse voltage is kept under 90 % of the forming voltage.",
        "<b>Voltage and current cost efficiency.</b> Above about −1 kV, designs need more stages and bigger tanks. The "
        "coil's share then grows again to 50–84 %.",
        "<b>Silicon is still better</b> by 0.01–0.09 in η. The aluminium stack works, but its film and its onset "
        "voltage are costs, not a catalyst.",
    ])

    s.append(P("The parametric route (2f pumping)", H2))
    s.append(P(
        "Driving the tank at 2f through the film's C(v) grows an oscillation at f only if the capacitance swing "
        "m = ΔC/C exceeds 2/Q (the Mathieu threshold). The formed oxide gives m ≤ 2 % even when swung to its "
        "forming voltage, so it would need Q > 100. The coil's Q of 50 needs m > 4 %. The oxide is too linear to "
        "pump parametrically; a varactor diode, whose capacitance changes several-fold, could."))

    s.append(PageBreak())
    s.append(P("One tank per rectifier: concentric pipes", H1))
    s.append(P(
        "In this layout an aluminium rod sits at the centre, pipes surround it and the tank wall is outside. Each "
        "annulus of water between two neighbouring pipes is one cell. Two wirings were evolved:"))
    s += bullets([
        "<b>Parallel (bus bars).</b> The pipes alternate aluminium and iron, and each metal goes to its own bus bar, "
        "so every annulus is the same valve between the same two buses. Neighbouring annuli sit at the same voltage, "
        "so the pipe ends need no seals. The tank holds one forming voltage.",
        "<b>Series (clad pipes).</b> Each pipe is aluminium outside and iron inside, and the current crosses every annulus "
        "in turn, so the tank holds N forming voltages. Neighbouring annuli differ by one cell voltage, so every pipe "
        "bottom must be sealed. In <b>graded series</b>, the outer pipes' aluminium is masked so that every annulus "
        "has the same area.",
    ])
    s.append(figure("pipes.png", W,
                    "Figure 2. (a) Pipe fronts over the plate-stack front. (b) Seal leakage on the best series tank "
                    "that needs its annuli. (c) Parallel annuli for current. (d) The two tank layouts."))
    s.append(table([
        ["tank", "bank", "load", "input", "η"],
        ["single cell: Al rod Ø3.9 mm in an iron tank Ø7 mm, 5 mm wetted, formed at 414 V", "−197 V", "4.0 mA", "0.92 W", "0.869"],
        ["graded series, 2 annuli, tank Ø18 mm", "−282 V", "6.7 mA", "2.2 W", "0.856"],
        ["graded series, 2 annuli, tank Ø10 mm", "−959 V", "31 mA", "37 W", "0.801"],
        ["graded series, 8 annuli, tank Ø42 mm, 3 stages", "−2.5 kV", "13 mA", "53 W", "0.606"],
        ["graded series, 5 annuli, tank Ø57 mm, 4 stages", "−4.8 kV", "19 mA", "228 W", "0.405"],
    ], [0.56 * W, 0.11 * W, 0.11 * W, 0.11 * W, 0.11 * W], bold_rows=(1,)))
    s.append(Spacer(1, 3 * mm))
    s += bullets([
        "<b>Pipes match plates.</b> The pipe front lies on the plate front: η 0.87 against 0.86 at a few hundred volts, "
        "and 0.6 against 0.6 at 2.5 kV. One tank of up to 10 annuli replaces up to 10 separate cells.",
        "<b>The less aluminium, the better.</b> The wetted length stays at its 5 mm minimum on almost every design. The "
        "best rectifier is a thimble: one aluminium rod in an iron cup.",
        "<b>Parallel annuli do not help at these currents.</b> Rewiring the best single cell as 2, 4 or 6 parallel "
        "annuli lowers η: 0.869, 0.854, 0.803 and 0.715. The current is set by the ladder and the drive, and the "
        "forward loss is the 2 V onset, not the water. Parallel annuli pay off only near 1 A.",
        "<b>Plain series tanks load the inner annulus.</b> The AC voltage divides as 1/C, so in a test design the annuli "
        "take 92, 16 and 9 % of their forming voltage. In the graded tank they share it equally (41/41/41 %), and "
        "graded series tanks fill 33 of the 37 places on the front.",
        "<b>Seals matter only when the tank needs its annuli.</b> The −959 V tank (1.38 forming voltages) keeps 99 % of "
        "its η only with at least about 0.3 MΩ between annuli. Under a 5 mm-wide seal, that means a water film "
        "thinner than about 70 nm: a bonded, potted seal, not a gasket. Pipes merely standing on the floor (about "
        "13 Ω) spark the inner annulus, and η falls to 0.57.",
    ])

    s.append(P("Limits of the model", H1))
    s += bullets([
        "The cell constants are representative literature values (water resistivity, Güntherschulze–Betz β, oxide "
        "growth, leakage at forming), not measurements on a particular cell. A bench cell should be characterised for "
        "forward onset, leakage against voltage and film capacitance, and these constants replaced.",
        "The film does not age or re-form during a run. A real oxide loses its holding voltage when idle and "
        "re-forms with current.",
        "The coil has a fixed Q of 50. A better coil, or a lower-loss tank, would move the whole front up.",
        "Electrolysis gas, heating and changes in the electrolyte are not modelled. Foam over the pipe tops could "
        "bridge annuli like a leaky seal.",
        "Forward current density is not limited. The high-current designs reach about 1 kA/m² on small aluminium "
        "faces, where gassing would be heavy.",
    ])

    s.append(P("Reproducing the results", H1))
    s.append(P(
        "<font name='Courier' size='8'>go run ./examples/al_cell_pump</font> runs both studies, two 40 × 25 evolutions "
        "taking about 15 min on 2 CPUs. <font name='Courier' size='8'>-study pipes</font> runs only the pipe study, and "
        "<font name='Courier' size='8'>-data DIR</font> writes the CSVs for the figures. "
        "<font name='Courier' size='8'>go test ./examples/al_cell_pump</font> takes about 6 s and checks, among other "
        "things, that electrons per joule never exceed the e·|V| ideal. The full output is in "
        "<font name='Courier' size='8'>run_output.txt</font>."))
    s.append(P("Status in the register: <b>bounded</b>. The standard prediction was computed, and nothing goes under e·|V|.", SMALL))

    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "al_cell_pump_report.pdf"))
