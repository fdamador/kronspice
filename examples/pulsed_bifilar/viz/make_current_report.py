# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""pulsed_bifilar current report: volts, amps and power of the bifilar solenoid with its
strands in series and in parallel, currents the same way and opposing (CLAUDE.md §3e).

python examples/pulsed_bifilar/viz/make_current_report.py examples/pulsed_bifilar/pulsed_bifilar_current_report.pdf

Every number is read from run_output.txt (coils, pulse sweeps, ledgers) and
viz/vi_sweep.csv (the AC V–I sweep) and ngspice_twin_output.txt (kronspice
against ngspice, ends open and joined: viz/ngspice_twin.py); figures are viz/*.png made by
viz/render.py. Layout from internal/pdfreport/kronreport.py.
"""
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "internal", "pdfreport"))
from kronreport import Report, esc  # noqa: E402,F401

RUN = os.path.join(HERE, "..", "run_output.txt")
VI = os.path.join(HERE, "vi_sweep.csv")
TWIN = os.path.join(HERE, "..", "ngspice_twin_output.txt")

# tag, label, heading printed by the Go program
COILS = [
    ("ordinary", "Ordinary solenoid", "ordinary cylindrical solenoid"),
    ("bifilar", "Series bifilar, same direction (Tesla)", "Tesla bifilar cylindrical solenoid"),
    ("counter", "Series bifilar, opposite currents", "opposite-current bifilar cylindrical solenoid"),
    ("parallel", "Parallel bifilar, same direction", "parallel bifilar (same direction)"),
    ("parallel_counter", "Parallel bifilar, opposite currents", "parallel bifilar (opposite currents)"),
]


def eng(x, unit, d=3):
    """x with an SI prefix."""
    if x == 0:
        return f"0 {unit}"
    for p, s in ((1e9, "G"), (1e6, "M"), (1e3, "k"), (1, ""), (1e-3, "m"), (1e-6, "µ"), (1e-9, "n"), (1e-12, "p")):
        if abs(x) >= p:
            return f"{x / p:.{d}g} {s}{unit}"
    return f"{x:.{d}g} {unit}"


def parse_run():
    txt = open(RUN, encoding="utf-8").read()
    out = {}
    for block in re.split(r"\n== ", txt)[1:]:
        head = block.split(" winding ==")[0]
        tag = next((c[0] for c in COILS if head.startswith(c[2])), None)
        if tag is None:
            continue
        d = {"weak": "no resonance with gain above 2" in block}
        d["Ctt"] = float(re.search(r"neighbouring turns ([0-9.e+-]+) nF", block).group(1)) * 1e-9
        m = re.search(r"f1 = ([0-9.e+]+) Hz \(wire = ([0-9.]+) wavelengths.*?×([0-9.]+), Q ≈ ([0-9]+)", block)
        d["f1"], d["wl"], d["gain"], d["Q"] = float(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))
        d["M"] = float(re.search(r"mutual inductance to the primary ([0-9.e+-]+) µH", block).group(1)) * 1e-6
        ml = re.search(r"; ([0-9.e+-]+) µH round the loop", block)
        d["Mloop"] = float(ml.group(1)) * 1e-6 if ml else 0.0
        li = re.search(r"inductance ([0-9.e+-]+) (mH|µH)", block)
        d["L"] = float(li.group(1)) * (1e-3 if li.group(2) == "mH" else 1e-6)
        num = r"\s+([0-9.e+-]+)"
        rows = re.findall(r"^" + num * 11 + r"\s*$", block, re.M)
        d["rows"] = [[float(x) for x in r] for r in rows]
        lg = re.search(r"ledger at f1: in ([0-9.e+-]+) J = heat ([0-9.e+-]+) J \+ stored ([0-9.e+-]+) J "
                       r"\(\+ time-stepping error ([0-9.e+-]+) J\)", block)
        d["ledger"] = [float(lg.group(i)) for i in range(1, 5)]
        out[tag] = d
    return out


def at_f1(d):
    return next(r for r in d["rows"] if abs(r[0] - 1) < 1e-9)


def vi_peaks():
    rows = list(csv.reader(open(VI, encoding="utf-8")))
    head, data = rows[0], [[float(x) for x in r] for r in rows[1:]]
    col = {h: [r[i] for r in data] for i, h in enumerate(head)}
    out = {}
    for tag, _, _ in COILS:
        if tag + "_v_ends" not in col:
            continue
        v, i = col[tag + "_v_ends"], col[tag + "_i_wire"]
        k = max(range(len(v)), key=lambda j: v[j])
        ki = max((j for j in range(len(i)) if col["freq"][j] < 1e6), key=lambda j: i[j])  # below 1 MHz, well inside the model
        out[tag] = {"f": col["freq"][k], "v": v[k], "i": i[k], "fi": col["freq"][ki], "imax": i[ki],
                    "v1k": v[0], "i1k": i[0]}
    return out


def parse_twin():
    """ngspice_twin_output.txt (viz/ngspice_twin.py): {(tag, run): {quantity: (kronspice, ngspice, rel)}}."""
    out = {}
    for line in open(TWIN, encoding="utf-8").read().split("\n")[1:]:
        p = line.split()
        if len(p) != 6:
            continue
        out.setdefault((p[0], p[1]), {})[p[2]] = (float(p[3]), float(p[4]), float(p[5]))
    return out


def build(out):
    run, vi = parse_run(), vi_peaks()
    o, b, c, p, pc = (run[t] for t in ("ordinary", "bifilar", "counter", "parallel", "parallel_counter"))
    r = Report("examples/pulsed_bifilar", "Bifilar coil currents: series, parallel, same and opposing",
               "Volts, amps and power of one 400-turn bifilar solenoid wired four ways — CLAUDE.md §3e",
               "CLAUDE.md §3e", "100 V, 50 Ω one-sign pulses through a 5-turn primary; 20 sections per coil")

    r.h1("Summary")
    r.p("One 400-turn bifilar solenoid (0.5 mm copper, 126 m of wire, 10 cm across) is wired four ways: its two "
        "strands in <b>series</b> or in <b>parallel</b>, with their currents circulating the <b>same way</b> or "
        "<b>opposing</b>. The ordinary 400-turn solenoid of the same wire is the reference. Each coil is driven by "
        "one-sign 100 V pulses at its own first resonance f<sub>1</sub> through the same 5-turn primary, and swept by AC "
        "from 1 kHz to the model's λ/10 limit.")
    ratio_c = (o["f1"] / b["f1"]) ** 2
    r.bullets([
        f"<b>Series, same direction (Tesla, US 512,340):</b> L = {eng(b['L'], 'H')}, the same as the ordinary coil; "
        f"f<sub>1</sub> falls {o['f1'] / b['f1']:.0f}× to {eng(b['f1'], 'Hz')} because the capacitance between strands rises "
        f"about {ratio_c:,.0f}×. Under AC at resonance the wire current is about the same as the ordinary coil's "
        f"({eng(vi['bifilar']['i'], 'A')} against {eng(vi['ordinary']['i'], 'A')} per volt of drive) while the coil voltage "
        f"is {vi['ordinary']['v'] / vi['bifilar']['v']:.0f}× lower: same L, so for the same current the voltage scales as "
        "√(L/C).",
        f"<b>Series, opposite currents</b> (the wire folded and wound doubled): the fields cancel, L = {eng(c['L'], 'H')} "
        f"({b['L'] / c['L']:,.0f}× less), the primary's coupling falls {b['M'] / c['M']:.0f}×, and there is no resonance "
        f"with gain above 2; its strongest peak ({eng(c['f1'], 'Hz')}) gives ×{c['gain']:.3g} at the terminals.",
        f"<b>Parallel, same direction:</b> L = {eng(p['L'], 'H')}, exactly a quarter of the series winding's; the strands "
        f"sit at the same potential, so the strand capacitance drops out and f<sub>1</sub> = {eng(p['f1'], 'Hz')}, "
        f"{p['f1'] / o['f1']:.2f}× the ordinary coil's and {p['f1'] / b['f1']:.0f}× the series bifilar's.",
        f"<b>Parallel, opposite currents:</b> L = {eng(pc['L'], 'H')}; no resonance at the terminals (×{pc['gain']:.3g}). "
        "Inside it, current can run up one strand and down the other, a loop whose turns all circulate the same way: "
        f"a shorted coil coupled to the primary by {eng(pc['Mloop'], 'H')} (against {eng(pc['M'], 'H')} at the terminals), "
        f"which the terminals do not see. Pulsed at its peak, its neighbouring turns reach {eng(at_f1(pc)[3], 'V')} apart "
        f"while the terminals show {eng(at_f1(pc)[2], 'V')}: the strands meet at the same point of the wire, but their "
        "induced EMFs run opposite ways, so the loop's voltage sits across the enamel.",
        "<b>Opposing currents move the current inside the winding.</b> Both opposite-current windings show a comb of "
        f"wire-current peaks under AC below 1 MHz (largest {eng(vi['counter']['imax'], 'A')} and {eng(vi['parallel_counter']['imax'], 'A')} "
        f"per volt of drive, against {eng(vi['bifilar']['imax'], 'A')} and {eng(vi['parallel']['imax'], 'A')} for the same-direction windings) with almost no voltage at the terminals: modes of the "
        "strand pair, whose fields cancel outside.",
        "<b>Every joule comes from the pulse source</b> in every run; the ledgers close to the integration error.",
        "<b>Joining the two ends lowers voltage and current</b> (the coil becomes a shorted secondary with no "
        "resonance), in kronspice and in ngspice alike: AC agrees to 10<super>−5</super>, pulses to a few per cent "
        "(section \"The two ends joined, checked against ngspice\").",
    ])

    r.h1("The four connections")
    r.p("All four use the same turn positions: strand 1 on the even positions from the base up, strand 2 on the odd "
        "positions between them. Only the connection and the sense of circulation differ.")
    r.table([
        ["Winding", "Connection", "Strand 2 circulates", "Neighbouring turns sit", "Inductance"],
        ["Series, same direction (Tesla)", "end of 1 to start of 2", "same way", "half the wire apart (V/2)", eng(b["L"], "H")],
        ["Series, opposite currents", "wire folded at the top", "the other way", "from the whole wire apart (base) to 0 (fold)", eng(c["L"], "H")],
        ["Parallel, same direction", "both strands base → top", "same way", "at the same point of their strands (≈ 0 V)", eng(p["L"], "H")],
        ["Parallel, opposite currents", "both strands base → top", "the other way", "at the same point, but with opposite induced EMFs", eng(pc["L"], "H")],
        ["Ordinary (reference)", "one strand", "—", "one turn apart (V/N)", eng(o["L"], "H")],
    ], [0.25, 0.18, 0.14, 0.28, 0.15])
    r.p("The parallel inductance follows from the strands' block sums L<sub>1</sub>, L<sub>2</sub> and their mutual M: "
        "(L<sub>1</sub>L<sub>2</sub> − M²)/(L<sub>1</sub> + L<sub>2</sub> − 2M). For identical strands this is (L<sub>1</sub> + M)/2, "
        "a quarter of the series value 2(L<sub>1</sub> + M). With opposing currents M is negative and nearly −L<sub>1</sub>, "
        "so both connections are left with only the leakage between neighbouring strands. kronspice's K element takes "
        "only positive coupling, so each returning section's inductor is entered with its terminals reversed (the dot "
        "convention); the wire current is unchanged.")

    r.h1("Resonance and coupling")
    r.table([
        ["Winding", "f<sub>1</sub>", "Wire / λ", "AC gain", "Q", "M to primary", "Strand capacitance"],
    ] + [[lab, eng(run[t]["f1"], "Hz") + (" (peak)" if run[t]["weak"] else ""), f"{run[t]['wl']:.3f}",
          f"×{run[t]['gain']:.3g}", f"{run[t]['Q']:.0f}", eng(run[t]["M"], "H"), eng(run[t]["Ctt"], "F")]
         for t, lab, _ in COILS], [0.30, 0.14, 0.09, 0.10, 0.07, 0.15, 0.15])
    r.p("\"Strand capacitance\" is the two-wire capacitance summed over every pair of neighbouring turns; what it stores "
        "depends on the voltage between them, which the connection sets (previous table). \"(peak)\": no resonance with "
        "gain above 2 below the λ/10 limit; f<sub>1</sub> is then the strongest peak, used as the pulse rate.")

    r.h1("V and I against frequency")
    rows = [["Winding", "Peak f", "Coil ends", "Wire current", "V / I", "Largest wire current below 1 MHz (at f)"]]
    for t, lab, _ in COILS:
        v = vi[t]
        rows.append([lab, eng(v["f"], "Hz"), eng(v["v"], "V"), eng(v["i"], "A"), eng(v["v"] / v["i"], "Ω"),
                     f"{eng(v['imax'], 'A')} ({eng(v['fi'], 'Hz')})"])
    r.table(rows, [0.30, 0.12, 0.12, 0.13, 0.12, 0.21])
    r.p("Per volt of AC drive on the primary. Lowering f<sub>1</sub> by adding capacitance at fixed L leaves the current "
        "at resonance about the same and lowers the voltage: at equal stored energy ½LI² = ½CV², so I is unchanged and "
        "V = I·√(L/C) falls by √(C ratio). Compared at equal voltage instead, the bifilar carries √(C ratio) more current "
        "and stores the C ratio more energy.")
    r.figure(os.path.join(HERE, "vi_frequency.png"),
             "Figure 1. AC sweep per volt of drive on the 5-turn primary. Top: voltage between the coil's ends. Middle: "
             "largest wire current in any section. Bottom: their ratio. Dotted lines mark each winding's peak.", 0.95,
             max_height_mm=230)

    r.h1("Pulse drive at f<sub>1</sub>: volts, amps, power")
    rows = [["Winding", "Coil ends", "Neighbouring turns", "Wire current", "Primary current", "Power in", "Into heat"]]
    for t, lab, _ in COILS:
        x = at_f1(run[t])
        rows.append([lab, eng(x[2], "V"), eng(x[3], "V"), eng(x[4], "A"), eng(x[5], "A"), eng(x[6], "W"), eng(x[7], "W")])
    r.table(rows, [0.30, 0.11, 0.12, 0.12, 0.12, 0.11, 0.12])
    r.p("Peaks in the last 10 % of each run; power is the mean over the run. The pulse width is 0.1/f<sub>1</sub>, so "
        "low-f<sub>1</sub> windings take wider pulses and more power, most of it as heat in the source's 50 Ω and the primary.")
    r.figure(os.path.join(HERE, "coil_vip.png"),
             "Figure 2. Volts, amps and power with pulses at each coil's own f<sub>1</sub>, all seven coils of the example "
             "(log scales for the voltages, linear for currents and power).", 0.98)
    r.figure(os.path.join(HERE, "coil_amps.png"),
             "Figure 3. Wire current (largest section) and primary current building up pulse by pulse at f<sub>1</sub>.", 0.98)
    r.figure(os.path.join(HERE, "coil_current.png"),
             "Figure 4. Strands in series. Top: every turn coloured by the current circulating in it at the moment of "
             "largest current, red one way round and blue the other. Bottom: each section's current along the wire and in "
             "the sense of its turns. In the opposite-current winding the two halves circulate against each other.", 0.98)
    r.figure(os.path.join(HERE, "coil_current_parallel.png"),
             "Figure 5. Strands in parallel, drawn as Figure 4. Sections 0–9 are strand 1, 10–19 strand 2.", 0.98)
    r.figure(os.path.join(HERE, "coil_voltage_counter.png"),
             "Figure 6. Voltage pattern at f<sub>1</sub>: series bifilar with currents the same way and opposing.", 0.95)
    r.figure(os.path.join(HERE, "coil_voltage_parallel.png"),
             "Figure 7. Voltage pattern at f<sub>1</sub>: parallel bifilar with currents the same way and opposing.", 0.95)

    r.h1("Energy ledger at f<sub>1</sub>")
    rows = [["Winding", "In", "Heat", "Stored", "Time-stepping error"]]
    for t, lab, _ in COILS:
        e = run[t]["ledger"]
        rows.append([lab, eng(e[0], "J", 5), eng(e[1], "J", 5), eng(e[2], "J"), f"{e[3]:.2g} J"])
    r.table(rows, [0.34, 0.16, 0.16, 0.16, 0.18])

    tw = parse_twin()
    oa, ja = tw[("open", "ac1")], tw[("joined", "ac1")]
    ot, jt = tw[("open", "tran")], tw[("joined", "tran")]
    zth = oa["ends"][0] / ja["join"][0]
    wl = 2 * math.pi * b["f1"] * b["L"]
    r.h1("The two ends joined, checked against ngspice")
    r.p("The Tesla bifilar solenoid (series, same direction) with its two ends joined through 1 mΩ, against the same "
        "coil with its ends open, driven through the primary as before. kronspice writes both circuits as SPICE "
        "netlists and ngspice-39 (the kronspice-tools image) runs them: AC at 10 kHz, 41.61 kHz and 100 kHz per volt "
        f"of drive, and one-sign 100 V pulses at the open coil's f<sub>1</sub> = 41.61 kHz for 151 periods (peaks over the "
        "last 10 %). \"Wire\" is the current in the middle section (section 9), \"join\" the current through the 1 mΩ.")
    rows = [["Case", "Quantity", "kronspice", "ngspice", "Rel. difference"]]
    names = {"ends": "coil-end voltage", "join": "current through the join", "wire": "mid-coil wire current",
             "prim": "primary current"}
    units = {"ends": "V", "join": "A", "wire": "A", "prim": "A"}
    for (tag, run_), vals in tw.items():
        if run_ not in ("ac1", "tran"):
            continue
        case = f"{tag}, " + ("AC 41.61 kHz (per V)" if run_ == "ac1" else "pulses at 41.61 kHz")
        for q in ("ends", "join", "wire", "prim"):
            if tag == "open" and q == "join":
                continue
            kv, nv, rel = vals[q]
            rows.append([case, names[q], eng(kv, units[q], 5), eng(nv, units[q], 5), f"{rel:.1e}"])
            case = ""
    r.table(rows, [0.27, 0.25, 0.17, 0.17, 0.14])
    worst_ac = max(v[q][2] for (_, rn), v in tw.items() if rn.startswith("ac") for q in v)
    r.bullets([
        f"<b>Both simulators agree.</b> Every AC value at the three frequencies matches to {worst_ac:.0e} or better "
        "(the precision ngspice prints). The pulse runs agree to "
        f"{100 * max(ot['ends'][2], ot['wire'][2]):.1f} % with the ends open and "
        f"{100 * max(jt['join'][2], jt['wire'][2]):.0f} % joined, where the join carries microvolts and the two solvers "
        "take different time steps.",
        f"<b>Joining the ends lowers the voltage and the current.</b> Under pulses the wire current falls from "
        f"{eng(ot['wire'][1], 'A')} to {eng(jt['wire'][1], 'A')} (ngspice), and the coil-end voltage from "
        f"{eng(ot['ends'][1], 'V')} to {eng(jt['ends'][1], 'V')}. The joined coil is a shorted transformer secondary: "
        "its loop current is (M/L)·I<sub>primary</sub>, set by flux linkage, and the 41.6 kHz resonance disappears.",
        f"<b>Why:</b> at resonance the open coil is a high-impedance source. Its Thevenin impedance is "
        f"V<sub>open</sub>/I<sub>joined</sub> = {eng(oa['ends'][0], 'V')} / {eng(ja['join'][0], 'A')} = {eng(zth, 'Ω')}, "
        f"the tank's Q·ωL = {b['Q']:.0f} × {eng(wl, 'Ω')} = {eng(b['Q'] * wl, 'Ω')}. Shorting a source of that "
        f"impedance draws little current. The {eng(oa['wire'][0], 'A')} per volt that circulates inside the open coil "
        "at resonance is reactive current between its L and its strand C, not available at the terminals.",
        "<b>Volts traded for amps</b> needs equal power or energy (P = V·I; ½LI² = ½CV²) and an impedance "
        "transformation: a matched load or a step-down winding, such as the example's tuned pickup (3.8 A at "
        "0.77 V from a coil ringing at 3.6 kV), which still delivers less than the source.",
    ])

    r.h1("What this means for use")
    r.bullets([
        "<b>Low-frequency resonance without a capacitor:</b> series, same direction (Tesla). Same L, large built-in C, "
        "lower Q, lower voltage at the same current, V/2 across the enamel.",
        "<b>Higher resonance from the same wire:</b> parallel, same direction. A quarter of the inductance and no strand "
        "capacitance; it rings above the ordinary coil.",
        "<b>No magnetic field, no resonance at the terminals:</b> either opposite-current connection (non-inductive "
        "resistors, heaters, low-field windings). The series fold puts the whole coil voltage across the enamel at "
        "the base. The parallel connection carries an internal loop that an outside field drives, and the loop's "
        "voltage appears between neighbouring turns; it shows nothing at the terminals, so it has to be checked "
        "inside the winding.",
        "No connection adds energy: each changes how L and C share the stored energy and where the current flows.",
    ])

    r.h1("Limits of the model")
    r.bullets([
        "Twenty sections per coil, each equipotential inside; the opposite-current windings' internal modes above a "
        "few hundred kHz are resolved only roughly, and nothing above the λ/10 limit is used.",
        "The bifilar junction and parallel terminal leads, radiation and proximity effect are not modelled.",
        "The AC sweep uses the wire resistance at 1 MHz throughout, so Q away from 1 MHz is approximate; the pulse "
        "runs use the resistance at each coil's f<sub>1</sub>.",
    ])
    r.h1("Reproducing the results")
    r.code("go run ./examples/pulsed_bifilar               # all seven coils, the V–I sweep and the pickup (writes viz/*.csv)\n"
           "go run ./examples/pulsed_bifilar -only-vi      # the V–I sweep alone\n"
           "python3 examples/pulsed_bifilar/viz/render.py\n"
           "python3 examples/pulsed_bifilar/viz/make_current_report.py\n"
           "go test ./examples/pulsed_bifilar\n"
           "# ngspice twin, ends open and joined:\n"
           "KRONSPICE_NGSPICE_TWIN=/tmp/twin go test ./examples/pulsed_bifilar -run TestNgspiceTwinShort\n"
           "docker run --rm -v /tmp/twin:/w -w /w kronspice-tools sh -c \\\n"
           "  'for f in open_ac joined_ac open_tran joined_tran; do ngspice -b $f.cir > $f.log; done'\n"
           "python3 examples/pulsed_bifilar/viz/ngspice_twin.py /tmp/twin > examples/pulsed_bifilar/ngspice_twin_output.txt")
    r.small("Sources: N. Tesla, US 512,340 (1894); J. C. Maxwell, <i>Treatise</i> Vol. II art. 701; F. W. Grover, "
            "<i>Inductance Calculations</i> (1946); Ramo, Whinnery &amp; Van Duzer, ch. 4.")
    r.build(out)


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "pulsed_bifilar_current_report.pdf"))
