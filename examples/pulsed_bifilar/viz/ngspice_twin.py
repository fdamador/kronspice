# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Compare kronspice and ngspice on the Tesla bifilar solenoid with its ends open and joined.

    KRONSPICE_NGSPICE_TWIN=<dir> go test ./examples/pulsed_bifilar -run TestNgspiceTwinShort
    docker run --rm -v <dir>:/w -w /w kronspice-tools sh -c \\
        'for f in open_ac joined_ac open_tran joined_tran; do ngspice -b $f.cir > $f.log; done'
    python3 examples/pulsed_bifilar/viz/ngspice_twin.py <dir>

Reads kronspice.txt (written by the test) and ngspice's wrdata files, and
prints both simulators' coil-end voltage, current through the join (V/1 mΩ),
mid-coil wire current (section 9, V/R) and primary current (V across the
source's 50 Ω): per volt of AC drive at 10 kHz, 41.61 kHz and 100 kHz, and the
peaks over the last 10 % of the pulse run at 41.61 kHz.
"""
import sys

RS = 50.0
SHORT = {"open": 0.0, "joined": 1e-3}


def read(path):
    lines = [l for l in open(path).read().split("\n") if l.strip()]
    return [[float(x) for x in l.split()] for l in lines[1:]]


def main(d):
    k, rsec = {}, {}
    for line in open(f"{d}/kronspice.txt"):
        p = line.split()
        if p[0] == "kronspice":
            k[(p[1], p[2])] = {a.split("=")[0]: float(a.split("=")[1]) for a in p[3:] if not a.startswith("f=")}
        elif p[0] == "rsection":
            rsec[p[1]] = {a.split("=")[0]: float(a.split("=")[1]) for a in p[2:]}
    print(f"{'case':16}{'quantity':10}{'kronspice':>14}{'ngspice':>14}{'rel diff':>10}")
    for tag in ("open", "joined"):
        short = SHORT[tag]
        runs = []
        for i in range(3):
            x = read(f"{d}/{tag}_ac{i}.txt")[0]
            c = [complex(x[j], x[j + 1]) for j in range(1, 13, 2)]  # top, n0, m9, n10, s, p
            ends = abs(c[0] - c[1])
            runs.append((f"ac{i}", {"ends": ends, "join": ends / short if short else 0.0,
                                    "wire": abs(c[2] - c[3]) / rsec[tag]["ac"], "prim": abs(c[4] - c[5]) / RS}))
        rows = read(f"{d}/{tag}_tran.txt")
        stop = rows[-1][0]
        ng = {"ends": 0.0, "wire": 0.0, "prim": 0.0}
        for x in rows:
            if x[0] < 0.9 * stop:
                continue
            ng["ends"] = max(ng["ends"], abs(x[1] - x[2]))
            ng["wire"] = max(ng["wire"], abs(x[3] - x[4]) / rsec[tag]["tran"])
            ng["prim"] = max(ng["prim"], abs(x[5] - x[6]) / RS)
        ng["join"] = ng["ends"] / short if short else 0.0
        runs.append(("tran", ng))
        for run, vals in runs:
            for q in ("ends", "join", "wire", "prim"):
                a, b = k[(tag, run)][q], vals[q]
                rel = abs(a - b) / abs(b) if b else 0.0
                print(f"{tag + ' ' + run:16}{q:10}{a:14.6g}{b:14.6g}{rel:10.2e}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
