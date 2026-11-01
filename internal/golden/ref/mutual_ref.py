# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent references for the K (mutual inductance) golden files.

mutual_ac:   nodal analysis written out by hand for the two loops (primary
             V1-R1-L1, secondary L2 || C2 || R2), solved per frequency.
mutual_tran: exact solution of the state equations
             [L1 M; M L2] d/dt [i1; i2] = [v - R1 i1; -R2 i2]
             by the matrix exponential (ideal step at t = 0).

Run from the repository root:  python3 internal/golden/ref/mutual_ref.py
"""
import math
import os

import numpy as np
from scipy.linalg import expm

OUT = os.path.join("testdata", "golden")


def mutual_ac():
    R1, L1, L2, k, C2, R2 = 50.0, 10e-6, 40e-6, 0.9, 1e-9, 2e3
    M = k * math.sqrt(L1 * L2)
    rows = []
    for j in range(0, 20 * 4 + 1):  # .ac dec 20 10k 100meg
        f = 10 ** (4 + j / 20)
        w = 2 * math.pi * f
        ys = 1j * w * C2 + 1 / R2  # secondary load admittance
        # Loop equations: 1 = (R1 + jwL1) i1 + jwM i2 ; v_s = jwM i1 + jwL2 i2 = -i2 / ys
        A = np.array([[R1 + 1j * w * L1, 1j * w * M], [1j * w * M, 1j * w * L2 + 1 / ys]])
        i1, i2 = np.linalg.solve(A, np.array([1.0, 0.0]))
        vs = -i2 / ys
        vp = 1.0 - R1 * i1
        rows.append((f, abs(vp), math.degrees(math.atan2(vp.imag, vp.real)),
                     abs(vs), math.degrees(math.atan2(vs.imag, vs.real))))
    write("mutual_ac.csv", ["freq", "vm(p)", "vp(p)", "vm(s)", "vp(s)"], rows,
          "loop equations of the coupled pair solved at each frequency")


def mutual_tran():
    R1, L1, L2, k, R2, V = 10.0, 1e-3, 4e-3, 0.8, 100.0, 1.0
    M = k * math.sqrt(L1 * L2)
    Lm = np.array([[L1, M], [M, L2]])
    Li = np.linalg.inv(Lm)
    A = -Li @ np.diag([R1, R2])
    b = Li @ np.array([V, 0.0])
    # x' = A x + b, x(0) = 0  ->  x(t) = A^-1 (e^{At} - I) b
    Ainv = np.linalg.inv(A)
    rows = []
    for j in range(0, 501):
        t = j * 2e-6
        x = Ainv @ (expm(A * t) - np.eye(2)) @ b
        rows.append((t, x[0], R2 * -x[1]))
    write("mutual_tran.csv", ["time", "i(l1)", "v(s)"], rows,
          "exact solution of the coupled state equations (matrix exponential, ideal step)")


def write(name, header, rows, note):
    path = os.path.join(OUT, name)
    with open(path, "w") as fh:
        fh.write("# reference: " + note + " (internal/golden/ref/mutual_ref.py)\n")
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(repr(float(v)) for v in r) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    mutual_ac()
    mutual_tran()
