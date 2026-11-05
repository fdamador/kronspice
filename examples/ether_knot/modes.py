"""Lowest vibration frequencies of the ether knot (run `go run ./examples/ether_knot -modes DIR` first).

    python3 examples/ether_knot/modes.py DIR/ops_b<box>_h<h>.txt ...

Each file holds the knot's Hessian H (triplets) and kinetic form K (a 2 x 2 block
per free node, in the order j = 1 .. NZ-2, i = 0 .. NR-2) from modes.go. The
frequencies are omega^2 = eigenvalues of (H, 2K), in c/l, by shift-invert Lanczos
(scipy.sparse.linalg.eigsh) near zero. For each mode it also prints the share of its
kinetic norm within 3 l of the knot's centre: the knot's own vibrations sit there,
the box's standing waves spread through the whole box.
"""
import re
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def load(path):
    with open(path) as f:
        n, nnz = map(int, f.readline().split())
        tri = np.loadtxt(f, max_rows=nnz)
        kb = np.loadtxt(f)
    H = sp.csc_matrix((tri[:, 2], (tri[:, 0].astype(int), tri[:, 1].astype(int))), shape=(n, n))
    rows, cols, vals = [], [], []
    for q, b in enumerate(kb):
        for x in range(2):
            for y in range(2):
                rows.append(2 * q + x)
                cols.append(2 * q + y)
                vals.append(2 * b[2 * x + y])
    M = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
    return H, M


def positions(box, h):
    NR, NZ = int(round(box / h)), 2 * int(round(box / h)) + 1
    r = []
    for j in range(1, NZ - 1):
        for i in range(NR - 1):
            rho, z = (i + 0.5) * h, -(NZ - 1) / 2 * h + j * h
            r.append(np.hypot(rho, z))
    return np.repeat(np.array(r), 2)


def spectrum(path, k=14):
    m = re.search(r"ops_b([\d.]+)_h([\d.]+)\.txt", path)
    box, h = float(m.group(1)), float(m.group(2))
    H, M = load(path)
    w2, V = sla.eigsh(H, k=k, M=M, sigma=-0.02, which="LM")
    order = np.argsort(w2)
    r = positions(box, h)
    out = []
    for i in order:
        v = V[:, i]
        mv = M @ v
        share = float((v * mv)[r < 3].sum() / (v * mv).sum())
        out.append((np.sign(w2[i]) * np.sqrt(abs(w2[i])), share))
    return box, h, H.shape[0], out


if __name__ == "__main__":
    for path in sys.argv[1:]:
        box, h, n, out = spectrum(path)
        print(f"box {box:g} l, h {h:g} l, {n} dofs")
        print("  omega (c/l):     ", " ".join(f"{w:7.4f}" for w, _ in out))
        print("  share within 3 l:", " ".join(f"{s:7.2f}" for _, s in out))
