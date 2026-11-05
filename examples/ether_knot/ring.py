"""Ringing frequencies of the struck ether knot (run `go run ./examples/ether_knot -ring DIR` first).

    python3 examples/ether_knot/ring.py DIR/ring_b<box>_h<h>.json ...

Each file holds the overlap o(t) of the knot's motion with its starting push, for
three pushes (breathing, squashing, twisting), from modes.go (Ring). The record up
to the time the first echo from the box wall can return is fitted to a sum of
damped cosines (matrix pencil method; Hua & Sarkar, IEEE Trans. ASSP 38, 814, 1990),
giving each component's frequency omega (c/l), amplitude decay rate g (c/l) and
quality factor Q = omega/(2g): the energy falls by e in Q/omega, about Q/(2 pi) cycles.
"""
import json
import re
import sys

import numpy as np


def pencil(y, dt, M=10):
    L = len(y) // 2
    Y = np.array([y[i:i + L + 1] for i in range(len(y) - L)])
    _, _, Vh = np.linalg.svd(Y, full_matrices=False)
    V = Vh[:M].conj().T
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1]) @ V[1:])
    lam = np.log(z.astype(complex)) / dt
    Z = np.vander(z, len(y), increasing=True).T
    amp = np.linalg.lstsq(Z, y, rcond=None)[0]
    rows = [(abs(l.imag), -l.real, abs(a)) for l, a in zip(lam, amp) if l.imag >= -1e-9]
    return sorted(rows, key=lambda r: -r[2])


def analyse(path):
    m = re.search(r"ring_b([\d.]+)_h([\d.]+)\.json", path)
    box = float(m.group(1))
    d = json.load(open(path))
    t = np.array(d["t"])
    dt = t[1] - t[0]
    tmax = 2 * (box - 3) - 6  # before the first echo from the wall comes back to the knot
    out = {}
    for name in ("breathing", "squashing", "twisting"):
        o = np.array(d[name]) / d[name][0]
        n = int(tmax / dt)
        stride = max(1, int(0.12 / dt))
        out[name] = pencil(o[:n:stride], stride * dt)
    return box, tmax, out


if __name__ == "__main__":
    for path in sys.argv[1:]:
        box, tmax, out = analyse(path)
        print(f"{path}: box {box:g} l, record to t = {tmax:g} l/c")
        for name, rows in out.items():
            top = [r for r in rows if r[2] > 0.05][:3]
            print(f"  {name:9s}", "; ".join(f"omega {w:.3f}, decay {g:.3f}, Q {w / (2 * g):.1f} (amp {a:.2f})" for w, g, a in top))
