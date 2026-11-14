#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""
quantum_circuit_demo.py -- the two simulations from the phone call, self-contained.

Part A  Hydrogen as a resonant ladder circuit (Kron, Phys. Rev. 67, 39, 1945).
        The radial Schrodinger equation  -1/2 u'' + V_l(r) u = E u  (atomic units)
        is built as an LC network in its dual form (as in examples/kron_schrodinger):
          * series capacitor 1/(2h^2) between neighbouring nodes   -> kinetic energy
          * shunt capacitor  V_l(r_i) + V0 at each node              -> potential energy
          * shunt inductor   1 H at each node                        -> the energy
          * short circuit (u = 0) at r = 0 and at the outer edge
        Node voltage = wavefunction u(r) = r*psi(r).
        Each energy level is a resonance at omega^2 = 1/(E + V0).  An AC sweep of the
        driving-point reactance finds it; a matrix eigen-solve is the twin check.

Part B  "Quantumtricity" (H22): is circuit electricity the quark model in disguise?
        Step 1  U(1): on a two-mesh RLC circuit, Kirchhoff's laws are the identities of
                an abelian lattice gauge theory (KVL = zero loop holonomy, KCL = Gauss).
        Step 2  SU(3): put a 3x3 colour matrix on every branch of the same graph and
                compare loop composition order (W1 W2 vs W2 W1) with what every circuit
                quantity does.

Only numpy and matplotlib are needed.   python3 quantum_circuit_demo.py [outdir]
"""
import os
import sys
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
HARTREE_EV = 27.211386245988          # CODATA 2018
HC_EV_NM = 1239.841984332
MU = 1.0 / (1.0 + 1.0 / 1836.15267343)  # reduced-mass factor for hydrogen

# ----------------------------------------------------------------------------------
# Part A: the ladder network
# ----------------------------------------------------------------------------------

class Ladder:
    """Kron's radial network for angular momentum l on a grid of spacing h out to Rmax."""

    def __init__(self, h, Rmax, l=0, Z=1.0):
        self.h, self.l = h, l
        self.N = int(round(Rmax / h)) - 1
        self.r = h * np.arange(1, self.N + 1)
        self.V = -Z / self.r + l * (l + 1) / (2 * self.r ** 2)
        self.V0 = max(-self.V.min(), 1.0) + 0.1      # keeps every shunt capacitor positive and E + V0 > 0
        self.Cs = 1.0 / (2 * h * h)                    # series capacitors (F)
        self.Csh = self.V + self.V0                    # shunt capacitors (F)
        self.L = 1.0                                   # shunt inductors (H)
        # nodal capacitance matrix K (tridiagonal): diag, off-diag
        self.d = 2 * self.Cs + self.Csh
        self.e = -self.Cs * np.ones(self.N - 1)

    def omega_of_E(self, E):
        return 1.0 / np.sqrt(np.asarray(E) + self.V0)

    def E_of_omega(self, w):
        return 1.0 / np.asarray(w) ** 2 - self.V0

    def solve(self, w, node, current=1.0):
        """AC nodal solve Y(w) v = I for a 1 A source into `node`.  w may be a vector:
        the Thomas algorithm then runs for every frequency at once.
        Y = j w K + 1/(j w L) I = (j/w)(w^2 K - I/L)."""
        w = np.atleast_1d(np.asarray(w, float))
        a = w[None, :] ** 2 * self.d[:, None] - 1.0 / self.L      # diag of (w^2K - I/L)
        b = w ** 2 * self.e[0] if self.N > 1 else None              # off-diag (constant)
        rhs = np.zeros((self.N, w.size))
        rhs[node] = current
        # forward sweep
        cp = np.empty((self.N, w.size)); dp = np.empty((self.N, w.size))
        cp[0] = b / a[0]; dp[0] = rhs[0] / a[0]
        for i in range(1, self.N):
            den = a[i] - b * cp[i - 1]
            cp[i] = b / den
            dp[i] = (rhs[i] - b * dp[i - 1]) / den
        x = np.empty((self.N, w.size)); x[-1] = dp[-1]
        for i in range(self.N - 2, -1, -1):
            x[i] = dp[i] - cp[i] * x[i + 1]
        # v = Y^-1 I = (w/j) x = -j w x  -> reactance X = -w x[node]
        return -1j * w[None, :] * x

    def reactance(self, w, node):
        return (self.solve(w, node)[node]).imag

    def resonances(self, Emin, Emax, node, n_sweep=4000):
        """AC sweep in energy; a resonance is a pole of the driving-point reactance,
        i.e. a zero of 1/X.  Bracket on the sweep, then bisect on 1/X."""
        Es = np.linspace(Emin, Emax, n_sweep)
        ws = self.omega_of_E(Es)
        X = self.reactance(ws, node)
        g = 1.0 / X
        idx = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        if idx.size == 0:
            return []
        # all brackets bisected together (one vectorized AC solve per step)
        lo, hi, glo = ws[idx].copy(), ws[idx + 1].copy(), g[idx].copy()
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            gm = 1.0 / self.reactance(mid, node)
            same = np.sign(gm) == np.sign(glo)
            lo = np.where(same, mid, lo); glo = np.where(same, gm, glo)
            hi = np.where(same, hi, mid)
        wr = 0.5 * (lo + hi)
        # a pole of X (resonance) leaves |X| huge; a zero of X (anti-resonance) leaves it small
        Xr = np.abs(self.reactance(wr, node))
        keep = Xr > 1e6 * np.median(np.abs(X))
        return sorted(float(e) for e in self.E_of_omega(wr[keep]))

    def eigen_twin(self, k):
        """k-th lowest eigenvalue of K by Sturm-sequence bisection (independent of AC).
        k may be a list: all are bisected together."""
        d, e2 = self.d, self.e ** 2
        ks = np.atleast_1d(k)

        def count_below(x):
            q = d[0] - x; c = (q < 0).astype(int)
            for i in range(1, self.N):
                q = d[i] - x - e2[i - 1] / np.where(q == 0, 1e-300, q)
                c += q < 0
            return c
        lo = np.full(ks.size, self.d.min() - 2 * self.Cs)
        hi = np.full(ks.size, self.d.max() + 2 * self.Cs)
        for _ in range(64):
            mid = 0.5 * (lo + hi)
            up = count_below(mid) > ks
            hi = np.where(up, mid, hi); lo = np.where(up, lo, mid)
        out = 0.5 * (lo + hi) - self.V0                  # E = kappa - V0
        return out if np.ndim(k) else float(out[0])

    def mode_shape(self, E, node):
        """Node voltages driven just off resonance: the voltage pattern is u(r)."""
        v = self.solve(self.omega_of_E(E + 1e-9), node)[:, 0].imag
        v /= np.sqrt(np.sum(v ** 2) * self.h)
        return v * np.sign(v[np.argmax(np.abs(v))])


def richardson(h, E):
    """Three-grid Richardson extrapolation with the order estimated from the data."""
    (h1, h2, h3), (E1, E2, E3) = h, E
    p = np.log(abs((E1 - E2) / (E2 - E3))) / np.log(h1 / h2)
    return E3 + (E3 - E2) / ((h2 / h3) ** p - 1), p


def part_a(lines):
    P = lines.append
    P("=" * 78)
    P("PART A  Hydrogen as a resonant ladder circuit (Kron 1945)")
    P("=" * 78)
    hs = (0.04, 0.02, 0.01)
    Rmax = 60.0
    targets = {0: [(1, "1s"), (2, "2s"), (3, "3s")], 1: [(2, "2p"), (3, "3p")]}
    results = {}
    for l, states in targets.items():
        per_h = []
        for h in hs:
            lad = Ladder(h, Rmax, l)
            node = int(round(1.0 / h)) - 1                  # drive at r = 1 a0
            res = lad.resonances(-0.55, -0.045, node)
            twins = list(lad.eigen_twin(list(range(len(states)))))
            per_h.append((lad, res, twins))
            P(f"  l={l}  h={h:5.3f} a0  {lad.N:5d} nodes  V0={lad.V0:7.2f}  "
              f"AC resonances (Ha): " + ", ".join(f"{x:.8f}" for x in res[:len(states)]))
            P(f"  {'':19s}eigen twin   (Ha): " + ", ".join(f"{x:.8f}" for x in twins)
              + f"   max |AC - twin| = {max(abs(a - b) for a, b in zip(res, twins)):.1e}")
        for k, (n, name) in enumerate(states):
            Es = [per_h[j][1][k] for j in range(3)]
            Ex, p = richardson(hs, Es)
            results[name] = Ex
            exact = -0.5 / n ** 2
            P(f"  {name}: h->0 (order {p:.2f})  E = {Ex:.9f} Ha = {Ex * HARTREE_EV:.5f} eV"
              f"   exact {exact:.9f} Ha   error {Ex - exact:+.1e}")
        P("")
    dE = (results["2p"] - results["1s"]) * HARTREE_EV * MU
    lam = HC_EV_NM / dE
    dEb = (results["3p"] - results["2s"]) * HARTREE_EV * MU
    P(f"  Lyman alpha (2p -> 1s), reduced mass applied: {lam:.4f} nm   (measured 121.5670 nm)")
    P(f"  Balmer alpha (3p -> 2s), reduced mass applied: {HC_EV_NM / dEb:.3f} nm"
      f"   (measured 656.46 nm in vacuum)")
    P("  Remaining gap to measurement: fine structure, Lamb shift (not in this network).")
    return hs, results


def write_spice(path):
    """Coarse ladder as a plain SPICE netlist (normalized units: farads, henries, rad/s)."""
    h, Rmax = 0.25, 20.0
    lad = Ladder(h, Rmax, 0)
    E1 = lad.eigen_twin(0)
    w1 = lad.omega_of_E(E1)
    node = int(round(1.0 / h)) - 1
    with open(path, "w") as f:
        f.write(f"* Kron 1945 hydrogen s-wave ladder, dual LC form, h={h} a0, R={Rmax} a0\n")
        f.write(f"* node n<i> voltage = u(r_i), r_i = i*h.  V0 = {lad.V0:.6f}\n")
        f.write(f"* 1s resonance on this grid: E = {E1:.8f} Ha -> f = {w1 / 2 / np.pi:.9f} Hz\n")
        f.write(f"* (E = 1/omega^2 - V0).  Run the .ac sweep and look for the driving-point pole.\n")
        f.write(f"I1 0 n{node + 1} AC 1\n")
        f.write(f"CB0 n1 0 {lad.Cs:.10g}\n")
        for i in range(lad.N):
            f.write(f"CS{i + 1} n{i + 1} 0 {lad.Csh[i]:.10g}\n")
            f.write(f"L{i + 1} n{i + 1} 0 1\n")
            if i < lad.N - 1:
                f.write(f"C{i + 1}_{i + 2} n{i + 1} n{i + 2} {lad.Cs:.10g}\n")
        f.write(f"CB1 n{lad.N} 0 {lad.Cs:.10g}\n")
        fl, fh = 0.995 * w1 / 2 / np.pi, 1.005 * w1 / 2 / np.pi
        f.write(f".ac lin 2001 {fl:.9g} {fh:.9g}\n.print ac vm(n{node + 1})\n.end\n")
    return E1, w1


# ----------------------------------------------------------------------------------
# Part B: U(1) circuit identities and the SU(3) composition test
# ----------------------------------------------------------------------------------

def haar_su3(rng):
    z = (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    q = q @ np.diag(np.diag(r) / np.abs(np.diag(r)))
    return q / np.linalg.det(q) ** (1 / 3)


def part_b(lines, rng):
    P = lines.append
    P("")
    P("=" * 78)
    P("PART B  Quantumtricity (H22): circuit electricity vs the quark model")
    P("=" * 78)
    # graph: nodes g=0 a=1 b=2 c=3 d=4 ; branches (from, to)
    names = ["V1", "R1", "L1", "C1", "L2", "R2"]
    br = [(0, 1), (1, 2), (2, 3), (3, 0), (3, 4), (4, 0)]
    nN, nB = 5, len(br)
    A = np.zeros((nN, nB))
    for k, (u, v) in enumerate(br):
        A[u, k], A[v, k] = 1, -1
    Cm = np.array([[1, 0], [1, 0], [1, 0], [1, -1], [0, 1], [0, 1]], float)  # meshes
    f = 1e3; w = 2 * np.pi * f
    R1, L1, C1, L2, R2, Vs = 100.0, 10e-3, 1e-6, 20e-3, 50.0, 1.0
    Zb = np.array([0, R1, 1j * w * L1, 1 / (1j * w * C1), 1j * w * L2, R2], complex)
    e = np.array([-Vs, 0, 0, 0, 0, 0], complex)       # v_k = Z_k i_k + e_k
    Im = np.linalg.solve(Cm.T @ np.diag(Zb) @ Cm, -Cm.T @ e)
    i = Cm @ Im
    v = Zb * i + e
    phi = np.zeros(nN, complex)
    phi[1:] = np.linalg.lstsq(A[1:].T, v, rcond=None)[0]   # ground phi_0 = 0

    P("Step 1 (U(1), standard): two-mesh RLC at 1 kHz "
      "(V1 1 V, R1 100, L1 10 mH, C1 1 uF, L2 20 mH, R2 50)")
    rows = [
        ("Curl o Grad = 0         (A . C = 0)", np.abs(A @ Cm).max()),
        ("KVL = zero holonomy      (C^T v = 0)", np.abs(Cm.T @ v).max()),
        ("voltage is pure gauge    (v = A^T phi)", np.abs(A.T @ phi - v).max()),
        ("KCL = Gauss's law        (A i = 0)", np.abs(A @ i).max()),
        ("gauge shift phi+const leaves v", np.abs(A.T @ (phi + 7.3) - v).max()),
    ]
    for name, val in rows:
        P(f"  {name:42s} residual {val:.1e}")
    q = C1 * v[3]; flux1, flux2 = L1 * i[2], L2 * i[4]
    P(f"  mesh currents |I1| = {abs(Im[0]) * 1e3:.4f} mA, |I2| = {abs(Im[1]) * 1e3:.4f} mA;"
      f" capacitor charge |Q| = {abs(q) * 1e6:.4f} uC; inductor flux |Phi_L1| = {abs(flux1) * 1e6:.3f} uWb")
    P(f"  stored energy: C1 {0.5 * C1 * abs(v[3]) ** 2 * 1e6:.4f} uJ, "
      f"L1 {0.5 * L1 * abs(i[2]) ** 2 * 1e6:.4f} uJ, L2 {0.5 * L2 * abs(i[4]) ** 2 * 1e6:.4f} uJ (peak)")
    pin = 0.5 * (-(e[0]) * np.conj(i[0])).real
    pdiss = 0.5 * (R1 * abs(i[1]) ** 2 + R2 * abs(i[5]) ** 2)
    P(f"  power ledger: source {pin * 1e3:.6f} mW, resistors {pdiss * 1e3:.6f} mW, "
      f"difference {abs(pin - pdiss):.1e} W")

    # loop holonomies, based at node c (3) where both meshes meet
    # mesh1 from c: C1 (c->0), V1 (0->a), R1, L1 (b->c) ; mesh2 from c: L2, R2, C1 reversed
    path1 = [(3, +1), (0, +1), (1, +1), (2, +1)]
    path2 = [(4, +1), (5, +1), (3, -1)]

    def wilson(U, path):
        W = np.eye(U[0].shape[0], dtype=complex)
        for k, s in path:
            W = W @ (U[k] if s > 0 else U[k].conj().T)
        return W

    def defect(U):
        W1, W2 = wilson(U, path1), wilson(U, path2)
        return np.linalg.norm(W1 @ W2 - W2 @ W1) / U[0].shape[0]

    P("")
    P("Every circuit quantity as a link variable exp(i x_k) on the same graph:")
    circ = {}
    qty = {"voltage": v, "current": i, "charge (C v)": np.array([0, 0, 0, C1 * v[3], 0, 0]) + 0 * v,
           "flux (L i)": np.array([0, 0, L1 * i[2], 0, L2 * i[4], 0]),
           "complex power v i*": v * np.conj(i)}
    for name, x in qty.items():
        xs = x / (np.abs(x).max() or 1)               # scale to order one so phases are not tiny
        U = [np.array([[np.exp(1j * xx)]]) for xx in xs]
        circ[name] = defect(U)
        P(f"  {name:22s} composition defect |W1W2 - W2W1| = {circ[name]:.1e}")

    P("")
    P("Step 2 (SU(3), the hypothesis): random colour matrices on every branch, 2000 fields")
    d_su3, d_cartan, gauge_tr, gauge_def = [], [], [], []
    for _ in range(2000):
        U = [haar_su3(rng) for _ in range(nB)]
        d_su3.append(defect(U))
        th = rng.uniform(0, 2 * np.pi, (nB, 2))
        Uc = [np.diag(np.exp(1j * np.array([a, b, -a - b]))) for a, b in th]
        d_cartan.append(defect(Uc))
        g = [haar_su3(rng) for _ in range(nN)]
        Ug = [g[u] @ U[k] @ g[vv].conj().T for k, (u, vv) in enumerate(br)]
        gauge_tr.append(max(abs(np.trace(wilson(U, p)) - np.trace(wilson(Ug, p))) for p in (path1, path2)))
        gauge_def.append(abs(defect(U) - defect(Ug)))
    d_su3 = np.array(d_su3)
    P(f"  SU(3) composition defect: mean {d_su3.mean():.3f}, min {d_su3.min():.3f}, max {d_su3.max():.3f}")
    P(f"  diagonal (Cartan, abelian part) SU(3): max defect {max(d_cartan):.1e}")
    P(f"  gauge invariance: Wilson-loop traces change by <= {max(gauge_tr):.1e}, "
      f"defect by <= {max(gauge_def):.1e}")
    P("")
    P("Verdict: every circuit quantity composes like plain numbers (defect exactly 0),")
    P("the colour field does not (defect ~0.75).  Circuits are the abelian U(1) limit of")
    P("a lattice gauge theory; the quark model is the non-abelian SU(3) one on the same")
    P("kind of graph.  H22 stays 'bounded' until some circuit quantity carries an index")
    P("along its branches whose loop products do not commute.")
    return d_su3, circ


# ----------------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------------

def figures(d_su3, circ):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.4))
    lad = Ladder(0.02, 60.0, 0)
    node = int(round(1.0 / 0.02)) - 1
    Es = np.linspace(-0.56, -0.04, 6000)
    X = lad.reactance(lad.omega_of_E(Es), node)
    X[np.abs(X) > 60] = np.nan
    ax1.plot(Es * HARTREE_EV, X, color="#2a6fdb", lw=1.2)
    for n in (1, 2, 3):
        ax1.axvline(-0.5 / n ** 2 * HARTREE_EV, color="#d1495b", ls="--", lw=1)
        ax1.text(-0.5 / n ** 2 * HARTREE_EV, 52, f" n={n}", color="#d1495b")
    ax1.set_xlabel("drive frequency, read as energy E = 1/ω² − V0  (eV)")
    ax1.set_ylabel("driving-point reactance X (Ω, normalized)")
    ax1.set_title("AC sweep of the hydrogen ladder: poles = energy levels")
    ax1.set_ylim(-60, 60)

    r = np.linspace(0, 25, 500)
    exact = {"1s": 2 * r * np.exp(-r),
             "2s": (1 / np.sqrt(2)) * r * (1 - r / 2) * np.exp(-r / 2),
             "2p": (1 / (2 * np.sqrt(6))) * r ** 2 * np.exp(-r / 2)}
    cols = {"1s": "#2a6fdb", "2s": "#e08a1e", "2p": "#2a9d8f"}
    for name, l, k in (("1s", 0, 0), ("2s", 0, 1), ("2p", 1, 0)):
        L = Ladder(0.05, 60.0, l)
        E = L.eigen_twin(k)
        u = L.mode_shape(E, int(round(1.0 / 0.05)) - 1)
        ue = exact[name] * np.sign(exact[name][np.argmax(np.abs(exact[name]))])
        ax2.plot(r, ue, color=cols[name], lw=1.2, label=f"{name} exact")
        sel = L.r <= 25
        ax2.plot(L.r[sel][::6], u[sel][::6], "o", ms=3.5, color=cols[name], label=f"{name} node voltages")
    ax2.set_xlabel("r (Bohr radii)")
    ax2.set_ylabel("u(r) = r ψ(r)")
    ax2.set_title("Voltage pattern at each resonance = the wavefunction")
    ax2.legend(fontsize=8, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "hydrogen_ladder.png"), dpi=150)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(d_su3, bins=50, color="#7b5ea7", alpha=0.85, label="SU(3) colour field (2000 random)")
    for name, val in circ.items():
        pass
    ax.axvline(max(max(circ.values()), 1e-17), color="#d1495b", lw=2,
               label="every circuit quantity (V, I, Q, Φ, P): ≈ 0")
    ax.set_xlabel("loop composition defect  ‖W₁W₂ − W₂W₁‖ / 3")
    ax.set_ylabel("count")
    ax.set_title("Does order matter going round two loops?")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "quark_test.png"), dpi=150)


def main():
    rng = np.random.default_rng(20261112)
    lines = []
    part_a(lines)
    E1, w1 = write_spice(os.path.join(OUT, "hydrogen_ladder.cir"))
    lines.append(f"  SPICE netlist hydrogen_ladder.cir written (h = 0.25 a0): 1s resonance "
                 f"E = {E1:.6f} Ha at f = {w1 / 2 / np.pi:.6f} Hz (normalized units)")
    d_su3, circ = part_b(lines, rng)
    figures(d_su3, circ)
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(OUT, "run_output.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
