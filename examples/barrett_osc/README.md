# Barrett's oscillator-shuttle-circuits (Tesla OSC), Figs. 2–5

T. W. Barrett, "Tesla's nonlinear oscillator-shuttle-circuit (OSC) theory
compared with linear, nonlinear-feedback and nonlinear-element electrical
engineering circuit theory", *Annales de la Fondation Louis de Broglie* 16(1),
23–41 (1991); the same disclosure is US Patent 5,493,691 (20 Feb 1996).
Barrett says Tesla's circuits, with energy shuttled to and from several
floating "grounds", are nonlinear in a way that Kirchhoff's and Ohm's laws
cannot describe, mix pump and signal into an idler that obeys the
Manley–Rowe relations, phase-conjugate the signal, and carry the algebra
of quaternions (SU(2)) and higher Cayley–Dickson numbers. Register entry
H36 (CLAUDE.md §1).

```
go run ./examples/barrett_osc
go test ./examples/barrett_osc
```

The circuits are read from the paper's drawings (Fig. 1 is Tesla's own
coil forms; Figs. 2–5 are Barrett's). Neither the paper nor the patent
gives a component value, so every value is declared in the code. Spark
gaps (Fig. 3) use the gap model and solver of
`examples/sparkgap_diodebridge` (copied in `gap.go`, `sim.go`, `cell.go`);
everything else runs in kronspice.

## Fig. 2 — the OSCs A–G as drawn

Pump ~ (10 V, 1 MHz, 10 Ω source), pump coil α 10 µH, shuttle coil β 30 µH,
load box 50 Ω with its coil γ 5 µH, 2 Ω per coil, A's store 1 nF, each
floating terminal (G2, G3) a 10 pF top-load with 1 GΩ insulation leak,
coupled coils k = 0.3 (0.6 for G's pancakes), F's capacitor 1 nF. The
shuttle's terminal receives the signal (1 V, 630 kHz), as the patent says
("oscillator shuttles … are also receive and transmit antennas").

| Fig. | load at fp | load at fs | transient vs AC (fp, fs) | largest mixing product | ledger |
|---|---|---|---|---|---|
| A | 58.7 mA | 2.4 µA | 1e-4, 1e-5 | 1e-14 | 7e-8 |
| B | 62.1 mA | 1.04 mA | 6e-5, 2e-5 | 6e-15 | 1e-8 |
| C | 58.7 mA | 5.0 µA | 1e-4, 2e-5 | 7e-15 | 1e-8 |
| D | 58.3 mA | 14.5 µA | 8e-5, 4e-5 | 1e-14 | 1e-8 |
| E | 54.6 mA | 1.5 µA | 1e-4, 2e-5 | 3e-15 | 7e-9 |
| F | 35.1 mA | 1.5 µA | 2e-4, 4e-5 | 1e-14 | 1e-7 |
| G | 58.5 mA | 9.9 µA | 7e-5, 5e-5 | 1e-14 | 1e-8 |

Every circuit is linear: the transient load current holds exactly the
lines AC analysis gives for each source alone (superposition), and no
mixing product rises above round-off. Coils and capacitors, however their
grounds float, make a linear network, and a linear network makes no idler.
**Bench test:** drive the pump, feed a signal to G2, and look for fp − fs
in the load. The standard prediction is only what the parts' own
nonlinearity gives (ferrite cores, ceramic dielectrics).

**With a nonlinear store.** Fig. 2A with the terminal top-load a varactor,
q = Ct·(v + v²/(2·20 V)), pump 2 V, signal 0.2 V:

- an idler at fp − fs = 370 kHz appears at the terminal (the shuttle's
  transmitting end), 4.1×10⁻⁵ of the pump line;
- Manley–Rowe sums for the power into the varactor: 4×10⁻⁷ and 5×10⁻⁷
  (zero for a lossless nonlinear reactance; harmonics to 4th);
- idler ∝ Vs^1.01 · Vp^1.03: standard χ⁽²⁾ mixing, Es*·Ep. The paper's
  Eq. 24 (patent Eq. 15) has |Es|²·Ep, which would give Vs².
- signal phase +60° moves the idler by −62.7°: phase conjugation.

So the mixing, Manley–Rowe and phase conjugation Barrett describes are
those of any parametric mixer, and they come from the nonlinear element,
not from the floating grounds.

Spectra use a Hann window: a floating terminal charges slowly when a
nonlinear store rectifies (its leak time constant is about a second), and
the β–top-load ring has not quite died out after the warm-up; neither is
periodic in the window.

## Fig. 3 — the switched −, +, − capacitors (spark gaps)

10 kV on 1 nF each, 1 µH into 50 Ω (G), 1 cm air gaps (Paschen 35.5 kV,
so they hold the 20 kV between neighbours), triggered at 0, 100 and
200 ns, hot channels recovering over 1 µs:

- load voltage lobes −6.71, +3.42, −1.37 kV: the drawn − + − pulse;
- ledger: 150 mJ stored = 66.3 mJ to the load + 76.5 mJ in the arcs +
  3.8 mJ left + 3.4 mJ solver damping (residual 1×10⁻¹⁰ J).

This is a lumped frozen-wave generator, the impulse source of later
ultrawideband radar. When a gap fires, the next capacitor shares charge
with the one before it with no inductance between them. That loss goes to
the arcs, which is why the + and second − lobes are smaller.
Transmission-line segments in place of capacitors avoid it.

## Fig. 4 — the coil networks (Q, H, J)

α…ψ = 10…70 µH (1 mΩ each), 1 pF per node, 10 pF per open terminal,
5 MHz.

- **One-terminal networks (A, B(i), C(i)).** Driven at the terminal
  against G1, the admittance is 1.21, 1.36 and 1.62 × that of the
  terminal's own capacitance, and only through the 1 pF at the inner
  nodes. With no inner-node capacitance every coil current is zero. A
  network with one terminal is only its terminal's capacitance, so the
  coils, and with them the α, β, γ "products", act only through stray
  capacitance, or as antennas in an outside field.
- **Two-ports (B(ii) = X, C(ii) = Y, B(ii) with values reversed = Z):**
  - kronspice cascade XY = ABCD(X)·ABCD(Y) to 2×10⁻¹⁶;
  - XY ≠ YX by 41 %, so order matters ("αβ ≠ βα");
  - XYZ = (XY)Z = X(YZ) to 1×10⁻¹⁰, so cascades are associative;
  - det ABCD = 1 to 1×10⁻¹⁶ (reciprocal: always invertible);
  - S†S = I to 2×10⁻⁶ (lossless S is unitary, an element of U(2)).

Order dependence is real, and it is ordinary two-port theory: cascades
multiply as 2×2 matrices, which do not commute. The SU(2) Barrett names is
the group of a lossless two-port's scattering matrix, up to phase.
Cascades are, however, associative and invertible. They cannot carry the
Cayley numbers' non-associativity (box H) or the zero divisors of box J.
Box H's "(αβ)γ ≠ γ(αβ)" is in any case a commutator, not an associator.

## Fig. 5 — the quaternion "idler"

x = x0 + x1·i + x2·j + x3·k with x0 = Es*·e^{−iωs t}, x1 = Es·e^{−iβs t},
x2 = √Ep·e^{iωp t}, x3 = √Ep·e^{−iβp t}. This is a sum: each component holds
only the signal or the pump frequency, with nothing at fp − fs (all
< 10⁻¹⁶). An idler needs a product, that is a nonlinear element, which is
what the Fig. 2A varactor shows. x1, x3 and Eq. 24 also write β·t. β is
a propagation constant (rad/m, the paper's own wave equation), so β·t is
not a phase; the travelling-wave phase is β·z.

## Verdict

- The circuits as drawn are linear and obey Kirchhoff's laws (kronspice
  is built on them, and transient, AC and the energy ledger agree).
- Barrett's mixing and phase conjugation are real effects of a nonlinear
  element, with standard scaling and the Manley–Rowe relations.
- The SU(2) and order effects are ordinary two-port matrix algebra.
- Fig. 3 is a working pulse former.
- The claims beyond this (noise rejection, lower transmission loss from
  "higher symmetry") are untested here: the paper gives no measurement
  and no equation to simulate.
