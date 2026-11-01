# Milestone M9a — status

M9a from CLAUDE.md §9: `kron` — primitive network, connection tensor, mesh
and node forms (linear DC/AC). Done when they agree with MNA to 1e−9 on all
linear golden circuits.

## Done-when checklist

| Criterion | Status |
|---|---|
| Primitive network | ✅ `kron.Build`: one primitive branch per element (R, C, L, V, I, E, F, G, H), plus a control coil for each voltage-controlled source; nonlinear elements return `ErrNonlinear` |
| Connection tensors | ✅ mesh tensor `C` (fundamental loops of a breadth-first spanning tree from ground) and node-pair tensor `A` (reduced incidence); `Aᵀ·C = 0`, `m = b − N + S`, `p = N − S` checked on 200 random graphs, including several sub-networks |
| Mesh form | ✅ `Z′ = Cₜ·Z·C`, `e′ = Cₜ·e`, solved for mesh currents (Kron 1939 Ch. IV pp. 101–104, Ch. V pp. 119–121) |
| Node-pair form | ✅ `Y′ = Aᵀ·Y·A`, solved for node potentials (Kron 1939 Ch. XIV p. 360) |
| Linear DC and AC | ✅ `kron.OP`, `kron.AC`, returning `output.Result` / `output.ACResult` |
| Agrees with MNA to 1e−9 on all linear golden circuits | ✅ 11 circuits (bridge, controlled_op, divider, isource, ladder, rc_lowpass_ac, rc_step, rl_step, rlc_bandpass_ac, rlc_step, two_sources), both forms, OP node voltages and element currents and every AC point: worst difference 9.8e−16 |
| Beyond the golden set | ✅ 300 random linear circuits with every element kind, at DC and 10 Hz, 10 kHz, 10 MHz, both forms: worst 2.1e−11 against MNA |
| Typed frames (§3a rule 1) | ✅ `Mat[R, C]`, `Vec[F]` with frames Primitive, Mesh, Node, Junction, Apparent: a transform between the wrong frames does not compile |
| No inverse (rule 2) | ✅ every system is factored and solved (`solver`), with two steps of iterative refinement |
| Robustness | ✅ singular networks (voltage-source loops, a controlled source across its own control) return `ErrSingular`; 60 s of fuzzing, no panic; `go test -race` clean |

## What is computed

```
primitive    one linear relation per branch: Σ a·E + Σ b·i = s
             (Kron's E + e = z·(i + I) with no infinite entries)
mesh         i = C·i′            Z′·i′ + Cₜ·Mⱼ·Eⱼ = e′      Z′ = Cₜ·Z·C, e′ = Cₜ·e
             + one relation per junction axis
node-pair    E = A·v             Y′·v + Aᵀ·Bₐ·iₐ = −Aᵀ·J    Y′ = Aᵀ·Y·A
             + one relation per apparent coil
```

A pure mesh formulation cannot carry a current source (its impedance is
infinite) and a pure node formulation cannot carry a voltage source (its
admittance is infinite). Kron's own answer is to keep such axes in the
other kind: along a current source E is unknown and i known, so it is a
junction pair, not a mesh (Kron 1939 p. 411); a branch with no impedance and
a known voltage is an apparent coil (pp. 407–408). The two forms here do
exactly that, with two typed frames for the extra unknowns:

- **Junction** (mesh form): voltages across current sources, controlled
  current sources and control coils.
- **Apparent** (node form): currents of voltage sources and controlled
  voltage sources.

Controlled sources are mutual terms between primitive axes. A
voltage-controlled source gets a control coil (zero admittance) across its
controlling node pair, so its gain couples two primitive axes; a
current-controlled source couples to the axis of the voltage source it
senses.

**Reactive state variables as axes.** At first the inductors were in the
node form's admittance tensor (`Y = 1/(jωL)`) and the capacitors in the mesh
form's impedance tensor (`Z = 1/(jωC)`). The random twin test found the
cost: a 69 nH inductor at 10 Hz has |jωL| = 4.3 µΩ, and its current
recovered from two nearly equal node potentials carried an error of
ε·|v|/|ωL| ≈ 3e−11 A, 6e−9 relative — the node form's conditioning, not a
solver error (iterative refinement did not change it). So the node form
keeps inductor currents as Apparent axes and the mesh form keeps capacitor
voltages as Junction axes, which is also what they are at DC (a shorted
inductor is an apparent coil, an open capacitor a junction pair). After
that change the worst random-circuit difference dropped to 2.1e−11.

**Invariants checked on every random solution:** Kirchhoff's current law
`Aᵀ·i = 0` and voltage law `Cₜ·E = 0`, Tellegen's theorem `Σ E·i = 0`,
and Kron's power invariance between frames, `e·i = e′·i′` (Kron 1939
eq 4.24).

## Speed

RC ladders, operating point, dense complex arithmetic throughout:

| Sections | mesh | node | MNA (`analysis.OP`) |
|---|---|---|---|
| 10 | 70 µs | 46 µs | 62 µs |
| 40 | 512 µs | 316 µs | 356 µs |

The primitive tensors are dense b × b; that is fine for evolved circuits of
tens of branches. The point of `kron` is not to beat MNA on one solve but
to make tearing possible (M9b), where a mutation re-solves only one torn
piece.

## Sign conventions

Each branch runs n+ → n− of its element: E = V(n+) − V(n−), current
n+ → n− through it, as the SPICE branch current. Kron's impressed emf is
`e = −s` for a voltage source of value s (E + e = 0). The power invariants
are contractions without conjugation (Kron's e·i); complex power is ½·E·i*.

## Not done here

- The energy ledger in Kron's frame (open paths to the exterior, the
  `E·I` output invariant): M9a.1.
- Tearing (M9b), nonlinear circuits (M9c), the GA frame (M9e).
- Sparse primitive tensors.
