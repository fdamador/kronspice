# kronspice and ECE2: similarities and differences

2026-11-14 · Franklin Amador. Exported from the Claude doc
[kronspice and ECE2: similarities and differences](https://claude.ai/artifact/Qpd6FnK5kFBnGokADxzuGS);
the per-row tracking table is `doc/ece2_comparison.csv`.

Kron's machine theory and ECE2 use the same geometry (torsion, connection,
curvature), and in a circuit their terms take the same form, but they read it
differently: in Kron every torsion term comes from a change of reference frame
and moves energy between declared ports, while ECE2 treats the spin connection
as a source of charge, current and energy. kronspice implements both and books
the difference on its energy ledger.

## Side by side

Every geometric object in ECE2 has a counterpart in Kron's machine theory as
Bewley sets it out, and kronspice now has code for each.

| Kron / Bewley 1961 | Where in Bewley | ECE / ECE2 counterpart | kronspice code |
| --- | --- | --- | --- |
| Metric a = L(θ) on the windings, J on the shaft | Eq. 15‑8 | Tetrad-built metric of Cartan geometry | `circuit.MachineSpec`, `device.Machine` |
| Frame transformation C(θ), rotation matrix γ | Eqs. 15‑10, 15‑12 | Tetrad q | `kron.MovingFrame`, `RotatingFrame`, `ParkFrame` |
| Christoffel symbols [αβ, γ] | Eq. 16‑13 | Symmetric (Levi-Civita) part of the connection | `kron.Christoffel` |
| Torsion T = −G of a rotating frame | Eqs. 16‑3, 16‑12 | Cartan torsion, F = A⁽⁰⁾T | `kron.Torsion` (Hoffmann Ω) |
| Affine connection Γ = Christoffel + torsion | Eq. 16‑16 | Spin connection ω | `kron.NonRiemannian` |
| Cutting voltage G·i·pθ | Table 16‑1 | Torsion term s·cω₀L·i in a coil | `device.Machine`; `device/ece.go`, `kron/ece.go` |
| Curvature tensor K | Eq. 19‑44 | Cartan curvature | `kron.Curvature` (new) |
| Invariance of power under every frame change | Eq. 16‑31 = Eq. 14‑9c | Spin connection as a source of charge, current and energy | Energy ledger, `Ledger.SpinWork` |

## Where they agree

In a circuit, both theories add a voltage that is a rate times an inductance
times a current, with no d/dt on the current.

```
Kron (cutting voltage, Table 16‑1):   v = G · i · pθ
ECE torsion (kron/ece.go):            v = s · cω₀ · L · i
```

- **Same form.** pθ is a rotor speed and cω₀ is a rate in 1/s; G and L are
  both inductance matrices. kronspice stamps both as a speed-dependent
  resistance-like term (r_t = s·cω₀L for ECE; ω·G for the machine element).
- **Same geometric origin.** In Kron the G term is the torsion of a
  non-holonomic frame (Bewley Eqs. 16‑3, 16‑12; Hoffmann 1949). In ECE the
  cω₀A term comes from Cartan torsion through the spin connection.
- **Same machinery for the frame.** A tetrad in ECE and Kron's connection
  tensor C(θ) play the same role: both carry quantities between a fixed frame
  and a moving one, and both produce a connection with a torsion part.
- **Both reduce to Maxwell when the extra terms vanish.** A machine at
  standstill has no cutting voltage; ECE with ω = 0 (or the book's free-space
  choice κ = κ₀ = 0) gives Maxwell's fields exactly. kronspice checks both
  cases bit for bit.

## Where they differ

The difference is where the energy goes, and it follows from the symmetry of
the matrix in the extra term.

| Question | Kron / Bewley | ECE / ECE2 as implemented in kronspice |
| --- | --- | --- |
| Symmetry of the extra term | G is not symmetric (for d–q coils G_qd = L_d, G_dq = −L_q), so the cutting term exchanges power between windings | s·cω₀L is symmetric like L, so it acts like a resistance (s·cω₀ > 0) or a negative resistance (< 0) |
| Power of the term | iᵀG i·pθ = shaft torque × speed: it leaves through the mechanical port | iᵀ(s·cω₀L)i = 2·s·cω₀ × stored magnetic energy, with no port to supply or absorb it |
| Ledger | Closes to 1e−9 with the shaft as a declared port (Ch. 14, 19, 21 tests) | Booked as model-supplied, labelled hypothetical (`Ledger.SpinWork`) |
| What the connection is | A frame quantity: change the frame and the torsion changes, power and torque do not (Eq. 16‑31 = Eq. 14‑9c) | A physical field: ECE2 turns κ₀, κ into charge and current, ρ = −2ε₀κ·E, μ₀J = 2(κ₀E/c + κ×B) |
| Curvature | A tensor of the L(θ), J metric; zero for constant inductances, nonzero for salient or slotted iron; it changes the small-oscillation equations (Art. 19‑6) but adds no energy | Part of the field equations through the Cartan identities |

Standard prediction, carried alongside: ω and ω₀ are frame quantities with no
field of their own, E and B obey Maxwell's equations, and any energy gain
needs a declared source.

## What kronspice has computed

The Kron side matches Bewley's results to the precision of the book; the ECE
side is implemented and its energy is booked, with no source found for it.

| Check | Result | Code |
| --- | --- | --- |
| Induction motor torque at s = 0.05 (Eq. 21‑11) | 3.23517 N·m simulated vs 3.23519 formula; ripple 1.6e−8 | `TestBewleyInductionTorque` |
| Hunting frequency (Eq. 19‑54) | 17.204 rad/s, decay 2 s⁻¹; error 7.4e−5 of amplitude | `TestBewleyHunting` |
| Reluctance torque (Problem 15‑3) | Exact to 1e−9, ledger closed | `TestBewleyReluctanceTorque` |
| Curvature of a salient winding | Equals the Gaussian curvature of L(θ)dq² + J dθ² to 1e−12 | `TestCurvatureSalientWinding` |
| Kron frame torsion (Hoffmann 1949 motor) | Frame I, frame II(b) and eq. 5′ agree (376.587 rad/s); without Ω the motor only coasts down | `examples/hoffmann_motor` |
| ECE / ECE2 terms in Kron's frame | Mesh and node forms equal MNA on 200 random circuits to 1e−11; ledger closed | `kron/ece.go` |
| ECE2 κ₀ in vacuum | Acts as a conductivity: energy ∝ e^(−2κ₀ct) to 0.2 %; growth for κ₀ < 0 booked to the spin connection (residual 3e−14) | `em/scalar/spin.go` |
| ECE2 Coulomb vs Ampère law | Not consistent: defect 0.55 on a TM110 mode | `CoulombDefect` |

## Next tests

These would show whether ECE2's spin connection can be read as a Kron frame
effect or needs a source of its own.

1. **Antisymmetric part only.** Split an ECE torsion term into symmetric and
   antisymmetric parts. If only the antisymmetric part survives a declared
   mechanical or frame coordinate, it is Kron's cutting voltage; the symmetric
   part is what the ledger books as model-supplied.
2. **Frame change test.** Put the same machine in frame I and frame II(b) with
   ECE terms on. Kron's torsion changes and power does not; check whether the
   ECE power changes with the frame.
3. **Curvature side by side.** Compute `kron.Curvature` for the Hoffmann motor
   and the curvature implied by ECE's spin connection for the same windings,
   and compare where each is zero.
4. **Bench reading.** For any ECE circuit that shows gain, name the port that
   supplies it (shaft, pump, source); the H48 standard prediction says one
   always exists.

## Sources

- L. V. Bewley, *Tensor Analysis of Electric Circuits and Machines*, Ronald
  Press 1961, Chs. 3–8, 14–16, 19–21 (photographs in
  `reference/bewley_tensor_analysis/`).
- B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Rev. Mod. Phys.* 21,
  535 (1949).
- H. Eckardt, *ECE Unified Field Theory* Vol. 1 (2022), eqs. 4.211–4.214,
  5.14, 6.164–6.190.
- H. Eckardt and Goel, UFT 456 (2024), eqs. 28–31.
- kronspice: `CLAUDE.md` H48, `doc/bewley.md`, `kron/ece.go`,
  `kron/hoffmann.go`, `kron/curvature.go`, `em/scalar/ece.go`,
  `em/scalar/spin.go`.
