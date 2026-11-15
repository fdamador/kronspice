# Bewley 1961 in kronspice

L. V. Bewley, *Tensor Analysis of Electric Circuits and Machines*, Ronald
Press 1961. Photographs of the chapters used are in
`reference/bewley_tensor_analysis/` (untracked, with its own README).
This page maps each chapter to the code and tests that implement or check it.

| Ch. | Topic | Code | Tests |
|---|---|---|---|
| 3 | Transformation of coordinates, Z′ = Cₜ·Z·C, partitioned elimination (Eq. 3‑18) | `kron.Transform` (existing) | `kron/bewley_test.go`: `TestBewleyFig31` (exact), `TestBewleyFig36` (exact fractions 3817/506 and −3817/517; the book prints 7.53 and −7.37) |
| 4 | Symmetrical components, power-invariant Fortescue (Eq. 4‑17) | `kron/sequence.go`: `FortescueMatrix`, `ToSequence`, `FromSequence`, `SequenceImpedance`, `PhaseImpedance` | `kron/sequence_test.go`: unitarity and power invariance, Eq. 4‑24 (cyclic → diagonal), Eq. 4‑26 (unbalanced load), Fig. 4‑2 / Eq. 4‑35, Problem 4‑4 line-to-ground fault (3E/(Z₀+Z₁+Z₂)) |
| 5 | Multiwinding transformers from pairwise leakage impedances (Eqs. 5‑4, 5‑9) | `kron/leakage.go`: `PairLeakage`, `LeakageMatrix` | `kron/leakage_test.go`: Eq. 5‑9 = Cₜ·z·C with the magnetizing current eliminated, Eq. 5‑21 (series secondaries), Eq. 5‑42 (forked autotransformer) |
| 6 | Multivelocity traveling waves: modes |L·K·v² − 1| = 0, surge impedance, transition points | `em/mtl.go`: `LosslessModes`, `Reflection` | `em/mtl_test.go`: single line = `TLine`, homogeneous medium (all modes at c, Eq. 6‑15), inhomogeneous pair, Fig. 6‑2 reflected and transmitted waves |
| 7 | Junction-pair (node) networks, Y′ = Cᵀ·Y·C | `kron` node form (existing) | `TestBewleyFig72` |
| 8 | Complete (orthogonal) networks | `kron.OrthogonalFrame` (existing) | `TestBewleyFig83` |
| 14–16 | Lagrange, Maxwell and Maxwell–Lorentz equations; Christoffel symbols, torsion, affine connection | `kron.Christoffel`, `kron.Geodesic`, `kron.RotatingFrame`, `kron.Torsion`, `kron.NonRiemannian` (existing); `device.Machine` | `analysis/bewley_machine_test.go`: `TestBewleyReluctanceTorque` (Problem 15‑3) |
| 19 | Small oscillations and hunting; curvature tensor (Art. 19‑6, Eq. 19‑44) | `kron/curvature.go`: `Curvature` | `kron/curvature_test.go` (flat machine, salient winding = Gaussian curvature of L(θ)dq² + J dθ²); `TestBewleyHunting` (Eq. 19‑54 frequency and decay) |
| 20 | Synchronous machines, phase-frame salient-pole inductances (Eqs. 20‑46, 20‑47) | `circuit.NewSalientPoleSpec` | `circuit/salient_test.go` (Park returns diag(L_d, L_q, L₀) at every θ) |
| 21 | Induction machines, torque at slip s (Eqs. 21‑7, 21‑11) | `device.Machine` (existing) | `TestBewleyInductionTorque` (two-phase motor at s = 0.05: 3.23517 vs 3.23519 N·m) |

## Still open

- **SPICE `T` element in the circuit engine.** `em.LosslessModes` gives the
  modal decoupling and surge impedance a multiconductor `T` card needs;
  the device itself (method of characteristics, one delay per mode, netlist
  card, ledger) is not written yet.
- **Lossy lines** (R, G ≠ 0): Bewley sets them to zero; Paul 2008 Ch. 7–8
  covers them.
- **Synchronous machine with field and dampers** (Art. 20‑1, Eq. 20‑21
  steady state on an infinite bus) and the two-machine hunting of Art. 20‑14.
- **Problem 7‑6** (Fig. 7‑7) has values but no printed answer; a kron
  vs MNA twin on it would add a node-form case.

## Comparison with ECE2 (H48)

Bewley's torsion (Eqs. 16‑3, 16‑12), affine connection (Eq. 16‑16) and
curvature (Eq. 19‑44) are the objects to compare with Cartan torsion, the
spin connection and Cartan curvature in ECE / ECE2 (`em/scalar/ece.go`,
`circuit/ece.go`). The torsion of a rotating frame changes no invariant
(power, torque); the curvature is a tensor and does not vanish by a frame
change. See `doc/ece2_comparison.md` (report) and `doc/ece2_comparison.csv`
(tracking table, one row per pairing).
