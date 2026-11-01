# Kron's moving inductors: frames I and II, and a parametric pump

The milestone M9c.1 example (CLAUDE.md §3a, §3l item 3; results in
[doc/M9c1.md](../../doc/M9c1.md)). One two-phase induction motor is
started from rest twice, in Kron's holonomic frame I (rotor windings on the
rotor, inductances trigonometric in θ) and in his non-holonomic rotor frame
II(b) (constant inductances, speed voltages G·i·ω), and the two runs must
agree. Then the same N machine element, with its coordinate driven at
constant speed, is a time-varying inductance pumping an LC tank, and its
growth is checked against Mathieu's equation. The frames are the ones
[kron_electron](../kron_electron) and [kron_nmr](../kron_nmr) build on, and
the pumped tank is the standard counterpart of the parametric circuits in
[parametric_lcr](../parametric_lcr) (H16) and
[parametric_pump](../parametric_pump).

```
go run ./examples/kron_machine
go test ./examples/kron_machine
KRONSPICE_EXPORT=1 go test ./examples/kron_machine -run TestExport   # netlists, KiCad schematics, results sheets in export/
```

## The construction

**The primitive machine.** Kron's two-phase primitive machine: stator
windings α, β and rotor windings a, b at angle θ, self-inductances
Ls = Lr = 0.12 H and stator–rotor mutuals M·R(θ) with M = 0.11 H, so the
inductance tensor is `[[Ls·I, M·R(θ)], [M·R(θ)ᵀ, Lr·I]]` with R the 2×2
rotation. In kronspice this is a `circuit.MachineSpec` with a first
harmonic in θ (`Primitive` in `main.go`).

**Frame I, holonomic.** The element is solved with L(θ) as written:
Kron's eq. 3, `e_α = R_αβ·i^β + L_αβ·di^β/dt + [βγ, α]·i^β·i^γ`, with the
Christoffel symbols of eq. 4 built from the metric (L(θ), J), and the
torque ½·iᵀ·(∂L/∂θ)·i (*Equivalent Circuits of Electric Machinery*, 1951,
Epilogue pp. 254–255). `kron.Christoffel` and `kron.Geodesic` reproduce
the element's voltages and torque from those symbols to 1e-12 (`kron`
tests, doc/M9c1.md).

**Frame II(b), the rotor frame.** `kron.RotatingFrame` applies the
connection tensor C(θ) of `kron.RotorFrame` to the rotor axes (power
invariant, e′ = Cᵀ·e): L′ = Cᵀ·L·C becomes constant and the rotor rows
gain the speed matrix G = Cᵀ·L·∂C/∂θ, Kron's non-holonomic object, so that
`e = R·i + p(L·i) + G·i·pθ` with torque iᵀ·G·i (Epilogue eq. 9, p. 257).
The program prints the transformed tensors:

```
L′ = [[0.12 0 0.11 0], [0 0.12 0 0.11], [0.11 0 0.12 0], [0 0.11 0 0.12]]
G′ = [[0 0 0 0], [0 0 0 0], [0 0.11 0 0.12], [−0.11 0 −0.12 0]]
```

Speed terms appear only in the rotor rows, as Kron's frame II requires. A
frame that leaves θ in L is refused by `RotatingFrame`.

**The motor circuit** (`Motor`): a balanced two-phase 100 V, 50 Hz supply
(sin and cos sources) through 0.5 Ω per phase, rotor windings shorted
through 0.8 Ω, rotor inertia 0.002 kg·m² as a capacitor on the mechanical
node w, and a friction load of 0.2 mN·m per rad/s as a 5 kΩ resistor on
the same node (the mechanical domain: angular velocity across, torque
through). Synchronous speed is 314.16 rad/s. Transient 1.5 s, 0.5 ms
step, RelTol 1e-5.

**The parametric pump** (`parametric`): one winding with
L(t) = L0·(1 + m·cos 2ω0t), L0 = 1 mH, m = 0.05, made from a second
harmonic in θ with the shaft turned at ω0 by a voltage source on the
mechanical node, so θ = ω0t and the inductance varies at 2ω0. It is in
series with C = 1/(ω0²L0) (1 V initial charge, f0 = 1 kHz) and a resistor
R = damp·(m·ω0/4)·2L0, where damp is the damping R/(2L0) as a multiple of
the pump rate. The tank runs 60 periods at 100 steps per period; the
growth rate is the slope of log(peak) over periods 30 to 59.

## Results

**Start-up in the two frames** (synchronous speed 314.16 rad/s):

| Frame | Final speed | Points | Supply | Heat | Stored (magnetic + kinetic) | Ledger mismatch |
|---|---|---|---|---|---|---|
| I, holonomic, L(θ) | 313.568 rad/s | 7953 | 326.10 J | 227.29 J | 98.75 J | 9.8e-5 |
| II(b), rotor frame, G·i·ω | 313.696 rad/s | 8439 | 326.21 J | 227.38 J | 98.83 J | −1.4e-6 |

Both ledgers close: supply = heat + stored. The test holds the two final
speeds within 1e-3 of each other and the motor between 300 rad/s and
synchronous; over the whole start-up the frames agree on speed within
1.1e-4 and on stator current within 5.8e-4 of peak
(`analysis.TestMachineFrameTwin`, doc/M9c1.md). The mismatch is the
trapezoidal integration error over the sources' breakpoints.

**The pumped tank against Mathieu** (ω0 = 2π·1 kHz, m = 0.05, pump rate
m·ω0/4 = 78.54 /s):

| Damping R/2L0 (× pump rate) | Growth, run | Mathieu m·ω0/4 − R/2L0 | Pump work | Heat |
|---|---|---|---|---|
| 0 | 78.550 /s | 78.540 /s | 0.07876 J | 5.1e-10 J |
| 0.5 | 39.273 /s | 39.270 /s | 1.396 mJ | 0.707 mJ |
| 2 | −78.571 /s | −78.540 /s | 4.38 µJ | 17.1 µJ |

The test holds the growth within 1 % of Mathieu's rate for the first two
rows. The energy the tank gains enters through the pump's mechanical
port: pump work = tank gain + heat, mismatch ≤ 6e-4
(`analysis.TestMachineParametric`). This is the standard account of
parametric excitation (Mandelstam & Papalexi 1932; Kouřil & Vrba §10.2.2,
whose averaged Mathieu exponent σ = −α + √((mω0/4)² − Δω²) reduces to the
column above at zero detuning, see
[doc/kouril-vrba-verification.md](../../doc/kouril-vrba-verification.md)),
and the agent that varies L pays for everything the tank gains. A COP
computed from the capacitor's source alone, as in H16, omits that work.

**Exports.** `export/` holds the frame-I and frame-II motors
(`.tran 0.5m 1.5 uic`) and the damping-0.5 parametric tank
(`.tran 10u 60m uic`) as netlists with the `N` element and
`.model … machine`, KiCad schematics and results sheets. The sheets show
the DC operating point by MNA only: the machine is nonlinear, so Kron's
linear torn solution is not drawn, and the torn Newton solve is checked in
`kron`'s own tests.

## Files

| File | What it holds |
|---|---|
| `main.go` | the primitive machine, the rotor-frame transform, the motor circuit, the parametric pump and the growth-rate fit |
| `main_test.go` | frame twin (speeds agree to 1e-3, near synchronous) and Mathieu growth within 1 % |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |
| `export/` | `kron_machine_motor_frame1`, `kron_machine_motor_frame2`, `kron_machine_parametric` |

## Sources

- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley, 1951,
  Epilogue pp. 253–265: eqs. 3–4 (frame I, Christoffel symbols), eq. 9
  (frame II(b), speed matrix G), part 2 (absolute time derivative). Scan
  in `reference/`.
- G. Kron, "Non-Riemannian Dynamics of Rotating Electrical Machinery,"
  *J. Math. Phys.* 13, 103–194, 1934 — the non-holonomic form the Epilogue
  cites.
- F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits*, Ellis Horwood,
  1988, §10.2.2 (pp. 377–393): Mathieu's equation and parametric
  self-excitation; `reference/kouril-vrba/`.
- L. Mandelstam, N. Papalexi, *Z. Phys.* 73, 223, 1932 — parametric
  excitation.
- [doc/M9c1.md](../../doc/M9c1.md) — the milestone record with the full
  test table.
