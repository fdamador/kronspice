# Ether cavitation: can the nonlinear ether tear? (proposed H47)

A liquid under tension tears into its own vapour, and the cavity later
collapses ([bubble_collapse](../bubble_collapse)). The ether of
[ether_knot](../ether_knot) (H29) cannot do this: its field is a unit
vector everywhere, so there is nothing that could drop to zero. That is
also why its knot holds, because a knot can untie only where the field
vanishes.

This example gives the ether an amplitude, φ = a·n, with an equation of
state V(a), and asks how and when it cavitates. Plan, rules and the
particle-model uses: [CLAUDE.md](CLAUDE.md).

```
go test -short ./examples/ether_cavitation    # everything but the ring runs (~4 s)
go test ./examples/ether_cavitation           # + ring relaxations (~5 s)
go run ./examples/ether_cavitation            # step 1 (~15 s); -part 1…5
```

## The model

```
E = κ2 ∫ [ (∂a)² + a²(∂n)² ] + κ4 ∫ a^p |B_n|² + ∫ V(a)
```

There are two equations of state, both written in u = a²:

| | V | The cavity a = 0 |
|---|---|---|
| Q (quartic) | (λ/4)(u − 1)² | a maximum, not a phase; it exists only where something forces a = 0 |
| S (sextic) | λu(1 − u)² + ε(1 − u)²(1 + 2u) | a second phase, costing ε per volume (ε < 0: the ether itself is metastable) |

Bending n costs κ2a²(∂n)². The cost falls on the ether only, never on a
cavity, so **bending is the ether's tension**: it favours the cavity, as
negative pressure favours vapour.

**Units.** κ2 = 1 and λ = 1 in the run, so lengths are in healing lengths
ξ = √(κ2/λ). For any β = λκ4/κ2², where ξ = ℓ/√β, the results scale as
follows:
- lengths ∝ ξ;
- wall tensions ∝ √β;
- twist rates ∝ √β;
- 3-D barriers ∝ 1/√β at fixed ε/λ.

Step 1 is one-dimensional and spherical. The 1-D amplitude problem is
solved as Kron's 1945 ladder network: node voltages are amplitudes,
resistors carry the stiffness, nonlinear shunts carry V′. Its stationary
points are exactly the operating points of that network. The same network
is also built as a kronspice circuit, the twin.

## Results (step 1, `run_output.txt`)

**1. Surface tension.** The Q cavity wall (a forced to 0 on a sheet) is
tanh(x/2ξ) with σ = ⅔√(κ2λ). The S coexistence wall is
a² = 1/(1 + e^(−2x/ξ)) with σ = ½√(κ2λ). The ladder reproduces both, the
error falling as h²: 5×10⁻⁶ and 1.3×10⁻⁵ at h = ξ/40.

**2. Tensile strength.** Twist n at rate k. The ether thins to a0(k) and
carries the twist stress j = 2κ2a0²k. The uniform twist is stable only
while ½V″(a0) > 3κ2k² (Eckhaus; Langer–Ambegaokar), which is the same
point as the maximum of j.

| Ether | Spinodal k_c (1/ξ) | Tensile strength j_c (√(κ2λ)) | Binodal k_B (1/ξ) |
|---|---|---|---|
| Q | 0.408248 = 1/√6 (exact) | 0.544331 = 4/(3√6) | none |
| S, ε = 0.01 λ | 0.5472 | 0.8620 | 0.1001 |
| S, ε = 0.05 λ | 0.5855 | 0.9181 | 0.2249 |
| S, ε = 0.20 λ | 0.7123 | 1.1041 | 0.4551 |

- In the S ether, cavitation becomes possible at the binodal (≈ √ε) and
  barrierless at the spinodal. Between the two it needs nucleation: the
  binodal and spinodal of a liquid.
- **A ring relaxed past k_c tears.** A 400-node ring with a fixed noise
  seed, at k = 1.08 k_c:
  - the amplitude falls to 0.008 at a point at t = 429;
  - the winding drops 14 → 6 in phase slips;
  - the energy goes 31.1 → 6.85, with the dissipation booked to 1×10⁻⁴.
- At 0.92 k_c the ring holds (winding 12, min a = 0.85).
- **An S² ether does not tear this way.** With n free to leave its plane,
  as in Faddeev's ether, the twist unwinds over the pole with the
  amplitude never below 0.96 (winding 6 → 0, no cavity). Tension in an S²
  ether needs topology that forbids the escape: antipodal ends, or a Hopf
  knot.

**3. A reversal squeezed into a slab.** n is held at +n̂ at one end and −n̂
at the other, a length L apart. The ether can bend through half a turn, or
pass straight through zero (a cavity sheet).

| Ether | Bending ether tears below L_t | Cavity sheet heals above L_c |
|---|---|---|
| Q | 4.9337 ξ | 4.9337 ξ (continuous: a pitchfork) |
| S, ε = 0.05 λ | 5.08 ξ | > 30 ξ |
| S, ε = 0.20 λ | 4.12 ξ | 12.65 ξ |
| S, ε = 0.50 λ | 3.16 ξ | 6.33 ξ |

- The S ether shows **hysteresis**: once torn, the sheet stays torn until
  the slab is opened past L_c, as a cavity in a liquid does.
- This is the case that applies to the knot: an S² ether cannot escape
  antipodal ends.

**4. Critical cavity bubble (S, 3-D).** The bubble is found by shooting
and then Newton on the radial ladder. It is a saddle with exactly one
unstable mode, and it tends to the thin-wall forms of classical nucleation
(Gibbs) and Coleman's static bubble:

| ε | R*/R_thin | E*/E_thin |
|---|---|---|
| −0.10 λ | 0.943 | 0.697 |
| −0.05 λ | 0.973 | 0.849 |
| −0.02 λ | 0.990 | 0.940 |
| −0.01 λ | 0.995 | 0.970 |

Under twist (ε = 0.05 λ, taking the twist as a uniform local term — a
labelled approximation; the full 3-D twisted bubble is step 3):
- the barrier falls from 1413 at 20 % of the way from binodal to spinodal
  to 7.9 at 95 %;
- R* falls from 27 ξ to 3.1 ξ.

**5. Kron twin.** The ladder built as a kronspice circuit (resistors and
PWL B elements), solved by OP:
- the Q wall: nodes to 1.1×10⁻⁷;
- the critical bubble, a saddle with 336 nodes: nodes to 7×10⁻⁶, barrier
  to 7 digits;
- the kronspice ledger balances in both. The negative-slope shunts are
  booked as model-supplied power.

**What this says about the knot.** The two 1-D tearing criteria agree on
the gradient at which the ether tears:
- the local one: (∂n)²ξ² > ½;
- the slab one: a half turn over less than 4.93 ξ, i.e. (∂n)²ξ² > 0.41.

`ether_knot`'s relaxed Q = 1 knot reaches (∂n)² = 16.7/ℓ² at its centre
on the axis (one coarse run). So its core should cavitate for β below
about **33–41**, opening first on the axis. Step 2 computes this directly
with the knot in the soft ether.

## Standard counterparts

- **Liquids:** Blake threshold, classical nucleation, Rayleigh collapse.
- **Superfluids and superconductors:** vortex cores at the healing length;
  phase slips past the depairing current (Langer & Ambegaokar, *Phys. Rev.*
  164, 498, 1967).
- **The QED vacuum:** Schwinger pair production at E_s = 1.32×10¹⁸ V/m,
  not yet reached.
- **False-vacuum decay:** Coleman, *PRD* 15, 2929, 1977.

The ether parameters κ2, κ4, λ and ε are not fixed by any measurement, so
every result is a ratio.

## Files

| File | What |
|---|---|
| `eos.go` | Q and S equations of state, twist term, ether amplitude a0(k) |
| `ladder.go` | Kron ladder (planar and radial): residual, energy, Newton, Sturm inertia |
| `vec.go` | two-component ladder: block Newton and block LDLᵀ inertia |
| `wall.go` | walls and the σ quadrature |
| `twist.go` | uniform twist, Eckhaus spinodal, binodal, ring relaxation with ledger |
| `reversal.go` | reversal slab: straight and turned branches, L_c and L_t |
| `radial.go` | critical cavity bubble: shooting, radial Newton, inertia |
| `kron.go` | kronspice circuit twin |
| `main.go`, `main_test.go`, `run_output.txt` | |
