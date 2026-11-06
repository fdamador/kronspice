# Mechanical and acoustic coupling (energy domains)

CLAUDE.md §3a.1 items 3 and 7: non-electrical channels with Kron's (and
Branin's) network analogies, joined to the circuit and booked on the ledger.
Built after M9a.1, so the acoustic drive of a single-bubble resonator can
be modeled end to end up to the pressure the bubble sees.

## Domains

Every node belongs to one energy domain, the exchange channel whose across
variable is its potential. Ground is the common reference of all of them
(zero velocity, ambient pressure).

| Domain | Across ("voltage") | Through ("current") | R | C | L |
|---|---|---|---|---|---|
| electrical | voltage [V] | current [A] | resistance | capacitance | inductance |
| mechanical | velocity [m/s] | force [N] | 1/damping [m/(N·s)] | mass [kg] | 1/stiffness [m/N] |
| acoustic | pressure [Pa] | volume velocity [m³/s] | acoustic resistance [Pa·s/m³] | compliance [m³/Pa] | inertance [kg/m⁴] |

The mechanical domain uses the mobility analogy (F. A. Firestone, J. Acoust.
Soc. Am. 4, 1933), the one the exchange channel table already declares:
velocity across, force through, so a mass is a capacitor to ground and
Kirchhoff's current law is d'Alembert's force balance at a point. The
acoustic domain is pressure across, volume velocity through.

Netlist (kronspice extensions):

```
.domain mechanical v
.domain acoustic p
Y1 e 0 v 0 0.6          * transducer: i1 = g·V2, i2 = −g·V1
RA p 0 6.9e8 EXT=acoustic * absorbed power leaves to the acoustic exterior
```

Builder: `Domain(exchange.Mechanical, "v")`, `Gyrator(name, p1, n1, p2, n2, g)`,
`ResistorTo(name, a, b, R, exchange.Acoustic)`.

Validation: an element other than a transducer must not join two domains;
diodes, transistors and B elements are electrical only; a resistor's
exterior channel must carry energy (not neutrino, gravitational or
hypothetical). A transducer port is a controlled current source, so it
gives no DC path: a free mass under a constant force has no operating
point, and the engine says so (a spring to ground supplies the path).

## Transducers

`Y` is an ideal, lossless gyrator, `i1 = g·V2`, `i2 = −g·V1`, so
`V1·i1 + V2·i2 = 0`: power in at one port leaves at the other, inside the
network. In the mobility analogy the classic transducers are gyrators:

- piezoelectric (or electrostatic) element of force factor φ [N/V]:
  motional current φ·v, blocked force φ·V — `g = φ`, with the clamped
  capacitance C0 as an ordinary capacitor across the electrical port;
- piston of area S between a mechanical and an acoustic port: U = S·v,
  F = S·p — `g = S`.

In Kron's primitive network a transducer is two axes coupled by an
antisymmetric mutual admittance (y₁₂ = g, y₂₁ = −g), a non-reciprocal
network; the mesh and node-pair forms carry it like any other element.

## Exterior and ledger

Sources draw on their domain's exterior by default (a force source on the
mechanical exterior, a pressure source on the acoustic one). A resistor's
absorbed power leaves as heat unless `EXT=` sends it elsewhere: the piston's
radiation resistance sends sound out through the acoustic open path. Kron's
exterior junctions then show, per channel, where the energy went: electrical
in at the source, thermal out at dampers and the liquid's losses, acoustic
out as radiated sound.

## Liquids (`fluid/`)

| Medium | ρ [kg/m³] | c [m/s] | Source |
|---|---|---|---|
| H₂O, 300 K | 996.556 | 1501.519 | IAPWS-95 verification values (Wagner & Pruß, JPCRD 31, 387, 2002, Table 7) |
| D₂O, 300 K | 1104.04 | 1403.75 | Herrig et al., JPCRD 47, 043102 (2018), Table 7 |

- `PistonRadiation(a, medium)`: a baffled piston's low-ka radiation load as
  an acoustic inertance 8ρ/(3π²a) in parallel with 128ρc/(9π³a²), matching
  the exact impedance (ρc/S)[(ka)²/2 + j·8ka/(3π)] (Kinsler et al. §7.5)
  to O((ka)²); `PistonImpedance` is the exact form (Bessel J₁ and Struve H₁,
  checked against SciPy).
- `SphericalResonator(R, medium, Q)`: the breathing mode of a liquid sphere
  with a pressure-release wall, f0 = c/(2R) (25.0 kHz for a 3 cm water
  sphere, 23.4 kHz in D₂O), as a parallel acoustic resonator whose across
  variable is the center pressure: Ca = 2R³/(πρc²) from the mode's potential
  energy (checked by numerical integration), Ma = 1/(ω0²Ca), Rpar = Q/(ω0·Ca).

## Checks

| Test | Result |
|---|---|
| Mass–spring–damper driven by 1 N: v = F/(jωm + b + k/jω), peak F/b at √(k/m) | exact to 1e−12; damper heat = F·v/2 |
| Piezo on a resonator: electrical admittance = Butterworth–Van Dyke, L1 = m/φ², C1 = φ²/k, R1 = b/φ² | to 1e−12 across the resonance; electrical in = mechanical heat; transducer absorbs 0 |
| Piezo → piston → water: radiated power through the acoustic path = ½·Re(Z)·\|U\|² (exact Z) | within the network's O((ka)²); ledger balances; H₂O and D₂O |
| Kron mesh and node forms on the electrical–mechanical–acoustic chain | equal to MNA to 1e−9 (each node against its own scale), ledger totals equal, acoustic junction = outflow |
| Netlist round trip of `.domain`, `Y`, `EXT=` | identical |

**Solver change.** The first acoustic circuits were declared singular: the
LU's pivot test compared each pivot with the largest entry of the whole
matrix, and an acoustic inertance of 10⁸ next to a transducer gain of 2e−4
made well-determined pivots look like round-off. The test is now relative
to each column's largest entry (dense, complex dense and sparse), which a
change of an unknown's unit does not affect; singular systems are still
reported (`solver.TestColumnScaledSystems`). Kron's augmented systems are
also equilibrated (rows and columns scaled to unit maximum) before LU. All
golden and GA tests are unchanged by this.

## Example: `examples/acoustic_drive`

Generator (50 Ω) → optional matching inductor → piezo ring (C0 = 4 nF,
φ = 0.6 N/V, its mechanical resonance at the mode) → coupling into the
flask's breathing mode → the mode (3 cm sphere, Q = 1500), in H₂O and D₂O.
For a 1.3 atm center amplitude (the SBSL range, Brenner, Hilgenfeldt &
Lohse, Rev. Mod. Phys. 74, 425, 2002):

| | f0 | Pa per volt | generator EMF | power into the liquid |
|---|---|---|---|---|
| H₂O, no matching | 25.03 kHz (peak 24.27) | 3.2 kPa | 41 V | 7.0 mW |
| H₂O, matched | 25.03 kHz | 69 kPa | 1.9 V | 7.0 mW |
| D₂O, no matching | 23.40 kHz (peak 22.69) | 3.2 kPa | 41 V | 6.7 mW |
| D₂O, matched | 23.40 kHz | 72 kPa | 1.8 V | 6.7 mW |

The power the liquid takes at a given pressure is fixed by the mode
(P²/2Rpar); what matching changes is the drive voltage and where the rest of
the power goes. **The Q, the coupling area and the transducer values are
declared assumptions,** not measurements: a hydrophone calibration of a real
flask sets them, after which the chain is a design tool (and an evolvable
one: matching network, transducer and drive frequency for the least power
at a target pressure).

## Not done here

- The bubble in the engine: Rayleigh–Plesset / Keller–Miksis dynamics,
  collapse heating and light emission stay out of `fluid/` and `device/`.
  The bubble enters as a cited equivalent-circuit model in
  `examples/bubble_collapse` (inertance ρ/4πR ∥ radiation resistance
  ρc/4πR², viscous resistance, gas compliance), coupled to the chain at
  node p by its fundamental admittance.
- Distributed acoustic fields (a resonator's higher modes, the transducer's
  position on the wall): one lumped mode per resonator node for now.
- Thermal networks with temperature as the across variable.
