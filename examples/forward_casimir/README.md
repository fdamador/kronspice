# Forward's vacuum-fluctuation battery as a Kron circuit

R. L. Forward, "Extracting electrical energy from the vacuum by cohesion
of charged foliated conductors," *Physical Review B* 30, 1700 (1984).

Two metal plates very close together are pulled toward each other by the
**Casimir force**, which comes from the electromagnetic zero-point
fluctuations of empty space. Forward's battery works like this:

- a stack of thin aluminium leaves, each given the same charge, so that
  their electric repulsion pushes them apart;
- a fast feedback ("servo") supply holds that repulsion just below the
  Casimir pull, so the leaves creep together;
- the work the Casimir force does on the closing leaves goes into the
  supply as electrical energy;
- pushing the leaves apart again "recharges" it.

He compares it to a hydroelectric dam: the force is a catalyst that
releases stored energy, here the surface energy of the separated leaves.

```
go run ./examples/forward_casimir -viz examples/forward_casimir/viz   # about 3 seconds
python3 examples/forward_casimir/viz/render.py
go test ./examples/forward_casimir
```

## 1. How strong is the pull?

Forward uses the formulas for perfect mirrors (his eqs. 1–2). Real
aluminium stops reflecting above its plasma frequency. Here the
**Lifshitz** formula is used with aluminium's plasma wavelength of 107 nm
(Lambrecht & Reynaud 2000).

| Gap | Energy, perfect mirrors | Energy, aluminium | Aluminium ÷ perfect | Pull, aluminium |
|---|---|---|---|---|
| 10 nm | 4.3×10⁻⁴ J/m² | 5.9×10⁻⁵ J/m² | 0.14 | 12,800 Pa |
| 20 nm | 5.4×10⁻⁵ J/m² | 1.2×10⁻⁵ J/m² | 0.23 | 1,400 Pa |
| 100 nm | 4.3×10⁻⁷ J/m² | 2.5×10⁻⁷ J/m² | 0.59 | 6.6 Pa |
| 1 µm | 4.3×10⁻¹⁰ J/m² | 4.1×10⁻¹⁰ J/m² | 0.94 | 0.0012 Pa |

With the plasma wavelength set to zero, the calculation reproduces
Forward's eq. 1 to 10⁻⁵.

## 2. One leaf pair as a Kron circuit

Kron's coordinates are the **charge q** on the leaves and the **gap a**,
as in his machines (Hoffmann 1949, [hoffmann_motor](../hoffmann_motor)).
The parts:

- **Capacitance.** Both leaves sit at the supply voltage, so there is no
  field in the gap. Each faces a grounded neighbour, giving
  C(a) = 4ε0/(L − a) per m². The field outside pushes the leaves apart:
  Forward's repulsion.
- **Equations.** Lagrange–Maxwell with Casimir's energy as a potential:
  - r·q̇ + q/C(a) = V
  - ½m·ä = (electric push) − (Casimir pull)
- **Servo.** At every instant it sets V so the push equals the pull, plus
  just enough to move the leaves at a steady 1 mm/s.

Values (declared): 100 nm aluminium leaves, 2 µm to the grounded
neighbours, closing from 1 µm to 20 nm.

| | Perfect mirrors | Aluminium |
|---|---|---|
| Work done by the Casimir force | 54.2 µJ/m² | 12.5 µJ/m² |
| Servo voltage at 20 nm | 64 V | 27 V |
| Energy the supply gives during the stroke | 24.2 mJ/m² | 4.22 mJ/m² |
| Energy stored in the field at 20 nm | 24.2 mJ/m² | 4.23 mJ/m² |
| **Net gain after discharging at the stop** | **54.0 µJ/m²** | **12.3 µJ/m²** |
| Energy per electron returned | 0.071 eV | 0.039 eV |

What happens:

- **During the stroke the supply pays.** The Casimir pull grows as 1/a⁴,
  so the servo must raise the voltage steeply. During the stroke the
  supply puts in hundreds of times the Casimir work, and it is stored in
  the field.
- **The gain appears only at the end.** The leaves must be stopped and
  discharged. All the charge comes back, each electron carrying a few
  hundredths of an eV more than it went in with. Net, the supply gains
  the Casimir work less tiny servo losses.
- **The ledger closes** to about 1 % of the Casimir work and 3×10⁻⁵ of the
  energy flows.

## 3. The same electrical side in kronspice

kronspice's varying capacitor follows the gap of the stroke and is driven
by the servo voltage through the supply resistance.

| | kronspice | own solution |
|---|---|---|
| Supply energy (Al) | 4.2198 mJ | 4.2208 mJ |
| Charge sent (Al) | 3.1715×10⁻⁴ C | 3.1721×10⁻⁴ C |
| Work of the moving leaves on the capacitor | 12.454 µJ | Casimir work 12.454 µJ |

## 4. A full cycle: close, then reopen

| | Closing gain | Reopening cost | Cycle |
|---|---|---|---|
| Perfect mirrors | +54.0 µJ/m² | −54.3 µJ/m² | **−0.29 µJ/m²** |
| Aluminium | +12.3 µJ/m² | −12.6 µJ/m² | **−0.28 µJ/m²** |

The Casimir force is conservative: reopening costs what closing gave, and
the servo's losses make every cycle a small net loss. Forward's battery is
a one-shot release of stored surface energy, like a spring — not a
continuous source. This agrees with Cole & Puthoff (*Phys. Rev. E* 48,
1562, 1993), who found the idea consistent with thermodynamics in
principle.

## 5. Scale, and electron pumps

- **Energy per kilogram.** One stroke of real aluminium yields 0.023 J
  per kg of leaf metal (100 nm leaves, two per gap). A lithium-ion cell
  stores about 900,000 J/kg.
- **The upper limit.** Pressing the leaves into full contact cannot
  release more than their surface energy, about 1.1 J/m² per aluminium
  surface (Tyson & Miller 1977). That is the same energy that had to be
  spent to make the separate leaves.
- **As an electron pump.** One stroke is a variable-capacitor charge pump
  driven by the Casimir force: it lifts about 2×10¹⁵ electrons per m² by
  0.04 eV each, once. A continuous electron pump needs a continuous source
  of work, as the pumps in [photon_electron](../photon_electron) have.
- **What it is good for.** The useful idea for pump circuits is the
  servo-controlled variable capacitor itself — the circuit in parts 2–3 —
  and the Casimir force as a nanoscale actuator. It is not a source of
  energy.
- **Measurements.** The Casimir force between parallel metal plates has
  been measured (Bressi et al., *Phys. Rev. Lett.* 88, 041804, 2002,
  0.5–3 µm, 15 %). Forward noted that no measurements yet existed for
  metal plates at the 1–20 nm spacings his battery needs.

Picture: `viz/forward_casimir.png`.

## Files

| File | What it holds |
|---|---|
| `casimir.go` | Forward's eqs. 1–2; Lifshitz energy with the plasma model; tables |
| `leaf.go` | the leaf pair as a Kron (q, a) system under Forward's servo, with its ledger |
| `circuit.go` | the electrical side in kronspice |
| `main.go` | the report |
| `main_test.go` | perfect-mirror limit, stroke ledger and twin, cycle |
| `viz/render.py` | `forward_casimir.png` |

## Sources

- R. L. Forward, *Phys. Rev. B* 30, 1700 (1984).
- H. B. G. Casimir, *Proc. K. Ned. Akad. Wet.* 51, 793 (1948); E. M. Lifshitz, *Sov. Phys. JETP* 2, 73 (1956).
- A. Lambrecht, S. Reynaud, "Casimir force between metallic mirrors," *Eur. Phys. J. D* 8, 309 (2000), arXiv:quant-ph/9907105 — plasma wavelengths (Al 107 nm).
- D. C. Cole, H. E. Puthoff, "Extracting energy and heat from the vacuum," *Phys. Rev. E* 48, 1562 (1993).
- G. Bressi, G. Carugno, R. Onofrio, G. Ruoso, *Phys. Rev. Lett.* 88, 041804 (2002), arXiv:quant-ph/0203002.
- W. R. Tyson, W. A. Miller, "Surface free energies of solid metals," *Surface Science* 62, 267 (1977).
