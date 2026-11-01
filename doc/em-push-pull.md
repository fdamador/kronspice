# Design notes: electromagnetic push and pull in space

Recorded at the user's request (2026-10-25) for future circuit designs.
The conversation covered pushing or pulling against fields in space: the
scalar parts of Maxwell's quaternion equations, the Sun's magnetic field,
and a charged body used for buoyancy. Every design built on these notes
must carry both the energy ledger and the momentum ledger (M9n): any net
force has to show where its momentum goes, as any net energy has to show
where it comes from.

## 1. What the scalar parts are

In Maxwell's quaternion form the scalar parts are S.∇E = −div E = −ρ/ε₀
(charge density) and S.∇B = −div B = 0. Creating "positive scalar over
ambient and negative under" in standard physics means separating charge.
A net charge q feels qE; a dipole feels only a torque in a uniform field
and a force along a gradient. Forces need a field that is actually
present, and momentum conservation forbids an isolated device pushing on
its own fields.

A toroidal or twisted (torsion) magnetic field, static or travelling, has
div B = 0 at every point and instant (∂(div B)/∂t = −div curl E = 0), so
it is not a monopole. It does carry a toroidal (anapole) moment and
magnetic helicity, both measurable (planned in `em/topology.go`, M9o).

The extended-electrodynamics hypothesis H8 (a physical gauge scalar,
van Vlaenderen, Hively; §3k) would add a scalar-field force term. It is
modelled only in hypothetical mode, with its own ledger columns, against
the standard prediction.

## 2. The Sun's magnetic field

The field is vast but weak where a craft would fly: about 10⁻⁴ T at the
photosphere (kG in sunspots), about 5 nT in the interplanetary field at
1 AU, reversing with the 11-year cycle and divided into sectors. A
magnetic dipole m feels F = ∇(m·B); at 1 AU the field varies over AU
scales (≈ 3e−20 T/m), so even m = 10⁶ A·m² feels ≈ 3e−14 N. Pushing on
the field means pushing on the plasma it is frozen into:

- **Electrodynamic tether**: EMF (v×B)·L along the conductor, force
  I·L×B, circuit closed through the plasma by contactors. Useful in
  planetary orbits (Earth ≈ 20–50 µT in LEO; NASA ProSEDS was designed to
  test it); small in the solar wind. Momentum goes to the planet's field
  and ionosphere.
- **Magnetic sail** (R. Zubrin, D. Andrews, 1991): a current loop inflates
  a magnetic bubble that deflects the solar wind (dynamic pressure
  ≈ 1–2 nPa at 1 AU). Momentum from the wind's particles.
- **Electric solar-wind sail** (P. Janhunen, 2004): long positively
  charged tethers (≈ 20 kV) repel solar-wind protons; an electron gun
  holds the potential against the electron current the plasma delivers.

The first is a circuit problem (motional EMF, conductor resistance,
contactor I–V characteristics as nonlinear elements, load) and is the
first propulsion example planned for the engine (M9n).

## 3. A charged body for buoyancy

The idea: a charge held in vacuum can reach far higher voltage than in
air, and could be used as buoyancy against the background charge in
space. What the physics says, with numbers to be checked by the engine
when the plasma current model exists:

- **Vacuum does allow much higher fields than air.** Air breaks down near
  3 MV/m; a clean positively charged surface in vacuum can hold tens of
  MV/m or more (electron field emission is a cathode effect). A 10 m
  sphere carrying 1 C sits at ≈ 900 MV with ≈ 90 MV/m at its surface.
- **Space is a plasma, not a vacuum, and it is neutral overall.** There
  is no background static charge to be buoyant against. The solar wind's
  electric field is motional, E = −v×B ≈ 2 mV/m at 1 AU in the craft's
  frame; in the plasma's own frame it is zero. A charged body gathers a
  sheath of opposite charge (Debye length ≈ 10 m in the solar wind), and
  the force it feels comes from momentum exchanged with deflected
  particles, as in the E-sail, not from qE on an unscreened charge.
- **The plasma discharges the body, and holding the charge costs
  power.** The thermal electron current density at 1 AU (n ≈ 5 cm⁻³,
  T ≈ 10 eV) is ≈ 4e−7 A/m². For a positively charged sphere,
  orbit-motion-limited collection grows as (1 + eV/kT); at 900 MV this
  is an enormous current (OML overestimates at such large sheaths, but
  even orders of magnitude less is impractical). That is why E-sail
  designs use thin wires at ≈ 20 kV.
- **For scale**: 1 C in 2 mV/m gives 2 mN, against the Sun's gravity of
  5.9 mN on 1 kg at 1 AU — comparable on paper, which is what makes the
  idea worth modelling, with the screening and the power bill included.
- **In the atmosphere**, the fair-weather field (≈ 100 V/m, pointing
  down) pushes a negatively charged body up. Floating 1 kg needs
  ≈ 0.1 C, which at the 3 MV/m breakdown limit needs a sphere of about
  17 m radius (whose own mass is far larger). The air's conductivity
  (≈ 1e−14 S/m near the ground) bleeds the charge with a time constant
  of about 15 minutes.

**How the engine will model it**: the charged body is a circuit — its
capacitance to the surroundings, the plasma (or air) current–voltage
characteristic as a nonlinear conductance (Langmuir-probe / OML model,
cited), the charging supply and its power, and the force from the
declared external field and particle flux, with energy and momentum
ledgers. The GA can then search for the geometry and potential that give
the most force per watt, and the ledgers will show whether any candidate
is exploiting a modelling gap.

## Sources to cite when these become models

- R. Zubrin, D. Andrews, "Magnetic sails and interplanetary travel,"
  J. Spacecraft and Rockets 28(2), 1991.
- P. Janhunen, "Electric sail for spacecraft propulsion," J. Propulsion
  and Power 20(4), 2004.
- Electrodynamic tethers: NASA ProSEDS mission documentation; L. Johnson
  et al., tether propulsion reviews.
- Orbit-motion-limited current collection: H. M. Mott-Smith,
  I. Langmuir, Phys. Rev. 28, 1926.
- Solar-wind and interplanetary-field values: to be taken from a cited
  dataset (e.g. OMNI) before use.
