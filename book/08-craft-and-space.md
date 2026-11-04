# 8. Craft and space

A circuit in orbit can push against a planet's magnetic field, and a
craft's hull can be a waveguide. The design notes of
[doc/em-push-pull.md](../doc/em-push-pull.md) set one rule for every such
design: beside the energy ledger it carries a **momentum ledger**, so that
any net force shows where its momentum goes, as any net energy shows where
it comes from (CLAUDE.md §3f). This chapter runs five examples under that
rule: the electrodynamic tether as generator and thruster
([tether](../examples/tether)), the same wire coiled up
([coiled_tether](../examples/coiled_tether)), towers driving standing waves
on a closed ring ([ring_towers](../examples/ring_towers)), the ring moved
inside a craft's double wall ([double_hull](../examples/double_hull)), and
one spinning rotor placed in three gravitational fields, where rotation
and curvature are kept apart ([machine_spacetime](../examples/machine_spacetime),
§3l).

## 8.1 Run it: the tether as a generator

A vertical aluminium wire 5 km long moves east at 7.7 km/s through a
northward field of 30 µT. The Lorentz force on its charges gives an EMF
(v×B)·L along it, and the circuit closes through the ionosphere by plasma
contactors at the two ends. The exported netlist of the generator case is
the circuit with that EMF already computed:

```
electrodynamic tether
VEMF top m DC 1155.0000000000002
RT m bot 168.70423967740905
RCT top x 100
RCB bot y 100
RL x y 500
RREF y 0 1e+09
.op
.end
```

`VEMF` is the motional EMF, 1155 V, equal to v·B·L. `RT` is the 5 km of
1 mm wire, 169 Ω. `RCT` and `RCB` stand in for the two contactors at
100 Ω each, a linear substitute for their nonlinear current–voltage
characteristics. `RL` is the load, and `RREF` gives the floating loop a
reference so the operating point exists. Run it:

```
go run ./cmd/kronspice examples/tether/export/tether_generator.cir
```

```
Elements
  name  kind            nodes  voltage     current     power absorbed
  VEMF  voltage source  top→m  1.155 kV    -1.32957 A  -1.53565 kW
  RT    resistor        m→bot  -224.303 V  -1.32957 A  298.226 W
  RCT   resistor        top→x  132.957 V   1.32957 A   176.775 W
  RCB   resistor        bot→y  -132.957 V  -1.32957 A  176.775 W
  RL    resistor        x→y    664.783 V   1.32957 A   883.873 W

  inflow        1.53565 kW
  outflow       0 W
  dissipation   1.53565 kW
  residual      227.374 fW (relative 0)
```

The load receives 884 W; the wire and the two contactors burn the other
652 W; the source delivers 1.536 kW. Where does that come from? The
netlist cannot say, and this is where the Go program begins. The current
feels ∫I·dl×B, which for this loop is 0.199 N of drag along the orbit,
and F·v = −1536 W: the orbit pays for the load, to the watt, since
F·v = −EMF·I. `em.MomentumLedger` books the force on the tether together
with its reaction, +0.199 N on the geomagnetic field and the ionosphere,
which is the exterior source. The test holds all of this to 1e-9.

This is the circuit NASA's *Tethers in Space Handbook* writes for the
electrodynamic tether (Cosmo & Lorenzini 1997, §4): the EMF
`V = ∫(v×B)·dl`, the force `F = ∫I·dl×B`, and the loop equation
`V_ind = I·R + ΔV_low + ΔV_up + ΔV_ion + ΔV_load`, where the ΔV terms at
the ends are the drops across the contactors' space-charge sheaths or
plasma clouds, the largest impedances in the loop, and the ionospheric
drop is of order 1 to 20 Ω. The example's two 100 Ω resistors stand in for
the sheaths with the ionosphere folded in.

## 8.2 The thruster and NASA's record

Drive the loop from a supply larger than the EMF and the current reverses:
with 3000 V the tether carries 5.0 A the other way, the force becomes
0.751 N of thrust, and F·v = +5780 W goes into the orbit while the supply
delivers 15.0 kW. The rest, 9.2 kW, is I²R heat in the 369 Ω of wire and
contactors. The handbook states the same reversal: fed against the induced
EMF, the same tether is a thruster at the expense of on-board power. NASA
flew this circuit and designed it twice more, and the tether README sets
the example beside those numbers:

| | This example | PMG (1993) | TSS-1R (1996) | ProSEDS (designed) | ISS reboost study (1998) |
|---|---|---|---|---|---|
| Tether | 5 km, 1 mm Al, 169 Ω | 500 m, #18 AWG Teflon-insulated Cu | 19.7 km deployed of 20.7 km; 10 strands of 34 AWG Cu, 2.54 mm overall | 5 km of 7-strand 28 AWG Al, mostly bare, plus 10 km Dyneema | 10 km, bare lower section, < 200 kg |
| EMF | 1155 V (231 V/km) | bias ±65 to −130 V in series with the EMF | 3482 V (177 V/km) | — | about 100 V/km |
| Current | 1.33 A generator, 5.0 A thruster | 100–300 mA by day, 10–50 mA at night; 0.3 A at 130 V bias | 0.97 A before the break; 1.1 A for 75 s after | — | — |
| Force | 0.199 N drag, 0.751 N thrust | — | — | about 0.4 N drag; 14 km/day descent | 0.5 N from 5 kW, 0.8 N from 10 kW |
| Ends | two 100 Ω linear contactors | hollow-cathode xenon plasma contactors, current fully reversible | 1.6 m conducting sphere; electron guns and the Orbiter ground | bare wire collects, hollow cathode emits | bare wire collects, plasma contactor emits |

Four things follow. The example's 231 V/km is the upper end, because it
takes v perpendicular to B at 30 µT; TSS-1R measured 177 V/km near the
equator at 297 km and the ISS study quotes about 100 V/km for the
station's inclined orbit, which is why the code says a real orbit must
take B from a field model. The example's thruster gives 0.05 N/kW against
the study's 0.1 N/kW at 5 kW, because NASA's bare tether collects
electrons along its length with no end sheath to cross (Sanmartín,
Martínez-Sánchez & Ahedo 1993) and a hollow cathode emits them at tens of
volts; the linear contactors are the model's weakest part, since TSS-1R's
sphere collected up to three times what the Parker–Murphy model
predicted. PMG demonstrated in flight the reversibility the two columns
assume. And TSS-1R's end is a design lesson: on 25 February 1996, with
19.7 km deployed and the EMF at 3482 V, an arc through a flaw in the
insulation burned the tether through at the Orbiter end, and the current
rose to 1.1 A as the bare broken strands entered the plasma (Gilchrist et
al. 1998; Szalai 1996). Every design here is checked against breakdown
at its peak voltage.

## 8.3 The same wire coiled

Would "multiplication between windings" make a coiled tether an energy
source? [coiled_tether](../examples/coiled_tether) takes a 5-mile
(8,049 m) copper tether of 2 mm wire, 43 Ω, and compares it straight, as an
ordinary coil 2 m across with 1,281 turns, and as the same coil wound
bifilar (Tesla, US 512,340). Energy may come only from the sources named:
Earth's magnetic field through the motion, and the Sun through a 2 m²
array at 30 %.

The motional EMF is computed along the actual wire path with
`em.MotionalEMF`:

| Winding | EMF |
|---|---|
| Straight 5-mile tether | 1,859 V |
| Coiled, ordinary | 0.888 V, the coil acting as a rod of its own 3.84 m length |
| Coiled, Tesla bifilar | 0.889 V, both strands running the same way through the field |
| Straight, folded out and back | 1.1×10⁻¹³ V, out and back cancel |

Coiling throws the tether's voltage away: in a uniform field the turns
cancel, and what is left is the coil's axial length. In the craft's frame
the same coupling is an electric field E′ = v×B = 0.231 V/m, the Lorentz
transformation of chapter 5, not a second source. The field changes only
over an orbit, far below any LC resonance, so Earth's field cannot drive
the coil; driving its 565.7 Hz resonance from the field would need the
coil to spin at 566 rev/s, a rim speed of 3,554 m/s against carbon
fibre's 1,881 m/s burst speed.

The Sun can. 816 W from the array drive the coil at resonance (L 1.36 H,
C 58 nF, Q 112) in a transient run with the ledger on:

| | Ordinary coil + tank capacitor | Tesla bifilar | Bifilar, matched load |
|---|---|---|---|
| Coil voltage | 29,784 V (×112) | 29,784 V (×112) | 21,064 V (×56) |
| Between neighbouring turns | 23 V | 14,892 V | 10,532 V |
| Energy stored between turns | 6×10⁻⁵ J | 25.8 J (all of it) | 12.9 J |
| Sun in = heat (+ load) | 815 W = 814 W | 815 W = 814 W | 815 W = 407 W + 407 W load |

Multiplication is real for voltage, by Q, and the bifilar winding does
put the stored energy between its turns, Tesla's N²/4; but every joule is
the array's, and the whole-run ledger closes: 309.9 J in = 284.1 J heat +
25.8 J stored, with 6 mJ of integration error. The bifilar turn voltage,
1.5×10⁷ V/m across a 1 mm gap, exceeds the declared 10 MV/m limit and is
flagged: it needs insulation or wider spacing. A multipactor scan between
the turns shows no avalanche at 566 Hz, 1 MHz or 100 MHz; it appears at
1 to 3 GHz, where the electrons take their energy from the drive and the
8 km wire is many wavelengths long.

## 8.4 Towers on a ring

[ring_towers](../examples/ring_towers) is a thought experiment as a
circuit: towers at points around the Earth–ionosphere ring drive one
frequency into it, each sending a wave east and a wave west, with phases
chosen so that the waves add, cancel or circulate. The ring is a lumped
transmission line scaled to MHz so that each lap takes microseconds: 96
sections of 1 µH and 0.132 Ω with 1 nF shunt, closed on itself, Z0 =
31.6 Ω, 0.20 Np loss per lap, run at its fourth resonance, 1.314 MHz, four
wavelengths around. Each tower is two 10 mA current sources a quarter
wavelength apart, an end-fire pair whose phasing launches east only or
west only; tower B sits at the antipode. Each state is run to steady state
by transient analysis.

| State | Peak | Envelope min/max | Eastward power flow | Towers deliver = ring loss | Stored |
|---|---|---|---|---|---|
| traveling | 3.294 V | 0.925 | 0.155 W | 0.06255 W = 0.06254 W | 0.475 µJ |
| standing | 3.328 V | 0.018 | 0.000163 W | 0.03136 W = 0.03135 W | 0.238 µJ |
| cancel | 0.167 V | 0.904 | 0.0003 W | 0.000161 W = 0.0001609 W | 0.00122 µJ |

With four wavelengths around the ring every quarter point is an antinode
of the A–B pattern, so a second pair C–D at 90° and 270°, phased to the
same reference waves, puts its pattern on top of the first and the two
add: 6.656 V peak from 0.1252 W, twice the amplitude and four times the
stored energy of one pair. C–D in antiphase cancels everything to 0.323 V;
shifted by an eighth of a wavelength it moves the antinodes off the
towers. Eight synchronized towers give 11.95 to 13.22 V around the ring
from 0.5006 W, sixteen times one pair's stored energy. The ledger
mismatch is at most 1.1×10⁻⁹ of the energy in; every joule in the ring is
the towers'.

A DC bias per pair through 1 MΩ feeds shows what a conductor can and
cannot hold. On a continuous ring the pairs cannot keep different levels:
with A–B at 100 V and C–D at 50 V every quarter point sits at 75.0 V and
50 µA circulates between the supplies. A sectored ring with 10 µF blocking
capacitors at 45°, 135°, 225° and 315° holds 100, 50, 100 and 50 V with no
DC current, and the RF rides unchanged on top.

## 8.5 Inside a double hull

[double_hull](../examples/double_hull) moves the ring inside a craft: the
gap between an inner and an outer conducting hull is a closed
parallel-plate waveguide around the craft, driven at a ring resonance by
radiators in the gap, with a DC bias between the hulls added to the RF.
The gap is the TEM line of the telegrapher's equations (Paul §4.2) with
L′ = μ0·d/w, C′ = ε0·w/d, Z0 = η0·d/w and the copper walls' skin-effect
loss, 24 sections per wavelength closing the ring. The example hull is
5 m in radius and 10 m high with a 0.10 m gap, driven at mode 4 with 1 kW:

| Run | f | Q | RF gap amplitude | Formula | Peak gap (RF + bias) |
|---|---|---|---|---|---|
| traveling, 0 V | 37.685 MHz | 9410 | 1681–1683 V | 1680 V | 1683 V |
| traveling, 2500 V | 37.685 MHz | 9410 | 1681–1683 V | 1680 V | 4183 V |
| standing, 0 V | 37.685 MHz | 9410 | 0–2379 V | 1680·√2 V | 2379 V |
| standing, 2500 V | 37.685 MHz | 9410 | 0–2379 V | 1680·√2 V | 4879 V |

The traveling-wave amplitude at mode n is √(Z0·P·Q/(π·n)), and a standing
wave doubles it at the antinodes. The ledger in every run reads radiators
1000.0 W = wall loss 999.9 W plus the bias feed; 39.7 mJ of RF energy sits
in the gap, and the 2.5 kV bias stores 87.8 mJ once and draws no steady
power.

Two limits bound the gap. A vacuum gap under RF can multipactor: electrons
crossing it in resonance with the field multiply at each impact when the
wall's secondary yield exceeds one (CLAUDE.md §3i). `em/multipactor.go`
tracks electrons launched from either wall at every phase through the
combined RF and DC field, impact after impact, with Furman's yield fit to
the measured copper values of Baglin et al. 2002. At this hull's
f·d = 3.77 MHz·m the worst-phase mean yield per impact is above one at RF
amplitudes of 300 and 1000 V with no bias, so ring-up through those
amplitudes meets a growth of 2.02 at about 760 V; a bias of 500 V clears
the whole ring-up, which is the cure measured on the LEP2 and KEKB
couplers. Kilpatrick's sparking limit at 37.7 MHz is 8.1 MV/m against
0.049 MV/m for the biased standing wave. An evolved design for a 20 kV
peak with at least 5 kV of it RF, least RF power, no multipactor on the
way up and under Kilpatrick, is a 0.470 m gap at mode 1 (9.09 MHz,
Q 21,729), 451 W of RF and 12,496 V of bias, giving 19,963 V peak; its
ledger reads 451.0 W = wall loss 448.3 W + feed 2.7 W.

## 8.6 One rotor in three gravitational fields

The last example leaves circuits for the frame they sit in. Three
different "curvatures" are easily confused (CLAUDE.md §3l): the
gravitational curvature of spacetime, rotating frames in flat spacetime,
and Kron's curved configuration space of a machine, which is the same
mathematics as the first two and none of their physics.
[machine_spacetime](../examples/machine_spacetime) takes one machine, a
0.5 m carbon-fibre disk spinning at half its burst speed and carrying a
1 m² loop, and puts it in an Earth laboratory at 45°, in Gravity Probe B's
642 km polar orbit, and 1000 km above the centre of a neutron star, as an
illustration. By the equivalence principle the hub runs at the
Schwarzschild rate for its position and speed, the rim is special
relativity in the hub's local frame, the spin axis precesses by Thomas,
geodetic and Lense–Thirring terms (Schiff 1960), and the loop sees the
Sagnac delay relative to the local gyroscopes.

| | Earth lab, 45° | GP-B orbit, 642 km | Neutron star, 1000 km |
|---|---|---|---|
| Hub clock rate | 0.999999999304 | 0.999999999052 | 0.99793 |
| Rim loses against hub | 0.43 µs/day | 0.43 µs/day | 0.43 µs/day |
| Geodetic precession | 350 mas/yr | 6,621 mas/yr (orbit average) | indicative only |
| Frame dragging | 172 mas/yr | 40.9 (39.1 against the guide star) | 4.9×10¹¹ mas/yr |
| Burst margin | 2.0 | 2.0 | 2.0 |
| Tidal stress on the disk | 0.003 Pa | 0.002 Pa | 3.3×10⁵ Pa |

The orbit column is the validation: Gravity Probe B measured 6601.8 ±
18.3 mas/yr geodetic and 37.2 ± 7.2 mas/yr frame dragging (Everitt et al.
2011), and the model gives 6621 and 39.1 with a spherical Earth. The rim's
0.43 µs per day is the same at all three sites because it is the rotor's
own speed, 941 m/s, and the burst margin of 2.0 is the material's, long
before any relativistic effect is large. Turned around, the machine's own
mass and spin curve spacetime 1 m away by a clock shift of 2×10⁻²⁶, a
tidal curvature 0.0012 of Earth's, which torsion balances can weigh, and
a frame dragging 3.7×10⁻¹⁰ of Earth's own, far below any instrument. No
energy enters a circuit from spacetime curvature by itself; the
gravitational channel of chapter 1 carries no power.

## 8.7 The Earth as a machine

[earth_machine](../examples/earth_machine) treats the planet itself as
three Kron networks. Layer 1 is the core dynamo as Rikitake's two-disk
machine: two N elements, each disk excited by the other's coil. It
reverses its field on its own, 71 times in t = 50 to 1000 at μ = 1, A = 1,
with books closed to 2×10⁻⁷.

Layer 2 is Wilson's global circuit. The fair-weather air is a stack of
60 leaky shells, each a resistor in parallel with a capacitor. The
thunderstorms are a current source that lifts charge from the ground to
the ionosphere. Calibrated on three measured values (V_I = 250 kV,
J = 2 pA/m², E0 = 120 V/m), the ladder gives a global current of 1,020 A
and a resistance of 244.9 Ω. The Earth's charge comes out at −5.42×10⁵ C,
Gauss's law at the surface. The balancing positive charge is in the air,
half of it below 1.42 km, and almost none of it at the ionosphere. So the
Earth is not a spherical capacitor of 0.076 F: its effective capacitance
is Q/V_I = 2.17 F. The relaxation time R·C_eff is exactly ε0/σ0 = 531 s,
and the storms keep 255 MW flowing. Against a daily ±15 % swing in storm
activity (the Carnegie curve), V_I follows 4.4 minutes behind.

Layer 3 is the Earth–ionosphere cavity as a transmission line in
colatitude, one band per degree. With perfect walls the ladder's modes
are Legendre's c/(2πa)·√(n(n+1)) to 4×10⁻⁴, and the engine's AC sweep
peaks on them. With v/c and the damping fitted to the first Schumann mode
alone (7.8 Hz, Q = 5), modes 2 to 5 fall at 13.6, 19.2, 24.8 and 30.3 Hz,
against the measured 14.3, 20.8, 27.3 and 33.8 Hz. The measured modes
imply a phase velocity that rises from 0.736c to 0.824c with frequency:
one constant ionosphere height cannot give them, which is the reason the
literature uses frequency-dependent electric and magnetic heights. With
those heights the ladder is rebuilt at every frequency, and all five
modes come within 2 %, with a predicted Q of 4.2 to 4.7.

The three layers then link. Each band's air column, the leaky shells of
layer 2, becomes the shunt branch of the cavity ladder. One fixed network
then gives a single 250 kV ionosphere at DC and Schumann peaks at ELF,
8 % low, with an electric height that comes from the conductivity
profile instead of a fit. The dynamo reaches the air through cosmic-ray
shielding of its conductivity, and every layer keeps its own channel in
one ledger.

## 8.8 What to keep

- A net force needs a momentum ledger as a net energy needs an energy
  ledger. The tether's drag is paid by the orbit, its thrust by the
  supply, and the reaction lands on the field and the ionosphere.
- Declared inputs are declared: the field, the speed and the contactor
  resistances are stated, and the NASA flights show where each is
  optimistic.
- Coiling a wire in a uniform field cancels its motional EMF. Resonant
  multiplication raises voltage by Q; it does not add energy, and the
  ledger says whose energy it is.
- A closed ring stores what its drivers deliver, and synchronized pairs
  add amplitude; a conductor cannot hold two DC levels unless it is
  sectored.
- A vacuum gap has two limits, multipactor and sparking, and a DC bias is
  the measured cure for the first.
- Gravity, rotation and Kron's configuration space are three curvatures.
  Keep them apart, convert by clock rates, and reject a rotor past its
  burst speed.

## Run it yourself

```
go run ./cmd/kronspice examples/tether/export/tether_generator.cir
go run ./cmd/kronspice examples/tether/export/tether_thruster.cir
go run ./examples/tether
go run ./examples/coiled_tether
go run ./examples/ring_towers
go run ./examples/double_hull
go run ./examples/machine_spacetime
go run ./examples/earth_machine
go test ./examples/tether ./examples/coiled_tether ./examples/machine_spacetime ./spacetime
```

The ring-tower and double-hull runs take minutes; their `run_output.txt`
files hold the output quoted above, and `examples/ring_towers/viz/` has
browser animations of the four- and eight-tower patterns.

## Sources

- M. L. Cosmo, E. C. Lorenzini (eds.), *Tethers in Space Handbook*, 3rd
  ed., NASA Marshall Space Flight Center, December 1997, NTRS 19980018321.
- L. Johnson, M. Herrmann, *International Space Station Electrodynamic
  Tether Reboost Study*, NASA/TM-1998-208538, July 1998.
- B. E. Gilchrist et al., "Enhanced electrodynamic tether currents due to
  electron emission from a neutral gas discharge: Results from the TSS-1R
  mission," *Geophys. Res. Lett.* 25, 437–440, 1998; K. Szalai, *TSS-1R
  Mission Failure Investigation Board: Final Report*, NASA, 1996.
- J. Ballance, L. Johnson, "Propulsive Small Expendable Deployer System
  (ProSEDS)," NASA MSFC, NTRS 20010020385; I. Katz et al., *J. Geophys.
  Res.* 100, 1995 (PMG).
- J. R. Sanmartín, M. Martínez-Sánchez, E. Ahedo, "Bare wire anodes for
  electrodynamic tethers," *J. Propulsion and Power* 9, 353–360, 1993.
- N. Tesla, "Coil for Electro-Magnets," US 512,340, 1894; G. Kopp,
  J. Lean, *Geophys. Res. Lett.* 38, L01706, 2011 (solar irradiance).
- C. P. Steinmetz, *Elementary Lectures on Electric Discharges, Waves and
  Impulses*, 1911, Lecture VI; C. R. Paul, *Analysis of Multiconductor
  Transmission Lines*, 2nd ed., 2008, §2 and §4.2.
- A. J. Hatch, H. B. Williams, 1954; J. R. M. Vaughan, *IEEE Trans.
  Electron Devices* 35, 1172, 1988; V. Baglin et al., CERN LHC Project
  Report 472, 2002; W. D. Kilpatrick, *Rev. Sci. Instrum.* 28, 824, 1957.
- C. W. F. Everitt et al., *Phys. Rev. Lett.* 106, 221101, 2011 (Gravity
  Probe B); L. I. Schiff, *Phys. Rev. Lett.* 4, 215, 1960; N. Ashby,
  *Living Rev. Relativ.* 6, 1, 2003; Toray T1000G data sheet, 2024.
- [doc/em-push-pull.md](../doc/em-push-pull.md), [doc/M9n.md](../doc/M9n.md),
  CLAUDE.md §3f and §3l.
