# 10. Hypotheses under test

CLAUDE.md keeps a register of ideas that go beyond established physics.
They are not removed for being unconventional and they are not promoted
without a test; each is stated plainly, modeled where a prediction can
be stated, and run with the standard prediction computed beside it. This
chapter explains how an entry is written and then walks through four of
them, where the engine's energy ledger does most of the work:
[eckardt_homopolar](../examples/eckardt_homopolar) (H15),
[parametric_lcr](../examples/parametric_lcr) (H16),
[pais_cavity](../examples/pais_cavity) (H18) and
[space_power_generator](../examples/space_power_generator) (H13), with
[bearing_motor](../examples/bearing_motor) supplying one bound. Later
sections add the H47 scalar bench, the H48 spin connection as a circuit
term, and two more ledgers (H10, H25).

## 10.1 How an entry is written

Every row of the register has four parts: the hypothesis in one sentence,
with its source; where it is modeled and tested; what would count as
evidence; and a status. The status values are fixed: *untested*,
*bounded* (a null result with the limit recorded), *deviation found*
(with data and conditions) and *reproduced*. An entry leaves the list
only after a test, and even then it stays with its bound written in.

Three rules govern the modeling (CLAUDE.md §1 and §3f):

- **The standard prediction is always computed alongside.** A bench
  result can then be read against both.
- **Hypothetical models are off by default**, labelled in every report,
  and keep their own ledger column, which is never counted as real inflow.
- **A surplus is never credited to the hypothesis.** A ledger residual is
  investigated as an undeclared input, a model error or a numerical error.
  Evolution exploits simulator flaws, and a flaw must never look like
  progress.

The last rule is why the chapter-1 ledger matters so much here. Each of
the four examples below is a published claim of energy from somewhere
other than the declared sources, and in each case the question the
engine asks is the same: when every input is on the books, where did the
energy come from?

## 10.2 Run it: Eckardt's self-excited homopolar (H15)

Eckardt (2012) describes a twin-disk homopolar machine whose loop current
excites its own field coil, and reports an energy surplus at a "giant
resonance" where the effective resistance vanishes. The example first
integrates the paper's equations as printed, then runs the same machine
as a kronspice circuit built on the N machine element. The exported
netlist is the circuit form:

```
Eckardt resonant homopolar
.model hg MACHINE(G_1_4=0.0005 G_2_3=-0.010053 L_3_3=0.001 L_4_4=0.001 N=4 ...)
.domain mechanical w
NH a b b c c d f 0 w 0 hg
R d a 0.3
CJ w 0 0.01 IC=0          ; the shaft inertia J
ITD 0 w DC 1              ; the drive torque, 1 N·m
.tran 0.0015 3 0 0.0015 uic
.end
```

```
go run ./cmd/kronspice examples/eckardt_homopolar/export/eckardt_homopolar_driven.cir
```

```
Node voltages at the end
  v(w)  30.3302 V          ; shaft speed, rad/s

Energy ledger over the interval
  energy in            84.7516 J
  energy dissipated    80.1052 J
  stored: start → end  500 µJ → 4.64692 J
  integration mismatch -27.013 µJ (time-stepping error, not a physical channel)
  worst point residual 0 (relative)
```

The shaft settles at 30.3 rad/s, which is the paper's ω_res = 29.84 rad/s
(its eq. 9), and the books close: 84.75 J of drive work became 80.1 J of
heat and 4.6 J of spin-up and field. The paper's equations as printed
tell a different story:

| | Paper's equations (J·ω̇ = α) | Same machine in kronspice |
|---|---|---|
| Current at 0.30 s | 0.56 A | 0.58 A |
| Current at 0.45 s | 1.2×10⁵ A, past 10⁶ A at 0.463 s | bounded, 9.68 A at 3 s |
| Energy into R and the coil | 1.42×10⁹ J | 80.1 J |
| Drive work | 10.7 J, all of it spin-up | 84.75 J |

The missing 1.423×10⁹ J equals, to one part in a thousand, the work the
omitted back torque ½r²(Bm + B_I)·I would have taken from the shaft. The
Lorentz force that gives the EMF ½ωr²B gives that torque too, so that
ω times torque equals EMF times current exactly; eq. (12) of the paper
leaves it out. In the machine element the back torque is not an added
term but part of the element, and with it the current stays bounded, the
speed settles at the resonance, and the steady state is α = K1·I + K2·I²,
giving 9.95 A and 29.84 W from the drive equal to R·I². The classical
resonance, R_eff = 0, is the self-excitation threshold of a series-wound
dynamo (Bullard 1955; Moffatt 1978), and it is reproduced. The "energy
from spacetime" in these equations is the shaft work they leave out. The
paper's two further resonances are not modeled because no equations are
supplied for them, so the status is *untested*.

## 10.3 The parametric LCR circuit (H16)

Amador and Eckardt (2012) drive a series R, L, C(t) loop with
C(t) = C(1 + a₀ sin(ω₁t + φ)) and report a coefficient of performance of
about 200 against the source's input, read as energy from the time
domain. The exported netlists are the paper's Fig. 4 and Fig. 5:

```
parametric LCR (Amador & Eckardt 2012)
VU s 0 AC 5 0 SIN(0 5 2346.6116 0 0 0)
L s a 0.02 IC=0
R a c 20
C c 0 2.3e-07 SCALE=SIN(1 0.2 4693.2232 0 0 0)     ; C(t) = C(1 + 0.2 sin 2ω₀t)
.tran 4.26e-06 0.03 0 4.26e-06 uic
.end
```

```
go run ./cmd/kronspice examples/parametric_lcr/export/parametric_lcr_fig4_constant_c.cir
go run ./cmd/kronspice examples/parametric_lcr/export/parametric_lcr_fig5_pumped_c.cir
```

| | Fig. 4, constant C | Fig. 5, C pumped at 2ω₀ |
|---|---|---|
| Energy in | 17.52 mJ | 286.6 J |
| Energy out, returned to sources | 0 | 34.1 J |
| Energy dissipated | 16.89 mJ | 162.8 J |
| Stored at 30 ms | 625 µJ | 89.7 J |
| Integration mismatch | 617 nJ | 12 mJ |
| Worst point residual | 0 | 0 |

With a constant capacitor the circuit does what a 5 V source into 20 Ω
at resonance does: 0.25 A peak, and the paper's COP of 0.964. With the
capacitor pumped at twice resonance the current is 87 A after 30 ms and
still growing, and the CLI's ledger shows 287 J entering a circuit whose
voltage source is 5 V. The source is not where it comes from. The
engine's time-varying capacitor is q = s(t)·C·v, the paper's eq. (25)
exactly, and the agent that moves the plates does work on the circuit at
the rate −ṡ·½Cv²: that is the pump of Mandelstam and Papalexi (1932),
and it is on the ledger as inflow because it is a declared mechanical
input. The example reads it out separately:

| Paper case | Current | Source | Pump | Heat in R + stored |
|---|---|---|---|---|
| Fig. 4, constant C | 0.25 A | 0.0175 J | 0 | 0.0169 + 0.0006 J |
| Fig. 5, pumped at 2ω₀, at f₀ | 87.3 A, unbounded | −0.085 J | 256 J | 165 + 91 J |
| Figs. 6–7, pump and drive at 0.97 f₀, φ = −90° | 16.1 A, bounded | 6.75 J | 1139 J | 1143 + 2 J |
| Fig. 8, a₀ = 0.2 at the drive frequency | 0.123 A, no growth | 1.8 mJ | 0.1 mJ | 1.9 mJ |

The books close to 0.04 % to 0.4 % in every bounded case, and the a₀ =
0.8 square-wave runs close to a few per cent, an integration error that
falls as the step is refined. The paper's COP divides the heat in R by the
source energy alone: 169 for the φ = −90° case it highlights. Divided by
source plus pump it is 0.998. The standard theory is Kouřil and Vrba's:
a₀ = 0.2 gives an elastance depth of 0.202 against a self-excitation
threshold of 2/Q = 0.136, so pumping at 2ω₀ self-excites, with the
envelope growing at 245 1/s as they predict; detuned to 0.97 f₀ the
threshold rises above the depth and the circuit is a parametric amplifier
of the drive. Every case in the paper is reproduced and the status is
*untested*, because the decisive measurement has not been made: the power
that drives C(t), at the motor shaft or the varactor's bias supply.

## 10.4 The Pais cavity (H18)

Pais's patent (US 10,144,532 B2, 2018) claims an annular microwave cavity
between two hull walls, charged, vibrated and rotated, that reduces the
craft's inertial mass. The patent gives no gap, height, voltage or drive
power, so the geometry is declared: radius 2 m from its example, a 5 cm
gap, 0.5 m walls, copper, 640 line sections, the n = 16 ring mode at
381.3 MHz inside its stated band, and 1 MW of drive. The exported netlist
is the ring as a closed parallel-plate line, 1,925 elements:

```
pais annular cavity ring
R0 n0 m0 0.00040011
L0 m0 n1 2.4674e-09
C0 n0 0 1.7385e-12
R1 n1 m1 0.00040011
...
```

Its results sheet is the one place in this chapter where Kron's tearing
of chapter 3 carries the load: Happ's two-level diakoptics tears the ring
into 16 sectors hinged on the inner wall and joins them by 16 links, so
the coupled solve has order 16 against 1,280 MNA unknowns. The torn and
full solutions differ by 1.1×10⁻⁸ at the mode, with the power taken by
the branches, 177,358 W, identical in both, as Roth's power invariance
requires.

| Case at 1 MW | Q | Stored | V_max | Energy in the emitter's quadrant |
|---|---|---|---|---|
| A. four emitters in phase | 1.48×10⁴ | 6.17 J | 148.9 kV | 0.250 |
| B. one emitter, on resonance | 1.48×10⁴ | 0.385 J for 62.5 kW | 37.2 kV | 0.250 |
| C. one emitter, between modes | 1.5×10⁴ | 1.1 µJ for 0.18 W | 63 V | 0.250 |
| D. one emitter, walls 2000× lossier | 7.4 | 3.3×10⁻⁷ J for 107 W | 64 V | 0.817 |

Input power equals wall loss to 2.5×10⁻¹⁰. The patent's sectional
steering does not work on a high-Q ring: a wave decays over about 3.7 km
and the ring is 12.6 m around, so energy fed into one sector fills the
ring evenly, and only a ring of single-digit Q keeps energy near its
emitter. Then the patent's figures against physical ceilings:

| Figure in the patent | Implies | Limit |
|---|---|---|
| σ = 1 C/m² | E = 1.13×10¹¹ V/m | Kilpatrick spark limit 1.9×10⁷ V/m; field emission drains it in 4×10⁻¹⁶ s; wall pressure 5.7×10¹⁰ Pa against copper's 2.1×10⁸ Pa |
| 10³³ W/m² from wall vibration | wall speed 8.9×10²¹ m/s | 3×10¹³ c |
| wall vibration at 10⁹ to 10¹⁸ Hz | | lattice vibrations end at copper's Debye frequency, 7.1×10¹² Hz |
| 30,000 rpm at 2 m | rim 6,283 m/s | 3.3 times a carbon-fibre hoop's burst speed |
| thrust | | zero net time-averaged force on a closed cavity; beaming all 1 MW out gives P/c = 3.3 mN |

The standard prediction is no net force, no mass change, and at most P/c
of thrust from radiation that escapes. The gyroscope weight loss the
patent cites (Hayasaka and Takeuchi 1989) was not reproduced by three
groups. The bench test is to weigh a sealed, driven, charged cavity
against a dummy load dissipating the same power. Status: *untested*.

## 10.5 Tewari's Space Power Generator (H13)

Tewari's book reports one loaded test of his homopolar generator and
reads it as an efficiency above one. The example runs the published
numbers through the ledger.

```
go run ./examples/space_power_generator
```

```
1. The book's arithmetic
   W_r = 8620 − 2943 = 5677 W (6.20); W_E = 3.41 V × 3438 A = 11724 W (6.21)
   incremental efficiency 2.065 (6.22); ÷ η_DM 0.83 = 2.49 (6.23)
2. All inputs
   DM 8620 W + excitation 650 W = 9270 W in; claimed electrical 11724 W → ratio 1.26
3. The assumption: E under load = E at no load
   extra shaft power into the SPG under load: 0.83 × 5677 W = 4712 W → 16.7 N·m
   standard electrodynamics: E_load·I = shaft power → E_load ≤ 1.37 V
```

The claim rests on multiplying the load current by the EMF measured at
no load. Standard electrodynamics requires the mechanical power the rotor
takes from the shaft to equal the loaded EMF times the current, and the
book's own numbers give that shaft power: 4712 W, a torque of 16.7 N·m at
282.7 rad/s. So the loaded EMF can be at most 1.37 V, a fall of 60 % from
no load at 3,438 A, which armature reaction, saturation and speed sag at
the rotor would produce; the claim needs it unchanged at 3.41 V and a
back torque of 41.5 N·m that nothing delivered. The book reports no
loaded-EMF measurement.

[bearing_motor](../examples/bearing_motor) asks whether the missing
24.8 N·m could come from the machine's own ball bearings acting as Huber
motors, with the load current through them. Its registry holds only
measured values, and the answer is a bound either way. In the measured
mode, where a back-EMF exists, the bearings would drop 2.04 V of the
3.41 V and the rest of the loop would get 4.71 kW, not 11.72 kW: the
1.37 V left is the same bound again, because a bearing motor moves shaft
work around and adds none. In the thermal mode, at the measured
efficiency of 0.10, the bearings would need 70 kW and 20.4 V across them,
more than the machine makes. The decisive measurements are the voltage
drops around the loaded loop, or a torque transducer on the shaft, or
insulated bearings. Status: *untested*; the claim rests on an unmeasured
loaded EMF.

## 10.6 Russell's octave bar magnet (H33)

Walter Russell painted the bar magnet as two spheres of lens-shaped rings
strung on its axis, and wrote beneath it that "a bar magnet has four
pairs of poles — not one pair".
[octave_magnet](../examples/octave_magnet) builds the standard field
exactly. A magnetized rod is a sheet of current on its side
(J = ∇×M), its poles are its magnetic charge −∇·M, and its field lines
are the contours of the flux function ρA_φ. The sheet matches the exact
solenoid formula to 1.3×10⁻⁶, and a search coil modelled in kronspice
reads the field model's flux to six digits.

The two readings separate cleanly. A uniform bar carries charge only on
its two end faces, so its outside lines run end to end in one family.
Four pairs of poles need the magnetization to reverse three times inside
the bar; its lines then bow from pole to pole along the axis, the chain
of lenses Russell drew. The bench test is a Hall probe run beside a real
bar: B_z changes sign twice for one pair, and five or more times for
four. Until that scan is made the entry stays *untested*.

The same example closes the bar into a ring of twelve sections, each
magnetized opposite to its neighbours: either radially (N S N S … round the
outside) or along the ring, so that like poles meet at every joint as in
the four-segment bar. Both give twelve poles round the outside, at the
section centres or at the joints, and a crown of twelve field arches
between them. A ring magnetized outward all round gives none. A pickup
coil modelled in kronspice sees six cycles per revolution, and its ledger
closes.

Read as hydrogen, Russell's words become entry H34,
[hydrogen_poles](../examples/hydrogen_poles). Poles come in powers of two:
one pair is a dipole, four pairs an octupole. An atom has no magnetic
quadrupole, and a 2^k-pole needs angular momentum j ≥ k/2, so ground-state
hydrogen has exactly one pair. The computed field confirms it, with the
dipole equal to the Landé value and the octupoles in the Wigner–Eckart
ratio. Four pairs first appear in the excited 2p3⁄2 level, whose spin
magnetization makes four N S N S zones round the atom out to 15.6 Bohr
radii. Read as two atoms, the painting's spheres are H₂ with its four small
magnets. The Heitler–London bond forms only with the electron spins
opposite, so those two magnets cancel, and what is left are the protons'
dipoles, 658 times weaker.

The same painting's "high pressure, short waves" and its two spheres
meeting at a plane hold in [light_hydrogen](../examples/light_hydrogen)
as the confined atom and the 2p_z state (chapter 9, section 9.7). The
painting and all three comparisons are in
[doc/russell-octave-report.md](../doc/russell-octave-report.md).

## 10.7 Konstapel's toroidal coherence propulsion (H35)

J. Konstapel (2025) proposes that heated fiber-optic toroidal coils,
phase-locked to the Earth's field, invert a scalar "Sn field" and make a
craft lighter. The paper gives no equation with a value for its coupling,
so [toroidal_coherence](../examples/toroidal_coherence) computes the
standard prediction for each of its three Tier 1 measurements, and the
ordinary effects a bench must remove first.

The predicted (3n ± 1) harmonics are what any balanced three-coil array
radiates. With three coils at 120° in kronspice, each closed through an
ordinary rectifying contact, multiples of 3 cancel to 10⁻¹⁰ by symmetry,
3n + 1 rotates forward and 3n − 1 backward, and the ledger closes. The
discrete infrared lines are a closed loop's longitudinal modes,
c/(n_g·ℓ) = 20.5 MHz apart for 10 m of silica. On the balance, ordinary
effects reach the claimed band of 500 mg to 50 g: air buoyancy lost round
the hot coil gives 16.7 g, and a feed lead in the geomagnetic field
1.3 g per metre. The paper's own figures fail on its arithmetic: 1 kN on
1000 t gives 0.001 m/s², not 1 m/s². Its fifth part, Hall's reading of
lensed galaxy images as total reflection at hydrogen clouds, turns light
by at most 10⁻³ arcsec in a dense core, against 1–30 arcsec observed.
That reflection would also be chromatic, so image positions that agree
across bands rule it out. Status: *untested*. The decisive bench is a
vacuum balance with a dummy coil at the same power.

## 10.8 Barrett's oscillator-shuttle-circuits (H36)

T. W. Barrett (1991; US 5,493,691) reads Tesla's coils, with energy
shuttled among floating grounds, as nonlinear circuits beyond Kirchhoff
that mix pump and signal into an idler, phase-conjugate it and carry
quaternion (SU(2)) algebra. [barrett_osc](../examples/barrett_osc)
builds his Figs. 2–5 from the drawings with declared values. Every
circuit of Fig. 2 is linear: transient and AC agree, and no mixing
product rises above round-off, however the grounds float. Replace a
top-load with a varactor and the idler appears at fp − fs, obeying
Manley–Rowe and phase-conjugating the signal, with the standard
Es*·Ep scaling rather than the paper's |Es|²·Ep. The mixing comes from
the nonlinear element, not the topology. Fig. 3 is a working frozen-wave
pulse former (−6.7, +3.4, −1.4 kV lobes; ledger to 10⁻¹⁰ J). Fig. 4's
"αβ ≠ βα" is two-port cascade algebra: ABCD matrices do not commute, a
lossless S matrix is unitary, and cascades are associative, so they
cannot carry the Cayley numbers Barrett also names. Status: *untested*
for the claims beyond these (noise rejection, lower transmission loss),
which the paper gives no measurement or equation for. Bench test: look
for fp − fs in a Fig. 2 build with air-core coils.

## 10.9 The magnetostrictive transformer with two transducers (H37)

The proposal follows Barrett: give the transformer a core whose strain is
coupled to its flux (Terfenol-D or Galfenol, as a rod, a ring or a helical
spring), put an ultrasonic transducer at each end, and look for the drive
and the sound mixing, with the swap test (drive the secondary, watch the
primary) as the sign of the nonlinearity.
[magnetostrictive_transformer](../examples/magnetostrictive_transformer)
builds it from one energy potential for B and stress, with declared Etrema
material values, and computes the standard prediction. The core does mix:
products of drive and sound (fa ± fe, fb ± fe) reach 2.9 mA on the
secondary of the Terfenol ring and vanish without sound, while sound ×
sound is 100× weaker. At small drive the transformer is reciprocal to
2×10⁻¹⁰; at full drive the products are not (×0.56 for the Galfenol ring,
×0.75 for the rod when swapped). In a ring only the breathing mode reaches
uniform windings, and the helix's spring shape does not amplify the
strain. Status: *bounded*: the mixing is ordinary magnetoelastic
nonlinearity, of a computed size.

## 10.10 The aluminium-cell ground pump (H38)

The negative ground pump of photon_electron is rebuilt with
aluminium–water valve cells (the Pollak and Nodon electrolytic rectifier
of sparkgap_diode) as its rectifiers, and a genetic search designs it for
"catalytic efficiency": pumping electrons from ground at an energy per
electron as close to e·|V| as possible. [al_cell_pump](../examples/al_cell_pump)
treats e·|V| as a floor that no circuit goes under, and the ledger checks
it. The baseline loses 98–99 % of its drive in the tank coil (η 0.0006
with silicon diodes, 0.00002 with aluminium cells). Evolution reaches
η 0.857 at −321 V (1.66×10¹⁶ electrons per joule against the ideal
1.94×10¹⁶), and a single aluminium rod in an iron cup reaches 0.869 at
−197 V. The oxide film's capacitance is the cell's main cost, silicon is
still better by 0.01–0.09, and the oxide is too linear (κ = 2 %) for
parametric pumping. Status: *bounded*: nothing goes under e·|V|.

## 10.11 The Rodin coil in quadrature (H43)

The Rodin coil winds its wire over a torus and through the hole, advancing
150° per pass between 12 poles: 1→6→11→4→9→2→7→12→5→10→3→8→1, a {12/5}
star from above and a torus-knot path of the (5, 12) class, five times
round the axis and twelve times through the hole. Claims of a "vortex"
field that accumulates energy, and of over-unity, come only from
promotional sources. [rodin_resonance](../examples/rodin_resonance) builds
the coil characterized by Ramírez Galindo et al. (*Rev. Mex. Fís.* 70,
030901, 2024) as a counter-wound pair; one fitted turn count (12 circuits
per winding) reproduces their centre field at 1.60 and 1.25 A within their
own spread. A sine on one winding and a cosine on the other turn the field
vector at each point (ellipticity about 0.1), but the field pattern stands
still round the axis; the mutual inductance passes 0.44–0.45 W from one source
to the other, ω·M·I_A·I_B·sin φ to four figures. Energy accumulates only by
resonance, to the Q-limited level (Q 76.4 against ω0L/R = 77.3), and a
resonant pair locks its currents in phase. Near 2 MHz the winding is a line
laid on a knot; closed in series on itself it is photon_electron's crossed
ring (H23), with 720° modes for the (5, 12) and the (2, 1) path alike.
Output coils evolved together with the drive frequency become a
transformer secondary: one coaxial coil at 4.1 kHz, or two coils under the
aperture at 29 kHz, charge a capacitor bank through Schottky bridges with
18–32 % of the source energy, and the ledger closes to the integration
error. Status: *bounded*: a coupled-inductor pair, output ≤ input.

## 10.12 Nachamkin's force-free plasmoids: harmonics and size (H26)

Nachamkin's 1992 report (PL-TR-92-3044) proposes ball lightning as a
plasma sphere holding a standing wave whose E and H are parallel
everywhere. [plasmoid_rings](../examples/plasmoid_rings) checked his field
and his worked example; [Force-Free_Harmonic_Plasmoids](../examples/Force-Free_Harmonic_Plasmoids)
studies the paper as written, for every size and energy reported for ball
lightning (diameter 10–50 cm, mean 24 cm; lifetime 2–5 s, mean about
10 s).

His resonance condition, j1(αR) = j1(kR) = 0 with k > α, allows each size
a comb of frequencies f_n = x_n·c/2πR, and his charge condition (eq. 49)
is the cold-plasma dispersion relation; a kronspice plasma-loaded line
resonates at it to the lattice error. His Appendix B has an arithmetic
slip: E0 = 3.566×10⁹ V/m follows from ∫|M|²dV = 35.87/α², while his own
B3, confirmed numerically, is 35.87/α³, so the field should be √α = 6.7
times larger. With his recipe the electrons move at v/c = ξ·kR, between
1.4 and 12 c over the whole reported range. Keeping them below 0.1 c moves
the plasmoid into the terahertz.

Free, it radiates. Its two halves (the TE and TM l = 1 modes, whose
spectrum kronspice returns to 10⁻¹⁵ in frequency and 2×10⁻⁶ in Q) lose
their energy in τ_TE = R(kR)²/(2c·x1²), a closed form checked against the
root finder. Lasting the seconds that ball lightning lasts sets kR,
independent of the energy, and with it the electron density: 7–74 times
sea-level molecular density for a 12 cm radius, where fully stripped air
holds only 14.5. Status unchanged: *untested* as ball lightning, its
equations sound, its example and its free lifetime bounded.

The standard counterpart stores the energy in chemistry instead.
[silicon_ball_lightning](../examples/silicon_ball_lightning) tests
Abrahamson and Dinniss's model, in which lightning reduces soil silica to
silicon nanoparticles that burn slowly in air, against the first spectrum
of a natural ball (Cen, Yuan and Xue 2014: silicon, iron and calcium
lines for 1.64 s). Thirty-one milligrams of silicon hold the 1 kJ of an
average ball; Deal-Grove kinetics burn 10 nm particles in seconds only at
about 2000-3000 K, which the laboratory estimate of 3140 K supports. A
welding arc at the laboratory operating point (kronspice: 120 A, 23 V,
2.5 kW) makes enough vapor for seconds of glow, and the 35 kV line beside
the Qinghai ball can flicker it at 100 Hz but delivers microwatts, not the
100 W it shines.

## 10.13 The scalar S on a rotating disc (H47)

Section 5.5 gives the field equations with Maxwell's scalar part kept.
H47 goes one step further: it proposes that S is physical and couples to
gravity, so a circuit that makes S would change the weight of nearby
masses. In Minotti and Modanese's theory S enters through a term β_AB S²
beside two scalar fields of the Brans–Dicke kind. One parameter,
C_φ ≈ 4–9×10²¹, fits the Podkletnov and Poher reports. Those reports have
not been replicated, and the authors call the fit phenomenological. S
has a source only where local charge conservation fails, so a bench has
to make such a place. Iadicicco and Verolino propose one: a fibre ring
on a spinning aluminium disc lights one sector and induces a charge
there, and the lit sector turns with the disc. In the rotating frame the
current splits so that I_G = Ω∂_φρ, and □S = μ₀I_G.

[scalar_disc](../examples/scalar_disc) builds that bench with declared
values and runs the scalar on it:

| Quantity | Value | Standard prediction |
|---|---|---|
| Peak S at the sensor | 1.10×10⁻¹² T (cS = 0.33 mV/m) | 0 |
| Reversed rotation | S → −S exactly | 0 |
| Whole ring lit | S = 0 exactly | 0 |
| S² force per mass (MM fit) | 5.4×10⁻³⁷ m/s², the same for both senses of rotation | 0 |
| At the bench's limits (air breakdown 1.5 nC; Al burst 100,000 rpm) | S 5.5×10⁻¹¹ T; force 1.3×10⁻³³ m/s² | 0 |
| Ordinary systematic: Coulomb field of the gated charge | 2.1×10⁴ V/m; 2.1×10⁻⁵ m/s² on a 1 g sensor carrying 1 pC | the same |

The table shows three things.

- **The two papers disagree on symmetry, and the bench can decide.** S is
  exactly odd in Ω. Minotti and Modanese's force goes as S², so it is
  even in Ω. Iadicicco, Modanese and Verolino's linear response is odd.
  A run clockwise and then anticlockwise tells them apart.
- **The predicted force is far too small to see.** The junction S that
  fits the gravity reports is 4×10¹³ times the disc's S, and the ordinary
  Coulomb force of the gated charge is 10²⁸ times the scalar force.
- **A null result still bounds the theory.** A null at 10⁻⁹ m/s² bounds
  C_φ < 5×10³³. The fit's ω₀ ≈ 3.5 already conflicts with Cassini's
  Brans–Dicke bound ω > 4×10⁴.

Status: untested. The bench measures Iadicicco, Modanese and Verolino's
Λ_opt, which has no published value. Step 3 of
[doc/h47-scalar-gravity.md](../doc/h47-scalar-gravity.md), the junction
numbers, is next.

## 10.14 The ECE spin connection in a circuit (H48)

Section 5.5 puts ECE's spin connection into the field equations. To run
it in a circuit, the field relation has to be integrated along a
component. This mapping is proposed here; it is not in the papers.
Integrating the ω₀ term along a coil's wire gives a torsion voltage.
Integrating the ω term along a lead gives a scalar factor on each
terminal potential:

```
torsion (ω₀A):   v = L di/dt + s·cω₀L·i        coupled coils: + s·cω₀M·i_other
scalar  (ωφ):    the body sees (1 + s·κ_a)φ_a − (1 + s·κ_b)φ_b,   κ = ∫ω·dl
```

The factor s is +1 for ECE (eq. 4.211) and −2 for ECE2 (eq. 6.189). So
ECE2 reverses every ECE effect and doubles it. UFT 456's field-built
currents become a shunt 2cκ₀·C across every capacitor (`k0=`) and a
series −2cκ_Λ0·L in every coil (`kl0=`). The series term is negative for
κ_Λ0 > 0, which makes it a gain. Ideal R and C have no series
inductance, so they carry only the scalar term.

The terms are declared with one line, `.ece form=ece|ece2 w0=… kappa=…
k0=… kl0=…`. Every element can override the defaults (`W0=`, `KAPPA=`,
`KAPPA_A=`, `KAPPA_B=`). The power the terms deliver is booked per
element as model-supplied and labelled hypothetical. Zero terms are
Maxwell bit for bit. Kron's primitive tensors refuse an `.ece` circuit,
so a torn or mesh solve cannot hide the terms. A 1 V source drives a
coil through 10 Ω, with ω₀ = 10⁻⁶ m⁻¹:

```
ece2 series RL
.ece form=ece2 w0=1e-6
V1 in 0 DC 1
R1 in a 10
L1 a 0 10m
.op
.end
```

```
  V1    voltage source  in→0   1 V         -249.741 mA  -249.741 mW
  R1    resistor        in→a   2.49741 V   249.741 mA   623.705 mW
  L1    inductor        a→0    -1.49741 V  249.741 mA   -373.964 mW
  ...
  model-supplied 373.964 mW (idealized active elements; not inflow from the exterior)
    L1  373.964 mW  [ECE2 spin connection (hypothetical, H48): w0=1e-06 1/m]
  residual      1.11022e-16 W (relative 0)
```

Under Maxwell the current is 100 mA. With ECE2 the torsion term is
s·cω₀L = −2 × 2.998×10⁸ × 10⁻⁶ × 0.01 = −5.996 Ω in series with the
coil. The current rises to 249.7 mA, and the coil delivers 374 mW at DC.
The resistor heats by 624 mW against 250 mW from the source. The ledger
closes only because the coil's 374 mW is booked as model-supplied, the
way a negative resistance is (chapter 1). Change the form to `ece` and
the same ω₀ gives +2.998 Ω. The current falls to 76.9 mA and the coil
absorbs 17.7 mW. The tests (`analysis/ece_test.go`, `kron/ece_test.go`)
check several things:

- Each case equals its ordinary-element twin, a resistor of s·cω₀L in
  series with the coil, to 10⁻⁹ in OP and AC.
- At resonance a series RLC carries 1/(R + r_t).
- The model-supplied power is −r_t|I|²/2.
- The ledger closes in OP, AC and transient.

[ece_rlc](../examples/ece_rlc) runs the two terms one at a time on a
series RLC: R 1 Ω, L 10 mH, f0 173 Hz, Q 10.87. Each run is set beside
its ordinary-element twin, in AC and in transient. The torsion term is
set to cω₀L = 0.25 Ω. It leaves f0 where it is and changes Q and the
current. With ECE, Q is 8.70, the peak current 0.80 A, and the term
absorbs 0.08 W. With ECE2, Q is 21.74, the peak 2.00 A, and the term
supplies 1.00 W while the source gives 1.00 W and the resistor heats
2.00 W. The scalar term is set to κ = 0.05 on the capacitor's top
terminal. It moves the resonance to 168.84 Hz (ECE) or 182.37 Hz (ECE2),
as the closed form predicts, and it supplies almost no power on average.
In transient the books close to a few 10⁻⁶ J for every case. The
netlists are in `examples/ece_rlc/export/`.

What this shows is the size of a claim, not a source of energy. Any ω₀
that makes a measurable gain is a negative resistance whose energy has
no declared source. The standard prediction is that ω₀ is a frame
quantity, so a coil has no torsion resistance. A bench test is the I–V
line of an air-core coil at DC: its slope is R_wire, and the torsion
term would add s·cω₀L.

[examples/AIAS](../examples/AIAS) is the running evaluation of the AIAS
UFT papers. It checks each paper's equations and, where it can,
simulates the paper's device with the ledger on. The circuit terms
above came out of UFT 454–456. There, `form=ece` equals UFT 455's
mapping case by case to 10⁻⁹. UFT 456's Table 2 speeds turn out to be
level points of exponential profiles, and the front is at c in every
medium: the Yee solver, the 3-D field solver and the kronspice ladder
all agree. The folder's `INDEX.md` gives each paper's verdict in a
fixed vocabulary. "Derivation holds" means the algebra checks.
"Derivation fails at eq. n" names the first equation that does not
follow from the ones before it. "Assumption added" marks a step that
needs something the paper does not state.

## 10.15 Two more ledgers: the MEG and the Casimir battery (H10, H25)

Two more entries are claims of energy from a source that is never
declared.

[bearden_meg](../examples/bearden_meg) builds the Motionless
Electromagnetic Generator (US 6,362,718) as a magnetic circuit with the
full ledger. Beside it sits a bucking-magnet "knife edge", two magnets N
to N with a split drive coil at the null. In 60 runs the MEG's COP
ranges from 0.002 to 0.16 and the knife edge's from 0.06 to 0.990. The
knife edge is a good transformer through all-iron paths. The COP is flat
in drive level, spacing and mismatch. Over whole cycles the magnets
supply no net power, at most 10⁻³ of the drive, and that goes into
stored energy while the circuit settles. Ngspice twins of both circuits
are included. The bench test still needed is a replica with true power
measured over whole cycles.

[forward_casimir](../examples/forward_casimir) runs Forward's 1984
vacuum-fluctuation battery. The leaves are charged, and the Casimir pull
closes them against a servo. One leaf pair is a Kron (q, a) system with
its ledger closed, and the Lifshitz energy for real aluminium is 0.23 of
the ideal at 20 nm. Each stroke releases 12 µJ/m², which is 0.02 J/kg.
This is the plates' surface energy, released once, as Cole and Puthoff
showed. Reopening the leaves costs more than closing them gave, so a
close–reopen cycle is a net loss, and the battery is not a continuous
source.

## 10.16 The bounds recorded so far

| Entry | What was computed | Standard reading | Status |
|---|---|---|---|
| H13, SPG | loaded EMF ≤ 1.37 V against 3.41 V assumed; bearings cannot close a 24.8 N·m gap | back torque k·I; Poynting | untested, unmeasured EMF |
| H13, electron | g = 15/8 against 2.0023; proton charge 12e; α = 15π²/16 | Dirac + QED | untested, Zeeman prediction needed |
| H15 | surplus = omitted back-torque work to 10⁻³; with it, 9.9 A and 30 W from the drive | series-dynamo self-excitation | untested, two resonances unmodeled |
| H16 | every case reproduced; COP 169 by source alone, 0.998 with the pump | Mandelstam–Papalexi; Kouřil & Vrba §10.2 | untested, pump power unmeasured |
| H18 | no steering on a high-Q ring; σ, flux, speed and vibration figures beyond material limits | no force on a closed cavity; thrust ≤ P/c | untested |
| H21 (chapter 9) | g ≤ 1 for a mass ring; massless ring needs 43.6 c | Dirac + QED; ESR 28 GHz/T | bounded |
| H23 (chapter 9) | light leaves the loop in under half a period | linear Maxwell | bounded |
| H22, [quantumtricity](../examples/quantumtricity) | circuit theory is U(1) on Kron's network (KVL zero holonomy, KCL the Gauss generator, Yee grid = FDTD to 7e-13); SU(3) loops on the same graph compose with defect 0.75, every circuit quantity with 3e-16 | circuit electricity is the abelian classical limit; V, I, R quantized by e and h | bounded |
| H35, [toroidal_coherence](../examples/toroidal_coherence) | 3n ± 1 harmonics from any balanced three-coil array; buoyancy 16.7 g and lead forces in the claimed 0.5–50 g band; Hall's reflection ≤ 10⁻³ arcsec against 1–30 observed | three-phase symmetry; loop modes c/(n_g·ℓ); E/c² makes the hot coil heavier | untested |
| H36, [barrett_osc](../examples/barrett_osc) | Figs. 2–5 linear (no idler); a varactor gives the idler, Manley–Rowe and phase conjugation; order effects are ABCD cascades | parametric mixing; two-port matrix algebra | untested |
| H37, [magnetostrictive_transformer](../examples/magnetostrictive_transformer) | drive × sound products up to 2.9 mA, sound × sound 100× weaker; reciprocal to 2e-10 at small drive, not at full drive | magnetoelastic nonlinearity of the core | bounded |
| H38, [al_cell_pump](../examples/al_cell_pump) | evolved η 0.857 at −321 V (0.869 for one rod in an iron cup); coil and oxide film the losses; Si still better | e·|V| per electron is the floor | bounded |
| H43, [rodin_resonance](../examples/rodin_resonance) | sin/cos drive turns the field vector, the pattern stands still; 0.44–0.45 W passed between the sources through M; Q-limited build-up; output coils charge a bank with 18–32 % of the input | coupled inductors / transformer; near 2 MHz a line on a (5, 12) knot | bounded |
| H10, [bearden_meg](../examples/bearden_meg) | MEG COP 0.002–0.16, knife edge 0.06–0.990; magnets supply ≤ 10⁻³ of the drive over whole cycles | transformer with permanent-magnet bias | untested, bench needed |
| H25, [forward_casimir](../examples/forward_casimir) | 12 µJ/m² per Al stroke (Lifshitz 0.23 of ideal); close–reopen is a net loss | one-shot surface energy (Cole & Puthoff 1993) | one-shot, not continuous |
| H47, [scalar_disc](../examples/scalar_disc) | S 1.1×10⁻¹² T, odd in Ω; S² force 10⁻³⁷–10⁻³³ m/s²; a 10⁻⁹ null bounds C_φ < 5×10³³ | S = 0 under charge conservation; GR weight | untested |
| H48, `.ece` and [AIAS](../examples/AIAS) | torsion s·cω₀L and k0, kl0 equal ordinary resistors to 10⁻⁹; their power model-supplied; fronts at c | ω a frame quantity; a gain needs a declared source | untested |
| H26, [Force-Free_Harmonic_Plasmoids](../examples/Force-Free_Harmonic_Plasmoids) | comb f_n = x_n·c/2πR; electrons at 1.4–12 c by his recipe; free lifetime τ_TE = R(kR)²/2c·x1²; seconds need 7–74 × sea-level density at 12 cm | cold-plasma standing wave; leaky sphere modes | untested; example and lifetime bounded |

In every case the engine found the energy, and in every case it was
already on the books: shaft work in H15, the pump in H16, the loaded EMF
in H13. H22 is the one entry about structure rather than energy: the
quark model and circuit theory are both gauge theories on a graph, and
the test found the difference between them to be exactly the difference
between SU(3) and U(1), with no circuit quantity on the non-abelian side. None of these is a verdict on the idea behind the entry. Each is
a statement of what the published equations or numbers contain, with the
measurement that would change the status written next to it.

## 10.17 What to keep

- A register entry is a model, a standard twin, an evidence criterion and
  a status, and it stays after a null result with its bound.
- The ledger is the instrument. A claimed surplus is traced to the input
  the claim left out before anything else is considered.
- Hypothetical physics runs only when switched on, is labelled, and never
  reaches the real ledger columns.
- The output of a test is the next measurement, stated so that a bench
  can make it.

## Run it yourself

```
go run ./examples/Force-Free_Harmonic_Plasmoids   # about five seconds
go run ./examples/eckardt_homopolar
go run ./examples/parametric_lcr
go run ./examples/space_power_generator
go run ./examples/bearing_motor
go run ./examples/pais_cavity
go run ./examples/quantumtricity
go run ./examples/toroidal_coherence
go run ./examples/barrett_osc
go run ./examples/magnetostrictive_transformer
go run ./examples/al_cell_pump -study pipes
go run ./examples/bearden_meg
go run ./examples/forward_casimir
go run ./examples/scalar_disc
go run ./examples/ece_rlc
go run ./cmd/kronspice examples/ece_rlc/export/ece_rlc_ece2_torsion.cir
go test ./analysis ./kron -run ECE                 # the .ece terms against their twins
go run ./cmd/kronspice examples/eckardt_homopolar/export/eckardt_homopolar_driven.cir
go run ./cmd/kronspice examples/parametric_lcr/export/parametric_lcr_fig5_pumped_c.cir
```

## Sources

- H. Eckardt, "Simulation of resonances in the homopolar generator,"
  A.I.A.S./UPITEC, 13 February 2012. E. C. Bullard, *Proc. Camb. Phil.
  Soc.* 51, 744 (1955); H. K. Moffatt, *Magnetic Field Generation in
  Electrically Conducting Fluids*, 1978, ch. 6.
- F. D. Amador, H. Eckardt, "Simulation of a Parametric Resonance
  Circuit," A.I.A.S./UPITEC, 30 March 2012. L. Mandelstam, N. Papalexi,
  *Z. Phys.* 73, 223 (1932); F. Kouřil, K. Vrba, *Non-linear and
  Parametric Circuits*, 1988, §10.2; J. M. Manley, H. E. Rowe, *Proc. IRE*
  44, 904 (1956).
- S. C. Pais, "Craft using an inertial mass reduction device," U.S.
  Patent 10,144,532 B2 (2018). H. H. Happ, *Diakoptics and Networks*,
  1971, Ch. XII §4. W. D. Kilpatrick, *Rev. Sci. Instrum.* 28, 824 (1957);
  R. G. Forbes, *J. Vac. Sci. Technol. B* 26, 788 (2008). H. Hayasaka,
  S. Takeuchi, *Phys. Rev. Lett.* 63, 2701 (1989); J. E. Faller et al.,
  *Phys. Rev. Lett.* 64, 825 (1990); T. J. Quinn, A. Picard, *Nature* 343,
  732 (1990); J. M. Nitschke, P. A. Wilmarth, *Phys. Rev. Lett.* 64, 2115
  (1990).
- P. Tewari, *Universal Principles of Space and Matter*, 2017, §6.7,
  pp. 149–161. The bearing-motor sources in the bearing_motor README:
  Huber 1959, Milroy 1967, Gruenberg 1978, Moyssides and
  Hatzikonstantinou 1990 and 1997, Watson 1993 and 1999, Choo, Soong and
  Abbott 2005, Marinov 1989, McDonald 2020.
- J. P. Roth, *Quart. Appl. Math.* 17 (1959), power invariance under
  tearing.
- J. Konstapel, "Toroidal Coherence Engineering: From Hall's Photon
  Theory to Practical Propulsion Systems," Leiden, 1 November 2025.
  P. Touboul et al., *Phys. Rev. Lett.* 129, 121102 (2022); E. R. Peck,
  S. Huang, *J. Opt. Soc. Am.* 67, 1550 (1977).
- T. W. Barrett, *Ann. Fond. Louis de Broglie* 16, 23–41 (1991); U.S.
  Patent 5,493,691 (1996).
- Etrema Products, Terfenol-D and Galfenol data sheets (material values
  declared in [magnetostrictive_transformer](../examples/magnetostrictive_transformer)).
- *Encyclopaedia Britannica*, 11th ed. (1911), "Transformers" (the Pollak
  and Nodon aluminium rectifier); A. Güntherschulze and H. Betz, the
  high-field ionic conduction law of anodic oxides (leakage in
  [al_cell_pump](../examples/al_cell_pump)).
- S. L. Patrick, T. E. Bearden et al., U.S. Patent 6,362,718 (2002).
- R. L. Forward, *Phys. Rev. B* 30, 1700 (1984); D. C. Cole,
  H. E. Puthoff, *Phys. Rev. E* 48, 1562 (1993); A. Lambrecht,
  S. Reynaud, *Eur. Phys. J. D* 8, 309 (2000).
- M. Minotti, G. Modanese, arXiv:2408.05230 and arXiv:2603.26986;
  Iadicicco, Modanese, Verolino, *Gen. Rel. Grav.* 58, 106 (2026);
  Iadicicco, Verolino, arXiv:2608.21403; B. Bertotti, L. Iess,
  P. Tortora, *Nature* 425, 374 (2003).
- H. Eckardt, *ECE Unified Field Theory* Vol. 1 (2022), eqs. 4.211,
  6.189; H. Eckardt, D. Goel, UFT 456 (2024); the AIAS papers evaluated
  in [examples/AIAS](../examples/AIAS).
- [doc/pais-us10144532-comparison.md](../doc/pais-us10144532-comparison.md),
  [doc/tewari-references.md](../doc/tewari-references.md),
  [doc/kouril-vrba-verification.md](../doc/kouril-vrba-verification.md).
