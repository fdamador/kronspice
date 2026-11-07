# 7. Networks beyond circuits

Kron wrote that separate networks joined by a connection tensor may be
"electrical, mechanical, acoustical, optical and so on" (*Tensor Analysis
of Networks*, 1939, p. 133), and Branin later made the analogy exact:
every energy domain has an across variable and a through variable whose
product is power, so a mechanical or acoustic system is a network the
same engine solves (CLAUDE.md §3a.1 item 3). This chapter follows that
idea outward from the circuit in four steps: a piezoelectric transducer
driving a flask of water ([acoustic_drive](../examples/acoustic_drive)),
the vibrations of a water molecule
([water_vibration](../examples/water_vibration)), the four DNA bases and
their pairs ([dna_vibration](../examples/dna_vibration)), and the
Schrödinger equation of the simplest chemical bond solved as Kron's 1945
network ([kron_schrodinger](../examples/kron_schrodinger)). The first is
a design problem; the other three are the chemistry plan of
[doc/chemistry-plan.md](../doc/chemistry-plan.md), each checked against a
measured spectrum or an exact solution.

## 7.1 Across, through and the transducer

Every node of a kronspice circuit belongs to one domain, and ground is the
common reference of all of them: zero velocity, ambient pressure
([doc/domains.md](../doc/domains.md)).

| Domain | Across ("voltage") | Through ("current") | R | C | L |
|---|---|---|---|---|---|
| electrical | voltage [V] | current [A] | resistance | capacitance | inductance |
| mechanical | velocity [m/s] | force [N] | 1/damping | mass [kg] | 1/stiffness [m/N] |
| acoustic | pressure [Pa] | volume velocity [m³/s] | acoustic resistance | compliance [m³/Pa] | inertance [kg/m⁴] |

The mechanical row is the mobility analogy of Firestone (1933): velocity
across and force through, so a mass is a capacitor to ground and
Kirchhoff's current law at a node is the force balance at a point. The
acoustic row is pressure across and volume velocity through.

What joins two domains is a transducer, and in this analogy the classic
transducers are gyrators. kronspice's `Y` element is an ideal, lossless
gyrator, `i1 = g·V2` and `i2 = −g·V1`, so `V1·i1 + V2·i2 = 0`: power in at
one port leaves at the other, inside the network. A piezoelectric element
of force factor φ has `g = φ`, with its clamped capacitance C0 as an
ordinary capacitor across the electrical port; a piston of area S between
mechanical and acoustic ports has `g = S`. In Kron's primitive network a
gyrator is two axes coupled by an antisymmetric mutual admittance, a
non-reciprocal element the mesh and node-pair forms carry like any other.
Two netlist extensions declare the domains and send a resistor's power to
a named exterior:

```
.domain mechanical v
.domain acoustic p
Y1 e 0 v 0 0.6              * transducer: i1 = g·V2, i2 = −g·V1
RA p 0 6.9e8 EXT=acoustic   * absorbed power leaves to the acoustic exterior
```

The ledger then shows, per channel, where the energy went: electrical in
at the source, thermal out at dampers and the liquid's loss, acoustic out
as radiated sound. That is Kron's exterior junction per channel of
chapter 1, now with three channels carrying power.

## 7.2 Run it: a piezo ring driving a flask

The exported netlist of the matched water run in
[acoustic_drive](../examples/acoustic_drive) is one chain across three
domains:

```
acoustic drive H2O
.domain mechanical v
.domain acoustic p
VGEN g 0 DC 0 AC 1 0
RGEN g a 50
LMATCH a e 0.010111626690817235
RLMATCH e x 0.5
C0 x 0 4e-09
YPZT x 0 v 0 0.6
CMASS v 0 0.004
LSPRING v 0 1.0111626690817235e-08
RLOSS v 0 0.3333333333333333
YCOUPLE v 0 p 0 0.0002
CMODE p 0 7.650334663231506e-15
LMODE p 0 5286.893782262894
RMODE p 0 1.2469565900682e+12
.ac lin 2001 24274.559430000001 25776.078570000005
.end
```

Read it stage by stage. A 1 V generator with its 50 Ω drives, through a
matching inductor that tunes out the ring's 4 nF clamped capacitance, the
electrical port of the piezo. `YPZT` with φ = 0.6 N/V turns volts into
force and velocity into current. On the mechanical side the ring is a
4 g mass, a spring placed so the ring resonates at the flask's mode, and
a loss of 3 N·s/m written as its reciprocal. `YCOUPLE`, with an effective
area of 2 cm², turns velocity into volume velocity and pressure into
force. The last three elements are the breathing mode of a 3 cm sphere of
water: compliance, inertance and the parallel resistance that sets its
Q of 1500, from `fluid.SphericalResonator`. Run it:

```
go run ./cmd/kronspice -json examples/acoustic_drive/export/acoustic_drive_water_matched.cir
```

In the document, `ac.v.p` is the centre pressure per volt of generator
EMF. Its peak is at 25,024.6 Hz, where it reaches 69.4 kPa per volt;
`ac.v.v` at the same point is 0.279 mm/s per volt, the ring's velocity.
The ledger at that frequency reads inflow 7.27 mW, dissipation 7.27 mW,
outflow 0 and a relative residual of 0: one volt in, 7.27 mW out as
heat, booked where it is burned. The text report prints the same sweep as
a Bode table; the pressure node shows 68.9 dB at the first point.

The program adds what a netlist cannot: it builds the chain for H₂O and
D₂O with and without the matching inductor, sweeps each, finds the peak,
scales the drive so the centre pressure reaches 1.3 atm, the range in
which single-bubble sonoluminescence is observed (Brenner, Hilgenfeldt &
Lohse 2002), and reads the ledger at that drive.

| Run | f0 (peak) | Pressure per volt | Generator EMF | Source power | Into the liquid | Transducer loss | Generator's own 50 Ω |
|---|---|---|---|---|---|---|---|
| H₂O, no matching | 25.03 kHz (24.27) | 3.22 kPa/V | 40.9 V | 10.9 mW | 6.96 mW | 3.50 mW | 0.42 mW |
| H₂O, matched | 25.03 kHz (25.02) | 69.4 kPa/V | 1.90 V | 26.2 mW | 6.96 mW | 0.42 µW | 19.0 mW |
| D₂O, no matching | 23.40 kHz (22.69) | 3.20 kPa/V | 41.2 V | 10.5 mW | 6.72 mW | 3.26 mW | 0.51 mW |
| D₂O, matched | 23.40 kHz (23.40) | 72.2 kPa/V | 1.82 V | 23.5 mW | 6.72 mW | 0.39 µW | 16.6 mW |

The 1.90 V matched drive is the 1 V run scaled: 1.90² × 7.27 mW is the
26.2 mW in the table. The power the liquid takes at a given pressure is
fixed by the mode, P²/(2·R_par), so matching does not change it; what
matching changes is the drive voltage, 21 times more pressure per volt,
and where the rest of the power goes. Without the inductor the
transducer's own loss takes a third of the input. With it the matched
chain is a low impedance at resonance and most of the input is burned in
the generator's 50 Ω. Heavy water resonates 6.5 % lower than light water
because its sound speed is 1404 m/s against 1502 m/s, both from the
IAPWS-class equations of state in `fluid/`.

The test holds, for both liquids, that the liquid's heat equals
P²/(2·R_par) to 1e-9, that the residual is below 1e-9, that matching puts
the peak back at f0 within 0.1 %, and that the exported chain torn
automatically at `LMATCH` into two pieces agrees with MNA to 5e-12. The
flask Q, the coupling area and the transducer values are declared
assumptions, to be replaced by a hydrophone calibration of a real flask.

The bubble the chain drives is
[bubble_collapse](../examples/bubble_collapse), and it is a circuit in the
same acoustic domain. The liquid around it is an inertance ρ/4πR in
parallel with the pulsating sphere's radiation resistance ρc/4πR², the
wall's viscosity a resistance, and the argon a nonlinear compliance with a
hard core. Because the inertance depends on the volume, the circuit's own
charge, its equation carries Kron's Christoffel term ½M′UU_M. Without the
losses it is the Rayleigh–Plesset equation, and Keller–Miksis runs beside
it as the twin. At 1.3 atm the 4.5 µm argon bubble swells to 32.6 µm and
collapses to 0.73 µm. Per cycle it takes 14.3 nJ from the sound field and
returns 13.3 nJ as a sound pulse, 91 % of it within ±10 ns of the
minimum, and 1.0 nJ as viscous heat; the ledger residual is 6×10⁻¹⁰. Put
back on node p as its fundamental admittance, it draws 0.36 mW, 5 % of
what the flask's own losses take, and raises the generator drive from
1.898 V to 1.924 V. The same example compares the collapse with the
knotted light of Chapter 9: lossless, the collapse retraces itself when
reversed, like the knot; the radiated sound makes it irreversible; and in
the field solver its thermal flash shows no net circular polarization,
where the knot shows +0.9996.

## 7.3 Water's bonds as a circuit

The same mobility analogy, taken down to a molecule. In
[water_vibration](../examples/water_vibration) the six in-plane Cartesian
coordinates of O, H1 and H2 are nodes whose capacitance is the atomic mass,
and the three bond coordinates, two O–H stretches and the H–O–H bend, are
springs, inductors of value 1/k. What connects the two is Wilson's B
matrix: each spring's rate of change is a linear combination of the atomic
velocities, built as a chain of voltage-controlled voltage sources (E),
the spring force flows in the inductor, and current-controlled current
sources (F) hand that force back to every atom. The E/F pair is an ideal
multiport transformer, Kron's connection tensor between the bond frame and
the atom frame, and the stretch–stretch interaction f_rr is a mutual
inductance between the two stretch inductors, a `K` element. In the
netlist:

```
CM_ox ox 0 15.9949146196                  * oxygen's mass, x coordinate
CM_h1x h1x 0 1.00782503223                * hydrogen's
E_r1_ox r1_s0 0 ox 0 -0.790569884316532   * B row of stretch r1, summed by a chain of E sources
F_r1_ox ox 0 VS_r1 -0.790569884316532     * the stretch force back on the oxygen
K_rr L_r1 L_r2 0.0122798366009168         * the stretch–stretch interaction
I_fx 0 h1x AC 1                           * a force on H1 to excite the modes
.ac lin 4001 0.0610822823 0.549740541
```

Capacitance is in daltons and inductance in Å/mdyn, so the circuit's
ω = 1/√(LC) maps onto the molecule's through √(100 N/m / 1 u), and the
netlist's header gives the factor: circuit frequency in hertz × 8185.68 is
wavenumber in cm⁻¹. Run it:

```
go run ./cmd/kronspice -json examples/water_vibration/netlists/h2o.cir
```

The sweep has 4001 points. Summing the atom nodes' response magnitudes and
taking the maxima gives three peaks:

| Circuit frequency | × 8185.68 | Mode |
|---|---|---|
| 0.194852 Hz | 1595.0 cm⁻¹ | bend |
| 0.446756 Hz | 3657.0 cm⁻¹ | symmetric stretch |
| 0.458850 Hz | 3756.0 cm⁻¹ | antisymmetric stretch |

Those are the three measured gas-phase fundamentals of H₂O, and they are
exact because the force constants were fitted to them: f_r 7.6712, f_rr
−0.0942 and f_a 0.7030 mdyn/Å by Newton's method, three constants for
three lines. The program's AC peaks equal the GF-matrix eigenvalues to
1e-6: the circuit is the molecule's equations, not an approximation of
them.

The test of the model is the isotopes. Deuterium changes masses, not
bonds, so the same three constants with heavier capacitors must predict
D₂O, HDO and H₂¹⁸O with no refit:

| Molecule | Mode | Circuit AC peak (cm⁻¹) | Measured (cm⁻¹) | Error |
|---|---|---|---|---|
| D₂O | bend | 1162.51 | 1178.0 | −1.32 % |
| D₂O | stretch | 2647.26 | 2671.0 | −0.89 % |
| D₂O | stretch | 2752.05 | 2788.0 | −1.29 % |
| HDO | bend | 1395.34 | 1402.0 | −0.48 % |
| HDO | stretch | 2697.36 | 2727.0 | −1.09 % |
| HDO | stretch | 3708.98 | 3707.0 | +0.05 % |
| H₂¹⁸O | bend, stretch, stretch | 1589.15, 3647.51, 3740.70 | not in CCCBDB | prediction only |

The remaining error is anharmonicity. A measured fundamental sits below
its harmonic frequency by a shift that shrinks with heavier atoms, so a
harmonic field fitted to H₂O's fundamentals absorbs hydrogen's larger
shift, comes out soft, and under-predicts the deuterium modes by about
1 %; HDO's O–H stretch, closest to H₂O, is within 0.1 %. The plan's next
step is Morse inductors, nonlinear springs that give fundamentals and
overtones directly. The H₂¹⁸O lines are shown as a prediction because no
primary measurement has been checked, which is the chemistry plan's rule:
numbers that cannot be traced to a source are left out, not guessed.

## 7.4 DNA bases and base pairs

[dna_vibration](../examples/dna_vibration) carries the same construction
to an all-atom force field. The molecules are 9-methyladenine,
1-methylthymine, 9-methylguanine and 1-methylcytosine (the methyl stands
in for the sugar, as in gas-phase spectroscopy) and the A·T and G·C pairs,
cut from the B-DNA dodecamer d(CGCGAATTCGCG) and minimized in the AMBER
OL15 force field; `prep/prep.py` writes `molecules.json` from OpenMM's
parameter and structure files, and the Go code evaluates every term.

Atoms are capacitors; every AMBER bond and angle is a spring branch
through Wilson's B, as for water. Everything else in the force field,
torsions, out-of-plane terms, electrostatics and van der Waals, and the
pre-stress at the minimum, is the residual stiffness R = K − Σ k·BᵀB,
applied through displacement integrators and one voltage-controlled
current source per non-zero entry, so the circuit's stiffness is the
force field's K exactly up to entries below 1e-9 of the largest. In a pair,
every force between the two bases is its own group G_HB; in a stack the
forces between bases of different pairs are a third group G_ST. Light
enters as a current source at every atom, the field's force on the atom's
charge, and the absorption spectrum is the power the field puts in during
the AC sweep from 0 to 3700 cm⁻¹ in 2 cm⁻¹ steps. Release per mode comes
by two routes that must agree: the Einstein A coefficient from the IR
intensity, and Maxwell's radiated power of the vibrating dipole from
`em.LarmorPower`.

| Molecule | Atoms | Vibrations | Range (cm⁻¹) | d/λ | Thermal release at 300 K |
|---|---|---|---|---|---|
| A (9-methyladenine) | 18 | 48 | 49–3367 | 2.5e-4 | 2.58e-20 W |
| T (1-methylthymine) | 18 | 48 | 32–3314 | 2.1e-4 | 1.28e-20 W |
| G (9-methylguanine) | 19 | 51 | 30–3363 | 2.5e-4 | 2.90e-20 W |
| C (1-methylcytosine) | 16 | 42 | 13–3366 | 2.2e-4 | 2.65e-20 W |
| A·T pair | 36 | 102 | 13–3365 | 4.4e-4 | 3.68e-20 W |
| G·C pair | 35 | 99 | 18–3370 | 4.2e-4 | 5.20e-20 W |

Every molecule is at a true minimum, with all 3N−6 vibrations real, and
every one is lumped by the λ/10 rule by four orders of magnitude. Against
measurement the picture is honest rather than flattering:

| Band | Measured | Model | Difference |
|---|---|---|---|
| A NH₂ antisymmetric stretch | 3569 | 3367.1 | −5.7 % |
| T N3–H stretch | 3435 | 3314.1 | −3.5 % |
| T C4=O stretch | 1720 | 1745.7 | +1.5 % |
| G N1–H stretch | 3440 | 3308.7 | −3.8 % |
| G C6=O stretch | 1752 | 1544.0 | −11.9 % |
| C C2=O stretch | 1740 | 1685.0 | −3.2 % |
| A·T inter-base H-bond stretch | 94 (duplex in water) | 123.5 (one pair in vacuum) | +31 % |

The N–H stretches sit 3 to 6 % below the data, the harmonic AMBER value
against anharmonic measurements, and the measured red shifts of those
stretches on pairing are not reproduced: with fixed charges the model's
H-bonded stretches move up. Pairing itself comes out as it should. A·T
has two hydrogen bonds and G·C three, with donor–acceptor distances of
2.82 to 2.92 Å near the crystal values, and the modes with a quarter or
more of their energy in the inter-base forces are the low ones: A·T at 24,
45, 54, 72 and 123 cm⁻¹, G·C at 62, 63, 116 and 154 cm⁻¹. Stacking two to
four pairs turns each line into a band and adds rigid-pair modes below
200 cm⁻¹; the lowest mode falls from 12.7 cm⁻¹ for one A·T pair to 3.4 cm⁻¹
for four. As a side note, `sound/` holds the spectra rendered as audio at
3 Hz per cm⁻¹, loudness by the square root of the IR intensity, which is
a way of hearing that G·C is the stiffer link.

[rna_vibration](../examples/rna_vibration) repeats every one of these
experiments for RNA: 9-methyladenine, 1-methyluracil, 9-methylguanine and
1-methylcytosine, the A·U and G·C pairs and RNA's own wobble pair G·U,
cut from the NMR structure of the duplex r(GAGUGCUC)₂ (PDB 1EKA, which
carries its hydrogens) and minimized in the AMBER RNA force field OL3.
The bases come out as they did for DNA, with N–H stretches 3 to 6 % below
the data, and uracil's N3–H stretch 3.4 % below its argon-matrix value;
the C4=O stretch is spread over several ring modes and has no clean line.
The wobble pair keeps exactly its two hydrogen bonds, G N1–H···O2 U and
U N3–H···O6 G, and its stiffest pairing mode, at 131 cm⁻¹, sits beside
A·U's 134 cm⁻¹ and below G·C's 152 cm⁻¹: two hydrogen bonds, and about
the strength of A·U. Its stacks of two to four pairs repeat one real
A-form step of the duplex, GA/UC, since 1EKA has no AA/UU or CG step to
cut; they split lines into bands as the DNA stacks do, and the lowest
mode falls from 12.5 cm⁻¹ for one A·U pair to 6.4 cm⁻¹ for four.

## 7.5 Kron's Schrödinger network

In 1945 Kron showed that the Schrödinger equation is a ladder network:
series inductors carry the kinetic energy, shunt inductors to ground the
potential energy, shunt capacitors the energy eigenvalue, ψ is the
junction voltage to ground, and an eigenstate is a resonance
(CLAUDE.md §3a.1). [kron_schrodinger](../examples/kron_schrodinger) builds
that network for H₂⁺, two protons and one electron, with the hydrogen atom
as its check; both have exact solutions. The ground state is axially
symmetric, so the network lies on the (ρ, z) half-plane with the nuclei on
the axis, and the equation multiplied by ρ, so that the network is
reciprocal, is

```
−½ [∂ρ(ρ ∂ρψ) + ρ ∂z²ψ] + ρ V ψ = E ρ ψ,    V = −Σ Z_a / |r − R_a|
```

On a grid of spacing h with junctions at ρ_i = (i + ½)h, so that no
junction sits on the axis, this is A ψ = E B ψ with B = diag(ρ_i): series
branches of admittance ρ_face/(2h²) between neighbours, a shunt branch
ρ_i·V_i at each junction, a shunt capacitor ρ_i. The outer boundary is a
short circuit far from the nuclei, Kron's termination "by a short circuit
at a point of zero voltage".

The example finds the eigenvalue two ways. Kron's own driven network is
excited at a trial energy and its response iterated to the nearest
resonance, with energies at h, h/2 and h/4 extrapolated to h → 0. The
second way is a kronspice netlist. Kron's all-inductor network has loops of
shorts at DC, so the netlist is written in its dual form, capacitors and
inductors exchanged: series capacitors between junctions, shunt capacitors
ρ(V + V0) and shunt inductors 1/ρ, with V0 shifting the potential so every
element is positive. Then every eigenstate is a resonance at
ω² = 1/(E + V0), the ground state the highest. The coarse netlist has 1,950
junctions on an h = 0.4 a₀ grid and one AC current source on node n960:

```
go run ./cmd/kronspice -json examples/kron_schrodinger/netlists/h2plus_coarse.cir
```

It runs in three seconds, 241 frequency points through the sparse solver.
Summing the node magnitudes, the largest response is at 0.0800748 Hz, and
with V0 = 5.08769 Ha from the netlist's header, E = 1/ω² − V0 = −1.1372 Ha
on that 241-point grid. The example refines the peak by a parabola through
its top three points to −1.13848455 Ha, and the driven-network eigenvalue
on the same grid is −1.13848455 Ha; the test requires the two to agree to
1e-7. The exported results sheet also solves this 1,950-node circuit by
automatic tearing into 4 pieces joined by 179 links, equal to MNA.

Finite-difference energies lie below the exact ones, because a difference
scheme is not variational, and converge with an order set by the Coulomb
cusps, about 1.45:

| Problem | h = 0.2 a₀ | h = 0.1 | h = 0.05 | h → 0 | Exact |
|---|---|---|---|---|---|
| H atom | −0.51087970 | −0.50382944 | −0.50123643 | −0.499728 Ha (−13.5983 eV) | −0.5 Ha (−13.6057 eV) |
| H₂⁺ at R = 2 a₀ | −1.11745397 | −1.10780517 | −1.10429632 | −1.102291 | −1.1026342 |

Adding the nuclear repulsion 1/R and sweeping R gives the bond curve, and
a least-squares cubic through its minimum gives the bond:

| | Model | Exact |
|---|---|---|
| Bond length R_e | 1.998 a₀ = 1.057 Å | 1.997 a₀ = 1.057 Å |
| Bond energy D_e | 0.1024 Ha = 2.786 eV | 0.1026 Ha = 2.79 eV |

Water's ten interacting electrons need a self-consistent field on a 3-D
network, with a Poisson network for the repulsion, which the plan lists as
its next step and which is not done here.

## 7.6 What to keep

- A domain is an across variable and a through variable whose product is
  power. Declare the domain of every node, and let ground be the common
  reference.
- A transducer is a gyrator. It passes power between domains without
  loss, so the ledger can place every watt: source, damper, liquid,
  radiated sound.
- The mobility analogy makes a molecule a circuit whose resonances are its
  spectrum, and the connection tensor between bond and atom frames is an
  E/F transformer. Fit to one molecule, then predict its isotopes with no
  refit.
- Report the model against measurement as it comes out, errors included,
  and leave out numbers with no traceable source.
- Kron's Schrödinger network is a netlist: in its dual LC form the ground
  state is the highest resonance, and the circuit solver finds it to the
  same digits as the eigenvalue iteration.

## Run it yourself

```
go run ./cmd/kronspice -json examples/acoustic_drive/export/acoustic_drive_water_matched.cir
go run ./examples/acoustic_drive
go run ./cmd/kronspice -json examples/water_vibration/netlists/h2o.cir
go run ./examples/water_vibration
go run ./examples/dna_vibration -stacks=false
go run ./examples/rna_vibration -stacks=false
go run ./cmd/kronspice -json examples/kron_schrodinger/netlists/h2plus_coarse.cir
go run ./examples/kron_schrodinger
```

The water netlists `d2o.cir`, `hdo.cir` and `h218o.cir` sit beside
`h2o.cir`, and `examples/dna_vibration/netlists/` holds each base and pair
as a netlist; all of them also open in ngspice.

## Kron's polyhedron, checked against his numbers

Kron's late work extended the field network to a sequence of networks of rising
dimension. [kron_polyhedron](../examples/kron_polyhedron) builds his 1945
cavity and the 1962 curve-fitting polyhedron and compares them with every number
he printed. A planar reference frame turns out to be a vector of whole numbers,
which the genetic algorithm recovers from a printed table column. The
"oscillatory" column, which fits any function without adjustment, returns its
input: all 13 characteristics depend on only four values, and a frame orthogonal
to them fits anything, random numbers included. That is an inference from the
tables, and it says nothing either way about the physical claims of register
entry H20. Two columns and the full three-dimensional cell remain unreproduced.

## Sources

- G. Kron, *Tensor Analysis of Networks*, Wiley, 1939, Ch. V p. 133 —
  networks of any kind joined by a connection tensor.
- G. Kron, "Electric circuit models of the Schrödinger equation," *Phys.
  Rev.* 67, 39–43, 1945.
- F. A. Firestone, "A new analogy between mechanical and electrical
  systems," *J. Acoust. Soc. Am.* 4, 249, 1933 — the mobility analogy.
- F. H. Branin Jr., "The algebraic-topological basis for network analogies
  and the vector calculus," 1966.
- L. E. Kinsler, A. R. Frey, A. B. Coppens, J. V. Sanders, *Fundamentals
  of Acoustics*, 4th ed., Wiley, 2000, §7.5 and §9.3.
- M. P. Brenner, S. Hilgenfeldt, D. Lohse, "Single-bubble
  sonoluminescence," *Rev. Mod. Phys.* 74, 425, 2002.
- W. Wagner, A. Pruß, *J. Phys. Chem. Ref. Data* 31, 387, 2002 (H₂O);
  S. Herrig, M. Thol, A. H. Harvey, E. W. Lemmon, *J. Phys. Chem. Ref.
  Data* 47, 043102, 2018 (D₂O).
- E. B. Wilson, J. C. Decius, P. C. Cross, *Molecular Vibrations*,
  McGraw-Hill, 1955, Ch. 4 — the B matrix and the GF method.
- A. R. Hoy, P. R. Bunker, *J. Mol. Spectrosc.* 74, 1, 1979 — water's
  geometry; NIST CCCBDB (Shimanouchi; Huber & Herzberg) — measured
  fundamentals.
- AMBER DNA force field OL15: W. D. Cornell et al., *JACS* 117, 5179,
  1995; A. Pérez et al., *Biophys. J.* 92, 3817, 2007; M. Zgarbová et al.,
  *JCTC* 11, 5723, 2015.
- C. Plützer et al., *PCCP* 3, 5466, 2001 and *ChemPhysChem* 4, 838, 2003;
  M. J. Nowak et al., *PCCP* 21, 22857, 2019; A. Min et al., *Bull. Korean
  Chem. Soc.* 30, 3039, 2009; M. González-Jiménez et al., *Nat. Commun.*
  7, 11799, 2016 — the measured bands in §7.4.
- D. R. Bates, K. Ledsham, A. L. Stewart, *Phil. Trans. R. Soc. A* 246,
  215, 1953; J. M. Peek, *J. Chem. Phys.* 43, 3004, 1965 — exact H₂⁺.
- [doc/domains.md](../doc/domains.md), [doc/chemistry-plan.md](../doc/chemistry-plan.md).
