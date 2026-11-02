# 3. Kron's machines

Gabriel Kron treated a circuit as a geometric object: branches on their
own form a primitive network, connecting them is a change of reference
frame, and a rotating machine is a circuit whose inductances depend on a
coordinate, the rotor angle, so that its voltage equation is the equation
of motion of a point in a curved space (CLAUDE.md §3a; §3l, curvature 3).
kronspice carries this as the N machine element, and this chapter runs it
through seven examples: an induction motor started in two of Kron's
frames ([kron_machine](../examples/kron_machine)), parametric pumping
three ways ([kron_machine](../examples/kron_machine),
[parametric_pump](../examples/parametric_pump),
[parametric_lcr](../examples/parametric_lcr)), Tesla's homopolar dynamo
([homopolar](../examples/homopolar)), the same machine with its current
through its bearings ([homopolar_bearing](../examples/homopolar_bearing),
[bearing_motor](../examples/bearing_motor)), and magnetic resonance as a
rotating frame ([kron_nmr](../examples/kron_nmr)). The last four are also
where three register hypotheses, H13, H14 and H16, meet their standard
predictions.

## 3.1 A circuit as a geometric object

Kron starts from the **primitive network**: every branch by itself,
unconnected, with its impedance tensor Z and its sources, so that
`e = Z·i`. Connecting the branches is a **connection tensor** C that
expresses the branch currents in terms of a smaller set, mesh currents or
node-pair voltages, `i = C·i′`. The connected network's equations follow
by transformation, `Z′ = Cᵀ·Z·C` and `e′ = Cᵀ·e`, and power is invariant
under the change, `e′·i′ = e·i`. kronspice builds both the mesh and the
node-pair forms (`kron/`), and on every linear circuit in the test suite
they agree with modified nodal analysis to 1e-9 or better (milestone
M9a). A machine is a network whose connection tensor depends on time.

## 3.2 The machine element

Kron's primitive machine has windings whose self and mutual inductances
depend on the rotor angle θ. With the inductances L(θ) and the rotor
inertia J together as the metric tensor, the voltage equation in his
**frame I**, axes fixed to the iron, is (Kron 1951, Epilogue eq. 3)

```
e_α = R_αβ·i^β + L_αβ·di^β/dt + [βγ, α]·i^β·i^γ
```

with the Christoffel symbol of eq. 4 built from the derivatives of L,
and the torque ½·iᵀ·(∂L/∂θ)·i on the shaft. In his **frame II(b)**,
axes that turn with the rotor, the inductances become constant and the
rotation reappears as a speed-voltage term (eq. 9):

```
e = R·i + d(L·i)/dt + G·i·dθ/dt        torque = iᵀ·G·i
```

G is the speed matrix, and G·i is the flux density of the armature.
kronspice's N element implements both: a `MachineSpec` with inductances
given as harmonics in θ, solved as written (frame I), or transformed by
`kron.RotatingFrame` into constant L′ and a G′ (frame II). The element
has one electrical port per winding and one **mechanical port**, a node
in the mechanical domain whose across variable is angular velocity and
whose through variable is torque, so shaft power is on the ledger like
any other channel (CLAUDE.md §3.1).

## 3.3 Run it: the induction motor in the rotor frame

[kron_machine](../examples/kron_machine) exports its circuits, and the
frame-II motor is a netlist the command-line tool runs as it is:

```
two-phase induction motor
.model mm MACHINE(G_3_2=0.11 G_3_4=0.12 G_4_1=-0.11 G_4_3=-0.12 L_1_1=0.12 L_1_3=0.11 L_2_2=0.12 L_2_4=0.11 L_3_3=0.12 L_4_4=0.12 N=4)
.domain mechanical w
VA sa 0 SIN(0 100 50 0 0 90)
VB sb 0 SIN(0 100 50 0 0 0)
RSA sa ia 0.5
RSB sb ib 0.5
RRA ra 0 0.8
RRB rb 0 0.8
N1 ia 0 ib 0 ra 0 rb 0 w 0 mm
CJ w 0 0.002 IC=0
RF w 0 5000
.tran 0.5m 1.5 uic
.end
```

Read the lines. The `.model … MACHINE` card holds the transformed
tensors: four windings (`N=4`), self-inductances 0.12 H, stator–rotor
mutuals 0.11 H, and the speed matrix entries `G_3_2`, `G_3_4`, `G_4_1`,
`G_4_3`, which appear only in the rotor rows 3 and 4, as frame II
requires (the program prints the full L′ and G′ matrices when it runs).
`.domain mechanical w` declares node `w` as a shaft. The `N1` line lists
the two terminals of each of the four windings, then the mechanical node
and its reference, then the model. `CJ`, a 0.002 F capacitor on the
shaft node, is the rotor inertia of 0.002 kg·m², and `RF`, 5 kΩ to
ground, is a friction load of 0.2 mN·m per rad/s: in the mechanical
domain a capacitor stores kinetic energy and a resistor dissipates
torque times speed. The two sources are a balanced two-phase 100 V, 50 Hz
supply, and the rotor windings are shorted through 0.8 Ω.

Both frames start the motor from rest over 1.5 s. The program's own
output:

| Frame | Final speed | Time points | Supply | Heat | Stored (magnetic + kinetic) | Ledger mismatch |
|---|---|---|---|---|---|---|
| I, holonomic, L(θ) | 313.568 rad/s | 7953 | 326.10 J | 227.29 J | 98.75 J | 9.8e-5 |
| II(b), rotor frame, G·i·ω | 313.696 rad/s | 8439 | 326.21 J | 227.38 J | 98.83 J | −1.4e-6 |

Synchronous speed is 314.16 rad/s. The two frames, which solve different
equations, agree on the final speed within 1e-3 and over the whole
start-up on speed within 1.1e-4 and on stator current within 5.8e-4 of
its peak ([doc/M9c1.md](../doc/M9c1.md)). Each ledger closes, supply
equal to heat plus stored energy, and the Christoffel form of eq. 3,
built by `kron.Christoffel` from the metric, reproduces the element's
voltages and torque to 1e-12 in the `kron` tests.

## 3.4 Parametric pumping, three ways

A machine whose coordinate is driven at constant speed is a time-varying
inductance. The second export of kron_machine is that circuit:

```
parametric
.model pump MACHINE(L_1_1=0.001 LC_2_1_1=5e-05 N=1)
.domain mechanical w
N1 a b w 0 pump
C1 a 0 2.533029591058444e-05 IC=1
VW w 0 DC 6283.185307179587
R1 b 0 0.07853981633974484
.tran 10u 60m uic
.end
```

One winding of L0 = 1 mH with a second-harmonic cosine term of
0.05 mH (`LC_2_1_1`), its shaft held at ω0 = 2π·1 kHz by the voltage
source `VW`, so the inductance is L0·(1 + 0.05·cos 2ω0t). In series with
C = 1/(ω0²L0) and a resistor it is an LC tank pumped at twice its
resonance, the textbook parametric oscillator. Its envelope grows at the
rate Mathieu's equation gives, m·ω0/4 less the damping R/2L0 (Kouřil &
Vrba §10.2.2):

| Damping, × pump rate | Growth, run | Mathieu | Pump work | Heat |
|---|---|---|---|---|
| 0 | 78.550 /s | 78.540 /s | 0.07876 J | 5.1e-10 J |
| 0.5 | 39.273 /s | 39.270 /s | 1.396 mJ | 0.707 mJ |
| 2 | −78.571 /s | −78.540 /s | 4.38 µJ | 17.1 µJ |

Every joule the tank gains enters through the pump's mechanical port:
pump work equals the tank's gain plus the heat.

The capacitive twin is [parametric_pump](../examples/parametric_pump),
one netlist with no Go program at all, using the `SCALE=` extension that
multiplies a capacitor's charge curve by a waveform:

```
C1 a 0 1u SCALE=SIN(1 0.1 10.066k 0 0 180) IC=1
L1 a b 1m
R1 b 0 0.1
.tran 1u 2m UIC
```

The capacitance falls while the voltage peaks, and the work of whatever
squeezes the plates is booked on the capacitor's mechanical port. The
tool prints:

```
  energy in            11.6963 µJ
  energy out           1.33306 µJ (returned to sources and ports)
  energy dissipated    611.799 nJ
  stored: start → end  500 nJ → 10.2503 µJ
  integration mismatch 1.1565 nJ (time-stepping error, not a physical channel)
```

The net work of the pump, 10.36 µJ, is the stored-energy gain of
9.75 µJ plus 0.61 µJ of heat; the envelope grows at 740.8 /s against
Mathieu's 740.6 /s. Nothing comes from nowhere.

That is the standard prediction beside which register entry **H16** is
read. Amador and Eckardt (2012) drive a series R, L, C(t) loop with
C(t) = C·(1 + a0·sin(ω1t + φ)) and report a coefficient of performance of
about 200 against the source's input, read as energy from the time
domain. [parametric_lcr](../examples/parametric_lcr) reproduces every
case of the paper with the same time-varying capacitor and puts the one
agent the paper's energy sum omits on the ledger, the pump that varies
C. From the run:

| Paper case | kronspice | COP by source alone | By source + pump |
|---|---|---|---|
| Fig. 4, constant C at f0 | peak 0.25 A = U0/R | 0.964 | 0.964 |
| Fig. 5, C(1 + 0.2 sin 2ω0t) at f0, 30 ms | 87.3 A and growing, σ = 245 1/s as Mathieu predicts; pump 256.2 J → R 165.2 J + stored 90.7 J | −1950 | 0.645 |
| Figs. 6–7, 0.97 f0, φ = −90° | bounded 16.1 A, 2.43 kW; source 6.75 J, pump 1139 J → R 1143 J | 169 | 0.998 |
| Fig. 8, sin⁷ pump at 1.07 ω0, a0 = 0.2 | 0.123 A, no growth | 0.987 | 0.926 |

The paper's highlighted case gives its COP of 169 only because the
source is credited with the heat while the pump that supplied 1139 J is
left out; with the pump counted, the ratio is 0.998 and below 1 in every
case. Kouřil and Vrba's threshold explains which cases grow: a0 = 0.2 is
an elastance depth of 0.202 against the tuned self-excitation threshold
2/Q = 0.136, so pumping at 2ω0 self-excites, while the 4.2 % detuning of
the 0.97 f0 cases raises the threshold to 0.214 and makes the circuit a
bounded parametric amplifier. The register keeps H16 as *untested*,
reproduced, with the bench test stated: measure the power that drives
C(t), the motor shaft or the varactor's bias supply.

## 3.5 Tesla's homopolar dynamo as machine elements

Tesla's U.S. Patent 406,968 (1889) puts two conducting disks on one
shaft in fields of opposite polarity, joins their rims by a belt and takes
the current from the two axles, so no brush bears on a fast rim.
[homopolar](../examples/homopolar) derives the EMF of a spinning disk from
the Lorentz force by grade selection (CLAUDE.md §3f, `em.MotionalEMF`),
which must equal ½·B·ω·(r_o² − r_i²), and then builds the machine from two
N elements. The exported netlist shows how:

```
.model da MACHINE(G_1_2=0.008960000000000001 L_1_1=1e-09 L_2_2=0.001 N=2)
.model db MACHINE(G_1_2=-0.008960000000000001 L_1_1=1e-09 L_2_2=0.001 N=2)
.domain mechanical w1 w2
IFA 0 fa DC 1
VW1 w1 0 DC 314.1592653589793
NA rim1 c1 fa 0 w1 0 da
RBELT rim1 rim2 0.0001
NB rim2 c2 fb 0 w2 0 db
RL t1 t2 0.01
```

Each disk is a two-winding machine: winding 1 is the disk's radial path
and winding 2 the field, driven at 1 A by `IFA`, with the one speed
entry `G_1_2 = ½·B·(r_o² − r_i²)` = 0.00896 V·s. The element then gives
the speed voltage ω·G·i_f in the disk and the torque G·i_disk·i_f on the
shaft; the second disk's G has the opposite sign for the opposite field.
Each shaft node is held at 314.16 rad/s (3000 r/min) by a mechanical
source. The tool's operating point for the reference design, 15 cm
copper disks at 0.8 T into a 10 mΩ load:

```
  VW1    voltage source  w1→0       314.159 V    -4.94393 A   -1.55318 kW
  NA     machine         rim1→c1    2.81487 V    -551.778 A   -227.374 fW
  RL     resistor        t1→t2      -5.51778 V   -551.778 A   3.04459 kW
  inflow 3.10636 kW   dissipation 3.10636 kW   residual -1.62049 nW (relative 0)
```

In the mechanical domain the shaft source's "current" of 4.944 A is a
torque of 4.944 N·m and its 1.553 kW is shaft power, booked as inflow on
the mechanical channel; the machine elements convert it inside the
network, so their own rows are 0 W. The program's comparison:

| | Two disks, belt (Tesla) | Single disk, rim brush |
|---|---|---|
| EMF per disk, Lorentz force | 2.8149 V | 2.8149 V |
| ½Bω(r_o² − r_i²) | 2.8149 V | 2.8149 V |
| Load | 5.5178 V, 3044.6 W | 2.6678 V, 711.7 W |
| Shaft (ledger inflow) | 3106.4 W | 750.9 W |
| Efficiency | 0.9801 | 0.9477 |

The two disks' EMFs add, and the belt costs far less than a brush on a
rim moving at 47 m/s. Evolved for at least 2 kW with the least loss
(seed 1, 562 evaluations, per the README), a 0.361 m disk at 1.22 T into
1 Ω puts 2485 W into the load from 2486 W at the shaft. The contact
resistances are declared, not measured; the tests hold the EMF to the
formula to 1e-12 and the ledger residual below 1e-12.

## 3.6 Current through the bearings: the P element and the registry

A current through a rotating ball bearing drives the shaft with a torque
that scales as I², the Huber effect, and kronspice carries it as the P
element ([doc/bearing-motor.md](../doc/bearing-motor.md)). Its
coefficients come only from fits to measured tables, and
[bearing_motor](../examples/bearing_motor) holds the registry of what has
been measured, each entry with its source and a status of *measured*,
*estimate* or *claim*: steady currents of 43.5–70.15 A with self-rotation
only at high current (Moyssides & Hatzikonstantinou 1990), torque
proportional to I² (Choo, Soong & Abbott 2005), mechanical output over
electrical input of 0.10 (Marinov 1989, calorimeter, data unpublished).
Marinov's claim that the heat is 8 % higher with rotation than without is
register entry **H14**, recorded as a claim and not used. No open table
fixes κ and ω0, so the element refuses to run without fitted values and
the program lists the nine papers whose tables would supply them.

The one calculation the registry allows is a bound. Tewari's Space Power
Generator test (**H13**, [space_power_generator](../examples/space_power_generator))
reports 3438 A at a no-load EMF of 3.41 V from a drive motor supplying
16.7 N·m. If all of E0·I were delivered, the shaft would need 41.5 N·m;
the gap of 24.8 N·m, 7.01 kW, is what a bearing torque would have to
supply. The exported netlist puts that bound into one loop:

```
.domain mechanical w
VE0 e 0 DC 3.41
PB e x w 0 1e-05 KAPPA=2.7548232420067764e-06 W0=282.7433388230814
RLOOP x 0 0.00038864400917643246
VW w 0 DC 282.7433388230814
```

The `PB` line is the bearing: its two electrical terminals, its shaft
node and reference, a contact resistance, and the fitted coefficients
κ and ω0 (here set to exactly the values that would close the gap). The
tool's operating point:

```
  VE0    voltage source  e→0    3.41 V     -3.438 kA  -11.7236 kW
  PB     bearing motor   e→x    2.07384 V  3.438 kA   118.198 W
  RLOOP  resistor        x→0    1.33616 V  3.438 kA   4.59371 kW
  VW     voltage source  w→0    282.743 V  24.7987 A  7.01167 kW
  inflow 11.7236 kW   outflow 7.01167 kW   dissipation 4.71191 kW   residual 0 W
```

A bearing motor that supplies 24.8 N·m to the shaft takes 2.07 V of the
3.41 V as back-EMF, so the rest of the loop gets 1.34 V and 4.59 kW, not
11.72 kW. That 1.34 V is the same loaded-EMF bound the SPG ledger gives
without any bearing: a bearing motor moves shaft work around and adds
none. In thermal mode at the measured efficiency of 0.10 it would need
70 kW of electrical input, more than the machine makes, and is excluded.

[homopolar_bearing](../examples/homopolar_bearing) then asks the design
question: Tesla's machine with its load current through the shaft
bearings, each bearing a P element. The bore removes magnetic area, but
larger bearings carry more current before damage (0.1–0.3 A/mm² of
Hertzian contact, 0.73 V contact drop, Safdarzadeh et al. 2024). From
its run_output.txt, load matched to the bearing bank's current limit:

| Bore | Magnetic area kept | EMF | Current limit | Load | Drive torque per watt |
|---|---|---|---|---|---|
| 10 mm | 99.9 % | 10.788 V | 0.238 A | 2.22 W | 0.0073 N·m/W |
| 40 mm | 98.2 % | 10.608 V | 0.600 A | 5.49 W | 0.0107 N·m/W |
| 120 mm | 84.0 % | 9.072 V | 1.248 A | 9.50 W | 0.0171 N·m/W |

A through-bearing unit is a watts-to-tens-of-watts machine. Evolved for
20 W with the least drive torque per watt, the design with insulated
bearings and a shaft collector wins by 7.4× (0.00135 against
0.0101 N·m/W), and in every case the drive torque per load watt is at
least 1/ω: `TestNoFreeTorque` holds the shaft power to the load plus the
measured losses.

## 3.7 Magnetic resonance as a rotating frame

Bloch's equations, dM/dt = γ·M × B with T1 and T2 relaxation, are a
rotation of the magnetization about the field, and a rotation of currents
at a fixed rate is exactly frame II(b): constant inductances and a
skew-symmetric G, so that L·dM/dt = −G·M. [kron_nmr](../examples/kron_nmr)
builds a proton at 1 T as a three-winding machine whose currents are Mx,
My and Mz, with G = L·γ·[B]ₓ and the shaft turning at ω = 1. A skew G has
zero torque iᵀ·G·i, so the precession takes no shaft power. Relaxation
is resistance on the windings, T2 = L/R₂ and T1 = L/R₁, with the heat on
the ledger, and a fourth open winding coupled to x is the pickup coil.
The NMR rotating frame is Kron's rotating reference frame, Park's
transform (CLAUDE.md §3b): at the RF frequency the RF field is a constant
x-field and the static field becomes B₀ − ω_rf/γ. The run:

| Check | Machine | Expected |
|---|---|---|
| Free precession | 42.5766 MHz | γB₀/2π = 42.5775 MHz (−2×10⁻⁵, trapezoidal phase error at 400 points per period) |
| Transverse decay at 4T2 | 0.01832 | e⁻⁴ = 0.01832 |
| Rabi nutation, B₁ = 1 mT | period 23.487 µs; shaft power 0 | 2π/γB₁ |
| Pickup peak after a π/2 pulse | 0.02676 V | M_p·γB₀·M₀ = 0.02675 V |
| Continuous-wave line, B₁ = 10 µT | absorption 0.00535 at resonance, 0.00054 at Δω·T2 = ±3 | Bloch's steady state, to 5 digits |
| Ledger, free precession | source 8.02×10⁻¹⁰ J, heat 9.83×10⁻¹⁰ J | mismatch −3×10⁻¹⁴ J |

The frequencies are measured inputs, so their match is by construction;
what the machine adds is Bloch's dynamics, confirmed since 1946 in every
spectrometer. The reason the example matters to the register is its last
table. A ring of charge whose mass circulates with it has γ = e/2m, the
g = 1 of chapter 9's [kron_electron](../examples/kron_electron) (**H21**):

| At 1 T | Measured | Ring, g = 1 |
|---|---|---|
| Proton NMR | 42.578 MHz (g_p = 5.5857) | 7.623 MHz |
| Electron ESR | 28.025 GHz (g_e = 2.0023) | 13.996 GHz |

A ring electron would resonate at half the ESR frequency. The standard
prediction beside it is the Dirac electron with QED,
g = 2.002 319 304 36, and the measured proton g.

## 3.8 What to keep

- A machine is a circuit with a coordinate: inductances L(θ) in frame I,
  or constant L′ and a speed matrix G in frame II. The two frames must
  agree, and in kronspice they do to 1e-3 on speed over a full start-up.
- The mechanical port is a channel like any other. Shaft power is inflow
  or outflow on the ledger, with torque as the through variable and
  angular velocity across.
- Parametric growth is paid for by whatever varies the element. A COP
  that omits the pump is not a measurement of anything.
- Coefficients for a measured effect come only from measured tables. When
  no table exists the element refuses to run, and the registry says so.
- A frame change can solve a problem from another field: Bloch's
  equations are a machine in frame II(b), and the NMR rotating frame is
  Park's.

## Run it yourself

```
go run ./cmd/kronspice examples/kron_machine/export/kron_machine_motor_frame2.cir
go run ./cmd/kronspice examples/kron_machine/export/kron_machine_parametric.cir
go run ./cmd/kronspice examples/parametric_pump/pump.cir
go run ./cmd/kronspice examples/homopolar/export/homopolar_two_disk_belt.cir
go run ./cmd/kronspice examples/bearing_motor/export/bearing_motor_spg_loop.cir
go run ./examples/kron_machine
go run ./examples/parametric_lcr
go run ./examples/homopolar -gens 0        # reference designs only; omit -gens 0 to evolve
go run ./examples/homopolar_bearing        # the sweeps and GA above (minutes)
go run ./examples/bearing_motor
go run ./examples/kron_nmr -gif viz/nmr.gif
```

## Sources

- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley, 1951,
  Epilogue pp. 253–265: eqs. 3–4 (frame I, Christoffel symbols), eq. 9
  (frame II(b), the speed matrix G), part 2 (the absolute time
  derivative and the equivalent circuits). Scan in `reference/`.
- G. Kron, "Non-Riemannian Dynamics of Rotating Electrical Machinery,"
  *J. Math. Phys.* 13, 103–194, 1934.
- G. Kron, *Tensor Analysis of Networks*, Wiley, 1939; *Tensors for
  Circuits*, Dover, 1959 (the primitive network and the connection tensor).
- F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits*, Ellis
  Horwood, 1988, §10.2.2; [doc/kouril-vrba-verification.md](../doc/kouril-vrba-verification.md).
- L. Mandelstam, N. Papalexi, *Z. Phys.* 73, 223, 1932.
- F. D. Amador, H. Eckardt, "Simulation of a Parametric Resonance
  Circuit," A.I.A.S./UPITEC, 30 March 2012 (hypothesis source for H16).
- N. Tesla, "Dynamo-Electric Machine," U.S. Patent 406,968, 1889.
- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16
  (reference only; CLAUDE.md §3f).
- J. Huber, *Elektrotechnik und Maschinenbau* 76(8), 169–174, 1959;
  P. G. Moyssides, P. Hatzikonstantinou, *IEEE Trans. Magn.* 26, 1274,
  1990; J. L. Choo, W. L. Soong, D. Abbott, *Proc. SPIE* 5649, 700, 2005;
  S. Marinov, *Electronics & Wireless World*, April 1989; K. T. McDonald,
  "Ball-Bearing Motor," Princeton, 2020; [doc/bearing-motor.md](../doc/bearing-motor.md).
- O. Safdarzadeh, A. Farahi, A. Binder et al., *Lubricants* 12(5), 148,
  2024 (bearing voltage and damage limit).
- P. Tewari, *Universal Principles of Space and Matter*, 2017, §6.7
  (hypothesis source for H13).
- F. Bloch, *Phys. Rev.* 70, 460, 1946; E. M. Purcell, H. C. Torrey,
  R. V. Pound, *Phys. Rev.* 69, 37, 1946; E. Zavoisky, *J. Phys. USSR* 9,
  211, 1945; CODATA 2018 (γ_p, γ_e, g_p, g_e).
- [doc/M9c1.md](../doc/M9c1.md), [doc/M9n.md](../doc/M9n.md): the
  milestone records with the full test tables.
