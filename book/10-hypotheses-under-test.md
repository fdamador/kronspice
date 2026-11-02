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
[bearing_motor](../examples/bearing_motor) supplying one bound.

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

## 10.6 The bounds recorded so far

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

In every case the engine found the energy, and in every case it was
already on the books: shaft work in H15, the pump in H16, the loaded EMF
in H13. H22 is the one entry about structure rather than energy: the
quark model and circuit theory are both gauge theories on a graph, and
the test found the difference between them to be exactly the difference
between SU(3) and U(1), with no circuit quantity on the non-abelian side. None of these is a verdict on the idea behind the entry. Each is
a statement of what the published equations or numbers contain, with the
measurement that would change the status written next to it.

## 10.7 What to keep

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
go run ./examples/eckardt_homopolar
go run ./examples/parametric_lcr
go run ./examples/space_power_generator
go run ./examples/bearing_motor
go run ./examples/pais_cavity
go run ./examples/quantumtricity
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
- [doc/pais-us10144532-comparison.md](../doc/pais-us10144532-comparison.md),
  [doc/tewari-references.md](../doc/tewari-references.md),
  [doc/kouril-vrba-verification.md](../doc/kouril-vrba-verification.md).
