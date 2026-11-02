# 4. Power and polyphase

Power in a circuit with one sinusoidal source is a scalar and a phase
angle. Power in a three-phase system with a rectifier on it is neither:
the three phases carry one quantity with a direction, and the harmonics
the rectifier draws carry power that is neither active nor reactive in the
textbook sense. kronspice handles both with the hypercomplex algebras of
CLAUDE.md §3b, quaternions for the three phases and geometric algebra for
the harmonics, and checks each against the conventional definitions. This
chapter builds the two tools and then uses them as the fitness of an
evolution, in [power_filter](../examples/power_filter): harmonic filters
for a six-pulse rectifier scored by geometric-algebra power.

## 4.1 Three phases as one quaternion

A set of three-phase values is written as the pure quaternion
x = x_a·i + x_b·j + x_c·k. The power-invariant Clarke transformation
(Clarke 1943), which turns abc into α, β and the zero sequence, is an
orthogonal matrix of determinant +1, so it is a rotation of that space:

```
√(2/3) · [ 1    −½     −½   ]
         [ 0    √3/2  −√3/2 ]
         [ 1/√2  1/√2  1/√2 ]
```

Park's transformation (Park 1929) adds a rotation by −θ about the
zero-sequence axis, so abc → dq0 is one unit quaternion r(θ) acting as
r·x·r*. That is Kron's rotating reference frame of the machine (chapter 3)
written as a rotor, and `kron/polyphase.go` keeps both forms, the rotor
and the matrix, as twins: on 500 random angles and vectors they agree to
1e-14.

Instantaneous power follows from the quaternion product. For pure
quaternions u·v = −u·v + u×v, and the conjugate of i is −i, so

```
s = v·i* = v·i − v×i
```

whose scalar part is p = v_a·i_a + v_b·i_b + v_c·i_c, the active power,
and whose vector part −v×i has the length of Akagi's instantaneous
imaginary power (Akagi, Kanazawa and Nabae 1984), in the vector form of
Peng and Lai (1996). Under any frame rotation s transforms as r·s·r*, so
p is invariant and the non-active vector turns with the frame. For a
balanced set the checks give p = 3·V·I·cos φ and |q| = 3·V·I·sin φ, to
1e-13 and 1e-10 (doc/M9f.md).

## 4.2 Power with harmonics: the geometric product

A periodic voltage with harmonics up to H is a vector in a Euclidean
space of 2H dimensions: each harmonic k contributes two orthonormal
functions over a period, √2·cos kωt and √2·sin kωt, and the RMS
components on them are the coordinates. The mean of v·i is then the dot
product V·I. In the geometric algebra Cl(2H) the full product

```
M = V·I = V⌋I + V∧I
```

has a scalar part, the active power P, and a bivector part. The
components of the bivector on the plane of one harmonic, cos k with
sin k, are the reactive powers Q_k = Im(V_k·I_k*), positive for an
inductive load. The components between two different harmonics h and k
are **distortion power**, each tied to its pair. Because |M| = |V|·|I| for
vectors,

```
S² = (V_rms·I_rms)² = P² + Σ_k Q_k² + Σ_{h≠k} D_hk²
```

and grouping V and I into fundamental and harmonic parts reproduces the
IEEE 1459 decomposition S² = S₁² + D_I² + D_V² + S_H². What the algebra
adds to the standard is the attribution: the bivector says which pairs
of harmonics carry the distortion, which is what a filter designer wants
to know. `analysis/power.go` builds M from the harmonics of a simulated
waveform (`HarmonicsOf`, a DFT over the last whole periods) and computes
the IEEE 1459 quantities conventionally beside it (`IEEE1459Of`). On 50
random band-limited waveforms with nine harmonics and a DC offset the two
agree on P, S, S₁, D_I, D_V, S_H, P₁ and Q₁ to 1e-9, and the identities
above hold to the same figure (Castro-Núñez and Castro-Puche 2012;
Montoya et al. 2021; doc/M9f.md).

## 4.3 Run it: the six-pulse rectifier

The example exports its circuits, and the unfiltered rectifier is a
netlist the tool runs as it is: a 60 Hz, 120 V (phase) supply with
0.05 Ω and 0.5 mH per line, a diode bridge, and a 2200 µF, 10 Ω DC load.

```
.model dr D(IS=1e-12 N=1 RS=0.001)
Va sa 0 SIN(0 169.7056274847714 60 0 0 -0)
RSa sa ma 0.05
LSa ma pa 0.0005
DUa pa pos dr
DDa neg pa dr
Vb sb 0 SIN(0 169.7056274847714 60 0 0 -120)
...
CDC pos neg 0.0022
RDC pos neg 10
RREF neg 0 100000
.tran 2.7777777777777779e-05 0.23333333333333334 0 2.7777777777777779e-05
```

The sources are peak amplitudes (120·√2), the phases 0°, −120° and
−240°, and the transient runs 14 cycles with a step of 1/600 of a
period. The report:

```
Time points  10968 recorded, 10967 steps accepted, 1943 rejected
Step size    26.5551 ps to 27.7778 µs
Interval     0 s to 233.333 ms

Energy ledger over the interval
  energy in            1.73106 kJ
  energy out           9.62753 pJ (returned to sources and ports)
  energy dissipated    1.74374 kJ
  stored: start → end  92.5359 J → 79.8255 J
  integration mismatch 32.4977 mJ (time-stepping error, not a physical channel)
  worst point residual 0 (relative)
```

The diodes' commutations are what reject 1,943 steps and drive the step
down to picoseconds; the ledger closes at every accepted point, and the
integration mismatch is 1.9e-5 of the energy delivered. The second
export adds, at each phase's point of common coupling, two series L–C
branches to a grounded star point, 1 mH with 270 µF and 150 µF, tuned to
300 Hz and 420 Hz, the 5th and 7th harmonics:

```
LF5a pa f5a 0.001
CF5a f5a r5a 0.00027
RF5a r5a 0 0.05
LF7a pa f7a 0.001
CF7a f7a r7a 0.00015
RF7a r7a 0 0.05
```

With the filters the run takes 8,637 points with 218 rejected, delivers
2.138 kJ, dissipates 2.060 kJ, and now returns 73.3 J to the sources:
the capacitors hand energy back to the supply during parts of each
cycle, which the ledger books as outflow, and the mismatch is 1.1e-5 of
the input.

## 4.4 Beyond the netlist: the analysis at the point of common coupling

The Go program is what turns those waveforms into a fitness. For each
phase it takes the voltage at the coupling point and the current from
the supply, resolves the last two of the 14 cycles into harmonics up to
the 25th, forms the geometric power and the IEEE 1459 quantities, and
sums P and S over the three phases:

```
PF = ΣP / ΣS
```

It also averages the quaternion instantaneous power s = v·i* of the three
phases over the same window, as an independent check that the scalar
part equals ΣP. The genetic algorithm (chapter 2) then evolves the four
branch values, L5, L7 from 0.1 to 10 mH and C5, C7 from 10 to 1000 µF,
all E12, with the single objective 1 − PF.

## 4.5 Results

The run on this machine, population 24 for 10 generations, 222
evaluations, took 9.4 s:

| | No filter | Evolved |
|---|---|---|
| Power factor | 0.9033 | 0.9559 |
| ΣP | 7389 W | 7673 W |
| Quaternion mean power | 7389 W | 7673 W |
| Phase a: S / P / Q₁ | 2731 VA / 2463 W / +515 var | 2701 VA / 2559 W / −73 var |
| Phase a: D_I | 1018 VA | 784 VA |
| THD_I / THD_V | 40.4 % / 9.9 % | 30.6 % / 12.8 % |
| Largest distortion pairs | 1–5 961 VA, 1–7 333 VA, 1–11 192 VA | 1–5 578 VA, 1–7 461 VA, 1–11 156 VA |

The GA and IEEE 1459 columns agree per phase (phase a's S by the two
routes, 2731 VA and the GA's P² + Q² + D², differ by the commutation
notches above the 25th harmonic), and the quaternion mean equals ΣP in
both cases. The evolved branches are 3.9 mH with 82 µF, resonant at
281 Hz, and 0.15 mH with 10 µF, resonant at 4109 Hz. They do not sit at
300 Hz and 420 Hz, and the bivector explains why: the 5th-harmonic branch
tuned a little low removes most of the 1–5 distortion power, from 961 to
578 VA, while the capacitors' leading fundamental reactive power swings
Q₁ from +515 var to −73 var, and that trade is what maximizes the power
factor; the second branch here is mostly a capacitor. The milestone's
longer run, 30 generations of 20, reached 0.957 with branches at 284 Hz
and 673 Hz (doc/M9f.md). The program exits non-zero if the evolved filter
fails to improve on the bare rectifier.

## 4.6 What to keep

- Three phases are one pure quaternion; Clarke and Park are rotors, and
  the quaternion product v·i* gives active power as its scalar part and
  Akagi's imaginary power as its vector part, invariant and co-rotating
  under any frame change.
- With harmonics, the geometric product V·I splits apparent power into
  active, reactive per harmonic and distortion per pair of harmonics,
  and reproduces IEEE 1459 exactly.
- The pair attribution is a design tool: it names the harmonics a filter
  must take, and shows when the best power factor is not at the nominal
  tuning.
- Every transient behind these numbers carries its ledger, and a filter
  that returns energy to the supply shows it as outflow, not as a gain.

## Run it yourself

```
go run ./cmd/kronspice examples/power_filter/export/power_filter_unfiltered.cir
go run ./cmd/kronspice examples/power_filter/export/power_filter_tuned_5_7.cir
go run ./examples/power_filter                       # ≈ 10 s; -gens 30 -pop 20 for the milestone's run
go test ./examples/power_filter                      # GA against IEEE 1459; -short skips the GA
go test ./kron -run 'Clarke|Park|Quaternion'         # the rotor twins
```

## Sources

- IEEE Std 1459-2010, *Definitions for the Measurement of Electric Power
  Quantities Under Sinusoidal, Nonsinusoidal, Balanced, or Unbalanced
  Conditions*.
- M. Castro-Núñez, R. Castro-Puche, "The IEEE Standard 1459, the CPC
  power theory, and geometric algebra in circuits with nonsinusoidal
  sources and linear loads," *IEEE Trans. Circuits Syst. I* 59(12), 2012.
- F. G. Montoya et al., *Mathematics* 9(11):1259, 2021 (geometric-algebra
  power).
- E. Clarke, *Circuit Analysis of A-C Power Systems*, Vol. I, Wiley, 1943,
  ch. X; R. H. Park, "Two-reaction theory of synchronous machines," *AIEE
  Trans.* 48, 1929.
- H. Akagi, Y. Kanazawa, A. Nabae, *IEEE Trans. Ind. Appl.* IA-20(3),
  1984; F. Z. Peng, J.-S. Lai, *IEEE Trans. Instrum. Meas.* 45(1), 1996
  (instantaneous reactive power and its vector form).
- A. S. F. Komeno et al., "Three-phase quaternion power in three-wire
  systems," *Renewable Energy and Power Quality Journal*.
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley, 1951,
  ch. 4 (the rotating reference frame).
- [doc/M9f.md](../doc/M9f.md): the quaternion polyphase and GA power
  checks and the milestone's run.
