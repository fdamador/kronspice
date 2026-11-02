# 5. Fields from Maxwell

A circuit simulator is the lumped limit of Maxwell's equations: parts much
smaller than a wavelength, fields confined to components, nothing
radiated. kronspice is built to stay honest about that limit (CLAUDE.md
§3c). The `em/` package holds the field equations themselves, in
Maxwell's own quaternion form and two others that must agree with it;
extracts R, L, C and K for real geometries from field solves, so that
layout becomes a gene; checks every circuit against the quasi-static
assumption; and carries a time-domain field solver for the cases where
the lumped picture ends. This chapter runs three examples through those
pieces: [pcb_layout](../examples/pcb_layout) evolves a microstrip pair
from a cross-section, [rotating_efield](../examples/rotating_efield)
asks whether a rotating electric field makes a magnetic one, and
[rate_of_change](../examples/rate_of_change) asks where higher time
derivatives enter, in the material and in the field equations.

## 5.1 One physics, three notations

Maxwell's *Treatise* (1873, Vol. II, Art. 618 and 619) states the field
equations in components and in Hamilton's quaternions. With
∇ = i ∂/∂x + j ∂/∂y + k ∂/∂z and Hamilton's signs, one product carries
both the divergence and the curl of a vector field:

```
∇v = −div v + curl v
```

so the induction is 𝔅 = V.∇𝔄 and the currents are 4π ℭ = V.∇ℌ, with the
displacement current ∂𝔇/∂t inside ℭ. The engine implements that form
(`em/maxwell.go`), the modern single equation ∇F = J/(ε₀c) in the
spacetime algebra Cl(1,3) with F = E + IcB (`em/biquat.go`), and plain
vector calculus (`em/vcalc.go`), and insists on a rule: the three are
the same physics, so any difference between them is a bug. The tests of
milestone M9h ([doc/M9h.md](../doc/M9h.md)) apply all three to
arbitrary non-Maxwellian fields and sources and require the residuals to
agree to 10⁻⁹; they check Hamilton's ∇ on polynomials to 10⁻¹², the
fields inside a charged sphere and a current cylinder, a plane wave with
∇F = 0 converging at fourth order, and two frame tests. Under a parity
inversion E flips and B does not, because B is an axial vector, carried
in the algebra as the bivector I·B; the vector rule, the algebra and a
recomputed Biot–Savart field agree to 10⁻¹² for inversion, a mirror and
fifty random rotations. A moving probe sees the Lorentz-transformed
field, obtained in the algebra by a rotor acting on F, and that equals
the textbook transformation for speeds to 0.9c. The point of the rule is
stated in CLAUDE.md §3c: using Maxwell's notation does not change any
prediction, and no "scalar wave" or other extra physics is to be smuggled
in through a choice of algebra.

## 5.2 From a cross-section to a circuit: the microstrip pair

Two parallel traces on a 1.6 mm FR-4 board, 10 cm long, are a circuit
whose values nobody looks up. `pcb_layout` solves the cross-section as
div(ε grad φ) = 0 with each conductor driven in turn, which gives the
2×2 capacitance matrix C; the same solve with ε = 1 gives C₀, and the
inductance matrix is L = μ₀ε₀C₀⁻¹ (`em.Section`, milestone M9i). The
pair then becomes twenty sections of the coupled telegrapher's equations
(`em.CoupledLadder`): per section a self-inductance on each trace, a K
element for the mutual, a capacitor from each trace to ground and one
between the traces. The export folder holds that netlist for a 3 mm
trace and a 2 mm gap. Its first section reads:

```
LX0_a0 a_in x0_m0 7.240277495153368e-10
LX1_a0 v_in x1_m0 7.240277495153364e-10
KX01_a0 LX0_a0 LX1_a0 0.1276995674610439
...
CX0_g0 x0_m0 0 5.942475582044419e-13
CX01_0 x0_m0 x1_m0 3.947220550682456e-14
CX1_g0 x1_m0 0 5.942475582044421e-13
```

The aggressor is driven through its own 47.8 Ω and both traces are
terminated in it at both ends. Run the netlist with `-json`; its single
AC point is 100 MHz:

| Node | Magnitude |
|---|---|
| a_in (aggressor, near end) | 0.5002 V |
| a_out (aggressor, far end) | 0.4997 V |
| v_in (victim, near end) | 0.01764 V |
| v_out (victim, far end) | 0.00622 V |

Near-end crosstalk is v_in over a_in: −29.1 dB at 100 MHz, which is the
value the example's test records for a 2 mm gap. The AC ledger's worst
relative residual over the sweep is 0, and the results sheet notes that
the circuit torn automatically into two pieces joined at the source
resistor agrees with the full solve to 1.2×10⁻¹⁵ with the same power.

What the Go program adds is the extraction and the search. The extraction
is validated in `em/extract_test.go` against closed forms:

| Case | Reference | Result |
|---|---|---|
| Coax, εr = 2.1 | 2πε/ln(b/a) | +0.07 % |
| Stripline, zero thickness | Cohn 1954 | −0.006 % after Richardson extrapolation |
| Two-wire line | πε₀/acosh(D/2a) | 0.1 % |
| Microstrip, FR-4, 3 mm | Hammerstad–Jensen 1980 | Z₀ −1.8 %, ε_eff −0.4 % |

The genome is the trace width (0.3 to 6 mm) and the gap (0.2 to 10 mm),
with a weighted sum of |Z₀ − 50|/50, crosstalk above −30 dB and a small
weight on board width. The default run, population 20 for 12
generations with seed 1, reports after 218 evaluations:

```
w 2.766 mm, s 3.020 mm, board 8.55 mm
Z₀ 49.82 Ω, ε_eff 3.257, L₁₁ 301 nH/m, C₁₁ 121.3 pF/m, K = L₁₂/L₁₁ 0.0857, NEXT -33.1 dB
```

Coupling falls with the gap, as the test checks for a 3 mm trace: K is
0.287 at 0.5 mm, 0.128 at 2 mm and 0.030 at 6 mm, with near-end
crosstalk of −20.8, −29.1 and −43.3 dB. The limits are declared with the
result: zero-thickness strips, the quasi-static TEM solution only, so
dispersion and radiation are left out, and the ladder's section count set
by the λ/10 rule.

That rule is the validity check of `em/limits.go`. Kirchhoff's laws
follow from Maxwell's when every part of a circuit is small against the
shortest wavelength present; a conductor of length d at frequency f has a
phase delay 2πd/λ that a lumped model sets to zero. `em.CheckQuasiStatic`
compares declared lengths with the wavelength at the analysis's highest
frequency, and for a transient takes that frequency from the fastest edge
as the knee 0.5/t_r (Johnson and Graham). Past λ/10 it returns
`ErrNotQuasiStatic`, and the genetic algorithm penalizes such a
candidate as invalid; the coil examples of chapter 6 use the same check.

## 5.3 Does a rotating electric field make a magnetic field?

The homopolar dynamo of chapter 3 turns a conductor in B and gets a
voltage. `rotating_efield` asks the inverse question and keeps two
answers apart. The real physics first: B comes from moving charge
(Rowland 1876) and from a changing electric field, the displacement
current in Maxwell's 4π ℭ = V.∇ℌ. Five cases are computed and checked
with the quaternion operators of §5.1; the program prints:

| Case | Result | Expected |
|---|---|---|
| 1a. Rowland disk, σ = 10⁻⁵ C/m², 10 cm, 20,000 rpm | B_z 1.316×10⁻⁹ T on the axis, 4.496×10⁻¹⁰ T at 5 cm up | analytic (μ₀σω/2)[(R² + 2z²)/√(R² + z²) − 2\|z\|], same values |
| 1b. Uniform 1 MV/m rotating at 1 MHz | \|B\| = 1.75×10⁻⁶ T at r = 5 cm | (μ₀ε₀/2)·E·ω·r; Ampère–Maxwell residual 3.3×10⁻¹⁰ |
| 1c. 1 µC ring, 5 cm, at rest and at 1000 rev/s | B at the centre 0 T at rest, 1.26×10⁻⁸ T spinning | the same static E pattern in both; B follows the moving charge |
| 1d. ±1 nC on a 10 cm rod at 1000 rev/s | \|B\| = 8.99×10⁻¹² T on the axis 2.5 cm up | the E pattern really rotates; residual 1.3×10⁻⁴ of the displacement term |
| 1e. εr = 3 disk, 10 cm, 1 mm, 50 rev/s in a fixed 1 MV/m | B_r = −3.46×10⁻¹⁰ T at r = 5 cm | thin-disk −(εr − 1)·E·v/c² = −3.5×10⁻¹⁰ T |

The answer is that a rotating E pattern by itself makes no B; B comes
from the charge that moves, or from ∂E/∂t where the field really changes
in time. Case 1e is the electric mirror of the homopolar disk: a
dielectric moving in E is magnetized, M = P × v (Röntgen 1888,
Eichenwald 1903). The ring sum uses only the bound face charge
ε₀(εr − 1)E carried round as two opposite surface currents and
Biot–Savart's μ₀, so the 1/c² = μ₀ε₀ comes out by itself; its deficit
against the thin-disk value falls linearly in thickness over radius, and
the test checks that halving the thickness halves it. That 1/c² is why
the mirror effect is nanotesla where the homopolar one is volts: 1 T at
100 m/s gives 100 V/m, while 1 MV/m at 100 m/s gives about 1 nT.

The second answer is hypothetical and labelled so. With
`em.Hypothetical{MagneticCharge: true}` the program computes the force a
magnetic charge of 1 A·m would feel, the dual pair of Joot's Theorem 3.16
(CLAUDE.md §3f, register entry H2): 1.32×10⁻⁹ N at rest at the Rowland
disk's centre, 1.11×10⁻⁹ N moving at 100 m/s through 1 MV/m. No magnetic
monopole has ever been observed, so these go to the hypothetical column
of the momentum ledger only, and the real net force stays zero. The
same output line that gives the number says so.

## 5.4 Rates of change: in the material and in the equations

Where do second, third and higher time derivatives enter electromagnetism?
`rate_of_change` was written for a premise recorded in
[doc/vassilatos-references.md](../doc/vassilatos-references.md): that a
fast pulse acts on space the way a fast slap acts on water, which answers
like a solid above some rate. The engine's answer lives in two places.

**In the material (standard physics).** Maxwell's constitutive relation
𝔇 = (1/4π)K𝔈 with a permittivity that depends on how fast the field
changes. Water's molecules need about 8 ps to turn with a field, so a
single Debye term, τ·dP/dt + P = ε₀Δε·E, describes it below about
100 GHz, with ε(0) = 78.6, ε∞ = 5.4 and τ = 8.3 ps (`em/dispersion.go`;
Debye 1929). Through the network analogy of CLAUDE.md §3a.1 each such
term is exactly a circuit: a plate cell filled with the material is
C∞ = ε₀ε∞A/d in parallel with a series R–C per Debye term. The export
folder holds that netlist for a 1 cm², 1 mm cell:

```
plate cell filled with liquid water, 25 °C (single Debye)
V1 top 0 DC 0 AC 1 0
CINF top 0 4.78126141891221e-12
RD0 top d0 0.128061410644749
CD0 d0 0 6.48126547896988e-11
.ac dec 10 1meg 100g
.end
```

Run it. The report shows the Debye node falling to −14.5 dB at −79° by
100 GHz, and the ledger's worst relative residual is 0. Reading the
source current from the JSON document and dividing its admittance by
jω·ε₀A/d gives the permittivity the cell presents:

| Frequency | ε′ | ε″ |
|---|---|---|
| 1 MHz | 78.60 | 0.00 |
| 100 MHz | 78.60 | 0.38 |
| 1 GHz | 78.40 | 3.81 |
| 10 GHz | 62.95 | 30.01 |
| 100 GHz | 8.00 | 13.54 |

which is the Debye function to the digits shown; the example's test
holds the worst difference over the sweep at 7.9×10⁻¹³. The same model
in the time domain, with an exponential integrator that is exact for a
field varying linearly between grid points, gives the slap test: a field
ramped from 0 to 1 in 1 µs is answered with ε = 78.6, in 100 ps with
72.5, in 10 ps with 36.1, in 1 ps with 9.6, and in every case 78.6 once
the ramp has settled for twenty time constants. Slow, the molecules
follow; fast, they cannot turn in time and the water answers like a stiff
solid with ε → 5.4. So the slap-on-water picture is real, and it is Debye
relaxation. Causality is checked too: the Kramers–Kronig integral of ε″
recovers ε′(0) − ε∞ = 73.2. The results sheet tears the cell into two
pieces joined by the Debye branch, and the torn and full solves agree to
2.2×10⁻¹⁶ with the same power.

**In the field equations (hypothesis H19).** Bopp (1940) and Podolsky
(1942) add a term to Maxwell's Lagrangian that makes the field equation
fourth order, (1 + a²□)□A = μ₀J, with one new length a; a = 0 is
Maxwell. `em/podolsky.go` holds the static point-charge potential
φ = q/(4πε₀r)·(1 − e^(−r/a)), finite at the centre, the field, the
self-energy and the free modes, and the example prints each beside the
standard value, as the register requires. For an electron with a = 1 fm,
chosen for illustration:

| r / a | Maxwell φ | Podolsky φ | Fraction below Coulomb |
|---|---|---|---|
| 0 | infinite | 1.440 MV | 1 |
| 1 | 1.440 MV | 0.910 MV | 0.368 |
| 3 | 480.0 kV | 456.1 kV | 0.050 |
| 10 | 144.0 kV | 144.0 kV | 4.5×10⁻⁵ |

The field energy of a point charge becomes finite, q²/(8πε₀a), and
equals m_e c² when a = 1.409×10⁻¹⁵ m, which was Podolsky's motivation.
The free field gains a massive mode starting at ω = c/a, 4.8×10²² Hz for
this a. A Coulomb-law test good to one part in 10⁶ at 1 cm bounds only
a < 0.72 mm; the Bhabha-scattering bound of Bufalo, Pimentel and Soto
(2014) is far tighter and is to be entered from the paper, not assumed.
The effect is gone a few a away, so only very short-range experiments
can see it. Status: untested, with a = 0 beside every line.

## 5.5 Where the lumped picture ends: the field solver

When a circuit is electrically large, or when the question is about the
field itself, `em/fdtd.go` solves Maxwell's two curl equations on a Yee
grid: E on cell edges, B on cell faces half a cell and half a step apart,
conducting walls, a graded matched absorbing layer whose electric and
magnetic loss rates are equal so it has the impedance of vacuum,
materials by εr and σ, and current sources on edges. Because the discrete
divergence of a discrete curl is exactly zero, div B stays at round-off,
and the scheme satisfies an energy ledger exactly, step by step: source
work in equals field energy plus heat in materials plus energy absorbed
by the layer, with a residual of about 10⁻¹³. The solver is validated
against the TM110 cavity mode (Yee's frequency to 10⁻⁷), a Hertzian
dipole against Larmor's formula (+0.5 %) and the exact Hellwarth–Nouchi
flying-doughnut pulse (`em/doughnut.go`, 3.8 % at 8 cells per q1). It
is the instrument of chapter 9, where it tests whether light on a double
loop stays confined, and it is the local field solve that §3g's
arc-surge experiments call for.

## 5.6 What to keep

- The quaternion, spacetime-algebra and vector forms of Maxwell's
  equations are one physics. Every field computation is checked three
  ways, and no notation adds a prediction.
- Layout is a circuit: R, L, C and K come from a field solve of the
  geometry, validated against closed forms, so a genome can carry
  widths and gaps.
- The lumped model has a limit, λ/10 of the shortest wavelength, and the
  engine refuses or penalizes a circuit past it.
- B comes from moving charge or a changing E, never from a static
  pattern that merely rotates; the mirror of the homopolar disk is
  nanotesla because of μ₀ε₀.
- Rate dependence in a material is a circuit (Debye as R–C), and a
  hypothesis that changes the field equations is printed beside Maxwell
  with its bound.

## Run it yourself

```
go run ./cmd/kronspice -json examples/pcb_layout/export/pcb_layout_pair_w3_s2.cir
go run ./examples/pcb_layout                      # the GA, about a minute
go run ./examples/rotating_efield
go run ./cmd/kronspice examples/rate_of_change/export/water_cell_debye.cir
go run ./examples/rate_of_change
go test ./em                                      # triple twins, extraction, FDTD, dispersion, Podolsky
```

## Sources

- J. C. Maxwell, *A Treatise on Electricity and Magnetism*, 1873,
  Vol. II, Art. 618–619 (field equations in quaternion form), Art. 701
  (coaxial circles).
- D. Hestenes, *Space-Time Algebra*, 1966, for ∇F = J; C. Doran,
  A. Lasenby, *Geometric Algebra for Physicists*, 2003.
- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed.,
  2008; E. Hammerstad, Ø. Jensen, "Accurate models for microstrip
  computer-aided design," *IEEE MTT-S*, 1980; S. B. Cohn, stripline,
  1954.
- D. M. Pozar, *Microwave Engineering*, 4th ed., 2012, §2.1;
  H. W. Johnson, M. Graham, *High-Speed Digital Design*, 1993, §1.1 and
  §4.1, for the lumped limit and the knee frequency.
- H. A. Rowland, *Am. J. Sci.* 15, 30, 1878; W. C. Röntgen, *Ann. Phys.*
  35, 264, 1888; A. Eichenwald, *Ann. Phys.* 11, 1, 1903; J. D. Jackson,
  *Classical Electrodynamics*, 3rd ed., §11.10.
- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16
  (reference only; the units erratum is recorded in CLAUDE.md §3f).
- P. Debye, *Polar Molecules*, 1929; H. A. Lorentz, *The Theory of
  Electrons*, 1909; A. Taflove, S. C. Hagness, *Computational
  Electrodynamics*, 3rd ed., 2005, ch. 3 and ch. 9; K. S. Yee, *IEEE
  Trans. Antennas Propag.* 14, 302, 1966.
- D. C. Elton, arXiv:1704.01667; J. Barthel et al., *Chem. Phys. Lett.*
  165, 369, 1990; R. Buchner, J. Barthel, J. Stauber, *Chem. Phys. Lett.*
  306, 57, 1999, for water's Debye parameters.
- F. Bopp, *Ann. Phys.* 430, 345, 1940; B. Podolsky, *Phys. Rev.* 62,
  68, 1942; R. Bufalo, B. M. Pimentel, D. E. Soto, *Phys. Rev. D* 90,
  085012, 2014.
- R. W. Hellwarth, P. Nouchi, *Phys. Rev. E*, 1996, the flying-doughnut
  pulse.
- [doc/M9h.md](../doc/M9h.md), [doc/M9i.md](../doc/M9i.md),
  [doc/M9n.md](../doc/M9n.md).
