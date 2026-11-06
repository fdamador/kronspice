# Kron's field and polyhedron networks: figure catalogue

A transcription of every drawing and table in the sources we hold that
bears on Kron's equivalent circuits of Maxwell's equations and on his
multidimensional (polyhedron) networks, read from page images rendered with
PyMuPDF at 110–130 dpi (2026-11-01), not from text layers. It is the
specification for register entry H20 (CLAUDE.md §1) and for the planned
`kron/polyhedron.go`. Element formulas and numbers are transcribed; the
drawings themselves stay in `reference/` (copyright), cited by page.

Sources and page offsets (printed page → PDF page):

| Tag | Source | Offset |
|---|---|---|
| K62 | G. Kron, "Multidimensional curve-fitting with self-organizing automata," *J. Math. Anal. Appl.* 5, 46–69, 1962. `reference/PII_ 0022-247X(62)90005-7.pdf` | printed p = PDF p − 45 |
| K45a | G. Kron, "Numerical solution of ordinary and partial differential equations by means of equivalent circuits," *J. Appl. Phys.* 16, 172–186, 1945. `reference/Gabriel Kron_ Num. Sol'n of ODE's & PDE's …pdf` (website transcription, 32 pp.) | figures only; cited by PDF page |
| K45b | G. Kron, "Electric circuit models of the Schrödinger equation," *Phys. Rev.* 67, 39–43, 1945. `reference/Gabriel Kron, Electric Circuit Model of the Schrödinger Equation …pdf` (11 pp.) | cited by PDF page |
| B68 | N. V. Balasubramanian, *Differential structures in electromagnetic field and circuit theory*, PhD thesis, Liverpool, 1968. `reference/642931.pdf` | printed p = PDF p − 8 |
| N74 | H. Nicholson, "System concepts in Kron's polyhedron model and the scattering problem," ACSE Report 26, Sheffield, 1974. `reference/acse research report 26.pdf` | printed p = PDF p − 2 |
| T01 | E. Tonti, "Finite formulation of the electromagnetic field," *PIER* 32, 1–44, 2001. `reference/tonti_2001_pier32_finite_formulation.pdf` | printed = PDF |

---

## 1. K62 Fig. 1 (p. 49): the three-dimensional Maxwell network

"Linear electric-network model of the field equations of Maxwell
(orthogonal, curvilinear reference frame)." One cell of the lattice in
coordinates u¹, u², u³ with scale factors h₁, h₂, h₃ and cell sides Δu¹,
Δu², Δu³; p = d/dt; Gaussian-type units with c written explicitly. Every
element is labelled on the drawing:

| Element | Along axis 1 | Axis 2 | Axis 3 |
|---|---|---|---|
| Dielectric branch (capacitor) | Z = (h₁/h₂h₃)(Δu¹/Δu²Δu³)(c/εp) | (h₂/h₁h₃)(Δu²/Δu¹Δu³)(c/εp) | (h₃/h₁h₂)(Δu³/Δu¹Δu²)(c/εp) |
| Conductive branch (resistor, in parallel with the capacitor) | Z = (h₁/h₂h₃)(Δu¹/Δu²Δu³)(c/σ) | (h₂/h₁h₃)(Δu²/Δu¹Δu³)(c/σ) | (h₃/h₁h₂)(Δu³/Δu¹Δu²)(c/σ) |
| Magnetic branch (inductor) | Z = (h₂h₃/h₁)(Δu²Δu³/Δu¹)(μp/4c) | (h₁h₃/h₂)(Δu¹Δu³/Δu²)(μp/4c) | (h₁h₂/h₃)(Δu¹Δu²/Δu³)(μp/4c) |
| Ideal transformers (the ring of coils drawn at every node) | L = M = ∞ | | |

Variables on the drawing: mesh currents H₍₁₎h₁Δu¹, H₍₂₎h₂Δu², H₍₃₎h₃Δu³
(line integrals of the magnetic force round the magnetic meshes); branch
voltages E₍₁₎h₁Δu¹, E₍₂₎h₂Δu², E₍₃₎h₃Δu³ (line integrals of the electric
force along the dielectric branches); nodal charge ρh₁h₂h₃Δu¹Δu²Δu³.

The 1/4 in the inductors: each magnetic mesh (a "large" mesh, B68 p. 155)
carries its flux ∂B/∂t·ΔzΔx as four voltage drops, one quarter each,
across the four inductors round it (B68 p. 157). SI form of the same
values (B68 eqs. 5-3, 5-4 and p. 156):

```
C₁ = (h₂h₃/h₁)(ΔyΔz/Δx)·ε        G₁ = (h₂h₃/h₁)(ΔyΔz/Δx)·σ
L₂ = ¼ (h₃h₁/h₂)(ΔzΔx/Δy)·μ        (and cyclically)
```

What the ideal transformers do (B68 pp. 153–157, Figs. 10–11): each
"large" mesh ABCDEFGH holds impedance elements plus the two primaries of
two ideal transformers and their secondaries placed diagonally opposite
in the same mesh, poled so that round the large mesh each primary voltage
cancels its secondary; the "small" meshes (EFKL) of inductors and
transformers only carry zero mesh current because the transformers offer
infinite impedance to them. Hence the mesh currents are of the type H₂Δy,
the voltage round a large mesh sums to ∂B²/∂t·ΔzΔx = 0 (Faraday), and
`*(J + Ḋ) = dH` (Ampère) is the branch–mesh relation i = C_c·i′. The
transformers at the boundary are energised from outside the region, which
is how a scalar magnetic potential from external sources enters (p. 158).
div B = 0 follows from the transformer flux linkages (B68 App. 2).

**Why this matters.** Fig. 1 is the 1944 *Proc. IRE* network complete,
with every value. It settles H20 step 1 at k = 3 without the 1944 paper.

## 2. K62 Fig. 2 (p. 50): primal and dual polyhedra

Axes X₍₁₎, X₍₂₎, X₍₃₎; the primal polyhedron drawn solid with its 0-space
(vertices), 1-space (edges), 2-space (a hatched triangle), 3-space (a
hatched tetrahedron); the dual polyhedron dashed, its vertices inside the
primal cells. Text (p. 54): all even-dimensional q-networks are magnetic
(inductors and resistances), all odd-dimensional ones dielectric
(capacitors and conductances); to introduce ideal transformers a dual
polyhedron is needed in which the roles of even and odd are exchanged; the
two are interconnected both conductively and inductively; one complete
wave (e, h, b, d) occupies the configuration of Fig. 1 spread over a q- and
a (q+1)-network, so a k-dimensional polyhedron carries k/2 full waves.

This is T01 Fig. 4 (p. 12) in Tonti's language: electric potential V on
primal points, voltage U on primal lines, magnetic flux Φ on primal
surfaces, magnetic charge content on primal volumes; electric charge Qᶜ on
dual volumes, electric flux Ψ on dual surfaces, magnetic voltage U_m on
dual lines, magnetic potential V_m on dual points. "A single cell complex
is not sufficient but it is necessary to introduce a dual complex"
(T01 p. 12).

## 3. K62 Fig. 3 (p. 56) and the regression algorithm (pp. 54–55)

The worked example: four points on a plane, data matrix

```
X = [ x  y ]        point 1: ( 1,  1)
    [ 1  1 ]        point 2: ( 1, −2)
    [ 1 −2 ]        point 3: (−2,  3)
    [−2  3 ]        point 4: ( 4,  2)
    [ 4  2 ]
```

drawn as a polyhedron on the (x, y) grid with its networks labelled:
"−1 network (electric), least-square model"; "0-network (magnetic), 4
junction points"; "1-network (dielectric), 6 branches"; "2-network
(magnetic), 3 planes" (the three hatched triangles 1-2-3, 1-3-4, 1-2-4
with the inductor symbols on their sides). Counts: 4 + 6 + 3 = 13
characteristics, 9 of them derived (p. 57); one residual degree of freedom
(n − k − 1 = 1), a 12 × 12 moment matrix (p. 57).

Six trial functions (p. 56):

```
z₁ = 3x + 4y + 2                  z₄ = x³ + 2y²x − 2x² − 2xy + 4y + 8
z₂ = 5xy + 3x + 4y + 2            z₅ = 10[1 + log(x + 4) − 1/y]
z₃ = x² + 2y² − 3xy + 5x − 3y + 6 z₆ = 10[2 + e^(x/10) − y·log(x + 10)]
```

The eight steps (pp. 54–55), in Kron's words condensed:

1. The dependent variables y are known currents I impressed on the
   vertices of the untorn polyhedron (vertices belong to the 0-network).
2. The polyhedron is torn into isolated q-networks and the open-circuit
   voltages E′ on each q-network are the "divided differences" to be
   fitted along with the function; they extend the y vector.
3. The open-circuit voltages are impressed as e = Xₜy on the
   "intersection" polyhedron whose impedance matrix is XₜX = S.
4. Invert it: (XₜX)⁻¹ = S⁻¹.
5. The resultant currents i are the regression coefficients,
   β̂ = (XₜX)⁻¹Xₜy.
6. The additional currents I′ on the isolated subdivisions are ŷ = Xβ̂,
   the estimated function and all its divided differences.
7. E′ − I′ = y − ŷ = ε is the error.
8. The unknown potentials E on the vertices of the untorn polyhedron are
   not used.

"Thus the physical tearing apart of the original polyhedron is an
absolutely necessary procedure, since it is the hypothetical 'intersection'
polyhedron that acts as an estimating model" (p. 55). The curve fitting
does not reduce the degrees of freedom: X is enlarged by an equal number of
rows and columns (p. 55).

### Table I (p. 58), z = 3x + 4y + 2, transcribed

Columns: the three functions to be fitted (given) · least-square ·
least-square + linear · least-square + linear + planar (frame 1) · the
same (frame 2, "two nonoscillatory polyhedra") · oscillatory polyhedron.

Function at the four points:

| given | LS | LS+lin | LS+lin+pl (1) | LS+lin+pl (2) | oscillatory |
|---|---|---|---|---|---|
| 9.0000000 | 9.0000002 | 9.0796617 | 8.8787335 | 8.9632293 | 8.9999997 |
| −3.0000000 | −2.9999995 | −3.0265549 | −2.9595791 | −2.9877437 | −2.9999993 |
| 8.0000000 | 8.0000000 | 7.9734455 | 8.0404218 | 8.0122575 | 8.0000011 |
| 22.000000 | 22.000000 | 21.973445 | 22.040421 | 22.012257 | 22.000001 |

First divided differences (six branches):

| given | LS+lin | LS+lin+pl (1) | LS+lin+pl (2) | oscillatory |
|---|---|---|---|---|
| 4.0000000 | 5.9057718 | 4.0538963 | 4.0649511 | 3.9999975 |
| 1.8864844 | −1.8864846 | 1.8795529 | 2.0501086 | 1.8864869 |
| 5.0000000 | 5.0000005 | 4.9757475 | 4.8752506 | 4.9999983 |
| −0.27735011 | −6.4341227 | −0.32219431 | −0.41067224 | −0.27734973 |
| 4.1109605 | 6.5376262 | 4.0598305 | 4.0416890 | 4.1109617 |
| 2.3015856 | −2.3015861 | 2.3148763 | 2.4084983 | 2.3015888 |

Second divided differences (three planes):

| given | LS+lin+pl (1) | LS+lin+pl (2) | oscillatory |
|---|---|---|---|
| −4.5000000 | −4.4910178 | −4.3879722 | −4.5000036 |
| −1.4558823 | −1.5007952 | −1.5821587 | −1.4558866 |
| −1.7027027 | −1.7206683 | −1.7817965 | −1.7027075 |

Tables II–VI (pp. 59–64) have the same layout for z₂ … z₆, transcribed in
`examples/kron_polyhedron/tables.go`. The "LS + linear" fourth entry above is
as printed; §11 gives the likely misprint.

### The oscillatory polyhedron (pp. 65–67)

To fit any number of functions without adjustment, "the polyhedron
itself (without the regression network or −1 network) should be made
oscillatory, or self-excited. Such a stage could be expected to be
reached by enlarging the purely electromagnetic wave (b, d, e, h) with
electric and magnetic charges and currents (ρᵉ, Jᵉ, ρₘ, Jₘ). This next
step in the generalization was suggested by the unused portions in the
configuration of Fig. 1, which still had just enough room left to
incorporate the added electrical source variables" (p. 65). In the
oscillatory runs "the previous linear network configuration was retained
… and only the planar network configuration was revamped" (p. 65); "no
exact fit is possible, since the polyhedral model still has one residual
degree of freedom" (p. 65). The choice of open-circuit and short-circuit
reference frames of minimum stored energy "is basically a nonlinear
programming problem … arrived at by physical reasoning … not been reduced
as yet to a routine computer procedure" (p. 67).

So the oscillatory frame is a search over the 2-network's reference frame
(which planes are open, which short, how they interconnect) with the
objective Kron states (minimum stored energy = minimum error-square,
p. 67) and Tables I–VI as the check. That is the GA step of the plan below.

## 4. K45a: the 1945 J. Appl. Phys. figures

- **Figs. 3–7 (PDF pp. 4–8), time networks.** The nonlinear wave equation
  ∂/∂u¹(f₁∂φ/∂u¹) + ∂/∂u²(f₂∂φ/∂u²) + ∂/∂u³(f₃∂φ/∂u³) = f₄ + f₅φ + f₆∂φ/∂t
  as a lattice in which time is an extra spatial axis: layers Δt apart,
  joined by impedances. Element labels: f₁Δu²Δu³Δt/Δu¹, f₂Δu¹Δu³Δt/Δu²,
  f₃Δu¹Δu²Δt/Δu³ (series, along the space axes), f₄Δu¹Δu²Δu³Δt (source),
  f₅Δu¹Δu²Δu³Δt (shunt), f₆Δu¹Δu²Δu³ (the time link, Fig. 3); Fig. 4 the
  three-dimensional time network (two space axes plus t, layers Δt/2
  apart); Fig. 5a two-dimensional (one space axis plus t), 5b
  one-dimensional (t alone: f₁∂φ/∂t + f₂φ + f₃ = 0 as a ladder of Δt/2
  sections); Figs. 6–7 the same when time enters as a second derivative
  (inductors in the time links). Curvilinear coordinates: f₁ = (h₂h₃/h₁)f₁′,
  f₂ = (h₁h₃/h₂)f₂′, f₃ = (h₁h₂/h₃)f₃′, f₄ = h₁h₂h₃f₄′, etc. (PDF p. 6).
  Note "the algebraic symbols may be replaced by a mechanical model (masses
  and springs) or an acoustical model (pipes and chambers)" (PDF p. 4).
- **Figs. 9–10 (PDF pp. 11, 13).** Worked step-by-step solutions:
  di/dt = 10 − 2i with Δt = 0.1 (Fig. 9: values 0.475, 0.905, 1.294, 1.646,
  1.965, 2.253, 2.514, 2.75, 2.964 against Doherty & Keller's 0.91, 1.65,
  2.26, 2.76, …), and a transient heat-flow lattice (Fig. 10,
  Δt = 0.0005/a, y = Δx = 0.1, y = aΔt/Δx = 0.005, columns 10.00 / 181.55 /
  19.61 / 1.272 / 0.056 at t = 0.0025/a …). Both are reproduction targets
  for a Kron time network solved by the engine.
- **Fig. 13 (PDF p. 17), a junction**, and the method of unbalanced
  currents and voltages (eqs. 1–14, PDF pp. 17–19): i_u = Σ(e_k − e₀)y_k + I,
  weighted average E = (Σe_k y_k + I)/Σy_k, e_u = E − e₀ = i_u/Σy; for
  Laplace's equation y_k = 1: i_u = e₁+e₂+e₃+e₄ − 4e₀; for Maxwell's
  equation y_k = y_L, y_g = y_c (a capacitance, negative): the 4 becomes
  k = 4 + y_c/y_L. Iterate: assume potentials, compute i_u and e_u at every
  junction, correct, repeat (PDF p. 19).
- **Fig. 14 (PDF p. 19), the square waveguide / cavity.** "Let a square
  wave guide or cavity resonator, by a conducting surface that extends to
  infinity in the y direction, be divided into sixteen squares. The net of
  coils in each square consists of two inductances (z = 0.8 or y = 1.25)
  and a ground capacitor (z = 1.3 or y = 0.77). The currents in the
  inductance represent E_x and E_z, while the voltages across the
  capacitors represent H_y." Bounded in the x–z plane on three sides; at
  the open (upper) end a sinusoidal E_x is impressed: coil currents 50,
  35.35, 35.35, 50 across the top. Junction potentials H_y found on the
  network analyzer (rows, top to bottom, left column then across):
  84.8 | 53.68, 0.24, −53.28, −84; 25.66 | 16.34, 0.24, −15.46, −25;
  −44.72 | −27.86, 0.44, 28.46, 46.2; −100.2 | −60.7, 0.48, 63.4, 101.8;
  −130.2 | −80.2, 0.7, 83.7, 133.7.
- **Fig. 15 (PDF p. 20).** The unbalanced currents and voltages that
  result (i_u, e_u per junction; largest i_u 7.93 A against a largest coil
  current of 102.8 A, "8 percent inaccuracy"). Rows: (2.4, 1.39), (2.15,
  0.72), (−0.02, −0.01), (−0.75, −0.25), (−2.2, −1.27); (2.12, 0.5),
  (−4.36, −1.03), (0.63, 0.15), (3.41, 0.81), (7.93, 2.66); (5.26, 1.24),
  (6.97, 1.65), (−0.21, 0.05), (2.15, −0.51), (−6.1, −2.04); (4.1, 0.97),
  (−2.87, −0.68), (2.69, 0.64), (−0.1, −0.02), (0.9, 0.3); (−0.2, −0.12),
  (1.17, 0.39), (2.97, 1.0), (−2.2, −0.74), (0.4, 0.23).
- **Fig. 16 (PDF p. 20).** An elastic beam on a triangular lattice with
  displacements and unbalanced currents (0, 7.676, 19.375, 34.856, 52.865,
  72.154 along the top; bracketed residuals).

Fig. 14's cavity is the first twin target: the same 4 × 4 cell network
built from Fig. 1's values, driven by the same impressed E_x, must give
Kron's H_y potentials to the accuracy his analyzer allowed (the Fig. 15
residuals), and its modes must match `em/fdtd.go`'s TM₁₁₀ cavity test.

## 5. K45b: the Schrödinger networks

Figs. 1–7 (PDF pp. 3–10): kinetic-energy operator as series coils of
impedance (2m/ħ²)Δx; potential energy as shunt coils of admittance VΔx;
total energy −E as shunt capacitors of admittance EΔx (Fig. 1); the sum of
operators as an interconnection (Fig. 2); ψ as junction voltages, αψ as
currents (Fig. 3); the potential well with the generator inserted and E
varied until the generator current is zero (Fig. 4); the variable-frequency
model L₁ = (2m/ħ²)Δx, L₂ = 1/(VΔx), C = ħΔx (Fig. 5); the free particle as a
transmission line (Fig. 6); the three-dimensional curvilinear model with
y = (ħ²/2m)(h₂h₃/h₁)(Δu²Δu³/Δu¹) series admittances and
V·h₁h₂h₃Δu¹Δu²Δu³, E·h₁h₂h₃Δu¹Δu²Δu³ shunts (Fig. 7). Already implemented
in `examples/kron_schrodinger`; listed here because Fig. 7 is the
same cell as K62 Fig. 1 with the magnetic branches removed.

## 6. B68: Balasubramanian's figures

- **Fig. 5a (p. 124), interconnected network model for fields.** One
  layer of the lattice: nodes M, N, P, Q, R, S, T; each branch an R ∥ C
  pair in series with a generator ∂A_k/∂t·Δx_k; branch currents
  (J + ∂D/∂t)·ΔyΔz; mesh currents i²Δy, i³Δz; meshes A and B with fluxes
  B²′ΔzΔx; E₁Δx, E₂Δy across the pairs. **Fig. 5b (p. 125), the primitive
  1-network**: a branch with impressed V₁ (= grad ψ) and generator v₁ (=
  ∂A/∂t), response currents I′ and i′, "E = −∂A/∂t + grad ψ"; the R ∥ C
  branch carries J + ∂D/∂t. Element values (eqs. 5-3, 5-4):
  C = (h₂h₃/h₁)(ΔyΔz/Δx)·κ, G = (h₂h₃/h₁)(ΔyΔz/Δx)·σ for branch PQ.
- **Fig. 6 (p. 127), elemental blocks and regions**: the area Σ normal to
  a branch over which the 2-form *J is integrated; the block with B¹′ΔyΔz,
  B²′ΔzΔx, B³′ΔxΔy on its faces.
- **Fig. 7 (p. 139), algebraic diagram for the field network.** (a) the
  1-circuit (branches: A₀₁, C_c¹, currents i¹, voltages V₁ + v₁) joined to
  the 2-circuit (planes: A₀₂, C_c², currents I² + i², voltage V₂) through
  the mesh/open-path square (i^{c1′}, I^{o2′}, v_{c1}′, V₀₂′). (b) the same
  with field quantities: J + Ḋ = *di (Ampère), E = dψ − ∂A/∂t,
  −*∂B/∂t (Faraday), *H, Ω; the impedances 1/(σ + jωκ) between E and
  J + Ḋ, and μη₀ between *H and *B; d*(J + Ḋ) = 0 and d*B = 0 at the ends;
  the 4-dimensional tensors H^{αβ} and F_{αβ} as the two "currents" and
  "voltages".
- **Fig. 8 (p. 150), extended algebraic diagram for the field network.**
  Four columns, Nodes · Branches · Planes · Volume, joined by the incidence
  matrices M⁰₁, M¹₂, M²₃ (and their transposes along the bottom), each
  column with its C_c, A₀ pair, the impedances 1/(σ + jωκε₀) on the
  branches and jωμη₀ on the planes; "*d*" (generalised divergence) and
  "d" (generalised curl) marked on the arrows. Voltages form a
  cohomology sequence, currents a homology sequence (p. 151).
- **Fig. 9 (p. 152), "Kron's algebraic diagram for space filters",
  reproduced from Kron's 1958/59 space-filter papers (B68 ref. 6).** Four
  columns Junction · Branch · Plane · Cube, each an orthogonal-network
  diamond, with the material tensor of each column written in its middle:
  ε¹¹ on the branch column, μ₂₂ on the plane column, ε³³ on the cube
  column; the operators between columns: DIV = (M⁰₁)ₜ from junction to
  branch, CURL⁻¹ = (M¹₂)⁻¹ from branch to plane, DIV = (M²₃)ₜ from plane
  to cube; the variables at top and bottom: ḃ₀, Ḃ₀ at the junctions; Ḋ¹,
  ḋ¹ and E₁, e₁ at the branches; h², H² and ḃ₂, Ḃ₂ at the planes; Ḋ³, ḋ³
  and E₃, e₃ at the cubes. "Kron's analysis, however, includes abstract
  quantities like magnetic poles and magnetic conduction currents"
  (p. 153).

  **This is the q-network rule I reported missing.** Each q-network is an
  orthogonal network (Nicholson's eq. 3 below) whose primitive impedance
  tensor is the permittivity on odd q and the permeability on even q;
  neighbouring networks are joined by the incidence matrices M^q_{q+1}
  and their transposes, playing div and curl⁻¹ alternately.
- **Figs. 10–11 (pp. 154–155), Kron's network model and its x–z planes.**
  Fig. 1 of K62 redrawn: the rings of coils at each node are the ideal
  transformers L = M = ∞; Fig. 11 labels the large mesh ABCDEFGH, the
  small mesh EFKL, inductor currents H₂Δy and voltage drops ¼Ḃ²′ΔzΔx.
- **Figs. 12–13 (pp. 163, 165), two-dimensional field models.** When
  nothing varies along y, Maxwell's equations split into two independent
  sets (eqs. 6-2a, b) and so does the network: Fig. 12 holds B¹′, B³′, E₂,
  H₁, H₃, J²′, D²′ (the TE set), and simplifies to Fig. 13, a plain lattice
  of series inductors (H₁Δx vertical, H₃Δz horizontal) with an R ∥ C shunt
  at every node carrying E₂Δy and J²′ΔzΔx + Ḋ²′ΔzΔx, and no ideal
  transformers at all (p. 164). Fig. 13 is K45a Fig. 14's lattice.
- **Figs. 14–16 (pp. 167, 173, 174), the rectangular waveguide.** WR-430
  (4.3 in × 2.15 in, TE₁₀ range 1.70–2.60 GHz; width taken as 11 cm);
  E₂ = Ê₂ sin(πx/a) sin(βz − ωt), β = √((ω/c)² − (π/a)²) (eq. 6-3); wave
  impedance 377/(βc/ω) (eq. 6-4) and the matching impedance per section
  (eq. 6-5). At f = 2 GHz: βc/ω = 0.731955; with Δx = Δz = 2.2 cm (five
  subdivisions of the width): X_L = j1.579×10⁴·Δy Ω per section,
  X_C = −j1.857×10⁴·Δy Ω per section, theoretical matching impedance
  2.3395×10⁴·Δy Ω per section (pp. 170–171). Side walls: the shunt
  admittances short-circuited (E_y = 0). The computed matching impedance
  is the termination for which the voltage along the axis changes in
  phase only (p. 172). Results: Fig. 15 (fixed f = 2 GHz, λ = 0.15 m, Δz
  varying), per-cent error against (Δz/λ)² nearly linear: ≈ 0 at 3 %,
  ≈ 1.2 at 7 %, ≈ 3 at 11 %, ≈ 6 at 14.5 %, ≈ 10 at 18.4 %, for n = 5 and
  n = 8 subdivisions of the width; Fig. 16 (Δz = 8 mm, f = 3–7 GHz): error
  ≈ 0.5 at 5 %, 4 at 8 %, 7 at 11 %, 12 at 15 %, 20 at 18.5 %. "A
  subdivision along the axis of transmission equal to 10 % of the wave
  length gives an error of about 4 % whereas a subdivision of 5 % of the
  wave length gives an error of about one per cent", agreeing with
  Whinnery and Ramo's network-analyzer studies (p. 175).
- **Fig. 17 (p. 177), transmission-line network for electric machines.**
  Polar coordinates (r, θ, z): μ_θ·ΔrΔz/(rΔθ) series, μ_r·rΔθΔz/Δr shunt,
  R = (Δz/(Δr·rΔθ))(1/σ); (b) the one-dimensional form under sinusoidal
  distribution, μ_r(r)²Δz/(p²Δr·rΔθ) shunt, agreeing with Barton & Cullen.

## 7. N74: Nicholson's figures and equations

- **Eqs. (1)–(3), pp. 2–4.** The orthogonal network: E + e = Z(I + i) or
  V = ZJ; square non-singular connection matrices C = [C_c C_o],
  A = [Aᶜ A°] = (Cᵀ)⁻¹; with trees and links C = [[C_T, B_T],[C_L, 0]],
  A = [[0, A_T],[δ_m, A_L]]; Kirchhoff: C_cᵀE = 0, (A°)ᵀi = 0, Eᵀi = 0; the
  solution

  ```
  [ i^{c′} ]   [ −Z₁⁻¹Z₂        Z₁⁻¹      ] [ I^{o′} ]        [ I^{o′} ]
  [ E_o′   ] = [ Z₄ − Z₃Z₁⁻¹Z₂  Z₃Z₁⁻¹    ] [ e_c′   ]  =  N  [ e_c′   ]     (3)
  ```
  with Z̄ = CᵀZC = [[Z₁,Z₂],[Z₃,Z₄]], Ȳ = AᵀYA, Z₁⁻¹ = (Z_L + A_LZ_TA_Lᵀ)⁻¹,
  Y₄⁻¹ = (Y_T + A_LᵀY_LA_L)⁻¹; and the mesh and node solutions i = L(e − ZĨ),
  E = M(I − Yẽ), L = C_c(C_cᵀZC_c)⁻¹C_cᵀ, M = A°[(A°)ᵀYA°]⁻¹(A°)ᵀ = Z − ZLZ.
- **Fig. 1 (p. 3).** A five-branch example (nodes A, B, C, D; branches
  Z₁ … Z₅ with sources e₁ … e₅; two closed paths 4′, 5′ and three open
  paths A′, B′, C′ carrying I^{A′}, I^{B′}, I^{C′}) with its C and A
  matrices written out:

  ```
         4′ 5′ | A′ B′ C′                4′ 5′ | A′ B′ C′
  C = 1 [−1  . |  1  .  . ]        A = 1 [ .  . |  1  .  . ]
      2 [ 1  1 |  .  1  . ]            2 [ .  . |  .  1  . ]
      3 [ . −1 |  .  .  1 ]            3 [ .  . |  .  .  1 ]
      4 [ 1  . |  .  .  . ]            4 [ 1  . |  1 −1  . ]
      5 [ .  1 |  .  .  . ]            5 [ .  1 |  . −1  1 ]
  ```
  A small, fully specified check case for `kron/orthogonal.go`.
- **Fig. 2 (p. 5), transformation diagram for the sequence of orthogonal
  networks.** Columns 0-network (nodes), 1-network (branches), 2-network
  (planes), 3-network; each a diamond with its C_c^{(i)}, A°_{(i)}, Z_{(i)},
  Y_{(i)}, the source/response pairs (D₍ᵢ₎, I₍ᵢ₎; d₍ᵢ₎, i₍ᵢ₎; E₍ᵢ₎, e₍ᵢ₎,
  B₍ᵢ₎, b₍ᵢ₎; H₍₂₎, h₍₂₎ on the planes), the operators U_i, W_i, S_i, R_i,
  and the incidence matrices M⁰₁, M¹₂, M²₃ joining neighbours; "scattering
  structure" marked across the 2–3 link. The incidence relations (p. 6):

  ```
  C_c^{(i)}(A°_{(i+1)})ᵀ = M^i_{i+1},   (C_c^{(i)})ᵀA°_{(i)} = 0,
  M^{i−1}_i M^i_{i+1} = 0,              L_i(M^{i−1}_i)ᵀ = 0
  ```
  and the wave-propagation recursion (p. 6), a cycle every two dimensions:

  ```
  ḃ₍ᵢ₊₁₎ + Ḃ₍ᵢ₊₁₎ = Z₍ᵢ₊₁₎(M^i_{i+1})† Y₍ᵢ₎ (M^{i−1}_i)ᵀ [ḃ₍ᵢ₋₁₎ + Ḃ₍ᵢ₋₁₎],  ḃ₍₀₎ = 0, i = 1, 3, …
  (M^i_{i+1})† = [C_c^{(i)}(A°_{(i+1)})ᵀ]† ≡ C_o^{(i+1)}(Aᶜ_{(i)})ᵀ,  "identified with curl⁻¹"
  ```
  Transverse wave on the closed-path variables (d, e; h, b), longitudinal
  on the open-path variables (D, E; H, B); propagation from the 0-points
  up through the dimensions (p. 6).
- **Figs. 3–5, 7–9 (pp. 7–8, 14, 17).** The scattering (star-product)
  reading: two obstacles T₁, T₂ with incident and reflected waves
  V₁ … V₆ (eq. 4), the combined matrix T₁*T₂ (eq. 5), and N as a star
  product of tree and link parts (eq. 6); least squares with a priori
  information (eqs. 8–10): x̂ = PHᵀR⁻¹y, P = (S⁻¹ + HᵀR⁻¹H)⁻¹, the
  return-difference matrix F = δ − MR⁻¹; the orthogonal network matrix of
  eq. 3 identified with the least-squares, multimachine and control-system
  problems (Fig. 7); the multistage optimal control sweep (eq. 11–12).
  Nicholson's thesis: Kron's "vertical propagation within each polyhedron
  component introduces effectively a least squares operation, with the
  measurement and covariance matrices defined by the properties of the
  higher-dimensional elements" (p. 14).

## 8. T01: Tonti's figures used here

Fig. 3 (p. 10), outer orientation of points, lines, surfaces and volumes
in 1, 2 and 3 dimensions; Fig. 4 (p. 12), the primal/dual association of
the six global variables (§2 above); Fig. 5 (p. 13), their operational
definitions; Fig. 6 (p. 17), the space complex and the space-time complex
with Φ, Ψ, U, U_m on surfaces; Fig. 9 (p. 21), cell complexes in 1, 2, 3
dimensions and in space-time; Fig. 10 (p. 22), simplicial complexes
adapted to a boundary. Tonti's rule for which variable lives on which
element of which complex is the modern statement of K62's odd/even
assignment and dual polyhedron, and his incidence matrices are Kron's
M^q_{q+1}.

---

## 9. What the catalogue settles

1. **The three-dimensional network is fully specified** by K62 Fig. 1
   (values), B68 Figs. 10–11 (how the ideal transformers work) and B68
   eqs. 5-3/5-4 (SI). Its two-dimensional TE reduction (B68 Fig. 13, K45a
   Fig. 14) has no transformers and is a plain L–C lattice.
2. **The q-network rule for q ≥ 2 is in hand** after all: B68 Fig. 9
   reproduces Kron's 1959 algebraic diagram (orthogonal network per
   dimension, ε on odd and μ on even dimensions, incidence matrices as
   div and curl⁻¹ between them), and N74 pp. 4–6 give the equations
   (eq. 3 per network, the four incidence relations, the propagation
   recursion). Kron 1959 is still wanted to confirm the 1/4 and sign
   conventions at q = 2, 3, but the build can start.
3. **The oscillatory frame is not in any source** (K62 p. 67; B68 p. 20
   confirms), but the tables determine what it does: see §11.
4. **Reproduction targets with numbers**: K45a Fig. 14 (cavity potentials
   and Fig. 15 residuals); K45a Figs. 9–10 (ODE and heat lattices); B68
   Figs. 15–16 (waveguide matching-impedance error against Δz/λ, with the
   1 % at 5 % and 4 % at 10 % rule); K62 Tables I–VI (regression and
   polyhedral fits, Table I transcribed above); N74 Fig. 1 (C and A
   matrices).

## 10. Build plan and status

| Stage | Build | Check | Status |
|---|---|---|---|
| A | `kron/polyhedron.go`: `FieldCell` (K62 Fig. 1 in SI), the TE lattice (B68 Fig. 13), `Kron1945Cavity`, `Waveguide` | K45a Figs. 14–15; TE box against the dispersion FDTD uses; B68 Figs. 15–16 | **Done.** Cavity worst 5.9 % of the largest potential; Fig. 15 recomputed at 23 of 25 junctions within 0.33 A; TE₁₁ to 2e-14 of the lattice value and 0.36 % of the continuum; waveguide +0.69, +1.50, +3.47, +6.43 % against about 1.0, 1.6, 3.4, 6.2 % |
| B | `kron/maxwell3d.go`: full 3-D cell with ideal transformers (K62 Fig. 1, B68 Figs. 10–11) | straight-wire B, coaxial TEM pulse, cavity modes against FDTD and TLM | **Partly done** (2026-11-05). Spheres of quarter inductors, R∥C branches, Kron's two diagonal transformers per large mesh as k = 1 coils of L_w = Ideal·L/4. MNA solution = space-discrete Yee to 8e-8 (E·Δ) and 1.5e-7 (every arc's H·Δ) at Ideal = 10⁷, falling as 1/Ideal; three ties per mesh give the same, so Kron's two suffice. TM110 pole = lattice to 2.5e-8; em/fdtd.go rings at the same ω through its leapfrog to 3.8e-7; AC ledger closed. Straight wire: Ampère round every loop to 5e-8, H against the line current in its tube with error ≈ 0.32·(Δ/r)² (discrete Laplacian), under 1 % from 5.5 cells. Coaxial TEM pulse (square line, kronspice transient, Ideal 1e3): V/I = lattice Z₀ 55.13 Ω to 0.2–0.7 %, inverted V at both shorts, group delay = lattice dispersion 1/(8σ²), ledger closed, em/fdtd.go (Conductor) V(t) within 1.5 %. Transients need Ideal ≲ 1e3: the attachment points are tied together only through the windings' magnetizing path h/(2L_w). Open: TLM twin |
| C | `kron/curvefit.go`: the K62 Fig. 3 polyhedron, frames as residual directions, the one-mesh circuit | Tables I–VI | **Partly done.** Given columns 1.9e-6; frame 1 8.9e-6; mesh circuit dissipates the error-square. Second planar frame and LS + linear not reproduced (§11). The general orthogonal network per q (N74 eq. 3) is not built |
| D | GA (`evo/`) over the frame | oscillatory columns | **Superseded in part.** The GA recovers frame 1 from the printed column exactly, but does not reach zero power on Kron's objective; the closed form does (§11) |
| E | Time networks (K45a Figs. 3–7) | K45a Figs. 9–10 numbers | Not started |

Stages A–D are numerics on standard physics and carry no hypothesis. Nothing
here tests the physical claims of H20 (plasma, crystals), which stay
*untested*.

## 11. What the tables say about Kron's frames

Results of `examples/kron_polyhedron`; inferences from internal consistency
of the tables, not statements in Kron's text. Kron 1959 and the General
Electric report 60GL164, not held, may bear on them.

- **Every column should be linear in z.** All 13 characteristics are linear in
  the four z values, and six tables overdetermine a 4-variable map, so each
  column can be tested for being one fixed map without knowing the frame.
- **First planar frame: an orthogonal projection.** Error rank one to the
  tables' rounding, with residual direction w = (−3, 1, 1, 1 | 4, −1, −3, −4, −4,
  2 over the branch lengths | 2, −10, −4 over 2A = 9). Integers reproduce the
  column to 8.9e-6.
- **Oscillatory column: the identity.** Equal to the given column to 2.8e-5
  for all six functions. A frame with w orthogonal to the 4-dimensional range of
  the characteristics does this exactly for any values at the four points,
  which accounts for Kron's random-number remark (p. 69). N = (9, M₂/2, M₁/2,
  −13, −10, M₃/2) with any even M. A residual below about 1e-5 cannot be
  excluded, and the printed column cannot identify the member of the family.
- **Second planar frame: unresolved.** Function rows are frame 1's residual
  times 0.3032; divided differences are not linear in z.
- **LS + linear: linear with two corrections, otherwise unresolved.** Rank-3
  error map on the inner branches; two outer branches reversed in sign.
- **Second differences.** S = −(Δz/L)(A/L) on the plane's branch opposite the
  interior vertex, read from the tables to 1.9e-6; the paper does not write it.
- **Apparent misprints:** Table II LS + linear row 2 (−12.339281 for
  −12.330283); Table I LS + linear fourth divided difference (−6.4341227 for
  −0.64341227); K45a Fig. 15 row 2 col 3 (sign), row 1 col 0 (2.12, computed
  −5.94, unresolved).
