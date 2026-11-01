# Verification of CLAUDE.md §3a and §3a.1 against Kron

Sources checked (both in `reference/`):

- **Kron 1939** — *Tensor Analysis of Networks*, Wiley 1939 (2nd printing 1940),
  scan `dli.ernet.8037`. Printed page = PDF page − 27.
- **TfC** — Kron, *Tensors for Circuits*, Dover 1959 (reprint of *A Short Course
  in Tensor Analysis for Electrical Engineers*, 1942). Image-only scan, read by
  OCR. Printed page = PDF page − 25.

Page numbers below are printed page numbers. Equation numbers are given where
the scan shows them legibly.

## Summary

| # | §3a / §3a.1 statement | Verdict |
|---|---|---|
| 1 | `i = C·i'`, `e' = Cᵀ·e`, `Z' = Cᵀ·Z·C` | **Confirmed** (Kron 1939 pp. 101–104, eqs 4.23–4.37; TfC p. 43, eqs 6.6–6.7) |
| 2 | Primitive equation `e = Z·i` | **Incomplete** — true for mesh networks only; Kron's general (orthogonal) form is `E + e = z·(i + I)` |
| 3 | Node-pair dual `Y' = Aᵀ·Y·A` | **Confirmed**, with `A = Cₜ⁻¹` in Kron's notation (Kron 1939 p. 360; TfC p. 43, eq 6.8) |
| 4 | `m = b − N + 1` closed paths, `N − 1` open paths | **Correct only for one connected piece.** Kron's formula is junction-pairs = junctions − sub-networks (Kron 1939 p. 75, eq 3.1; TfC p. 45, eq 6.12). Torn networks need the sub-network term |
| 5 | Power invariant `P = eᵢ iⁱ` | **Incomplete** — in orthogonal networks Kron keeps two invariants: input `e·i` and output `E·I` |
| 6 | Inactive junction pairs are in Ch. V | **Wrong chapter.** Ch. V has potentials across open-circuited coils but not the term "inactive" |
| 7 | Opening of circuits as transformations, Ch. VI | **Confirmed** (Ch. VI §XIV, p. 162) |
| 8 | Impressed voltages, Ch. X | **Confirmed but narrow** — Ch. X covers impressed voltages *in eliminated axes* (§V p. 248, §IX p. 254); the concept is introduced in Ch. III (p. 76) |
| 9 | "Verified" page cites pp. 75, 76, 407, 408, 411, 429 | **All confirmed** |
| 10 | Each element port is an apparent coil | **Conflicts with Kron's method** — see below |
| 11 | Kron field networks (1944): region boundary as open paths | **Not verifiable from these books.** TfC Ch. 14 is not an equivalent-circuit treatment |

## Details and suggested edits

### 1. Transformation laws — confirmed

- Kron 1939 p. 101–102: `i = C·i'` (4.23); power invariance `e·i = e'·i'` (4.24),
  stated as a consequence of the first relation, not a separate assumption;
  `e' = Cₜ·e` (4.26–4.27).
- Kron 1939 p. 103–104: `e = z·i` (4.33) keeps its form `e' = z'·i'` (4.34), giving
  `z' = Cₜ·z·C` (p. 104, following 4.37).
- TfC p. 43: same result as eqs 6.6–6.7, and the note that in most problems C
  is rectangular (singular); laws that do not need C⁻¹ stay valid.

Caveat: Kron 1939 Ch. IV derives these with `Cₜ⁻¹`, i.e. for square C. The
singular-C case (the normal mesh case) is handled in Ch. V §§III–IV
(pp. 119–121). Cite Ch. V for the code path that uses a rectangular C.

**Suggested citation for §3a:** Kron 1939 Ch. IV pp. 101–104 (eqs 4.23–4.37) and
Ch. V pp. 119–121; TfC Ch. 6 p. 43.

### 2. The primitive equation

`e = Z·i` is the mesh-network form. Kron's general form, which §3a.1 needs
because it uses junction-pair (open-path) variables, is the orthogonal
equation `E + e = z·(i + I)`:

- Kron 1939 p. 417 (final set of equations of an orthogonal network) and p. 480
  (eq 19.3, split into mesh and junction-pair blocks).
- TfC p. 49, footnote: the primitive branches obey `e + E = z·(i + I)`, not
  `e = z·i`, once the dual axes are included.

Here `e`, `I` are impressed (known) and `i`, `E` are responses. §3a should
show both forms and say the mesh form is the special case `E = 0, I = 0` on
the axes used.

### 3. Junction (node-pair) form

Kron 1939 p. 360 (Ch. XIV): `Y' = C⁻¹·Y·Cₜ⁻¹`. TfC p. 43, eq 6.8: same. Writing
`A = Cₜ⁻¹` gives the spec's `Y' = Aᵀ·Y·A`. Suggest stating that identity in
§3a so the `kron/` code names match Kron's.

### 4. Counting closed and open paths — affects tearing tests

- Kron 1939 p. 75, eq 3.1: junction-pairs = junctions − sub-networks; eq 3.2:
  coils = meshes + junction-pairs.
- TfC p. 45, eqs 6.12–6.13: `M = C − (J − S)`, `C = M + P`.

So with S sub-networks: open paths `N − S`, closed paths `b − N + S`. The spec's
`m = b − N + 1` is the S = 1 case. **A torn network is several sub-networks**,
so the §3a.1 test "closed-path and open-path counts satisfy m = b − N + 1" will
fail on torn pieces as written.

**Suggested test wording:** `m = b − N + S` and `p = N − S`, with S the number
of connected sub-networks (1 before tearing).

### 5. Power invariance — the ledger uses the output invariant

Kron 1939 p. 409 (Ch. XVI §II(d)): the orthogonal method rests on keeping the
number of variables equal to the number of coils, and on two invariants —
power input `e·i` and power output `E·I`. Kron 1939 pp. 412–413 (§VI) explains why
both stay invariant under any interconnection once all four quantities are
kept.

For §3a.1 this matters: the power through the open paths (junction pairs) is
the `E·I` term, not `e·i`. Suggested wording for "Open paths carry the
exchange with the outside": "Power is invariant under every Kron
transformation; in the orthogonal frame both the input `e·i` and the output
`E·I` are invariant (Kron 1939 p. 409, pp. 412–413). The exchange through open
paths is the `E·I` term."

### 6. Where "inactive" junction pairs are actually defined

Ch. V does not use the term. The sources are:

- **Kron 1939 Ch. XIX p. 485 (§IV, "The elimination of variables")** — the
  definition: junction-pairs through which no currents I flow (and meshes
  with no impressed e) are inactive for analysis. In setting up z, inactive
  junction-pairs are ignored; inactive meshes need the Ch. X reduction
  formulas.
- **Kron 1939 Ch. XVI pp. 428–429** — the same rule for orthogonal networks (already
  cited in the spec).
- **Kron 1939 Ch. XXII pp. 566–567** — junction-pairs classed as **input, output and
  inactive**, with the inactive ones eliminated so the analysis runs on active
  axes only. This is the closest match in Kron to the KronSpice port model and
  is worth citing directly. Ch. XXII is not in your chapter list; consider
  adding it.

Ch. V is still a good cite for **potentials across open-circuited coils**
(pp. 128–129) and for **interconnection of several networks** (§X p. 133,
below).

**Suggested replacement** in §3a.1, first bullet: "…inactive junction pairs
(Ch. XIX p. 485; Ch. XVI pp. 428–429; Ch. XXII pp. 566–567), potentials across
open-circuited coils (Ch. V pp. 128–129), opening of circuits as
transformations (Ch. VI p. 162), impressed voltages (Ch. III p. 76; in
eliminated axes, Ch. X pp. 248, 254)…"

### 7–9. Verified cites — all confirmed

- p. 75 eqs 3.1–3.2 ✔︎ (see item 4).
- p. 76 open mesh, impressed vs response quantities ✔︎.
- p. 407 apparent coil defined for an impedanceless branch with known
  impressed voltage ✔︎. The admittanceless case (known impressed current)
  continues onto **p. 408**; cite pp. 407–408.
- p. 408 zero-valued currents impressed across junction-pairs "for the sake
  of completeness" ✔︎.
- p. 411 closed ↔ open mesh by an apparent coil; axis is a junction-pair when
  E is unknown, a mesh when i is unknown ✔︎.
- p. 429 inactive junction-pairs ignored in z; rectangular C as part of a
  square C ✔︎. Note the same page also says **none of the meshes may be
  ignored** when setting up z — worth adding to the rule, since a port
  modelled as a closed mesh cannot be skipped the same way.

### 10. Ports as apparent coils

The spec says each element's external port **is** an apparent coil. Kron
introduces apparent coils only to show the equivalence, then removes them:

- Kron 1939 p. 408 (§I(c)): one purpose of the orthogonal analysis is to remove
  zero-impedance/zero-admittance coils, so apparent coils and their
  superfluous axes disappear.
- Kron 1939 p. 412: apparent coils do not appear in the primitive network.
- Kron 1939 p. 417: they no longer appear as coils; their effect is known
  components of E' and i'.

**Suggested wording:** "Each element port is a junction-pair axis (open path)
from the element to `Ω_channel`. A driven port has a known impressed current
I (or a known E); in Kron's terms it could be drawn as an apparent coil, but
it enters the orthogonal equation only as known components of E' and i', not
as an extra coil in the primitive z (Kron 1939 pp. 407–408, 412, 417)."

This also matches the implementation advice: no zero-impedance rows in z.

### 11. Field regions as networks

TfC Ch. 14 ("The Field Equations of Maxwell", pp. 109–112) gives the 3-D
vector equations and Minkowski's 4-D tensor form (`F_αβ`, `H^αβ`, Heaviside–
Lorentz units). It is **not** an equivalent-circuit treatment and does not
treat region boundaries as open paths. It may be a useful tensor cross-check
for §3c, but §3a.1 item 8 still needs the 1944 Proc. IRE paper.

## Additional supporting passages worth citing

- **Multi-domain interconnection (§3a.1 items 3 and 7):** Kron 1939 Ch. V §X p. 133
  — separate networks interconnected by a C; the apparatus may be electrical,
  mechanical, acoustical, optical and so on. It also describes dividing a
  system by opening branches, analysing parts separately and recombining by
  a C — the 1939 seed of tearing.
- **Physical subdivision (diakoptics precursor):** Kron 1939 Ch. XIX §IX p. 492.
- **Computing C from Kirchhoff constraints:** TfC Ch. 7 pp. 46–50 — writing
  the junction constraints as `B·i = 0` and turning them into `i = C·i'`
  (dependent and independent currents). Directly usable for building C from
  the netlist incidence matrix (C spans the null space of B).
- **Open/closed exchange of variables under constraint:** Kron 1939 p. 409–410 —
  opening a branch replaces a contravariant variable (mesh current) by a
  covariant one (junction-pair voltage), and vice versa; the total stays n.

## Not changed

No edits were made to CLAUDE.md. The suggested wording above can be applied
section by section.

## Addendum (24 Oct 2026): TAN papers, tearing, Maxwell

**TAN claims in §3a.1 — confirmed.** arXiv:1412.0187 states currents as
1-chains and potentials as 0-cochains (§5), and that the natural space for
the electromotive force is the cochain of degree 2. WCE 2013 defines chords
as off-diagonal impedance terms replaced by functions of any complexity,
and the affine interaction as a delayed source emf `E·e^(−x p / v)`. The
validations are as stated: cavity with aperture against FDTD (arXiv §13),
10 GHz horn and metal wall or absorbers against measurement, 1.2% (§14),
and cavity monopoles against measured S21, 100 MHz–2 GHz (Leman et al.
2009). Added: the radiation resistance is intrinsic to the radiating edge
but modified by the environment (walls), so the radiation edge's value is
environment-dependent.

**Tearing — algorithm added.** Aitchison (1983) gives an explicit,
implementable form, `z = Q(I − B⁻¹K(LB⁻¹K)⁻¹L)B⁻¹Pc`, with `LB⁻¹K`
non-singular iff `A` is (Lemma 1), and requires each torn block of `B` to
be non-singular (Lemma 2). Roth (1959) gives ohmicness as the necessary and
sufficient existence condition and proves power invariance under tearing.
§3a now carries the algorithm, the singular-block caveat, the cost
estimate, GA caching rules and tests for torn = untorn and power
invariance.

**§3c against Maxwell, Treatise Vol. II, Arts. 618–619 (pp. 236–238).**
Equations (A), (B), (E), Ohm's law, displacement and total current match.
Corrections made: (A) carries Maxwell's condition `S.∇𝔄 = 0` (Coulomb
gauge); the general magnetization relation is (D) `𝔅 = ℌ + 4π𝔍`, with
`𝔅 = μℌ` only for induced magnetization; `ℌ = −∇Ω` added; Maxwell's
"convergence" `S.∇σ = −div σ` cited (Vol. I, Art. 25). Open item: the signs
of the electric and magnetic density relations in Art. 619 must be read from
the page images, as the scan's text layer drops signs.
