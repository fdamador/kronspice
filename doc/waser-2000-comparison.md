# Waser (2000) against kronspice's Maxwell equations

Source: André Waser, *On the Notation of Maxwell's Field Equations*,
AW-Verlag, issued 28.06.2000, 14 pp. (copy:
`reference/waser_2000_notation_maxwell.pdf`). Equation numbers below are
Waser's. Engine side: `em/maxwell.go`, `em/biquat.go`, `em/vcalc.go`,
conventions in `em/doc.go`, CLAUDE.md §3c. Tests: `em/waser_test.go`.

Summary: every physical relation in the paper that belongs to standard
electromagnetism is already in the engine. The paper has sign and labeling
errors in its modern-vector section that lead it to a wrong result
(∂ρ/∂t = 0). Its extensions are either rewrites of the same equations or
hypotheses already in, or suited to, the register (CLAUDE.md §1).

## 1. Maxwell 1865, eqs. (1.1)–(1.8) (Waser p. 2–3)

| Waser | Content | Engine |
|---|---|---|
| (1.1) | total current = conduction + displacement | `MaxwellQuaternion` / `MaxwellVector` Ampère residual |
| (1.2) | μH = ∇×A | `FieldsFromPotentials` |
| (1.3) | ∇×H = J (total, Gaussian 4π) | Ampère–Maxwell residual (SI) |
| (1.4) | E = μ(v×H) − ∂A/∂t − ∇φ | `FieldsFromPotentials` for conductors at rest (no v×B term); moving conductors via `em/lorentz.go` and the machine element |
| (1.5) | εE = D | SI constitutive relation |
| (1.6) | σE = j | Ohm's law, resistor stamp |
| (1.7) | ρ + ∇·D = 0, i.e. −ρ = ∇·D | **sign differs**: engine uses ∇·D = +ρ (§4 below) |
| (1.8) | ∂ρ/∂t + ∇·j = 0 | KCL and the ledger's charge check; `TestWaserContinuity` |

## 2. Maxwell 1873 quaternion form, eqs. (1.9)–(1.20) (Waser p. 4)

| Waser | Treatise (Art. 619 letter) | Engine | Verdict |
|---|---|---|---|
| (1.9) B = V.∇A | (A) | `FieldsFromPotentials` (`B = V.∇𝔄`) | same |
| (1.10) E = V.vB − Ȧ − ∇Ψ | (B), with 𝔊 the velocity | same, at rest | same |
| (1.11) F = V.vB + eE − m∇Ω | (C) mechanical force | `em/lorentz.go` (per charge) | **check**: as recalled the Treatise has V.ℭ𝔅 (current density), not velocity; Waser's "v" looks like a transcription slip — to be read from the page image |
| (1.12) B = H + 4πM | (D) | SI B = μ0(H + M), `device/core.go` | same |
| (1.13) 4πJ_tot = V.∇H | (E) | Ampère–Maxwell (vector part of `R2`) | same |
| (1.14) J = CE | (F) | Ohm's law | same |
| (1.15) D = KE/4π | (G) | D = εE | same |
| (1.16) J_tot = J + Ḋ | (H) | total current | same |
| (1.17) B = μH | (L), induced magnetization only | linear media | same |
| (1.18) e = S.∇D | (J) | ρ = −S.∇D = +∇·D | **sign question**, §4 |
| (1.19) m = S.∇M | (K) | pole density −∇·M | same; not a monopole (§5) |
| (1.20) H = −∇Ω | (I) | magnetostatic potential | same |

Treatise letters (F)–(L) follow the 1873 order as recalled and are to be
confirmed from the page image with (C), (J) and (K).

## 3. Errors in Waser's modern-vector section (pp. 5–6)

1. **Labels swapped, (1.26)–(1.27).** ∇×H = ∂D/∂t + j is the
   Ampère–Maxwell law; −∇×E = ∂B/∂t is Faraday's law.
2. **Sign, (1.24).** ∇×E = −μ ∂H/∂t, not +. A vacuum plane wave fails
   Waser's form by 2·∂B/∂t (`TestWaserFaradaySign`).
3. **Sign, (1.28).** Modern Gauss law is ∇·D = +ρ: a positive charge has
   outward D (Coulomb's law). Waser carries over the sign of 1865 (1.7)
   without the convention that goes with it.
4. **Wrong conclusion, (1.32).** From (1.28) and (1.8) Waser obtains
   ∂ρ/∂t = 0 — no charge could ever accumulate, so no capacitor could
   charge. With ∇·D = +ρ the divergence of Ampère–Maxwell is the identity
   ∂ρ/∂t + ∇·J = 0 for any fields. `TestWaserContinuity` builds a field
   whose charge grows (∂ρ/∂t up to 6 C/m³/s) with the current Ampère–Maxwell
   requires: the engine's continuity residual is 1.4e−14 (round-off), while
   Waser's (1.28) + (1.8) leave a residual of exactly 2·∂ρ/∂t.
5. **Typos.** (1.52)–(1.53) both start "E ="; the second is the magnetic
   field. (1.74) writes ∂/∂x three times; it is ∂/∂x, ∂/∂y, ∂/∂z.

## 4. The sign of e = S.∇D

Hamilton's operator gives S.∇D = −∇·D (Maxwell's "convergence", Treatise
Vol. I Art. 25; `em/doc.go`). So (1.18) read literally is e = −∇·D, the
same sign as the 1865 equation (1.7) `e + df/dx + dg/dy + dh/dz = 0`. The
engine uses ρ = +∇·D, the sign every measurement of Coulomb's law gives.
The difference is a sign convention in Maxwell's text (or a slip in it),
not different physics; `TestWaserGaussSign` pins the engine's sign: a
positively charged ball has S.∇D = −ρ0, the engine's Gauss residual
vanishes with ρ = +ρ0, and Waser's sign leaves −2ρ0.

**Open item (CLAUDE.md §3c):** read the Treatise Art. 619 page image
(scan in `reference/` on the project machine) for (C), (J) and (K), and
record Maxwell's printed sign and convention in `em/doc.go`.

## 5. "Magnetic mass" is not a monopole

Waser calls (1.19) Maxwell's first formulation of a magnetic charge
density. It is the pole density of magnetization — the field-theory form
of a bar magnet's poles. With B = H + 4πM (SI B = μ0(H + M)) and B = ∇×A,
∇·B = 0 everywhere while ∇·H = −∇·M is non-zero where magnetization
varies (`TestWaserMagneticMass`). Free magnetic charge, with ∇·B ≠ 0, is
hypothesis H2 (`em/duality.go`), separate from this.

## 6. Extensions (pp. 7–12)

| Waser section | What it is | Engine | Status |
|---|---|---|---|
| Hertz–Phipps, (1.40)–(1.47): ∂/∂t → d/dt = ∂/∂t + v·∇ | Galilean-invariant equations for moving media | `em/frames.go` uses Lorentz transformations (inertial frames, flat spacetime); rotating and other accelerated frames use Kron's Riemannian / non-holonomic machine frames (`kron/machine.go`, Kron 1934); gravitational curvature is out of scope (CLAUDE.md §1 non-goals) | Contradicted for moving media: Fizeau (1851) partial light drag, Eichenwald (1903) and H. A. Wilson (1904) moving-dielectric experiments. Candidate for the register with these as its recorded bound |
| Dirac, (1.48)–(1.53) | dual-symmetric equations, magnetic charge and current | `em/duality.go` | H2, hypothetical, off by default |
| Harmuth, (1.54)–(1.58): magnetic conductivity s | magnetic current J_m = sH, no sources | can be a J_m term in `em/duality.go`, own ledger column | not yet modeled; Barrett's critique (IEEE Trans. EMC 30, 1988) to be recorded with it |
| Meyl, same equations for all cases | no electric or magnetic charges | — | H11; ∇·E = 0 everywhere is contradicted by measured charge |
| Múnera–Guzmán, (1.59)–(1.68): N = B − E, P = B + E | linear recombination (units c = 1) | equivalent to `em/biquat.go` F = E + IcB | same physics; their extra scalar field needs assumptions beyond (1.59)–(1.62) |
| Minkowski tensor, p. 10–11 | covariant form | `em/biquat.go` (∇F = J) | same physics |
| Inomata, Rauscher, (1.69)–(1.72) | magnetic sources labeled with complex i | `em/duality.go` | same as Dirac; the "i" labels, adds no physics |
| Honig, (1.73)–(1.75) | biquaternion form | `em/biquat.go` | same physics |

On imaginary numbers: Maxwell used only the scalar and vector parts (S.
and V.) of Hamilton's quaternions, whose units i, j, k are not the complex
i (Waser p. 3). In the engine the complex i appears only in AC phasors
(`complex128`).

## 7. Next steps

- Page-image check of Treatise Art. 619 (C), (J), (K) and the letter
  order; record in `em/doc.go`.
- Optional hypothetical modes: Harmuth's magnetic conductivity in
  `em/duality.go`; Hertz's total-derivative equations with the
  Fizeau/Eichenwald/Wilson bounds as validation data.
