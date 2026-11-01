# M9i — lumped parameters and transmission lines from geometry

Files: `em/extract.go`, `em/inductance.go`, `em/tline.go`, tests
`em/extract_test.go`, `em/tline_test.go`; example `examples/pcb_layout`.

## Field extraction (`em/extract.go`)

`Section` is a 2-D cross-section: conductors (`Disk`, `Rect`, `Strip`,
`Outside`; names starting with `gnd` are grounded, as is the box edge) and a
permittivity map `EpsR(x, y)`. `Solve` runs the quasi-static equation
`div(ε grad φ) = 0` — the scalar part of Maxwell's `∇𝔇 = 4πρ` with
`𝔇 = K𝔈/4π`, `𝔈 = −∇ψ` (Treatise Art. 619; §3c table row "Capacitor") —
by finite volumes with Shortley–Weller cut edges at curved boundaries,
conjugate gradients with an SSOR preconditioner (ω = 1.6), and the charge
from the flux through each conductor. It returns the Maxwell capacitance
matrix `C` and its vacuum twin `C0`; the TEM inductance matrix is
`L = μ₀ε₀ C0⁻¹`.

| Case | Reference | Result |
|---|---|---|
| Coax, εr = 2.1 | 2πε/ln(b/a) | +0.07 % (C, L, Z₀) |
| Stripline, zero thickness | Cohn 1954, exact K(k)/K(k′) | −0.006 % after Richardson extrapolation (b/40, b/80); field ∝ r^−½ at the edge makes the raw error first order |
| Two-wire line | πε₀/acosh(D/2a) | 0.1 % in a 120 mm box |
| Microstrip, FR-4, w = 3 mm | Hammerstad–Jensen 1980 | Z₀ −1.8 %, ε_eff −0.4 % |

## Inductance (`em/inductance.go`)

`Neumann` integrates `M = (μ₀/4π)∮∮ dl₁·dl₂/r` over polygon loops with
Gauss–Legendre points; `MutualCircles` is Maxwell's elliptic-integral form
for coaxial circles (Treatise Art. 701), with K and E by the AGM;
`SolenoidSheet` integrates it over two current sheets and equals
`Nagaoka`'s closed form to 1e−6. Neumann vs Maxwell: 0.007 %. Wheeler's
formula (`em.Coil`) lies within 0.6 % of Nagaoka.

## Lines in the circuit (`em/tline.go`)

`TLine{R, L, G, C, Length}` per unit length (`FromSection` builds one from
an extraction); `Z0`, `Delay`; `Ladder` adds n T-sections to a builder;
`CoupledLadder` adds an N-conductor line with mutual inductance as K
elements and mutual capacitance as capacitors between conductors. The
telegrapher test: a 2 m microstrip-like line gives Z₀ 51.85 Ω and a
transient delay of 9.67 ns (0.690 c), matching √(LC) and √(L/C).

## Example: `examples/pcb_layout`

Genes: width w and gap s of two microstrips on 1.6 mm FR-4, 10 cm long.
Each candidate is extracted, built as a 20-section coupled ladder with
both traces terminated in Z₀, and solved by `.AC` at 100 MHz. Objectives:
|Z₀ − 50 Ω|, near-end crosstalk above −30 dB, and a small weight on board
width.

    go run ./examples/pcb_layout -gens 8 -pop 16
    evolved (114 evaluations): w 2.769 mm, s 8.913 mm, board 14.45 mm
      Z₀ 49.86 Ω, ε_eff 3.281, K = L₁₂/L₁₁ 0.0130, NEXT −50.8 dB

The test also checks that K and NEXT fall monotonically as the gap grows
(s = 0.5, 2, 6 mm: K 0.287, 0.128, 0.030; NEXT −20.8, −29.1, −43.3 dB).

Limits: zero-thickness strips; quasi-static (TEM) only, so dispersion and
radiation are not included — `em/limits.go` applies; the ladder's own
λ/10 rule sets the section count.
