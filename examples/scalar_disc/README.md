# The optically gated rotating disc: a bench test of the scalar S (H47 step 4)

Iadicicco & Verolino (arXiv:2608.21403) propose a falsifiable test of the
Aharonov–Bohm scalar S = ∂_μA^μ, the scalar part of Maxwell's quaternion
equations that standard theory sets to zero (CLAUDE.md §3c, "With the
scalar part kept"). A fibre ring on a rotating aluminium disc lights one
angular sector, which induces a charge in the metal under it; the sector
turns with the disc. In a rotating frame the split of the conserved current
gives a source for S (Iadicicco, Modanese & Verolino, *Gen. Rel. Grav.* 58,
106, 2026, eqs. 13, 19):

```
I_G = Ω ∂_φ ρ,      □S = μ₀ I_G
```

This example computes S from `em/scalar`'s field equations (quasi-static:
the rotation is 10⁸ times slower than light crossing the disc) and the
gravitational effect of S² from Minotti & Modanese's weak-field coupling
(*Eur. Phys. J. C* 86, 593, 2026, eqs. 11–13: −∇²χ = (K_φC_φ/2c²)S²,
f = −∇χ), with their fitted constants K_φ = 10⁻¹⁶ A² N⁻², C_φ = 6×10²¹.
The standard prediction is S = 0 and no force beyond Newton and Coulomb.

```
go run ./examples/scalar_disc
go test ./examples/scalar_disc
```

## Declared bench

The paper gives no dimensions, speeds or charge. Declared: 5 cm Al disc,
lit ring at r = 4 cm, 5 mm wide, sector 15° FWHM holding 1 nC, 3000 rpm,
sensor 20 mm above the track.

## Results (8 November 2026)

| Quantity | Value | Standard prediction |
|---|---|---|
| Peak source \|I_G\| | 0.0324 C m⁻² s⁻¹ (= IMV eq. B3 to 3e-4); ring integral 0 (eq. 24) | — |
| Peak S at the sensor | 1.10×10⁻¹² T (cS = 0.33 mV/m) | 0 |
| Reversed rotation | S → −S exactly | 0 |
| Whole ring lit | S = 0 exactly | 0 |
| S² force per mass (MM fit) | 5.4×10⁻³⁷ m/s² peak, the same for both senses of rotation | 0 |
| At the bench's limits (air breakdown: 1.5 nC; Al burst: 1.05×10⁴ rad/s, 100,000 rpm) | S 5.5×10⁻¹¹ T; force 1.3×10⁻³³ m/s² | 0 |
| Ordinary systematic: Coulomb field of the gated charge at the sensor | 2.1×10⁴ V/m; on a 1 g sensor carrying 1 pC, 2.1×10⁻⁵ m/s² | the same |

The full log is `run_output.txt`.

## What it shows

1. **The two papers predict different symmetries.** Iadicicco & Verolino
   assume a response linear in the source, so it reverses with the rotation
   and that reversal is their decisive control. Minotti & Modanese's
   gravity coupling is through S², which does not reverse. A CW/CCW run
   tells the two apart.
2. **Under the 2026 fit the disc's force is far below any measurement.**
   The superconducting junction S that the fit relies on is ≈ 2×10³ T, about
   4×10¹³ times the disc's best. To reach 10⁻⁹ m/s² at the bench's limits,
   C_φ would have to be 5×10³³, 9×10¹¹ times the fitted value; a null at
   that sensitivity sets that bound.
3. **The fit implies a small Brans–Dicke parameter.** 2ω₀ + 3 =
   8πG₀ε₀C_φ/(c²K_φ) ≈ 9.9, so ω₀ ≈ 3.5, against Cassini's ω > 4×10⁴; the
   authors need the matter coupling suppressed (their Appendix B).
4. **The linear response is what the disc can measure.** Neither paper
   gives a value for Λ_opt. The disc's S is 3.7×10⁻¹³ T per nC per
   1000 rpm, so a measured response reversing with Ω and vanishing for the
   whole ring lit, divided by this S, is Λ_opt; a null bounds it.
5. **The controls separate it from the ordinary systematic.** The Coulomb
   force is 10³¹ times larger than the S² force, but it follows the charge,
   does not reverse with rotation and stays when the whole ring is lit;
   the S signal reverses and vanishes.

## Checks (`go test`)

- `TestSource`: total charge Q; peak source = IMV eq. B3; ring integral 0.
- `TestSmoothing`: widening the sector at fixed charge lowers the peak source
  by σ²/(σ² + λ²) (IMV eq. B7).
- `TestSignatures`: S odd in Ω, linear in QΩ, zero for the whole ring; the
  S² force even in Ω.
- `TestFarField`: S falls as 1/R² far away (the source has no net charge).
