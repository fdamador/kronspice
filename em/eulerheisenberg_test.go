// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"math"
	"testing"
)

// The Euler–Heisenberg vacuum: constants, constitutive law as the
// Lagrangian's derivatives, small-signal tensors as the constitutive law's
// derivatives, Lorentz invariance, energy density, and vacuum
// birefringence n∥ = 1 + 7a, n⊥ = 1 + 4a with Δn/B² = 4.0×10⁻²⁴ T⁻².
func TestEulerHeisenberg(t *testing.T) {
	// Critical fields and the coefficient.
	if r := math.Abs(SchwingerField/1.3232855e18 - 1); r > 1e-6 {
		t.Errorf("E_s = %.8g V/m", SchwingerField)
	}
	if r := math.Abs(1/FineStructure/137.035999 - 1); r > 1e-8 {
		t.Errorf("1/α = %.10g", 1/FineStructure)
	}
	// ξc²B²/ε₀ = (α/90π)(B/B_c)²: the two forms of a agree.
	if r := math.Abs(EHXi*SpeedOfLight*SpeedOfLight/Eps0/EHBirefringenceParameter(1) - 1); r > 1e-12 {
		t.Errorf("ξc²/ε₀ against α/(90π B_c²): %.3g", r)
	}
	// PVLAS's target: Δn/B² = 3a/B² ≈ 4.0×10⁻²⁴ T⁻² (Ejlli et al. 2020).
	// (n∥ − n⊥ = 3a; formed as a difference of indices it is lost in 1.)
	if dn := 3 * EHBirefringenceParameter(1); math.Abs(dn/3.99e-24-1) > 1e-2 {
		t.Errorf("Δn at 1 T = %.4g", dn)
	}
	t.Logf("E_s %.6g V/m, B_c %.6g T, ξ %.4g, a(1 T) %.4g, Δn(1 T) = 3a = %.4g", SchwingerField, CriticalB, EHXi, EHBirefringenceParameter(1), 3*EHBirefringenceParameter(1))

	// A field at 3 % of critical, all components and E·B ≠ 0.
	e := [3]float64{0.02, -0.011, 0.017}
	b := [3]float64{-0.009, 0.015, 0.012}
	for i := range e {
		e[i] *= SchwingerField
		b[i] *= CriticalB
	}
	d, h := EHConstitutive(e, b)
	// D = ∂L/∂E, H = −∂L/∂B by central differences of the nonlinear part
	// (the linear part is exact and would drown the 1e-6 correction).
	nl := func(e, b [3]float64) float64 {
		c2 := SpeedOfLight * SpeedOfLight
		f, g := dot(e, e)-c2*dot(b, b), dot(e, b)
		return EHXi * (f*f + 7*c2*g*g)
	}
	// Relative to the largest component of each quartic part (single
	// components can nearly cancel between the F and G terms).
	dScale, hScale := 0.0, 0.0
	for i := 0; i < 3; i++ {
		dScale = math.Max(dScale, math.Abs(d[i]-Eps0*e[i]))
		hScale = math.Max(hScale, math.Abs(h[i]-b[i]/Mu0))
	}
	worst := 0.0
	for i := 0; i < 3; i++ {
		he, hb := 1e-6*SchwingerField, 1e-6*CriticalB
		ep, em := e, e
		ep[i] += he
		em[i] -= he
		dn := (nl(ep, b) - nl(em, b)) / (2 * he)
		bp, bm := b, b
		bp[i] += hb
		bm[i] -= hb
		hn := -(nl(e, bp) - nl(e, bm)) / (2 * hb)
		worst = math.Max(worst, math.Abs(dn-(d[i]-Eps0*e[i]))/dScale)
		worst = math.Max(worst, math.Abs(hn-(h[i]-b[i]/Mu0))/hScale)
	}
	if worst > 1e-6 {
		t.Errorf("D, H against ∂L: %.3g", worst)
	}
	// Small-signal tensors against differences of D, H.
	s, err := EHLinearize(e, b)
	if err != nil {
		t.Fatal(err)
	}
	worstT := 0.0
	for j := 0; j < 3; j++ {
		he, hb := 1e-5*SchwingerField, 1e-5*CriticalB
		ep, em := e, e
		ep[j] += he
		em[j] -= he
		dp, hp := EHConstitutive(ep, b)
		dm, hm := EHConstitutive(em, b)
		bp, bm := b, b
		bp[j] += hb
		bm[j] -= hb
		dbp, hbp := EHConstitutive(e, bp)
		dbm, hbm := EHConstitutive(e, bm)
		for i := 0; i < 3; i++ {
			lin := func(ii, jj int, eps float64) float64 { // the linear part of a diagonal
				if ii == jj {
					return eps
				}
				return 0
			}
			check := func(num, ana, lin, scale float64) {
				worstT = math.Max(worstT, math.Abs((num-lin)-(ana-lin))/scale)
			}
			nlScale := EHXi * SchwingerField * SchwingerField // size of the quartic part of ∂D/∂E
			check((dp[i]-dm[i])/(2*he), s.Eps[i][j], lin(i, j, Eps0), nlScale)
			check((hbp[i]-hbm[i])/(2*hb), s.Nu[i][j], lin(i, j, 1/Mu0), nlScale*SpeedOfLight*SpeedOfLight*SpeedOfLight*SpeedOfLight)
			check((dbp[i]-dbm[i])/(2*hb), s.Xeb[i][j], 0, nlScale*SpeedOfLight*SpeedOfLight)
			check((hp[i]-hm[i])/(2*he), s.Xhe[i][j], 0, nlScale*SpeedOfLight*SpeedOfLight)
		}
	}
	if worstT > 1e-6 {
		t.Errorf("small-signal tensors against differences: %.3g of the quartic scale", worstT)
	}
	// Lorentz invariance: L depends only on E² − c²B² and E·B.
	l0 := EHLagrangian(e, b)
	eb, bb := BoostVector(e, b, [3]float64{0.3 * SpeedOfLight, -0.5 * SpeedOfLight, 0.2 * SpeedOfLight})
	if r := math.Abs(nl(eb, bb)-nl(e, b)) / math.Abs(nl(e, b)); r > 1e-9 {
		t.Errorf("quartic part of L changes under a boost by %.3g", r)
	}
	if r := math.Abs(EHLagrangian(eb, bb)-l0) / (0.5 * Eps0 * dot(e, e)); r > 1e-12 {
		t.Errorf("L changes under a boost by %.3g", r)
	}
	// Energy density of a pure electric field: ε₀E²/2 + 3ξE⁴.
	e1 := [3]float64{0.05 * SchwingerField, 0, 0}
	u := EHEnergyDensity(e1, [3]float64{})
	want := 0.5*Eps0*e1[0]*e1[0] + 3*EHXi*math.Pow(e1[0], 4)
	if r := math.Abs(u-want) / (3 * EHXi * math.Pow(e1[0], 4)); r > 1e-6 {
		t.Errorf("energy density %.12g, want %.12g", u, want)
	}

	// Birefringence from the tensors: light along x across B₀ ẑ.
	b0 := 0.1 * CriticalB
	s0, err := EHLinearize([3]float64{}, [3]float64{0, 0, b0})
	if err != nil {
		t.Fatal(err)
	}
	a := EHBirefringenceParameter(b0)
	nPar := math.Sqrt(s0.Eps[2][2] / Eps0 / (Mu0 * s0.Nu[1][1]))  // e ∥ z, b ∥ y
	nPerp := math.Sqrt(s0.Eps[1][1] / Eps0 / (Mu0 * s0.Nu[2][2])) // e ∥ y, b ∥ z
	wPar, wPerp := EHIndices(b0)
	if math.Abs(nPar-wPar) > 100*a*a || math.Abs(nPerp-wPerp) > 100*a*a {
		t.Errorf("indices from the tensors n∥ − 1 = %.6g, n⊥ − 1 = %.6g; Adler %.6g, %.6g", nPar-1, nPerp-1, wPar-1, wPerp-1)
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if s0.Xeb[i][j] != 0 || s0.Xhe[i][j] != 0 || (i != j && (s0.Eps[i][j] != 0 || s0.Nu[i][j] != 0)) {
				t.Errorf("magnetic background along z: off-diagonal or magnetoelectric term at %d%d", i, j)
			}
		}
	}
	t.Logf("B₀ = 0.1 B_c (%.4g T): a = %.4g; from the tensors n∥ − 1 = %.6g (7a = %.6g), n⊥ − 1 = %.6g (4a = %.6g)",
		b0, a, nPar-1, 7*a, nPerp-1, 4*a)

	if _, err := EHLinearize([3]float64{0.2 * SchwingerField}, [3]float64{}); !errors.Is(err, ErrEHStrongField) {
		t.Errorf("0.2 E_s accepted: %v", err)
	}
}
