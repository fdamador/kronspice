// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/hyper/quat"
)

// Checks of A. Waser, "On the Notation of Maxwell's Field Equations"
// (AW-Verlag, 2000) against the engine's equations. Comparison and
// verdicts: doc/waser-2000-comparison.md. Each test states which of
// Waser's equations it checks; the engine's form must hold and, where
// Waser's form differs, the difference must show up as a non-zero
// residual on a field for which the engine's form is exact.

// Waser (1.32): from his Gauss law (1.28) ∇·D = −ρ and continuity (1.8)
// −∂ρ/∂t = ∇·j he derives ∂ρ/∂t = 0. With the engine's ∇·D = +ρ the
// divergence of Ampère–Maxwell is the identity ∂ρ/∂t + ∇·J = 0 for any
// fields, and charge can change (a charging capacitor). Here D(t) has a
// growing charge density, J is the current Ampère–Maxwell requires, and:
// the engine's continuity holds to round-off; with Waser's ρ_W = −∇·D his
// continuity (1.8) fails by 2·∂ρ_W/∂t, i.e. it holds only if charge never
// changes.
func TestWaserContinuity(t *testing.T) {
	g := cubeGrid(13, 0.05)
	// D = (1 + s·t)·D0(r) with D0 = r·exp(−r²/0.04) C/m², so ∂D/∂t = s·D0
	// and the charge density ∇·D grows at s·∇·D0. Amplitudes are chosen so
	// ∂ρ/∂t and ∇×H are of the same order (round-off comparable).
	const s = 2.0
	envelope := func(x, y, z float64) float64 { return math.Exp(-(x*x + y*y + z*z) / 0.04) }
	dDdt := SampleV(g, func(x, y, z float64) [3]float64 {
		e := s * envelope(x, y, z)
		return [3]float64{x * e, y * e, z * e}
	})
	H := SampleV(g, func(x, y, z float64) [3]float64 {
		return [3]float64{math.Sin(3 * y), math.Cos(2*z) * x, math.Sin(x*z + y)}
	})
	cH := Curl(H)
	J := VField{G: g, V: make([][3]float64, g.Len())}
	for n := range J.V {
		for a := 0; a < 3; a++ {
			J.V[n][a] = cH.V[n][a] - dDdt.V[n][a] // Ampère–Maxwell (1.26) solved for J
		}
	}
	divJ := Div(J)
	dRhoDt := Div(dDdt) // ∂(∇·D)/∂t, engine's ρ = +∇·D
	var engine, waser, scale float64
	for n := 0; n < g.Len(); n++ {
		if math.IsNaN(divJ[n]) || math.IsNaN(dRhoDt[n]) {
			continue // stencil does not fit (two-stage derivative)
		}
		scale = math.Max(scale, math.Abs(dRhoDt[n]))
		engine = math.Max(engine, math.Abs(dRhoDt[n]+divJ[n]))
		// Waser: ρ_W = −∇·D, continuity (1.8) reads ∂ρ_W/∂t + ∇·j = 0.
		waser = math.Max(waser, math.Abs(-dRhoDt[n]+divJ[n]))
	}
	if scale == 0 {
		t.Fatal("no interior points")
	}
	if engine > 1e-9*scale {
		t.Errorf("engine continuity residual %g of %g", engine, scale)
	}
	if waser < 1.9*scale {
		t.Errorf("Waser's signs should leave a residual of 2·∂ρ/∂t (%g), got %g", 2*scale, waser)
	}
	t.Logf("∂ρ/∂t up to %.3g C/m³/s: engine continuity residual %.2g; Waser (1.28)+(1.8) residual %.3g = 2·∂ρ/∂t",
		scale, engine, waser)
}

// Waser (1.24) and the labels of (1.26)–(1.27): a vacuum plane wave
// E = x̂ E0 cos(kz − ωt), B = ŷ (E0/c) cos(kz − ωt) obeys Faraday's law
// ∇×E = −∂B/∂t (engine; Treatise (B) with 𝔅 = V.∇𝔄); Waser's (1.24)
// ∇×E = +μ ∂H/∂t fails by 2·∂B/∂t.
func TestWaserFaradaySign(t *testing.T) {
	g := cubeGrid(11, 0.02)
	const e0, k = 1.0, 7.0
	w := k * SpeedOfLight
	t0 := 0.1 / w
	E := SampleV(g, func(x, y, z float64) [3]float64 { return [3]float64{e0 * math.Cos(k*z-w*t0), 0, 0} })
	dBdt := SampleV(g, func(x, y, z float64) [3]float64 {
		return [3]float64{0, e0 / SpeedOfLight * w * math.Sin(k*z-w*t0), 0}
	})
	cE := Curl(E)
	var engine, waser, scale float64
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		for a := 0; a < 3; a++ {
			scale = math.Max(scale, math.Abs(dBdt.V[n][a]))
			engine = math.Max(engine, math.Abs(cE.V[n][a]+dBdt.V[n][a]))
			waser = math.Max(waser, math.Abs(cE.V[n][a]-dBdt.V[n][a]))
		}
	}
	// Fourth-order stencil error ~ (kh)⁴/30 relative.
	if engine > 1e-4*scale {
		t.Errorf("Faraday residual %g of %g", engine, scale)
	}
	if waser < 1.9*scale {
		t.Errorf("Waser (1.24) should fail by 2·∂B/∂t (%g), got %g", 2*scale, waser)
	}
}

// Waser (1.18) e = S.∇D and (1.28) ∇·D = −ρ against the measured sign of
// Coulomb's law. A uniformly charged ball of positive charge density ρ0
// has D = ρ0·r/3 inside, pointing outward. Hamilton's ∇ gives
// S.∇D = −∇·D = −ρ0 (the engine's convention, em/doc.go), so e = S.∇D
// read literally gives e = −ρ0, the sign of Maxwell's 1865 equation (G);
// the engine's quaternion residual uses ρ = −S.∇D = +∇·D and vanishes
// with ρ = +ρ0. Settling which sign Maxwell printed in Treatise art. 619
// needs the page image (CLAUDE.md §3c); this test pins the engine's sign.
func TestWaserGaussSign(t *testing.T) {
	g := cubeGrid(9, 0.05)
	const rho0 = 1e-6 // C/m³
	E := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(x, y, z).Scale(rho0 / (3 * Eps0)) })
	zero := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Q{} })
	nE := Nabla(E)
	rho := make([]float64, g.Len())
	for n := range rho {
		rho[n] = rho0
	}
	r := MaxwellQuaternion(E, zero, zero, zero, rho, nil)
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		if s := nE.Q[n].W * Eps0; math.Abs(s+rho0) > 1e-9*rho0 {
			t.Fatalf("S.∇D = %g, want −ρ0 = %g", s, -rho0)
		}
		if res := r.R1.Q[n].W * Eps0; math.Abs(res) > 1e-9*rho0 {
			t.Fatalf("engine Gauss residual %g with ρ = +ρ0", res)
		}
	}
	rho2 := make([]float64, g.Len())
	for n := range rho2 {
		rho2[n] = -rho0 // Waser (1.28)
	}
	r2 := MaxwellQuaternion(E, zero, zero, zero, rho2, nil)
	for n := 0; n < g.Len(); n++ {
		if g.Interior(n) && math.Abs(r2.R1.Q[n].W*Eps0+2*rho0) > 1e-9*rho0 {
			t.Fatalf("Waser sign: residual %g, want −2ρ0", r2.R1.Q[n].W*Eps0)
		}
	}
}

// Waser (1.12) B = H + 4πM and (1.19) m = S.∇M: the "magnetic mass" is
// the pole density of magnetization, not a free monopole. In SI,
// B = μ0(H + M): with B = ∇×A (Treatise (A)) and any magnetization M,
// ∇·B = 0 while ∇·H = −∇·M ≠ 0 — poles in H, none in B.
func TestWaserMagneticMass(t *testing.T) {
	g := cubeGrid(13, 0.05)
	A := SampleV(g, func(x, y, z float64) [3]float64 {
		return [3]float64{math.Sin(2 * y), math.Cos(3*z) * x, math.Sin(x + y*z)}
	})
	B := Curl(A)
	M := SampleV(g, func(x, y, z float64) [3]float64 {
		return [3]float64{1e5 * math.Exp(-(x*x+y*y+z*z)/0.03), 0, 0} // a magnetized blob along x [A/m]
	})
	H := VField{G: g, V: make([][3]float64, g.Len())}
	for n := range H.V {
		for a := 0; a < 3; a++ {
			H.V[n][a] = B.V[n][a]/Mu0 - M.V[n][a]
		}
	}
	divB, divH, divM := Div(B), Div(H), Div(M)
	var maxB, maxPole float64
	for n := 0; n < g.Len(); n++ {
		if math.IsNaN(divB[n]) || math.IsNaN(divH[n]) {
			continue
		}
		maxB = math.Max(maxB, math.Abs(divB[n]))
		maxPole = math.Max(maxPole, math.Abs(divM[n]))
		if math.Abs(divH[n]+divM[n]) > 1e-6*math.Max(1, math.Abs(divM[n])) {
			t.Fatalf("∇·H = %g, −∇·M = %g", divH[n], -divM[n])
		}
	}
	if maxB > 1e-9 {
		t.Errorf("∇·B = %g, want 0 (no free monopoles)", maxB)
	}
	if maxPole < 1e5 {
		t.Errorf("pole density only %g; test field too weak", maxPole)
	}
}
