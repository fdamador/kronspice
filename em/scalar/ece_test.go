// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/em"
)

func vdiff(a, b [3]float64) float64 {
	return math.Max(math.Abs(a[0]-b[0]), math.Max(math.Abs(a[1]-b[1]), math.Abs(a[2]-b[2])))
}

func vnorm(a [3]float64) float64 { return math.Sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) }

// Quadratic potentials (central differences are exact on them).
var quadPot = Potentials{
	Phi: func(x [3]float64, t float64) float64 { return 5 + 3*x[0]*x[0] - x[1]*x[2] + 2*t*x[0] },
	A: func(x [3]float64, t float64) [3]float64 {
		return [3]float64{2*t + x[0]*x[1], x[2]*x[2] - t*x[0], 0.5*x[0]*x[1] + 1}
	},
}

var pt = [3]float64{0.3, -0.2, 0.7}

// With no spin connection all three conventions are Maxwell's, and Maxwell's
// is the analytic −∇φ − ∂A/∂t, ∇×A.
func TestECEZeroSpinIsMaxwell(t *testing.T) {
	p := quadPot.At(pt, 0.4, 1e-3, 1e-3)
	x, y, z, tt := pt[0], pt[1], pt[2], 0.4
	wantE := [3]float64{-(6*x + 2*tt) - 2, z + x, y}
	wantB := [3]float64{0.5*x - 2*z, -0.5 * y, -tt - x}
	if d := vdiff(p.EStd, wantE); d > 1e-9 {
		t.Errorf("E_std off by %.3g: %v vs %v", d, p.EStd, wantE)
	}
	if d := vdiff(p.BStd, wantB); d > 1e-9 {
		t.Errorf("B_std off by %.3g: %v vs %v", d, p.BStd, wantB)
	}
	if vdiff(p.EECE, p.EStd) != 0 || vdiff(p.EECE2, p.EStd) != 0 || vdiff(p.BECE, p.BStd) != 0 || vdiff(p.BECE2, p.BStd) != 0 {
		t.Errorf("ω = 0 but ECE fields differ from Maxwell's")
	}
	if math.Abs(p.DivBECE) > 1e-6 || vnorm(p.FaradayECE) > 1e-6 {
		t.Errorf("ω = 0 but ECE fields carry magnetic sources: div B %.3g, Faraday %.3g", p.DivBECE, vnorm(p.FaradayECE))
	}
}

// The book's free space (6.178–6.179): ω = A/W⁽⁰⁾, ω₀ = φ/(cW⁽⁰⁾), so κ₀ = κ = 0.
// Then ω ∥ A and cω₀A = ωφ, both spin-connection terms vanish identically,
// and the ECE fields are Maxwell's for any W⁽⁰⁾.
func TestECEFreeSpaceIsMaxwell(t *testing.T) {
	const w0 = 0.37
	p := quadPot
	p.Omega = func(x [3]float64, t float64) [3]float64 {
		a := quadPot.A(x, t)
		return [3]float64{a[0] / w0, a[1] / w0, a[2] / w0}
	}
	p.Omega0 = func(x [3]float64, t float64) float64 { return quadPot.Phi(x, t) / (em.SpeedOfLight * w0) }
	r := p.At(pt, 0.4, 1e-3, 1e-3)
	k0, k := p.Kappa(pt, 0.4, w0)
	scale := vnorm(r.EStd)
	t.Logf("κ₀ = %.2g, |κ| = %.2g; E_ECE − E = %.2g, E_ECE2 − E = %.2g of %.3g", k0, vnorm(k),
		vdiff(r.EECE, r.EStd), vdiff(r.EECE2, r.EStd), scale)
	if math.Abs(k0) > 1e-15 || vnorm(k) > 1e-12 {
		t.Errorf("free-space spin connection has κ₀ = %.3g, κ = %v", k0, k)
	}
	for _, d := range []float64{vdiff(r.EECE, r.EStd), vdiff(r.EECE2, r.EStd), vdiff(r.BECE, r.BStd), vdiff(r.BECE2, r.BStd)} {
		if d > 1e-12*scale {
			t.Errorf("free-space ECE field differs from Maxwell's by %.3g", d)
		}
	}
}

// A uniform vector spin connection ω₀ = 0, ω = const: B_ECE = ∇×A − ω×A has
// ∇·B = ω·(∇×A) (the identity ∇·(ω×A) = A·∇×ω − ω·∇×A), magnetic charge
// unless ω ⊥ ∇×A. ECE's Gauss constraint (4.213) is what forbids it.
func TestECEUniformSpinGivesMagneticCharge(t *testing.T) {
	om := [3]float64{0.8, -0.5, 1.1}
	p := quadPot
	p.Omega = func([3]float64, float64) [3]float64 { return om }
	r := p.At(pt, 0.4, 1e-3, 1e-3)
	want := om[0]*r.BStd[0] + om[1]*r.BStd[1] + om[2]*r.BStd[2]
	t.Logf("∇·B_ECE = %.6g, ω·∇×A = %.6g (ρ_m = %.3g A·m/m³ per unit)", r.DivBECE, want, want/em.Mu0)
	if math.Abs(r.DivBECE-want) > 1e-6*math.Abs(want) {
		t.Errorf("∇·B_ECE = %.9g, want ω·∇×A = %.9g", r.DivBECE, want)
	}
	// Faraday (4.214 with ω₀ = 0): ∇×E + ∂B/∂t = ∇×(ωφ) − ∂(ω×A)/∂t, and for
	// uniform ω that is ∇φ×ω − ω×∂A/∂t.
	x, y, z, tt := pt[0], pt[1], pt[2], 0.4
	gphi := [3]float64{6*x + 2*tt, -z, -y}
	dAdt := [3]float64{2, -x, 0}
	f := cross(gphi, om)
	g := cross(om, dAdt)
	wantF := [3]float64{f[0] - g[0], f[1] - g[1], f[2] - g[2]}
	if d := vdiff(r.FaradayECE, wantF); d > 1e-6*vnorm(wantF) {
		t.Errorf("Faraday residual %v, want %v (4.214)", r.FaradayECE, wantF)
	}
}

// The slides' split (slide 17) cω₀A = ωφ = ½(∇φ − ∂A/∂t) satisfies the
// electric antisymmetry law (5.14) only where ∇φ − ∂A/∂t ∥ A; there the
// spin-connection terms of E cancel and E is Maxwell's.
func TestECEAntisymmetrySplit(t *testing.T) {
	cl := em.SpeedOfLight
	phi := func(x [3]float64, t float64) float64 { return 4 + 3*x[0]*x[0] }
	g := func(x [3]float64, t float64) [3]float64 { return [3]float64{6*x[0] - 2, 0, 0} } // ∇φ − ∂A/∂t
	mk := func(axis int) Potentials {
		a := func(x [3]float64, t float64) [3]float64 {
			v := [3]float64{}
			v[axis] = 2*t + 1 + x[1]*x[1]
			return v
		}
		// A along x gives ∂A/∂t = (2,0,0); along y the gradient is unchanged
		// and ∂A/∂t moves to y, so g is recomputed for that case below.
		p := Potentials{Phi: phi, A: a}
		gf := g
		if axis != 0 {
			gf = func(x [3]float64, t float64) [3]float64 { return [3]float64{6 * x[0], -2, 0} }
		}
		p.Omega = func(x [3]float64, t float64) [3]float64 {
			v, f := gf(x, t), phi(x, t)
			return [3]float64{v[0] / (2 * f), v[1] / (2 * f), v[2] / (2 * f)}
		}
		p.Omega0 = func(x [3]float64, t float64) float64 {
			v, av := gf(x, t), a(x, t)
			return (v[0]*av[0] + v[1]*av[1] + v[2]*av[2]) / (2 * cl * (av[0]*av[0] + av[1]*av[1] + av[2]*av[2]))
		}
		return p
	}
	par := mk(0).At(pt, 0.4, 1e-3, 1e-3)
	perp := mk(1).At(pt, 0.4, 1e-3, 1e-3)
	t.Logf("∥: antisymmetry residual %.2g, E_ECE − E %.2g;  ⊥: residual %.3g, E_ECE − E %.3g (|E| %.3g)",
		vnorm(par.ElecAnti), vdiff(par.EECE, par.EStd), vnorm(perp.ElecAnti), vdiff(perp.EECE, perp.EStd), vnorm(perp.EStd))
	if vnorm(par.ElecAnti) > 1e-9 || vdiff(par.EECE, par.EStd) > 1e-9 {
		t.Errorf("parallel case: antisymmetry %.3g, E differs by %.3g", vnorm(par.ElecAnti), vdiff(par.EECE, par.EStd))
	}
	if vnorm(perp.ElecAnti) < 0.1 {
		t.Errorf("perpendicular case should break the antisymmetry law, residual %.3g", vnorm(perp.ElecAnti))
	}
	// The spin connection read back from the potentials is the one declared.
	if d := vdiff(par.SpinFromPot, mk(0).Omega(pt, 0.4)); d > 1e-9 {
		t.Errorf("slide-17 ω read back off by %.3g", d)
	}
}
