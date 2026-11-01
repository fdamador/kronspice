// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/hyper/quat"
)

func cubeGrid(n int, h float64) Grid {
	o := -float64(n-1) * h / 2
	return Grid{N: [3]int{n, n, n}, H: [3]float64{h, h, h}, Origin: [3]float64{o, o, o}}
}

// fieldAt samples a time-dependent vector field at t and its fourth-order
// time derivative from snapshots at t ± δ, t ± 2δ.
func fieldAt(g Grid, f func(x, y, z, t float64) quat.Q, t, dt float64) (QField, QField) {
	at := func(tt float64) QField {
		return SampleQ(g, func(x, y, z float64) quat.Q { return f(x, y, z, tt) })
	}
	d, _ := TimeDerivative(at(t-2*dt), at(t-dt), at(t+dt), at(t+2*dt), dt)
	return at(t), d
}

// maxInterior returns the largest |x| over interior points of a residual
// accessor.
func maxInterior(g Grid, val func(n int) float64) float64 {
	m := 0.0
	for n := 0; n < g.Len(); n++ {
		if g.Interior(n) {
			m = math.Max(m, math.Abs(val(n)))
		}
	}
	return m
}

// Triple twin (CLAUDE.md §3c rule 2): on arbitrary smooth fields and
// sources — not solutions, so every residual is non-zero — the Treatise
// quaternion form, the spacetime-algebra form ∇F = J/(ε₀c) and plain
// vector calculus give the same four residuals to 1e−9 relative (they
// differ only in round-off).
func TestTripleTwin(t *testing.T) {
	r := rand.New(rand.NewSource(51))
	g := cubeGrid(9, 0.07)
	sta := NewSTA()
	for trial := 0; trial < 5; trial++ {
		var a [12]float64
		for k := range a {
			a[k] = r.NormFloat64()
		}
		field := func(s float64) func(x, y, z float64) quat.Q {
			return func(x, y, z float64) quat.Q {
				return quat.Vec(s*math.Sin(a[0]*x+a[1]*y+a[2]), s*math.Cos(a[3]*y+a[4]*z)*(1+a[5]*x), s*math.Sin(a[6]*x*z+a[7]*y))
			}
		}
		E, B := SampleQ(g, field(100)), SampleQ(g, field(3e-7))
		dE, dB := SampleQ(g, field(2e9)), SampleQ(g, field(-50))
		rho := SampleS(g, func(x, y, z float64) float64 { return 1e-9 * (a[8] + x*y + a[9]*z) })
		J := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(a[10]+x, y*z, a[11]*z*z) })
		Jv := ToV(J)
		vr := MaxwellVector(ToV(E), ToV(B), ToV(dE), ToV(dB), rho, &Jv)
		qr := MaxwellQuaternion(E, B, dE, dB, rho, &J)
		sr, ok := sta.MaxwellSpacetime(E, B, dE, dB, rho, &J)
		scale := [4]float64{}
		for n := 0; n < g.Len(); n++ {
			if !g.Interior(n) {
				continue
			}
			scale[0] = math.Max(scale[0], math.Abs(vr.Gauss[n]))
			scale[2] = math.Max(scale[2], math.Abs(vr.DivB[n]))
			for k := 0; k < 3; k++ {
				scale[1] = math.Max(scale[1], math.Abs(vr.Faraday[n][k]))
				scale[3] = math.Max(scale[3], math.Abs(vr.Ampere[n][k]))
			}
		}
		for n := 0; n < g.Len(); n++ {
			if !g.Interior(n) {
				if ok[n] {
					t.Fatal("spacetime residual at an edge point")
				}
				continue
			}
			q1, q2 := qr.R1.Q[n], qr.R2.Q[n]
			gs, fs, ds, as := sta.SpacetimeParts(sr[n])
			faraday := [3]float64{q1.X, q1.Y, q1.Z}
			ampere := [3]float64{q2.X, q2.Y, q2.Z}
			check := func(name string, i int, vec, quat, st float64) {
				if math.Abs(vec-quat) > 1e-9*scale[i] || math.Abs(vec-st) > 1e-9*scale[i] {
					t.Fatalf("%s at %d: vector %.12g, quaternion %.12g, spacetime %.12g", name, n, vec, quat, st)
				}
			}
			check("Gauss", 0, vr.Gauss[n], -q1.W, gs)
			check("div B", 2, vr.DivB[n], -q2.W, ds)
			for k := 0; k < 3; k++ {
				check("Faraday", 1, vr.Faraday[n][k], faraday[k], fs[k])
				check("Ampère", 3, vr.Ampere[n][k], ampere[k], as[k])
			}
		}
		if scale[0] == 0 || scale[1] == 0 || scale[2] == 0 || scale[3] == 0 {
			t.Fatalf("a residual is identically zero: %v", scale)
		}
	}
}

// Hamilton's ∇ on a vector field: scalar part −div, vector part curl;
// on a scalar field, grad — exact (to round-off) on polynomials of degree
// ≤ 4, which the fourth-order stencil differentiates exactly.
func TestNablaPolynomials(t *testing.T) {
	g := cubeGrid(7, 0.3)
	v := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(x*x*y, y*z*z*z, x*y*z+x*x*x*x) })
	nv := Nabla(v)
	psi := SampleS(g, func(x, y, z float64) float64 { return x*x*x*y + z*z })
	gp := Nabla(Scalars(g, psi))
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		p := g.Point(g.IJK(n))
		x, y, z := p[0], p[1], p[2]
		div := 2*x*y + z*z*z + x*y
		curl := [3]float64{x*z - 3*y*z*z, -(y*z + 4*x*x*x), -x * x}
		got := nv.Q[n]
		if math.Abs(-got.W-div) > 1e-12 || math.Abs(got.X-curl[0]) > 1e-12 || math.Abs(got.Y-curl[1]) > 1e-12 || math.Abs(got.Z-curl[2]) > 1e-12 {
			t.Fatalf("∇v at (%g,%g,%g) = %v, want −div %g, curl %v", x, y, z, got, div, curl)
		}
		gr := gp.Q[n]
		if math.Abs(gr.X-3*x*x*y) > 1e-12 || math.Abs(gr.Y-x*x*x) > 1e-12 || math.Abs(gr.Z-2*z) > 1e-12 || gr.W != 0 {
			t.Fatalf("∇ψ = %v", gr)
		}
	}
}

// Analytic cases with sources, exact for the stencil: inside a uniformly
// charged sphere E = ρr/(3ε₀) (div E = ρ/ε₀, curl E = 0); inside a wire of
// uniform current density J along z, B = μ₀J(−y, x, 0)/2 (curl B = μ₀J,
// div B = 0).
func TestStaticSourcesExact(t *testing.T) {
	g := cubeGrid(9, 1e-3)
	const rho0, j0 = 1e-6, 1e6
	E := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(x, y, z).Scale(rho0 / (3 * Eps0)) })
	B := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(-y, x, 0).Scale(Mu0 * j0 / 2) })
	zero := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Q{} })
	rho := SampleS(g, func(x, y, z float64) float64 { return rho0 })
	J := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(0, 0, j0) })
	r := MaxwellQuaternion(E, B, zero, zero, rho, &J)
	if m := maxInterior(g, func(n int) float64 { return r.R1.Q[n].Abs() }); m > 1e-9*rho0/Eps0 {
		t.Errorf("charged sphere: |R1| = %g (scale %g)", m, rho0/Eps0)
	}
	if m := maxInterior(g, func(n int) float64 { return r.R2.Q[n].Abs() }); m > 1e-9*Mu0*j0 {
		t.Errorf("current cylinder: |R2| = %g (scale %g)", m, Mu0*j0)
	}
}

// A plane wave in free space satisfies ∇F = 0: the residuals are the
// stencil's truncation error, falling as h⁴ (16× per halving), below 1e−6
// of k·|E| on a grid of 20 points per wavelength. Direction and
// polarization are general.
func TestPlaneWave(t *testing.T) {
	kdir := [3]float64{1, 2, 2}
	pol := [3]float64{2, -1, 0} // ⟂ kdir
	norm := func(v [3]float64) [3]float64 {
		l := math.Sqrt(dot(v, v))
		return [3]float64{v[0] / l, v[1] / l, v[2] / l}
	}
	kh, eh := norm(kdir), norm(pol)
	bh := cross(kh, eh)
	const lambda, e0 = 1.0, 100.0
	k := 2 * math.Pi / lambda
	w := k * SpeedOfLight
	wave := func(amp float64, dir [3]float64) func(x, y, z, t float64) quat.Q {
		return func(x, y, z, t float64) quat.Q {
			c := amp * math.Cos(k*(kh[0]*x+kh[1]*y+kh[2]*z)-w*t)
			return quat.Vec(c*dir[0], c*dir[1], c*dir[2])
		}
	}
	sta := NewSTA()
	var prev float64
	for _, ppl := range []float64{10, 20, 40} {
		h := lambda / ppl
		g := cubeGrid(7, h)
		dt := h / SpeedOfLight
		E, dE := fieldAt(g, wave(e0, eh), 0.1/w, dt)
		B, dB := fieldAt(g, wave(e0/SpeedOfLight, bh), 0.1/w, dt)
		qr := MaxwellQuaternion(E, B, dE, dB, nil, nil)
		m1 := maxInterior(g, func(n int) float64 { return qr.R1.Q[n].Abs() }) / (k * e0)
		m2 := maxInterior(g, func(n int) float64 { return qr.R2.Q[n].Abs() }) / (k * e0 / SpeedOfLight)
		sr, ok := sta.MaxwellSpacetime(E, B, dE, dB, nil, nil)
		ms := 0.0
		for n := range sr {
			if ok[n] {
				ms = math.Max(ms, sr[n].MaxAbs())
			}
		}
		ms /= k * e0
		rel := math.Max(m1, m2)
		t.Logf("%g points/λ: quaternion residual %.3g, spacetime %.3g (relative to k·E₀)", ppl, rel, ms)
		if prev > 0 && prev/rel < 12 {
			t.Errorf("convergence ratio %g at %g points/λ, want ≈16 (fourth order)", prev/rel, ppl)
		}
		if ppl == 20 && rel > 1e-3 {
			t.Errorf("plane wave residual %g at 20 points/λ", rel)
		}
		prev = rel
	}
	if prev > 1e-5 {
		t.Errorf("plane wave residual %g at 40 points/λ", prev)
	}
}

// The Treatise relations 𝔅 = V.∇𝔄 and 𝔈 = −∂𝔄/∂t − ∇ψ: the Lorenz-gauge
// potentials of the plane wave give its fields, and a gauge change
// A → A + ∇χ, ψ → ψ − ∂χ/∂t (χ a polynomial, differentiated exactly)
// leaves E and B unchanged.
func TestPotentials(t *testing.T) {
	const lambda, e0 = 1.0, 100.0
	k := 2 * math.Pi / lambda
	w := k * SpeedOfLight
	g := cubeGrid(9, lambda/40)
	dt := g.H[0] / SpeedOfLight
	t0 := 0.3 / w
	// E = E₀x̂ cos(kz − ωt), B = (E₀/c)ŷ cos(kz − ωt) from
	// A = (E₀/ω)x̂ sin(kz − ωt), ψ = 0.
	A, dA := fieldAt(g, func(x, y, z, tt float64) quat.Q {
		return quat.Vec(e0/w*math.Sin(k*z-w*tt), 0, 0)
	}, t0, dt)
	psi := make([]float64, g.Len())
	E, B := FieldsFromPotentials(psi, A, dA)
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		p := g.Point(g.IJK(n))
		c := math.Cos(k*p[2] - w*t0)
		// Truncation error of the fourth-order stencils at 40 points/λ is
		// (2π/40)⁴/30 ≈ 2e−5 of the amplitude.
		if math.Abs(E.Q[n].X-e0*c) > 1e-4*e0 || math.Abs(E.Q[n].Y)+math.Abs(E.Q[n].Z) > 1e-9*e0 ||
			math.Abs(B.Q[n].Y-e0/SpeedOfLight*c) > 1e-4*e0/SpeedOfLight {
			t.Fatalf("fields from potentials at z = %g: E %v, B %v", p[2], E.Q[n], B.Q[n])
		}
	}
	// Gauge change with χ = a·x²y + b·z³ − c·t·x: ∇χ and ∂χ/∂t exact.
	chi := func(x, y, z, tt float64) float64 { return 2*x*x*y + 0.5*z*z*z - 3*tt*x }
	A2, dA2 := fieldAt(g, func(x, y, z, tt float64) quat.Q {
		a := quat.Vec(e0/w*math.Sin(k*z-w*tt), 0, 0)
		return a.Add(quat.Vec(4*x*y-3*tt, 2*x*x, 1.5*z*z))
	}, t0, dt)
	psi2 := SampleS(g, func(x, y, z float64) float64 { return -(chi(x, y, z, t0+1e-9) - chi(x, y, z, t0-1e-9)) / 2e-9 })
	E2, B2 := FieldsFromPotentials(psi2, A2, dA2)
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		if E2.Q[n].Sub(E.Q[n]).Abs() > 1e-6*e0 || B2.Q[n].Sub(B.Q[n]).Abs() > 1e-9*e0/SpeedOfLight {
			t.Fatalf("gauge change altered the fields: ΔE %v, ΔB %v", E2.Q[n].Sub(E.Q[n]), B2.Q[n].Sub(B.Q[n]))
		}
	}
}

// Parity and rotations (CLAUDE.md §3c rule 3): the field of a current
// segment (Biot–Savart, closed form) recomputed for the transformed
// geometry equals both the vector rule (E → ME, B → det(M)·MB) and the GA
// rule (B as the bivector I·B, sandwiched by reflections), for proper and
// improper orthogonal M, including the inversion −1 (E flips, B does not).
func TestParityAndRotations(t *testing.T) {
	segB := func(a, b, p [3]float64, cur float64) [3]float64 {
		// B = μ₀I/(4πd)·(cos θ₂ − cos θ₁)·ŝ×d̂ for the segment a→b.
		l := [3]float64{b[0] - a[0], b[1] - a[1], b[2] - a[2]}
		ll := math.Sqrt(dot(l, l))
		u := [3]float64{l[0] / ll, l[1] / ll, l[2] / ll}
		ap := [3]float64{p[0] - a[0], p[1] - a[1], p[2] - a[2]}
		bp := [3]float64{p[0] - b[0], p[1] - b[1], p[2] - b[2]}
		perp := [3]float64{ap[0] - dot(ap, u)*u[0], ap[1] - dot(ap, u)*u[1], ap[2] - dot(ap, u)*u[2]}
		d := math.Sqrt(dot(perp, perp))
		c1 := dot(ap, u) / math.Sqrt(dot(ap, ap))
		c2 := dot(bp, u) / math.Sqrt(dot(bp, bp))
		dir := cross(u, [3]float64{perp[0] / d, perp[1] / d, perp[2] / d})
		mag := Mu0 * cur / (4 * math.Pi * d) * (c1 - c2)
		return [3]float64{mag * dir[0], mag * dir[1], mag * dir[2]}
	}
	pointE := func(q, src, p [3]float64) [3]float64 {
		r := [3]float64{p[0] - src[0], p[1] - src[1], p[2] - src[2]}
		r3 := math.Pow(dot(r, r), 1.5)
		s := q[0] / (4 * math.Pi * Eps0 * r3)
		return [3]float64{s * r[0], s * r[1], s * r[2]}
	}
	rng := rand.New(rand.NewSource(52))
	ortho := func(improper bool) [3][3]float64 {
		// Random rotation via a random unit quaternion, optionally × −1.
		q, _ := quat.Q{W: rng.NormFloat64(), X: rng.NormFloat64(), Y: rng.NormFloat64(), Z: rng.NormFloat64()}.Unit()
		m := q.Matrix()
		if improper {
			for i := range m {
				for j := range m[i] {
					m[i][j] = -m[i][j]
				}
			}
		}
		return m
	}
	mats := [][3][3]float64{{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}, {{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}
	for n := 0; n < 50; n++ {
		mats = append(mats, ortho(n%2 == 1))
	}
	for _, m := range mats {
		a, b := [3]float64{0.1, -0.2, 0.3}, [3]float64{-0.4, 0.5, 0.2}
		src := [3]float64{0.3, 0.1, -0.2}
		p := [3]float64{1.1, 0.7, -0.5}
		e := pointE([3]float64{1e-9}, src, p)
		bf := segB(a, b, p, 2)
		// Recompute for the transformed geometry.
		eT := pointE([3]float64{1e-9}, mulVec(m, src), mulVec(m, p))
		bT := segB(mulVec(m, a), mulVec(m, b), mulVec(m, p), 2)
		ev, bv := TransformVector(m, e, bf)
		eg, bg := TransformGA(m, e, bf)
		for k := 0; k < 3; k++ {
			if math.Abs(ev[k]-eT[k]) > 1e-12*norm3(e) || math.Abs(eg[k]-eT[k]) > 1e-12*norm3(e) ||
				math.Abs(bv[k]-bT[k]) > 1e-12*norm3(bf) || math.Abs(bg[k]-bT[k]) > 1e-12*norm3(bf) {
				t.Fatalf("det %g: E recomputed %v, vector %v, GA %v; B recomputed %v, vector %v, GA %v",
					det3(m), eT, ev, eg, bT, bv, bg)
			}
		}
	}
}

func norm3(v [3]float64) float64 { return math.Sqrt(dot(v, v)) }

// Moving observer (CLAUDE.md §3c rule 3): the spacetime-algebra boost
// R̃FR equals the textbook transformation for random fields and velocities
// up to 0.9c; a probe moving past a static charge sees B′ = −γ·v×E/c²
// (perpendicular part), which for v ≪ c is B − v×E/c²; and E² − c²B² and
// E·B are invariant.
func TestMovingObserver(t *testing.T) {
	sta := NewSTA()
	r := rand.New(rand.NewSource(53))
	c := SpeedOfLight
	for n := 0; n < 200; n++ {
		e := [3]float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		b := [3]float64{r.NormFloat64() / c, r.NormFloat64() / c, r.NormFloat64() / c}
		vv := [3]float64{r.NormFloat64(), r.NormFloat64(), r.NormFloat64()}
		s := 0.9 * c * r.Float64() / norm3(vv)
		v := [3]float64{vv[0] * s, vv[1] * s, vv[2] * s}
		et, bt := BoostVector(e, b, v)
		eg, bg := sta.Boost(e, b, v)
		for k := 0; k < 3; k++ {
			if math.Abs(et[k]-eg[k]) > 1e-12*10 || math.Abs(bt[k]-bg[k]) > 1e-12*10/c {
				t.Fatalf("boost by %v: textbook E′ %v B′ %v, GA E′ %v B′ %v", v, et, bt, eg, bg)
			}
		}
		inv1 := func(e, b [3]float64) float64 { return dot(e, e) - c*c*dot(b, b) }
		if math.Abs(inv1(eg, bg)-inv1(e, b)) > 1e-10*(1+math.Abs(inv1(e, b))) || math.Abs(dot(eg, bg)-dot(e, b)) > 1e-10/c {
			t.Fatalf("invariants changed")
		}
	}
	// Static charge, slow probe: B′ ≈ −v×E/c².
	e := [3]float64{0, 1e3, 0}
	v := [3]float64{30, 0, 0}
	_, bg := sta.Boost(e, [3]float64{}, v)
	want := cross(v, e)
	for k := range want {
		want[k] *= -1 / (c * c)
	}
	if norm3([3]float64{bg[0] - want[0], bg[1] - want[1], bg[2] - want[2]}) > 1e-12*norm3(want) {
		t.Errorf("slow probe past a charge: B′ = %v, want −v×E/c² = %v", bg, want)
	}
}

func BenchmarkNabla32(b *testing.B) {
	g := cubeGrid(32, 0.01)
	f := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(math.Sin(x), y*z, math.Cos(z)) })
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		Nabla(f)
	}
}

func BenchmarkMaxwellSpacetime16(b *testing.B) {
	g := cubeGrid(16, 0.01)
	f := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(math.Sin(x), y*z, math.Cos(z)) })
	sta := NewSTA()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		sta.MaxwellSpacetime(f, f, f, f, nil, nil)
	}
}
