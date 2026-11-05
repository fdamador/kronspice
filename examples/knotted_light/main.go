// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command knotted_light tests the first step of hypothesis H29 (CLAUDE.md
// §1): particles as knots of field energy. If a particle were a knotted
// field, its stability could come from topology, since a knot cannot untie
// without being cut. The question here is whether topology alone, in
// ordinary (linear) Maxwell electrodynamics, holds a knot of light
// together.
//
// The field is Rañada's knotted light, the hopfion: an exact, finite-energy
// solution of Maxwell's equations in vacuum whose electric and magnetic
// field lines are closed loops, every pair of them linked once (the Hopf
// fibration). A. F. Rañada, Lett. Math. Phys. 18, 97 (1989); W. T. M.
// Irvine, D. Bouwmeester, Nature Phys. 4, 716 (2008). It is written in the
// Bateman form (H. Kedia, I. Bialynicki-Birula, D. Peralta-Salas,
// W. T. M. Irvine, PRL 111, 150404, 2013), with lengths in units of L0 and
// times in L0/c:
//
//	α = (r² − t² − 1 + 2iz) / (r² − (t − i)²)
//	β = 2(x − iy) / (r² − (t − i)²)
//	F = E + icB = ∇α × ∇β = ∇ × (α∇β)
//
// so E = curl Re(α∇β) and cB = curl Im(α∇β): the field is set in the solver
// as discrete curls of these potentials, which makes div E = div B = 0
// exactly on the grid.
//
// What is measured, each against the exact solution at the same time:
//
//  1. Energy: the fraction within 2 L0 of the centre, from em/fdtd.go and
//     from the exact field, with the solver's energy ledger.
//  2. Topology: two magnetic field lines are traced (from the solver's
//     field and from the exact field) and their Gauss linking number is
//     computed. The Hopf fibration has linking number 1 for every pair.
//
// Standard prediction: the linking is conserved (Maxwell's equations carry
// the topology along, Rañada 1989; helicity is conserved), but nothing in
// linear Maxwell holds the energy in place: the knot expands at the speed
// of light while staying linked. A stable knot needs a nonlinear term
// (Skyrme 1961; Faddeev & Niemi, Nature 387, 58, 1997), which Derrick's
// theorem (J. Math. Phys. 5, 1252, 1964) shows is required for any static
// lump in three dimensions.
package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// Hopfion is the exact field in units of L0 (length) and L0/c (time).
// F returns F = E + icB (dimensionless); P returns the potential α∇β.
type Hopfion struct{}

func (Hopfion) parts(x, y, z, t float64) (alpha complex128, gA, gB [3]complex128) {
	r2 := x*x + y*y + z*z
	tm := complex(t, -1)
	D := complex(r2, 0) - tm*tm
	N := complex(r2-t*t-1, 2*z)
	alpha = N / D
	M := complex(2*x, -2*y)
	beta := M / D
	_ = beta
	p := [3]float64{x, y, z}
	for j := 0; j < 3; j++ {
		dD := complex(2*p[j], 0)
		dN := dD
		if j == 2 {
			dN += 2i
		}
		var dM complex128
		switch j {
		case 0:
			dM = 2
		case 1:
			dM = -2i
		}
		gA[j] = dN/D - N*dD/(D*D)
		gB[j] = dM/D - M*dD/(D*D)
	}
	return
}

// F is ∇α × ∇β at (x, y, z, t).
func (h Hopfion) F(x, y, z, t float64) [3]complex128 {
	_, a, b := h.parts(x, y, z, t)
	return [3]complex128{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

// P is the potential α∇β, with F = ∇ × P.
func (h Hopfion) P(x, y, z, t float64) [3]complex128 {
	al, _, b := h.parts(x, y, z, t)
	return [3]complex128{al * b[0], al * b[1], al * b[2]}
}

// Sign is +1 when E = Re F, cB = Im F satisfy Maxwell's equations
// (i∂F/∂t = ∇ × F) and −1 when the conjugate does (then cB = −Im F). It is
// checked numerically at a few points, so the helicity convention of the
// formula cannot be taken on trust.
func (h Hopfion) Sign() float64 {
	const e = 1e-5
	var plus, minus float64
	for _, p := range [][4]float64{{0.3, 0.2, -0.1, 0.4}, {-0.5, 0.7, 0.3, -0.2}, {0.1, -0.4, 0.8, 0.9}} {
		x, y, z, t := p[0], p[1], p[2], p[3]
		dt := h.F(x, y, z, t+e)
		dm := h.F(x, y, z, t-e)
		// curl F by central differences
		d := func(c int, ax int) complex128 {
			q1, q2 := [3]float64{x, y, z}, [3]float64{x, y, z}
			q1[ax] += e
			q2[ax] -= e
			return (h.F(q1[0], q1[1], q1[2], t)[c] - h.F(q2[0], q2[1], q2[2], t)[c]) / complex(2*e, 0)
		}
		curl := [3]complex128{d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)}
		for c := 0; c < 3; c++ {
			ft := (dt[c] - dm[c]) / complex(2*e, 0)
			plus += cmplx.Abs(1i*ft - curl[c])
			minus += cmplx.Abs(-1i*ft - curl[c])
		}
	}
	if plus <= minus {
		return 1
	}
	return -1
}

// Knot sets up the run.
type Knot struct {
	L0      float64   // length scale [m]
	Cells   float64   // grid cells per L0
	Half    float64   // half-width of the region inside the absorber, in L0
	Times   []float64 // sample times, in L0/c
	Inside  float64   // radius of the "still here" sphere, in L0
	Absorb  int       // absorber thickness, cells
	Segment float64   // field-line tracing step, in L0
}

// DefaultKnot: L0 = 1 m, 10 cells per L0, region ±6 L0, samples to t = 4,
// sphere of 2 L0.
func DefaultKnot() Knot {
	return Knot{L0: 1, Cells: 10, Half: 6, Times: []float64{0, 0.5, 1, 1.5, 2, 3, 4}, Inside: 2, Absorb: 12, Segment: 0.01}
}

// Sample is one time.
type Sample struct {
	T                 float64 // in L0/c
	InsideFDTD        float64 // fraction of the starting energy within Inside·L0, solver
	InsideExact       float64 // the same from the exact field on the same points
	Absorbed          float64 // fraction that left through the walls
	LinkFDTD, LinkExt float64 // Gauss linking number of two B lines, solver and exact
	ClosedFDTD        bool    // both solver lines closed on themselves
	ClosedExact       bool    // both exact lines closed on themselves
	RMSRadius         float64 // exact √⟨r²⟩ of the energy, in L0
	HelicityFDTD      float64 // solver magnetic helicity ∫A·B dV / its value at t = 0
	HelicityExact     float64 // exact magnetic helicity on the same points / its value at t = 0
	CurlAErr          float64 // max |curl A − B| / max|B| inside the absorber (bookkeeping check)
}

// Result is the outcome.
type Result struct {
	Sign      float64
	Samples   []Sample
	W0        float64 // starting energy [J]
	Residual  float64 // ledger residual / W0
	DivB      float64 // max |div B| relative
	Seeds     [2][3]float64
	LinkStart float64 // exact linking at t = 0
	H0Ratio   float64 // solver helicity at t = 0 / exact helicity at t = 0
}

// Seeds of the two traced magnetic field lines, in L0.
var Seeds = [2][3]float64{{0.6, 0, 0}, {-0.6, 0, 0}}

// trace follows the field line of b from p0 with RK4 on the unit direction
// and returns the loop, closed reports whether it came back within 2·step
// of its start after leaving it.
func trace(b func(p [3]float64) [3]float64, p0 [3]float64, step float64, maxLen float64) (line [][3]float64, closed bool) {
	dir := func(p [3]float64) [3]float64 {
		v := b(p)
		n := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
		if n == 0 {
			return [3]float64{}
		}
		return [3]float64{v[0] / n, v[1] / n, v[2] / n}
	}
	add := func(p, d [3]float64, s float64) [3]float64 {
		return [3]float64{p[0] + s*d[0], p[1] + s*d[1], p[2] + s*d[2]}
	}
	p := p0
	line = append(line, p)
	left := false
	for s := 0.0; s < maxLen; s += step {
		k1 := dir(p)
		k2 := dir(add(p, k1, step/2))
		k3 := dir(add(p, k2, step/2))
		k4 := dir(add(p, k3, step))
		p = [3]float64{
			p[0] + step/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]),
			p[1] + step/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]),
			p[2] + step/6*(k1[2]+2*k2[2]+2*k3[2]+k4[2]),
		}
		line = append(line, p)
		d := math.Sqrt((p[0]-p0[0])*(p[0]-p0[0]) + (p[1]-p0[1])*(p[1]-p0[1]) + (p[2]-p0[2])*(p[2]-p0[2]))
		if d > 10*step {
			left = true
		}
		if left && d < 2*step {
			return line, true
		}
	}
	return line, false
}

// Linking is the Gauss linking number of two closed curves:
// (1/4π) ∮∮ (r1 − r2)·(dr1 × dr2) / |r1 − r2|³.
func Linking(a, b [][3]float64) float64 {
	var s float64
	for i := 0; i < len(a); i++ {
		a0, a1 := a[i], a[(i+1)%len(a)]
		ma := [3]float64{(a0[0] + a1[0]) / 2, (a0[1] + a1[1]) / 2, (a0[2] + a1[2]) / 2}
		da := [3]float64{a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]}
		for j := 0; j < len(b); j++ {
			b0, b1 := b[j], b[(j+1)%len(b)]
			mb := [3]float64{(b0[0] + b1[0]) / 2, (b0[1] + b1[1]) / 2, (b0[2] + b1[2]) / 2}
			db := [3]float64{b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]}
			r := [3]float64{ma[0] - mb[0], ma[1] - mb[1], ma[2] - mb[2]}
			cr := [3]float64{da[1]*db[2] - da[2]*db[1], da[2]*db[0] - da[0]*db[2], da[0]*db[1] - da[1]*db[0]}
			d := math.Sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])
			s += (r[0]*cr[0] + r[1]*cr[1] + r[2]*cr[2]) / (d * d * d)
		}
	}
	return s / (4 * math.Pi)
}

// potentials keeps the vector potential A on the E edges beside the
// solver, in the temporal gauge: the solver advances E first, then
// B ← B − Δt·curl E, so A ← A − Δt·E after each step keeps curl A = B
// (outside the absorber) at every step.
type potentials struct {
	f  *em.FDTD
	A  [3][]float64
	lo int // first interior index on each axis (inside the absorber)
	hi int // last interior index
	n  int
	dx float64
	sy int
	sz int
}

func newPotentials(f *em.FDTD, n, absorb int, dx float64, A func(x, y, z float64) [3]float64) *potentials {
	p := &potentials{f: f, lo: absorb + 1, hi: n - absorb - 2, n: n, dx: dx, sy: n, sz: n * n}
	for c := 0; c < 3; c++ {
		p.A[c] = make([]float64, len(f.E[c]))
	}
	for c := 0; c < 3; c++ {
		for k := 0; k < n; k++ {
			for j := 0; j < n; j++ {
				for i := 0; i < n; i++ {
					q := [3]float64{float64(i), float64(j), float64(k)}
					q[c] += 0.5
					x, y, z := f.Origin[0]+q[0]*dx, f.Origin[1]+q[1]*dx, f.Origin[2]+q[2]*dx
					p.A[c][i+n*(j+n*k)] = A(x, y, z)[c]
				}
			}
		}
	}
	return p
}

// step advances the solver one step and A with it.
func (p *potentials) step() {
	p.f.Step()
	dt := p.f.Dt
	for c := 0; c < 3; c++ {
		for i, e := range p.f.E[c] {
			p.A[c][i] -= dt * e
		}
	}
}

// curlAt is the discrete curl of A at B component c's face (i, j, k).
func (p *potentials) curlAt(c, idx int) float64 {
	c1, c2 := (c+1)%3, (c+2)%3
	st := [3]int{1, p.sy, p.sz}
	d := func(a, b int) float64 { return (p.A[b][idx+st[a]] - p.A[b][idx]) / p.dx }
	return d(c1, c2) - d(c2, c1)
}

// helicity is ∫A·B dV over the interior, with B averaged onto each A edge
// from its four neighbouring faces of the same component (the faces of
// component c around edge c lie at ±½ in the two other directions).
func (p *potentials) helicity() (h, curlErr float64) {
	B := p.f.B
	st := [3]int{1, p.sy, p.sz}
	var bmax, emax float64
	dV := p.dx * p.dx * p.dx
	for c := 0; c < 3; c++ {
		c1, c2 := (c+1)%3, (c+2)%3
		for k := p.lo; k <= p.hi; k++ {
			for j := p.lo; j <= p.hi; j++ {
				for i := p.lo; i <= p.hi; i++ {
					idx := i + p.n*(j+p.n*k)
					// A_c at (i+½δc): B_c faces at (i, j±½, k±½) for c = x → indices
					// idx, idx−st[c1], idx−st[c2], idx−st[c1]−st[c2], shifted +½ along c
					// by averaging the two along c as well.
					var bs float64
					for _, o := range []int{0, st[c]} {
						bs += B[c][idx+o] + B[c][idx+o-st[c1]] + B[c][idx+o-st[c2]] + B[c][idx+o-st[c1]-st[c2]]
					}
					h += p.A[c][idx] * bs / 8 * dV
					b := B[c][idx]
					bmax = math.Max(bmax, math.Abs(b))
					emax = math.Max(emax, math.Abs(p.curlAt(c, idx)-b))
				}
			}
		}
	}
	return h, emax / bmax
}

// Run builds the region, sets the hopfion, follows it and measures.
func (k Knot) Run() (Result, error) {
	h := Hopfion{}
	c := em.SpeedOfLight
	L0 := k.L0
	dx := L0 / k.Cells
	n := 2*int(math.Ceil(k.Half*k.Cells)) + 2*k.Absorb + 1
	half := float64(n-1) / 2 * dx
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx},
		Origin: [3]float64{-half, -half, -half}, Absorber: k.Absorb, Ledger: true})
	if err != nil {
		return Result{}, err
	}
	res := Result{Sign: h.Sign(), Seeds: Seeds}
	const E0 = 1.0 // V/m scale of the field
	// E = curl C with C = E0·L0·Re(α∇̃β) [V]; B = curl A with
	// A = sign·(E0·L0/c)·Im(α∇̃β) [T·m], B taken at t = Dt/2 (leapfrog).
	f.SetEFromPotential(func(x, y, z float64) [3]float64 {
		p := h.P(x/L0, y/L0, z/L0, 0)
		return [3]float64{E0 * L0 * real(p[0]), E0 * L0 * real(p[1]), E0 * L0 * real(p[2])}
	})
	tb := c * f.Dt / 2 / L0
	aFn := func(x, y, z float64) [3]float64 {
		p := h.P(x/L0, y/L0, z/L0, tb)
		s := res.Sign * E0 * L0 / c
		return [3]float64{s * imag(p[0]), s * imag(p[1]), s * imag(p[2])}
	}
	f.SetBFromPotential(aFn)
	pot := newPotentials(f, n, k.Absorb, dx, aFn)
	all := func(x, y, z float64) bool { return true }
	ins := func(x, y, z float64) bool { return x*x+y*y+z*z < k.Inside*k.Inside*L0*L0 }
	res.W0 = f.RegionEnergy(all)
	// Exact energy on the interior grid points, |F|² ∝ energy density.
	exact := func(t float64) (inside, total, r2, hel float64) {
		lim := int(math.Ceil(k.Half * k.Cells))
		for i := -lim; i <= lim; i++ {
			for j := -lim; j <= lim; j++ {
				for l := -lim; l <= lim; l++ {
					x, y, z := float64(i)/k.Cells, float64(j)/k.Cells, float64(l)/k.Cells
					F := h.F(x, y, z, t)
					P := h.P(x, y, z, t)
					hel += imag(P[0])*imag(F[0]) + imag(P[1])*imag(F[1]) + imag(P[2])*imag(F[2])
					u := real(F[0]*cmplx.Conj(F[0]) + F[1]*cmplx.Conj(F[1]) + F[2]*cmplx.Conj(F[2]))
					total += u
					rr := x*x + y*y + z*z
					r2 += u * rr
					if rr < k.Inside*k.Inside {
						inside += u
					}
				}
			}
		}
		return
	}
	_, total0, _, hel0 := exact(0)
	// exact helicity in SI on the grid: A·B = (E0·L0/c)(E0/c)·Im P·Im F, dV = dx³.
	helScale := E0 * L0 / c * E0 / c * dx * dx * dx
	var hFDTD0 float64
	maxLen := 40.0
	bExact := func(t float64) func(p [3]float64) [3]float64 {
		return func(p [3]float64) [3]float64 {
			F := h.F(p[0], p[1], p[2], t)
			return [3]float64{res.Sign * imag(F[0]), res.Sign * imag(F[1]), res.Sign * imag(F[2])}
		}
	}
	bFDTD := func(p [3]float64) [3]float64 {
		_, B := f.FieldAt(p[0]*L0, p[1]*L0, p[2]*L0)
		return B
	}
	for _, t := range k.Times {
		for f.Time() < t*L0/c-1e-9*f.Dt {
			pot.step()
		}
		s := Sample{T: t}
		hf, ce := pot.helicity()
		if t == 0 {
			hFDTD0 = hf
			res.H0Ratio = hf / (hel0 * helScale)
		}
		s.HelicityFDTD = hf / hFDTD0
		s.CurlAErr = ce
		s.InsideFDTD = f.RegionEnergy(ins) / res.W0
		s.Absorbed = f.Ledger().Absorbed / res.W0
		in, tot, r2, hl := exact(t)
		s.HelicityExact = hl / hel0
		s.InsideExact = in / total0
		s.RMSRadius = math.Sqrt(r2 / tot)
		le1, e1 := trace(bExact(t), Seeds[0], k.Segment, maxLen)
		le2, e2 := trace(bExact(t), Seeds[1], k.Segment, maxLen)
		s.ClosedExact = e1 && e2
		s.LinkExt = Linking(le1, le2)
		lf1, c1 := trace(bFDTD, Seeds[0], k.Segment, maxLen)
		lf2, c2 := trace(bFDTD, Seeds[1], k.Segment, maxLen)
		s.ClosedFDTD = c1 && c2
		s.LinkFDTD = Linking(lf1, lf2)
		if t == 0 {
			res.LinkStart = s.LinkExt
		}
		res.Samples = append(res.Samples, s)
	}
	res.Residual = f.Ledger().Residual() / res.W0
	wd, sc := f.DivB()
	res.DivB = wd / sc
	return res, nil
}

func main() {
	fmt.Println("knotted_light (H29, step 1): Rañada's knotted light (hopfion) in the field solver")
	fmt.Println("Question: does topology alone, under linear Maxwell equations, hold a knot of light together?")
	k := DefaultKnot()
	r, err := k.Run()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("\nGrid: %g cells per L0, region ±%g L0, absorber %d cells. Helicity sign of the formula: %+g\n", k.Cells, k.Half, k.Absorb, r.Sign)
	fmt.Printf("Field-line seeds (L0): %v and %v; linking number at t = 0 (exact): %.4f\n", r.Seeds[0], r.Seeds[1], r.LinkStart)
	fmt.Printf("Solver helicity at t = 0 / exact: %.4f\n", r.H0Ratio)
	fmt.Println("\n  t (L0/c) | energy within 2 L0: solver / exact | absorbed | rms radius | helicity: solver / exact | curl A − B | linking of two B lines: exact (closed) / solver (closed)")
	for _, s := range r.Samples {
		fmt.Printf("  %8.2f | %17.4f / %.4f | %8.4f | %10.3f | %9.4f / %.4f | %10.1e | %+.4f (%v) / %+.4f (%v)\n",
			s.T, s.InsideFDTD, s.InsideExact, s.Absorbed, s.RMSRadius, s.HelicityFDTD, s.HelicityExact, s.CurlAErr,
			s.LinkExt, s.ClosedExact, s.LinkFDTD, s.ClosedFDTD)
	}
	fmt.Printf("\nLedger residual %.1e of the starting energy; max |div B| %.1e (relative)\n", r.Residual, r.DivB)
	fmt.Println("\nStandard prediction: linking conserved (topology carried by Maxwell's equations), energy not held:")
	fmt.Println("the knot expands at c while staying linked. Holding it needs a nonlinear term (Skyrme; Faddeev–Niemi;")
	fmt.Println("Derrick's theorem: no static lump of a linear field in three dimensions).")
}
