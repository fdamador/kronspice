// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// H30 (CLAUDE.md §1): implosion and explosion. Rañada's knot is an exact
// solution for all times, not only t ≥ 0: before t = 0 the same knot comes
// in from far away and contracts, at t = 0 it is tightest, after it expands
// (energy within 2 L0: 0.02 at t = −3, 0.96 at 0, 0.02 at +3). Step 1 of
// H29 watched only the expanding half. Here the incoming half is fed in
// through a closed Huygens surface around the region (total-field /
// scattered-field: inside it the solver carries the whole field, outside
// only what departs from the declared incoming wave; Taflove & Hagness,
// Computational Electrodynamics, 3rd ed., ch. 5). The surface is a row of
// equivalent edge currents and magnetic face currents (em.EdgeSource,
// em.FaceSource), so the solver's ledger books every joule the surface
// puts in and takes out.
//
// The declared incoming wave is the stand-in for what an ether, an
// absorber or any other medium would have to supply. Its power is read off
// the ledger; where it would come from is the open question, and the run
// does not answer it.
//
// Part 1: one implosion then explosion, t = −4 … 4 L0/c.
// Part 2: a pulsing knot, a train of knots one period apart fed in from
// the boundary: per cycle the surface puts in what it takes out, and the
// knot is rebuilt at the centre every cycle.

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// Cycle sets up a run with a Huygens surface.
type Cycle struct {
	L0     float64   // length scale [m]
	Cells  float64   // cells per L0
	Half   float64   // half-width of the region inside the absorber, L0
	Absorb int       // absorber cells
	Inside float64   // radius of the energy sphere, L0
	Start  float64   // first time, L0/c
	Times  []float64 // sample times, L0/c (ascending, ≥ Start)
	Pulses []float64 // times at which a knot is tightest, L0/c
}

// DefaultImplosion is part 1: one knot, tightest at t = 0, followed from
// t = −10 (when 0.4 % of it is inside the surface) to t = 10.
func DefaultImplosion() Cycle {
	return Cycle{L0: 1, Cells: 10, Half: 6, Absorb: 12, Inside: 2, Start: -10,
		Times: []float64{-10, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 10}, Pulses: []float64{0}}
}

// DefaultPulsing is part 2: knots tightest at t = 0, 8, 16, followed from
// t = −4 to t = 20 (three cycles of period 8 L0/c).
func DefaultPulsing() Cycle {
	var ts []float64
	for t := -4.0; t <= 20; t += 1 {
		ts = append(ts, t)
	}
	return Cycle{L0: 1, Cells: 10, Half: 6, Absorb: 12, Inside: 2, Start: -4, Times: ts, Pulses: []float64{0, 8, 16}}
}

// CycleSample is one time.
type CycleSample struct {
	T                       float64 // L0/c
	InsideFDTD, InsideExact float64 // energy within Inside·L0 / one knot's energy
	Box                     float64 // energy inside the surface / one knot's energy
	In, Out                 float64 // energy the surface put in (steps with net inflow) and took out (net outflow), cumulative / one knot's energy
	Leak                    float64 // energy outside the surface plus absorbed / one knot's energy (numerical error)
	FieldErr                float64 // ‖E_solver − E_exact‖/‖E_exact‖ within Inside·L0
	Link                    float64 // Gauss linking of two solver B lines, each closed at its nearest return
	Gap                     float64 // the larger nearest return of the two, L0 (0 for an exact closure)
	Traced                  bool    // lines traced at this time (the times a knot is tightest)
}

// CycleResult is the outcome.
type CycleResult struct {
	Samples  []CycleSample
	Knot     float64 // one knot's energy [J]
	Residual float64 // ledger residual / one knot's energy
	NE, NB   int     // number of edge and face sources on the surface
}

// knotEnergy is one knot's energy in units of ε0E0²L0³/2: ∫|F|² d³u at
// t = 0, where the energy density is spherically symmetric (checked in
// the tests), as a radial integral.
func knotEnergy() float64 {
	h := Hopfion{}
	var s float64
	const dr, rmax = 1e-3, 200.0
	for r := dr / 2; r < rmax; r += dr {
		F := h.F(r, 0, 0, 0)
		u := real(F[0]*cmplx.Conj(F[0]) + F[1]*cmplx.Conj(F[1]) + F[2]*cmplx.Conj(F[2]))
		s += 4 * math.Pi * r * r * u * dr
	}
	return s
}

// Run builds the region, sets the field at Start inside the surface, puts
// the incoming train on the surface and follows it.
func (cy Cycle) Run() (CycleResult, error) {
	h := Hopfion{}
	sign := h.Sign()
	c := em.SpeedOfLight
	L0 := cy.L0
	dx := L0 / cy.Cells
	n := 2*int(math.Ceil(cy.Half*cy.Cells)) + 2*cy.Absorb + 1
	half := float64(n-1) / 2 * dx
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx},
		Origin: [3]float64{-half, -half, -half}, Absorber: cy.Absorb, Ledger: true})
	if err != nil {
		return CycleResult{}, err
	}
	f.SetStartTime(cy.Start * L0 / c)
	const E0 = 1.0
	// The declared incoming wave: the train of knots, in SI.
	F := func(x, y, z, t float64) [3]complex128 {
		var s [3]complex128
		for _, tk := range cy.Pulses {
			q := h.F(x/L0, y/L0, z/L0, c*t/L0-tk)
			s[0] += q[0]
			s[1] += q[1]
			s[2] += q[2]
		}
		return s
	}
	P := func(x, y, z, t float64) [3]complex128 {
		var s [3]complex128
		for _, tk := range cy.Pulses {
			q := h.P(x/L0, y/L0, z/L0, c*t/L0-tk)
			s[0] += q[0]
			s[1] += q[1]
			s[2] += q[2]
		}
		return s
	}
	eInc := func(comp int, x, y, z, t float64) float64 { return E0 * real(F(x, y, z, t)[comp]) }
	bInc := func(comp int, x, y, z, t float64) float64 { return sign * E0 / c * imag(F(x, y, z, t)[comp]) }

	// The surface: a cube of half-width b, between grid planes so that no
	// field component lies on it.
	b := (math.Floor(cy.Half*cy.Cells) - 2.25) * dx
	pos := func(p [3]float64) [3]float64 {
		return [3]float64{-half + p[0]*dx, -half + p[1]*dx, -half + p[2]*dx}
	}
	ePos := func(cc, i, j, k int) [3]float64 {
		p := [3]float64{float64(i), float64(j), float64(k)}
		p[cc] += 0.5
		return pos(p)
	}
	bPos := func(cc, i, j, k int) [3]float64 {
		p := [3]float64{float64(i) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
		p[cc] -= 0.5
		return pos(p)
	}
	in := func(p [3]float64) float64 {
		if math.Abs(p[0]) < b && math.Abs(p[1]) < b && math.Abs(p[2]) < b {
			return 1
		}
		return 0
	}
	type term struct {
		comp  int
		p     [3]float64
		kappa float64 // curl coefficient × (χ_self − χ_neighbour)
	}
	// Curl stencils. For E edge c at (i,j,k): curl_h B = Σ κ B over faces
	// (Yee, as em.FDTD's ampere); for B face c: curl_e E = Σ κ E over edges.
	type nb struct {
		comp       int
		di, dj, dk int
		kappa      float64
	}
	eSt := [3][]nb{
		{{2, 0, 0, 0, 1 / dx}, {2, 0, -1, 0, -1 / dx}, {1, 0, 0, 0, -1 / dx}, {1, 0, 0, -1, 1 / dx}},
		{{0, 0, 0, 0, 1 / dx}, {0, 0, 0, -1, -1 / dx}, {2, 0, 0, 0, -1 / dx}, {2, -1, 0, 0, 1 / dx}},
		{{1, 0, 0, 0, 1 / dx}, {1, -1, 0, 0, -1 / dx}, {0, 0, 0, 0, -1 / dx}, {0, 0, -1, 0, 1 / dx}},
	}
	bSt := [3][]nb{
		{{2, 0, 1, 0, 1 / dx}, {2, 0, 0, 0, -1 / dx}, {1, 0, 0, 1, -1 / dx}, {1, 0, 0, 0, 1 / dx}},
		{{0, 0, 0, 1, 1 / dx}, {0, 0, 0, 0, -1 / dx}, {2, 1, 0, 0, -1 / dx}, {2, 0, 0, 0, 1 / dx}},
		{{1, 1, 0, 0, 1 / dx}, {1, 0, 0, 0, -1 / dx}, {0, 0, 1, 0, -1 / dx}, {0, 0, 0, 0, 1 / dx}},
	}
	lo, hi := cy.Absorb+1, n-cy.Absorb-2
	res := CycleResult{}
	area := dx * dx
	for cc := 0; cc < 3; cc++ {
		for k := lo; k <= hi; k++ {
			for j := lo; j <= hi; j++ {
				for i := lo; i <= hi; i++ {
					// E edge
					self := in(ePos(cc, i, j, k))
					var ts []term
					for _, s := range eSt[cc] {
						d := self - in(bPos(s.comp, i+s.di, j+s.dj, k+s.dk))
						if d != 0 {
							ts = append(ts, term{s.comp, bPos(s.comp, i+s.di, j+s.dj, k+s.dk), s.kappa * d})
						}
					}
					if len(ts) > 0 {
						ts := ts
						// E ← … + cb·(curl_h B_used)/μ0 − cb·J: J = −Σκ B_inc/μ0.
						cur := func(t float64) float64 {
							var corr float64
							for _, q := range ts {
								corr += q.kappa * bInc(q.comp, q.p[0], q.p[1], q.p[2], t)
							}
							return -corr / em.Mu0 * area
						}
						if err := f.AddSource(em.EdgeSource{Comp: cc, I: i, J: j, K: k, Current: cur}); err != nil {
							return res, err
						}
						res.NE++
					}
					// B face
					self = in(bPos(cc, i, j, k))
					ts = nil
					for _, s := range bSt[cc] {
						d := self - in(ePos(s.comp, i+s.di, j+s.dj, k+s.dk))
						if d != 0 {
							ts = append(ts, term{s.comp, ePos(s.comp, i+s.di, j+s.dj, k+s.dk), s.kappa * d})
						}
					}
					if len(ts) > 0 {
						ts := ts
						// B ← B − dt·(curl_e E_used): M = Σκ E_inc.
						m := func(t float64) float64 {
							var corr float64
							for _, q := range ts {
								corr += q.kappa * eInc(q.comp, q.p[0], q.p[1], q.p[2], t)
							}
							return corr
						}
						if err := f.AddMagneticSource(em.FaceSource{Comp: cc, I: i, J: j, K: k, M: m}); err != nil {
							return res, err
						}
						res.NB++
					}
				}
			}
		}
	}
	// Starting field: the train at Start inside the surface, zero outside
	// (the departure from the incoming wave starts at zero). Set as curls
	// of the potentials (div-free inside), then cut at the surface.
	t0 := cy.Start * L0 / c
	f.SetEFromPotential(func(x, y, z float64) [3]float64 {
		p := P(x, y, z, t0)
		return [3]float64{E0 * L0 * real(p[0]), E0 * L0 * real(p[1]), E0 * L0 * real(p[2])}
	})
	f.SetBFromPotential(func(x, y, z float64) [3]float64 {
		p := P(x, y, z, t0+f.Dt/2)
		s := sign * E0 * L0 / c
		return [3]float64{s * imag(p[0]), s * imag(p[1]), s * imag(p[2])}
	})
	for cc := 0; cc < 3; cc++ {
		for k := 0; k < n; k++ {
			for j := 0; j < n; j++ {
				for i := 0; i < n; i++ {
					idx := i + n*(j+n*k)
					if in(ePos(cc, i, j, k)) == 0 {
						f.E[cc][idx] = 0
					}
					if in(bPos(cc, i, j, k)) == 0 {
						f.B[cc][idx] = 0
					}
				}
			}
		}
	}
	res.Knot = knotEnergy() * 0.5 * em.Eps0 * E0 * E0 * L0 * L0 * L0
	U := res.Knot
	insFn := func(x, y, z float64) bool { return x*x+y*y+z*z < cy.Inside*cy.Inside*L0*L0 }
	boxFn := func(x, y, z float64) bool { return math.Abs(x) < b && math.Abs(y) < b && math.Abs(z) < b }
	outFn := func(x, y, z float64) bool { return !boxFn(x, y, z) }
	var cumIn, cumOut, prevWork float64
	lim := int(math.Ceil(cy.Inside * cy.Cells))
	exactInside := func(t float64) float64 {
		var s float64
		for i := -lim; i <= lim; i++ {
			for j := -lim; j <= lim; j++ {
				for l := -lim; l <= lim; l++ {
					x, y, z := float64(i)*dx, float64(j)*dx, float64(l)*dx
					if x*x+y*y+z*z >= cy.Inside*cy.Inside*L0*L0 {
						continue
					}
					q := F(x, y, z, t)
					s += real(q[0]*cmplx.Conj(q[0]) + q[1]*cmplx.Conj(q[1]) + q[2]*cmplx.Conj(q[2]))
				}
			}
		}
		return s * 0.5 * em.Eps0 * E0 * E0 * dx * dx * dx
	}
	fieldErr := func(t float64) float64 {
		var num, den float64
		for i := -lim; i <= lim; i += 2 {
			for j := -lim; j <= lim; j += 2 {
				for l := -lim; l <= lim; l += 2 {
					x, y, z := float64(i)*dx, float64(j)*dx, float64(l)*dx
					if x*x+y*y+z*z >= cy.Inside*cy.Inside*L0*L0 {
						continue
					}
					E, _ := f.FieldAt(x, y, z)
					q := F(x, y, z, t)
					for a := 0; a < 3; a++ {
						d := E[a] - E0*real(q[a])
						num += d * d
						den += E0 * E0 * real(q[a]) * real(q[a])
					}
				}
			}
		}
		return math.Sqrt(num / den)
	}
	bFDTD := func(p [3]float64) [3]float64 {
		_, B := f.FieldAt(p[0]*L0, p[1]*L0, p[2]*L0)
		return B
	}
	step := func() {
		f.Step()
		w := f.Ledger().SourceWork
		if d := w - prevWork; d > 0 {
			cumIn += d
		} else {
			cumOut -= d
		}
		prevWork = w
	}
	for _, ts := range cy.Times {
		for f.Time() < ts*L0/c-1e-9*f.Dt {
			step()
		}
		t := f.Time()
		s := CycleSample{T: ts}
		s.InsideFDTD = f.RegionEnergy(insFn) / U
		s.InsideExact = exactInside(t) / U
		s.Box = f.RegionEnergy(boxFn) / U
		s.In, s.Out = cumIn/U, cumOut/U
		s.Leak = (f.RegionEnergy(outFn) + f.Ledger().Absorbed) / U
		s.FieldErr = fieldErr(t)
		for _, tk := range cy.Pulses {
			if ts == tk {
				l1, g1 := traceGap(bFDTD, Seeds[0], 0.01, 40)
				l2, g2 := traceGap(bFDTD, Seeds[1], 0.01, 40)
				s.Traced, s.Gap, s.Link = true, math.Max(g1, g2), Linking(l1, l2)
			}
		}
		res.Samples = append(res.Samples, s)
	}
	res.Residual = f.Ledger().Residual() / U
	return res, nil
}
