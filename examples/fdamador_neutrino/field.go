// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 5b: implosion and explosion of the E-ring in the field solver
// (em/fdtd.go), the field twin of Part 5a.
//
// The ring is a dielectric torus (εr = 100, major radius a = 10 mm, tube
// radius r = 2.5 mm, lossless) on a Yee grid of r/Cells cells, staircased.
// Two runs per grid:
//
//  1. Ring-down. The ring starts with an azimuthal displacement field in
//     the tube, D_φ = D0·cos²(πs/2r) (s the distance from the tube's centre
//     line), B = 0: E round the ring, no charge (D = curl of a potential,
//     so the discrete Gauss law holds exactly). It rings and radiates
//     through an absorbing layer. From it: f0 (zero crossings of E_φ in
//     the tube), Q (decay of the energy inside the Huygens surface) and the
//     ring's outgoing field on the surface, as a phasor at f0 (a discrete
//     Fourier transform of samples taken 16 times a period, after the
//     start-up transient has left).
//
//  2. Implosion then explosion. A closed Huygens (total-field /
//     scattered-field) surface round the ring, as in examples/knotted_light
//     (cycle.go), injects the declared incoming wave: the time reverse of
//     that outgoing mode field, Ê → conj(Ê), B̂ → −conj(B̂), growing as
//     e^{+γt} (γ = ω0/2Q). It is switched off at the peak, chosen at an E
//     maximum in the ring (interior B zero); the ring then explodes and its
//     outgoing wave crosses the surface into the absorber. The surface's
//     edge and face sources book every joule they put in; energy appearing
//     outside the surface or in the absorber before the peak is what the
//     ring reflected (or the declared wave failed to match), after it is
//     the explosion. Switching off loses the declared field outside the
//     surface at that instant (the mode's weak outer E), counted with it.
//
//  3. A constant burst of the same incoming pattern and the same incident
//     energy (raised-cosine edges one period long, 2 τ_E), not continued
//     as an outgoing wave: what the ring does not take is reflected through
//     the surface into the absorber.
//
// The incident energy of a declared wave with envelope g(t) (relative to
// the mode pattern Ê at t = 0) is 2γ·W0·∫g² dt: the flux of the time-
// reversed pattern is the power the mode radiates, 2γ·W0 at unit
// amplitude, with W0 the ring-down's energy extrapolated to t = 0.
//
// Approximations stated: staircased torus, lumped-free (no tan δ), the
// absorber a graded lossy layer (reflection ~1e-6 at normal incidence) at
// 3a from the centre, where the ring's near field is small but not zero;
// convergence is checked on two grids.

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// FieldRing is the field-solver ring.
type FieldRing struct {
	A, R, EpsR float64
	Cells      float64 // grid cells per tube radius
	HalfA      float64 // half-width of the region inside the absorber, in a
	SurfA      float64 // half-width of the Huygens surface, in a
	Absorb     int     // absorbing-layer cells
	FilterT    float64 // length of the filtering run, in lumped periods
	RingDownTE float64 // ring-down length, in energy decay times τ_E (lumped estimate)
	ImplodeTE  float64 // incoming growth window, in τ_E (measured)
	ExplodeTE  float64 // outgoing window after the peak, in τ_E
	// Log, when set, receives progress lines.
	Log func(format string, args ...any)
}

// DefaultFieldRing is the bench ring of Part 1 at r/4 cells.
func DefaultFieldRing() FieldRing {
	return FieldRing{A: 10e-3, R: 2.5e-3, EpsR: 100, Cells: 4, HalfA: 3, SurfA: 2, Absorb: 12,
		FilterT: 15, RingDownTE: 4, ImplodeTE: 6, ExplodeTE: 3}
}

// FieldSample is one time of a run.
type FieldSample struct {
	T                float64 // time, periods (from the run's start)
	Inside           float64 // energy inside the surface [J]
	Outside          float64 // energy outside the surface, not yet absorbed [J]
	Absorbed         float64 // absorbed by the layer [J]
	SurfIn, SurfOut  float64 // cumulative work the surface put in / took out [J]
	Probe            float64 // E_φ at (a, 0, 0) [V/m]
	ExactInc         float64 // declared incident energy delivered so far [J]
	Residual         float64 // ledger residual [J]
	DeclaredIncoming bool
}

// FieldResult is everything Part 5b reports for one grid.
type FieldResult struct {
	Cells, Dx           float64
	Nodes               int
	F0, Q               float64 // ring-down: frequency [Hz], Q from the energy decay
	QProbe              float64 // Q from the decay of the probe's peaks
	W0                  float64 // ring-down energy inside the surface extrapolated to t = 0 [J]
	RingDown            []FieldSample
	RDResidual          float64 // ring-down ledger residual / initial energy
	EInc                float64 // incident energy of the declared incoming half [J]
	Peak                float64 // energy inside the surface at the peak [J]
	PeakT               float64 // periods from the drive's start
	LeakAtPeak          float64 // outside + absorbed at the peak [J] (reflection + mismatch)
	Delivered           float64 // net work of the surface up to the peak [J]
	LeakEnd             float64 // outside + absorbed at the end [J]
	SurfInTotal         float64 // surface work in, whole run [J]
	SurfOutTotal        float64 // surface work out, whole run [J]
	InsideEnd           float64
	Residual            float64 // implosion-run ledger residual / EInc
	Implosion           []FieldSample
	EPhiFrac, BPhiFrac  float64 // azimuthal share of electric energy in the tube at the peak; of magnetic energy a quarter period later
	BzCentre, BzOutside float64 // B_z at the centre and at (2a, 0, 0) a quarter period after the peak [T]
	Corr                float64 // |correlation| of E at the peak with the ring-down's E at an E maximum
	Burst               []FieldSample
	BurstEInc           float64
	BurstPeak           float64
	BurstLeakDrive      float64 // outside + absorbed when the burst ends [J]
	BurstResidual       float64
	Plane               []PlanePoint // y = 0 plane at the peak and a quarter period later
	Seconds             float64
}

// PlanePoint is the field on the y = 0 plane (x ≥ 0).
type PlanePoint struct {
	X, Z               float64 // m
	Ex, Ey, Ez         float64 // at the peak
	Bx, By, Bz         float64 // a quarter period later
	InRing             bool
	ExQ, EyQ, EzQ      float64 // E a quarter period later
	BxPeak, ByP, BzPea float64 // B at the peak
}

func (fr FieldRing) inTorus(x, y, z float64) bool {
	s := math.Hypot(math.Hypot(x, y)-fr.A, z)
	return s < fr.R
}

// grid is a built region with its Huygens surface terms.
type grid struct {
	fr      FieldRing
	n       int
	dx      float64
	half    float64
	b       float64 // surface half-width [m]
	pts     []spt   // field samples the surface needs
	eTerms  [][]sterm
	bTerms  [][]sterm
	eEdges  [][4]int // comp, i, j, k of each E source edge
	bFaces  [][4]int
	matFunc func(x, y, z float64) (float64, float64)
}

type spt struct {
	isB     bool
	c, i, j int
	k       int
}

type sterm struct {
	pt    int
	kappa float64
}

func (fr FieldRing) newGrid() *grid {
	dx := fr.R / fr.Cells
	m := int(math.Ceil(fr.HalfA * fr.A / dx))
	n := 2*m + 2*fr.Absorb + 1
	g := &grid{fr: fr, n: n, dx: dx, half: float64(n-1) / 2 * dx}
	g.b = (math.Round(fr.SurfA*fr.A/dx) + 0.25) * dx
	g.matFunc = func(x, y, z float64) (float64, float64) {
		if fr.inTorus(x, y, z) {
			return fr.EpsR, 0
		}
		return 1, 0
	}
	// Surface stencils (as examples/knotted_light/cycle.go): an E edge
	// inside whose curl uses a B face outside (or the reverse) gets an
	// edge current from the declared B there; a B face whose curl crosses
	// the surface gets a magnetic current from the declared E.
	pos := func(p [3]float64) [3]float64 {
		return [3]float64{-g.half + p[0]*dx, -g.half + p[1]*dx, -g.half + p[2]*dx}
	}
	ePos := func(c, i, j, k int) [3]float64 {
		p := [3]float64{float64(i), float64(j), float64(k)}
		p[c] += 0.5
		return pos(p)
	}
	bPos := func(c, i, j, k int) [3]float64 {
		p := [3]float64{float64(i) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
		p[c] -= 0.5
		return pos(p)
	}
	in := func(p [3]float64) float64 {
		if math.Abs(p[0]) < g.b && math.Abs(p[1]) < g.b && math.Abs(p[2]) < g.b {
			return 1
		}
		return 0
	}
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
	index := map[spt]int{}
	ptOf := func(p spt) int {
		if q, ok := index[p]; ok {
			return q
		}
		index[p] = len(g.pts)
		g.pts = append(g.pts, p)
		return index[p]
	}
	lo, hi := fr.Absorb+1, n-fr.Absorb-2
	for c := 0; c < 3; c++ {
		for k := lo; k <= hi; k++ {
			for j := lo; j <= hi; j++ {
				for i := lo; i <= hi; i++ {
					self := in(ePos(c, i, j, k))
					var ts []sterm
					for _, s := range eSt[c] {
						if d := self - in(bPos(s.comp, i+s.di, j+s.dj, k+s.dk)); d != 0 {
							ts = append(ts, sterm{ptOf(spt{true, s.comp, i + s.di, j + s.dj, k + s.dk}), s.kappa * d})
						}
					}
					if len(ts) > 0 {
						g.eTerms = append(g.eTerms, ts)
						g.eEdges = append(g.eEdges, [4]int{c, i, j, k})
					}
					self = in(bPos(c, i, j, k))
					ts = nil
					for _, s := range bSt[c] {
						if d := self - in(ePos(s.comp, i+s.di, j+s.dj, k+s.dk)); d != 0 {
							ts = append(ts, sterm{ptOf(spt{false, s.comp, i + s.di, j + s.dj, k + s.dk}), s.kappa * d})
						}
					}
					if len(ts) > 0 {
						g.bTerms = append(g.bTerms, ts)
						g.bFaces = append(g.bFaces, [4]int{c, i, j, k})
					}
				}
			}
		}
	}
	return g
}

func (g *grid) fdtd() (*em.FDTD, error) {
	return em.NewFDTD(em.FDTDConfig{N: [3]int{g.n, g.n, g.n}, D: [3]float64{g.dx, g.dx, g.dx},
		Origin: [3]float64{-g.half, -g.half, -g.half}, Absorber: g.fr.Absorb, Ledger: true, Material: g.matFunc})
}

func (g *grid) idx(i, j, k int) int { return i + g.n*(j+g.n*k) }

// regionFns are the inside-surface and outside-surface tests.
func (g *grid) regionFns() (in, out func(x, y, z float64) bool) {
	b := g.b
	in = func(x, y, z float64) bool { return math.Abs(x) < b && math.Abs(y) < b && math.Abs(z) < b }
	out = func(x, y, z float64) bool { return !in(x, y, z) }
	return
}

// Run does the ring-down, the implosion–explosion and the burst on one grid.
func (fr FieldRing) Run() (FieldResult, error) {
	g := fr.newGrid()
	res := FieldResult{Cells: fr.Cells, Dx: g.dx, Nodes: g.n * g.n * g.n}
	// --- 1. ring-down ---
	f, err := g.fdtd()
	if err != nil {
		return res, err
	}
	// D_φ = −∂C_z/∂ρ with C_z(ρ, z) = −∫₀^ρ D_φ dρ'; E = D/ε per edge.
	bump := func(rho, z float64) float64 {
		s := math.Hypot(rho-fr.A, z)
		if s >= fr.R {
			return 0
		}
		c := math.Cos(math.Pi * s / (2 * fr.R))
		return c * c
	}
	f.SetEFromPotential(func(x, y, z float64) [3]float64 {
		rho := math.Hypot(x, y)
		lo := math.Max(0, fr.A-fr.R)
		if rho <= lo || math.Abs(z) >= fr.R {
			return [3]float64{}
		}
		hi := math.Min(rho, fr.A+fr.R)
		const m = 64
		h := (hi - lo) / m
		var s float64
		for q := 0; q < m; q++ {
			s += bump(lo+(float64(q)+0.5)*h, z) * h
		}
		return [3]float64{0, 0, -s}
	})
	for c := 0; c < 3; c++ {
		for k := 0; k < g.n; k++ {
			for j := 0; j < g.n; j++ {
				for i := 0; i < g.n; i++ {
					p := [3]float64{float64(i), float64(j), float64(k)}
					p[c] += 0.5
					er, _ := g.matFunc(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx)
					f.E[c][g.idx(i, j, k)] /= er
				}
			}
		}
	}
	inFn, outFn := g.regionFns()
	e := DefaultERing()
	e.A, e.R, e.EpsR = fr.A, fr.R, fr.EpsR
	Tl := 2 * math.Pi / e.Omega0()
	tauL := e.QRad() / e.Omega0()
	stopRD := fr.RingDownTE * tauL
	every := int(math.Max(1, math.Round(Tl/16/f.Dt)))
	nPts := len(g.pts)
	var sampT []float64
	var sampE, sampB [][]float32 // per sample: values at every surface point (E at E time, B at B time)
	var probeT, probe []float64
	var times []float64
	var winside []float64
	// The bump excites the TE01δ-like mode near 2.05 GHz and, more weakly,
	// higher n = 0 modes of the tube's cross-section (near 7 GHz, higher
	// Q). A first ring-down of FilterT lumped periods filters the field:
	// the discrete Fourier transform of every E and B component at the
	// lumped f0 over [2T, FilterT·T] is the low mode's field (the others are
	// 5 GHz away, rejected by the window). The ring-down proper starts
	// from that filtered field.
	{
		wl := e.Omega0()
		var pE, pB [3][]complex128
		for c := 0; c < 3; c++ {
			pE[c] = make([]complex128, len(f.E[c]))
			pB[c] = make([]complex128, len(f.B[c]))
		}
		for f.Time() < fr.FilterT*Tl {
			f.Step()
			if f.Steps%every != 0 || f.Time() < 2*Tl {
				continue
			}
			ce := cmplx.Exp(complex(0, -wl*f.Time()))
			cb := cmplx.Exp(complex(0, -wl*(f.Time()+f.Dt/2)))
			for c := 0; c < 3; c++ {
				ec, bc, pe, pb := f.E[c], f.B[c], pE[c], pB[c]
				for n := range ec {
					pe[n] += complex(ec[n], 0) * ce
					pb[n] += complex(bc[n], 0) * cb
				}
			}
		}
		f, err = g.fdtd()
		if err != nil {
			return res, err
		}
		cb := cmplx.Exp(complex(0, wl*f.Dt/2))
		for c := 0; c < 3; c++ {
			for n := range f.E[c] {
				f.E[c][n] = real(pE[c][n])
				f.B[c][n] = real(pB[c][n] * cb)
			}
		}
		if fr.Log != nil {
			fr.Log("filter run done (%d steps)", int(fr.FilterT*Tl/f.Dt))
		}
	}
	W00 := f.RegionEnergy(func(x, y, z float64) bool { return true })
	var snap [3][]float64
	snapT := -1.0
	for f.Time() < stopRD {
		f.Step()
		E, _ := f.FieldAt(fr.A, 0, 0)
		probeT, probe = append(probeT, f.Time()), append(probe, E[1])
		if f.Steps%every == 0 {
			sampT = append(sampT, f.Time())
			se := make([]float32, nPts)
			sb := make([]float32, nPts)
			for q, p := range g.pts {
				if p.isB {
					sb[q] = float32(f.B[p.c][g.idx(p.i, p.j, p.k)])
				} else {
					se[q] = float32(f.E[p.c][g.idx(p.i, p.j, p.k)])
				}
			}
			sampE, sampB = append(sampE, se), append(sampB, sb)
			if f.Steps%(every*2) == 0 {
				times = append(times, f.Time())
				winside = append(winside, f.RegionEnergy(inFn))
				l := f.Ledger()
				res.RingDown = append(res.RingDown, FieldSample{T: f.Time() / Tl, Inside: winside[len(winside)-1],
					Absorbed: l.Absorbed, Probe: E[1], Residual: l.Residual()})
			}
		}
	}
	res.RDResidual = f.Ledger().Residual() / W00
	// f0 from the zero crossings of the probe after 2 lumped periods.
	var zc []float64
	for k := 1; k < len(probe); k++ {
		if probeT[k] > 0.5*Tl && probe[k-1] < 0 && probe[k] >= 0 {
			zc = append(zc, probeT[k-1]+(probeT[k]-probeT[k-1])*(-probe[k-1])/(probe[k]-probe[k-1]))
		}
	}
	if len(zc) < 3 {
		return res, errFieldRun("ring-down too short to find f0")
	}
	res.F0 = float64(len(zc)-1) / (zc[len(zc)-1] - zc[0])
	w := 2 * math.Pi * res.F0
	T := 1 / res.F0
	// Q from ln W(t) inside the surface after 0.5 τ_L.
	var xs, ys []float64
	for k := range times {
		if times[k] > 0.5*tauL {
			xs, ys = append(xs, times[k]), append(ys, math.Log(winside[k]))
		}
	}
	slope, icpt := fitLine(xs, ys)
	res.Q = -w / slope
	res.W0 = math.Exp(icpt)
	gam := w / (2 * res.Q)
	// Q from the probe's peaks.
	var pt, py []float64
	for k := 1; k+1 < len(probe); k++ {
		a := math.Abs(probe[k])
		if probeT[k] > 0.5*tauL && a > math.Abs(probe[k-1]) && a >= math.Abs(probe[k+1]) {
			pt, py = append(pt, probeT[k]), append(py, math.Log(a))
		}
	}
	ps, _ := fitLine(pt, py)
	res.QProbe = -w / (2 * ps)
	if fr.Log != nil {
		fr.Log("ring-down: f0 %.5g Hz, Q %.4g (probe %.4g), W0 %.4g J, %d crossings", res.F0, res.Q, res.QProbe, res.W0, len(zc))
	}
	// Snapshot of E at an E maximum (probe peak) late in the ring-down, for the
	// correlation: re-run would be costly, so take the last stored state at
	// a probe peak by stepping a few more steps.
	for k := 0; k < int(1.5*T/f.Dt); k++ {
		Ea, _ := f.FieldAt(fr.A, 0, 0)
		f.Step()
		Eb, _ := f.FieldAt(fr.A, 0, 0)
		f.Step()
		Ec, _ := f.FieldAt(fr.A, 0, 0)
		if math.Abs(Eb[1]) > math.Abs(Ea[1]) && math.Abs(Eb[1]) >= math.Abs(Ec[1]) && snapT < 0 {
			// Eb was at the previous step; Ec is current. Close enough (one
			// step is 1/400 of a period): take the current field.
			for c := 0; c < 3; c++ {
				snap[c] = append([]float64(nil), f.E[c]...)
			}
			snapT = f.Time()
			break
		}
	}
	// Phasors of the outgoing field on the surface over [0.5 τ_L, end]:
	// Ê = 2·Σ v e^{−iωt} Δt / Σ e^{−γt} Δt.
	pe := make([]complex128, nPts)
	var norm float64
	dtS := float64(every) * f.Dt
	for s, ts := range sampT {
		if ts < 0.5*tauL {
			continue
		}
		norm += math.Exp(-gam*ts) * dtS
		ce := cmplx.Exp(complex(0, -w*ts)) * complex(dtS, 0)
		cb := cmplx.Exp(complex(0, -w*(ts+f.Dt/2))) * complex(dtS, 0)
		for q := range g.pts {
			if g.pts[q].isB {
				pe[q] += complex(float64(sampB[s][q]), 0) * cb
			} else {
				pe[q] += complex(float64(sampE[s][q]), 0) * ce
			}
		}
	}
	sampE, sampB = nil, nil
	for q := range pe {
		pe[q] *= complex(2/norm, 0)
	}
	// Phase: the drive stops at an E maximum inside the ring (the probe's
	// phasor real and positive), where the interior B of the mode is zero.
	var pp complex128
	for k := range probe {
		if probeT[k] >= 0.5*tauL {
			pp += complex(probe[k], 0) * cmplx.Exp(complex(0, -w*probeT[k]))
		}
	}
	rot := cmplx.Exp(complex(0, -cmplx.Phase(pp)))
	for q := range pe {
		pe[q] *= rot
	}

	// --- 2. implosion then explosion ---
	tp := math.Ceil(fr.ImplodeTE/(2*gam)/T) * T
	stop := tp + fr.ExplodeTE/(2*gam)
	// The declared wave: the time-reversed mode, conj(Ê)·e^{(iω+γ)(t−tp)},
	// −conj(B̂)·…, until tp; zero after (the drive switched off). The
	// factor is cached per time value (every source asks at the same two
	// times each step).
	var lastT = math.NaN()
	var fac complex128
	declared := func(q int, t float64) float64 {
		x := t - tp
		if x > 0 {
			return 0
		}
		if t != lastT {
			lastT, fac = t, cmplx.Exp(complex(gam*x, w*x))
		}
		v := real(cmplx.Conj(pe[q]) * fac)
		if g.pts[q].isB {
			return -v
		}
		return v
	}
	run2 := func(inc func(q int, t float64) float64, start, stop float64, onSample func(f *em.FDTD, s *FieldSample), onStep func(f *em.FDTD)) ([]FieldSample, *em.FDTD, error) {
		f, err := g.fdtd()
		if err != nil {
			return nil, nil, err
		}
		f.SetStartTime(start)
		area := g.dx * g.dx
		for m, ts := range g.eTerms {
			ts := ts
			ed := g.eEdges[m]
			cur := func(t float64) float64 {
				var corr float64
				for _, q := range ts {
					corr += q.kappa * inc(q.pt, t)
				}
				return -corr / em.Mu0 * area
			}
			if err := f.AddSource(em.EdgeSource{Comp: ed[0], I: ed[1], J: ed[2], K: ed[3], Current: cur}); err != nil {
				return nil, nil, err
			}
		}
		for m, ts := range g.bTerms {
			ts := ts
			fc := g.bFaces[m]
			mm := func(t float64) float64 {
				var corr float64
				for _, q := range ts {
					corr += q.kappa * inc(q.pt, t)
				}
				return corr
			}
			if err := f.AddMagneticSource(em.FaceSource{Comp: fc[0], I: fc[1], J: fc[2], K: fc[3], M: mm}); err != nil {
				return nil, nil, err
			}
		}
		var out []FieldSample
		var cin, cout, prev float64
		n := 0
		for f.Time() < stop {
			f.Step()
			wk := f.Ledger().SourceWork
			if d := wk - prev; d > 0 {
				cin += d
			} else {
				cout -= d
			}
			prev = wk
			n++
			if n%(2*every) == 0 || f.Time() >= stop {
				l := f.Ledger()
				Ep, _ := f.FieldAt(fr.A, 0, 0)
				s := FieldSample{T: (f.Time() - start) / T, Inside: f.RegionEnergy(inFn), Outside: f.RegionEnergy(outFn),
					Absorbed: l.Absorbed, SurfIn: cin, SurfOut: cout, Probe: Ep[1], Residual: l.Residual()}
				if onSample != nil {
					onSample(f, &s)
				}
				out = append(out, s)
			}
			if onStep != nil {
				onStep(f)
			}
		}
		return out, f, nil
	}
	start := tp - fr.ImplodeTE/(2*gam)
	start = tp - math.Floor((tp-start)/f.Dt)*f.Dt // whole steps to tp
	res.EInc = res.W0 * (1 - math.Exp(-2*gam*(tp-start)))
	peakDone, quarterDone := false, false
	planeX := []float64{}
	for x := 0.0; x <= 2.4*fr.A; x += g.dx / 2 {
		planeX = append(planeX, x)
	}
	planeZ := []float64{}
	for z := -1.2 * fr.A; z <= 1.2*fr.A+1e-12; z += g.dx / 2 {
		planeZ = append(planeZ, z)
	}
	capture := func(f *em.FDTD, quarter bool) {
		k := 0
		for _, x := range planeX {
			for _, z := range planeZ {
				E, B := f.FieldAt(x, 0, z)
				if !quarter {
					res.Plane = append(res.Plane, PlanePoint{X: x, Z: z, Ex: E[0], Ey: E[1], Ez: E[2], BxPeak: B[0], ByP: B[1], BzPea: B[2], InRing: fr.inTorus(x, 0, z)})
				} else {
					p := &res.Plane[k]
					p.Bx, p.By, p.Bz = B[0], B[1], B[2]
					p.ExQ, p.EyQ, p.EzQ = E[0], E[1], E[2]
				}
				k++
			}
		}
	}
	samples, f2, err := run2(declared, start, stop, func(f *em.FDTD, s *FieldSample) {
		x := f.Time() - tp
		if x <= 0 {
			s.ExactInc = res.W0 * (math.Exp(2*gam*x) - math.Exp(-2*gam*(tp-start)))
		} else {
			s.ExactInc = res.EInc
		}
		s.DeclaredIncoming = x <= 0
	}, func(f *em.FDTD) {
		if !quarterDone && peakDone && f.Time() >= tp+T/4-f.Dt/2 {
			quarterDone = true
			fr.quarter(g, f, &res, capture)
		}
		if !peakDone && f.Time() >= tp-f.Dt/2 {
			peakDone = true
			l := f.Ledger()
			res.Peak = f.RegionEnergy(inFn)
			res.PeakT = (f.Time() - start) / T
			res.LeakAtPeak = f.RegionEnergy(outFn) + l.Absorbed
			res.Delivered = l.SourceWork
			// E_φ share of the electric energy in the tube.
			var eph, eall float64
			for x := -fr.A - fr.R; x <= fr.A+fr.R; x += g.dx {
				for y := -fr.A - fr.R; y <= fr.A+fr.R; y += g.dx {
					for z := -fr.R; z <= fr.R; z += g.dx {
						if !fr.inTorus(x, y, z) {
							continue
						}
						E, _ := f.FieldAt(x, y, z)
						rho := math.Hypot(x, y)
						ephi := (-E[0]*y + E[1]*x) / rho
						eph += ephi * ephi
						eall += E[0]*E[0] + E[1]*E[1] + E[2]*E[2]
					}
				}
			}
			res.EPhiFrac = eph / eall
			// correlation with the ring-down's E at an E maximum (ε-weighted)
			if snapT > 0 {
				var s12, s11, s22 float64
				for c := 0; c < 3; c++ {
					for k := 0; k < g.n; k++ {
						for j := 0; j < g.n; j++ {
							for i := 0; i < g.n; i++ {
								p := [3]float64{float64(i), float64(j), float64(k)}
								p[c] += 0.5
								x, y, z := -g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx
								if !inFn(x, y, z) {
									continue
								}
								er, _ := g.matFunc(x, y, z)
								n := g.idx(i, j, k)
								a, b := snap[c][n], f.E[c][n]
								s12 += er * a * b
								s11 += er * a * a
								s22 += er * b * b
							}
						}
					}
				}
				res.Corr = math.Abs(s12) / math.Sqrt(s11*s22)
			}
			capture(f, false)
		}
	})
	if err != nil {
		return res, err
	}
	res.Implosion = samples
	if fr.Log != nil {
		fr.Log("implosion and explosion done")
	}
	l := f2.Ledger()
	res.Residual = l.Residual() / res.EInc
	last := samples[len(samples)-1]
	res.LeakEnd = last.Outside + last.Absorbed
	res.SurfInTotal, res.SurfOutTotal = last.SurfIn, last.SurfOut
	res.InsideEnd = last.Inside

	// --- 3. constant burst, same incident energy, not continued outgoing ---
	ramp := T
	length := math.Round(2/(2*gam)/T) * T // 2 τ_E
	bStart := start
	bEnd := bStart + length
	shape := func(t float64) float64 {
		x := t - bStart
		switch {
		case x < 0 || x > length:
			return 0
		case x < ramp:
			return 0.5 * (1 - math.Cos(math.Pi*x/ramp))
		case x > length-ramp:
			return 0.5 * (1 - math.Cos(math.Pi*(length-x)/ramp))
		}
		return 1
	}
	var s2 float64
	h := T / 200
	for t := bStart; t < bEnd; t += h {
		v := shape(t + h/2)
		s2 += v * v * h
	}
	amp := math.Sqrt(res.EInc / (2 * gam * res.W0 * s2))
	var lastB = math.NaN()
	var facB complex128
	burst := func(q int, t float64) float64 {
		if t != lastB {
			lastB, facB = t, cmplx.Exp(complex(0, w*(t-bEnd)))*complex(amp*shape(t), 0)
		}
		if facB == 0 {
			return 0
		}
		v := real(cmplx.Conj(pe[q]) * facB)
		if g.pts[q].isB {
			return -v
		}
		return v
	}
	res.BurstEInc = res.EInc
	bs, f3, err := run2(burst, bStart, bEnd+2/(2*gam), func(f *em.FDTD, s *FieldSample) {
		if f.Time() <= bEnd {
			res.BurstLeakDrive = s.Outside + s.Absorbed
		}
		res.BurstPeak = math.Max(res.BurstPeak, s.Inside)
	}, nil)
	if err != nil {
		return res, err
	}
	res.Burst = bs
	res.BurstResidual = f3.Ledger().Residual() / res.EInc
	return res, nil
}

// quarter takes the snapshot a quarter period after the peak (B threads
// the hole a quarter period after the E maximum).
func (fr FieldRing) quarter(g *grid, f *em.FDTD, res *FieldResult, capture func(*em.FDTD, bool)) {
	capture(f, true)
	var bph, ball float64
	for x := -fr.A - fr.R; x <= fr.A+fr.R; x += g.dx {
		for y := -fr.A - fr.R; y <= fr.A+fr.R; y += g.dx {
			for z := -fr.R; z <= fr.R; z += g.dx {
				if !fr.inTorus(x, y, z) {
					continue
				}
				_, B := f.FieldAt(x, y, z)
				rho := math.Hypot(x, y)
				b := (-B[0]*y + B[1]*x) / rho
				bph += b * b
				ball += B[0]*B[0] + B[1]*B[1] + B[2]*B[2]
			}
		}
	}
	res.BPhiFrac = bph / ball
	_, B := f.FieldAt(0, 0, 0)
	res.BzCentre = B[2]
	_, B = f.FieldAt(2*fr.A, 0, 0)
	res.BzOutside = B[2]
}

type errFieldRun string

func (e errFieldRun) Error() string { return "field run: " + string(e) }

// fitLine returns slope and intercept of the least-squares line.
func fitLine(x, y []float64) (slope, icpt float64) {
	n := float64(len(x))
	var sx, sy, sxx, sxy float64
	for k := range x {
		sx += x[k]
		sy += y[k]
		sxx += x[k] * x[k]
		sxy += x[k] * y[k]
	}
	slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
	return slope, (sy - slope*sx) / n
}
