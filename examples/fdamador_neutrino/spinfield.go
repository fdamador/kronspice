// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 7b: the spinning field (H40) in the field solver (em/fdtd.go), on
// the torus of Part 5b (εr = 100, a = 10 mm, r = 2.5 mm, lossless).
//
// In vacuum a field "spins" only by carrying circulating energy flow,
// S = E×H. Angular momentum about the axis (Abraham form, the vacuum
// momentum density ε0μ0·E×H = ε0·E×B for μ = μ0) is
//
//	J_z = ε0·∫ ρ·(E×B)_φ dV ,
//
// beside the Minkowski form ∫ ρ·(D×B)_φ dV (the two differ inside a
// medium; the Abraham–Minkowski question is not settled here). For a
// field made of modes with azimuthal orders m_k and energies U_k at
// frequencies ω_k, J_z = Σ m_k·U_k/ω_k (the angular momentum of a mode is
// m/ω per unit energy, ħm per quantum ħω; J. D. Jackson, Classical
// Electrodynamics, 3rd ed., §9.7 and problem 7.27). So J_z·ω/U is the
// energy-weighted mean of whole numbers m: ½ per quantum is not available
// to one coherent classical mode, only to a mixture.
//
// Helicity (electric–magnetic, J. L. Trueba, A. F. Rañada, Eur. J. Phys.
// 17, 141, 1996; R. P. Cameron, S. M. Barnett, A. M. Yao, New J. Phys. 14,
// 053050, 2012), in SI with B = curl A, D = −curl C, ∂A/∂t = −E,
// ∂C/∂t = −H (temporal gauges, both potentials carried beside the
// solver):
//
//	h = (1/2c)·∫ (A·B/μ0 − C·D/ε0) dV,
//
// normalized so that h·ω/U = ±1 for a circularly polarized wave (tested).
// It is conserved in vacuum and, in a medium, only when ε/μ is uniform
// (duality symmetry; I. Fernandez-Corbaton et al., Phys. Rev. Lett. 111,
// 060401, 2013). The dielectric ring breaks that symmetry, so its h
// need not be conserved. ∫E·B is reported beside it.
//
// The fields compared:
//
//	E-ring:   TE01δ, E round the ring, B poloidal (H28) — m = 0;
//	B-ring:   the dual pattern, B round the ring, E poloidal (H23's
//	          ring-current picture in the same dielectric torus) — m = 0;
//	E + B:    the two in quadrature, equal energies: azimuthal E×B flow
//	          appears beside the E-ring (the reading "spin about the axis");
//	m = 1:    the ring's first azimuthal harmonic, standing (cos φ) and
//	          travelling (cos φ ± i·sin φ): a pattern that turns.
//
// Each mode is made as in Part 5b: a start field excites it, a first run
// finds its frequency (peak of the probe's spectrum), a second run's
// discrete Fourier transform at that frequency is the mode's field. The
// sin φ copy of the m = 1 mode is the cos φ field turned by 90° on the
// grid (exact: the grid and the staircased torus have that symmetry).
//
// Approximations: staircased torus on r/Cells cells (r/2 by default), the
// absorber at 3a, single-frequency phasors (other modes filtered out),
// the Abraham/Minkowski split stated, not resolved.

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
)

// SpinField runs Part 7b on one grid.
type SpinField struct {
	Ring    FieldRing
	ProbeT  float64 // probe run, lumped periods
	FilterT float64 // filter run, lumped periods (DFT over [2T, FilterT·T])
	// ProbeLong, FilterLong replace them for the m = 1 family, whose two
	// lowest odd lattice modes lie close together (high Q, so the longer
	// window resolves them).
	ProbeLong, FilterLong float64
	RunT                  float64 // ring-down of each combination, lumped periods
	Every                 int     // diagnostics every this many steps (0: a sixteenth of a lumped period)
	Log                   func(format string, args ...any)
}

// DefaultSpinField is the r/2 grid with Part 5b's region.
func DefaultSpinField() SpinField {
	fr := DefaultFieldRing()
	fr.Cells = 2
	return SpinField{Ring: fr, ProbeT: 10, FilterT: 14, ProbeLong: 40, FilterLong: 70, RunT: 20}
}

// Mode is a filtered mode field: E(t) = Re(Ê·e^{iωt}), B likewise.
type Mode struct {
	Name    string
	F       float64 // Hz
	E, B    [3][]complex128
	Energy  float64 // field energy of Re(Ê), Re(B̂) set at t = 0 [J]
	Turn90  float64 // energy-weighted overlap with itself turned 90° about z: cos(m·90°) for a pure order m
	Turn180 float64 // with itself turned 180°: cos(m·180°)
	Inside  float64 // share of the mode's energy inside the torus
}

// SpinRun is one ring-down with diagnostics.
type SpinRun struct {
	Name                         string
	Omega                        float64   // reference ω (energy-weighted: U/Σ(U_k/ω_k) at the start)
	T                            []float64 // lumped periods
	U, UOut, Abs                 []float64 // energy inside the surface; outside it; absorbed [J]
	Jz, JzMink                   []float64 // Abraham and Minkowski angular momentum inside [J·s]
	JzAbs                        []float64 // ε0∫ρ|(E×B)_φ| inside: the circulation, whatever its sign [J·s]
	EB                           []float64 // ∫E·B inside [V·T·m²]
	H                            []float64 // helicity inside [J·s]
	JzW, HW                      []float64 // J_z·ω/U and h·ω/U
	JzMW                         []float64 // Minkowski J_z·ω/U
	Q                            float64   // from the decay of U after 2 periods
	Residual                     float64   // field ledger residual / starting energy
	CurlARes                     float64   // max |curl A − B| / max |B| inside at the start (potential consistency)
	JzW0, HW0, JzWEnd, HWEnd     float64   // "end": the last sample with U ≥ 1e-3 of the start
	JzWMin, JzWMax, HWMin, HWMax float64   // over the samples up to Last
	JzMW0, JzMWMin, JzMWMax      float64   // Minkowski: start and range
	Last                         int
}

// spinGrid is the region of Part 5b without the Huygens surface.
type spinGrid struct {
	fr   FieldRing
	n    int
	dx   float64
	half float64
	b    float64
}

func (sf SpinField) grid() spinGrid {
	fr := sf.Ring
	dx := fr.R / fr.Cells
	m := int(math.Ceil(fr.HalfA * fr.A / dx))
	n := 2*m + 2*fr.Absorb + 1
	return spinGrid{fr: fr, n: n, dx: dx, half: float64(n-1) / 2 * dx, b: (math.Round(fr.SurfA*fr.A/dx) + 0.25) * dx}
}

func (g spinGrid) mat(x, y, z float64) (float64, float64) {
	if g.fr.inTorus(x, y, z) {
		return g.fr.EpsR, 0
	}
	return 1, 0
}

func (g spinGrid) fdtd() (*em.FDTD, error) {
	return em.NewFDTD(em.FDTDConfig{N: [3]int{g.n, g.n, g.n}, D: [3]float64{g.dx, g.dx, g.dx},
		Origin: [3]float64{-g.half, -g.half, -g.half}, Absorber: g.fr.Absorb, Ledger: true, Material: g.mat})
}

func (g spinGrid) idx(i, j, k int) int { return i + g.n*(j+g.n*k) }

func (g spinGrid) inside(x, y, z float64) bool {
	return math.Abs(x) < g.b && math.Abs(y) < g.b && math.Abs(z) < g.b
}

// zInt is −∫_{a−r}^{ρ} bump(ρ', z) dρ' for the start fields.
func (g spinGrid) zInt(rho, z float64) float64 {
	fr := g.fr
	lo := math.Max(0, fr.A-fr.R)
	if rho <= lo || math.Abs(z) >= fr.R {
		return 0
	}
	hi := math.Min(rho, fr.A+fr.R)
	const m = 64
	h := (hi - lo) / m
	var s float64
	for q := 0; q < m; q++ {
		r := lo + (float64(q)+0.5)*h
		sd := math.Hypot(r-fr.A, z)
		if sd < fr.R {
			c := math.Cos(math.Pi * sd / (2 * fr.R))
			s += c * c * h
		}
	}
	return -s
}

// start sets the start field of a mode family: "E" (D_φ bump), "B" (B_φ
// bump), "m1" (D = curl(C_z·cos φ·ẑ): D_φ ∝ bump·cos φ, divergence-free).
func (g spinGrid) start(f *em.FDTD, kind string) {
	switch kind {
	case "B":
		f.SetBFromPotential(func(x, y, z float64) [3]float64 { return [3]float64{0, 0, g.zInt(math.Hypot(x, y), z)} })
	default:
		f.SetEFromPotential(func(x, y, z float64) [3]float64 {
			rho := math.Hypot(x, y)
			c := g.zInt(rho, z)
			if kind == "m1" && rho > 0 {
				c *= x / rho
			}
			return [3]float64{0, 0, c}
		})
		g.divideEps(f)
	}
}

func (g spinGrid) divideEps(f *em.FDTD) {
	for c := 0; c < 3; c++ {
		for k := 0; k < g.n; k++ {
			for j := 0; j < g.n; j++ {
				for i := 0; i < g.n; i++ {
					p := [3]float64{float64(i), float64(j), float64(k)}
					p[c] += 0.5
					er, _ := g.mat(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx)
					f.E[c][g.idx(i, j, k)] /= er
				}
			}
		}
	}
}

// probe is E_φ (or B_φ for the B family) at (a, 0, 0).
func (g spinGrid) probe(f *em.FDTD, kind string) float64 {
	E, B := f.FieldAt(g.fr.A, 0, 0)
	if kind == "B" {
		return B[1]
	}
	return E[1]
}

// FindMode runs the start field, finds the strongest spectral peak of the
// probe between fLo and fHi, and filters the field at it.
func (sf SpinField) FindMode(name, kind string, fLo, fHi float64) (Mode, error) {
	g := sf.grid()
	e := DefaultERing()
	Tl := 2 * math.Pi / e.Omega0()
	f, err := g.fdtd()
	if err != nil {
		return Mode{}, err
	}
	g.start(f, kind)
	pT, fT := sf.ProbeT, sf.FilterT
	if kind == "m1" && sf.ProbeLong > 0 {
		pT, fT = sf.ProbeLong, sf.FilterLong
	}
	var ts, ps []float64
	for f.Time() < pT*Tl {
		f.Step()
		if f.Time() > 1*Tl {
			ts, ps = append(ts, f.Time()), append(ps, g.probe(f, kind))
		}
	}
	spec := func(fr float64) float64 {
		var s complex128
		for k := range ts {
			// Hann window over the record
			w := 0.5 - 0.5*math.Cos(2*math.Pi*float64(k)/float64(len(ts)-1))
			s += complex(w*ps[k], 0) * cmplx.Exp(complex(0, -2*math.Pi*fr*ts[k]))
		}
		return cmplx.Abs(s)
	}
	best, bf := -1.0, fLo
	for fr := fLo; fr <= fHi; fr += (fHi - fLo) / 600 {
		if v := spec(fr); v > best {
			best, bf = v, fr
		}
	}
	lo, hi := bf-(fHi-fLo)/600, bf+(fHi-fLo)/600
	gr := (math.Sqrt(5) - 1) / 2
	for it := 0; it < 60; it++ {
		m1, m2 := hi-gr*(hi-lo), lo+gr*(hi-lo)
		if spec(m1) > spec(m2) {
			hi = m2
		} else {
			lo = m1
		}
	}
	md := Mode{Name: name, F: (lo + hi) / 2}
	if sf.Log != nil {
		sf.Log("%s: probe peak at %.5g GHz", name, md.F/1e9)
	}
	// filter run
	f, err = g.fdtd()
	if err != nil {
		return md, err
	}
	g.start(f, kind)
	w := 2 * math.Pi * md.F
	every := int(math.Max(1, math.Round(Tl/16/f.Dt)))
	for c := 0; c < 3; c++ {
		md.E[c] = make([]complex128, len(f.E[c]))
		md.B[c] = make([]complex128, len(f.B[c]))
	}
	for f.Time() < fT*Tl {
		f.Step()
		if f.Steps%every != 0 || f.Time() < 2*Tl {
			continue
		}
		ce := cmplx.Exp(complex(0, -w*f.Time()))
		cb := cmplx.Exp(complex(0, -w*(f.Time()+f.Dt/2)))
		for c := 0; c < 3; c++ {
			ec, bc, pe, pb := f.E[c], f.B[c], md.E[c], md.B[c]
			for n := range ec {
				pe[n] += complex(ec[n], 0) * ce
				pb[n] += complex(bc[n], 0) * cb
			}
		}
	}
	// phase: the E-type modes start at an E maximum (Σε Ê·Ê real and
	// positive: B = 0 inside), the B-ring at a B maximum (Σ B̂·B̂ real and
	// positive: E = 0), so "in phase" and "in quadrature" have a meaning.
	var ss complex128
	for c := 0; c < 3; c++ {
		if kind == "B" {
			for _, z := range md.B[c] {
				ss += z * z
			}
		} else {
			for _, z := range md.E[c] {
				ss += z * z // the ε weight matters little for the phase; inside the ring ε is uniform
			}
		}
	}
	rot := cmplx.Exp(complex(0, -cmplx.Phase(ss)/2))
	for c := 0; c < 3; c++ {
		for n := range md.E[c] {
			md.E[c][n] *= rot
			md.B[c][n] *= rot
		}
	}
	// normalize: unit energy when set alone at t = 0
	f, err = g.fdtd()
	if err != nil {
		return md, err
	}
	g.set(f, []Mode{md}, []complex128{1})
	md.Energy = f.RegionEnergy(func(x, y, z float64) bool { return true })
	s := complex(1/math.Sqrt(md.Energy), 0)
	for c := 0; c < 3; c++ {
		for n := range md.E[c] {
			md.E[c][n] *= s
			md.B[c][n] *= s
		}
	}
	md.Energy = 1
	md.Inside = f.RegionEnergy(func(x, y, z float64) bool { return g.fr.inTorus(x, y, z) }) / f.RegionEnergy(func(x, y, z float64) bool { return true })
	r90 := g.Rot90(md)
	md.Turn90 = g.overlap(md, r90)
	md.Turn180 = g.overlap(md, g.Rot90(r90))
	if sf.Log != nil {
		sf.Log("%s: filtered at %.5g GHz; overlap with itself turned 90° %.4f, 180° %.4f", name, md.F/1e9, md.Turn90, md.Turn180)
	}
	return md, nil
}

// overlap is Re Σ(ε Ê*·Ê' + B̂*·B̂'/μ0) / Σ(ε|Ê|² + |B̂|²/μ0).
func (g spinGrid) overlap(a, b Mode) float64 {
	var num, den float64
	for c := 0; c < 3; c++ {
		for k := 0; k < g.n; k++ {
			for j := 0; j < g.n; j++ {
				for i := 0; i < g.n; i++ {
					q := g.idx(i, j, k)
					p := [3]float64{float64(i), float64(j), float64(k)}
					p[c] += 0.5
					er, _ := g.mat(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx)
					num += Eps0*er*real(cmplx.Conj(a.E[c][q])*b.E[c][q]) + real(cmplx.Conj(a.B[c][q])*b.B[c][q])/Mu0
					den += Eps0*er*sq(cmplx.Abs(a.E[c][q])) + sq(cmplx.Abs(a.B[c][q]))/Mu0
				}
			}
		}
	}
	return num / den
}

// Rot90 is the mode turned by +90° about z: F'(r) = R·F(R⁻¹r).
func (g spinGrid) Rot90(m Mode) Mode {
	n := g.n
	out := Mode{Name: m.Name + " (turned 90°)", F: m.F, Energy: m.Energy}
	for c := 0; c < 3; c++ {
		out.E[c] = make([]complex128, len(m.E[c]))
		out.B[c] = make([]complex128, len(m.B[c]))
	}
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				q := g.idx(i, j, k)
				if i <= n-2 {
					out.E[0][q] = -m.E[1][g.idx(j, n-2-i, k)]
					out.B[1][q] = m.B[0][g.idx(j, n-2-i, k)]
					out.B[2][q] = m.B[2][g.idx(j, n-2-i, k)]
				}
				out.E[1][q] = m.E[0][g.idx(j, n-1-i, k)]
				out.E[2][q] = m.E[2][g.idx(j, n-1-i, k)]
				out.B[0][q] = -m.B[1][g.idx(j, n-1-i, k)]
			}
		}
	}
	return out
}

// set puts Σ a_k·mode_k at t = 0: E = Re Σ a_k Ê_k, B(dt/2) = Re Σ a_k B̂_k e^{iω_k dt/2}.
func (g spinGrid) set(f *em.FDTD, ms []Mode, a []complex128) {
	for c := 0; c < 3; c++ {
		for n := range f.E[c] {
			var e, b complex128
			for k, m := range ms {
				e += a[k] * m.E[c][n]
				b += a[k] * m.B[c][n] * cmplx.Exp(complex(0, 2*math.Pi*m.F*f.Dt/2))
			}
			f.E[c][n] = real(e)
			f.B[c][n] = real(b)
		}
	}
}

// RingDown runs Σ a_k·mode_k for RunT lumped periods with the diagnostics.
func (sf SpinField) RingDown(name string, ms []Mode, a []complex128) (SpinRun, error) {
	g := sf.grid()
	e := DefaultERing()
	Tl := 2 * math.Pi / e.Omega0()
	f, err := g.fdtd()
	if err != nil {
		return SpinRun{}, err
	}
	g.set(f, ms, a)
	// potentials: A at the B time (dt/2), C at the E time (0)
	var A, C, D [3][]float64
	for c := 0; c < 3; c++ {
		A[c] = make([]float64, len(f.E[c]))
		C[c] = make([]float64, len(f.B[c]))
		D[c] = make([]float64, len(f.E[c]))
		for n := range A[c] {
			var av, cv complex128
			for k, m := range ms {
				w := 2 * math.Pi * m.F
				av += a[k] * complex(0, 1/w) * m.E[c][n] * cmplx.Exp(complex(0, w*f.Dt/2))
				cv += a[k] * complex(0, 1/(Mu0*w)) * m.B[c][n]
			}
			A[c][n] = real(av)
			if n < len(C[c]) {
				C[c][n] = real(cv)
			}
		}
	}
	eps := make([][]float64, 3) // ε per E edge
	for c := 0; c < 3; c++ {
		eps[c] = make([]float64, len(f.E[c]))
		for k := 0; k < g.n; k++ {
			for j := 0; j < g.n; j++ {
				for i := 0; i < g.n; i++ {
					p := [3]float64{float64(i), float64(j), float64(k)}
					p[c] += 0.5
					er, _ := g.mat(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx)
					eps[c][g.idx(i, j, k)] = Eps0 * er
				}
			}
		}
	}
	run0CurlA := g.curlARes(f, A)
	all := func(x, y, z float64) bool { return true }
	out := func(x, y, z float64) bool { return !g.inside(x, y, z) }
	W0 := f.RegionEnergy(all)
	run := SpinRun{Name: name, CurlARes: run0CurlA}
	var den float64
	for k, m := range ms {
		den += sq(cmplx.Abs(a[k])) * m.Energy / (2 * math.Pi * m.F)
	}
	var num float64
	for k := range ms {
		num += sq(cmplx.Abs(a[k])) * ms[k].Energy
	}
	run.Omega = num / den
	every := sf.Every
	if every <= 0 {
		every = int(math.Max(1, math.Round(Tl/16/f.Dt)))
	}
	bOld := [3][]float64{make([]float64, len(f.B[0])), make([]float64, len(f.B[1])), make([]float64, len(f.B[2]))}
	diag := func() {
		// fields at the E time n+1: E, B̄ = (B^{n+½} + B^{n+3/2})/2,
		// Ā = A^{n+3/2} + dt·E/2, C at n+1, D = εE.
		var bm, am [3][]float64
		for c := 0; c < 3; c++ {
			bm[c] = make([]float64, len(f.B[c]))
			am[c] = make([]float64, len(f.E[c]))
			for n := range bm[c] {
				bm[c][n] = 0.5 * (bOld[c][n] + f.B[c][n])
			}
			for n := range am[c] {
				am[c][n] = A[c][n] + 0.5*f.Dt*f.E[c][n]
				D[c][n] = eps[c][n] * f.E[c][n]
			}
		}
		jz, jzm, jza, eb, hh := g.integrals(f.E, bm, D, am, C)
		u := f.RegionEnergy(g.inside)
		run.T = append(run.T, f.Time()/Tl)
		run.U = append(run.U, u)
		run.UOut = append(run.UOut, f.RegionEnergy(out))
		run.Abs = append(run.Abs, f.Ledger().Absorbed)
		run.Jz, run.JzMink, run.JzAbs, run.EB, run.H = append(run.Jz, jz), append(run.JzMink, jzm), append(run.JzAbs, jza), append(run.EB, eb), append(run.H, hh)
		run.JzW = append(run.JzW, jz*run.Omega/u)
		run.JzMW = append(run.JzMW, jzm*run.Omega/u)
		run.HW = append(run.HW, hh*run.Omega/u)
	}
	for f.Time() < sf.RunT*Tl {
		isDiag := (f.Steps+1)%every == 0
		if isDiag {
			for c := 0; c < 3; c++ {
				copy(bOld[c], f.B[c])
			}
		}
		// C^{n+1} = C^n − dt·H^{n+½}
		for c := 0; c < 3; c++ {
			bc, cc := f.B[c], C[c]
			for n := range cc {
				cc[n] -= f.Dt * bc[n] / Mu0
			}
		}
		f.Step()
		// A^{n+3/2} = A^{n+½} − dt·E^{n+1}
		for c := 0; c < 3; c++ {
			ec, ac := f.E[c], A[c]
			for n := range ac {
				ac[n] -= f.Dt * ec[n]
			}
		}
		if isDiag {
			diag()
		}
	}
	run.Residual = f.Ledger().Residual() / W0
	// Q from ln U after 1.5 lumped periods, while U ≥ 1e-6 of the start
	// (a low-Q mode is gone long before the run ends; what is left then is
	// the box's slow residue, not the mode).
	var xs, ys []float64
	last := 0
	for k := range run.T {
		if run.U[k] >= 1e-3*run.U[0] {
			last = k
		}
		if run.T[k] > 1.5 && run.U[k] >= 1e-6*run.U[0] {
			xs, ys = append(xs, run.T[k]*Tl), append(ys, math.Log(run.U[k]))
		}
	}
	if len(xs) >= 3 {
		sl, _ := fitLine(xs, ys)
		run.Q = -run.Omega / sl
	}
	run.Last = last
	run.JzW0, run.HW0 = run.JzW[0], run.HW[0]
	run.JzWEnd, run.HWEnd = run.JzW[last], run.HW[last]
	run.JzWMin, run.JzWMax, run.HWMin, run.HWMax = math.Inf(1), math.Inf(-1), math.Inf(1), math.Inf(-1)
	run.JzMW0, run.JzMWMin, run.JzMWMax = run.JzMW[0], math.Inf(1), math.Inf(-1)
	for k := 0; k <= last; k++ {
		run.JzMWMin, run.JzMWMax = math.Min(run.JzMWMin, run.JzMW[k]), math.Max(run.JzMWMax, run.JzMW[k])
		run.JzWMin, run.JzWMax = math.Min(run.JzWMin, run.JzW[k]), math.Max(run.JzWMax, run.JzW[k])
		run.HWMin, run.HWMax = math.Min(run.HWMin, run.HW[k]), math.Max(run.HWMax, run.HW[k])
	}
	return run, nil
}

// integrals returns, over the cells inside the surface, J_z (Abraham and
// Minkowski), ε0∫ρ|(E×B)_φ|, ∫E·B and the helicity, with every component
// taken to the cell centre (E, A, D: mean of the four parallel edges; B,
// C: mean of the two parallel faces).
func (g spinGrid) integrals(E, B, D, A, C [3][]float64) (jz, jzm, jza, eb, hh float64) {
	n := g.n
	dv := g.dx * g.dx * g.dx
	sy, sz := n, n*n
	edge := func(arr [3][]float64, c, q int) float64 {
		var o1, o2 int
		switch c {
		case 0:
			o1, o2 = sy, sz
		case 1:
			o1, o2 = 1, sz
		default:
			o1, o2 = 1, sy
		}
		return 0.25 * (arr[c][q] + arr[c][q+o1] + arr[c][q+o2] + arr[c][q+o1+o2])
	}
	face := func(arr [3][]float64, c, q int) float64 {
		o := [3]int{1, sy, sz}[c]
		return 0.5 * (arr[c][q] + arr[c][q+o])
	}
	for k := 0; k < n-1; k++ {
		z := -g.half + (float64(k)+0.5)*g.dx
		if math.Abs(z) >= g.b {
			continue
		}
		for j := 0; j < n-1; j++ {
			y := -g.half + (float64(j)+0.5)*g.dx
			if math.Abs(y) >= g.b {
				continue
			}
			for i := 0; i < n-1; i++ {
				x := -g.half + (float64(i)+0.5)*g.dx
				if math.Abs(x) >= g.b {
					continue
				}
				q := g.idx(i, j, k)
				var e, b, d, av, cv [3]float64
				for c := 0; c < 3; c++ {
					e[c], d[c], av[c] = edge(E, c, q), edge(D, c, q), edge(A, c, q)
					b[c], cv[c] = face(B, c, q), face(C, c, q)
				}
				rho := math.Hypot(x, y)
				if rho > 0 {
					// (E×B)_φ·ρ = x·(E×B)_y − y·(E×B)_x
					sx := e[1]*b[2] - e[2]*b[1]
					syy := e[2]*b[0] - e[0]*b[2]
					l := x*syy - y*sx
					jz += Eps0 * l
					jza += Eps0 * math.Abs(l)
					mx := d[1]*b[2] - d[2]*b[1]
					my := d[2]*b[0] - d[0]*b[2]
					jzm += x*my - y*mx
				}
				eb += e[0]*b[0] + e[1]*b[1] + e[2]*b[2]
				hh += (av[0]*b[0]+av[1]*b[1]+av[2]*b[2])/Mu0 - (cv[0]*d[0]+cv[1]*d[1]+cv[2]*d[2])/Eps0
			}
		}
	}
	return jz * dv, jzm * dv, jza * dv, eb * dv, hh * dv / (2 * C0)
}

// curlARes is max |curl A − B| / max |B| over faces inside the surface,
// with A at the B time (they evolve together; the start phasors fix it).
func (g spinGrid) curlARes(f *em.FDTD, A [3][]float64) float64 {
	n := g.n
	sy, sz := n, n*n
	var worst, bmax float64
	for c := 0; c < 3; c++ {
		for k := 1; k < n-2; k++ {
			for j := 1; j < n-2; j++ {
				for i := 1; i < n-2; i++ {
					p := [3]float64{float64(i) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
					p[c] -= 0.5
					if !g.inside(-g.half+p[0]*g.dx, -g.half+p[1]*g.dx, -g.half+p[2]*g.dx) {
						continue
					}
					q := g.idx(i, j, k)
					var cu float64
					switch c {
					case 0:
						cu = (A[2][q+sy]-A[2][q])/g.dx - (A[1][q+sz]-A[1][q])/g.dx
					case 1:
						cu = (A[0][q+sz]-A[0][q])/g.dx - (A[2][q+1]-A[2][q])/g.dx
					default:
						cu = (A[1][q+1]-A[1][q])/g.dx - (A[0][q+sy]-A[0][q])/g.dx
					}
					worst = math.Max(worst, math.Abs(cu-f.B[c][q]))
					bmax = math.Max(bmax, math.Abs(f.B[c][q]))
				}
			}
		}
	}
	return worst / bmax
}

// HelicityPlaneWave is h·ω/u for a plane wave E = E0(x̂ cos ωt ∓ ŷ sin ωt)
// travelling along z, with A, C from the temporal gauges: the check of the
// normalization (±1 for circular, 0 for linear when lin is true).
func HelicityPlaneWave(sign float64, lin bool, t float64) float64 {
	const e0, w = 1.0, 1.0e9
	c := C0
	ey := sign
	if lin {
		ey = 0
	}
	// at z = 0: E = (cos ωt, −ey·sin ωt, 0)·E0; B = ẑ×E/c.
	E := [3]float64{e0 * math.Cos(w*t), -ey * e0 * math.Sin(w*t), 0}
	B := [3]float64{-E[1] / c, E[0] / c, 0}
	// A = −∫E dt
	A := [3]float64{-e0 * math.Sin(w*t) / w, -ey * e0 * math.Cos(w*t) / w, 0}
	// C = −∫H dt: ∫B dt = ẑ×∫E dt/c = −ẑ×A/c, so C = ẑ×A/(cμ0).
	C := [3]float64{-A[1] / (c * Mu0), A[0] / (c * Mu0), 0}
	D := [3]float64{Eps0 * E[0], Eps0 * E[1], 0}
	h := ((A[0]*B[0]+A[1]*B[1])/Mu0 - (C[0]*D[0]+C[1]*D[1])/Eps0) / (2 * c)
	u := 0.5*Eps0*(E[0]*E[0]+E[1]*E[1]) + (B[0]*B[0]+B[1]*B[1])/(2*Mu0)
	return h * w / u
}

// AzimuthalContent decomposes the mode's E_φ phasor sampled on the ring's
// centre circle (ρ = a, z = 0, 256 points) into azimuthal orders: the
// share of Σ|E_φ|² in each m from −8 to 8.
func (g spinGrid) AzimuthalContent(ms []Mode, a []complex128) map[int]float64 {
	f, err := g.fdtd()
	if err != nil {
		return nil
	}
	fi, _ := g.fdtd()
	g.set(f, ms, a)
	ai := make([]complex128, len(a))
	for k := range a {
		ai[k] = a[k] * complex(0, -1) // Re(−i·z) = Im z
	}
	g.set(fi, ms, ai)
	const np = 256
	ephi := make([]complex128, np)
	for k := 0; k < np; k++ {
		ph := 2 * math.Pi * float64(k) / np
		x, y := g.fr.A*math.Cos(ph), g.fr.A*math.Sin(ph)
		er, _ := f.FieldAt(x, y, 0)
		ei, _ := fi.FieldAt(x, y, 0)
		ephi[k] = complex(-er[0]*math.Sin(ph)+er[1]*math.Cos(ph), -ei[0]*math.Sin(ph)+ei[1]*math.Cos(ph))
	}
	out := map[int]float64{}
	var tot float64
	for m := -8; m <= 8; m++ {
		var c complex128
		for k := range ephi {
			c += ephi[k] * cmplx.Exp(complex(0, -float64(m)*2*math.Pi*float64(k)/np))
		}
		out[m] = sq(cmplx.Abs(c))
		tot += out[m]
	}
	for m := range out {
		out[m] /= tot
	}
	return out
}
