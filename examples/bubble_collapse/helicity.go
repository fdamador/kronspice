// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"

	"github.com/fdamador/kronspice/em"
)

// Part 5: does the flash carry the knot's topology?
//
// If the light of a collapse were a knot of field energy (H29, H30), it
// would carry net optical helicity: more photons of one handedness than of
// the other. Thermal light from random currents carries none on average.
// The measurement a bench can make is the degree of circular polarization
// of the whole flash, integrated over directions and frequencies:
//
//	P_c = Σ V / Σ I,  V = −2·Im(Ê_θ*·Ê_φ),  I = |Ê_θ|² + |Ê_φ|²
//
// summed over the positive-frequency spectrum of the transverse field at
// every probe of a sphere around the source (θ̂ × φ̂ = r̂, so V > 0 is
// spin along the direction of travel: positive helicity). The same
// measurement is made in em/fdtd.go on four sources:
//
//   - the hopfion of examples/knotted_light, with the same P_c computed
//     from its exact field at the same probes (the measurement's twin);
//   - a chiral dipole pair, electric dipole p along z and magnetic dipole
//     m = c·p along z in quadrature, the simplest source of one helicity
//     (|P_c| near 1);
//   - a rotating electric dipole, circular on its axis but of zero net
//     helicity (the two hemispheres cancel): an electric dipole alone
//     cannot carry helicity, which needs electric and magnetic parts
//     together (helicity is the charge of electric–magnetic duality;
//     M. G. Calkin, Am. J. Phys. 33, 958, 1965; R. P. Cameron, S. M. Barnett,
//     A. M. Yao, New J. Phys. 14, 053050, 2012);
//   - the thermal core: random, independent current elements filling a
//     sphere of the collapsed bubble's size relative to the light's
//     wavelength (R_min ≈ 0.79 µm against λ ≈ 450 nm, so 1.75 λ), with
//     random directions, frequencies (0.5–1.5 ω₀) and phases, under one
//     Gaussian envelope: the fluctuating currents of fluctuational
//     electrodynamics (S. M. Rytov, 1953) standing for a hot gas.
//
// Maxwell's equations have no length scale in vacuum, so the run is made
// in units of the central wavelength λ₀ = 8 cells.

// HelicitySetup sets the field runs.
type HelicitySetup struct {
	Lambda0  float64 // central wavelength [m]
	Cells    int     // cells per λ₀
	Core     float64 // thermal core radius [λ₀]
	Probe    float64 // probe-sphere radius [λ₀]
	Half     float64 // interior half-width [λ₀]
	Absorb   int     // absorber cells
	Probes   int     // probes on the sphere
	Dipoles  int     // thermal current elements
	Seeds    int     // thermal realizations
	HopfL0   float64 // hopfion scale [λ₀]
	Envelope float64 // Gaussian envelope σ_t [periods]
}

// DefaultHelicity: λ₀ = 450 nm at 8 cells, core 1.75 λ₀, probes at 6 λ₀.
func DefaultHelicity() HelicitySetup {
	return HelicitySetup{Lambda0: 450e-9, Cells: 8, Core: 1.75, Probe: 6, Half: 6.75, Absorb: 10,
		Probes: 600, Dipoles: 300, Seeds: 6, HopfL0: 0.625, Envelope: 1.5}
}

// HelicityRun is one source's measurement.
type HelicityRun struct {
	Name     string
	Pc       float64 // solver
	PcExact  float64 // exact field at the probes (hopfion only; NaN otherwise)
	Residual float64 // field ledger residual / energy scale
	Captured float64 // energy that crossed the probe sphere / energy put in (≈ fraction reaching the absorber)
	// For the report (-data): per-probe V and I, and |E|² slices on z = 0.
	ProbeV, ProbeI []float64
	Snaps          [][]float32
	SnapT          []float64
}

func (r *HelicityRun) take(x recorded) {
	r.Pc, r.PcExact, r.ProbeV, r.ProbeI, r.Snaps, r.SnapT = x.pc, x.pcExact, x.v, x.in, x.snaps, x.snapT
}

// HelicityResult collects the runs.
type HelicityResult struct {
	Runs            []HelicityRun
	ThermalMean     float64
	ThermalRMS      float64
	ThermalExpected float64 // 1/√(number of independent radiators × 2 polarizations)
}

type probe struct {
	x, y, z float64
	th, ph  [3]float64
}

func (s HelicitySetup) probes() []probe {
	r := s.Probe * s.Lambda0
	ps := make([]probe, s.Probes)
	g := math.Pi * (3 - math.Sqrt(5))
	for i := range ps {
		ct := 1 - 2*(float64(i)+0.5)/float64(s.Probes)
		st := math.Sqrt(1 - ct*ct)
		phi := g * float64(i)
		cp, sp := math.Cos(phi), math.Sin(phi)
		ps[i] = probe{x: r * st * cp, y: r * st * sp, z: r * ct,
			th: [3]float64{ct * cp, ct * sp, -st}, ph: [3]float64{-sp, cp, 0}}
	}
	return ps
}

func (s HelicitySetup) grid() (*em.FDTD, int, error) {
	dx := s.Lambda0 / float64(s.Cells)
	n := 2*int(math.Ceil(s.Half*float64(s.Cells))) + 2*s.Absorb + 1
	half := float64(n-1) / 2 * dx
	f, err := em.NewFDTD(em.FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx},
		Origin: [3]float64{-half, -half, -half}, Absorber: s.Absorb, Ledger: true})
	return f, n, err
}

// fft returns the transform Σ x_j e^{−2πijk/n} of x zero-padded to a power
// of two, bins 0 … n/2 (iterative radix-2, Cooley–Tukey).
func fft(x []float64) []complex128 {
	n := 1
	for n < len(x) {
		n <<= 1
	}
	a := make([]complex128, n)
	for i, v := range x {
		a[i] = complex(v, 0)
	}
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
	for size := 2; size <= n; size <<= 1 {
		w := cmplx.Exp(complex(0, -2*math.Pi/float64(size)))
		for st := 0; st < n; st += size {
			wk := complex(1, 0)
			for k := 0; k < size/2; k++ {
				u, v := a[st+k], a[st+k+size/2]*wk
				a[st+k], a[st+k+size/2] = u+v, u-v
				wk *= w
			}
		}
	}
	return a[:n/2+1]
}

// stokesPer returns V and I at each probe from transverse time series.
func stokesPer(eth, eph [][]float64) (v, in []float64) {
	v, in = make([]float64, len(eth)), make([]float64, len(eth))
	for p := range eth {
		a := fft(eth[p])
		b := fft(eph[p])
		for k := 1; k < len(a); k++ {
			v[p] += -2 * imag(cmplx.Conj(a[k])*b[k])
			in[p] += real(a[k]*cmplx.Conj(a[k])) + real(b[k]*cmplx.Conj(b[k]))
		}
	}
	return
}

// stokes returns ΣV and ΣI over the probes.
func stokes(eth, eph [][]float64) (V, I float64) {
	v, in := stokesPer(eth, eph)
	for p := range v {
		V += v[p]
		I += in[p]
	}
	return
}

// slice returns |E|² on the z = 0 plane at the interior nodes, row-major
// in y then x.
func slice(f *em.FDTD, absorb int) []float32 {
	n := f.N[0]
	out := make([]float32, 0, (n-2*absorb)*(n-2*absorb))
	for j := absorb; j < n-absorb; j++ {
		for i := absorb; i < n-absorb; i++ {
			x := f.Origin[0] + float64(i)*f.D[0]
			y := f.Origin[1] + float64(j)*f.D[1]
			E, _ := f.FieldAt(x, y, 0)
			out = append(out, float32(E[0]*E[0]+E[1]*E[1]+E[2]*E[2]))
		}
	}
	return out
}

// record steps f to tEnd, sampling the transverse E at the probes; exact,
// if not nil, gives the exact E at a point and time for the twin.
func record(f *em.FDTD, ps []probe, tEnd float64, exact func(x, y, z, t float64) [3]float64) (pc, pcExact float64) {
	r := recordFull(f, ps, tEnd, exact, nil, 0)
	return r.pc, r.pcExact
}

type recorded struct {
	pc, pcExact float64
	v, in       []float64   // per-probe V and I (solver)
	snaps       [][]float32 // |E|² on z = 0 at the requested times
	snapT       []float64
}

// recordFull is record with per-probe results and field slices at snapAt.
func recordFull(f *em.FDTD, ps []probe, tEnd float64, exact func(x, y, z, t float64) [3]float64, snapAt []float64, absorb int) (out recorded) {
	takeSnaps := func() {
		for len(out.snaps) < len(snapAt) && f.Time() >= snapAt[len(out.snaps)]-0.5*f.Dt {
			out.snaps = append(out.snaps, slice(f, absorb))
			out.snapT = append(out.snapT, f.Time())
		}
	}
	takeSnaps()
	var eth, eph, xth, xph [][]float64
	eth = make([][]float64, len(ps))
	eph = make([][]float64, len(ps))
	if exact != nil {
		xth = make([][]float64, len(ps))
		xph = make([][]float64, len(ps))
	}
	dot := func(a, b [3]float64) float64 { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
	for f.Time() < tEnd {
		f.Step()
		takeSnaps()
		t := f.Time()
		for i, p := range ps {
			E, _ := f.FieldAt(p.x, p.y, p.z)
			eth[i] = append(eth[i], dot(E, p.th))
			eph[i] = append(eph[i], dot(E, p.ph))
			if exact != nil {
				X := exact(p.x, p.y, p.z, t)
				xth[i] = append(xth[i], dot(X, p.th))
				xph[i] = append(xph[i], dot(X, p.ph))
			}
		}
	}
	out.v, out.in = stokesPer(eth, eph)
	var V, I float64
	for p := range out.v {
		V += out.v[p]
		I += out.in[p]
	}
	out.pc = V / I
	out.pcExact = math.NaN()
	if exact != nil {
		V, I = stokes(xth, xph)
		out.pcExact = V / I
	}
	return
}

func ledgerOf(f *em.FDTD) (res, captured float64) {
	l := f.Ledger()
	in := l.W0 + l.SourceWork
	return math.Abs(l.Residual()) / in, l.Absorbed / in
}

// hopfion runs knotted_light's knot through the probe sphere.
func (s HelicitySetup) hopfion() (HelicityRun, error) {
	f, _, err := s.grid()
	if err != nil {
		return HelicityRun{}, err
	}
	h := Hopfion{}
	sign := h.Sign()
	c := em.SpeedOfLight
	L0 := s.HopfL0 * s.Lambda0
	f.SetEFromPotential(func(x, y, z float64) [3]float64 {
		p := h.P(x/L0, y/L0, z/L0, 0)
		return [3]float64{L0 * real(p[0]), L0 * real(p[1]), L0 * real(p[2])}
	})
	tb := c * f.Dt / 2 / L0
	f.SetBFromPotential(func(x, y, z float64) [3]float64 {
		p := h.P(x/L0, y/L0, z/L0, tb)
		k := sign * L0 / c
		return [3]float64{k * imag(p[0]), k * imag(p[1]), k * imag(p[2])}
	})
	exact := func(x, y, z, t float64) [3]float64 {
		F := h.F(x/L0, y/L0, z/L0, c*t/L0)
		return [3]float64{real(F[0]), real(F[1]), real(F[2])}
	}
	tEnd := (s.Probe*s.Lambda0 + 8*L0) / c
	r := HelicityRun{Name: "hopfion (knotted_light)"}
	r.take(recordFull(f, s.probes(), tEnd, exact, []float64{0, 2.5 * s.Lambda0 / c, 5 * s.Lambda0 / c}, s.Absorb))
	r.Residual, r.Captured = ledgerOf(f)
	return r, nil
}

func (s HelicitySetup) envelope() (func(t float64) float64, float64) {
	T0 := s.Lambda0 / em.SpeedOfLight
	sig := s.Envelope * T0
	tc := 4 * sig
	return func(t float64) float64 { d := (t - tc) / sig; return math.Exp(-0.5 * d * d) }, tc + 4*sig
}

// dipoles runs the chiral pair (chiral = true) or the rotating dipole.
func (s HelicitySetup) dipoles(chiral bool) (HelicityRun, error) {
	f, n, err := s.grid()
	if err != nil {
		return HelicityRun{}, err
	}
	w := 2 * math.Pi * em.SpeedOfLight / s.Lambda0
	env, tOff := s.envelope()
	m := n / 2
	sin := func(a float64) func(float64) float64 {
		return func(t float64) float64 { return a * env(t) * math.Sin(w*t) }
	}
	cos := func(a float64) func(float64) float64 {
		return func(t float64) float64 { return a * env(t) * math.Cos(w*t) }
	}
	var src []em.EdgeSource
	name := "rotating electric dipole"
	if chiral {
		name = "chiral dipole pair (p ∥ m, quadrature)"
		// p: z edge at the centre; m: a 2 × 2-cell square loop in the
		// z = m plane, |m| = c|p| ⇒ I_loop = I_p·D/(k·A) = I_p/(4kD).
		k := w / em.SpeedOfLight
		D := s.Lambda0 / float64(s.Cells)
		il := 1 / (4 * k * D)
		src = append(src, em.EdgeSource{Comp: 2, I: m, J: m, K: m, Current: sin(1)})
		for _, i := range []int{m - 1, m} {
			src = append(src, em.EdgeSource{Comp: 0, I: i, J: m - 1, K: m, Current: sin(il)})
			src = append(src, em.EdgeSource{Comp: 0, I: i, J: m + 1, K: m, Current: sin(-il)})
		}
		for _, j := range []int{m - 1, m} {
			src = append(src, em.EdgeSource{Comp: 1, I: m + 1, J: j, K: m, Current: sin(il)})
			src = append(src, em.EdgeSource{Comp: 1, I: m - 1, J: j, K: m, Current: sin(-il)})
		}
	} else {
		src = append(src, em.EdgeSource{Comp: 0, I: m, J: m, K: m, Current: cos(1)},
			em.EdgeSource{Comp: 1, I: m, J: m, K: m, Current: sin(1)})
	}
	for _, e := range src {
		if err := f.AddSource(e); err != nil {
			return HelicityRun{}, err
		}
	}
	tEnd := tOff + (s.Probe+1)*s.Lambda0/em.SpeedOfLight
	r := HelicityRun{Name: name}
	r.take(recordFull(f, s.probes(), tEnd, nil, []float64{tOff - s.Envelope*4*s.Lambda0/em.SpeedOfLight + 2.5*s.Lambda0/em.SpeedOfLight}, s.Absorb))
	r.Residual, r.Captured = ledgerOf(f)
	return r, nil
}

// thermal runs one realization of the random core.
func (s HelicitySetup) thermal(seed uint64) (HelicityRun, error) {
	f, n, err := s.grid()
	if err != nil {
		return HelicityRun{}, err
	}
	rng := rand.New(rand.NewPCG(seed, 0x50b0b1e))
	w0 := 2 * math.Pi * em.SpeedOfLight / s.Lambda0
	env, tOff := s.envelope()
	m := n / 2
	rc := s.Core * float64(s.Cells)
	for added := 0; added < s.Dipoles; {
		i := m + int(math.Round((2*rng.Float64()-1)*rc))
		j := m + int(math.Round((2*rng.Float64()-1)*rc))
		k := m + int(math.Round((2*rng.Float64()-1)*rc))
		di, dj, dk := float64(i-m), float64(j-m), float64(k-m)
		if di*di+dj*dj+dk*dk > rc*rc {
			continue
		}
		w := w0 * (0.5 + rng.Float64())
		ph := 2 * math.Pi * rng.Float64()
		if err := f.AddSource(em.EdgeSource{Comp: rng.IntN(3), I: i, J: j, K: k,
			Current: func(t float64) float64 { return env(t) * math.Cos(w*t+ph) }}); err != nil {
			return HelicityRun{}, err
		}
		added++
	}
	tEnd := tOff + (s.Probe+s.Core+1)*s.Lambda0/em.SpeedOfLight
	lc := s.Lambda0 / em.SpeedOfLight
	tc := tOff - s.Envelope*4*lc
	r := HelicityRun{Name: fmt.Sprintf("thermal core, seed %d", seed)}
	r.take(recordFull(f, s.probes(), tEnd, nil, []float64{tc, tc + 2.5*lc, tc + 5*lc}, s.Absorb))
	r.Residual, r.Captured = ledgerOf(f)
	return r, nil
}

// RunAll makes every measurement.
func (s HelicitySetup) RunAll() (HelicityResult, error) {
	var res HelicityResult
	r, err := s.hopfion()
	if err != nil {
		return res, err
	}
	res.Runs = append(res.Runs, r)
	for _, chiral := range []bool{true, false} {
		r, err := s.dipoles(chiral)
		if err != nil {
			return res, err
		}
		res.Runs = append(res.Runs, r)
	}
	var sum, sq float64
	for k := 1; k <= s.Seeds; k++ {
		r, err := s.thermal(uint64(k))
		if err != nil {
			return res, err
		}
		res.Runs = append(res.Runs, r)
		sum += r.Pc
		sq += r.Pc * r.Pc
	}
	res.ThermalMean = sum / float64(s.Seeds)
	res.ThermalRMS = math.Sqrt(sq / float64(s.Seeds))
	res.ThermalExpected = 1 / math.Sqrt(2*float64(s.Dipoles))
	return res, nil
}

func printHelicity(h HelicityResult) {
	for _, r := range h.Runs {
		ex := ""
		if !math.IsNaN(r.PcExact) {
			ex = fmt.Sprintf(" (exact field at the probes %+.4f)", r.PcExact)
		}
		fmt.Printf("  %-42s P_c = %+.4f%s; ledger residual %.1e, %.3f of the energy absorbed\n",
			r.Name, r.Pc, ex, r.Residual, r.Captured)
	}
	fmt.Printf("  thermal core: mean P_c %+.4f, rms %.4f over the seeds (independent radiators: ~1/√(2N) = %.4f)\n",
		h.ThermalMean, h.ThermalRMS, h.ThermalExpected)
}
