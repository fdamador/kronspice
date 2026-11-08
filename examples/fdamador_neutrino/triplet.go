// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 8: the user's hypothesis H41 (CLAUDE.md §1, 2026-11-07): the three
// neutrino flavours are three spin states of the ring at detection —
// turning one way with full magnetic field, not turning (no magnetic
// field), turning the other way with full magnetic field — and, if the
// whole unit also turns because of outside conditions, the Coriolis and
// centrifugal terms add energy that could account for some of the mass.
//
// Testable form: the three flavours are the azimuthal states m = +1, 0, −1
// of the ring (a spin-1-like triplet), labelled by the angular momentum
// J_z about the axis, the quantity a free ring conserves; the
// instantaneous B passes through zero twice a period in every standing
// pattern, so it cannot label a state.
//
// On the CL ring line of Parts 6 and 7 (z_j = u_j + iωL_g·g_j, as in
// flavour.go) the triplet is
//
//	m = +1: z_j ∝ e^{−iκ₁j} (CW, turning toward increasing j),
//	m =  0: z_j ∝ 1        (the uniform state: every shunt inductor in
//	                        phase, no voltage across the series C),
//	m = −1: z_j ∝ e^{+iκ₁j} (CCW),
//
// with κ₁ = 2π/N. In the field solver (tripletfield.go) m = 0 is the
// E-ring (TE01δ) and m = ±1 are Part 7b's turning patterns.
//
// The angular momentum of a field made of harmonics of order m with
// actions I_m = U_m/ω_m is J = Σ m·I_m: the generator of rotations, so
// that a rotation by α changes the phase of harmonic m by m·α. A defect
// that breaks the rotation symmetry exerts a torque on the field,
// dJ/dt = −∂H/∂α (H the energy with the field turned by α against the
// defect); the defect takes up −dJ/dt. Landau & Lifshitz, Mechanics, §9
// and §39; Kron's torque ½iᵀ(∂L/∂θ)i (Equivalent Circuits of Electric
// Machinery, 1951, Epilogue) is the same generalized force.

import (
	"math"
	"math/cmplx"
	"runtime"
	"sort"
	"sync"
)

// tripletVecs returns the unit z-patterns of m = +1, 0, −1 at κ₁ = 2π/N.
func (r FRing) tripletVecs() [3][]complex128 {
	n := r.N
	k := r.Kappa(1)
	var v [3][]complex128
	for a := range v {
		v[a] = make([]complex128, n)
	}
	s := complex(1/math.Sqrt(float64(n)), 0)
	for j := 0; j < n; j++ {
		v[0][j] = s * cmplx.Exp(complex(0, -k*float64(j)))
		v[1][j] = s
		v[2][j] = s * cmplx.Exp(complex(0, k*float64(j)))
	}
	return v
}

// typeVecs returns Part 6's C, S, D at κ₁ (reference c = ½) as unit z-patterns.
func (r FRing) typeVecs() [3][]complex128 {
	var v [3][]complex128
	k := r.Kappa(1)
	for a, name := range []string{"C", "S", "D"} {
		c := r.Type(name, k, 0.5)
		v[a] = make([]complex128, r.N)
		var s float64
		for j := 0; j < r.N; j++ {
			v[a][j] = c.A*cmplx.Exp(complex(0, -k*float64(j))) + c.B*cmplx.Exp(complex(0, k*float64(j)))
			s += sq(cmplx.Abs(v[a][j]))
		}
		for j := range v[a] {
			v[a][j] /= complex(math.Sqrt(s), 0)
		}
	}
	return v
}

// gram returns ⟨v_a|v_b⟩.
func gram(vs [3][]complex128) [3][3]complex128 {
	var g [3][3]complex128
	for a := 0; a < 3; a++ {
		for b := 0; b < 3; b++ {
			for j := range vs[a] {
				g[a][b] += cmplx.Conj(vs[a][j]) * vs[b][j]
			}
		}
	}
	return g
}

// hermEig3 returns the eigenvalues (ascending) of a 3×3 Hermitian matrix,
// from its real symmetric embedding [[Re, −Im], [Im, Re]] (each eigenvalue
// appears twice there).
func hermEig3(h [3][3]complex128) [3]float64 {
	m := make([][]float64, 6)
	for i := range m {
		m[i] = make([]float64, 6)
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			m[i][j], m[i+3][j+3] = real(h[i][j]), real(h[i][j])
			m[i][j+3], m[i+3][j] = -imag(h[i][j]), imag(h[i][j])
		}
	}
	vals, _ := eigSym(m)
	sort.Float64s(vals)
	return [3]float64{vals[0], vals[2], vals[4]}
}

// TripletIndep is Part 8a on the ring line.
type TripletIndep struct {
	GramTriplet, GramTypes [3]float64 // Gram eigenvalues
	MaxOffTriplet          float64    // largest |⟨a|b⟩|, a ≠ b
	TypeCS, TypeCD         float64    // |⟨C|S⟩|, |⟨C|D⟩|
	Omega                  [3]float64 // exact ω of m = +1, 0, −1 on the uniform ring [rad/s]
	Persist                [3]Persist8
}

// Persist8 is one triplet state left alone on the perfect ring (kronspice).
type Persist8 struct {
	Name      string
	MinOwn    float64 // smallest share of Σ|z|² in the state's own pattern
	Periods   float64
	LedgerRel float64
	Jw0, JwT  float64 // J·ω/U at the start and the end (+1, 0, −1 expected)
}

// share is |⟨f|z⟩|²/Σ|z|² at sample t with z = u + iωL_g·g.
func (r FRing) share(fr FRun, t int, f []complex128, omega float64) float64 {
	var p complex128
	var tot float64
	for j := 0; j < r.N; j++ {
		z := complex(fr.U[t][j], omega*r.Lg*fr.G[t][j])
		p += cmplx.Conj(f[j]) * z
		tot += sq(cmplx.Abs(z))
	}
	return sq(cmplx.Abs(p)) / tot
}

// ringJ is J·ω₁/U of the harmonic-1 and uniform content at sample t: the
// energies U_CW = N·c_κ|A|²/4, U_CCW = N·c_κ|B|²/4 (and U_0 = N·Cg|A₀|²/4
// for the uniform state) give J = (U_CW − U_CCW)/ω₁; U is their sum.
func (r FRing) ringJ(fr FRun, t int) (j, u float64) {
	k := r.Kappa(1)
	w1, w0 := r.OmegaK(k), r.OmegaK(0)
	ck := 4*r.C*sq(math.Sin(k/2)) + r.Cg
	A, B := r.Amp(fr, t, k, w1)
	A0, _ := r.Amp(fr, t, 0, w0)
	n := float64(r.N)
	ucw, uccw := n*ck*sq(cmplx.Abs(A))/4, n*ck*sq(cmplx.Abs(B))/4
	u0 := n * r.Cg * sq(cmplx.Abs(A0)) / 4
	return (ucw - uccw) / w1, ucw + uccw + u0
}

// TripletIndependence runs Part 8a on the plain ring line.
func TripletIndependence(periods float64) (TripletIndep, error) {
	r := DefaultFRing(false)
	var ti TripletIndep
	tv, yv := r.tripletVecs(), r.typeVecs()
	gt, gy := gram(tv), gram(yv)
	ti.GramTriplet, ti.GramTypes = hermEig3(gt), hermEig3(gy)
	for a := 0; a < 3; a++ {
		for b := 0; b < 3; b++ {
			if a != b {
				ti.MaxOffTriplet = math.Max(ti.MaxOffTriplet, cmplx.Abs(gt[a][b]))
			}
		}
	}
	ti.TypeCS, ti.TypeCD = cmplx.Abs(gy[0][1]), cmplx.Abs(gy[0][2])
	k := r.Kappa(1)
	ti.Omega = [3]float64{r.OmegaK(k), r.OmegaK(0), r.OmegaK(k)}
	starts := [3]Comp{r.Type("C", k, 0.5), {Kappa: 0, Omega: r.OmegaK(0), A: 1}, r.Type("CCW", k, 0.5)}
	names := [3]string{"m = +1 (CW)", "m = 0 (uniform)", "m = −1 (CCW)"}
	for a := 0; a < 3; a++ {
		w := ti.Omega[a]
		fr, err := r.Run([]Comp{starts[a]}, w, periods*2*math.Pi/w, 60)
		if err != nil {
			return ti, err
		}
		p := Persist8{Name: names[a], MinOwn: 1, Periods: periods, LedgerRel: fr.LedgerRel}
		for t := range fr.Time {
			p.MinOwn = math.Min(p.MinOwn, r.share(fr, t, tv[a], w))
		}
		j0, u0 := r.ringJ(fr, 0)
		jt, ut := r.ringJ(fr, len(fr.Time)-1)
		w1 := ti.Omega[0]
		p.Jw0, p.JwT = j0*w1/u0, jt*w1/ut
		ti.Persist[a] = p
	}
	return ti, nil
}

// ---------------------------------------------------------------- 8b mixing

// TripletKnobs are the ring's tuning parameters for the triplet cluster.
// The shunt (node-to-ground) capacitance varies smoothly round the ring,
//
//	ΔC_j = C·(M1·cos κ₁j + M2·cos(2κ₁j + 2ψ)),
//
// the m = 1 deformation (an off-centre ring) coupling m = 0 to m = ±1 and
// the m = 2 deformation (an elliptical ring, turned by ψ) coupling +1 to
// −1; Cg is the stray capacitance of every node (it maps the capacitance
// spectrum onto ω), Nu the rotation on the κ₁ pair (Part 7's gyrators).
// D on section S is a series-capacitor defect as in Part 6, and N1, N2
// point shunt defects on nodes 0 and P2 (used for the series-only check).
type TripletKnobs struct {
	M1, M2, Psi float64
	Cg          float64
	Nu          float64
	D           float64
	S           int
	N1, N2      float64
	P2          int
}

// ring builds the FRing for the knobs.
func (k TripletKnobs) ring() FRing {
	r := DefaultFRing(false)
	r.Cg = k.Cg
	r.Def = map[int]float64{k.S: k.D}
	r.NodeC = map[int]float64{}
	k1 := r.Kappa(1)
	for j := 0; j < r.N; j++ {
		x := k1 * float64(j)
		r.NodeC[j] = r.C * (k.M1*math.Cos(x) + k.M2*math.Cos(2*x+2*k.Psi))
	}
	r.NodeC[0] += k.N1 * r.C
	r.NodeC[k.P2] += k.N2 * r.C
	r.Gyro = r.GyroFor(k.Nu, k1)
	return r
}

// valid reports whether every node keeps a positive capacitance to ground.
func (k TripletKnobs) valid() bool {
	r := k.ring()
	for j := 0; j < r.N; j++ {
		if r.Cg+r.NodeC[j] <= 0 {
			return false
		}
	}
	return true
}

// TripletCl is a triplet cluster: the three exact modes with the most
// m = +1, 0, −1 content, the mixing in Part 7's conventions (RotCluster:
// angles from |U|, J and δ_CP from the Löwdin-unitarized U), and the
// splitting ratio in ω² (mass²) beside Part 6's ratio in ω.
type TripletCl struct {
	RotCluster
	Knobs   TripletKnobs
	Ratio2  float64 // |ω3² − ω1²| / |ω2² − ω1²|
	Cost2   float64 // angles + ratio in ω² (no δ_CP)
	CostCP2 float64 // with δ_CP against 212° (+26/−41)
}

// TripletMix computes the cluster for the knobs, flavours e, μ, τ =
// triplet index perm (0 = +1, 1 = 0, 2 = −1).
func TripletMix(k TripletKnobs, perm [3]int) TripletCl {
	r := k.ring()
	n := r.N
	var us [][]complex128
	var om []float64
	if r.Gyro == 0 {
		w, vec := r.Modes()
		for i := range w {
			u := make([]complex128, n)
			for j := 0; j < n; j++ {
				u[j] = complex(vec[j][i], 0)
			}
			us = append(us, u)
		}
		om = w
	} else {
		rm := r.RotModes()
		us, om = rm.U, rm.Omega
	}
	f := r.tripletVecs()
	proj := func(a, i int) complex128 {
		var d complex128
		for j := 0; j < n; j++ {
			d += cmplx.Conj(f[a][j]) * us[i][j]
		}
		return d
	}
	type cand struct {
		i int
		w float64
	}
	var cs []cand
	for i := range us {
		var s float64
		for a := 0; a < 3; a++ {
			s += sq(cmplx.Abs(proj(a, i)))
		}
		cs = append(cs, cand{i, s})
	}
	sort.SliceStable(cs, func(a, b int) bool { return cs[a].w > cs[b].w })
	ids := []int{cs[0].i, cs[1].i, cs[2].i}
	sort.Slice(ids, func(a, b int) bool { return om[ids[a]] < om[ids[b]] })
	tc := TripletCl{Knobs: k}
	tc.D1, tc.P1, tc.Nu, tc.Gyro = k.D, k.S, k.Nu, r.Gyro
	for m, i := range ids {
		tc.Omega[m] = om[i]
		var s float64
		for a := 0; a < 3; a++ {
			tc.U[a][m] = proj(a, i)
			s += sq(cmplx.Abs(tc.U[a][m]))
		}
		tc.Leak = math.Max(tc.Leak, 1-s)
	}
	tc.angles(perm, NuFIT6())
	w := tc.Omega
	i1, i2, i3 := tc.Mass[0], tc.Mass[1], tc.Mass[2]
	tc.Ratio2 = math.Abs(w[i3]*w[i3]-w[i1]*w[i1]) / math.Abs(w[i2]*w[i2]-w[i1]*w[i1])
	t := NuFIT6()
	tc.Cost2 = sq((tc.Th12-t.Th12)/t.S12) + sq((tc.Th23-t.Th23)/t.S23) + sq((tc.Th13-t.Th13)/t.S13) + sq(math.Log(tc.Ratio2/t.Ratio)/t.SRatio)
	dd := tc.DeltaCP - t.DeltaCP
	sig := 26.0
	if dd < 0 {
		sig = 41
	}
	tc.CostCP2 = tc.Cost2 + sq(dd/sig)
	return tc
}

// tripletPerms are the six assignments of m = +1, 0, −1 to e, μ, τ.
var tripletPerms = [][3]int{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}

// TripletFit is the best cluster found for one setting of the search.
type TripletFit struct {
	Best      TripletCl
	PerAssign [6]float64 // best objective (cost + (leak/leakScale)²) for each assignment (tripletPerms order)
	Evals     int
}

// FitTriplet searches the knobs M1, M2, ψ and log Cg for the measured
// angles and the splitting ratio in ω², for each assignment of
// m = +1, 0, −1 to e, μ, τ, by Nelder–Mead from a grid of starts (the
// pair m = 0 and one standing ±1 combination is near-degenerate when
// M2 ≈ 2·4 sin²(κ₁/2) ≈ 0.14, so M2 is started finely there). With rot
// (from), the rotation ν is added, δ_CP enters the cost, and the search
// starts from the reciprocal result `from` with a set of ν. In parallel,
// results kept by index (deterministic).
func FitTriplet(rot bool, leakScale float64, from ...TripletCl) TripletFit {
	type job struct {
		perm  int
		start []float64
	}
	var jobs []job
	if rot && len(from) > 0 {
		k := from[0].Knobs
		pi := 0
		for q, pm := range tripletPerms {
			if pm == from[0].Assign {
				pi = q
			}
		}
		for _, nu := range []float64{-3, -1, -0.3, -0.1, 0.1, 0.3, 1, 3} {
			jobs = append(jobs, job{pi, []float64{k.M1, k.M2, k.Psi, math.Log10(k.Cg / 1e-9), nu}})
		}
	} else {
		for pi := range tripletPerms {
			for i := 0; i <= 20; i++ {
				for _, psi := range []float64{0.04, 0.07, 0.1, 0.14, 0.2} {
					for _, m1 := range []float64{0.02, 0.05} {
						for _, lc := range []float64{-1, -0.3} {
							jobs = append(jobs, job{pi, []float64{m1, 0.10 + 0.004*float64(i), psi, lc}})
						}
					}
				}
			}
		}
	}
	res := make([]TripletCl, len(jobs))
	evals := make([]int, len(jobs))
	sem := make(chan struct{}, runtime.NumCPU())
	var wg sync.WaitGroup
	// objective: the cost plus (leak/leakScale)², which can keep the three
	// mass states inside the triplet (large deformations pull in m = ±2)
	// (leakScale < 0: a hard cap, leak ≤ |leakScale|)
	cost := func(c TripletCl) float64 {
		base := c.Cost2
		if rot {
			base = c.CostCP2
		}
		if leakScale < 0 {
			if c.Leak > -leakScale {
				return math.Inf(1)
			}
			return base
		}
		return base + sq(c.Leak/leakScale)
	}
	knobs := func(x []float64) (TripletKnobs, bool) {
		k := TripletKnobs{M1: x[0], M2: x[1], Psi: x[2], Cg: 1e-9 * math.Pow(10, x[3])}
		if rot {
			k.Nu = x[4] * 1e-3
		}
		ok := math.Abs(k.M1) <= 1 && math.Abs(k.M2) <= 1 && k.Cg >= 1e-13 && k.Cg <= 1e-8 && math.Abs(k.Nu) <= 0.05 && k.valid()
		return k, ok
	}
	for q, jb := range jobs {
		wg.Add(1)
		sem <- struct{}{}
		go func(q int, jb job) {
			defer wg.Done()
			defer func() { <-sem }()
			perm := tripletPerms[jb.perm]
			ev := 0
			eval := func(x []float64) float64 {
				ev++
				k, ok := knobs(x)
				if !ok {
					return math.Inf(1)
				}
				return cost(TripletMix(k, perm))
			}
			step := []float64{0.2 * jb.start[0], 0.01, 0.02, 0.1}
			if rot {
				step = append(step, 0.3*math.Abs(jb.start[4]))
			}
			x, _ := nelderMead(eval, jb.start, step, 300)
			x, _ = nelderMead(eval, x, step, 300) // restart once from the end point
			best := TripletCl{Cost2: math.Inf(1), CostCP2: math.Inf(1)}
			if k, ok := knobs(x); ok {
				best = TripletMix(k, perm)
			}
			res[q], evals[q] = best, ev
		}(q, jb)
	}
	wg.Wait()
	tf := TripletFit{Best: TripletCl{Cost2: math.Inf(1), CostCP2: math.Inf(1)}}
	for pi := range tf.PerAssign {
		tf.PerAssign[pi] = math.Inf(1)
	}
	for q, c := range res {
		tf.Evals += evals[q]
		pi := jobs[q].perm
		tf.PerAssign[pi] = math.Min(tf.PerAssign[pi], cost(c))
		if cost(c) < cost(tf.Best) {
			tf.Best = c
		}
	}
	return tf
}

// TripletRun is the kronspice run of a reciprocal triplet cluster started
// in the e flavour, against the exact mode expansion (same readout) and
// the 3×3 formula.
type TripletRun struct {
	T              []float64    // seconds, thinned
	P              [3][]float64 // e, μ, τ readout / total, kronspice
	PFull          [3][]float64 // the same from all N exact modes
	P3             [3][]float64 // 3×3 formula with the unitarized U
	MaxDevAll      float64
	MaxDev3        float64
	LedgerRel      float64
	Periods        float64
	MinM0, MaxM0   float64 // the m = 0 share's range (how far the start leaves its state)
	ExchangeJ      float64 // largest |J|·ω/U reached (angular momentum the defects hand to the field)
	StartFlavourM  int     // triplet index of the start state
	AssignTriplets [3]int
}

// RunTriplet runs the cluster tc (reciprocal, Nu = 0) in kronspice from
// its e flavour for `slow` periods of the slowest beat.
func RunTriplet(tc TripletCl, slow float64) (TripletRun, error) {
	k := tc.Knobs
	r := k.ring()
	n := r.N
	f := r.tripletVecs()
	from := tc.Assign[0]
	kap := r.Kappa(1)
	// readout frequency of each flavour: the uniform ring's (with this Cg)
	wf := [3]float64{r.OmegaK(kap), r.OmegaK(0), r.OmegaK(kap)}
	// start: the e flavour's pattern with g = 0 (u = Re z at t = 0)
	var init Comp
	switch from {
	case 0:
		init = Comp{Kappa: kap, Omega: wf[0], A: 1}
	case 1:
		init = Comp{Kappa: 0, Omega: wf[1], A: 1}
	default:
		init = Comp{Kappa: kap, Omega: wf[2], B: 1}
	}
	// a start with g = 0 (u = Re z at t = 0); for m = 0 it is already so,
	// for ±1 it is the standing cos pattern, half of each sense
	if from != 1 {
		init.A, init.B = 0.5*init.A+0.5*init.B, 0.5*init.A+0.5*init.B
	}
	u0, _ := r.state([]Comp{init})
	dw := math.Inf(1)
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			dw = math.Min(dw, math.Abs(tc.Omega[i]-tc.Omega[j]))
		}
	}
	wref := (tc.Omega[0] + tc.Omega[1] + tc.Omega[2]) / 3
	stop := slow * 2 * math.Pi / dw
	const per = 120
	fr, err := r.Run([]Comp{init}, wref, stop, per)
	if err != nil {
		return TripletRun{}, err
	}
	tr := TripletRun{LedgerRel: fr.LedgerRel, Periods: stop * wref / (2 * math.Pi), StartFlavourM: from, AssignTriplets: tc.Assign,
		MinM0: math.Inf(1), MaxM0: math.Inf(-1)}
	read := func(u, g []float64) [3]float64 {
		var out [3]float64
		var tot float64
		for j := 0; j < n; j++ {
			tot += sq(u[j]) + sq(wref*r.Lg*g[j])
		}
		for a := 0; a < 3; a++ {
			var p complex128
			for j := 0; j < n; j++ {
				p += cmplx.Conj(f[a][j]) * complex(u[j], wf[a]*r.Lg*g[j])
			}
			out[a] = sq(cmplx.Abs(p)) / tot
		}
		return out
	}
	// exact expansion at the trapezoidal frequencies (as ThreeFlavour)
	om, vec := r.Modes()
	h := 2 * math.Pi / wref / per
	for i := range om {
		om[i] = 2 / h * math.Atan(om[i]*h/2)
	}
	c := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			c[i] += vec[j][i] * u0[j]
		}
	}
	U := tc.UU
	ue, ug := make([]float64, n), make([]float64, n)
	for t := range fr.Time {
		tt := fr.Time[t]
		pm := read(fr.U[t], fr.G[t])
		for j := 0; j < n; j++ {
			var u, g float64
			for i := 0; i < n; i++ {
				u += c[i] * vec[j][i] * math.Cos(om[i]*tt)
				g += c[i] * vec[j][i] * math.Sin(om[i]*tt) / (om[i] * r.Lg)
			}
			ue[j], ug[j] = u, g
		}
		pf := read(ue, ug)
		var p3 [3]float64
		for b := 0; b < 3; b++ {
			var amp complex128
			for i := 0; i < 3; i++ {
				amp += U[0][i] * cmplx.Conj(U[b][i]) * cmplx.Exp(complex(0, tc.Omega[tc.Mass[i]]*tt))
			}
			p3[b] = sq(cmplx.Abs(amp))
		}
		var pe [3]float64 // by e, μ, τ
		var pfe [3]float64
		for b := 0; b < 3; b++ {
			pe[b], pfe[b] = pm[tc.Assign[b]], pf[tc.Assign[b]]
			tr.MaxDevAll = math.Max(tr.MaxDevAll, math.Abs(pe[b]-pfe[b]))
			tr.MaxDev3 = math.Max(tr.MaxDev3, math.Abs(pe[b]-p3[b]))
		}
		tr.MinM0, tr.MaxM0 = math.Min(tr.MinM0, pm[1]), math.Max(tr.MaxM0, pm[1])
		jj, uu := r.ringJ(fr, t)
		tr.ExchangeJ = math.Max(tr.ExchangeJ, math.Abs(jj*wf[0]/uu))
		if t%10 == 0 {
			tr.T = append(tr.T, tt)
			for b := 0; b < 3; b++ {
				tr.P[b] = append(tr.P[b], pe[b])
				tr.PFull[b] = append(tr.PFull[b], pfe[b])
				tr.P3[b] = append(tr.P3[b], p3[b])
			}
		}
	}
	return tr, nil
}

// SymmetryCheck reports how far a reciprocal cluster is from μ–τ
// reflection symmetry when e is the m = 0 state: max_i ||U_μi| − |U_τi||.
func (tc TripletCl) MuTauAsym() float64 {
	var m float64
	for i := 0; i < 3; i++ {
		m = math.Max(m, math.Abs(cmplx.Abs(tc.UU[1][i])-cmplx.Abs(tc.UU[2][i])))
	}
	return m
}

// ---------------------------------------------------------------- 8c torque

// TorqueRun is +1 → −1 under one capacitor defect, with the angular
// momentum of the field and the torque the defect exerts.
type TorqueRun struct {
	Delta      float64
	T          []float64 // seconds, thinned
	Jw         []float64 // J·ω₁/U₀ of the field (harmonic 1 and the uniform state)
	TorqueInt  []float64 // −∫(∂H/∂α) dt · ω₁/U₀: the angular momentum the defect has given the field
	U          []float64 // field energy (harmonic 1 + uniform) / U₀
	J0, JMin   float64   // J·ω₁/U₀ at the start and its minimum
	DeltaJ     float64   // J(t_min) − J(0), in U₀/ω₁
	TorqueAt   float64   // ∫ −∂H/∂α dt to t_min, in U₀/ω₁
	MaxDiff    float64   // largest |ΔJ(t) − ∫τ dt| over the run, in U₀/ω₁
	Df         float64   // beat frequency [Hz]
	LedgerRel  float64
	EnergyLoss float64 // 1 − U(t_min)/U₀ (the 1 mΩ shunt loss)
	DefectWork float64 // ∫ τ·Ω_defect dt = 0: the defect does not move, so it takes angular momentum and no energy
}

// Torque runs the +1 (CW) state on the plain ring with defect δ on
// section 0 for `beats` beat periods. The torque on the field is
// −∂H/∂α: H = Σ q_j²/c_j with the section charges q_j = c_j(u_j − u_{j+1})/2
// per wire turned by α against the defect (charges trigonometrically
// interpolated round the ring), so only the defect section contributes:
// ∂H/∂x = −2·q_s·q′(s)·(1/c_s − 1/C) per section of turn, α = 2πx/N.
func Torque(delta, beats float64) (TorqueRun, error) {
	r := DefaultFRing(false)
	r.Def = map[int]float64{0: delta}
	n := r.N
	k := r.Kappa(1)
	w1 := r.OmegaK(k)
	lo, hi := r.pairNear(w1)
	df := (hi - lo) / (2 * math.Pi)
	const per = 120
	stop := beats / df
	fr, err := r.Run([]Comp{r.Type("C", k, 0.5)}, w1, stop, per)
	if err != nil {
		return TorqueRun{}, err
	}
	tq := TorqueRun{Delta: delta, Df: df, LedgerRel: fr.LedgerRel}
	j0, u0 := r.ringJ(fr, 0)
	tq.J0 = j0 * w1 / u0
	cs := r.C * (1 + delta)
	dHdA := func(t int) float64 {
		q := make([]float64, n)
		for j := 0; j < n; j++ {
			c := r.C * (1 + r.Def[j])
			q[j] = c * (fr.U[t][j] - fr.U[t][(j+1)%n]) / 2
		}
		// q′ at section 0 from the trigonometric interpolant (Nyquist term dropped)
		var dq float64
		for m := -n/2 + 1; m < n/2; m++ {
			var qm complex128
			for j := 0; j < n; j++ {
				qm += complex(q[j], 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(m*j)/float64(n)))
			}
			qm /= complex(float64(n), 0)
			dq += real(complex(0, 2*math.Pi*float64(m)/float64(n)) * qm) // at s = 0: e^{0} = 1
		}
		dHdx := -2 * q[0] * dq * (1/cs - 1/r.C)
		return dHdx * float64(n) / (2 * math.Pi)
	}
	var integ float64
	prev := dHdA(0)
	tq.JMin = tq.J0
	imin := 0
	var jw, ti, uu []float64
	for t := range fr.Time {
		cur := dHdA(t)
		if t > 0 {
			integ += -0.5 * (prev + cur) * (fr.Time[t] - fr.Time[t-1])
		}
		prev = cur
		jj, ut := r.ringJ(fr, t)
		jw = append(jw, jj*w1/u0)
		ti = append(ti, integ*w1/u0)
		uu = append(uu, ut/u0)
		if jw[t] < tq.JMin {
			tq.JMin, imin = jw[t], t
		}
		tq.MaxDiff = math.Max(tq.MaxDiff, math.Abs((jw[t]-tq.J0)-ti[t]))
	}
	tq.DeltaJ = jw[imin] - tq.J0
	tq.TorqueAt = ti[imin]
	tq.EnergyLoss = 1 - uu[imin]
	for t := 0; t < len(jw); t += 20 {
		tq.T, tq.Jw, tq.TorqueInt, tq.U = append(tq.T, fr.Time[t]), append(tq.Jw, jw[t]), append(tq.TorqueInt, ti[t]), append(tq.U, uu[t])
	}
	return tq, nil
}

// FreeRingJ shows that on the ring without a defect the +1 state keeps
// its J: the largest |J·ω/U − 1| over `periods` periods.
func FreeRingJ(periods float64) (float64, error) {
	r := DefaultFRing(false)
	k := r.Kappa(1)
	w1 := r.OmegaK(k)
	fr, err := r.Run([]Comp{r.Type("C", k, 0.5)}, w1, periods*2*math.Pi/w1, 60)
	if err != nil {
		return 0, err
	}
	var m float64
	for t := range fr.Time {
		j, u := r.ringJ(fr, t)
		m = math.Max(m, math.Abs(j*w1/u-1))
	}
	return m, nil
}

// ---------------------------------------------------------------- 8e outside rotation

// RotMass is the lab-frame energy a whole-unit rotation Ω adds to one
// quantum of a ring state of order m and frequency ω (in the co-rotating
// frame), radius a (⟨ρ²⟩ = rho2·a²): Coriolis (first order, exact form
// E = E′ + Ω·J, Landau & Lifshitz, Mechanics, §39) and centrifugal (second
// order, ½IΩ² with the field's moment of inertia I = (ħω/c²)⟨ρ²⟩).
type RotMass struct {
	Name           string
	A, Omega, Rho2 float64
	M              int
	// Ω at which each contribution equals the target energy (0.1 eV declared):
	OmCor, OmCen   float64 // [rad/s] (OmCor = +Inf for m = 0)
	RimCor, RimCen float64 // aΩ [m/s]
	CorOverOmega   float64 // Ω_cor/ω (the first-order form needs ≪ 1)
	// at the burst speed of T1000G (Ω_b = v_b/a):
	OmBurst            float64
	CorAtBurst         float64 // eV per quantum
	CenAtBurst         float64 // eV per quantum
	QuantumEV          float64 // ħω [eV]
	CenCeilingOverHbar float64 // ½(v/c)² at v = c: the centrifugal term can never exceed ½ħω·(v/c)² < ½ħω
}

// OutsideRotation evaluates RotMass for a target energy (eV).
func OutsideRotation(name string, m int, omega, a, rho2, targetEV, burst float64) RotMass {
	rm := RotMass{Name: name, A: a, Omega: omega, Rho2: rho2, M: m}
	e := targetEV * E
	q := Hbar * omega
	rm.QuantumEV = q / E
	if m != 0 {
		rm.OmCor = e / (Hbar * math.Abs(float64(m)))
	} else {
		rm.OmCor = math.Inf(1)
	}
	rm.RimCor = rm.OmCor * a
	rm.CorOverOmega = rm.OmCor / omega
	// ½·(q/c²)·rho2·a²·Ω² = e
	rm.OmCen = math.Sqrt(2 * e * C0 * C0 / (q * rho2 * a * a))
	rm.RimCen = rm.OmCen * a
	rm.OmBurst = burst / a
	rm.CorAtBurst = math.Abs(float64(m)) * Hbar * rm.OmBurst / E
	rm.CenAtBurst = 0.5 * q * rho2 * a * a * rm.OmBurst * rm.OmBurst / (C0 * C0) / E
	rm.CenCeilingOverHbar = 0.5 * rho2
	return rm
}

// RatioOmega is the rotation that makes the triplet's mass² splittings
// (m = 0 at ω0, m = ±1 at ω1 ± Ω) have the measured ratio
// Δm²31/Δm²21 = R: the ±1 pair is the close pair (split 4ω1Ω in ω²), the
// m = 0 state the third. Solved by bisection on the exact levels.
func RatioOmega(w0, w1, R float64) float64 {
	ratio := func(om float64) float64 {
		lp, lm := sq(w1+om), sq(w1-om)
		d21 := lp - lm
		// Δ31 from the pair member nearer m = 0 (ν1), as the labelling
		// of angles(): ν1, ν2 the close pair
		d31 := math.Min(math.Abs(sq(w0)-lp), math.Abs(sq(w0)-lm))
		return d31 / d21
	}
	lo, hi := 1e-12*w1, 0.5*math.Abs(w1-w0)
	for it := 0; it < 200; it++ {
		mid := math.Sqrt(lo * hi)
		if ratio(mid) > R {
			lo = mid
		} else {
			hi = mid
		}
	}
	return math.Sqrt(lo * hi)
}
