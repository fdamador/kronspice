// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 9 (H42): how two of these rings interact. Two (or N) identical
// E-rings, each the lumped loop of Part 1 (C, L, R_rad + R_d), coupled
// through the field of the other. The rings carry no charge, so their
// whole interaction is the vector potential of one ring's displacement
// current acting along the other (E = −jωA, no scalar potential): for two
// filament loops with uniform current the mutual impedance is the retarded
// Neumann integral
//
//	Z12 = jω·(μ0/4π)·∮∮ dl1·dl2·e^{−jkR}/R,  R = |r1 − r2|
//
// (R. W. P. King, The Theory of Linear Antennas, 1956, ch. III; for one
// loop with itself its real part is the exact radiation resistance of a
// uniform-current loop, η0(π/2)·ka·∫0^{2ka} J2(y) dy, Balanis, Antenna
// Theory, 3rd ed., eq. 5-63, which → (η0π/6)(ka)⁴ of Part 1 as ka → 0).
// It splits into
//
//	the static mutual inductance (μ0/4π)∮∮ dl1·dl2/R: coaxial rings (axes
//	  shared, centres d apart on the axis) by Maxwell's elliptic formula
//	  (em.MutualCircles, Treatise Vol. II art. 701), coplanar side-by-side
//	  rings by the Neumann integral of two 256-gons (em.Neumann);
//	the retardation correction ΔM = (μ0/4π)∮∮ dl1·dl2·(cos kR − 1)/R
//	  (a reactance, added to M);
//	the mutual radiation resistance R_m = ω(μ0/4π)∮∮ dl1·dl2·sin(kR)/R.
//
// The self terms use the same kernel: R11(ω) and ΔL(ω) from the d = 0
// integral, beside Part 1's thick-ring static L. One kernel for self and
// mutual terms matters: subradiance is the small difference R11 − R_m
// (≈ R11·(kd)²/10 for a close coaxial pair), smaller than the (ka)²/5 by
// which the exact R11 differs from Part 1's small-loop formula. So the
// single ring of Part 9 is Part 1's with its first retardation corrections
// (a few per cent in f0 and Q; reported), and every pair and chain is
// compared with that ring.
//
// For far rings the result tends to the dipole–dipole form (Jackson,
// Classical Electrodynamics, 3rd ed., eq. 9.36):
//
//	Z12/R11 → (3/2)·j·e^{−jx}·[ sin²θ/x + (3cos²θ − 1)(1/x³ + j/x²) ],  x = kd,
//
// whose real part is the classical collective decay rate of two dipoles
// (R. H. Lehmberg, Phys. Rev. A 2, 883, 1970; Z. Ficek, R. Tanaś, Phys.
// Rep. 372, 369, 2002): coaxial (θ = 0) R_m/R11 = 3(sin x/x³ − cos x/x²),
// side by side (θ = 90°) (3/2)(sin x/x + cos x/x² − sin x/x³). It is
// printed beside the exact values; close rings need the exact integral.
//
// So the impedance matrix of a pair is
//
//	Z = [[R11 + R_d + jω(L + ΔL) + 1/(jωC), R_m + jω(M + ΔM)],[…, …]],
//
// and of N rings in a row the same with R_m(d_ij), M(d_ij). R11, ΔL, R_d =
// tan δ·ωL, R_m and ΔM depend on ω; each mode is solved self-consistently
// at its own frequency (the matrices frozen there, the frequency updated,
// until it repeats to 1e-12), and the kronspice run of a mode uses the
// values at that frequency.
//
// In kronspice each ring is one series loop C–L–R11–R_d–(0 V sense
// source), the rings coupled by K elements (k = (M + ΔM)/(L + ΔL)) and the
// mutual radiation resistances by current-controlled voltage sources in
// series (H elements, V = R_m·I_j). The H elements' power is the
// interference term of the radiated power, Σ_{i≠j} R_m·I_i·I_j: negative
// (they hand power back) for an antiphase pair, which is subradiance. The
// ledger books it as model-supplied, beside the self radiation on the
// magnetic channel, and the total radiated power Iᵀ·R·I is checked to stay
// ≥ 0 at every step.

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"sort"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
)

// Geom is the arrangement of the rings.
type Geom int

const (
	Coaxial  Geom = iota // axes shared, rings stacked along the axis
	Coplanar             // axes parallel, rings side by side in one plane
)

func (g Geom) String() string {
	return map[Geom]string{Coaxial: "coaxial", Coplanar: "side by side"}[g]
}

// cos2 is cos²θ between the dipole axis and the separation.
func (g Geom) cos2() float64 {
	if g == Coaxial {
		return 1
	}
	return 0
}

// Array is N identical rings in a row, neighbours D apart.
type Array struct {
	Ring ERing
	Geom Geom
	D    float64
	N    int
	mst  map[int]float64 // static mutual inductance by |i − j|
}

// NewArray builds an array and its static mutual inductances.
func NewArray(e ERing, g Geom, d float64, n int) *Array {
	a := &Array{Ring: e, Geom: g, D: d, N: n, mst: map[int]float64{}}
	for k := 1; k < n; k++ {
		a.mst[k] = StaticM(e.A, g, float64(k)*d)
	}
	return a
}

// StaticM is the quasi-static mutual inductance of two filament rings of
// radius a at centre distance d [H].
func StaticM(a float64, g Geom, d float64) float64 {
	if g == Coaxial {
		return em.MutualCircles(a, a, d)
	}
	// n-gons inscribed in the circles: the error falls as 1/n², so the
	// 128- and 256-gon values are extrapolated (Richardson)
	m := func(n int) float64 {
		l1 := em.Loop{P: em.Circle([3]float64{0, 0, 0}, [3]float64{0, 0, 1}, a, n), I: 1}
		l2 := em.Loop{P: em.Circle([3]float64{d, 0, 0}, [3]float64{0, 0, 1}, a, n), I: 1}
		return em.Neumann(l1, l2)
	}
	return (4*m(256) - m(128)) / 3
}

// DipoleZ12 is Z12/R11 of two parallel small loops at x = kd.
func DipoleZ12(g Geom, x float64) complex128 {
	c2 := g.cos2()
	s2 := 1 - c2
	br := complex(s2/x+(3*c2-1)/(x*x*x), (3*c2-1)/(x*x))
	return 1.5 * complex(0, 1) * cmplx.Exp(complex(0, -x)) * br
}

// dipoleStatic is the x → 0 limit of Im Z12/R11, (3/2)(3cos²θ − 1)/x³.
func dipoleStatic(g Geom, x float64) float64 { return 1.5 * (3*g.cos2() - 1) / (x * x * x) }

// retarded is (μ0/4π)∮∮ dl1·dl2·(cos kR − 1)/R and (μ0/4π)∮∮ dl1·dl2·sin(kR)/R
// for two filament circles of radius a, centres d apart (d = 0: one ring
// with itself). Coaxial (and self): one integral over Δ = φ1 − φ2,
// 2πa²∫ cos Δ·f(R) dΔ, R² = d² + 4a²sin²(Δ/2), composite Gauss–Legendre;
// coplanar: the double integral by the trapezoidal rule (periodic and
// smooth, so it converges geometrically).
func retarded(a float64, g Geom, d, k float64) (gc, gs float64) {
	f := func(r float64) (float64, float64) {
		if r < 1e-300 {
			return 0, k
		}
		kr := k * r
		if kr < 1e-4 { // series, to keep (cos − 1)/R accurate
			return -k * kr / 2 * (1 - kr*kr/12), k * (1 - kr*kr/6)
		}
		return (math.Cos(kr) - 1) / r, math.Sin(kr) / r
	}
	if g == Coaxial || d == 0 {
		const panels = 128
		h := math.Pi / panels
		for p := 0; p < panels; p++ {
			for q := range gl8x {
				dl := (float64(p) + 0.5*(1+gl8x[q])) * h
				sn := math.Sin(dl / 2)
				r := math.Sqrt(d*d + 4*a*a*sn*sn)
				c, s := f(r)
				w := 0.5 * h * gl8w[q] * math.Cos(dl)
				gc += w * c
				gs += w * s
			}
		}
		// ∫0^{2π} = 2∫0^π (even about π); ∮∮ = 2πa²∫dΔ
		sc := Mu0 / (4 * math.Pi) * 2 * math.Pi * a * a * 2
		return gc * sc, gs * sc
	}
	const n = 256
	h := 2 * math.Pi / n
	for i := 0; i < n; i++ {
		s1, c1 := math.Sincos(float64(i) * h)
		for j := 0; j < n; j++ {
			s2, c2 := math.Sincos(float64(j) * h)
			dx, dy := a*c1-(d+a*c2), a*s1-a*s2
			c, s := f(math.Hypot(dx, dy))
			dot := c1*c2 + s1*s2
			gc += dot * c
			gs += dot * s
		}
	}
	sc := Mu0 / (4 * math.Pi) * a * a * h * h
	return gc * sc, gs * sc
}

var gl8x = [8]float64{-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363}
var gl8w = [8]float64{0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
	0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763}

// Self returns one ring's radiation resistance R11 and retardation
// inductance ΔL at ω.
func (a *Array) Self(w float64) (r11, dl float64) {
	gc, gs := retarded(a.Ring.A, a.Geom, 0, w/C0)
	return w * gs, gc
}

// R11 is one ring's radiation resistance at ω (exact uniform loop).
func (a *Array) R11(w float64) float64 { r, _ := a.Self(w); return r }

// LSelf is L + ΔL at ω.
func (a *Array) LSelf(w float64) float64 { _, dl := a.Self(w); return a.Ring.L() + dl }

// RD is the dielectric loss resistance tan δ·ωL at ω.
func (a *Array) RD(w float64) float64 { return a.Ring.TanD * w * a.Ring.L() }

// Coupling returns M + ΔM and R_m between rings i, j at ω.
func (a *Array) Coupling(i, j int, w float64) (m, rm float64) {
	k := i - j
	if k < 0 {
		k = -k
	}
	gc, gs := retarded(a.Ring.A, a.Geom, float64(k)*a.D, w/C0)
	return a.mst[k] + gc, w * gs
}

// Mats returns the inductance and resistance matrices at ω.
func (a *Array) Mats(w float64) (L, R [][]float64) {
	n := a.N
	r11, dl := a.Self(w)
	rd := a.RD(w)
	cm, cr := map[int]float64{}, map[int]float64{}
	for k := 1; k < n; k++ {
		cm[k], cr[k] = a.Coupling(0, k, w)
	}
	L, R = make([][]float64, n), make([][]float64, n)
	for i := range L {
		L[i], R[i] = make([]float64, n), make([]float64, n)
		for j := 0; j < n; j++ {
			if i == j {
				L[i][j] = a.Ring.L() + dl
				R[i][j] = r11 + rd
				continue
			}
			k := i - j
			if k < 0 {
				k = -k
			}
			L[i][j], R[i][j] = cm[k], cr[k]
		}
	}
	return
}

// AMode is one mode of the array: currents I(t) = Re(V·e^{st}).
type AMode struct {
	S      complex128
	F, Q   float64
	QRad   float64 // radiation alone: ω·(vᵀLv)/(vᵀR_rad v) at the same vector
	W      float64 // frequency the matrices were evaluated at [rad/s]
	V      []complex128
	Dipole float64 // |Σ V_i| / (√N·|V|): the share of the net magnetic moment
}

// csolve solves A·x = b (complex, partial pivoting); A is overwritten.
func csolve(A [][]complex128, b []complex128) []complex128 {
	n := len(b)
	for c := 0; c < n; c++ {
		p := c
		for r := c + 1; r < n; r++ {
			if cmplx.Abs(A[r][c]) > cmplx.Abs(A[p][c]) {
				p = r
			}
		}
		A[c], A[p] = A[p], A[c]
		b[c], b[p] = b[p], b[c]
		if A[c][c] == 0 {
			A[c][c] = 1e-300
		}
		for r := c + 1; r < n; r++ {
			f := A[r][c] / A[c][c]
			if f == 0 {
				continue
			}
			for k := c; k < n; k++ {
				A[r][k] -= f * A[c][k]
			}
			b[r] -= f * b[c]
		}
	}
	x := make([]complex128, n)
	for r := n - 1; r >= 0; r-- {
		s := b[r]
		for k := r + 1; k < n; k++ {
			s -= A[r][k] * x[k]
		}
		x[r] = s / A[r][r]
	}
	return x
}

// quadForm is vᵀ·M·v (no conjugate: Z is complex symmetric).
func quadForm(M [][]float64, v []complex128) complex128 {
	var s complex128
	for i := range v {
		var r complex128
		for j := range v {
			r += complex(M[i][j], 0) * v[j]
		}
		s += v[i] * r
	}
	return s
}

// refine finds the eigenpair of s²L + sR + 1/C near (s, v) with the
// matrices at ω fixed: Rayleigh-functional root (the quadratic vᵀT(s)v = 0
// solved for the root nearest s) and inverse iteration, alternately.
func refine(L, R [][]float64, invC float64, s complex128, v []complex128) (complex128, []complex128) {
	n := len(v)
	for it := 0; it < 40; it++ {
		a, b := quadForm(L, v), quadForm(R, v)
		var c complex128
		for i := range v {
			c += complex(invC, 0) * v[i] * v[i]
		}
		disc := cmplx.Sqrt(b*b - 4*a*c)
		r1, r2 := (-b+disc)/(2*a), (-b-disc)/(2*a)
		sn := r1
		if cmplx.Abs(r2-s) < cmplx.Abs(r1-s) {
			sn = r2
		}
		done := cmplx.Abs(sn-s) < 1e-15*cmplx.Abs(s) && it > 2
		s = sn
		if done {
			break
		}
		// inverse iteration with a shift slightly off s (T(s) is singular at the root)
		sh := s * complex(1+1e-13, 1e-13)
		T := make([][]complex128, n)
		rhs := make([]complex128, n)
		for i := range T {
			T[i] = make([]complex128, n)
			for j := range T[i] {
				T[i][j] = sh*sh*complex(L[i][j], 0) + sh*complex(R[i][j], 0)
			}
			T[i][i] += complex(invC, 0)
			for j := 0; j < n; j++ {
				rhs[i] += (2*sh*complex(L[i][j], 0) + complex(R[i][j], 0)) * v[j]
			}
		}
		u := csolve(T, rhs)
		var nn float64
		for _, z := range u {
			nn += real(z * cmplx.Conj(z))
		}
		nn = math.Sqrt(nn)
		// fix the phase: largest component real and positive
		big := 0
		for i := range u {
			if cmplx.Abs(u[i]) > cmplx.Abs(u[big]) {
				big = i
			}
		}
		ph := cmplx.Exp(complex(0, -cmplx.Phase(u[big])))
		for i := range u {
			v[i] = u[i] * ph / complex(nn, 0)
		}
	}
	return s, v
}

// SolveMode follows a mode from a start vector, re-evaluating the
// matrices at the mode's own frequency until it repeats to 1e-12.
func (a *Array) SolveMode(v0 []float64) AMode {
	invC := 1 / a.Ring.C()
	w := a.Ring.Omega0()
	v := make([]complex128, a.N)
	var nn float64
	for _, x := range v0 {
		nn += x * x
	}
	for i := range v {
		v[i] = complex(v0[i]/math.Sqrt(nn), 0)
	}
	L, R := a.Mats(w)
	lv, rv := real(quadForm(L, v)), real(quadForm(R, v))
	s := complex(-rv/(2*lv), 1/math.Sqrt(lv*a.Ring.C()))
	for it := 0; it < 30; it++ {
		w = imag(s)
		L, R = a.Mats(w)
		sn, vn := refine(L, R, invC, s, v)
		done := math.Abs(imag(sn)-imag(s)) < 1e-12*imag(s) && it > 0
		s, v = sn, vn
		if done {
			break
		}
	}
	m := AMode{S: s, W: w, V: v, F: imag(s) / (2 * math.Pi), Q: imag(s) / (-2 * real(s))}
	// radiation-only Q at the same vector (first order): replace R by R_rad
	Rr := make([][]float64, a.N)
	for i := range Rr {
		Rr[i] = append([]float64(nil), R[i]...)
		Rr[i][i] -= a.RD(w)
	}
	vr := make([]complex128, a.N)
	var sum complex128
	var vn float64
	for i := range v {
		vr[i] = complex(real(v[i]), 0)
		sum += v[i]
		vn += real(v[i] * cmplx.Conj(v[i]))
	}
	m.QRad = w * real(quadForm(L, vr)) / real(quadForm(Rr, vr))
	m.Dipole = cmplx.Abs(sum) / math.Sqrt(float64(a.N)*vn)
	return m
}

// PairModes returns the in-phase (I2 = I1) and antiphase (I2 = −I1)
// modes of a pair.
func (a *Array) PairModes() (in, anti AMode) {
	return a.SolveMode([]float64{1, 1}), a.SolveMode([]float64{1, -1})
}

// PairAnalytic is the closed form for a pair at ω (both modes share the
// eigenvectors (1, ±1)): ω± = 1/√((L ± M)C), Q± = ω±(L ± M)/(R11 + R_d ± R_m).
func (a *Array) PairAnalytic(sign float64, w float64) (f, q float64) {
	m, rm := a.Coupling(0, 1, w)
	l := a.LSelf(w) + sign*m
	ww := 1 / math.Sqrt(l*a.Ring.C())
	return ww / (2 * math.Pi), ww * l / (a.R11(w) + a.RD(w) + sign*rm)
}

// ChainModes returns the brightest and darkest modes of the array: the
// eigenvectors of L at ω0 ranked by first-order Q, the darkest three and
// the brightest refined, the darkest of the refined kept.
func (a *Array) ChainModes() (bright, dark AMode, all []AMode) {
	L, R := a.Mats(a.Ring.Omega0())
	vals, vec := eigSym(L)
	type cand struct {
		k int
		q float64
	}
	var cs []cand
	for k := range vals {
		v := make([]float64, a.N)
		var rv float64
		for i := range v {
			v[i] = vec[i][k]
		}
		for i := range v {
			for j := range v {
				rv += v[i] * R[i][j] * v[j]
			}
		}
		cs = append(cs, cand{k, math.Sqrt(vals[k]/a.Ring.C()) / rv})
	}
	sort.Slice(cs, func(i, j int) bool { return cs[i].q > cs[j].q })
	pick := []int{}
	for i := 0; i < len(cs) && i < 3; i++ {
		pick = append(pick, cs[i].k)
	}
	pick = append(pick, cs[len(cs)-1].k)
	dark.Q = -1
	for idx, k := range pick {
		v := make([]float64, a.N)
		for i := range v {
			v[i] = vec[i][k]
		}
		m := a.SolveMode(v)
		all = append(all, m)
		if idx == len(pick)-1 {
			bright = m
		} else if m.Q > dark.Q {
			dark = m
		}
	}
	return
}

// ArrayRun is a kronspice transient of the array.
type ArrayRun struct {
	Time       []float64
	I, VC      [][]float64 // [ring][t] loop current, capacitor voltage
	U          []float64   // ½ΣLI² + ΣMII + ½ΣCV² [J]
	UKron      []float64   // kronspice's own stored energy [J]
	Q          float64     // from ln U, after the first period
	SelfRad    float64     // energy radiated by the self radiation resistances [J]
	Heat       float64     // dielectric heat [J]
	Mutual     float64     // energy absorbed by the H elements (interference term) [J]
	MinRadP    float64     // smallest total radiated power Iᵀ·R_rad·I seen / its largest
	Mismatch   float64     // (stored start − stored end − radiated − heat) / stored start
	LedgerRel  float64     // kronspice's worst relative ledger residual
	StoredFrom float64     // kronspice stored energy at the start vs U (relative difference)
	W          float64
	L, R       [][]float64
}

// Run builds the array with its matrices at ω and initial loop currents
// i0 and capacitor voltages v0, and integrates for `periods` periods of ω
// at `per` steps a period.
func (a *Array) Run(w float64, i0, v0 []float64, periods float64, per int) (ArrayRun, error) {
	L, R := a.Mats(w)
	b := circuit.NewBuilder(fmt.Sprintf("%d E-rings, %s, d = %.4g m", a.N, a.Geom, a.D))
	s := strconv.Itoa
	r11, rd := a.R11(w), a.RD(w)
	// kronspice's K takes 0 < k ≤ 1, so a negative mutual inductance M_ij
	// is made with a helper inductor X in ring j, wound backwards (the loop
	// current runs through it from its second node to its first) and
	// coupled to L_i's main inductor with k = ½: L_X = 4M²/L_i(main),
	// taken off L_j's main inductor, so ring j's total self-inductance is
	// unchanged and the mutual term is −½√(L_i(main)·L_X) = M exactly
	// (rings in order, so L_i(main) is final before it is used).
	main := make([]float64, a.N)
	for i := range main {
		main[i] = L[i][i]
	}
	type helper struct {
		i, j int
		lx   float64
	}
	var helpers []helper
	for j := 0; j < a.N; j++ {
		for i := 0; i < j; i++ {
			if L[i][j] < 0 {
				lx := 4 * L[i][j] * L[i][j] / main[i]
				helpers = append(helpers, helper{i, j, lx})
				main[j] -= lx
			}
		}
	}
	for i := 0; i < a.N; i++ {
		x := s(i)
		b.CapacitorIC("C"+x, "a"+x, "b"+x, a.Ring.C(), v0[i])
		b.InductorIC("L"+x, "b"+x, "c"+x, main[i], i0[i])
		from := "c" + x
		for _, h := range helpers {
			if h.j == i {
				nx := "x" + x + "_" + s(h.i)
				b.InductorIC("X"+x+"_"+s(h.i), nx, from, h.lx, -i0[i])
				// a shunt so the UIC start (every inductor a pinned current) is
				// not two current sources in series. X carries the coupling EMF
				// ωM·I_i (up to ~1 Ω·I), so the shunt must be large: at 1e8 Ω it
				// takes ≤ 1e-8·I² W, under 1e-5 of even the darkest chain mode's
				// radiation (booked as heat; 1e4 Ω cost the N = 16 chain a third
				// of its Q)
				b.Resistor("RX"+x+"_"+s(h.i), nx, from, 1e8)
				from = nx
			}
		}
		if rd > 0 {
			b.ResistorTo("RR"+x, from, "d"+x, r11, exchange.Magnetic)
			b.Resistor("RD"+x, "d"+x, "e"+x, rd)
		} else {
			b.ResistorTo("RR"+x, from, "e"+x, r11, exchange.Magnetic)
		}
		b.VoltageSource("VS"+x, "e"+x, "h"+x+"_0", 0)
		node := "h" + x + "_0"
		k := 0
		for j := 0; j < a.N; j++ {
			if j == i {
				continue
			}
			k++
			next := "h" + x + "_" + s(k)
			if k == a.N-1 {
				next = "a" + x
			}
			b.CCVS("H"+x+"_"+s(j), node, next, "VS"+s(j), R[i][j])
			node = next
		}
		if a.N == 1 {
			b.Resistor("RZ"+x, node, "a"+x, 1e-12)
		}
		// one tie to ground per loop: each ring is a closed loop touching
		// ground at one node, so no current flows through the tie
		b.Resistor("RB"+x, "a"+x, "0", 1e-3)
	}
	for i := 0; i < a.N; i++ {
		for j := i + 1; j < a.N; j++ {
			if L[i][j] > 0 {
				b.Mutual("K"+s(i)+"_"+s(j), "L"+s(i), "L"+s(j), L[i][j]/math.Sqrt(main[i]*main[j]))
			}
		}
	}
	for _, h := range helpers {
		b.Mutual("KX"+s(h.j)+"_"+s(h.i), "L"+s(h.i), "X"+s(h.j)+"_"+s(h.i), 0.5)
	}
	ckt, err := b.BuildTransient()
	if err != nil {
		return ArrayRun{}, err
	}
	h := 2 * math.Pi / w / float64(per)
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: h, Stop: periods * 2 * math.Pi / w, MaxStep: h, UIC: true}, analysis.TranOptions{MaxSteps: 20_000_000})
	if err != nil {
		return ArrayRun{}, err
	}
	run := ArrayRun{Time: wf.Time, W: w, L: L, R: R, UKron: wf.Ledger.StoredEnergy, LedgerRel: wf.Ledger.MaxRelativeResidual}
	nt := len(wf.Time)
	run.I, run.VC = make([][]float64, a.N), make([][]float64, a.N)
	for i := 0; i < a.N; i++ {
		x := s(i)
		run.I[i], _ = wf.ElementCurrent("L" + x)
		va, _ := wf.NodeVoltage("a" + x)
		vb, _ := wf.NodeVoltage("b" + x)
		run.VC[i] = make([]float64, nt)
		for t := range va {
			run.VC[i][t] = va[t] - vb[t]
		}
		run.SelfRad += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RR"+x)
		if rd > 0 {
			run.Heat += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RD"+x)
		}
		for _, h := range helpers {
			if h.j == i {
				run.Heat += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RX"+x+"_"+s(h.i))
			}
		}
		for j := 0; j < a.N; j++ {
			if j != i {
				run.Mutual += elementEnergy(wf.Time, wf.ElementNames, wf.P, "H"+x+"_"+s(j))
			}
		}
	}
	run.U = make([]float64, nt)
	maxP, minP := 0.0, math.Inf(1)
	for t := 0; t < nt; t++ {
		var u, p float64
		for i := 0; i < a.N; i++ {
			u += 0.5 * a.Ring.C() * run.VC[i][t] * run.VC[i][t]
			for j := 0; j < a.N; j++ {
				u += 0.5 * L[i][j] * run.I[i][t] * run.I[j][t]
				rr := R[i][j]
				if i == j {
					rr -= rd
				}
				p += rr * run.I[i][t] * run.I[j][t]
			}
		}
		run.U[t] = u
		maxP, minP = math.Max(maxP, p), math.Min(minP, p)
	}
	if maxP > 0 {
		run.MinRadP = minP / maxP
	}
	run.Mismatch = (run.U[0] - run.U[nt-1] - run.SelfRad - run.Heat - run.Mutual) / run.U[0]
	run.StoredFrom = (run.UKron[0] - run.U[0]) / run.U[0]
	// Q from ln U after the first period
	T := 2 * math.Pi / w
	var xs, ys []float64
	for t := 0; t < nt; t++ {
		if wf.Time[t] > T && run.U[t] > 0 {
			xs, ys = append(xs, wf.Time[t]), append(ys, math.Log(run.U[t]))
		}
	}
	if len(xs) > 3 {
		run.Q = -w / fitSlope(xs, ys)
	}
	return run, nil
}

// RunMode starts the array in a mode (I = Re V, V_C = Re(V/(sC)), scaled
// to 1 J) and runs it.
func (a *Array) RunMode(m AMode, periods float64, per int) (ArrayRun, error) {
	i0, v0 := make([]float64, a.N), make([]float64, a.N)
	for i := range m.V {
		i0[i] = real(m.V[i])
		v0[i] = real(m.V[i] / (m.S * complex(a.Ring.C(), 0)))
	}
	L, _ := a.Mats(m.W)
	var u float64
	for i := 0; i < a.N; i++ {
		u += 0.5 * a.Ring.C() * v0[i] * v0[i]
		for j := 0; j < a.N; j++ {
			u += 0.5 * L[i][j] * i0[i] * i0[j]
		}
	}
	sc := 1 / math.Sqrt(u)
	for i := range i0 {
		i0[i] *= sc
		v0[i] *= sc
	}
	return a.Run(m.W, i0, v0, periods, per)
}

// ringEnergy is ring i's own energy ½LI² + ½CV² at sample t.
func (r ArrayRun) ringEnergy(a *Array, i, t int) float64 {
	return 0.5*r.L[i][i]*r.I[i][t]*r.I[i][t] + 0.5*a.Ring.C()*r.VC[i][t]*r.VC[i][t]
}

// Beats is ring 0 started alone (capacitor charged), ring 1 at rest.
type Beats struct {
	D, KD          float64
	DfMatrix       float64 // f(antiphase) − f(in-phase) from the self-consistent modes [Hz]
	DfMeasured     float64 // from the zero crossings of E0 − E1 [Hz]
	DfFrozen       float64 // the same pair with its matrices frozen at the run's ω (what kronspice integrates) [Hz]
	Crossings      int
	MaxTransfer    float64 // largest E1/(E0 + E1)
	TMaxTransfer   float64 // when [periods]
	DarkShare      float64 // at the end: share of the energy in the antiphase pattern
	EnergyLeft     float64 // at the end, of the start
	Run            ArrayRun
	FIn, FAnti     float64
	QIn, QAnti     float64
	PeriodsRun     float64
	CrossingTimesT []float64
}

// RunBeats runs the beat test on a pair at distance d.
func RunBeats(e ERing, g Geom, d, periods float64) (Beats, error) {
	a := NewArray(e, g, d, 2)
	in, anti := a.PairModes()
	w := 0.5 * (imag(in.S) + imag(anti.S))
	bt := Beats{D: d, KD: e.Omega0() * d / C0, DfMatrix: anti.F - in.F, FIn: in.F, FAnti: anti.F, QIn: in.Q, QAnti: anti.Q, PeriodsRun: periods}
	run, err := a.Run(w, []float64{0, 0}, []float64{1, 0}, periods, 200)
	if err != nil {
		return bt, err
	}
	bt.Run = run
	{
		m, _ := a.Coupling(0, 1, w)
		ls := a.LSelf(w)
		bt.DfFrozen = (1/math.Sqrt((ls-m)*e.C()) - 1/math.Sqrt((ls+m)*e.C())) / (2 * math.Pi)
	}
	nt := len(run.Time)
	T := 2 * math.Pi / w
	diff := make([]float64, nt)
	for t := 0; t < nt; t++ {
		e0, e1 := run.ringEnergy(a, 0, t), run.ringEnergy(a, 1, t)
		diff[t] = e0 - e1
		if f := e1 / (e0 + e1); f > bt.MaxTransfer {
			bt.MaxTransfer, bt.TMaxTransfer = f, run.Time[t]/T
		}
	}
	// carrier-average E0 − E1 over one period (moving window), then zero crossings
	win := int(math.Round(T / (run.Time[1] - run.Time[0])))
	var avg []float64
	var tav []float64
	var acc float64
	for t := 0; t < nt; t++ {
		acc += diff[t]
		if t >= win {
			acc -= diff[t-win]
		}
		if t >= win-1 {
			avg = append(avg, acc/float64(win))
			tav = append(tav, run.Time[t]-0.5*float64(win-1)*(run.Time[1]-run.Time[0]))
		}
	}
	var zc []float64
	for k := 1; k < len(avg); k++ {
		if (avg[k-1] > 0) != (avg[k] > 0) && avg[k-1] != 0 {
			fr := avg[k-1] / (avg[k-1] - avg[k])
			tz := tav[k-1] + fr*(tav[k]-tav[k-1])
			// stop where the signal is down to round-off
			if math.Abs(avg[k-1]) < 1e-9*math.Abs(avg[0]) {
				break
			}
			zc = append(zc, tz)
			bt.CrossingTimesT = append(bt.CrossingTimesT, tz/T)
		}
	}
	bt.Crossings = len(zc)
	if len(zc) >= 2 {
		// spacing 1/(2Δf): least-squares slope of crossing time against index
		idx := make([]float64, len(zc))
		for k := range idx {
			idx[k] = float64(k)
		}
		bt.DfMeasured = 1 / (2 * fitSlope(idx, zc))
	}
	// end state: projection on (1, −1)/√2 of current and voltage
	last := nt - 1
	var ua float64
	for _, sg := range []float64{1} {
		ia := (run.I[0][last] - sg*run.I[1][last]) / math.Sqrt2
		va := (run.VC[0][last] - sg*run.VC[1][last]) / math.Sqrt2
		m, _ := a.Coupling(0, 1, w)
		ua = 0.5*(run.L[0][0]-m)*ia*ia + 0.5*e.C()*va*va
	}
	bt.DarkShare = ua / run.U[last]
	bt.EnergyLeft = run.U[last] / run.U[0]
	return bt, nil
}

// ForceAt is the time-averaged force between two rings of a pair in a
// mode, F = ⟨I1·I2⟩·∂M/∂d (energy method at constant currents; M is the
// full mutual inductance, static + retardation), per stored energy U and
// in units of U/a. Positive = repulsion (pushes d up).
type ForceAt struct {
	D, DMdd          float64
	FInA, FAntiA     float64 // analytic: ±(U/2(L ± M))·∂M/∂d·a/U
	FInK, FAntiK     float64 // from the kronspice mode runs (⟨I1I2⟩ over the first 3 periods / ⟨U⟩)
	FRandEarly       float64 // random relative phases (32 pairs): ensemble ⟨F⟩·a/⟨U⟩, first period
	FRandEarlySigma  float64 // its standard error
	FRandLate        float64 // the same over the last period of the run
	FRandLateSigma   float64
	LateEnergyLeft   float64
	FRandPeriods     float64
	CoherentInRatio  float64 // FRandLate / FAntiA
	InResidual       float64
	AntiResidual     float64
	InQ, AntiQ       float64
	InMinRad         float64
	AntiMinRad       float64
	InMismatch       float64
	AntiMismatch     float64
	RandomMaxMismach float64
}

// dMdd is ∂M/∂d at ω by a central difference (step d·1e-4).
func dMdd(e ERing, g Geom, d, w float64) float64 {
	h := d * 1e-4
	mp, _ := NewArray(e, g, d+h, 2).Coupling(0, 1, w)
	mm, _ := NewArray(e, g, d-h, 2).Coupling(0, 1, w)
	return (mp - mm) / (2 * h)
}

// meanI1I2 is ⟨I0·I1⟩ and ⟨U⟩ over the samples with time in [t0, t1].
func meanI1I2(r ArrayRun, t0, t1 float64) (p, u float64) {
	var n float64
	for t := range r.Time {
		if r.Time[t] >= t0 && r.Time[t] <= t1 {
			p += r.I[0][t] * r.I[1][t]
			u += r.U[t]
			n++
		}
	}
	return p / n, u / n
}

// Force computes the force at distance d.
func Force(e ERing, g Geom, d float64, randPeriods float64, seed uint64) (ForceAt, error) {
	a := NewArray(e, g, d, 2)
	in, anti := a.PairModes()
	f := ForceAt{D: d, InQ: in.Q, AntiQ: anti.Q, FRandPeriods: randPeriods}
	for _, md := range []struct {
		m    AMode
		sg   float64
		a, k *float64
		res  *float64
		mr   *float64
		mm   *float64
	}{{in, 1, &f.FInA, &f.FInK, &f.InResidual, &f.InMinRad, &f.InMismatch}, {anti, -1, &f.FAntiA, &f.FAntiK, &f.AntiResidual, &f.AntiMinRad, &f.AntiMismatch}} {
		dm := dMdd(e, g, d, md.m.W)
		m, _ := a.Coupling(0, 1, md.m.W)
		*md.a = md.sg * dm / (2 * (a.LSelf(md.m.W) + md.sg*m)) * e.A
		run, err := a.RunMode(md.m, 3, 200)
		if err != nil {
			return f, err
		}
		T := 2 * math.Pi / md.m.W
		p, u := meanI1I2(run, 0, 3*T)
		*md.k = p * dm / u * e.A
		*md.res, *md.mr, *md.mm = run.LedgerRel, run.MinRadP, run.Mismatch
		if md.sg > 0 {
			f.DMdd = dm
		}
	}
	// random relative phases: each ring started as a free oscillation of
	// unit amplitude with phase φ_i (V_C = cos φ, I = −sin φ·√(C/L)·...)
	w := 0.5 * (imag(in.S) + imag(anti.S))
	T := 2 * math.Pi / w
	dm := dMdd(e, g, d, w)
	rng := newPCG(seed)
	var early, late, early2, late2, uE, uL, left float64
	const nr = 32
	z := math.Sqrt(e.C() / a.LSelf(e.Omega0()))
	for k := 0; k < nr; k++ {
		ph := [2]float64{2 * math.Pi * rng.Float64(), 2 * math.Pi * rng.Float64()}
		i0 := []float64{-z * math.Sin(ph[0]), -z * math.Sin(ph[1])}
		v0 := []float64{math.Cos(ph[0]), math.Cos(ph[1])}
		run, err := a.Run(w, i0, v0, randPeriods, 120)
		if err != nil {
			return f, err
		}
		pe, ue := meanI1I2(run, 0, T)
		tl := randPeriods * T
		pl, ul := meanI1I2(run, tl-T, tl)
		early += pe
		early2 += pe * pe
		late += pl
		late2 += pl * pl
		uE += ue
		uL += ul
		left += run.U[len(run.U)-1] / run.U[0]
		f.RandomMaxMismach = math.Max(f.RandomMaxMismach, math.Abs(run.Mismatch))
	}
	mE, mL := early/nr, late/nr
	sE := math.Sqrt(math.Max(0, early2/nr-mE*mE) / nr)
	sL := math.Sqrt(math.Max(0, late2/nr-mL*mL) / nr)
	f.FRandEarly = mE * dm / (uE / nr) * e.A
	f.FRandEarlySigma = sE * math.Abs(dm) / (uE / nr) * e.A
	f.FRandLate = mL * dm / (uL / nr) * e.A
	f.FRandLateSigma = sL * math.Abs(dm) / (uL / nr) * e.A
	f.LateEnergyLeft = left / nr
	f.CoherentInRatio = f.FRandLate / f.FAntiA
	return f, nil
}

// newPCG is a seeded generator (CLAUDE.md §4: never the global source).
func newPCG(seed uint64) *rand.Rand { return rand.New(rand.NewPCG(seed, 0x9e3779b97f4a7c15)) }
