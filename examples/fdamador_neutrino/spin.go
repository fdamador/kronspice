// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 7a: a rotating ring (H40, user hypothesis 2026-11-07): "If the
// torus spins about its axis and spins around the E axis (the ring), would
// this also add to the missing information of neutrino data?"
//
// Reading 1, a material ring that turns. A ring turning at Ω makes its
// co- and counter-rotating waves non-degenerate: the Sagnac effect
// (G. Sagnac, C. R. Acad. Sci. 157, 708, 1913; E. J. Post, Rev. Mod. Phys.
// 39, 475, 1967). In the ring's own frame a wave going round once is
// delayed by Δt = 4AΩ/c² against the one going the other way, for any
// medium (Fresnel drag of a co-moving medium makes the delay universal,
// Post §IV), so a resonance of a ring of radius a whose guided wave has
// group velocity v_g splits by
//
//	Δω/ω = Δt/T_rt = 2aΩ·v_g/c²,      T_rt = 2πa/v_g,
//
// which is spacetime.SagnacSplitting·(v_g/c) (the vacuum ring-laser value
// for v_g = c; a co-moving dielectric of index n gives 2mΩ/n²).
//
// In the ring's frame the rotation enters the circuit as a gyroscopic
// (Coriolis) term: a velocity-proportional force that is antisymmetric and
// so does no work (L. D. Landau, E. M. Lifshitz, Mechanics, §39; the same
// structure as Kron's speed-voltage tensor G, Equivalent Circuits of
// Electric Machinery, 1951, Epilogue eq. 9). For the CL ring line with
// node voltages u and shunt currents g it is an antisymmetric current
// 2·Gy·Q·u added to Kirchhoff's current law,
//
//	Cmat·du/dt = −g + 2·Gy·Q·u,   Lg·dg/dt = u,   Qᵀ = −Q,
//
// and in kronspice it is a gyrator (device/gyrator.go: i1 = Gy·V2,
// i2 = −Gy·V1, lossless) across every section. On a uniform ring the CW
// wave e^{−iκj} and the CCW wave e^{+iκj} then sit at
// c_κω² ± 4·Gy·sin κ·ω − 1/Lg = 0, split exactly by
//
//	Δω = 4·|Gy|·sin κ / c_κ,    c_κ = 4C sin²(κ/2) + Cg.
//
// Gy is set from the non-reciprocity ν = Δω/ω wanted on a chosen pair;
// ν is converted to the rim speed aΩ = ν·c²/(2v_g) and checked against the
// burst speed of the strongest rotor material in spacetime/ (T1000G carbon
// fibre, 1,881 m/s). Ledger: the gyrators' total power is zero at every
// instant (Σ u·Q·u = 0); booked as the rotation's (mechanical) work, it
// is reported, never counted as inflow. The ring's kinetic energy, ½v² per
// kilogram, is paid by whatever turns it.
//
// With rotation the modes are complex (travelling components), so the
// 3×3 mixing of Part 6 can carry a CP phase: J = Im(U_e1 U_μ2 U*_e2 U*_μ1)
// (C. Jarlskog, Phys. Rev. Lett. 55, 1039, 1985). The question is how much
// rotation that needs.

import (
	"math"
	"math/cmplx"
	"runtime"
	"sort"
	"sync"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
	"github.com/fdamador/kronspice/spacetime"
)

// RotModes are the exact modes of a ring with gyrators (no loss).
type RotModes struct {
	N         int
	Omega     []float64      // positive frequencies, ascending
	V         [][]complex128 // [mode][2N]: unit eigenvectors of A in y = (Cmat^½·u, √Lg·g), A·v = iω·v
	U         [][]complex128 // [mode][N]: node-voltage pattern û (z(t) = û·e^{iωt}), Euclidean unit
	CHalf     [][]float64    // Cmat^½
	CHalfInv  [][]float64    // Cmat^−½
	Lg        float64
	ResidualA float64 // largest |A·v − iω·v|
}

// cmatQ builds the capacitance matrix of Modes() and the antisymmetric
// gyrator incidence Q (KCL: Cmat·u̇ = −g + 2·Gy·Q·u).
func (r FRing) cmatQ() (cm, q [][]float64) {
	n := r.N
	cm, q = make([][]float64, n), make([][]float64, n)
	for i := range cm {
		cm[i], q[i] = make([]float64, n), make([]float64, n)
	}
	for j := 0; j < n; j++ {
		c := r.C * (1 + r.Def[j])
		m := (j + 1) % n
		s := 1.0
		if r.Crossed && j == n-1 {
			s = -1
		}
		cm[j][j] += c + r.Cg + r.NodeC[j]
		cm[m][m] += c
		cm[j][m] -= s * c
		cm[m][j] -= s * c
		// gyrator j: port 1 = pair j, port 2 = pair m (sign s): current
		// Gy·s·u_m leaves a_j, Gy·s·u_j enters a_m.
		q[j][m] -= s
		q[m][j] += s
	}
	return
}

// symFunc returns f(M) for a symmetric matrix M.
func symFunc(m [][]float64, f func(float64) float64) [][]float64 {
	vals, v := eigSym(m)
	n := len(m)
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, n)
		for j := range out[i] {
			var s float64
			for k := 0; k < n; k++ {
				s += v[i][k] * f(vals[k]) * v[j][k]
			}
			out[i][j] = s
		}
	}
	return out
}

func matMul(a, b [][]float64) [][]float64 {
	n, m, p := len(a), len(b), len(b[0])
	out := make([][]float64, n)
	for i := range out {
		out[i] = make([]float64, p)
		for k := 0; k < m; k++ {
			if a[i][k] == 0 {
				continue
			}
			for j := 0; j < p; j++ {
				out[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return out
}

// SystemA is the real antisymmetric generator ẏ = A·y of the lossless ring.
func (r FRing) SystemA() (a, ch, chi [][]float64) {
	n := r.N
	cm, q := r.cmatQ()
	ch = symFunc(cm, math.Sqrt)
	chi = symFunc(cm, func(x float64) float64 { return 1 / math.Sqrt(x) })
	s := matMul(matMul(chi, q), chi)
	a = make([][]float64, 2*n)
	for i := range a {
		a[i] = make([]float64, 2*n)
	}
	sl := 1 / math.Sqrt(r.Lg)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			a[i][j] = 2 * r.Gyro * s[i][j]
			a[i][n+j] = -chi[i][j] * sl
			a[n+i][j] = chi[i][j] * sl
		}
	}
	// exact antisymmetry (round-off of the products)
	for i := 0; i < 2*n; i++ {
		for j := i + 1; j < 2*n; j++ {
			x := 0.5 * (a[i][j] - a[j][i])
			a[i][j], a[j][i] = x, -x
		}
		a[i][i] = 0
	}
	return
}

// RotModes solves the lossless ring with gyrators: the eigenvalues of
// −A² = AᵀA are ω² (each twice); for an eigenvector p, v = p − iAp/ω
// satisfies A·v = iω·v. Duplicates (v and i·v) are removed by complex
// Gram–Schmidt, which also keeps exactly degenerate pairs orthogonal.
func (r FRing) RotModes() RotModes {
	n := r.N
	a, ch, chi := r.SystemA()
	ata := make([][]float64, 2*n)
	for i := range ata {
		ata[i] = make([]float64, 2*n)
		for j := range ata[i] {
			var s float64
			for k := 0; k < 2*n; k++ {
				s += a[k][i] * a[k][j]
			}
			ata[i][j] = s
		}
	}
	vals, p := eigSym(ata)
	idx := make([]int, 2*n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(x, y int) bool { return vals[idx[x]] < vals[idx[y]] })
	rm := RotModes{N: n, CHalf: ch, CHalfInv: chi, Lg: r.Lg}
	for _, k := range idx {
		w := math.Sqrt(math.Max(vals[k], 0))
		v := make([]complex128, 2*n)
		for i := 0; i < 2*n; i++ {
			var ap float64
			for j := 0; j < 2*n; j++ {
				ap += a[i][j] * p[j][k]
			}
			v[i] = complex(p[i][k], -ap/w)
		}
		for _, x := range rm.V {
			var d complex128
			for i := range v {
				d += cmplx.Conj(x[i]) * v[i]
			}
			for i := range v {
				v[i] -= d * x[i]
			}
		}
		var nn float64
		for _, z := range v {
			nn += real(z)*real(z) + imag(z)*imag(z)
		}
		if nn < 0.5 { // norm² of v is 2 before projection
			continue
		}
		for i := range v {
			v[i] /= complex(math.Sqrt(nn), 0)
		}
		rm.V = append(rm.V, v)
		rm.Omega = append(rm.Omega, w)
		if len(rm.V) == n {
			break
		}
	}
	// re-order by ω (Gram–Schmidt kept ascending order already) and build û
	for m, v := range rm.V {
		u := make([]complex128, n)
		var s float64
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				u[i] += complex(chi[i][j], 0) * v[j]
			}
			s += real(u[i])*real(u[i]) + imag(u[i])*imag(u[i])
		}
		for i := range u {
			u[i] /= complex(math.Sqrt(s), 0)
		}
		rm.U = append(rm.U, u)
		for i := 0; i < 2*n; i++ {
			var av complex128
			for j := 0; j < 2*n; j++ {
				av += complex(a[i][j], 0) * v[j]
			}
			rm.ResidualA = math.Max(rm.ResidualA, cmplx.Abs(av-complex(0, rm.Omega[m])*v[i])/rm.Omega[m])
		}
	}
	return rm
}

// GyroFor is the gyrator transconductance that splits the uniform ring's
// pair at κ by Δω = ν·ω(κ): Gy = ν·ω·c_κ/(4 sin κ). Positive ν lowers the
// CCW wave (e^{+iκj}) and raises the CW wave.
func (r FRing) GyroFor(nu, kappa float64) float64 {
	ck := 4*r.C*sq(math.Sin(kappa/2)) + r.Cg
	return nu * r.OmegaK(kappa) * ck / (4 * math.Sin(kappa))
}

// PairSplitExact is the exact uniform-ring splitting ω_CW − ω_CCW at κ for
// transconductance gy: the positive roots of c_κω² ± 4·gy·sin κ·ω − 1/Lg = 0.
func (r FRing) PairSplitExact(gy, kappa float64) (wCW, wCCW float64) {
	ck := 4*r.C*sq(math.Sin(kappa/2)) + r.Cg
	b := 4 * gy * math.Sin(kappa)
	d := math.Sqrt(b*b + 4*ck/r.Lg)
	return (b + d) / (2 * ck), (-b + d) / (2 * ck)
}

// GroupVelocityRing is the circuit line's own group velocity at κ if the N
// sections are laid round a ring of radius a: v_g = (2πa/N)·|dω/dκ|, with
// dω/dκ = −ω·C·sin κ/c_κ.
func (r FRing) GroupVelocityRing(kappa, a float64) float64 {
	ck := 4*r.C*sq(math.Sin(kappa/2)) + r.Cg
	return 2 * math.Pi * a / float64(r.N) * r.OmegaK(kappa) * r.C * math.Abs(math.Sin(kappa)) / ck
}

// SagnacRim converts a non-reciprocity ν = Δω/ω into the rim speed aΩ of a
// ring whose wave has group velocity vg: aΩ = ν·c²/(2·vg), and reports the
// rotor against the burst speed of T1000G (spacetime.AnalyzeRotor). Speeds
// at or above c are returned with Possible = false.
type SagnacRim struct {
	Nu, Vg      float64
	RimSpeed    float64 // aΩ [m/s]
	Omega       float64 // Ω for the ring radius a [rad/s]
	Beta2       float64 // (aΩ/c)²
	BurstMargin float64 // rim speed / T1000G hoop burst speed
	KinPerKg    float64 // ½(aΩ)² [J/kg], paid by the drive
	Possible    bool    // below c and below the burst speed
	BelowC      bool
	TwinRel     float64 // |ν from spacetime.SagnacSplitting·vg/c − ν| / ν
}

// Rim maps ν onto a ring of radius a with group velocity vg at frequency f.
func Rim(nu, vg, a, f float64) SagnacRim {
	s := SagnacRim{Nu: nu, Vg: vg}
	s.RimSpeed = math.Abs(nu) * C0 * C0 / (2 * vg)
	s.Omega = s.RimSpeed / a
	s.Beta2 = sq(s.RimSpeed / C0)
	s.BurstMargin = s.RimSpeed / spacetime.BurstSpeed(spacetime.T1000G)
	s.KinPerKg = 0.5 * s.RimSpeed * s.RimSpeed
	s.BelowC = s.RimSpeed < C0
	s.Possible = s.BelowC && s.BurstMargin < 1
	if s.BelowC {
		// twin: the ring-laser splitting (vacuum) times v_g/c, as Δω/ω.
		df := spacetime.SagnacSplitting(math.Pi*a*a, 2*math.Pi*a, C0/f, s.Omega) * vg / C0
		s.TwinRel = math.Abs(df/f-math.Abs(nu)) / math.Abs(nu)
		if ro, err := spacetime.AnalyzeRotor(s.Omega, a, spacetime.T1000G); err == nil {
			s.BurstMargin = ro.BurstMargin
			s.Beta2 = ro.Beta2
		}
	}
	return s
}

// ---------------------------------------------------------------- three states with rotation

// RotCluster is Part 6's band-edge cluster with rotation.
type RotCluster struct {
	D1, D2       float64
	P1, P2       int
	Nu           float64 // non-reciprocity on the uniform ring's pair at κ = 11π/12
	Gyro         float64 // gyrator transconductance [S]
	Omega        [3]float64
	U            [3][3]complex128 // rows f1, f2, f3; columns the three modes by ω (raw projections)
	UU           [3][3]complex128 // rows e, μ, τ; columns ν1, ν2, ν3: Löwdin-unitarized
	Leak         float64
	Assign, Mass [3]int
	Th12, Th23   float64 // degrees, from the column-normalized |U| (as Part 6)
	Th13         float64
	Ratio        float64
	J            float64 // Jarlskog invariant of UU
	JMax         float64 // c12 s12 c23 s23 c13² s13 (|J| at δ = ±90° for these angles)
	DeltaCP      float64 // degrees in [0, 360)
	Cost         float64 // Part 6 cost (angles and ratio)
	CostCP       float64 // with δ_CP against 212° (+26/−41)
	Modes        RotModes
}

// RotClusterMix computes the cluster for defects (d1 on p1, d2 on p2) and
// non-reciprocity nu, flavours assigned e, μ, τ = f(perm).
func RotClusterMix(d1, d2 float64, p1, p2 int, nu float64, perm [3]int) RotCluster {
	r := DefaultFRing(false)
	r.Def = map[int]float64{p1: d1}
	r.Def[p2] += d2
	k := r.Kappa(r.N/2 - 1)
	r.Gyro = r.GyroFor(nu, k)
	rm := r.RotModes()
	f := r.clusterFlavours()
	n := r.N
	type cand struct {
		i int
		w float64
	}
	var cs []cand
	proj := func(a, i int) complex128 {
		var d complex128
		for j := 0; j < n; j++ {
			d += complex(f[a][j], 0) * rm.U[i][j]
		}
		return d
	}
	for i := range rm.U {
		var s float64
		for a := 0; a < 3; a++ {
			s += sq(cmplx.Abs(proj(a, i)))
		}
		cs = append(cs, cand{i, s})
	}
	sort.Slice(cs, func(a, b int) bool { return cs[a].w > cs[b].w })
	ids := []int{cs[0].i, cs[1].i, cs[2].i}
	sort.Slice(ids, func(a, b int) bool { return rm.Omega[ids[a]] < rm.Omega[ids[b]] })
	cl := RotCluster{D1: d1, D2: d2, P1: p1, P2: p2, Nu: nu, Gyro: r.Gyro, Modes: rm}
	for m, i := range ids {
		cl.Omega[m] = rm.Omega[i]
		var s float64
		for a := 0; a < 3; a++ {
			cl.U[a][m] = proj(a, i)
			s += sq(cmplx.Abs(cl.U[a][m]))
		}
		cl.Leak = math.Max(cl.Leak, 1-s)
	}
	cl.angles(perm, NuFIT6())
	return cl
}

// angles: mass labels as Part 6 (closest pair = 1, 2; ν1 the larger e
// content), angles from the column-normalized |U|, J and δ_CP from the
// Löwdin-unitarized matrix (U·(U†U)^−½, the unitary closest to U).
func (cl *RotCluster) angles(perm [3]int, target PMNS) {
	w := cl.Omega
	pairs := [][2]int{{0, 1}, {1, 2}, {0, 2}}
	best := 0
	for q := range pairs {
		if math.Abs(w[pairs[q][0]]-w[pairs[q][1]]) < math.Abs(w[pairs[best][0]]-w[pairs[best][1]]) {
			best = q
		}
	}
	i1, i2 := pairs[best][0], pairs[best][1]
	i3 := 3 - i1 - i2
	col := func(i int) [3]complex128 {
		v := [3]complex128{cl.U[perm[0]][i], cl.U[perm[1]][i], cl.U[perm[2]][i]}
		n := math.Sqrt(sq(cmplx.Abs(v[0])) + sq(cmplx.Abs(v[1])) + sq(cmplx.Abs(v[2])))
		return [3]complex128{v[0] / complex(n, 0), v[1] / complex(n, 0), v[2] / complex(n, 0)}
	}
	if cmplx.Abs(col(i1)[0]) < cmplx.Abs(col(i2)[0]) {
		i1, i2 = i2, i1
	}
	c1, c2, c3 := col(i1), col(i2), col(i3)
	cl.Assign, cl.Mass = perm, [3]int{i1, i2, i3}
	cl.Th13 = math.Asin(math.Min(1, cmplx.Abs(c3[0]))) * 180 / math.Pi
	cl.Th12 = math.Atan2(cmplx.Abs(c2[0]), cmplx.Abs(c1[0])) * 180 / math.Pi
	cl.Th23 = math.Atan2(cmplx.Abs(c3[1]), cmplx.Abs(c3[2])) * 180 / math.Pi
	cl.Ratio = math.Abs(w[i3]-w[i1]) / math.Abs(w[i2]-w[i1])
	cl.Cost = sq((cl.Th12-target.Th12)/target.S12) + sq((cl.Th23-target.Th23)/target.S23) + sq((cl.Th13-target.Th13)/target.S13) +
		sq(math.Log(cl.Ratio/target.Ratio)/target.SRatio)
	var m [3][3]complex128
	for a := 0; a < 3; a++ {
		for i, c := range [3][3]complex128{c1, c2, c3} {
			m[a][i] = c[a]
		}
	}
	cl.UU = lowdin(m)
	u := cl.UU
	cl.J = imag(u[0][0] * u[1][1] * cmplx.Conj(u[0][1]) * cmplx.Conj(u[1][0]))
	s13 := cmplx.Abs(u[0][2])
	c13 := math.Sqrt(math.Max(0, 1-s13*s13))
	s12 := cmplx.Abs(u[0][1]) / c13
	c12 := math.Sqrt(math.Max(0, 1-s12*s12))
	s23 := cmplx.Abs(u[1][2]) / c13
	c23 := math.Sqrt(math.Max(0, 1-s23*s23))
	cl.JMax = c12 * s12 * c23 * s23 * c13 * c13 * s13
	sd := cl.J / cl.JMax
	cd := (sq(cmplx.Abs(u[1][0])) - s12*s12*c23*c23 - c12*c12*s23*s23*s13*s13) / (2 * s12 * c12 * s23 * c23 * s13)
	d := math.Atan2(sd, cd) * 180 / math.Pi
	if d < 0 {
		d += 360
	}
	cl.DeltaCP = d
	dd := d - target.DeltaCP
	sig := 26.0
	if dd < 0 {
		sig = 41
	}
	cl.CostCP = cl.Cost + sq(dd/sig)
}

// lowdin returns U·(U†U)^−½ for a 3×3 complex matrix, using the real
// symmetric embedding [[Re, −Im], [Im, Re]] of the Hermitian U†U.
func lowdin(u [3][3]complex128) [3][3]complex128 {
	var s [3][3]complex128
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			for a := 0; a < 3; a++ {
				s[i][j] += cmplx.Conj(u[a][i]) * u[a][j]
			}
		}
	}
	m := make([][]float64, 6)
	for i := range m {
		m[i] = make([]float64, 6)
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			m[i][j], m[i+3][j+3] = real(s[i][j]), real(s[i][j])
			m[i][j+3], m[i+3][j] = -imag(s[i][j]), imag(s[i][j])
		}
	}
	f := symFunc(m, func(x float64) float64 { return 1 / math.Sqrt(x) })
	var si [3][3]complex128
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			si[i][j] = complex(f[i][j], f[i+3][j])
		}
	}
	var out [3][3]complex128
	for a := 0; a < 3; a++ {
		for j := 0; j < 3; j++ {
			for i := 0; i < 3; i++ {
				out[a][j] += u[a][i] * si[i][j]
			}
		}
	}
	return out
}

// nelderMead minimizes f from x0 with initial steps; deterministic.
func nelderMead(f func([]float64) float64, x0, step []float64, iters int) ([]float64, float64) {
	n := len(x0)
	xs := make([][]float64, n+1)
	fs := make([]float64, n+1)
	for i := range xs {
		xs[i] = append([]float64(nil), x0...)
		if i > 0 {
			xs[i][i-1] += step[i-1]
		}
		fs[i] = f(xs[i])
	}
	order := func() {
		sort.Sort(simplex{xs, fs})
	}
	for it := 0; it < iters; it++ {
		order()
		c := make([]float64, n)
		for i := 0; i < n; i++ {
			for d := 0; d < n; d++ {
				c[d] += xs[i][d] / float64(n)
			}
		}
		at := func(t float64) []float64 {
			p := make([]float64, n)
			for d := range p {
				p[d] = c[d] + t*(xs[n][d]-c[d])
			}
			return p
		}
		xr := at(-1)
		fr := f(xr)
		switch {
		case fr < fs[0]:
			xe := at(-2)
			if fe := f(xe); fe < fr {
				xs[n], fs[n] = xe, fe
			} else {
				xs[n], fs[n] = xr, fr
			}
		case fr < fs[n-1]:
			xs[n], fs[n] = xr, fr
		default:
			xc := at(0.5)
			if fc := f(xc); fc < fs[n] {
				xs[n], fs[n] = xc, fc
			} else {
				for i := 1; i <= n; i++ {
					for d := range xs[i] {
						xs[i][d] = xs[0][d] + 0.5*(xs[i][d]-xs[0][d])
					}
					fs[i] = f(xs[i])
				}
			}
		}
	}
	order()
	return xs[0], fs[0]
}

type simplex struct {
	x [][]float64
	f []float64
}

func (s simplex) Len() int           { return len(s.f) }
func (s simplex) Less(i, j int) bool { return s.f[i] < s.f[j] }
func (s simplex) Swap(i, j int) {
	s.f[i], s.f[j] = s.f[j], s.f[i]
	s.x[i], s.x[j] = s.x[j], s.x[i]
}

// FitCP searches the two defect sizes (|δ| ≤ 0.5) on sections p1, p2 and
// the non-reciprocity ν for the measured angles, splitting ratio and
// δ_CP = 212° (+26/−41), flavours e, μ, τ = f1, f2, f3, by Nelder–Mead
// from Part 6's defect sizes and each of the given ν starts. With pin > 0
// the δ_CP term is instead ((δ − 212°)/pin)², pinning the phase at the
// measured central value to see what the angles then do. Deterministic.
func FitCP(p1, p2 int, nuStarts []float64, pin float64) RotCluster {
	perm := [3]int{0, 1, 2}
	cost := func(c RotCluster) float64 {
		if pin > 0 {
			return c.Cost + sq((c.DeltaCP-NuFIT6().DeltaCP)/pin)
		}
		return c.CostCP
	}
	eval := func(x []float64) float64 {
		if math.Abs(x[0]) > 0.5 || math.Abs(x[1]) > 0.5 {
			return math.Inf(1)
		}
		return cost(RotClusterMix(x[0], x[1], p1, p2, x[2], perm))
	}
	best, bc := RotCluster{CostCP: math.Inf(1)}, math.Inf(1)
	for _, n0 := range nuStarts {
		for _, d0 := range [][2]float64{{0.2064, 0.1628}, {0.2, -0.2}, {-0.2, 0.2}} {
			x, _ := nelderMead(eval, []float64{d0[0], d0[1], n0}, []float64{0.03, 0.03, 0.3 * n0}, 200)
			if c := RotClusterMix(x[0], x[1], p1, p2, x[2], perm); cost(c) < bc {
				best, bc = c, cost(c)
			}
		}
	}
	return best
}

// FitCPAll runs FitCP on every section pair (p1 < N/2, the reflection maps
// the rest onto these) in parallel, results kept by pair index.
func FitCPAll(nuStarts []float64, pin float64) RotCluster {
	n := DefaultFRing(false).N
	var pairs [][2]int
	for p1 := 0; p1 < n/2; p1++ {
		for p2 := p1 + 1; p2 < n; p2++ {
			pairs = append(pairs, [2]int{p1, p2})
		}
	}
	res := make([]RotCluster, len(pairs))
	sem := make(chan struct{}, runtime.NumCPU())
	var wg sync.WaitGroup
	for q, pp := range pairs {
		wg.Add(1)
		sem <- struct{}{}
		go func(q, p1, p2 int) {
			defer wg.Done()
			defer func() { <-sem }()
			res[q] = FitCP(p1, p2, nuStarts, pin)
		}(q, pp[0], pp[1])
	}
	wg.Wait()
	cost := func(c RotCluster) float64 {
		if pin > 0 {
			return c.Cost + sq((c.DeltaCP-NuFIT6().DeltaCP)/pin)
		}
		return c.CostCP
	}
	best := res[0]
	for _, c := range res[1:] {
		if cost(c) < cost(best) {
			best = c
		}
	}
	return best
}

// ---------------------------------------------------------------- kronspice runs

// SagnacRun is the harmonic splitting of a rotating ring measured in a
// kronspice transient.
type SagnacRun struct {
	Ring                string
	Nu, Gyro            float64
	Kappa, F0           float64 // pair phase per section, unperturbed f [Hz]
	DwExact             float64 // ω_CW − ω_CCW from the quadratic [rad/s]
	DwModes             float64 // from RotModes
	DwWarped            float64 // the exact value at the trapezoidal frequencies
	DwKron              float64 // phase slopes of A and B in kronspice
	DwFirst             float64 // first order 4·Gy·sin κ/c_κ
	ConvAmp             float64 // S → D conversion amplitude (fitted)
	GyroPowerMax        float64 // largest power absorbed by any gyrator / (ω·W)
	GyroPass            float64 // largest power one gyrator passes between sections / (ω·W)
	GyroWork            float64 // |∫ Σ gyrator power dt| / stored energy
	LedgerRel, Mismatch float64
	Periods             float64
}

// Sagnac runs harmonic 1 with non-reciprocity nu, started in S (c = ½).
func Sagnac(crossed bool, nu, periods float64) (SagnacRun, error) {
	r := DefaultFRing(crossed)
	k := r.Kappa(1)
	w := r.OmegaK(k)
	r.Gyro = r.GyroFor(nu, k)
	s := SagnacRun{Ring: map[bool]string{false: "plain", true: "crossed"}[crossed], Nu: nu, Gyro: r.Gyro, Kappa: k, F0: w / (2 * math.Pi), Periods: periods}
	wcw, wccw := r.PairSplitExact(r.Gyro, k)
	s.DwExact = wcw - wccw
	ck := 4*r.C*sq(math.Sin(k/2)) + r.Cg
	s.DwFirst = 4 * r.Gyro * math.Sin(k) / ck
	rm := r.RotModes()
	// the two modes nearest w: CW has the larger projection on e^{−iκj}
	type mw struct {
		w, cw float64
	}
	var near []mw
	for i, om := range rm.Omega {
		if math.Abs(om-w) < 0.05*w {
			var p complex128
			for j := 0; j < r.N; j++ {
				p += rm.U[i][j] * cmplx.Exp(complex(0, k*float64(j)))
			}
			near = append(near, mw{om, cmplx.Abs(p)})
		}
	}
	if len(near) == 2 {
		if near[0].cw < near[1].cw {
			near[0], near[1] = near[1], near[0]
		}
		s.DwModes = near[0].w - near[1].w
	}
	const per = 60
	h := 2 * math.Pi / w / per
	warp := func(x float64) float64 { return 2 / h * math.Atan(x*h/2) }
	s.DwWarped = warp(wcw) - warp(wccw)
	ckt, err := r.Build([]Comp{r.Type("S", k, 0.5)})
	if err != nil {
		return s, err
	}
	wf, err := tranRun(ckt, h, periods*2*math.Pi/w)
	if err != nil {
		return s, err
	}
	fr, err := r.read(wf)
	if err != nil {
		return s, err
	}
	var phA, phB []float64
	var ts, pd []float64
	for t := range fr.Time {
		A, B := r.Amp(fr, t, k, w)
		phA, phB = append(phA, cmplx.Phase(A)), append(phB, cmplx.Phase(B))
		ts = append(ts, fr.Time[t])
		pd = append(pd, Frac(A, B, k, 0.5).D)
	}
	// z ∝ e^{iωt}: A (CW) phase advances at ω_CW, B at ω_CCW.
	sa, _ := fitLine(ts, unwrap(phA))
	sb, _ := fitLine(ts, unwrap(phB))
	s.DwKron = sa - sb
	_, s.ConvAmp, _ = fitBeat(ts, pd, math.Abs(s.DwWarped)/(2*math.Pi))
	s.LedgerRel = fr.LedgerRel
	s.Mismatch = math.Abs(fr.Mismatch) / fr.StoredStart
	// The rotation's power. Each gyrator absorbs V1·i1 + V2·i2 = 0 (power
	// in at one section is passed to the next); the power it passes,
	// V1·i1 = u_j·Gy·u_{j+1}, is compared with ω·W.
	scale := w * fr.StoredStart
	var mx, pass, work float64
	prev := 0.0
	for t := range wf.Time {
		var p float64
		for i, n := range wf.ElementNames {
			if len(n) > 2 && n[:2] == "YR" {
				p += wf.P[i][t]
				mx = math.Max(mx, math.Abs(wf.P[i][t]))
			}
		}
		for j := 0; j < r.N; j++ {
			m := (j + 1) % r.N
			sg := 1.0
			if r.Crossed && j == r.N-1 {
				sg = -1
			}
			pass = math.Max(pass, math.Abs(fr.U[t][j]*r.Gyro*sg*fr.U[t][m]))
		}
		if t > 0 {
			work += 0.5 * (p + prev) * (wf.Time[t] - wf.Time[t-1])
		}
		prev = p
	}
	s.GyroPowerMax = mx / scale
	s.GyroPass = pass / scale
	s.GyroWork = math.Abs(work) / fr.StoredStart
	return s, nil
}

func unwrap(ph []float64) []float64 {
	out := make([]float64, len(ph))
	var acc float64
	for i := range ph {
		if i > 0 {
			d := ph[i] - ph[i-1]
			for d > math.Pi {
				d -= 2 * math.Pi
			}
			for d < -math.Pi {
				d += 2 * math.Pi
			}
			acc += d
		}
		out[i] = acc
	}
	return out
}

// TwoState is the CW → CCW conversion with a defect and rotation together.
type TwoState struct {
	Delta, Nu               float64
	DfDefect, DfRot         float64 // splittings each alone (exact modes) [Hz]
	DfExact                 float64 // both together (exact modes) [Hz]
	DfQuad                  float64 // √(Δd² + Δr²)
	AmpPred, AmpFit, DfFit  float64
	ThetaPred, ThetaFit     float64 // degrees
	MaxResid, LedgerRel     float64
	T, P                    []float64
	AmpPredNoRot, DfNoRotHz float64
}

// TwoStateRot: defect δ on section 0 and the rotation that splits the
// uniform pair by the same amount; started CW. sin²2θ = Δd²/(Δd² + Δr²).
func TwoStateRot(delta, rotOverDefect, beats float64) (TwoState, error) {
	r := DefaultFRing(false)
	k := r.Kappa(1)
	w := r.OmegaK(k)
	ts := TwoState{Delta: delta}
	pair := func(x FRing) float64 {
		rm := x.RotModes()
		var near []float64
		for _, om := range rm.Omega {
			if math.Abs(om-w) < 0.05*w {
				near = append(near, om)
			}
		}
		if len(near) < 2 {
			return math.NaN()
		}
		return math.Abs(near[1]-near[0]) / (2 * math.Pi)
	}
	rd := r
	rd.Def = map[int]float64{0: delta}
	ts.DfDefect = pair(rd)
	ts.Nu = rotOverDefect * ts.DfDefect * 2 * math.Pi / w
	rr := r
	rr.Gyro = r.GyroFor(ts.Nu, k)
	ts.DfRot = pair(rr)
	rb := rd
	rb.Gyro = rr.Gyro
	ts.DfExact = pair(rb)
	ts.DfQuad = math.Hypot(ts.DfDefect, ts.DfRot)
	ts.AmpPred = sq(ts.DfDefect) / (sq(ts.DfDefect) + sq(ts.DfRot))
	ts.ThetaPred = 0.5 * math.Asin(math.Sqrt(ts.AmpPred)) * 180 / math.Pi
	stop := beats / ts.DfExact
	fr, err := rb.Run([]Comp{rb.Type("C", k, 0.5)}, w, stop, 60)
	if err != nil {
		return ts, err
	}
	var tt, pp []float64
	for t := range fr.Time {
		A, B := rb.Amp(fr, t, k, w)
		tt, pp = append(tt, fr.Time[t]), append(pp, Frac(A, B, k, 0.5).CCW)
	}
	ts.DfFit, ts.AmpFit, ts.MaxResid = fitBeat(tt, pp, ts.DfExact)
	ts.ThetaFit = 0.5 * math.Asin(math.Sqrt(math.Min(1, ts.AmpFit))) * 180 / math.Pi
	for t := 0; t < len(tt); t += 20 {
		ts.T, ts.P = append(ts.T, tt[t]), append(ts.P, pp[t])
	}
	ts.LedgerRel = fr.LedgerRel
	return ts, nil
}

// RotThree is a kronspice run of the rotating cluster from one flavour.
type RotThree struct {
	From      int
	TAll      []float64
	PAll      [3][]float64 // measured (energy in f1, f2, f3 / the ring's)
	PredAll   [3][]float64 // full complex mode expansion, same readout
	Pred3     [3][]float64 // 3×3 formula with the unitarized U
	MaxDevAll float64
	MaxDev3   float64
	LedgerRel float64
	GyroRel   float64 // |∫ Σ gyrator power dt| / ring energy at the start
	Periods   float64
}

// RotThreeFlavour runs the rotating cluster from flavour `from` for `slow`
// periods of its slow beat and compares it with the exact complex mode
// expansion (same readout as Part 6, trapezoidal frequencies).
func RotThreeFlavour(cl RotCluster, from int, slow float64) (RotThree, error) {
	r := DefaultFRing(false)
	r.Def = map[int]float64{cl.P1: cl.D1}
	r.Def[cl.P2] += cl.D2
	r.Gyro = cl.Gyro
	n := r.N
	k11 := r.Kappa(n/2 - 1)
	w11, w12 := r.OmegaK(k11), r.OmegaK(math.Pi)
	init := map[int]Comp{0: r.Type("S", k11, 0.5), 1: r.Type("D", k11, 0.5), 2: {Kappa: math.Pi, Omega: w12, A: 1}}[from]
	dw := math.Inf(1)
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			dw = math.Min(dw, math.Abs(cl.Omega[i]-cl.Omega[j]))
		}
	}
	stop := slow * 2 * math.Pi / dw
	const per = 120
	h := 2 * math.Pi / w11 / per
	ckt, err := r.Build([]Comp{init})
	if err != nil {
		return RotThree{}, err
	}
	wf, err := tranRun(ckt, h, stop)
	if err != nil {
		return RotThree{}, err
	}
	fr, err := r.read(wf)
	if err != nil {
		return RotThree{}, err
	}
	tr := RotThree{From: from, LedgerRel: fr.LedgerRel, Periods: stop * w11 / (2 * math.Pi)}
	var work float64
	prev := 0.0
	for t := range wf.Time {
		var p float64
		for i, nm := range wf.ElementNames {
			if len(nm) > 2 && nm[:2] == "YR" {
				p += wf.P[i][t]
			}
		}
		if t > 0 {
			work += 0.5 * (p + prev) * (wf.Time[t] - wf.Time[t-1])
		}
		prev = p
	}
	tr.GyroRel = math.Abs(work) / fr.StoredStart
	wk := func(k float64) float64 { return sq(math.Sin(k/2)) + r.Cg/(4*r.C) }
	energies := func(fr FRun, t int) [3]float64 {
		A, B := r.Amp(fr, t, k11, w11)
		P, _ := r.Amp(fr, t, math.Pi, w12)
		fs := Frac(A, B, k11, 0.5)
		s11 := sq(cmplx.Abs(A)) + sq(cmplx.Abs(B))
		return [3]float64{wk(k11) * fs.S * s11, wk(k11) * fs.D * s11, wk(math.Pi) * sq(cmplx.Abs(P))}
	}
	total := func(fr FRun, t int) float64 {
		var s float64
		for m := 1; m <= n/2; m++ {
			k := r.Kappa(m)
			A, B := r.Amp(fr, t, k, r.OmegaK(k))
			if m == n/2 {
				s += wk(k) * sq(cmplx.Abs(A))
			} else {
				s += wk(k) * (sq(cmplx.Abs(A)) + sq(cmplx.Abs(B)))
			}
		}
		return s
	}
	// exact expansion: y(t) = 2·Re Σ_k (v_k†·y0)·v_k·e^{iω_k t}, ω warped.
	rm := cl.Modes
	u0, _ := r.state([]Comp{init})
	y0 := make([]float64, 2*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			y0[i] += rm.CHalf[i][j] * u0[j]
		}
	}
	cf := make([]complex128, len(rm.V))
	wW := make([]float64, len(rm.V))
	for m, v := range rm.V {
		for i := range y0 {
			cf[m] += cmplx.Conj(v[i]) * complex(y0[i], 0)
		}
		wW[m] = 2 / h * math.Atan(rm.Omega[m]*h/2)
	}
	sl := 1 / math.Sqrt(r.Lg)
	one := FRun{U: [][]float64{make([]float64, n)}, G: [][]float64{make([]float64, n)}}
	y := make([]complex128, 2*n)
	for t := range fr.Time {
		tt := fr.Time[t]
		for i := range y {
			y[i] = 0
		}
		for m, v := range rm.V {
			ph := cf[m] * cmplx.Exp(complex(0, wW[m]*tt))
			for i := range y {
				y[i] += ph * v[i]
			}
		}
		for i := 0; i < n; i++ {
			var u float64
			for j := 0; j < n; j++ {
				u += rm.CHalfInv[i][j] * 2 * real(y[j])
			}
			one.U[0][i] = u
			one.G[0][i] = 2 * real(y[n+i]) * sl
		}
		e, ef := energies(fr, t), energies(one, 0)
		tm, tf := total(fr, t), total(one, 0)
		tr.TAll = append(tr.TAll, tt)
		for b := 0; b < 3; b++ {
			var amp complex128
			// 3×3: flavour a = Assign^{-1}; UU rows are e, μ, τ.
			ra, rb := rowOf(cl.Assign, from), rowOf(cl.Assign, b)
			for i := 0; i < 3; i++ {
				amp += cmplx.Conj(cl.UU[ra][i]) * cl.UU[rb][i] * cmplx.Exp(complex(0, cl.Omega[cl.Mass[i]]*tt))
			}
			p3 := sq(cmplx.Abs(amp))
			tr.PAll[b] = append(tr.PAll[b], e[b]/tm)
			tr.PredAll[b] = append(tr.PredAll[b], ef[b]/tf)
			tr.Pred3[b] = append(tr.Pred3[b], p3)
			tr.MaxDevAll = math.Max(tr.MaxDevAll, math.Abs(e[b]/tm-ef[b]/tf))
			tr.MaxDev3 = math.Max(tr.MaxDev3, math.Abs(e[b]/tm-p3))
		}
	}
	return tr, nil
}

// rowOf returns the UU row (e, μ, τ = 0, 1, 2) of ring flavour f.
func rowOf(assign [3]int, f int) int {
	for i, a := range assign {
		if a == f {
			return i
		}
	}
	return -1
}

// TAsym compares P(e→μ)(t) from run a (started in e) with P(μ→e)(t) from
// run b (started in μ): largest difference, measured and predicted.
func TAsym(a, b RotThree, e, mu int) (meas, pred float64) {
	for i := range a.TAll {
		pm := interp(b.TAll, b.PAll[e], a.TAll[i])
		pp := interp(b.TAll, b.PredAll[e], a.TAll[i])
		meas = math.Max(meas, math.Abs(a.PAll[mu][i]-pm))
		pred = math.Max(pred, math.Abs(a.PredAll[mu][i]-pp))
	}
	return
}

// tranRun is the transient Part 6 uses (UIC, fixed step h).
func tranRun(ckt *circuit.Circuit, h, stop float64) (*output.Waveforms, error) {
	return analysis.Tran(ckt, analysis.TranSpec{Step: h, Stop: stop, MaxStep: h, UIC: true}, analysis.TranOptions{MaxSteps: 5_000_000})
}
