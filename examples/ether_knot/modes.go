// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Vibrations of the ether knot. The Faddeev Lagrangian (Battye &
// Sutcliffe's eq. 1, signature +−−−) is
//
//	L = ∂_μn·∂^μn − ½(∂_μn × ∂_νn)·(∂^μn × ∂^νn) = T − E
//	T = ∫ ( |ṅ|² + Σ_i |ṅ × ∂_i n|² ) d³x
//
// with E the static energy. Small motions x(t) about the relaxed knot obey
// 2K ẍ = −H x, so the frequencies are ω² = eigenvalues of (H, 2K): H the
// Hessian of E on the sphere (the Euclidean Hessian projected on the
// tangent plane, minus (n·∇E) at each node, the sphere's curvature term),
// K the kinetic form. Units: c = 1, lengths in ℓ = √(κ4/κ2), so ω is in
// c/ℓ. Here only motions that keep the knot's axial symmetry are taken
// (n = R_z(φ) m, as in axi.go): for these ṅ = R_z(φ) ṁ and
//
//	T = 2π ∫∫ ( |ṁ|² + |ṁ × m_ρ|² + |ṁ × m_z|² + |ṁ × (ẑ × m)|²/ρ² ) ρ dρ dz,
//
// so K is a 3 × 3 block at each node. Two motions cost no energy and must
// come out at ω = 0: sliding the knot along z, and turning the field about
// the vacuum direction ẑ (the ether's internal symmetry). Small waves of
// this ether travel at c with no lowest frequency, so a vibration of the
// knot can leak away as waves: in a finite box the knot's own resonances
// keep their frequency when the box grows, while the box's modes fall as
// 1/(box size).

import (
	"bufio"
	"fmt"
	"math"
	"os"
)

// Operators holds H and K on the tangent space: two dofs per free node,
// along the unit vectors E1, E2 ⊥ m.
type Operators struct {
	Free   []int // node index of each free node
	E      [][2][3]float64
	Hi, Hj []int
	Hv     []float64
	K      [][2][2]float64
	Lambda []float64 // m·∇E per free node (equilibrium check: ∇E ∥ m)
	Resid  float64   // largest tangent gradient / largest normal one
}

func (a *Axi) free(i, j int) bool { return i < a.NR-1 && j > 0 && j < a.NZ-1 }

// tangent basis ⊥ m
func tangent(m [3]float64) (e1, e2 [3]float64) {
	ref := [3]float64{1, 0, 0}
	if math.Abs(m[0]) > 0.8 {
		ref = [3]float64{0, 1, 0}
	}
	e1 = cross(m, ref)
	l := math.Sqrt(dot(e1, e1))
	e1 = scale(e1, 1/l)
	e2 = cross(m, e1)
	return
}

// Operators builds H (by central differences of the analytic gradient,
// eps, over a 5 × 5 colouring: the gradient at a node depends on nodes at
// most two steps away) and K.
func (a *Axi) Operators(eps float64) *Operators {
	N := a.NR * a.NZ
	o := &Operators{}
	slot := make([]int, N)
	for n := range slot {
		slot[n] = -1
	}
	for j := 0; j < a.NZ; j++ {
		for i := 0; i < a.NR; i++ {
			if a.free(i, j) {
				n := a.at(i, j)
				slot[n] = len(o.Free)
				o.Free = append(o.Free, n)
				e1, e2 := tangent(a.vec(n))
				o.E = append(o.E, [2][3]float64{e1, e2})
			}
		}
	}
	grad := func() [3][]float64 {
		var g [3][]float64
		for c := 0; c < 3; c++ {
			g[c] = make([]float64, N)
		}
		a.Energy(&g)
		return g
	}
	g0 := grad()
	var tmax, nmax float64
	for q, n := range o.Free {
		g := [3]float64{g0[0][n], g0[1][n], g0[2][n]}
		l := dot(g, a.vec(n))
		o.Lambda = append(o.Lambda, l)
		nmax = math.Max(nmax, math.Abs(l))
		for b := 0; b < 2; b++ {
			tmax = math.Max(tmax, math.Abs(dot(g, o.E[q][b])))
		}
	}
	o.Resid = tmax / nmax
	type key struct{ r, c int }
	H := map[key]float64{}
	for ci := 0; ci < 5; ci++ {
		for cj := 0; cj < 5; cj++ {
			for d := 0; d < 2; d++ {
				var moved []int
				for q, n := range o.Free {
					i, j := n%a.NR, n/a.NR
					if i%5 == ci && j%5 == cj {
						moved = append(moved, q)
					}
				}
				if len(moved) == 0 {
					continue
				}
				shift := func(s float64) {
					for _, q := range moved {
						n := o.Free[q]
						for c := 0; c < 3; c++ {
							a.m[c][n] += s * o.E[q][d][c]
						}
					}
				}
				shift(eps)
				gp := grad()
				shift(-2 * eps)
				gm := grad()
				shift(eps)
				for _, q := range moved {
					n := o.Free[q]
					i0, j0 := n%a.NR, n/a.NR
					for dj := -2; dj <= 2; dj++ {
						for di := -2; di <= 2; di++ {
							i, j := i0+di, j0+dj
							if i < 0 || j < 0 || i >= a.NR || j >= a.NZ {
								continue
							}
							p := slot[a.at(i, j)]
							if p < 0 {
								continue
							}
							m := a.at(i, j)
							dg := [3]float64{(gp[0][m] - gm[0][m]) / (2 * eps), (gp[1][m] - gm[1][m]) / (2 * eps), (gp[2][m] - gm[2][m]) / (2 * eps)}
							for b := 0; b < 2; b++ {
								if v := dot(o.E[p][b], dg); v != 0 {
									H[key{2*p + b, 2*q + d}] += v
								}
							}
						}
					}
				}
			}
		}
	}
	// symmetrize, add the curvature term
	for k, v := range H {
		if k.r <= k.c {
			w := 0.5 * (v + H[key{k.c, k.r}])
			o.Hi, o.Hj, o.Hv = append(o.Hi, k.r), append(o.Hj, k.c), append(o.Hv, w)
			if k.r != k.c {
				o.Hi, o.Hj, o.Hv = append(o.Hi, k.c), append(o.Hj, k.r), append(o.Hv, w)
			}
		}
	}
	for q := range o.Free {
		for b := 0; b < 2; b++ {
			o.Hi, o.Hj, o.Hv = append(o.Hi, 2*q+b), append(o.Hj, 2*q+b), append(o.Hv, -o.Lambda[q])
		}
	}
	// K: 2π ρ h² (I + Σ_a (|a|² I − a aᵀ)), a ∈ {m_ρ, m_z, (ẑ×m)/ρ}
	h := a.H
	for q, n := range o.Free {
		i, j := n%a.NR, n/a.NR
		p := a.vec(n)
		r := a.rho(i)
		var pm [3]float64
		if i == 0 {
			pm = [3]float64{-p[0], -p[1], p[2]}
		} else {
			pm = a.vec(a.at(i-1, j))
		}
		vs := [][3]float64{
			scale(sub(a.vec(a.at(i+1, j)), pm), 1/(2*h)),
			scale(sub(a.vec(a.at(i, j+1)), a.vec(a.at(i, j-1))), 1/(2*h)),
			scale(zc(p), 1/r),
		}
		var K3 [3][3]float64
		for x := 0; x < 3; x++ {
			K3[x][x] = 1
		}
		for _, v := range vs {
			vv := dot(v, v)
			for x := 0; x < 3; x++ {
				for y := 0; y < 3; y++ {
					d := 0.0
					if x == y {
						d = vv
					}
					K3[x][y] += d - v[x]*v[y]
				}
			}
		}
		w := 2 * math.Pi * r * h * h
		var K2 [2][2]float64
		for b := 0; b < 2; b++ {
			for c := 0; c < 2; c++ {
				var s float64
				for x := 0; x < 3; x++ {
					for y := 0; y < 3; y++ {
						s += o.E[q][b][x] * K3[x][y] * o.E[q][c][y]
					}
				}
				K2[b][c] = w * s
			}
		}
		o.K = append(o.K, K2)
	}
	return o
}

// Rayleigh is vᵀHv / vᵀ(2K)v for a tangent vector given as a 3-vector per
// node (projected on the basis).
func (o *Operators) Rayleigh(v func(n int) [3]float64) float64 {
	x := make([]float64, 2*len(o.Free))
	for q, n := range o.Free {
		w := v(n)
		x[2*q], x[2*q+1] = dot(w, o.E[q][0]), dot(w, o.E[q][1])
	}
	var num, den float64
	for k := range o.Hv {
		num += x[o.Hi[k]] * o.Hv[k] * x[o.Hj[k]]
	}
	for q := range o.Free {
		for b := 0; b < 2; b++ {
			for c := 0; c < 2; c++ {
				den += 2 * x[2*q+b] * o.K[q][b][c] * x[2*q+c]
			}
		}
	}
	return num / den
}

// ZeroModes returns ω² of the two exact symmetries, by Rayleigh quotient:
// sliding along z (v = ∂m/∂z) and turning the field about ẑ (v = ẑ × m).
func (a *Axi) ZeroModes(o *Operators) (slide, turn float64) {
	h := a.H
	slide = o.Rayleigh(func(n int) [3]float64 {
		return scale(sub(a.vec(n+a.NR), a.vec(n-a.NR)), 1/(2*h))
	})
	turn = o.Rayleigh(func(n int) [3]float64 { return zc(a.vec(n)) })
	return
}

// Write saves the operators for modes.py: H as triplets, K as 2 × 2 blocks.
func (o *Operators) Write(path string) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	fmt.Fprintf(w, "%d %d\n", 2*len(o.Free), len(o.Hv))
	for k := range o.Hv {
		fmt.Fprintf(w, "%d %d %.17g\n", o.Hi[k], o.Hj[k], o.Hv[k])
	}
	for _, b := range o.K {
		fmt.Fprintf(w, "%.17g %.17g %.17g %.17g\n", b[0][0], b[0][1], b[1][0], b[1][1])
	}
	return w.Flush()
}

// Ring strikes the knot with a tangent push v (a 3-vector per node,
// projected on the basis, with the two zero modes removed), lets the
// linearized motion 2K ẍ = −H x run for time T (leapfrog, step dt), and
// returns the overlap o(t) = x(t)·2K·x(0) at each step: the knot "ringing".
// Waves that leave the knot are lost to the record until they come back
// from the box wall, about 2·(box − 3ℓ)/c later.
func (a *Axi) Ring(o *Operators, v func(n int) [3]float64, T, dt float64) (ts, ov []float64) {
	nd := 2 * len(o.Free)
	project := func(f func(n int) [3]float64) []float64 {
		x := make([]float64, nd)
		for q, n := range o.Free {
			w := f(n)
			x[2*q], x[2*q+1] = dot(w, o.E[q][0]), dot(w, o.E[q][1])
		}
		return x
	}
	mdot := func(x, y []float64) float64 {
		var s float64
		for q := range o.Free {
			for b := 0; b < 2; b++ {
				for c := 0; c < 2; c++ {
					s += 2 * x[2*q+b] * o.K[q][b][c] * y[2*q+c]
				}
			}
		}
		return s
	}
	x := project(v)
	h := a.H
	zero := [][]float64{
		project(func(n int) [3]float64 { return scale(sub(a.vec(n+a.NR), a.vec(n-a.NR)), 1/(2*h)) }),
		project(func(n int) [3]float64 { return zc(a.vec(n)) }),
	}
	// Gram–Schmidt the zero modes, then remove them from x (2K inner product)
	for i := range zero {
		for j := 0; j < i; j++ {
			c := mdot(zero[i], zero[j]) / mdot(zero[j], zero[j])
			for k := range zero[i] {
				zero[i][k] -= c * zero[j][k]
			}
		}
	}
	for _, z := range zero {
		c := mdot(x, z) / mdot(z, z)
		for k := range x {
			x[k] -= c * z[k]
		}
	}
	x0 := append([]float64(nil), x...)
	// (2K)^{-1} per node
	inv := make([][2][2]float64, len(o.Free))
	for q, b := range o.K {
		a11, a12, a21, a22 := 2*b[0][0], 2*b[0][1], 2*b[1][0], 2*b[1][1]
		d := a11*a22 - a12*a21
		inv[q] = [2][2]float64{{a22 / d, -a12 / d}, {-a21 / d, a11 / d}}
	}
	acc := func(x []float64) []float64 {
		hx := make([]float64, nd)
		for k := range o.Hv {
			hx[o.Hi[k]] += o.Hv[k] * x[o.Hj[k]]
		}
		ac := make([]float64, nd)
		for q := range o.Free {
			ac[2*q] = -(inv[q][0][0]*hx[2*q] + inv[q][0][1]*hx[2*q+1])
			ac[2*q+1] = -(inv[q][1][0]*hx[2*q] + inv[q][1][1]*hx[2*q+1])
		}
		return ac
	}
	vel := make([]float64, nd)
	ac := acc(x)
	for t := 0.0; t <= T; t += dt {
		ts = append(ts, t)
		ov = append(ov, mdot(x, x0))
		for k := range x {
			vel[k] += 0.5 * dt * ac[k]
			x[k] += dt * vel[k]
		}
		ac = acc(x)
		for k := range x {
			vel[k] += 0.5 * dt * ac[k]
		}
	}
	return
}

// Strikes are the three localized pushes: breathing (the knot's size, the
// direction Derrick's balance fixes), squashing (taller and thinner), and
// twisting (turning the field by an angle that changes along z), each
// under an envelope exp(−r²/8ℓ²) that keeps it on the knot.
func (a *Axi) Strikes() map[string]func(n int) [3]float64 {
	h := a.H
	d := func(n int) (mr, mz [3]float64, rho, z float64) {
		i, j := n%a.NR, n/a.NR
		var pm [3]float64
		if i == 0 {
			p := a.vec(n)
			pm = [3]float64{-p[0], -p[1], p[2]}
		} else {
			pm = a.vec(n - 1)
		}
		mr = scale(sub(a.vec(n+1), pm), 1/(2*h))
		mz = scale(sub(a.vec(n+a.NR), a.vec(n-a.NR)), 1/(2*h))
		return mr, mz, a.rho(i), a.z(j)
	}
	env := func(rho, z float64) float64 { return math.Exp(-(rho*rho + z*z) / 8) }
	return map[string]func(n int) [3]float64{
		"breathing": func(n int) [3]float64 {
			mr, mz, rho, z := d(n)
			e := env(rho, z)
			return [3]float64{e * (rho*mr[0] + z*mz[0]), e * (rho*mr[1] + z*mz[1]), e * (rho*mr[2] + z*mz[2])}
		},
		"squashing": func(n int) [3]float64 {
			mr, mz, rho, z := d(n)
			e := env(rho, z)
			return [3]float64{e * (z*mz[0] - 0.5*rho*mr[0]), e * (z*mz[1] - 0.5*rho*mr[1]), e * (z*mz[2] - 0.5*rho*mr[2])}
		},
		"twisting": func(n int) [3]float64 {
			_, _, rho, z := d(n)
			return scale(zc(a.vec(n)), env(rho, z)*z)
		},
	}
}
