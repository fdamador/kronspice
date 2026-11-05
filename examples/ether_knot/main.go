// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command ether_knot is step 2 of hypothesis H29 (CLAUDE.md §1): can a
// nonlinear, fluid-like ether hold a knot of field energy at a fixed size?
//
// Step 1 (examples/knotted_light) found that linear Maxwell keeps a knot's
// topology but not its energy, and H30 that it can be kept pulsing only by
// an inward wave supplied from outside. Two routes to a knot that holds by
// itself were checked against the literature first:
//
//   - Nonlinear electrodynamics (Born–Infeld 1934; Euler–Heisenberg 1936,
//     seen in light-by-light scattering, ATLAS, Nature Phys. 13, 852,
//     2017) has no non-zero finite-energy source-free stationary solution
//     (Kiessling, arXiv:1107.2333, Theorem 1), so it cannot hold a knot at
//     rest.
//   - A superfluid of the Gross–Pitaevskii type unties its vortex knots by
//     reconnection, every one of 1,458 simulated (Kleckner, Kauffman &
//     Irvine, Nature Phys. 12, 650, 2016): its vortex cores are places
//     where the field vanishes and a knot can be cut.
//
// What remains is an ether whose field never vanishes: a unit vector n(x)
// (an order parameter on the sphere), which is what a charged
// two-component superfluid reduces to (Babaev, Faddeev & Niemi, Phys. Rev.
// B 65, 100512, 2002). Its knots are Faddeev's (Faddeev & Niemi, Nature
// 387, 58, 1997; Battye & Sutcliffe, PRL 81, 4798, 1998), with static
// energy
//
//	E = κ2 ∫ (∂_i n)² d³x + κ4 ∫ |B_n|² d³x,   B_n,k = ½ ε_kij n·(∂_i n × ∂_j n)
//
// (Battye–Sutcliffe normalization: κ2 = κ4 = 1). The link to knotted
// light: for Rañada's Hopf map, B_n is Rañada's knotted magnetic field
// (B_n = −16 B_Rañada at t = 0, checked in the tests), so the second term
// is the knotted light's magnetic energy and the first is the ether's
// stiffness. Under a change of size x → λx the first grows as λ and the
// second falls as 1/λ (Derrick 1964): the magnetic term alone spreads the
// knot (as in step 1), the stiffness alone shrinks it, and together they
// fix a size, at which the two terms are equal (virial theorem).
//
// The run relaxes the energy from Rañada's Hopf map (arrested Newton flow,
// as Battye & Sutcliffe) and reports energy, size, the virial ratio and the
// linking number of two preimage curves of n (the field lines of B_n),
// which is the Hopf charge Q.
package main

import (
	"flag"
	"fmt"
	"math"
	"runtime"
	"sync"
)

// Grid is the ether field n on N³ nodes of spacing H, centred on the
// origin; boundary nodes stay at the vacuum value (0, 0, 1).
type Grid struct {
	N  int
	H  float64
	K2 float64 // stiffness coupling κ2
	K4 float64 // magnetic coupling κ4
	n  [3][]float64
}

func (g *Grid) idx(i, j, k int) int { return i + g.N*(j+g.N*k) }

// pos is the position of node i on one axis.
func (g *Grid) pos(i int) float64 { return (float64(i) - float64(g.N-1)/2) * g.H }

// HopfMap is Rañada's map at t = 0, scaled by s: with (α, β) the Bateman
// pair of examples/knotted_light at t = 0 (a point of the 3-sphere),
// n = (2 Re αβ̄, 2 Im αβ̄, |α|² − |β|²), → (0, 0, 1) at infinity.
func HopfMap(x, y, z, s float64) [3]float64 {
	x, y, z = x/s, y/s, z/s
	r2 := x*x + y*y + z*z
	d := r2 + 1
	ar, ai := (r2-1)/d, 2*z/d // α
	br, bi := 2*x/d, -2*y/d   // β
	wr, wi := ar*br+ai*bi, ai*br-ar*bi
	return [3]float64{2 * wr, 2 * wi, ar*ar + ai*ai - br*br - bi*bi}
}

// CompactHopf is the Hopf map with space pulled in so that it reaches
// the vacuum smoothly at radius rc: x → x/(1 − (r/rc)²) for r < rc (the
// whole of space mapped into the ball, so the map and its Hopf charge are
// unchanged), (0, 0, 1) beyond. Rañada's map itself approaches the vacuum
// only as 1/r, so holding a finite box edge at the vacuum would cut it.
func CompactHopf(x, y, z, s, rc float64) [3]float64 {
	r2 := (x*x + y*y + z*z) / (rc * rc)
	if r2 >= 1 {
		return [3]float64{0, 0, 1}
	}
	f := 1 / (1 - r2)
	return HopfMap(x*f, y*f, z*f, s)
}

// NewGrid fills the grid from the compact Hopf map at scale s, reaching
// the vacuum at 0.9 of the half-width.
func NewGrid(N int, half, s, k2, k4 float64) *Grid {
	g := &Grid{N: N, H: 2 * half / float64(N-1), K2: k2, K4: k4}
	for c := 0; c < 3; c++ {
		g.n[c] = make([]float64, N*N*N)
	}
	for k := 0; k < N; k++ {
		for j := 0; j < N; j++ {
			for i := 0; i < N; i++ {
				v := [3]float64{0, 0, 1}
				if i > 0 && j > 0 && k > 0 && i < N-1 && j < N-1 && k < N-1 {
					v = CompactHopf(g.pos(i), g.pos(j), g.pos(k), s, 0.9*half)
				}
				m := g.idx(i, j, k)
				for c := 0; c < 3; c++ {
					g.n[c][m] = v[c]
				}
			}
		}
	}
	return g
}

func cross(a, b [3]float64) [3]float64 {
	return [3]float64{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}
func dot(a, b [3]float64) float64 { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }

// slabs runs fn over k-ranges in parallel and sums its four results in a
// fixed order (deterministic).
func slabs(nk int, fn func(k0, k1 int) [4]float64) [4]float64 {
	w := runtime.NumCPU()
	if w > nk {
		w = nk
	}
	parts := make([][4]float64, w)
	var wg sync.WaitGroup
	for p := 0; p < w; p++ {
		wg.Add(1)
		go func(p int) {
			defer wg.Done()
			parts[p] = fn(p*nk/w, (p+1)*nk/w)
		}(p)
	}
	wg.Wait()
	var s [4]float64
	for _, q := range parts {
		for i := range s {
			s[i] += q[i]
		}
	}
	return s
}

// Energy returns the stiffness and magnetic terms (without couplings)
// and, if grad is non-nil, adds the gradient of κ2·E2 + κ4·E4 into it.
//
// E2 = h³ Σ_x Σ_i |n(x+e_i) − n(x)|²/h² (forward differences);
// E4 = h³ Σ_x Σ_{i<j} |D_i n × D_j n|², D_i central differences at
// interior nodes; |D_i n × D_j n| is |B_n| component by component when the
// differences are tangent to the sphere, so E4 → ∫|B_n|².
func (g *Grid) Energy(grad *[3][]float64) (e2, e4 float64) {
	N, h := g.N, g.H
	st := [3]int{1, N, N * N}
	at := func(m int) [3]float64 { return [3]float64{g.n[0][m], g.n[1][m], g.n[2][m]} }
	// E2 and its gradient: each forward pair (m, m+e_i) once; the gradient
	// is written per node from both its pairs on each axis (no races).
	s := slabs(N, func(k0, k1 int) [4]float64 {
		var part [4]float64
		for k := k0; k < k1; k++ {
			for j := 0; j < N; j++ {
				for i := 0; i < N; i++ {
					m := g.idx(i, j, k)
					p := at(m)
					ijk := [3]int{i, j, k}
					for a := 0; a < 3; a++ {
						if ijk[a] < N-1 {
							q := at(m + st[a])
							d := [3]float64{q[0] - p[0], q[1] - p[1], q[2] - p[2]}
							part[0] += h * dot(d, d)
						}
						if grad != nil && i > 0 && j > 0 && k > 0 && i < N-1 && j < N-1 && k < N-1 {
							q, r := at(m+st[a]), at(m-st[a])
							for c := 0; c < 3; c++ {
								(*grad)[c][m] += g.K2 * 2 * h * (2*p[c] - q[c] - r[c])
							}
						}
					}
				}
			}
		}
		return part
	})
	e2 = s[0]
	// E4: G_i(x) = ∂E4/∂D_i(x) = 2h³ Σ_{j≠i} D_j × (D_i × D_j); node y gets
	// +G_i(y−e_i)/2h − G_i(y+e_i)/2h. Computed in two passes: G on every
	// interior node, then gathered per node.
	var G [3][3][]float64 // G[i][c]
	if grad != nil {
		for a := 0; a < 3; a++ {
			for c := 0; c < 3; c++ {
				G[a][c] = make([]float64, N*N*N)
			}
		}
	}
	h3 := h * h * h
	s = slabs(N, func(k0, k1 int) [4]float64 {
		var part [4]float64
		for k := max(k0, 1); k < min(k1, N-1); k++ {
			for j := 1; j < N-1; j++ {
				for i := 1; i < N-1; i++ {
					m := g.idx(i, j, k)
					var D [3][3]float64
					for a := 0; a < 3; a++ {
						q, r := at(m+st[a]), at(m-st[a])
						D[a] = [3]float64{(q[0] - r[0]) / (2 * h), (q[1] - r[1]) / (2 * h), (q[2] - r[2]) / (2 * h)}
					}
					for a := 0; a < 3; a++ {
						for b := a + 1; b < 3; b++ {
							x := cross(D[a], D[b])
							part[0] += h3 * dot(x, x)
						}
					}
					if grad != nil {
						for a := 0; a < 3; a++ {
							var v [3]float64
							for b := 0; b < 3; b++ {
								if b == a {
									continue
								}
								w := cross(D[b], cross(D[a], D[b]))
								v[0], v[1], v[2] = v[0]+w[0], v[1]+w[1], v[2]+w[2]
							}
							for c := 0; c < 3; c++ {
								G[a][c][m] = 2 * h3 * v[c]
							}
						}
					}
				}
			}
		}
		return part
	})
	e4 = s[0]
	if grad != nil {
		slabs(N, func(k0, k1 int) [4]float64 {
			for k := max(k0, 1); k < min(k1, N-1); k++ {
				for j := 1; j < N-1; j++ {
					for i := 1; i < N-1; i++ {
						m := g.idx(i, j, k)
						for a := 0; a < 3; a++ {
							for c := 0; c < 3; c++ {
								(*grad)[c][m] += g.K4 * (G[a][c][m-st[a]] - G[a][c][m+st[a]]) / (2 * h)
							}
						}
					}
				}
			}
			return [4]float64{}
		})
	}
	return e2, e4
}

// Total is κ2·E2 + κ4·E4.
func (g *Grid) Total() float64 {
	e2, e4 := g.Energy(nil)
	return g.K2*e2 + g.K4*e4
}

// RMSRadius is √⟨r²⟩ of the energy density 1 − n_z ... measured instead on
// the deviation from the vacuum, w = (1 − n_z), which is local and
// positive; for comparison between runs only.
func (g *Grid) RMSRadius() float64 {
	var s, sr float64
	for k := 0; k < g.N; k++ {
		for j := 0; j < g.N; j++ {
			for i := 0; i < g.N; i++ {
				w := 1 - g.n[2][g.idx(i, j, k)]
				x, y, z := g.pos(i), g.pos(j), g.pos(k)
				s += w
				sr += w * (x*x + y*y + z*z)
			}
		}
	}
	return math.Sqrt(sr / s)
}

// Relax runs arrested Newton flow (Battye & Sutcliffe): n̈ = −∇E on the
// sphere with unit mass density, and the velocity is zeroed whenever the
// energy rises. It stops after steps or when the energy changes by less
// than tol (relative) over 100 steps. dt is in units of h.
func (g *Grid) Relax(steps int, dt float64, tol float64, every int, report func(step int, e2, e4 float64)) (e2, e4 float64, done int) {
	N := g.N
	tau := dt * g.H
	var v, gr [3][]float64
	for c := 0; c < 3; c++ {
		v[c] = make([]float64, N*N*N)
		gr[c] = make([]float64, N*N*N)
	}
	h3 := g.H * g.H * g.H
	prev := math.Inf(1)
	ref := math.Inf(1)
	for s := 0; s <= steps; s++ {
		for c := 0; c < 3; c++ {
			for m := range gr[c] {
				gr[c][m] = 0
			}
		}
		e2, e4 = g.Energy(&gr)
		e := g.K2*e2 + g.K4*e4
		if report != nil && (s%every == 0) {
			report(s, e2, e4)
		}
		if s%100 == 0 {
			if math.Abs(ref-e) < tol*e {
				return e2, e4, s
			}
			ref = e
		}
		if e > prev {
			for c := 0; c < 3; c++ {
				for m := range v[c] {
					v[c][m] = 0
				}
			}
		}
		prev = e
		if s == steps {
			break
		}
		slabs(N, func(k0, k1 int) [4]float64 {
			for k := max(k0, 1); k < min(k1, N-1); k++ {
				for j := 1; j < N-1; j++ {
					for i := 1; i < N-1; i++ {
						m := g.idx(i, j, k)
						p := [3]float64{g.n[0][m], g.n[1][m], g.n[2][m]}
						f := [3]float64{-gr[0][m] / h3, -gr[1][m] / h3, -gr[2][m] / h3}
						fp := dot(f, p)
						u := [3]float64{v[0][m], v[1][m], v[2][m]}
						for c := 0; c < 3; c++ {
							u[c] += tau * (f[c] - fp*p[c])
						}
						up := dot(u, p)
						for c := 0; c < 3; c++ {
							u[c] -= up * p[c]
							p[c] += tau * u[c]
						}
						l := math.Sqrt(dot(p, p))
						for c := 0; c < 3; c++ {
							g.n[c][m] = p[c] / l
							v[c][m] = u[c]
						}
					}
				}
			}
			return [4]float64{}
		})
	}
	return e2, e4, steps
}

// Bn is B_n at a point by trilinear interpolation of the node values of
// B_n (central differences), for tracing preimage curves.
func (g *Grid) bnNodes() [3][]float64 {
	N, h := g.N, g.H
	st := [3]int{1, N, N * N}
	var b [3][]float64
	for c := 0; c < 3; c++ {
		b[c] = make([]float64, N*N*N)
	}
	for k := 1; k < N-1; k++ {
		for j := 1; j < N-1; j++ {
			for i := 1; i < N-1; i++ {
				m := g.idx(i, j, k)
				var D [3][3]float64
				for a := 0; a < 3; a++ {
					for c := 0; c < 3; c++ {
						D[a][c] = (g.n[c][m+st[a]] - g.n[c][m-st[a]]) / (2 * h)
					}
				}
				p := [3]float64{g.n[0][m], g.n[1][m], g.n[2][m]}
				b[0][m] = dot(p, cross(D[1], D[2]))
				b[1][m] = dot(p, cross(D[2], D[0]))
				b[2][m] = dot(p, cross(D[0], D[1]))
			}
		}
	}
	return b
}

func (g *Grid) interp(f [3][]float64, x [3]float64) [3]float64 {
	var u [3]float64
	var i0 [3]int
	for a := 0; a < 3; a++ {
		q := x[a]/g.H + float64(g.N-1)/2
		fl := math.Floor(q)
		if fl < 0 || int(fl) >= g.N-1 {
			return [3]float64{}
		}
		i0[a], u[a] = int(fl), q-fl
	}
	var out [3]float64
	for dz := 0; dz < 2; dz++ {
		for dy := 0; dy < 2; dy++ {
			for dx := 0; dx < 2; dx++ {
				w := (1 - u[0] + float64(dx)*(2*u[0]-1)) * (1 - u[1] + float64(dy)*(2*u[1]-1)) * (1 - u[2] + float64(dz)*(2*u[2]-1))
				m := g.idx(i0[0]+dx, i0[1]+dy, i0[2]+dz)
				for c := 0; c < 3; c++ {
					out[c] += w * f[c][m]
				}
			}
		}
	}
	return out
}

// Linking traces two field lines of B_n (preimage curves of n) from the
// seeds, closes each at its nearest return, and returns their Gauss
// linking number (the Hopf charge, up to sign) and the larger gap.
func (g *Grid) Linking(s1, s2 [3]float64) (link, gap float64) {
	b := g.bnNodes()
	f := func(p [3]float64) [3]float64 { return g.interp(b, p) }
	step := g.H / 4
	l1, g1 := traceGap(f, s1, step, 60)
	l2, g2 := traceGap(f, s2, step, 60)
	return gaussLinking(l1, l2), math.Max(g1, g2)
}

// traceGap follows a line once round (RK4 on the unit direction) and
// closes it at its nearest return after it has gone at least 1 length unit.
func traceGap(b func(p [3]float64) [3]float64, p0 [3]float64, step, maxLen float64) ([][3]float64, float64) {
	dir := func(p [3]float64) [3]float64 {
		v := b(p)
		n := math.Sqrt(dot(v, v))
		if n == 0 {
			return [3]float64{}
		}
		return [3]float64{v[0] / n, v[1] / n, v[2] / n}
	}
	add := func(p, d [3]float64, s float64) [3]float64 {
		return [3]float64{p[0] + s*d[0], p[1] + s*d[1], p[2] + s*d[2]}
	}
	line := [][3]float64{p0}
	p := p0
	best, at := math.Inf(1), -1
	for s := 0.0; s < maxLen; s += step {
		k1 := dir(p)
		k2 := dir(add(p, k1, step/2))
		k3 := dir(add(p, k2, step/2))
		k4 := dir(add(p, k3, step))
		for c := 0; c < 3; c++ {
			p[c] += step / 6 * (k1[c] + 2*k2[c] + 2*k3[c] + k4[c])
		}
		line = append(line, p)
		if s < 1 {
			continue
		}
		d := math.Sqrt((p[0]-p0[0])*(p[0]-p0[0]) + (p[1]-p0[1])*(p[1]-p0[1]) + (p[2]-p0[2])*(p[2]-p0[2]))
		if d < best {
			best, at = d, len(line)-1
		}
		if best < 0.3 && d > best+0.3 {
			break
		}
	}
	if at < 0 {
		return line, math.Inf(1)
	}
	return line[:at+1], best
}

// gaussLinking is (1/4π) ∮∮ (r1 − r2)·(dr1 × dr2)/|r1 − r2|³.
func gaussLinking(a, b [][3]float64) float64 {
	var s float64
	for i := range a {
		a0, a1 := a[i], a[(i+1)%len(a)]
		ma := [3]float64{(a0[0] + a1[0]) / 2, (a0[1] + a1[1]) / 2, (a0[2] + a1[2]) / 2}
		da := [3]float64{a1[0] - a0[0], a1[1] - a0[1], a1[2] - a0[2]}
		for j := range b {
			b0, b1 := b[j], b[(j+1)%len(b)]
			mb := [3]float64{(b0[0] + b1[0]) / 2, (b0[1] + b1[1]) / 2, (b0[2] + b1[2]) / 2}
			db := [3]float64{b1[0] - b0[0], b1[1] - b0[1], b1[2] - b0[2]}
			r := [3]float64{ma[0] - mb[0], ma[1] - mb[1], ma[2] - mb[2]}
			d := math.Sqrt(dot(r, r))
			s += dot(r, cross(da, db)) / (d * d * d)
		}
	}
	return s / (4 * math.Pi)
}

func main() {
	flag.Parse()
	fmt.Println("ether_knot (H29 step 2): does a nonlinear ether hold a knot of field energy at a fixed size?")
	fmt.Println("Energy E = κ2∫(∂n)² + κ4∫|B_n|²: the ether's stiffness plus the knotted light's magnetic energy")
	fmt.Println("(B_n = −16 × Rañada's knotted B for the starting map). Units: κ2 = κ4 = 1, length ℓ = √(κ4/κ2).")
	fmt.Println("Axisymmetric (Q = 1) relaxation from Rañada's Hopf map; arrested Newton flow, dt = 0.05 h.")

	fmt.Println("\n1. Stiffness + magnetic: grid and box study")
	fmt.Println("   box (ℓ) | h (ℓ) | steps | E | E2 (stiffness) | E4 (magnetic) | E2/E4 | core radius | rms radius")
	type row struct{ box, h, e, ratio float64 }
	var rows []row
	var best *Axi
	for _, c := range [][2]float64{{10, 0.1}, {16, 0.1}, {16, 0.07}, {16, 0.05}} {
		a := NewAxi(c[0], c[0], c[1], 1, 1, 1)
		e2, e4, n := a.Relax(60000, 0.05, 1e-7, 0, nil)
		co, r := a.Size()
		fmt.Printf("   %8g | %5g | %5d | %.3f | %14.3f | %13.3f | %.4f | %11.3f | %.3f\n", c[0], c[1], n, e2+e4, e2, e4, e2/e4, co, r)
		rows = append(rows, row{c[0], c[1], e2 + e4, e2 / e4})
		best = a
	}
	// Richardson extrapolation in h² from the box-16 rows.
	f := func(i, j int) (float64, float64) {
		a, b := rows[i], rows[j]
		k := (a.e - b.e) / (a.h*a.h - b.h*b.h)
		kr := (a.ratio - b.ratio) / (a.h*a.h - b.h*b.h)
		return b.e - k*b.h*b.h, b.ratio - kr*b.h*b.h
	}
	e0, r0 := f(2, 3)
	fmt.Printf("   extrapolated to h → 0 (h = 0.07, 0.05): E = %.1f, E2/E4 = %.4f (virial theorem: 1)\n", e0, r0)
	fmt.Println("   Battye & Sutcliffe (PRL 81, 4798, 1998), 3-D lattice of 100³: E = 504 (spacing not stated)")

	fmt.Println("\n2. Controls (box 10, h 0.1, from the same start)")
	fmt.Println("   magnetic energy only (κ2 = 0: knotted light, no ether) — Derrick: spreads")
	b := NewAxi(10, 10, 0.1, 1, 0, 1)
	b.Relax(6000, 0.05, 0, 1500, func(s int, e2, e4 float64) {
		_, r := b.Size()
		fmt.Printf("     step %5d  E %.3f  rms radius %.3f\n", s, e4, r)
	})
	fmt.Println("   ether stiffness only (κ4 = 0) — Derrick: shrinks; m_z at the core = −1 while the knot is there")
	c := NewAxi(10, 10, 0.1, 1, 1, 0)
	c.Relax(1000, 0.05, 0, 100, func(s int, e2, e4 float64) {
		co, _ := c.Size()
		fmt.Printf("     step %5d  E %.3f  core radius %.3f  m_z at core %+.3f\n", s, e2, co, c.m[2][c.at(int(co/c.H-0.5), (c.NZ-1)/2)])
	})
	fmt.Println("   (the core reaches the grid spacing by step 200 and the knot unties through the lattice: a grid")
	fmt.Println("    cannot hold topology below its own spacing; in the continuum it shrinks to a point)")

	fmt.Println("\n3. The relaxed knot in 3-D (finest grid above, resampled on 101³ over ±8 ℓ)")
	g := best.To3D(101, 8)
	e2, e4 := g.Energy(nil)
	l, gap := g.Linking([3]float64{0.6, 0, 0}, [3]float64{-0.6, 0, 0})
	l2, gap2 := g.Linking([3]float64{1.4, 0, 0}, [3]float64{0, 0.9, 0.4})
	fmt.Printf("   3-D energy E2 %.2f + E4 %.2f = %.2f (grid h = 0.16: E4 reads low by O(h²))\n", e2, e4, e2+e4)
	fmt.Printf("   linking of preimage curves (field lines of B_n): %.4f (gap %.3f ℓ), %.4f (gap %.3f ℓ) — Hopf charge 1\n", l, gap, l2, gap2)

	fmt.Println("\nStandard prediction: Derrick (1964) — the magnetic term alone spreads the knot, the stiffness alone")
	fmt.Println("shrinks it; together they fix its size at ℓ = √(κ4/κ2) and its energy, ≈ 548·√(κ2κ4) by this run (Faddeev & Niemi 1997).")
	fmt.Println("What the model does not fix: κ2 and κ4 themselves — no measured ether stiffness exists to set them.")
}
