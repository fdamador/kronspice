// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"sort"
)

// Where the void's charge could sit (all normalised to one charge e):
//
//	Deficit: where the ether is missing, ∝ 1 − |ψ|²;
//	Wall:    on the ether's surface (Tewari's interface), ∝ |∇|ψ|²|;
//	Light:   riding with the trapped light, ∝ Σ_q χ_q².
type Placement int

const (
	Deficit Placement = iota
	Wall
	Light
)

func (p Placement) String() string {
	return [...]string{"where the ether is missing", "on the wall (Tewari's interface)", "with the light"}[p]
}

// Charge returns the charge density per cell (∫ρ dV = 1) for a placement,
// counting only cells within rmax of the origin (rmax ≤ 0: all).
func (s *State) Charge(p Placement, rmax float64) []float64 {
	g := s.G
	s.updateDens()
	w := make([]float64, len(s.dens))
	tot := 0.0
	for i := 0; i < g.NR; i++ {
		for j := 0; j < g.NZ; j++ {
			k := g.at(i, j)
			if rmax > 0 && math.Hypot(g.Rho[i], g.Z[j]) > rmax {
				continue
			}
			var v float64
			switch p {
			case Deficit:
				v = 1 - s.dens[k]
			case Wall:
				var gr, gz float64
				if i > 0 && i+1 < g.NR {
					gr = (s.dens[k+g.NZ] - s.dens[k-g.NZ]) / (2 * g.H)
				}
				if j > 0 && j+1 < g.NZ {
					gz = (s.dens[k+1] - s.dens[k-1]) / (2 * g.H)
				}
				v = math.Hypot(gr, gz)
			case Light:
				for _, q := range s.Occ {
					v += s.Chi[q][k] * s.Chi[q][k]
				}
			}
			w[k] = v
			tot += v * g.DV(i)
		}
	}
	for k := range w {
		w[k] /= tot
	}
	return w
}

// R2 is ⟨r²⟩ of a normalised density about the origin.
func (s *State) R2(rho []float64) float64 {
	g := s.G
	m := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			m += rho[g.at(i, j)] * dv * (g.Rho[i]*g.Rho[i] + g.Z[j]*g.Z[j])
		}
	}
	return m
}

// Enclosed is the share of the charge within each radius r (the Gauss's-law
// charge a probe at r would see).
func (s *State) Enclosed(rho []float64, rs []float64) []float64 {
	g := s.G
	type cell struct{ r, q float64 }
	var cs []cell
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			cs = append(cs, cell{math.Hypot(g.Rho[i], g.Z[j]), rho[g.at(i, j)] * dv})
		}
	}
	sort.Slice(cs, func(a, b int) bool { return cs[a].r < cs[b].r })
	out := make([]float64, len(rs))
	acc, k := 0.0, 0
	for n, r := range rs {
		for k < len(cs) && cs[k].r <= r {
			acc += cs[k].q
			k++
		}
		out[n] = acc
	}
	return out
}

// FormFactor is the orientation-averaged form factor ∫ρ sin(qr)/(qr) dV.
func (s *State) FormFactor(rho []float64, q float64) float64 {
	g := s.G
	f := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			x := q * math.Hypot(g.Rho[i], g.Z[j])
			sx := 1.0
			if x > 1e-9 {
				sx = math.Sin(x) / x
			}
			f += rho[g.at(i, j)] * dv * sx
		}
	}
	return f
}

// LocalDeficit is ∫(1 − |ψ|²) dV within radius R (the ether missing near
// the particle) and over the whole box.
func (s *State) LocalDeficit(R float64) (local, total float64) {
	g := s.G
	s.updateDens()
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			d := (1 - s.dens[g.at(i, j)]) * dv
			total += d
			if math.Hypot(g.Rho[i], g.Z[j]) <= R {
				local += d
			}
		}
	}
	return
}

// Quadrupole is ∫(1 − |ψ|²)(2z² − ρ²) dV within R: the void's elongation.
func (s *State) Quadrupole(R float64) float64 {
	g := s.G
	s.updateDens()
	q := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			if math.Hypot(g.Rho[i], g.Z[j]) <= R {
				q += (1 - s.dens[g.at(i, j)]) * dv * (2*g.Z[j]*g.Z[j] - g.Rho[i]*g.Rho[i])
			}
		}
	}
	return q
}

// EtherOutside is the ether energy beyond radius R from the point (0, zc):
// what has left as sound.
func (s *State) EtherOutside(R, zc float64) float64 {
	g := s.G
	e := 0.0
	for i := 0; i < g.NR; i++ {
		dv := g.DV(i)
		for j := 0; j < g.NZ; j++ {
			if math.Hypot(g.Rho[i], g.Z[j]-zc) <= R {
				continue
			}
			k := g.at(i, j)
			p := s.Psi[k]
			var gr, gz complex128
			if i+1 < g.NR {
				gr = (s.Psi[k+g.NZ] - p) / complex(g.H, 0)
			}
			if j+1 < g.NZ {
				gz = (s.Psi[k+1] - p) / complex(g.H, 0)
			}
			d := sq(p)
			e += (0.5*(sq(gr)+sq(gz)) + 0.25*(d-1)*(d-1)) * dv
		}
	}
	return e
}

// Profile samples a relaxed state's ether amplitude |ψ|(ρ, z) by bilinear
// interpolation; outside the grid it is 1 (undisturbed ether).
func (s *State) Profile() func(r, z float64) float64 {
	g := s.G
	a := make([]float64, len(s.Psi))
	for k, p := range s.Psi {
		a[k] = math.Sqrt(sq(p))
	}
	return func(r, z float64) float64 {
		x := r/g.H - 0.5
		y := (z+g.LZ)/g.H - 0.5
		if x < 0 {
			x = 0
		}
		i, j := int(math.Floor(x)), int(math.Floor(y))
		if i+1 >= g.NR || j < 0 || j+1 >= g.NZ {
			return 1
		}
		fx, fy := x-float64(i), y-float64(j)
		return (1-fx)*(1-fy)*a[g.at(i, j)] + fx*(1-fy)*a[g.at(i+1, j)] + (1-fx)*fy*a[g.at(i, j+1)] + fx*fy*a[g.at(i+1, j+1)]
	}
}

// Pair places two copies of a relaxed bubble on the axis at z = ±d.
func Pair(f func(r, z float64) float64, d float64) func(r, z float64) float64 {
	return func(r, z float64) float64 { return f(r, z-d) * f(r, z+d) }
}

// Scaled stretches a relaxed bubble by sr across and sz along the axis.
func Scaled(f func(r, z float64) float64, sr, sz float64) func(r, z float64) float64 {
	return func(r, z float64) float64 { return f(r/sr, z/sz) }
}

// smooth averages a signal over a window of w samples.
func smooth(x []float64, w int) []float64 {
	y := make([]float64, len(x))
	for i := range x {
		a, n := 0.0, 0
		for j := i - w/2; j <= i+w/2; j++ {
			if j >= 0 && j < len(x) {
				a += x[j]
				n++
			}
		}
		y[i] = a / float64(n)
	}
	return y
}

// turning returns the times of the turning points of a signal smoothed
// over w samples, after the first skip samples.
func turning(x []float64, dt float64, skip, w int) []float64 {
	y := smooth(x, w)
	var tp []float64
	for i := skip + 1; i+1 < len(y)-w/2; i++ {
		if (y[i]-y[i-1])*(y[i+1]-y[i]) < 0 {
			tp = append(tp, float64(i)*dt)
		}
	}
	return tp
}

// RingFrequency reads the angular frequency of a damped ringing from the
// spacing of its turning points (smoothed over w samples, which removes the
// sound's ripple): ω = π / median half-period. It also returns how many
// turning points it used.
func RingFrequency(x []float64, dt float64, skip, w int) (float64, int) {
	tp := turning(x, dt, skip, w)
	if len(tp) < 2 {
		return 0, len(tp)
	}
	var gaps []float64
	for i := 1; i < len(tp); i++ {
		gaps = append(gaps, tp[i]-tp[i-1])
	}
	sort.Float64s(gaps)
	return math.Pi / gaps[len(gaps)/2], len(tp)
}

// FirstTurn is the time of the first turning point after the start (a
// disturbed start is itself a turning point, so this is half a swing).
func FirstTurn(x []float64, dt float64, w int) float64 {
	tp := turning(x, dt, w, w)
	if len(tp) == 0 {
		return math.Inf(1)
	}
	return tp[0]
}
