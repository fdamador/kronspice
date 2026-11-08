// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/material"
)

// Model holds the inductances of a pair.
type Model struct {
	P *Pair
	// Pass is the inductance matrix of the sections of one circuit of A
	// (0..Sec−1) and one of B (Sec..2Sec−1), each a single bundle path [H].
	Pass [][]float64
	// Delta is the leakage between the wires of a bundle, added on the
	// diagonal for each wire's own section: the pass matrix uses the
	// bundle radius r_b for the self terms, a wire's own self term uses
	// r_w, so Δ = μ0ℓ/2π·ln(r_b/r_w) per section (makes the full section
	// matrix positive definite).
	Delta      []float64
	SecLen     []float64 // length of each pass section of A then B [m]
	LA, LB, M  float64   // whole windings (Nc circuits), M for path-sense currents [H]
	RDC        float64   // wire resistance of winding A [Ω]
	LenA, LenB float64   // wire length of A and of B (B lies one bundle further out) [m]
	FieldScale float64   // ampere-turns per amp of terminal current (Nc)
}

// NewModel computes the pass matrix and the winding inductances. A
// section of circuit c and one of circuit c′ at the same pass are the same
// bundle path, so the whole winding's section matrix is the pass matrix
// repeated (Pass ⊗ ones(Nc)); the self terms use the bundle radius.
func NewModel(p *Pair) *Model {
	var s Segs
	s.AddPolyline(p.A, 0, p.Sec, p.Radius)
	s.AddPolyline(p.B, p.Sec, p.Sec, p.Radius)
	m := &Model{P: p, Pass: s.SectionMatrix(2 * p.Sec), FieldScale: float64(p.Nc)}
	n2 := float64(p.Nc * p.Nc)
	m.SecLen = make([]float64, 2*p.Sec)
	for i := range s.A {
		m.SecLen[s.Sec[i]] += norm(sub(s.B[i], s.A[i]))
	}
	m.Delta = make([]float64, 2*p.Sec)
	for i, l := range m.SecLen {
		m.Delta[i] = em.Mu0 / (2 * math.Pi) * l * math.Log(p.Radius/WireRadius)
		if i < p.Sec {
			m.LA += float64(p.Nc) * m.Delta[i]
		} else {
			m.LB += float64(p.Nc) * m.Delta[i]
		}
	}
	for i := 0; i < p.Sec; i++ {
		for j := 0; j < p.Sec; j++ {
			m.LA += n2 * m.Pass[i][j]
			m.LB += n2 * m.Pass[p.Sec+i][p.Sec+j]
			m.M += n2 * m.Pass[i][p.Sec+j]
		}
	}
	m.LenA = float64(p.Nc) * polyLen(p.A)
	m.LenB = float64(p.Nc) * polyLen(p.B)
	m.RDC = WireResistance(m.LenA, 0)
	return m
}

// K is the coupling coefficient M/√(LA·LB).
func (m *Model) K() float64 { return m.M / math.Sqrt(m.LA*m.LB) }

// WireResistance of a length of the 0.6 mm copper wire at f: DC value
// raised by the skin effect, R_dc·d/(4δ) when δ is small (Terman 1943, as
// em.Coil.Resistance).
func WireResistance(length, f float64) float64 {
	area := math.Pi * WireRadius * WireRadius
	rdc := material.Copper.Resistivity * length / area
	if f <= 0 {
		return rdc
	}
	delta := math.Sqrt(material.Copper.Resistivity / (math.Pi * f * em.Mu0))
	return rdc * math.Max(1, 2*WireRadius/(4*delta))
}

// FieldA, FieldB are the fields [T] of windings A and B per ampere of
// path-sense terminal current.
func (m *Model) FieldA(p P3) P3 { return FieldOf(m.P.A, m.FieldScale, p) }
func (m *Model) FieldB(p P3) P3 { return FieldOf(m.P.B, m.FieldScale, p) }

// Reciprocity computes M_AB as Σ segMutual(a_i, b_j) and M_BA as
// Σ segMutual(b_j, a_i): the recursive splitting differs between the two
// orders, so they agree only to quadrature error.
func (m *Model) Reciprocity() (mab, mba float64) {
	pa, pb := m.P.A, m.P.B
	var mu sync.Mutex
	var wg sync.WaitGroup
	ch := make(chan int, 64)
	go func() {
		for i := 0; i+1 < len(pa); i++ {
			ch <- i
		}
		close(ch)
	}()
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			x, y := 0.0, 0.0
			for i := range ch {
				for j := 0; j+1 < len(pb); j++ {
					x += segMutual(pa[i], pa[i+1], pb[j], pb[j+1], 0)
					y += segMutual(pb[j], pb[j+1], pa[i], pa[i+1], 0)
				}
			}
			mu.Lock()
			mab += x
			mba += y
			mu.Unlock()
		}()
	}
	wg.Wait()
	n2 := float64(m.P.Nc*m.P.Nc) * em.Mu0 / (4 * math.Pi)
	return mab * n2, mba * n2
}

// FluxCheck compares the mutual inductance of winding A and a test circle
// two ways: Neumann's double line integral, and the flux of A's
// Biot–Savart field through the circle's disc (polar Gauss quadrature).
func (m *Model) FluxCheck(centre, axis P3, r float64) (neumann, flux float64) {
	circ := em.Circle(centre, axis, r, 256)
	loop := make([]P3, len(circ))
	for i := range circ {
		loop[i] = circ[i]
	}
	var s Segs
	s.AddPolyline(m.P.A, 0, 1, m.P.Radius)
	neumann = s.MutualTo(loop, 1, 1)[0] * m.FieldScale
	// disc quadrature: 24 Gauss radii × 96 angles
	u := unit(sub(loop[0], centre))
	w := cross(unit(axis), u)
	nr := 24
	var mu sync.Mutex
	var wg sync.WaitGroup
	for ir := 0; ir < nr; ir++ {
		wg.Add(1)
		go func(ir int) {
			defer wg.Done()
			x, wt := gaussPt(ir, nr)
			rr := r * x
			acc := 0.0
			for ia := 0; ia < 96; ia++ {
				th := 2 * math.Pi * (float64(ia) + 0.5) / 96
				p := add(centre, add(scl(u, rr*math.Cos(th)), scl(w, rr*math.Sin(th))))
				acc += dot(m.FieldA(p), unit(axis))
			}
			mu.Lock()
			flux += acc * (2 * math.Pi / 96) * rr * r * wt
			mu.Unlock()
		}(ir)
	}
	wg.Wait()
	return
}

// gaussPt returns the i-th of n Gauss–Legendre points mapped to [0, 1]
// with its weight (by Newton iteration on P_n).
func gaussPt(i, n int) (x, w float64) {
	z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
	var pp float64
	for it := 0; it < 100; it++ {
		p1, p2 := 1.0, 0.0
		for j := 1; j <= n; j++ {
			p3 := p2
			p2 = p1
			p1 = ((2*float64(j)-1)*z*p2 - (float64(j)-1)*p3) / float64(j)
		}
		pp = float64(n) * (z*p1 - p2) / (z*z - 1)
		dz := p1 / pp
		z -= dz
		if math.Abs(dz) < 1e-15 {
			break
		}
	}
	return 0.5 * (1 - z), 1 / ((1 - z*z) * pp * pp)
}

// EnergyCheck compares ½L_ext·I² of one closed filament path (external
// inductance outside a tube of radius rt: self term 2ℓ(ln(2ℓ/rt) − 1)) with
// the field energy ∫B²/2μ0 dV outside the same tube, on a fine grid of
// spacing h around the coil and a coarse one (4h) out to ±25 cm.
func EnergyCheck(path []P3, rt, h float64) (circuit, field float64) {
	var s Segs
	s.AddPolyline(path, 0, 1, rt)
	// replace the internal-inclusive self terms by the external ones
	l := 0.0
	for i := range s.A {
		li := norm(sub(s.B[i], s.A[i]))
		l += 2*li*(math.Log(2*li/rt)-1) - selfSeg(li, rt)
	}
	lext := s.SectionMatrix(1)[0][0] + l*em.Mu0/(4*math.Pi)
	circuit = 0.5 * lext
	inner := [3]float64{0.068, 0.068, 0.03}
	field = gridEnergy(path, rt, h, inner, nil) + gridEnergy(path, rt, 4*h, [3]float64{0.25, 0.25, 0.25}, &inner)
	return
}

// gridEnergy integrates B²/2μ0 (1 A) over the cells of spacing h inside
// the box ±half, outside the box ±skip (when set) and outside the tube.
func gridEnergy(path []P3, rt, h float64, half [3]float64, skip *[3]float64) float64 {
	nx, ny, nz := int(2*half[0]/h), int(2*half[1]/h), int(2*half[2]/h)
	var mu sync.Mutex
	var wg sync.WaitGroup
	tot := 0.0
	ch := make(chan int, nx)
	for i := 0; i < nx; i++ {
		ch <- i
	}
	close(ch)
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			acc := 0.0
			for i := range ch {
				x := -half[0] + (float64(i)+0.5)*h
				for j := 0; j < ny; j++ {
					y := -half[1] + (float64(j)+0.5)*h
					for k := 0; k < nz; k++ {
						z := -half[2] + (float64(k)+0.5)*h
						if skip != nil && math.Abs(x) < skip[0] && math.Abs(y) < skip[1] && math.Abs(z) < skip[2] {
							continue
						}
						p := P3{x, y, z}
						if distToPath(path, p) < rt {
							continue
						}
						b := FieldOf(path, 1, p)
						acc += dot(b, b)
					}
				}
			}
			mu.Lock()
			tot += acc
			mu.Unlock()
		}()
	}
	wg.Wait()
	return tot * h * h * h / (2 * em.Mu0)
}

func distToPath(path []P3, p P3) float64 {
	best := math.Inf(1)
	for k := 0; k+1 < len(path); k++ {
		a, b := path[k], path[k+1]
		d := sub(b, a)
		t := dot(sub(p, a), d) / dot(d, d)
		t = math.Max(0, math.Min(1, t))
		if r := norm(sub(p, add(a, scl(d, t)))); r < best {
			best = r
		}
	}
	return best
}
