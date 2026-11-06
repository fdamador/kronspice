// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The octave bar bent into a ring (user request, 2026-11-05): a ring magnet
// cut into twelve 30° sections, each magnetized opposite to its neighbours,
// and the field it makes outside the ring.
//
// Two ways to make neighbouring sections "opposite" are built and compared:
//
//   - radial: each section is magnetized along its own centre line, outward
//     and inward in turn (N S N S … on the outer face). This is the standard
//     multipole ring of motors, couplings and encoders.
//   - tangential: each section is magnetized along the ring, in turn
//     clockwise and anticlockwise, so like poles meet at every joint
//     (N|N, S|S). This is the four-segment bar of H33 (S N | N S | …) closed
//     into a ring: its poles are consequent poles at the joints.
//
// Both put twelve poles, six N and six S, round the outside. Each section is
// magnetized uniformly ("parallel" magnetization, how segmented magnets are
// made), so its magnetic charge is all on its surfaces, σ = M·n̂ (Jackson
// §5.9): the outer and inner arcs and the two cut faces; the flat top and
// bottom carry none, since M lies in the plane. H outside is the Coulomb
// field of that charge, and B = μ0H there. The twin is the same volume
// written as point dipoles m = M dV (Jackson eq. 5.56), summed directly.

import (
	"math"
	"runtime"
	"sync"
)

// Ring is a ring magnet of N equal sections: inner radius Ri, outer Ro,
// height H (centred on z = 0), section k spanning [k, k+1]·2π/N with
// magnetization M0·Sign[k] along its centre line (Radial) or along the ring
// (tangential).
type Ring struct {
	Name      string
	Ri, Ro, H float64
	M0        float64
	Radial    bool
	Sign      []float64
	Short     string // short label for figures
}

func twelve(name, short string, radial bool, alternate bool) Ring {
	s := make([]float64, 12)
	for k := range s {
		s[k] = 1
		if alternate && k%2 == 1 {
			s[k] = -1
		}
	}
	// A ring 50 mm across, 30 mm bore, 10 mm tall, of the same declared
	// material as the bar (Br 1.25 T, demagnetization ignored).
	return Ring{Name: name, Short: short, Ri: 15e-3, Ro: 25e-3, H: 10e-3, M0: Br / Mu0, Radial: radial, Sign: s}
}

// RadialTwelve: twelve sections magnetized radially, outward and inward in turn.
func RadialTwelve() Ring {
	return twelve("12 sections, radial, alternating (N S N S … outside)", "radial", true, true)
}

// TangentialTwelve: twelve sections magnetized along the ring, in turn
// clockwise and anticlockwise — like poles meet at every joint.
func TangentialTwelve() Ring {
	return twelve("12 sections, along the ring, alternating (N|N, S|S joints)", "tangential", false, true)
}

// RadialUniform is the control: all twelve sections magnetized outward.
func RadialUniform() Ring {
	return twelve("12 sections, radial, all outward (control)", "all outward", true, false)
}

// N is the number of sections.
func (r Ring) N() int { return len(r.Sign) }

// Span is section k's angular range.
func (r Ring) Span(k int) (phi0, phi1 float64) {
	w := 2 * math.Pi / float64(r.N())
	return float64(k) * w, float64(k+1) * w
}

// Mag is section k's magnetization vector [A/m].
func (r Ring) Mag(k int) [3]float64 {
	p0, p1 := r.Span(k)
	c := 0.5 * (p0 + p1)
	m := r.M0 * r.Sign[k]
	if r.Radial {
		return [3]float64{m * math.Cos(c), m * math.Sin(c), 0}
	}
	return [3]float64{-m * math.Sin(c), m * math.Cos(c), 0}
}

// Inside reports whether p is in the magnet's volume.
func (r Ring) Inside(p [3]float64) bool { return r.Near(p, 0) }

// Near reports whether p is in the magnet or within skin of its surface.
func (r Ring) Near(p [3]float64, skin float64) bool {
	rho := math.Hypot(p[0], p[1])
	return rho >= r.Ri-skin && rho <= r.Ro+skin && math.Abs(p[2]) <= r.H/2+skin
}

// Skin is where a traced line counts as having reached a face: the surface
// charge is a grid of points (≤ 0.5 mm apart), whose individual Coulomb
// fields dominate within about half a spacing and would trap a line there.
const Skin = 0.25e-3

// Charge is a point of magnetic charge Q [A·m] at P.
type Charge struct {
	P [3]float64
	Q float64
}

// gauss returns the n-point Gauss–Legendre nodes and weights on [0, 1]
// (n = 1 is the midpoint rule).
func gauss(n int) (x, w []float64) {
	switch n {
	case 1:
		return []float64{0.5}, []float64{1}
	case 2:
		d := 0.5 / math.Sqrt(3)
		return []float64{0.5 - d, 0.5 + d}, []float64{0.5, 0.5}
	default: // 3
		d := 0.5 * math.Sqrt(0.6)
		return []float64{0.5 - d, 0.5, 0.5 + d}, []float64{5.0 / 18, 8.0 / 18, 5.0 / 18}
	}
}

// Charges is the surface charge σ = M·n̂ of every section, on panels no
// larger than h with an n×n Gauss–Legendre rule in each.
func (r Ring) Charges(h float64, n int) []Charge {
	gx, gw := gauss(n)
	panels := func(length float64) int { return int(math.Ceil(length/h - 1e-9)) }
	var out []Charge
	nz := panels(r.H)
	for k := 0; k < r.N(); k++ {
		M := r.Mag(k)
		p0, p1 := r.Span(k)
		// outer (n̂ = r̂) and inner (n̂ = −r̂) arcs: dS = ρ dφ dz
		for _, face := range []struct{ rho, sgn float64 }{{r.Ro, 1}, {r.Ri, -1}} {
			np := panels(face.rho * (p1 - p0))
			dp, dz := (p1-p0)/float64(np), r.H/float64(nz)
			for i := 0; i < np; i++ {
				for j := 0; j < nz; j++ {
					for a := range gx {
						for b := range gx {
							phi := p0 + (float64(i)+gx[a])*dp
							z := -r.H/2 + (float64(j)+gx[b])*dz
							sig := face.sgn * (M[0]*math.Cos(phi) + M[1]*math.Sin(phi))
							q := sig * face.rho * dp * dz * gw[a] * gw[b]
							out = append(out, Charge{[3]float64{face.rho * math.Cos(phi), face.rho * math.Sin(phi), z}, q})
						}
					}
				}
			}
		}
		// cut faces at φ0 (n̂ = −φ̂) and φ1 (n̂ = +φ̂): dS = dρ dz
		nr := panels(r.Ro - r.Ri)
		dr, dz := (r.Ro-r.Ri)/float64(nr), r.H/float64(nz)
		for _, face := range []struct{ phi, sgn float64 }{{p0, -1}, {p1, 1}} {
			sig := face.sgn * (-M[0]*math.Sin(face.phi) + M[1]*math.Cos(face.phi))
			if sig == 0 {
				continue
			}
			for i := 0; i < nr; i++ {
				for j := 0; j < nz; j++ {
					for a := range gx {
						for b := range gx {
							rho := r.Ri + (float64(i)+gx[a])*dr
							z := -r.H/2 + (float64(j)+gx[b])*dz
							out = append(out, Charge{[3]float64{rho * math.Cos(face.phi), rho * math.Sin(face.phi), z}, sig * dr * dz * gw[a] * gw[b]})
						}
					}
				}
			}
		}
	}
	return out
}

// BCharges is B = μ0H at p (outside the magnet) from point charges:
// H = Σ Q (p − P)/(4π|p − P|³).
func BCharges(cs []Charge, p [3]float64) [3]float64 {
	var b [3]float64
	for _, c := range cs {
		dx, dy, dz := p[0]-c.P[0], p[1]-c.P[1], p[2]-c.P[2]
		r2 := dx*dx + dy*dy + dz*dz
		f := c.Q / (r2 * math.Sqrt(r2))
		b[0] += f * dx
		b[1] += f * dy
		b[2] += f * dz
	}
	k := Mu0 / (4 * math.Pi)
	return [3]float64{k * b[0], k * b[1], k * b[2]}
}

// Dipole is a point dipole m [A·m²] at P.
type Dipole struct{ P, M [3]float64 }

// Dipoles is the twin: the volume cut into cells no larger than h, each a
// point dipole m = M·dV at its centre.
func (r Ring) Dipoles(h float64) []Dipole {
	var out []Dipole
	nr := int(math.Ceil((r.Ro - r.Ri) / h))
	nz := int(math.Ceil(r.H / h))
	dr, dz := (r.Ro-r.Ri)/float64(nr), r.H/float64(nz)
	for k := 0; k < r.N(); k++ {
		M := r.Mag(k)
		p0, p1 := r.Span(k)
		for i := 0; i < nr; i++ {
			rho := r.Ri + (float64(i)+0.5)*dr
			np := int(math.Ceil(rho * (p1 - p0) / h))
			dp := (p1 - p0) / float64(np)
			for j := 0; j < np; j++ {
				phi := p0 + (float64(j)+0.5)*dp
				for l := 0; l < nz; l++ {
					z := -r.H/2 + (float64(l)+0.5)*dz
					dv := rho * dr * dp * dz
					out = append(out, Dipole{[3]float64{rho * math.Cos(phi), rho * math.Sin(phi), z}, [3]float64{M[0] * dv, M[1] * dv, M[2] * dv}})
				}
			}
		}
	}
	return out
}

// BDipoles is B at p from point dipoles: μ0/4π [3(m·r̂)r̂ − m]/r³.
func BDipoles(ds []Dipole, p [3]float64) [3]float64 {
	var b [3]float64
	for _, d := range ds {
		x, y, z := p[0]-d.P[0], p[1]-d.P[1], p[2]-d.P[2]
		r2 := x*x + y*y + z*z
		r := math.Sqrt(r2)
		mr := (d.M[0]*x + d.M[1]*y + d.M[2]*z) / r2
		f := 1 / (r2 * r)
		b[0] += f * (3*mr*x - d.M[0])
		b[1] += f * (3*mr*y - d.M[1])
		b[2] += f * (3*mr*z - d.M[2])
	}
	k := Mu0 / (4 * math.Pi)
	return [3]float64{k * b[0], k * b[1], k * b[2]}
}

// Cyl returns B_r, B_φ, B_z of b at azimuth phi.
func Cyl(b [3]float64, phi float64) (br, bp, bz float64) {
	c, s := math.Cos(phi), math.Sin(phi)
	return b[0]*c + b[1]*s, -b[0]*s + b[1]*c, b[2]
}

// RingScan is B_r, B_φ, B_z round the ring on the circle (rho, z), at n
// angles from 0 to 2π (exclusive).
func RingScan(cs []Charge, rho, z float64, n int) (phi, br, bp, bz []float64) {
	phi, br, bp, bz = make([]float64, n), make([]float64, n), make([]float64, n), make([]float64, n)
	parallel(n, func(i int) {
		f := 2 * math.Pi * float64(i) / float64(n)
		b := BCharges(cs, [3]float64{rho * math.Cos(f), rho * math.Sin(f), z})
		phi[i] = f
		br[i], bp[i], bz[i] = Cyl(b, f)
	})
	return
}

// CyclicSignChanges counts sign changes of a periodic sequence.
func CyclicSignChanges(v []float64) int {
	n := SignChanges(v)
	if len(v) > 1 && v[len(v)-1]*v[0] < 0 {
		n++
	}
	return n
}

// PeakAngles returns the angles of the positive local maxima of a periodic
// sequence (the N poles seen from outside, when v is B_r).
func PeakAngles(phi, v []float64) []float64 {
	n := len(v)
	var out []float64
	for i := range v {
		a, c := v[(i+n-1)%n], v[(i+1)%n]
		if v[i] > 0 && v[i] >= a && v[i] > c {
			out = append(out, phi[i])
		}
	}
	return out
}

// TraceLine follows a field line from p0 along +B (dir = 1) or −B (dir = −1)
// with fixed-length RK4 steps ds, until it reaches the magnet (within Skin), leaves the box
// (ρ > rmax or |z| > zmax), or runs nmax steps.
func (r Ring) TraceLine(cs []Charge, p0 [3]float64, dir, ds float64, nmax int, rmax, zmax float64) [][3]float64 {
	unit := func(p [3]float64) [3]float64 {
		b := BCharges(cs, p)
		n := math.Sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2])
		if n == 0 {
			return [3]float64{}
		}
		return [3]float64{dir * b[0] / n, dir * b[1] / n, dir * b[2] / n}
	}
	add := func(p, v [3]float64, s float64) [3]float64 {
		return [3]float64{p[0] + s*v[0], p[1] + s*v[1], p[2] + s*v[2]}
	}
	line := [][3]float64{p0}
	p := p0
	for i := 0; i < nmax; i++ {
		k1 := unit(p)
		k2 := unit(add(p, k1, ds/2))
		k3 := unit(add(p, k2, ds/2))
		k4 := unit(add(p, k3, ds))
		for c := 0; c < 3; c++ {
			p[c] += ds / 6 * (k1[c] + 2*k2[c] + 2*k3[c] + k4[c])
		}
		line = append(line, p)
		if r.Near(p, Skin) || math.Hypot(p[0], p[1]) > rmax || math.Abs(p[2]) > zmax {
			break
		}
	}
	return line
}

// Seeds places field-line starts on the circle (rho, z) where B_r has the
// sign sgn, k per pole, at equal steps of the flux through the circle, so
// line density follows field strength.
func Seeds(cs []Charge, rho, z float64, sgn float64, k int) [][3]float64 {
	const n = 1440
	phi, br, _, _ := RingScan(cs, rho, z, n)
	// start at a sign change so no interval wraps
	start := 0
	for i := range br {
		if br[i]*br[(i+1)%n] < 0 {
			start = (i + 1) % n
			break
		}
	}
	var out [][3]float64
	for i := 0; i < n; {
		j := (start + i) % n
		if sgn*br[j] <= 0 {
			i++
			continue
		}
		// collect the interval
		var idx []int
		for i < n && sgn*br[(start+i)%n] > 0 {
			idx = append(idx, (start+i)%n)
			i++
		}
		var total float64
		for _, m := range idx {
			total += sgn * br[m]
		}
		var acc float64
		next := 0
		for _, m := range idx {
			acc += sgn * br[m]
			for next < k && acc >= total*(float64(next)+0.5)/float64(k) {
				out = append(out, [3]float64{rho * math.Cos(phi[m]), rho * math.Sin(phi[m]), z})
				next++
			}
		}
	}
	return out
}

// parallel runs f(0..n−1) on all CPUs; results go into caller-owned slots,
// so the outcome does not depend on scheduling.
func parallel(n int, f func(i int)) {
	w := runtime.NumCPU()
	var wg sync.WaitGroup
	ch := make(chan int)
	for g := 0; g < w; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range ch {
				f(i)
			}
		}()
	}
	for i := 0; i < n; i++ {
		ch <- i
	}
	close(ch)
	wg.Wait()
}

// Lines traces a field line from every seed (forward, along +B, from the N
// poles to the S poles).
func (r Ring) Lines(cs []Charge, seeds [][3]float64, ds float64) [][][3]float64 {
	out := make([][][3]float64, len(seeds))
	parallel(len(seeds), func(i int) {
		out[i] = r.TraceLine(cs, seeds[i], 1, ds, 2000, 2.4*r.Ro, 1.6*r.Ro)
	})
	return out
}
