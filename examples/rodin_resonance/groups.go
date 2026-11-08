// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/material"
)

// Output coil groups (Part 3). Every group is a round coil of Turns turns
// of 0.5 mm copper (declared), radius Radius, centred at angle θ round the
// axis at radius Rho and height Z, its axis turned by β from the azimuthal
// direction ê_φ (β = 0: the coil links the toroidal flux inside the tube)
// towards the vertical ẑ (β = 90°: it links flux through the hole's
// direction). The groups sit at θ0, θ0 + Δθ, …; each is tuned to the drive
// frequency by a series capacitor (declared) and feeds the load(s).
const (
	OutWireRadius = 0.25e-3
	Clearance     = 2e-3 // m between a group coil and any bundle surface
	GroupPts      = 64
)

// GroupShape is the shape and placement shared by the groups.
type GroupShape struct {
	Rho, Z, Beta, Radius float64 // m, m, rad, m
	Turns                int
}

// Loop returns the coil's polygon at angle theta.
func (s GroupShape) Loop(theta float64) []P3 {
	er := P3{math.Cos(theta), math.Sin(theta), 0}
	ephi := P3{-math.Sin(theta), math.Cos(theta), 0}
	ez := P3{0, 0, 1}
	axis := add(scl(ephi, math.Cos(s.Beta)), scl(ez, math.Sin(s.Beta)))
	centre := add(scl(er, s.Rho), scl(ez, s.Z))
	// in-plane basis: ê_ρ and axis × ê_ρ
	u := er
	w := cross(axis, u)
	pts := make([]P3, GroupPts)
	for k := range pts {
		a := 2 * math.Pi * float64(k) / GroupPts
		pts[k] = add(centre, add(scl(u, s.Radius*math.Cos(a)), scl(w, s.Radius*math.Sin(a))))
	}
	return pts
}

// SelfL is the coil's inductance: N²·μ0R(ln(8R/r_b) − 7/4), r_b the radius
// of its round bundle.
func (s GroupShape) SelfL() float64 {
	rb := OutWireRadius * math.Sqrt(float64(s.Turns)/Fill)
	return float64(s.Turns*s.Turns) * em.Mu0 * s.Radius * (math.Log(8*s.Radius/rb) - 1.75)
}

// Resistance of its wire at f (skin effect as WireResistance).
func (s GroupShape) Resistance(f float64) float64 {
	l := float64(s.Turns) * 2 * math.Pi * s.Radius
	area := math.Pi * OutWireRadius * OutWireRadius
	rdc := material.Copper.Resistivity * l / area
	delta := math.Sqrt(material.Copper.Resistivity / (math.Pi * f * em.Mu0))
	return rdc * math.Max(1, 2*OutWireRadius/(4*delta))
}

// Clear returns how far (m) the coil at theta stays from the bundles
// beyond the required clearance (negative = violation), and from the
// former's inside.
func (m *Model) Clear(s GroupShape, theta float64) float64 {
	worst := math.Inf(1)
	for _, p := range s.Loop(theta) {
		d := math.Min(distToPath(m.P.A, p), distToPath(m.P.B, p)) - m.P.Radius - Clearance
		dt := math.Hypot(math.Hypot(p[0], p[1])-TorusR0, p[2]) // distance from the tube's centre circle
		if dt < TorusA {
			d = math.Min(d, dt-TorusA)
		}
		worst = math.Min(worst, d)
	}
	return worst
}

// PassMutual is the mutual inductance between the coil (all turns) and
// each pass section of one circuit of A then of B [H].
func (m *Model) PassMutual(s GroupShape, theta float64) []float64 {
	var sg Segs
	sg.AddPolyline(m.P.A, 0, m.P.Sec, m.P.Radius)
	sg.AddPolyline(m.P.B, m.P.Sec, m.P.Sec, m.P.Radius)
	return sg.MutualTo(s.Loop(theta), float64(s.Turns), 2*m.P.Sec)
}

// SectionMutual expands pass mutuals onto the sections of a net (every
// circuit of a bundle shares its pass's mutual).
func (n *Net) SectionMutual(m *Model, pass []float64) []float64 {
	out := make([]float64, len(n.L))
	for g := range out {
		w := g / n.NS
		q := g % m.P.Sec
		out[g] = pass[w*m.P.Sec+q]
	}
	return out
}

// GroupMutual is the mutual inductance of two group coils (Neumann, all
// turns).
func GroupMutual(s GroupShape, t1, t2 float64) float64 {
	return float64(s.Turns*s.Turns) * em.Neumann(em.Loop{P: toEm(s.Loop(t1))}, em.Loop{P: toEm(s.Loop(t2))})
}

func toEm(p []P3) [][3]float64 {
	out := make([][3]float64, len(p))
	copy(out, p)
	return out
}

// Layout is a full Part 3 design.
type Layout struct {
	F           float64
	Sense       int
	NG          int
	Theta0, DTh float64 // rad
	Shape       GroupShape
	Comb        Combine
	RLoad       float64
}

// Thetas are the group angles.
func (l Layout) Thetas() []float64 {
	t := make([]float64, l.NG)
	for k := range t {
		t[k] = l.Theta0 + float64(k)*l.DTh
	}
	return t
}

// Part3Volts is the amplitude of each drive source in Part 3 [V].
var Part3Volts = 3.4

// Groups builds the group coils of a layout for a net.
func (m *Model) Groups(n *Net, l Layout) ([]Group, [][]float64) {
	th := l.Thetas()
	gs := make([]Group, l.NG)
	lg := l.Shape.SelfL()
	w := 2 * math.Pi * l.F
	for k := range gs {
		gs[k] = Group{L: lg, R: l.Shape.Resistance(l.F), Ctune: 1 / (w * w * lg),
			MSec: n.SectionMutual(m, m.PassMutual(l.Shape, th[k]))}
	}
	mg := newMat(l.NG)
	for i := range mg {
		mg[i][i] = lg
		for j := i + 1; j < l.NG; j++ {
			v := GroupMutual(l.Shape, th[i], th[j])
			mg[i][j], mg[j][i] = v, v
		}
	}
	return gs, mg
}

// DriveSources returns the sin (A) and cos (B) sources of a driven pair.
func (n *Net) DriveSources(sense int, v float64) []Source {
	bn := n.B[0]
	if sense == Dipole {
		bn = n.B[n.NS]
	}
	return []Source{{Node: n.A[0], V: v, PhaseDeg: -90, R: SourceOhms}, {Node: bn, V: v, PhaseDeg: 0, R: SourceOhms}}
}

// Evaluate solves a layout on the driven net (sense must match the net).
func (m *Model) Evaluate(n *Net, l Layout) (*ACResult, error) {
	n.SetFreq(l.F)
	gs, mg := m.Groups(n, l)
	return n.SolveAC(l.F, n.DriveSources(l.Sense, Part3Volts), gs, mg, l.Comb, l.RLoad)
}

// EMFMap returns, for a reference coil at each angle, the amplitude [V]
// of its open-circuit EMF at each frequency, from the section currents of
// the unloaded driven pair (an open coil draws no current, so this is
// exact). Rows: frequencies; columns: angles.
func (m *Model) EMFMap(n *Net, s GroupShape, freqs, thetas []float64, sense int) ([][]float64, error) {
	pass := make([][]float64, len(thetas))
	parallelFor(len(thetas), func(k int) { pass[k] = m.PassMutual(s, thetas[k]) })
	out := make([][]float64, len(freqs))
	for i, f := range freqs {
		n.SetFreq(f)
		r, err := n.SolveAC(f, n.DriveSources(sense, Part3Volts), nil, nil, Separate, 1)
		if err != nil {
			return nil, err
		}
		ip := make([]complex128, 2*m.P.Sec) // current per pass, summed over the circuits of each bundle
		for g, I := range r.I {
			ip[(g/n.NS)*m.P.Sec+g%m.P.Sec] += I
		}
		out[i] = make([]float64, len(thetas))
		for k := range thetas {
			var e complex128
			for q, I := range ip {
				e += complex(0, 2*math.Pi*f*pass[k][q]) * I
			}
			out[i][k] = cmplx.Abs(e)
		}
	}
	return out, nil
}

// PowerMap returns, for one reference coil (tuned, loaded by R_L equal to
// its wire resistance) at each angle, the power into its load from the
// full AC solution with the coil in the network (it loads the drive), and
// the power the two sources give. Rows: frequencies; columns: angles.
func (m *Model) PowerMap(n *Net, s GroupShape, freqs, thetas []float64, sense int) (pl, ps [][]float64, err error) {
	pass := make([][]float64, len(thetas))
	parallelFor(len(thetas), func(k int) { pass[k] = m.PassMutual(s, thetas[k]) })
	pl = make([][]float64, len(freqs))
	ps = make([][]float64, len(freqs))
	lg := s.SelfL()
	for i, f := range freqs {
		pl[i] = make([]float64, len(thetas))
		ps[i] = make([]float64, len(thetas))
		nn := *n
		nn.R = append([]float64(nil), n.R...)
		nn.SetFreq(f)
		rg := s.Resistance(f)
		w := 2 * math.Pi * f
		for k := range thetas {
			g := []Group{{L: lg, R: rg, Ctune: 1 / (w * w * lg), MSec: nn.groupMSec(m, pass[k])}}
			r, e := nn.SolveAC(f, nn.DriveSources(sense, Part3Volts), g, [][]float64{{lg}}, Separate, rg)
			if e != nil {
				return nil, nil, e
			}
			pl[i][k] = r.PLoad
			ps[i][k] = r.PSrc[0] + r.PSrc[1]
		}
	}
	return
}
