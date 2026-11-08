// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
)

// Geometry of the Rodin coil and of the comparison windings.
//
// The former is the coil of Ramírez Galindo et al., Rev. Mex. Fís. 70,
// 030901 (2024): outer diameter 9.95 cm, aperture 2.9 cm, 12 poles 30°
// apart, each pass advancing 150° (5 poles) through the hole, the sequence
// 1→6→11→4→9→2→7→12→5→10→3→8→1 (a {12/5} star). One circuit of 12 passes
// goes 5 times round the axis and 12 times through the hole: the (5, 12)
// torus knot class. The photon electron of examples/photon_electron (H23)
// is the (2, 1) class.
//
// Wire model ("taut wire"): seen from above, a pass is the straight chord
// from pole P_k to pole P_{k+5}; its first half lies on the top surface of
// the torus, its second half on the bottom, and it changes side through
// the hole at the chord's midpoint (where the chord is nearest the axis).
// Inside the aperture radius the wire crosses at mid-height (z = 0); when a
// layer's chord never reaches the aperture, the wire drops vertically down
// the inner surface at the midpoint.
//
// Bundles (declared): the paper gives no turn count. A winding of Nc
// circuits lays the same 12-pass path Nc times, so its wires form a round
// bundle along one path (radius BundleRadius(Nc), fill 0.7); winding A's
// bundle lies on the former (layer 0), the counter-wound B's bundle on top
// of it (layer 1, one bundle pitch out). Nc is fitted to the paper's
// centre field (Table II) by least squares: Circuits = 12.
const (
	OuterRadius = 0.0995 / 2 // m, Rp of the paper (layer 0)
	HoleRadius  = 0.029 / 2  // m, aperture radius of the paper (layer 0)
	Poles       = 12
	Skip        = 5      // poles advanced per pass (150°)
	ChordPts    = 40     // points per pass
	WireRadius  = 0.3e-3 // m: 0.6 mm copper (declared; the paper gives no gauge)
	Fill        = 0.7    // copper fill of a bundle (declared)
	Circuits    = 12     // circuits per winding (fitted to Table II, FitCircuits)
)

// BundleRadius is the radius of a round bundle of n wires.
func BundleRadius(n int) float64 { return WireRadius * math.Sqrt(float64(n)/Fill) }

// Pitch is the distance between the centres of the A and B bundles of n
// circuits each (two bundle radii plus 0.1 mm).
func Pitch(n int) float64 { return 2*BundleRadius(n) + 0.1e-3 }

// TorusR0, TorusA are the former's major and minor radii.
var (
	TorusR0 = (OuterRadius + HoleRadius) / 2
	TorusA  = (OuterRadius - HoleRadius) / 2
)

// P3 is a point.
type P3 = [3]float64

// rodinCircuit returns one closed circuit (12 passes) of the Rodin winding
// on layer `layer`; mirror reflects it (y → −y), which gives the
// counter-wound partner (P_k → P_{k−5}); rot turns it about the axis.
func rodinCircuit(layer int, mirror bool, rot, pitch float64) []P3 {
	a := TorusA + (float64(layer)+0.5)*pitch
	rp := TorusR0 + a
	rin := TorusR0 - a
	var pts []P3
	k := 0
	for pass := 0; pass < Poles; pass++ {
		phi0 := 2 * math.Pi * float64(k) / Poles
		phi1 := 2 * math.Pi * float64(k+Skip) / Poles
		x0, y0 := rp*math.Cos(phi0), rp*math.Sin(phi0)
		x1, y1 := rp*math.Cos(phi1), rp*math.Sin(phi1)
		for i := 0; i < ChordPts; i++ {
			s := float64(i) / ChordPts
			x, y := x0+(x1-x0)*s, y0+(y1-y0)*s
			rho := math.Hypot(x, y)
			h := a*a - (rho-TorusR0)*(rho-TorusR0)
			z := 0.0
			if h > 0 && rho > rin {
				z = math.Sqrt(h)
			}
			if s >= 0.5 {
				z = -z
			}
			if i == ChordPts/2 && z != 0 {
				// vertical drop down the inner surface at the midpoint
				pts = append(pts, P3{x, y, -z})
			}
			pts = append(pts, P3{x, y, z})
		}
		k += Skip
	}
	return orient(pts, mirror, rot)
}

// torusKnot returns the (p, q) torus knot on layer `layer`: p times round
// the axis, q times through the hole, n points.
func torusKnot(p, q, layer, n int, mirror bool, rot, pitch float64) []P3 {
	a := TorusA + (float64(layer)+0.5)*pitch
	pts := make([]P3, n)
	for i := range pts {
		u := 2 * math.Pi * float64(i) / float64(n)
		ph, th := float64(p)*u, float64(q)*u
		r := TorusR0 + a*math.Cos(th)
		pts[i] = P3{r * math.Cos(ph), r * math.Sin(ph), a * math.Sin(th)}
	}
	return orient(pts, mirror, rot)
}

func orient(pts []P3, mirror bool, rot float64) []P3 {
	c, s := math.Cos(rot), math.Sin(rot)
	out := make([]P3, len(pts))
	for i, p := range pts {
		if mirror {
			p[1] = -p[1]
		}
		out[i] = P3{c*p[0] - s*p[1], s*p[0] + c*p[1], p[2]}
	}
	return out
}

// Pair is the counter-wound pair: closed paths of A (layer 0) and of its
// mirror B (layer 1), each laid Nc times, and the sections of one circuit
// (Sec per circuit, 12 = one per pass for the Rodin path).
type Pair struct {
	Kind   Kind
	Nc     int
	A, B   []P3 // closed polylines (last point = first)
	Sec    int
	Radius float64 // bundle radius used for the self terms
	Pitch  float64 // A to B bundle centres
}

// Kind of winding.
type Kind int

const (
	Rodin  Kind = iota // (5, 12) taut-wire Rodin path
	Knot21             // (2, 1) torus knot
)

func (k Kind) String() string {
	if k == Rodin {
		return "Rodin (5, 12)"
	}
	return "(2, 1) torus knot"
}

// NewPair builds the pair of nc circuits per winding, one section per
// pass (12 per circuit).
func NewPair(kind Kind, nc int) *Pair { return NewPairSec(kind, nc, Poles) }

// NewPairSec is NewPair with sec sections per circuit.
func NewPairSec(kind Kind, nc, sec int) *Pair {
	p := &Pair{Kind: kind, Nc: nc, Sec: sec, Radius: BundleRadius(nc), Pitch: Pitch(nc)}
	switch kind {
	case Rodin:
		p.A, p.B = rodinCircuit(0, false, 0, p.Pitch), rodinCircuit(1, true, 0, p.Pitch)
	case Knot21:
		p.A, p.B = torusKnot(2, 1, 0, Knot21Pts, false, 0, p.Pitch), torusKnot(2, 1, 1, Knot21Pts, true, 0, p.Pitch)
	}
	p.A = append(p.A, p.A[0])
	p.B = append(p.B, p.B[0])
	return p
}

// PathLength is the length of one circuit of A [m].
func (p *Pair) PathLength() float64 { return polyLen(p.A) }

// WireLength is the wire of one winding [m].
func (p *Pair) WireLength() float64 { return float64(p.Nc) * polyLen(p.A) }

func polyLen(pl []P3) float64 {
	s := 0.0
	for i := 0; i+1 < len(pl); i++ {
		s += dist(pl[i], pl[i+1])
	}
	return s
}

// Knot21Pts is the number of points of one (2, 1) circuit.
const Knot21Pts = 480

// OuterShell is the minor radius of the outside of the B bundle.
func (p *Pair) OuterShell() float64 { return TorusA + 1.5*p.Pitch + p.Radius }

func dist(a, b P3) float64 {
	return math.Sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]))
}
