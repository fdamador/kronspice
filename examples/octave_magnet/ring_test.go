// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The surface-charge field equals the point-dipole twin outside the ring.
func TestRingTwin(t *testing.T) {
	for _, r := range []Ring{RadialTwelve(), TangentialTwelve()} {
		if e := TwinError(r); e > 2e-4 {
			t.Errorf("%s: charge vs dipole twin differ by %.2g of |B|", r.Short, e)
		}
	}
}

// Every section's surface charge sums to zero (∇·M integrates to zero over a
// closed body), and B outside is divergence- and curl-free.
func TestRingChargeDivCurl(t *testing.T) {
	for _, r := range []Ring{RadialTwelve(), TangentialTwelve()} {
		cs := r.Charges(0.25e-3, 2)
		var net, abs float64
		for _, c := range cs {
			net, abs = net+c.Q, abs+math.Abs(c.Q)
		}
		if math.Abs(net)/abs > 1e-12 {
			t.Errorf("%s: net charge %.2g of Σ|Q|", r.Short, net/abs)
		}
		for _, p := range ringTestPoints(r)[:3] {
			if d, c := DivCurl(cs, p); d > 1e-5 || c > 1e-5 {
				t.Errorf("%s at %v: ∇·B %.2g, ∇×B %.2g of |∇B|", r.Short, p, d, c)
			}
		}
	}
}

// The bench scan: 12 sign changes round both alternating rings, none round
// the all-outward control; N poles at the section centres (radial) or at
// the joints (along the ring).
func TestRingScan(t *testing.T) {
	for _, c := range []struct {
		r       Ring
		changes int
		first   float64 // first N pole [deg]
	}{{RadialTwelve(), 12, 15}, {TangentialTwelve(), 12, 30}, {RadialUniform(), 0, -1}} {
		phi, br, _, _ := RingScan(c.r.Charges(0.5e-3, 1), c.r.Ro+3e-3, 0, 720)
		if n := CyclicSignChanges(br); n != c.changes {
			t.Errorf("%s: %d sign changes, want %d", c.r.Short, n, c.changes)
		}
		if c.first < 0 {
			continue
		}
		pk := PeakAngles(phi, br)
		if len(pk) != 6 {
			t.Fatalf("%s: %d N poles", c.r.Short, len(pk))
		}
		for i, p := range pk {
			if want := c.first + 60*float64(i); math.Abs(p*180/math.Pi-want) > 0.6 {
				t.Errorf("%s: N pole %d at %.2f°, want %.0f°", c.r.Short, i, p*180/math.Pi, want)
			}
		}
	}
}

// A field line from an outer N face of the alternating radial ring returns
// to the magnet (on an S face) without leaving the box.
func TestRingLinesClose(t *testing.T) {
	r := RadialTwelve()
	cs := r.Charges(0.5e-3, 1)
	seeds := Seeds(cs, r.Ro+0.4e-3, 0, 1, 4)
	if len(seeds) != 24 {
		t.Fatalf("%d seeds, want 24", len(seeds))
	}
	for i, l := range r.Lines(cs, seeds, 0.3e-3) {
		end := l[len(l)-1]
		if !r.Near(end, Skin) {
			t.Errorf("line %d ends outside the magnet at %v", i, end)
			continue
		}
		phi := math.Atan2(end[1], end[0])
		k := int(math.Floor(math.Mod(phi+2*math.Pi, 2*math.Pi) / (math.Pi / 6)))
		if r.Sign[k] > 0 {
			t.Errorf("line %d ends in section %d, an N section", i, k)
		}
	}
}

// The kronspice pickup coil: 6 electrical cycles per revolution (6 pole
// pairs) and the ledger closes.
func TestRingPickup(t *testing.T) {
	r := RadialTwelve()
	p := DefaultPickup(r)
	res, err := p.Run(r.Charges(0.5e-3, 1), 360)
	if err != nil {
		t.Fatal(err)
	}
	if res.Cycles != 6 {
		t.Errorf("%d cycles per revolution, want 6", res.Cycles)
	}
	if e := math.Abs(res.Ein-res.Ediss-res.dE) / res.Ein; e > 1e-4 {
		t.Errorf("ledger: in %.6g, heat %.6g, Δstored %.3g (residual %.2g)", res.Ein, res.Ediss, res.dE, e)
	}
}
