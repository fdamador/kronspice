// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"sort"
	"testing"

	"github.com/fdamador/kronspice/em"
)

// TestCoilNetwork: both windings use every physical position once; the
// turn-by-turn inductance sum agrees with Nagaoka's current-sheet formula
// (filament vs sheet, within 1 %) and is the same for both windings, as
// the strands of a series bifilar coil carry the same current the same
// way (em.Coil.Inductance doc); the inter-section capacitance matrix is
// symmetric; the bifilar coil holds (N − 1) turn gaps between its strands.
func TestCoilNetwork(t *testing.T) {
	ord, bif := NewCoil(Solenoid(), false), NewCoil(Solenoid(), true)
	for _, c := range []*Coil{ord, bif, NewCoil(Cone(), true), NewCoil(Flat(), true)} {
		if c.Geom.N()%2 != 0 {
			t.Fatalf("%s: odd turn count", c.Name)
		}
		if l := c.Geom.WireLength(); math.Abs(l-ord.Geom.WireLength())/ord.Geom.WireLength() > 0.02 { // whole even turns: within 2 %
			t.Errorf("%s: wire %.2f m, solenoid %.2f m", c.Name, l, ord.Geom.WireLength())
		}
		p := append([]int(nil), c.Pos...)
		sort.Ints(p)
		for i := range p {
			if p[i] != i {
				t.Fatalf("%s: positions not a permutation", c.Name)
			}
		}
		for i := range c.CAdj {
			for j := range c.CAdj {
				if c.CAdj[i][j] != c.CAdj[j][i] {
					t.Fatalf("%s: CAdj not symmetric", c.Name)
				}
			}
		}
	}
	if l, n := ord.TotalL(), em.Nagaoka(Radius, 400*Pitch, 400); math.Abs(l-n)/n > 0.01 {
		t.Errorf("L %.4g vs Nagaoka %.4g", l, n)
	}
	if math.Abs(ord.TotalL()-bif.TotalL())/ord.TotalL() > 1e-12 {
		t.Errorf("bifilar L %.6g vs ordinary %.6g", bif.TotalL(), ord.TotalL())
	}
	sum := 0.0
	for i := range bif.CAdj {
		for j := i + 1; j < len(bif.CAdj); j++ {
			sum += bif.CAdj[i][j]
		}
	}
	if want := float64(bif.Geom.N()-1) * bif.Ctt; math.Abs(sum-want)/want > 1e-12 {
		t.Errorf("bifilar strand capacitance %.4g, want (N−1)·C_tt = %.4g", sum, want)
	}
}

// TestFlatSpiralInductance: the flat spiral's turn-by-turn inductance
// against Wheeler's flat-spiral formula L[µH] = a²N²/(8a + 11c), a the mean
// radius and c the radial width in inches (H. A. Wheeler, Proc. IRE 16,
// 1398, 1928; about 5 % for c > 0.2a), held to 10 %.
func TestFlatSpiralInductance(t *testing.T) {
	c := NewCoil(Flat(), false)
	g := c.Geom
	const inch = 0.0254
	a := (g.R[0] + g.R[g.N()-1]) / 2 / inch
	w := (g.R[g.N()-1] - g.R[0]) / inch
	n := float64(g.N())
	want := a * a * n * n / (8*a + 11*w) * 1e-6
	if l := c.TotalL(); math.Abs(l-want)/want > 0.10 {
		t.Errorf("flat spiral L %.4g H, Wheeler %.4g H", l, want)
	}
	t.Logf("flat spiral %.3f mH, Wheeler %.3f mH", c.TotalL()*1e3, want*1e3)
}

// TestPulseBuildUp: on the ordinary coil, pulses at f1 build more than
// three times higher than pulses at 1.5·f1 within the same time, and the
// ledger closes: energy in = heat + stored to 1 % of the input.
func TestPulseBuildUp(t *testing.T) {
	c := NewCoil(Solenoid(), false)
	f1, _, _, err := c.Resonance(1e4, 1e8)
	if err != nil {
		t.Fatal(err)
	}
	on, err := c.Pulse(f1, f1, 0.1/f1, 60, false)
	if err != nil {
		t.Fatal(err)
	}
	off, err := c.Pulse(f1, 1.5*f1, 0.1/f1, 60, false)
	if err != nil {
		t.Fatal(err)
	}
	if on.PeakTop < 3*off.PeakTop {
		t.Errorf("in step %.1f V, off step %.1f V", on.PeakTop, off.PeakTop)
	}
	for _, r := range []PulseRun{on, off} {
		if math.Abs(r.EIn-r.EDiss-r.EStore) > 0.01*r.EIn {
			t.Errorf("rate %.3g: in %.4g, heat %.4g, stored %.4g", r.Rate, r.EIn, r.EDiss, r.EStore)
		}
	}
}

// TestPickup: the tuning capacitor resonates the pickup at f1
// (1/(2π√(LC)) = f1); with the pickup tuned and loaded, the ledger closes
// to 1 % and the load never gets more than the power that reaches the
// primary past the source's own resistance.
func TestPickup(t *testing.T) {
	c := NewCoil(Solenoid(), false)
	f1, _, _, err := c.Resonance(1e4, 1e8)
	if err != nil {
		t.Fatal(err)
	}
	if f := 1 / (2 * math.Pi * math.Sqrt(PickupL()*TuningC(f1))); math.Abs(f-f1)/f1 > 1e-12 {
		t.Errorf("tuned pickup resonates at %.6g Hz, f1 %.6g Hz", f, f1)
	}
	pk := Pickup{Centre: 0.5 * c.Geom.Z[c.Geom.N()-1], Tuned: true, CFactor: 1, RLoad: 0.2}
	f, _, err := c.PickupAC(f1, pk)
	if err != nil {
		t.Fatal(err)
	}
	r, err := c.PulsePickup(f1, f, 60, pk)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.EIn-r.EDiss-r.EStore) > 0.01*r.EIn {
		t.Errorf("ledger: in %.4g, heat %.4g, stored %.4g", r.EIn, r.EDiss, r.EStore)
	}
	if r.PLoad <= 0 || r.PLoad >= r.PIn-r.PSource {
		t.Errorf("load %.4g W, past the source resistor %.4g W", r.PLoad, r.PIn-r.PSource)
	}
}
