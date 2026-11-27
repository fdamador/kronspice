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
	for _, c := range []*Coil{ord, bif, NewCoil(Cone(), true), NewCoil(Flat(), true), NewCounterCoil(Solenoid())} {
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

// TestCounterCoil: the opposite-current (folded) bifilar winding of the
// solenoid. Its strands' fields cancel: terminal inductance below 1 % of
// the same-direction (Tesla) winding's and its coupling to the primary
// below 2 %; every section lies on one strand (none spans the fold), so
// all K couplings are positive in Build's reference directions; it still
// holds N − 1 strand-to-strand gaps; pulses at its strongest peak close
// the ledger to 1 % of the input.
func TestCounterCoil(t *testing.T) {
	ctr, tes := NewCounterCoil(Solenoid()), NewCoil(Solenoid(), true)
	if r := ctr.TotalL() / tes.TotalL(); r > 0.01 || r <= 0 {
		t.Errorf("L counter / Tesla = %.4g", r)
	}
	sum := func(v []float64) (s float64) {
		for _, x := range v {
			s += x
		}
		return
	}
	if r := math.Abs(sum(ctr.Mp)) / sum(tes.Mp); r > 0.02 {
		t.Errorf("primary coupling counter / Tesla = %.4g", r)
	}
	n := ctr.Geom.N()
	for s := 0; s < NumSections(); s++ {
		for tt := sectionStart(s, n); tt < sectionStart(s+1, n); tt++ {
			if ctr.Sign[tt] != ctr.SectionSign(s) {
				t.Fatalf("section %d spans the fold", s)
			}
		}
	}
	for i := 0; i < NumSections(); i++ {
		for j := i + 1; j < NumSections(); j++ {
			if k := ctr.SectionSign(i) * ctr.SectionSign(j) * ctr.L[i][j]; k < 0 {
				t.Fatalf("K(%d,%d) negative in Build's frame", i, j)
			}
		}
	}
	cs := 0.0
	for i := range ctr.CAdj {
		for j := i + 1; j < len(ctr.CAdj); j++ {
			cs += ctr.CAdj[i][j]
		}
	}
	if want := float64(n-1) * ctr.Ctt; math.Abs(cs-want)/want > 1e-12 {
		t.Errorf("strand capacitance %.4g, want %.4g", cs, want)
	}
	f1, gain, _, weak, err := ctr.ResonanceWeak(1e3, 1e8)
	if err != nil {
		t.Fatal(err)
	}
	if !weak || gain >= 2 {
		t.Errorf("counter coil peak gain %.3g (weak %v): expected no resonance above 2", gain, weak)
	}
	r, err := ctr.Pulse(f1, f1, 0.1/f1, 60, false)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(r.EIn-r.EDiss-r.EStore) > 0.01*r.EIn {
		t.Errorf("ledger: in %.4g, heat %.4g, stored %.4g", r.EIn, r.EDiss, r.EStore)
	}
	t.Logf("counter: L %.4g µH, f1 %.4g Hz, gain %.3g, pulse peak %.3g V", ctr.TotalL()*1e6, f1, gain, r.PeakTop)
}

// TestParallelCoils: the bifilar solenoid with its strands in parallel.
// With currents the same way the terminal inductance is a quarter of the
// series winding's: for strands of equal self inductance L1 and mutual M,
// (L1² − M²)/(2(L1 − M)) = (L1 + M)/2 against 2(L1 + M) in series (the two
// interleaved strands differ only by half a pitch, so held to 1e-3); with
// currents opposing it is below 1 %. The parallel strands sit at nearly the same potential, so
// no gain-2 resonance is lost to strand capacitance: the same-direction
// winding resonates above the ordinary coil. Pulses at f1 close the ledger
// to 1 %.
func TestParallelCoils(t *testing.T) {
	ser := NewCoil(Solenoid(), true)
	same, opp := NewParallelCoil(Solenoid(), false), NewParallelCoil(Solenoid(), true)
	if r := same.TotalL() / ser.TotalL(); math.Abs(r-0.25) > 1e-3 {
		t.Errorf("parallel same-direction L / series = %.4f, want ≈ 0.25", r)
	}
	if r := opp.TotalL() / ser.TotalL(); r <= 0 || r > 0.01 {
		t.Errorf("parallel opposite L / series = %.4g", r)
	}
	ord := NewCoil(Solenoid(), false)
	fo, _, _, err := ord.Resonance(1e4, 1e8)
	if err != nil {
		t.Fatal(err)
	}
	for _, c := range []*Coil{same, opp} {
		f1, gain, _, weak, err := c.ResonanceWeak(1e3, 1e8)
		if err != nil {
			t.Fatal(err)
		}
		if c == same && (weak || f1 <= fo) {
			t.Errorf("parallel same-direction f1 %.4g Hz (weak %v), ordinary %.4g Hz", f1, weak, fo)
		}
		if c == opp && !weak {
			t.Errorf("parallel opposite: resonance with gain %.3g", gain)
		}
		r, err := c.Pulse(f1, f1, 0.1/f1, 60, false)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.EIn-r.EDiss-r.EStore) > 0.01*r.EIn {
			t.Errorf("%s ledger: in %.4g, heat %.4g, stored %.4g", c.Name, r.EIn, r.EDiss, r.EStore)
		}
		t.Logf("%s: L %.4g µH, f1 %.4g Hz, gain %.3g, pulse peak %.3g V, wire %.3g A", c.Name, c.TotalL()*1e6, f1, gain, r.PeakTop, r.PeakWire)
	}
}
