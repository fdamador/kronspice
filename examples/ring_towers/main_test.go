// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Traveling: flat envelope and power circulating east; standing: nodes and
// no net flow; cancel: almost nothing. In each state the towers supply
// exactly the ring loss, and the ledger closes.
func TestStates(t *testing.T) {
	for _, s := range States() {
		r, err := Run(s)
		if err != nil {
			t.Fatal(err)
		}
		var flow, peak float64
		for k, e := range r.Envelope {
			flow += r.Flow[k] / N
			peak = math.Max(peak, e)
		}
		pt := r.PTowers[0] + r.PTowers[1] + r.PTowers[2] + r.PTowers[3]
		switch s.Name {
		case "traveling":
			if r.Uniformity() < 0.85 || flow < 0.1 {
				t.Errorf("traveling: envelope min/max %.3f, flow %.3g W", r.Uniformity(), flow)
			}
		case "standing":
			if r.Uniformity() > 0.05 || math.Abs(flow) > 1e-3 {
				t.Errorf("standing: envelope min/max %.3f, flow %.3g W", r.Uniformity(), flow)
			}
		case "cancel":
			if peak > 0.2 {
				t.Errorf("cancel: peak %.3f V", peak)
			}
			continue
		}
		if math.Abs(pt-r.PLoss) > 0.02*r.PLoss {
			t.Errorf("%s: towers %.5g W, ring loss %.5g W", s.Name, pt, r.PLoss)
		}
		if math.Abs(r.Mismatch) > 1e-3*r.EnergyIn {
			t.Errorf("%s: ledger mismatch %g of %g J", s.Name, r.Mismatch, r.EnergyIn)
		}
	}
}

// Four towers: the C–D standing wave is A–B's pattern rotated by 90°; with
// both pairs on the same reference waves the two patterns coincide and
// add (twice the amplitude, four times the stored energy and the tower
// power), all four quarter points swinging in step; with C–D in antiphase
// they cancel. The towers always supply exactly the ring loss.
func TestFourTowers(t *testing.T) {
	res := map[string]*Result{}
	for _, s := range FourStates() {
		r, err := Run(s)
		if err != nil {
			t.Fatal(err)
		}
		res[s.Name] = r
		var pt float64
		for _, p := range r.PTowers {
			pt += p
		}
		if math.Abs(pt-r.PLoss) > 0.02*r.PLoss+1e-6 {
			t.Errorf("%s: towers %.5g W, ring loss %.5g W", s.Name, pt, r.PLoss)
		}
		if math.Abs(r.Mismatch) > 1e-3*r.EnergyIn {
			t.Errorf("%s: ledger mismatch %g of %g J", s.Name, r.Mismatch, r.EnergyIn)
		}
	}
	one, both, anti := res["A–B standing"], res["both, synchronized"], res["both, C–D antiphase"]
	for _, at := range []int{nodeA, nodeC, nodeB, nodeD} {
		if ratio := both.Envelope[at] / one.Envelope[at]; ratio < 1.9 {
			t.Errorf("node %d: synchronized/single %.3f, want ≈ 2", at, ratio)
		}
		if anti.Envelope[at] > 0.15*one.Envelope[at] {
			t.Errorf("node %d: antiphase leaves %.3f V", at, anti.Envelope[at])
		}
	}
	if q := both.Stored / one.Stored; math.Abs(q-4) > 0.05 {
		t.Errorf("stored energy ratio %.3f, want 4", q)
	}
	// Same sign at all four quarter points when V(0°) peaks.
	best, bf := 0.0, 0
	for i, fr := range both.Frames {
		if math.Abs(fr[nodeA]) > best {
			best, bf = math.Abs(fr[nodeA]), i
		}
	}
	for _, at := range []int{nodeC, nodeB, nodeD} {
		if both.Frames[bf][at]*both.Frames[bf][nodeA] <= 0 {
			t.Errorf("node %d not in step with 0°", at)
		}
	}
}

// DC bias per pair: on a continuous ring both pairs' supplies see one
// conductor, so the ring sits at the average and a steady current
// (VAB − VCD)/RFeed runs between the supplies; on a sectored ring
// (blocking capacitors at 45°, 135°, 225°, 315°) each sector holds its own
// pair's bias and no DC flows. In both, the RF standing wave rides
// unchanged on the DC level and the towers still supply exactly the ring
// loss.
func TestBias(t *testing.T) {
	unb, err := Run(FourStates()[2])
	if err != nil {
		t.Fatal(err)
	}
	for _, sec := range []bool{false, true} {
		b := DefaultBias(100, 50, sec)
		lv, sup, err := BiasOP(b)
		if err != nil {
			t.Fatal(err)
		}
		want := DCLevels(b)
		for k := range lv {
			if math.Abs(lv[k]-want[k]) > 1e-3 {
				t.Fatalf("sectored=%v node %d: OP %g V, expected %g V", sec, k, lv[k], want[k])
			}
		}
		iab := sup["A"] + sup["B"]
		if sec {
			if math.Abs(iab) > 1e-9 {
				t.Errorf("sectored: DC current %g A, want 0", iab)
			}
		} else if math.Abs(iab-(b.VAB-b.VCD)/b.RFeed) > 1e-9 {
			t.Errorf("continuous: A–B supplies deliver %g A, want %g", iab, (b.VAB-b.VCD)/b.RFeed)
		}
		r, err := RunBiased(FourStates()[2], b)
		if err != nil {
			t.Fatal(err)
		}
		var pt float64
		for _, p := range r.PTowers {
			pt += p
		}
		if math.Abs(pt-r.PLoss) > 0.02*r.PLoss {
			t.Errorf("sectored=%v: towers %g W, loss %g W", sec, pt, r.PLoss)
		}
		for _, at := range []int{nodeA, nodeC, nodeB, nodeD} {
			a := math.Hypot(r.Phasor[at][0], r.Phasor[at][1])
			a0 := math.Hypot(unb.Phasor[at][0], unb.Phasor[at][1])
			if math.Abs(a/a0-1) > 0.03 {
				t.Errorf("sectored=%v node %d: RF %g V vs unbiased %g V", sec, at, a, a0)
			}
			if math.Abs(r.MeanV[at]-want[at]) > 0.2 {
				t.Errorf("sectored=%v node %d: mean %g V, bias %g V", sec, at, r.MeanV[at], want[at])
			}
		}
	}
}

// Eight towers: the 45°-rotated pairs make the same standing pattern, so
// all four pairs add — four times one pair's amplitude at 0°, 45°, 90°, …
// (the 45° points swinging opposite to the quarter points), sixteen times
// the stored energy, towers' input equal to the ring loss. With the 45°
// set in antiphase the two sets cancel.
func TestEightTowers(t *testing.T) {
	one, err := Run(FourStates()[0])
	if err != nil {
		t.Fatal(err)
	}
	res := map[string]*Result{}
	for _, s := range EightStates() {
		r, err := Run(s)
		if err != nil {
			t.Fatal(err)
		}
		res[s.Name] = r
		var pt float64
		for _, p := range r.PTowers {
			pt += p
		}
		if math.Abs(pt-r.PLoss) > 0.02*r.PLoss+1e-6 {
			t.Errorf("%s: towers %g W, loss %g W", s.Name, pt, r.PLoss)
		}
	}
	all, anti := res["eight towers, synchronized"], res["eight towers, 45° set in antiphase"]
	amp := func(r *Result, k int) float64 { return math.Hypot(r.Phasor[k][0], r.Phasor[k][1]) }
	for k := 0; k < N; k += N / 8 {
		// ±10 %: one pair's amplitude varies around the ring with the loss
		// per lap (0.2 Np), the four-pair sum much less.
		if q := amp(all, k) / amp(one, k); q < 3.5 || q > 4.5 {
			t.Errorf("node %d: eight/two = %.3f, want ≈ 4", k, q)
		}
		if amp(anti, k) > 0.15*amp(one, k) {
			t.Errorf("node %d: antiphase leaves %.3f V", k, amp(anti, k))
		}
	}
	if q := all.Stored / one.Stored; math.Abs(q-16) > 0.5 {
		t.Errorf("stored energy ratio %.2f, want 16", q)
	}
	// Opposite sign at 45° to 0° at the instant 0° peaks.
	v0 := all.Phasor[nodeA][0]
	if v0*all.Phasor[nodeE][0] >= 0 {
		t.Errorf("45° point not in antiphase with 0°: %g, %g", v0, all.Phasor[nodeE][0])
	}
}
