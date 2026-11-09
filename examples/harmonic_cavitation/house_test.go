// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// The heater design of part 4, rebuilt from its genes (run_output.txt).
func heaterForTest() Scored {
	return Scored{M: DimpleMachine(0.0843, 2, 0.006, 0.5, 3.33, 1.20, 0.8e-3, 1828, 1.18e5)}
}

// Part 11: the books close, the matched speed delivers the shaft power, the
// rows carry the load at the loop's best speed, and the condenser and
// emitter formulas hold.
func TestHouseSizing(t *testing.T) {
	p := DefaultParams()
	best := heaterForTest()
	for _, em := range emitters {
		load := 25 * houseFt2 * btuhToW
		s := SizeHouse(p, best, load, em)
		if math.Abs(s.Loop.PHome-load) > 1e-6*load {
			t.Errorf("%s: condenser %g W, load %g W", em.Name, s.Loop.PHome, load)
		}
		if math.Abs(s.Loop.Residual) > 1e-9*s.Loop.PElec {
			t.Errorf("%s: books residual %g W", em.Name, s.Loop.Residual)
		}
		if got := shaftAt(p, best.M, s.Rows, s.RPM, 3e5); math.Abs(got/s.Loop.PShaft-1) > 1e-6 {
			t.Errorf("%s: matched speed gives %g W, want %g", em.Name, got, s.Loop.PShaft)
		}
		if s.RPM > s.RPMFull*1.0001 {
			t.Errorf("%s: %d rows need %0.f rpm, above the loop's best %0.f", em.Name, s.Rows, s.RPM, s.RPMFull)
		}
		d1, d2 := em.Cond.T-em.TRet, em.Cond.T-em.TSup
		if math.Abs(s.CondArea*uCondenser*(d1-d2)/math.Log(d1/d2)-load) > 1e-6*load {
			t.Errorf("%s: condenser UA·LMTD ≠ load", em.Name)
		}
	}
	// EN 442: a radiator at mean water 60 °C, room 21.1 °C gives (38.9/50)^1.3 of its rating
	s := SizeHouse(p, best, 10000, emitters[1])
	if f := 10000 / s.RadRatedW; math.Abs(f-math.Pow((60-fToC(70))/50, 1.3)) > 1e-12 {
		t.Errorf("radiator factor %g", f)
	}
}

// The turbine's efficiency moves heat around the loop but cannot change the
// electricity a given home heat needs (all turbine work returns as heat).
func TestTurbineCannotLowerElectricity(t *testing.T) {
	p := DefaultParams()
	best := heaterForTest()
	em := emitters[1]
	load := 20000.0
	var first float64
	for i, et := range []float64{0.3, 0.7, 1.0} {
		u := LoopB(p, best.M, 3e5, em.Top, em.Flash, em.Cond, et, 1000)
		r := LoopB(p, best.M, 3e5, em.Top, em.Flash, em.Cond, et, 1000*load/u.PHome)
		if i == 0 {
			first = r.PElec
		} else if math.Abs(r.PElec/first-1) > 1e-3 {
			t.Errorf("turbine η %.1f: electricity %g W against %g W", et, r.PElec, first)
		}
		if r.PHome+r.PRoom > r.PElec*(1+1e-9) {
			t.Errorf("turbine η %.1f: heat %g W above electricity %g W", et, r.PHome+r.PRoom, r.PElec)
		}
	}
}
