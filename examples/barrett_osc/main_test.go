// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"
)

// Fig. 2: coils and capacitors with floating terminals are a linear
// network: the transient load current holds the pump and signal lines AC
// analysis gives for each source alone, and no mixing product.
func TestFig2Linear(t *testing.T) {
	p := DefaultFig2
	for _, v := range []byte("ACEG") {
		s, err := p.Run(v, true, true, 3)
		if err != nil {
			t.Fatal(err)
		}
		if s.MixMax > 1e-10 {
			t.Errorf("Fig. 2%c: mixing product %.2g at %v", v, s.MixMax, s.MixAt)
		}
		for _, pump := range []bool{true, false} {
			a, err := p.ACLine(v, pump)
			if err != nil {
				t.Fatal(err)
			}
			got := s.Line[[2]int{1, 0}]
			if !pump {
				got = s.Line[[2]int{0, 1}]
			}
			if e := cmplx.Abs(got-a) / cmplx.Abs(a); e > 1e-3 {
				t.Errorf("Fig. 2%c pump=%v: transient vs AC %.2g", v, pump, e)
			}
		}
		if s.Ledger > 1e-6 {
			t.Errorf("Fig. 2%c: ledger mismatch %.2g", v, s.Ledger)
		}
	}
}

// Fig. 2A with a varactor: an idler appears, scales as Es·Ep (not |Es|²),
// is phase-conjugate to the signal, and the varactor obeys Manley–Rowe.
func TestFig2Varactor(t *testing.T) {
	q := VaractorFig2()
	s, err := q.Run('A', true, true, 4)
	if err != nil {
		t.Fatal(err)
	}
	if r := cmplx.Abs(dftReal(s.TerminalV, q.MP-q.MS)) / cmplx.Abs(dftReal(s.TerminalV, q.MP)); r < 1e-5 {
		t.Errorf("idler %.2g of the pump", r)
	}
	if a, b := q.ManleyRowe(s, 4); math.Abs(a) > 1e-4 || math.Abs(b) > 1e-4 {
		t.Errorf("Manley–Rowe sums %.2g, %.2g", a, b)
	}
	ex, err := q.idlerExponents()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(ex[0]-1) > 0.05 || math.Abs(ex[1]-1) > 0.08 {
		t.Errorf("idler exponents Vs %.3f Vp %.3f, want 1, 1", ex[0], ex[1])
	}
	d, err := q.idlerPhaseShift(60)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(d+60) > 5 {
		t.Errorf("idler phase moved %.2f° for +60° of signal, want −60°", d)
	}
}

// Fig. 3: the spark-gap chain gives a − + − pulse and the ledger closes.
func TestFig3FrozenWave(t *testing.T) {
	r, err := DefaultFrozenWave.Run()
	if err != nil {
		t.Fatal(err)
	}
	if len(r.Lobes) != 3 || r.Lobes[0] >= 0 || r.Lobes[1] <= 0 || r.Lobes[2] >= 0 {
		t.Errorf("lobes %v, want − + −", r.Lobes)
	}
	if math.Abs(r.Residual) > 1e-6*r.Stored0 {
		t.Errorf("ledger residual %.2g J of %.2g J", r.Residual, r.Stored0)
	}
	if r.Load >= r.Stored0 {
		t.Errorf("load %.3g J ≥ stored %.3g J", r.Load, r.Stored0)
	}
	for k, n := range r.Strikes {
		if n == 0 {
			t.Errorf("gap %d never struck", k+1)
		}
	}
}

// Fig. 4: one-terminal networks pass no coil current without stray
// capacitance; two-port cascades multiply as ABCD matrices (non-commuting,
// associative, det 1, unitary when lossless).
func TestFig4(t *testing.T) {
	q := DefaultFig4
	q.Cnode = 0
	for _, v := range []string{"A", "B(i)", "C(i)"} {
		y, cm, err := q.OneTerminal(v)
		if err != nil {
			t.Fatal(err)
		}
		if cm > 1e-12 || cmplx.Abs(y-1) > 1e-9 {
			t.Errorf("Fig. 4%s: coil current %.2g, Y/(jωC) %v", v, cm, y)
		}
	}
	r, err := TwoPorts()
	if err != nil {
		t.Fatal(err)
	}
	if r.CascadeVsProduct > 1e-9 || r.Associator > 1e-8 || r.DetErr > 1e-9 || r.UnitarityErr > 1e-4 {
		t.Errorf("two-port checks: %+v", r)
	}
	if r.Commutator < 0.05 {
		t.Errorf("XY vs YX only %.3g apart", r.Commutator)
	}
}

// Fig. 5: a sum of pump and signal has no line at fp − fs.
func TestFig5(t *testing.T) {
	for c, v := range FigFiveLines() {
		if v > 1e-12 {
			t.Errorf("component %d: %.2g at fp − fs", c, v)
		}
	}
}

func TestMainRuns(t *testing.T) {
	if testing.Short() {
		t.Skip("full run")
	}
	if err := run(); err != nil {
		t.Fatal(err)
	}
}
