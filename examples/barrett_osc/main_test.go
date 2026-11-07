// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
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

// Fig. 2G bench: air-core pancakes (ordinary or bifilar) make no mixing
// product; an unbiased Metglas core makes odd orders, a biased one the
// three-wave idler too; Manley–Rowe holds for the core's flux curve.
// (A coarser step than the printed run: these are yes/no checks.)
func TestBench(t *testing.T) {
	b := DefaultBench()
	b.Points = 20000
	cs, err := benchCases(b, []int{1, 2, 3}) // bifilar air; Metglas; Metglas + bias
	if err != nil {
		t.Fatal(err)
	}
	if c := cs[0]; c.Idler > 1e-12 || c.IM3 > 1e-12 {
		t.Errorf("%s: idler %.2g, 2fp−fs %.2g", c.Name, c.Idler, c.IM3)
	}
	if c := cs[1]; c.IM3 < 1e-4 || c.Idler > 1e-2*c.IM3 {
		t.Errorf("unbiased Metglas: 2fp−fs %.2g, fp−fs %.2g", c.IM3, c.Idler)
	}
	if c := cs[2]; c.Idler < 1e-4 {
		t.Errorf("biased Metglas: idler %.2g", c.Idler)
	}
	for _, c := range cs[1:] {
		if math.Abs(c.MR1) > 2e-2 || math.Abs(c.MR2) > 2e-2 || c.Run.Ledger > 1e-2 {
			t.Errorf("%s: Manley–Rowe %.2g %.2g, ledger %.2g", c.Name, c.MR1, c.MR2, c.Run.Ledger)
		}
	}
}

// The Tesla-coil solver against kronspice: the coil with the gap held
// closed (a 0.2 Ω resistor) ringing from a charged Cp.
func TestTeslaTwin(t *testing.T) {
	tc := NewTeslaCoil(false)
	lp := tc.Lp()
	const dt, steps = 2e-9, 15000
	n := newTC()
	n.C("CP", "hv", "0", tc.Cp, tc.Vdc)
	n.R("RG", "hv", "pa", 0.2)
	l1 := n.L("LP", "pa", "pb", lp)
	n.R("RP", "pb", "0", 0.05)
	l2 := n.L("LS", "0", "sm", tc.Ls)
	n.R("RSEC", "sm", "top", tc.Rsec)
	n.C("CSELF", "0", "top", tc.Cself+tc.Ctop, 0)
	n.Couple(l1, l2, tc.K)
	if err := n.Prepare(); err != nil {
		t.Fatal(err)
	}
	var mine []float64
	for k := 0; k < steps; k++ {
		if err := n.Step(dt); err != nil {
			t.Fatal(err)
		}
		mine = append(mine, n.Volt("top"))
	}
	b := circuit.NewBuilder("tesla twin")
	b.CapacitorIC("CP", "hv", "0", tc.Cp, tc.Vdc)
	b.Resistor("RG", "hv", "pa", 0.2)
	b.Inductor("LP", "pa", "pb", lp)
	b.Resistor("RP", "pb", "0", 0.05)
	b.Inductor("LS", "0", "sm", tc.Ls)
	b.Resistor("RSEC", "sm", "top", tc.Rsec)
	b.Capacitor("CSELF", "0", "top", tc.Cself+tc.Ctop)
	b.Mutual("K1", "LP", "LS", tc.K)
	c, err := b.Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: dt * steps, MaxStep: dt, UIC: true}, analysis.TranOptions{MaxSteps: 1_000_000})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("top")
	ref := resample(wf.Time, v, dt, dt, steps)
	var peak, diff float64
	for j := range ref {
		peak = math.Max(peak, math.Abs(ref[j]))
		diff = math.Max(diff, math.Abs(ref[j]-mine[j]))
	}
	if diff > 1e-3*peak {
		t.Errorf("top-load voltage: max difference %.3g V of peak %.3g V", diff, peak)
	}
	if math.Abs(n.Residual()) > 1e-9*n.stored0 {
		t.Errorf("ledger residual %.2g J of %.3g J", n.Residual(), n.stored0)
	}
}

// Spark-gap Tesla coil: the timed gap puts sidebands on the signal at
// fs ± k·1 kHz, for both windings; the ledger closes.
func TestTeslaMixing(t *testing.T) {
	for _, bif := range []bool{false, true} {
		r, err := RunTesla(bif, true)
		if err != nil {
			t.Fatal(err)
		}
		for _, k := range []int{-1, 1} {
			if a := 2 * cmplx.Abs(r.MixLines[k]); a < 1e-3*r.SigLine {
				t.Errorf("bifilar=%v: sideband %+d at %.2g V of signal %.2g V", bif, k, a, r.SigLine)
			}
		}
		if math.Abs(r.Both.Residual) > 1e-4*r.Both.Supplied {
			t.Errorf("bifilar=%v: ledger %.2g of %.3g J", bif, r.Both.Residual, r.Both.Supplied)
		}
	}
}
