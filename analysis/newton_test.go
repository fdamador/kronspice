// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"
	"sync"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/netlist"
)

func mustDeck(t testing.TB, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// A stack of diodes driven hard from a stiff source, plus a zener in
// breakdown: plain Newton needs limiting, and the fallbacks have work to do.
const hardDiodes = `hard
V1 in 0 100
R1 in a 10
D1 a b dd
D2 b c dd
D3 c 0 dd
R2 a z 100
DZ 0 z dz
.model dd D(IS=1e-15 N=1 RS=0.1)
.model dz D(IS=1e-14 BV=5.6 IBV=5m RS=1)
.op
`

func TestDiodeOPLedger(t *testing.T) {
	c := mustDeck(t, hardDiodes)
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-12 {
		t.Fatalf("ledger residual %g", rr)
	}
	d1, _ := r.Element("D1")
	if d1.Current < 1 || d1.ExteriorPower[exchange.Thermal] != -d1.PowerAbsorbed {
		t.Fatalf("D1 = %+v", d1)
	}
	if vz, _ := r.Voltage("z"); vz > 5.6 || vz < 5 {
		t.Fatalf("zener clamps at %g V", vz)
	}
}

func TestNoConvergenceIsTyped(t *testing.T) {
	c := mustDeck(t, hardDiodes)
	_, err := OP(c, Options{ITL1: 2, NoGminStepping: true, NoSourceStepping: true, NoPseudoTransient: true})
	if !errors.Is(err, ErrNoConvergence) {
		t.Fatalf("err = %v, want ErrNoConvergence", err)
	}
}

// Gmin stepping and source stepping, run on their own, reach the same
// operating point as plain Newton.
func TestSteppingAgrees(t *testing.T) {
	c := mustDeck(t, hardDiodes)
	ref, err := New(c)
	if err != nil {
		t.Fatal(err)
	}
	o := Options{}
	o.defaults()
	ref.lc.Mode = 0
	if err := ref.dcSolve(&o); err != nil {
		t.Fatal(err)
	}
	for name, run := range map[string]func(*Sim) error{
		"gmin":   func(s *Sim) error { return s.gminStepping(&o) },
		"source": func(s *Sim) error { return s.sourceStepping(&o) },
		"ptran":  func(s *Sim) error { return s.pseudoTransient(&o) },
	} {
		s, _ := New(c)
		if err := run(s); err != nil {
			t.Fatalf("%s stepping: %v", name, err)
		}
		for i := range s.x {
			if d := math.Abs(s.x[i] - ref.x[i]); d > 1e-3*math.Abs(ref.x[i])+1e-5 {
				t.Errorf("%s stepping: x[%d] = %g, Newton %g", name, i, s.x[i], ref.x[i])
			}
		}
	}
}

// Small-signal AC of a forward-biased diode with RS and junction charge:
// H = Zd / (R + Zd), Zd = RS + 1/(gd + jωC), with gd and C from the model
// equations at the solved junction voltage.
func TestDiodeAC(t *testing.T) {
	c := mustDeck(t, `ac
V1 in 0 DC 5 AC 1
R1 in a 1k
D1 a 0 dm
.model dm D(IS=1e-14 RS=20 CJO=4p VJ=0.75 M=0.4 TT=10n)
.ac dec 5 1k 1g
`)
	// AC linearizes at the operating point polished to round-off
	// (Sim.polish), so the expected response uses the exact junction
	// voltage, from a solve with tolerances far below the defaults.
	op, err := OP(c, Options{RelTol: 1e-10, VnTol: 1e-13, AbsTol: 1e-18})
	if err != nil {
		t.Fatal(err)
	}
	d1, _ := op.Element("D1")
	vd := d1.Voltage - 20*d1.Current
	vt := 1.380649e-23 * 300.15 / 1.602176634e-19
	gd := 1e-14*math.Exp(vd/vt)/vt + 1e-12
	fcv := 0.5 * 0.75
	f2 := math.Pow(0.5, 1.4)
	cap := 10e-9*(gd-1e-12) + 4e-12/f2*(1-0.5*1.4+0.4*vd/0.75)
	if vd < fcv {
		t.Fatalf("expected forward bias above FC·VJ, vd = %g", vd)
	}
	r, err := AC(c, ACSpec{Variation: "dec", Points: 5, FStart: 1e3, FStop: 1e9}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	va, _ := r.NodeVoltage("a")
	for k, f := range r.Freq {
		zd := complex(20, 0) + 1/complex(gd, 2*math.Pi*f*cap)
		want := zd / (1000 + zd)
		if cmplx.Abs(va[k]-want) > 1e-6*cmplx.Abs(want) {
			t.Errorf("f = %g: H = %v, want %v", f, va[k], want)
		}
	}
}

// Diode transients in parallel goroutines give bit-identical results
// (CLAUDE.md §1 goals 3 and 5). Run with -race.
func TestParallelDiodeTransients(t *testing.T) {
	src := "r\nV1 in 0 SIN(0 10 1k)\nD1 in out dr\nR1 out 0 1k\nC1 out 0 10u\n.model dr D(IS=1e-14 CJO=2p TT=5n)\n"
	c := mustDeck(t, src)
	run := func() []float64 {
		wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 2e-3}, TranOptions{})
		if err != nil {
			t.Error(err)
			return nil
		}
		v, _ := wf.NodeVoltage("out")
		return v
	}
	want := run()
	var wg sync.WaitGroup
	for g := 0; g < 8; g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			got := run()
			if len(got) != len(want) {
				t.Errorf("length %d, want %d", len(got), len(want))
				return
			}
			for i := range got {
				if got[i] != want[i] {
					t.Errorf("point %d differs: %v vs %v", i, got[i], want[i])
					return
				}
			}
		}()
	}
	wg.Wait()
}

func BenchmarkDiodeOP(b *testing.B) {
	c := mustDeck(b, hardDiodes)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := OP(c, Options{}); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkRectifierTran(b *testing.B) {
	c := mustDeck(b, "r\nV1 in 0 SIN(0 10 1k)\nD1 in out dr\nR1 out 0 1k\nC1 out 0 10u\n.model dr D(IS=1e-14)\n")
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := Tran(c, TranSpec{Step: 10e-6, Stop: 3e-3}, TranOptions{}); err != nil {
			b.Fatal(err)
		}
	}
}
