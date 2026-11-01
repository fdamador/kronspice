// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/cmplx"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/netlist"
)

// N-type (tunnel-diode-like) curve, i = f(v), with a negative-slope region
// from 0.1 V to 0.3 V. With 0.4 V behind 50 Ω the load line crosses it three
// times: v = 0.4/6 (rising part), 0.26 (negative slope), 0.34 (valley).
const tunnel = `tunnel
V1 in 0 0.4
R1 in a 50
B1 a 0 I=pwl(V(a), 0,0, 0.1,10m, 0.3,1m, 0.5,2m, 0.6,20m)
.op
`

func TestVoltageControlledMultipleOperatingPoints(t *testing.T) {
	c := mustDeck(t, tunnel)
	sols := []float64{0.4 / 6, 0.26, 0.34}
	isSol := func(v float64) int {
		for k, s := range sols {
			if math.Abs(v-s) < 1e-6 {
				return k
			}
		}
		return -1
	}
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("a"); isSol(v) < 0 {
		t.Fatalf("cold start found v(a) = %g, not a solution", v)
	}
	for want, guess := range map[int]float64{0: 0.05, 1: 0.26, 2: 0.36} {
		r, err := OP(c, Options{NodeSet: map[string]float64{"a": guess, "in": 0.4}})
		if err != nil {
			t.Fatal(err)
		}
		v, _ := r.Voltage("a")
		if isSol(v) != want {
			t.Errorf("nodeset %g: v(a) = %.9g, want %g", guess, v, sols[want])
		}
		if rr := r.Ledger.RelativeResidual(); rr > 1e-12 {
			t.Errorf("ledger residual %g", rr)
		}
	}
	if _, err := OP(c, Options{NodeSet: map[string]float64{"zz": 1}}); err == nil {
		t.Error("unknown nodeset node accepted")
	}
}

// Small-signal gain on the negative-slope branch: g = −0.045 S, so
// v(a)/v(in) = (1/R)/(1/R + g) = 0.02/(0.02 − 0.045) = −0.8.
func TestNegativeSlopeACGain(t *testing.T) {
	c := mustDeck(t, strings.Replace(tunnel, "V1 in 0 0.4", "V1 in 0 DC 0.4 AC 1", 1))
	r, err := AC(c, ACSpec{Variation: "lin", Points: 1, FStart: 1e3, FStop: 1e3},
		Options{NodeSet: map[string]float64{"a": 0.26, "in": 0.4}})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := r.NodeVoltage("a")
	if cmplx.Abs(v[0]-(-0.8)) > 1e-12 {
		t.Fatalf("gain %v, want -0.8", v[0])
	}
}

// S-type (glow-lamp-like) curve, v = g(i), driven from 100 V through 10 kΩ.
// The only intersection is on the last segment: v = 60 + (10/9m)(i − 1m),
// so i = (100 − 60 + 10/9)/(10k + 10/9m). The drive is above the 90 V
// striking point, so the low-current branch has no solution and Newton,
// gmin stepping and source stepping all stall at the turning point; only
// pseudo-transient continuation reaches the glow branch.
func TestCurrentControlledConverges(t *testing.T) {
	c, err := circuit.NewBuilder("lamp").VoltageSource("V1", "in", "0", 100).
		Resistor("R1", "in", "a", 1e4).
		NonlinearCC("B1", "a", "0", []float64{0, 0, 1e-6, 90, 1e-3, 60, 1e-2, 70}).Build()
	if err != nil {
		t.Fatal(err)
	}
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := (100 - 60 + 10.0/9) / (1e4 + 10/9e-3)
	b1, _ := r.Element("B1")
	if math.Abs(b1.Current-want) > 1e-9*want {
		t.Fatalf("i = %.12g, want %.12g", b1.Current, want)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-12 {
		t.Fatalf("ledger residual %g", rr)
	}
	if !c.SelfControlled(c.Elements[c.Index("B1")]) {
		t.Fatal("series-sensed B1 not recognized as a nonlinear resistor")
	}
	if b1.ExteriorPower[exchange.Thermal] != -b1.PowerAbsorbed {
		t.Fatalf("B1 heat %g, absorbed %g", b1.ExteriorPower[exchange.Thermal], b1.PowerAbsorbed)
	}
}

// A curve through the 4th quadrant (i = −1 mA at v = 0) delivers power:
// with 1 kΩ across it, v = 1/3 V, and both the model-supplied power and the
// resistor's heat are v²/R = 1/9 mW.
func TestNonlinearResistorSupplyingIsModelSupplied(t *testing.T) {
	c := mustDeck(t, "q\nB1 a 0 I=pwl(V(a,0), -1,-3m, 0,-1m, 1,1m)\nR1 a 0 1k\n.op\n")
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("a"); !near(v, 1.0/3, 1e-12) {
		t.Fatalf("v(a) = %g", v)
	}
	l := r.Ledger
	if !near(l.ModelSupplied, 1e-3/9, 1e-15) || !near(l.Dissipation, 1e-3/9, 1e-15) || l.Inflow != 0 {
		t.Fatalf("ledger %+v", l)
	}
}

// A behavioral element controlled by another node is a controlled source:
// all of its power is model-supplied. Here a transconductance of 1 mS into
// 1 kΩ inverts the input.
func TestBehavioralControlledSource(t *testing.T) {
	c := mustDeck(t, "g\nV1 in 0 0.5\nRin in 0 1k\nB1 out 0 I=pwl(V(in), -1,-1m, 1,1m)\nRL out 0 1k\n.op\n")
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := r.Voltage("out"); !near(v, -0.5, 1e-12) {
		t.Fatalf("v(out) = %g", v)
	}
	if b, _ := r.Element("B1"); b.ModelSupplied != 0.25e-3 || b.ExteriorPower[exchange.Thermal] != 0 {
		t.Fatalf("B1 = %+v", b)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-12 {
		t.Fatalf("residual %g", rr)
	}
}

func TestOpenShort(t *testing.T) {
	c := mustDeck(t, "th\nV1 in 0 10\nR1 in out 1k\nR2 out 0 1k\n.op\n")
	p, err := OpenShort(c, "out", "0", Options{})
	if err != nil {
		t.Fatal(err)
	}
	if !near(p.Voc, 5, 1e-12) || !near(p.Isc, 10e-3, 1e-15) {
		t.Fatalf("%+v", p)
	}
	// Nonlinear: a diode clamps the open-circuit voltage, but shorted it
	// carries nothing, so Isc = 5 V / 1 kΩ exactly.
	c = mustDeck(t, hardDiodes)
	c = mustDeck(t, "d\nV1 in 0 5\nR1 in out 1k\nD1 out 0 dm\n.model dm D(IS=1e-14)\n.op\n")
	p, err = OpenShort(c, "out", "gnd", Options{})
	if err != nil {
		t.Fatal(err)
	}
	if p.Voc < 0.6 || p.Voc > 0.8 || !near(p.Isc, 5e-3, 1e-12) {
		t.Fatalf("%+v", p)
	}
	if _, err := OpenShort(c, "out", "out", Options{}); err == nil {
		t.Fatal("same node accepted")
	}
}

func TestBehavioralNetlistRoundTrip(t *testing.T) {
	src := "rt\nV1 in 0 1\nVs_b2 in m 0\nB1 in 0 I=pwl(V(in), 0,0, 1,1m)\nB2 m 0 V=pwl(I(Vs_b2), 0,0, 1m,1)\nR1 in 0 1k\n.op\n.end\n"
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb, ".op"); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(sb.String(), "B1 in 0 I=pwl(V(in), 0, 0, 1, 0.001)") ||
		!strings.Contains(sb.String(), "B2 m 0 V=pwl(I(vs_b2), 0, 0, 0.001, 1)") {
		t.Fatalf("netlist:\n%s", sb.String())
	}
	c2 := mustDeck(t, sb.String())
	r1, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	r2, _ := OP(c2, Options{})
	for i := range r1.NodeVoltages {
		if r1.NodeVoltages[i] != r2.NodeVoltages[i] {
			t.Fatalf("round trip changed the solution")
		}
	}
}

func BenchmarkTunnelDiodeOP(b *testing.B) {
	c := mustDeck(b, tunnel)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := OP(c, Options{}); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkGlowLampPseudoTransient(b *testing.B) {
	c, _ := circuit.NewBuilder("lamp").VoltageSource("V1", "in", "0", 100).Resistor("R1", "in", "a", 1e4).
		NonlinearCC("B1", "a", "0", []float64{0, 0, 1e-6, 90, 1e-3, 60, 1e-2, 70}).Build()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := OP(c, Options{}); err != nil {
			b.Fatal(err)
		}
	}
}
