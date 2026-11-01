// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

func newTestDiode(params map[string]float64, area float64) *Diode {
	return NewDiode(circuit.Element{Name: "D1", Kind: circuit.Diode, Nodes: [2]int{1, 0}, Model: "d", Area: area},
		circuit.Model{Name: "d", Type: "d", Params: params})
}

func TestDiodeDefaultsAndArea(t *testing.T) {
	d := newTestDiode(map[string]float64{"rs": 10, "cjo": 2e-12}, 4)
	if d.P.IS != 4e-14 || d.P.RS != 2.5 || d.P.CJO != 8e-12 || d.P.N != 1 || d.P.VJ != 1 || d.P.FC != 0.5 {
		t.Fatalf("params = %+v", d.P)
	}
	// kT/q at 27 °C.
	if vt := d.vte; math.Abs(vt-0.0258649258) > 1e-10 {
		t.Fatalf("vt = %.9g", vt)
	}
}

// The current and its derivative are continuous where the equations switch
// (−3·vte and −BV), and gd is the derivative of id.
func TestDiodeCurrentContinuity(t *testing.T) {
	d := newTestDiode(map[string]float64{"is": 1e-12, "n": 1.5, "bv": 6}, 0)
	// At −XBV the SPICE3 equations meet to within the small cubic term of
	// the reverse formula, (3·vte/(e·XBV))³ ≈ 1e-7 here.
	for _, v := range []float64{-3 * d.vte, -d.xbv} {
		a, _ := d.current(v + 1e-12)
		b, _ := d.current(v - 1e-12)
		if math.Abs(a-b) > 1e-6*math.Max(math.Abs(a), d.P.IS) {
			t.Errorf("id jumps at %g V: %g vs %g", v, a, b)
		}
	}
	for _, v := range []float64{0.6, 0.1, -0.01, -0.5, -5.9, -6.2, -7} {
		h := 1e-6
		ip, _ := d.current(v + h)
		im, _ := d.current(v - h)
		_, gd := d.current(v)
		num := (ip - im) / (2 * h)
		if math.Abs(num-gd) > 1e-4*math.Abs(gd)+1e-20 {
			t.Errorf("gd(%g) = %g, numeric %g", v, gd, num)
		}
	}
}

// The capacitance is the derivative of the charge, on both sides of FC·VJ,
// and the charge is continuous there.
func TestDiodeChargeConsistency(t *testing.T) {
	d := newTestDiode(map[string]float64{"cjo": 5e-12, "vj": 0.8, "m": 0.33, "fc": 0.5, "tt": 5e-9}, 0)
	fcv := d.P.FC * d.P.VJ
	q := func(v float64) float64 {
		id, gd := d.current(v)
		qq, _ := d.charge(v, id, gd)
		return qq
	}
	if a, b := q(fcv-1e-12), q(fcv+1e-12); math.Abs(a-b) > 1e-9*math.Abs(a) {
		t.Errorf("charge jumps at FC·VJ: %g vs %g", a, b)
	}
	for _, v := range []float64{-3, -0.2, 0.2, fcv + 0.05, 0.65} {
		h := 1e-6
		id, gd := d.current(v)
		_, c := d.charge(v, id, gd)
		num := (q(v+h) - q(v-h)) / (2 * h)
		if math.Abs(num-c) > 1e-5*c {
			t.Errorf("C(%g) = %g, numeric %g", v, c, num)
		}
	}
}

func TestPnjlim(t *testing.T) {
	vt, vcrit := 0.025864186, 0.6
	// Small steps and steps below vcrit pass through unchanged.
	if v, lim := pnjlim(0.5, 0.49, vt, vcrit); v != 0.5 || lim {
		t.Errorf("below vcrit: %g %v", v, lim)
	}
	if v, lim := pnjlim(0.7, 0.69, vt, vcrit); v != 0.7 || lim {
		t.Errorf("small step: %g %v", v, lim)
	}
	// A large step from a forward-biased point is compressed logarithmically.
	v, lim := pnjlim(5, 0.7, vt, vcrit)
	if !lim || math.Abs(v-(0.7+vt*math.Log(1+(5-0.7)/vt))) > 1e-15 {
		t.Errorf("forward jump: %g %v", v, lim)
	}
	// From reverse bias the new value is vt·ln(v/vt).
	v, lim = pnjlim(5, -1, vt, vcrit)
	if !lim || math.Abs(v-vt*math.Log(5/vt)) > 1e-15 {
		t.Errorf("from reverse: %g %v", v, lim)
	}
}

// Load stamps the Newton companion: g = gd and ieq = id − gd·vd, so the
// linearized current at the expansion point equals the model current.
func TestDiodeLoadCompanion(t *testing.T) {
	d := newTestDiode(map[string]float64{"is": 1e-14}, 0)
	ctx := NewSetupContext(2)
	_ = d.Setup(ctx)
	m := &denseStub{n: 1, a: make([]float64, 1)}
	rhs := make([]float64, 1)
	x := []float64{0.62}
	_ = d.Load(&LoadContext{M: m, RHS: rhs, X: x, Gmin: 1e-12})
	id, gd := d.current(0.62)
	id += 1e-12 * 0.62
	gd += 1e-12
	if math.Abs(m.a[0]-gd) > 1e-15*gd {
		t.Errorf("g = %g, want %g", m.a[0], gd)
	}
	if ieq := -rhs[0]; math.Abs(ieq-(id-gd*0.62)) > 1e-12*math.Abs(ieq) {
		t.Errorf("ieq = %g, want %g", ieq, id-gd*0.62)
	}
	op := d.Operating(x)
	if math.Abs(op.Current-id) > 1e-12*id {
		t.Errorf("current = %g, want %g", op.Current, id)
	}
	if !d.ConvTest(x, 1e-3, 1e-12) {
		t.Error("ConvTest fails at the expansion point")
	}
	if d.ConvTest([]float64{0.70}, 1e-3, 1e-12) {
		t.Error("ConvTest passes far from the expansion point")
	}
	// Heat leaves through the thermal port; nothing is stored at DC.
	if got := d.ExteriorPower(exchange.Thermal, op); math.Abs(got+op.PowerAbsorbed) > 0 {
		t.Errorf("thermal = %g, want %g", got, -op.PowerAbsorbed)
	}
}

// XBV makes the breakdown current at −BV equal to IBV (for N = 1 the
// exponential term dominates, so |id(−BV)| ≈ IBV).
func TestDiodeBreakdownCurrent(t *testing.T) {
	d := newTestDiode(map[string]float64{"bv": 5.1, "ibv": 1e-3}, 0)
	id, _ := d.current(-5.1)
	if math.Abs(-id-1e-3) > 1e-3*1e-3 {
		t.Fatalf("id(−BV) = %g, want −1e-3", id)
	}
	// IBV below IS·BV/vt: XBV = BV.
	d = newTestDiode(map[string]float64{"bv": 5, "ibv": 1e-15}, 0)
	if d.xbv != 5 {
		t.Fatalf("xbv = %g, want 5", d.xbv)
	}
}

// The closed-form stored energy matches a fine numerical ∫₀^vd v·C(v) dv
// in every region: reverse, breakdown, depletion below and above FC·VJ,
// and strong forward bias where diffusion charge dominates.
func TestDiodeStoredEnergyAllRegions(t *testing.T) {
	d := newTestDiode(map[string]float64{"cjo": 2e-12, "vj": 0.7, "m": 0.4, "fc": 0.5, "tt": 5e-9, "bv": 6}, 0)
	numeric := func(vd float64) float64 {
		const n = 200000
		h := vd / n
		f := func(v float64) float64 {
			id, gd := d.current(v)
			_, c := d.charge(v, id, gd)
			return v * c
		}
		sum := f(0) + f(vd)
		for i := 1; i < n; i++ {
			w := 2.0
			if i%2 == 1 {
				w = 4
			}
			sum += w * f(float64(i)*h)
		}
		return sum * h / 3
	}
	for _, vd := range []float64{-8, -5, -0.05, 0.2, 0.5, 0.75} {
		d.vdSol = vd
		got, want := d.StoredEnergy(OperatingPoint{}), numeric(vd)
		if math.Abs(got-want) > 1e-8*math.Abs(want) {
			t.Errorf("vd = %g: stored %.12g, numeric %.12g", vd, got, want)
		}
	}
}

// StoredEnergy is ∫ v dq from 0 to vd.
func TestDiodeStoredEnergy(t *testing.T) {
	d := newTestDiode(map[string]float64{"cjo": 1e-12}, 0)
	ctx := NewSetupContext(2)
	_ = d.Setup(ctx)
	vd := -4.0
	d.Operating([]float64{vd})
	got := d.StoredEnergy(OperatingPoint{})
	// M = 0.5, VJ = 1: C = CJO/√(1 − v); ∫₀^vd v C dv in closed form.
	f := func(v float64) float64 { s := math.Sqrt(1 - v); return 1e-12 * (2.0/3*s*s*s - 2*s) }
	want := f(vd) - f(0)
	if math.Abs(got-want) > 1e-9*math.Abs(want) {
		t.Fatalf("stored = %g, want %g", got, want)
	}
}
