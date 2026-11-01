// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// inverseCurve solves y = curve(x) for x on a rising piecewise-linear curve
// (end segments extended).
func inverseCurve(c []float64, y float64) float64 {
	n := len(c) / 2
	k := 1
	for k < n-1 && y > c[2*k+1] {
		k++
	}
	x0, y0, x1, y1 := c[2*(k-1)], c[2*(k-1)+1], c[2*k], c[2*k+1]
	return x0 + (y-y0)*(x1-x0)/(y1-y0)
}

func nearRel(t *testing.T, what string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol*math.Abs(want) {
		t.Errorf("%s = %.6g, want %.6g (relative tolerance %g)", what, got, want, tol)
	}
}

// A charge curve that is a straight line, Q = C·v, gives the same RC step
// response as the linear capacitor, and a straight flux curve the same RL
// response as the linear inductor.
func TestNonlinearReactiveLinearCurves(t *testing.T) {
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-9, 1e-9, 1, 2}}
	c, err := circuit.NewBuilder("rc").VoltageWave("V1", "in", "0", w).Resistor("R1", "in", "out", 1e3).
		NonlinearCapacitor("C1", "out", "0", []float64{-10, -1e-5, 10, 1e-5}, circuit.Waveform{}).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 5e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("out")
	checkSeries(t, "v(out)", wf, v, func(tt float64) float64 { return 1 - math.Exp(-tt/1e-3) }, 1)
	checkLedger(t, wf)

	w = circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 2, 0, 1e-9, 1e-9, 1, 2}}
	c, err = circuit.NewBuilder("rl").VoltageWave("V1", "in", "0", w).Resistor("R1", "in", "a", 10).
		NonlinearInductor("L1", "a", "0", []float64{-10, -0.1, 10, 0.1}, circuit.Waveform{}).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err = Tran(c, TranSpec{Step: 10e-6, Stop: 5e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	i, _ := wf.ElementCurrent("L1")
	checkSeries(t, "i(L1)", wf, i, func(tt float64) float64 { return 0.2 * (1 - math.Exp(-tt/1e-3)) }, 0.2)
	checkLedger(t, wf)
}

// C(t) = C0·(1 + ½ sin ωt) held at V by a source: i = V·dC/dt exactly. Over
// the first quarter period the source delivers V²·ΔC, half of it is stored
// (½V²·ΔC) and half leaves as mechanical work — the v·dC/dt term at work.
func TestVaryingCapacitorAtFixedVoltage(t *testing.T) {
	const V, C0, f = 2.0, 1e-6, 1e3
	scale := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.5, f}}
	c, err := circuit.NewBuilder("cv").VoltageSource("V1", "a", "0", V).VaryingCapacitor("C1", "a", "0", C0, scale).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 1e-6, Stop: 0.25 / f}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	i, _ := wf.ElementCurrent("C1")
	peak := V * C0 * 0.5 * 2 * math.Pi * f
	// t = 0 is the DC operating point, where the capacitor is open (as in
	// SPICE); the analytic current applies from the first step on.
	for k := 1; k < len(wf.Time); k++ {
		if want := peak * math.Cos(2*math.Pi*f*wf.Time[k]); math.Abs(i[k]-want) > tranFullScaleTol*peak {
			t.Fatalf("i(C1) at t = %g s: %g, want %g", wf.Time[k], i[k], want)
		}
	}
	checkLedger(t, wf)
	l := wf.Ledger
	dC := 0.5 * C0
	nearRel(t, "source energy", l.EnergyIn, V*V*dC, 2e-3)
	nearRel(t, "mechanical work out", l.EnergyOut, 0.5*V*V*dC, 2e-3)
	nearRel(t, "stored energy gain", l.StoredEnd-l.StoredStart, 0.5*V*V*dC, 1e-9)
}

// An isolated charged capacitor whose capacitance halves keeps its charge,
// so its voltage doubles (v = q/C) and its energy q²/2C doubles: the work
// of pulling the plates apart enters through the mechanical port.
func TestVaryingCapacitorKeepsCharge(t *testing.T) {
	const C0 = 1e-6
	scale := circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 1, 1e-3, 0.5}}
	b := circuit.NewBuilder("q")
	b.Add(circuit.Element{Name: "C1", Kind: circuit.NonlinearCapacitor, Value: C0, Scale: scale, IC: 1, HasIC: true}, "a", "0")
	c, err := b.BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 1e-3, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("a")
	checkSeries(t, "v(a)", wf, v, func(tt float64) float64 { return 1 / (1 - 0.5*tt/1e-3) }, 2)
	checkLedger(t, wf)
	l := wf.Ledger
	nearRel(t, "stored energy at the end", l.StoredEnd, C0, 1e-9)
	nearRel(t, "mechanical work in", l.EnergyIn, 0.5*C0, 2e-3)
}

// A constant current into a capacitor with a saturating charge curve: the
// charge rises linearly, q = I·t, and the voltage follows the inverse of the
// curve, v = Q⁻¹(I·t). The stored energy is ∫ v dq, all of it from the source.
func TestNonlinearCapacitorChargeCurve(t *testing.T) {
	curve := []float64{-1, -1e-6, 0, 0, 1, 1e-6, 2, 1.5e-6, 4, 2e-6}
	b := circuit.NewBuilder("nlc").CurrentSource("I1", "0", "a", 1e-3)
	b.Add(circuit.Element{Name: "C1", Kind: circuit.NonlinearCapacitor, Curve: curve, HasIC: true}, "a", "0")
	c, err := b.BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 2e-3, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("a")
	checkSeries(t, "v(a)", wf, v, func(tt float64) float64 { return inverseCurve(curve, 1e-3*tt) }, 4)
	checkLedger(t, wf)
	// ∫₀⁴ v dq = 4·Q(4) − ∫₀⁴ Q dv = 8e-6 − 5.25e-6.
	nearRel(t, "stored energy", wf.Ledger.StoredEnd, 2.75e-6, 1e-6)
	nearRel(t, "source energy", wf.Ledger.EnergyIn, 2.75e-6, 2e-3)
}

// A constant voltage across an inductor with a saturating flux curve: the
// flux rises linearly, λ = V·t, and the current is i = Λ⁻¹(V·t), climbing
// steeply once the core saturates.
func TestNonlinearInductorFluxCurve(t *testing.T) {
	curve := []float64{-1, -1e-3, 0, 0, 1, 1e-3, 2, 1.2e-3, 10, 2e-3}
	b := circuit.NewBuilder("nll").VoltageSource("V1", "a", "0", 1)
	b.Add(circuit.Element{Name: "L1", Kind: circuit.NonlinearInductor, Curve: curve, HasIC: true}, "a", "0")
	c, err := b.BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 2e-3, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	i, _ := wf.ElementCurrent("L1")
	checkSeries(t, "i(L1)", wf, i, func(tt float64) float64 { return inverseCurve(curve, tt) }, 10)
	checkLedger(t, wf)
	// ∫₀¹⁰ i dλ = 10·Λ(10) − ∫₀¹⁰ Λ di = 20e-3 − 14.4e-3.
	nearRel(t, "stored energy", wf.Ledger.StoredEnd, 5.6e-3, 1e-6)
	nearRel(t, "source energy", wf.Ledger.EnergyIn, 5.6e-3, 2e-3)
}

// A current circulating in a near-lossless loop keeps its flux when the
// inductance doubles, so the current halves (i = λ/L) and the energy λ²/2L
// halves: the difference leaves as mechanical work (the core is pulled in).
func TestVaryingInductorKeepsFlux(t *testing.T) {
	const L0 = 1e-3
	scale := circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 1, 1e-3, 2}}
	b := circuit.NewBuilder("flux").Resistor("R1", "a", "0", 1e-6)
	b.Add(circuit.Element{Name: "L1", Kind: circuit.NonlinearInductor, Value: L0, Scale: scale, IC: 1, HasIC: true}, "a", "0")
	c, err := b.BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 1e-3, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	i, _ := wf.ElementCurrent("L1")
	checkSeries(t, "i(L1)", wf, i, func(tt float64) float64 { return 1 / (1 + tt/1e-3) }, 1)
	checkLedger(t, wf)
	nearRel(t, "mechanical work out", wf.Ledger.EnergyOut, 0.25*L0, 2e-3)
}

// Parametric pumping: an LC tank whose capacitance is modulated at twice the
// resonant frequency, decreasing while the voltage peaks, grows in energy.
// Every joule of the growth, and of the loss in R, enters through the
// mechanical port; with the pump phase reversed the tank decays.
func TestParametricPumping(t *testing.T) {
	const L, C0 = 1e-3, 1e-6
	f0 := 1 / (2 * math.Pi * math.Sqrt(L*C0))
	run := func(phase float64) (stored0, stored1, mech float64) {
		scale := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 0.1, 2 * f0, 0, 0, phase}}
		b := circuit.NewBuilder("pump").Inductor("L1", "a", "b", L).Resistor("R1", "b", "0", 0.1)
		b.Add(circuit.Element{Name: "C1", Kind: circuit.NonlinearCapacitor, Value: C0, Scale: scale, IC: 1, HasIC: true}, "a", "0")
		c, err := b.BuildTransient()
		if err != nil {
			t.Fatal(err)
		}
		wf, err := Tran(c, TranSpec{Step: 1e-6, Stop: 10 / f0, UIC: true}, TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		checkLedger(t, wf)
		l := wf.Ledger
		nearRel(t, "mechanical energy against storage + loss", l.EnergyIn-l.EnergyOut,
			l.StoredEnd-l.StoredStart+l.EnergyDissipated, 1e-2)
		t.Logf("pump phase %g°: stored %.4g → %.4g J, mechanical net %.4g J, dissipated %.4g J, %d steps",
			phase, l.StoredStart, l.StoredEnd, l.EnergyIn-l.EnergyOut, l.EnergyDissipated, wf.Stats.Accepted)
		return l.StoredStart, l.StoredEnd, l.EnergyIn - l.EnergyOut
	}
	w0, w1, mech := run(180) // s = 1 − 0.1·sin(2ω₀t): C falls as |v| peaks
	if w1 < 5*w0 || mech <= 0 {
		t.Errorf("pumped tank: energy %.3g → %.3g J, mechanical %.3g J; want growth fed from the port", w0, w1, mech)
	}
	w0, w1, mech = run(0)
	if w1 > w0/5 || mech >= 0 {
		t.Errorf("reverse-phase tank: energy %.3g → %.3g J, mechanical %.3g J; want decay into the port", w0, w1, mech)
	}
}

// The small-signal capacitance of a biased nonlinear capacitor is the slope
// of its charge curve at the operating point.
func TestNonlinearCapacitorAC(t *testing.T) {
	curve := []float64{-1, -1e-6, 0, 0, 1, 1e-6, 2, 1.5e-6}
	mk := func(add func(*circuit.Builder)) *circuit.Circuit {
		b := circuit.NewBuilder("ac")
		b.Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, Value: 1.5, HasDC: true, ACMag: 1}, "in", "0")
		b.Resistor("R1", "in", "a", 1e3)
		add(b)
		c, err := b.Build()
		if err != nil {
			t.Fatal(err)
		}
		return c
	}
	nl := mk(func(b *circuit.Builder) { b.NonlinearCapacitor("C1", "a", "0", curve, circuit.Waveform{}) })
	lin := mk(func(b *circuit.Builder) { b.Capacitor("C1", "a", "0", 0.5e-6) })
	spec := ACSpec{Variation: "dec", Points: 5, FStart: 10, FStop: 1e5}
	r1, err := AC(nl, spec, Options{})
	if err != nil {
		t.Fatal(err)
	}
	r2, err := AC(lin, spec, Options{})
	if err != nil {
		t.Fatal(err)
	}
	v1, _ := r1.NodeVoltage("a")
	v2, _ := r2.NodeVoltage("a")
	for k := range v1 {
		if cmplx.Abs(v1[k]-v2[k]) > 1e-9*cmplx.Abs(v2[k]) {
			t.Fatalf("f = %g: v = %v, want %v", r1.Freq[k], v1[k], v2[k])
		}
	}
}

// A pumped tank with a saturating flux curve: a charge-curve capacitor,
// a flux-curve inductor and a varying scale in one transient.
func BenchmarkNonlinearReactiveTran(b *testing.B) {
	c := mustDeck(b, "pump\nV1 in 0 SIN(0 1 5k)\nR1 in a 10\n"+
		"CN a 0 Q=pwl(-2 -1.5u -1 -1u 0 0 1 1u 2 1.5u) SCALE=SIN(1 0.1 10k)\n"+
		"LN a 0 FLUX=pwl(-10 -2m -1 -1m 0 0 1 1m 10 2m)\n")
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := Tran(c, TranSpec{Step: 1e-6, Stop: 2e-3}, TranOptions{}); err != nil {
			b.Fatal(err)
		}
	}
}
