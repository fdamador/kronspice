// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/output"
)

// Tolerances for transient comparisons against exact analytic solutions,
// stated per CLAUDE.md §8 rule 5. With SPICE default tolerances (reltol 1e-3,
// trtol 7) the time-stepping error is a few parts in 10⁴ of full scale; 2e-3
// of full scale is the acceptance limit.
const tranFullScaleTol = 2e-3

func checkSeries(t *testing.T, name string, wf *output.Waveforms, y []float64, exact func(float64) float64, fullScale float64) {
	t.Helper()
	worst, at := 0.0, 0.0
	for k, tt := range wf.Time {
		if d := math.Abs(y[k] - exact(tt)); d > worst {
			worst, at = d, tt
		}
	}
	if worst > tranFullScaleTol*fullScale {
		t.Errorf("%s: max error %.3g at t = %.4g s exceeds %.3g", name, worst, at, tranFullScaleTol*fullScale)
	}
}

func checkLedger(t *testing.T, wf *output.Waveforms) {
	t.Helper()
	l := wf.Ledger
	if l.MaxRelativeResidual > 1e-9 {
		t.Errorf("per-point ledger residual %.3g", l.MaxRelativeResidual)
	}
	scale := math.Max(l.EnergyIn, math.Max(math.Abs(l.StoredEnd), math.Abs(l.StoredStart)))
	if scale > 0 && math.Abs(l.IntegrationMismatch) > 1e-2*scale {
		t.Errorf("integrated energy mismatch %.3g J of %.3g J", l.IntegrationMismatch, scale)
	}
}

func TestRCStep(t *testing.T) {
	// 1 V step (1 ns edge) into R = 1 kΩ, C = 1 µF: v_c = 1 − e^(−t/τ), τ = 1 ms.
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-9, 1e-9, 1, 2}}
	c, err := circuit.NewBuilder("rc").VoltageWave("V1", "in", "0", w).
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1e-6).Build()
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
	// Final stored energy ½CV² with V ≈ 1 − e^(−5).
	vf := 1 - math.Exp(-5)
	if got, want := wf.Ledger.StoredEnd, 0.5*1e-6*vf*vf; math.Abs(got-want) > 1e-3*want {
		t.Errorf("stored energy %.6g J, want %.6g J", got, want)
	}
	// Energy into an RC from a step: half stored, half dissipated (as t → ∞).
	if r := wf.Ledger.EnergyDissipated / wf.Ledger.StoredEnd; math.Abs(r-1) > 0.02 {
		t.Errorf("dissipated/stored = %.4f, want ≈ 1", r)
	}
	if wf.Stats.Accepted > 2000 {
		t.Errorf("too many steps: %+v", wf.Stats)
	}
}

func TestRLStep(t *testing.T) {
	// i_L = (V/R)(1 − e^(−tR/L)), V = 2 V, R = 10 Ω, L = 10 mH: τ = 1 ms.
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 2, 0, 1e-9, 1e-9, 1, 2}}
	c, err := circuit.NewBuilder("rl").VoltageWave("V1", "in", "0", w).
		Resistor("R1", "in", "a", 10).Inductor("L1", "a", "0", 10e-3).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 5e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	i, _ := wf.ElementCurrent("L1")
	checkSeries(t, "i(L1)", wf, i, func(tt float64) float64 { return 0.2 * (1 - math.Exp(-tt/1e-3)) }, 0.2)
	checkLedger(t, wf)
}

func TestRLCUnderdamped(t *testing.T) {
	// Series RLC step response, R = 10 Ω, L = 1 mH, C = 1 µF:
	// α = R/2L = 5000 s⁻¹, ω0 = 1/√(LC) = 31623 rad/s, ωd = √(ω0² − α²).
	// v_c = V[1 − e^(−αt)(cos ωd t + (α/ωd) sin ωd t)].
	R, L, C, V := 10.0, 1e-3, 1e-6, 1.0
	alpha := R / (2 * L)
	w0 := 1 / math.Sqrt(L*C)
	wd := math.Sqrt(w0*w0 - alpha*alpha)
	wave := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, V, 0, 1e-10, 1e-10, 1, 2}}
	c, err := circuit.NewBuilder("rlc").VoltageWave("V1", "in", "0", wave).
		Resistor("R1", "in", "a", R).Inductor("L1", "a", "b", L).Capacitor("C1", "b", "0", C).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 1e-6, Stop: 1e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("b")
	exact := func(tt float64) float64 {
		return V * (1 - math.Exp(-alpha*tt)*(math.Cos(wd*tt)+alpha/wd*math.Sin(wd*tt)))
	}
	checkSeries(t, "v(b)", wf, v, exact, 2*V) // peak ≈ 1.6 V
	checkLedger(t, wf)
}

func TestLCOscillatorConservesEnergy(t *testing.T) {
	// Lossless LC with C charged to 1 V (UIC): energy ½CV² must stay constant.
	// The trapezoidal rule has no numerical damping, so amplitude is preserved.
	L, C := 1e-3, 1e-6
	c, err := circuit.NewBuilder("lc").CapacitorIC("C1", "a", "0", C, 1).Inductor("L1", "a", "0", L).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	period := 2 * math.Pi * math.Sqrt(L*C)
	wf, err := Tran(c, TranSpec{Step: period / 200, Stop: 10 * period, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	e0 := 0.5 * C
	for k, e := range wf.Ledger.StoredEnergy {
		if math.Abs(e-e0) > 1e-3*e0 {
			t.Fatalf("stored energy %.6g J at t = %.4g s, want %.6g J", e, wf.Time[k], e0)
		}
	}
	if wf.Ledger.EnergyIn != 0 || wf.Ledger.EnergyDissipated != 0 {
		t.Errorf("closed LC should exchange nothing: %+v", wf.Ledger)
	}
	// The trapezoidal rule preserves amplitude exactly but has a small phase
	// error that accumulates: about 5e-4 of a period per period at SPICE's
	// default tolerances. Over the first period the waveform must meet the
	// usual 2e-3 full-scale limit; over 10 periods the stated limit is 1e-2.
	v, _ := wf.NodeVoltage("a")
	cosine := func(tt float64) float64 { return math.Cos(2 * math.Pi * tt / period) }
	worst1, worst10 := 0.0, 0.0
	for k, tt := range wf.Time {
		d := math.Abs(v[k] - cosine(tt))
		if tt <= period && d > worst1 {
			worst1 = d
		}
		if d > worst10 {
			worst10 = d
		}
	}
	if worst1 > tranFullScaleTol || worst10 > 1e-2 {
		t.Errorf("v(a) error: %.3g over 1 period, %.3g over 10 periods", worst1, worst10)
	}
}

func TestBackwardEulerDampsLC(t *testing.T) {
	// Backward Euler is known to add numerical damping; the ledger must show
	// the lost energy as integration mismatch, not as a physical loss channel.
	L, C := 1e-3, 1e-6
	c, _ := circuit.NewBuilder("lc").CapacitorIC("C1", "a", "0", C, 1).Inductor("L1", "a", "0", L).BuildTransient()
	period := 2 * math.Pi * math.Sqrt(L*C)
	wf, err := Tran(c, TranSpec{Step: period / 100, Stop: 5 * period, UIC: true}, TranOptions{Method: device.BackwardEuler})
	if err != nil {
		t.Fatal(err)
	}
	l := wf.Ledger
	if !(l.StoredEnd < l.StoredStart) {
		t.Fatalf("backward Euler should lose energy: %g → %g", l.StoredStart, l.StoredEnd)
	}
	if l.EnergyDissipated != 0 {
		t.Fatal("numerical damping must not appear as dissipation")
	}
	if l.MaxRelativeResidual > 1e-9 {
		t.Fatalf("per-point residual %g", l.MaxRelativeResidual)
	}
}

func TestSourceWaveforms(t *testing.T) {
	// Sources driving resistors: node voltages equal the waveform exactly.
	sin := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0.5, 2, 1e3, 0.2e-3, 0, 30}}
	pwl := circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-3, 1, 1.5e-3, 1, 2e-3, -1}}
	exp := circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{0, 1, 0.1e-3, 0.2e-3, 1e-3, 0.3e-3}}
	c, err := circuit.NewBuilder("src").
		VoltageWave("V1", "s", "0", sin).Resistor("R1", "s", "0", 1e3).
		VoltageWave("V2", "p", "0", pwl).Resistor("R2", "p", "0", 1e3).
		CurrentWave("I1", "0", "e", exp).Resistor("R3", "e", "0", 2).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 10e-6, Stop: 3e-3}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	vs, _ := wf.NodeVoltage("s")
	vp, _ := wf.NodeVoltage("p")
	ve, _ := wf.NodeVoltage("e")
	sinV := func(tt float64) float64 {
		if tt < 0.2e-3 {
			return 0.5 + 2*math.Sin(30*math.Pi/180)
		}
		return 0.5 + 2*math.Sin(2*math.Pi*1e3*(tt-0.2e-3)+30*math.Pi/180)
	}
	pwlV := func(tt float64) float64 {
		switch {
		case tt <= 1e-3:
			return tt / 1e-3
		case tt <= 1.5e-3:
			return 1
		case tt <= 2e-3:
			return 1 - 2*(tt-1.5e-3)/0.5e-3
		default:
			return -1
		}
	}
	expV := func(tt float64) float64 {
		i := 0.0
		if tt >= 0.1e-3 {
			i = 1 - math.Exp(-(tt-0.1e-3)/0.2e-3)
		}
		if tt >= 1e-3 {
			i += -(1 - math.Exp(-(tt-1e-3)/0.3e-3))
		}
		return 2 * i
	}
	for k, tt := range wf.Time {
		if math.Abs(vs[k]-sinV(tt)) > 1e-9 || math.Abs(vp[k]-pwlV(tt)) > 1e-9 || math.Abs(ve[k]-expV(tt)) > 1e-9 {
			t.Fatalf("t = %g: sin %g/%g pwl %g/%g exp %g/%g", tt, vs[k], sinV(tt), vp[k], pwlV(tt), ve[k], expV(tt))
		}
	}
	// PWL corners must be time points (breakpoints).
	for _, b := range []float64{1e-3, 1.5e-3, 2e-3} {
		found := false
		for _, tt := range wf.Time {
			if math.Abs(tt-b) < 1e-15 {
				found = true
			}
		}
		if !found {
			t.Errorf("breakpoint %g s not landed on", b)
		}
	}
	checkLedger(t, wf)
}

func TestPeriodicPulseBreakpoints(t *testing.T) {
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 1e-6, 1e-7, 1e-7, 4e-7, 1e-6}}
	c, _ := circuit.NewBuilder("p").VoltageWave("V1", "a", "0", w).Resistor("R1", "a", "0", 1).Build()
	wf, err := Tran(c, TranSpec{Step: 1e-8, Stop: 5e-6}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("a")
	// Value at the middle of the 3rd pulse's flat top: td + 2·per + tr + pw/2.
	if got := output.Interp(wf.Time, v, 1e-6+2e-6+1e-7+2e-7); math.Abs(got-1) > 1e-9 {
		t.Fatalf("pulse top = %g", got)
	}
	if got := output.Interp(wf.Time, v, 1e-6+2e-6+0.8e-6); math.Abs(got) > 1e-9 {
		t.Fatalf("pulse low = %g", got)
	}
}

func TestTranStartAndErrors(t *testing.T) {
	c, _ := circuit.NewBuilder("rc").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "b", 1).Capacitor("C1", "b", "0", 1).Build()
	wf, err := Tran(c, TranSpec{Step: 0.1, Stop: 2, Start: 1}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if wf.Time[0] < 1 {
		t.Fatalf("recorded before tstart: %g", wf.Time[0])
	}
	if _, err := Tran(c, TranSpec{Step: 0, Stop: 1}, TranOptions{}); !errors.Is(err, ErrBadTranSpec) {
		t.Fatalf("err = %v", err)
	}
	if _, err := Tran(c, TranSpec{Step: 1e-3, Stop: 1}, TranOptions{MaxSteps: 5}); !errors.Is(err, ErrTooManySteps) {
		t.Fatalf("err = %v, want ErrTooManySteps", err)
	}
	// A capacitor-only node needs UIC.
	fl := &circuit.Circuit{NodeNames: []string{"0", "a", "b"}, Elements: []circuit.Element{
		{Name: "V1", Kind: circuit.VoltageSource, Nodes: [2]int{1, 0}, Value: 1, HasDC: true},
		{Name: "R1", Kind: circuit.Resistor, Nodes: [2]int{1, 0}, Value: 1},
		{Name: "C1", Kind: circuit.Capacitor, Nodes: [2]int{2, 0}, Value: 1},
	}}
	if _, err := Tran(fl, TranSpec{Step: 0.1, Stop: 1}, TranOptions{}); !errors.Is(err, circuit.ErrFloatingNode) {
		t.Fatalf("err = %v, want ErrFloatingNode", err)
	}
}

func BenchmarkTranRC(b *testing.B) {
	w := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 1, 0, 1e-9, 1e-9, 1, 2}}
	c, _ := circuit.NewBuilder("rc").VoltageWave("V1", "in", "0", w).
		Resistor("R1", "in", "out", 1e3).Capacitor("C1", "out", "0", 1e-6).Build()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := Tran(c, TranSpec{Step: 10e-6, Stop: 5e-3}, TranOptions{}); err != nil {
			b.Fatal(err)
		}
	}
}

// Coupled inductors (K): a lossless pair, L1 charged with a current and L2
// open-circuited through a capacitor, exchanges energy with C and keeps
// ½L1i1² + ½L2i2² + M·i1·i2 + ½Cv² constant; with losses the ledger closes.
func TestMutualInductanceEnergy(t *testing.T) {
	L1, L2, k, C := 1e-3, 4e-3, 0.6, 1e-6
	c, err := circuit.NewBuilder("k").InductorIC("L1", "a", "0", L1, 0.1).Inductor("L2", "b", "0", L2).
		Capacitor("C1", "b", "0", C).CapacitorIC("C0", "a", "0", C, 0).Mutual("K1", "L1", "L2", k).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 1e-6, Stop: 2e-3, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	e0 := 0.5 * L1 * 0.1 * 0.1
	for kk, e := range wf.Ledger.StoredEnergy {
		if math.Abs(e-e0) > 2e-3*e0 {
			t.Fatalf("stored energy %.6g J at t = %.4g s, want %.6g J", e, wf.Time[kk], e0)
		}
	}
	// The secondary is driven only through M: without K it stays at zero.
	vb, _ := wf.NodeVoltage("b")
	peak := 0.0
	for _, v := range vb {
		peak = math.Max(peak, math.Abs(v))
	}
	if peak < 1 {
		t.Fatalf("secondary peak %g V: no coupling", peak)
	}
	checkLedger(t, wf)
}
