// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// Power relations of non-linear and parametric elements, from F. Kouřil
// and K. Vrba, Non-linear and Parametric Circuits, Ellis Horwood 1988,
// §4.5 pp. 133–142 (after Manley & Rowe and Pantell). The engine knows
// nothing of these relations; they must emerge from its transient
// solutions.

// spectrum returns the complex amplitudes at multiples h·f0 (h = 0…hmax)
// of a waveform over the last `periods` periods of 1/f0, by trapezoidal
// integration on a fine uniform resampling.
func spectrum(tm, x []float64, f0 float64, periods, hmax int) []complex128 {
	T := float64(periods) / f0
	t1 := tm[len(tm)-1]
	t0 := t1 - T
	const n = 200000
	dt := T / n
	out := make([]complex128, hmax+1)
	j := 0
	for s := 0; s <= n; s++ {
		t := t0 + float64(s)*dt
		for j < len(tm)-2 && tm[j+1] < t {
			j++
		}
		a := (t - tm[j]) / (tm[j+1] - tm[j])
		v := x[j] + a*(x[j+1]-x[j])
		w := 1.0
		if s == 0 || s == n {
			w = 0.5
		}
		for h := 0; h <= hmax; h++ {
			out[h] += complex(w*v*dt, 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(h)*f0*(t-t0)))
		}
	}
	for h := range out {
		f := 2 / T
		if h == 0 {
			f = 1 / T
		}
		out[h] *= complex(f, 0)
	}
	return out
}

// harmonicPower is the average power ½·Re(V·I*) at each harmonic (V·I at
// dc).
func harmonicPower(v, i []complex128) []float64 {
	p := make([]float64, len(v))
	for h := range v {
		if h == 0 {
			p[h] = real(v[0]) * real(i[0])
			continue
		}
		p[h] = 0.5 * real(v[h]*cmplx.Conj(i[h]))
	}
	return p
}

func mustTran(t *testing.T, c *circuit.Circuit, spec TranSpec) *output.Waveforms {
	t.Helper()
	wf, err := Tran(c, spec, TranOptions{Options: Options{RelTol: 1e-6}})
	if err != nil {
		t.Fatal(err)
	}
	return wf
}

// Manley–Rowe for a parametric (pumped) inductance (Kouřil & Vrba eqs.
// 4.95 and 4.101, Fig. 81b): an inductance L(t) = L0·(1 + m·cos ωp·t),
// pumped through its mechanical port, a signal source at ωs and an idler
// load at ωi = ωp − ωs. Over the combination frequencies k1·ωs + k2·ωp,
//
//	Σ k1·P(k1,k2)/(k1·ωs + k2·ωp) = 0,   Σ k2·P(k1,k2)/(…) + Pp/ωp = 0
//
// with P the power delivered to the element at each frequency and Pp the
// pump power; with only ωs and ωi present this is K&V eq. 4.101,
// −Ps : Pp : Pi = ωs : ωp : (ωp − ωs). The pumped element amplifies: it
// returns power to the signal source while drawing it from the pump.
func TestManleyRowePumpedInductor(t *testing.T) {
	const (
		fs, fp = 3e3, 11e3 // signal and pump; idler 8 kHz; common period 1 ms
		l0, m  = 1e-3, 0.2
		lf     = 10e-3 // filter inductance of the signal and idler branches
	)
	fi := fp - fs
	cf := func(f float64) float64 { return 1 / (math.Pow(2*math.Pi*f, 2) * lf) }
	s := circuit.NewMachineSpec(1)
	s.L0[0][0] = l0
	s.Harmonic(2).Cos[0][0] = m * l0 // L(θ) with θ = Ωt, Ω = ωp/2
	c, err := circuit.NewBuilder("pumped inductor").Domain(exchange.Mechanical, "w").
		MachineModel("pump", s).
		Machine("N1", "pump", [][2]string{{"a", "0"}}, "w", "0").
		VoltageSource("VW", "w", "0", math.Pi*fp).
		VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 1, fs}}).
		Resistor("RS", "s", "s1", 20).Inductor("LS", "s1", "s2", lf).Capacitor("CS", "s2", "a", cf(fs)).
		Inductor("LI", "a", "i1", lf).Capacitor("CI", "i1", "i2", cf(fi)).Resistor("RI", "i2", "0", 20).
		Build()
	if err != nil {
		t.Fatal(err)
	}
	wf := mustTran(t, c, TranSpec{Step: 1e-6, Stop: 40e-3, MaxStep: 1e-6})
	tm := wf.Time
	va, _ := wf.NodeVoltage("a")
	ia, _ := wf.ElementCurrent("N1")
	ivw, _ := wf.ElementCurrent("VW")
	const f0, periods, hmax = 1e3, 4, 60
	V, I := spectrum(tm, va, f0, periods, hmax), spectrum(tm, ia, f0, periods, hmax)
	P := harmonicPower(V, I)
	// Pump power into the element: the velocity source delivers −V·i.
	vw := make([]float64, len(tm))
	for k := range vw {
		vw[k] = -math.Pi * fp * ivw[k]
	}
	Pp := real(spectrum(tm, vw, f0, periods, 0)[0])
	// Assign each harmonic to its lowest-order combination k1·fs + k2·fp.
	var s1, s2, pmax float64
	for h := 1; h <= hmax; h++ {
		best, k1b, k2b := 99, 0, 0
		for k2 := -6; k2 <= 6; k2++ {
			r := float64(h)*f0 - float64(k2)*fp
			k1 := int(math.Round(r / fs))
			if math.Abs(float64(k1)*fs-r) > 1e-6 {
				continue
			}
			if o := abs(k1) + abs(k2); o < best {
				best, k1b, k2b = o, k1, k2
			}
		}
		if best == 99 {
			continue
		}
		w := 2 * math.Pi * float64(h) * f0
		s1 += float64(k1b) * P[h] / w
		s2 += float64(k2b) * P[h] / w
		pmax = math.Max(pmax, math.Abs(P[h]))
	}
	s2 += Pp / (2 * math.Pi * fp)
	hs, hi := int(fs/f0), int(fi/f0)
	ws, wi, wp := 2*math.Pi*fs, 2*math.Pi*fi, 2*math.Pi*fp
	ref := math.Abs(P[hs]) / ws
	if math.Abs(s1) > 0.01*ref || math.Abs(s2) > 0.01*ref {
		t.Errorf("Manley–Rowe sums %.3g and %.3g W·s, scale %.3g", s1, s2, ref)
	}
	// Eq. 4.101: the element returns power at ωs (P < 0 delivered to it),
	// takes Pp from the pump and delivers Pi at ωi, in ratio ωs : ωp : ωi.
	if !(P[hs] < 0 && P[hi] < 0 && Pp > 0) {
		t.Errorf("signs: P(ωs) %.3g, P(ωi) %.3g, Pp %.3g", P[hs], P[hi], Pp)
	}
	// The three-frequency proportion is exact only when no other
	// combination frequency carries power; here ωp ± ωs and the higher
	// products carry a few per cent, which the sums above include.
	a, b, d := -P[hs]/ws, Pp/wp, -P[hi]/wi
	if math.Abs(a-b) > 0.05*b || math.Abs(d-b) > 0.05*b {
		t.Errorf("−Ps/ωs = %.4g, Pp/ωp = %.4g, Pi/ωi = %.4g (K&V eq. 4.101: equal)", a, b, d)
	}
	t.Logf("Ps %.4g W at %g kHz, Pi %.4g W at %g kHz, Pp %.4g W; per rad/s: %.4g, %.4g, %.4g; sums %.2g, %.2g of scale",
		P[hs], fs/1e3, P[hi], fi/1e3, Pp, a, d, b, s1/ref, s2/ref)
}

// Pantell's relations for a non-linear resistor (Kouřil & Vrba eqs.
// 4.85–4.91): for one exciting frequency, Σ k²·P_k = ⟨G_d(v)·(dv/dα)²⟩ ≥ 0,
// where P_k is the power delivered to the resistor at k·ω and G_d its
// differential conductance. So a resistive frequency multiplier passes at
// most 1/c² of its input power to the c-th harmonic (eq. 4.91): 25 % to
// the second. A piecewise-linear nonlinear resistor (B element) driven
// through a source resistance is checked on both.
func TestPantellResistiveMultiplier(t *testing.T) {
	const f1 = 1e3
	// i = f(v): 1 mS below 0.5 V, 21 mS above (a hard knee).
	curve := []float64{-10, -10e-3, 0.5, 0.5e-3, 10, 0.5e-3 + 9.5*21e-3}
	c, err := circuit.NewBuilder("resistive doubler").
		VoltageWave("V1", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0.3, 1, f1}}).
		Resistor("RS", "s", "a", 50).
		NonlinearVC("B1", "a", "0", curve).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf := mustTran(t, c, TranSpec{Step: 1e-7, Stop: 3e-3, MaxStep: 1e-7})
	va, _ := wf.NodeVoltage("a")
	ib, _ := wf.ElementCurrent("B1")
	const hmax = 200
	V, I := spectrum(wf.Time, va, f1, 2, hmax), spectrum(wf.Time, ib, f1, 2, hmax)
	P := harmonicPower(V, I)
	lhs := 0.0
	for k := 1; k <= hmax; k++ {
		lhs += float64(k*k) * P[k]
	}
	// ⟨G_d(v)·(dv/dα)²⟩ with α = ωt, dv/dα = V′(t)/ω, over the same window.
	tm := wf.Time
	T := 2 / f1
	t0 := tm[len(tm)-1] - T
	rhs := 0.0
	for k := 1; k < len(tm); k++ {
		if tm[k-1] < t0 {
			continue
		}
		dt := tm[k] - tm[k-1]
		vm := 0.5 * (va[k] + va[k-1])
		g := 1e-3
		if vm > 0.5 {
			g = 21e-3
		}
		dv := (va[k] - va[k-1]) / dt / (2 * math.Pi * f1)
		rhs += g * dv * dv * dt
	}
	rhs /= T
	if math.Abs(lhs-rhs) > 0.01*rhs {
		t.Errorf("Σ k²·P_k = %.5g W, ⟨G_d·(dv/dα)²⟩ = %.5g W (K&V eq. 4.87)", lhs, rhs)
	}
	if !(P[1] > 0 && -P[2] <= P[1]/4 && -P[3] <= P[1]/9) {
		t.Errorf("P1 %.4g, P2 %.4g, P3 %.4g: multiplier bound 1/c² (eq. 4.91) violated", P[1], P[2], P[3])
	}
	t.Logf("Σ k²P_k %.5g W vs ⟨G_d(dv/dα)²⟩ %.5g W; P1 %.4g W in, P2 %.4g (%.1f %% ≤ 25 %%), P3 %.4g (%.1f %% ≤ 11 %%)",
		lhs, rhs, P[1], P[2], 100*-P[2]/P[1], P[3], 100*-P[3]/P[1])
}

func abs(k int) int {
	if k < 0 {
		return -k
	}
	return k
}
