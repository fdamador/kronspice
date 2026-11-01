// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// Resonance in non-linear and parametric circuits and parametric
// amplifiers, from F. Kouřil and K. Vrba, Non-linear and Parametric
// Circuits, Ellis Horwood 1988, §10.2 pp. 371–393 and §6.4.4–6.4.5
// pp. 222–228. As in manley_rowe_test.go, the engine knows nothing of these
// results; they must emerge from its transient solutions.
//
// The engine's time-varying capacitor scales the charge, q = s(t)·Q(v), so
// its elastance is 1/s(t). With s = c0 + m·sin θ and c0 = √(1+m²), the mean
// elastance is exactly D0 and its first harmonic is −mEff·sin θ, with
// mEff = 2m/(1+c0) (the book's m in D = D0(1 + m·sin 2ω1t), eq. 10.71).

func pumpScale(m, f, phaseDeg float64) (circuit.Waveform, float64) {
	c0 := math.Sqrt(1 + m*m)
	return circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{c0, m, f, 0, 0, phaseDeg}}, 2 * m / (1 + c0)
}

// growthRate fits ln|peak| of x against time over [t0, t1] by least
// squares (the envelope's exponential rate, 1/s).
func growthRate(tm, x []float64, t0, t1 float64) float64 {
	var n, sx, sy, sxx, sxy float64
	for k := 1; k < len(x)-1; k++ {
		if tm[k] < t0 || tm[k] > t1 {
			continue
		}
		a := math.Abs(x[k])
		if a >= math.Abs(x[k-1]) && a > math.Abs(x[k+1]) && a > 0 {
			y := math.Log(a)
			n++
			sx += tm[k]
			sy += y
			sxx += tm[k] * tm[k]
			sxy += tm[k] * y
		}
	}
	return (n*sxy - sx*sy) / (n*sxx - sx*sx)
}

// Self-excitation of a linear parametric resonant circuit (K&V §10.2.2,
// eqs. 10.55–10.70, 10.101 and 10.103): a series R, L, C(t) circuit whose
// elastance is pumped at 2ω1 near twice its natural frequency ω0. The
// envelope grows at σ = −α + √((m·ω0/4)² − Δω²), α = R/2L, Δω = ω1 − ω0, so
// oscillation starts above m_c = 2·√(d² + (2Δω/ω0)²), d = 1/Q (10.103; 2/Q
// when tuned, 10.101).
func TestParametricSelfExcitation(t *testing.T) {
	const (
		l, c0 = 1e-3, 1e-9
		q     = 50.0
	)
	w0 := 1 / math.Sqrt(l*c0)
	r := 1 / (w0 * c0 * q)
	alpha := r / (2 * l)
	d := 1 / q
	for _, tc := range []struct {
		detune float64 // Δω/ω0
		ratio  float64 // m/m_c
		rate   bool    // check σ quantitatively (else its sign only)
	}{
		{0, 0.9, false}, {0, 1.1, false}, {0, 2, true},
		{d, 0.9, false}, {d, 1.1, false}, {d, 1.5, true},
	} {
		mc := 2 * math.Sqrt(d*d+4*tc.detune*tc.detune)
		w1 := w0 * (1 + tc.detune)
		f1 := w1 / (2 * math.Pi)
		// Find the SIN amplitude whose effective elastance index is m.
		mWant := tc.ratio * mc
		mSin := mWant
		for k := 0; k < 4; k++ {
			_, me := pumpScale(mSin, 0, 0)
			mSin *= mWant / me
		}
		sc, m := pumpScale(mSin, 2*f1, 0)
		c, err := circuit.NewBuilder("parametric RLC").
			InductorIC("L1", "a", "0", l, 1e-3).Resistor("R1", "a", "c", r).
			VaryingCapacitor("C1", "c", "0", c0, sc).BuildTransient()
		if err != nil {
			t.Fatal(err)
		}
		T := 2 * math.Pi / w0
		wf, err := Tran(c, TranSpec{Step: T / 80, Stop: 600 * T, MaxStep: T / 80, UIC: true}, TranOptions{Options: Options{RelTol: 1e-6}})
		if err != nil {
			t.Fatal(err)
		}
		il, _ := wf.ElementCurrent("L1")
		s := growthRate(wf.Time, il, 100*T, 600*T)
		g := (m*w0/4)*(m*w0/4) - (tc.detune*w0)*(tc.detune*w0)
		want := -alpha
		if g > 0 {
			want += math.Sqrt(g)
		}
		if tc.rate {
			if math.Abs(s-want) > 0.03*math.Sqrt(g) {
				t.Errorf("Δω/ω0 %g, m %.4g: σ = %.5g /s, want %.5g", tc.detune, m, s, want)
			}
		} else if (s > 0) != (tc.ratio > 1) {
			t.Errorf("Δω/ω0 %g, m/m_c %g: σ = %.4g /s has the wrong sign (m_c = %.4g, eq. 10.103)", tc.detune, tc.ratio, s, mc)
		}
		t.Logf("Δω/ω0 %g, m/m_c %g (m_c %.4g): σ/ω0 %.5f, predicted %.5f", tc.detune, tc.ratio, mc, s/w0, want/w0)
	}
}

// Resonance of a circuit with a capacitive diode (K&V §10.2.1, eqs.
// 10.31–10.38): C(v) = C0·(1 + √a·v)^(−1/2) around a reverse bias. The
// natural frequency falls with amplitude. The book's equivalent
// capacitance (10.37) gives ω_r/ω0 = (1 + 3aV1²/16)^(−1/2) ≈ 1 − (18/192)aV1².
// This curve's charge is q = (2C0/√a)(√(1+√a·v) − 1), so v(q) is exactly
// quadratic in q, v = q/C0 + √a·q²/(4C0²), and the Lindstedt–Poincaré
// expansion of L·q̈ + v(q) = 0 gives ω/ω0 = 1 − (5/192)·a·V1² + O(V1⁴).
// The engine must give the latter. (The book's fundamental of eq. 10.35
// keeps a2V1²/2 but drops the −a2V1²/4 from the cos 2ωt·cos ωt product,
// and its equivalent linearization neglects the second-harmonic voltage,
// which the circulating current of a free or high-Q tank makes ~√a·V1²/6,
// not small.)
func TestNonlinearTankBackbone(t *testing.T) {
	const (
		l, c0  = 1e-3, 1e-9
		sqrtA  = -1 / 4.7 // 1/(V0 − Vφ): V0 = −4 V, Vφ = 0.7 V
		i0     = 1e-3     // initial inductor current: V1 ≈ 1 V
		cycles = 20
	)
	a := sqrtA * sqrtA
	w0 := 1 / math.Sqrt(l*c0)
	T := 2 * math.Pi / w0
	var curve []float64
	for k := -500; k <= 400; k++ {
		v := float64(k) * 0.01
		curve = append(curve, v, 2*c0/sqrtA*(math.Sqrt(1+sqrtA*v)-1))
	}
	run := func(nonlinear bool) (w, v1 float64) {
		b := circuit.NewBuilder("varactor tank").InductorIC("L1", "a", "0", l, i0)
		if nonlinear {
			b.NonlinearCapacitor("C1", "a", "0", curve, circuit.Waveform{})
		} else {
			b.Capacitor("C1", "a", "0", c0)
		}
		c, err := b.Build()
		if err != nil {
			t.Fatal(err)
		}
		wf, err := Tran(c, TranSpec{Step: T / 2000, Stop: cycles * T, MaxStep: T / 2000, UIC: true},
			TranOptions{Options: Options{RelTol: 1e-7}})
		if err != nil {
			t.Fatal(err)
		}
		v, _ := wf.NodeVoltage("a")
		var up []float64
		for k := 1; k < len(v); k++ {
			if v[k-1] < 0 && v[k] >= 0 {
				up = append(up, wf.Time[k-1]+(wf.Time[k]-wf.Time[k-1])*(-v[k-1])/(v[k]-v[k-1]))
			}
		}
		n := len(up) - 1
		per := (up[n] - up[1]) / float64(n-1)
		// Fundamental over the last whole periods.
		var acc complex128
		var span float64
		for k := 1; k < len(v); k++ {
			tt := wf.Time[k]
			if tt < up[1] || tt > up[n] {
				continue
			}
			dt := tt - wf.Time[k-1]
			acc += complex(v[k]*dt, 0) * cmplx.Exp(complex(0, -2*math.Pi*tt/per))
			span += dt
		}
		return 2 * math.Pi / per, 2 * cmplx.Abs(acc) / span
	}
	wl, _ := run(false) // removes the integrator's own phase error
	wn, v1 := run(true)
	ratio := wn / wl
	x := a * v1 * v1
	coef := (1 - ratio) / x * 192
	lind := 1 - 5.0/192*x
	book := 1 / math.Sqrt(1+3.0/16*x)
	if math.Abs(coef-5) > 0.25 {
		t.Errorf("ω/ω0 = %.6f at V1 = %.3f V: shift %.3f/192·aV1², Lindstedt 5/192", ratio, v1, coef)
	}
	t.Logf("V1 %.3f V, aV1² %.4f: engine ω/ω0 %.6f (%.3f/192·aV1²); Lindstedt %.6f; K&V eq. 10.38 %.6f (18/192)",
		v1, x, ratio, coef, lind, book)
}

// Amplitude of parametric oscillation with a non-linear capacitor (K&V
// eqs. 10.71–10.92, Fig. 234): elastance D(q, t) = D0(1 + m·sin 2ω1t) +
// a·D0·q², tuned (x = ω1/ω0 = 1). Eq. 10.92 at x = 1 gives the steady
// normalized amplitude A² = a·Q1² from (3A²/4)² + d² = m²/4. Built as the
// pumped linear capacitor in series with a capacitor of cubic elastance
// v = a·D0·q³.
func TestParametricOscillatorAmplitude(t *testing.T) {
	const (
		l, c0 = 1e-3, 1e-9
		q     = 50.0
		a     = 1e18 // 1/C²
		mSin  = 0.1
	)
	d0 := 1 / c0
	w0 := 1 / math.Sqrt(l*c0)
	r := 1 / (w0 * c0 * q)
	d := 1 / q
	var curve []float64
	for k := -80; k <= 80; k++ {
		qq := float64(k) * 12.5e-12
		curve = append(curve, a*d0*qq*qq*qq, qq)
	}
	sc, m := pumpScale(mSin, w0/math.Pi, 0)
	c, err := circuit.NewBuilder("parametric oscillator").
		InductorIC("L1", "a", "0", l, 1e-5).Resistor("R1", "a", "c", r).
		VaryingCapacitor("CP", "c", "n", c0, sc).
		NonlinearCapacitor("CN", "n", "0", curve, circuit.Waveform{}).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	T := 2 * math.Pi / w0
	wf, err := Tran(c, TranSpec{Step: T / 100, Stop: 800 * T, MaxStep: T / 100, UIC: true}, TranOptions{Options: Options{RelTol: 1e-6}})
	if err != nil {
		t.Fatal(err)
	}
	il, _ := wf.ElementCurrent("L1")
	I := spectrum(wf.Time, il, w0/(2*math.Pi), 50, 1)
	q1 := cmplx.Abs(I[1]) / w0
	got := a * q1 * q1
	want := 4.0 / 3 * 0.5 * math.Sqrt(m*m-4*d*d)
	if math.Abs(got-want) > 0.03*want {
		t.Errorf("A² = %.5g, K&V eq. 10.92 at x = 1: %.5g", got, want)
	}
	t.Logf("m %.4g, d %.3g: A² = %.5g (Q1 %.4g C), eq. 10.92 %.5g", m, d, got, q1, want)
}

// Synchronous parametric amplifier with one resonant circuit (K&V §6.4.4,
// eqs. 6.23–6.27): a series R, L, C(t) circuit tuned to the signal ω1 and
// pumped at 2ω1. The pumped capacitor adds R_e = m·cos 2φ/(2ω1C0) (eq.
// 6.24), so the voltage gain V_L/V_s = Q/(1 ∓ m/m_c), m_c = 2/Q, depending
// on the signal's phase relative to the pump (eq. 6.27).
func TestSynchronousParametricAmplifier(t *testing.T) {
	const (
		l, c0 = 1e-3, 1e-9
		q     = 50.0
	)
	w0 := 1 / math.Sqrt(l*c0)
	f0 := w0 / (2 * math.Pi)
	r := 1 / (w0 * c0 * q)
	mc := 2 / q
	gain := func(phaseDeg float64) (float64, float64) {
		sc, m := pumpScale(0.5*mc, 2*f0, 0)
		c, err := circuit.NewBuilder("parametric amplifier").
			VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 1e-3, f0, 0, 0, phaseDeg}}).
			Resistor("R1", "s", "a", r).Inductor("L1", "a", "b", l).
			VaryingCapacitor("C1", "b", "0", c0, sc).Build()
		if err != nil {
			t.Fatal(err)
		}
		T := 1 / f0
		wf, err := Tran(c, TranSpec{Step: T / 80, Stop: 400 * T, MaxStep: T / 80}, TranOptions{Options: Options{RelTol: 1e-6}})
		if err != nil {
			t.Fatal(err)
		}
		va, _ := wf.NodeVoltage("a")
		vb, _ := wf.NodeVoltage("b")
		vl := make([]float64, len(va))
		for k := range vl {
			vl[k] = va[k] - vb[k]
		}
		return cmplx.Abs(spectrum(wf.Time, vl, f0, 40, 1)[1]) / 1e-3, m
	}
	// Engine elastance −m·sin 2ω0t, signal sin(ω0t + ψ): R_e ∝ cos 2ψ.
	gLow, m := gain(0)
	gHigh, _ := gain(90)
	wantHigh, wantLow := q/(1-m/mc), q/(1+m/mc)
	if math.Abs(gHigh-wantHigh) > 0.02*wantHigh || math.Abs(gLow-wantLow) > 0.02*wantLow {
		t.Errorf("gains %.4g and %.4g, K&V eq. 6.27: %.4g and %.4g", gHigh, gLow, wantHigh, wantLow)
	}
	t.Logf("m/m_c %.3f, Q %g: gain %.4g (eq. 6.27: %.4g) and %.4g (%.4g)", m/mc, q, gHigh, wantHigh, gLow, wantLow)
}

// Parametric amplifier with two resonant circuits (K&V §6.4.5, eqs.
// 6.30–6.40, Fig. 131): parallel tanks at ω1 and ω2 = ωs − ω1 coupled by
// C(t) = C0 − ΔC·cos ωs·t. Seen from the signal tank, the coupling and the
// idler tank are C0 in parallel with the negative conductance
// G_e = −¼·ω1·ω2·ΔC²·R2 (eq. 6.39), so the signal voltage rises by
// G1/(G1 + G_e); the idler voltage is V2 = ½·V1·ω2·ΔC·R2 (eq. 6.34); and the
// amplifier is stable while ¼·ω1·ω2·ΔC²·R2 < 1/R1 (eq. 6.40).
//
// The book assumes each tank is a short circuit at the other's frequency.
// Here an off-tune tank is about 5 Ω against the coupling capacitor's
// 145–330 Ω, which lowers the gain 2 of eq. 6.39 to 1.76. The exact
// steady state of this linear periodically varying circuit is found
// independently by harmonic balance over the sidebands ω1 + n·ωs
// (convMatrixGain); the engine is held to that within 2.5 %, the error of
// trapezoidal integration at its default tolerance on two high-Q tanks
// whose gain is amplified by 1/(1 − k), and to the book within 15 %.
func TestTwoCircuitParametricAmplifier(t *testing.T) {
	const (
		f1, f2 = 100e3, 230e3
		ctank  = 100e-9 // total tank capacitance, including C0
		c0     = 3e-9
		q      = 200.0
		k      = 0.5 // ¼·ω1·ω2·ΔC²·R1·R2: gain 1/(1 − k) = 2
	)
	w1, w2 := 2*math.Pi*f1, 2*math.Pi*f2
	r1, r2 := q/(w1*ctank), q/(w2*ctank)
	l1, l2 := 1/(w1*w1*ctank), 1/(w2*w2*ctank)
	dc := math.Sqrt(4 * k / (w1 * w2 * r1 * r2))
	run := func(pump bool) (v1, v2 float64) {
		mod := 0.0
		if pump {
			mod = dc / c0
		}
		c, err := circuit.NewBuilder("two-circuit paramp").
			CurrentWave("IS", "0", "a", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 1e-3, f1}}).
			Inductor("L1", "a", "0", l1).Capacitor("C1", "a", "0", ctank-c0).Resistor("R1", "a", "0", r1).
			VaryingCapacitor("C", "a", "b", c0, circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, mod, f1 + f2, 0, 0, 270}}).
			Inductor("L2", "b", "0", l2).Capacitor("C2", "b", "0", ctank-c0).Resistor("R2", "b", "0", r2).Build()
		if err != nil {
			t.Fatal(err)
		}
		const fc = 10e3 // common period of f1, f2 and f1 + f2
		h := 1 / (40 * (f1 + f2))
		wf, err := Tran(c, TranSpec{Step: h, Stop: 8e-3, MaxStep: h}, TranOptions{Options: Options{RelTol: 1e-6}})
		if err != nil {
			t.Fatal(err)
		}
		va, _ := wf.NodeVoltage("a")
		vb, _ := wf.NodeVoltage("b")
		A := spectrum(wf.Time, va, fc, 10, int(f1/fc))
		B := spectrum(wf.Time, vb, fc, 10, int(f2/fc))
		return cmplx.Abs(A[int(f1/fc)]), cmplx.Abs(B[int(f2/fc)])
	}
	v10, _ := run(false)
	v1, v2 := run(true)
	gain := v1 / v10
	tank := func(l, r float64) func(w float64) complex128 {
		return func(w float64) complex128 {
			return complex(1/r, w*(ctank-c0)-1/(w*l))
		}
	}
	hbGain, hbV2 := convMatrixGain(w1, w1+w2, c0, dc, tank(l1, r1), tank(l2, r2), 1e-3)
	book, bookV2 := 1/(1-k), 0.5*v1*w2*dc*r2
	if math.Abs(gain-hbGain) > 0.025*hbGain || math.Abs(v2-hbV2) > 0.025*hbV2 {
		t.Errorf("gain %.4g, V2 %.4g V; harmonic balance %.4g, %.4g V", gain, v2, hbGain, hbV2)
	}
	if math.Abs(gain-book) > 0.15*book || math.Abs(v2-bookV2) > 0.15*bookV2 {
		t.Errorf("gain %.4g (eq. 6.39: %.4g), V2 %.4g V (eq. 6.34: %.4g V)", gain, book, v2, bookV2)
	}
	t.Logf("ΔC %.4g F: gain %.4g (harmonic balance %.4g, eq. 6.39 %.4g); V2 %.4g V (harmonic balance %.4g, eq. 6.34 %.4g)",
		dc, gain, hbGain, book, v2, hbV2, bookV2)
}

// convMatrixGain solves the periodic steady state of two tanks (admittances
// y1, y2 from each node to ground) coupled by C(t) = C0 − ΔC·cos ωs·t, with
// a current is at ω1 into node 1, by harmonic balance over the frequencies
// ω1 + n·ωs, |n| ≤ 6 (negative frequencies are conjugate components). It
// returns the signal voltage's gain over ΔC = 0 and the idler amplitude at
// ωs − ω1.
func convMatrixGain(w1, ws, c0, dc float64, y1, y2 func(float64) complex128, is float64) (float64, float64) {
	const nh = 6
	m := 2*nh + 1
	solve := func(dc float64) []complex128 {
		n := 2 * m
		a := make([][]complex128, n)
		for i := range a {
			a[i] = make([]complex128, n+1)
		}
		for i := 0; i < m; i++ {
			w := w1 + float64(i-nh)*ws
			a[i][i] += y1(w)
			a[m+i][m+i] += y2(w)
			for j := 0; j < m; j++ {
				var c float64
				switch i - j {
				case 0:
					c = c0
				case 1, -1:
					c = -dc / 2
				default:
					continue
				}
				// Coupling charge c·(u_j) at ω_i, u = v1 − v2, current j·ω_i·charge.
				y := complex(0, w*c)
				a[i][j] += y
				a[i][m+j] -= y
				a[m+i][j] -= y
				a[m+i][m+j] += y
			}
		}
		a[nh][n] = complex(is, 0)
		for c := 0; c < n; c++ { // Gaussian elimination, partial pivoting
			p := c
			for r := c + 1; r < n; r++ {
				if cmplx.Abs(a[r][c]) > cmplx.Abs(a[p][c]) {
					p = r
				}
			}
			a[c], a[p] = a[p], a[c]
			for r := c + 1; r < n; r++ {
				f := a[r][c] / a[c][c]
				for k := c; k <= n; k++ {
					a[r][k] -= f * a[c][k]
				}
			}
		}
		x := make([]complex128, n)
		for r := n - 1; r >= 0; r-- {
			s := a[r][n]
			for k := r + 1; k < n; k++ {
				s -= a[r][k] * x[k]
			}
			x[r] = s / a[r][r]
		}
		return x
	}
	x0, x := solve(0), solve(dc)
	return cmplx.Abs(x[nh]) / cmplx.Abs(x0[nh]), cmplx.Abs(x[m+nh-1])
}
