// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Layer 3 refinement: an ionosphere whose heights depend on frequency.
//
// In the two-scale-height picture (P. Greifinger & C. Greifinger, Radio Sci.
// 13, 831, 1978; D. D. Sentman, J. Atmos. Terr. Phys. 52, 35, 1990) the
// cavity's electric height h_e is where the ionospheric conductivity reaches
// ωε0, and its magnetic height h_m where it reaches 1/(4μ0ωζm²). With
// conductivity rising exponentially (scale heights ζe below, ζm above),
//
//	h_e(f) = h_e0 + ζe·ln(f/f0) − iπζe/2
//	h_m(f) = h_m0 − ζm·ln(f/f0) + iπζm/2,
//
// so h_e rises and h_m falls with frequency, and the imaginary parts are the
// losses. The line then has, per unit arc length, at each frequency
//
//	series:  R = ω·μ0·Im h_m / (2π sinθ),   L = μ0·Re h_m / (2π sinθ)
//	shunt:   C = ε0·w·Re h_e / |h_e|²,       G = ω·ε0·w·|Im h_e| / |h_e|²   (w = 2πa sinθ),
//
// from −iωμ0h_m and −iωε0/h_e (e^(−iωt)). A fixed RLC network cannot carry
// frequency-dependent heights, so the ladder is rebuilt at every frequency
// of the sweep and solved by the engine's AC analysis at that frequency.
// The modes solve n(n+1) = (ωa/c)²·h_m(ω)/h_e(ω) for complex ω.

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// TwoHeight is the frequency-dependent ionosphere.
type TwoHeight struct {
	He0, Hm0, ZetaE, ZetaM float64 // [m] at F0
	F0                     float64 // reference frequency [Hz]
}

// Heights returns the complex electric and magnetic heights at a
// (complex) frequency f.
func (t TwoHeight) Heights(f complex128) (he, hm complex128) {
	l := cmplx.Log(f / complex(t.F0, 0))
	he = complex(t.He0, 0) + complex(t.ZetaE, 0)*l - complex(0, math.Pi*t.ZetaE/2)
	hm = complex(t.Hm0, 0) - complex(t.ZetaM, 0)*l + complex(0, math.Pi*t.ZetaM/2)
	return
}

// Mode solves the zonal mode n for its complex frequency (e^(−iωt): a
// decaying mode has Im f < 0) and returns its frequency Re f and Q.
func (t TwoHeight) Mode(n int) (f, q float64) {
	nn := complex(float64(n*(n+1)), 0)
	F := func(f complex128) complex128 {
		he, hm := t.Heights(f)
		k := 2 * math.Pi * f * complex(rEarth/cLight, 0)
		return k*k*hm/he - nn
	}
	z := complex(0.75*LegendreFreq(n), -0.5)
	for it := 0; it < 100; it++ {
		d := 1e-6 * cmplx.Abs(z)
		g := (F(z+complex(d, 0)) - F(z-complex(d, 0))) / complex(2*d, 0)
		step := F(z) / g
		z -= step
		if cmplx.Abs(step) < 1e-12*cmplx.Abs(z) {
			break
		}
	}
	return real(z), real(z) / (2 * math.Abs(imag(z)))
}

// fitTwoHeight fits (h_e0, h_m0, ζe, ζm) at f0 = 8 Hz to the five measured
// frequencies and a first-mode Q, by Nelder–Mead on the relative errors.
func fitTwoHeight(q1 float64) (TwoHeight, float64) {
	model := func(x []float64) TwoHeight {
		return TwoHeight{He0: x[0] * 1e3, Hm0: x[1] * 1e3, ZetaE: x[2] * 1e3, ZetaM: x[3] * 1e3, F0: 8}
	}
	cost := func(x []float64) float64 {
		for _, v := range x {
			if v <= 0.1 {
				return 1e9
			}
		}
		t := model(x)
		var c float64
		for n := 1; n <= 5; n++ {
			f, q := t.Mode(n)
			if math.IsNaN(f) {
				return 1e9
			}
			e := (f - MeasuredSchumann[n-1]) / MeasuredSchumann[n-1]
			c += e * e
			if n == 1 {
				e = (q - q1) / q1
				c += e * e
			}
		}
		return c
	}
	x := nelderMead(cost, []float64{50, 90, 3, 4}, 0.5, 4000)
	return model(x), cost(x)
}

// nelderMead minimizes f from x0 (standard coefficients 1, 2, ½, ½).
func nelderMead(f func([]float64) float64, x0 []float64, step float64, iters int) []float64 {
	n := len(x0)
	pts := make([][]float64, n+1)
	val := make([]float64, n+1)
	for i := range pts {
		pts[i] = append([]float64(nil), x0...)
		if i > 0 {
			pts[i][i-1] += step * math.Max(1, math.Abs(x0[i-1]))
		}
		val[i] = f(pts[i])
	}
	for it := 0; it < iters; it++ {
		// order
		for i := 1; i <= n; i++ {
			for j := i; j > 0 && val[j] < val[j-1]; j-- {
				val[j], val[j-1] = val[j-1], val[j]
				pts[j], pts[j-1] = pts[j-1], pts[j]
			}
		}
		cen := make([]float64, n)
		for i := 0; i < n; i++ {
			for j := range cen {
				cen[j] += pts[i][j] / float64(n)
			}
		}
		at := func(a float64) []float64 {
			p := make([]float64, n)
			for j := range p {
				p[j] = cen[j] + a*(pts[n][j]-cen[j])
			}
			return p
		}
		xr := at(-1)
		fr := f(xr)
		switch {
		case fr < val[0]:
			xe := at(-2)
			if fe := f(xe); fe < fr {
				pts[n], val[n] = xe, fe
			} else {
				pts[n], val[n] = xr, fr
			}
		case fr < val[n-1]:
			pts[n], val[n] = xr, fr
		default:
			xc := at(0.5)
			if fc := f(xc); fc < val[n] {
				pts[n], val[n] = xc, fc
			} else {
				for i := 1; i <= n; i++ {
					for j := range pts[i] {
						pts[i][j] = pts[0][j] + 0.5*(pts[i][j]-pts[0][j])
					}
					val[i] = f(pts[i])
				}
			}
		}
	}
	best := 0
	for i := range val {
		if val[i] < val[best] {
			best = i
		}
	}
	return pts[best]
}

// circuitAt builds the cavity ladder for one frequency.
func (t TwoHeight) circuitAt(bands int, f float64) (*circuit.Circuit, error) {
	he, hm := t.Heights(complex(f, 0))
	w := 2 * math.Pi * f
	d := math.Pi / float64(bands)
	b := circuit.NewBuilder("Earth-ionosphere cavity, two-height")
	nd := func(k int) string { return fmt.Sprintf("b%d", k) }
	he2 := real(he)*real(he) + imag(he)*imag(he)
	for k := 0; k < bands; k++ {
		th0, th1 := float64(k)*d, float64(k+1)*d
		area := 2 * math.Pi * rEarth * rEarth * (math.Cos(th0) - math.Cos(th1))
		b.Capacitor(fmt.Sprintf("C%d", k), nd(k), "0", eps0*area*real(he)/he2)
		b.Resistor(fmt.Sprintf("G%d", k), nd(k), "0", he2/(w*eps0*area*math.Abs(imag(he))))
		if k < bands-1 {
			per := d / (2 * math.Pi * math.Sin(th1))
			mid := fmt.Sprintf("m%d", k)
			b.Resistor(fmt.Sprintf("RS%d", k), nd(k), mid, w*mu0*imag(hm)*per)
			b.Inductor(fmt.Sprintf("L%d", k), mid, nd(k+1), mu0*real(hm)*per)
		}
	}
	b.Add(circuit.Element{Name: "ISTROKE", Kind: circuit.CurrentSource, ACMag: 1}, "0", nd(0))
	return b.Build()
}

// Spectrum sweeps the rebuilt ladder: |V| at the pole and at the equator.
func (t TwoHeight) Spectrum(bands int, f0, f1, df float64) (Spectrum, error) {
	var s Spectrum
	for f := f0; f <= f1+1e-9; f += df {
		c, err := t.circuitAt(bands, f)
		if err != nil {
			return s, err
		}
		ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return s, err
		}
		vp, _ := ac.NodeVoltage("b0")
		ve, _ := ac.NodeVoltage(fmt.Sprintf("b%d", bands/2))
		s.F = append(s.F, f)
		s.Pole = append(s.Pole, cmplx.Abs(vp[0]))
		s.Equator = append(s.Equator, cmplx.Abs(ve[0]))
		s.LedgerRes = math.Max(s.LedgerRes, ac.Ledger.MaxRelativeResidual)
	}
	return s, nil
}

func runCavityTwoHeight() error {
	fmt.Println("Layer 3 refinement — frequency-dependent ionospheric heights (two-scale-height model)")
	const q1 = 5.0
	t, cost := fitTwoHeight(q1)
	fmt.Printf("  fitted to the five measured frequencies and Q1 = %.0f (rms error %.2g):\n", q1, math.Sqrt(cost/6))
	fmt.Printf("    at 8 Hz: h_e = %.1f km, h_m = %.1f km; scale heights ζe = %.2f km, ζm = %.2f km\n",
		t.He0/1e3, t.Hm0/1e3, t.ZetaE/1e3, t.ZetaM/1e3)
	sp, err := t.Spectrum(181, 3, 42, 0.05)
	if err != nil {
		return err
	}
	pp := Peaks(sp.F, sp.Pole)
	fmt.Println("    n   mode f [Hz]  Q     measured  model/measured  pole spectrum peak (rebuilt ladder, engine AC)")
	for n := 1; n <= 5; n++ {
		f, q := t.Mode(n)
		pk := math.NaN()
		if n-1 < len(pp) {
			pk = pp[n-1]
		}
		fmt.Printf("    %d   %7.2f      %4.1f   %5.1f     %.3f           %.2f\n", n, f, q, MeasuredSchumann[n-1], f/MeasuredSchumann[n-1], pk)
	}
	fmt.Printf("    engine AC ledger residual ≤ %.1g over %d rebuilt ladders\n", sp.LedgerRes, len(sp.F))
	fmt.Println("  the heights' frequency dependence lifts modes 2–5 onto the measurements and keeps Q from")
	fmt.Println("  growing ∝ f; Q2–Q5 are predictions (reported values stay near 4–6, Sentman 1995).")
	fmt.Println()
	return nil
}
