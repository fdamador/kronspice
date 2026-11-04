// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Layer 3 of the Earth machine: the Earth–ionosphere cavity and its
// Schumann resonances (W. O. Schumann, Z. Naturforsch. 7a, 149, 1952;
// observed by M. Balser & C. A. Wagner, Nature 188, 638, 1960). The thin
// shell of air between the ground and the ionosphere is a transmission line
// along the colatitude θ, the network form of the cavity used by
// T. Madden & W. Thompson, "Low-frequency electromagnetic oscillations of the
// Earth–ionosphere cavity," Rev. Geophys. 3, 211, 1965.
//
// For the zonal (axisymmetric) TM₀ modes the voltage is V = E_r·h and the
// current I = H_φ·2πa·sinθ, the total current crossing the circle of
// colatitude θ. Per unit arc length a·dθ the line has
//
//	L' = μ0·h/(2πa·sinθ)      C' = ε0·2πa·sinθ/h,
//
// and the telegrapher's equations reduce to Legendre's equation with
// (ωa/c)² = n(n+1): f_n = c/(2πa)·√(n(n+1)) = 10.6, 18.4, 25.9 … Hz for a
// perfectly conducting cavity, whatever its height h.
//
// The Kron ladder: one band of the shell per node (capacitor to ground,
// the band's exact area), one inductor between neighbouring bands. The
// ionosphere is not a perfect wall: losses and the difference between the
// electric height h_e (C') and the magnetic height h_m (L') lower the
// frequencies (the two-height picture, e.g. V. C. Mushtak & E. R. Williams,
// J. Atmos. Sol.-Terr. Phys. 64, 1989, 2002). Here they enter as two fitted
// constants: the phase velocity v/c = √(h_e/h_m), and a uniform damping α
// (a shunt conductance G = 2α·C at every band), both fitted to the first mode
// alone. Modes 2–5 are then predictions.

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// MeasuredSchumann are commonly quoted average peak frequencies [Hz] of the
// first five modes (Balser & Wagner 1960; D. D. Sentman, "Schumann
// resonances," in Handbook of Atmospheric Electrodynamics, Vol. I, CRC, 1995).
var MeasuredSchumann = []float64{7.8, 14.3, 20.8, 27.3, 33.8}

// Cavity is the shell ladder.
type Cavity struct {
	Bands  int     // bands from pole to pole (odd: a band centred on the equator)
	Hm, He float64 // magnetic and electric heights [m]
	Alpha  float64 // uniform damping [1/s]; shunt G = 2α·C per band
}

// IdealCavity is the perfectly conducting cavity at 1° bands (h_e = h_m).
func IdealCavity() Cavity { return Cavity{Bands: 181, Hm: 90e3, He: 90e3} }

func (c Cavity) dtheta() float64 { return math.Pi / float64(c.Bands) }

// Elements returns each band's capacitance and the inductance between
// band k and k+1.
func (c Cavity) Elements() (C, L []float64) {
	d := c.dtheta()
	for k := 0; k < c.Bands; k++ {
		th0, th1 := float64(k)*d, float64(k+1)*d
		C = append(C, eps0*2*math.Pi*rEarth*rEarth*(math.Cos(th0)-math.Cos(th1))/c.He)
		if k < c.Bands-1 {
			L = append(L, mu0*c.Hm*d/(2*math.Pi*math.Sin(th1)))
		}
	}
	return
}

// LegendreFreq is the perfect-cavity mode frequency c/(2πa)·√(n(n+1)).
func LegendreFreq(n int) float64 {
	return cLight / (2 * math.Pi * rEarth) * math.Sqrt(float64(n*(n+1)))
}

// Modes returns the ladder's undamped mode frequencies [Hz], n = 1 … count,
// from the eigenvalues of C^(-1/2)·K·C^(-1/2) (K the inductor Laplacian),
// a symmetric tridiagonal matrix, by Sturm-sequence bisection. The n = 0
// eigenvalue (uniform charge, ω = 0) is skipped.
func (c Cavity) Modes(count int) []float64 {
	C, L := c.Elements()
	n := len(C)
	d := make([]float64, n)
	e := make([]float64, n-1)
	for k := 0; k < n; k++ {
		if k > 0 {
			d[k] += 1 / L[k-1] / C[k]
		}
		if k < n-1 {
			d[k] += 1 / L[k] / C[k]
			e[k] = -1 / (L[k] * math.Sqrt(C[k]*C[k+1]))
		}
	}
	below := func(x float64) int { // eigenvalues < x (Sturm count)
		cnt, q := 0, 1.0
		for k := 0; k < n; k++ {
			if k == 0 {
				q = d[0] - x
			} else {
				if q == 0 {
					q = 1e-300
				}
				q = d[k] - x - e[k-1]*e[k-1]/q
			}
			if q < 0 {
				cnt++
			}
		}
		return cnt
	}
	hi := 0.0
	for k := range d {
		hi = math.Max(hi, math.Abs(d[k])+2*math.Abs(e[min(k, n-2)]))
	}
	var out []float64
	for m := 1; m <= count; m++ { // the (m+1)-th smallest eigenvalue
		lo, up := 0.0, hi
		for it := 0; it < 200; it++ {
			mid := 0.5 * (lo + up)
			if below(mid) > m {
				up = mid
			} else {
				lo = mid
			}
		}
		out = append(out, math.Sqrt(0.5*(lo+up))/(2*math.Pi))
	}
	return out
}

// Circuit builds the ladder driven by a 1 A (AC) current source into the
// north-pole band, a vertical lightning stroke between the ground and the
// ionosphere.
func (c Cavity) Circuit() (*circuit.Circuit, error) {
	C, L := c.Elements()
	b := circuit.NewBuilder("Earth-ionosphere cavity")
	nd := func(k int) string { return fmt.Sprintf("b%d", k) }
	for k := range C {
		b.Capacitor(fmt.Sprintf("C%d", k), nd(k), "0", C[k])
		b.Resistor(fmt.Sprintf("G%d", k), nd(k), "0", 1/(2*c.Alpha*C[k]))
		if k < len(L) {
			b.Inductor(fmt.Sprintf("L%d", k), nd(k), nd(k+1), L[k])
		}
	}
	b.Add(circuit.Element{Name: "ISTROKE", Kind: circuit.CurrentSource, ACMag: 1}, "0", nd(0))
	return b.Build()
}

// Spectrum is |V| at the source band and at the equator from the engine's
// AC sweep, with the worst relative residual of its power ledger.
type Spectrum struct {
	F, Pole, Equator []float64
	LedgerRes        float64
}

// Sweep runs the AC analysis.
func (c Cavity) Sweep(f0, f1 float64, points int) (Spectrum, error) {
	var s Spectrum
	ck, err := c.Circuit()
	if err != nil {
		return s, err
	}
	ac, err := analysis.AC(ck, analysis.ACSpec{Variation: "lin", Points: points, FStart: f0, FStop: f1}, analysis.Options{})
	if err != nil {
		return s, err
	}
	vp, _ := ac.NodeVoltage("b0")
	ve, _ := ac.NodeVoltage(fmt.Sprintf("b%d", c.Bands/2))
	s.F = ac.Freq
	for k := range s.F {
		s.Pole = append(s.Pole, cmplx.Abs(vp[k]))
		s.Equator = append(s.Equator, cmplx.Abs(ve[k]))
	}
	s.LedgerRes = ac.Ledger.MaxRelativeResidual
	return s, nil
}

// Peaks returns the frequencies of the local maxima of y, refined by a
// parabola through the three samples.
func Peaks(f, y []float64) []float64 {
	var out []float64
	for k := 1; k+1 < len(y); k++ {
		if y[k] > y[k-1] && y[k] >= y[k+1] {
			den := y[k-1] - 2*y[k] + y[k+1]
			x := 0.0
			if den != 0 {
				x = 0.5 * (y[k-1] - y[k+1]) / den
			}
			out = append(out, f[k]+x*(f[k+1]-f[k]))
		}
	}
	return out
}

// FitLossy fits v/c and α to the first mode's damped frequency fd and Q
// (Q = ω0/2α): ω0 = 2π·fd/√(1 − 1/4Q²).
func FitLossy(fd, Q float64) Cavity {
	ideal := IdealCavity()
	w1 := 2 * math.Pi * ideal.Modes(1)[0]
	w0 := 2 * math.Pi * fd / math.Sqrt(1-1/(4*Q*Q))
	v := w0 / w1
	c := ideal
	c.He = v * v * c.Hm
	c.Alpha = w0 / (2 * Q)
	return c
}

func runCavity() error {
	fmt.Println("Layer 3 — the Earth–ionosphere cavity as a Kron transmission-line ladder (Schumann resonances)")
	ideal := IdealCavity()
	m := ideal.Modes(5)
	fmt.Printf("  perfect cavity, %d bands of 1°: ladder modes vs Legendre c/(2πa)·√(n(n+1))\n", ideal.Bands)
	for n := 1; n <= 5; n++ {
		fmt.Printf("    n=%d  ladder %.4f Hz  Legendre %.4f Hz  (rel. %.1g)\n", n, m[n-1], LegendreFreq(n), m[n-1]/LegendreFreq(n)-1)
	}
	// Twin: the engine's AC sweep of the nearly lossless ladder peaks at its eigenfrequencies.
	hq := ideal
	hq.Alpha = 2 * math.Pi * m[0] / (2 * 300) // Q1 = 300
	sp, err := hq.Sweep(5, 45, 4001)
	if err != nil {
		return err
	}
	pk := Peaks(sp.F, sp.Pole)
	fmt.Printf("  engine AC sweep, Q1 = 300, driven at the pole: peaks")
	for _, p := range pk {
		fmt.Printf(" %.3f", p)
	}
	fmt.Printf(" Hz; AC ledger residual ≤ %.1g\n", sp.LedgerRes)

	const q1 = 5.0 // first-mode Q, about 4–6 (Sentman 1995)
	lossy := FitLossy(MeasuredSchumann[0], q1)
	fmt.Printf("  lossy cavity fitted to mode 1 only (f = %.1f Hz, Q = %.0f): v/c = √(h_e/h_m) = %.4f (h_e = %.1f km at h_m = %.0f km), α = %.2f 1/s\n",
		MeasuredSchumann[0], q1, math.Sqrt(lossy.He/lossy.Hm), lossy.He/1e3, lossy.Hm/1e3, lossy.Alpha)
	sp, err = lossy.Sweep(2, 45, 4301)
	if err != nil {
		return err
	}
	pp, pe := Peaks(sp.F, sp.Pole), Peaks(sp.F, sp.Equator)
	um := lossy.Modes(5)
	fmt.Println("    n  model mode f [Hz]  Q    spectrum peak at pole  measured   model/measured   implied v/c")
	for n := 1; n <= 5; n++ {
		w0 := 2 * math.Pi * um[n-1]
		fd := math.Sqrt(w0*w0-lossy.Alpha*lossy.Alpha) / (2 * math.Pi)
		peak := math.NaN()
		if n-1 < len(pp) {
			peak = pp[n-1]
		}
		meas := MeasuredSchumann[n-1]
		fmt.Printf("    %d  %8.2f        %5.1f   %8.2f              %5.1f      %.3f            %.3f\n",
			n, fd, w0/(2*lossy.Alpha), peak, meas, fd/meas, meas/LegendreFreq(n))
	}
	fmt.Printf("  equator station (odd modes have a node there, P_n(0) = 0): peaks")
	for _, p := range pe {
		fmt.Printf(" %.2f", p)
	}
	fmt.Printf(" Hz; AC ledger residual ≤ %.1g\n", sp.LedgerRes)
	fmt.Println("  one constant velocity puts modes 2–5 increasingly low (the measured implied v/c rises with")
	fmt.Println("  frequency) and gives Q ∝ f; both need ionospheric heights that depend on frequency.")
	fmt.Println()
	return nil
}
