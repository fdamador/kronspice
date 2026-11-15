// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The measurement as circuits. A trapped electron's axial oscillation
// drives a current e·ż/d through the endcap electrodes (d an effective
// distance set by the trap's geometry), and seen from those electrodes it
// is a series LC: L_e = m d²/e², C_e = 1/(ω_z² L_e) (Wineland & Dehmelt,
// J. Appl. Phys. 46, 919 (1975)). The detector is a parallel tank tuned to
// ω_z; the electron shorts it at ω_z, a dip of width γ_z = R_t e²/(m d²).

// SI constants (CODATA 2022).
const (
	Qe   = 1.602176634e-19
	Me   = 9.1093837139e-31
	MeC2 = 8.1871057880e-14 // J
	Hpl  = 6.62607015e-34
)

// Trap holds the run's trap (Fan et al. 2023: B = 5.3 T, ν̄_z ≈ 114 MHz;
// the effective electrode distance D is illustrative).
type Trap struct {
	B  float64 // T
	Fz float64 // axial frequency [Hz]
	D  float64 // effective electrode distance [m]
	Lt float64 // detector tank inductance [H]
	Qt float64 // detector tank quality factor
}

// DefaultTrap is the run's trap.
func DefaultTrap() Trap { return Trap{B: 5.3, Fz: 114e6, D: 3.5e-3, Lt: 100e-9, Qt: 1000} }

// ElectronLC is the electron's series equivalent.
func (t Trap) ElectronLC() (L, C float64) {
	L = Me * t.D * t.D / (Qe * Qe)
	w := 2 * math.Pi * t.Fz
	return L, 1 / (w * w * L)
}

// TankRC is the detector's capacitance and parallel resistance.
func (t Trap) TankRC() (C, R float64) {
	w := 2 * math.Pi * t.Fz
	return 1 / (w * w * t.Lt), t.Qt * w * t.Lt
}

// TransformerN is the ideal transformer ratio used to bring the electron's
// equivalent circuit to the detector's impedance scale.
const TransformerN = 1e5

// DipWidth is the expected full width of the dip [Hz]: γ_z/2π = R_t/(2π L_e).
func (t Trap) DipWidth() float64 {
	L, _ := t.ElectronLC()
	_, R := t.TankRC()
	return R / L / (2 * math.Pi)
}

// AxialDip sweeps the detector with and without the electron in kronspice
// AC and returns the dip's centre and full width at half depth [Hz], and
// the depth (|V| at the dip over |V| without the electron).
func (t Trap) AxialDip(span float64, points int) (fc, width, depth float64, freqs, with, without []float64, err error) {
	build := func(electron bool) (*circuit.Circuit, error) {
		b := circuit.NewBuilder("Penning trap axial detection")
		ct, rt := t.TankRC()
		b.Inductor("LT", "n", "0", t.Lt)
		b.Capacitor("CT", "n", "0", ct)
		b.Resistor("RT", "n", "0", rt)
		if electron {
			// The electron's series LC seen through an ideal 1:N transformer
			// (VCVS + CCCS): the circuit carries L_e/N², C_e·N², and the tank
			// sees N²·Z = Z_e. Same physics, with element values on the
			// tank's scale so the matrix stays well conditioned.
			le, ce := t.ElectronLC()
			n := TransformerN
			b.VCVS("ES", "s", "0", "n", "0", 1/n)
			b.Inductor("LE", "s", "m", le/(n*n))
			b.Capacitor("CE", "m", "k", ce*n*n)
			b.VoltageSource("VSENSE", "k", "0", 0)
			b.CCCS("FP", "n", "0", "VSENSE", 1/n)
		}
		b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "0", "n")
		return b.Build()
	}
	sweep := func(electron bool) ([]float64, []float64, error) {
		c, err := build(electron)
		if err != nil {
			return nil, nil, err
		}
		res, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: points, FStart: t.Fz - span, FStop: t.Fz + 1.013*span}, analysis.Options{})
		if err != nil {
			return nil, nil, err
		}
		v, _ := res.NodeVoltage("n")
		m := make([]float64, len(v))
		for i := range v {
			m[i] = cmplx.Abs(v[i])
		}
		return res.Freq, m, nil
	}
	freqs, with, err = sweep(true)
	if err != nil {
		return
	}
	_, without, err = sweep(false)
	if err != nil {
		return
	}
	k := 0
	for i := range with {
		if with[i]/without[i] < with[k]/without[k] {
			k = i
		}
	}
	depth = with[k] / without[k]
	half := 0.5 * (1 + depth)
	edge := func(dir int) float64 {
		for i := k; i+dir >= 0 && i+dir < len(with); i += dir {
			a, b := with[i]/without[i], with[i+dir]/without[i+dir]
			if b > half {
				return freqs[i] + (freqs[i+dir]-freqs[i])*(half-a)/(b-a)
			}
		}
		return math.NaN()
	}
	lo, hi := edge(-1), edge(1)
	return 0.5 * (lo + hi), hi - lo, depth, freqs, with, without, nil
}

// ---------------------------------------------------------------- the two clocks

// TrapFrequencies are the free cyclotron frequency and the three trap
// eigenfrequencies plus the anomaly, for an electron of given g in the
// trap (ideal quadrupole: ν̄_c ν̄_m = ν_z²/2, ν̄_c + ν̄_m = ν_c).
type TrapFrequencies struct {
	Nuc, NucBar, NuZ, NuM, NuA float64
}

// Frequencies computes them from B, ν_z and g.
func (t Trap) Frequencies(g float64) TrapFrequencies {
	nuc := Qe * t.B / (2 * math.Pi * Me)
	disc := math.Sqrt(nuc*nuc/4 - t.Fz*t.Fz/2)
	nucb := nuc/2 + disc
	num := nuc/2 - disc
	nus := g / 2 * nuc
	return TrapFrequencies{Nuc: nuc, NucBar: nucb, NuZ: t.Fz, NuM: num, NuA: nus - nucb}
}

// ReadTank finds a tank's resonance in kronspice AC: a coarse sweep of
// ±2 %, then a fine sweep of ±¼ linewidth round the coarse peak, with a
// parabola through the top three points.
func ReadTank(f0, q float64) (float64, error) {
	L := 1e-9
	w := 2 * math.Pi * f0
	C := 1 / (w * w * L)
	R := q * w * L
	b := circuit.NewBuilder("clock tank")
	b.Inductor("L", "n", "0", L)
	b.Capacitor("C", "n", "0", C)
	b.Resistor("R", "n", "0", R)
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "0", "n")
	c, err := b.Build()
	if err != nil {
		return 0, err
	}
	peak := func(fa, fb float64, n int) (float64, error) {
		res, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: n, FStart: fa, FStop: fb}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		v, _ := res.NodeVoltage("n")
		k := 1
		for i := 1; i+1 < len(v); i++ {
			if cmplx.Abs(v[i]) > cmplx.Abs(v[k]) {
				k = i
			}
		}
		// fit |V|² near the top: a Lorentzian's |V|⁻² is a parabola in f
		y0, y1, y2 := 1/sq(cmplx.Abs(v[k-1])), 1/sq(cmplx.Abs(v[k])), 1/sq(cmplx.Abs(v[k+1]))
		h := res.Freq[k+1] - res.Freq[k]
		return res.Freq[k] + h*(y0-y2)/(2*(y0-2*y1+y2)), nil
	}
	fc, err := peak(f0*0.98, f0*1.0201, 4001)
	if err != nil {
		return 0, err
	}
	bw := f0 / q
	return peak(fc-bw/4, fc+bw/4*1.0003, 401)
}

func sq(x float64) float64 { return x * x }

// GFromClocks is the measured g/2 from the trap's frequencies by the
// Brown–Gabrielse invariance theorem, ν_c² = ν̄_c² + ν̄_z² + ν̄_m², with
// g/2 = ν_s/ν_c = (ν̄_c + ν̄_a)/ν_c.
func GFromClocks(nucb, nuz, num, nua float64) float64 {
	nuc := math.Sqrt(nucb*nucb + nuz*nuz + num*num)
	return (nucb + nua) / nuc
}
