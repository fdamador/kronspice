// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"math/rand"
)

// The one-dimensional charged oscillator in zero-point radiation (Puthoff
// eqs. 4–16), in scaled units m = ω0 = ħ = 1. The radiation damping
// e²/(6πε0c³) is m·Γ; Γ is chosen far larger than the electron's
// (6.3×10⁻²⁴ s) so a run covers many relaxation times, and every result
// is compared with the formula for that Γ. The radiation reaction Γ·x⃛ is
// used in its reduced (Landau–Lifshitz) form −Γω0²·ẋ, equal at resonance:
//
//	ẍ + Γω0²·ẋ + ω0²·x = F(t)/m
//
// The force F = e·E_x of the zero-point field has the one-sided spectrum
// (eq. 1 with ⟨E_x²⟩ a third of ⟨E²⟩ = ρ/ε0)
//
//	S_F(ω) = e²·ħω³/(6π²ε0c³) = m·ħ·Γ·ω³/π
//
// generated with random phases (eq. 2) on a grid by an inverse FFT, kept
// below Band·ω0 (the full ω³ spectrum makes the kinetic energy diverge in
// any classical model; the cut is declared).

// Osc is one set of oscillator runs.
type Osc struct {
	Gamma float64 // radiation damping (scaled)
	Band  float64 // highest field frequency ÷ ω0
	N     int     // grid points (power of two)
	Dt    float64 // grid step
	Runs  int     // independent random phases
	Seed  int64
}

// OscResult holds the time averages over all runs.
type OscResult struct {
	Energy   float64 // ⟨½mẋ² + ½mω0²x²⟩
	Absorbed float64 // ⟨F·ẋ⟩
	Radiated float64 // ⟨mΓω0²·ẋ²⟩
	Spread   float64 // standard error of Energy across runs
	// Expected values for the same band-limited spectrum and damping
	// (exact sums over the grid's modes), and the narrowband formulas.
	ExpEnergy, ExpPower float64
	HalfHbarW, Eq16     float64 // ½ħω0 and ħΓω0³/2 (eq. 16 scaled)
	EnergyTrace         []float64
}

// noise returns one random-phase realization of F(t) on the grid, and the
// mode frequencies and amplitudes used.
func (o Osc) noise(rng *rand.Rand) []float64 {
	n := o.N
	T := float64(n) * o.Dt
	dw := 2 * math.Pi / T
	spec := make([]complex128, n)
	for k := 1; k < n/2; k++ {
		w := float64(k) * dw
		if w > o.Band {
			break
		}
		amp := math.Sqrt(2 * o.Gamma * w * w * w / math.Pi * dw) // S_F·dω, m = ħ = 1
		spec[k] = cmplx.Rect(amp/2, 2*math.Pi*rng.Float64())
		spec[n-k] = cmplx.Conj(spec[k])
	}
	fft(spec, true)
	out := make([]float64, n)
	for i := range out {
		out[i] = real(spec[i])
	}
	return out
}

// Run integrates the oscillator (RK4, force linearly interpolated) over
// each realization, discarding the first 10 relaxation times.
func (o Osc) Run() OscResult {
	rng := rand.New(rand.NewSource(o.Seed))
	var res OscResult
	tau := 1 / o.Gamma
	skip := int(10 * tau / o.Dt)
	var energies []float64
	for run := 0; run < o.Runs; run++ {
		f := o.noise(rng)
		x, v := 0.0, 0.0
		acc := func(x, v, F float64) float64 { return F - o.Gamma*v - x }
		var sE, sA, sR float64
		count := 0
		for i := 0; i+1 < len(f); i++ {
			F0, F1 := f[i], f[i+1]
			Fm := 0.5 * (F0 + F1)
			h := o.Dt
			k1x, k1v := v, acc(x, v, F0)
			k2x, k2v := v+h/2*k1v, acc(x+h/2*k1x, v+h/2*k1v, Fm)
			k3x, k3v := v+h/2*k2v, acc(x+h/2*k2x, v+h/2*k2v, Fm)
			k4x, k4v := v+h*k3v, acc(x+h*k3x, v+h*k3v, F1)
			x += h / 6 * (k1x + 2*k2x + 2*k3x + k4x)
			v += h / 6 * (k1v + 2*k2v + 2*k3v + k4v)
			if i < skip {
				continue
			}
			e := 0.5*v*v + 0.5*x*x
			sE += e
			sA += F1 * v
			sR += o.Gamma * v * v
			count++
			if run == 0 && i%20 == 0 {
				res.EnergyTrace = append(res.EnergyTrace, e)
			}
		}
		energies = append(energies, sE/float64(count))
		res.Absorbed += sA / float64(count) / float64(o.Runs)
		res.Radiated += sR / float64(count) / float64(o.Runs)
	}
	for _, e := range energies {
		res.Energy += e / float64(len(energies))
	}
	for _, e := range energies {
		res.Spread += (e - res.Energy) * (e - res.Energy)
	}
	res.Spread = math.Sqrt(res.Spread/float64(len(energies)-1)) / math.Sqrt(float64(len(energies)))
	// expected values over the same modes
	T := float64(o.N) * o.Dt
	dw := 2 * math.Pi / T
	for k := 1; k < o.N/2; k++ {
		w := float64(k) * dw
		if w > o.Band {
			break
		}
		sf := o.Gamma * w * w * w / math.Pi * dw
		d := (1-w*w)*(1-w*w) + o.Gamma*o.Gamma*w*w
		res.ExpEnergy += sf * (w*w + 1) / 2 / d
		res.ExpPower += sf * o.Gamma * w * w / d
	}
	res.HalfHbarW, res.Eq16 = 0.5, o.Gamma/2
	return res
}

// fft is an in-place radix-2 FFT; inverse=true computes Σ X_k e^{+2πikn/N}
// (no 1/N), so a spectrum of mode amplitudes becomes the time series.
func fft(a []complex128, inverse bool) {
	n := len(a)
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j ^= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
	sign := -1.0
	if inverse {
		sign = 1
	}
	for size := 2; size <= n; size <<= 1 {
		w := cmplx.Rect(1, sign*2*math.Pi/float64(size))
		for start := 0; start < n; start += size {
			wk := complex(1, 0)
			for k := 0; k < size/2; k++ {
				u := a[start+k]
				t := wk * a[start+k+size/2]
				a[start+k] = u + t
				a[start+k+size/2] = u - t
				wk *= w
			}
		}
	}
}
