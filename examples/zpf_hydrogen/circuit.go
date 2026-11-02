// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The same oscillator as a circuit in kronspice: a series L–C tank whose
// resistance R plays the radiation resistance, driven by that resistor's
// own zero-point voltage noise. By the fluctuation–dissipation theorem
// (Callen & Welton, Phys. Rev. 83, 34, 1951 — Puthoff's footnote 13) a
// resistor at zero temperature has the one-sided voltage spectrum
//
//	S_V(ω) = ħ·ω·R/π     (per unit angular frequency; 2hfR per hertz)
//
// The charge q plays x, L plays m, R/L plays Γω0². Run in scaled units
// (L = 1 H, C = 1 F so ω0 = 1 rad/s, and ħ replaced by 1 J·s): the circuit
// is linear, so every energy simply scales with ħ, and the result is read
// as a multiple of ½ħω0.

// Tank is the circuit test.
type Tank struct {
	R    float64 // [Ω]
	N    int     // PWL points (≤ circuit.MaxPWL, power of two)
	Dt   float64 // PWL spacing [s]
	Step float64 // kronspice time step [s]
	Band float64 // highest noise frequency [rad/s]
	Runs int
	Seed int64
}

// TankResult holds the averages over the runs.
type TankResult struct {
	Energy, Spread float64 // ⟨½Li² + ½Cv²⟩ ÷ ħ, its standard error
	Supplied, Heat float64 // energy from the noise source and into R over all runs [J, scaled]
	Stored         float64 // change of stored energy, summed over runs
	Mismatch       float64 // worst ledger mismatch ÷ energy supplied
	ExpEnergy      float64 // exact expectation for the same modes
}

func (k Tank) noise(rng *rand.Rand) []float64 {
	n := k.N
	T := float64(n) * k.Dt
	dw := 2 * math.Pi / T
	spec := make([]complex128, n)
	for j := 1; j < n/2; j++ {
		w := float64(j) * dw
		if w > k.Band {
			break
		}
		amp := math.Sqrt(2 * w * k.R / math.Pi * dw)
		spec[j] = cmplx.Rect(amp/2, 2*math.Pi*rng.Float64())
		spec[n-j] = cmplx.Conj(spec[j])
	}
	fft(spec, true)
	out := make([]float64, n)
	for i := range out {
		out[i] = real(spec[i])
	}
	return out
}

// Run simulates each realization in kronspice and averages the stored
// energy after 10 relaxation times L/R.
func (k Tank) Run() (TankResult, error) {
	rng := rand.New(rand.NewSource(k.Seed))
	var res TankResult
	var es []float64
	T := float64(k.N-1) * k.Dt
	skip := 10 / k.R
	for run := 0; run < k.Runs; run++ {
		v := k.noise(rng)
		pwl := make([]float64, 0, 2*len(v))
		for i, x := range v {
			pwl = append(pwl, float64(i)*k.Dt, x)
		}
		b := circuit.NewBuilder("tank in zero-point noise").
			VoltageWave("VN", "a", "0", circuit.Waveform{Kind: circuit.WavePWL, Args: pwl}).
			Resistor("R", "a", "b", k.R).
			Inductor("L", "b", "c", 1).
			Capacitor("C", "c", "0", 1)
		ckt, err := b.BuildTransient()
		if err != nil {
			return res, err
		}
		wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: k.Step, Stop: T, MaxStep: k.Step, UIC: true}, analysis.TranOptions{MaxSteps: 5_000_000})
		if err != nil {
			return res, fmt.Errorf("run %d: %w", run, err)
		}
		il, _ := wf.ElementCurrent("L")
		vc, _ := wf.NodeVoltage("c")
		var s, tt float64
		for i := 1; i < len(wf.Time); i++ {
			if wf.Time[i] < skip {
				continue
			}
			dt := wf.Time[i] - wf.Time[i-1]
			s += (0.5*il[i]*il[i] + 0.5*vc[i]*vc[i]) * dt
			tt += dt
		}
		es = append(es, s/tt)
		l := wf.Ledger
		res.Supplied += l.EnergyIn - l.EnergyOut // the noise source both gives and takes
		res.Heat += l.EnergyDissipated
		res.Stored += l.StoredEnd - l.StoredStart
		res.Mismatch = math.Max(res.Mismatch, math.Abs(l.IntegrationMismatch)/l.EnergyDissipated)
	}
	for _, e := range es {
		res.Energy += e / float64(len(es))
	}
	for _, e := range es {
		res.Spread += (e - res.Energy) * (e - res.Energy)
	}
	res.Spread = math.Sqrt(res.Spread/float64(len(es)-1)) / math.Sqrt(float64(len(es)))
	dw := 2 * math.Pi / (float64(k.N) * k.Dt)
	for j := 1; j < k.N/2; j++ {
		w := float64(j) * dw
		if w > k.Band {
			break
		}
		sv := w * k.R / math.Pi * dw
		// kronspice joins the PWL points by straight lines, which passes
		// a mode of frequency ω with amplitude sinc²(ω·Dt/2).
		x := w * k.Dt / 2
		sinc := math.Sin(x) / x
		sv *= math.Pow(sinc, 4)
		d := (1-w*w)*(1-w*w) + k.R*k.R*w*w
		res.ExpEnergy += sv * (w*w + 1) / 2 / d
	}
	return res, nil
}
