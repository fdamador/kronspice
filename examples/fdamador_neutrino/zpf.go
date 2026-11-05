// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 4: the zero-point field as the E-ring's inflow.
//
// The ring loses energy through its radiation resistance R_rad. By the
// fluctuation–dissipation theorem (H. B. Callen, T. A. Welton, Phys. Rev.
// 83, 34, 1951) the same resistance is a noise source: the field the ring
// radiates into drives it back. At ħω ≫ kT only the zero-point part is
// left, with one-sided voltage spectrum
//
//	S_V(ω) = ħ·ω·R_rad/π     (per unit angular frequency)
//
// For the 0.1 eV ring of Part 2 in the 2.7 K cosmic microwave background
// the thermal part is negligible (photon occupation e^(−ħω/kT) ≈ e^(−430)),
// so this is the whole inflow. It is booked on the magnetic-field channel
// (the ring couples to the radiation field as a magnetic dipole), and
// R_rad is a resistor whose power leaves on that same channel, so inflow
// and outflow sit in one ledger column.
//
// The loop is the E-ring of Part 1 in scaled units: total L = 1 H,
// C = 1 F (ω0 = 1 rad/s), R_rad = 1/Q, ħ replaced by 1 J·s. The circuit
// is linear, so every energy scales with ħ and is read as a multiple of
// ½ħω0 (as in examples/zpf_hydrogen, H24, whose noise construction this
// follows).

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// ZPFRing is the scaled E-ring in its radiation resistance's zero-point noise.
type ZPFRing struct {
	Q    float64 // radiation Q (R_rad = 1/Q)
	N    int     // ring sections
	Pts  int     // PWL points (power of two, ≤ circuit.MaxPWL)
	Dt   float64 // PWL spacing [s]
	Step float64 // kronspice time step [s]
	Band float64 // highest noise frequency [rad/s]
	Runs int
	Seed int64
	// Ramp is the raised-cosine switch-on of the noise, in units of Q
	// (seconds, with ω0 = 1). 0 starts the noise abruptly, which leaves a
	// frozen DC charge on the series capacitors (see Run); the default 5
	// ends inside the 10Q settling time that is discarded.
	Ramp float64
}

// DefaultZPFRing: the bench ring's radiation Q (40.68), 4 sections, 2¹⁵
// points 0.25 s apart (8,192 s, about 200 relaxation times), noise to 2ω0,
// 4 realizations.
func DefaultZPFRing() ZPFRing {
	return ZPFRing{Q: DefaultERing().QRad(), N: 4, Pts: 1 << 15, Dt: 0.25, Step: 0.02, Band: 2, Runs: 4, Seed: 1951, Ramp: 5}
}

// ZPFResult holds averages over the realizations (energies in units of ħ,
// so ½ħω0 = 0.5).
type ZPFResult struct {
	Energy, Spread    float64 // ⟨stored energy⟩ and its standard error
	ExpEnergy         float64 // exact expectation for the same noise modes
	Fluctuation       float64 // σ(E)/⟨E⟩ within a run, averaged over runs
	MeanMomentRatio   float64 // |⟨i⟩| / rms(i): the moment's time average
	ZPIn, RadOut      float64 // net energy from the zero-point source, out through R_rad (all runs)
	Stored            float64 // change of stored energy (all runs)
	Mismatch          float64 // worst |time-stepping mismatch| / R_rad energy
	MaxResidual       float64 // worst per-point ledger residual
	InflowOverOutflow float64 // ZPIn / RadOut
}

func (z ZPFRing) noise(rng *rand.Rand, r float64) []float64 {
	n := z.Pts
	dw := 2 * math.Pi / (float64(n) * z.Dt)
	spec := make([]complex128, n)
	for j := 1; j < n/2; j++ {
		w := float64(j) * dw
		if w > z.Band {
			break
		}
		amp := math.Sqrt(2 * w * r / math.Pi * dw) // ħ = 1
		spec[j] = cmplx.Rect(amp/2, 2*math.Pi*rng.Float64())
		spec[n-j] = cmplx.Conj(spec[j])
	}
	fft(spec)
	out := make([]float64, n)
	for i := range out {
		out[i] = real(spec[i])
	}
	return out
}

// Build makes the scaled ring with the noise waveform v as its source, in
// section 0 between p0 and q0 (the same layout as ERing.Build).
func (z ZPFRing) Build(v []float64) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("E-ring in zero-point noise (fdamador_neutrino)")
	s := strconv.Itoa
	nn := func(k int) string { return "n" + s(k%z.N) }
	cs, ls, rs := float64(z.N), 1/float64(z.N), 1/z.Q/float64(z.N)
	for k := 0; k < z.N; k++ {
		b.Capacitor("CS"+s(k), nn(k), "m"+s(k), cs)
		b.Inductor("LS"+s(k), "m"+s(k), "p"+s(k), ls)
		from := "p" + s(k)
		if k == 0 {
			from = "q0"
		}
		b.ResistorTo("RR"+s(k), from, nn(k+1), rs, exchange.Magnetic)
		b.Resistor("RB"+s(k), "m"+s(k), "0", 1e9)
		b.Resistor("RN"+s(k), nn(k), "0", 1e9)
	}
	pwl := make([]float64, 0, 2*len(v))
	for i, x := range v {
		pwl = append(pwl, float64(i)*z.Dt, x)
	}
	b.Add(circuit.Element{Name: "VZP", Kind: circuit.VoltageSource, Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: pwl},
		Channel: exchange.Magnetic, ChannelNote: "zero-point radiation field (Callen–Welton noise of R_rad)"}, "p0", "q0")
	return b.BuildTransient()
}

// elementEnergy is ∫ (absorbed power) dt of a named element.
func elementEnergy(time []float64, names []string, p [][]float64, name string) float64 {
	for i, n := range names {
		if n == name {
			var e float64
			for k := 1; k < len(time); k++ {
				e += 0.5 * (p[i][k] + p[i][k-1]) * (time[k] - time[k-1])
			}
			return e
		}
	}
	return math.NaN()
}

// Run simulates each realization in kronspice and averages the stored
// energy after 10 energy-relaxation times (Q/ω0).
func (z ZPFRing) Run() (ZPFResult, error) {
	rng := rand.New(rand.NewSource(z.Seed))
	var res ZPFResult
	var es []float64
	stop := float64(z.Pts-1) * z.Dt
	skip := 10 * z.Q
	r := 1 / z.Q
	for run := 0; run < z.Runs; run++ {
		v := z.noise(rng, r)
		// Ramp the noise on from zero over the first Ramp·Q seconds (raised
		// cosine), inside the discarded settling time. Without it the
		// t = 0 operating point charges the series capacitors to the
		// noise's first value through the 1 GΩ bleeders, an offset that
		// relaxes only over R_B·C ≈ 4×10⁹ s and sits in the stored energy
		// as a bias proportional to the noise power, i.e. to R_rad.
		ramp := z.Ramp * z.Q
		for i := range v {
			if t := float64(i) * z.Dt; t < ramp {
				v[i] *= 0.5 * (1 - math.Cos(math.Pi*t/ramp))
			}
		}
		ckt, err := z.Build(v)
		if err != nil {
			return res, err
		}
		wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: z.Step, Stop: stop, MaxStep: z.Step}, analysis.TranOptions{MaxSteps: 5_000_000})
		if err != nil {
			return res, fmt.Errorf("run %d: %w", run, err)
		}
		l := wf.Ledger
		var s, s2, tt, si, si2 float64
		il, _ := wf.ElementCurrent("LS0")
		for i := 1; i < len(wf.Time); i++ {
			if wf.Time[i] < skip {
				continue
			}
			dt := wf.Time[i] - wf.Time[i-1]
			e := l.StoredEnergy[i]
			s += e * dt
			s2 += e * e * dt
			si += il[i] * dt
			si2 += il[i] * il[i] * dt
			tt += dt
		}
		mean := s / tt
		es = append(es, mean)
		res.Fluctuation += math.Sqrt(s2/tt-mean*mean) / mean / float64(z.Runs)
		res.MeanMomentRatio = math.Max(res.MeanMomentRatio, math.Abs(si/tt)/math.Sqrt(si2/tt))
		zp := -elementEnergy(wf.Time, wf.ElementNames, wf.P, "VZP")
		var rad float64
		for k := 0; k < z.N; k++ {
			rad += elementEnergy(wf.Time, wf.ElementNames, wf.P, "RR"+strconv.Itoa(k))
		}
		res.ZPIn += zp
		res.RadOut += rad
		res.Stored += l.StoredEnd - l.StoredStart
		res.Mismatch = math.Max(res.Mismatch, math.Abs(l.IntegrationMismatch)/rad)
		res.MaxResidual = math.Max(res.MaxResidual, l.MaxRelativeResidual)
	}
	for _, e := range es {
		res.Energy += e / float64(len(es))
	}
	for _, e := range es {
		res.Spread += (e - res.Energy) * (e - res.Energy)
	}
	res.Spread = math.Sqrt(res.Spread/float64(len(es)-1)) / math.Sqrt(float64(len(es)))
	res.InflowOverOutflow = res.ZPIn / res.RadOut
	// Exact expectation for the same discrete modes. kronspice joins PWL
	// points by straight lines, which passes a mode ω with amplitude
	// sinc²(ω·Dt/2) (as in zpf_hydrogen).
	dw := 2 * math.Pi / (float64(z.Pts) * z.Dt)
	for j := 1; j < z.Pts/2; j++ {
		w := float64(j) * dw
		if w > z.Band {
			break
		}
		x := w * z.Dt / 2
		sv := w * r / math.Pi * dw * math.Pow(math.Sin(x)/x, 4)
		d := (1-w*w)*(1-w*w) + r*r*w*w
		res.ExpEnergy += sv * (w*w + 1) / 2 / d
	}
	return res, nil
}

// fft is an in-place radix-2 inverse FFT, Σ X_k e^{+2πikn/N} (no 1/N), so a
// spectrum of mode amplitudes becomes the time series (as in zpf_hydrogen).
func fft(a []complex128) {
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
	for size := 2; size <= n; size <<= 1 {
		w := cmplx.Rect(1, 2*math.Pi/float64(size))
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
