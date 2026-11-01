// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Pulse test: the same 10 J dumped from a charged capacitor through the
// device, a 1 µH loop inductance and a 0.1 Ω load (declared). It answers
// "which is the most fast-acting for amps": peak current, rise time and
// how much of the pulse reaches the load.

const (
	pulseJ = 10.0
	pulseL = 1e-6
	pulseR = 0.1
)

// PulseResult is one device's pulse.
type PulseResult struct {
	Name       string
	V0, C      float64 // charge voltage [V], capacitance [F]
	Peak       float64 // peak current [A]
	Rise       float64 // 10–90 % rise time of the current [s]
	ToLoad     float64 // energy into the load [J]
	InDevice   float64 // energy lost in the device [J]
	Left       float64 // energy left in the capacitor [J]
	Mismatch   float64 // ledger mismatch (kronspice runs)
	Time, Amps []float64
}

func riseTime(t, i []float64, peak float64) float64 {
	t10, t90 := math.NaN(), math.NaN()
	for k := range i {
		a := math.Abs(i[k])
		if math.IsNaN(t10) && a >= 0.1*peak {
			t10 = t[k]
		}
		if math.IsNaN(t90) && a >= 0.9*peak {
			t90 = t[k]
			break
		}
	}
	return t90 - t10
}

// PulseCell discharges C (charged to v0) through the cell in kronspice.
func PulseCell(c Cell, v0 float64) (PulseResult, error) {
	cap := 2 * pulseJ / (v0 * v0)
	b := circuit.NewBuilder("cell pulse")
	b.CapacitorIC("CS", "a", "0", cap, v0)
	b.Inductor("LP", "a", "b", pulseL)
	c.Add(b, "X", "b", "c")
	b.Resistor("RL", "c", "0", pulseR)
	ckt, err := b.BuildTransient()
	if err != nil {
		return PulseResult{}, err
	}
	stop, step := 4e-3, 20e-9
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step, UIC: true}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return PulseResult{}, err
	}
	vc, _ := wv.NodeVoltage("c")
	va, _ := wv.NodeVoltage("a")
	r := PulseResult{Name: c.Name(), V0: v0, C: cap}
	for k := 1; k < len(wv.Time); k++ {
		i := vc[k] / pulseR
		r.Time = append(r.Time, wv.Time[k])
		r.Amps = append(r.Amps, i)
		r.Peak = math.Max(r.Peak, math.Abs(i))
		r.ToLoad += i * i * pulseR * (wv.Time[k] - wv.Time[k-1])
	}
	r.Rise = riseTime(r.Time, r.Amps, r.Peak)
	r.Left = 0.5 * cap * va[len(va)-1] * va[len(va)-1]
	r.InDevice = wv.Ledger.EnergyDissipated - r.ToLoad
	r.Mismatch = wv.Ledger.IntegrationMismatch
	return r, nil
}

// PulseGap discharges C (charged to v0, at or above the gap's forward
// breakdown) through a spark gap: an arc of Varc + Rarc that goes out at
// the first current zero after Hold, and re-strikes if the capacitor's
// reversed voltage exceeds the reverse breakdown (RK4, own loop).
func PulseGap(name string, g PolarGap, v0 float64) PulseResult {
	cap := 2 * pulseJ / (v0 * v0)
	r := PulseResult{Name: name, V0: v0, C: cap}

	tq := math.Pi * math.Sqrt(pulseL*cap) // half-period of the ring
	dt := tq / 4000
	vcap, i, t := v0, 0.0, 0.0
	on, sgn, tOn := v0 >= g.Vpos(), 1.0, 0.0

	deriv := func(vc, i, s float64) (float64, float64) {
		return -i / cap, (vc - s*g.Varc - (g.Rarc+pulseR)*i) / pulseL
	}
	for t < 20*tq {
		if on {
			k1v, k1i := deriv(vcap, i, sgn)
			k2v, k2i := deriv(vcap+dt/2*k1v, i+dt/2*k1i, sgn)
			k3v, k3i := deriv(vcap+dt/2*k2v, i+dt/2*k2i, sgn)
			k4v, k4i := deriv(vcap+dt*k3v, i+dt*k3i, sgn)
			ni := i + dt/6*(k1i+2*k2i+2*k3i+k4i)
			vcap += dt / 6 * (k1v + 2*k2v + 2*k3v + k4v)
			im := (i + ni) / 2
			r.ToLoad += im * im * pulseR * dt
			r.InDevice += (g.Varc*math.Abs(im) + g.Rarc*im*im) * dt
			if t-tOn >= g.Hold && ni*sgn <= g.Iext {
				on, ni = false, 0
			}
			i = ni
		} else {
			// off: the capacitor sits across the gap (leak through the water ignored over µs)
			if vcap >= g.Vpos() {
				on, sgn, tOn = true, 1, t
			} else if -vcap >= g.Vneg() {
				on, sgn, tOn = true, -1, t
			} else {
				break
			}
		}
		t += dt
		r.Time = append(r.Time, t)
		r.Amps = append(r.Amps, i)
		r.Peak = math.Max(r.Peak, math.Abs(i))
	}
	r.Rise = riseTime(r.Time, r.Amps, r.Peak)
	r.Left = 0.5 * cap * vcap * vcap
	r.Mismatch = math.Abs(pulseJ-r.ToLoad-r.InDevice-r.Left) / pulseJ
	return r
}

// MistGap is a 1 mm air gap in fog: PaschenAir(760 Torr, 1 mm) raised by
// the fog factor 29.4/21.7 (Liu et al. 2018, EXTRAPOLATED from 5 cm), as
// in examples/water_lamp; symmetric (no polarity), arc 30 V + 0.5 Ω.
func MistGap() PolarGap {
	return PolarGap{Fixed: 6822, NegRatio: 1, Varc: 30, Rarc: 0.5, Iext: 0.1, Hold: 0}
}

// Pulses runs every device with 10 J.
func Pulses() ([]PulseResult, error) {
	var out []PulseResult
	for _, a := range []Electrode{Pb, FeAgCoat, CuAgTips} {
		p, err := PulseCell(DefaultCell(a, Al), 100)
		if err != nil {
			return nil, err
		}
		out = append(out, p)
	}
	wg := DefaultGap()
	wg.Hold = 0
	wg.Tau = 1e-6                        // the capacitor's voltage is applied in about a microsecond
	v0 := math.Ceil(wg.Vpos()/1e3) * 1e3 // just above its breakdown
	out = append(out, PulseGap("water polarity gap, 3 mm", wg, v0))
	out = append(out, PulseGap("mist-filled air gap, 1 mm", MistGap(), 8000))
	return out, nil
}
