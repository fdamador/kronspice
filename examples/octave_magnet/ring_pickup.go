// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Pickup is a small flat coil facing the ring's outer side: radius A, Turns
// turns, axis radial, centre at radius Rc in the mid-plane. Its resistance
// and inductance are declared: 200 turns of 0.1 mm copper on a 6 mm former,
// R = ρ_Cu·N·2πA/(πd²/4); L from Wheeler's multilayer formula (Proc. IRE
// 16, 1398, 1928) for a 2 mm long, 1 mm deep winding.
type Pickup struct {
	A, Rc float64
	Turns int
	R, L  float64
	Load  float64
	RPM   float64
}

// DefaultPickup sits 1.5 mm outside the ring; the ring turns at 3000 rpm.
func DefaultPickup(r Ring) Pickup {
	const rhoCu, d = 1.72e-8, 0.1e-3
	n := 200
	a := 3e-3
	res := rhoCu * float64(n) * 2 * math.Pi * a / (math.Pi * d * d / 4)
	// Wheeler multilayer: L[µH] = 0.8 a²N²/(6a + 9b + 10c), lengths in inches.
	in := 1 / 0.0254
	l := 0.8 * math.Pow(a*in, 2) * float64(n*n) / (6*a*in + 9*2e-3*in + 10*1e-3*in) * 1e-6
	return Pickup{A: a, Rc: r.Ro + 1.5e-3 + 0, Turns: n, R: res, L: l, Load: 8.2, RPM: 3000}
}

// Flux is the flux through one turn of the coil when it faces the ring at
// azimuth theta (polar Gauss rule over the disc).
func (p Pickup) Flux(cs []Charge, theta float64) float64 {
	gx, gw := gauss(3)
	nr, na := 4, 24
	c, s := math.Cos(theta), math.Sin(theta)
	var phi float64
	for i := 0; i < nr; i++ {
		for a := range gx {
			u := (float64(i) + gx[a]) / float64(nr) // r/A
			rr := u * p.A
			w := gw[a] / float64(nr) * p.A * rr // dr · r
			for j := 0; j < na; j++ {
				t := 2 * math.Pi * (float64(j) + 0.5) / float64(na)
				y, z := rr*math.Cos(t), rr*math.Sin(t)
				// local frame: x along r̂(θ), y along φ̂(θ)
				pt := [3]float64{p.Rc*c - y*s, p.Rc*s + y*c, z}
				b := BCharges(cs, pt)
				phi += (b[0]*c + b[1]*s) * w * 2 * math.Pi / float64(na)
			}
		}
	}
	return phi
}

// PickupResult is the coil's run in kronspice.
type PickupResult struct {
	Theta, Flux, EMF []float64 // over one revolution
	Cycles           int       // electrical cycles per revolution
	PeakEMF, Vrms    float64   // open-circuit peak, load rms [V]
	Power            float64   // mean load power [W]
	Ein, Ediss, dE   float64   // ledger over the run [J]
	Residual         float64   // kronspice ledger, max relative residual
}

// Run computes the flux round one revolution, turns it into the EMF
// e = N·Ω·dΦ/dθ (the ring turning at Ω past the fixed coil), and drives
// coil resistance, coil inductance and the load with it in kronspice for two
// revolutions. The EMF source stands in for the shaft: every joule the load
// gets is work done turning the ring.
func (p Pickup) Run(cs []Charge, n int) (PickupResult, error) {
	var res PickupResult
	res.Theta, res.Flux, res.EMF = make([]float64, n), make([]float64, n), make([]float64, n)
	parallel(n, func(i int) {
		res.Theta[i] = 2 * math.Pi * float64(i) / float64(n)
		res.Flux[i] = p.Flux(cs, res.Theta[i])
	})
	omega := p.RPM / 60 * 2 * math.Pi
	dth := 2 * math.Pi / float64(n)
	for i := range res.EMF {
		d := (res.Flux[(i+1)%n] - res.Flux[(i+n-1)%n]) / (2 * dth)
		res.EMF[i] = float64(p.Turns) * omega * d
		res.PeakEMF = math.Max(res.PeakEMF, math.Abs(res.EMF[i]))
	}
	res.Cycles = CyclicSignChanges(res.EMF) / 2
	trev := 2 * math.Pi / omega
	var pwl []float64
	for rev := 0; rev < 2; rev++ {
		for i := 0; i < n; i++ {
			// the ring turns by +Ωt, so the coil sees angle −Ωt
			k := (n - i) % n
			pwl = append(pwl, float64(rev)*trev+float64(i)*trev/float64(n), res.EMF[k])
		}
	}
	pwl = append(pwl, 2*trev, res.EMF[0])
	b := circuit.NewBuilder("pickup coil beside the 12-section ring")
	b.VoltageWave("EP", "e", "0", circuit.Waveform{Kind: circuit.WavePWL, Args: pwl})
	b.Resistor("RC", "e", "a", p.R)
	b.Inductor("LC", "a", "b", p.L)
	b.Resistor("RL", "b", "0", p.Load)
	ckt, err := b.BuildTransient()
	if err != nil {
		return res, err
	}
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: trev / float64(4*n), Stop: 2 * trev}, analysis.TranOptions{})
	if err != nil {
		return res, err
	}
	v, ok := w.NodeVoltage("b")
	if !ok {
		return res, fmt.Errorf("no node b")
	}
	// rms and power over the second revolution
	var e2, tt float64
	for i := 1; i < len(w.Time); i++ {
		if w.Time[i-1] < trev {
			continue
		}
		dt := w.Time[i] - w.Time[i-1]
		e2 += 0.5 * (v[i]*v[i] + v[i-1]*v[i-1]) * dt
		tt += dt
	}
	res.Vrms = math.Sqrt(e2 / tt)
	res.Power = res.Vrms * res.Vrms / p.Load
	L := w.Ledger
	res.Ein, res.Ediss = L.EnergyIn, L.EnergyDissipated
	res.dE = L.StoredEnd - L.StoredStart
	res.Residual = L.MaxRelativeResidual
	return res, nil
}
