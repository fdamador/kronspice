// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// A collapsing electron ring: relativistic electrons circling the axis of
// a solenoid are compressed as the solenoid's field rises — the electron
// ring compressor (V. I. Veksler, CERN Symposium 1956; the Dubna and Berkeley
// electron-ring accelerators). The rising field B(t) induces the azimuthal
// field E_φ = −(r/2)·dB/dt (Faraday), which speeds the electrons up; the
// stronger field bends them tighter. Two invariants: the electrons'
// azimuthal angle is a cyclic coordinate, so the canonical angular
// momentum P_φ = r·(p_φ − e·A_φ) is conserved exactly (Kron's and
// Kaluza's cyclic coordinate again), and for a ring centred on the axis
// this means
//
//	p²/B = constant:   p ∝ B^{1/2},   R = p/(eB) ∝ B^{−1/2}
//
// The field comes from a capacitor bank discharged into the solenoid,
// solved in kronspice; the electron is pushed in that field (relativistic
// Boris integrator), and the energy it gains is the work of the induced
// electric field — paid by the bank.

// Compressor is the bank, coil and electron.
type Compressor struct {
	C, V0, L, R float64 // bank [F], charge [V], coil [H], resistance [Ω]
	N           float64 // solenoid turns per metre
	Bias        float64 // steady field before the discharge [T]
	E0          float64 // starting electron kinetic energy [eV]
	Electrons   float64 // electrons in the ring (for the ring current)
}

// DefaultCompressor (declared): 1 mF at 10 kV into 10 µH, 16 turns/m
// (≈ 2 T peak), 0.05 T bias, 2 MeV electrons, 10¹³ of them.
func DefaultCompressor() Compressor {
	return Compressor{C: 1e-3, V0: 1e4, L: 10e-6, R: 5e-3, N: 16, Bias: 0.05, E0: 2e6, Electrons: 1e13}
}

// Compression is the result.
type Compression struct {
	Time, Radius, Gamma, Bfield []float64
	B0, B1, R0, R1, P0, P1      float64 // start and end field, radius, momentum [kg·m/s]
	Invariant                   float64 // largest relative change of r·(p_φ − eA_φ)
	WorkE, DeltaKE              float64 // work of the induced field, kinetic energy gain [J] per electron
	BankEnergy, ToCoil          float64 // bank's starting energy, energy in the coil at peak [J]
	RingCurrent0, RingCurrent1  float64 // N·e·c·β/(2πR) at start and end [A]
	Mismatch                    float64 // kronspice ledger
}

// Run discharges the bank in kronspice and pushes the electron through
// the rise of the coil current (a quarter period).
func (c Compressor) Run() (Compression, error) {
	var out Compression
	ckt, err := circuit.NewBuilder("ring compressor").
		CapacitorIC("CB", "a", "0", c.C, c.V0).
		Resistor("RC", "a", "b", c.R).
		Inductor("LC", "b", "0", c.L).BuildTransient()
	if err != nil {
		return out, err
	}
	quarter := math.Pi / 2 * math.Sqrt(c.L*c.C)
	step := quarter / 4000
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: quarter, MaxStep: step, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return out, fmt.Errorf("kronspice: %w", err)
	}
	il, _ := wf.ElementCurrent("LC")
	vb, _ := wf.NodeVoltage("b") // coil voltage = L·dI/dt
	l := wf.Ledger
	out.Mismatch = l.IntegrationMismatch / l.StoredStart
	out.BankEnergy = 0.5 * c.C * c.V0 * c.V0
	ipk := 0.0
	for _, i := range il {
		ipk = math.Max(ipk, math.Abs(i))
	}
	out.ToCoil = 0.5 * c.L * ipk * ipk
	k := Mu0 * c.N
	last := 0
	at := func(t float64) (B, dB float64) { // field and its rate from the kronspice run
		n := len(wf.Time)
		j := last
		for j > 0 && wf.Time[j] > t {
			j--
		}
		for j < n-2 && wf.Time[j+1] < t {
			j++
		}
		last = j
		f := (t - wf.Time[j]) / (wf.Time[j+1] - wf.Time[j])
		I := il[j] + f*(il[j+1]-il[j])
		v := vb[j] + f*(vb[j+1]-vb[j])
		// the inductor current flows b → 0 and v(b) = L·dI/dt
		return c.Bias + k*math.Abs(I), k * math.Abs(v) / c.L
	}
	// electron starting on the circle centred on the axis
	mc := Me * C0
	gam := 1 + c.E0*Qe/(Me*C0*C0)
	p := mc * math.Sqrt(gam*gam-1)
	B, _ := at(0)
	r := p / (Qe * B)
	out.B0, out.R0, out.P0 = B, r, p
	// electron (charge −e) moving so that it circles the axis: position
	// (r, 0, 0), momentum along −y for B along +z
	x, y := r, 0.0
	px, py := 0.0, p
	q := -Qe
	if q*B > 0 {
		py = -p
	}
	inv0 := x*py - y*px + 0.5*q*B*(x*x+y*y) // r·(p_φ + q·A_φ), A_φ = Br/2
	t := 0.0
	tc := 2 * math.Pi * gam * Me / (Qe * B) // starting cyclotron period
	dt := tc / 400
	out.DeltaKE = 0
	ke0 := (gam - 1) * Me * C0 * C0
	for t < quarter {
		Bm, dBm := at(t + dt/2)
		// fields at the half step (uniform B_z, E_φ = −(r/2)·dB/dt)
		ex, ey := 0.5*dBm*y, -0.5*dBm*x // E = −(dB/2)(−y, x)
		// Boris push
		pxm := px + q*ex*dt/2
		pym := py + q*ey*dt/2
		g := math.Sqrt(1 + (pxm*pxm+pym*pym)/(mc*mc))
		tz := q * Bm * dt / (2 * g * Me)
		sz := 2 * tz / (1 + tz*tz)
		pxp := pxm + pym*tz
		pyp := pym - pxm*tz
		pxm2 := pxm + pyp*sz
		pym2 := pym - pxp*sz
		px = pxm2 + q*ex*dt/2
		py = pym2 + q*ey*dt/2
		g = math.Sqrt(1 + (px*px+py*py)/(mc*mc))
		vx, vy := px/(g*Me), py/(g*Me)
		out.WorkE += q * (ex*vx + ey*vy) * dt
		x += vx * dt
		y += vy * dt
		t += dt
		// adapt the step to the cyclotron period
		dt = 2 * math.Pi * g * Me / (Qe * Bm) / 400
		Bn, _ := at(t)
		inv := x*py - y*px + 0.5*q*Bn*(x*x+y*y)
		out.Invariant = math.Max(out.Invariant, math.Abs(inv-inv0)/math.Abs(inv0))
		if len(out.Time) == 0 || t-out.Time[len(out.Time)-1] > quarter/400 {
			pp := math.Hypot(px, py)
			out.Time = append(out.Time, t)
			out.Radius = append(out.Radius, pp/(Qe*Bn))
			out.Gamma = append(out.Gamma, g)
			out.Bfield = append(out.Bfield, Bn)
		}
	}
	gEnd := math.Sqrt(1 + (px*px+py*py)/(mc*mc))
	out.DeltaKE = (gEnd-1)*Me*C0*C0 - ke0
	out.P1 = math.Hypot(px, py)
	out.B1, _ = at(quarter)
	out.R1 = out.P1 / (Qe * out.B1)
	beta := func(g float64) float64 { return math.Sqrt(1 - 1/(g*g)) }
	out.RingCurrent0 = c.Electrons * Qe * C0 * beta(gam) / (2 * math.Pi * out.R0)
	out.RingCurrent1 = c.Electrons * Qe * C0 * beta(gEnd) / (2 * math.Pi * out.R1)
	return out, nil
}
