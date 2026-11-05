// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func rel(a, b float64) float64 { return math.Abs(a-b) / math.Abs(b) }

// The loop in kronspice resonates at 1/√(LC) with Q = ω0·L/(R_rad + R_d),
// and the AC ledger balances.
func TestERingAC(t *testing.T) {
	e := DefaultERing()
	ac, err := e.Resonance()
	if err != nil {
		t.Fatal(err)
	}
	f0 := e.Omega0() / (2 * math.Pi)
	if d := rel(ac.F0, f0); d > 2e-4 {
		t.Errorf("f0 %.6g Hz, want %.6g (rel %.1e)", ac.F0, f0, d)
	}
	if d := rel(ac.Q, e.Q()); d > 2e-3 {
		t.Errorf("Q %.5g, want %.5g (rel %.1e)", ac.Q, e.Q(), d)
	}
	if ac.MaxResidual > 1e-9 {
		t.Errorf("AC ledger residual %.1e", ac.MaxResidual)
	}
}

// Kicked and left to ring: the decay gives the same Q, every joule of the
// kick is radiated, heated or still stored, no charge leaves the ring, and
// the magnetic moment averages to zero.
func TestRingDown(t *testing.T) {
	e := DefaultERing()
	rd, err := e.RingDown(40)
	if err != nil {
		t.Fatal(err)
	}
	if d := rel(rd.Q, e.Q()); d > 0.01 {
		t.Errorf("ring-down Q %.5g, want %.5g (rel %.1e)", rd.Q, e.Q(), d)
	}
	if bal := math.Abs(rd.EnergyIn-rd.Dissipated-rd.StoredEnd) / rd.EnergyIn; bal > 1e-3 {
		t.Errorf("ledger: in %.4g J, out %.4g + stored %.4g J (rel %.1e)", rd.EnergyIn, rd.Dissipated, rd.StoredEnd, bal)
	}
	if rd.MaxResidual > 1e-6 {
		t.Errorf("tran ledger residual %.1e", rd.MaxResidual)
	}
	if rd.GroundCharge > 1e-9 {
		t.Errorf("net charge lost %.1e of a section's charge", rd.GroundCharge)
	}
	if m := math.Abs(rd.MeanMoment) / rd.PeakMomentLate; m > 0.02 {
		t.Errorf("mean moment %.3g of peak, want ≈ 0", m)
	}
}

// The permittivity is the confinement: Q_rad ∝ εr^(3/2), and kronspice
// agrees with the formula across the sweep.
func TestQScaling(t *testing.T) {
	e := DefaultERing()
	q := func(er float64) float64 { x := e; x.EpsR = er; return x.QRad() }
	if d := rel(q(1000)/q(100), math.Pow(10, 1.5)); d > 1e-12 {
		t.Errorf("Q_rad(1000)/Q_rad(100) off ε^1.5 by %.1e", d)
	}
	if d := rel(e.QRad(), e.Omega0()*e.L()/e.RRad()); d > 1e-12 {
		t.Errorf("closed form 6Λ/(π(ka)³) differs from ω0L/R_rad by %.1e", d)
	}
	x := e
	x.EpsR = 300
	ac, err := x.Resonance()
	if err != nil {
		t.Fatal(err)
	}
	if d := rel(ac.Q, x.Q()); d > 2e-3 {
		t.Errorf("εr 300: AC Q %.5g, want %.5g", ac.Q, x.Q())
	}
}

// The bench's required ring, built as an ERing, resonates at mc²/ħ and has
// exactly the Q the SN1987A proper time demands; the vacuum ring is
// electrically large and lives under one period.
func TestBench(t *testing.T) {
	b := DefaultBench()
	r := b.Run()
	if d := rel(r.ProperTime, b.TravelS*b.MassEV/b.EnergyEV); d > 1e-12 {
		t.Errorf("proper time off by %.1e", d)
	}
	e := ERing{A: r.RadiusNeeded, R: r.RadiusNeeded / b.AoverR, EpsR: r.EpsRNeeded, N: 1}
	if d := rel(e.Omega0(), r.Omega0); d > 1e-9 {
		t.Errorf("required ring resonates at %.6g rad/s, want %.6g", e.Omega0(), r.Omega0)
	}
	if d := rel(e.QRad(), r.QNeeded); d > 1e-9 {
		t.Errorf("required ring Q %.6g, want %.6g", e.QRad(), r.QNeeded)
	}
	if r.VacuumKA < 1 || r.VacuumQ > 1 {
		t.Errorf("vacuum ring ka %.3g, Q %.3g: expected electrically large and Q < 1", r.VacuumKA, r.VacuumQ)
	}
}

// The inverted ring line: AC peaks follow the CL ladder's dispersion, whole
// wavelengths on the plain ring and half-integer ones on the crossed ring.
func TestDualRingModes(t *testing.T) {
	if testing.Short() {
		t.Skip("AC sweeps of two 100-section rings")
	}
	for _, crossed := range []bool{false, true} {
		d := DefaultDualRing(crossed)
		peaks, resid, err := d.Modes(4, 10e6, 200e6)
		if err != nil {
			t.Fatal(err)
		}
		if len(peaks) < 4 {
			t.Fatalf("crossed=%v: found %d peaks", crossed, len(peaks))
		}
		for m := 1; m <= 4; m++ {
			if dd := rel(peaks[m-1], d.ModeFrequency(m)); dd > 1e-4 {
				t.Errorf("crossed=%v mode %d: %.6g Hz, theory %.6g (rel %.1e)", crossed, m, peaks[m-1], d.ModeFrequency(m), dd)
			}
		}
		if resid > 1e-9 {
			t.Errorf("crossed=%v: AC ledger residual %.1e", crossed, resid)
		}
	}
}

// The zero-point field as inflow: the ring settles at the expected energy
// (½ħω0 for its noise modes), the magnetic channel's books close (zero-point
// in = radiated out + stored), the energy fluctuates like random noise, and
// the magnetic moment still averages to zero.
func TestZeroPointInflow(t *testing.T) {
	if testing.Short() {
		t.Skip("zero-point noise runs")
	}
	z := DefaultZPFRing()
	z.Runs = 3
	r, err := z.Run()
	if err != nil {
		t.Fatal(err)
	}
	if d := math.Abs(r.Energy - r.ExpEnergy); d > math.Max(4*r.Spread, 0.03) {
		t.Errorf("stored %.4f ± %.4f ħ, expected %.4f", r.Energy, r.Spread, r.ExpEnergy)
	}
	if bal := math.Abs(r.ZPIn-r.RadOut-r.Stored) / r.RadOut; bal > 2e-3 {
		t.Errorf("magnetic channel: in %.4g = out %.4g + stored %.4g? (rel %.1e)", r.ZPIn, r.RadOut, r.Stored, bal)
	}
	if r.Fluctuation < 0.7 || r.Fluctuation > 1.3 {
		t.Errorf("energy fluctuation σ/⟨E⟩ = %.2f, expected ≈ 1 for a random field", r.Fluctuation)
	}
	if r.MeanMomentRatio > 0.01 {
		t.Errorf("mean moment %.2g of its rms, expected ≈ 0", r.MeanMomentRatio)
	}
}
