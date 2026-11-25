// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// TestShortedEnds: the coil's two ends joined through 1 mΩ (Coil.Short).
// It is then a shorted transformer secondary: below its own loop
// resonances the current round it is (M/L)·I_primary (the flux linkage
// M·I_p + L·I = 0 of a lossless closed loop), checked at 10 kHz to 1 % for
// the ordinary and the Tesla bifilar solenoid. Pulsed at the bifilar's
// open-ends f1, the joined coil carries less current (wire peak below a
// fifth of the open coil's) and stores less than a thousandth of its
// energy: joining the ends removes the resonance instead of trading volts
// for amps. Ledgers close to 1 %.
func TestShortedEnds(t *testing.T) {
	for _, c := range []*Coil{NewCoil(Solenoid(), false), NewCoil(Solenoid(), true)} {
		c.Short = 1e-3
		ckt, err := c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, c.SectionResistance(1e4))
		if err != nil {
			t.Fatal(err)
		}
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: 1e4, FStop: 1e4}, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		is, _ := res.ElementCurrent("RSHORT")
		ip, _ := res.ElementCurrent("RS")
		m, _ := c.PrimaryM()
		if got, want := cabs(is[0])/cabs(ip[0]), m/c.TotalL(); math.Abs(got-want)/want > 0.01 {
			t.Errorf("%s: I_loop/I_primary %.5g, M/L %.5g", c.Name, got, want)
		}
	}
	const f1 = 41.61e3 // the Tesla bifilar solenoid's open-ends resonance
	open, shorted := NewCoil(Solenoid(), true), NewCoil(Solenoid(), true)
	shorted.Short = 1e-3
	var runs [2]PulseRun
	for k, c := range []*Coil{open, shorted} {
		r, err := c.Pulse(f1, f1, 0.1/f1, 60, false)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r.EIn-r.EDiss-r.EStore) > 0.01*r.EIn {
			t.Errorf("%v ledger: in %.4g, heat %.4g, stored %.4g", k, r.EIn, r.EDiss, r.EStore)
		}
		runs[k] = r
	}
	if runs[1].PeakWire > runs[0].PeakWire/5 {
		t.Errorf("wire current joined %.4g A, open %.4g A", runs[1].PeakWire, runs[0].PeakWire)
	}
	if runs[1].EStore > runs[0].EStore/1000 {
		t.Errorf("stored joined %.4g J, open %.4g J", runs[1].EStore, runs[0].EStore)
	}
}
