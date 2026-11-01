// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Coefficients for these tests are arbitrary (they test the element, not
// a bearing); real ones come from examples/bearing_motor's measurement fit.
const (
	bbR, bbKappa, bbW0 = 2e-3, 4e-6, 20.0 // Ω, N·m/A², rad/s
	bbJ, bbRd          = 1e-5, 5000.0     // kg m², rad/s per N·m (damper)
)

func bearingRig(t *testing.T, amps, w0 float64, thermal bool, eta float64) *circuit.Circuit {
	t.Helper()
	c, err := circuit.NewBuilder("bearing").Domain(exchange.Mechanical, "w").
		CurrentSource("I1", "0", "a", amps).
		BearingMotor("P1", "a", "0", "w", "0", bbR, bbKappa, bbW0, thermal, eta).
		CapacitorIC("J", "w", "0", bbJ, w0).
		Resistor("RD", "w", "0", bbRd).Build()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

// steadyOmega solves κ·I²·ψ(ω) = ω/Rd for the nonzero root by bisection.
func steadyOmega(amps float64, thermal bool, eta float64) float64 {
	psi := func(w float64) float64 {
		s := math.Tanh(w / bbW0)
		if thermal {
			s /= 1 + bbKappa*s*w/(eta*bbR)
		}
		return s
	}
	f := func(w float64) float64 { return bbKappa*amps*amps*psi(w) - w/bbRd }
	lo, hi := 1e-6, 1e6
	for k := 0; k < 200; k++ {
		m := 0.5 * (lo + hi)
		if f(m) > 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi)
}

// From rest the bearing does not start, whatever the current. Given a small
// spin it runs up, in the sense of the spin, to the speed where its
// I²-torque meets the damper — if the current is above the threshold
// I_c = √(ω0/(κ·Rd)) (31.6 A here) at which dT/dω at rest exceeds the
// damping; below it the spin decays (measured: no self-rotation at low
// currents, Moyssides & Hatzikonstantinou 1990). The ledger closes: the source's work equals the
// heat in bearing and damper plus the shaft's kinetic energy.
func TestBearingMotorTran(t *testing.T) {
	for _, mode := range []struct {
		thermal bool
		eta     float64
	}{{false, 0}, {true, 0.1}} {
		for _, amps := range []float64{20, 60, 120} {
			for _, w0 := range []float64{0, 0.5, -0.5} {
				c := bearingRig(t, amps, w0, mode.thermal, mode.eta)
				wf, err := Tran(c, TranSpec{Step: 2e-3, Stop: 2, UIC: true}, TranOptions{})
				if err != nil {
					t.Fatal(err)
				}
				w, _ := wf.NodeVoltage("w")
				end := w[len(w)-1]
				switch {
				case w0 == 0:
					if end != 0 {
						t.Errorf("thermal=%v I=%g: started from rest, ω = %g", mode.thermal, amps, end)
					}
				case amps < math.Sqrt(bbW0/(bbKappa*bbRd)):
					if math.Abs(end) > 1e-3*math.Abs(w0) {
						t.Errorf("thermal=%v I=%g below threshold: ω_end %g, want decay", mode.thermal, amps, end)
					}
				default:
					want := math.Copysign(steadyOmega(amps, mode.thermal, mode.eta), w0)
					if math.Abs(end-want) > 1e-3*math.Abs(want) {
						t.Errorf("thermal=%v I=%g ω0=%g: ω_end %g, want %g", mode.thermal, amps, w0, end, want)
					}
				}
				l := wf.Ledger
				if s := l.EnergyIn + l.EnergyDissipated; s > 0 && math.Abs(l.IntegrationMismatch) > 1e-4*s {
					t.Errorf("thermal=%v I=%g ω0=%g: integration mismatch %g of %g J", mode.thermal, amps, w0, l.IntegrationMismatch, s)
				}
				if l.MaxRelativeResidual > 1e-9 {
					t.Errorf("thermal=%v I=%g ω0=%g: ledger residual %g", mode.thermal, amps, w0, l.MaxRelativeResidual)
				}
				// EM mode: V = R·I + back-EMF; thermal: V = R·I.
				va, _ := wf.NodeVoltage("a")
				v := va[len(va)-1]
				back := 0.0
				if !mode.thermal {
					back = bbKappa * math.Tanh(end/bbW0) * end * amps
				}
				if want := bbR*amps + back; math.Abs(v-want) > 1e-6*want {
					t.Errorf("thermal=%v I=%g: V = %g, want %g", mode.thermal, amps, v, want)
				}
			}
		}
	}
}

// At a running operating point (held by a speed source) the shaft power is
// below η·R·I² in thermal mode and equals the back-EMF power in EM mode.
func TestBearingMotorOPLedger(t *testing.T) {
	for _, thermal := range []bool{false, true} {
		eta := 0.0
		if thermal {
			eta = 0.1
		}
		c, err := circuit.NewBuilder("bearing op").Domain(exchange.Mechanical, "w").
			CurrentSource("I1", "0", "a", 100).
			BearingMotor("P1", "a", "0", "w", "0", bbR, bbKappa, bbW0, thermal, eta).
			VoltageSource("VW", "w", "0", 300).Build()
		if err != nil {
			t.Fatal(err)
		}
		r, err := OP(c, Options{})
		if err != nil {
			t.Fatal(err)
		}
		if rr := r.Ledger.RelativeResidual(); rr > 1e-9 {
			t.Errorf("thermal=%v: ledger residual %g", thermal, rr)
		}
		va := r.NodeVoltages[c.NodeIndex("a")]
		pe := va * 100
		heat := 0.0
		for _, row := range r.Ledger.Rows {
			if row.Element == "P1" && row.Channel == exchange.Thermal {
				heat = -row.PowerW
			}
		}
		shaft := pe - heat
		if thermal {
			if !(shaft > 0 && shaft < eta*bbR*100*100) {
				t.Errorf("thermal: shaft %g W, must be in (0, η·R·I² = %g)", shaft, eta*bbR*1e4)
			}
		} else if math.Abs(heat-bbR*1e4) > 1e-9*heat {
			t.Errorf("EM: heat %g W, want R·I² = %g", heat, bbR*1e4)
		}
	}
}
