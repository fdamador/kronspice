// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// Three coils at 120° with an ordinary diode loss: the array field holds
// only 3n ± 1 (3n+1 forward, 3n−1 backward), multiples of 3 cancel to the
// numerical floor, while one coil carries every harmonic. The ledger closes.
func TestThreeCoilSignature(t *testing.T) {
	s, err := DefaultCoils.Run(9)
	if err != nil {
		t.Fatal(err)
	}
	const floor = 1e-7
	for n := 1; n <= 9; n++ {
		switch n % 3 {
		case 0:
			if s.Fwd[n] > floor || s.Bwd[n] > floor {
				t.Errorf("n=%d (multiple of 3): fwd %.2g bwd %.2g, want < %.0g", n, s.Fwd[n], s.Bwd[n], floor)
			}
			if s.Single[n] < 1e-5 {
				t.Errorf("n=%d: one coil %.2g, want a visible harmonic", n, s.Single[n])
			}
		case 1:
			if s.Fwd[n] < 1e-5 || s.Bwd[n] > floor {
				t.Errorf("n=%d (3k+1): fwd %.2g bwd %.2g, want forward only", n, s.Fwd[n], s.Bwd[n])
			}
		case 2:
			if s.Bwd[n] < 1e-5 || s.Fwd[n] > floor {
				t.Errorf("n=%d (3k−1): fwd %.2g bwd %.2g, want backward only", n, s.Fwd[n], s.Bwd[n])
			}
		}
		// The rotating parts of a real three-phase set equal the one-coil
		// harmonics when every coil is identical.
		if n%3 != 0 && n > 1 {
			got := math.Max(s.Fwd[n], s.Bwd[n])
			if math.Abs(got-s.Single[n]) > 1e-3*s.Single[n] {
				t.Errorf("n=%d: array %.6g vs one coil %.6g", n, got, s.Single[n])
			}
		}
	}
	bal := s.EnergyIn - s.EnergyOut - s.Dissipated - (s.StoredEnd - s.StoredStart)
	if math.Abs(bal) > 1e-6*s.EnergyIn || s.Residual > 1e-9 {
		t.Errorf("ledger: imbalance %.3g J of %.3g J in, residual %.2g", bal, s.EnergyIn, s.Residual)
	}
}

// Without the nonlinearity a coil carries no harmonics at all.
func TestLinearControl(t *testing.T) {
	k := DefaultCoils
	k.Diode = false
	k.N = 1
	s, err := k.Run(6)
	if err != nil {
		t.Fatal(err)
	}
	for n := 2; n <= 6; n++ {
		if s.Single[n] > 1e-9 {
			t.Errorf("n=%d: %.2g in a linear coil", n, s.Single[n])
		}
	}
}

func TestLoop(t *testing.T) {
	l := Loop{Length: 10, Ng: 1.46, T: 1000}
	if f := l.FSR(); math.Abs(f-20.533e6) > 1e3 {
		t.Errorf("FSR %.6g Hz", f)
	}
	if w := l.WienPeak(); math.Abs(w-2.8978e-6) > 1e-9 {
		t.Errorf("Wien peak %.6g m", w)
	}
}

// The budget's entries against hand values.
func TestBudget(t *testing.T) {
	b := DefaultRig.Budget()
	checks := []struct {
		name      string
		got, want float64
	}{
		{"claim low [kg]", b.ClaimLow, 5e-4},
		{"claim high [kg]", b.ClaimHigh, 0.05},
		{"E/c² [kg]", b.HeatMass, 50 * 1000 * (1000 - 293.15) / (c0 * c0)},
		{"P/c [kg]", b.RecoilMax, 1e5 / c0 / gn},
		{"lead [kg/m]", b.LeadPerMetre, 250 * 50e-6 / gn},
	}
	for _, c := range checks {
		if math.Abs(c.got-c.want) > 1e-12*math.Abs(c.want) {
			t.Errorf("%s: %.6g, want %.6g", c.name, c.got, c.want)
		}
	}
	// The ordinary effects reach the claimed band: in air, buoyancy alone
	// exceeds the low end; recoil stays below it.
	if b.AirBuoyancy < b.ClaimLow || b.RecoilMax > b.ClaimLow {
		t.Errorf("budget ordering changed: buoyancy %.3g, recoil %.3g, claim %.3g kg", b.AirBuoyancy, b.RecoilMax, b.ClaimLow)
	}
	// Steady state: the coil radiates less than the heater supplies.
	if b.RadiatedPower > DefaultRig.Power {
		t.Errorf("radiated %.3g W above heater %.3g W", b.RadiatedPower, DefaultRig.Power)
	}
}

func TestPaperArithmetic(t *testing.T) {
	if a := Accel(1e3, 1e6); a != 1e-3 {
		t.Errorf("1 kN on 1000 t: %g m/s²", a)
	}
	if f := ResidualMassFraction(1e3, 1e6, 6000*gn); math.Abs(f-1.6996e-8) > 1e-11 {
		t.Errorf("residual mass fraction %.5g", f)
	}
}

func TestMain_runs(t *testing.T) {
	if err := run(); err != nil {
		t.Fatal(err)
	}
}

// Peck & Huang's H₂ value at the sodium D line, and the TIR ceiling.
func TestHydrogenTIR(t *testing.T) {
	if d := H2Refractivity(589.3e-9); math.Abs(d-1.388e-4) > 0.005e-4 {
		t.Errorf("H2 n−1 at 589 nm: %.5g", d)
	}
	got := TIRMaxDeflection(1e12, 589e-9) // 10⁶ cm⁻³
	want := 2 * math.Sqrt(2*H2Refractivity(589e-9)*1e12/loschmidt)
	if math.Abs(got-want) > 1e-12*want || got*radToArcsec > 2e-3 {
		t.Errorf("deflection %.4g rad = %.3g arcsec", got, got*radToArcsec)
	}
	// Chromatic: blue turned more than near-IR.
	if r := TIRMaxDeflection(1e12, 400e-9) / TIRMaxDeflection(1e12, 1e-6); r < 1.01 || r > 1.05 {
		t.Errorf("400 nm / 1 µm ratio %.4f", r)
	}
}
