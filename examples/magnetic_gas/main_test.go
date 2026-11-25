// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

func near(t *testing.T, what string, got, want, rel float64) {
	t.Helper()
	if math.Abs(got-want) > rel*math.Abs(want) {
		t.Errorf("%s: %.6g, want %.6g (±%.2g)", what, got, want, rel)
	}
}

// The Paschen rule passes through each gas's measured minimum and is
// lowest there; the air check and the 500 psi argon value are the closed
// form V_min·x/(1 + ln x).
func TestPaschen(t *testing.T) {
	for _, g := range append(Gases, Air) {
		d := g.PDmin * torr / atm / 100 // gap [m] that puts 1 atm at the minimum
		near(t, g.Formula+" at its minimum", g.Breakdown(atm, d), g.Vmin, 1e-12)
		if g.Breakdown(atm, 0.9*d) <= g.Vmin || g.Breakdown(atm, 1.1*d) <= g.Vmin {
			t.Errorf("%s: not a minimum", g.Formula)
		}
	}
	x := 760.0 / 0.567
	near(t, "air 1 cm, 1 atm", Air.Breakdown(atm, 0.01), 327*x/(1+math.Log(x)), 1e-12)
	x = 500 * psi / torr * 0.02 / 0.9
	near(t, "argon 0.2 mm, 500 psi", Argon.Breakdown(500*psi, 0.2e-3), 137*x/(1+math.Log(x)), 1e-12)
	if !math.IsInf(Argon.Breakdown(torr, 0.1*0.9/math.E/100), 1) {
		t.Error("no breakdown expected left of x = 1/e")
	}
	// Raizer's form for argon, 0.2 mm at 500 psi, by hand: 12.9 kV.
	near(t, "Raizer argon", Argon.BreakdownRaizer(500*psi, 0.2e-3), 12.9e3, 0.01)
}

// Neutral susceptibilities: O₂ at 1 atm is the handbook's 1.9e-6 (CGS
// molar 3449e-6 cm³/mol), and its Curie law with S = 1 (μ_eff = 2√2)
// gives the same within 2 %; argon and hydrogen are diamagnetic; HHO is
// paramagnetic through its oxygen.
func TestMagnetism(t *testing.T) {
	near(t, "O₂ χ at 1 atm", Oxygen.ChiNeutral(atm, Troom), 1.80e-6, 0.01)
	n := Density(atm, Troom)
	near(t, "O₂ Curie S = 1", ChiCurie(n, 2*math.Sqrt2, Troom), Oxygen.ChiNeutral(atm, Troom), 0.02)
	if Argon.ChiNeutral(atm, Troom) >= 0 || Hydrogen.ChiNeutral(atm, Troom) >= 0 || HHO.ChiNeutral(atm, Troom) <= 0 {
		t.Error("signs: Ar, H₂ diamagnetic; HHO paramagnetic")
	}
	// Removing an electron turns argon paramagnetic and makes O₂ less so.
	if z := Argon.Ionized(500*psi, 1, Troom); z.Total <= 0 {
		t.Error("ionized argon should be paramagnetic")
	}
	if z := Oxygen.Ionized(500*psi, 1, Troom); z.Total >= Oxygen.ChiNeutral(500*psi, Troom) {
		t.Error("O₂⁺ should be less paramagnetic than O₂")
	}
	// Electrons: spin 1, Landau −1/3: net 2/3 of the spin Curie value.
	near(t, "free electrons", ChiElectrons(n, Troom), 2.0/3*ChiCurie(n, math.Sqrt(3), Troom)*3/3, 1e-12)
}

// Compression numbers checked by hand.
func TestCompression(t *testing.T) {
	near(t, "1 T magnetic pressure", MagneticPressure(1), 397887, 1e-5)
	near(t, "field to hold 500 psi", FieldToHold(500*psi), 2.9434, 1e-4)
	near(t, "argon adiabatic", Argon.AdiabaticT(Troom, atm, 500*psi), 1201.7, 1e-3)
	near(t, "HHO chemical, 1 cm³", HHOChemical(500*psi, 1e-6), 228.1, 2e-3)
	near(t, "ionize 1 cm³ argon", Argon.IonizeEnergy(500*psi, 1e-6), 2148, 2e-3)
}

// The Marx generator charges every stage to Vc, erects to 3.06·Vc into
// 50 Ω through 1 µH (3.85·Vc would need no inductance), and the ledger
// closes in both phases: firing heat = stored energy.
func TestMarx(t *testing.T) {
	for _, g := range Gases {
		r, err := DefaultMarx(g).Run()
		if err != nil {
			t.Fatal(err)
		}
		for k, v := range r.StageV {
			near(t, g.Formula+" stage voltage "+string(rune('1'+k)), v, r.Vc, 1e-3)
		}
		near(t, g.Formula+" gain", r.Gain, 3.06, 0.01)
		near(t, g.Formula+" charging heat = stored", r.ChargeHeat, r.Stored, 2e-3)
		near(t, g.Formula+" firing heat = stored", r.EHeat, r.Stored, 1e-3)
		if math.Abs(r.Mismatch) > 1e-3 || math.Abs(r.ChargeMismatch) > 1e-3 {
			t.Errorf("%s: books %.2g, %.2g", g.Formula, r.ChargeMismatch, r.Mismatch)
		}
	}
}

// ECE2 (H48): the torsion term moves the shot by < 1e-5 here; the scalar
// term scales every terminal voltage by 1/(1 − 2κ) and delivers no energy.
func TestMarxECE2(t *testing.T) {
	var rs []MarxResult
	for _, m := range MarxECE(Argon) {
		r, err := m.Run()
		if err != nil {
			t.Fatal(err)
		}
		rs = append(rs, r)
	}
	near(t, "torsion peak", rs[1].Peak, rs[0].Peak, 1e-5)
	near(t, "scalar peak", rs[2].Peak, rs[0].Peak/(1-2*0.05), 1e-4)
	near(t, "scalar heat", rs[2].EHeat, rs[0].EHeat, 1e-4)
	if math.Abs(rs[2].EModel) > 1e-9 {
		t.Errorf("scalar model-supplied %.3g J", rs[2].EModel)
	}
}

// The pump reaches close to r·Vin; the 100 V supply pays the load's current
// at 100 V and the plasma supply pays the rest; books close.
func TestPump(t *testing.T) {
	r, err := DefaultPump(Argon).Run()
	if err != nil {
		t.Fatal(err)
	}
	if r.Vout > r.Ideal || r.Vout < 0.95*r.Ideal {
		t.Errorf("pump output %.4g V, ideal %.4g V", r.Vout, r.Ideal)
	}
	near(t, "supply power = I_load·Vin", r.PIn, r.Vout/r.Rload*r.Vin, 0.03)
	if r.PPlasma < 0.5*r.PLoad {
		t.Errorf("plasma supply pays %.3g W of %.3g W", r.PPlasma, r.PLoad)
	}
	if math.Abs(r.Mismatch) > 5e-3 {
		t.Errorf("books %.2g", r.Mismatch)
	}
	if r.NeutralRatio > 1.02 {
		t.Errorf("neutral compression ratio %.4f", r.NeutralRatio)
	}
}
