// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command space_power_generator runs P. Tewari's published test of his
// Space Power Generator (SPG, hypothesis H13) through the energy ledger.
// Source: P. Tewari, Universal Principles of Space and Matter, Grain and
// Chaff, Toronto, 2017 (1st ed. 2009), §6.7, pp. 149–161
// (reference/tewari_upsm_2017.pdf). The SPG is a homopolar machine — an
// electromagnet mounted on a rotating iron core, a static yoke, output
// brushes on slip rings (Fig. 6-10, p. 155) — whose EMF the book gives as
// V = ½ωBR² (eq. 6.16, p. 154), the same as Tesla's and Faraday's disk.
//
// The test (pp. 159–160): a dc drive motor (DM) at 2700 r/min;
//
//	no load:  DM input W0 = 2943 W (windage, friction), excitation 650 W,
//	          EMF E = 3.41 V;
//	loaded:   I = 3438 A, DM input W_L = 8620 W, speed unchanged;
//	          W_r = W_L − W0 = 5677 W                        (6.20)
//	          W_E = E·I = 3.41 V × 3438 A = 11 724 W        (6.21)
//	          η_incremental = W_E / W_r = 2.066             (6.22)
//	          η_SPG = 2.066 / 0.83 = 2.49 (DM efficiency 0.83)  (6.23)
//
// The ledger checks three things:
//
//  1. Arithmetic, as the book states it.
//
//  2. Accounting: all inputs, including the 650 W excitation, against the
//     claimed output.
//
//  3. The assumption behind W_E. The book multiplies the load current by
//     the EMF measured at no load. Under load, standard electrodynamics
//     requires the mechanical power the rotor takes from the shaft to equal
//     E_load·I (Poynting's theorem; the torque on a homopolar rotor is
//     k·I when E = k·ω, wherever the reaction lands — on the static yoke
//     and the static conductor A–A₁). With the shaft power the book's own
//     numbers allow, that fixes the largest loaded EMF standard physics
//     permits. The book reports no loaded-EMF measurement; it says the
//     drops "can also be" measured (p. 160). That measurement — the sum of
//     the voltage drops around the loaded circuit, or a torque transducer
//     on the shaft — decides between the two.
//
//     go run ./examples/space_power_generator
package main

import (
	"fmt"
	"math"
)

// Test is the published test, p. 159–160.
type Test struct {
	RPM        float64 // r/min
	W0         float64 // DM input at no load [W]
	Excitation float64 // SPG excitation power [W]
	E0         float64 // no-load EMF [V]
	I          float64 // load current [A]
	WL         float64 // DM input loaded [W]
	EtaDM      float64 // DM efficiency used by the book
}

// Published is the test as printed.
func Published() Test {
	return Test{RPM: 2700, W0: 2943, Excitation: 650, E0: 3.41, I: 3438, WL: 8620, EtaDM: 0.83}
}

// Ledger is the analysis.
type Ledger struct {
	Wr, WE                  float64 // eqs. 6.20, 6.21
	EtaInc, EtaSPG          float64 // eqs. 6.22, 6.23
	TotalIn                 float64 // DM loaded input + excitation
	TotalRatio              float64 // W_E / total input
	ShaftIncrease           float64 // η_DM · W_r: extra shaft power into the SPG
	ELoadMax, ELoadMaxIdeal float64 // largest loaded EMF allowed (η_DM as stated; η_DM = 1) [V]
	FluxDropMin             float64 // smallest loss of EMF under load the ledger requires (fraction)
	TorqueClaimed           float64 // k·I with k = E0/ω: torque the book's W_E implies [N m]
	TorqueMeasured          float64 // shaft torque increase from η_DM·W_r [N m]
	RCircuit                float64 // E0/I: total loop resistance implied [Ω]
	Omega                   float64 // [rad/s]
}

// Analyse builds the ledger.
func (t Test) Analyse() Ledger {
	var l Ledger
	l.Omega = t.RPM * 2 * math.Pi / 60
	l.Wr = t.WL - t.W0
	l.WE = t.E0 * t.I
	l.EtaInc = l.WE / l.Wr
	l.EtaSPG = l.EtaInc / t.EtaDM
	l.TotalIn = t.WL + t.Excitation
	l.TotalRatio = l.WE / l.TotalIn
	l.ShaftIncrease = t.EtaDM * l.Wr
	// Standard: E_load·I = shaft power into the electrical conversion ≤
	// the shaft-power increase (any extra mechanical loss under load only
	// lowers it).
	l.ELoadMax = l.ShaftIncrease / t.I
	l.ELoadMaxIdeal = l.Wr / t.I
	l.FluxDropMin = 1 - l.ELoadMax/t.E0
	l.TorqueClaimed = t.E0 * t.I / l.Omega
	l.TorqueMeasured = l.ShaftIncrease / l.Omega
	l.RCircuit = t.E0 / t.I
	return l
}

func main() {
	t := Published()
	l := t.Analyse()
	fmt.Println("Tewari, Space Power Generator test — Universal Principles of Space and Matter (2017), pp. 159–160")
	fmt.Println()
	fmt.Println("1. The book's arithmetic")
	fmt.Printf("   W_r = %.0f − %.0f = %.0f W (6.20); W_E = %.2f V × %.0f A = %.0f W (6.21)\n", t.WL, t.W0, l.Wr, t.E0, t.I, l.WE)
	fmt.Printf("   incremental efficiency %.3f (6.22, printed 2.066 from 11727 W); ÷ η_DM %.2f = %.2f (6.23)\n", l.EtaInc, t.EtaDM, l.EtaSPG)
	fmt.Println()
	fmt.Println("2. All inputs")
	fmt.Printf("   DM %.0f W + excitation %.0f W = %.0f W in; claimed electrical %.0f W → ratio %.2f\n", t.WL, t.Excitation, l.TotalIn, l.WE, l.TotalRatio)
	fmt.Println()
	fmt.Println("3. The assumption: E under load = E at no load")
	fmt.Printf("   ω = %.1f rad/s; implied loop resistance E0/I = %.3g mΩ (brushes, rings, SPG, load)\n", l.Omega, l.RCircuit*1e3)
	fmt.Printf("   extra shaft power into the SPG under load: %.2f × %.0f W = %.0f W → %.1f N·m\n", t.EtaDM, l.Wr, l.ShaftIncrease, l.TorqueMeasured)
	fmt.Printf("   standard electrodynamics: E_load·I = shaft power (less extra mechanical loss) → E_load ≤ %.2f V (≤ %.2f V even if the DM were lossless)\n", l.ELoadMax, l.ELoadMaxIdeal)
	fmt.Printf("   i.e. the EMF must fall by at least %.0f %% under 3438 A (armature reaction, saturation, speed sag at the rotor)\n", 100*l.FluxDropMin)
	fmt.Printf("   the book's W_E needs the EMF unchanged and a back torque of only %.1f N·m against the %.1f N·m that k·I gives\n", l.TorqueMeasured, l.TorqueClaimed)
	fmt.Println()
	fmt.Println("   Decisive measurement (not reported in the book): the sum of the voltage drops around the loaded loop")
	fmt.Println("   (brushes, rings, conductor, shunt, load) at 3438 A, or the shaft torque with a torque transducer.")
	fmt.Printf("   Standard: loaded EMF ≤ %.2f V, consistent with the %.1f N·m the drive motor delivered. Claim: %.2f V, which would need %.1f N·m.\n", l.ELoadMax, l.TorqueMeasured, t.E0, l.TorqueClaimed)
}
