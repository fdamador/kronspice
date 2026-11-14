// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
)

func part11() {
	fmt.Println()
	fmt.Println("Part 11a. g of the photon ring: g = M c² / E_c (rest energy over the part that circulates)")
	mu, g := WvdMMoment()
	fmt.Printf("  H23 as written (loop λ_C/4π, all energy circulating): μ = %.3f μ_B, g = %.3f   (measured 2.00232)\n", mu, g)
	fmt.Println("  light held by a squeeze k·R^n: equilibrium puts 1/(1+1/n) of the mass in the light, so g = 1 + 1/n")
	fmt.Println("  squeeze                      n      light share   g")
	for _, s := range []struct {
		name string
		n    float64
	}{{"ether pressure (Part 9)", 3}, {"ether tension (Part 9)", 2}, {"constant-force string", 1}, {"softer, R^½", 0.5}} {
		gg := SplitG(s.n, 0)
		fmt.Printf("  %-27s  %-5.2g  %-12.4f  %.4f\n", s.name, s.n, 1/gg, gg)
	}
	R, k, ratio := StringRing(Me)
	fmt.Println("  Only a constant-force string (n = 1) gives g = 2. Then the light carries exactly half the rest energy,")
	fmt.Printf("  its spin (mc/2)·R = ħ/2 puts it on R = ħ/mc = %.4g m (twice H23's loop), the string pulls with\n", R)
	fmt.Printf("  %.4g N (%.3g keV/fm; QCD's string is ~0.9 GeV/fm), and the half-energy photon's wavelength is\n", k, k*1e-15/1.602176634e-16)
	fmt.Printf("  %.6f circumferences: it closes on itself after two turns, the double loop. (That part holds for any\n", ratio)
	fmt.Println("  photon carrying ħ/2 on a circle, λ = h/p = 4πR; what n = 1 adds is the energy split that makes g = 2.)")

	fmt.Println()
	fmt.Println("Part 11b. kronspice: the crossed ring holds the half-energy photon")
	plain := DefaultRing()
	plain.Crossed = false
	mp, err := plain.Measure(1, 0.02, 4001)
	if err != nil {
		log.Fatal(err)
	}
	mc, err := DefaultRing().Measure(0, 0.02, 4001)
	if err != nil {
		log.Fatal(err)
	}
	fp, _ := plain.Mode(1)
	fc, _ := DefaultRing().Mode(0)
	fmt.Printf("  plain ring, lowest wave (one wavelength per turn)     %.6g Hz (ladder %.6g)\n", mp.F, fp)
	fmt.Printf("  crossed ring, lowest wave (one wavelength per 2 turns) %.6g Hz (ladder %.6g)\n", mc.F, fc)
	fmt.Printf("  ratio %.6f: on the same ring the double loop's photon has half the energy, as 11a needs.\n", mc.F/mp.F)

	fmt.Println()
	fmt.Println("Part 11c. The anomaly: the charge's static field does not go round")
	fmt.Println("  self-energy U = α·ħc·β/R sits outside the light, so g/2 = 1 + 2αβ; Schwinger's α/2π needs β = 1/4π = 0.0796")
	fmt.Println("  charge spread over                    β        (g − 2)/2    / measured")
	for _, a := range []float64{0.01, 0.1, 0.5, 1} {
		b := Beta(TorusInvD(a, 800))
		fmt.Printf("  the light's torus, tube a = %-4g R     %-7.4f  %-11.5f  %.2f\n", a, b, Anomaly(b), Anomaly(b)/AeMeasured)
	}
	for _, s := range []float64{1, 2} {
		b := 1 / (2 * s)
		fmt.Printf("  a shell of radius %g R                  %-7.4f  %-11.5f  %.2f\n", s, b, Anomaly(b), Anomaly(b)/AeMeasured)
	}
	se, sm := ShellForSchwinger(false), ShellForSchwinger(true)
	fmt.Printf("  a shell of radius 2π R = λ_C           %-7.4f  %-11.5f  %.4f\n", BetaSchwinger, Anomaly(BetaSchwinger), Anomaly(BetaSchwinger)/AeMeasured)
	fmt.Printf("  Charge on the light's own torus gives 3.7–13× too much. The anomaly lands on α/2π if the charge's\n")
	fmt.Printf("  field is spread over a shell of radius %.4f R = λ_C (the ring's circumference), or %.3f R if the\n", se, sm)
	fmt.Println("  moving charge's magnetic self-energy counts as non-circulating too. That radius is one fitted number.")
	fmt.Printf("  α/2π = %.8f against measured %.8f: 0.15%% high, the size of QED's next term, −0.328(α/π)².\n", AeSchwinger, AeMeasured)
	fmt.Printf("  Alternatively, no self-energy and a string slightly softer than constant force, n = %.6f.\n", NForAnomaly())
	fmt.Println("  Standard prediction: g − 2 from QED loops, to 12 digits, with a pointlike electron.")
}
