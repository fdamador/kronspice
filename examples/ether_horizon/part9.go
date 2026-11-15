// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
	"math"
	"path/filepath"
)

// BagB is the MIT bag constant, B^{1/4} ≈ 145 MeV ⇒ B ≈ 57 MeV/fm³ [Pa]
// (DeGrand, Jaffe, Johnson & Kiskis, Phys. Rev. D 12, 2060, 1975).
const BagB = 57.5 * 1.602176634e-13 / 1e-45

func part9(vizDir string) {
	fmt.Println()
	fmt.Println("Part 9a. A round void: the light inside pushes out, the nonlinear ether pushes in")
	fmt.Println("  E(R) = ħc·x/R + (4π/3)·p∞·R³ + 4π·σ·R²  (one photon in a cavity mode of root x; the MIT bag's structure)")
	fmt.Println("  mode     held by        radius R [m]   R / ring radius   R / ƛ_C   needs                       stiffness d²E/dR²")
	for _, m := range []struct {
		name string
		x    float64
	}{{"TM₁₁", XTM11}, {"TE₁₁", XTE11}} {
		R, p := PressureForElectron(m.x)
		Rs, sg := TensionForElectron(m.x)
		lc := Hbar / (Me * C0)
		fmt.Printf("  %-7s  pressure p∞     %-13.4g  %-16.2f  %-8.3f  p∞ = %.3g Pa (%.0e B)   %.3g N/m\n", m.name, R, R/RingRadius(), R/lc, p, p/BagB, VoidStiffness(R, m.x, p, 0))
		fmt.Printf("  %-7s  tension σ       %-13.4g  %-16.2f  %-8.3f  σ = %.3g N/m               %.3g N/m\n", m.name, Rs, Rs/RingRadius(), Rs/lc, sg, VoidStiffness(Rs, m.x, 0, sg))
	}
	fmt.Println("  Both are stable (d²E/dR² > 0): squeeze the void and the light pushes back harder (1/R²), let it grow")
	fmt.Println("  and the ether wins. Sized for m_e c², the void is 7–13 ring radii across: the photon ring fits inside it.")
	fmt.Println("  Standard prediction: the electron shows no size down to < 10⁻¹⁸ m, a million times smaller.")

	fmt.Println()
	fmt.Println("Part 9b. kronspice twin: the void's light as an LC tank whose L and C grow with the void (s = R/R0)")
	bt := BreathingTank{L0: 1e-6, C0: 1e-9, I0: 1, SEnd: 2, Periods: 50, StepsPerPeriod: 800}
	r, err := bt.Run()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  void doubles (τ = 50 periods): frequency %.6g → %.6g Hz (ratio %.6f), energy %.6g → %.6g J\n",
		r.Freq0, r.Freq1, r.Freq1/r.Freq0, r.W0, r.W1)
	fmt.Printf("  W·R constant to %.1e (the photon's ħω ∝ 1/R); work on the moving wall %.6g J = the energy lost to %.1e;\n",
		r.Invariant-1, r.WallWork, r.BooksResidual)
	fmt.Printf("  kronspice's own integration mismatch %.1e J. That wall work is the light's radiation pressure on the void.\n", r.EngineMismatch)

	fmt.Println()
	fmt.Println("Part 9c. The fountain round the void: does it stay round?")
	fmt.Println("  At the centre of the two-pole flow the ether comes in round the equator and leaves along the axis: a")
	fmt.Println("  stagnation point, where the ether's pressure is highest. A round void there makes the flow slide over it,")
	fmt.Println("  |u| = (5/2)·E·b·|sin 2θ|: still at the poles and the equator, fastest at 45°. Linearized Young–Laplace:")
	fmt.Println("  Weber We = ρ(Eb)²b/σ   ε₂        ε₄         pole r/b   45° r/b   equator r/b   pole/equator")
	for _, we := range []float64{0, 0.1, 0.3, 1} {
		e := ShapeModes(we)
		fmt.Printf("  %-22.2f %-9.4f %-10.4f %-10.4f %-9.4f %-13.4f %.3f\n", we, e[2], e[4], WallRadius(0, we), WallRadius(math.Pi/4, we), WallRadius(math.Pi/2, we), Elongation(we))
	}
	fmt.Printf("  The flow stretches the void along the axis, pinching its waist: 10 %% longer than wide at We = %.3f,\n", WeberFor(1.1))
	fmt.Println("  and on towards the line void of Parts 1 and 8 as the flow wins (beyond the linear theory).")
	fmt.Printf("  So the void stays round only while its surface tension beats the flow: ½ρ(Eb)² ≲ %.2f·σ/b.\n", WeberFor(1.1)/2)
	fmt.Println("  Held by pressure alone (σ = 0) it has nothing to restore its shape in this model: any flow stretches it.")
	fmt.Println("  The axis vortex adds hollow tails out to the poles, of radius C0/√(2p∞/ρ) (Rankine's hollow core):")
	fmt.Println("  thin while the swirl's ½ρU_φ² at the wall stays below p∞ — a round bead on a thread.")
	if vizDir != "" {
		var rows [][]float64
		for i := 0; i <= 180; i++ {
			th := math.Pi * float64(i) / 180
			rows = append(rows, []float64{th, WallRadius(th, 0), WallRadius(th, 0.1), WallRadius(th, 0.3), WallRadius(th, 1)})
		}
		writeCSV(filepath.Join(vizDir, "round_void_shapes.csv"), []string{"theta", "we0", "we0.1", "we0.3", "we1"}, rows)
	}
}
