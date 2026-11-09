// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command ether_cavitation asks whether a nonlinear ether can cavitate
// (proposed H47; plan in CLAUDE.md). Step 1: the ether in one dimension and
// the spherical critical cavity. The soft ether is φ = a·n with an
// amplitude a and an equation of state V(a) (eos.go); κ2 = 1 and λ = 1 here,
// so lengths are healing lengths ξ and energies κ2-units: every result
// scales to any β = λκ4/κ2² by lengths ∝ ξ = ℓ/√β, wall tensions ∝ √β
// (E0/ℓ²), twist rates ∝ √β (1/ℓ), 3-D barriers ∝ 1/√β (E0) at fixed ε/λ.
//
//	go run ./examples/ether_cavitation            # all of step 1 (~1 min)
//	go run ./examples/ether_cavitation -part 3    # one part
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

func main() {
	part := flag.Int("part", 0, "1 walls, 2 tension, 3 reversal, 4 bubble, 5 Kron twin (0 = all)")
	flag.Parse()
	run := func(p int) bool { return *part == 0 || *part == p }
	fail := func(err error) {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
	fmt.Println("ether_cavitation, step 1 (κ2 = 1, λ = 1: lengths in ξ, energies in κ2 units)")
	if run(1) {
		if err := part1(); err != nil {
			fail(err)
		}
	}
	if run(2) {
		part2()
	}
	if run(3) {
		if err := part3(); err != nil {
			fail(err)
		}
	}
	if run(4) {
		if err := part4(); err != nil {
			fail(err)
		}
	}
	if run(5) {
		if err := part5(); err != nil {
			fail(err)
		}
	}
}

func part1() error {
	fmt.Println("\n== 1. Surface tension of the ether (1-D walls) ==")
	for _, h := range []float64{0.1, 0.05, 0.025} {
		n := int(math.Round(20/h)) + 1
		w, err := QCavityWall(1, 20, n)
		if err != nil {
			return err
		}
		s, x0, err := SCoexistenceWall(1, 12, int(math.Round(24/h))+1)
		if err != nil {
			return err
		}
		fmt.Printf("  h = %.3f ξ:  Q cavity wall σ = %.7f (exact 2/3 = %.7f, rel %+.1e, profile %.1e) | S coexistence σ = %.7f (exact 1/2, rel %+.1e, profile %.1e, centre %.1e ξ)\n",
			h, w.Sigma, w.SigmaExact, w.Sigma/w.SigmaExact-1, w.ProfileErr, s.Sigma, s.Sigma/s.SigmaExact-1, s.ProfileErr, x0)
	}
	fmt.Printf("  quadrature check: Q %.12f, S(ε=0) %.12f\n", SigmaQuad(Q(1), 0, 1, 0), SigmaQuad(S(1, 0), 0, 1, 0))
	return nil
}

func part2() {
	fmt.Println("\n== 2. Tension: n twisted at rate k (tensile strength) ==")
	kc, kj, jc := Spinodal(Q(1))
	fmt.Printf("  Q: spinodal k_c = %.6f/ξ (Eckhaus) , max j at k = %.6f/ξ, exact 1/√6 = %.6f; j_c = %.6f (exact 4/(3√6) = %.6f), a0² = %.4f (2/3)\n",
		kc, kj, 1/math.Sqrt(6), jc, 4/(3*math.Sqrt(6)), 1-2*kc*kc)
	fmt.Println("  S (sextic): binodal k_B (cavity phase as cheap as the twisted ether) and spinodal k_c")
	for _, eps := range []float64{0.01, 0.05, 0.1, 0.2} {
		e := S(1, eps)
		kc, kj, jc := Spinodal(e)
		kb := Binodal(e)
		a0, _ := e.Twist(kc).Ether()
		fmt.Printf("    ε = %.2f λ: k_B = %.4f/ξ (≈ √ε = %.4f), k_c = %.4f/ξ (j-max %.4f), j_c = %.4f, a0(k_c) = %.4f\n",
			eps, kb, math.Sqrt(eps), kc, kj, jc, a0)
	}
	fmt.Println("  Ring relaxation (U(1): n turns in a plane; 400 nodes, h = 0.5 ξ, L = 200 ξ, noise 1e-3, seed 1)")
	for _, w := range []int{12, 14} {
		k := 2 * math.Pi * float64(w) / 200
		r := Ring(Q(1), 400, 0.5, w, 2, 1e-3, 3000, 1)
		fmt.Printf("    Q, winding %d (k = %.4f = %.2f k_c): winding %d → %d, min a %.4f, first a < 0.05 at t = %.0f, E %.3f → %.3f, ledger %.1e\n",
			w, k, k/kc, r.W0, r.W1, r.MinA, r.SlipT, r.E0, r.E1, r.LedgerRel)
	}
	fmt.Println("  Ring relaxation (S²: n free to leave the plane, as Faddeev's ether)")
	for _, w := range []int{6, 14} {
		r := Ring(Q(1), 400, 0.5, w, 3, 1e-3, 3000, 1)
		fmt.Printf("    Q, winding %d: winding %d → %d, min a %.4f, E %.3f → %.3f, ledger %.1e\n",
			w, r.W0, r.W1, r.MinA, r.E0, r.E1, r.LedgerRel)
	}
}

func part3() error {
	fmt.Println("\n== 3. Reversal squeezed into a slab (n: +n̂ → −n̂ over L) ==")
	for _, e := range []EOS{Q(1), S(1, 0.05), S(1, 0.2), S(1, 0.5)} {
		r, err := CriticalLengths(e, 0.02)
		if err != nil {
			return err
		}
		name := "Q"
		if e.Kind == 'S' {
			name = fmt.Sprintf("S ε=%.2f", e.Eps)
		}
		lc := fmt.Sprintf("%.4f ξ", r.Lc)
		if r.Lc > 29.9 {
			lc = "> 30 ξ (search limit)"
		}
		fmt.Printf("  %-9s L_c (cavity sheet unstable above) = %s, L_t (bending ether tears below) = %.4f ξ; at L_t: E_turned %.4f, E_straight %.4f\n",
			name, lc, r.Lt, r.ETurned, r.EStraight)
	}
	return nil
}

func part4() error {
	fmt.Println("\n== 4. Critical cavity bubble in 3-D (sextic ether) ==")
	fmt.Println("  ether metastable (ε < 0), no twist: against thin wall R* = 2σ/Δ, E* = 16πσ³/3Δ²")
	for _, eps := range []float64{-0.1, -0.05, -0.02, -0.01} {
		b, err := CriticalBubble(S(1, eps), 0.05)
		if err != nil {
			return err
		}
		fmt.Printf("    ε = %+.2f: R* = %8.3f ξ (thin wall %8.3f, ratio %.4f), E* = %10.3f (thin wall %10.3f, ratio %.4f), negative modes %d\n",
			eps, b.R, b.RTW, b.R/b.RTW, b.E, b.ETW, b.E/b.ETW, b.Neg)
	}
	fmt.Println("  ether stable (ε = 0.05 λ) under twist k (local approximation): binodal → spinodal")
	e := S(1, 0.05)
	kb, kc := Binodal(e), 0.0
	kc, _, _ = Spinodal(e)
	for _, f := range []float64{0.2, 0.4, 0.6, 0.8, 0.95} {
		k := kb + f*(kc-kb)
		b, err := CriticalBubble(e.Twist(k), 0.05)
		if err != nil {
			fmt.Printf("    k = %.4f: %v\n", k, err)
			continue
		}
		fmt.Printf("    k = %.4f/ξ (%.0f %% of the way): a0 = %.4f, Δ = %.5f, R* = %7.3f ξ, E* = %9.3f (thin wall %9.3f), negative modes %d\n",
			k, 100*f, b.A0, b.Delta, b.R, b.E, b.ETW, b.Neg)
	}
	return nil
}

func part5() error {
	fmt.Println("\n== 5. Kron twin: the ladder as a kronspice circuit (OP) ==")
	// Q cavity wall
	{
		n := 201
		h := 20.0 / float64(n-1)
		l := Planar(n, h)
		e := Q(1)
		a := make([]float64, n)
		for i := range a {
			a[i] = math.Tanh(float64(i) * h / 2)
		}
		a[0], a[n-1] = 0, 1
		if _, err := l.Solve(e, a, 1e-12, 200); err != nil {
			return err
		}
		k, err := KronTwin(l, e, a, 0)
		if err != nil {
			return err
		}
		fmt.Printf("  Q cavity wall (%d nodes): max |V − a| = %.1e, σ circuit %.7f vs Go %.7f, ledger %.1e, model-supplied %.3e W\n",
			k.Nodes, k.MaxDiff, k.Energy, l.Energy(e, a, 0), k.LedgerResidual, k.ModelSupplied)
	}
	// S critical bubble (a saddle)
	{
		e := S(1, -0.1)
		b, err := CriticalBubble(e, 0.1)
		if err != nil {
			return err
		}
		l := Radial(len(b.Profile), b.H)
		k, err := KronTwin(l, e, b.Profile, e.V(b.A0))
		if err != nil {
			return err
		}
		fmt.Printf("  S critical bubble ε = −0.1 (%d nodes, a saddle): max |V − a| = %.1e, E* circuit %.5f vs Go %.5f, ledger %.1e, model-supplied %.3e W\n",
			k.Nodes, k.MaxDiff, k.Energy, b.E, k.LedgerResidual, k.ModelSupplied)
	}
	return nil
}
