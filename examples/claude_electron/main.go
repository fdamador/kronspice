// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command claude_electron tests H44: Williamson & van der Mark's
// toroidal-photon electron (H23) held together by the knotted ether of
// H29, with a Hestenes zitterbewegung electron as its standard twin.
//
// Step 1 relaxes the knots (Rañada's single loop A_{1,1}, the double
// loop A_{1,2} and the Q = 2 minimum A_{2,1}) and follows their degree;
// step 2 quantizes them as rotors (Finkelstein–Rubinstein: spin ½ only
// for odd Q); step 3 sets out the charge options; step 4 fixes κ2, κ4 by
// the electron's mass and a size and compares every prediction with
// measurement; step 5 prints the verdict for the register.
//
//	go run ./examples/claude_electron            # h = 0.1 (about five minutes)
//	go run ./examples/claude_electron -h 0.15    # quick
package main

import (
	"flag"
	"fmt"
	"math"
	"time"
)

type spec struct {
	name    string
	M, K    int
	box     float64
	fermion bool
}

func main() {
	h := flag.Float64("h", 0.1, "grid spacing in ℓ")
	steps := flag.Int("steps", 80000, "relaxation step limit")
	flag.Parse()

	tw := StandardTwin()
	fmt.Println("== Step 0: standard twin (Dirac, read by Hestenes) ==")
	fmt.Printf("zitter radius ħ/2mc = %.3f fm (H23's loop radius λ/4π = %.3f fm: the same)\n", tw.ZitterRadiusFm, WvdMRadiusFm())
	fmt.Printf("zitter frequency 2mc²/ħ = %.4e rad/s, speed c\n", tw.ZitterOmega)
	fmt.Printf("g: Dirac %.0f, measured %.11f; classical charge at c on ħ/2mc: g = %.4f (μ = μ_B/2)\n", tw.GDirac, tw.GQED, tw.GClassicalAtZitter)
	fmt.Printf("g = 2 by classical circulation at c needs radius %.3f fm = ħ/mc, where an orbiting mass m would carry ħ, not ħ/2\n\n", tw.RadiusForG2Fm)

	specs := []spec{
		{"A11 Rañada single loop", 1, 1, 10, true},
		{"A21 (Q = 2 minimum)", 2, 1, 12, false},
	}
	fmt.Printf("== Step 1: confinement (h = %.2f ℓ) ==\n", *h)
	fmt.Println("the double loop A12 (twice round the ring, once round the tube; Q = 2), every 1000 steps:")
	d := NewAxi(10, 10, *h, 1, 1, 1, 1, 2)
	fmt.Printf("  start: degree %.3f, E %.1f, core %.3f ℓ\n", d.Degree(), d.Total(), d.Core())
	for s := 1; s <= 8; s++ {
		d.Relax(1000, 0.05, 0)
		fmt.Printf("  %5d: degree %.3f, E %.1f, core %.3f ℓ\n", 1000*s, d.Degree(), d.Total(), d.Core())
	}
	fmt.Println("  → no minimum in the axial ansatz: the core falls onto the axis and the knot unties through the grid")
	{
		b := NewAxi(10, 10, *h, 1, 0, 1, 1, 2)
		r0 := math.Sqrt(b.Inertia(8).R2)
		b.Relax(3000, 0.05, 0)
		r1 := math.Sqrt(b.Inertia(8).R2)
		fmt.Printf("  control, double loop as knotted light alone (κ2 = 0): rms %.2f → %.2f ℓ in 3000 steps\n", r0, r1)
	}
	fmt.Println("knot                     | Q | degree start→end | E/E0    | E2/E4  | core/ℓ | rms r/ℓ | steps | s")
	var knots []Knot
	for _, sp := range specs {
		t0 := time.Now()
		a := NewAxi(sp.box, sp.box, *h, 1, 1, 1, sp.M, sp.K)
		d0 := a.Degree()
		e2, e4, n := a.Relax(*steps, 0.05, 1e-7)
		d1 := a.Degree()
		mo := a.Inertia(32)
		core := a.Core()
		q := sp.M * int(math.Round(d1))
		fmt.Printf("%-24s | %d | %.3f→%.3f      | %7.2f | %.4f | %.3f  | %.3f   | %5d | %.0f\n",
			sp.name, q, d0, d1, e2+e4, e2/e4, core, math.Sqrt(mo.R2), n, time.Since(t0).Seconds())
		if sp.M == 1 {
			fmt.Printf("  inertia check: I∥ from φ-average %.3f vs direct %.3f\n", mo.IParVia, mo.IParFull)
		}
		knots = append(knots, Knot{
			Name: sp.name, Q: q, Fermion: q%2 == 1,
			E: e2 + e4, Core: core, RmsRho: math.Sqrt(mo.Rho2), RmsR: math.Sqrt(mo.R2),
			IPar: mo.IPar, IPer: mo.IPerp,
		})
	}
	fmt.Printf("E(A21)/E(A11) = %.3f (Vakulenko–Kapitanski scaling Q^¾ = %.3f)\n", knots[1].E/knots[0].E, math.Pow(2, 0.75))

	fmt.Println("\n== Step 2: spin (Finkelstein–Rubinstein; Krusch & Speight 2006: fermion only for odd Q) ==")
	fmt.Println("knot                     | Q | fermion | I∥     | I⊥     | c_r      | ℓ_min (spin ½ in m_e) | g (rigid)")
	for _, k := range knots {
		lmin, c := "—", "—"
		if k.Fermion {
			lmin = fmt.Sprintf("%.1f fm = %.3f λ̄", k.MinSizeSpin(), k.MinSizeSpin()/lbarC)
			c = fmt.Sprintf("%.3e", k.SpinC())
		}
		fmt.Printf("%-24s | %d | %-7v | %6.1f | %6.1f | %-8s | %-21s | %.3f\n", k.Name, k.Q, k.Fermion, k.IPar, k.IPer, c, lmin, k.GRigid())
	}

	fmt.Println("\n== Step 3: charge ==")
	for _, o := range ChargeOptions() {
		fmt.Printf("%-36s q = %.6f e | anomaly %.9f (measured %.11f) | Lyman-α %.2f nm | atom %+.3f e | derived: %v\n",
			o.Name, o.Q, o.Anomaly, aMeas, o.LymanNm, o.NetAtomCharge, o.Derived)
	}
	fmt.Println("gauged ether (Babaev–Faddeev–Niemi): the pair charge is an input, not derived")

	k := knots[0] // the only knot that can carry spin ½
	fmt.Println("\n== Step 4: fix κ2, κ4 (A11) and predict ==")
	lA := WvdMRadiusFm() / k.Core
	for _, r := range []Reading{k.At("A: knot ring = λ̄/2", lA), k.At("B: ℓ = 10⁻¹⁸ m", 1e-3), k.At("ℓ_min (spin)", k.MinSizeSpin())} {
		fmt.Printf("%-22s ℓ = %.4g fm, ring %.4g fm, rms %.4g fm, E0 = %.4g keV\n", r.Name, r.LFm, r.RingFm, r.RmsFm, r.E0MeV*1e3)
		fmt.Printf("   spin energy %.3g m_e c² (real E0 with spin: %v), rim speed %.3g c\n", r.ERotOverM, r.E0SpinOK, r.RimSpeedC)
		fmt.Printf("   tone ħω %.4g MeV (excited electron), Coulomb %.3g m_e c², 1S–2S shift %.4g MHz, g %.3f, ESR %.2f GHz/T\n",
			r.ToneMeV, r.CoulombOverM, r.H1S2SMHz, r.G, r.ESRGHzT)
	}
	fmt.Printf("H23 loop with charge circulating at c: g = %.4f (eq. 12's g = 2(1 + α′/2π) is asserted, not obtained by circulation)\n", GCirculation(WvdMRadiusFm()))
	fmt.Println("measured: no excited electron below 5.6 TeV, size < 10⁻¹⁸ m, g = 2.00232, ESR 28.025 GHz/T")

	fmt.Println("\nwindows on ℓ for A11:")
	ws := k.Windows()
	for _, w := range ws {
		fmt.Printf("  %-46s ℓ ∈ [%.3g, %.3g] fm\n", w.Name, w.Lo, w.Hi)
	}
	lo, hi := Overlap(ws)
	fmt.Printf("  overlap: [%.3g, %.3g] fm → ", lo, hi)
	if lo > hi {
		fmt.Printf("empty; the lower bounds (spin, Coulomb) exceed the upper bounds (scattering, excited states, hydrogen) by %.1e\n", lo/hi)
	} else {
		fmt.Println("non-empty")
	}

	fmt.Println("\n== Step 5: verdict ==")
	fmt.Println("double loop (A12, Q = 2): not held, it unties; the held Q = 2 knot (A21, windings swapped) is bosonic under FR: no spin ½")
	fmt.Println("single loop (A11) can be a fermion, but no size satisfies spin, Coulomb, scattering, excited-state and hydrogen bounds together")
	fmt.Println("charge: not derived (H23 0.910 e, Rañada 3.30 e); status: bounded")
}
