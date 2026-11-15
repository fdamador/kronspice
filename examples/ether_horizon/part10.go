// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
	"math"
)

func part10() {
	fmt.Println()
	fmt.Println("Part 10a. Meyl beside H49: which contraction law holds a void at one size?")
	fmt.Println("  E(R) = ħc·x/R + S(R), units with ħc·x = 1; stable needs a minimum of E")
	fmt.Println("  contraction S(R)                    stationary R      kind")
	for _, c := range Contractions() {
		roots, kinds, neutral := Balance(c)
		switch {
		case neutral:
			fmt.Printf("  %-34s  every R           neutral: no size is preferred\n", c.Name)
		case len(roots) == 0:
			fmt.Printf("  %-34s  none              the void runs away\n", c.Name)
		default:
			for i, r := range roots {
				k := "stable"
				if kinds[i] < 0 {
					k = "unstable"
				}
				fmt.Printf("  %-34s  %-16.4g  %s\n", c.Name, r, k)
			}
		}
	}
	fmt.Println("  Meyl's self-similar field (E, H ∝ 1/L², eq. 2.20) has field energy ∝ 1/R, the same power as the light:")
	fmt.Println("  the two cancel (k = 1) or one wins outright. It sets no size; something that grows with R must.")
	fmt.Println("  A squeeze that weakens faster (−1/R²) gives a balance, but a maximum: the void collapses or flies apart.")

	fmt.Println()
	fmt.Println("Part 10b. Proton and neutron: the 9a void at the nucleon mass vs measurement")
	fmt.Println("  particle  mode   held by    void R [fm]   uniform rms [fm]   measured [fm]")
	for _, p := range []struct {
		name string
		m    float64
		meas float64
		lbl  string
	}{{"proton", Mp, RpCODATA, "0.84075 charge (CODATA 2022)"}, {"neutron", Mn, RnMag, "0.864 magnetic (PDG)"}} {
		for _, md := range []struct {
			n string
			x float64
		}{{"TM₁₁", XTM11}, {"TE₁₁", XTE11}} {
			for _, t := range []bool{false, true} {
				h := "pressure"
				if t {
					h = "tension"
				}
				R := VoidAtMass(p.m, md.x, t)
				fmt.Printf("  %-8s  %-5s  %-9s  %-12.4f  %-17.4f  %s\n", p.name, md.n, h, R*1e15, UniformRMS(R)*1e15, p.lbl)
			}
		}
	}
	fmt.Println("  The TM₁₁ void's edge lands within 10% of the proton's rms radius with no fitted constant — the bag")
	fmt.Println("  model's old success (a cavity mode sized by ħ/mc). Read as an rms of a filled ball it is 20–30% small,")
	fmt.Println("  so the charge would have to sit near the wall. The neutron comes out the same size to 0.1%, as measured.")
	fmt.Println()
	fmt.Println("  Magnetic moments: a void of neutral light has none. Measured μ_p = 2.793 μ_N, μ_n = −1.913 μ_N;")
	fmt.Printf("  μ_n/μ_p measured %.4f, quark/bag model %.4f (%.1f%% off). The void needs charged content to match.\n",
		MuN/MuP, BagMomentRatio, 100*(BagMomentRatio/(MuN/MuP)-1))

	fmt.Println()
	fmt.Println("Part 10c. Meyl's probe-dependent radius vs the 2010–2026 data")
	d, b := ProbeBound()
	fmt.Printf("  Meyl (§3.5): (H_p + H_e) ∝ 1/R², so a muon (heavier, smaller) should see a smaller proton than an electron.\n")
	fmt.Printf("  2010: e–p %.4f fm vs μp %.5f fm, a %.1f%% gap — the proton radius puzzle, which his rule explains in sign.\n",
		RpOld*1e15, RpMuonic*1e15, 100*(RpOld-RpMuonic)/RpMuonic)
	fmt.Printf("  2026: hydrogen 2s–6p %.4f(15) fm vs μp %.5f(39) fm: Δr/r = %+.3f%%, so |Δr/r| < %.2f%% (2σ).\n",
		RpH2026*1e15, RpMuonic*1e15, 100*d, 100*b)
	gap := (RpOld - RpMuonic) / RpMuonic
	fmt.Printf("  The 2010 gap his rule explained is now at least %.0f× smaller than it was. The puzzle closed on the\n", gap/b)
	fmt.Println("  electronic side (new hydrogen data agree with muonic), so any probe dependence is below 0.4%. H49's")
	fmt.Println("  void radius is fixed by the target's mass alone, which is what the data now show.")

	fmt.Println()
	fmt.Println("Part 10d. Onion of light–void shells as coupled LC tanks (kronspice AC)")
	o := Onion{N: 3, Q: 1 / MuOverE, K: 0.01, L: 1e-3, F0: 1e3}
	peaks, err := o.Modes()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("  self-similar shells, radius ratio fitted to m_μ/m_e = 206.77; mass ∝ mode frequency")
	fs := o.Freqs()
	for k, f := range fs {
		p := math.NaN()
		if k < len(peaks) {
			p = peaks[k]
		}
		fmt.Printf("  shell %d   uncoupled %-11.5g Hz   coupled peak %-11.5g Hz   ratio to shell 0  %.2f\n", k, f, p, p/peaks[0])
	}
	pred := MuOverE * MuOverE * 0.51099895
	fmt.Printf("  third shell predicts m_τ ≈ %.1f GeV; measured 1.777 GeV (%.0f× too heavy).\n", pred/1000, pred/1776.93)
	fmt.Printf("  Lepton steps are 206.77 and %.2f, not equal: the onion is not self-similar. Shells would need\n", TauOverM)
	fmt.Println("  their own rule (a Koide-like relation fits all three; a geometric ladder does not).")
}
