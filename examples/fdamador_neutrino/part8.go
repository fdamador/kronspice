// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The report for Part 8 (H41): the flavours as three spin states.

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"time"

	"github.com/fdamador/kronspice/spacetime"
)

// Part8 holds what the report computed (also for viz8.go).
type Part8 struct {
	Indep     TripletIndep
	Recip     TripletFit
	Rot       TripletFit
	RecipRun  TripletRun
	RotRim    SagnacRim
	RecipCap  TripletFit // the same with the leak held ≤ 10 %
	Torque    TorqueRun
	FreeJ     float64
	Field     []TripletFieldResult
	RotMasses []RotMass
	OmRatio   float64 // rotation that sets Δm²31/Δm²21 (field frequencies, finest grid; or the ring line's)
	OmRatioW  [2]float64
	Seconds   float64
}

func runPart8(field bool, cellsList []float64, viz string) error {
	t0 := time.Now()
	var p Part8
	fmt.Println("\n8. The flavours as three spin states (H41, user hypothesis 2026-11-07): \"turning one way with full magnetic field,")
	fmt.Println("   not turning (no magnetic field), turning the other way with full magnetic field\"; and if the whole unit also turns")
	fmt.Println("   from outside conditions, Coriolis and centrifugal terms add energy that could account for some of the mass.")
	fmt.Println("   Testable form: the flavours are the azimuthal states m = +1, 0, −1 of the ring (a spin-1-like triplet), labelled by")
	fmt.Println("   the angular momentum J_z about the axis (what a free ring conserves), not by the instantaneous B (which passes")
	fmt.Println("   through zero twice a period in every standing pattern). Ring line (Parts 6, 7): m = ±1 = the CW/CCW harmonic-1 waves")
	fmt.Println("   z_j ∝ e^{∓iκ₁j}, m = 0 = the uniform state z_j ∝ 1; field solver: m = 0 = the E-ring (TE01δ), m = ±1 = Part 7b's")
	fmt.Println("   turning patterns. J = Σ m·U_m/ω_m (the generator of rotations: harmonic m's phase moves by m·α).")

	// 8a
	ti, err := TripletIndependence(50)
	if err != nil {
		return err
	}
	p.Indep = ti
	fmt.Println("8a. Independence. Gram matrices ⟨a|b⟩ of the unit z-patterns on the 24-section ring line (harmonic 1, κ₁ = 2π/24):")
	fmt.Printf("    triplet m = +1, 0, −1: eigenvalues %.4f, %.4f, %.4f (largest off-diagonal |⟨a|b⟩| %.1e): three dimensions\n",
		ti.GramTriplet[0], ti.GramTriplet[1], ti.GramTriplet[2], ti.MaxOffTriplet)
	fmt.Printf("    Part 6's types C, S, D:  eigenvalues %.4f, %.4f, %.4f (|⟨C|S⟩| = |⟨C|D⟩| = %.4f): two dimensions, C = (S + iD)/√2\n",
		math.Abs(ti.GramTypes[0]), ti.GramTypes[1], ti.GramTypes[2], ti.TypeCS)
	fmt.Println("    The triplet differs because m = 0 is a new pattern (another harmonic), not a superposition of the ±1 waves.")
	fmt.Println("    Persistence on the perfect ring (kronspice, 50 periods each at the state's own ω, 60 steps a period):")
	for a, q := range ti.Persist {
		fmt.Printf("      %-18s ω = %.5g rad/s: own share ≥ %.10f; J·ω₁/U start %+.6f, end %+.6f; ledger %.1e\n",
			q.Name, ti.Omega[a], q.MinOwn, q.Jw0, q.JwT, q.LedgerRel)
	}
	fmt.Printf("    Frequencies: m = ±1 share one (degenerate, %.5g MHz); m = 0 is the uniform state at 1/√(L_g·C_g) = %.5g MHz,\n",
		ti.Omega[0]/(2*math.Pi)/1e6, ti.Omega[1]/(2*math.Pi)/1e6)
	fmt.Printf("    %.3f× higher: on the CL line the uniform state does not load the series C, so only the 1 pF stray sets it.\n", ti.Omega[1]/ti.Omega[0])

	// 8b
	fmt.Println("8b. Mixing. With series-capacitor defects alone (Parts 6, 7) m = 0 never mixes: the uniform pattern is an exact")
	r0 := DefaultFRing(false)
	r0.Def = map[int]float64{2: 0.2064, 23: 0.1628}
	_, v0 := r0.Modes()
	var unif float64 // largest overlap of a mode with the uniform pattern
	for i := range v0[0] {
		var s float64
		for j := range v0 {
			s += v0[j][i] / math.Sqrt(float64(r0.N))
		}
		unif = math.Max(unif, s*s)
	}
	fmt.Printf("    eigenvector of the capacitance matrix for any series defects (Cmat·1 = C_g·1); with Part 6's two defects one mode is\n")
	fmt.Printf("    the uniform pattern to %.1e. Coupling it needs the shunt capacitance (node to ground, FRing.NodeC) to vary round the\n", 1-unif)
	fmt.Println("    ring: ΔC_j = C·(M1·cos κ₁j + M2·cos(2κ₁j + 2ψ)), the m = 1 deformation (an off-centre ring: couples 0 to ±1) and the")
	fmt.Println("    m = 2 deformation (an elliptical ring turned by ψ: couples +1 to −1). Knobs M1, M2, ψ and the stray C_g (maps the")
	fmt.Println("    capacitance spectrum onto ω); each of the six assignments of m = +1, 0, −1 to e, μ, τ; cost = squared deviations in")
	fmt.Println("    measured σ (NuFIT 6.0, NO) of the three angles and of Δm²31/Δm²21 in ω² (mass²); leak = largest share of a mass state")
	fmt.Println("    outside the triplet.")
	lab := []string{"+1", "0", "−1"}
	t1 := time.Now()
	p.Recip = FitTriplet(false, 1e3)
	fmt.Printf("    Reciprocal ring (no rotation), %d evaluations, %.0f s: best cost per assignment (e, μ, τ as m):\n", p.Recip.Evals, time.Since(t1).Seconds())
	for pi, pm := range tripletPerms {
		fmt.Printf("      e = %-2s μ = %-2s τ = %-2s : %.4g\n", lab[pm[0]], lab[pm[1]], lab[pm[2]], p.Recip.PerAssign[pi])
	}
	printTriplet("best (reciprocal)", p.Recip.Best)
	fmt.Println("    Without rotation the ring is reciprocal: every mode is real, so its m = +1 and m = −1 projections are complex")
	fmt.Println("    conjugates. With e = m = 0 this is μ–τ reflection symmetry (U_τi = U*_μi; Harrison & Scott, Phys. Lett. B 547, 219,")
	fmt.Printf("    2002), which forces θ23 = 45° and δ_CP = 90° or 270° whatever the knobs (here ||U_μi| − |U_τi|| ≤ %.1e). Any other\n", p.Recip.Best.MuTauAsym())
	fmt.Println("    assignment gives |U_e i| = |U_μ i| or |U_e i| = |U_τ i| for every i, which cannot hold θ13 small. Measured θ23 = 43.3 ± 0.9°")
	fmt.Println("    (45° at 1.9σ); δ_CP = 212° (+26/−41): 270° is 2.2σ above, 90° 3.0σ below. Knobs: 4 continuous (M1, M2, ψ, C_g)")
	fmt.Println("    against 3 free targets (θ12, θ13, ratio): a fit, not a prediction, except θ23 = 45° and δ_CP = ±90°, which come")
	fmt.Println("    from the symmetry. (Angles from the column-normalized projections as in Parts 6, 7; the U printed is the Löwdin-")
	fmt.Println("    unitarized one, which with a large leak differs from them, e.g. |U_e3| above.)")
	fmt.Printf("    But the fit needs deformations of half the series capacitance (M1 = %.2f, M2 = %.2f), and then %.0f %% of a mass state\n",
		p.Recip.Best.Knobs.M1, p.Recip.Best.Knobs.M2, 100*p.Recip.Best.Leak)
	fmt.Println("    lies outside the triplet (in m = ±2 and above): the mass states are no longer mainly the three spin states.")
	t1 = time.Now()
	p.RecipCap = FitTriplet(false, -0.10)
	fmt.Printf("    The same search with the leak held ≤ 10 %% (%.0f s): best cost per assignment\n", time.Since(t1).Seconds())
	for pi, pm := range tripletPerms {
		fmt.Printf("      e = %-2s μ = %-2s τ = %-2s : %.4g\n", lab[pm[0]], lab[pm[1]], lab[pm[2]], p.RecipCap.PerAssign[pi])
	}
	printTriplet("best (leak ≤ 10 %)", p.RecipCap.Best)
	fmt.Println("    The ring line's m = ±2 states lie only about four times farther from m = 0 than m = ±1 (dispersion sin²(κ/2)), so a")
	fmt.Println("    deformation strong enough to bring m = 0 near a ±1 combination (needed for θ12 with a close pair) also pulls in ±2.")
	rr, err := RunTriplet(p.Recip.Best, 1)
	if err != nil {
		return err
	}
	p.RecipRun = rr
	fmt.Printf("    kronspice at the reciprocal fit, started in e (m = %s, u uniform, g = 0), %.0f periods (one slow beat), 120 steps a period:\n",
		lab[rr.StartFlavourM], rr.Periods)
	fmt.Printf("      against the exact expansion over all %d modes (same readout): largest deviation %.2e; against the 3×3 formula %.3f\n",
		r0.N, rr.MaxDevAll, rr.MaxDev3)
	fmt.Printf("      (the formula leaves out the %.0f %% leak); the m = 0 share swings between %.3f and %.3f; |J|·ω₁/U reaches %.3f; ledger %.1e\n",
		100*p.Recip.Best.Leak, rr.MinM0, rr.MaxM0, rr.ExchangeJ, rr.LedgerRel)
	t1 = time.Now()
	p.Rot = FitTriplet(true, 1e3, p.Recip.Best)
	fmt.Printf("    With the ring turning (Part 7's gyrators on the κ₁ pair, ν = Δω/ω) and δ_CP in the cost, from the reciprocal fit\n")
	fmt.Printf("    with eight starting ν (%d evaluations, %.0f s):\n", p.Rot.Evals, time.Since(t1).Seconds())
	printTriplet("best (rotating)", p.Rot.Best)
	r1 := DefaultFRing(false)
	f1 := r1.OmegaK(r1.Kappa(1)) / (2 * math.Pi)
	p.RotRim = Rim(p.Rot.Best.Nu, C0, bench7a, f1)
	fmt.Printf("    Rotation breaks μ–τ symmetry and moves θ23 off 45°; δ_CP stays by −90° (%.0f°). Knobs 5 against 5 targets. The ν needed,\n", p.Rot.Best.DeltaCP)
	fmt.Printf("    %.3g, is a rim speed aΩ = νc²/2v_g = %.3g m/s on the bench ring", math.Abs(p.Rot.Best.Nu), p.RotRim.RimSpeed)
	if p.RotRim.BelowC {
		fmt.Printf(" (a = 10 mm, v_g = c, the most favourable guide): %.3g× the T1000G burst speed.\n", p.RotRim.BurstMargin)
	} else {
		fmt.Println(" (a = 10 mm, v_g = c): at or above c, impossible.")
	}
	fmt.Println("    Standard prediction: flavour states are defined by the charged lepton of the weak interaction; the PMNS matrix")
	fmt.Println("    is measured, not derived; mixing between states of different J_z needs a torque from outside (8c).")

	// 8c
	tq, err := Torque(0.1, 0.5)
	if err != nil {
		return err
	}
	p.Torque = tq
	fr, err := FreeRingJ(200)
	if err != nil {
		return err
	}
	p.FreeJ = fr
	fmt.Println("8c. Angular-momentum ledger. +1 started alone (CW) on the ring with one defect δ = 0.1 on section 0 (half a beat, the")
	fmt.Println("    full conversion to −1). J of the field read from the run; the torque the defect exerts on the field computed")
	fmt.Println("    independently as −∂H/∂α (H = Σ q_j²/c_j with the charges turned by α against the defect), integrated over time:")
	fmt.Printf("      J·ω₁/U₀: start %+.5f, at full conversion %+.5f, ΔJ = %+.5f; ∫ torque dt = %+.5f (same units); largest |ΔJ(t) − ∫τ dt| %.1e\n",
		tq.J0, tq.JMin, tq.DeltaJ, tq.TorqueAt, tq.MaxDiff)
	fmt.Printf("      beat %.5g Hz; energy lost meanwhile %.2e (the 1 mΩ shunt resistance); ledger %.1e\n", tq.Df, tq.EnergyLoss, tq.LedgerRel)
	fmt.Println("    So the 2 units of J (2ħ per quantum) go to the defect, which takes the angular momentum and, being fixed, no energy")
	fmt.Println("    (torque × its angular velocity = 0). Without the defect the +1 state keeps its J:")
	fmt.Printf("      perfect ring, 200 periods: largest |J·ω₁/U − 1| = %.1e.\n", fr)
	fmt.Println("    A free ring in vacuum has nothing to push against: rotation symmetry makes J_z exactly conserved, so +1 → −1 (ΔJ_z = 2ħ)")
	fmt.Println("    and +1 → 0 (ħ) are forbidden unless something outside takes the difference. Measured neutrino oscillation keeps the")
	fmt.Println("    helicity (−½ for every flavour) in vacuum, with no partner to recoil. Standard prediction: flavour oscillation is")
	fmt.Println("    interference of mass states with the same spin; it exchanges no angular momentum.")

	// 8d
	burst := spacetime.BurstSpeed(spacetime.T1000G)
	e := DefaultERing()
	b := DefaultBench().Run()
	if field {
		fmt.Println("8d. The states in the field solver (em/fdtd.go, the torus of 5b/7b), from filtered phasors normalized to Ū = 1 J:")
		for _, c := range cellsList {
			st := time.Now()
			sf := DefaultSpinField()
			sf.Ring.Cells = c
			res, err := RunTripletField(sf, 48, func(f string, a ...any) {
				fmt.Fprintf(os.Stderr, "[8d r/%g %s] "+f+"\n", append([]any{c, time.Since(st).Round(time.Second)}, a...)...)
			})
			if err != nil {
				return err
			}
			res.Seconds = time.Since(st).Seconds()
			p.Field = append(p.Field, res)
			reportTripletField(res, e, b)
		}
		fmt.Println("    Every state's moment averages to zero over a period: a single-frequency linear field has no static (rectified)")
		fmt.Println("    current, so Borexino's bound on a static moment (< 2.8e-11 μB) is met trivially. Against the user's description: the")
		fmt.Println("    turning states do keep a constant magnetic energy (U_B = Ū/2 at every instant: B never vanishes as a whole), but they")
		fmt.Println("    carry no magnetic dipole; what turns with them is an electric dipole of constant size. The m = 0 state is the one")
		fmt.Println("    with a magnetic moment, along the axis, and it is not without B: its B and its moment pass through zero twice a")
		fmt.Println("    period, half the energy magnetic on average. Standard prediction: Dirac-neutrino moment 3eG_F m_ν/(8√2π²) ≈ 3e-20 μB")
		fmt.Println("    (m/0.1 eV), static; no electric dipole without CP violation beyond the PMNS phase.")
	} else {
		fmt.Println("8d. Field solver: run with -field; see run_output.txt for the recorded run.")
	}

	// 8e
	fmt.Printf("8e. Outside rotation of the whole unit (the user's extension). In the lab a state of order m in a ring turning at Ω has\n")
	fmt.Println("    E = E′ + Ω·J_z (Landau & Lifshitz, Mechanics, §39): per quantum the Coriolis (first-order) shift is m·ħΩ, i.e. mΩ/ω of")
	fmt.Println("    the energy (Part 7c); the centrifugal (second-order) term is ½IΩ² with the field's moment of inertia I = (ħω/c²)⟨ρ²⟩.")
	fmt.Println("    Target 0.1 eV (declared; KATRIN m < 0.45 eV). Burst limit: T1000G hoop, 1,881 m/s (spacetime.BurstSpeed).")
	rho2 := [3]float64{1, 1, 1}
	if len(p.Field) > 0 {
		fl := p.Field[len(p.Field)-1]
		for a := 0; a < 3; a++ {
			rho2[a] = fl.States[a].Rho2 / (e.A * e.A)
		}
		fmt.Printf("    ⟨ρ²⟩/a² from the field solver (r/%g): m = +1 %.3f, m = 0 %.3f, m = −1 %.3f\n", fl.Cells, rho2[0], rho2[1], rho2[2])
	} else {
		fmt.Println("    ⟨ρ²⟩ = a² (thin ring; the field solver's values with -field)")
	}
	fmt.Printf("    %-42s %3s %11s %11s %11s %9s %11s %11s %12s %12s\n", "ring", "m", "ħω [eV]", "Ω_cor", "aΩ_cor", "Ω_cor/ω", "Ω_cen", "aΩ_cen",
		"cor @ burst", "cen @ burst")
	rings := []struct {
		name  string
		a, w  float64
		wPlus float64
	}{
		{"bench ring (a = 10 mm, εr = 100)", e.A, e.Omega0(), 0},
		{fmt.Sprintf("Part 2 ring (a = %.3g m, εr_eff %.2g)", b.RadiusNeeded, b.EpsRNeeded), b.RadiusNeeded, b.Omega0, 0},
	}
	if len(p.Field) > 0 {
		fl := p.Field[len(p.Field)-1]
		rings[0].w = 2 * math.Pi * fl.FE
		rings[0].wPlus = 2 * math.Pi * fl.F1
	}
	for _, rg := range rings {
		for a, m := range []int{1, 0, -1} {
			w := rg.w
			if m != 0 && rg.wPlus > 0 {
				w = rg.wPlus
			}
			rm := OutsideRotation(rg.name, m, w, rg.a, rho2[a], 0.1, burst)
			p.RotMasses = append(p.RotMasses, rm)
			cor, corR, corW := "—", "—", "—"
			if m != 0 {
				cor, corR, corW = fmt.Sprintf("%.3g", rm.OmCor), fmt.Sprintf("%.3g", rm.RimCor), fmt.Sprintf("%.3g", rm.CorOverOmega)
			}
			fmt.Printf("    %-42s %+3d %11.3g %11s %11s %9s %11.3g %11.3g %12.3g %12.3g\n", rg.name, m, rm.QuantumEV, cor, corR, corW,
				rm.OmCen, rm.RimCen, rm.CorAtBurst, rm.CenAtBurst)
		}
	}
	fmt.Println("    (Ω in rad/s, aΩ in m/s; 'cor/cen @ burst' = eV per quantum at the burst speed v_b = aΩ_b. The bench ring's m = ±1")
	fmt.Println("    quantum is the m = 1 mode's when the field solver ran.)")
	fmt.Printf("    Coriolis: Ω_cor = 0.1 eV/ħ = %.4g rad/s for every ring size. On the bench ring that is %.3g× its own ω and a rim speed\n",
		0.1*E/Hbar, 0.1*E/Hbar/e.Omega0())
	fmt.Printf("    of %.2g c. On the Part 2 ring Ω_cor equals its own ω₀ (ħω₀ = 0.1 eV by construction): the rim moves at %.3g m/s, below the\n",
		0.1*E/Hbar*e.A/C0, 0.1*E/Hbar*b.RadiusNeeded)
	fmt.Printf("    burst speed, but Ω/ω = 1 is outside the first-order form and the ring is %.0f× smaller than a hydrogen atom: no material\n", 5.29177e-11/b.RadiusNeeded)
	fmt.Println("    ring of that size exists to turn, and any material ring carries a rest mass far above 0.1 eV.")
	fmt.Println("    Centrifugal: ½IΩ² = ½ħω·⟨ρ²⟩Ω²/c² < ½ħω·(v/c)² < ½ħω for any rim speed below c. It can add at most half the quantum's")
	fmt.Println("    own energy, so it cannot supply a mass the field does not already have; on the bench ring the rim would need")
	fmt.Printf("    %.0f c. The m = 0 state has J_z = 0 and gets no first-order shift; only its tiny centrifugal term.\n", p.RotMasses[1].RimCen/C0)
	// pattern
	w0, w1 := e.Omega0(), e.Omega0()*2.155
	src := "Part 1's lumped E-ring and an assumed ratio 2.155"
	if len(p.Field) > 0 {
		fl := p.Field[len(p.Field)-1]
		w0, w1 = 2*math.Pi*fl.FE, 2*math.Pi*fl.F1
		src = fmt.Sprintf("the field solver's m = 0 and m = 1 frequencies (r/%g)", fl.Cells)
	}
	m := NuFIT6()
	p.OmRatioW = [2]float64{w0, w1}
	p.OmRatio = RatioOmega(w0, w1, m.Ratio)
	fmt.Println("    Mass pattern under rotation: E(±1) = ħ(ω₁ ± Ω), E(0) = ħω₀ (+ centrifugal): the ±1 pair splits symmetrically, m = 0")
	fmt.Printf("    stays. With %s, %.4g and %.4g GHz,\n", src, w0/(2*math.Pi)/1e9, w1/(2*math.Pi)/1e9)
	fmt.Printf("    the measured Δm²31/Δm²21 = %.2f needs Ω/ω₁ = %.3g. The close pair is then the turning pair, heavier than m = 0 by\n", m.Ratio, p.OmRatio/w1)
	fmt.Printf("    ω₁/ω₀ = %.3f: an inverted ordering (NuFIT 6.0 favours the normal one). Scaling the ring so that Δm²31 = 2.513e-3 eV²\n", w1/w0)
	mm0 := math.Sqrt(2.513e-3 / (sq(w1/w0) - 1))
	fmt.Printf("    gives m(0) = %.3g eV and m(±1) ≈ %.3g eV. The rim speed aΩ = (Ω/ω₁)·(aω₁) is the same for every size at εr = 100:\n", mm0, mm0*w1/w0)
	vr := p.OmRatio / w1 * w1 * e.A
	fmt.Printf("    %.3g m/s on the bench ring, %.0f× the burst speed, (v/c)² = %.2g.\n", vr, vr/burst, sq(vr/C0))
	fmt.Println("    Ledger: the rotation's energy (Ω·J_z of the wave, plus ½I_matΩ² of any material ring) is work done by the outside")
	fmt.Println("    agent that turns the unit, a declared mechanical input column, never inflow from nowhere; at constant Ω the")
	fmt.Println("    Coriolis terms do no work (Part 7a: gyrators' net work 6e-19 of the stored energy). Rest mass: a ring turning")
	fmt.Println("    about its own centre has p = 0, so its invariant mass Mc² = √(E² − p²c²) = E includes the rotation energy, as a")
	fmt.Println("    spinning flywheel weighs more by E_rot/c². It is not intrinsic: it is set by Ω, which is frame energy put in from")
	fmt.Println("    outside and held only while J is held. Standard prediction: neutrino masses are not explained by the Standard Model")
	fmt.Println("    (Dirac masses via the Higgs, or a seesaw with heavy partners are the leading proposals); measured sum and splittings")
	fmt.Println("    constrain them; no rotation of the neutrino enters.")
	fmt.Println("    Status H41: bounded. The triplet is three-dimensional (unlike Part 6's types) and its m = 0 state is B-free only at")
	fmt.Println("    two instants a period. Mixing it needs a torque from a structure (shunt defect); a free ring conserves J_z. With e = 0")
	fmt.Println("    the reciprocal ring gives θ23 = 45°, δ_CP = ±90° (μ–τ symmetry); the rest is a fit that needs deformations so large")
	fmt.Println("    that about half of each mass state lies outside the triplet. Its J_z = ±ħ, 0 are whole units,")
	fmt.Println("    not ½; spin sign separates neutrino from antineutrino in measurement (all flavours helicity −½), not flavours.")
	fmt.Println("    Rotation adds mΩħ per quantum (none to m = 0) and at most ½ħω centrifugally; the splitting ratio needs rim speeds")
	fmt.Println("    hundreds of times the burst speed. Every joule of it is booked to the outside agent.")
	p.Seconds = time.Since(t0).Seconds()
	fmt.Printf("   (Part 8: %s)\n", time.Since(t0).Round(time.Second))
	if viz != "" {
		if err := writeViz8(viz, p); err != nil {
			return err
		}
		fmt.Println("\nviz data for Part 8 written to", viz, "— run python3 viz/render_part8.py and python3 viz/plot_triplet.py")
	}
	return nil
}

func printTriplet(name string, c TripletCl) {
	m := NuFIT6()
	lab := []string{"+1", "0", "−1"}
	k := c.Knobs
	fmt.Printf("    %s: e, μ, τ = m %s, %s, %s; knobs M1 = %.4f, M2 = %.4f, ψ = %.4f rad, C_g = %.4g F, ν = %.4g\n",
		name, lab[c.Assign[0]], lab[c.Assign[1]], lab[c.Assign[2]], k.M1, k.M2, k.Psi, k.Cg, k.Nu)
	fmt.Printf("      θ12 %.2f° (%.2f ± %.2f), θ23 %.2f° (%.1f ± %.1f), θ13 %.2f° (%.2f ± %.2f), Δm²31/Δm²21 %.2f (%.2f; in ω: %.2f),\n",
		c.Th12, m.Th12, m.S12, c.Th23, m.Th23, m.S23, c.Th13, m.Th13, m.S13, c.Ratio2, m.Ratio, c.Ratio)
	fmt.Printf("      δ_CP %.1f° (212 +26/−41), J %.5f (J_max %.5f); leak %.3f; cost %.3g (with δ_CP %.3g)\n", c.DeltaCP, c.J, c.JMax, c.Leak, c.Cost2, c.CostCP2)
	fmt.Println("      U (rows e, μ, τ; columns ν1, ν2, ν3), |U|∠arg U:")
	for a := 0; a < 3; a++ {
		fmt.Printf("        %.4f∠%+7.1f°  %.4f∠%+7.1f°  %.4f∠%+7.1f°\n", cmplx.Abs(c.UU[a][0]), cmplx.Phase(c.UU[a][0])*180/math.Pi,
			cmplx.Abs(c.UU[a][1]), cmplx.Phase(c.UU[a][1])*180/math.Pi, cmplx.Abs(c.UU[a][2]), cmplx.Phase(c.UU[a][2])*180/math.Pi)
	}
}

func reportTripletField(res TripletFieldResult, e ERing, b BenchResult) {
	fmt.Printf("    grid r/%g (dx = %.3g mm): m = 0 at %.5g GHz, m = ±1 at %.5g GHz (ratio %.4f); %.0f s\n", res.Cells, res.Dx*1e3, res.FE/1e9, res.F1/1e9, res.F1/res.FE, res.Seconds)
	fmt.Printf("      Gram eigenvalues of (+1, 0, −1): %.4f, %.4f, %.4f; |⟨+1|−1⟩| %.1e, |⟨+1|0⟩| %.1e, |⟨0|−1⟩| %.1e\n",
		res.GramEig[0], res.GramEig[1], res.GramEig[2], cmplx.Abs(res.Gram[0][2]), cmplx.Abs(res.Gram[0][1]), cmplx.Abs(res.Gram[1][2]))
	fmt.Printf("      %-28s %3s %9s %9s %15s %11s %11s %11s %16s %16s\n", "state", "m", "J_zω/Ū A", "J_zω/Ū M", "U_B/Ū range",
		"|m_z| [μB]", "|m_⊥| [μB]", "|p|c [μB]", "|m(t)|/max", "|p(t)|/max")
	for _, s := range res.States {
		w := 2 * math.Pi * s.F
		ax, pp := s.MomentPerQuantum()
		fmt.Printf("      %-28s %+3d %+9.4f %+9.4f  [%.3f, %.3f] %11.3g %11.3g %11.3g   [%.3f, %.3f]   [%.3f, %.3f]\n", s.Name, s.M, s.JzAMean*w, s.JzMMean*w,
			s.UBMin, s.UBMax, ax/MuB, pp/MuB, s.DipolePerQuantum()*C0/MuB, s.MRange[0], s.MRange[1], s.PRange[0], s.PRange[1])
	}
	fmt.Println("      (J_z averaged over a period; U_B/Ū and the dipoles' size over a period; moments are amplitudes with one quantum ħω")
	fmt.Println("      stored; the electric dipole p as p·c in μB; a moment below 1e-12 of the other is round-off, and its range is noise)")
	l := e.L()
	i0 := math.Sqrt(2 * Hbar * e.Omega0() / l)
	lump := i0 * math.Pi * e.A * e.A
	ax, _ := res.States[1].MomentPerQuantum()
	fmt.Printf("      twin: the lumped E-ring of Part 1 with one quantum, I₀·πa² (½LI₀² = ħω₀), gives %.3g μB; the field's m = 0 axial moment\n", lump/MuB)
	fmt.Printf("      %.3g μB (ratio %.3f). Per quantum the moment scales ∝ a at fixed shape and εr; the Part 2 ring's lumped value\n", ax/MuB, ax/lump)
	fmt.Printf("      (%.3g μB, Part 2) times this ratio is %.3g μB for m = 0 (oscillating; εr differs there, so a scale, not a result).\n",
		b.MomentAmplitude, b.MomentAmplitude*ax/lump)
}
