// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command void_charge tests proposed hypothesis H51: give the void of a
// nonlinear ether electrical properties (a charge, a pressure) and test them
// against electrical measurements of the electron. It starts from
// Tewari's charged void (H13, re-read: Universal Principles of Space and
// Matter, §2.3, §3.1, §8) and the round light bubble of examples/void_shells
// (H50), and asks four questions:
//
//  1. Tewari's charged void, as he wrote it: field, energy, annihilation;
//  2. charge radius: where can the charge sit, and what do hydrogen and LEP
//     allow (Darwin smearing 0.866 ƛ; form factor at √s = 200 GeV)?
//  3. exact charge: does a void's charge stay fixed and come in units?
//  4. electron and positron: how two voids interact, and what happens when
//     they meet (annihilation into light?);
//  5. stiffness: the void's own vibrations and what they predict;
//  6. the scorecard.
//
// Standard prediction beside it: a pointlike charge −e (< 10⁻¹⁸ m) whose
// Coulomb field is conservative and isotropic, Darwin smearing from Dirac's
// equation, QED vacuum polarization making α grow at short distance, e+e−
// annihilating into two 511 keV photons, no internal excitations.
package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
)

// The run's grid and constants, those of examples/void_shells.
const (
	GridH  = 0.4
	GridLR = 24.0
	GridLZ = 24.0
	V0     = 20.0
	EpsRef = 300.0
	Relax  = 8000
	RLocal = 9.0 // "near the particle": within 9 ξ (the bubble reaches 4.4 ξ)
)

func grid() *Grid { return NewGrid(GridH, GridLR, GridLZ) }

func relaxedBall(eps float64, occ []int, r0 float64) *State {
	s := NewState(grid(), eps, V0, occ, Ball(r0))
	s.Relax(Relax)
	return s
}

var csvOut = map[string][]string{}

func csvAdd(name, header, row string) {
	if len(csvOut[name]) == 0 {
		csvOut[name] = append(csvOut[name], header)
	}
	csvOut[name] = append(csvOut[name], row)
}

func main() {
	viz := flag.String("viz", "", "write CSVs and animation frames to DIR")
	flag.Parse()
	var frameDir string
	if *viz != "" {
		frameDir = filepath.Join(*viz, "frames")
		if err := os.MkdirAll(frameDir, 0o755); err != nil {
			log.Fatal(err)
		}
	}

	fmt.Println("void_charge (proposed H51): electrical properties for the void, tested against the electron")
	fmt.Printf("Ether and light as examples/void_shells: grid %.1f ξ, box ρ < %.0f, |z| < %.0f, V0 = %.0f, ε = %.0f (plays ħc).\n", GridH, GridLR, GridLZ, V0, EpsRef)
	fmt.Println("Standard prediction: a pointlike charge (< 1e-18 m), conservative isotropic Coulomb field, Darwin smearing")
	fmt.Println("0.866 ƛ from Dirac's equation, α growing at short distance (QED), e+e− → two 511 keV photons, no excited electron.")

	part1()
	ball := relaxedBall(EpsRef, []int{0}, 5)
	part2(ball)
	part3(ball, frameDir)
	part4(ball, frameDir)
	part5(ball, frameDir)
	part6()

	if *viz != "" {
		for name, rows := range csvOut {
			if err := os.WriteFile(filepath.Join(*viz, name+".csv"), []byte(strings.Join(rows, "\n")+"\n"), 0o644); err != nil {
				log.Fatal(err)
			}
		}
	}
}

func part1() {
	t := TewariVoid()
	fmt.Println()
	fmt.Println("Part 1. Tewari's charged void, as written (UPSM 2017 §2.3, §2.6, §3.1, §8)")
	fmt.Println("  He puts no charge inside the void: it is fieldless. The charge is the surface integral of the space")
	fmt.Println("  fluid's speed on the void's surface (eq. 2.4, q = π² r_e² c), its field is E = −c² r_e² sin²θ / 2r²")
	fmt.Println("  (eq. 3.1.3), zero along the spin axis; a positron is the same void spinning the other way (§2.3);")
	fmt.Println("  in annihilation the voids collapse and a single shell of light of width r_e leaves (§8.0–8.1).")
	eUnit := Alpha * MC2keV / t.ReLbar // e²/r_e in keV
	fmt.Printf("  void radius r_e = %.3e m = %.3f ƛ (fixed by the measured charge, eq. 2.8)\n", t.ReCm*1e-2, t.ReLbar)
	fmt.Println("  test                               his void                         measured                         verdict")
	fmt.Printf("  field is conservative              work %.2f e²/r_e = %.2f keV per    zero (static fields conserve      fails: a test charge gains\n", t.LoopWork, t.LoopWork*eUnit)
	fmt.Println("    (∮E·dl = 0)                       loop (equator out, axis back)     energy)                           energy every lap")
	fmt.Printf("  hydrogen 2p by spin orientation    levels spread by %.2f eV           degenerate (fine structure        fails by ~1e5\n", t.Ly2pSplit)
	fmt.Println("                                                                       45 µeV); Lyman-α one line")
	r2 := t.ReLbar * t.ReLbar
	fmt.Printf("  charge radius vs hydrogen 2S       ⟨r²⟩ = %.3f ƛ² vs Darwin 0.750:   2S shift known to ±9 kHz          fails by ~1e6\n", r2)
	fmt.Printf("                                     extra 2S shift %.1f GHz (%.1f GHz\n", ContactShiftHz(r2-DarwinR2, 2)/1e9, ContactShiftHz(r2, 2)/1e9)
	fmt.Println("                                     if added to Dirac's)")
	q := LEPq()
	fmt.Printf("  form factor at LEP2 (√s = 200 GeV) |F| ≤ %.1e (shell: sin x / x,       ≥ %.3f                          fails\n", 1/(q*t.ReLbar), LEPFmin)
	fmt.Printf("                                     x = q r_e = %.1e)\n", q*t.ReLbar)
	fmt.Printf("  energy freed per annihilating pair %.2f × 2mc² (eq. 2.14)             2mc² (two 511 keV photons)        fails by %.0f %%\n", t.EnergyMC2, 100*(1-t.EnergyMC2))
	fmt.Printf("  annihilation light's wavelength    r_e = %.3f × h/mc                  h/mc = 2.426e-12 m (511 keV)      fails by ×%.1f\n", t.WaveOverL, 1/t.WaveOverL)
	fmt.Println("  annihilation light's shape         one spherical shell, no frequency  two back-to-back 511 keV photons, fails")
	fmt.Println("                                                                       polarisations crossed (Wu &")
	fmt.Println("                                                                       Shaknov 1950)")
	fmt.Println("  What survives: his field-free void, the charge on its surface, the positron as opposite rotation and")
	fmt.Println("  annihilation as the voids filling in are pictures the void model can keep and test. His formulas for")
	fmt.Println("  them fail. The tests below use our ether (void_shells), where the void and its light are computed.")
}

func part2(s *State) {
	ee, el := s.Energies()
	M := ee + el
	lbar := EpsRef / M
	fmt.Println()
	fmt.Println("Part 2. Charge radius: where can the charge sit on the computed bubble (ε = 300, one photon)?")
	fmt.Printf("  bubble: M = %.4f, ƛ = ε/M = %.4f ξ; void reaches %.1f ξ = %.2f ƛ across\n", M, lbar, s.Measure().REq, s.Measure().REq/lbar)
	fmt.Println("  Hydrogen: any charge spread shifts nS levels by (2/3)α⁴mc²⟨r²⟩/(n³ƛ²). Dirac's Darwin term is exactly such")
	fmt.Println("  a spread, ⟨r²⟩ = ¾ƛ², and it is confirmed; the 2S Lamb shift (1057.845 ± 0.009 MHz) leaves room for")
	fmt.Printf("  |Δ⟨r²⟩| ≤ %.1e ƛ² more. LEP2 needs rms ≤ %.1e ƛ (%.1e m).\n", LambBound/ContactShiftHz(1, 2), LEPrmsBound(), LEPrmsBound()*LbarM)
	fmt.Println("  placement                          rms radius   ⟨r²⟩/ƛ²   extra 2S shift if it replaces   if added to   F at q = 1/ƛ")
	fmt.Println("                                     (ƛ)                    Darwin (GHz)                    Dirac (GHz)")
	rs := []float64{}
	for r := 0.05; r <= 12.0001; r += 0.05 {
		rs = append(rs, r*lbar)
	}
	encl := map[Placement][]float64{}
	for _, p := range []Placement{Deficit, Wall, Light} {
		rho := s.Charge(p, 0)
		r2 := s.R2(rho) / (lbar * lbar)
		encl[p] = s.Enclosed(rho, rs)
		fmt.Printf("  %-34s %-12.3f %-9.2f %-31.1f %-13.1f %.3f\n", p, math.Sqrt(r2), r2, ContactShiftHz(r2-DarwinR2, 2)/1e9, ContactShiftHz(r2, 2)/1e9, s.FormFactor(rho, 1/lbar))
	}
	fmt.Println("  Every placement is 3–7 ƛ across, so it shifts hydrogen's 2S by hundreds of GHz against a ±9 kHz margin,")
	fmt.Println("  and its form factor is already far below 1 at q = 1/ƛ (about 0.5 MeV), let alone at LEP's 200 GeV.")
	fmt.Println()
	fmt.Println("  The charge a probe sees inside radius r (Gauss's law), against QED's running charge:")
	fmt.Println("  r / ƛ    void: missing ether   wall     light     QED (Uehling): α(r)/α")
	for _, want := range []float64{8, 4, 2, 1, 0.5, 0.25, 0.1} {
		n := int(math.Round(want/0.05)) - 1
		fmt.Printf("  %-8.2f %-21.4f %-8.4f %-9.4f %.5f\n", want, encl[Deficit][n], encl[Wall][n], encl[Light][n], UehlingRatio(want))
	}
	for n, r := range rs {
		csvAdd("enclosed", "r_over_lbar,deficit,wall,light,qed", fmt.Sprintf("%.3f,%.6f,%.6f,%.6f,%.6f", r/lbar, encl[Deficit][n], encl[Wall][n], encl[Light][n], UehlingRatio(r/lbar)))
	}
	fmt.Println("  Opposite directions: inside the void the seen charge falls toward zero; in QED it rises above e (measured:")
	fmt.Println("  α = 1/128 at the Z mass). A charged void must be smaller than ~1e-6 ƛ to look like the measured electron.")
}

func part3(ball *State, frameDir string) {
	fmt.Println()
	fmt.Println("Part 3. Exact charge: if the charge is the missing ether, is it fixed and does it come in units?")
	fmt.Println("  Measured: every electron, muon and proton carries the same |e| to ~1e-21 (atoms are neutral) whatever")
	fmt.Println("  its mass, speed or surroundings.")
	_, t1 := ball.LocalDeficit(1e9)
	fmt.Println("  state                 missing ether (ξ³)   / one photon at ε = 300")
	for _, c := range []struct {
		name string
		eps  float64
		occ  []int
		r0   float64
	}{{"ε = 100", 100, []int{0}, 3.5}, {"ε = 300", 300, []int{0}, 5}, {"ε = 1000", 1000, []int{0}, 7}, {"2 photons, ε = 300", 300, []int{0, 1}, 6.4}, {"3 photons, ε = 300", 300, []int{0, 1, 2}, 7.6}} {
		s := ball
		if c.name != "ε = 300" {
			s = relaxedBall(c.eps, c.occ, c.r0)
		}
		_, tot := s.LocalDeficit(1e9)
		fmt.Printf("  %-21s %-20.3f %.4f\n", c.name, tot, tot/t1)
	}
	fmt.Println("  Not in units: two photons do not make twice the missing ether, and the amount tracks the light's strength.")
	fmt.Println()
	fmt.Println("  In real time, a bubble started 3 % too large across breathes (no friction, ε = 300); the missing")
	fmt.Println("  ether within 7 ξ, against the relaxed bubble's:")
	f := ball.Profile()
	ref, _ := ball.LocalDeficit(7)
	s := NewState(grid(), EpsRef, V0, []int{0}, Scaled(f, 1.03, 1.03))
	e0, l0 := s.Energies()
	_, tot0 := s.LocalDeficit(1e9)
	dt, T := 0.01, 60.0
	lo, hi := math.Inf(1), math.Inf(-1)
	var series []float64
	fr := openFrames(frameDir, "rt_breathe.bin", s)
	n := int(math.Round(T / dt))
	fmt.Println("    t     missing ether within 7 ξ   whole box     energy drift")
	for k := 0; k <= n; k++ {
		if k%10 == 0 {
			loc, _ := s.LocalDeficit(7)
			series = append(series, loc/ref)
			if k >= 500 { // after t = 5, once the start has passed through
				lo, hi = math.Min(lo, loc/ref), math.Max(hi, loc/ref)
			}
			csvAdd("breathe", "t,local_deficit", fmt.Sprintf("%.2f,%.6f", float64(k)*dt, loc/ref))
		}
		if fr != nil && k%50 == 0 {
			writeFrame(fr, float64(k)*dt, s)
		}
		if k%1000 == 0 {
			loc, tot := s.LocalDeficit(7)
			ee, el := s.Energies()
			fmt.Printf("    %-5.0f %-26.4f %-13.6f %+.2e\n", float64(k)*dt, loc/ref, tot/tot0, (ee+el-e0-l0)/(e0+l0))
		}
		if k < n {
			s.Step(dt)
		}
	}
	closeFrames(fr)
	fmt.Printf("  A 3 %% disturbance makes the missing ether near the particle swing from %.3f to %.3f of its rest value.\n", lo, hi)
	fmt.Println("  The whole box keeps it exactly (the ether is conserved), but near the particle it flickers as the void")
	fmt.Println("  breathes and sheds sound. A charge carried by missing ether would change by percent with every jolt;")
	fmt.Println("  measured, the charge holds to 1e-21. It would also differ between bubbles of different light (above).")
	part5breath = FirstTurn(series, 0.1, 21)
	fmt.Println()
	fmt.Println("  What does come in units: the light's winding. Each photon circulates once around the axis in two turns")
	fmt.Println("  (m = ±½); m cannot change continuously, and its sign is the handedness — Tewari's opposite rotation for")
	fmt.Println("  the positron. If the charge is the winding (e for m = −½, −e for m = +½ by convention), it is exact and")
	fmt.Println("  conserved. But here m enters the light only as m² = ¼: nothing in these equations feels its sign, so the")
	fmt.Println("  winding labels the charge without producing any electric force (Part 4).")
}

func part4(ball *State, frameDir string) {
	ee1, el1 := ball.Energies()
	M1 := ee1 + el1
	lbar := EpsRef / M1
	fmt.Println()
	fmt.Println("Part 4. Electron and positron: two voids on the axis, centres 2d apart (relaxed bubbles, light solved)")
	fmt.Println("  e+e−: opposite windings may share the lowest light mode; e−e−: same winding, one photon per mode.")
	fmt.Println("  Coulomb for comparison: ∓α ε / 2d (ε = ħc). Light: photon energy ε√λ against the escape threshold ε√V0.")
	fmt.Printf("  one bubble: M = %.4f, photon %.2f, escape needs %.2f (λ = %.4f against V0 = %.0f)\n", M1, EpsRef*math.Sqrt(ball.Lam[0]), EpsRef*math.Sqrt(V0), ball.Lam[0], V0)
	fmt.Println("  2d (ξ)   2d (ƛ)   wall gap   e+e− energy − 2M   e−e− energy − 2M   e+e− minus e−e−   Coulomb αε/2d")
	f := ball.Profile()
	for _, d := range []float64{4.4, 4.8, 5.2, 5.6, 6.0, 6.8, 8.0, 10.0} {
		s := NewState(grid(), EpsRef, V0, []int{0, 1}, Pair(f, d))
		s.refreshLight(40)
		ee, _ := s.Energies()
		epm := ee + 2*EpsRef*math.Sqrt(s.Lam[0]) - 2*M1
		emm := ee + EpsRef*(math.Sqrt(s.Lam[0])+math.Sqrt(s.Lam[1])) - 2*M1
		coul := Alpha * EpsRef / (2 * d)
		gap := 2*d - 2*ball.Measure().RAxis
		fmt.Printf("  %-8.1f %-8.2f %-10.1f %-18.5f %-18.5f %-17.2e %.5f\n", 2*d, 2*d/lbar, gap, epm, emm, epm-emm, coul)
		csvAdd("pair", "D,D_over_lbar,gap,e_pm,e_mm,coulomb", fmt.Sprintf("%.2f,%.4f,%.2f,%.8f,%.8f,%.8f", 2*d, 2*d/lbar, gap, epm, emm, coul))
	}
	fmt.Println("  Both pairs feel the same short-range pull once their walls overlap, and it dies within a healing length")
	fmt.Println("  (the ether's static response is screened: perturbations fall as e^{−√2 r}/r). The only difference")
	fmt.Println("  between e+e− and e−e− is the light's exchange splitting, also exponential. No 1/r force, no sign:")
	fmt.Println("  like charges do not repel. Coulomb's law needs a field that reaches out, which this ether lacks.")
	fmt.Println()
	fmt.Println("  In real time (no friction), two bubbles started with walls 1.6 ξ apart (2d = 8.8 ξ):")
	fmt.Println("  pair    t     voids   centre density   light energy   light in voids   ether energy beyond 14 ξ   drift")
	for _, c := range []struct {
		tag string
		occ []int
	}{{"e+e−", []int{0, 0}}, {"e−e−", []int{0, 1}}} {
		s := NewState(grid(), EpsRef, V0, c.occ, Pair(f, 4.4))
		s.refreshLight(40)
		e0, l0 := s.Energies()
		name := "rt_pm.bin"
		if c.tag == "e−e−" {
			name = "rt_mm.bin"
		}
		fr := openFrames(frameDir, name, s)
		dt, T := 0.01, 40.0
		n := int(math.Round(T / dt))
		for k := 0; k <= n; k++ {
			if fr != nil && k%50 == 0 {
				writeFrame(fr, float64(k)*dt, s)
			}
			if k%1000 == 0 {
				ee, el := s.Energies()
				l := s.Measure()
				fmt.Printf("  %-7s %-5.0f %-7d %-16.3f %-14.3f %-16.4f %-26.4f %+.1e\n", c.tag, float64(k)*dt, l.Voids, l.CoreDens, el, l.LightIn, s.EtherOutside(14, 0), (ee+el-e0-l0)/(e0+l0))
				csvAdd("collide", "pair,t,voids,core,light_energy,light_in,ether_out", fmt.Sprintf("%s,%.1f,%d,%.5f,%.5f,%.5f,%.5f", c.tag, float64(k)*dt, l.Voids, l.CoreDens, el, l.LightIn, s.EtherOutside(14, 0)))
			}
			if k < n {
				s.Step(dt)
			}
		}
		closeFrames(fr)
	}
	two := relaxedBall(EpsRef, []int{0, 0}, 6.4)
	e2, l2 := two.Energies()
	like := relaxedBall(EpsRef, []int{0, 1}, 6.4)
	e3, l3 := like.Energies()
	fmt.Printf("  Relaxed, the merged e+e− bubble has M = %.4f against 2M = %.4f (e−e−: %.4f): %.1f %% (%.1f %%) is freed,\n", e2+l2, 2*M1, e3+l3, 100*(1-(e2+l2)/(2*M1)), 100*(1-(e3+l3)/(2*M1)))
	fmt.Println("  and it leaves as")
	fmt.Println("  sound in the ether, not light. How they meet differs (see the frames): e+e− may share one light mode, so")
	fmt.Println("  within t ≈ 5 both photons gather in one void; the emptied void collapses and sends out a spherical shell")
	fmt.Println("  of sound — Tewari's collapsing void and outgoing shell, but the shell is sound and the light stays.")
	fmt.Println("  e−e− cannot share a mode, so each void keeps its photon until the two coalesce, by t ≈ 20.")
	fmt.Println("  Both photons stay trapped (their energy is ~1/5 of what escape needs), so")
	fmt.Println("  the pair never annihilates: it becomes one neutral (windings cancel) bubble of light that lasts forever.")
	fmt.Println("  Measured: positronium decays in 0.125 ns (para) or 142 ns (ortho) into 2 or 3 photons carrying all 2mc².")
	fmt.Println("  Like pairs (e−e−) merge just the same. Tewari's picture — the voids fill in and the light leaves — needs")
	fmt.Println("  the light to escape when the windings cancel; in this ether nothing releases it.")
}

var part5breath float64 // time of the breathing's first swing back

func part5(ball *State, frameDir string) {
	ee, el := ball.Energies()
	M := ee + el
	fmt.Println()
	fmt.Println("Part 5. Stiffness: the void's own vibrations (ħ = 1 in these units, so ħω / mc² = ω / M)")
	f := ball.Profile()
	s := NewState(grid(), EpsRef, V0, []int{0}, Scaled(f, 1/math.Sqrt(1.03), 1.03))
	q0 := ball.Quadrupole(RLocal)
	var q []float64
	dt, T := 0.01, 60.0
	n := int(math.Round(T / dt))
	for k := 0; k <= n; k++ {
		if k%10 == 0 {
			q = append(q, s.Quadrupole(RLocal)-q0)
			csvAdd("quad", "t,quadrupole", fmt.Sprintf("%.2f,%.6f", float64(k)*dt, q[len(q)-1]))
		}
		if k < n {
			s.Step(dt)
		}
	}
	quadW, quadN := RingFrequency(q, 0.1, 30, 21)
	bw := math.Pi / part5breath
	fmt.Println("  Read from the motion after a 3 % disturbance (signals smoothed over 2 time units to remove the sound's")
	fmt.Println("  ripple). Both modes lose their energy to sound within a few swings: broad resonances, good to ~20 %.")
	fmt.Println("  mode                       ω (μ/ħ)        ħω / mc²     energy      how read")
	fmt.Printf("  breathing (l = 0)          %-14.2f %-12.1e %-11s does not ring: one slow swing back, first turn at t = %.1f (ω ≈ π/t)\n", bw, bw/M, fmt.Sprintf("%.2f keV", bw/M*MC2keV), part5breath)
	fmt.Printf("  stretching (l = 2)         %-14.2f %-12.1e %-11s rings: median half-period over %d turning points\n", quadW, quadW/M, fmt.Sprintf("%.2f keV", quadW/M*MC2keV), quadN)
	fmt.Println("  sliding (l = 1)            0 (it is a free translation: no restoring force, no dipole)")
	fmt.Println("  A round void has no electric dipole moment (measured |d_e| < 4.1e-30 e·cm: passes), but it does have")
	fmt.Println("  soft internal modes below a keV: an excited electron. They sit near ħ × (sound speed / void size),")
	fmt.Println("  and the ether's sound (1/√2) is ~400 times slower than its light (ε = 300), so the void rings at about")
	fmt.Println("  1/1000 of the Compton frequency. Measured: no excited state of the electron has ever been seen — light")
	fmt.Println("  scatters off free electrons as off a structureless charge (Thomson, Klein–Nishina), and collider")
	fmt.Println("  searches exclude excited electrons up to TeV masses.")
}

func part6() {
	fmt.Println()
	fmt.Println("Part 6. Scorecard: electrical properties assigned to the void")
	fmt.Println("  property                    Tewari's void (Part 1)        computed bubble (Parts 2–5)      measured")
	fmt.Println("  charge −e                   input (fixes r_e)             input (assigned)                 suggestion: put in by hand")
	fmt.Println("  field conservative          fails (∮E·dl ≠ 0)             no field at all                  yes")
	fmt.Println("  field isotropic             fails (sin²θ; Lyman-α split)  —                                yes (spin-independent Coulomb)")
	fmt.Println("  charge radius / hydrogen    fails (1.04 ƛ shell)          fails (3–7 ƛ)                    0.866 ƛ Darwin, extra ≤ 6e-4 ƛ")
	fmt.Println("  form factor at LEP2         fails                         fails                            pointlike to ~3e-19 m")
	fmt.Println("  running charge              (none; Coulomb stops at r_e)  opposite sign (falls inside)     α rises at short distance")
	fmt.Println("  charge quantized            claims (no smaller void)      winding yes; missing ether no    exact to 1e-21")
	fmt.Println("  charge constant in time     —                             missing ether flickers by %      yes")
	fmt.Println("  like charges repel          via flow (not computed)       fails (merge)                    yes, 1/r²")
	fmt.Println("  e+e− → 2mc² of light        fails (0.8, one shell)        fails (sound only, light kept)   2 × 511 keV photons")
	fmt.Println("  EDM                         0 (axial)                     0 (round)                        < 4.1e-30 e·cm: passes")
	fmt.Println("  no internal excitations     —                             fails (keV modes)                none seen")
	fmt.Println("  Verdict: assigning electrical properties to the void makes the model testable, and the tests are decided.")
	fmt.Println("  The void is too big for hydrogen and LEP by six orders, and the ether as built has no long-range field to")
	fmt.Println("  carry a charge's force or to let light out when voids meet. The one property that survives is the light's")
	fmt.Println("  winding as an exact, signed charge label — which still needs a field that feels its sign.")
}

func openFrames(dir, name string, s *State) *os.File {
	if dir == "" {
		return nil
	}
	f, err := os.Create(filepath.Join(dir, name))
	if err != nil {
		log.Fatal(err)
	}
	binary.Write(f, binary.LittleEndian, []int32{int32(s.G.NR), int32(s.G.NZ)})
	return f
}

func closeFrames(f *os.File) {
	if f != nil {
		f.Close()
	}
}

func writeFrame(f *os.File, t float64, s *State) {
	s.updateDens()
	binary.Write(f, binary.LittleEndian, float32(t))
	d := make([]float32, len(s.dens))
	x := make([]float32, len(s.dens))
	for i := range d {
		d[i] = float32(s.dens[i])
		for _, q := range s.Occ {
			x[i] += float32(s.Chi[q][i] * s.Chi[q][i])
		}
	}
	binary.Write(f, binary.LittleEndian, d)
	binary.Write(f, binary.LittleEndian, x)
}
