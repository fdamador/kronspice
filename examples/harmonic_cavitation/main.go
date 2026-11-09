// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command harmonic_cavitation designs an efficient cavitation generator: two
// dimpled metal drums, one at each outlet, spinning on one shaft (driven by
// a motor or a Tesla turbine) inside a fixed housing, sized from the
// published rotational cavitation reactors (CLAUDE.md, LITERATURE.md). The
// measure is cavitation per joule: collapse energy over shaft work.
//
//	go run ./examples/harmonic_cavitation          # about 1–2 minutes
//	go run ./examples/harmonic_cavitation -part 4  # one part
//	go test ./examples/harmonic_cavitation
//
// Parts:
//  1. calibration on Xue et al. 2024;
//  2. sizing a drum: speeds for inception and the best σ, dimples per row;
//  3. the reduced cavitation model on reference drums, with the
//     kronspice shaft circuit as its ledger twin;
//  4. NSGA-II over the dimple drums, and the same search under other
//     values of the declared constants;
//  5. the gap as an acoustic ring in kronspice: azimuthal modes and the
//     sonic match of the dimple pattern;
//  6. one cavity's collapse (Rayleigh–Plesset, ledger) and many collapses
//     gathering on the axis;
//  7. the energy books, socket to radiator, motor against pump + turbine;
//  8. the Mandelbulb challenger (H46, proposed): geometry checks, n = 8 on
//     the winning machine, NSGA-II over n, best efficiency for each n, and
//     the bulb without its fine detail;
//  9. summary against the controls;
//  10. the closed steam loop (user's layout): water in at the drums' outer
//     top, steam from the drums drives the Tesla turbine on the same shaft,
//     condensed at the outlet (the home's heat) and pumped back; boiling in
//     the drums against subcooled drums with a flash tank;
//  11. a 3,000 ft² home at −10 °F, 24/7: rows, drum speed matched to the
//     demand, condenser, emitters, the winter, hole size against the
//     Tesla turbine's speed.
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"sort"
	"time"
)

func main() {
	part := flag.Int("part", 0, "run one part (1–11); 0 runs all")
	seed := flag.Uint64("seed", 1, "GA seed")
	viz := flag.String("viz", "", "write the report's data (CSV) to this directory and exit")
	flag.Parse()
	run := func(k int) bool { return *part == 0 || *part == k }
	t0 := time.Now()
	p := DefaultParams()
	if *viz != "" {
		if err := WriteViz(*viz, p, *seed); err != nil {
			fmt.Println(err)
			os.Exit(1)
		}
		return
	}

	if run(1) || run(2) {
		stage1(p.Cal, run(1), run(2))
	}

	ref := DimpleMachine(0.05, 10, 0.010, 0.5, 1.5, 2.2, 1.5e-3, 3050, pAtm)
	if run(3) {
		part3(p, ref)
	}

	var best Scored
	var front []Scored
	if run(4) || run(5) || run(6) || run(7) || run(8) || run(9) || run(10) || run(11) {
		fmt.Println("4. NSGA-II over the dimple drums (population 80, 80 generations), shaft ≤ 3 kW")
		var err error
		front, err = Evolve(p, dimpleSpace, DimpleFrom, dimpleLand, 80, 80, *seed, nil)
		if err != nil || len(front) == 0 {
			fmt.Println("   GA failed:", err)
			os.Exit(1)
		}
		best = Heater(front)
		if run(4) {
			part4(p, front, *seed)
		}
	}
	if run(5) {
		part5(best)
	}
	if run(6) {
		part6(p, best)
	}
	if run(7) {
		part7(best)
	}
	var bulbBest Scored
	if run(8) || run(9) {
		bulbBest = part8(p, best, *seed, run(8))
	}
	if run(9) {
		part9(p, best, bulbBest)
	}
	if run(10) {
		part10(p, best)
	}
	if run(11) {
		part11(p, best)
	}
	fmt.Printf("\n(run time %.0f s)\n", time.Since(t0).Seconds())
}

func stage1(cal Calibration, one, two bool) {
	x := Xue2024
	if one {
		fmt.Println("1. Calibration on", x.Name, "(rotor Ø 80 mm, 20 grooves 6 mm wide; 1 atm at the rotor declared)")
		fmt.Printf("   surface speed: %.2f m/s none, %.2f m/s cavitating, %.2f m/s best\n",
			omega(x.RPMNone)*x.R, omega(x.RPMOn)*x.R, omega(x.RPMBest)*x.R)
		fmt.Printf("   inception σ between %.3f and %.3f (σ_i = %.3f taken); best σ = %.3f\n",
			cal.SigmaOn, cal.SigmaNone, 0.5*(cal.SigmaOn+cal.SigmaNone), cal.SigmaBest)
		fmt.Printf("   pulsation St = fW/U = %.3f; dominant f / groove passing = %.3f\n", cal.StLock, cal.PassRatio)
		fmt.Printf("   resonance rule on Xue's rotor: N = %.1f (built 20)\n", 2*math.Pi*x.R*cal.StLock/x.W)
		p := DefaultParams()
		fmt.Printf("   vapour law g(σ) = x·e^(1−x), x = (σ_i − σ)/Δ; Δ = %.3f fitted so Xue's collapse rate peaks at 3800 rpm\n\n", p.Delta)
	}
	if !two {
		return
	}
	fmt.Println("2. Sizing a drum from the calibration (outlet 1 atm, water at 300 K)")
	si := 0.5 * (cal.SigmaOn + cal.SigmaNone)
	for _, r := range []float64{0.04, 0.05, 0.084} {
		ui := math.Sqrt(2 * (pAtm - pVap) / (rho * si))
		ub := math.Sqrt(2 * (pAtm - pVap) / (rho * cal.SigmaBest))
		n := 2 * math.Pi * r * cal.StLock / 0.010
		fmt.Printf("   drum R %.0f mm: cavitation from %.0f rpm (U %.1f m/s), strongest pulsation at %.0f rpm (U %.1f m/s); resonance rule for Ø 10 mm dimples: %.0f per row\n",
			r*1e3, ui/(2*math.Pi*r)*60, ui, ub/(2*math.Pi*r)*60, ub, n)
	}
	fmt.Println("   every dimple on a drum runs at the same speed and σ, so one speed suits them all")
	fmt.Println()
}

func printEval(e Eval) {
	fmt.Printf("   U %.2f m/s, p %.1f kPa, σ %.3f, g(σ) %.2f, %d cavities\n", e.U, e.P/1e3, e.Sigma, e.G, e.Cavities)
	fmt.Printf("   shaft %.0f W = friction %.0f W + cavity drag %.0f W; collapse energy %.0f W of the drag → η = %.4f\n",
		e.PShaft, e.PFriction, e.PDrag, e.PCollapse, e.Eta)
}

func part3(p Params, m Machine) {
	fmt.Println("3. Reference drums: R 50 mm, 10 rows of Ø 10 mm hemispheres, rotor pitch 1.5 W, housing pitch 2.2 W, gap 1.5 mm, 3050 rpm")
	e, err := p.Evaluate(m)
	if err != nil {
		fmt.Println("   ", err)
		return
	}
	printEval(e)
	pin, res, err := LedgerTwin(e, m.RPM)
	fmt.Printf("   kronspice shaft circuit: source gives %.2f W against the model's %.2f W, ledger residual %.1e (err %v)\n", pin, e.PShaft, res, err)
	m.Housing = nil
	e2, _ := p.Evaluate(m)
	fmt.Printf("   same with a plain housing: η = %.4f (collapse %.0f W of %.0f W)\n\n", e2.Eta, e2.PCollapse, e2.PShaft)
}

func describe(s Scored) string {
	x := s.X
	h := "plain housing"
	if math.Round(x[5]) == 2 {
		h = fmt.Sprintf("housing pitch %.2f W", x[6])
	}
	nr := 0
	if len(s.M.Rotor) > 0 {
		nr = s.M.Rotor[0][0].Count
	}
	return fmt.Sprintf("R %.1f mm, %d row(s), Ø %.1f mm, depth %.2f W, %d per row (pitch %.2f W), %s, gap %.2f mm, %.0f rpm, outlet %.2f bar",
		1e3*x[0], int(math.Round(x[1])), 1e3*x[2], x[3], nr, x[4], h, 1e3*x[7], x[8], x[9]/1e5)
}

func part4(p Params, front []Scored, seed uint64) {
	fmt.Printf("   %d feasible designs on the front; by efficiency:\n", len(front))
	fmt.Println("        η     collapse W   shaft W     σ")
	for i, s := range front {
		if i >= 8 {
			break
		}
		fmt.Printf("   %7.4f   %8.0f   %8.0f   %.3f   %s\n", s.E.Eta, s.E.PCollapse, s.E.PShaft, s.E.Sigma, describe(s))
	}
	best := Heater(front)
	fmt.Printf("   the front is flat in η (%.4f–%.4f), so the drums can be sized to the load; heater design (most collapse power within 0.5 %% of the best η):\n   %s\n",
		front[len(front)-1].E.Eta, front[0].E.Eta, describe(best))
	pin, res, _ := LedgerTwin(best.E, best.M.RPM)
	fmt.Printf("   heater design: kronspice shaft circuit %.2f W (model %.2f W), ledger residual %.1e\n", pin, best.E.PShaft, res)

	fmt.Println("   The same search with the declared constants changed (population 60, 60 generations):")
	for _, r := range Sensitivity(p, seed) {
		if r.Err != nil {
			fmt.Printf("   %-20s failed: %v\n", r.Name, r.Err)
			continue
		}
		fmt.Printf("   %-20s η %.4f  %s\n", r.Name, r.Best.E.Eta, describe(r.Best))
	}
	fmt.Println()
}

func part5(best Scored) {
	m := best.M
	e := best.E
	nr := m.Rotor[0][0].Count
	bandLen := float64(m.Rows) * m.RowLen
	fmt.Printf("5. The gap as an acoustic ring (kronspice AC): best design, R %.1f mm, gap %.2f mm, band %.1f mm\n", 1e3*m.R, 1e3*m.Gap, 1e3*bandLen)
	sections := 8 * nr
	if sections < 64 {
		sections = 64
	}
	pk, ex, err := RingModes(m.R, m.Gap, bandLen, cLiq, sections, 3)
	if err == nil && len(pk) >= 3 {
		fmt.Printf("   clear water: modes 1–3 at %.0f, %.0f, %.0f Hz (lattice exact %.0f, %.0f, %.0f; m·c/2πR = %.0f, %.0f, %.0f)\n",
			pk[0], pk[1], pk[2], ex[0], ex[1], ex[2], cLiq/(2*math.Pi*m.R), 2*cLiq/(2*math.Pi*m.R), 3*cLiq/(2*math.Pi*m.R))
	}
	fpass := float64(nr) * m.RPM / 60
	fmt.Printf("   rotor pattern: %d cavities per ring → passing %.0f Hz; mode %d in clear water is at %.0f Hz (×%.0f)\n",
		nr, fpass, nr, float64(nr)*cLiq/(2*math.Pi*m.R), float64(nr)*cLiq/(2*math.Pi*m.R)/fpass)
	a := SonicVoid(e.U, e.P)
	if math.IsNaN(a) {
		fmt.Printf("   the pattern meets mode %d only if the mixture's sound speed falls to U = %.1f m/s; Wood's minimum at %.0f kPa is %.1f m/s, so it never does here\n",
			nr, e.U, e.P/1e3, WoodSpeed(0.5, e.P))
	} else {
		fmt.Printf("   it meets mode %d when the bubbly mixture's sound speed falls to U = %.1f m/s: void fraction α = %.3f (Wood)\n", nr, e.U, a)
		pk2, ex2, err := RingModes(m.R, m.Gap, bandLen, e.U, 8*nr, nr)
		if err == nil && len(pk2) >= nr {
			fmt.Printf("   kronspice ring at c = U, %d sections: mode %d at %.1f Hz (lattice exact %.1f) against passing %.1f Hz\n", 8*nr, nr, pk2[nr-1], ex2[nr-1], fpass)
		}
	}
	fmt.Println()
}

func part6(p Params, best Scored) {
	m, e := best.M, best.E
	c0 := m.Rotor[0][0]
	pulse := 0.0
	if len(m.Housing) > 0 {
		pulse = p.Cpulse * 0.5 * rho * e.U * e.U * math.Pow(p.G0/m.Gap, p.GapExp)
	}
	dp := e.P - pVap + pulse
	vc := p.KappaV * e.G * c0.V
	r0 := math.Cbrt(3 * vc / (4 * math.Pi))
	fmt.Printf("6. One cavity's collapse: vapour %.1f mm³ (sphere R0 %.2f mm) under %.1f kPa (local %.1f + housing pulse %.1f)\n",
		vc*1e9, r0*1e3, dp/1e3, (e.P-pVap)/1e3, pulse/1e3)
	cv := RunCollapse(r0, dp, 0, 0, 0)
	fmt.Printf("   empty cavity, no viscosity: collapse %.2f µs against Rayleigh 0.915·R0·√(ρ/Δp) = %.2f µs\n", cv.T*1e6, 0.915*r0*math.Sqrt(rho/dp)*1e6)
	c := RunCollapse(r0, dp, 1e-3, 0.0717, p.Mu)
	fmt.Printf("   with gas (1e-3 of Δp at R0, declared), surface tension, viscosity: collapse %.2f µs to R_min %.1f µm, gas peak %.2g Pa\n",
		c.T*1e6, c.Rmin*1e6, c.PGasMax)
	fmt.Println("   (incompressible Rayleigh–Plesset: the final peaks are overstated; liquid compressibility, as in Keller–Miksis, limits them)")
	fmt.Printf("   energy: (p∞ − p_v)·V0 = %.3g J; peak liquid kinetic energy %.3g J; ledger residual %.1e\n", c.EPot, c.KMax, c.LedgerRelResidual)
	ns := 0
	if len(m.Housing) > 0 {
		ns = m.Housing[0].Count
	}
	nr := c0.Count
	k := Simultaneous(nr, ns)
	p1 := c.PRadPeakAt1m / m.R
	fmt.Printf("   radiated peak at the axis (distance R = %.1f mm) from one collapse: %.3g Pa\n", 1e3*m.R, p1)
	fmt.Printf("   rotor %d and housing %d cavities per ring line up %d at a time → %d collapses together: focal peak ≤ %.3g Pa\n", nr, ns, k, k, Focus(c, k, m.R))
	if ns > 0 {
		fmt.Printf("   with housing count = rotor count (%d), all %d collapse together: focal peak ≤ %.3g Pa (×%d), the energy still only the %d collapses'\n",
			nr, nr, Focus(c, nr, m.R), nr/max(k, 1), nr)
	}
	fmt.Println("   the focal values are upper bounds: linear superposition of incompressible near fields, at amplitudes where neither holds")
	fmt.Println("   the axis of a solid drum is metal: the focused pulse goes into the rotor, not the water, unless the drum is hollow with a water core")
	fmt.Println()
}

func part7(best Scored) {
	fmt.Println("7. Energy books for the heater design (declared: motor 0.90, pump 0.70, Tesla turbine 0.30)")
	fmt.Println("   drive                    electric W  shaft W  heat to water W  heat to room W  collapse W  collapse/J   steam kg/h  L/min +45 K  heat pump (COP 3.5) W")
	for _, d := range Drives(best.E) {
		fmt.Printf("   %-24s %9.0f  %7.0f  %15.0f  %14.0f  %10.0f  %10.4f  %10.2f  %11.2f  %10.0f\n",
			d.Name, d.PElec, d.PShaft, d.PWaterHeat, d.PRoomHeat, d.PCollapse, d.CavPerJoule, d.SteamKgH, d.WaterLmin45K, d.HeatPumpW)
	}
	fmt.Println("   all heat ≤ the electricity; the turbine drive costs cavitation per joule but its losses still warm the water")
	fmt.Println()
}

func bulbDescribe(s Scored) string {
	x := s.X
	h := "plain housing"
	if math.Round(x[3]) == 2 {
		h = fmt.Sprintf("housing pitch %.2f × 10 mm", x[4])
	}
	return fmt.Sprintf("n %d, R %.1f mm, %d ring(s), %s, gap %.2f mm, %.0f rpm, outlet %.2f bar",
		int(math.Round(x[0])), 1e3*x[1], int(math.Round(x[2])), h, 1e3*x[5], x[6], x[7]/1e5)
}

func part8(p Params, best Scored, seed uint64, show bool) Scored {
	say := func(f string, a ...any) {
		if show {
			fmt.Printf(f, a...)
		}
	}
	say("8. The Mandelbulb challenger (H46, proposed): polar convention, MaxIt %d, detail ≥ %.1f mm, rings within ±%.1f rad\n", bulbMaxIt, bulbRes*1e3, bulbLatMax)
	say("   tips round the equator within 1 %% of the outer radius (n − 1 main lobes, plus smaller bulbs reaching as far), polar / latitude:\n  ")
	for _, n := range []float64{3, 4, 6, 8, 12, 16} {
		pp := Bulb{N: n, MaxIt: bulbMaxIt}.Profile(0, 720)
		pl := Bulb{N: n, Latitude: true, MaxIt: bulbMaxIt}.Profile(0, 720)
		say(" n=%g: %d/%d", n, Lobes(pp, 0.01), Lobes(pl, 0.01))
	}
	say("\n")
	x := best.X
	ph := 0.0
	if math.Round(x[5]) == 2 {
		ph = x[6] * x[2] / 0.01
	}
	m8 := BulbMachine(Bulb{N: 8, MaxIt: bulbMaxIt}, x[0], 1, bulbLatMax, bulbRes, ph, x[7], x[8], x[9])
	if e8, err := p.Evaluate(m8); err == nil {
		say("   n = 8 bulb ring on the heater design (same R, gap, speed, pressure, housing): %d cavities, η = %.4f (dimples %.4f), shaft %.0f W\n",
			e8.Cavities, e8.Eta, best.E.Eta, e8.PShaft)
	}
	front, err := Evolve(p, bulbSpace, BulbFrom, bulbLand, 60, 60, seed, nil)
	if err != nil || len(front) == 0 {
		say("   bulb GA failed: %v\n", err)
		return Scored{}
	}
	bb := front[0]
	say("   NSGA-II over n and the machine (60 × 60): best η %.4f, collapse %.0f W, shaft %.0f W: %s\n", bb.E.Eta, bb.E.PCollapse, bb.E.PShaft, bulbDescribe(bb))
	if !show {
		return bb
	}
	say("   best η for each n (n fixed, NSGA-II 50 × 50 each):\n")
	rows := PerN(p, seed)
	for _, r := range rows {
		say("     n %2d  η %.4f  (%d cavities)\n", r.N, r.Eta, r.Cavities)
	}
	sort.Slice(rows, func(i, j int) bool { return rows[i].Eta > rows[j].Eta })
	say("   best n: %d\n", rows[0].N)
	// fine detail removed: MaxIt 3 against 10, same machine
	xb := bb.X
	phb := 0.0
	if math.Round(xb[3]) == 2 {
		phb = xb[4]
	}
	coarse := BulbMachine(Bulb{N: math.Round(xb[0]), MaxIt: 3}, xb[1], int(math.Round(xb[2])), bulbLatMax, bulbRes, phb, xb[5], xb[6], xb[7])
	if ec, err := p.Evaluate(coarse); err == nil {
		say("   the winning bulb with its fine detail removed (MaxIt 3): η %.4f against %.4f with detail (MaxIt %d)\n", ec.Eta, bb.E.Eta, bulbMaxIt)
	}
	say("\n")
	return bb
}

func part9(p Params, best, bulb Scored) {
	fmt.Println("9. Summary (cavitation per joule of shaft work; declared constants in CLAUDE.md, ranges in part 4)")
	fmt.Printf("   heater dimple drums: η %.4f, collapse %.0f W of %.0f W shaft — %s\n", best.E.Eta, best.E.PCollapse, best.E.PShaft, describe(best))
	m := best.M
	m.Housing = nil
	if e, err := p.Evaluate(m); err == nil {
		fmt.Printf("   control, same drums with a plain housing: η %.4f\n", e.Eta)
	}
	if bulb.E.PShaft > 0 {
		fmt.Printf("   best Mandelbulb drums: η %.4f — %s\n", bulb.E.Eta, bulbDescribe(bulb))
	}
	ref := DimpleMachine(0.05, 10, 0.010, 0.5, 1.5, 2.2, 1.5e-3, 3050, pAtm)
	if e, err := p.Evaluate(ref); err == nil {
		fmt.Printf("   reference drums of part 3: η %.4f at %.0f W shaft (over the 3 kW limit)\n", e.Eta, e.PShaft)
	}
}

// Loops returns the two closed-loop cycles for a design at 3 kW of shaft
// work (declared: Tesla turbine isentropic efficiency 0.30, condenser at
// 0.2 bar / 60 °C for home heating, motor 0.90, pump 0.70).
func Loops(p Params, best Scored) []LoopResult {
	m := best.M
	return []LoopResult{
		LoopA(p, m, sat150k, sat20k, 0.30, 3000),
		LoopB(p, m, 3e5, sat200k, sat101k, sat20k, 0.30, 3000),
	}
}

func part10(p Params, best Scored) {
	fmt.Println("10. Closed steam loop (user's layout): drums → steam → Tesla turbine on the shaft → condenser (home heat) → pump back")
	fmt.Println("   declared: turbine isentropic efficiency 0.30, condenser 0.2 bar (60 °C), motor 0.90, pump 0.70; 3 kW of shaft work in the drums")
	m := best.M
	fmt.Println("   η against water temperature, heater drums at the speed best for 60 °C water at each pressure:")
	for _, pd := range []float64{1.18e5, 1.5e5, 3e5} {
		_, rpm := pathEta(p, m, pd, 60, 60)
		fmt.Printf("     drum %.2f bar, %5.0f rpm: ", pd/1e5, rpm)
		mm := m
		mm.POut, mm.RPM = pd, rpm
		for _, t := range []float64{25, 60, 80, 100, 110, 120, 130} {
			q := p
			q.PVap = SatPressure(t + 273.15)
			if q.PVap >= pd {
				fmt.Printf(" %3.0f°C boils", t)
				continue
			}
			e, err := q.Evaluate(mm)
			if err != nil {
				fmt.Printf(" %3.0f°C  —   ", t)
				continue
			}
			fmt.Printf(" %3.0f°C %.3f", t, e.Eta)
		}
		fmt.Println()
	}
	carnot := 1 - (sat20k.T+273.15)/(sat150k.T+273.15)
	fmt.Printf("   Carnot limit between %.0f °C and %.0f °C: %.3f of the heat can become work\n", sat150k.T, sat20k.T, carnot)
	fmt.Println("   cycle                                         steam kg/h  turbine W  motor shaft W  electric W  home heat W  room W  collapse W  η loop  collapse/elec J  rpm")
	for _, r := range Loops(p, best) {
		mot, res, err := LoopCircuit(r, r.RPM)
		fmt.Printf("   %-45s %9.2f  %9.0f  %13.0f  %10.0f  %11.0f  %6.0f  %10.0f  %6.3f  %15.3f  %4.0f\n",
			r.Name, r.SteamKgH, r.PTurbine, r.PMotorShaft, r.PElec, r.PHome, r.PRoom, r.PCollapse, r.EtaLoop, r.CavPerElecJ, r.RPM)
		fmt.Printf("     books: electricity − (home + room) = %.1e W; kronspice shaft circuit: motor %.1f W, ledger residual %.1e (err %v)\n",
			r.Residual, mot, res, err)
		fmt.Printf("     sensible (subcooled) share of the heat %.2f, η on it %.3f\n", r.FracSensible, r.EtaSensible)
	}
	fmt.Println("   the loop cannot run itself: the turbine returns a few per cent of the shaft work; the motor supplies the rest,")
	fmt.Println("   and the home receives about the electricity, as with no turbine. Boiling in the drums wastes the cavitation:")
	fmt.Println("   steam that leaves has not collapsed, and p − p_v → 0 as the water nears boiling. Subcooled drums + flash keep it.")
}
