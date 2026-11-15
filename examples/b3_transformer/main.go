// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command b3_transformer: a transformer whose core is itself an N-turn
// spring (a helix closed by a return leg), so the magnetic flux goes round
// the spring's axis N times: the "B(3)" helical field of the user's design.
// A winding wrapped round the spring's wall links the flux once for every
// core pass it encloses, so the voltage ratio is
//
//	v_s / v_p = (N_s · n_s) / (N_p · n_p)
//
// with N the winding turns and n the core passes each turn encloses; the
// rate of change dΦ/dt is the same in every pass when no flux leaks across
// the gaps. The magnetic network (core segments, gap leakage, return leg)
// gives the inductances; kronspice runs the transformer with its ledger.
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

// Design windings: primary round all 10 passes at azimuth 0, secondary
// round 2 passes (turns 5 and 6) at 180°, 100 turns each.
func designWindings(n int) (pri, sec Winding) {
	pri = Winding{Name: "primary", Azimuth: 0, First: 0, Passes: n, Turns: 100}
	sec = Winding{Name: "secondary", Azimuth: math.Pi, First: n/2 - 1, Passes: 2, Turns: 100}
	if sec.First < 0 {
		sec.First = 0
	}
	return
}

func spread(v []float64) float64 {
	mn, mx := v[0], v[0]
	for _, x := range v {
		mn, mx = math.Min(mn, x), math.Max(mx, x)
	}
	return (mx - mn) / mx
}

func main() {
	data := flag.String("data", "", "also write the report's data (CSV) to this directory")
	flag.Parse()
	err := run()
	if err == nil && *data != "" {
		err = WriteData(*data)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "b3_transformer:", err)
		os.Exit(1)
	}
}

func run() error {
	fmt.Println("B(3) transformer — 10-turn spring core, primary round 10 core passes, secondary round 2")
	fmt.Println()

	// 1. Keeping flux off the gaps.
	fmt.Println("1. Gap leakage. Flux stays on the wire only if the air between neighbouring turns is a poor path next")
	fmt.Println("   to the wire itself: helix number H = (2πR)²/(μr·a²·acosh(p/2a)) ≪ 1. Two leaks are measured:")
	fmt.Println("   spread = (max−min)/max of dΦ/dt over the passes with the 10-pass primary driven, and the")
	fmt.Println("   2-pass winding's edge leakage (flux round its edges through the unwound neighbours), seen as")
	fmt.Println("   the coupling k and the swap ratio v(10-pass)/v(2-pass) with the 2-pass winding driven (ideal 5).")
	fmt.Println()
	fmt.Printf("   %-22s %-30s %9s %9s %9s %8s\n", "core", "spring", "H", "spread", "k", "swap")
	plain := func(m Material) Spring {
		s := DefaultSpring(m)
		s.A, s.P, s.RetArea = 1e-3, 3e-3, 1
		return s
	}
	type lever struct {
		name string
		mod  func(*Spring)
	}
	levers := []lever{
		{"plain: a 1 mm, gap 1 mm, ret 1×", func(*Spring) {}},
		{"gap 1 → 8 mm", func(s *Spring) { s.P = 2*s.A + 8e-3 }},
		{"gap 1 → 18 mm", func(s *Spring) { s.P = 2*s.A + 18e-3 }},
		{"return 1× → 16×", func(s *Spring) { s.RetArea = 16 }},
		{"core wire a 1 → 2 mm", func(s *Spring) { s.A = 2e-3; s.P = 2*s.A + 1e-3 }},
		{"design: a 2, gap 8, ret 16×", func(s *Spring) { *s = DefaultSpring(s.Mat) }},
	}
	for _, m := range Materials {
		for _, l := range levers {
			s := plain(m)
			l.mod(&s)
			lk, err := s.Leakage()
			if err != nil {
				return err
			}
			fmt.Printf("   %-22s %-30s %9.3g %9.2g %9.4f %8.3f\n", m.Name, l.name, s.HelixNumber(), lk.Spread, lk.K, lk.Swap)
		}
	}
	fmt.Println("   Wider gaps help only logarithmically (gap permeance ∝ 1/acosh(p/2a)); the thick return leg removes")
	fmt.Println("   the spread; the edge leakage falls with μr and (a/R)². With the nanocrystalline design the flux")
	fmt.Println("   stays on the wire: spread 1e-4, k 0.989. Silicon steel and Galfenol cannot reach that at this size.")
	fmt.Println()

	// 2. The ratio and the rate of change, design spring.
	fmt.Println("2. Ratio and rate of change (design spring, 100-turn windings). Ideal: v_s/v_p = (100·2)/(100·10) = 0.2;")
	fmt.Println("   per volt on the primary, dΦ/dt = 1/(N_p·n_p) = 1 mWb/s in every pass, and the bore (axis) wire gives N·dΦ/dt.")
	fmt.Printf("   %-22s %9s %10s %10s %11s %11s\n", "core", "H", "ratio", "ideal", "spread", "axis/v_p")
	for _, m := range Materials {
		s := DefaultSpring(m)
		pri, sec := designWindings(s.N)
		axis := Winding{Name: "axis", Azimuth: -math.Pi / 2, First: 0, Passes: s.N, Turns: 1}
		r, err := s.Ratio(pri, sec)
		if err != nil {
			return err
		}
		ra, err := s.Ratio(pri, axis)
		if err != nil {
			return err
		}
		f, _ := s.unitField(pri)
		fmt.Printf("   %-22s %9.3g %10.6f %10.6f %11.2g %11.6f (ideal %.3f)\n", m.Name, s.HelixNumber(), r, IdealRatio(pri, sec), spread(s.PassFlux(f)), ra, IdealRatio(pri, axis))
	}
	fmt.Println()

	// 3. Core turns N.
	fmt.Println("3. Core turns N (primary round all N passes, secondary round 2, 100 turns each; nanocrystalline design):")
	fmt.Printf("   %4s %10s %10s %14s %12s\n", "N", "ratio", "2/N", "dΦ/dt per V", "spread")
	for _, n := range []int{2, 4, 6, 10, 16, 20, 30, 40} {
		s := DefaultSpring(Materials[2])
		s.N = n
		pri, sec := designWindings(n)
		r, err := s.Ratio(pri, sec)
		if err != nil {
			return err
		}
		f, _ := s.unitField(pri)
		fmt.Printf("   %4d %10.6f %10.6f %11.3g Wb/s %12.2g\n", n, r, IdealRatio(pri, sec), 1/float64(pri.Turns*pri.Passes), spread(s.PassFlux(f)))
	}
	fmt.Println("   The core's turn count enters the ratio exactly like winding turns: 2/N. Per primary volt,")
	fmt.Println("   dΦ/dt falls as 1/N (the primary links the flux N_p·N times), so each pass carries less flux.")
	fmt.Println()

	// 4. Passes enclosed by the secondary, and the primary round fewer passes.
	fmt.Println("4. Passes enclosed (nanocrystalline design):")
	fmt.Printf("   %-40s %10s %10s %9s\n", "windings", "ratio", "ideal", "spread")
	s := DefaultSpring(Materials[2])
	cases := []struct {
		name     string
		pri, sec Winding
	}{}
	for _, ns := range []int{1, 2, 5, 10} {
		p, q := designWindings(10)
		q.Passes, q.First = ns, 0
		cases = append(cases, struct {
			name     string
			pri, sec Winding
		}{fmt.Sprintf("primary 10 passes, secondary %d", ns), p, q})
	}
	for _, np := range []int{2, 5} {
		p, q := designWindings(10)
		p.Passes, p.First = np, 0
		q.First = 8
		cases = append(cases, struct {
			name     string
			pri, sec Winding
		}{fmt.Sprintf("primary %d passes (1–%d), secondary 9–10", np, np), p, q})
	}
	for _, c := range cases {
		r, err := s.Ratio(c.pri, c.sec)
		if err != nil {
			return err
		}
		f, _ := s.unitField(c.pri)
		fmt.Printf("   %-40s %10.6f %10.6f %9.2g\n", c.name, r, IdealRatio(c.pri, c.sec), spread(s.PassFlux(f)))
	}
	fmt.Println("   A primary round only some passes drives the others' flux through the core and the return leg,")
	fmt.Println("   so each part sees the potential difference that pushes flux across the gaps: dΦ/dt is no longer the same in every pass.")
	fmt.Println()

	// 5. kronspice.
	r := DefaultRig()
	fmt.Printf("5. kronspice: %.0f V peak at %.0f Hz, %.1f Ω source, %.0f %% declared winding air leakage, copper windings.\n", r.Vpk, r.Freq, r.Rsrc, 100*r.AirLeak)
	fmt.Printf("   %-22s %9s %10s %9s %9s %9s %9s %8s %11s %9s\n", "core", "k(p,s)", "ratio AC", "axis AC", "Rload", "Pin", "Pout", "η", "ledger res", "B peak")
	for _, m := range Materials {
		s := DefaultSpring(m)
		pri, sec := designWindings(s.N)
		cr, err := s.RunCircuit(r, pri, sec)
		if err != nil {
			return fmt.Errorf("%s: %w", m.Name, err)
		}
		fmt.Printf("   %-22s %9.6f %10.6f %9.6f %7.3g Ω %7.3g W %7.3g W %8.3f %11.2g %7.3f T\n", m.Name, cr.K[0], cr.RatioAC, cr.RatioAxis, cr.Rload, cr.Pin, cr.Pout, cr.Pout/cr.Pin, cr.Mismatch/cr.EIn, cr.Bpk)
	}
	s = DefaultSpring(Materials[2])
	pri, sec := designWindings(s.N)
	rs := r
	rs.Vpk = r.Vpk * IdealRatio(pri, sec) // same core flux as the forward run
	sw, err := s.RunCircuit(rs, sec, pri)
	if err != nil {
		return err
	}
	fmt.Printf("   swap (drive the 2-pass winding at %.0f V peak, nanocrystalline): v_10-pass/v_2-pass = %.4f (ideal 5), η %.3f, ledger %.1g, B %.3f T\n", rs.Vpk, sw.RatioAC, sw.Pout/sw.Pin, sw.Mismatch/sw.EIn, sw.Bpk)
	fmt.Println("   Eddy currents and hysteresis are not modelled. Skin depth and saturation at this drive:")
	for _, m := range Materials {
		sm := DefaultSpring(m)
		ok := "eddy currents negligible"
		if !sm.EddyOK(r.Freq) {
			ok = "eddy currents NOT negligible (solid wire thicker than δ): laminate, strand or use tape"
		}
		fmt.Printf("     %-22s δ = %.3g mm, B peak / Bs = %.2f — %s\n", m.Name, m.SkinDepth(r.Freq)*1e3, sw.Bpk/m.Bsat, ok)
	}
	fmt.Println()

	// 6. The helical field (B(3)).
	cr, err := s.RunCircuit(r, pri, sec)
	if err != nil {
		return err
	}
	phi := cr.Bpk * s.Area()
	h := s.HelixField(phi)
	fmt.Println("6. The helical field (nanocrystalline design, design drive), one pass at its peak:")
	fmt.Printf("   flux Φ %.3g Wb, B along the wire %.3f T, pitch angle %.2f°\n", phi, h.B, h.Alpha*180/math.Pi)
	fmt.Printf("   transverse part B(0) = B·cos α = %.4f T, turning once per %.0f mm of height (κ = %.0f 1/m)\n", h.Trans, s.P*1e3, h.Kappa)
	fmt.Printf("   B Cyclic theorem B(1)×B(2) = i·B(0)·B(3)*:  |B(3)| = B(0) = %.4f T along the axis (Evans, AIAS Essay 36)\n", h.B3Cyclic)
	fmt.Printf("   standard: B(1)×B(2) = i·B(0)²·k, |.| = %.4f T² (an axial vector, not a field); field along the axis in the wire B·sin α = %.4f T;\n", h.Product, h.Axial)
	fmt.Println("             net axial flux through any cut of the spring = Φ (returned by the leg)")
	e, bo, wl := s.PickupPerDPhi(2)
	fmt.Printf("   pickups per dΦ/dt: a turn round the whole spring %.0f×, the wire up the bore %.0f×, a wall turn round 2 passes %.0f×\n", e, bo, wl)
	fmt.Println("   The bore wire reads the circulation (transverse part, N times); a coil round the outside reads the axial part once.")
	fmt.Println()

	// 7. Energy ledger against core turns.
	fmt.Println("7. Energy ledger against core turns N (nanocrystalline; primary 100 turns round all N passes, secondary")
	fmt.Println("   100 turns round 2; 20 V peak, 1 kHz, 20 cycles, one fixed 632 Ω load). Energies over the whole run [mJ]:")
	fmt.Printf("   %4s %8s %10s %7s %9s %9s %9s %9s %9s %9s %10s %9s %9s %6s\n", "N", "ratio", "dΦ/dt pk", "B pk", "E in", "E back", "net in", "E load", "E copper", "E stored", "residual", "Pin", "Pout", "η")
	rows, err := LedgerSweep([]int{2, 4, 6, 10, 16, 20, 30, 40})
	if err != nil {
		return err
	}
	for _, l := range rows {
		fmt.Printf("   %4d %8.4f %7.4f Wb/s %5.3f T %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %10.1e %7.3g W %7.3g W %6.3f\n", l.N, l.Ratio, l.DPhi, l.Bpk,
			l.EIn*1e3, l.EBack*1e3, (l.EIn-l.EBack)*1e3, l.ELoad*1e3, (l.ECu+l.ESrcR)*1e3, l.EStore*1e3, l.Mismatch/l.EIn, l.Pin, l.Pout, l.Pout/l.Pin)
		if l.Bpk > Materials[2].Bsat {
			fmt.Printf("        ↑ B above the tape's %.1f T saturation: outside the linear model (a real core would saturate)\n", Materials[2].Bsat)
		}
	}
	fmt.Println("   E in is what the source delivers, E back what the reactive (magnetizing) current returns to it each cycle;")
	fmt.Println("   net in = load + copper (incl. the source's 0.5 Ω) + stored, to the residual. At a fixed drive and load, v_s = 2·v_p/N,")
	fmt.Println("   so the power delivered falls as 1/N²; the input falls with it (the core turns raise the primary's")
	fmt.Println("   inductance as N², so its magnetizing current falls too). Efficiency stays near 1: the core turns change")
	fmt.Println("   the ratio and the rate of change per pass, never the energy balance.")
	fmt.Println()

	// 8. As a magnetic amplifier.
	fmt.Println("8. As a magnetic amplifier (first-order estimate; the core here is linear, saturation not simulated).")
	fmt.Println("   A DC control winding saturates the core; the gate (load) winding blocks AC until the core saturates.")
	fmt.Println("   Ampere-turn balance: I_load·N_g·n_g = I_c·N_c·n_c, so the current gain is (N_c·n_c)/(N_g·n_g).")
	fmt.Printf("   %-44s %9s %9s %9s %8s %8s %10s\n", "windings (nanocrystalline design, 60 Hz, 100 Ω load)", "I_c sat", "V gate", "Ki", "Ki plain", "Rc", "Kp")
	sm := DefaultSpring(Materials[2])
	for _, c := range []struct {
		name string
		c, g Winding
	}{
		{"control 100 t × 10 passes, gate 100 t × 2", Winding{"control", 0, 0, 10, 100}, Winding{"gate", math.Pi, 4, 2, 100}},
		{"control 100 t × 10 passes, gate 100 t × 10", Winding{"control", 0, 0, 10, 100}, Winding{"gate", math.Pi, 0, 10, 100}},
		{"control 100 t × 2 passes, gate 100 t × 10", Winding{"control", 0, 4, 2, 100}, Winding{"gate", math.Pi, 0, 10, 100}},
		{"control 500 t × 10 passes, gate 50 t × 2", Winding{"control", 0, 0, 10, 500}, Winding{"gate", math.Pi, 4, 2, 50}},
	} {
		m := sm.MagAmpEstimate(c.c, c.g, 60, 100, r.CuDia)
		fmt.Printf("   %-44s %6.2f mA %7.1f V %9.2f %8.2f %6.2f Ω %10.3g\n", c.name, m.IcSat*1e3, m.Vgate, m.Ki, m.KiPlain, m.Rc, m.Kp)
	}
	fmt.Println("   The core passes multiply each winding's effective turns, so a control winding round all the passes and a")
	fmt.Println("   gate round a few gains (n_c/n_g)× more current gain than the same turns on an ordinary core. It is the")
	fmt.Println("   gain of an ordinary core wound N·n turns: the spring sets those ratios with one winding each, it adds no power.")
	return nil
}
