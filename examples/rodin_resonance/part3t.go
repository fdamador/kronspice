// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"strconv"
	"time"

	"github.com/fdamador/kronspice/exchange"
)

// part3t runs Part 3t: the output winding wound round the torus as a
// toroidal inductor winding (through the hole and over the outside), in
// 12 sections 30° apart, the drive frequency and the sections evolved
// together, the winner charging the bank in a kronspice transient.
func part3t(m *Model, v float64) {
	fmt.Println("\nPart 3t — the output winding as a toroidal inductor winding in 12 sections round the torus")
	if v > 0 {
		Part3Volts = v
	}
	kit, err := NewTorKit(m.P)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  drive: sine on A, cosine on B, %.4g V peak each through %.0f Ω, as Part 3\n", Part3Volts, SourceOhms)
	fmt.Printf("  winding: 0.5 mm copper (enamel %.2f mm), every turn round the tube's cross-section, through the hole and over the outside,\n", TorEnamel*1e3)
	fmt.Printf("  laid on winding B's outer surface (first layer's turn radius a0 = %.3f mm round the tube's centre circle R0 = %.4f mm);\n", kit.A0*1e3, TorusR0*1e3)
	fmt.Printf("  12 sections centred at θ0 + k·30° (pole 1 at 0°, the wire's hole crossings at 15° + k·30°), arc %.0f°–%.0f° each, up to %d layers;\n",
		TorMinArc*180/math.Pi, TorMaxArc*180/math.Pi, TorMaxLayers)
	fmt.Print("  turns that fit on the inside of the hole per section (arc → turns in 1 … 6 layers):")
	for _, deg := range []float64{10, 20, 29} {
		n := 0
		for kit.Layers(n+1, deg*math.Pi/180) > 0 {
			n++
		}
		fmt.Printf(" %.0f° → %d;", deg, n)
	}
	fmt.Println()

	// 1. checks of the tables
	c := kit.Checks(m, 114, 8, 20)
	fmt.Println("  checks of the inductance tables:")
	fmt.Printf("    turn–pass mutuals, interpolated and carried round by the 30° symmetry vs direct, at 7.3°, 133.71°, 251.06°, 359.9°: largest error %.1e of the largest\n", c.TableErr)
	fmt.Printf("    closely wound uniform toroid of %d turns (one layer all round): turns summed %.5g µH vs current sheet μ0N²(R0 − √(R0² − a²)) %.5g µH (%+.2e)\n",
		c.NTor, c.TorRings*1e6, c.TorSheet*1e6, c.TorRings/c.TorSheet-1)
	fmt.Printf("    one section of 8 turns over 20°: as rings vs as a helix closed by its return lead along the outer equator —\n")
	fmt.Printf("      mutual to one circuit of A %.5g vs %.5g nH (%+.2e), of B %.5g vs %.5g nH (%+.2e), self %.5g vs %.5g µH (%+.2e; the helix adds its lead loop)\n",
		c.RingMA*1e9, c.HelMA*1e9, c.RingMA/c.HelMA-1, c.RingMB*1e9, c.HelMB*1e9, c.RingMB/c.HelMB-1, c.RingL*1e6, c.HelL*1e6, c.RingL/c.HelL-1)

	// 2. one turn round the torus at every angle: its EMF, full sectional net
	nphi := 720
	nf := 90
	if *short {
		nphi, nf = 180, 24
	}
	phis := make([]float64, nphi)
	for i := range phis {
		phis[i] = 2 * math.Pi * float64(i) / float64(nphi)
	}
	open := m.NewNet(OpenA, 0, 1e6)
	fo, _, _ := open.Modes()
	fmt.Printf("  one turn round the tube at every azimuth (%d angles), its open-circuit EMF from the full sectional net (12 sections per circuit):\n", nphi)
	rows := [][]string{{"phi_deg", "tor_sincos_amp", "tor_sincos_ph", "tor_A_amp", "tor_B_amp", "tor_inphase_amp", "dip_sincos_amp", "dip_sincos_ph", "dip_A_amp", "dip_B_amp", "dip_inphase_amp"}}
	cols := map[string][]complex128{}
	for _, f := range []float64{1e3, 1e4, 1e5, 1e6, fo[0], 5e6} {
		for _, sense := range []int{Toroidal, Dipole} {
			sn := map[int]string{Toroidal: "toroidal", Dipole: "dipole"}[sense]
			net := m.NewNet(DrivenPair, sense, f)
			sc, err1 := kit.TurnEMF(m, net, f, sense, 0, 0, phis)
			ea, err2 := kit.TurnEMF(m, net, f, sense, 0, 1, phis)
			eb, err3 := kit.TurnEMF(m, net, f, sense, 0, 2, phis)
			ip, err4 := kit.TurnEMF(m, net, f, sense, -90, 0, phis)
			if err := firstErr(err1, err2, err3, err4); err != nil {
				log.Fatal(err)
			}
			st := func(z []complex128) (mean, lo, hi float64) {
				lo = math.Inf(1)
				for _, x := range z {
					a := cmplx.Abs(x)
					mean += a
					lo, hi = math.Min(lo, a), math.Max(hi, a)
				}
				return mean / float64(len(z)), lo, hi
			}
			ms, los, his := st(sc)
			ma, _, _ := st(ea)
			mb, _, _ := st(eb)
			mi, _, _ := st(ip)
			var sum complex128
			for _, x := range sc {
				sum += x
			}
			amps := make([]float64, len(sc))
			plo, phi := math.Inf(1), math.Inf(-1)
			ref := cmplx.Phase(sum)
			for i, x := range sc {
				amps[i] = cmplx.Abs(x)
				d := math.Remainder(cmplx.Phase(x)-ref, 2*math.Pi) * 180 / math.Pi
				plo, phi = math.Min(plo, d), math.Max(phi, d)
			}
			hm, ha := dominantHarmonic(amps)
			fmt.Printf("    f = %8.4g Hz, %-8s sense: |EMF| per turn mean %.4g mV (over θ %.4g … %.4g, m = %d at %.2g of the mean), phase over θ %+.1f … %+.1f°;\n",
				f, sn, ms*1e3, los*1e3, his*1e3, hm, ha, plo, phi)
			fmt.Printf("      A alone %.4g mV, B alone %.4g mV, √(A² + B²) %.4g mV; A and B in phase %.4g mV; the whole ring of turns summed (all alike) %.4g mV per turn\n",
				ma*1e3, mb*1e3, math.Hypot(ma, mb)*1e3, mi*1e3, cmplx.Abs(sum)/float64(len(sc))*1e3)
			if f == 1e3 {
				cols[sn+"sc"], cols[sn+"a"], cols[sn+"b"], cols[sn+"ip"] = sc, ea, eb, ip
			}
		}
	}
	for i, p := range phis {
		r := []string{f6(p * 180 / math.Pi)}
		for _, sn := range []string{"toroidal", "dipole"} {
			r = append(r, f6(cmplx.Abs(cols[sn+"sc"][i])), f6(cmplx.Phase(cols[sn+"sc"][i])*180/math.Pi),
				f6(cmplx.Abs(cols[sn+"a"][i])), f6(cmplx.Abs(cols[sn+"b"][i])), f6(cmplx.Abs(cols[sn+"ip"][i])))
		}
		rows = append(rows, r)
	}
	writeCSV("tor_turn_emf.csv", rows)
	// frequency × angle map of |EMF| per turn, both senses
	var freqs []float64
	for i := 0; i < nf; i++ {
		freqs = append(freqs, 1e3*math.Pow(3e7/1e3, float64(i)/float64(nf-1)))
	}
	for _, sense := range []int{Toroidal, Dipole} {
		sn := map[int]string{Toroidal: "toroidal", Dipole: "dipole"}[sense]
		mp := [][]string{{"f"}}
		for k := 0; k < len(phis); k += 4 {
			mp[0] = append(mp[0], f6(phis[k]*180/math.Pi))
		}
		net := m.NewNet(DrivenPair, sense, 1e3)
		for _, f := range freqs {
			e, err := kit.TurnEMF(m, net, f, sense, 0, 0, phis)
			if err != nil {
				log.Fatal(err)
			}
			r := []string{f6(f)}
			for k := 0; k < len(phis); k += 4 {
				r = append(r, f6(cmplx.Abs(e[k])))
			}
			mp = append(mp, r)
		}
		writeCSV("tor_map_"+sn+".csv", mp)
	}

	// 3. the genetic search on the reduced net (2 sections per circuit):
	// 48 sections + 12 output sections = 60 windings, within the N element
	mr := NewModel(NewPairSec(Rodin, Circuits, 2))
	tp := NewTorPart(kit, mr)
	pop, gens := 64, 50
	if *short {
		pop, gens = 16, 5
	}
	fmt.Println("\n  model for the searches: the sectional net with 2 sections per circuit (48 + 12 sections = 60 windings, within the kronspice N element)")
	fmt.Println("  constraints: turns fit in ≤ 6 layers on the inside of the hole; wire ≤ λ/20 (each section, the whole winding in series);")
	fmt.Println("  tuning capacitor (equivalent single capacitor for the series string) ≥ 30 pF; fitness: P_load × bridge factor, as Part 3")
	torSearch(kit, mr, tp, "sine on A, cosine on B (the drive of Parts 2 and 3)", "", GAOptions{Population: pop, Generations: gens, Seed: 45}, false)
	torSearch(kit, mr, tp, "the phase between the two drives evolved too", "_phase", GAOptions{Population: pop, Generations: gens, Seed: 46}, true)
}

// torSearch runs one genetic search for the toroidal winding, reports the
// winner, the variations round it, and its kronspice transient.
func torSearch(kit *TorKit, mr *Model, tp *TorPart, title, tag string, o GAOptions, phase bool) {
	fmt.Printf("\n  genetic search, %s (evo, seed %d, %d × %d): genes", title, o.Seed, o.Population, o.Generations)
	for _, g := range TorSpace(phase) {
		fmt.Printf(" | %s ∈ [%g, %g]", g.Name, g.Min, g.Max)
	}
	fmt.Println()
	t0 := time.Now()
	ga, err := tp.Evolve(o, phase)
	if err != nil {
		log.Fatal(err)
	}
	hist := [][]string{{"generation", "best", "mean", "failed", "evaluations"}}
	for _, h := range ga.History {
		hist = append(hist, []string{strconv.Itoa(h.Generation), f6(h.Best), f6(h.Mean), strconv.Itoa(h.Failed), strconv.Itoa(h.Evaluations)})
	}
	writeCSV("tor_ga_history"+tag+".csv", hist)
	l := ga.Best
	ac := ga.AC
	fmt.Printf("  fitness P_load × bridge factor: best of generation 0 %.4g W; evolved after %d evaluations (%.0f s) %.4g W (P_load %.4g W, bridge factor %.3f)\n",
		math.Pow(10, 6-ga.History[0].Best), ga.Evals, time.Since(t0).Seconds(), ga.Power*BridgeFactor(ac), ga.Power, BridgeFactor(ac))
	fmt.Printf("  winner: %v\n", l)
	reportTor(kit, mr, tp, l, ac, tag)

	// 4. around the winner
	fmt.Print("    sections moved round the torus (θ0, all else fixed), P_load [W]:")
	th := [][]string{{"theta0_deg", "pload", "eff"}}
	for d := 0.0; d <= 30.0001; d += 1.5 {
		ll := l
		ll.Theta0 = d * math.Pi / 180
		pw, r, err := tp.Score(ll, 0)
		if err != nil {
			continue
		}
		th = append(th, []string{f6(d), f6(pw), f6(pw / (r.PSrc[0] + r.PSrc[1]))})
		if math.Mod(d, 7.5) < 1e-9 {
			fmt.Printf(" %.1f°: %.4g;", d, pw)
		}
	}
	fmt.Println()
	writeCSV("tor_theta0"+tag+".csv", th)
	fmt.Print("    section arc (turns kept, layers recomputed), P_load [W]:")
	for _, deg := range []float64{5, 10, 15, 20, 25, 29} {
		ll := l
		ll.Arc = deg * math.Pi / 180
		if pw, _, err := tp.Score(ll, 0); err == nil {
			fmt.Printf(" %.0f°: %.4g;", deg, pw)
		} else {
			fmt.Printf(" %.0f°: infeasible;", deg)
		}
	}
	fmt.Println()
	fmt.Print("    turns per section (arc kept), P_load [W]:")
	for _, n := range []int{1, 2, 4, 8, 12, 16, 24, 32, 40} {
		ll := l
		ll.Turns = n
		if pw, _, err := tp.Score(ll, 0); err == nil {
			fmt.Printf(" %d: %.4g;", n, pw)
		} else {
			fmt.Printf(" %d: infeasible;", n)
		}
	}
	fmt.Println()
	for _, cmb := range []Combine{Separate, Series, Parallel} {
		for _, alt := range []bool{false, true} {
			ll := l
			ll.Comb, ll.Alt = cmb, alt
			pol := "all alike"
			if alt {
				pol = "alternating"
			}
			if pw, r, err := tp.Score(ll, 0); err == nil {
				fmt.Printf("    combined %-32s sections %-11s %.4g W (efficiency %.4f, bridge factor %.3f)\n", cmb.String()+",", pol+":", pw, pw/(r.PSrc[0]+r.PSrc[1]), BridgeFactor(r))
			} else {
				fmt.Printf("    combined %-32s sections %-11s %v\n", cmb.String()+",", pol+":", err)
			}
		}
	}
	for _, d := range []struct {
		name string
		ph   float64
	}{{"sine on A, cosine on B", 0}, {"A and B in phase", -90}, {"A and B in antiphase", 90}, {"cosine on A, sine on B", 180}} {
		ll := l
		ll.PhaseB = 0
		if pw, r, err := tp.Score(ll, d.ph); err == nil {
			fmt.Printf("    drive %-36s P_load %.4g W, sources %.4g W, efficiency %.4f\n", d.name+":", pw, r.PSrc[0]+r.PSrc[1], pw/(r.PSrc[0]+r.PSrc[1]))
		}
	}
	other := l
	other.Sense = -l.Sense
	if pw, r, err := tp.Score(other, 0); err == nil {
		fmt.Printf("    the primaries connected in the other sense: P_load %.4g W (efficiency %.4f)\n", pw, pw/(r.PSrc[0]+r.PSrc[1]))
	}
	fmt.Print("    the same winding at other drive frequencies (sections retuned), P_load [W]:")
	for _, f := range []float64{1e3, 2e3, 4e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6} {
		ll := l
		ll.F = f
		if pw, _, err := tp.Score(ll, 0); err == nil {
			fmt.Printf(" %.3g Hz: %.3g;", f, pw)
		} else {
			fmt.Printf(" %.3g Hz: infeasible;", f)
		}
	}
	fmt.Println()
	for _, sec := range []int{4, 6, 12} {
		mm := NewModel(NewPairSec(Rodin, Circuits, sec))
		pw, _, err := NewTorPart(kit, mm).Score(l, 0)
		fmt.Printf("    discretization check: the same winding on %2d sections per circuit gives %.4g W (%v)\n", sec, pw, errStr(err))
	}

	// the winner charging the bank (kronspice transient)
	torTransient(kit, mr, tp, l, ac, tag)
}

func firstErr(es ...error) error {
	for _, e := range es {
		if e != nil {
			return e
		}
	}
	return nil
}

// torLumped returns the winding's string inductance (all sections in
// series with the layout's polarity), resistance, and its mutual to the
// pair (A + sense·B, path sense) on net n of model m.
func torLumped(kit *TorKit, m *Model, n *Net, l TorLayout) (ls, rs, msp, ma, mb float64) {
	gs, mg, err := kit.Groups(m, n, l)
	if err != nil {
		return
	}
	for i := range mg {
		for j := range mg[i] {
			ls += mg[i][j]
		}
		rs += gs[i].R
		for g, v := range gs[i].MSec {
			if g < n.NS {
				ma += v
			} else {
				mb += v
			}
		}
	}
	msp = ma + float64(l.Sense)*mb
	return
}

func reportTor(kit *TorKit, mr *Model, tp *TorPart, l TorLayout, ac *ACResult, tag string) {
	n := tp.net(l.Sense)
	n.SetFreq(l.F)
	gs, mg, _ := kit.Groups(mr, n, l)
	ly := kit.Layers(l.Turns, l.Arc)
	ls, rs, msp, ma, mb := torLumped(kit, mr, n, l)
	lp := mr.LA + mr.LB + 2*float64(l.Sense)*mr.M
	fmt.Printf("    %d layer(s), mean turn radius %.3f mm; wire %.4g m per section, %.4g m in all (λ/20 at f: %.4g m)\n",
		ly, kit.LayerRadius(ly)*1e3, kit.SectionWire(l), TorSections*kit.SectionWire(l), 3e8/(20*l.F))
	fmt.Printf("    section L %.4g µH, R %.4g Ω at f, tuning C %.4g nF; neighbouring sections' mutual %.4g µH (k %.3f), across the torus %.3g µH\n",
		mg[0][0]*1e6, gs[0].R, gs[0].Ctune*1e9, mg[0][1]*1e6, mg[0][1]/mg[0][0], mg[0][TorSections/2]*1e6)
	fmt.Printf("    whole winding in series: L %.4g µH, R %.4g Ω; mutual to winding A %.4g µH, to B %.4g µH (path sense); to the series pair in its sense %.4g µH,\n",
		ls*1e6, rs, ma*1e6, mb*1e6, msp*1e6)
	fmt.Printf("    coupling to each winding k = M/√(L_winding·L): to A %.3f, to B %.3f; to the series pair in its sense %.3f (L_pair = %.4g µH)\n",
		ma/math.Sqrt(ls*mr.LA), mb/math.Sqrt(ls*mr.LB), msp/math.Sqrt(ls*lp), lp*1e6)
	if math.Abs(math.Remainder(l.PhaseB, 180)) < 45 {
		fmt.Println("    (with the two sources 90° apart the two couplings neither add nor cancel, so the pair figure says little)")
	}
	fmt.Printf("    AC books: sources %.4g + %.4g W = load %.4g W + copper %.4g W + source resistors %.4g W (difference %.1e W)\n",
		ac.PSrc[0], ac.PSrc[1], ac.PLoad, ac.PCopper, ac.PSrcR, ac.PSrc[0]+ac.PSrc[1]-ac.PLoad-ac.PCopper-ac.PSrcR)
	fmt.Printf("    section current phases [°]: %.1f\n    load voltage amplitudes [V]: %.4g\n", GroupPhases(ac), absAll(ac.VLoad))
	src0 := n.DriveSources(l.Sense, Part3Volts)
	src0[1].PhaseDeg += l.PhaseB
	r0, _ := n.SolveAC(l.F, src0, nil, nil, Separate, 1)
	fmt.Printf("    reaction on the drives: without the winding the sources give %.4g + %.4g W (|I| %.4g, %.4g A); with it %.4g + %.4g W (|I| %.4g, %.4g A)\n",
		r0.PSrc[0], r0.PSrc[1], cmplx.Abs(r0.ISrc[0]), cmplx.Abs(r0.ISrc[1]), ac.PSrc[0], ac.PSrc[1], cmplx.Abs(ac.ISrc[0]), cmplx.Abs(ac.ISrc[1]))
	fmt.Printf("    efficiency at AC (load / sources) %.4f\n", ac.PLoad/(ac.PSrc[0]+ac.PSrc[1]))
	rows := [][]string{{"section", "centre_deg", "arc_deg", "turns", "layers", "radius", "sign"}}
	for s := 0; s < TorSections; s++ {
		rows = append(rows, []string{strconv.Itoa(s), f6((l.Theta0 + float64(s)*math.Pi/6) * 180 / math.Pi), f6(l.Arc * 180 / math.Pi),
			strconv.Itoa(l.Turns), strconv.Itoa(ly), f6(kit.LayerRadius(ly)), f6(l.sign(s))})
	}
	writeCSV("tor_sections"+tag+".csv", rows)
}

func torTransient(kit *TorKit, mr *Model, tp *TorPart, l TorLayout, ac *ACResult, tag string) {
	n := tp.net(l.Sense)
	n.SetFreq(l.F)
	gs, mg, err := kit.Groups(mr, n, l)
	if err != nil {
		log.Fatal(err)
	}
	ls, rs, _, _, _ := torLumped(kit, mr, n, l)
	wl := 2 * math.Pi * l.F
	var ql float64
	switch l.Comb {
	case Series:
		ql = wl * ls / (rs + l.RLoad)
	case Parallel:
		ql = wl * (ls / (TorSections * TorSections)) / (rs/(TorSections*TorSections) + l.RLoad)
	default:
		ql = wl * (ls / TorSections) / (rs/TorSections + l.RLoad)
	}
	tauDrive := mr.LA / (SourceOhms + WireResistance(mr.LenA, l.F))
	cyc := int(math.Max(10*ql/math.Pi, 10*tauDrive*l.F))
	cyc = max(400, min(cyc, 8000))
	steps := 60
	if *short {
		cyc, steps = min(cyc, 300), 40
	}
	run := float64(cyc) / l.F
	vpk := 0.0
	for _, v := range ac.VLoad {
		vpk = math.Max(vpk, cmplx.Abs(v))
	}
	bank := e6Down(ac.PLoad * run / (vpk * vpk))
	t0 := time.Now()
	tr, err := TranGroups(n, l.F, l.Sense, l.PhaseB, gs, mg, l.Comb, l.RLoad, cyc, steps, bank)
	if err != nil {
		log.Fatal("transient: ", err)
	}
	fmt.Printf("\n  kronspice transient of the winner: %d cycles (%.4g ms) at %d steps per cycle (winding loaded Q %.3g, drive L/R %.3g ms),\n", cyc, run*1e3, steps, ql, tauDrive*1e3)
	fmt.Printf("  48 primary sections and the 12 output sections as the 60 windings of one N element, Schottky bridges (as Part 3), bank %.4g µF\n", bank*1e6)
	fmt.Printf("  sized to fill within the run (AC load voltage %.4g V peak); %.0f s\n", vpk, time.Since(t0).Seconds())
	esrc := tr.ESrc[0] + tr.ESrc[1]
	nT := len(tr.Time)
	fmt.Printf("    bank: %.4g V, %.4g µJ at the end (%.4g ms)\n", tr.VBank[nT-1], tr.EBankEnd*1e6, tr.Time[nT-1]*1e3)
	for _, frac := range []float64{0.1, 0.25, 0.5, 1} {
		k := int(frac*float64(nT-1) + 0.5)
		fmt.Printf("      at %.3g ms: %.4g V, %.4g µJ\n", tr.Time[k]*1e3, tr.VBank[k], tr.EBank[k]*1e6)
	}
	fmt.Printf("    energy from the sources: A %.4g µJ, B %.4g µJ (total %.4g µJ)\n", tr.ESrc[0]*1e6, tr.ESrc[1]*1e6, esrc*1e6)
	fmt.Printf("    split: bank %.4g µJ (%.2f %%), copper %.4g µJ (%.2f %%), source resistors %.4g µJ (%.2f %%), diodes %.4g µJ (%.2f %%), bleeders and joints %.3g µJ, left in coils and capacitors %.4g µJ\n",
		tr.EBankEnd*1e6, 100*tr.EBankEnd/esrc, tr.ECopper*1e6, 100*tr.ECopper/esrc, tr.ESrcR*1e6, 100*tr.ESrcR/esrc, tr.EDiode*1e6, 100*tr.EDiode/esrc, tr.EBleed*1e6, tr.ECoilsEnd*1e6)
	fmt.Printf("    charging efficiency (bank energy / source energy) %.4f; books: sources − (bank + copper + source R + diodes + bleed + stored) = %.3g µJ; ledger integration mismatch %.3g µJ (%.1e of the energy in), largest point residual %.1e\n",
		tr.EBankEnd/esrc, (esrc-tr.EBankEnd-tr.ECopper-tr.ESrcR-tr.EDiode-tr.EBleed-tr.ECoilsEnd)*1e6, tr.Mismatch*1e6, tr.Mismatch/tr.EIn, tr.MaxRel)
	fmt.Printf("    external channels other than electrical and thermal (%v among them) carry nothing: every joule is booked against the two drive sources\n", exchange.Neutrino)
	rows := [][]string{{"t", "vbank", "ebank"}}
	step := max(1, nT/3000)
	for k := 0; k < nT; k += step {
		rows = append(rows, []string{f6(tr.Time[k]), f6(tr.VBank[k]), f6(tr.EBank[k])})
	}
	writeCSV("tor_bank"+tag+".csv", rows)
}
