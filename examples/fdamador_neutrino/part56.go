// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The report for Parts 5 and 6, and the data the charts are drawn from.

import (
	"fmt"
	"math"
	"os"
	"time"
)

// Part5a holds the circuit implosion results.
type Part5a struct {
	Set    ImplosionSet
	Locked TrainRun
	Anti   TrainRun
}

func reportPart5a() (Part5a, error) {
	var p Part5a
	e := DefaultERing()
	fmt.Println("\n5. Implosion and explosion of the E-ring (H28 × H30)")
	fmt.Println("5a. Circuit: R_rad as a port to the outside (a line of impedance R_rad); an incoming wave on it is a source")
	fmt.Println("    2·v_inc behind R_rad; v_ref = v_inc − R_rad·i is the outgoing wave (reflection while driven, explosion after)")
	t0 := time.Now()
	set, err := e.Implosions()
	if err != nil {
		return p, err
	}
	p.Set = set
	fmt.Printf("    bench ring: f0 = %.5g GHz, Q = %.4g, τ_E = Q/ω0 = %.4g ns; growth rates: exact reversal (R_rad+R_d)/2L = %.5g /s,\n",
		e.Omega0()/2/math.Pi/1e9, e.Q(), set.TauE*1e9, set.AlphaTot)
	fmt.Printf("    no-reflection (CPA) (R_rad−R_d)/2L = %.5g /s; predicted reflection of the exact reversal R_d/(R_rad+R_d) = %.4g in amplitude (%.2g in power)\n",
		set.AlphaCPA, set.RefCPA, set.RefCPA*set.RefCPA)
	fmt.Printf("    every drive carries the same incident energy %.4g J over %.0f τ_E (growing drives start at e^-8 of their final power), then off;\n",
		set.Runs[0].EInc, set.WindowTauE)
	fmt.Println("    explosion followed for as long again. Fractions of the incident energy:")
	fmt.Printf("    %-52s %9s %9s %9s %9s %9s %9s\n", "drive", "absorbed", "peak W", "explosion", "heat R_d", "ledger", "waves")
	for _, r := range set.Runs {
		fmt.Printf("    %-52s %9.6f %9.5f %9.5f %9.2e %9.1e %9.1e\n", r.Name, r.Absorbed, r.PeakStored/r.EInc, r.EOutAfter/r.EInc,
			r.HeatDiel/r.EInc, r.LedgerRel, r.WaveRel)
	}
	c := set.Runs[0]
	fmt.Printf("    CPA run, ledger: source work %.6g J = radiated (R_rad) %.6g J + heat (R_d) %.4g J + stored at end %.3g J; time-stepping %.1e\n",
		c.SourceWork, c.RadHeat, c.HeatDiel, c.StoredEnd, c.Mismatch)
	fmt.Printf("    wave books: incident %.6g J = reflected %.4g J + explosion %.6g J + heat %.4g J + stored %.3g J\n",
		c.EInc, c.ERefDrive, c.EOutAfter, c.HeatDiel, c.StoredEnd)
	fmt.Printf("    the CPA reflection %.2g of the incident energy is the drive's switch-on: e^-8 = %.2g\n", 1-c.Absorbed, math.Exp(-8))
	locked, err := e.Train("pulsing, locked to the free ring", 40, 6, 1)
	if err != nil {
		return p, err
	}
	anti, err := e.Train("pulsing, antiphase to the free ring", 40, 6, -1)
	if err != nil {
		return p, err
	}
	p.Locked, p.Anti = locked, anti
	fmt.Println("    Pulsing ring: an implosion (3 τ_E of CPA growth) every 40 periods, 6 cycles:")
	for _, tr := range []TrainRun{locked, anti} {
		fmt.Printf("    %s (ledger %.1e):\n", tr.Name, tr.LedgerRel)
		for i, cy := range tr.Cycles[:len(tr.Cycles)-1] {
			fmt.Printf("      cycle %d: in %.4g J, out %.4g J, heat %.3g J, stored %.3g → %.3g J, peak %.4g J\n",
				i+1, cy.EInc, cy.EOut, cy.Heat, cy.StoredStart, cy.StoredEnd, cy.Peak)
		}
	}
	l := locked.Cycles[len(locked.Cycles)-2]
	fmt.Printf("    steady cycle: the outside supplies %.4g J and takes back %.4g J; the difference %.4g J is the dielectric heat %.4g J\n",
		l.EInc, l.EOut, l.EInc-l.EOut, l.Heat)
	fmt.Printf("    (%.2f %% of the incoming wave); nothing is made inside. (%s)\n", 100*(l.EInc-l.EOut)/l.EInc, time.Since(t0).Round(time.Second))
	fmt.Println("    Standard prediction: the implosion is the explosion reversed; its energy is the declared incoming wave's, the")
	fmt.Println("    time-reversed drive is absorbed without reflection (coherent perfect absorption, Chong et al., PRL 105, 053901, 2010),")
	fmt.Println("    any other drive is partly reflected, and a pulsing ring needs one implosion's energy from outside every cycle.")
	return p, nil
}

func reportPart5b(cells []float64) ([]FieldResult, error) {
	fmt.Println("\n5b. Field: the dielectric torus (a = 10 mm, r = 2.5 mm, εr = 100, lossless) in em/fdtd.go, Huygens surface at ±2a,")
	fmt.Println("    absorber from ±3a; incoming wave = the time reverse of the ring's own outgoing radiation (phasor at f0 on the surface,")
	fmt.Println("    from a ring-down started from the filtered low mode), growing as e^{+γt} over 6 τ_E and switched off at an E maximum;")
	fmt.Println("    the explosion then crosses the surface into the absorber (followed 3 τ_E). Compared with a constant burst of the same")
	fmt.Println("    pattern and incident energy. 'Delivered' = the surface's net work up to the peak.")
	var out []FieldResult
	for _, c := range cells {
		fr := DefaultFieldRing()
		fr.Cells = c
		st := time.Now()
		fr.Log = func(f string, a ...any) {
			fmt.Fprintf(os.Stderr, "[5b r/%g %s] "+f+"\n", append([]any{c, time.Since(st).Round(time.Second)}, a...)...)
		}
		t0 := time.Now()
		r, err := fr.Run()
		if err != nil {
			return out, err
		}
		r.Seconds = time.Since(t0).Seconds()
		out = append(out, r)
		fmt.Printf("    grid r/%g (dx = %.3g mm, %d nodes, %.0f s):\n", c, r.Dx*1e3, r.Nodes, r.Seconds)
		fmt.Printf("      ring-down: f0 = %.5g GHz, Q = %.4g (energy decay), %.4g (field peaks); ledger residual %.1e of the start\n",
			r.F0/1e9, r.Q, r.QProbe, r.RDResidual)
		fmt.Printf("      implosion: incident %.4g J (2γ·W0·∫g²dt), delivered by the surface %.4f of it; at the peak (%.1f periods)\n",
			r.EInc, r.Delivered/r.EInc, r.PeakT)
		fmt.Printf("        the ring region holds %.4f of the delivered energy and %.2e of it is outside the surface or absorbed\n",
			r.Peak/r.Delivered, r.LeakAtPeak/r.Delivered)
		fmt.Printf("        (reflection + mismatch of the declared wave); explosion afterwards %.4f of the delivered energy;\n",
			(r.LeakEnd-r.LeakAtPeak)/r.Delivered)
		fmt.Printf("        at the end inside %.2e; surface work in %.4f, out %.4f; ledger residual %.1e of the incident energy\n",
			r.InsideEnd/r.Delivered, r.SurfInTotal/r.EInc, r.SurfOutTotal/r.EInc, r.Residual)
		fmt.Printf("      pattern at the peak: E_φ carries %.4f of the electric energy in the tube; a quarter period later B_φ carries %.2e\n",
			r.EPhiFrac, r.BPhiFrac)
		fmt.Printf("        of the magnetic energy there, B_z at the centre %.3g T, at (2a, 0, 0) %.3g T (opposite: the flux closes outside);\n",
			r.BzCentre, r.BzOutside)
		fmt.Printf("        correlation of E with the free ring-down's E at an E maximum: %.4f\n", r.Corr)
		fmt.Printf("      constant burst (2 τ_E, same energy and pattern): peak inside %.4f, outside + absorbed when it ends %.4f; ledger %.1e\n",
			r.BurstPeak/r.BurstEInc, r.BurstLeakDrive/r.BurstEInc, r.BurstResidual)
	}
	e := DefaultERing()
	fmt.Printf("    lumped model (Part 1, lossless): f0 = %.5g GHz, Q_rad = %.4g\n", e.Omega0()/2/math.Pi/1e9, e.QRad())
	return out, nil
}

// Part6 holds the flavour results.
type Part6 struct {
	Persist     []Persist
	Beats       []Beat
	Thetas      []Beat
	Crossed     Beat
	Cluster     Cluster
	ClusterZero Cluster
	Three       ThreeRun
	ThreeRev    ThreeRun
	Explosions  []Explosion
	Rebuilds    []Rebuild
	PortR       float64
	FExplode    Explosion
	FRebuild    Rebuild
	FGrowth     Rebuild
}

func reportPart6() (Part6, error) {
	var p Part6
	t0 := time.Now()
	r := DefaultFRing(false)
	k := r.Kappa(1)
	fmt.Println("\n6. Flavours as wave types (H39, user hypothesis 2026-11-07): constructive, destructive, standing")
	fmt.Println("   The uniform E-ring (TE01δ, azimuthal order 0) is one state; the types live on harmonics n ≥ 1 of the inverted")
	fmt.Println("   ring lines (Part 3 values, 24 sections, 1 pF stray to ground, 1 mΩ shunt loss). z_j = u_j + iωL_g·g_j = A e^{−iκj} + B e^{iκj}.")
	fmt.Println("   Fixed definition, relative to a reference point c: C (constructive) = travelling, B = 0; S (standing) = |A| = |B|,")
	fmt.Println("   antinode at c; D (destructive) = |A| = |B|, node at c.")
	fmt.Printf("   Full cancellation: Σ|z_j|² = N(|A|²+|B|²) (the two waves are orthogonal; smallest singular value √N = %.4f),\n", math.Sqrt(float64(r.N)))
	fmt.Println("   so no non-zero state cancels everywhere; C = (S + iD)/√2, so the three types span two dimensions, not three.")
	ps, err := Persistence(100)
	if err != nil {
		return p, err
	}
	p.Persist = ps
	fmt.Println("   Test 1, persistence on the perfect ring (harmonic 1, c = ½, 100 periods, kronspice transient from the imposed state):")
	fmt.Printf("   %-8s %-4s %12s %12s %12s %10s %10s %9s\n", "ring", "type", "own min", "other max", "purity min", "turns/T", "energy end", "ledger")
	for _, x := range ps {
		fmt.Printf("   %-8s %-4s %12.10f %12.2e %12.10f %10.5f %10.5f %9.1e\n", x.Ring, x.Type, x.MinOwn, x.MaxOther, x.MinPurity, x.Rotation, x.EnergyEnd, x.LedgerRel)
	}
	fmt.Println("   All three are stationary: C rotates once a period (pattern turning, nothing converts), S and D keep their nodes;")
	fmt.Println("   energy falls only by the 1 mΩ shunt loss. Standard prediction: every superposition of a degenerate pair is stationary.")

	fmt.Println("   Test 2, flavour oscillation: one section's C made C(1+δ) (both wires). The pair splits into standing waves with")
	fmt.Println("   node / antinode at the defect; CW → CCW → CW with P = sin²2θ·sin²(πΔf·t):")
	fmt.Printf("   %-7s %6s %6s %12s %12s %12s %8s %8s %9s %8s\n", "ring", "δ", "start", "Δf exact", "Δf 1st ord.", "Δf kronspice", "sin²2θ", "θ fit", "max |res|", "ledger")
	for _, d := range []float64{0.05, 0.1, 0.2, 0.4} {
		b, err := FlavourBeat(false, d, "C", 0, 1.25)
		if err != nil {
			return p, err
		}
		p.Beats = append(p.Beats, b)
	}
	cb, err := FlavourBeat(true, 0.2, "C", 0, 1.25)
	if err != nil {
		return p, err
	}
	p.Crossed = cb
	for _, b := range append(append([]Beat{}, p.Beats...), cb) {
		fmt.Printf("   %-7s %6.2f %6s %9.5g Hz %9.5g Hz %9.5g Hz %8.4f %7.2f° %9.2e %8.1e\n", b.Ring, b.Delta, b.Start, b.DfExact, b.DfFirst, b.DfFit, b.AmpFit, b.ThetaFit, b.MaxResid, b.LedgerRel)
	}
	fmt.Printf("   (f = %.5g MHz; first order Δf = δ·f/(N(1+δ))·4C sin²(κ/2)/(4C sin²(κ/2) + Cg); the plain ring's kronspice Δf is %.2f %% below the exact, trapezoidal frequency warping at 60 steps a period)\n",
		p.Beats[0].F0/1e6, 100*(1-p.Beats[0].DfFit/p.Beats[0].DfExact))
	fmt.Println("   Mixing angle set by where the state is made: start in the standing wave with its antinode s sections from the")
	fmt.Println("   defect's centre (θ = κ·s), δ = 0.2:")
	for _, s := range []float64{1, 2, 3} {
		b, err := FlavourBeat(false, 0.2, "S", s, 1.25)
		if err != nil {
			return p, err
		}
		p.Thetas = append(p.Thetas, b)
		fmt.Printf("     s = %g: θ predicted %.1f°, sin²2θ predicted %.4f; kronspice sin²2θ %.4f (θ %.2f°), Δf %.5g Hz (exact %.5g)\n",
			s, b.ThetaPred, sq(math.Sin(2*b.ThetaPred*math.Pi/180)), b.AmpFit, b.ThetaFit, b.DfFit, b.DfExact)
	}
	fmt.Println("   In the travelling-wave basis θ is 45° for any reciprocal defect (the eigenmodes are real standing waves); a")
	fmt.Println("   different θ needs a basis tied to a production point, or broken time-reversal symmetry (rotation, a gyrator).")

	fmt.Println("   Three states: one harmonic gives two, so three flavours need a third near-degenerate state. The plain ring's band")
	fmt.Println("   edge has one: the π-mode (κ = π, neighbouring nodes in antiphase, a standing wave that cannot travel) lies 0.86 %")
	fmt.Println("   below the pair at κ = 11π/12 (the CL dispersion is flat there). Ring flavours f1, f2 = standing cos / sin at c = ½")
	fmt.Println("   on κ = 11π/12, f3 = the π-mode; two capacitor defects (sizes ≤ 0.5, any two sections) set the mixing.")
	cz := ClusterMix(0, 0, 0, 1)
	p.ClusterZero = cz
	fmt.Printf("   unperturbed: ω = %.6g, %.6g, %.6g rad/s\n", cz.Omega[0], cz.Omega[1], cz.Omega[2])
	cl := FitPMNS()
	p.Cluster = cl
	t := NuFIT6()
	fmt.Printf("   best match to NuFIT 6.0 (normal ordering, with SK-atm): δ = %.4f on section %d, %.4f on section %d; flavours e, μ, τ = f%d, f%d, f%d\n",
		cl.D1, cl.P1, cl.D2, cl.P2, cl.Assign[0]+1, cl.Assign[1]+1, cl.Assign[2]+1)
	fmt.Printf("     θ12 %.2f° (measured %.2f ± %.2f), θ23 %.2f° (%.1f ± %.1f), θ13 %.2f° (%.2f ± %.2f), Δ31/Δ21 %.2f (Δm²31/Δm²21 = %.2f)\n",
		cl.Th12, t.Th12, t.S12, cl.Th23, t.Th23, t.S23, cl.Th13, t.Th13, t.S13, cl.Ratio, t.Ratio)
	fmt.Printf("     cost %.3g (sum of squared deviations in measured σ); largest share of a mode outside the three ring states %.3f\n", cl.Cost, cl.Leak)
	fmt.Printf("     modes ω = %.6g, %.6g, %.6g rad/s; U (rows f1, f2, f3; columns modes by ω):\n", cl.Omega[0], cl.Omega[1], cl.Omega[2])
	for a := 0; a < 3; a++ {
		fmt.Printf("       %+8.4f %+8.4f %+8.4f\n", cl.U[a][0], cl.U[a][1], cl.U[a][2])
	}
	fmt.Println("     U is real (a reciprocal, lossless ring): δ_CP is 0 or 180°, no CP violation is possible; measured δ_CP = 212° (+26/−41),")
	fmt.Println("     consistent with CP conservation at about 1σ in normal ordering (NuFIT 6.0).")
	from := cl.Assign[0]
	three, err := ThreeFlavour(cl, from, 1)
	if err != nil {
		return p, err
	}
	rev, err := ThreeFlavour(cl, cl.Assign[1], 1)
	if err != nil {
		return p, err
	}
	p.Three, p.ThreeRev = three, rev
	// T symmetry: P(e→μ)(t) against P(μ→e)(t).
	var dmax float64
	for i := range three.TAll {
		pr := interp(rev.TAll, rev.PAll[from], three.TAll[i])
		dmax = math.Max(dmax, math.Abs(three.PAll[cl.Assign[1]][i]-pr))
	}
	three.Reverse = dmax
	p.Three.Reverse = dmax
	fmt.Printf("   kronspice, started in e (f%d), %.0f periods (one slow beat); energy in each ring flavour / the whole ring's:\n", from+1, three.Periods)
	fmt.Printf("     against the ring's full mode expansion (all %d modes, same readout, trapezoidal frequencies): largest deviation %.2e\n", r.N, three.MaxDevAll)
	fmt.Printf("     against the 3×3 formula P(α→β) = |Σ_i U_αi U_βi e^{iω_i t}|²: largest deviation %.3f (up to %.1f %% of a mode lies\n", three.MaxDev, 100*cl.Leak)
	fmt.Println("       outside the three states: the defects these angles need mix other harmonics in)")
	fmt.Printf("     ledger %.1e; T symmetry: |P(e→μ)(t) − P(μ→e)(t)| at most %.1e\n", three.LedgerRel, dmax)
	fmt.Println("   Standard prediction: this is the generic mathematics of coupled modes (beats), which the measured neutrino formula")
	fmt.Println("   shares; the neutrino's content is that flavour is defined by the weak interaction (the charged lepton made), the")
	fmt.Println("   splittings are mass differences (Δm²L/4E), and matter changes them (MSW). The ring reproduces the beat form; its")
	fmt.Println("   flavours are positions on a ring, its splittings are defect sizes, and its numbers match only as far as tuned.")

	fmt.Println("   Test 3, implosion and explosion of each type through a port (the ring's radiation port, R set for Q_ext = 40.5 at an")
	pr := r.PortR(k, 40.5)
	p.PortR = pr
	fmt.Printf("   antinode: R = %.1f Ω, differential, at node 0; a second port at node 6 = a quarter wavelength on). Reference c = node 0.\n", pr)
	conjD := func(f Fractions) float64 { return f.D }
	conjCCW := func(f Fractions) float64 { return f.CCW }
	for _, ty := range []string{"S", "D", "C"} {
		ex, _, err := r.Explode(ty, []Comp{r.Type(ty, k, 0)}, []int{0}, pr, 80, 0, conjD)
		if err != nil {
			return p, err
		}
		p.Explosions = append(p.Explosions, ex)
		fmt.Printf("     free explosion of %s through port 0 (80 periods): energy left %.4f; end state S %.4f, D %.4f, CW %.4f; ledger %.1e\n",
			ty, ex.Remaining, ex.End.S, ex.End.D, ex.End.CW, ex.LedgerRel)
	}
	fmt.Println("     S leaves completely, D is dark to the port and stays, C loses its S half and is left as D: the explosion converts")
	fmt.Println("     a travelling wave into a standing wave with its node at the port.")
	type rc struct {
		name   string
		target Comp
		nodes  []int
		conj   func(Fractions) float64
	}
	for _, c := range []rc{
		{"S, one port", r.Type("S", k, 0), []int{0}, conjD},
		{"C, ports 0 and 6", r.Type("C", k, 0), []int{0, 6}, conjCCW},
		{"C, port 0 only", r.Type("C", k, 0), []int{0}, conjCCW},
	} {
		rb, err := r.RebuildByReversal(c.name, []Comp{c.target}, c.nodes, pr, 80, 40, 0, c.conj)
		if err != nil {
			return p, err
		}
		p.Rebuilds = append(p.Rebuilds, rb)
		fmt.Printf("     rebuild %-17s: absorbed %.5f, ring at the peak %.4f of the incident energy, fidelity %.5f (S %.4f D %.4f CW %.4f),\n",
			c.name, rb.Absorbed, rb.StoredPeak, rb.Fidelity, rb.PeakFrac.S, rb.PeakFrac.D, rb.PeakFrac.CW)
		fmt.Printf("       explosion %.4f; recording left %.1e in the ring; ledger %.1e, time-stepping %.1e\n", rb.EOutAfter, rb.RecRemaining, rb.LedgerRel, rb.Mismatch)
	}
	fmt.Println("     (rebuild = the recorded explosion of the time-reversed target, played back reversed as the incident wave on each port)")
	fmt.Println("     D cannot be built through port 0 at all (its explosion there is zero); C needs both ports, in quadrature.")
	rd := r
	rd.Def = map[int]float64{0: 0.2}
	F := []Comp{rd.Type("S", k, 2)}
	fex, _, err := rd.Explode("F", F, []int{2}, pr, 120, 2, conjD)
	if err != nil {
		return p, err
	}
	p.FExplode = fex
	var mx float64
	for _, v := range fex.PConj {
		mx = math.Max(mx, v)
	}
	fmt.Printf("     flavour superposition: defect δ = 0.2 on section 0, port at node 2, F = standing wave with antinode at the port\n")
	fmt.Printf("       (θ = 22.5° to the defect's modes). Free explosion (120 periods): energy left %.4f, conversion to F⊥ reaches %.4f\n", fex.Remaining, mx)
	fmt.Printf("       while it decays (the beat goes on through the explosion: the dark part beats into F and leaves).\n")
	frb, err := rd.RebuildByReversal("F", F, []int{2}, pr, 120, 60, 2, conjD)
	if err != nil {
		return p, err
	}
	fg, err := rd.RebuildByGrowth("F single tone", F, 2, pr, rd.OmegaK(k)/(2*40.5), 80, 60, 2, conjD)
	if err != nil {
		return p, err
	}
	p.FRebuild, p.FGrowth = frb, fg
	fmt.Printf("       rebuild by the reversed explosion: absorbed %.5f, ring %.4f, fidelity %.5f (F %.4f, F⊥ %.4f); the recording left %.4f\n",
		frb.Absorbed, frb.StoredPeak, frb.Fidelity, frb.PeakFrac.S, frb.PeakFrac.D, frb.RecRemaining)
	fmt.Printf("       rebuild by a single-tone growth at the port: absorbed %.5f, ring %.4f, fidelity %.5f (F %.4f, F⊥ %.4f)\n",
		fg.Absorbed, fg.StoredPeak, fg.Fidelity, fg.PeakFrac.S, fg.PeakFrac.D)
	fmt.Println("       To rebuild a flavour state the incoming wave must carry both mode frequencies with the right amplitudes and phases:")
	fmt.Println("       the time-reversed beat. After the peak the beat resumes as in the free explosion.")
	fmt.Printf("   (Part 6: %s)\n", time.Since(t0).Round(time.Second))
	return p, nil
}
