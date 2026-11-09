// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command rodin_resonance builds the Rodin coil of Ramírez Galindo et al.
// (Rev. Mex. Fís. 70, 030901, 2024) as a counter-wound bifilar pair on its
// toroidal former, drives the two windings with a sine and a cosine, maps
// the magnetic field they make, builds it up by resonance, and adds groups
// of output coils round the toroid that charge a capacitor bank through
// rectifiers, the drive frequency and the groups evolved together. Every
// joule is booked against the two drive sources (register entry H43).
//
//	Part 1  geometry and validation: the {12/5} star path, the centre field
//	        against the paper's Table II, inductances by the segment method,
//	        reciprocity, flux and field-energy checks.
//	Part 2  quadrature drive: A only, B only, A + B in phase, A + B in
//	        quadrature, in both connection senses; field maps, B(t),
//	        magnetic energy, the power the sources pass to each other; the
//	        tuned (resonant) build-up and its Q.
//	Part 2b the winding as a line on a knot: self-resonances, and the
//	        crossed (series, counter-wound) ring against the plain one for
//	        the (5, 12) Rodin path and a (2, 1) winding on the same former —
//	        the 720° modes of examples/photon_electron (H23).
//	Part 3  output groups: power against (frequency, angle), the genetic
//	        search over frequency and groups, the winner in a kronspice
//	        transient with Schottky bridges and the bank, and the ledger.
//
//	go run ./examples/rodin_resonance            # full run (several minutes)
//	go run ./examples/rodin_resonance -short     # quick, coarse
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
)

// TableII is the paper's centre field (current A; theory, FEM, measured mT).
var TableII = [][4]float64{{1.60, 3.90, 4.08, 4.13}, {1.25, 3.01, 3.00, 3.05}, {0.60, 1.30, 1.10, 1.16}}

var (
	short  = flag.Bool("short", false, "coarse grids and a small GA")
	vizDir = flag.String("viz", filepath.Join("examples", "rodin_resonance", "viz"), "directory for the plot data")
	part   = flag.String("part", "all", "all, 1, 2, 2b, 3 or 3t")
)

func main() {
	flag.Parse()
	if err := os.MkdirAll(*vizDir, 0o755); err != nil {
		log.Fatal(err)
	}
	start := time.Now()
	p := NewPair(Rodin, Circuits)
	m := NewModel(p)
	if *part == "all" || *part == "1" {
		part1(m)
	}
	var lowV float64
	if *part == "all" || *part == "2" {
		lowV = part2(m)
	}
	if lowV == 0 {
		// Part 2's drive amplitude: A alone carries 1.25 A peak at 1 kHz
		lowV = 1.25 * cmplx.Abs(complex(SourceOhms+m.RDC, 2*math.Pi*1e3*m.LA))
	}
	if *part == "all" || *part == "2b" {
		part2b()
	}
	if *part == "all" || *part == "3" {
		part3(m, lowV)
	}
	if *part == "all" || *part == "3t" {
		part3t(m, lowV)
	}
	fmt.Printf("\n(run time %.0f s)\n", time.Since(start).Seconds())
}

// FitCircuits returns the number of circuits per winding that fits the
// paper's measured centre field best (least squares over Table II, the
// pair in the dipole sense, the paper's series connection), with each
// candidate's field per ampere.
func FitCircuits(lo, hi int) (best int, perAmp map[int]float64) {
	perAmp = map[int]float64{}
	bestSS := math.Inf(1)
	for n := lo; n <= hi; n++ {
		q := NewPair(Rodin, n)
		b := sub(FieldOf(q.A, float64(n), P3{}), FieldOf(q.B, float64(n), P3{}))[2]
		perAmp[n] = b
		ss := 0.0
		for _, r := range TableII {
			d := b*r[0]*1e3 - r[3]
			ss += d * d
		}
		if ss < bestSS {
			best, bestSS = n, ss
		}
	}
	return
}

func part1(m *Model) {
	p := m.P
	fmt.Println("Part 1 — the Rodin coil: geometry and validation")
	fmt.Printf("  former (Ramírez Galindo et al. 2024): outer diameter %.2f cm, aperture %.2f cm; torus R0 = %.4f cm, a = %.4f cm\n",
		2*OuterRadius*100, 2*HoleRadius*100, TorusR0*100, TorusA*100)
	fmt.Printf("  path: %d poles 30° apart, %d poles (150°) per pass: 1→6→11→4→9→2→7→12→5→10→3→8→1; one circuit goes 5 times round the axis and 12 times through the hole\n", Poles, Skip)
	fmt.Printf("  one circuit of A: %.4f m of wire; of B (mirror path, one bundle out): %.4f m\n", polyLen(p.A), polyLen(p.B))
	best, per := FitCircuits(6, 16)
	fmt.Println("  circuits per winding fitted to the measured centre field (A and B in series, dipole sense):")
	for n := 8; n <= 14; n++ {
		mark := ""
		if n == best {
			mark = "  ← least squares"
		}
		fmt.Printf("    Nc = %2d: %.4f mT/A%s\n", n, per[n]*1e3, mark)
	}
	if best != Circuits {
		fmt.Printf("  WARNING: fit gives %d, the code uses %d\n", best, Circuits)
	}
	fmt.Printf("  wire 0.6 mm copper, %d circuits per winding in a round bundle of radius %.3f mm; A to B bundle pitch %.3f mm\n", p.Nc, p.Radius*1e3, p.Pitch*1e3)
	fmt.Println("  centre field, Biot–Savart (em.SegmentB) against Table II [mT]:")
	fmt.Printf("    %6s %10s %10s %10s %12s %14s %14s\n", "I [A]", "theory", "FEM", "measured", "kronspice", "vs measured", "vs theory")
	bc := sub(m.FieldA(P3{}), m.FieldB(P3{}))[2]
	for _, r := range TableII {
		b := bc * r[0] * 1e3
		fmt.Printf("    %6.2f %10.2f %10.2f %10.2f %12.3f %13.1f%% %13.1f%%\n", r[0], r[1], r[2], r[3], b, 100*(b/r[3]-1), 100*(b/r[1]-1))
	}
	fmt.Println("    (Biot–Savart is linear in I; the paper's own numbers are not: measured 2.58, 2.44, 1.93 mT/A at 1.60, 1.25, 0.60 A)")
	ring := em.Circle([3]float64{}, [3]float64{0, 0, 1}, 0.05, 200)
	pr := make([]P3, len(ring))
	copy(pr, ring)
	ls, lk := LoopInductance(pr, 3e-4), CircleRing(0.05, 3e-4)
	fmt.Printf("  segment method on a 5 cm ring of 0.3 mm wire: %.6g H against Kirchhoff's μ0R(ln(8R/r) − 7/4) = %.6g H (%.1e)\n", ls, lk, ls/lk-1)
	fmt.Printf("  inductances: L_A = %.5g µH, L_B = %.5g µH, M = %+.5g µH (path-sense currents), k = %+.4f\n", m.LA*1e6, m.LB*1e6, m.M*1e6, m.K())
	fmt.Printf("    in the dipole sense (B reversed) M = %+.5g µH; series inductance toroidal %.5g µH, dipole %.5g µH\n",
		-m.M*1e6, (m.LA+m.LB+2*m.M)*1e6, (m.LA+m.LB-2*m.M)*1e6)
	fmt.Printf("  wire resistance (DC): A %.4f Ω (%.3f m), B %.4f Ω (%.3f m)\n", WireResistance(m.LenA, 0), m.LenA, WireResistance(m.LenB, 0), m.LenB)
	mab, mba := m.Reciprocity()
	fmt.Printf("  reciprocity: M_AB = %.10g µH, M_BA = %.10g µH (relative difference %.1e)\n", mab*1e6, mba*1e6, math.Abs(mab-mba)/math.Abs(mab))
	n1, f1 := m.FluxCheck(P3{}, P3{0, 0, 1}, 0.01)
	fmt.Printf("  flux check, 1 cm circle in the aperture: Neumann %.6g nH, ∫B·dS of A's field %.6g nH (%.1e)\n", n1*1e9, f1*1e9, f1/n1-1)
	th := 15 * math.Pi / 180
	n2, f2 := m.FluxCheck(P3{TorusR0 * math.Cos(th), TorusR0 * math.Sin(th), 0}, P3{-math.Sin(th), math.Cos(th), 0}, p.OuterShell()+3e-3)
	fmt.Printf("  flux check, circle round the tube at 15°: Neumann %.6g nH, ∫B·dS %.6g nH (%.1e)\n", n2*1e9, f2*1e9, f2/n2-1)
	h := 0.75e-3
	if *short {
		h = 2e-3
	}
	q := NewPair(Knot21, 1)
	for _, c := range []struct {
		name string
		path []P3
	}{{"(2, 1) knot, one circuit", q.A}, {"Rodin path, one circuit", p.A}} {
		w0, w1 := EnergyCheck(c.path, 2.5e-3, h)
		fmt.Printf("  energy check, %s at 1 A, outside a 2.5 mm tube: ½L_ext·I² = %.5g nJ, ∫B²/2μ0 dV = %.5g nJ (%+.2f %%, grid %.2g mm)\n",
			c.name, w0*1e9, w1*1e9, 100*(w1/w0-1), h*1e3)
	}
	fmt.Println("    (the Rodin path crosses itself on the torus surface: the tubes overlap at the crossings, so the grid misses part of their volume)")
	writePaths(m)
}

func writePaths(m *Model) {
	rows := [][]string{{"winding", "x", "y", "z"}}
	add := func(name string, pl []P3) {
		for _, pt := range pl {
			rows = append(rows, []string{name, f6(pt[0]), f6(pt[1]), f6(pt[2])})
		}
	}
	add("A", m.P.A)
	add("B", m.P.B)
	q := NewPair(Knot21, 1)
	add("K21A", q.A)
	add("K21B", q.B)
	writeCSV("paths.csv", rows)
}

func f6(x float64) string { return strconv.FormatFloat(x, 'g', 7, 64) }

func writeCSV(name string, rows [][]string) {
	f, err := os.Create(filepath.Join(*vizDir, name))
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := csv.NewWriter(f)
	w.WriteAll(rows)
	if err := w.Error(); err != nil {
		log.Fatal(err)
	}
}

// part2 runs the drive cases and returns the low-frequency amplitude.
func part2(m *Model) float64 {
	fmt.Println("\nPart 2 — sine on A, cosine on B: the field they make, and the power they pass")
	f := 1e3
	z := cmplx.Abs(complex(SourceOhms+m.RDC, 2*math.Pi*f*m.LA))
	v := 1.25 * z
	cyc, steps := 20, 200
	if *short {
		cyc, steps = 12, 100
	}
	fmt.Printf("  lumped two-winding circuit in kronspice (valid far below the first self-resonance, Part 2b): sources %.4g V peak at %.0f Hz through %.0f Ω,\n", v, f, SourceOhms)
	fmt.Println("  chosen so that A alone carries 1.25 A peak (a current of the paper's Table II)")
	type c struct {
		name   string
		va, vb float64
		ph     float64
		sense  int
	}
	cases := []c{
		{"A only", v, 0, 0, Toroidal},
		{"B only", 0, v, 0, Toroidal},
		{"A + B in phase, toroidal sense", v, v, 0, Toroidal},
		{"A + B in phase, dipole sense", v, v, 0, Dipole},
		{"A + B quadrature (sin/cos), toroidal sense", v, v, 90, Toroidal},
		{"A + B quadrature (sin/cos), dipole sense", v, v, 90, Dipole},
	}
	n := 81
	if *short {
		n = 41
	}
	planes := []*Plane{m.NewPlane("xy", n, 0.07), m.NewPlane("xz", n, 0.07)}
	writePlanes(planes)
	pts := []struct {
		name string
		p    P3
	}{{"centre (0, 0, 0)", P3{}}, {"in the tube at a pole (R0, 0°)", P3{TorusR0, 0, 0}},
		{"in the tube between poles (R0, 15°)", P3{TorusR0 * math.Cos(math.Pi/12), TorusR0 * math.Sin(math.Pi/12), 0}},
		{"on the axis 3 cm up", P3{0, 0, 0.03}}}
	ptRows := [][]string{{"point", "x", "y", "z", "bAx", "bAy", "bAz", "bBx", "bBy", "bBz"}}
	pfa := make([]P3, len(pts))
	pfb := make([]P3, len(pts))
	for i, q := range pts {
		pfa[i], pfb[i] = m.FieldA(q.p), m.FieldB(q.p)
		ptRows = append(ptRows, []string{q.name, f6(q.p[0]), f6(q.p[1]), f6(q.p[2]), f6(pfa[i][0]), f6(pfa[i][1]), f6(pfa[i][2]), f6(pfb[i][0]), f6(pfb[i][1]), f6(pfb[i][2])})
	}
	writeCSV("points.csv", ptRows)
	drv := [][]string{{"case", "f", "sense", "IA_re", "IA_im", "IB_re", "IB_im"}}
	wm := [][]string{{"case", "t", "Wm"}}
	fmt.Printf("  %-44s %8s %8s %9s %9s %10s %10s %11s %11s %9s\n", "case", "|IA| A", "|IB| A", "P_A W", "P_B W", "A→B W", "ωMI_AI_Bsinφ", "Wm min µJ", "Wm max µJ", "ledger")
	type res struct {
		r *DriveRun
		c c
	}
	var all []res
	for _, cs := range cases {
		d := Drive{Name: cs.name, F: f, VA: cs.va, VB: cs.vb, PhaseB: cs.ph, Sense: cs.sense, Cycles: cyc, StepsPerCyc: steps}
		r, err := m.Run(d)
		if err != nil {
			log.Fatal(cs.name, ": ", err)
		}
		all = append(all, res{r, cs})
		fmt.Printf("  %-44s %8.4f %8.4f %9.4f %9.4f %+10.4f %+10.4f %11.4g %11.4g %9.1e\n", cs.name, cmplx.Abs(r.IA), cmplx.Abs(r.IB), r.PA, r.PB, r.Exchange, r.ExchangePhasor,
			r.WmMin*1e6, r.WmMax*1e6, r.Mismatch/r.EIn)
		drv = append(drv, []string{cs.name, f6(f), strconv.Itoa(cs.sense), f6(real(r.IA)), f6(imag(r.IA)), f6(real(r.IB)), f6(imag(r.IB))})
		for k := 0; k < len(r.Time); k += 4 {
			wm = append(wm, []string{cs.name, f6(r.Time[k]), f6(r.Wm[k])})
		}
	}
	fmt.Println("  ledger column: integration mismatch / energy in (sources in = copper + source resistors + Δ stored); per-source books over the last cycle:")
	for _, a := range all {
		r := a.r
		fmt.Printf("    %-44s sources %.4f W = copper %.4f W + source resistors %.4f W (difference %.1e W)\n", a.c.name, r.PA+r.PB, r.Copper, r.SrcR, r.PA+r.PB-r.Copper-r.SrcR)
	}
	fmt.Println("  field (|B| peak over a cycle, rms, ellipticity = minor/major axis of the ellipse B(t) traces; 0 = it only swings, 1 = it turns at constant size):")
	for _, a := range all {
		r := a.r
		fmt.Printf("    %s\n", a.c.name)
		for i, q := range pts {
			s := pointStats(r.IA, r.IB, pfa[i], pfb[i])
			fmt.Printf("      %-38s peak %7.4f mT, rms %7.4f mT, ellipticity %.3f\n", q.name, s.Peak*1e3, s.RMS*1e3, s.Ellipticity)
		}
		for _, pl := range planes {
			pk, at, en, el := pl.Stats(r.IA, r.IB)
			fmt.Printf("      plane %s (±7 cm, wires ±3 mm excluded): largest peak %.4f mT at (%.1f, %.1f, %.1f) mm; mean ⟨B²⟩/2μ0 %.4g J/m³; ⟨B²⟩-weighted ellipticity %.3f\n",
				pl.Name, pk*1e3, at[0]*1e3, at[1]*1e3, at[2]*1e3, en, el)
		}
		adv, tr, ar := m.RingTravel(r.IA, r.IB, 360, 96)
		fmt.Printf("      12-fold pattern of |B|² round the ring ρ = R0, z = 0: net phase advance over one cycle %.2g° (travel %.2g°), min/max amplitude %.3g\n", adv, tr, ar)
	}
	writeCSV("drives.csv", drv)
	writeCSV("wm_lowf.csv", wm)
	fmt.Println("  the quadrature pair: the toroidal and the dipole parts of the field are a quarter cycle apart, so the field vector turns at each point")
	fmt.Println("  (ellipticity > 0) while the pattern itself stands still round the axis (12-fold symmetric at every instant: zero net advance).")

	// Resonant build-up
	f0 := 1e5
	ra := WireResistance(m.LenA, f0)
	v0 := 1.25 * (SourceOhms + ra)
	cyc0 := 300
	if *short {
		cyc0 = 120
	}
	fmt.Printf("\n  resonance at f0 = %.0f kHz: series capacitors tuned to each winding's own L (C_A = %.4g nF, C_B = %.4g nF); sources %.4g V peak (1.25 A in A alone when tuned)\n",
		f0/1e3, 1e9/(math.Pow(2*math.Pi*f0, 2)*m.LA), 1e9/(math.Pow(2*math.Pi*f0, 2)*m.LB), v0)
	rcases := []Drive{
		{Name: "A only, tuned", F: f0, VA: v0, Tuned: true},
		{Name: "A only, not tuned", F: f0, VA: v0},
		{Name: "quadrature, toroidal sense, tuned", F: f0, VA: v0, VB: v0, PhaseB: 90, Sense: Toroidal, Tuned: true},
		{Name: "quadrature, dipole sense, tuned", F: f0, VA: v0, VB: v0, PhaseB: 90, Sense: Dipole, Tuned: true},
		{Name: "quadrature, toroidal sense, not tuned", F: f0, VA: v0, VB: v0, PhaseB: 90, Sense: Toroidal},
		{Name: "in phase, toroidal, tuned to L + M", F: f0, VA: v0, VB: v0, Sense: Toroidal, Tuned: true, TuneLA: m.LA + m.M, TuneLB: m.LB + m.M},
		{Name: "quadrature, toroidal, tuned to L + M", F: f0, VA: v0, VB: v0, PhaseB: 90, Sense: Toroidal, Tuned: true, TuneLA: m.LA + m.M, TuneLB: m.LB + m.M},
	}
	ru := [][]string{{"case", "t", "Wm"}}
	fmt.Printf("  %-40s %9s %9s %9s %11s %10s %9s %9s %9s %10s\n", "case", "|IA| A", "|IB| A", "φB−φA °", "Wm mean µJ", "P_in W", "τ µs", "Q=πf0τ", "Q=ω0L/R", "ledger")
	for _, d := range rcases {
		d.Cycles, d.StepsPerCyc = cyc0, 100
		r, err := m.Run(d)
		if err != nil {
			log.Fatal(d.Name, ": ", err)
		}
		tau, qc := math.NaN(), math.NaN()
		if d.Tuned {
			tau = riseTime(r, f0)
		}
		if d.Tuned && d.VB == 0 {
			qc = 2 * math.Pi * f0 * m.LA / (SourceOhms + ra)
		}
		dph := math.NaN()
		if d.VB != 0 {
			dph = math.Remainder(cmplx.Phase(r.IB)-cmplx.Phase(r.IA), 2*math.Pi) * 180 / math.Pi
		}
		fmt.Printf("  %-40s %9.4f %9.4f %9.1f %11.4g %10.4f %9.3g %9.3g %9.3g %10.1e\n", d.Name, cmplx.Abs(r.IA), cmplx.Abs(r.IB), dph, r.WmMean*1e6, r.PA+r.PB,
			tau*1e6, math.Pi*f0*tau, qc, r.Mismatch/r.EIn)
		// the envelope: largest magnetic energy in each cycle, at the cycle's centre
		T0 := 1 / f0
		env := make([]float64, d.Cycles)
		for k, t := range r.Time {
			if c := int(t / T0); c >= 0 && c < d.Cycles {
				env[c] = math.Max(env[c], r.Wm[k])
			}
		}
		for c, w := range env {
			ru = append(ru, []string{d.Name, f6((float64(c) + 0.5) * T0), f6(w)})
		}
	}
	writeCSV("ringup.csv", ru)
	fmt.Println("  τ: time for the current amplitude (√ of the cycle-maximum magnetic energy) to reach 1 − 1/e of its final value; a series")
	fmt.Println("  tank switched on at resonance grows as 1 − e^(−t/τ) with τ = 2L/R, so Q = π·f0·τ. The build-up is the resonance alone: the stored energy rises to")
	fmt.Println("  the level where the loss per cycle equals what the sources put in; the untuned windings at the same voltage hold far less.")
	return v
}

// riseTime is the amplitude time constant of the build-up: the time for
// the current amplitude, √(cycle-maximum magnetic energy), to reach
// (1 − 1/e) of its value over the last cycle. A series tank switched on at
// resonance from a voltage source grows as 1 − e^(−t/τ), τ = 2L/R, so
// Q = ω0·L/R = π·f0·τ.
func riseTime(r *DriveRun, f float64) float64 {
	T := 1 / f
	n := int(r.Time[len(r.Time)-1]/T + 0.5)
	env := make([]float64, n)
	for k, t := range r.Time {
		c := int(t / T)
		if c >= 0 && c < n {
			env[c] = math.Max(env[c], math.Sqrt(math.Max(0, r.Wm[k])))
		}
	}
	final := env[n-1]
	for c := 0; c < n; c++ {
		if env[c] >= (1-1/math.E)*final {
			// interpolate within the cycle
			if c == 0 {
				return T
			}
			fr := ((1-1/math.E)*final - env[c-1]) / (env[c] - env[c-1])
			return (float64(c) - 0.5 + fr) * T // envelope values sit at the cycle centres
		}
	}
	return math.NaN()
}

func writePlanes(planes []*Plane) {
	rows := [][]string{{"plane", "x", "y", "z", "near", "bAx", "bAy", "bAz", "bBx", "bBy", "bBz"}}
	for _, pl := range planes {
		for k, q := range pl.Pts {
			nr := "0"
			if pl.Near[k] {
				nr = "1"
			}
			a, b := pl.BA[k], pl.BB[k]
			rows = append(rows, []string{pl.Name, f6(q[0]), f6(q[1]), f6(q[2]), nr, f6(a[0]), f6(a[1]), f6(a[2]), f6(b[0]), f6(b[1]), f6(b[2])})
		}
	}
	writeCSV("field_planes.csv", rows)
}

func part2b() {
	fmt.Println("\nPart 2b — the winding as a line on a knot: self-resonances and the 720° modes")
	fmt.Println("  sectional model: one section per pass of every circuit (12·Nc per winding), inductance matrix from the geometry, declared")
	fmt.Printf("  capacitances (wire to wire in a bundle %.4g pF/m, A–B contact ×%.1f, %.3g pF to the surroundings in all); lossless eigenmodes\n", WireC()*1e12, ABContact, 4*math.Pi*em.Eps0*OuterRadius*1e12)
	rows := [][]string{{"kind", "topology", "mode", "f", "parity"}}
	for _, kind := range []Kind{Rodin, Knot21} {
		p := NewPair(kind, Circuits)
		m := NewModel(p)
		fmt.Printf("  %s winding, %d circuits per winding: %.3f m of wire per winding, L_A = %.4g µH, k = %+.3f\n", kind, p.Nc, m.LenA, m.LA*1e6, m.K())
		for _, top := range []Topology{OpenA, SeriesTor, SeriesDip, PlainRing, CrossedRing} {
			n := m.NewNet(top, 0, 1e6)
			f, v, err := n.Modes()
			if err != nil {
				log.Fatal(err)
			}
			var list []string
			cnt := 0
			for i := 0; i < len(f) && cnt < 6; i++ {
				if f[i] < 1e3 { // the zero mode (all nodes at one potential)
					continue
				}
				par := math.NaN()
				s := fmt.Sprintf("%.4g MHz", f[i]/1e6)
				if top == CrossedRing || top == PlainRing {
					par = n.TripParity(v, i)
					s += fmt.Sprintf(" (%+.2f)", par)
				}
				list = append(list, s)
				rows = append(rows, []string{kind.String(), top.String(), strconv.Itoa(cnt + 1), f6(f[i]), f6(par)})
				cnt++
			}
			fmt.Printf("    %-24s %v\n", top.String()+":", list)
		}
	}
	writeCSV("modes.csv", rows)
	fmt.Println("  (parity: correlation of the node voltages along A with those along B; for the crossed ring −1 means the wave comes back")
	fmt.Println("  inverted after one trip through one winding, so it needs two trips (720°) to come back the same way up, as on")
	fmt.Println("  photon_electron's crossed ring line; for the plain rings it says whether the two closed windings swing together or opposite)")
}

func part3(m *Model, v float64) {
	fmt.Println("\nPart 3 — groups of output coils round the toroid, the drive frequency evolved with them")
	if v > 0 {
		Part3Volts = v
	}
	fmt.Printf("  drive: sine on A, cosine on B, %.4g V peak each through %.0f Ω; groups tuned to the drive frequency by series capacitors\n", Part3Volts, SourceOhms)
	// 1. the map, on the full sectional net (12 sections per circuit)
	nf, nth := 90, 180
	if *short {
		nf, nth = 30, 72
	}
	var freqs, thetas []float64
	for i := 0; i < nf; i++ {
		freqs = append(freqs, 1e3*math.Pow(3e7/1e3, float64(i)/float64(nf-1)))
	}
	for k := 0; k < nth; k++ {
		thetas = append(thetas, 2*math.Pi*float64(k)/float64(nth))
	}
	shapes := []struct {
		name string
		s    GroupShape
	}{
		{"tube", GroupShape{Rho: TorusR0, Z: 0, Beta: 0, Radius: m.P.OuterShell() + 3e-3, Turns: 20}},
		{"top", GroupShape{Rho: TorusR0, Z: m.P.OuterShell() + 4e-3, Beta: math.Pi / 2, Radius: 0.01, Turns: 20}},
	}
	open := m.NewNet(OpenA, 0, 1e6)
	fo, _, _ := open.Modes()
	mr := NewModel(NewPairSec(Rodin, Circuits, 2))
	fmt.Printf("  maps, %d frequencies × %d angles, for a reference coil of 20 turns at each angle θ:\n", nf, nth)
	fmt.Println("    (a) its open-circuit EMF, from the full sectional net (12 sections per circuit); exact, since an open coil draws no current;")
	fmt.Println("    (b) the power into its load when tuned and loaded by R_L = its wire resistance, the coil in the network (it loads the drive),")
	fmt.Println("        on the reduced net (2 sections per circuit) that the genetic search and the kronspice transient use")
	header := func() []string {
		h := []string{"f"}
		for _, t := range thetas {
			h = append(h, f6(t*180/math.Pi))
		}
		return h
	}
	dump := func(name string, mp [][]float64) {
		rows := [][]string{header()}
		for i := range freqs {
			row := []string{f6(freqs[i])}
			for k := range thetas {
				row = append(row, f6(mp[i][k]))
			}
			rows = append(rows, row)
		}
		writeCSV(name, rows)
	}
	argmax := func(mp [][]float64) (int, int, float64) {
		bi, bk, bv := 0, 0, 0.0
		for i := range mp {
			for k, v := range mp[i] {
				if v > bv {
					bi, bk, bv = i, k, v
				}
			}
		}
		return bi, bk, bv
	}
	for _, sh := range shapes {
		fmt.Printf("  reference coil %q: radius %.1f mm, β = %.0f°, at ρ = %.1f mm, z = %.1f mm (clearance beyond the required 2 mm: %.2f mm)\n", sh.name, sh.s.Radius*1e3, sh.s.Beta*180/math.Pi, sh.s.Rho*1e3, sh.s.Z*1e3, m.Clear(sh.s, 0.1)*1e3)
		for _, sense := range []int{Toroidal, Dipole} {
			sname := map[int]string{Toroidal: "toroidal", Dipole: "dipole"}[sense]
			emf, err := m.EMFMap(m.NewNet(DrivenPair, sense, 1e3), sh.s, freqs, thetas, sense)
			if err != nil {
				log.Fatal(err)
			}
			pl, ps, err := mr.PowerMap(mr.NewNet(DrivenPair, sense, 1e3), sh.s, freqs, thetas, sense)
			if err != nil {
				log.Fatal(err)
			}
			dump(fmt.Sprintf("map_emf_%s_%s.csv", sh.name, sname), emf)
			dump(fmt.Sprintf("map_power_%s_%s.csv", sh.name, sname), pl)
			ei, ek, ev := argmax(emf)
			pi, pk, pv := argmax(pl)
			fmt.Printf("    %s sense: largest EMF %.4g V at f = %.4g Hz, θ = %.0f°; largest load power %.4g W at f = %.4g Hz, θ = %.0f° (sources give %.4g W there: %.3f)\n",
				sname, ev, freqs[ei], thetas[ek]*180/math.Pi, pv, freqs[pi], thetas[pk]*180/math.Pi, ps[pi][pk], pv/ps[pi][pk])
			for _, ft := range []float64{1e3, 1e4, 1e5, 1e6, fo[0], 5e6, 1e7, 2e7} {
				i := nearest(freqs, ft)
				lo, hi := math.Inf(1), 0.0
				for _, x := range emf[i] {
					lo, hi = math.Min(lo, x), math.Max(hi, x)
				}
				hm, ha := dominantHarmonic(emf[i])
				plo, phi := math.Inf(1), 0.0
				for _, x := range pl[i] {
					plo, phi = math.Min(plo, x), math.Max(phi, x)
				}
				fmt.Printf("      f = %9.4g Hz: EMF over θ %.4g … %.4g V (max/min %.3g), strongest angular harmonic m = %d (%.2g of the mean); load power %.4g … %.4g W\n",
					freqs[i], lo, hi, hi/lo, hm, ha, plo, phi)
			}
		}
	}
	fmt.Printf("  (first self-resonance of one winding, its far end open: %.4g MHz; below it every pass carries the same current and the pattern\n", fo[0]/1e6)
	fmt.Println("  repeats every 30° (harmonics m = 12, 24 …); above it the current changes along the wire, and since the wire goes")
	fmt.Println("  round the axis 5 times per circuit and 150° per pass, the field round the toroid takes the harmonics of that standing wave)")

	// 2. the GA, on the reduced sectional net (2 sections per circuit), which the kronspice transient can hold
	p3 := NewPart3(mr)
	pop, gens := 64, 50
	if *short {
		pop, gens = 16, 5
	}
	fmt.Println("\n  model for the searches: the sectional net with 2 sections per circuit (48 + the groups ≤ 64 windings, which the kronspice N element holds)")
	search(mr, p3, "free search", "", GAOptions{Population: pop, Generations: gens, Seed: 43})
	search(mr, p3, "distributed groups (at least 2 coils, radius ≤ 20 mm)", "_dist", GAOptions{Population: pop, Generations: gens, Seed: 44, MinGroups: 2, MaxRadius: 0.02})
}

// search runs one genetic search, reports the winner and verifies it in a
// kronspice transient; tag suffixes the CSV files.
func search(mr *Model, p3 *Part3, title, tag string, o GAOptions) {
	space := Part3Space()
	if o.MinGroups > 1 || o.MaxRadius > 0 {
		space = Part3SpaceLimited(max(1, o.MinGroups), o.MaxRadius)
	}
	fmt.Printf("\n  genetic search, %s (evo, seed %d, %d × %d): genes", title, o.Seed, o.Population, o.Generations)
	for _, g := range space {
		fmt.Printf(" | %s ∈ [%g, %g]", g.Name, g.Min, g.Max)
	}
	fmt.Println()
	ga, err := p3.Evolve(o)
	if err != nil {
		log.Fatal(err)
	}
	hist := [][]string{{"generation", "best", "mean", "failed", "evaluations"}}
	for _, h := range ga.History {
		hist = append(hist, []string{strconv.Itoa(h.Generation), f6(h.Best), f6(h.Mean), strconv.Itoa(h.Failed), strconv.Itoa(h.Evaluations)})
	}
	writeCSV("ga_history"+tag+".csv", hist)
	first := ga.History[0].Best
	fmt.Printf("  fitness P_load × bridge factor: best of generation 0 %.4g W; evolved after %d evaluations %.4g W (P_load %.4g W, bridge factor %.3f)\n",
		math.Pow(10, 6-first), ga.Evals, ga.Power*BridgeFactor(ga.AC), ga.Power, BridgeFactor(ga.AC))
	l := ga.Best
	fmt.Printf("  winner: %v\n", l)
	fmt.Printf("    coil L = %.4g µH, R = %.4g Ω at f, tuning C = %.4g nF; clearance %.2f mm\n", l.Shape.SelfL()*1e6, l.Shape.Resistance(l.F), 1e9/(math.Pow(2*math.Pi*l.F, 2)*l.Shape.SelfL()), mr.Feasible(l)*1e3)
	ac := ga.AC
	fmt.Printf("    AC books: sources %.4g + %.4g W = load %.4g W + copper %.4g W + source resistors %.4g W (difference %.1e W)\n",
		ac.PSrc[0], ac.PSrc[1], ac.PLoad, ac.PCopper, ac.PSrcR, ac.PSrc[0]+ac.PSrc[1]-ac.PLoad-ac.PCopper-ac.PSrcR)
	fmt.Printf("    group current phases [°]: %.1f; load voltage amplitudes [V]: %.4g\n", GroupPhases(ac), absAll(ac.VLoad))
	// loading: the same drive without the groups
	nb := p3.net(l.Sense)
	nb.SetFreq(l.F)
	r0, _ := nb.SolveAC(l.F, nb.DriveSources(l.Sense, Part3Volts), nil, nil, Separate, 1)
	fmt.Printf("    reaction on the drives: without the groups the sources give %.4g + %.4g W (|I| %.4g, %.4g A); with them %.4g + %.4g W (|I| %.4g, %.4g A)\n",
		r0.PSrc[0], r0.PSrc[1], cmplx.Abs(r0.ISrc[0]), cmplx.Abs(r0.ISrc[1]), ac.PSrc[0], ac.PSrc[1], cmplx.Abs(ac.ISrc[0]), cmplx.Abs(ac.ISrc[1]))
	fmt.Printf("    efficiency at AC (load / sources) %.4f\n", ac.PLoad/(ac.PSrc[0]+ac.PSrc[1]))
	for _, sec := range []int{4, 6, 12} {
		mm := NewModel(NewPairSec(Rodin, Circuits, sec))
		pw, _, err := NewPart3(mm).Score(l)
		fmt.Printf("    discretization check: the same layout on %2d sections per circuit gives %.4g W (%v)\n", sec, pw, errStr(err))
	}
	// combination check: the same groups combined the other ways
	for _, cmb := range []Combine{Separate, Series, Parallel} {
		ll := l
		ll.Comb = cmb
		pw, rr, err := p3.Score(ll)
		if err != nil {
			continue
		}
		fmt.Printf("    combined %-32s %.4g W (efficiency %.4f, bridge factor %.3f)\n", cmb.String()+":", pw, pw/(rr.PSrc[0]+rr.PSrc[1]), BridgeFactor(rr))
	}
	// frequency sweep of the winner's layout (the groups' tuning follows f)
	fmt.Print("    the same layout at other drive frequencies (groups retuned), P_load [W]:")
	for _, f := range []float64{1e3, 2e3, 4e3, 1e4, 3e4, 1e5, 1e6, 5e6} {
		ll := l
		ll.F = f
		if pw, _, err := p3.Score(ll); err == nil {
			fmt.Printf(" %.3g Hz: %.3g;", f, pw)
		} else {
			fmt.Printf(" %.3g Hz: infeasible;", f)
		}
	}
	fmt.Println()
	// 3. the winner in a kronspice transient with bridges and the bank
	vrun(mr, p3, l, ac, tag)
}

func absAll(z []complex128) []float64 {
	out := make([]float64, len(z))
	for i, v := range z {
		out[i] = cmplx.Abs(v)
	}
	return out
}

// e6Down rounds down to the E6 series.
func e6Down(x float64) float64 {
	d := math.Pow(10, math.Floor(math.Log10(x)))
	best := d
	for _, m := range []float64{1, 1.5, 2.2, 3.3, 4.7, 6.8} {
		if m*d <= x {
			best = m * d
		}
	}
	return best
}

func errStr(err error) string {
	if err == nil {
		return "ok"
	}
	return err.Error()
}

// vrun runs the transient verification.
func vrun(mr *Model, p3 *Part3, l Layout, ac *ACResult, tag string) {
	// run length: the groups ring up over their loaded Q, the drive over L/R
	wl := 2 * math.Pi * l.F
	ql := wl * l.Shape.SelfL() / (l.Shape.Resistance(l.F) + l.RLoad)
	if l.Comb != Separate {
		ql = wl * l.Shape.SelfL() * float64(l.NG) / (float64(l.NG)*l.Shape.Resistance(l.F) + l.RLoad)
	}
	tauDrive := mr.LA / (SourceOhms + WireResistance(mr.LenA, l.F))
	cyc := int(math.Max(10*ql/math.Pi, 10*tauDrive*l.F))
	cyc = max(400, min(cyc, 8000))
	steps := 60
	if *short {
		cyc, steps = min(cyc, 300), 40
	}
	run := float64(cyc) / l.F
	// bank sized to fill in about the run: C = P·t/V², V the load voltage
	// amplitude in the AC solution (rounded down to the E6 series)
	vpk := 0.0
	for _, v := range ac.VLoad {
		vpk = math.Max(vpk, cmplx.Abs(v))
	}
	bank := e6Down(ac.PLoad * run / (vpk * vpk))
	n := p3.net(l.Sense)
	t0 := time.Now()
	tr, err := mr.TranVerify(n, l, true, cyc, steps, bank)
	if err != nil {
		log.Fatal("transient: ", err)
	}
	fmt.Printf("\n  kronspice transient of the winner: %d cycles (%.4g ms) at %d steps per cycle (group loaded Q %.3g, drive L/R %.3g ms),\n", cyc, run*1e3, steps, ql, tauDrive*1e3)
	fmt.Printf("  Schottky bridges (is %.0e A, n %.2f, rs %.2f Ω), bank %.4g µF sized to fill within the run (AC load voltage %.4g V peak); %.0f s\n",
		DiodeModel["is"], DiodeModel["n"], DiodeModel["rs"], bank*1e6, vpk, time.Since(t0).Seconds())
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
	rad := RadiationEstimate(l)
	fmt.Printf("    radiation: each coil a magnetic dipole, R_rad = 20π²N²(ka)⁴ Ω = %.3g Ω against its wire %.3g Ω (not modelled)\n", rad, l.Shape.Resistance(l.F))
	fmt.Printf("    external channels other than electrical and thermal (%v among them) carry nothing in these runs: every joule is booked against the two drive sources\n", exchange.Neutrino)
	rows := [][]string{{"t", "vbank", "ebank"}}
	step := max(1, nT/3000)
	for k := 0; k < nT; k += step {
		rows = append(rows, []string{f6(tr.Time[k]), f6(tr.VBank[k]), f6(tr.EBank[k])})
	}
	writeCSV("bank"+tag+".csv", rows)
	g := [][]string{{"theta_deg", "rho", "z", "beta_deg", "radius", "turns"}}
	for _, t := range l.Thetas() {
		g = append(g, []string{f6(t * 180 / math.Pi), f6(l.Shape.Rho), f6(l.Shape.Z), f6(l.Shape.Beta * 180 / math.Pi), f6(l.Shape.Radius), strconv.Itoa(l.Shape.Turns)})
	}
	writeCSV("groups"+tag+".csv", g)
}

// RadiationEstimate is a small loop's radiation resistance (Balanis 5-24).
func RadiationEstimate(l Layout) float64 {
	k := 2 * math.Pi * l.F / em.SpeedOfLight
	ka := k * l.Shape.Radius
	n := float64(l.Shape.Turns)
	return 20 * math.Pi * math.Pi * n * n * math.Pow(ka, 4)
}

func nearest(xs []float64, x float64) int {
	b, bi := math.Inf(1), 0
	for i, v := range xs {
		if d := math.Abs(math.Log(v / x)); d < b {
			b, bi = d, i
		}
	}
	return bi
}

// dominantHarmonic returns the strongest non-zero angular harmonic of a
// periodic sequence and its amplitude relative to the mean.
func dominantHarmonic(y []float64) (int, float64) {
	n := len(y)
	mean := 0.0
	for _, v := range y {
		mean += v
	}
	mean /= float64(n)
	best, bm := 0.0, 0
	for mh := 1; mh <= n/2; mh++ {
		var c complex128
		for k, v := range y {
			c += complex(v, 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(mh*k)/float64(n)))
		}
		a := 2 * cmplx.Abs(c) / float64(n)
		if a > best {
			best, bm = a, mh
		}
	}
	return bm, best / mean
}
