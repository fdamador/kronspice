// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
)

// Ratings the evolved design must respect (declared from the parts in the
// bill of materials). A design beyond one is penalized in proportion to
// the excess, so the GA cannot buy output power by overrating a part.
const (
	limVds    = 150.0  // IRFP260N is a 200 V part; 25 % margin for switching spikes [V]
	limIsup   = 20.0   // 24 V, 20 A supply (480 W) [A]
	limI2     = 30.0   // coil 2 strand, Ø2 mm copper at tens of kHz, air-cooled [A rms]
	limI3     = 8.0    // coil 3, Ø1 mm copper [A rms]
	limVC3    = 2000.0 // C3: polypropylene film, 2 kV DC / ~700 V AC class; peak [V]
	limVHV    = 20e3   // HV terminal to ground, insulation of coil 1 [V]
	limPHV    = 400.0  // the 48 V HV generator's power [W]
	pRef      = 100.0  // fitness scale: error = pRef/(pRef + P_out)
	noRoyer   = 1.0    // penalty when the Royer does not oscillate (latched or dead)
	evalLimit = 90 * time.Second
)

// syncMode selects the synchronized search (flag -sync): the HV amplifier is
// clocked by the Royer's tank, so the genes are its phase and gain instead
// of a free frequency and starting phase.
var syncMode bool

// Space is the genome: every gene is something a builder chooses.
func Space() evo.Space {
	s := evo.Space{
		{Name: "N3", Min: 10, Max: 200, Integer: true},         // output coil turns
		{Name: "N2", Min: 2, Max: 20, Integer: true},           // Royer coil turns per strand
		{Name: "gap23", Min: 1e-3, Max: 20e-3},                 // radial clearance coil 3 → coil 2 [m]
		{Name: "N1s", Min: 50, Max: 1500, Integer: true},       // HV secondary turns
		{Name: "N1p", Min: 2, Max: 30, Integer: true},          // HV primary turns
		{Name: "C2", Min: 100e-9, Max: 4.7e-6, Series: "E12"},  // Royer tank capacitor
		{Name: "C3", Min: 100e-12, Max: 2.2e-6, Series: "E12"}, // output capacitor
		{Name: "RL", Min: 1, Max: 10e3, Series: "E24"},         // load
	}
	if syncMode {
		return append(s,
			evo.Gene{Name: "syncphase", Min: 1, Max: 359}, // HV drive phase relative to the tank [deg]
			evo.Gene{Name: "gain", Min: 0.05, Max: 3})     // HV amplifier voltage gain from the tank voltage
	}
	return append(s,
		evo.Gene{Name: "fratio", Min: 0.7, Max: 1.4}, // HV generator frequency / unloaded Royer tank frequency
		evo.Gene{Name: "phase", Min: 1, Max: 360})    // HV generator phase at t = 0 [deg]
}

// Decode turns gene values (Space order) into a design.
func Decode(v []float64) Design {
	d := Design{N3: int(v[0]), N2: int(v[1]), Gap23: v[2], N1s: int(v[3]), N1p: int(v[4]),
		C2: v[5], C3: v[6], RL: v[7]}
	if syncMode {
		d.Sync = true
		d.SetSyncAngle(v[8])
		d.SyncGain = v[9]
		return d
	}
	d.PhiHV = v[9]
	d.FHV = v[8] * d.RoyerFreq0()
	return d
}

// SetSyncAngle sets the slaved drive's phase from an angle in 0–360°:
// the all-pass gives 0–180° of lag, the sign the other half turn.
func (d *Design) SetSyncAngle(a float64) {
	a = math.Mod(math.Mod(a, 360)+360, 360)
	d.SyncSign = 1
	if a > 180 {
		d.SyncSign, a = -1, a-180
	}
	d.SyncPhase = a
}

// SyncAngle is the inverse of SetSyncAngle.
func (d Design) SyncAngle() float64 {
	if d.SyncSign < 0 {
		return d.SyncPhase + 180
	}
	return d.SyncPhase
}

// Score is the fitness error of a measured design (lower is better).
func Score(r Result) float64 {
	if r.Err != nil {
		p, _ := evo.Penalty(r.Err)
		return p
	}
	over := func(x, lim float64) float64 { return math.Max(0, x/lim-1) }
	e := pRef / (pRef + math.Max(r.POut, 0))
	e += over(r.VdsPk, limVds) + over(math.Abs(r.IChoke), limIsup) + over(r.I2rms, limI2) +
		over(r.I3rms, limI3) + over(r.VC3Pk, limVC3) + over(r.VHVPk, limVHV) + over(r.PInHV, limPHV)
	if r.FRoyer == 0 && !r.Design.RoyerOff {
		e += noRoyer
	}
	return e
}

func objective() evo.Objective {
	return evo.Objective{Name: "output power (with ratings)", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
		d := Decode(e.Values)
		r := Result{Design: d, FR0: d.RoyerFreq0()}
		wf, err := e.Tran(analysis.TranSpec{Step: maxStep, Stop: tStop, MaxStep: maxStep})
		if err != nil {
			return 0, err
		}
		if err := r.measure(wf); err != nil {
			return 0, err
		}
		return Score(r), nil
	}}
}

// Evolve runs the GA and returns its final state.
func Evolve(ctx context.Context, gens, pop int, seed uint64, log func(evo.GenStats)) (*evo.State, error) {
	p := evo.Problem{
		Space: Space(),
		Build: func(v []float64) (*circuit.Circuit, error) {
			c, _, err := Decode(v).Build()
			return c, err
		},
		Objectives: []evo.Objective{objective()},
	}
	opt := evo.Options{Population: pop, Generations: gens, Seed: seed, Elite: 2, Timeout: evalLimit, Sigma: 0.15}
	return evo.Run(ctx, p, opt, nil, func(st *evo.State) error {
		if log != nil {
			log(st.History[len(st.History)-1])
		}
		return nil
	})
}

// parallel simulates every design on all CPUs, results in input order.
func parallel(ctx context.Context, ds []Design, keep bool) []Result {
	out := make([]Result, len(ds))
	jobs := make(chan int)
	var wg sync.WaitGroup
	for w := 0; w < runtime.NumCPU(); w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				out[i] = Simulate(ctx, ds[i], keep)
			}
		}()
	}
	for i := range ds {
		jobs <- i
	}
	close(jobs)
	wg.Wait()
	return out
}

func run(ctx context.Context, dir string, gens, pop int, seed uint64) error {
	mode := "free-running HV generator"
	if syncMode {
		mode = "HV amplifier clocked by the Royer (synchronized)"
	}
	fmt.Printf("\n== Genetic evolution (%s): population %d, %d generations, seed %d ==\n", mode, pop, gens, seed)
	t0 := time.Now()
	var hist []evo.GenStats
	st, err := Evolve(ctx, gens, pop, seed, func(g evo.GenStats) {
		hist = append(hist, g)
		fmt.Printf("  gen %3d  best %.5f  mean %.4f  failed %d  timed out %d  (P_out ≥ %.1f W)\n",
			g.Generation, g.Best, g.Mean, g.Failed, g.TimedOut, pRef/g.Best-pRef)
	})
	if err != nil {
		return err
	}
	fmt.Printf("  %d evaluations in %s\n", st.Evaluations, time.Since(t0).Round(time.Second))
	best := Decode(st.Best.Values)
	fmt.Println("\n== Evolved design ==")
	rb := Simulate(ctx, best, true)
	report(rb)
	fmt.Printf("  fitness %.5f (GA %.5f)\n", Score(rb), st.Best.Fitness)
	printLosses(rb)

	// Comparisons: each drive alone, both together.
	hv, roy := best, best
	hv.RoyerOff, roy.HVOff = true, true
	if syncMode {
		// The slaved amplifier stops with the Royer; drive coil 1 alone with
		// a sine at the same frequency and amplitude it had when slaved.
		hv = best
		hv.Sync, hv.RoyerOff = false, true
		hv.FHV, hv.PhiHV = rb.FRoyer, 0
		hv.VgAmp = rb.VgPk
	}
	base := Baseline()
	bhv, broy := base, base
	bhv.RoyerOff, broy.HVOff = true, true
	cmp := parallel(ctx, []Design{hv, roy, base, bhv, broy}, false)
	fmt.Println("\n== Each drive alone against both (evolved design) ==")
	for i, name := range []string{"HV coil only", "Royer only"} {
		fmt.Printf("  %-13s P_out %8.3f W   in: HV %8.2f W, Royer %8.2f W   η %.3f\n", name, cmp[i].POut, cmp[i].PInHV, cmp[i].PInRoy, cmp[i].Eta)
	}
	fmt.Printf("  %-13s P_out %8.3f W   in: HV %8.2f W, Royer %8.2f W   η %.3f\n", "both", rb.POut, rb.PInHV, rb.PInRoy, rb.Eta)
	sumAlone := cmp[0].POut + cmp[1].POut
	coh := math.Pow(math.Sqrt(cmp[0].POut)+math.Sqrt(cmp[1].POut), 2)
	fmt.Printf("  sum of the two alone %.3f W; coherent (in-phase) sum (√P1 + √P2)² = %.3f W; both measured %.3f W\n", sumAlone, coh, rb.POut)
	fmt.Println("\n== Same comparison, hand-picked baseline ==")
	for i, name := range []string{"both", "HV coil only", "Royer only"} {
		r := cmp[2+i]
		fmt.Printf("  %-13s P_out %8.3f W   in: HV %8.2f W, Royer %8.2f W\n", name, r.POut, r.PInHV, r.PInRoy)
	}

	// Trends around the evolved design.
	fmt.Println("\n== Trends around the evolved design ==")
	type sweep struct {
		name, unit string
		xs         []float64
		set        func(d *Design, x float64)
	}
	logspace := func(a, b float64, n int) []float64 {
		v := make([]float64, n)
		for i := range v {
			v[i] = a * math.Pow(b/a, float64(i)/float64(n-1))
		}
		return v
	}
	linspace := func(a, b float64, n int) []float64 {
		v := make([]float64, n)
		for i := range v {
			v[i] = a + (b-a)*float64(i)/float64(n-1)
		}
		return v
	}
	sweeps := []sweep{
		{"detune", "f_HV / f_HV,best", linspace(0.85, 1.15, 31), func(d *Design, x float64) { d.FHV *= x }},
		{"phase", "deg", linspace(0, 330, 12), func(d *Design, x float64) { d.PhiHV = x }},
	}
	if syncMode {
		sweeps = []sweep{
			{"gain", "G / G,best", linspace(0, 1.5, 16), func(d *Design, x float64) { d.SyncGain *= x }},
			{"phase", "deg", linspace(0, 345, 24), func(d *Design, x float64) { d.SetSyncAngle(x) }},
		}
	}
	sweeps = append(sweeps, []sweep{
		{"load", "RL / RL,best", logspace(0.1, 10, 17), func(d *Design, x float64) { d.RL *= x }},
		{"c3", "C3 / C3,best", logspace(0.5, 2, 17), func(d *Design, x float64) { d.C3 *= x }},
		{"gap", "mm", linspace(1, 20, 12), func(d *Design, x float64) { d.Gap23 = x * 1e-3 }},
	}...)
	var trendRows []string
	for _, s := range sweeps {
		ds := make([]Design, len(s.xs))
		for i, x := range s.xs {
			ds[i] = best
			s.set(&ds[i], x)
		}
		rs := parallel(ctx, ds, false)
		fmt.Printf("  %s (%s):\n", s.name, s.unit)
		for i, r := range rs {
			k := 0.5 * (r.E.K2a3 + r.E.K2b3)
			fmt.Printf("    x %-9.4g P_out %9.3f W  HV in %8.2f  Royer in %8.2f  f_R %8.0f  locked %-5v beat %.3f  k23 %.3f  err %v\n",
				s.xs[i], r.POut, r.PInHV, r.PInRoy, r.FRoyer, r.Locked, r.Beat, k, r.Err != nil)
			trendRows = append(trendRows, fmt.Sprintf("%s,%g,%g,%g,%g,%g,%g,%t,%g,%g,%g", s.name, s.xs[i], r.POut, r.PInHV, r.PInRoy,
				r.FRoyer, r.Design.FHV, r.Locked, r.Beat, k, r.Eta))
		}
	}

	if dir == "" {
		return nil
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	if err := writeCSV(filepath.Join(dir, "trends.csv"), "sweep,x,p_out,p_hv,p_roy,f_royer,f_hv,locked,beat,k23,eta", trendRows); err != nil {
		return err
	}
	var gh []string
	for _, g := range hist {
		gh = append(gh, fmt.Sprintf("%d,%g,%g,%g,%d", g.Generation, g.Best, g.Mean, pRef/g.Best-pRef, g.Failed))
	}
	if err := writeCSV(filepath.Join(dir, "convergence.csv"), "gen,best,mean,p_lower_bound,failed", gh); err != nil {
		return err
	}
	var cr []string
	for _, x := range []struct {
		name string
		r    Result
	}{{"evolved both", rb}, {"evolved HV only", cmp[0]}, {"evolved Royer only", cmp[1]},
		{"baseline both", cmp[2]}, {"baseline HV only", cmp[3]}, {"baseline Royer only", cmp[4]}} {
		cr = append(cr, fmt.Sprintf("%s,%g,%g,%g,%g,%g,%t", x.name, x.r.POut, x.r.PInHV, x.r.PInRoy, x.r.Eta, x.r.FRoyer, x.r.Locked))
	}
	if err := writeCSV(filepath.Join(dir, "compare.csv"), "case,p_out,p_hv,p_roy,eta,f_royer,locked", cr); err != nil {
		return err
	}
	if err := writeDesign(filepath.Join(dir, "design.csv"), rb); err != nil {
		return err
	}
	return writeWaves(ctx, dir, best)
}

// writeWaves simulates a design and writes wave.csv, startup.csv and
// best.cir into dir.
func writeWaves(ctx context.Context, dir string, d Design) error {
	rb := Simulate(ctx, d, true)
	if rb.Err != nil {
		return rb.Err
	}
	if err := writeWave(filepath.Join(dir, "wave.csv"), rb); err != nil {
		return err
	}
	if err := writeStartup(filepath.Join(dir, "startup.csv"), rb); err != nil {
		return err
	}
	// Accuracy check: the same design at a 4× finer largest step.
	fine := SimulateStep(ctx, d, false, maxStep/4)
	if fine.Err != nil {
		return fine.Err
	}
	fmt.Printf("%s: P_out %.5g W at %.3g s steps, %.5g W at %.3g s (%.2g relative); η %.4f / %.4f; mismatch %.2g / %.2g\n",
		dir, rb.POut, maxStep, fine.POut, maxStep/4, (fine.POut-rb.POut)/rb.POut, rb.Eta, fine.Eta, rb.Mismatch, fine.Mismatch)
	if err := writeCSV(filepath.Join(dir, "check.csv"), "p_out,p_out_fine,rel,eta,eta_fine,mismatch,mismatch_fine",
		[]string{fmt.Sprintf("%g,%g,%g,%g,%g,%g,%g", rb.POut, fine.POut, (fine.POut-rb.POut)/rb.POut, rb.Eta, fine.Eta, rb.Mismatch, fine.Mismatch)}); err != nil {
		return err
	}
	c, _, err := d.Build()
	if err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, "best.cir"))
	if err != nil {
		return err
	}
	defer f.Close()
	return c.WriteNetlist(f, fmt.Sprintf(".tran %g %g 0 %g", maxStep, tStop, maxStep))
}

// printLosses lists where the input power goes.
func printLosses(r Result) {
	if r.Err != nil {
		return
	}
	keys := make([]string, 0, len(r.PDiss))
	for k := range r.PDiss {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return r.PDiss[keys[i]] > r.PDiss[keys[j]] })
	var b strings.Builder
	for _, k := range keys {
		fmt.Fprintf(&b, " %s %.2f W;", k, r.PDiss[k])
	}
	fmt.Printf("  losses:%s load %.2f W\n", b.String(), r.POut)
}

func writeCSV(path, head string, rows []string) error {
	return os.WriteFile(path, []byte(head+"\n"+strings.Join(rows, "\n")+"\n"), 0o644)
}

// writeDesign writes the evolved design, its electrical values and the
// steady-state measurements as name,value,unit rows.
func writeDesign(path string, r Result) error {
	d, e := r.Design, r.E
	c1p, c1s, c2a, _, c3 := d.Coils()
	rows := []string{}
	add := func(n string, v float64, u string) { rows = append(rows, fmt.Sprintf("%s,%g,%s", n, v, u)) }
	add("N3", float64(d.N3), "turns")
	add("N2", float64(d.N2), "turns per strand")
	add("gap23", d.Gap23, "m")
	add("N1s", float64(d.N1s), "turns")
	add("N1p", float64(d.N1p), "turns")
	add("C2", d.C2, "F")
	add("C3", d.C3, "F")
	add("RL", d.RL, "ohm")
	add("FHV", d.FHV, "Hz")
	add("PhiHV", d.PhiHV, "deg")
	add("Sync", map[bool]float64{false: 0, true: 1}[d.Sync], "")
	add("SyncAngle", d.SyncAngle(), "deg")
	add("SyncGain", d.SyncGain, "")
	add("VgPk", r.VgPk, "V")
	add("FR0", r.FR0, "Hz")
	add("L3", e.L3, "H")
	add("L2strand", e.L2a, "H")
	add("Ltank", e.L2a+e.L2b+2*e.M2ab, "H")
	add("L1s", e.L1s, "H")
	add("L1p", e.L1p, "H")
	add("K1", e.K1, "")
	add("K2ab", e.K2ab, "")
	add("K23", 0.5*(e.K2a3+e.K2b3), "")
	add("C1self", e.C1self, "F")
	add("R1p", e.R1p, "ohm")
	add("R1s", e.R1s, "ohm")
	add("R2", e.R2, "ohm")
	add("R3", e.R3, "ohm")
	for _, h := range []Helix{c1p, c1s, c2a, c3} {
		key := strings.ReplaceAll(h.Name, " ", "_")
		add(key+"_diameter", 2*h.R, "m")
		add(key+"_length", h.Length(), "m")
		add(key+"_wire_length", h.WireLength(), "m")
	}
	add("POut", r.POut, "W")
	add("PInHV", r.PInHV, "W")
	add("PInRoy", r.PInRoy, "W")
	add("Eta", r.Eta, "")
	add("FRoyer", r.FRoyer, "Hz")
	add("FOut", r.FOut, "Hz")
	add("Beat", r.Beat, "")
	add("VdsPk", r.VdsPk, "V")
	add("IdPk", r.IdPk, "A")
	add("VHVPk", r.VHVPk, "V")
	add("VC3Pk", r.VC3Pk, "V")
	add("VOutPk", r.VOutPk, "V")
	add("I3rms", r.I3rms, "A")
	add("I2rms", r.I2rms, "A")
	add("ISupply", r.IChoke, "A")
	add("Mismatch", r.Mismatch, "")
	for k, v := range r.PDiss {
		add("loss_"+k, v, "W")
	}
	sort.Strings(rows[len(rows)-len(r.PDiss):])
	return writeCSV(path, "name,value,unit", rows)
}

// writeWave writes the last 5 cycles of the steady state.
func writeWave(path string, r Result) error {
	wf := r.Wave
	if wf == nil {
		return fmt.Errorf("no waveforms")
	}
	t := wf.Time
	from := t[len(t)-1] - 5/r.Design.FHV
	i0 := window(t, from)
	hv, out, c1, c2 := node(wf, "hv"), node(wf, "out"), node(wf, "c1"), node(wf, "c2")
	b3 := node(wf, "b3")
	I := func(n string) []float64 { return wf.I[elem(wf, n)] }
	src := "Vg"
	if r.Design.Sync {
		src = "Eg"
	}
	i3, i2a, i2b, ig := I("L3"), I("L2a"), I("L2b"), I(src)
	var rows []string
	for k := i0; k < len(t); k++ {
		rows = append(rows, fmt.Sprintf("%g,%g,%g,%g,%g,%g,%g,%g,%g", (t[k]-t[i0])*1e6, hv[k], out[k], c1[k]-c2[k], b3[k]-out[k], i3[k], 0.5*(i2a[k]+i2b[k]), -ig[k], out[k]/r.Design.RL))
	}
	return writeCSV(path, "t_us,v_hv,v_out,v_tank,v_c3,i_coil3,i_tank,i_hvgen,i_load", rows)
}

// writeStartup writes the output power envelope from switch-on.
func writeStartup(path string, r Result) error {
	wf := r.Wave
	if wf == nil {
		return fmt.Errorf("no waveforms")
	}
	t, out := wf.Time, node(wf, "out")
	c1, c2 := node(wf, "c1"), node(wf, "c2")
	var rows []string
	step := 10e-6
	for a := 0.0; a+step <= t[len(t)-1]+1e-12; a += step {
		i0, i1 := window(t, a), window(t, a+step)
		pk, tk := 0.0, 0.0
		for k := i0; k < i1 && k < len(t); k++ {
			pk = math.Max(pk, math.Abs(out[k]))
			tk = math.Max(tk, math.Abs(c1[k]-c2[k]))
		}
		rows = append(rows, fmt.Sprintf("%g,%g,%g", (a+step/2)*1e3, pk, tk))
	}
	return writeCSV(path, "t_ms,v_out_peak,v_tank_peak", rows)
}
