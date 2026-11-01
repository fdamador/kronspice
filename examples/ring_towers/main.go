// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command ring_towers is a thought experiment as a circuit: two towers at
// opposite points of the Earth drive the same frequency into the
// Earth–ionosphere waveguide, a closed ring around the planet, and their
// phases are shifted so that the four waves they launch (each tower sends
// one wave east and one west) add or subtract.
//
// Circuit equivalent (lumped transmission line, telegrapher's equations;
// Steinmetz, Lecture VI; Paul §2): the ring is N sections of series L and R
// with a shunt C, the last section closing back to the first. Wave speed
// 1/√(LC) per section, characteristic impedance √(L/C), loss R per section.
// Each tower is two current sources a quarter wavelength apart, so its
// launch direction is set by their relative phase (a two-element end-fire
// array): the eastern element lagging 90° launches east only, leading 90°
// launches west only.
//
// States:
//
//	traveling — both towers launch east, tower B phased to the arrival of
//	            A's wave (φ_B = −k·(half circumference)): one wave circulates
//	            around the ring, the same amplitude everywhere, power flowing
//	            around the planet.
//	standing  — A launches east, B launches west: counter-propagating waves
//	            of equal amplitude form a standing wave with fixed nodes and
//	            antinodes; energy sloshes in place, no net flow.
//	cancel    — both launch east, B phased 180° from A's arrival: the waves
//	            subtract to nothing; neither tower delivers power.
//
// Four towers: a second pair C (90°) and D (270°) is added at right angles
// to A–B. Each pair alone makes a standing wave (C launches east, D west,
// like A and B). With n = 4 wavelengths around the ring every quarter
// point — 0°, 90°, 180°, 270° — is an antinode of the A–B pattern
// cos(kx), all four swinging in step (kx = 2π·m there). Phasing C and D to
// the same reference waves as A and B (east wave −kx, west wave +kx at
// their own positions) puts the C–D pattern exactly on top of A–B's, in
// space and in time, so the two add: double the amplitude at every
// antinode, four times the stored energy, from four towers each feeding
// twice the power of one tower of a single pair. Two contrasts show why
// the synchronization matters: C–D in antiphase cancels A–B, and C–D
// shifted by an eighth of a wavelength in phase moves the combined
// antinodes off the tower points.
//
// The ring runs at its n-th resonance (n wavelengths around), so waves
// returning after each lap add to the new ones; the build-up is limited by
// the loss per lap, and the ledger shows every joule in the ring supplied
// by the towers. Scale: Earth's ring is 40,000 km and its Schumann
// fundamental is 7.8 Hz; here the ring is scaled to MHz so each lap takes
// microseconds. The patterns depend only on n and the phases.
//
//	go run ./examples/ring_towers [-json out.json]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Ring parameters.
const (
	N    = 96   // sections around the ring
	Lsec = 1e-6 // H per section
	Csec = 1e-9 // F per section
	Nwav = 4    // wavelengths around the ring
	Amp  = 0.01 // A, each radiator
	Loss = 0.2  // amplitude loss per lap [Np]
)

// Rsec is the series resistance per section giving Loss nepers per lap:
// α = R/(2·Z0) per section.
func Rsec() float64 { return Loss * 2 * math.Sqrt(Lsec/Csec) / N }

// Freq is the n-th ring resonance of the LC ladder, ω = (2/√(LC))·sin(kΔ/2)
// with kΔ = 2πn/N (lattice dispersion; ≈ n·v/circumference).
func Freq() float64 {
	return 2 / math.Sqrt(Lsec*Csec) * math.Sin(math.Pi*Nwav/N) / (2 * math.Pi)
}

// Tower B sits at the antipode; each tower's second element is a quarter
// wavelength east of its first.
const (
	nodeA = 0
	nodeB = N / 2
	quart = N / (4 * Nwav)
)

// Radiator is one current source of a tower.
type Radiator struct {
	Name  string
	Node  int
	Phase float64 // deg
}

// State is a set of radiators with their phases.
type State struct {
	Name string
	Rad  []Radiator
}

func pair(tag string, at int, p1, p2 float64) []Radiator {
	return []Radiator{{tag + "1", at, p1}, {tag + "2", at + quart, p2}}
}

func fold(d float64) float64 { return math.Mod(math.Mod(d, 360)+360, 360) }

// k is the phase step per section [deg].
const kDeg = 360.0 * Nwav / N

// East launches an eastward wave from node x0 whose phase at node x is
// −k·x + off (reference east wave); West launches a westward wave whose
// phase is +k·x + off. Each is a two-element end-fire pair λ/4 apart.
func East(tag string, x0 int, off float64) []Radiator {
	p := fold(-kDeg*float64(x0) + off)
	return pair(tag, x0, p, fold(p-90))
}

// West: see East.
func West(tag string, x0 int, off float64) []Radiator {
	p := fold(kDeg*float64(x0) + off)
	return pair(tag, x0, p, fold(p+90))
}

func join(rs ...[]Radiator) []Radiator {
	var out []Radiator
	for _, r := range rs {
		out = append(out, r...)
	}
	return out
}

// Tower positions [nodes]: A 0°, C 90°, B 180°, D 270°, and the second
// set rotated by 45°: E 45°, G 135°, F 225°, H 315°.
const (
	nodeC = N / 4
	nodeD = 3 * N / 4
	nodeE = N / 8
	nodeG = 3 * N / 8
	nodeF = 5 * N / 8
	nodeH = 7 * N / 8
)

// EightStates adds the 45°-rotated pairs E–F (45°, 225°) and G–H (135°,
// 315°), phased to the same reference waves. With n = 4 wavelengths
// around the ring the 45° points are antinodes too, swinging opposite to
// the quarter points (cos 4θ = −1 there), so all four pairs' standing
// waves coincide and add: four times one pair's amplitude.
func EightStates() []State {
	ab := join(East("A", nodeA, 0), West("B", nodeB, 0))
	cd := join(East("C", nodeC, 0), West("D", nodeD, 0))
	ef := join(East("E", nodeE, 0), West("F", nodeF, 0))
	gh := join(East("G", nodeG, 0), West("H", nodeH, 0))
	return []State{
		{"E–F standing", ef},
		{"eight towers, synchronized", join(ab, cd, ef, gh)},
		{"eight towers, 45° set in antiphase", join(ab, cd, join(East("E", nodeE, 180), West("F", nodeF, 180)), join(East("G", nodeG, 180), West("H", nodeH, 180)))},
	}
}

// FourStates returns the two-pair settings: each pair alone, both
// synchronized (antinodes add at 0°, 90°, 180°, 270°), C–D in antiphase,
// and C–D shifted by 45° of phase (an eighth of a wavelength).
func FourStates() []State {
	ab := join(East("A", nodeA, 0), West("B", nodeB, 0))
	cd := func(off float64) []Radiator { return join(East("C", nodeC, off), West("D", nodeD, off)) }
	return []State{
		{"A–B standing", ab},
		{"C–D standing", cd(0)},
		{"both, synchronized", join(ab, cd(0))},
		{"both, C–D antiphase", join(ab, cd(180))},
		{"both, C–D shifted 45°", join(ab, join(East("C", nodeC, -45), West("D", nodeD, 45)))},
	}
}

// States returns the three phase settings. φ_B = −k·(N/2 sections) is the
// phase of A's eastward wave arriving at B.
func States() []State {
	phiB := -360.0 * Nwav * (N / 2) / N
	phiB = math.Mod(phiB, 360)
	st := func(name string, p [4]float64) State {
		return State{name, join(pair("A", nodeA, p[0], p[1]), pair("B", nodeB, p[2], p[3]))}
	}
	return []State{
		st("traveling", [4]float64{0, -90, phiB, phiB - 90}),
		st("standing", [4]float64{0, -90, phiB, phiB + 90}),
		st("cancel", [4]float64{0, -90, phiB + 180, phiB + 90}),
	}
}

func node(k int) string { return fmt.Sprintf("n%d", ((k%N)+N)%N) }

// Circuit builds the ring and the four radiators.
func Circuit(s State) (*circuit.Circuit, error) { return build(s, nil, true) }

// Bias puts a DC voltage on each tower pair: a supply of VAB at towers A
// and B and one of VCD at C and D, each through its own feed resistor
// RFeed to the tower's first node (high enough that the feeds carry
// almost no RF). The ring is one conductor at DC — its inductors connect
// every node — so with a continuous ring the two pairs cannot hold
// different DC levels: the whole ring settles at the average and a steady
// current (VAB − VCD)/RFeed flows from one pair's supplies to the other's.
// Sectored inserts DC-blocking capacitors CBlock in series in the ring at
// 45°, 135°, 225° and 315°, splitting it into four sectors (one per
// tower) that each hold their own tower's bias while the RF passes
// through the capacitors (reactance ≪ section impedance).
type Bias struct {
	VAB, VCD float64 // V
	RFeed    float64 // Ω
	Sectored bool
	CBlock   float64 // F
}

// DefaultBias returns RFeed 1 MΩ and CBlock 10 µF (reactance 0.012 Ω at
// the ring frequency, against ωL = 8.3 Ω per section).
func DefaultBias(vab, vcd float64, sectored bool) *Bias {
	return &Bias{VAB: vab, VCD: vcd, RFeed: 1e6, Sectored: sectored, CBlock: 10e-6}
}

// sectorEdge reports whether a blocking capacitor follows section k.
func sectorEdge(k int) bool { return k == N/8-1 || k == 3*N/8-1 || k == 5*N/8-1 || k == 7*N/8-1 }

// SectorOf returns the tower whose sector holds node k (sectored ring):
// 'A' 315°–45°, 'C' 45°–135°, 'B' 135°–225°, 'D' 225°–315°.
func SectorOf(k int) byte {
	return "ACBD"[((k+N/8)%N)/(N/4)]
}

// DCLevels returns each node's DC voltage for bias b (analytic): the
// average of the two supplies on a continuous ring (equal feeds), the
// sector's own supply on a sectored ring.
func DCLevels(b *Bias) []float64 {
	out := make([]float64, N)
	for k := range out {
		switch {
		case b == nil:
		case !b.Sectored:
			out[k] = (b.VAB + b.VCD) / 2
		case SectorOf(k) == 'A' || SectorOf(k) == 'B':
			out[k] = b.VAB
		default:
			out[k] = b.VCD
		}
	}
	return out
}

func build(s State, bias *Bias, withSources bool) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("ring towers: " + s.Name)
	r := Rsec()
	dc := DCLevels(bias)
	for k := 0; k < N; k++ {
		mid := fmt.Sprintf("m%d", k)
		b.Inductor(fmt.Sprintf("L%d", k), node(k), mid, Lsec)
		if bias != nil && bias.Sectored && sectorEdge(k) {
			blk := fmt.Sprintf("x%d", k)
			b.Resistor(fmt.Sprintf("R%d", k), mid, blk, r)
			b.CapacitorIC(fmt.Sprintf("CB%d", k), blk, node(k+1), bias.CBlock, dc[k]-dc[(k+1)%N])
		} else {
			b.Resistor(fmt.Sprintf("R%d", k), mid, node(k+1), r)
		}
		b.CapacitorIC(fmt.Sprintf("C%d", k), node(k), "0", Csec, dc[k])
	}
	if bias != nil {
		for _, t := range []struct {
			name string
			at   int
			v    float64
		}{{"A", nodeA, bias.VAB}, {"B", nodeB, bias.VAB}, {"C", nodeC, bias.VCD}, {"D", nodeD, bias.VCD}} {
			b.VoltageSource("VB"+t.name, "vb"+t.name, "0", t.v)
			b.Resistor("RF"+t.name, "vb"+t.name, node(t.at), bias.RFeed)
		}
	}
	if withSources {
		f := Freq()
		for _, r := range s.Rad {
			b.CurrentWave("I"+r.Name, "0", node(r.Node), circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, Amp, f, 0, 0, r.Phase}})
		}
		return b.BuildTransient()
	}
	return b.Build()
}

// BiasOP solves the DC operating point of the biased ring alone (no RF
// sources): node voltages and the steady current each supply delivers.
func BiasOP(bias *Bias) (levels []float64, supply map[string]float64, err error) {
	c, err := build(State{Name: "bias"}, bias, false)
	if err != nil {
		return nil, nil, err
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	levels = make([]float64, N)
	for k := range levels {
		levels[k] = op.NodeVoltages[c.NodeIndex(node(k))]
	}
	supply = map[string]float64{}
	for _, t := range []string{"A", "B", "C", "D"} {
		vb := op.NodeVoltages[c.NodeIndex("vb"+t)]
		var at int
		switch t {
		case "A":
			at = nodeA
		case "B":
			at = nodeB
		case "C":
			at = nodeC
		default:
			at = nodeD
		}
		supply[t] = (vb - levels[at]) / bias.RFeed // current out of the supply into the ring [A]
	}
	return levels, supply, nil
}

// Result is the steady state of one state.
type Result struct {
	Name     string
	Rad      []Radiator
	Envelope []float64   // peak |V| per node over the last period [V]
	Frames   [][]float64 // V per node at 16 instants of the last period
	IFrames  [][]float64 // I per section at the same instants
	MeanE    []float64   // cycle-mean electric energy ½CV² per node [J]
	MeanM    []float64   // cycle-mean magnetic energy ½LI² per section [J]
	PTowers  []float64   // mean power each radiator delivers [W]
	PLoss    float64     // mean power lost in the ring [W]
	Flow     []float64   // mean power flowing east through each section [W]
	Stored   float64     // mean stored energy [J]
	// Phasor is each node's voltage over the last period as
	// V(t) = Phasor[k][0]·cos ω(t−t0) + Phasor[k][1]·sin ω(t−t0), t0 a whole
	// number of periods (so source phases are measured from the same
	// origin as the SIN sources' phase argument).
	Phasor   [][2]float64
	MeanV    []float64 // cycle-mean (DC) voltage per node over the last period [V]
	Mismatch float64   // ledger integration mismatch over the run [J]
	EnergyIn float64   // [J]
}

// Run simulates a state to steady state (8 amplitude-decay times).
func Run(s State) (*Result, error) { return RunBiased(s, nil) }

// RunBiased is Run with a DC bias on the tower pairs (nil for none). The
// run starts from the bias's DC levels (capacitor initial conditions), so
// only the RF has to build up.
func RunBiased(s State, bias *Bias) (*Result, error) {
	c, err := build(s, bias, true)
	if err != nil {
		return nil, err
	}
	f := Freq()
	T := 1 / f
	lap := N * math.Sqrt(Lsec*Csec)
	tEnd := math.Ceil(8*lap/Loss/T) * T
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: T / 64, Stop: tEnd, MaxStep: T / 64, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return nil, err
	}
	res := &Result{Name: s.Name, Rad: s.Rad, PTowers: make([]float64, len(s.Rad)), Envelope: make([]float64, N), Flow: make([]float64, N),
		MeanE: make([]float64, N), MeanM: make([]float64, N),
		Mismatch: wf.Ledger.IntegrationMismatch, EnergyIn: wf.Ledger.EnergyIn}
	v := make([][]float64, N)
	il := make([][]float64, N)
	for k := 0; k < N; k++ {
		v[k], _ = wf.NodeVoltage(node(k))
		il[k], _ = wf.ElementCurrent(fmt.Sprintf("L%d", k))
	}
	n := len(wf.Time)
	t0 := wf.Time[n-1] - T
	k0 := 0
	for k0 < n && wf.Time[k0] < t0 {
		k0++
	}
	var span float64
	res.Phasor = make([][2]float64, N)
	res.MeanV = make([]float64, N)
	// Integrate over exactly one period [t0, tEnd]: the first segment is
	// clipped at t0 by linear interpolation, so DC levels do not leak into
	// the Fourier coefficients and every average covers a whole cycle.
	srcs := make([][]float64, len(s.Rad))
	for i, r := range s.Rad {
		srcs[i], _ = wf.ElementCurrent("I" + r.Name)
	}
	for j := max(k0, 1); j < n; j++ {
		ta, tb := wf.Time[j-1], wf.Time[j]
		frac := 0.0
		if ta < t0 {
			frac = (t0 - ta) / (tb - ta)
			ta = t0
		}
		prev := func(x []float64) float64 { return x[j-1] + (x[j]-x[j-1])*frac }
		dt := tb - ta
		span += dt
		wa, wb := 2*math.Pi*f*(ta-t0), 2*math.Pi*f*(tb-t0)
		for k := 0; k < N; k++ {
			va, vb := prev(v[k]), v[k][j]
			ia, ib := prev(il[k]), il[k][j]
			res.Envelope[k] = math.Max(res.Envelope[k], math.Abs(vb))
			res.Flow[k] += dt * 0.5 * (vb*ib + va*ia)
			res.PLoss += dt * 0.5 * Rsec() * (ib*ib + ia*ia)
			ee := dt * 0.5 * 0.5 * Csec * (vb*vb + va*va)
			em := dt * 0.5 * 0.5 * Lsec * (ib*ib + ia*ia)
			res.MeanE[k] += ee
			res.MeanM[k] += em
			res.Stored += ee + em
			res.MeanV[k] += dt * 0.5 * (va + vb)
			res.Phasor[k][0] += dt * (va*math.Cos(wa) + vb*math.Cos(wb)) / T
			res.Phasor[k][1] += dt * (va*math.Sin(wa) + vb*math.Sin(wb)) / T
		}
		for i, r := range s.Rad {
			at := r.Node % N
			// The source drives its current from ground into the node, so
			// it delivers V(node)·i.
			res.PTowers[i] += dt * 0.5 * (v[at][j]*srcs[i][j] + prev(v[at])*prev(srcs[i]))
		}
	}
	for k := range res.Flow {
		res.MeanV[k] /= span
		res.Flow[k] /= span
		res.MeanE[k] /= span
		res.MeanM[k] /= span
	}
	res.PLoss /= span
	res.Stored /= span
	for i := range res.PTowers {
		res.PTowers[i] /= span
	}
	for fr := 0; fr < 16; fr++ {
		tt := t0 + T*float64(fr)/16
		j := k0
		for j < n-1 && wf.Time[j] < tt {
			j++
		}
		frame, iframe := make([]float64, N), make([]float64, N)
		for k := 0; k < N; k++ {
			frame[k], iframe[k] = v[k][j], il[k][j]
		}
		res.Frames = append(res.Frames, frame)
		res.IFrames = append(res.IFrames, iframe)
	}
	return res, nil
}

// Uniformity is min/max of the envelope: 1 for a pure traveling wave, near
// 0 for a pure standing wave.
func (r *Result) Uniformity() float64 {
	lo, hi := math.Inf(1), 0.0
	for _, e := range r.Envelope {
		lo, hi = math.Min(lo, e), math.Max(hi, e)
	}
	return lo / hi
}

func main() {
	out := flag.String("json", "", "write frames and envelopes as JSON")
	ph := flag.String("phasors", "", "write each launcher's steady-state voltage phasors (the data behind viz/*.html) as JSON and exit")
	flag.Parse()
	if *ph != "" {
		if err := WritePhasors(*ph); err != nil {
			log.Fatal(err)
		}
		return
	}
	z0 := math.Sqrt(Lsec / Csec)
	fmt.Printf("Ring: %d sections, Z0 %.1f Ω, %d wavelengths around, f = %.4g MHz, loss %.2f Np per lap (R %.3g Ω/section)\n",
		N, z0, Nwav, Freq()/1e6, Loss, Rsec())
	fmt.Printf("Towers: A at 0°, B at 180°; each two radiators %d sections (λ/4) apart, %.0f mA each\n\n", quart, Amp*1e3)
	var all []*Result
	for _, s := range States() {
		r, err := Run(s)
		if err != nil {
			log.Fatalf("%s: %v", s.Name, err)
		}
		all = append(all, r)
		var hi, mean float64
		for k, e := range r.Envelope {
			hi = math.Max(hi, e)
			mean += r.Flow[k] / N
		}
		pt := r.PTowers[0] + r.PTowers[1] + r.PTowers[2] + r.PTowers[3]
		fmt.Printf("%-9s %s\n", r.Name, phases(r.Rad))
		fmt.Printf("          peak %.3f V, envelope min/max %.3f, mean eastward power flow %.3g W\n", hi, r.Uniformity(), mean)
		fmt.Printf("          towers deliver %.4g W (A %.3g W, B %.3g W) = ring loss %.4g W; stored %.3g µJ; ledger mismatch %.1e of %.3g J in\n\n",
			pt, r.PTowers[0]+r.PTowers[1], r.PTowers[2]+r.PTowers[3], r.PLoss, r.Stored*1e6, r.Mismatch, r.EnergyIn)
	}

	fmt.Printf("Four towers: A 0°, C 90°, B 180°, D 270° (each pair: east launcher + west launcher)\n\n")
	fmt.Printf("%-24s %8s %8s %8s %8s %9s %10s %10s %10s\n", "state", "|V| 0°", "90°", "180°", "270°", "peak", "towers W", "loss W", "stored µJ")
	for _, s := range FourStates() {
		r, err := Run(s)
		if err != nil {
			log.Fatalf("%s: %v", s.Name, err)
		}
		all = append(all, r)
		var hi, pt float64
		for _, e := range r.Envelope {
			hi = math.Max(hi, e)
		}
		for _, p := range r.PTowers {
			pt += p
		}
		fmt.Printf("%-24s %8.3f %8.3f %8.3f %8.3f %9.3f %10.4g %10.4g %10.3g\n", r.Name,
			r.Envelope[nodeA], r.Envelope[nodeC], r.Envelope[nodeB], r.Envelope[nodeD], hi, pt, r.PLoss, r.Stored*1e6)
		fmt.Printf("%-24s %s\n", "", phases(r.Rad))
		fmt.Printf("%-24s V at the quarter points in step: %s\n", "", inStep(r))
	}

	fmt.Printf("\nEight towers: A 0°, E 45°, C 90°, G 135°, B 180°, F 225°, D 270°, H 315°\n")
	for _, s := range EightStates() {
		r, err := Run(s)
		if err != nil {
			log.Fatal(err)
		}
		var pt float64
		for _, p := range r.PTowers {
			pt += p
		}
		fmt.Printf("  %-36s |V| at 0°…315°:", s.Name)
		for k := 0; k < N; k += N / 8 {
			fmt.Printf(" %5.2f", math.Hypot(r.Phasor[k][0], r.Phasor[k][1]))
		}
		fmt.Printf("  towers %.4g W = loss %.4g W, stored %.3g µJ\n", pt, r.PLoss, r.Stored*1e6)
	}
	fmt.Printf("\nDC bias per pair (synchronized RF): A–B supply 100 V, C–D supply 50 V, 1 MΩ feeds\n")
	for _, sec := range []bool{false, true} {
		b := DefaultBias(100, 50, sec)
		lv, sup, err := BiasOP(b)
		if err != nil {
			log.Fatal(err)
		}
		r, err := RunBiased(FourStates()[2], b)
		if err != nil {
			log.Fatal(err)
		}
		name := "continuous ring"
		if sec {
			name = "sectored ring (10 µF blocking caps at 45°, 135°, 225°, 315°)"
		}
		amp := func(k int) float64 { return math.Hypot(r.Phasor[k][0], r.Phasor[k][1]) }
		iab := sup["A"] + sup["B"]
		var hi float64
		for _, e := range r.Envelope {
			hi = math.Max(hi, e)
		}
		fmt.Printf("  %s\n", name)
		fmt.Printf("    DC at 0°/90°/180°/270°: %.1f / %.1f / %.1f / %.1f V; RF on top: %.2f / %.2f / %.2f / %.2f V; peak %.1f V\n",
			lv[nodeA], lv[nodeC], lv[nodeB], lv[nodeD], amp(nodeA), amp(nodeC), amp(nodeB), amp(nodeD), hi)
		fmt.Printf("    DC current between the pairs' supplies %.3g µA, dissipated in the feeds %.3g mW\n",
			iab*1e6, iab*(b.VAB-b.VCD)*1e3)
	}
	if *out != "" {
		bs, _ := json.Marshal(all)
		if err := os.WriteFile(*out, bs, 0o644); err != nil {
			log.Fatal(err)
		}
	}
}

func phases(rs []Radiator) string {
	out := ""
	for _, r := range rs {
		out += fmt.Sprintf("%s %+4.0f°  ", r.Name, r.Phase)
	}
	return out
}

// inStep reports whether the instantaneous voltages at 0°, 90°, 180° and
// 270° have the same sign at the frame where V(0°) is largest (they peak
// together), with their spread (max − min)/max.
func inStep(r *Result) string {
	best, bf := 0.0, 0
	for i, fr := range r.Frames {
		if math.Abs(fr[nodeA]) > best {
			best, bf = math.Abs(fr[nodeA]), i
		}
	}
	fr := r.Frames[bf]
	vs := []float64{fr[nodeA], fr[nodeC], fr[nodeB], fr[nodeD]}
	ok := best > 0.05
	lo, hi := math.Inf(1), 0.0
	for _, v := range vs {
		if v*vs[0] <= 0 {
			ok = false
		}
		lo, hi = math.Min(lo, math.Abs(v)), math.Max(hi, math.Abs(v))
	}
	return fmt.Sprintf("%v, spread %.0f %% (%+.3f %+.3f %+.3f %+.3f V)", ok, 100*(hi-lo)/hi, vs[0], vs[1], vs[2], vs[3])
}

// LauncherData is one launcher's steady state: the node voltage phasors
// (V(t) = re·cos ωt + im·sin ωt) it produces alone, and its radiators.
type LauncherData struct {
	V   [][2]float64 `json:"V"`
	Rad []struct {
		Node  int     `json:"node"`
		Phase float64 `json:"phase"`
	} `json:"rad"`
}

// Launchers simulates each of the eight launchers alone (A, C, E, G
// east; B, D, F, H west). The ring is linear, so any phasing of the four is the
// sum of these responses rotated by the chosen phase offsets — what the
// animations in viz/ compute in the browser.
func Launchers() (map[string]LauncherData, error) {
	ls := []struct {
		name string
		rad  []Radiator
	}{
		{"A east", East("A", nodeA, 0)}, {"B west", West("B", nodeB, 0)},
		{"C east", East("C", nodeC, 0)}, {"D west", West("D", nodeD, 0)},
		{"E east", East("E", nodeE, 0)}, {"F west", West("F", nodeF, 0)},
		{"G east", East("G", nodeG, 0)}, {"H west", West("H", nodeH, 0)},
	}
	out := map[string]LauncherData{}
	for _, l := range ls {
		r, err := Run(State{l.name, l.rad})
		if err != nil {
			return nil, err
		}
		var d LauncherData
		for _, p := range r.Phasor {
			d.V = append(d.V, [2]float64{math.Round(p[0]*1e4) / 1e4, math.Round(p[1]*1e4) / 1e4})
		}
		for _, x := range l.rad {
			d.Rad = append(d.Rad, struct {
				Node  int     `json:"node"`
				Phase float64 `json:"phase"`
			}{x.Node, x.Phase})
		}
		out[l.name] = d
	}
	return out, nil
}

// WritePhasors writes Launchers() as JSON to path.
func WritePhasors(path string) error {
	d, err := Launchers()
	if err != nil {
		return err
	}
	b, err := json.Marshal(d)
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}
