// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"sync"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/material"
)

// Part 3t: the output winding as a toroidal inductor winding. Every turn
// goes over the outside of the torus, down through the hole and back
// under it, round the tube's cross-section, wound on top of the two Rodin
// windings (over the outside of winding B's bundle). The winding is split
// into TorSections = 12 sections spaced 30° apart round the axis, section
// k centred at θ0 + k·30°, each spanning an arc Δφ with N turns of 0.5 mm
// copper (declared, as the Part 3 coils).
//
// Model (declared): a turn is a circle in the meridian plane at its
// azimuth φ, radius a_m round the tube's centre circle (R0, φ, 0). A
// section's N turns are spread evenly over its arc. When they do not fit
// in one layer on the inside of the hole (the crowded side, circumference
// Δφ·(R0 − a) per layer), further layers are wound over the first, up to
// TorMaxLayers; the turns then sit at the layers' mean radius, spread over
// the arc, and two turns are never closer than touching on the inside
// (the ring–ring mutual is clamped there). The helix's advance round the
// axis within a section is cancelled by its leads, which run back
// together to the section's terminals (checked against an explicit helix
// with its return lead, TorChecks). Leads between sections and to the
// bridges are twisted pairs with no loop area.
//
// Inductances by the segment method (induct.go): every turn's mutual to
// every segment of the Rodin windings is tabulated over 0–30° in 0.25°
// steps and carried round the axis by the windings' 30° symmetry (a 30°
// turn maps pass q of A onto pass q + 5 and pass q of B onto pass q + 7,
// both in the same position along the pass); the turn–turn mutual is
// tabulated against the azimuth difference.

const (
	TorSections  = 12
	TorMaxLayers = 6
	TorRingPts   = 96
	TorStepDeg   = 0.25               // table step for the turn–Rodin mutuals
	TorTable     = 120                // steps over 30°
	TorEnamel    = 0.03e-3            // enamel per wire (as line.go)
	TorMinArc    = 2 * math.Pi / 180  // smallest section arc
	TorMaxArc    = 29 * math.Pi / 180 // largest: 1° left between sections for the leads
	TorWirePitch = 2*OutWireRadius + 2*TorEnamel
)

// TorKit holds the inductance tables of the toroidal secondary on one
// Rodin pair geometry (the same polylines serve every sectioning of the
// pair; only the labels differ).
type TorKit struct {
	P0 *Pair // geometry (any Sec)
	A0 float64
	np int // segments per pass

	mu  sync.Mutex
	seg map[int][][]float64 // layer count → [segment][angle step] mutual of one turn [H]
	rr  map[int]*ringTable  // layer count → turn–turn table
	lab map[[2]int][][]float64
}

type ringTable struct {
	am, delta float64
	x, y      []float64 // log-spaced Δφ and the mutual [H]
	self      float64
}

// NewTorKit prepares the tables' geometry. The first layer's turns lie on
// winding B's outer surface: a0 = OuterShell + r_w + enamel.
func NewTorKit(p *Pair) (*TorKit, error) {
	nA, nB := len(p.A)-1, len(p.B)-1
	if nA%Poles != 0 || nB != nA {
		return nil, fmt.Errorf("pass segments not uniform (%d, %d)", nA, nB)
	}
	k := &TorKit{P0: p, A0: p.OuterShell() + OutWireRadius + TorEnamel, np: nA / Poles,
		seg: map[int][][]float64{}, rr: map[int]*ringTable{}, lab: map[[2]int][][]float64{}}
	// the symmetry the tables rely on: a 30° turn maps segment (q, j) of A
	// to (q+5, j) and of B to (q+7, j)
	rot := func(p P3) P3 {
		c, s := math.Cos(math.Pi/6), math.Sin(math.Pi/6)
		return P3{c*p[0] - s*p[1], s*p[0] + c*p[1], p[2]}
	}
	for q := 0; q < Poles; q++ {
		for j := 0; j < k.np; j++ {
			if d := dist(rot(p.A[q*k.np+j]), p.A[((q+5)%Poles)*k.np+j]); d > 1e-9 {
				return nil, fmt.Errorf("A not 30°-symmetric at pass %d point %d (%g m)", q, j, d)
			}
			if d := dist(rot(p.B[q*k.np+j]), p.B[((q+7)%Poles)*k.np+j]); d > 1e-9 {
				return nil, fmt.Errorf("B not 30°-symmetric at pass %d point %d (%g m)", q, j, d)
			}
		}
	}
	return k, nil
}

// LayerRadius is the mean turn radius of a section of `layers` layers.
func (k *TorKit) LayerRadius(layers int) float64 {
	return k.A0 + float64(layers-1)*TorWirePitch/2
}

// Layers returns how many layers a section of n turns over arc needs
// (filling the inside of the hole layer by layer), 0 when it does not fit
// in TorMaxLayers.
func (k *TorKit) Layers(n int, arc float64) int {
	left := n
	for l := 1; l <= TorMaxLayers; l++ {
		a := k.A0 + float64(l-1)*TorWirePitch
		fit := int(arc * (TorusR0 - a) / TorWirePitch)
		left -= fit
		if left <= 0 {
			return l
		}
	}
	return 0
}

// ring returns a turn at azimuth phi of radius a.
func torRing(phi, a float64) []P3 {
	pts := make([]P3, TorRingPts)
	c, s := math.Cos(phi), math.Sin(phi)
	for i := range pts {
		u := 2 * math.Pi * float64(i) / TorRingPts
		r := TorusR0 + a*math.Cos(u)
		pts[i] = P3{r * c, r * s, a * math.Sin(u)}
	}
	return pts
}

// segTable returns, for a layer count, the mutual of one turn at each
// table angle (0 … 30°) with every segment of A then B.
func (k *TorKit) segTable(layers int) [][]float64 {
	k.mu.Lock()
	if t, ok := k.seg[layers]; ok {
		k.mu.Unlock()
		return t
	}
	k.mu.Unlock()
	am := k.LayerRadius(layers)
	var sg Segs
	sg.AddPolyline(k.P0.A, 0, 1, k.P0.Radius)
	sg.AddPolyline(k.P0.B, 0, 1, k.P0.Radius)
	ns := len(sg.A)
	t := make([][]float64, ns)
	for i := range t {
		t[i] = make([]float64, TorTable)
	}
	parallelFor(TorTable, func(a int) {
		ring := torRing(float64(a)*TorStepDeg*math.Pi/180, am)
		for i := 0; i < ns; i++ {
			v := 0.0
			for q := range ring {
				v += segMutual(sg.A[i], sg.B[i], ring[q], ring[(q+1)%len(ring)], 0)
			}
			t[i][a] = v * em.Mu0 / (4 * math.Pi)
		}
	})
	k.mu.Lock()
	k.seg[layers] = t
	k.mu.Unlock()
	return t
}

// Labels returns, for a model's sectioning and a layer count, the mutual
// of one turn with every pass label of the model (0 … 2·Sec−1) at every
// angle round the axis in TorStepDeg steps (12·TorTable columns).
func (k *TorKit) Labels(m *Model, layers int) [][]float64 {
	key := [2]int{m.P.Sec, layers}
	k.mu.Lock()
	if t, ok := k.lab[key]; ok {
		k.mu.Unlock()
		return t
	}
	k.mu.Unlock()
	t := k.segTable(layers)
	sec := m.P.Sec
	n := Poles * k.np // segments per circuit
	out := make([][]float64, 2*sec)
	for i := range out {
		out[i] = make([]float64, Poles*TorTable)
	}
	for w := 0; w < 2; w++ {
		step := 5
		if w == 1 {
			step = 7
		}
		for i := 0; i < n; i++ {
			label := w*sec + i*sec/n
			q, j := i/k.np, i%k.np
			for m30 := 0; m30 < Poles; m30++ {
				src := w*n + ((q-step*m30)%Poles+Poles)%Poles*k.np + j
				row := t[src]
				dst := out[label][m30*TorTable : (m30+1)*TorTable]
				for a := range dst {
					dst[a] += row[a]
				}
			}
		}
	}
	k.mu.Lock()
	k.lab[key] = out
	k.mu.Unlock()
	return out
}

// interp reads a periodic label table at azimuth phi [rad].
func interpPhi(row []float64, phi float64) float64 {
	n := len(row)
	x := phi / (TorStepDeg * math.Pi / 180)
	x = math.Mod(x, float64(n))
	if x < 0 {
		x += float64(n)
	}
	i := int(x)
	f := x - float64(i)
	return row[i%n]*(1-f) + row[(i+1)%n]*f
}

// ringRR returns the turn–turn table for a layer count.
func (k *TorKit) ringRR(layers int) *ringTable {
	k.mu.Lock()
	if t, ok := k.rr[layers]; ok {
		k.mu.Unlock()
		return t
	}
	k.mu.Unlock()
	am := k.LayerRadius(layers)
	rt := &ringTable{am: am, delta: TorWirePitch / (TorusR0 - am), self: CircleRing(am, OutWireRadius)}
	const npt = 240
	rt.x = make([]float64, npt)
	rt.y = make([]float64, npt)
	r0 := torRing(0, am)
	parallelFor(npt, func(i int) {
		x := rt.delta * math.Pow(math.Pi/rt.delta, float64(i)/float64(npt-1))
		r1 := torRing(x, am)
		v := 0.0
		for p := range r0 {
			for q := range r1 {
				v += segMutual(r0[p], r0[(p+1)%len(r0)], r1[q], r1[(q+1)%len(r1)], 0)
			}
		}
		rt.x[i], rt.y[i] = x, v*em.Mu0/(4*math.Pi)
	})
	k.mu.Lock()
	k.rr[layers] = rt
	k.mu.Unlock()
	return rt
}

// M returns the mutual of two turns dphi apart (clamped at touching).
func (rt *ringTable) M(dphi float64) float64 {
	d := math.Abs(math.Remainder(dphi, 2*math.Pi))
	if d <= rt.x[0] {
		return rt.y[0]
	}
	if d >= rt.x[len(rt.x)-1] {
		return rt.y[len(rt.y)-1]
	}
	u := math.Log(d/rt.x[0]) / math.Log(rt.x[len(rt.x)-1]/rt.x[0]) * float64(len(rt.x)-1)
	i := min(int(u), len(rt.x)-2)
	f := u - float64(i)
	return rt.y[i]*(1-f) + rt.y[i+1]*f
}

// TorLayout is a toroidal secondary design.
type TorLayout struct {
	F      float64
	Sense  int     // primary connection (Toroidal / Dipole)
	Theta0 float64 // centre of section 0 [rad]; pole 1 is at φ = 0, the hole crossings at 15° + k·30°
	Arc    float64 // arc of each section [rad]
	Turns  int     // per section
	Alt    bool    // alternate the winding sense from section to section
	Comb   Combine
	RLoad  float64
	PhaseB float64 // shift of B's source from the cosine [deg]: 0 = sine on A, cosine on B; −90 = in phase
}

// Angles of the turns of section s.
func (l TorLayout) turnAngles(s int) []float64 {
	c := l.Theta0 + float64(s)*2*math.Pi/TorSections
	out := make([]float64, l.Turns)
	for j := range out {
		out[j] = c - l.Arc/2 + l.Arc*(float64(j)+0.5)/float64(l.Turns)
	}
	return out
}

func (l TorLayout) sign(s int) float64 {
	if l.Alt && s%2 == 1 {
		return -1
	}
	return 1
}

// SectionWire is the wire length of one section [m].
func (k *TorKit) SectionWire(l TorLayout) float64 {
	return float64(l.Turns) * 2 * math.Pi * k.LayerRadius(max(1, k.Layers(l.Turns, l.Arc)))
}

// torResistance of a length of the 0.5 mm wire at f (skin effect as
// GroupShape.Resistance).
func torResistance(length, f float64) float64 {
	area := math.Pi * OutWireRadius * OutWireRadius
	rdc := material.Copper.Resistivity * length / area
	delta := math.Sqrt(material.Copper.Resistivity / (math.Pi * f * em.Mu0))
	return rdc * math.Max(1, 2*OutWireRadius/(4*delta))
}

// Feasible returns 0 or a negative violation measure: layers, arc, lumped
// wire (λ/20: each section's wire, the whole winding's when in series),
// tuning capacitor ≥ MinTuneC.
func (k *TorKit) Feasible(l TorLayout) float64 {
	worst := 0.0
	if l.Arc < TorMinArc-1e-12 || l.Arc > TorMaxArc+1e-12 {
		worst = math.Min(worst, -1e-3)
	}
	ly := k.Layers(l.Turns, l.Arc)
	if ly == 0 {
		return math.Min(worst, -1e-3*float64(l.Turns)/float64(max(1, k.Layers(1, l.Arc))))
	}
	lw := k.SectionWire(l)
	if l.Comb == Series {
		lw *= TorSections
	}
	if d := em.SpeedOfLight/(20*l.F) - lw; d < 0 {
		worst = math.Min(worst, d*1e-2)
	}
	_, mg := k.inductances(nil, l, ly)
	leff := 0.0
	for j := range mg[0] {
		leff += mg[0][j]
	}
	if leff <= 0 {
		return math.Min(worst, -1e-3)
	}
	ceq := 1 / (math.Pow(2*math.Pi*l.F, 2) * leff)
	if l.Comb == Series {
		ceq /= TorSections
	}
	if ceq < MinTuneC {
		worst = math.Min(worst, -1e-3*math.Log10(MinTuneC/ceq)-1e-6)
	}
	return worst
}

// inductances returns each section's mutual to the pass labels of model
// m (nil m: skip) and the section matrix (signed for the polarity).
func (k *TorKit) inductances(m *Model, l TorLayout, layers int) ([][]float64, [][]float64) {
	rt := k.ringRR(layers)
	// turn–turn: row 0 against every section (the layout is 12-fold
	// symmetric apart from the polarity signs)
	a0 := l.turnAngles(0)
	row := make([]float64, TorSections)
	for s := 0; s < TorSections; s++ {
		as := l.turnAngles(s)
		v := 0.0
		for i, x := range a0 {
			for j, y := range as {
				if s == 0 && i == j {
					v += rt.self
					continue
				}
				v += rt.M(y - x)
			}
		}
		row[s] = v
	}
	mg := newMat(TorSections)
	for i := 0; i < TorSections; i++ {
		for j := 0; j < TorSections; j++ {
			mg[i][j] = row[(j-i+TorSections)%TorSections] * l.sign(i) * l.sign(j)
		}
	}
	if m == nil {
		return nil, mg
	}
	lab := k.Labels(m, layers)
	pass := make([][]float64, TorSections)
	for s := range pass {
		pass[s] = make([]float64, len(lab))
		for _, phi := range l.turnAngles(s) {
			for q := range lab {
				pass[s][q] += interpPhi(lab[q], phi)
			}
		}
		for q := range pass[s] {
			pass[s][q] *= l.sign(s)
		}
	}
	return pass, mg
}

// Groups builds the 12 sections as output groups on a net of model m:
// each tuned by a series capacitor to its effective inductance with all
// sections carrying equal currents (L_k + Σ M_kj), exact for the series
// string and for the symmetric layout.
func (k *TorKit) Groups(m *Model, n *Net, l TorLayout) ([]Group, [][]float64, error) {
	ly := k.Layers(l.Turns, l.Arc)
	if ly == 0 {
		return nil, nil, fmt.Errorf("%d turns do not fit in %d layers over %.1f°", l.Turns, TorMaxLayers, l.Arc*180/math.Pi)
	}
	pass, mg := k.inductances(m, l, ly)
	w := 2 * math.Pi * l.F
	r := torResistance(k.SectionWire(l), l.F)
	gs := make([]Group, TorSections)
	for s := range gs {
		leff := 0.0
		for j := range mg[s] {
			leff += mg[s][j]
		}
		gs[s] = Group{L: mg[s][s], R: r, Ctune: 1 / (w * w * leff), MSec: n.groupMSec(m, pass[s])}
	}
	return gs, mg, nil
}

// TorPart holds the driven nets of a model.
type TorPart struct {
	K    *TorKit
	M    *Model
	Nets map[int]*Net
	mu   sync.Mutex
}

// NewTorPart builds the driven nets (both senses) on model m.
func NewTorPart(k *TorKit, m *Model) *TorPart {
	return &TorPart{K: k, M: m, Nets: map[int]*Net{Toroidal: m.NewNet(DrivenPair, Toroidal, 1e5), Dipole: m.NewNet(DrivenPair, Dipole, 1e5)}}
}

func (p *TorPart) net(sense int) *Net {
	p.mu.Lock()
	defer p.mu.Unlock()
	src := p.Nets[sense]
	c := *src
	c.R = append([]float64(nil), src.R...)
	return &c
}

// Score solves a layout with the sin (A) / cos (B) drive; phaseB shifts
// B's source from the cosine (0 = cosine, −90 = in phase with A).
func (p *TorPart) Score(l TorLayout, phaseB float64) (float64, *ACResult, error) {
	if c := p.K.Feasible(l); c < 0 {
		return 0, nil, fmt.Errorf("infeasible (%.3g)", c)
	}
	n := p.net(l.Sense)
	n.SetFreq(l.F)
	gs, mg, err := p.K.Groups(p.M, n, l)
	if err != nil {
		return 0, nil, err
	}
	src := n.DriveSources(l.Sense, Part3Volts)
	src[1].PhaseDeg += l.PhaseB + phaseB
	r, err := n.SolveAC(l.F, src, gs, mg, l.Comb, l.RLoad)
	if err != nil {
		return 0, nil, err
	}
	return r.PLoad, r, nil
}

// TorSpace is the genome; with phase, the drive phase of B is a gene too.
func TorSpace(phase bool) evo.Space {
	sp := evo.Space{
		{Name: "f [Hz]", Min: 1e3, Max: 3e7},
		{Name: "θ0 (log10 x · 30°)", Min: 1, Max: 10},
		{Name: "arc (2° + log10 x · 27°)", Min: 1, Max: 10},
		{Name: "turns per section", Min: 1, Max: 40, Integer: true},
		{Name: "polarity (1 all alike, 2 alternating)", Min: 1, Max: 2, Integer: true},
		{Name: "combine (1 separate, 2 series, 3 parallel)", Min: 1, Max: 3, Integer: true},
		{Name: "R_L [Ω]", Min: 0.1, Max: 1e4},
		{Name: "sense (1 toroidal, 2 dipole)", Min: 1, Max: 2, Integer: true},
	}
	if phase {
		sp = append(sp, evo.Gene{Name: "B's shift from the cosine (log10 x · 360°)", Min: 1, Max: 10})
	}
	return sp
}

// DecodeTor turns gene values into a layout.
func DecodeTor(v []float64) TorLayout {
	lin := math.Log10
	l := TorLayout{F: v[0], Theta0: math.Pi / 6 * lin(v[1]), Arc: TorMinArc + (TorMaxArc-TorMinArc)*lin(v[2]),
		Turns: int(v[3]), Alt: int(v[4]) == 2, Comb: Combine(int(v[5])), RLoad: v[6], Sense: Toroidal}
	if int(v[7]) == 2 {
		l.Sense = Dipole
	}
	if len(v) > 8 {
		l.PhaseB = 360 * lin(v[8])
	}
	return l
}

// Evolve runs the genetic search (fitness as Part 3: 6 − log10 of the
// load power times the bridge factor; infeasible layouts 100 + 10³·violation).
func (p *TorPart) Evolve(o GAOptions, phase bool) (*TorGAResult, error) {
	dummy := func([]float64) (*circuit.Circuit, error) {
		return circuit.NewBuilder("placeholder").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "0", 1).Build()
	}
	prob := evo.Problem{Build: dummy, Space: TorSpace(phase), Objectives: []evo.Objective{{Name: "6 − log10(P_load/W)", Weight: 1,
		Eval: func(e *evo.Eval) (float64, error) {
			l := DecodeTor(e.Values)
			if c := p.K.Feasible(l); c < 0 {
				return 100 + 1e3*(-c), nil
			}
			pw, r, err := p.Score(l, 0)
			if err != nil {
				return 0, err
			}
			return 6 - math.Log10(pw*BridgeFactor(r)+1e-30), nil
		}}}}
	st, err := evo.Run(context.Background(), prob, evo.Options{Population: o.Population, Generations: o.Generations, Seed: o.Seed, Timeout: 600e9}, nil, nil)
	if err != nil {
		return nil, err
	}
	if st.Best.Failure != "" {
		return nil, fmt.Errorf("best failed: %s %s", st.Best.Failure, st.Best.Error)
	}
	best := DecodeTor(st.Best.Values)
	pw, ac, err := p.Score(best, 0)
	if err != nil {
		return nil, err
	}
	return &TorGAResult{Best: best, Power: pw, AC: ac, History: st.History, Evals: st.Evaluations}, nil
}

// TorGAResult is the outcome of a search.
type TorGAResult struct {
	Best    TorLayout
	Power   float64
	AC      *ACResult
	History []evo.GenStats
	Evals   int
}

func (l TorLayout) String() string {
	s := map[int]string{Toroidal: "toroidal", Dipole: "dipole"}[l.Sense]
	pol := "all sections wound alike"
	if l.Alt {
		pol = "alternate sections reversed"
	}
	return fmt.Sprintf("f = %.4g Hz, primary sense %s, B's source %.1f° ahead of A's (90° = sine/cosine), 12 sections centred at %.2f° + k·30°, arc %.2f°, %d turns each, %s, %s, R_L = %.3g Ω",
		l.F, s, math.Mod(90+l.PhaseB+360, 360), l.Theta0*180/math.Pi, l.Arc*180/math.Pi, l.Turns, pol, l.Comb, l.RLoad)
}

// TurnEMF returns, for one turn at each azimuth, its open-circuit EMF
// phasor [V] on the unloaded driven net of model m (sources as Score,
// phaseB shift on B; onlyA / onlyB switch the other source off).
func (k *TorKit) TurnEMF(m *Model, n *Net, f float64, sense int, phaseB float64, which int, phis []float64) ([]complex128, error) {
	n.SetFreq(f)
	src := n.DriveSources(sense, Part3Volts)
	src[1].PhaseDeg += phaseB
	switch which {
	case 1:
		src[1].V = 0
	case 2:
		src[0].V = 0
	}
	r, err := n.SolveAC(f, src, nil, nil, Separate, 1)
	if err != nil {
		return nil, err
	}
	ip := make([]complex128, 2*m.P.Sec)
	for g, I := range r.I {
		ip[(g/n.NS)*m.P.Sec+g%m.P.Sec] += I
	}
	lab := k.Labels(m, 1)
	out := make([]complex128, len(phis))
	for i, phi := range phis {
		var e complex128
		for q, I := range ip {
			e += complex(0, 2*math.Pi*f*interpPhi(lab[q], phi)) * I
		}
		out[i] = e
	}
	return out, nil
}

// TorChecks validates the tables: (1) interpolated, rotated turn–pass
// mutuals against a direct computation at angles off the table grid;
// (2) a closely wound uniform toroid of nTor turns against the current-
// sheet formula μ0N²(R0 − √(R0² − a²)); (3) one section of nh turns over
// arcDeg as rings against an explicit helix closed by a return lead along
// the outer equator: mutual to A and to B, and self inductance.
type TorCheck struct {
	TableErr                     float64 // largest |interp − direct| / largest |direct|
	TorRings, TorSheet           float64 // H
	NTor                         int
	RingMA, RingMB, HelMA, HelMB float64 // H (one circuit's path)
	RingL, HelL                  float64
}

func (k *TorKit) Checks(m *Model, nTor, nh int, arcDeg float64) TorCheck {
	var c TorCheck
	lab := k.Labels(m, 1)
	am := k.LayerRadius(1)
	worst, big := 0.0, 0.0
	for _, deg := range []float64{7.3, 133.71, 251.06, 359.9} {
		phi := deg * math.Pi / 180
		var sg Segs
		sg.AddPolyline(m.P.A, 0, m.P.Sec, m.P.Radius)
		sg.AddPolyline(m.P.B, m.P.Sec, m.P.Sec, m.P.Radius)
		direct := sg.MutualTo(torRing(phi, am), 1, 2*m.P.Sec)
		for q := range direct {
			worst = math.Max(worst, math.Abs(interpPhi(lab[q], phi)-direct[q]))
			big = math.Max(big, math.Abs(direct[q]))
		}
	}
	c.TableErr = worst / big
	// uniform toroid
	rt := k.ringRR(1)
	c.NTor = nTor
	for i := 0; i < nTor; i++ {
		for j := 0; j < nTor; j++ {
			if i == j {
				c.TorRings += rt.self
			} else {
				c.TorRings += rt.M(2 * math.Pi * float64(j-i) / float64(nTor))
			}
		}
	}
	c.TorSheet = em.Mu0 * float64(nTor*nTor) * (TorusR0 - math.Sqrt(TorusR0*TorusR0-am*am))
	// one section: rings vs helix + return lead
	arc := arcDeg * math.Pi / 180
	l := TorLayout{Theta0: 0, Arc: arc, Turns: nh}
	pass, mg := k.inductances(m, l, 1)
	for q := 0; q < m.P.Sec; q++ {
		c.RingMA += pass[0][q]
		c.RingMB += pass[0][m.P.Sec+q]
	}
	c.RingL = mg[0][0]
	hel := torHelix(-arc/2, arc, nh, am)
	var sg Segs
	sg.AddPolyline(m.P.A, 0, 1, m.P.Radius)
	sg.AddPolyline(m.P.B, 1, 1, m.P.Radius)
	mm := sg.MutualTo(hel, 1, 2)
	c.HelMA, c.HelMB = mm[0], mm[1]
	var sh Segs
	sh.AddPolyline(append(append([]P3(nil), hel...), hel[0]), 0, 1, OutWireRadius)
	c.HelL = sh.SectionMatrix(1)[0][0]
	return c
}

// torHelix is a section wound as a helix of n turns from azimuth phi0
// over arc, starting and ending on the outer equator, followed by its
// return lead back along the outer equator (closed polygon).
func torHelix(phi0, arc float64, n int, a float64) []P3 {
	pts := []P3{}
	tot := n * TorRingPts
	for i := 0; i <= tot; i++ {
		u := float64(i) / float64(tot)
		phi := phi0 + arc*u
		th := 2 * math.Pi * float64(n) * u
		r := TorusR0 + a*math.Cos(th)
		pts = append(pts, P3{r * math.Cos(phi), r * math.Sin(phi), a * math.Sin(th)})
	}
	// return lead along the outer equator, just outside the turns
	rl := TorusR0 + a + TorWirePitch
	for i := 16; i >= 0; i-- {
		phi := phi0 + arc*float64(i)/16
		pts = append(pts, P3{rl * math.Cos(phi), rl * math.Sin(phi), 0})
	}
	return pts
}
