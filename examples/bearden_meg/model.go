// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// Magnetic circuits are written in kronspice's magnetic domain (exchange
// channel "magnetic": across = magnetomotive force [A], through = flux rate
// [V]), the gyrator-capacitor model (D. C. Hamill, "Lumped equivalent
// circuits of magnetic components: the gyrator-capacitor approach," IEEE
// Trans. Power Electron. 8(2), 97–103, 1993):
//
//   - a core leg is a capacitor whose charge is the flux, Φ = A·B(MMF/l);
//   - a winding of N turns is a gyrator with g = 1/N (MMF = N·i, v = N·dΦ/dt);
//   - a permanent magnet is an MMF source Fc = Hc·lm in series with its
//     recoil permeance μ0·μr·A/lm.
//
// Every element is lossless except the resistors, so the ledger books the
// magnets' energy like any other source.

// Core material: the B-H curve of the ngspice drafts (soft iron, no
// hysteresis), extended with the slope of vacuum beyond 2 kA/m.
var (
	coreH = []float64{-1e5, -2000, -1000, -400, -200, -100, 0, 100, 200, 400, 1000, 2000, 1e5}
	coreB = []float64{-1.6 - mu0*98000, -1.6, -1.5, -1.3, -1.0, -0.6, 0, 0.6, 1.0, 1.3, 1.5, 1.6, 1.6 + mu0*98000}
)

const (
	mu0     = 4e-7 * math.Pi
	legArea = 1e-4 // core leg cross-section [m²]
	legLen  = 0.05 // core leg magnetic length [m]
	magHc   = 900e3
	magLen  = 5e-3 // NdFeB, 5 mm long, same area as the leg
	magMuR  = 1.05
	leak    = 1e12 // MMF-node leak to ground [Ω analogue]: defines the DC point
	tOn     = 1e-3
)

func magnetF() float64 { return magHc * magLen }
func magnetP() float64 { return mu0 * magMuR * legArea / magLen }

// Kind selects the experiment.
type Kind int

const (
	// MEG: magnet in the centre leg, two side legs, each with a control
	// and an output winding; the control windings are pulsed alternately
	// (US 6,362,718 arrangement).
	MEG Kind = iota
	// Knife: two matched magnets in opposition round the outer loop, so the
	// centre leg sits near zero flux; a drive and an output winding share
	// the centre leg.
	Knife
)

func (k Kind) String() string { return [...]string{"meg", "knife"}[k] }

// Spec is one run. Zero fields take the defaults in withDefaults.
type Spec struct {
	Name     string
	Kind     Kind
	V        float64 // drive amplitude [V] (MEG: pulse height)
	F        float64 // drive frequency [Hz]
	RL       float64 // load per output winding [Ω]
	Rs       float64 // drive source resistance [Ω]
	Nd, No   int     // drive turns (Knife: total of the two series coils) / output turns
	Mismatch float64 // Knife: magnet B strength relative to A
	Gap      float64 // Knife: iron between the magnets' N faces [m] (4 equal bar segments)
	NullLeak float64 // Knife: air-leakage permeance from the null to the yokes [Wb/A] (0: none)
	Cres     float64 // Knife: series-resonant drive capacitor [F] (0: none)
	Duty     float64 // MEG: pulse duty per coil
	Step     float64 // max time step [s]
	Cycles   int     // measured drive periods
}

func (s Spec) withDefaults() Spec {
	d := func(p *float64, v float64) {
		if *p == 0 {
			*p = v
		}
	}
	if s.Kind == MEG {
		d(&s.V, 5)
	} else {
		d(&s.V, 2)
	}
	d(&s.F, 10e3)
	d(&s.RL, 100)
	d(&s.Rs, 1)
	d(&s.Mismatch, 1)
	d(&s.Duty, 0.45)
	d(&s.Step, 0.1e-6)
	d(&s.Gap, 0.15)
	if s.Nd == 0 {
		s.Nd = 50
	}
	if s.No == 0 {
		s.No = 200
	}
	if s.Cycles == 0 {
		s.Cycles = 10
	}
	return s
}

// Window is the measured interval: whole drive periods after the magnets
// have ramped in and the drive has settled.
func (s Spec) Window() (t0, t1 float64) {
	s = s.withDefaults()
	settle := 1e-3
	if s.Kind == Knife {
		settle = 10e-3 // the all-iron drive path has L/R ≈ 3 ms
	}
	t1 = tOn + settle + float64(s.Cycles)/s.F
	return t1 - float64(s.Cycles)/s.F, t1
}

// element is one netlist element in a neutral form, written out in either
// dialect so the two simulators run the same circuit.
type element struct {
	kind  byte // 'V' 'R' 'C' 'Q' (flux-defined core) 'Y' (gyrator)
	name  string
	nodes []string
	value float64 // 'Q': magnetic length [m] (0: legLen)
	wave  string  // source waveform, SPICE syntax
}

type role struct {
	drives, magnets, loads []string // source and load element names
	probes                 []string // nodes compared between simulators
}

func (s Spec) elements() ([]element, role, []string) {
	s = s.withDefaults()
	F, P := magnetF(), magnetP()
	ramp := fmt.Sprintf("PWL(0 0 %g %g)", tOn/2, F)
	var e []element
	add := func(x element) { e = append(e, x) }
	var r role
	var mag []string
	switch s.Kind {
	case MEG:
		T := 1 / s.F
		mag = []string{"m1", "t", "l1", "l2", "r1", "r2"}
		add(element{'V', "VM", []string{"m1", "0"}, 0, ramp})
		add(element{'C', "CM", []string{"m1", "t"}, P, ""})
		add(element{'Y', "YCL", []string{"cl", "0", "t", "l1"}, 1 / float64(s.Nd), ""})
		add(element{'Y', "YOL", []string{"ol", "0", "l1", "l2"}, 1 / float64(s.No), ""})
		add(element{'Q', "CL", []string{"l2", "0"}, 0, ""})
		add(element{'Y', "YCR", []string{"cr", "0", "t", "r1"}, 1 / float64(s.Nd), ""})
		add(element{'Y', "YOR", []string{"or", "0", "r1", "r2"}, 1 / float64(s.No), ""})
		add(element{'Q', "CR", []string{"r2", "0"}, 0, ""})
		for _, n := range []string{"t", "l1", "l2", "r1", "r2"} {
			add(element{'R', "RK" + strings.ToUpper(n), []string{n, "0"}, leak, ""})
		}
		pulse := func(delay float64) string {
			return fmt.Sprintf("PULSE(0 %g %g 1u 1u %g %g)", s.V, delay, s.Duty*T, T)
		}
		add(element{'V', "VDL", []string{"dl", "0"}, 0, pulse(tOn)})
		add(element{'V', "VDR", []string{"dr", "0"}, 0, pulse(tOn + T/2)})
		add(element{'R', "RSL", []string{"dl", "cl"}, s.Rs, ""})
		add(element{'R', "RSR", []string{"dr", "cr"}, s.Rs, ""})
		add(element{'R', "RLL", []string{"ol", "0"}, s.RL, ""})
		add(element{'R', "RLR", []string{"or", "0"}, s.RL, ""})
		r = role{[]string{"VDL", "VDR"}, []string{"VM"}, []string{"RLL", "RLR"}, []string{"ol", "or", "l2", "r2"}}
	case Knife:
		// Bucking magnets on one bar, N facing N:
		//   yoke ─ [S A N] ─ jL ═ drive, output ═ jR ─ [N B S] ─ yoke
		// Each magnet's flux returns from its junction (jL, jR) to its own
		// S end through its own iron yoke (CRL, CRR); the bridge between
		// the junctions carries the drive and output windings and no flux
		// when the magnets match. The drive's own flux path (bridge, CRR,
		// CRL) is all iron.
		//
		// The bar between the N faces is four equal segments of Gap/4:
		//   A ─ CYA ─ jl ─ drive 1 ─ CB1 ─ n (null, output) ─ CB2 ─ drive 2 ─ jr ─ CYB ─ B
		// The two drive coils (Nd/2 turns each) are in series, wound to
		// aid along the bar. NullLeak, when set, is the air path from the
		// null to the yokes; without it every bridge segment carries the
		// same flux and the two coils act as one of Nd turns.
		seg := s.Gap / 4
		half := 1 / (float64(s.Nd) / 2)
		mag = []string{"a1", "al", "jl", "y1", "n", "y2", "y3", "jr", "br", "b1"}
		add(element{'V', "VMA", []string{"a1", "0"}, 0, ramp})
		add(element{'C', "CMA", []string{"a1", "al"}, P, ""})
		add(element{'Q', "CYA", []string{"al", "jl"}, seg, ""})
		add(element{'V', "VMB", []string{"b1", "0"}, 0, fmt.Sprintf("PWL(0 0 %g %g)", tOn/2, F*s.Mismatch)})
		add(element{'C', "CMB", []string{"b1", "br"}, P, ""})
		add(element{'Q', "CYB", []string{"br", "jr"}, seg, ""})
		add(element{'Y', "YDR1", []string{"d1", "dm", "jl", "y1"}, half, ""})
		add(element{'Q', "CB1", []string{"y1", "n"}, seg, ""})
		add(element{'Y', "YOUT", []string{"o1", "0", "n", "y2"}, 1 / float64(s.No), ""})
		add(element{'Q', "CB2", []string{"y2", "y3"}, seg, ""})
		add(element{'Y', "YDR2", []string{"dm", "0", "y3", "jr"}, half, ""})
		add(element{'R', "RDM", []string{"dm", "0"}, 1e6, ""}) // DC path for the coils' midpoint
		add(element{'Q', "CRL", []string{"jl", "0"}, 0, ""})
		add(element{'Q', "CRR", []string{"jr", "0"}, 0, ""})
		if s.NullLeak > 0 {
			add(element{'C', "CLK", []string{"n", "0"}, s.NullLeak, ""})
		}
		for _, n := range []string{"al", "jl", "y1", "n", "y2", "y3", "jr", "br"} {
			add(element{'R', "RK" + strings.ToUpper(n), []string{n, "0"}, leak, ""})
		}
		add(element{'V', "VD", []string{"ds", "0"}, 0, fmt.Sprintf("SIN(0 %g %g 0)", s.V, s.F)}) // starts with the magnets; 10 ms to settle
		if s.Cres > 0 {                                                                          // series resonance: source → Rs → C → winding
			add(element{'R', "RS", []string{"ds", "dx"}, s.Rs, ""})
			add(element{'C', "CRES", []string{"dx", "d1"}, s.Cres, ""})
			add(element{'R', "RDC", []string{"d1", "0"}, 1e6, ""})
		} else {
			add(element{'R', "RS", []string{"ds", "d1"}, s.Rs, ""})
		}
		add(element{'R', "RL", []string{"o1", "0"}, s.RL, ""})
		r = role{[]string{"VD"}, []string{"VMA", "VMB"}, []string{"RL"}, []string{"o1", "d1", "dm", "jl", "n", "jr"}}
	}
	return e, r, mag
}

// curve returns the flux-versus-MMF knots of a core segment of length l
// (0: legLen).
func curve(l float64) (mmf, flux []float64) {
	if l == 0 {
		l = legLen
	}
	for i := range coreH {
		mmf = append(mmf, coreH[i]*l)
		flux = append(flux, coreB[i]*legArea)
	}
	return
}

func pairs(sep string, l float64) string {
	x, y := curve(l)
	var b strings.Builder
	for i := range x {
		if i > 0 {
			b.WriteString(sep)
		}
		fmt.Fprintf(&b, "%.10g%s%.10g", x[i], sep, y[i])
	}
	return b.String()
}

func (s Spec) title() string {
	s = s.withDefaults()
	t := fmt.Sprintf("bearden_meg %s V=%g f=%g RL=%g", s.Kind, s.V, s.F, s.RL)
	if s.Kind == Knife {
		t += fmt.Sprintf(" mismatch=%g Cres=%g gap=%g leak=%g", s.Mismatch, s.Cres, s.Gap, s.NullLeak)
	}
	return t
}

// Kron writes the kronspice netlist.
func (s Spec) Kron() string {
	s = s.withDefaults()
	e, _, mag := s.elements()
	_, t1 := s.Window()
	var b strings.Builder
	fmt.Fprintln(&b, s.title())
	fmt.Fprintln(&b, "* magnetic domain: node voltage = MMF [A], branch current = flux rate [V]")
	fmt.Fprintln(&b, ".domain magnetic", strings.Join(mag, " "))
	for _, x := range e {
		switch x.kind {
		case 'V':
			fmt.Fprintf(&b, "%s %s %s\n", x.name, strings.Join(x.nodes, " "), x.wave)
		case 'R', 'C':
			fmt.Fprintf(&b, "%s %s %.10g\n", x.name, strings.Join(x.nodes, " "), x.value)
		case 'Q':
			fmt.Fprintf(&b, "%s %s Q=pwl(%s)\n", x.name, strings.Join(x.nodes, " "), pairs(" ", x.value))
		case 'Y':
			fmt.Fprintf(&b, "%s %s %.10g\n", x.name, strings.Join(x.nodes, " "), x.value)
		}
	}
	fmt.Fprintf(&b, ".tran %g %.10g\n.end\n", s.Step, t1)
	return b.String()
}

// Ngspice writes the ngspice twin of Kron: each gyrator as two VCCS
// (i1 = g·V2 out of p1+, i2 = −g·V1 out of p2+), each core as a B source
// carrying d/dt of its flux curve. Its .control block writes the columns
// named on the "* columns:" line, at a fixed step, to NAME.dat.
func (s Spec) Ngspice(name string) string {
	s = s.withDefaults()
	e, r, _ := s.elements()
	_, t1 := s.Window()
	var b strings.Builder
	fmt.Fprintln(&b, s.title()+" (ngspice twin)")
	cols := s.columns(r)
	fmt.Fprintln(&b, "* columns:", strings.Join(cols, " "))
	for _, x := range e {
		n := x.nodes
		switch x.kind {
		case 'V':
			fmt.Fprintf(&b, "%s %s %s %s\n", x.name, n[0], n[1], x.wave)
		case 'R', 'C':
			fmt.Fprintf(&b, "%s %s %s %.10g\n", x.name, n[0], n[1], x.value)
		case 'Q':
			fmt.Fprintf(&b, "B%s %s %s I=ddt(pwl(v(%s,%s), %s))\n", x.name, n[0], n[1], n[0], n[1], pairs(", ", x.value))
		case 'Y':
			fmt.Fprintf(&b, "G%s_1 %s %s %s %s %.10g\n", x.name, n[0], n[1], n[2], n[3], x.value)
			fmt.Fprintf(&b, "G%s_2 %s %s %s %s %.10g\n", x.name, n[2], n[3], n[0], n[1], -x.value)
		}
	}
	fmt.Fprintln(&b, ".options reltol=1e-4 abstol=1e-12 vntol=1e-9 chgtol=1e-16 method=gear")
	fmt.Fprintf(&b, ".tran %g %.10g 0 %g\n", s.Step, t1, s.Step)
	fmt.Fprintln(&b, ".control")
	fmt.Fprintln(&b, "run")
	fmt.Fprintf(&b, "linearize\n")
	var vec []string
	for _, c := range cols {
		vec = append(vec, c)
	}
	fmt.Fprintf(&b, "wrdata %s.dat %s\n", name, strings.Join(vec, " "))
	fmt.Fprintln(&b, "quit")
	fmt.Fprintln(&b, ".endc")
	fmt.Fprintln(&b, ".end")
	return b.String()
}

// columns lists what both simulators must report for the energy split:
// every source's voltage nodes and branch current, every load's node.
func (s Spec) columns(r role) []string {
	e, _, _ := s.elements()
	byName := map[string]element{}
	for _, x := range e {
		byName[x.name] = x
	}
	seen := map[string]bool{}
	var cols []string
	addV := func(n string) {
		if n != "0" && !seen[n] {
			seen[n] = true
			cols = append(cols, "v("+n+")")
		}
	}
	for _, g := range [][]string{r.drives, r.magnets} {
		for _, n := range g {
			addV(byName[n].nodes[0])
			addV(byName[n].nodes[1])
			cols = append(cols, "i("+strings.ToLower(n)+")")
		}
	}
	for _, n := range r.loads {
		addV(byName[n].nodes[0])
	}
	for _, n := range r.probes {
		addV(n)
	}
	return cols
}

// Result is the energy split over the window, as average powers.
type Result struct {
	Pdrive, Pmagnet, Pload float64 // [W]
	Pdiss                  float64 // all resistors [W] (kronspice ledger only)
	DeltaStored, Residual  float64 // [J] (kronspice ledger only)
	COP                    float64 // Pload / Pdrive
}

// RelResidual is the ledger residual against the drive energy.
func (r Result) RelResidual(s Spec) float64 {
	t0, t1 := s.Window()
	return r.Residual / (r.Pdrive * (t1 - t0))
}

// Run simulates s in kronspice and returns its energy split and waveforms.
func (s Spec) Run() (Result, *output.Waveforms, error) {
	s = s.withDefaults()
	deck, err := netlist.ParseString(s.Kron())
	if err != nil {
		return Result{}, nil, err
	}
	c, err := circuit.FromDeck(deck)
	if err != nil {
		return Result{}, nil, err
	}
	a := deck.Analyses[0]
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: a.Step, Stop: a.Stop, MaxStep: a.MaxStep}, analysis.TranOptions{})
	if err != nil {
		return Result{}, nil, err
	}
	_, r, _ := s.elements()
	t0, t1 := s.Window()
	idx := map[string]int{}
	for i, n := range wf.ElementNames {
		idx[strings.ToLower(n)] = i
	}
	power := func(names []string, sign float64) (float64, error) {
		var sum float64
		for _, n := range names {
			i, ok := idx[strings.ToLower(n)]
			if !ok {
				return 0, fmt.Errorf("element %s not in the results", n)
			}
			sum += sign * integrate(wf.Time, wf.P[i], t0, t1)
		}
		return sum / (t1 - t0), nil
	}
	var res Result
	if res.Pdrive, err = power(r.drives, -1); err != nil {
		return res, nil, err
	}
	if res.Pmagnet, err = power(r.magnets, -1); err != nil {
		return res, nil, err
	}
	if res.Pload, err = power(r.loads, 1); err != nil {
		return res, nil, err
	}
	res.Pdiss = integrate(wf.Time, wf.Ledger.Dissipation, t0, t1) / (t1 - t0)
	res.DeltaStored = at(wf.Time, wf.Ledger.StoredEnergy, t1) - at(wf.Time, wf.Ledger.StoredEnergy, t0)
	res.Residual = (res.Pdrive+res.Pmagnet-res.Pdiss)*(t1-t0) - res.DeltaStored
	res.COP = res.Pload / res.Pdrive
	return res, wf, nil
}

// FromColumns computes the same split from a table of the "* columns:"
// quantities (an ngspice twin's output), with sources delivering −V·I.
func (s Spec) FromColumns(time []float64, col map[string][]float64) (Result, error) {
	s = s.withDefaults()
	e, r, _ := s.elements()
	byName := map[string]element{}
	for _, x := range e {
		byName[x.name] = x
	}
	v := func(n string) ([]float64, error) {
		if n == "0" {
			return make([]float64, len(time)), nil
		}
		y, ok := col["v("+n+")"]
		if !ok {
			return nil, fmt.Errorf("column v(%s) missing", n)
		}
		return y, nil
	}
	t0, t1 := s.Window()
	src := func(names []string) (float64, error) {
		var sum float64
		for _, n := range names {
			x := byName[n]
			a, err := v(x.nodes[0])
			if err != nil {
				return 0, err
			}
			b, err := v(x.nodes[1])
			if err != nil {
				return 0, err
			}
			i, ok := col["i("+strings.ToLower(n)+")"]
			if !ok {
				return 0, fmt.Errorf("column i(%s) missing", strings.ToLower(n))
			}
			p := make([]float64, len(time))
			for k := range p {
				p[k] = -(a[k] - b[k]) * i[k]
			}
			sum += integrate(time, p, t0, t1)
		}
		return sum / (t1 - t0), nil
	}
	var res Result
	var err error
	if res.Pdrive, err = src(r.drives); err != nil {
		return res, err
	}
	if res.Pmagnet, err = src(r.magnets); err != nil {
		return res, err
	}
	for _, n := range r.loads {
		x := byName[n]
		a, err := v(x.nodes[0])
		if err != nil {
			return res, err
		}
		p := make([]float64, len(time))
		for k := range p {
			p[k] = a[k] * a[k] / x.value
		}
		res.Pload += integrate(time, p, t0, t1) / (t1 - t0)
	}
	res.COP = res.Pload / res.Pdrive
	return res, nil
}

// integrate is the trapezoidal integral of y over [t0, t1], with the ends
// interpolated.
func integrate(t, y []float64, t0, t1 float64) float64 {
	var sum float64
	for k := 1; k < len(t); k++ {
		a, b := math.Max(t[k-1], t0), math.Min(t[k], t1)
		if b <= a {
			continue
		}
		lerp := func(x float64) float64 {
			if t[k] == t[k-1] {
				return y[k]
			}
			return y[k-1] + (x-t[k-1])/(t[k]-t[k-1])*(y[k]-y[k-1])
		}
		sum += 0.5 * (lerp(a) + lerp(b)) * (b - a)
	}
	return sum
}

// at interpolates y linearly at time x.
func at(t, y []float64, x float64) float64 {
	k := 1
	for k < len(t)-1 && t[k] < x {
		k++
	}
	if t[k] == t[k-1] {
		return y[k]
	}
	f := (x - t[k-1]) / (t[k] - t[k-1])
	return y[k-1] + f*(y[k]-y[k-1])
}
