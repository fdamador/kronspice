// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command pulsed_bifilar drives a coil with a train of one-sign voltage
// pulses and finds the pulse rates at which the energy builds up into a
// large standing wave (CLAUDE.md §3e step 3), for four coils of the same
// 0.5 mm copper wire (about 125.7 m of it):
//
//   - a cylindrical solenoid wound ordinarily (400 turns, 10 cm across);
//   - the same solenoid wound as a Tesla bifilar coil (US 512,340);
//   - a conical bifilar coil (45° cone opening upward from a 4 cm base);
//   - a flat-spiral bifilar coil (Archimedean spiral from 4 cm across);
//   - the solenoid wound as an opposite-current ("non-inductive") bifilar:
//     the wire folded in half and the doubled wire wound side by side, so
//     strand 2 comes back down between the turns of strand 1 with its
//     current circulating the other way;
//   - the bifilar solenoid with its two strands in parallel between the same
//     two terminals, currents circulating the same way, and with currents
//     opposing.
//
// Each coil is a distributed network, not a lumped L and C: its turns in
// 20 sections. Every turn is a circle; the inductance matrix of the
// sections is summed turn by turn from Maxwell's formula for coaxial
// circles (em.MutualCircles, Treatise art. 701) with a single turn's self
// inductance μ0R(ln(8R/r_w) − 2) (high-frequency form, F. W. Grover,
// Inductance Calculations, 1946), and coupled in kronspice with K elements.
// Capacitance between physically adjacent turns is the two-wire formula
// (as em.Coil.TurnCapacitance, with each pair's mean circumference); in a
// bifilar coil adjacent turns belong to the two strands, half a wire apart.
// Capacitance to the surroundings treats the winding as a sheet facing the
// nearest grounded surface: the solenoid stands in a grounded cylindrical
// can of twice its radius (2πε0/ln(R_can/R) per metre of length); the cone
// and the spiral sit 5 cm above a grounded plane (ε0·area/height for each
// turn's ring of the sheet, a parallel-plate estimate). The coil floats;
// the pulses drive a short 5-turn primary just below its base (around it for
// the solenoid), and reach the coil only through the magnetic coupling
// (Maxwell's coaxial-circle mutuals again), as in a Tesla coil.
//
// For each winding the first resonance f1 is found by an AC sweep. Then
// one-sign pulses (100 V, 50 Ω source) are applied at a set of repetition
// rates. A one-sign pulse train at rate r contains only the frequencies r,
// 2r, 3r, …, so it builds up the coil only when one of them lands on f1:
// at r = f1 and at whole fractions f1/2, f1/3 (each pulse then catches the
// wave at the same phase); at 1.5·f1 or 2·f1 none does. The ledger shows
// that every joule in the coil came from the pulse source.
//
// Approximations (declared): sections are equipotential inside;
// the junction lead of the bifilar coil and radiation are left out; wire
// resistance with skin effect at f1 (as em.Coil.Resistance). Data for the
// coil pictures are written to viz/.
//
//	go run ./examples/pulsed_bifilar
package main

import (
	"encoding/csv"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime/debug"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/material"
	"github.com/fdamador/kronspice/output"
)

// Declared wire, geometry and drive.
const (
	Sections   = 20
	WireRadius = 0.25e-3 // m (0.5 mm copper)
	Pitch      = 0.7e-3  // m between neighbouring wires: 0.2 mm of enamel
	EnamelEps  = 3.0     // declared
	PulseVolts = 100.0
	SourceOhms = 50.0
	PulseEdge  = 0.02 // rise and fall time as a fraction of the period of f1
	// Primary: a short coil that the pulses drive; the coil under test
	// floats and is excited only through the magnetic coupling, as a
	// Tesla-coil secondary is.
	PrimTurns = 5
	PrimPitch = 2e-3 // m
	PrimWireR = 1e-3 // m
	BleedOhms = 1e9  // DC reference for the floating coil
	PlaneGap  = 0.05 // m from the base of the cone and spiral to their ground plane
	Radius    = 0.05 // m: the solenoid's radius
	CanRadius = 2 * Radius
)

// Geometry places the turns: R[p], Z[p] for each physical position p (0 =
// innermost / lowest), each turn's capacitance to the surroundings, and the
// primary.
type Geometry struct {
	Name, Tag   string
	R, Z        []float64
	CGround     []float64 // per turn [F]
	PrimRadius  float64
	PrimTop     float64 // axial position of the primary's top turn
	Surrounding string
}

// N is the number of turns.
func (g Geometry) N() int { return len(g.R) }

// WireLength is Σ 2πR.
func (g Geometry) WireLength() float64 {
	s := 0.0
	for _, r := range g.R {
		s += 2 * math.Pi * r
	}
	return s
}

// Solenoid is the 400-turn cylinder, 10 cm across, in its grounded can.
func Solenoid() Geometry {
	g := Geometry{Name: "cylindrical solenoid", Tag: "", PrimRadius: 0.06, PrimTop: -4e-3,
		Surrounding: "grounded can of twice its radius"}
	cCan := 2 * math.Pi * em.Eps0 / math.Log(CanRadius/Radius) * Pitch
	for p := 0; p < 400; p++ {
		g.R = append(g.R, Radius)
		g.Z = append(g.Z, float64(p)*Pitch)
		g.CGround = append(g.CGround, cCan)
	}
	return g
}

// spiralTurns returns the even number of turns of a winding that starts at
// radius r0 and grows by dr per turn, whose wire is as long as the
// solenoid's.
func spiralTurns(r0, dr float64) int {
	target, s, n := Solenoid().WireLength(), 0.0, 0
	for s < target {
		s += 2 * math.Pi * (r0 + float64(n)*dr)
		n++
	}
	return n + n%2
}

// Cone is a conical winding: a 45° cone opening upward from 2 cm radius,
// neighbouring wires Pitch apart along the cone's surface, 5 cm above a
// grounded plane; the primary sits just below its base.
func Cone() Geometry {
	const r0 = 0.02
	d := Pitch / math.Sqrt2
	g := Geometry{Name: "conical (45°)", Tag: "cone_", PrimRadius: r0 + 0.01, PrimTop: -4e-3,
		Surrounding: "grounded plane 5 cm below the base"}
	for p := 0; p < spiralTurns(r0, d); p++ {
		r, z := r0+float64(p)*d, float64(p)*d
		g.R = append(g.R, r)
		g.Z = append(g.Z, z)
		g.CGround = append(g.CGround, em.Eps0*2*math.Pi*r*d/(z+PlaneGap))
	}
	return g
}

// Flat is a flat (Archimedean) spiral from 2 cm radius, neighbouring wires
// Pitch apart, 5 cm above a grounded plane; the primary sits just below it.
func Flat() Geometry {
	const r0 = 0.02
	g := Geometry{Name: "flat spiral", Tag: "flat_", PrimRadius: r0 + 0.01, PrimTop: -4e-3,
		Surrounding: "grounded plane 5 cm below"}
	for p := 0; p < spiralTurns(r0, Pitch); p++ {
		r := r0 + float64(p)*Pitch
		g.R = append(g.R, r)
		g.Z = append(g.Z, 0)
		g.CGround = append(g.CGround, em.Eps0*2*math.Pi*r*Pitch/PlaneGap)
	}
	return g
}

// WireResistance is the resistance of a length of the wire at frequency
// f: the DC value raised by the skin effect, R_dc·d/(4δ) when the skin
// depth δ is small against the wire (Terman 1943, as em.Coil.Resistance).
func WireResistance(length, f float64) float64 {
	area := math.Pi * WireRadius * WireRadius
	rdc := material.Copper.Resistivity * length / area
	if f <= 0 {
		return rdc
	}
	delta := math.Sqrt(material.Copper.Resistivity / (math.Pi * f * em.Mu0))
	return rdc * math.Max(1, 2*WireRadius/(4*delta))
}

// Coil is the sectioned network of one winding.
type Coil struct {
	Name     string
	Geom     Geometry
	Bifilar  bool
	Counter  bool        // opposite currents: strand 2 circulates against strand 1
	Parallel bool        // the two strands in parallel between the terminals (else in series)
	From, To []int       // end nodes (n<k>) of each section, in the wire's direction
	Top      int         // the far terminal's node (the near one is n0)
	NNodes   int         // nodes n0 … n<NNodes−1>
	Short    float64     // > 0: the coil's two ends joined through this resistance [Ω] (RSHORT); 0: ends open
	Pos      []int       // Pos[t]: physical position (0 = innermost / lowest) of the t-th turn along the wire
	Sign     []float64   // Sign[t]: +1 or −1, the sense in which turn t's current circulates
	L        [][]float64 // section inductance matrix [H]
	CAdj     [][]float64 // capacitance between sections [F] (symmetric)
	CCan     []float64   // capacitance of each section to the surroundings [F]
	Lp       float64     // primary self inductance [H]
	Mp       []float64   // primary–section mutual inductance [H]
	Ctt      float64     // between the two innermost adjacent turns [F]
	CttSum   float64     // all adjacent-turn capacitances summed [F]
}

// NumSections is the number of sections of every coil.
func NumSections() int { return Sections }

// sectionOf is the section of wire turn t of n: turns are shared out evenly.
func sectionOf(t, n int) int { return t * Sections / n }

// sectionStart is the first wire turn of section s.
func sectionStart(s, n int) int { return (s*n + Sections - 1) / Sections }

// NewCoil builds the network for a geometry. In a bifilar coil the wire
// runs out along strand 1 (positions 0, 2, 4, …) and back to start strand
// 2 (1, 3, 5, …), the end of the first joined to the start of the second.
func NewCoil(g Geometry, bifilar bool) *Coil {
	n := g.N()
	c := &Coil{Name: "ordinary " + g.Name, Geom: g, Bifilar: bifilar}
	if bifilar {
		c.Name = "Tesla bifilar " + g.Name
	}
	c.Pos = make([]int, n)
	for t := range c.Pos {
		if !bifilar {
			c.Pos[t] = t
		} else if t < n/2 {
			c.Pos[t] = 2 * t
		} else {
			c.Pos[t] = 2*(t-n/2) + 1
		}
	}
	return c.network()
}

// NewCounterCoil builds the opposite-current ("non-inductive") bifilar
// winding of a geometry: the wire folded in half at its middle and the
// doubled wire wound side by side. Strand 1 (wire turns 0 … n/2 − 1) runs
// up the even positions 0, 2, 4, …; past the fold strand 2 runs back down
// the odd positions n − 1, n − 3, …, 1, its current circulating the other
// way. Neighbouring turns then sit n − 1 − 2k wire turns apart: almost the
// whole coil voltage between them at the base, none at the fold. This is
// the winding of wire-wound non-inductive resistors (Ayrton & Perry's
// bifilar), the counterpart of Tesla's same-direction series connection.
func NewCounterCoil(g Geometry) *Coil {
	n := g.N()
	c := &Coil{Name: "opposite-current bifilar " + g.Name, Geom: g, Bifilar: true, Counter: true}
	c.Pos = make([]int, n)
	c.Sign = make([]float64, n)
	for t := range c.Pos {
		if t < n/2 {
			c.Pos[t], c.Sign[t] = 2*t, 1
		} else {
			c.Pos[t], c.Sign[t] = 2*(n-1-t)+1, -1
		}
	}
	return c.network()
}

// NewParallelCoil builds the bifilar winding of a geometry with its two
// strands in parallel: both start at the base terminal n0 and run up the
// coil side by side (strand 1 on the even positions, strand 2 on the odd),
// both ending at the far terminal. With opposite = false their currents
// circulate the same way (inductance about a quarter of the series
// winding's, as two tightly coupled halves in parallel); with opposite =
// true strand 2 is wound the other way round, so the currents oppose and
// the fields cancel. Neighbouring turns sit at the same point of their
// strands, so in both the strands are at nearly the same potential.
func NewParallelCoil(g Geometry, opposite bool) *Coil {
	n := g.N()
	c := &Coil{Name: "parallel bifilar (same direction) " + g.Name, Geom: g, Bifilar: true, Parallel: true, Counter: opposite}
	if opposite {
		c.Name = "parallel bifilar (opposite currents) " + g.Name
	}
	c.Pos = make([]int, n)
	c.Sign = make([]float64, n)
	for t := range c.Pos {
		c.Sign[t] = 1
		if t < n/2 {
			c.Pos[t] = 2 * t
		} else {
			c.Pos[t] = 2*(t-n/2) + 1
			if opposite {
				c.Sign[t] = -1
			}
		}
	}
	return c.network()
}

// Kind is the winding's file tag: ordinary, bifilar, counter, parallel or
// parallel_counter.
func (c *Coil) Kind() string {
	switch {
	case c.Parallel && c.Counter:
		return "parallel_counter"
	case c.Parallel:
		return "parallel"
	case c.Counter:
		return "counter"
	case c.Bifilar:
		return "bifilar"
	}
	return "ordinary"
}

// network fills the inductance, capacitance and primary coupling of a coil
// whose turn positions (and, if set, circulation signs) are known.
func (c *Coil) network() *Coil {
	g := c.Geom
	n := g.N()
	if c.Sign == nil {
		c.Sign = make([]float64, n)
		for t := range c.Sign {
			c.Sign[t] = 1
		}
	}
	ns := NumSections()
	// Section end nodes: in series n0 → n1 → … → n<ns>; in parallel strand 1
	// is n0 → … → n<ns/2> (the far terminal) and strand 2 runs n0 →
	// n<ns/2+1> → … → n<ns−1> → n<ns/2>.
	c.From, c.To = make([]int, ns), make([]int, ns)
	if !c.Parallel {
		for s := 0; s < ns; s++ {
			c.From[s], c.To[s] = s, s+1
		}
		c.Top, c.NNodes = ns, ns+1
	} else {
		h := ns / 2
		for s := 0; s < h; s++ {
			c.From[s], c.To[s] = s, s+1
		}
		for j := 0; j < h; j++ {
			c.From[h+j], c.To[h+j] = h+j, h+j+1
			if j == 0 {
				c.From[h] = 0
			}
			if j == h-1 {
				c.To[ns-1] = h
			}
		}
		c.Top, c.NNodes = h, ns
	}
	sec := func(t int) int { return sectionOf(t, n) }
	c.L = make([][]float64, ns)
	c.CAdj = make([][]float64, ns)
	for i := range c.L {
		c.L[i] = make([]float64, ns)
		c.CAdj[i] = make([]float64, ns)
	}
	// Mutual inductance between every pair of turns (coaxial circles).
	for a := 0; a < n; a++ {
		pa := c.Pos[a]
		for b := 0; b < n; b++ {
			pb := c.Pos[b]
			var m float64
			if a == b {
				m = em.Mu0 * g.R[pa] * (math.Log(8*g.R[pa]/WireRadius) - 2)
			} else {
				m = em.MutualCircles(g.R[pa], g.R[pb], math.Abs(g.Z[pa]-g.Z[pb]))
			}
			c.L[sec(a)][sec(b)] += c.Sign[a] * c.Sign[b] * m
		}
	}
	// Capacitance between physically adjacent turns.
	ctt := func(p int) float64 {
		d := g.R[p] + g.R[p+1] // mean diameter
		return math.Pi * em.Eps0 * EnamelEps * math.Pi * d / math.Acosh(Pitch/(2*WireRadius))
	}
	c.Ctt = ctt(0)
	at := make([]int, n)
	for t, p := range c.Pos {
		at[p] = t
	}
	for p := 0; p+1 < n; p++ {
		c.CttSum += ctt(p)
		a, b := sec(at[p]), sec(at[p+1])
		if a != b {
			c.CAdj[a][b] += ctt(p)
			c.CAdj[b][a] += ctt(p)
		}
	}
	// Primary: its own inductance and its coupling to every section.
	pz := func(i int) float64 { return g.PrimTop - float64(i)*PrimPitch }
	for i := 0; i < PrimTurns; i++ {
		for j := 0; j < PrimTurns; j++ {
			if i == j {
				c.Lp += em.Mu0 * g.PrimRadius * (math.Log(8*g.PrimRadius/PrimWireR) - 2)
			} else {
				c.Lp += em.MutualCircles(g.PrimRadius, g.PrimRadius, pz(i)-pz(j))
			}
		}
	}
	c.Mp = make([]float64, ns)
	c.CCan = make([]float64, ns)
	for t := 0; t < n; t++ {
		p := c.Pos[t]
		for i := 0; i < PrimTurns; i++ {
			c.Mp[sec(t)] += c.Sign[t] * em.MutualCircles(g.PrimRadius, g.R[p], g.Z[p]-pz(i))
		}
		c.CCan[sec(t)] += g.CGround[p]
	}
	return c
}

// PrimaryM is the coil's mutual inductance to the primary seen at its
// terminals (the open-circuit EMF per unit of primary di/dt), and for a
// parallel winding also round the internal loop that runs up strand 1 and
// back down strand 2. In series the terminal value is the sum over
// sections. In parallel, with L1, L2, Mx the strands' self and mutual
// block sums (wire direction) and M1, M2 their couplings to the primary,
// the open-circuit loop current is di_loop/dt = (M1 − M2)/(L1 + L2 − 2Mx)
// per unit primary di/dt, so the terminals see
// M1 − (L1 − Mx)(M1 − M2)/(L1 + L2 − 2Mx), and the loop M1 − M2.
func (c *Coil) PrimaryM() (terminal, loop float64) {
	if !c.Parallel {
		for _, m := range c.Mp {
			terminal += m
		}
		return terminal, 0
	}
	h := NumSections() / 2
	var l1, l2, mx, m1, m2 float64
	for i, r := range c.L {
		for j, x := range r {
			switch {
			case i < h && j < h:
				l1 += x
			case i >= h && j >= h:
				l2 += x
			default:
				mx += x / 2
			}
		}
	}
	for s, m := range c.Mp {
		if s < h {
			m1 += m
		} else {
			m2 += m
		}
	}
	loop = m1 - m2
	return m1 - (l1-mx)*loop/(l1+l2-2*mx), loop
}

// SectionSign is the sense of circulation of section s's turns (+1 or −1;
// a section never spans the fold of an opposite-current winding). Build
// gives a −1 section's inductor the reversed reference direction, so its
// K couplings, written for the inductors' own reference currents, are all
// positive, as kronspice's K element requires (the dot convention).
func (c *Coil) SectionSign(s int) float64 {
	return c.Sign[sectionStart(s, c.Geom.N())]
}

// TotalL is the coil's low-frequency terminal inductance: in series the
// sum of the whole section matrix; in parallel, with L1, L2 the two
// strands' block sums and M their cross block (in the wire's direction),
// (L1·L2 − M²)/(L1 + L2 − 2M).
func (c *Coil) TotalL() float64 {
	if !c.Parallel {
		s := 0.0
		for _, r := range c.L {
			for _, x := range r {
				s += x
			}
		}
		return s
	}
	h := NumSections() / 2
	var l1, l2, m float64
	for i, r := range c.L {
		for j, x := range r {
			switch {
			case i < h && j < h:
				l1 += x
			case i >= h && j >= h:
				l2 += x
			default:
				m += x / 2
			}
		}
	}
	return (l1*l2 - m*m) / (l1 + l2 - 2*m)
}

// SectionResistance is the wire resistance of one section at f (the
// wire shared evenly).
func (c *Coil) SectionResistance(f float64) float64 {
	return WireResistance(c.Geom.WireLength(), f) / float64(NumSections())
}

// Build makes the circuit: source (wave) → 50 Ω → primary; the coil
// sections n_s → m_s (L) → n_{s+1} (R) float, coupled to the primary and
// to each other by K elements, with their capacitances split half to each
// end node of their sections and a 1 GΩ bleeder for a DC reference.
func (c *Coil) Build(src circuit.Element, rSection float64) (*circuit.Circuit, error) {
	return c.build(src, rSection, nil)
}

// build is Build with an optional hook that adds more elements (a pickup).
func (c *Coil) build(src circuit.Element, rSection float64, extra func(*circuit.Builder)) (*circuit.Circuit, error) {
	ns := NumSections()
	b := circuit.NewBuilder(c.Name + " coil, pulse driven through a primary")
	src.Name = "VP"
	b.Add(src, "s", "0")
	b.Resistor("RS", "s", "p", SourceOhms)
	b.Inductor("LP", "p", "0", c.Lp)
	b.Resistor("RB", "n0", "0", BleedOhms)
	if c.Short > 0 {
		b.Resistor("RSHORT", "n0", "n"+strconv.Itoa(c.Top), c.Short)
	}
	for s := 0; s < ns; s++ {
		b.Mutual(fmt.Sprintf("KP%d", s), "LP", fmt.Sprintf("L%d", s), c.SectionSign(s)*c.Mp[s]/math.Sqrt(c.Lp*c.L[s][s]))
	}
	node := func(k int) string { return "n" + strconv.Itoa(k) }
	for s := 0; s < ns; s++ {
		if c.SectionSign(s) > 0 {
			b.Inductor(fmt.Sprintf("L%d", s), node(c.From[s]), fmt.Sprintf("m%d", s), c.L[s][s])
		} else {
			b.Inductor(fmt.Sprintf("L%d", s), fmt.Sprintf("m%d", s), node(c.From[s]), c.L[s][s])
		}
		b.Resistor(fmt.Sprintf("R%d", s), fmt.Sprintf("m%d", s), node(c.To[s]), rSection)
	}
	k := 0
	for i := 0; i < ns; i++ {
		for j := i + 1; j < ns; j++ {
			kk := c.SectionSign(i) * c.SectionSign(j) * c.L[i][j] / math.Sqrt(c.L[i][i]*c.L[j][j])
			if math.Abs(kk) < 1e-6 {
				continue
			}
			b.Mutual(fmt.Sprintf("K%d", k), fmt.Sprintf("L%d", i), fmt.Sprintf("L%d", j), kk)
			k++
		}
	}
	// Capacitances: to the can, and between sections; half at each end.
	cn := make([]float64, c.NNodes) // to ground per node
	for s := 0; s < ns; s++ {
		cn[c.From[s]] += c.CCan[s] / 2
		cn[c.To[s]] += c.CCan[s] / 2
	}
	for n, v := range cn {
		b.Capacitor(fmt.Sprintf("CG%d", n), node(n), "0", v)
	}
	q := 0
	for i := 0; i < ns; i++ {
		for j := i + 1; j < ns; j++ {
			if c.CAdj[i][j] == 0 {
				continue
			}
			for _, e := range [][2]int{{c.From[i], c.From[j]}, {c.To[i], c.To[j]}} {
				if e[0] == e[1] {
					continue
				}
				b.Capacitor(fmt.Sprintf("CA%d", q), node(e[0]), node(e[1]), c.CAdj[i][j]/2)
				q++
			}
		}
	}
	if extra != nil {
		extra(b)
	}
	return b.BuildTransient()
}

// SectionLimit is the highest frequency at which one section's wire is
// still shorter than a tenth of a free-space wavelength: the edge of the
// sectioned model's validity (CLAUDE.md §3e, the λ/10 rule).
func (c *Coil) SectionLimit() float64 {
	return em.SpeedOfLight / (10 * c.Geom.WireLength() / float64(NumSections()))
}

// Resonance sweeps the coil's terminal voltage |V(top) − V(base)| per volt
// of drive and returns the first resonance and its peak gain. A winding
// with no peak above gain 2 (the opposite-current bifilar) gets its
// strongest peak below SectionLimit instead; ResonanceWeak reports that.
func (c *Coil) Resonance(fLo, fHi float64) (f1, gain, q float64, err error) {
	f1, gain, q, _, err = c.ResonanceWeak(fLo, fHi)
	return
}

// ResonanceWeak is Resonance that also reports whether the peak is weak
// (no gain above 2 anywhere in the model's range).
func (c *Coil) ResonanceWeak(fLo, fHi float64) (f1, gain, q float64, weak bool, err error) {
	rs := c.SectionResistance(1e6)
	ckt, err := c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, rs)
	if err != nil {
		return
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "dec", Points: 400, FStart: fLo, FStop: fHi}, analysis.Options{})
	if err != nil {
		return
	}
	top := c.terminal(res)
	// First local maximum above the low-frequency gain, below the sections' λ/10 limit.
	for i := 1; i+1 < len(res.Freq) && res.Freq[i] < c.SectionLimit(); i++ {
		g := cabs(top[i])
		if g > cabs(top[i-1]) && g >= cabs(top[i+1]) && g > 2 {
			f1, gain = res.Freq[i], g
			break
		}
	}
	if f1 == 0 {
		weak = true
		lim := c.SectionLimit()
		for i := 1; i+1 < len(res.Freq) && res.Freq[i] < lim; i++ {
			g := cabs(top[i])
			if g > cabs(top[i-1]) && g >= cabs(top[i+1]) && g > gain {
				f1, gain = res.Freq[i], g
			}
		}
	}
	if f1 == 0 {
		return 0, 0, 0, false, fmt.Errorf("no resonance between %g and %g Hz", fLo, fHi)
	}
	// Refine with a linear sweep ±5 % and the resistance at f1.
	rs = c.SectionResistance(f1)
	ckt, err = c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, rs)
	if err != nil {
		return
	}
	res, err = analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 2001, FStart: 0.95 * f1, FStop: 1.05 * f1}, analysis.Options{})
	if err != nil {
		return
	}
	top = c.terminal(res)
	gain = 0
	for i := range res.Freq {
		if g := cabs(top[i]); g > gain {
			f1, gain = res.Freq[i], g
		}
	}
	// Half-power width.
	lo, hi := f1, f1
	for i := range res.Freq {
		if cabs(top[i]) >= gain/math.Sqrt2 {
			lo = math.Min(lo, res.Freq[i])
			hi = math.Max(hi, res.Freq[i])
		}
	}
	q = f1 / (hi - lo)
	return f1, gain, q, weak, nil
}

func cabs(z complex128) float64 { return math.Hypot(real(z), imag(z)) }

// terminal is V(top) − V(base) at every frequency.
func (c *Coil) terminal(res *output.ACResult) []complex128 {
	top, _ := res.NodeVoltage("n" + strconv.Itoa(c.Top))
	base, _ := res.NodeVoltage("n0")
	out := make([]complex128, len(top))
	for i := range top {
		out[i] = top[i] - base[i]
	}
	return out
}

// PulseRun is the result of one pulse train.
type PulseRun struct {
	Rate, Ratio        float64 // pulse rate [Hz], rate / f1
	PeakTop            float64 // largest coil terminal voltage |V(top) − V(base)| in the last 10 % [V]
	PeakGround         float64 // largest voltage of any coil node to ground [V]
	PeakTurn           float64 // largest voltage between physically adjacent turns [V]
	PeakWire           float64 // largest wire current in any coil section, last 10 % [A]
	PeakPrim           float64 // largest primary current, last 10 % [A]
	PeakShort          float64 // largest current through the join of the two ends (Short > 0), last 10 % [A]
	PIn, PHeat         float64 // mean power from the source and into heat over the run [W]
	EIn, EDiss, EStore float64 // ledger over the run [J]
	Mismatch           float64
	Time               []float64
	Nodes              [][]float64 // node voltages n0..n_ns (only when keep)
	Wire               [][]float64 // wire current of each section, along the wire (only when keep)
	Prim               []float64   // primary current (only when keep)
}

// Pulse drives the coil with one-sign pulses of width `width` at `rate`
// for `cycles` periods of f1, with the wire resistance at f1.
func (c *Coil) Pulse(f1, rate, width float64, cycles int, keep bool) (PulseRun, error) {
	per := 1 / rate
	src := circuit.Element{Kind: circuit.VoltageSource, Wave: circuit.Waveform{Kind: circuit.WavePulse,
		Args: []float64{0, PulseVolts, 0, PulseEdge / f1, PulseEdge / f1, width, per}}}
	ckt, err := c.Build(src, c.SectionResistance(f1))
	if err != nil {
		return PulseRun{}, err
	}
	stop := float64(cycles) / f1
	step := 1 / f1 / 50
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, TranOpts)
	if err != nil {
		return PulseRun{}, err
	}
	ns := NumSections()
	r := PulseRun{Rate: rate, Ratio: rate / f1}
	nodes := make([][]float64, c.NNodes)
	for k := range nodes {
		nodes[k], _ = w.NodeVoltage("n" + strconv.Itoa(k))
	}
	n := c.Geom.N()
	at := make([]int, n)
	for t, p := range c.Pos {
		at[p] = t
	}
	wire := make([][]float64, ns)
	for s := range wire {
		wire[s], _ = w.ElementCurrent(fmt.Sprintf("R%d", s)) // the wire current, the same in both inductor orientations
	}
	prim, _ := w.ElementCurrent("LP")
	short, _ := w.ElementCurrent("RSHORT")
	for i, t := range w.Time {
		if t < 0.9*stop {
			continue
		}
		for s := range wire {
			if wire[s] != nil {
				r.PeakWire = math.Max(r.PeakWire, math.Abs(wire[s][i]))
			}
		}
		if prim != nil {
			r.PeakPrim = math.Max(r.PeakPrim, math.Abs(prim[i]))
		}
		if short != nil {
			r.PeakShort = math.Max(r.PeakShort, math.Abs(short[i]))
		}
		r.PeakTop = math.Max(r.PeakTop, math.Abs(nodes[c.Top][i]-nodes[0][i]))
		for k := range nodes {
			r.PeakGround = math.Max(r.PeakGround, math.Abs(nodes[k][i]))
		}
		for p := 0; p+1 < n; p++ {
			d := math.Abs(c.TurnVoltage(nodes, i, at[p]) - c.TurnVoltage(nodes, i, at[p+1]))
			r.PeakTurn = math.Max(r.PeakTurn, d)
		}
	}
	r.EIn, r.EDiss, r.EStore, r.Mismatch = w.Ledger.EnergyIn, w.Ledger.EnergyDissipated, w.Ledger.StoredEnd-w.Ledger.StoredStart, w.Ledger.IntegrationMismatch
	r.PIn, r.PHeat = r.EIn/stop, r.EDiss/stop
	if keep {
		r.Time, r.Nodes, r.Wire, r.Prim = w.Time, nodes, wire, prim
	}
	return r, nil
}

// TranOpts loosens the absolute tolerances to the scale of this circuit
// (volts to kilovolts, milliamps to amps): SPICE's defaults (1 pA, 1 µV,
// 0.01 pC) make the step control chase the sectioned network's
// high-frequency internal ringing at the picoamp level. Relative tolerance
// stays at the SPICE 1e-3.
var TranOpts = analysis.TranOptions{Options: analysis.Options{AbsTol: 1e-6, VnTol: 1e-3}, ChgTol: 1e-12}

// TurnVoltage interpolates the voltage of wire turn t at time index i
// between the end nodes of its section.
func (c *Coil) TurnVoltage(nodes [][]float64, i, t int) float64 {
	from, to, f := c.turnFrac(t)
	return nodes[from][i]*(1-f) + nodes[to][i]*f
}

// turnFrac gives turn t's section end nodes and its fraction along the
// section.
func (c *Coil) turnFrac(t int) (from, to int, f float64) {
	n := c.Geom.N()
	s := sectionOf(t, n)
	a, b := sectionStart(s, n), sectionStart(s+1, n)
	return c.From[s], c.To[s], (float64(t-a) + 0.5) / float64(b-a)
}

// WriteViz writes, from the full-resolution run at f1, the data for the
// pictures: <coil>_envelope.csv (largest coil-end voltage in each period of
// f1), <coil>_glow.csv (largest |voltage| of every turn, by physical
// position, in each of `frames` equal time windows) and <coil>_peak.csv
// (every turn's voltage at the instant of the largest coil-end voltage in
// the last 10 % of the run).
func WriteViz(dir string, c *Coil, r PulseRun, f1 float64, frames int) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	tag := c.Geom.Tag + c.Kind()
	Turns := c.Geom.N()
	at := make([]int, Turns)
	for t, p := range c.Pos {
		at[p] = t
	}
	n := len(r.Time)
	write := func(name string, head []string, rows [][]float64) error {
		f, err := os.Create(filepath.Join(dir, tag+"_"+name+".csv"))
		if err != nil {
			return err
		}
		defer f.Close()
		cw := csv.NewWriter(f)
		cw.Write(head)
		for _, row := range rows {
			rec := make([]string, len(row))
			for i, x := range row {
				rec[i] = strconv.FormatFloat(x, 'g', 7, 64)
			}
			cw.Write(rec)
		}
		cw.Flush()
		return cw.Error()
	}
	posHead := func(first string) []string {
		h := []string{first}
		for p := 0; p < Turns; p++ {
			h = append(h, fmt.Sprintf("pos%d", p))
		}
		return h
	}
	// Envelope per period of f1.
	var env [][]float64
	per := 1 / f1
	cur, start, k := 0.0, 0.0, 0
	for i, t := range r.Time {
		for t >= start+per {
			if k > 0 {
				env = append(env, []float64{start + per/2, cur})
			}
			start += per
			cur, k = 0, 0
		}
		cur = math.Max(cur, math.Abs(r.Nodes[c.Top][i]-r.Nodes[0][i]))
		k++
	}
	if err := write("envelope", []string{"time", "coil_end_volts"}, env); err != nil {
		return err
	}
	// Glow frames.
	var glow [][]float64
	for f := 0; f < frames; f++ {
		a, b := f*n/frames, (f+1)*n/frames
		row := make([]float64, Turns+1)
		row[0] = r.Time[b-1]
		for i := a; i < b; i++ {
			for p := 0; p < Turns; p++ {
				row[p+1] = math.Max(row[p+1], math.Abs(c.TurnVoltage(r.Nodes, i, at[p])))
			}
		}
		glow = append(glow, row)
	}
	if err := write("glow", posHead("time"), glow); err != nil {
		return err
	}
	// Snapshot at the largest coil-end voltage near the end.
	best, bi := -1.0, 0
	for i := 9 * n / 10; i < n; i++ {
		if v := math.Abs(r.Nodes[c.Top][i] - r.Nodes[0][i]); v > best {
			best, bi = v, i
		}
	}
	row := []float64{r.Time[bi]}
	for p := 0; p < Turns; p++ {
		row = append(row, c.TurnVoltage(r.Nodes, bi, at[p]))
	}
	if err := write("peak", posHead("time"), [][]float64{row}); err != nil {
		return err
	}
	// Geometry: radius and height of every position.
	var geom [][]float64
	for p := 0; p < Turns; p++ {
		geom = append(geom, []float64{float64(p), c.Geom.R[p], c.Geom.Z[p], float64(at[p]), c.Sign[at[p]]})
	}
	if err := write("geom", []string{"pos", "radius", "height", "turn", "sign"}, geom); err != nil {
		return err
	}
	return writeCurrents(write, c, r, f1)
}

// writeCurrents writes, from the run at f1, <coil>_ienv.csv (largest wire
// current in any section and largest primary current in each period of
// f1) and <coil>_current.csv (each section's wire current along the wire
// and the circulating current Sign·I of its turns, at the instant of the
// largest wire current in the last 10 % of the run, with each section's
// peak |I| over that window).
func writeCurrents(write func(string, []string, [][]float64) error, c *Coil, r PulseRun, f1 float64) error {
	ns := NumSections()
	n := len(r.Time)
	if len(r.Wire) != ns || r.Prim == nil {
		return fmt.Errorf("%s: currents not kept", c.Name)
	}
	var env [][]float64
	per := 1 / f1
	wmax, pmax, start, k := 0.0, 0.0, 0.0, 0
	for i, t := range r.Time {
		for t >= start+per {
			if k > 0 {
				env = append(env, []float64{start + per/2, wmax, pmax})
			}
			start += per
			wmax, pmax, k = 0, 0, 0
		}
		for s := 0; s < ns; s++ {
			wmax = math.Max(wmax, math.Abs(r.Wire[s][i]))
		}
		pmax = math.Max(pmax, math.Abs(r.Prim[i]))
		k++
	}
	if err := write("ienv", []string{"time", "wire_amps", "primary_amps"}, env); err != nil {
		return err
	}
	best, bi := -1.0, 0
	peak := make([]float64, ns)
	for i := 9 * n / 10; i < n; i++ {
		for s := 0; s < ns; s++ {
			a := math.Abs(r.Wire[s][i])
			peak[s] = math.Max(peak[s], a)
			if a > best {
				best, bi = a, i
			}
		}
	}
	var rows [][]float64
	for s := 0; s < ns; s++ {
		i := r.Wire[s][bi]
		rows = append(rows, []float64{float64(s), i, c.SectionSign(s) * i, peak[s]})
	}
	return write("current", []string{"section", "wire_amps", "circulating_amps", "peak_amps"}, rows)
}

// Cases are the seven coils, in report order.
func Cases() []*Coil {
	return []*Coil{NewCoil(Solenoid(), false), NewCoil(Solenoid(), true), NewCoil(Cone(), true), NewCoil(Flat(), true),
		NewCounterCoil(Solenoid()), NewParallelCoil(Solenoid(), false), NewParallelCoil(Solenoid(), true)}
}

func main() {
	// Each long transient holds every element's waveform until it returns;
	// keep the heap from growing across the 36 runs.
	debug.SetMemoryLimit(1 << 30)
	only := flag.Bool("only-pickup", false, "run only the pickup study")
	onlyVI := flag.Bool("only-vi", false, "run only the V–I sweep of the three solenoids (viz/vi_sweep.csv)")
	flag.Parse()
	if *onlyVI {
		if err := WriteVISweep("viz"); err != nil {
			log.Fatal(err)
		}
		return
	}
	if *only {
		if err := RunPickups(); err != nil {
			log.Fatal(err)
		}
		return
	}
	fmt.Println("Seven coils of 0.5 mm copper, each about 125.7 m of wire, driven through a 5-turn primary by one-sign 100 V pulses")
	for _, c := range Cases() {
		g := c.Geom
		f1, gain, q, weak, err := c.ResonanceWeak(1e3, 1e8)
		if err != nil {
			log.Fatal(c.Name, ": ", err)
		}
		lo, hi, zmax := g.R[0], g.R[g.N()-1], g.Z[g.N()-1]
		fmt.Printf("\n== %s winding ==\n", c.Name)
		fmt.Printf("   %d turns, radius %.1f to %.1f cm, height %.1f cm, wire %.1f m; %s\n", g.N(), lo*100, hi*100, zmax*100, g.WireLength(), g.Surrounding)
		if l := c.TotalL(); l >= 1e-3 {
			fmt.Printf("   inductance %.4g mH", l*1e3)
		} else {
			fmt.Printf("   inductance %.4g µH", l*1e6)
		}
		if g.Tag == "" && !c.Parallel && !c.Counter {
			fmt.Printf(" (Nagaoka sheet formula %.3f mH)", em.Nagaoka(Radius, float64(g.N())*Pitch, float64(g.N()))*1e3)
		}
		fmt.Printf("; capacitance between all neighbouring turns %.3g nF\n", c.CttSum*1e9)
		if weak {
			fmt.Printf("   no resonance with gain above 2 from 1 kHz to the sections' λ/10 limit (%.3g Hz): the strands' fields cancel\n", c.SectionLimit())
			fmt.Printf("   strongest peak f1 = %.4g Hz (wire = %.3f wavelengths in free space), AC gain at the terminals ×%.3f, Q ≈ %.0f\n",
				f1, g.WireLength()*f1/em.SpeedOfLight, gain, q)
		} else {
			fmt.Printf("   first resonance f1 = %.4g Hz (wire = %.3f wavelengths in free space), AC gain at top ×%.1f, Q ≈ %.0f\n",
				f1, g.WireLength()*f1/em.SpeedOfLight, gain, q)
		}
		mt, ml := c.PrimaryM()
		fmt.Printf("   mutual inductance to the primary %.4g µH", mt*1e6)
		if c.Parallel {
			fmt.Printf(" at the terminals; %.4g µH round the loop up one strand and down the other", ml*1e6)
		}
		fmt.Println()
		width := 0.1 / f1
		cycles := int(math.Min(600, math.Max(60, 5*q/math.Pi+20))) // ≈ 5 ring-up time constants Q/(πf1)
		fmt.Printf("   one-sign pulses: %.0f V, %.0f Ω, width %.3g s, %d cycles of f1 each run\n", PulseVolts, SourceOhms, width, cycles)
		fmt.Printf("   %10s %12s %14s %18s %12s %12s %12s %12s %12s %12s %12s\n", "rate/f1", "rate [Hz]", "coil ends [V]", "adjacent turns [V]",
			"wire [A]", "primary [A]", "in [W]", "heat [W]", "in [J]", "heat [J]", "stored [J]")
		var best PulseRun
		for _, ratio := range []float64{1.0 / 3, 0.5, 0.9, 0.97, 1, 1.03, 1.1, 1.5, 2} {
			keep := ratio == 1
			r, err := c.Pulse(f1, ratio*f1, width, cycles, keep)
			if err != nil {
				log.Fatal(err)
			}
			fmt.Printf("   %10.3f %12.4g %14.1f %18.1f %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g %12.4g\n", ratio, r.Rate, r.PeakTop, r.PeakTurn,
				r.PeakWire, r.PeakPrim, r.PIn, r.PHeat, r.EIn, r.EDiss, r.EStore)
			if keep {
				best = r
			}
			debug.FreeOSMemory()
		}
		fmt.Printf("   ledger at f1: in %.5g J = heat %.5g J + stored %.5g J (+ time-stepping error %.2g J)\n", best.EIn, best.EDiss, best.EStore, best.Mismatch)
		if err := WriteViz("viz", c, best, f1, 48); err != nil {
			log.Fatal(err)
		}
	}
	fmt.Println("\nPulses build each coil up only in step with its resonance (rate f1, or f1/2, f1/3: one of the train's")
	fmt.Println("frequencies r, 2r, 3r … lands on f1); at 0.9, 1.1, 1.5 or 2 × f1 none does. Every joule came from the pulse source.")
	if err := WriteVISweep("viz"); err != nil {
		log.Fatal(err)
	}
	if err := RunPickups(); err != nil {
		log.Fatal(err)
	}
}
