// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// Pipes is a rectifier built in one tank from concentric pipes standing
// in water: a rod at the centre, pipes around it, the tank wall outside.
// Each annulus of water between two neighbouring pipes is one cell. Two
// ways to wire them:
//
//   - Parallel: the pipes alternate aluminium and iron; every aluminium
//     pipe goes to one bus bar and every iron pipe to the other. Every
//     annulus is aluminium facing iron, so every cell is a valve in the
//     same direction between the two buses, and an inner pipe works on
//     both faces. Neighbouring annuli sit at the same voltage, so water
//     around the pipe ends carries no current: the ends need no seals.
//     The tank holds one forming voltage.
//   - Series: no bus bars; current crosses every annulus in turn. For
//     all cells to face the same way each pipe must be aluminium on its
//     outer face and iron (or stainless) on its inner face — a clad pipe
//     (plain alternating aluminium and iron would put neighbouring cells
//     back to back and block both ways). The tank holds N forming
//     voltages, but neighbouring annuli now differ by a cell's voltage,
//     so water around a pipe's end shorts a cell: each end must be
//     sealed, and Rseal is the leak resistance left between neighbouring
//     annuli (0: a perfect seal).
//
// Annulus j (j = 0 … N−1) runs from radius a_j = R0 + j·(Gap + Wall) to
// b_j = a_j + Gap over the wetted length Len: its water resistance is
// ρ·ln(b/a)/(2π·Len) and its aluminium face 2π·r·Len, at r = a_j for the
// series stack (aluminium outside each inner pipe) and on whichever side
// is aluminium for the parallel one (the rod is aluminium).
type Pipes struct {
	Series bool
	N      int     // annuli (cells)
	R0     float64 // rod radius [m]
	Gap    float64 // water annulus [m]
	Wall   float64 // pipe wall [m]
	Len    float64 // wetted length [m]
	Vform  float64 // forming voltage of every aluminium face [V]
	Kappa  float64 // varactor coefficient of the oxide film
	Rseal  float64 // series only: leak between neighbouring annuli around the sealed pipe ends [Ω], 0 = perfect
	// Graded (series only): the outer pipes' aluminium is masked
	// (epoxy-coated) above a height that gives every annulus the same
	// aluminium area. Ungraded, the inner annulus has the least aluminium,
	// so the least film capacitance, and takes most of the reverse voltage
	// (the AC voltage divides as 1/C): the outer annuli hardly help.
	Graded bool
}

const (
	wallDef = 1.5e-3 // 1.5 mm pipe walls (declared)
	rFloat  = 1e10   // Ω, DC shunt per series film
)

func (p Pipes) radii(j int) (a, b float64) {
	a = p.R0 + float64(j)*(p.Gap+p.Wall)
	return a, a + p.Gap
}

// alFace is the radius of annulus j's aluminium face.
func (p Pipes) alFace(j int) float64 {
	a, b := p.radii(j)
	if p.Series || j%2 == 0 {
		return a // the inner pipe of the annulus is aluminium (outside)
	}
	return b
}

// Cell is annulus j as a valve cell (one cell, one "plate" of its
// aluminium face's area), its water resistance and its area.
func (p Pipes) Cell(j int) (s Stack, rw, area float64) {
	a, b := p.radii(j)
	L := p.Len
	if p.Series && p.Graded {
		L *= p.alFace(0) / p.alFace(j)
	}
	rw = rhoWater * math.Log(b/a) / (2 * math.Pi * L)
	area = 2 * math.Pi * p.alFace(j) * L
	return Stack{Ns: 1, Np: 1, Area: area, Gap: rw * area / rhoWater, Vform: p.Vform, Kappa: p.Kappa}, rw, area
}

// Hold is the largest reverse voltage the tank holds [V].
func (p Pipes) Hold() float64 {
	if p.Series {
		return float64(p.N) * p.Vform
	}
	return p.Vform
}

// TankRadius is the tank wall's inner radius [m].
func (p Pipes) TankRadius() float64 { _, b := p.radii(p.N - 1); return b }

// Film is the oxide capacitance between the terminals [F].
func (p Pipes) Film() float64 {
	var inv, sum float64
	for j := 0; j < p.N; j++ {
		s, _, _ := p.Cell(j)
		inv += 1 / s.C0()
		sum += s.C0()
	}
	if p.Series {
		return 1 / inv
	}
	return sum
}

// Rw is the water resistance between the terminals [Ω].
func (p Pipes) Rw() float64 {
	var sum, inv float64
	for j := 0; j < p.N; j++ {
		_, rw, _ := p.Cell(j)
		sum += rw
		inv += 1 / rw
	}
	if p.Series {
		return sum
	}
	return 1 / inv
}

// Area is the total aluminium face [m²].
func (p Pipes) Area() float64 {
	var s float64
	for j := 0; j < p.N; j++ {
		_, _, a := p.Cell(j)
		s += a
	}
	return s
}

func (p Pipes) kind() string {
	if p.N == 1 {
		return "single cell"
	}
	switch {
	case p.Series && p.Graded:
		return "graded series"
	case p.Series:
		return "series"
	}
	return "parallel"
}

// Short names the geometry.
func (p Pipes) Short() string {
	return fmt.Sprintf("%s, %d annuli, rod %.3g mm, gap %.1f mm, %.3g mm long, tank Ø%.0f mm",
		p.kind(), p.N, p.R0*1e3, p.Gap*1e3, p.Len*1e3, 2e3*p.TankRadius())
}

// Describe prints the tank.
func (p Pipes) Describe() string {
	return fmt.Sprintf("%s, formed at %.0f V: holds %.0f V, aluminium %.3g mm², water %.3g Ω, film %.3g µF",
		p.Short(), p.Vform, p.Hold(), p.Area()*1e6, p.Rw(), p.Film()*1e6)
}

// Add puts the tank in like a diode from anode to cathode (conducting
// when V(anode) > V(cathode), aluminium on the cathode side).
func (p Pipes) Add(b *circuit.Builder, name, anode, cathode string) {
	if !p.Series {
		for j := 0; j < p.N; j++ {
			s, _, _ := p.Cell(j)
			s.Add(b, fmt.Sprintf("%sg%d", name, j), anode, cathode)
		}
		return
	}
	// metal k_0 = cathode (the aluminium rod) … k_N = anode (the tank
	// wall); annulus j: k_j —[film valve with half the water]— m_j —[the
	// other half]— k_{j+1}; m_j is mid-annulus water. The leak around
	// pipe j's ends joins m_{j−1} and m_j.
	metal := func(j int) string {
		switch j {
		case 0:
			return cathode
		case p.N:
			return anode
		}
		return fmt.Sprintf("%sk%d", name, j)
	}
	for j := 0; j < p.N; j++ {
		s, rw, area := p.Cell(j)
		s.Gap = 0.5 * rw * area / rhoWater
		mid := fmt.Sprintf("%sm%d", name, j)
		s.Add(b, fmt.Sprintf("%sV%d", name, j), mid, metal(j))
		b.Resistor(fmt.Sprintf("%sR%d", name, j), metal(j+1), mid, 0.5*rw)
		// The film's leakage at zero volts is ~1e-12 A: the pipes between
		// cells would float at DC. A 10 GΩ shunt per film (far below the
		// leakage at any working voltage) fixes their DC level.
		b.Resistor(fmt.Sprintf("%sG%d", name, j), mid, metal(j), rFloat)
		if p.Rseal > 0 && j > 0 {
			b.Resistor(fmt.Sprintf("%sS%d", name, j), fmt.Sprintf("%sm%d", name, j-1), mid, p.Rseal)
		}
	}
}

// meter reads a finished run over the measuring window.
type meter struct {
	node func(string) []float64
	cur  func(string) []float64
	avg  func(func(i int) float64) float64
	max  func(func(i int) float64) float64
}

// twoTerm splits an element's dissipation into conducting and reverse.
func twoTerm(m meter, elem, a, c string) (fwd, leak float64) {
	id := m.cur(elem)
	if id == nil {
		return 0, 0
	}
	va, vc := m.node(a), m.node(c)
	fwd = m.avg(func(i int) float64 { return math.Max(0, (va[i]-vc[i])*id[i]) * b2f(va[i] > vc[i]) })
	leak = m.avg(func(i int) float64 { return math.Max(0, (va[i]-vc[i])*id[i]) * b2f(va[i] <= vc[i]) })
	return
}

// measure: conducting loss (electrolysis and water), leakage, seal loss,
// and the largest reverse voltage on any one oxide over its forming
// voltage.
func (s Stack) measure(m meter, name, anode, cathode string) (fwd, leak, seal, rev float64) {
	fwd, leak = twoTerm(m, name, anode, cathode)
	va, vc := m.node(anode), m.node(cathode)
	rev = m.max(func(i int) float64 { return (vc[i] - va[i]) / s.Hold() })
	return
}

func (p Pipes) measure(m meter, name, anode, cathode string) (fwd, leak, seal, rev float64) {
	if !p.Series {
		for j := 0; j < p.N; j++ {
			f, l := twoTerm(m, fmt.Sprintf("%sg%d", name, j), anode, cathode)
			fwd += f
			leak += l
		}
		va, vc := m.node(anode), m.node(cathode)
		rev = m.max(func(i int) float64 { return (vc[i] - va[i]) / p.Vform })
		return
	}
	metal := func(j int) string {
		switch j {
		case 0:
			return cathode
		case p.N:
			return anode
		}
		return fmt.Sprintf("%sk%d", name, j)
	}
	ohmic := func(elem, a, c string) float64 {
		i := m.cur(elem)
		if i == nil {
			return 0
		}
		va, vc := m.node(a), m.node(c)
		return m.avg(func(k int) float64 { return (va[k] - vc[k]) * i[k] })
	}
	for j := 0; j < p.N; j++ {
		mid := fmt.Sprintf("%sm%d", name, j)
		f, l := twoTerm(m, fmt.Sprintf("%sV%d", name, j), mid, metal(j))
		fwd += f + ohmic(fmt.Sprintf("%sR%d", name, j), metal(j+1), mid)
		leak += l + ohmic(fmt.Sprintf("%sG%d", name, j), mid, metal(j))
		vm, vk := m.node(mid), m.node(metal(j))
		rev = math.Max(rev, m.max(func(i int) float64 { return (vk[i] - vm[i]) / p.Vform }))
		if p.Rseal > 0 && j > 0 {
			seal += ohmic(fmt.Sprintf("%sS%d", name, j), fmt.Sprintf("%sm%d", name, j-1), mid)
		}
	}
	return
}

// cellRev is the largest reverse voltage on each annulus's oxide over its
// forming voltage (series tanks: the annuli share the voltage unequally).
func (p Pipes) cellRev(m meter, name, anode, cathode string) []float64 {
	out := make([]float64, p.N)
	for j := range out {
		mid := fmt.Sprintf("%sm%d", name, j)
		k := cathode
		if j > 0 {
			k = fmt.Sprintf("%sk%d", name, j)
		}
		if !p.Series {
			mid, k = anode, cathode
		}
		vm, vk := m.node(mid), m.node(k)
		out[j] = m.max(func(i int) float64 { return (vk[i] - vm[i]) / p.Vform })
	}
	return out
}
