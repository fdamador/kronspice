// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command homopolar models Tesla's homopolar (unipolar) dynamo, U.S.
// Patent 406,968 ("Dynamo-Electric Machine," filed March 23, 1889, issued
// July 16, 1889), from the Lorentz force (CLAUDE.md §3f, milestone M9n).
//
// The physics: a conducting disk spinning at ω in an axial field B carries
// its charges at v = ω×r; the Lorentz force q(v×B), obtained by grade
// selection ⟨F·Q†⟩₁ (em.LorentzGA), pushes them radially, so the disk is a
// source of EMF ½Bω(r_o² − r_i²) between axle and rim (em.MotionalEMF).
//
// The circuit: each disk is kronspice's machine element (N) with two
// windings — the disk's radial path and the field winding — and a speed
// matrix G whose one entry K = ½·B_per_amp·(r_o² − r_i²) gives the speed
// voltage ω·K·i_f in the disk and the torque K·i_disk·i_f on the shaft.
// The shaft is a mechanical port driven at ω, so the ledger books the
// shaft power as inflow and shows it leaving as heat in the load and the
// losses.
//
// Tesla's arrangement (patent, Background and Specification): two disks in
// fields of opposite polarity, joined at their flanged rims by a flexible
// metal belt, current taken from the two shafts — the disks' EMFs add, and
// no brush has to bear on a fast-moving rim. The program compares it with a
// single disk with a rim brush, then evolves the two-disk machine for
// efficiency at a required output power.
//
// Contact and belt resistances are illustrative inputs (declared below),
// not measured values.
//
//	go run ./examples/homopolar [-gens 20] [-pop 30]
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
	"github.com/fdamador/kronspice/output"
)

// Design is one machine.
type Design struct {
	Ro, Ri    float64 // disk outer and inner (axle) radius [m]
	Thickness float64 // [m]
	B         float64 // axial field [T]
	Omega     float64 // shaft speed [rad/s]
	Load      float64 // load resistance [Ω]
	TwoDisk   bool    // Tesla's two disks and belt (else one disk, rim brush)
	// Declared contact resistances (illustrative): sliding contact on a
	// slow shaft end, on a fast rim, and the belt between two rims.
	RShaft, RRim, RBelt float64
}

// Defaults returns a reference design: 15 cm copper disks at 3000 rpm in
// 0.8 T.
func Defaults() Design {
	return Design{Ro: 0.15, Ri: 0.01, Thickness: 5e-3, B: 0.8, Omega: 2 * math.Pi * 50, Load: 0.01,
		TwoDisk: true, RShaft: 50e-6, RRim: 500e-6, RBelt: 100e-6}
}

// EMF is one disk's EMF from the Lorentz force, integrated along a radius.
func (d Design) EMF() float64 {
	var path [][3]float64
	for k := 0; k <= 100; k++ {
		r := d.Ri + (d.Ro-d.Ri)*float64(k)/100
		path = append(path, [3]float64{r, 0, 0})
	}
	return em.MotionalEMF(path,
		func(p [3]float64) [3]float64 { return [3]float64{-d.Omega * p[1], d.Omega * p[0], 0} },
		func([3]float64) [3]float64 { return [3]float64{0, 0, d.B} })
}

// DiskResistance is the radial spreading resistance ln(r_o/r_i)/(2πσt).
func (d Design) DiskResistance() float64 {
	return math.Log(d.Ro/d.Ri) * material.Copper.Resistivity / (2 * math.Pi * d.Thickness)
}

// disk returns the machine spec of one disk: winding 1 the disk's radial
// path (axle → rim), winding 2 the field winding carrying 1 A, with
// K = ½B(r_o² − r_i²) per ampere of field and the sign of the field.
func (d Design) disk(polarity float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0] = 1e-9 // the disk path's own inductance (negligible)
	s.L0[1][1] = 1e-3
	s.G[0][1] = polarity * 0.5 * d.B * (d.Ro*d.Ro - d.Ri*d.Ri)
	return s
}

// Circuit builds the machine with its load. Shaft node(s) "w1" ("w2") are
// mechanical, driven at Omega by velocity sources VW1 (VW2); the load is
// "RL".
func (d Design) Circuit() (*circuit.Circuit, error) {
	b := circuit.NewBuilder("homopolar dynamo").Domain(exchange.Mechanical, "w1", "w2")
	b.MachineModel("dA", d.disk(1)).MachineModel("dB", d.disk(-1))
	rd := d.DiskResistance()
	// Field windings, each driven at 1 A (the field strength is in K).
	b.CurrentSource("IFA", "0", "fa", 1).Resistor("RFA", "fa", "0", 1e3)
	b.VoltageSource("VW1", "w1", "0", d.Omega)
	// Terminal + is the first shaft; the disk path runs axle → rim.
	b.Resistor("RS1", "t1", "ax1", d.RShaft).Resistor("RD1", "ax1", "c1", rd)
	b.Machine("NA", "dA", [][2]string{{"rim1", "c1"}, {"fa", "0"}}, "w1", "0")
	if d.TwoDisk {
		b.CurrentSource("IFB", "0", "fb", 1).Resistor("RFB", "fb", "0", 1e3)
		b.VoltageSource("VW2", "w2", "0", d.Omega)
		// Belt from rim to rim, then the second disk (winding axle → rim
		// like the first, field reversed: its EMF points rim → axle, so
		// traversed rim → axle it adds to the first), then the second shaft
		// terminal.
		b.Resistor("RBELT", "rim1", "rim2", d.RBelt)
		b.Machine("NB", "dB", [][2]string{{"rim2", "c2"}, {"fb", "0"}}, "w2", "0")
		b.Resistor("RD2", "c2", "ax2", rd).Resistor("RS2", "ax2", "t2", d.RShaft)
	} else {
		// One disk: a brush on the fast rim to the second terminal.
		b.Resistor("RRIM", "rim1", "t2", d.RRim)
		b.VoltageSource("VW2", "w2", "0", 0).Resistor("RW2", "w2", "0", 1)
	}
	b.Resistor("RL", "t1", "t2", d.Load)
	// The armature circuit floats; a megohm references it to ground.
	b.Resistor("RREF", "t2", "0", 1e6)
	return b.Build()
}

// Result is a solved design.
type Result struct {
	EMF        float64 // Lorentz-force EMF of one disk [V]
	VLoad      float64 // load voltage [V]
	PLoad      float64 // [W]
	PShaft     float64 // mechanical input from the ledger [W]
	Efficiency float64
	Torque     [2]float64 // shaft torques [N·m]
}

// Solve finds the operating point.
func (d Design) Solve() (*Result, *output.Result, error) {
	c, err := d.Circuit()
	if err != nil {
		return nil, nil, err
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	return d.result(r)
}

func (d Design) result(r *output.Result) (*Result, *output.Result, error) {
	res := &Result{EMF: d.EMF()}
	rl, _ := r.Element("RL")
	res.VLoad, res.PLoad = rl.Voltage, rl.PowerAbsorbed
	res.PShaft = r.Ledger.Inflow - r.Ledger.Outflow
	for k, n := range []string{"VW1", "VW2"} {
		if el, ok := r.Element(n); ok {
			res.Torque[k] = -el.Current
		}
	}
	eta, err := em.OPEfficiency(r, "RL")
	if err != nil {
		return nil, nil, err
	}
	res.Efficiency = eta
	return res, r, nil
}

// Problem evolves the two-disk machine: outer radius, thickness, field and
// load, for the least loss at an output of at least pTarget watts.
func Problem(pTarget float64) evo.Problem {
	design := func(v []float64) Design {
		d := Defaults()
		d.Ro, d.Thickness, d.B, d.Load = v[0], v[1], v[2], v[3]
		return d
	}
	return evo.Problem{
		Space: evo.Space{
			{Name: "r_o", Min: 0.05, Max: 0.4},
			{Name: "thickness", Min: 1e-3, Max: 30e-3},
			{Name: "B", Min: 0.1, Max: 1.5}, // iron saturates near 1.5–2 T
			{Name: "R_load", Min: 1e-4, Max: 1},
		},
		Build: func(v []float64) (*circuit.Circuit, error) { return design(v).Circuit() },
		Objectives: []evo.Objective{
			evo.OPEfficiency("RL"),
			{Name: "power shortfall", Weight: 10, Eval: func(e *evo.Eval) (float64, error) {
				r, err := e.OP()
				if err != nil {
					return 0, err
				}
				el, _ := r.Element("RL")
				return math.Max(0, (pTarget-el.PowerAbsorbed)/pTarget), nil
			}},
		},
	}
}

func report(name string, d Design) {
	res, r, err := d.Solve()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("%s: r_o %.3f m, t %.1f mm, B %.2f T, ω %.0f rad/s, R_L %.3g Ω\n", name, d.Ro, d.Thickness*1e3, d.B, d.Omega, d.Load)
	fmt.Printf("  EMF per disk %.4f V (Lorentz) vs ½Bω(r_o²−r_i²) %.4f V; disk %.3g Ω\n",
		res.EMF, 0.5*d.B*d.Omega*(d.Ro*d.Ro-d.Ri*d.Ri), d.DiskResistance())
	fmt.Printf("  load %.4f V, %.1f W; shaft %.1f W (torques %.3g, %.3g N·m); efficiency %.4f; ledger residual %.2g\n",
		res.VLoad, res.PLoad, res.PShaft, res.Torque[0], res.Torque[1], res.Efficiency, r.Ledger.RelativeResidual())
}

func main() {
	gens := flag.Int("gens", 20, "generations")
	pop := flag.Int("pop", 30, "population")
	ptarget := flag.Float64("p", 2000, "required output power [W]")
	flag.Parse()
	d := Defaults()
	report("Tesla two-disk, belt", d)
	d1 := d
	d1.TwoDisk = false
	report("single disk, rim brush", d1)
	s, err := evo.Run(context.Background(), Problem(*ptarget), evo.Options{Population: *pop, Generations: *gens, Seed: 1}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	v := s.Best.Values
	best := Defaults()
	best.Ro, best.Thickness, best.B, best.Load = v[0], v[1], v[2], v[3]
	report(fmt.Sprintf("evolved for ≥ %.0f W (%d evaluations)", *ptarget, s.Evaluations), best)
}
