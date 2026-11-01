// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command pcb_layout evolves the geometry of a pair of parallel microstrip
// traces on FR-4 (milestone M9i). For every candidate the per-unit-length
// L and C matrices are extracted from the cross-section by the
// quasi-static field solver (em.Section: div(ε grad φ) = 0, L = μ₀ε₀C₀⁻¹),
// the pair becomes a coupled-line ladder in the circuit (em.CoupledLadder,
// mutual inductance as K, mutual capacitance as a capacitor), and an AC
// analysis measures the near-end crosstalk onto the quiet trace.
//
// Genes: trace width w and gap s (substrate h = 1.6 mm, εr = 4.4 fixed).
// Objectives: characteristic impedance 50 Ω; near-end crosstalk at
// 100 MHz no worse than −30 dB; and the least board width 2w + s.
//
//	go run ./examples/pcb_layout [-gens 12] [-pop 20]
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
)

const (
	subH     = 1.6e-3 // substrate thickness [m]
	epsR     = 4.4
	length   = 0.1 // trace length [m]
	sections = 20
	fNext    = 100e6
)

// Pair is a trace geometry.
type Pair struct{ W, S float64 }

// Extract returns the 2×2 L and C matrices per unit length.
func (p Pair) Extract() (*em.SectionResult, error) {
	step := math.Min(math.Min(p.W, p.S), subH) / 4
	rows := math.Ceil(subH / step)
	step = subH / rows
	box := 2*p.W + p.S + 12*subH
	x0 := box/2 - p.S/2 - p.W
	s := em.Section{W: box, H: 8 * subH, Step: step,
		Conductors: []em.Conductor2D{
			em.Strip("a", x0, x0+p.W, subH),
			em.Strip("v", x0+p.W+p.S, x0+2*p.W+p.S, subH),
		},
		EpsR: func(x, y float64) float64 {
			if y < subH {
				return epsR
			}
			return 1
		}}
	return s.Solve()
}

// Z0 is a trace's characteristic impedance with the other trace present
// (√(L₁₁/C₁₁)).
func Z0(r *em.SectionResult) float64 { return math.Sqrt(r.L[0][0] / r.C[0][0]) }

// Circuit is the coupled pair: aggressor driven through Z₀ at "a_in",
// both traces terminated in Z₀ at both ends; near-end crosstalk is
// v(v_in)/v(a_in).
func (p Pair) Circuit(r *em.SectionResult) (*circuit.Circuit, error) {
	z := Z0(r)
	b := circuit.NewBuilder("coupled microstrip pair")
	b.Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "src", "0")
	b.Resistor("RS", "src", "a_in", z).Resistor("RA", "a_out", "0", z)
	b.Resistor("RVN", "v_in", "0", z).Resistor("RVF", "v_out", "0", z)
	if err := em.CoupledLadder(b, "X", []string{"a_in", "v_in"}, []string{"a_out", "v_out"}, "0", r.L, r.C, length, sections); err != nil {
		return nil, err
	}
	return b.Build()
}

// NEXT is the near-end crosstalk |v(v_in)/v(a_in)| in dB at f.
func NEXT(c *circuit.Circuit, f float64, opt analysis.Options) (float64, error) {
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, opt)
	if err != nil {
		return 0, err
	}
	va, _ := ac.NodeVoltage("a_in")
	vv, _ := ac.NodeVoltage("v_in")
	return 20 * math.Log10(cmplx.Abs(vv[0])/cmplx.Abs(va[0])), nil
}

// Problem is the GA problem.
func Problem() evo.Problem {
	extract := func(e *evo.Eval) (*em.SectionResult, error) {
		v, err := e.Memo("extract", func() (any, error) { return Pair{W: e.Values[0], S: e.Values[1]}.Extract() })
		if err != nil {
			return nil, err
		}
		return v.(*em.SectionResult), nil
	}
	return evo.Problem{
		Space: evo.Space{{Name: "w", Min: 0.3e-3, Max: 6e-3}, {Name: "s", Min: 0.2e-3, Max: 10e-3}},
		Build: func(v []float64) (*circuit.Circuit, error) {
			p := Pair{W: v[0], S: v[1]}
			r, err := p.Extract()
			if err != nil {
				return nil, err
			}
			return p.Circuit(r)
		},
		Objectives: []evo.Objective{
			{Name: "Z0 error", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				r, err := extract(e)
				if err != nil {
					return 0, err
				}
				return math.Abs(Z0(r)-50) / 50, nil
			}},
			{Name: "NEXT above −30 dB", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				n, err := NEXT(e.Circuit, fNext, e.AnalysisOptions())
				if err != nil {
					return 0, err
				}
				return math.Max(0, (n+30)/30), nil
			}},
			{Name: "board width", Weight: 0.01, Eval: func(e *evo.Eval) (float64, error) {
				return (2*e.Values[0] + e.Values[1]) / 10e-3, nil
			}},
		},
	}
}

func main() {
	gens := flag.Int("gens", 12, "generations")
	pop := flag.Int("pop", 20, "population")
	flag.Parse()
	s, err := evo.Run(context.Background(), Problem(), evo.Options{Population: *pop, Generations: *gens, Seed: 1}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	p := Pair{W: s.Best.Values[0], S: s.Best.Values[1]}
	r, err := p.Extract()
	if err != nil {
		log.Fatal(err)
	}
	c, _ := p.Circuit(r)
	n, _ := NEXT(c, fNext, analysis.Options{})
	k := r.L[0][1] / math.Sqrt(r.L[0][0]*r.L[1][1])
	fmt.Printf("evolved (%d evaluations): w %.3f mm, s %.3f mm, board %.2f mm\n", s.Evaluations, p.W*1e3, p.S*1e3, (2*p.W+p.S)*1e3)
	fmt.Printf("  Z₀ %.2f Ω, ε_eff %.3f, L₁₁ %.4g nH/m, C₁₁ %.4g pF/m, K = L₁₂/L₁₁ %.4f, NEXT %.1f dB at %.0f MHz\n",
		Z0(r), r.C[0][0]/r.C0[0][0], r.L[0][0]*1e9, r.C[0][0]*1e12, k, n, fNext/1e6)
}
