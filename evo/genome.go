// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package evo is the genetic-algorithm layer (CLAUDE.md §4). It evolves
// circuits through circuit.Circuit values and the analysis entry points
// only, never through simulator internals.
//
// This milestone (M5) provides the value genome: the topology is fixed and
// each gene is one component value or model parameter, on a log scale,
// optionally snapped to an E-series. Genes may be negative (negative
// resistance) or free to take either sign; nonlinear circuits (diodes) are
// evaluated like any other.
package evo

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// ErrBadGene is returned for an invalid gene definition.
var ErrBadGene = errors.New("evo: invalid gene")

// SignMode says which signs a gene may take.
type SignMode uint8

const (
	// Positive genes take values in [Min, Max] (the default).
	Positive SignMode = iota
	// Negative genes take values in [−Max, −Min] (for example a negative
	// resistance).
	Negative
	// Either lets evolution choose the sign as well as the magnitude: the
	// sign is a separate gene, the magnitude stays on a log scale in
	// [Min, Max]. Zero is never produced, so a resistor cannot become a
	// short.
	Either
)

func (s SignMode) String() string {
	switch s {
	case Negative:
		return "negative"
	case Either:
		return "either"
	default:
		return "positive"
	}
}

// ParseSignMode reads "positive", "negative" or "either" ("" = positive).
func ParseSignMode(s string) (SignMode, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "positive", "+":
		return Positive, nil
	case "negative", "-":
		return Negative, nil
	case "either", "any", "±", "+-":
		return Either, nil
	}
	return 0, fmt.Errorf("%w: sign %q (want positive, negative or either)", ErrBadGene, s)
}

// Gene is one evolvable value of a fixed-topology circuit.
//
// Exactly one target is set: Element (the element's value: ohms, farads,
// henries, or a source's DC value), Element + Param (an instance parameter,
// such as a MOSFET's channel width W), Model + Param (a .model parameter, such
// as a diode's IS or RS), or Name alone: a free design parameter (a coil's
// turn count, a plate gap) passed to Problem.Build, which turns physical
// parameters into a circuit.
type Gene struct {
	Name    string
	Element string
	Model   string
	Param   string
	// Min and Max bound the magnitude (0 < Min ≤ Max). Evolution works on
	// log10 of the magnitude, so a range of several decades is explored
	// evenly.
	Min, Max float64
	Sign     SignMode
	// Series snaps the magnitude to a preferred-number series: "" (none),
	// "E6", "E12", "E24", "E48" or "E96".
	Series string
	// Integer rounds the magnitude to a whole number (a turn count).
	Integer bool
}

// Label is a readable name for reports: "R1" or "dmod.is".
func (g Gene) Label() string {
	if g.Name != "" {
		return g.Name
	}
	if g.Element != "" {
		if g.Param != "" {
			return g.Element + "." + g.Param
		}
		return g.Element
	}
	return g.Model + "." + g.Param
}

func (g Gene) logRange() (lo, hi float64) { return math.Log10(g.Min), math.Log10(g.Max) }

// Validate checks the gene on its own and against the template circuit.
func (g Gene) Validate(t *circuit.Circuit) error {
	targets := 0
	for _, s := range []string{g.Element, g.Model} {
		if s != "" {
			targets++
		}
	}
	if targets > 1 || (targets == 0 && g.Name == "") {
		return fmt.Errorf("%w: %s: set exactly one of element, model or (with a Build function) name", ErrBadGene, g.Label())
	}
	if g.Model != "" && g.Param == "" {
		return fmt.Errorf("%w: %s: model gene needs a parameter", ErrBadGene, g.Label())
	}
	if !(g.Min > 0) || !(g.Max >= g.Min) || math.IsInf(g.Max, 0) {
		return fmt.Errorf("%w: %s: need 0 < min ≤ max, got [%g, %g]", ErrBadGene, g.Label(), g.Min, g.Max)
	}
	if g.Sign > Either {
		return fmt.Errorf("%w: %s: sign mode %d", ErrBadGene, g.Label(), g.Sign)
	}
	if _, err := series(g.Series); err != nil {
		return fmt.Errorf("%w: %s: %v", ErrBadGene, g.Label(), err)
	}
	if g.Integer && math.Floor(g.Max) < math.Ceil(g.Min) {
		return fmt.Errorf("%w: %s: no whole number in [%g, %g]", ErrBadGene, g.Label(), g.Min, g.Max)
	}
	if g.Element == "" && g.Model == "" {
		if t != nil {
			return fmt.Errorf("%w: %s: a free parameter gene needs Problem.Build", ErrBadGene, g.Label())
		}
		return nil
	}
	if t == nil {
		return fmt.Errorf("%w: %s: element and model genes need a template circuit", ErrBadGene, g.Label())
	}
	if g.Element != "" {
		i := t.Index(g.Element)
		if i < 0 {
			return fmt.Errorf("%w: %s: no such element", ErrBadGene, g.Label())
		}
		k := t.Elements[i].Kind
		if g.Param != "" {
			if !circuit.IsInstanceParam(k, strings.ToLower(g.Param)) {
				return fmt.Errorf("%w: %s: not an instance parameter of a %s", ErrBadGene, g.Label(), k)
			}
			if g.Sign != Positive {
				return fmt.Errorf("%w: %s: instance parameters must be positive", ErrBadGene, g.Label())
			}
			return nil
		}
		if g.Sign != Positive && k != circuit.Resistor && k != circuit.VoltageSource && k != circuit.CurrentSource {
			return fmt.Errorf("%w: %s: a %s value must be positive", ErrBadGene, g.Label(), k)
		}
		if k == circuit.Diode || k == circuit.BJT || k == circuit.MOSFET {
			return fmt.Errorf("%w: %s: evolve a %s through its model or instance parameters", ErrBadGene, g.Label(), k)
		}
		return nil
	}
	m, ok := t.Models[strings.ToLower(g.Model)]
	if !ok {
		return fmt.Errorf("%w: %s: no such model", ErrBadGene, g.Label())
	}
	_ = m
	if g.Sign != Positive {
		return fmt.Errorf("%w: %s: model parameters must be positive", ErrBadGene, g.Label())
	}
	return nil
}

// Genome is one candidate: log10 magnitudes (one per gene) and signs (+1 or
// −1, one per gene; fixed for Positive and Negative genes).
type Genome struct {
	X    []float64 `json:"x"`
	Sign []int8    `json:"sign"`
	// Parts is the topology genome (Problem.Topology); empty for value genomes.
	Parts []Part `json:"parts,omitempty"`
}

// Clone returns an independent copy.
func (g Genome) Clone() Genome {
	return Genome{X: append([]float64(nil), g.X...), Sign: append([]int8(nil), g.Sign...),
		Parts: append([]Part(nil), g.Parts...)}
}

// Space is the set of genes of one experiment.
type Space []Gene

// Validate checks every gene against the template; duplicate targets are
// rejected.
func (s Space) Validate(t *circuit.Circuit) error {
	if len(s) == 0 {
		return fmt.Errorf("%w: no genes", ErrBadGene)
	}
	seen := map[string]bool{}
	for _, g := range s {
		if err := g.Validate(t); err != nil {
			return err
		}
		key := strings.ToLower(g.Label())
		if seen[key] {
			return fmt.Errorf("%w: %s appears twice", ErrBadGene, g.Label())
		}
		seen[key] = true
	}
	return nil
}

// fixSign returns the sign a gene must have, or 0 if it is free.
func (g Gene) fixSign() int8 {
	switch g.Sign {
	case Negative:
		return -1
	case Either:
		return 0
	default:
		return 1
	}
}

// repair clamps magnitudes into range and forces fixed signs, so every
// genome produced by any operator decodes to a valid value.
func (s Space) repair(gn *Genome) {
	for i, g := range s {
		lo, hi := g.logRange()
		x := gn.X[i]
		if math.IsNaN(x) {
			x = (lo + hi) / 2
		}
		gn.X[i] = math.Min(math.Max(x, lo), hi)
		if f := g.fixSign(); f != 0 {
			gn.Sign[i] = f
		} else if gn.Sign[i] != -1 {
			gn.Sign[i] = 1
		}
	}
}

// Values returns the gene values a genome decodes to, after E-series
// snapping and with signs applied.
func (s Space) Values(gn Genome) []float64 {
	out := make([]float64, len(s))
	for i, g := range s {
		mag := math.Pow(10, gn.X[i])
		// At a bound, the bound itself (10^log10(x) can miss x by an ulp).
		if lo, hi := g.logRange(); gn.X[i] >= hi {
			mag = g.Max
		} else if gn.X[i] <= lo {
			mag = g.Min
		}
		ser, _ := series(g.Series)
		mag = snap(mag, ser, g.Min, g.Max)
		if g.Integer {
			mag = math.Min(math.Max(math.Round(mag), math.Ceil(g.Min)), math.Floor(g.Max))
		}
		out[i] = float64(gn.Sign[i]) * mag
	}
	return out
}

// Decode applies a genome to a copy of the template. The template is never
// modified, so many goroutines can decode from it at once.
func (s Space) Decode(t *circuit.Circuit, gn Genome) (*circuit.Circuit, []float64) {
	c := t.Clone()
	vals := s.Values(gn)
	for i, g := range s {
		if g.Element != "" {
			el := &c.Elements[c.Index(g.Element)]
			if g.Param != "" {
				if el.Params == nil {
					el.Params = map[string]float64{}
				}
				el.Params[strings.ToLower(g.Param)] = vals[i]
				continue
			}
			el.Value = vals[i]
			continue
		}
		c.Models[strings.ToLower(g.Model)].Params[strings.ToLower(g.Param)] = vals[i]
	}
	return c, vals
}
