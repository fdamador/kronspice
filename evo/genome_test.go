// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"errors"
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func TestSeriesTables(t *testing.T) {
	if len(e6) != 6 || len(e12) != 12 || len(e24) != 24 || len(e48) != 48 || len(e96) != 96 {
		t.Fatal("series lengths")
	}
	// E48 and E96 follow round(10^(i/n), 3 digits) exactly.
	for i, m := range e96 {
		if want := roundSig(math.Pow(10, float64(i)/96), 3); m != want {
			t.Errorf("E96[%d] = %g, formula %g", i, m, want)
		}
	}
	in := func(x float64, s []float64) bool {
		for _, v := range s {
			if v == x {
				return true
			}
		}
		return false
	}
	for _, m := range e6 {
		if !in(m, e12) {
			t.Errorf("E6 %g not in E12", m)
		}
	}
	for _, m := range e12 {
		if !in(m, e24) {
			t.Errorf("E12 %g not in E24", m)
		}
	}
}

func TestSnap(t *testing.T) {
	cases := []struct {
		v, want  float64
		ser      []float64
		min, max float64
	}{
		{4600, 4700, e12, 1, 1e9},
		{4.9e-9, 4.7e-9, e12, 1e-12, 1},
		{9.7e3, 1e4, e24, 1, 1e9},     // rounds up into the next decade
		{1234, 1240, e96, 1, 1e9},     // 1.24
		{50, 100, e12, 100, 1e3},      // clamped into range first
		{2e3, 1.8e3, e12, 100, 1.9e3}, // nearest in-range value
		{1234, 1234, nil, 1, 1e9},
	}
	for _, c := range cases {
		if got := snap(c.v, c.ser, c.min, c.max); math.Abs(got-c.want) > 1e-12*c.want {
			t.Errorf("snap(%g) = %.15g, want %g", c.v, got, c.want)
		}
	}
}

func signedTemplate(t *testing.T) *circuit.Circuit {
	c, err := circuit.NewBuilder("t").VoltageSource("V1", "in", "0", 1).
		Resistor("R1", "in", "out", 1e3).Resistor("R2", "out", "0", 1e3).
		Model("dm", "d", map[string]float64{"is": 1e-14}).Diode("D1", "out", "0", "dm").Build()
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func TestDecodeSignsAndModels(t *testing.T) {
	tpl := signedTemplate(t)
	sp := Space{
		{Element: "R1", Min: 100, Max: 1e5, Series: "E24"},
		{Element: "R2", Min: 100, Max: 1e5, Sign: Negative},
		{Element: "V1", Min: 0.1, Max: 10, Sign: Either},
		{Model: "dm", Param: "is", Min: 1e-16, Max: 1e-12},
	}
	if err := sp.Validate(tpl); err != nil {
		t.Fatal(err)
	}
	g := Genome{X: []float64{math.Log10(4600), 3, 0, -13}, Sign: []int8{-1, 1, -1, 1}}
	sp.repair(&g) // forces R1 positive and R2 negative
	c, vals := sp.Decode(tpl, g)
	want := []float64{4700, -1000, -1, 1e-13}
	for i := range want {
		if math.Abs(vals[i]-want[i]) > 1e-12*math.Abs(want[i]) {
			t.Fatalf("values = %v, want %v", vals, want)
		}
	}
	if c.Elements[1].Value != 4700 || c.Elements[2].Value != -1000 || c.Elements[0].Value != -1 ||
		c.Models["dm"].Params["is"] != 1e-13 {
		t.Fatalf("decoded circuit = %+v", c)
	}
	if tpl.Elements[1].Value != 1e3 || tpl.Models["dm"].Params["is"] != 1e-14 {
		t.Fatal("template modified")
	}
	if err := c.Validate(); err != nil {
		t.Fatalf("decoded circuit invalid: %v", err)
	}
}

func TestGeneValidation(t *testing.T) {
	tpl := signedTemplate(t)
	bad := []Gene{
		{Element: "R9", Min: 1, Max: 2},
		{Element: "R1", Min: 0, Max: 2},
		{Element: "R1", Min: 3, Max: 2},
		{Element: "R1", Model: "dm", Param: "is", Min: 1, Max: 2},
		{Model: "dm", Min: 1, Max: 2},
		{Model: "zz", Param: "is", Min: 1, Max: 2},
		{Model: "dm", Param: "is", Min: 1e-15, Max: 1e-12, Sign: Either},
		{Element: "D1", Min: 1, Max: 2},
		{Element: "R1", Min: 1, Max: 2, Series: "E7"},
	}
	for _, g := range bad {
		if err := (Space{g}).Validate(tpl); err == nil {
			t.Errorf("%+v accepted", g)
		}
	}
	if err := (Space{{Element: "R1", Min: 1, Max: 2}, {Element: "r1", Min: 1, Max: 2}}).Validate(tpl); err == nil {
		t.Error("duplicate gene accepted")
	}
}

// Element + Param genes set a MOSFET instance parameter; the element's own
// value and non-instance parameters are rejected.
func TestInstanceParamGene(t *testing.T) {
	tpl := deck(t, "m\nVDD d 0 5\nVG g 0 2\nR1 d x 1k\nM1 x g 0 0 nm L=1u W=10u\n.model nm NMOS(LEVEL=1)\n.op\n")
	s := Space{{Element: "M1", Param: "W", Min: 1e-6, Max: 1e-4}}
	if err := s.Validate(tpl); err != nil {
		t.Fatal(err)
	}
	if s[0].Label() != "M1.W" {
		t.Errorf("label %q", s[0].Label())
	}
	c, vals := s.Decode(tpl, Genome{X: []float64{-5}, Sign: []int8{1}})
	if m := c.Elements[c.Index("M1")]; m.Params["w"] != vals[0] || m.Params["l"] != 1e-6 {
		t.Errorf("decoded params %v, value %g", m.Params, vals[0])
	}
	if tpl.Elements[tpl.Index("M1")].Params["w"] != 10e-6 {
		t.Error("template modified")
	}
	for _, g := range []Gene{
		{Element: "M1", Min: 1, Max: 2},
		{Element: "M1", Param: "vto", Min: 1, Max: 2},
		{Element: "R1", Param: "w", Min: 1, Max: 2},
		{Element: "M1", Param: "w", Min: 1e-6, Max: 1e-4, Sign: Either},
	} {
		if err := g.Validate(tpl); !errors.Is(err, ErrBadGene) {
			t.Errorf("%+v: %v", g, err)
		}
	}
}
