// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// LTspice schematic (.asc) writer for the experiment's two topologies.
//
// Pin positions of LTspice's built-in symbols in orientation R0, relative to
// the symbol origin (from the symbol files res.asy, cap.asy, ind.asy,
// voltage.asy, current.asy shipped with LTspice): res and ind (16,16) and
// (16,96); cap (16,0) and (16,64); voltage (0,16) "+" and (0,96) "−";
// current (0,0) and (0,80), current flowing from the first pin through the
// source to the second. Coordinates are on LTspice's 16-unit grid.

type ascPart struct {
	sym   string // res, cap, ind, voltage, current
	name  string
	value string
}

// span is the distance between a symbol's two pins.
func (p ascPart) span() int {
	if p.sym == "cap" {
		return 64
	}
	return 80
}

// asc collects the sections of an LTspice schematic; they are written in
// the order LTspice itself uses (header, wires, flags, symbols, text).
type asc struct {
	head, wires, flags, syms, text strings.Builder
	cur                            *strings.Builder
}

func (a *asc) f(format string, x ...any) {
	b := a.cur
	if b == nil {
		b = &a.syms
	}
	fmt.Fprintf(b, format+"\n", x...)
}

func (a *asc) to(b *strings.Builder) *asc { a.cur = b; return a }

func (a *asc) wire(x1, y1, x2, y2 int) {
	if x1 != x2 || y1 != y2 {
		fmt.Fprintf(&a.wires, "WIRE %d %d %d %d\n", x1, y1, x2, y2)
	}
}

func (a *asc) String() string {
	return a.head.String() + a.wires.String() + a.flags.String() + a.syms.String() + a.text.String()
}

// symbol places p with its upper pin at (px, top); for a current source the
// upper pin is its second terminal (n−) and the lower pin its first (n+).
func (a *asc) symbol(p ascPart, px, top int) {
	switch p.sym {
	case "res", "ind":
		a.f("SYMBOL %s %d %d R0", p.sym, px-16, top-16)
	case "cap":
		a.f("SYMBOL cap %d %d R0", px-16, top)
	case "voltage":
		a.f("SYMBOL voltage %d %d R0", px, top-16)
	case "current":
		a.f("SYMBOL current %d %d R180", px, top+80)
	}
	a.f("SYMATTR InstName %s", p.name)
	a.f("SYMATTR Value %s", p.value)
}

// column places a series string from (x, top) down to (x, bottom).
func (a *asc) column(x, top, bottom int, parts []ascPart) {
	total := 0
	for _, p := range parts {
		total += p.span()
	}
	gap := (bottom - top - total) / (len(parts) + 1) / 16 * 16
	y := top
	for _, p := range parts {
		a.wire(x, y, x, y+gap)
		y += gap
		a.symbol(p, x, y)
		y += p.span()
	}
	a.wire(x, y, x, bottom)
}

func ascValue(v float64) string { return num(v) }

// schematicASC returns an LTspice schematic of a measured design with the
// same element values as its netlist. The .tran directive is active and the
// .ac directive is a comment: LTspice runs one analysis at a time, so swap
// the leading "!" and ";" to run the other.
func schematicASC(r Result, label string) (string, error) {
	d := r.Design
	c, err := d.Build(r.Rcoil, true)
	if err != nil {
		return "", err
	}
	val := func(name string) circuit.Element { return c.Elements[c.Index(name)] }
	var a asc
	a.to(&a.head).f("Version 4")
	a.f("SHEET 1 1360 720")
	a.to(&a.syms)
	const top, bottom = 96, 576
	xs := map[string]int{}
	var cols []int
	switch d.Variant {
	case Dry, Thermopile:
		xs["elem"] = 128
		a.column(128, top, bottom, []ascPart{
			{"res", "Rint", ascValue(val("Rint").Value)},
			{"voltage", "Vbial", wave(val("Vbial"))},
			{"cap", "Cplate", ascValue(val("Cplate").Value)},
		})
		cols = append(cols, 128)
	case Wet:
		split := 384
		a.column(160, top, split, []ascPart{
			{"res", "Relec", ascValue(val("Relec").Value)},
			{"voltage", "Vbial", wave(val("Vbial"))},
		})
		a.wire(96, split, 224, split)
		a.column(96, split, bottom, []ascPart{{"res", "Rct", ascValue(val("Rct").Value)}})
		a.column(224, split, bottom, []ascPart{{"cap", "Cdl", ascValue(val("Cdl").Value)}})
		a.column(400, top, bottom, []ascPart{{"cap", "Cgeo", ascValue(val("Cgeo").Value)}})
		cols = append(cols, 96, 160, 224, 400)
	}
	a.column(624, top, bottom, []ascPart{{"cap", "Cself", ascValue(val("Cself").Value)}})
	a.column(816, top, bottom, []ascPart{
		{"res", "Rcoil", ascValue(val("Rcoil").Value)},
		{"ind", "Lcoil", ascValue(val("Lcoil").Value)},
	})
	a.column(1008, top, bottom, []ascPart{{"current", "Iprobe", "AC 1"}})
	cols = append(cols, 624, 816, 1008)
	// Rails: the top rail joins every column's top (node a); the bottom
	// rail joins their bottoms (ground). The wet element's Rct and Cdl
	// columns start at the split, not at the top rail.
	topCols := []int{}
	for _, x := range cols {
		if d.Variant == Wet && (x == 96 || x == 224) {
			continue
		}
		topCols = append(topCols, x)
	}
	a.wire(topCols[0], top, topCols[len(topCols)-1], top)
	a.wire(cols[0], bottom, cols[len(cols)-1], bottom)
	a.to(&a.flags).f("FLAG %d %d a", topCols[0], top)
	a.f("FLAG %d %d 0", 512, bottom)
	a.to(&a.text)
	periods := r.Q * 3
	if periods < 20 {
		periods = 20
	}
	if periods > 2000 {
		periods = 2000
	}
	a.f("TEXT 64 656 Left 2 !.tran 0 %s 0 %s", num(periods/r.F0), num(1/(40*r.F0)))
	a.f("TEXT 64 688 Left 2 ;.ac dec 200 %s %s", num(r.F0/100), num(r.F0*100))
	title := fmt.Sprintf("Bi-Al %s element with %s coil, %s. f0 = %s, Q = %.4g (kronspice examples/bial_bifilar)",
		strings.Fields(d.Variant.String())[0], coilKind(d.Coil.Bifilar), label,
		strings.ReplaceAll(eng(r.F0, "Hz"), "µ", "u"), r.Q)
	a.f("TEXT 64 24 Left 2 ;%s", title)
	return a.String(), nil
}

// wave formats a source's value for LTspice.
func wave(e circuit.Element) string {
	if e.Wave.Kind == circuit.WaveNone {
		return ascValue(e.Value)
	}
	args := make([]string, len(e.Wave.Args))
	for i, x := range e.Wave.Args {
		args[i] = num(x)
	}
	return e.Wave.Kind.String() + "(" + strings.Join(args, " ") + ")"
}
