// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"sort"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
)

// ascNets rebuilds the netlist an .asc file describes: wires join points,
// symbol pins sit at the offsets documented in asc.go, and FLAGs name nets.
// It returns, for each instance, its two nets (first pin, second pin).
func ascNets(t *testing.T, text string) map[string][2]string {
	t.Helper()
	parent := map[[2]int][2]int{}
	var find func(p [2]int) [2]int
	find = func(p [2]int) [2]int {
		if q, ok := parent[p]; ok && q != p {
			r := find(q)
			parent[p] = r
			return r
		}
		parent[p] = p
		return p
	}
	union := func(a, b [2]int) { parent[find(a)] = find(b) }
	type wire struct{ x1, y1, x2, y2 int }
	var wires []wire
	pins := map[string][2][2]int{}
	flags := map[[2]int]string{}
	var cur string
	var curPins [2][2]int
	for _, line := range strings.Split(text, "\n") {
		f := strings.Fields(line)
		if len(f) == 0 {
			continue
		}
		switch f[0] {
		case "WIRE":
			var w wire
			fmt.Sscan(strings.Join(f[1:5], " "), &w.x1, &w.y1, &w.x2, &w.y2)
			wires = append(wires, w)
		case "FLAG":
			var x, y int
			fmt.Sscan(f[1]+" "+f[2], &x, &y)
			flags[[2]int{x, y}] = f[3]
		case "SYMBOL":
			var x, y int
			fmt.Sscan(f[2]+" "+f[3], &x, &y)
			switch f[1] {
			case "res", "ind":
				curPins = [2][2]int{{x + 16, y + 16}, {x + 16, y + 96}}
			case "cap":
				curPins = [2][2]int{{x + 16, y}, {x + 16, y + 64}}
			case "voltage":
				curPins = [2][2]int{{x, y + 16}, {x, y + 96}}
			case "current": // R180: pins (0,0) and (0,−80)
				curPins = [2][2]int{{x, y}, {x, y - 80}}
			}
		case "SYMATTR":
			if f[1] == "InstName" {
				cur = f[2]
				pins[cur] = curPins
			}
		}
	}
	// A wire joins its two ends and every pin, flag or wire end lying on it.
	points := map[[2]int]bool{}
	for _, w := range wires {
		points[[2]int{w.x1, w.y1}] = true
		points[[2]int{w.x2, w.y2}] = true
	}
	for _, p := range pins {
		points[p[0]] = true
		points[p[1]] = true
	}
	for p := range flags {
		points[p] = true
	}
	on := func(w wire, p [2]int) bool {
		if w.x1 == w.x2 {
			lo, hi := min(w.y1, w.y2), max(w.y1, w.y2)
			return p[0] == w.x1 && p[1] >= lo && p[1] <= hi
		}
		lo, hi := min(w.x1, w.x2), max(w.x1, w.x2)
		return p[1] == w.y1 && p[0] >= lo && p[0] <= hi
	}
	for _, w := range wires {
		for p := range points {
			if on(w, p) {
				union(p, [2]int{w.x1, w.y1})
			}
		}
	}
	name := func(p [2]int) string {
		r := find(p)
		for fp, n := range flags {
			if find(fp) == r {
				return n
			}
		}
		return fmt.Sprintf("n%d_%d", r[0], r[1])
	}
	out := map[string][2]string{}
	for inst, p := range pins {
		out[inst] = [2]string{name(p[0]), name(p[1])}
	}
	return out
}

// The LTspice schematic connects every element to the same nodes as the
// SPICE netlist (up to the names of internal nodes).
func TestASCMatchesNetlist(t *testing.T) {
	for _, v := range variants {
		d := designFrom(v, DefaultAssumptions(), true, append(baseX(), 100, 0.5e-3, 10e-3))
		r, err := Measure(d, analysis.Options{})
		if err != nil && v != Wet {
			t.Fatal(err)
		}
		if r.F0 == 0 { // wet test design has no resonance: still draw it
			r.F0, r.Q, r.Rcoil = 1e3, 1, d.Coil.Resistance(0)
		}
		text, err := schematicASC(r, "test")
		if err != nil {
			t.Fatal(err)
		}
		got := ascNets(t, text)
		c, _ := d.Build(r.Rcoil, true)
		// Map netlist node names to schematic net names consistently.
		mapping := map[string]string{"0": "0", "a": "a"}
		for _, e := range c.Elements {
			nets, ok := got[e.Name]
			if !ok {
				t.Fatalf("%s: %s missing from the schematic", v, e.Name)
			}
			for k := 0; k < 2; k++ {
				n := c.NodeNames[e.Nodes[k]]
				if m, ok := mapping[n]; ok {
					if m != nets[k] {
						t.Errorf("%s: %s pin %d on %s, want %s (node %s)", v, e.Name, k+1, nets[k], m, n)
					}
				} else {
					mapping[n] = nets[k]
				}
			}
		}
		// Distinct netlist nodes must be distinct nets.
		seen := map[string]string{}
		keys := make([]string, 0, len(mapping))
		for k := range mapping {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		for _, n := range keys {
			if other, dup := seen[mapping[n]]; dup {
				t.Errorf("%s: nodes %s and %s are shorted in the schematic", v, other, n)
			}
			seen[mapping[n]] = n
		}
		if len(got) != len(c.Elements) {
			t.Errorf("%s: schematic has %d instances, netlist %d", v, len(got), len(c.Elements))
		}
	}
}
