// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kicad

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Every netlist in the repository converts to a schematic whose geometry
// alone (pins, wires, labels, GND symbols) gives back the netlist's
// connectivity, element by element and pin by pin; the S-expression is
// well formed; and conversion is deterministic (same bytes twice).
func TestRoundTrip(t *testing.T) {
	var files []string
	for _, g := range []string{"../testdata/netlists/*.cir", "../examples/*/*.cir", "../examples/*/netlists/*.cir"} {
		m, _ := filepath.Glob(g)
		files = append(files, m...)
	}
	if len(files) < 20 {
		t.Fatalf("only %d netlists found", len(files))
	}
	for _, f := range files {
		src, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		p, err := Split(string(src))
		if err != nil {
			t.Errorf("%s: %v", f, err)
			continue
		}
		var a, b bytes.Buffer
		if err := Write(&a, string(src), Options{}); err != nil {
			t.Errorf("%s: %v", f, err)
			continue
		}
		if err := Write(&b, string(src), Options{}); err != nil || !bytes.Equal(a.Bytes(), b.Bytes()) {
			t.Errorf("%s: output not deterministic", f)
		}
		got, err := Read(a.String())
		if err != nil {
			t.Errorf("%s: %v", f, err)
			continue
		}
		for _, e := range p.Elements {
			nets := got[e.Name]
			if len(nets) != len(e.Nodes) {
				t.Errorf("%s %s: %d pins in schematic, %d nodes", f, e.Name, len(nets), len(e.Nodes))
				continue
			}
			for i, n := range e.Nodes {
				want := n
				if strings.EqualFold(n, "gnd") {
					want = "0"
				}
				if nets[i] != want {
					t.Errorf("%s %s pin %d: net %q, node %q", f, e.Name, i+1, nets[i], want)
				}
			}
		}
		if len(got) != len(p.Elements) {
			t.Errorf("%s: %d symbols, %d elements", f, len(got), len(p.Elements))
		}
	}
	t.Logf("%d netlists", len(files))
}

// The raw-SPICE simulation fields carry each element's model text, and
// dot cards become text directives.
func TestSimFields(t *testing.T) {
	src := "t\nV1 in 0 DC 1 AC 1\nR1 in out 1k\nC1 out 0 1u\nE1 a 0 out 0 2\nRL a 0 1k\nK1 L1 L2 0.5\nL1 x 0 1m\nL2 y 0 1m\nRx x y 1\nRy y 0 1\n.ac dec 10 1 1meg\n.end\n"
	var b bytes.Buffer
	if err := Write(&b, src, Options{}); err != nil {
		t.Fatal(err)
	}
	s := b.String()
	for _, want := range []string{
		`(property "Sim.Params" "type=\"V\" model=\"DC 1 AC 1\" lib=\"\""`,
		`(property "Sim.Params" "type=\"E\" model=\"2\" lib=\"\""`,
		`(property "Sim.Params" "type=\"K\" model=\"L1 L2 0.5\" lib=\"\""`,
		`(property "Sim.Pins" "1=1 2=2 3=3 4=4"`,
		`(text ".ac dec 10 1 1meg"`,
		`(lib_id "power:GND")`,
	} {
		if !strings.Contains(s, want) {
			t.Errorf("missing %s", want)
		}
	}
}

// Pin counts by element kind, from the SPICE3 element syntax.
func TestPinCounts(t *testing.T) {
	src := "t\nV1 in 0 1\nR1 in a 1k\nE1 b 0 a 0 2\nG1 c 0 a 0 1m\nF1 d 0 V1 3\nH1 e 0 V1 100\n" +
		"RB b 0 1\nRC c 0 1\nRD d 0 1\nRE e 0 1\nRA a 0 1\n.end\n"
	p, err := Split(src)
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]int{"V1": 2, "R1": 2, "E1": 4, "G1": 4, "F1": 2, "H1": 2}
	for _, e := range p.Elements {
		if n, ok := want[e.Name]; ok && len(e.Nodes) != n {
			t.Errorf("%s: %d nodes %v, want %d (model %q)", e.Name, len(e.Nodes), e.Nodes, n, e.Model)
		}
	}
}
