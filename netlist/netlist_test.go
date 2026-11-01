// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package netlist

import (
	"errors"
	"math"
	"testing"
)

func TestParseValue(t *testing.T) {
	cases := []struct {
		in   string
		want float64
	}{
		{"1", 1}, {"-2.5", -2.5}, {"1k", 1e3}, {"1K", 1e3}, {"4.7kohm", 4.7e3},
		{"1meg", 1e6}, {"1MEG", 1e6}, {"1m", 1e-3}, {"1M", 1e-3}, {"2u", 2e-6},
		{"3n", 3e-9}, {"4p", 4e-12}, {"5f", 5e-15}, {"1g", 1e9}, {"1t", 1e12},
		{"1mil", 25.4e-6}, {"1e3", 1e3}, {"1.5e-3", 1.5e-3}, {"5V", 5}, {".5", 0.5},
		{"1e", 1}, // "e" without digits is a unit letter, as in SPICE
	}
	for _, c := range cases {
		got, err := ParseValue(c.in)
		if err != nil {
			t.Errorf("ParseValue(%q) error: %v", c.in, err)
			continue
		}
		if math.Abs(got-c.want) > 1e-12*math.Max(1, math.Abs(c.want)) {
			t.Errorf("ParseValue(%q) = %g, want %g", c.in, got, c.want)
		}
	}
	for _, bad := range []string{"", "abc", "k", "--1", "1e999", "1%"} {
		if _, err := ParseValue(bad); !errors.Is(err, ErrBadValue) {
			t.Errorf("ParseValue(%q) err = %v, want ErrBadValue", bad, err)
		}
	}
}

func TestParseDivider(t *testing.T) {
	src := `voltage divider
* a comment
V1 in 0 DC 10 ; inline comment
R1 in out 1k
R2 out 0
+ 1k
.op
.end
R3 ignored after end 5
`
	d, err := ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	if d.Title != "voltage divider" {
		t.Errorf("title = %q", d.Title)
	}
	if len(d.Elements) != 3 {
		t.Fatalf("got %d elements, want 3", len(d.Elements))
	}
	if r2 := d.Elements[2]; r2.Name != "R2" || r2.Value != 1000 || r2.Nodes[0] != "out" {
		t.Errorf("R2 = %+v", r2)
	}
	if len(d.Analyses) != 1 || d.Analyses[0].Kind != "op" {
		t.Errorf("analyses = %+v", d.Analyses)
	}
}

func TestParseErrors(t *testing.T) {
	bad := []string{
		"t\nR1 a 0\n",                        // missing value
		"t\nJ1 d g s model\n",                // unsupported element (for now)
		"t\nM1 d g s model\n",                // MOSFET needs four nodes and a model
		"t\nM1 d g s b nm L 1u\n",            // MOSFET parameter without =
		"t\nM1 d g s b nm L=1u L=2u\n",       // parameter given twice
		"t\nQ1 c b\n",                        // BJT needs c b e and a model
		"t\nE1 a 0 b\n",                      // E needs nc+ nc- gain
		"t\nF1 a 0 V1\n",                     // F needs a gain
		"t\nK1 L1 0.5\n",                     // K needs two inductors and k
		"t\nK1 L1 L2 tight\n",                // bad coupling value
		"t\nV1 a 0\n",                        // source without a value
		"t\n.tran 1u\n",                      // .tran without a stop time
		"t\n.tran 1u 0\n",                    // stop time must be positive
		"t\nC1 a 0 1u IC 2\n",                // IC needs "="
		"t\nR1 a 0 1k IC=1\n",                // resistors take no IC
		"t\nV1 a 0 PULSE(0 1) SIN(0 1 1k)\n", // two waveforms
		"t\n.ac dec 10 0 1k\n",               // log sweep needs fstart > 0
		"t\n.ac log 10 1 1k\n",               // unknown sweep type
		"t\n.dc V1 0 5 -1\n",                 // step moves away from stop
		"t\n.noise v(out) V1 dec 10 1 1k\n",  // unsupported analysis (for now)
		"t\nR1 a 0 notanumber\n",             // bad value
	}
	for _, src := range bad {
		if _, err := ParseString(src); err == nil {
			t.Errorf("expected error for %q", src)
		}
	}
}

func TestParseTransientDeck(t *testing.T) {
	src := `rc step
V1 in 0 PULSE(0 5 1u 10n 10n 1m 2m)
V2 a 0 DC 1 AC 1 90 SIN 0 1 1k
I1 0 b pwl(0 0, 1m 1m, 2m 0)
R1 in out 1k
C1 out 0 1u IC=2.5
L1 a b 10m ic = 1m
.tran 1u 5m 0 2u UIC
.end
`
	d, err := ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	v1 := d.Elements[0]
	if v1.WaveKind != "pulse" || len(v1.WaveArgs) != 7 || v1.WaveArgs[1] != 5 || v1.WaveArgs[3] != 10e-9 || v1.HasDC {
		t.Errorf("V1 = %+v", v1)
	}
	v2 := d.Elements[1]
	if !v2.HasDC || v2.Value != 1 || v2.ACMag != 1 || v2.ACPhase != 90 || v2.WaveKind != "sin" || len(v2.WaveArgs) != 3 {
		t.Errorf("V2 = %+v", v2)
	}
	if i1 := d.Elements[2]; i1.WaveKind != "pwl" || len(i1.WaveArgs) != 6 || i1.WaveArgs[3] != 1e-3 {
		t.Errorf("I1 = %+v", i1)
	}
	if c1 := d.Elements[4]; c1.Value != 1e-6 || !c1.HasIC || c1.IC != 2.5 {
		t.Errorf("C1 = %+v", c1)
	}
	if l1 := d.Elements[5]; l1.Value != 10e-3 || !l1.HasIC || l1.IC != 1e-3 {
		t.Errorf("L1 = %+v", l1)
	}
	a := d.Analyses[0]
	if a.Kind != "tran" || a.Step != 1e-6 || a.Stop != 5e-3 || a.Start != 0 || a.MaxStep != 2e-6 || !a.UIC {
		t.Errorf("tran = %+v", a)
	}
}

func TestParseACDC(t *testing.T) {
	d, err := ParseString("t\nV1 a 0 DC 1 AC 1\nR1 a 0 1k\n.ac dec 10 1 1meg\n.dc V1 0 5 0.5 I1 1m 2m 0.5m\n")
	if err != nil {
		t.Fatal(err)
	}
	ac, dc := d.Analyses[0], d.Analyses[1]
	if ac.Kind != "ac" || ac.Variation != "dec" || ac.Points != 10 || ac.FStart != 1 || ac.FStop != 1e6 {
		t.Errorf("ac = %+v", ac)
	}
	if dc.Kind != "dc" || len(dc.Sweeps) != 2 || dc.Sweeps[0].Source != "V1" || dc.Sweeps[0].Step != 0.5 || dc.Sweeps[1].Stop != 2e-3 {
		t.Errorf("dc = %+v", dc)
	}
}

func FuzzParse(f *testing.F) {
	f.Add("t\nV1 1 0 PULSE(0 1 0 1n 1n 1u 2u)\nC1 1 0 1p IC=1\n.tran 1n 1u\n")
	f.Add("t\nV1 1 0 5\nR1 1 0 1k\n.op\n.end\n")
	f.Add("t\n+ continuation first\nR1 a b 1meg ; c\n")
	f.Add("")
	f.Add("t\nC1 a 0 Q=pwl(-1 -1u 0 0 1 1u) SCALE=SIN(1 0.2 1k) EXT=thermal IC=0\nL1 a 0 FLUX=pwl(0 0 1 1m) SCALE=pwl(0 1 1m 2)\n")
	f.Fuzz(func(t *testing.T, src string) {
		d, err := ParseString(src)
		if err != nil {
			return
		}
		for _, e := range d.Elements {
			if math.IsNaN(e.Value) || math.IsInf(e.Value, 0) {
				t.Fatalf("non-finite value accepted: %+v", e)
			}
		}
	})
}
