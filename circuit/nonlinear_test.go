// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"errors"
	"reflect"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/netlist"
)

// Nonlinear and time-varying capacitors and inductors parse into the new
// kinds, keep plain C and L lines linear, and survive a write and re-parse.
func TestNonlinearReactiveRoundTrip(t *testing.T) {
	src := `nonlinear reactive
V1 in 0 PULSE(0 1 0 1u 1u 1m 2m)
R1 in a 1k
CN a 0 Q=pwl(-1 -1u, 0 0, 1 1u, 2 1.5u) IC=0.5
CV a b 1u SCALE=SIN(1 0.2 10k) EXT=thermal
LN b 0 FLUX=pwl(-1 -1m 0 0 1 1m 2 1.2m) SCALE=PWL(0 1 1m 2)
LV b c 10u SCALE=pulse(1 2 0 1u 1u 5u 10u) IC=0.1
L1 c 0 1m
C1 c 0 1n IC=0.25
.end
`
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]Kind{"CN": NonlinearCapacitor, "CV": NonlinearCapacitor, "LN": NonlinearInductor,
		"LV": NonlinearInductor, "L1": Inductor, "C1": Capacitor}
	for name, k := range want {
		if got := c.Elements[c.Index(name)].Kind; got != k {
			t.Errorf("%s: kind %v, want %v", name, got, k)
		}
	}
	cn := c.Elements[c.Index("CN")]
	if !reflect.DeepEqual(cn.Curve, []float64{-1, -1e-6, 0, 0, 1, 1e-6, 2, 1.5e-6}) || !cn.HasIC || cn.IC != 0.5 {
		t.Errorf("CN = %+v", cn)
	}
	cv := c.Elements[c.Index("CV")]
	if cv.Value != 1e-6 || cv.Scale.Kind != WaveSin || cv.Channel != exchange.Thermal {
		t.Errorf("CV = %+v", cv)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	d2, err := netlist.ParseString(sb.String())
	if err != nil {
		t.Fatalf("%v\n%s", err, sb.String())
	}
	c2, err := FromDeck(d2)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(c, c2) {
		t.Fatalf("round trip differs:\n%s\n%+v\n%+v", sb.String(), c, c2)
	}
	if cl := c.Clone(); &cl.Elements[c.Index("LN")].Scale.Args[0] == &c.Elements[c.Index("LN")].Scale.Args[0] {
		t.Error("Clone shares the scale waveform")
	}
}

func TestNonlinearReactiveValidation(t *testing.T) {
	sin := func(args ...float64) Waveform { return Waveform{Kind: WaveSin, Args: args} }
	base := func() *Builder { return NewBuilder("v").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "b", 1) }
	bad := map[string]*Builder{
		"curve misses the origin": base().NonlinearCapacitor("C1", "b", "0", []float64{-1, -1, 1, 2}, Waveform{}),
		"falling charge":          base().NonlinearCapacitor("C1", "b", "0", []float64{-1, 1, 0, 0, 1, -1}, Waveform{}),
		"flat flux":               base().NonlinearInductor("L1", "b", "0", []float64{-1, -1, 0, 0, 1, 1, 2, 1}, Waveform{}),
		"one point":               base().NonlinearCapacitor("C1", "b", "0", []float64{0, 0}, Waveform{}),
		"zero varying C":          base().VaryingCapacitor("C1", "b", "0", 0, sin(1, 0.5, 1e3)),
		"scale reaches 0":         base().VaryingCapacitor("C1", "b", "0", 1e-6, sin(1, 1, 1e3)),
		"growing scale":           base().VaryingInductor("L1", "b", "0", 1e-3, sin(2, 0.5, 1e3, 0, -10)),
		"scale jump":              base().VaryingInductor("L1", "b", "0", 1e-3, Waveform{Kind: WavePWL, Args: []float64{0, 1, 1e-3, 1, 1e-3, 2}}),
		"value and curve": base().Add(Element{Name: "C1", Kind: NonlinearCapacitor, Value: 1e-6,
			Curve: []float64{-1, -1, 1, 1}}, "b", "0"),
		"neutrino port": base().Add(Element{Name: "C1", Kind: NonlinearCapacitor, Value: 1e-6,
			Scale: sin(1, 0.5, 1e3), Channel: exchange.Neutrino}, "b", "0"),
		"K on a nonlinear inductor": base().VaryingInductor("L1", "b", "0", 1e-3, sin(1, 0.5, 1e3)).
			Inductor("L2", "b", "0", 1e-3).Mutual("K1", "L1", "L2", 0.5),
	}
	for name, b := range bad {
		_, err := b.Build()
		if !errors.Is(err, ErrBadValue) && !errors.Is(err, ErrBadWaveform) {
			t.Errorf("%s: err = %v, want ErrBadValue or ErrBadWaveform", name, err)
		}
	}
	// A nonlinear capacitor is no DC path; a nonlinear inductor is.
	if _, err := NewBuilder("p").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "0", 1).
		NonlinearCapacitor("C1", "b", "0", []float64{-1, -1, 1, 1}, Waveform{}).Build(); !errors.Is(err, ErrFloatingNode) {
		t.Errorf("capacitor-only node: err = %v, want ErrFloatingNode", err)
	}
	if _, err := NewBuilder("p").VoltageSource("V1", "a", "0", 1).VaryingInductor("L1", "a", "b", 1e-3, sin(1, 0.5, 1e3)).
		Resistor("R1", "b", "0", 1).Build(); err != nil {
		t.Errorf("inductor path: %v", err)
	}
}

func TestNonlinearReactiveParseErrors(t *testing.T) {
	for _, line := range []string{
		"C1 a 0 Q=1u",                        // Q= needs pwl(...)
		"C1 a 0 FLUX=pwl(0 0 1 1)",           // FLUX= is for inductors
		"L1 a 0 1m SCALE=1",                  // SCALE needs a waveform
		"C1 a 0 1u IC=1 IC=2",                // repeated option
		"C1 a 0 Q=pwl(0 0 1)",                // odd curve
		"R1 a 0 1k SCALE=SIN(1 0.5 1k)",      // resistors take EXT= only
		"L1 a 0 FLUX=pwl(-1 -1 1 1) BOGUS=2", // unknown option
	} {
		if _, err := netlist.ParseString("t\nV1 a 0 1\n" + line + "\n.end\n"); err == nil {
			t.Errorf("%q: no error", line)
		}
	}
}
