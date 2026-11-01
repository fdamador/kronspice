// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"bytes"
	"errors"
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/netlist"
)

const machineDeck = `machine round trip
.domain mechanical w
.model mm machine (n=2 theta0=0.25 l_1_1=10m l_2_2=20m lc_1_1_2=5m ls_2_1_1=1m g_1_2=0.5m)
I1 0 a 1
I2 0 b 2
N1 a 0 b 0 w 0 mm
RD w 0 4
.end
`

func TestMachineNetlistRoundTrip(t *testing.T) {
	d, err := netlist.ParseString(machineDeck)
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	s, err := c.MachineFor(c.Elements[2])
	if err != nil {
		t.Fatal(err)
	}
	if s.N != 2 || s.Theta0 != 0.25 || s.L0[1][1] != 20e-3 || s.G[0][1] != 0.5e-3 || s.G[1][0] != 0 {
		t.Errorf("spec %+v", s)
	}
	l, d1 := [][]float64{{0, 0}, {0, 0}}, [][]float64{{0, 0}, {0, 0}}
	s.Inductance(0.3, l, d1, nil)
	if want := 5e-3 * math.Cos(0.3); math.Abs(l[0][1]-want) > 1e-15 || l[1][0] != l[0][1] {
		t.Errorf("L12(0.3) = %g, want %g", l[0][1], want)
	}
	if want := 10e-3 + 1e-3*math.Sin(0.6); math.Abs(l[0][0]-want) > 1e-15 {
		t.Errorf("L11(0.3) = %g, want %g", l[0][0], want)
	}
	if want := 2e-3 * math.Cos(0.6); math.Abs(d1[0][0]-want) > 1e-15 {
		t.Errorf("L11′(0.3) = %g, want %g", d1[0][0], want)
	}
	var buf bytes.Buffer
	if err := c.WriteNetlist(&buf); err != nil {
		t.Fatal(err)
	}
	d2, err := netlist.ParseString(buf.String())
	if err != nil {
		t.Fatalf("%v\n%s", err, buf.String())
	}
	c2, err := FromDeck(d2)
	if err != nil {
		t.Fatalf("%v\n%s", err, buf.String())
	}
	s2, _ := c2.MachineFor(c2.Elements[2])
	if s2.Params()["lc_1_1_2"] != s.Params()["lc_1_1_2"] || len(s2.Params()) != len(s.Params()) {
		t.Errorf("round trip changed the spec:\n%s", buf.String())
	}
}

func TestMachineValidation(t *testing.T) {
	for _, tc := range []struct{ src, want string }{
		{strings.Replace(machineDeck, "l_2_2=20m", "l_2_2=1m", 1), "not positive definite"},
		{strings.Replace(machineDeck, "n=2", "n=3 l_3_3=1m", 1), "nodes"},
		{strings.Replace(machineDeck, "g_1_2", "g_1_7", 1), "out of range"},
		{strings.Replace(machineDeck, "g_1_2", "q_1_2", 1), "unknown parameter"},
		{strings.Replace(machineDeck, ".domain mechanical w\n", "", 1), "not mechanical"},
		{strings.Replace(machineDeck, "N1 a 0 b 0 w 0 mm", "N1 a 0 b 0 w 0 nomodel", 1), "undefined model"},
	} {
		d, err := netlist.ParseString(tc.src)
		if err == nil {
			_, err = FromDeck(d)
		}
		if err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("want %q, got %v", tc.want, err)
		}
		if err != nil && !errors.Is(err, ErrBadModel) && !errors.Is(err, ErrBadValue) {
			t.Errorf("%v is neither ErrBadModel nor ErrBadValue", err)
		}
	}
}
