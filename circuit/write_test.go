// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"reflect"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/netlist"
)

func TestWriteNetlistRoundTrip(t *testing.T) {
	src := `round trip
.model dz D(IS=1e-14 BV=5.1 RS=2)
V1 in 0 DC 1.5 AC 1 30 SIN(0 2 1k 1u)
I1 0 b PWL(0 0 1m 2e-3)
R1 in a 1.2345678901234k
RN a 0 -470
C1 a b 3.3n IC=0.25
L1 b 0 10u
L2 a 0 22u
K1 L1 L2 0.75
D1 0 a dz 2
.tran 1u 1m
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
	var sb strings.Builder
	if err := c.WriteNetlist(&sb, ".tran 1u 1m"); err != nil {
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
}

// Transistors (with and without substrate and area) and the four
// controlled sources survive a write and re-parse unchanged.
func TestWriteNetlistRoundTripM7(t *testing.T) {
	src := `m7 round trip
.model qn NPN(IS=2e-15 BF=150 VAF=80 RB=100 CJE=1p)
.model qp PNP(BF=60)
VCC vcc 0 12
VS a m 0
R1 vcc b 47k
R2 b 0 10k
RA m 0 1k
Q1 c b e qn
Q2 c2 b e2 sub qn 2
Q3 0 b vcc qp 0.5
RC vcc c 3k
RE e 0 1k
R3 vcc c2 1k
R4 e2 0 1k
RS sub 0 1
E1 o1 0 b 0 2.5
G1 0 o2 c e 1m
F1 0 o3 VS 3
H1 o4 0 VS 500
RO1 o1 0 1k
RO2 o2 0 1k
RO3 o3 0 1k
RO4 o4 0 1k
V2 vcc a 1
.op
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
	q2 := c.Elements[c.Index("q2")]
	if q2.Area != 2 || len(q2.Extra) != 2 || c.NodeNames[q2.Extra[1]] != "sub" || c.NodeNames[q2.Nodes[1]] != "e2" {
		t.Fatalf("Q2 parsed as %+v", q2)
	}
	if q3 := c.Elements[c.Index("q3")]; q3.Area != 0.5 || len(q3.Extra) != 1 {
		t.Fatalf("Q3 parsed as %+v", q3)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb, ".op"); err != nil {
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
		t.Fatalf("round trip differs:\n%s", sb.String())
	}
}

// MOSFETs with instance parameters survive a write and re-parse.
func TestWriteNetlistRoundTripMOS(t *testing.T) {
	src := "mos\n.model nm NMOS(LEVEL=1 VTO=0.7 KP=50u TOX=20n)\n.model pm PMOS(VTO=-0.7)\n" +
		"VDD vdd 0 5\nVIN in 0 1\nMN out in 0 0 nm L=1u W=10u AD=50p AS=50p PD=30u PS=30u NRD=0.5 NRS=0.5 M=2\n" +
		"MP out in vdd vdd pm L=1u W=25u\nRL out 0 1meg\n.op\n.end\n"
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	mn := c.Elements[c.Index("mn")]
	if mn.Kind != MOSFET || mn.Params["w"] != 10e-6 || mn.Params["m"] != 2 || c.NodeNames[mn.Extra[0]] != "in" || c.NodeNames[mn.Nodes[1]] != "0" {
		t.Fatalf("MN parsed as %+v", mn)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb, ".op"); err != nil {
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
		t.Fatalf("round trip differs:\n%s", sb.String())
	}
	for name, bad := range map[string]string{
		"unknown instance parameter": "t\nV1 a 0 1\nM1 a a 0 0 nm Q=1\n.model nm NMOS\n",
		"zero length":                "t\nV1 a 0 1\nM1 a a 0 0 nm L=0\n.model nm NMOS\n",
		"level 4":                    "t\nV1 a 0 1\nM1 a a 0 0 nm\n.model nm NMOS(LEVEL=4)\n",
		"bjt model":                  "t\nV1 a 0 1\nM1 a a 0 0 nm\n.model nm NPN\n",
	} {
		d, err := netlist.ParseString(bad)
		if err != nil {
			continue
		}
		if _, err := FromDeck(d); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

// Validation of the new kinds returns typed errors.
func TestValidateM7(t *testing.T) {
	cases := map[string]string{
		"F controlled by a resistor": "t\nV1 a 0 1\nR1 a 0 1k\nF1 0 b R1 2\nR2 b 0 1k\n",
		"E with equal control nodes": "t\nV1 a 0 1\nE1 b 0 a a 2\nR2 b 0 1k\n",
		"BJT with a diode model":     "t\nV1 a 0 1\nQ1 a a 0 d\n.model d D(IS=1e-14)\n",
		"BJT with excess phase":      "t\nV1 a 0 1\nQ1 a a 0 q\n.model q NPN(PTF=30)\n",
		"BJT with no model":          "t\nV1 a 0 1\nQ1 a a 0 q\n",
		"G output floating":          "t\nV1 a 0 1\nR1 a 0 1k\nG1 b c a 0 1m\n",
	}
	for name, src := range cases {
		d, err := netlist.ParseString(src)
		if err != nil {
			t.Fatalf("%s: parse: %v", name, err)
		}
		if _, err := FromDeck(d); err == nil {
			t.Errorf("%s: no error", name)
		}
	}
}

func TestCloneIsDeep(t *testing.T) {
	c, err := NewBuilder("x").VoltageSource("V1", "a", "0", 1).Resistor("R1", "a", "0", 1e3).
		Model("d", "d", map[string]float64{"is": 1e-14}).Build()
	if err != nil {
		t.Fatal(err)
	}
	k := c.Clone()
	k.Elements[c.Index("r1")].Value = 5
	k.Models["d"].Params["is"] = 1
	if c.Elements[1].Value != 1e3 || c.Models["d"].Params["is"] != 1e-14 {
		t.Fatal("clone shares state with the original")
	}
}

// A builder name without its SPICE kind letter is written with the letter
// prefixed, and references to it follow, so the netlist parses back.
func TestWriteNetlistKindLetter(t *testing.T) {
	d, err := netlist.ParseString("t\nV1 a 0 1\nR1 a b 1k\nF1 b 0 V1 2\nL1 b 0 1m\nL2 c 0 1m\nRc c 0 1\nK1 L1 L2 0.5\n.end\n")
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	ren := map[string]string{"V1": "U", "L1": "coil", "K1": "M"}
	for i := range c.Elements {
		e := &c.Elements[i]
		if n, ok := ren[strings.ToUpper(e.Name)]; ok {
			e.Name = n
		}
		if strings.EqualFold(e.CtrlSource, "V1") {
			e.CtrlSource = "U"
		}
		for k, l := range e.Couples {
			if strings.EqualFold(l, "L1") {
				e.Couples[k] = "coil"
			}
		}
	}
	var b strings.Builder
	if err := c.WriteNetlist(&b); err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"\nVU a 0", "F1 b 0 VU 2", "\nLcoil b 0", "KM Lcoil "} {
		if !strings.Contains(b.String(), want) {
			t.Errorf("missing %q in\n%s", want, b.String())
		}
	}
	if _, err := netlist.ParseString(b.String()); err != nil {
		t.Errorf("does not parse back: %v", err)
	}
}
