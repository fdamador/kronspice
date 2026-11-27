// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package netlist

import (
	"errors"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"testing/fstest"
)

func TestEvalExpr(t *testing.T) {
	cases := []struct {
		in   string
		want float64
	}{
		{"1+2*3", 7}, {"(1+2)*3", 9}, {"2^3^2", 512}, {"2**3", 8}, {"-2^2", -4},
		{"10k/2", 5000}, {"1meg*1u", 1}, {"4.7kohm", 4700}, {"sqrt(16)+abs(-1)", 5},
		{"max(1, 3) - min(2,5)", 1}, {"2*pi", 2 * math.Pi}, {"log10(1000)", 3},
		{"ln(exp(2))", 2}, {"1 < 2", 1}, {"1 >= 2", 0}, {"1 == 1 && 0 || 1", 1},
		{"!0", 1}, {"3 != 3", 0}, {"1 ? 10 : 20", 10}, {"0 ? 10 : 2>1 ? 30 : 40", 30},
		{"7 % 4", 3}, {"pwr(-2, 2)", -4}, {"nint(2.5)", 2}, {".5e1", 5}, {"sgn(-3)", -1},
	}
	for _, c := range cases {
		got, err := EvalExpr(c.in)
		if err != nil {
			t.Errorf("EvalExpr(%q): %v", c.in, err)
			continue
		}
		if math.Abs(got-c.want) > 1e-12*math.Max(1, math.Abs(c.want)) {
			t.Errorf("EvalExpr(%q) = %g, want %g", c.in, got, c.want)
		}
	}
	for _, bad := range []string{"", "1+", "(1", "1/0", "foo", "sqrt(1,2)", "nope(1)", "1 2",
		"ln(0)", "sqrt(-1)", strings.Repeat("(", 500) + "1" + strings.Repeat(")", 500)} {
		if _, err := EvalExpr(bad); !errors.Is(err, ErrParam) || !errors.Is(err, ErrSyntax) {
			t.Errorf("EvalExpr(%q) err = %v, want ErrParam", bad, err)
		}
	}
}

func TestParam(t *testing.T) {
	d, err := ParseString(`params
.param rtop = {2*rbot} rbot=1k
.param vin='5*gain', gain=2
V1 in 0 DC {vin}
R1 in out {rtop}
R2 out 0 'rbot'
C1 out 0 {1/(2*pi*rbot*1k)} IC={vin/2}
.tran {1u} {10*1u}
.end
`)
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]float64{"V1": 10, "R1": 2000, "R2": 1000, "C1": 1 / (2 * math.Pi * 1e6)}
	for _, e := range d.Elements {
		if math.Abs(e.Value-want[e.Name]) > 1e-12*want[e.Name] {
			t.Errorf("%s = %g, want %g", e.Name, e.Value, want[e.Name])
		}
	}
	if c := d.Elements[3]; !c.HasIC || c.IC != 5 {
		t.Errorf("C1 IC = %g", c.IC)
	}
	if a := d.Analyses[0]; a.Step != 1e-6 || math.Abs(a.Stop-1e-5) > 1e-20 {
		t.Errorf(".tran = %g %g", a.Step, a.Stop)
	}
	if d.Params["rtop"] != 2000 || d.Params["gain"] != 2 || len(d.Params) != 4 {
		t.Errorf("Params = %v", d.Params)
	}

	for _, bad := range []string{
		".param a={b} b={a}\nR1 1 0 {a}",
		".param a=1 a=2\nR1 1 0 {a}",
		"R1 1 0 {undefined}",
		".param x={1/0}\nR1 1 0 1",
		".param =3",
		"R1 1 0 {1",
		".param sqrt=2",
	} {
		if _, err := ParseString("t\n" + bad + "\n"); !errors.Is(err, ErrParam) {
			t.Errorf("%q: err = %v, want ErrParam", bad, err)
		}
	}
}

const subcktDeck = `subcircuits
.param rload=4k
.subckt divider top mid params: rt=1k rb={rt}
R1 top mid {rt}
R2 mid 0 {rb}
.ends divider
.subckt stage in out params: g=2
.subckt buf a b
E1 b 0 a 0 1
.ends
Vs in sense 0
X1 sense n divider rt={1k*g}
X2 n out buf
F1 0 out Vs {g}
D1 out 0 dloc
.model dloc D(IS=1e-14)
.ends
V1 vin 0 DC 10
XA vin a divider
XB vin b divider rt=2k rb={rload}
XS a s stage params: g=3
.op
.end
`

func TestSubckt(t *testing.T) {
	d, err := ParseString(subcktDeck)
	if err != nil {
		t.Fatal(err)
	}
	byName := map[string]ElementLine{}
	for _, e := range d.Elements {
		byName[e.Name] = e
	}
	check := func(name string, value float64, nodes ...string) ElementLine {
		t.Helper()
		e, ok := byName[name]
		if !ok {
			t.Fatalf("no element %s in %v", name, names(d))
		}
		if value != 0 && e.Value != value {
			t.Errorf("%s = %g, want %g", name, e.Value, value)
		}
		if strings.Join(e.Nodes, " ") != strings.Join(nodes, " ") {
			t.Errorf("%s nodes %v, want %v", name, e.Nodes, nodes)
		}
		return e
	}
	check("R.XA.R1", 1000, "vin", "a")
	check("R.XA.R2", 1000, "a", "0")
	check("R.XB.R1", 2000, "vin", "b")
	check("R.XB.R2", 4000, "b", "0")
	check("V.XS.Vs", 0, "a", "xs.sense")
	check("R.XS.X1.R1", 3000, "xs.sense", "xs.n")
	check("R.XS.X1.R2", 3000, "xs.n", "0") // rb defaults to rt
	e := check("E.XS.X2.E1", 1, "s", "0")
	if strings.Join(e.CtrlNodes, " ") != "xs.n 0" {
		t.Errorf("E1 control %v", e.CtrlNodes)
	}
	f := check("F.XS.F1", 3, "0", "s")
	if f.CtrlSource != "v.xs.vs" {
		t.Errorf("F1 controlled by %q", f.CtrlSource)
	}
	dl := check("D.XS.D1", 0, "s", "0")
	if dl.Model != "xs.dloc" || len(d.Models) != 1 || d.Models[0].Name != "xs.dloc" {
		t.Errorf("local model: element uses %q, models %v", dl.Model, d.Models)
	}
	if len(d.Elements) != 11 {
		t.Errorf("%d elements: %v", len(d.Elements), names(d))
	}
}

func names(d *Deck) []string {
	var out []string
	for _, e := range d.Elements {
		out = append(out, e.Name)
	}
	return out
}

func TestSubcktErrors(t *testing.T) {
	def := ".subckt r2 a b params: r=1\nR1 a b {r}\n.ends\n"
	for _, c := range []struct {
		src  string
		want error
	}{
		{".subckt r2 a b\nR1 a b 1\n", ErrSubckt},                    // no .ends
		{".ends\n", ErrSubckt},                                       // stray .ends
		{def + "X1 1 r2\n", ErrSubckt},                               // too few nodes
		{def + "X1 1 0 nosuch\n", ErrSubckt},                         // unknown
		{def + "X1 1 0 r2 q=3\n", ErrParam},                          // undeclared parameter
		{".subckt loop a\nX1 a loop\n.ends\nX1 1 loop\n", ErrSubckt}, // recursion
		{".subckt s a\n.tran 1 2\n.ends\nX1 1 s\n", ErrSubckt},       // analysis inside
		{".subckt s a a\n.ends\n", ErrSubckt},                        // repeated port
		{def + def, ErrSubckt},                                       // defined twice
		{".subckt s a\n.ends t\n", ErrSubckt},                        // wrong .ends name
	} {
		_, err := ParseString("t\n" + c.src)
		if !errors.Is(err, c.want) || !errors.Is(err, ErrSyntax) {
			t.Errorf("%q: err = %v, want %v", c.src, err, c.want)
		}
	}
}

func TestSubcktGlobalAndScope(t *testing.T) {
	d, err := ParseString(`scopes
.global vdd
.param r=10
.subckt cell a
.param r2={2*r}
R1 a vdd {r2}
R2 a mid {r}
.ends
X1 n cell
`)
	if err != nil {
		t.Fatal(err)
	}
	if e := d.Elements[0]; e.Value != 20 || e.Nodes[1] != "vdd" {
		t.Errorf("R1 = %g on %v", e.Value, e.Nodes)
	}
	if e := d.Elements[1]; e.Value != 10 || e.Nodes[1] != "x1.mid" {
		t.Errorf("R2 = %g on %v", e.Value, e.Nodes)
	}
}

func TestInclude(t *testing.T) {
	fs := fstest.MapFS{
		"lib/parts.inc":    {Data: []byte("* parts\n.include sub/more.inc\n.param rr=2k\nR2 out 0 {rr}\n.end\nR9 never 0 1\n")},
		"lib/sub/more.inc": {Data: []byte(".subckt pair a b\nR1 a b\n+ 1k\n.ends\n")},
		"loop.inc":         {Data: []byte(".include loop.inc\n")},
	}
	opt := Options{Open: func(n string) (io.ReadCloser, error) { return fs.Open(filepath.ToSlash(n)) }}
	d, err := ParseWith(strings.NewReader("inc\nV1 in 0 1\n.inc \"lib/parts.inc\"\nX1 in out pair\n.op\n"), opt)
	if err != nil {
		t.Fatal(err)
	}
	if got := strings.Join(names(d), " "); got != "V1 R2 R.X1.R1" {
		t.Errorf("elements %s", got)
	}
	if d.Elements[1].Value != 2000 || d.Elements[2].Value != 1000 {
		t.Errorf("values %g %g", d.Elements[1].Value, d.Elements[2].Value)
	}
	if _, err := ParseWith(strings.NewReader("t\n.include loop.inc\n"), opt); !errors.Is(err, ErrInclude) {
		t.Errorf("loop: %v", err)
	}
	if _, err := ParseWith(strings.NewReader("t\n.include missing.inc\n"), opt); !errors.Is(err, ErrInclude) {
		t.Errorf("missing: %v", err)
	}
	if _, err := ParseString("t\n.include x.inc\n"); !errors.Is(err, ErrInclude) {
		t.Errorf("no file system: %v", err)
	}
	// An error inside an included file names that file.
	bad := fstest.MapFS{"bad.inc": {Data: []byte("R1 a\n")}}
	_, err = ParseWith(strings.NewReader("t\n.include bad.inc\n"),
		Options{Open: func(n string) (io.ReadCloser, error) { return bad.Open(n) }})
	if err == nil || !strings.Contains(err.Error(), "bad.inc") {
		t.Errorf("error does not name the file: %v", err)
	}
}

func TestParseFile(t *testing.T) {
	dir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(dir, "models"), 0o755); err != nil {
		t.Fatal(err)
	}
	write := func(name, s string) {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(s), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("models/d.lib", ".model dd D(IS=1e-15)\n")
	write("top.cir", "top\n.include models/d.lib\nV1 a 0 1\nD1 a 0 dd\n.op\n.end\n")
	d, err := ParseFile(filepath.Join(dir, "top.cir"))
	if err != nil {
		t.Fatal(err)
	}
	if len(d.Models) != 1 || d.Models[0].Params["is"] != 1e-15 {
		t.Errorf("models %v", d.Models)
	}
}
