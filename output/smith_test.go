// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"bytes"
	"encoding/json"
	"encoding/xml"
	"io"
	"math"
	"math/cmplx"
	"strings"
	"testing"
)

func TestSmithGammaAndRegion(t *testing.T) {
	cases := []struct {
		z      complex128
		g      complex128
		region SmithRegion
	}{
		{50, 0, SmithPassive},
		{0, -1, SmithLossless},
		{complex(0, 50), complex(0, 1), SmithLossless},
		{-25, 3, SmithActive},  // r = −½: Γ = (−25−50)/(−25+50) = −3
		{-100, 3, SmithActive}, // r = −2: Γ = (−150)/(−50) = 3
	}
	cases[3].g = -3
	for _, c := range cases {
		if g := SmithGamma(c.z, 50); cmplx.Abs(g-c.g) > 1e-15 {
			t.Errorf("Γ(%v) = %v, want %v", c.z, g, c.g)
		}
		if r := ClassifySmith(c.z, 50, 0); r != c.region {
			t.Errorf("region(%v) = %v, want %v", c.z, r, c.region)
		}
	}
	if g := SmithGamma(-50, 50); !cmplx.IsInf(g) {
		t.Errorf("Γ(−Z0) = %v, want ∞", g)
	}
}

func activeSweep() *SmithResult {
	r := &SmithResult{Title: "t", Plus: "a", Minus: "0", Z0: 50}
	for i, z := range []complex128{complex(20, -40), complex(0, 10), complex(-20, 5), complex(-50, 0), complex(-20, -5), complex(20, 40)} {
		g := SmithGamma(z, 50)
		r.Points = append(r.Points, SmithPoint{Freq: float64(i + 1), Z: z, Gamma: g, Region: ClassifySmith(z, 50, 0), PortPower: real(z) / 2, ModelSupplied: math.Max(0, -real(z)/2)})
	}
	return r
}

func TestSmithSweepSummaries(t *testing.T) {
	r := activeSweep()
	if b := r.ActiveBands(); len(b) != 1 || b[0] != [2]float64{3, 5} {
		t.Fatalf("bands %v", b)
	}
	if m, inf := r.MaxGamma(); !inf || math.Abs(m-cmplx.Abs(SmithGamma(complex(-20, 5), 50))) > 1e-12 {
		t.Fatalf("max |Γ| %g inf %v", m, inf)
	}
	var b bytes.Buffer
	if err := WriteSmithReport(&b, r); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(b.String(), "Negative resistance") || !strings.Contains(b.String(), "Γ = ∞") {
		t.Fatalf("report:\n%s", b.String())
	}
	b.Reset()
	if err := WriteSmithCSV(&b, r); err != nil {
		t.Fatal(err)
	}
	if lines := strings.Count(b.String(), "\n"); lines != 7 || !strings.HasPrefix(b.String(), "freq,re_z,im_z,re_gamma") {
		t.Fatalf("csv:\n%s", b.String())
	}
	// JSON has no infinity: Γ at Z = −Z0 is null, and the document encodes.
	js, err := json.Marshal(SmithDocument(r))
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(js), `"re_gamma":[`) || !strings.Contains(string(js), "null") || !strings.Contains(string(js), `"active"`) {
		t.Fatalf("json %s", js)
	}
}

// The extended chart is well-formed SVG, covers |Γ| > 1, draws the
// negative-r circles and the r = −1 line, and dashes the active segments.
func TestWriteSmithSVG(t *testing.T) {
	r := activeSweep()
	var gs []complex128
	for _, p := range r.Points {
		gs = append(gs, p.Gamma)
	}
	var b bytes.Buffer
	if err := WriteSmithSVG(&b, []SmithTrace{{Name: "port <a>", Gamma: gs}}, SmithChartOptions{Title: "Extended & full"}); err != nil {
		t.Fatal(err)
	}
	s := b.String()
	d := xml.NewDecoder(strings.NewReader(s))
	for {
		_, err := d.Token()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("SVG is not well-formed: %v", err)
		}
	}
	for _, want := range []string{"r=−1", "r=-0.5", "r=-2", `stroke-dasharray="6 3"`, "port &lt;a&gt;", "Re Γ", "#fdecea"} {
		if !strings.Contains(s, want) {
			t.Errorf("SVG lacks %q", want)
		}
	}
	// Passive only: the default extent stays near the unit circle.
	b.Reset()
	if err := WriteSmithSVG(&b, []SmithTrace{{Name: "p", Gamma: []complex128{0, 0.5, complex(0, 0.9)}}}, SmithChartOptions{}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(b.String(), "r=0.5") {
		t.Errorf("passive chart should label r = 0.5")
	}
	if strings.Contains(b.String(), "r=-2") {
		t.Errorf("passive chart should not reach r = −2 (Γ = 3)")
	}
}
