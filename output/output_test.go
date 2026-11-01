// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"bytes"
	"strings"
	"testing"
)

func TestEng(t *testing.T) {
	cases := map[float64]string{
		0: "0 V", 5: "5 V", 0.0047: "4.7 mV", -2.2e-6: "-2.2 µV", 1500: "1.5 kV", 3.3e6: "3.3 MV", 1e-15: "1 fV",
	}
	for x, want := range cases {
		if got := Eng(x, "V"); got != want {
			t.Errorf("Eng(%g) = %q, want %q", x, got, want)
		}
	}
}

func TestWriteCSVAndReport(t *testing.T) {
	r := &Result{
		Title: "t", Analysis: "op",
		NodeNames: []string{"0", "a"}, NodeVoltages: []float64{0, 2},
		Elements: []ElementResult{{Name: "R1", Kind: "resistor", Nodes: [2]string{"a", "0"}, Voltage: 2, Current: 1, PowerAbsorbed: 2}},
		Ledger:   &Ledger{},
	}
	var b bytes.Buffer
	if err := WriteCSV(&b, r); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(b.String(), "v(a),2\n") || !strings.Contains(b.String(), "p(r1),2\n") {
		t.Fatalf("csv:\n%s", b.String())
	}
	b.Reset()
	if err := WriteReport(&b, r, ReportOptions{}); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(b.String(), "v(a)") || !strings.Contains(b.String(), "residual") {
		t.Fatalf("report:\n%s", b.String())
	}
}
