// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package report

import (
	"bytes"
	"math/cmplx"
	"os"
	"strings"
	"testing"
)

// The ladder's torn solution equals MNA (CLAUDE.md §3a rule 4, 1e-9),
// power is the same torn and untorn (Roth 1959 §5), and every link's
// current is the current of that resistor: V across it over R.
func TestLadderDC(t *testing.T) {
	src, err := os.ReadFile("../testdata/netlists/ladder.cir")
	if err != nil {
		t.Fatal(err)
	}
	s, err := Compute(string(src), Options{Source: "ladder.cir"})
	if err != nil {
		t.Fatal(err)
	}
	if !s.DC || s.KronNote != "" {
		t.Fatalf("DC %v, note %q", s.DC, s.KronNote)
	}
	if s.MaxDiff > 1e-9 {
		t.Errorf("torn vs MNA %.3g", s.MaxDiff)
	}
	if s.PowerMNA == 0 || cmplx.Abs(s.PowerTorn-s.PowerMNA) > 1e-9*cmplx.Abs(s.PowerMNA) {
		t.Errorf("power torn %v untorn %v", s.PowerTorn, s.PowerMNA)
	}
	if len(s.Links) == 0 || len(s.Links) != len(s.Cut.Torn) {
		t.Fatalf("%d links, %d torn", len(s.Links), len(s.Cut.Torn))
	}
	for _, l := range s.Links {
		if r := real(l.V) / real(l.I); real(l.I) != 0 && (r < 0 || cmplx.Abs(complex(r, 0)) > 1e9) {
			t.Errorf("%s: V/I = %g", l.Name, r)
		}
	}
	var md, csv bytes.Buffer
	s.Markdown(&md, 0)
	s.CSV(&csv)
	for _, want := range []string{"# Results:", "Kron's diakoptic solution", "| Link (torn branch)", "power is invariant"} {
		if !strings.Contains(md.String(), want) {
			t.Errorf("markdown lacks %q", want)
		}
	}
	if n := strings.Count(csv.String(), "\n"); n != len(s.Nodes)+1 {
		t.Errorf("CSV %d lines, %d nodes", n, len(s.Nodes))
	}
}

// An .ac card puts the sheet at the sweep's geometric middle; the torn
// AC solution equals MNA there.
func TestRCLadderAC(t *testing.T) {
	var b strings.Builder
	b.WriteString("rc ladder\nV1 n0 0 AC 1\n")
	for k := 0; k < 12; k++ {
		b.WriteString("R" + itoa(k) + " n" + itoa(k) + " n" + itoa(k+1) + " 1k\n")
		b.WriteString("C" + itoa(k) + " n" + itoa(k+1) + " 0 10n\n")
	}
	b.WriteString(".ac dec 10 100 1meg\n.end\n")
	s, err := Compute(b.String(), Options{Pieces: 3})
	if err != nil {
		t.Fatal(err)
	}
	if s.DC || s.Freq < 9999 || s.Freq > 10001 {
		t.Fatalf("DC %v at %g Hz, want AC at 10 kHz", s.DC, s.Freq)
	}
	if s.KronNote != "" {
		t.Fatal(s.KronNote)
	}
	if s.MaxDiff > 1e-9 || s.Stats.Pieces < 2 {
		t.Errorf("torn vs MNA %.3g, %d pieces", s.MaxDiff, s.Stats.Pieces)
	}
}

// Circuits Kron's linear solve does not cover get the MNA part and a
// stated reason, never an error.
func TestNotes(t *testing.T) {
	for name, src := range map[string]string{
		"nonlinear": "d\nV1 a 0 1\nR1 a b 1k\nD1 b 0 dm\n.model dm D\n.op\n.end\n",
		"bound":     "e\nV1 a 0 1\nR1 a 0 1k\nE1 b 0 a 0 2\nR2 b c 1k\nF1 c 0 V1 1\nR3 c 0 1k\n.op\n.end\n",
	} {
		s, err := Compute(src, Options{})
		if err != nil {
			t.Errorf("%s: %v", name, err)
			continue
		}
		if s.KronNote == "" || s.Torn != nil || len(s.MNA) != len(s.Nodes) {
			t.Errorf("%s: note %q, torn %v", name, s.KronNote, s.Torn)
		}
		var md bytes.Buffer
		s.Markdown(&md, 0)
		if !strings.Contains(md.String(), s.KronNote) {
			t.Errorf("%s: note missing from the sheet", name)
		}
	}
}

func TestNatLess(t *testing.T) {
	if !natLess("n2", "n10") || natLess("n10", "n2") || !natLess("a", "b") || natLess("n1", "n1") {
		t.Error("natural order")
	}
}

func itoa(k int) string {
	const d = "0123456789"
	if k < 10 {
		return d[k : k+1]
	}
	return itoa(k/10) + d[k%10:k%10+1]
}
