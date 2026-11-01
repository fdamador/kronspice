// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"math"
	"testing"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/evo"
)

// A widely spaced pair behaves like a single microstrip; coupling and
// crosstalk fall as the gap grows; the GA finds a 50 Ω pair under −30 dB.
func TestCoupledPair(t *testing.T) {
	var prevK, prevN float64 = 1, 0
	for i, s := range []float64{0.5e-3, 2e-3, 6e-3} {
		p := Pair{W: 3e-3, S: s}
		r, err := p.Extract()
		if err != nil {
			t.Fatal(err)
		}
		k := r.L[0][1] / r.L[0][0]
		c, err := p.Circuit(r)
		if err != nil {
			t.Fatal(err)
		}
		n, err := NEXT(c, fNext, analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("s = %.1f mm: Z₀ %.2f Ω, K %.4f, NEXT %.1f dB", s*1e3, Z0(r), k, n)
		if k <= 0 || k >= prevK || (i > 0 && n >= prevN) {
			t.Errorf("coupling should fall with the gap: K %.4f (was %.4f), NEXT %.1f (was %.1f)", k, prevK, n, prevN)
		}
		prevK, prevN = k, n
	}
}

func TestEvolvePair(t *testing.T) {
	if testing.Short() {
		t.Skip("GA run")
	}
	// A generous per-evaluation timeout: the default 10 s is wall-clock
	// time, and under -race some field extractions exceed it, which
	// penalizes them and changes the evolved result (seen on 2026-10-27:
	// s = 2.29 mm, NEXT −30.0 dB under -race against s = 10 mm, −53 dB
	// without). With no timeouts the run is the same on every machine.
	s, err := evo.Run(context.Background(), Problem(), evo.Options{Population: 12, Generations: 6, Seed: 1, Timeout: 30 * time.Minute}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	for _, g := range s.History {
		if g.TimedOut > 0 {
			t.Fatalf("generation %d: %d candidates timed out, so the result depends on machine speed", g.Generation, g.TimedOut)
		}
	}
	p := Pair{W: s.Best.Values[0], S: s.Best.Values[1]}
	r, err := p.Extract()
	if err != nil {
		t.Fatal(err)
	}
	c, _ := p.Circuit(r)
	n, _ := NEXT(c, fNext, analysis.Options{})
	t.Logf("w %.3f mm, s %.3f mm: Z₀ %.2f Ω, NEXT %.1f dB", p.W*1e3, p.S*1e3, Z0(r), n)
	if math.Abs(Z0(r)-50) > 2.5 || n > -30 {
		t.Errorf("evolved pair misses the spec: Z₀ %.2f Ω, NEXT %.1f dB", Z0(r), n)
	}
}
