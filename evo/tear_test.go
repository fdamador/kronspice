// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/kron"
)

const tearSrc = `two-stage RC with a fixed load ladder
V1 in 0 AC 1
R1 in a 10k
C1 a 0 10n
RT a b 1k
R3 b c 1k
C3 c 0 1n
R4 c d 1k
C4 d 0 1n
RL d 0 100k
.ac dec 20 10 1meg
`

func tearProblem(t testing.TB, tear *Tear) Problem {
	spec := analysis.ACSpec{Variation: "dec", Points: 20, FStart: 10, FStop: 1e6}
	return Problem{
		Template: deck(t, tearSrc),
		Space: Space{
			{Element: "R1", Min: 100, Max: 1e6, Series: "E24"},
			{Element: "C1", Min: 1e-10, Max: 1e-5, Series: "E12"},
		},
		Objectives: []Objective{Cutoff("d", "0", spec, 2000, true)},
		Tear:       tear,
	}
}

// Tearing changes the time a run takes, not its result: the torn run's
// best individual equals the MNA run's, and the shared cache is used.
func TestTearSameResults(t *testing.T) {
	opt := Options{Population: 20, Generations: 8, Seed: 5}
	a, err := Run(context.Background(), tearProblem(t, nil), opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	cache := kron.NewTornCache(0)
	b, err := Run(context.Background(), tearProblem(t, &Tear{Torn: []string{"RT"}, Cache: cache}), opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	for i := range a.Best.Values {
		if a.Best.Values[i] != b.Best.Values[i] {
			t.Fatalf("best %v torn, %v MNA", b.Best.Values, a.Best.Values)
		}
	}
	if math.Abs(a.Best.Fitness-b.Best.Fitness) > 1e-9*math.Max(1, a.Best.Fitness) {
		t.Errorf("fitness %g torn, %g MNA", b.Best.Fitness, a.Best.Fitness)
	}
	if h, m, _ := cache.Stats(); h == 0 || m == 0 {
		t.Errorf("cache hits %d misses %d", h, m)
	}
}

func TestTearValidation(t *testing.T) {
	opt := Options{Population: 4, Generations: 1, Seed: 1}
	_, err := Run(context.Background(), tearProblem(t, &Tear{Torn: []string{"V1"}}), opt, nil, nil)
	if !errors.Is(err, ErrBadOptions) {
		t.Errorf("tearing a source: %v", err)
	}
	if _, err := Run(context.Background(), tearProblem(t, &Tear{Torn: []string{"R9"}}), opt, nil, nil); !errors.Is(err, ErrBadOptions) {
		t.Errorf("tearing a missing element: %v", err)
	}
}

// Candidates the torn path cannot take are solved by MNA: a cut that
// leaves a floating piece gives the same fitness as no tearing.
func TestTearFallback(t *testing.T) {
	opt := Options{Population: 6, Generations: 2, Seed: 2}
	// Tearing R1, RT, R3 and R4 leaves nodes a and c held to ground by a
	// capacitor alone, nearly floating at the low end of the sweep; any
	// sweep the torn path rejects is solved by MNA.
	a, err := Run(context.Background(), tearProblem(t, nil), opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	b, err := Run(context.Background(), tearProblem(t, &Tear{Torn: []string{"R1", "RT", "R3", "R4"}}), opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(a.Best.Fitness-b.Best.Fitness) > 1e-9*math.Max(1, a.Best.Fitness) {
		t.Errorf("fitness %g torn, %g MNA", b.Best.Fitness, a.Best.Fitness)
	}
}

// The GA precheck (kron.Precheck, M9e) leaves a sound problem's result
// unchanged and penalizes a candidate with a loop of voltage sources as
// singular before simulating it, naming the loop.
func TestPrecheckInGA(t *testing.T) {
	opt := Options{Population: 20, Generations: 4, Seed: 6}
	a, err := Run(context.Background(), tearProblem(t, nil), opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	p := tearProblem(t, nil)
	p.Precheck = &kron.PrecheckOptions{DC: true}
	b, err := Run(context.Background(), p, opt, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if a.Best.Fitness != b.Best.Fitness {
		t.Errorf("precheck changed the best fitness: %g → %g", a.Best.Fitness, b.Best.Fitness)
	}
	bad := tearProblem(t, nil)
	bad.Template = deck(t, tearSrc+"V2 in 0 AC 1\n")
	bad.Precheck = &kron.PrecheckOptions{}
	s, err := Run(context.Background(), bad, Options{Population: 4, Generations: 1, Seed: 1}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if s.Best.Failure != "singular" || !containsAll(s.Best.Error, "loop", "V1", "V2") {
		t.Errorf("V1 ∥ V2 candidate: failure %q, error %q", s.Best.Failure, s.Best.Error)
	}
}

func containsAll(s string, parts ...string) bool {
	for _, p := range parts {
		if !strings.Contains(s, p) {
			return false
		}
	}
	return true
}
