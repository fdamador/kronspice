// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"testing"

	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/kron"
)

// The M9b acceptance check, small: the same seed gives the same best
// design by MNA, untorn Kron and torn Kron, and the torn run is faster
// than the untorn one.
func TestSameBestFasterTorn(t *testing.T) {
	if testing.Short() {
		t.Skip("GA run")
	}
	res, err := Run("line", 20, 16, 6, 3)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range res[1:] {
		for i, v := range r.Best.Values {
			if v != res[0].Best.Values[i] {
				t.Fatalf("%s: best %v, MNA best %v", r.Name, r.Best.Values, res[0].Best.Values)
			}
		}
	}
	if res[2].Time >= res[1].Time {
		t.Errorf("torn %v not faster than untorn %v", res[2].Time, res[1].Time)
	}
	for _, r := range res {
		t.Logf("%-30s %v", r.Name, r.Time)
	}
}

func benchEval(b *testing.B, plant string, n int, tear *evo.Tear) {
	p, err := Problem(plant, n, tear)
	if err != nil {
		b.Fatal(err)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		// A new seed each time: new candidates, as in a running GA (the
		// fixed line's pieces stay cached, the pre-filter's do not).
		opt := evo.Options{Population: 8, Generations: 1, Seed: uint64(i + 1), Workers: 1}
		if _, err := evo.Run(context.Background(), p, opt, nil, nil); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkGeneration(b *testing.B) {
	for _, pl := range []struct {
		name string
		n    int
	}{{"line", 40}, {"grid", 12}} {
		b.Run(pl.name+"/MNA", func(b *testing.B) { benchEval(b, pl.name, pl.n, nil) })
		b.Run(pl.name+"/torn", func(b *testing.B) {
			benchEval(b, pl.name, pl.n, &evo.Tear{Torn: []string{"RT"}, Cache: kron.NewTornCache(0)})
		})
	}
}
