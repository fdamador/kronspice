// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command torn_ga evolves a two-pole RC pre-filter that drives a fixed
// 40-section LC line, three ways: AC solved by MNA, by Kron's node-pair
// form untorn, and by Kron's diakoptics with the line torn off at RT and
// its factors cached (milestone M9b). Every run uses the same seed; the
// program prints the time of each and checks that all three find the same
// best design.
//
//	go run ./examples/torn_ga [-plant line|grid] [-n 40] [-gens 30] [-pop 40]
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"strings"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
)

// Netlist is the problem circuit: source, evolved RC pre-filter (R1 C1 R2
// C2), the torn coupling resistor RT, and a fixed plant from p0 to the
// output node "out":
//
//   - "line": an n-section LC line (1 Ω, 1 µH, 1 nF per section,
//     √(L/C) ≈ 31.6 Ω) matched at its far end;
//   - "grid": an n×n RC mesh (10 Ω between neighbours, 1 nF to ground at
//     every node), driven at one corner and loaded by 100 Ω at the other —
//     the kind of plant (a resistive plane, a substrate) whose LU
//     factorization fills in.
func Netlist(plant string, n int) string {
	var b strings.Builder
	b.WriteString("torn GA: RC pre-filter into a fixed " + plant + "\n")
	b.WriteString("V1 in 0 AC 1\nR1 in m 1k\nC1 m 0 10n\nR2 m f 1k\nC2 f 0 10n\nRT f p0 10\n")
	switch plant {
	case "grid":
		node := func(i, j int) string {
			if i == 0 && j == 0 {
				return "p0"
			}
			if i == n-1 && j == n-1 {
				return "out"
			}
			return fmt.Sprintf("g%d_%d", i, j)
		}
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				fmt.Fprintf(&b, "CG%d_%d %s 0 1n\n", i, j, node(i, j))
				if i+1 < n {
					fmt.Fprintf(&b, "RV%d_%d %s %s 10\n", i, j, node(i, j), node(i+1, j))
				}
				if j+1 < n {
					fmt.Fprintf(&b, "RH%d_%d %s %s 10\n", i, j, node(i, j), node(i, j+1))
				}
			}
		}
		b.WriteString("RL out 0 100\n.end\n")
	default:
		for i := 0; i < n; i++ {
			to := fmt.Sprintf("p%d", i+1)
			if i == n-1 {
				to = "out"
			}
			fmt.Fprintf(&b, "RP%d p%d q%d 1\nLP%d q%d %s 1u\nCP%d %s 0 1n\n", i, i, i, i, i, to, i, to)
		}
		b.WriteString("RL out 0 31.6\n.end\n")
	}
	return b.String()
}

// Problem is the GA problem: the cutoff at the plant's output, 5 kHz.
func Problem(plant string, n int, tear *evo.Tear) (evo.Problem, error) {
	d, err := netlist.ParseString(Netlist(plant, n))
	if err != nil {
		return evo.Problem{}, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return evo.Problem{}, err
	}
	spec := analysis.ACSpec{Variation: "dec", Points: 50, FStart: 100, FStop: 1e6}
	return evo.Problem{
		Template: c,
		Space: evo.Space{
			{Element: "R1", Min: 10, Max: 1e5, Series: "E24"},
			{Element: "C1", Min: 1e-10, Max: 1e-5, Series: "E12"},
			{Element: "R2", Min: 10, Max: 1e5, Series: "E24"},
			{Element: "C2", Min: 1e-10, Max: 1e-5, Series: "E12"},
		},
		Objectives: []evo.Objective{evo.Cutoff("out", "0", spec, 5000, true)},
		Tear:       tear,
	}, nil
}

// Result is one timed run.
type Result struct {
	Name    string
	Time    time.Duration
	Best    evo.Individual
	Hits    int64
	Misses  int64
	Entries int
}

// Run runs the three variants.
func Run(plant string, n, pop, gens int, seed uint64) ([]Result, error) {
	opt := evo.Options{Population: pop, Generations: gens, Seed: seed}
	var out []Result
	for _, v := range []struct {
		name string
		tear *evo.Tear
	}{
		{"MNA", nil},
		{"Kron node form, untorn", &evo.Tear{Torn: nil, Cache: kron.NewTornCache(1)}},
		{"Kron diakoptics, torn at RT", &evo.Tear{Torn: []string{"RT"}, Cache: kron.NewTornCache(0)}},
	} {
		p, err := Problem(plant, n, v.tear)
		if err != nil {
			return nil, err
		}
		t0 := time.Now()
		st, err := evo.Run(context.Background(), p, opt, nil, nil)
		if err != nil {
			return nil, err
		}
		r := Result{Name: v.name, Time: time.Since(t0), Best: st.Best}
		if v.tear != nil {
			r.Hits, r.Misses, r.Entries = v.tear.Cache.Stats()
		}
		out = append(out, r)
	}
	return out, nil
}

func main() {
	plant := flag.String("plant", "line", `fixed plant: "line" or "grid"`)
	size := flag.Int("n", 40, "line sections, or grid side")
	pop := flag.Int("pop", 40, "population")
	gens := flag.Int("gens", 30, "generations")
	seed := flag.Uint64("seed", 1, "seed")
	flag.Parse()
	res, err := Run(*plant, *size, *pop, *gens, *seed)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("plant %s n=%d, population %d, %d generations, seed %d\n\n", *plant, *size, *pop, *gens, *seed)
	for _, r := range res {
		fmt.Printf("%-30s %8.2f s   best fitness %.6g  R1 %s C1 %s R2 %s C2 %s",
			r.Name, r.Time.Seconds(), r.Best.Fitness, eng(r.Best.Values[0]), eng(r.Best.Values[1]), eng(r.Best.Values[2]), eng(r.Best.Values[3]))
		if r.Hits+r.Misses > 0 {
			fmt.Printf("   pieces cached %d / factored %d", r.Hits, r.Misses)
		}
		fmt.Println()
	}
	same := true
	for _, r := range res[1:] {
		for i, v := range r.Best.Values {
			same = same && v == res[0].Best.Values[i]
		}
		same = same && math.Abs(r.Best.Fitness-res[0].Best.Fitness) <= 1e-9*math.Max(1, res[0].Best.Fitness)
	}
	fmt.Printf("\nsame best design in all three: %v\n", same)
	fmt.Printf("torn vs untorn Kron: %.1f× faster; torn vs MNA: %.2f×\n",
		res[1].Time.Seconds()/res[2].Time.Seconds(), res[0].Time.Seconds()/res[2].Time.Seconds())
	if !same {
		os.Exit(1)
	}
}

func eng(v float64) string { return fmt.Sprintf("%.3g", v) }
