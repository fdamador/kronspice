// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"fmt"
	"math"
	"math/rand/v2"
	"sort"
	"strings"
)

// Algorithm selects how a population is ranked.
type Algorithm uint8

const (
	// WeightedSum ranks candidates by Σ weight·error (the default).
	WeightedSum Algorithm = iota
	// NSGA2 keeps every objective separate and evolves a Pareto front:
	// fast non-dominated sorting, crowding distance and the crowded
	// tournament of K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, "A fast and
	// elitist multiobjective genetic algorithm: NSGA-II," IEEE Trans. Evol.
	// Comput. 6(2), 2002 (§III-A sort, §III-B crowding, §III-C main loop),
	// with the constrained dominance of §VI: a failed candidate (singular
	// matrix, no convergence …) or one that violates a constraint
	// (Objective.Constraint) is dominated by every feasible candidate;
	// constraint violators compare by their total weighted violation, and
	// rank ahead of failed candidates, which compare by their penalty.
	NSGA2
)

func (a Algorithm) String() string {
	if a == NSGA2 {
		return "nsga2"
	}
	return "weighted"
}

// ParseAlgorithm reads "weighted" (or "") and "nsga2" (or "nsga-ii").
func ParseAlgorithm(s string) (Algorithm, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "", "weighted", "weighted_sum", "sum":
		return WeightedSum, nil
	case "nsga2", "nsga-ii", "nsgaii":
		return NSGA2, nil
	}
	return 0, fmt.Errorf("%w: algorithm %q (want weighted or nsga2)", ErrBadOptions, s)
}

// boundaryCrowding stands for an infinite crowding distance (the extreme
// points of a front); it is finite so checkpoints stay valid JSON.
const boundaryCrowding = math.MaxFloat64

// FastNonDominatedSort splits points (all objectives minimized) into
// Pareto fronts: front 0 is non-dominated, front 1 is non-dominated once
// front 0 is removed, and so on (Deb et al. 2002, §III-A). Indices within a
// front are in ascending order, so the result depends only on the input.
// Equal points share a front.
func FastNonDominatedSort(points [][]float64) [][]int {
	n := len(points)
	dominates := func(p, q []float64) bool {
		lt := false
		for k := range p {
			if p[k] > q[k] {
				return false
			}
			if p[k] < q[k] {
				lt = true
			}
		}
		return lt
	}
	count := make([]int, n)  // n_p: how many points dominate p
	sets := make([][]int, n) // S_p: the points p dominates
	var front []int
	for p := 0; p < n; p++ {
		for q := p + 1; q < n; q++ {
			switch {
			case dominates(points[p], points[q]):
				sets[p] = append(sets[p], q)
				count[q]++
			case dominates(points[q], points[p]):
				sets[q] = append(sets[q], p)
				count[p]++
			}
		}
	}
	for p := 0; p < n; p++ {
		if count[p] == 0 {
			front = append(front, p)
		}
	}
	var fronts [][]int
	for len(front) > 0 {
		fronts = append(fronts, front)
		var next []int
		for _, p := range front {
			for _, q := range sets[p] {
				if count[q]--; count[q] == 0 {
					next = append(next, q)
				}
			}
		}
		sort.Ints(next)
		front = next
	}
	return fronts
}

// CrowdingDistance returns the crowding distance of each member of one
// front (Deb et al. 2002, §III-B): for every objective, the normalized gap
// between a point's two neighbours, summed. The two extreme points of each
// objective get +Inf. Ties in an objective keep index order, so the
// distances are deterministic.
func CrowdingDistance(points [][]float64, front []int) []float64 {
	d := make([]float64, len(front))
	if len(front) == 0 {
		return d
	}
	if len(front) <= 2 {
		for i := range d {
			d[i] = math.Inf(1)
		}
		return d
	}
	order := make([]int, len(front)) // positions into front
	for k := range points[front[0]] {
		for i := range order {
			order[i] = i
		}
		sort.SliceStable(order, func(a, b int) bool {
			return points[front[order[a]]][k] < points[front[order[b]]][k]
		})
		lo, hi := points[front[order[0]]][k], points[front[order[len(order)-1]]][k]
		d[order[0]], d[order[len(order)-1]] = math.Inf(1), math.Inf(1)
		if !(hi > lo) {
			continue
		}
		for i := 1; i < len(order)-1; i++ {
			d[order[i]] += (points[front[order[i+1]]][k] - points[front[order[i-1]]][k]) / (hi - lo)
		}
	}
	return d
}

// objectivePoint is the vector NSGA-II minimizes for one candidate:
// weight·error for every non-constraint objective with a positive weight
// (a zero weight leaves a term reported but unranked).
func objectivePoint(ind Individual, obs []Objective) []float64 {
	out := make([]float64, 0, len(obs))
	for k, ob := range obs {
		if ob.Weight > 0 && !ob.Constraint {
			out = append(out, ob.Weight*ind.Terms[k])
		}
	}
	return out
}

// violation is a simulated candidate's total weighted constraint violation.
func violation(ind Individual, obs []Objective) float64 {
	v := 0.0
	for k, ob := range obs {
		if ob.Constraint {
			v += ob.Weight * ind.Terms[k]
		}
	}
	return v
}

// rankPopulation sets Rank, Crowding and Violation on every individual:
// feasible candidates by fast non-dominated sorting and crowding; then
// constraint violators, one front per violation value (smaller first);
// then failed candidates, one front per penalty value. Infeasible
// candidates get zero crowding.
func rankPopulation(pop []Individual, obs []Objective) {
	var feas, infeas []int
	for i := range pop {
		ind := &pop[i]
		ind.Violation = 0
		if ind.Failure == "" && len(ind.Terms) == len(obs) {
			if ind.Violation = violation(*ind, obs); ind.Violation == 0 {
				feas = append(feas, i)
				continue
			}
		}
		infeas = append(infeas, i)
	}
	pts := make([][]float64, len(feas))
	for j, i := range feas {
		pts[j] = objectivePoint(pop[i], obs)
	}
	fronts := FastNonDominatedSort(pts)
	for r, f := range fronts {
		cd := CrowdingDistance(pts, f)
		for j, m := range f {
			ind := &pop[feas[m]]
			ind.Rank = r
			ind.Crowding = math.Min(cd[j], boundaryCrowding)
		}
	}
	// key orders infeasible candidates: violators (class 0) by violation,
	// then failures (class 1) by penalty.
	key := func(i int) (int, float64) {
		if pop[i].Failure == "" && len(pop[i].Terms) == len(obs) {
			return 0, pop[i].Violation
		}
		return 1, pop[i].Fitness
	}
	sort.SliceStable(infeas, func(a, b int) bool {
		ca, va := key(infeas[a])
		cb, vb := key(infeas[b])
		if ca != cb {
			return ca < cb
		}
		return va < vb
	})
	r := len(fronts) - 1
	for j, i := range infeas {
		if j == 0 {
			r++
		} else {
			c0, v0 := key(infeas[j-1])
			c1, v1 := key(i)
			if c0 != c1 || v0 != v1 {
				r++
			}
		}
		pop[i].Rank = r
		pop[i].Crowding = 0
	}
}

// crowdedLess is Deb's crowded-comparison operator ≺ₙ: lower rank wins,
// then larger crowding distance, then (for determinism) lower index.
func crowdedLess(pop []Individual, a, b int) bool {
	if pop[a].Rank != pop[b].Rank {
		return pop[a].Rank < pop[b].Rank
	}
	if pop[a].Crowding != pop[b].Crowding {
		return pop[a].Crowding > pop[b].Crowding
	}
	return a < b
}

// crowdedTournament returns the winner of k individuals drawn with
// replacement under the crowded-comparison operator.
func crowdedTournament(r *rand.Rand, pop []Individual, k int) int {
	best := r.IntN(len(pop))
	for i := 1; i < k; i++ {
		if c := r.IntN(len(pop)); crowdedLess(pop, c, best) {
			best = c
		}
	}
	return best
}

// nsgaSurvivors ranks the combined parent + offspring population and keeps
// the n best under ≺ₙ (Deb et al. 2002, §III-C): whole fronts while they
// fit, then the most widely spread members of the front that does not.
// The survivors are returned in ≺ₙ order with their combined-population
// rank and crowding.
func nsgaSurvivors(all []Individual, obs []Objective, n int) []Individual {
	rankPopulation(all, obs)
	idx := make([]int, len(all))
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool { return crowdedLess(all, idx[a], idx[b]) })
	out := make([]Individual, n)
	for i := range out {
		out[i] = all[idx[i]]
	}
	return out
}

// ParetoFront returns the feasible rank-0 candidates of a ranked
// population (the current non-dominated set), in population order, each
// distinct design once (candidates with equal decoded values and errors
// are the same design).
func ParetoFront(pop []Individual) []Individual {
	var out []Individual
	seen := map[string]bool{}
	for _, ind := range pop {
		if ind.Rank != 0 || ind.Failure != "" || ind.Violation != 0 {
			continue
		}
		key := fmt.Sprint(ind.Values, ind.Terms, ind.Genome.Parts)
		if seen[key] {
			continue
		}
		seen[key] = true
		out = append(out, ind)
	}
	return out
}
