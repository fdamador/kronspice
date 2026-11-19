// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"fmt"
	"math"
	"sort"
)

// LocalSearch is a multi-part local search applied to the best individual
// (a memetic step; doc/M7.md, doc/M8.md). A GA on E-series values stalls
// where every single-part change makes the design worse but a coordinated
// change of two or three parts improves it — a bias pair that must move
// together, an L and a C that keep a resonance while the impedance moves.
// The search tries every set of up to Size genes, each moved one step on
// its own grid in either direction, takes the best neighbour if it
// improves the fitness, and repeats until none does or MaxEvals is spent.
//
// A step is to the next value of the gene's E-series, ±1 for an integer
// gene, and Step decades for a gene without a series. Topology genomes
// (direct encoding) move their part values the same way; the topology
// itself is left to the GA. Developmental genomes are not searched (their
// values come from arithmetic subtrees, not one gene per part).
//
// The improved individual replaces the worst one of the population
// (Lamarckian). Neighbours are evaluated on the worker pool and stored by
// index, so the search is deterministic and a checkpointed run resumes
// exactly. Weighted-sum runs only.
type LocalSearch struct {
	// Size is the largest number of genes moved together (default 2,
	// at most 4). The neighbourhood of n genes has Σ_{k≤Size} C(n,k)·2^k
	// members.
	Size int
	// Every runs the search after every Every-th generation; 0 runs it
	// once, after the last generation.
	Every int
	// MaxEvals bounds the evaluations of one search (default 5000).
	MaxEvals int
	// Step is the move in decades for genes without a series (default 0.02).
	Step float64
}

func (ls *LocalSearch) defaults() error {
	if ls.Size == 0 {
		ls.Size = 2
	}
	if ls.MaxEvals == 0 {
		ls.MaxEvals = 5000
	}
	if ls.Step == 0 {
		ls.Step = 0.02
	}
	if ls.Size < 1 || ls.Size > 4 || ls.Every < 0 || ls.MaxEvals < 1 || !(ls.Step > 0) || math.IsInf(ls.Step, 0) {
		return fmt.Errorf("%w: local search size %d (1–4), every %d, max_evals %d, step %g",
			ErrBadOptions, ls.Size, ls.Every, ls.MaxEvals, ls.Step)
	}
	return nil
}

// due reports whether the search runs after generation gen.
func (ls *LocalSearch) due(gen, last int) bool {
	if ls.Every > 0 {
		return gen%ls.Every == 0
	}
	return gen == last
}

// seriesStep returns the value one step from v along ser (dir ±1) within
// [min, max], or false at the end of the range. v is taken as already
// snapped; the step goes to the nearest series value strictly beyond it.
func seriesStep(v float64, ser []float64, min, max float64, dir int) (float64, bool) {
	lv := math.Log10(v)
	dec := math.Floor(lv)
	var cand []float64
	for d := dec - 2; d <= dec+2; d++ {
		p := math.Pow(10, d)
		for _, m := range ser {
			c := roundSig(m*p, 3)
			if c >= min*(1-1e-12) && c <= max*(1+1e-12) {
				cand = append(cand, c)
			}
		}
	}
	sort.Float64s(cand)
	const tol = 1e-9 // relative: the same series value
	if dir > 0 {
		for _, c := range cand {
			if c > v*(1+tol) {
				return c, true
			}
		}
		return 0, false
	}
	for i := len(cand) - 1; i >= 0; i-- {
		if cand[i] < v*(1-tol) {
			return cand[i], true
		}
	}
	return 0, false
}

// gridStep moves one log10 gene value x one step: along the series when
// there is one, by ±1 for integers, else by step decades, clamped to
// [min, max]. It returns false when the move leaves the value unchanged.
func gridStep(x float64, min, max float64, serName string, integer bool, step float64, dir int) (float64, bool) {
	ser, _ := series(serName)
	switch {
	case ser != nil:
		v := snap(math.Pow(10, x), ser, min, max)
		n, ok := seriesStep(v, ser, min, max, dir)
		if !ok {
			return x, false
		}
		return math.Log10(n), true
	case integer:
		v := math.Min(math.Max(math.Round(math.Pow(10, x)), math.Ceil(min)), math.Floor(max))
		n := v + float64(dir)
		if n < math.Ceil(min) || n > math.Floor(max) || n <= 0 {
			return x, false
		}
		return math.Log10(n), true
	}
	n := math.Min(math.Max(x+float64(dir)*step, math.Log10(min)), math.Log10(max))
	return n, n != x
}

// lsGenes returns, for genome g, the number of searchable genes and a
// function moving gene i one step (nil when the genome is not searched).
func (p Problem) lsGenes(g Genome, step float64) (int, func(g *Genome, i, dir int) bool) {
	if t := p.Topology; t != nil {
		if t.Encoding == Developmental {
			return 0, nil
		}
		return len(g.Parts), func(g *Genome, i, dir int) bool {
			pt := &g.Parts[i]
			pk, ok := t.kind(pt.Kind)
			if !ok {
				return false
			}
			x, moved := gridStep(pt.X, pk.Min, pk.Max, pk.Series, false, step, dir)
			pt.X = x
			return moved
		}
	}
	return len(p.Space), func(g *Genome, i, dir int) bool {
		gn := p.Space[i]
		x, moved := gridStep(g.X[i], gn.Min, gn.Max, gn.Series, gn.Integer, step, dir)
		g.X[i] = x
		return moved
	}
}

// neighbours lists every genome one multi-gene step from g, in a fixed
// order: subsets of size 1..size in lexicographic order, each with its
// direction patterns. Moves that change nothing are skipped.
func neighbours(g Genome, n, size int, move func(*Genome, int, int) bool) []Genome {
	var out []Genome
	idx := make([]int, 0, size)
	var subsets func(start, k int)
	subsets = func(start, k int) {
		if len(idx) == k {
			for mask := 0; mask < 1<<k; mask++ {
				h := g.Clone()
				ok := true
				for b, i := range idx {
					dir := 1
					if mask&(1<<b) != 0 {
						dir = -1
					}
					if !move(&h, i, dir) {
						ok = false
						break
					}
				}
				if ok {
					out = append(out, h)
				}
			}
			return
		}
		for i := start; i < n; i++ {
			idx = append(idx, i)
			subsets(i+1, k)
			idx = idx[:len(idx)-1]
		}
	}
	for k := 1; k <= size && k <= n; k++ {
		subsets(0, k)
	}
	return out
}

// improve runs the local search from ind. It returns the best individual
// found (ind itself if nothing improves) and the evaluations spent.
func (ev *evaluator) improve(ctx context.Context, ind Individual, ls LocalSearch) (Individual, int) {
	best, evals := ind, 0
	if best.Failure != "" {
		return best, 0
	}
	for evals < ls.MaxEvals {
		if ctx.Err() != nil {
			break
		}
		n, move := ev.p.lsGenes(best.Genome, ls.Step)
		if move == nil || n == 0 {
			break
		}
		nb := neighbours(best.Genome, n, ls.Size, move)
		if room := ls.MaxEvals - evals; len(nb) > room {
			nb = nb[:room]
		}
		if len(nb) == 0 {
			break
		}
		scored := ev.evaluate(ctx, nb)
		evals += len(nb)
		k := -1
		for i, s := range scored {
			if s.Failure == "" && s.Fitness < best.Fitness && (k < 0 || s.Fitness < scored[k].Fitness) {
				k = i
			}
		}
		if k < 0 {
			break
		}
		best = scored[k]
	}
	return best, evals
}

// localSearch improves the population's best individual and puts the
// result in place of the worst (ties: the later index). It returns the
// evaluations spent.
func (ev *evaluator) localSearch(ctx context.Context, pop []Individual, ls LocalSearch) int {
	bi, wi := 0, 0
	for i, ind := range pop {
		if ind.Fitness < pop[bi].Fitness {
			bi = i
		}
		if ind.Fitness >= pop[wi].Fitness {
			wi = i
		}
	}
	got, evals := ev.improve(ctx, pop[bi], ls)
	if got.Fitness < pop[bi].Fitness {
		pop[wi] = got
	}
	return evals
}
