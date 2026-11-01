// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"math"
	"math/rand/v2"
)

// rngFor returns the random stream for one stage of a run. Every stream is
// derived from the experiment seed and a stream number (0 = initial
// population, g = breeding of generation g), so a run resumed from a
// checkpoint draws exactly the numbers an uninterrupted run would
// (CLAUDE.md §4, determinism). The global math/rand source is never used.
func rngFor(seed, stream uint64) *rand.Rand {
	return rand.New(rand.NewPCG(seed, stream^0x9e3779b97f4a7c15))
}

// randomGenome draws magnitudes uniformly in log space and free signs with
// equal probability.
func (s Space) randomGenome(r *rand.Rand) Genome {
	g := Genome{X: make([]float64, len(s)), Sign: make([]int8, len(s))}
	for i, gene := range s {
		lo, hi := gene.logRange()
		g.X[i] = lo + r.Float64()*(hi-lo)
		g.Sign[i] = 1
		if r.IntN(2) == 0 {
			g.Sign[i] = -1
		}
	}
	s.repair(&g)
	return g
}

// tournament returns the index of the fittest of k individuals drawn with
// replacement; ties go to the lower index, so the result never depends on
// anything but the random stream.
func tournament(r *rand.Rand, fit []float64, k int) int {
	best := r.IntN(len(fit))
	for i := 1; i < k; i++ {
		c := r.IntN(len(fit))
		if fit[c] < fit[best] || (fit[c] == fit[best] && c < best) {
			best = c
		}
	}
	return best
}

// crossover is blend crossover (BLX-α, α = 0.5; Eshelman & Schaffer 1993)
// on the log magnitudes, and uniform crossover on the signs.
func crossover(r *rand.Rand, a, b Genome) (Genome, Genome) {
	const alpha = 0.5
	c1, c2 := a.Clone(), b.Clone()
	for i := range a.X {
		lo, hi := math.Min(a.X[i], b.X[i]), math.Max(a.X[i], b.X[i])
		d := hi - lo
		c1.X[i] = lo - alpha*d + r.Float64()*(1+2*alpha)*d
		c2.X[i] = lo - alpha*d + r.Float64()*(1+2*alpha)*d
		if r.IntN(2) == 0 {
			c1.Sign[i], c2.Sign[i] = b.Sign[i], a.Sign[i]
		}
	}
	return c1, c2
}

// mutate perturbs each gene with probability rate: a Gaussian step of sigma
// decades on the magnitude and, for free-sign genes, a sign flip with
// probability flip.
func mutate(r *rand.Rand, g *Genome, rate, sigma, flip float64) {
	for i := range g.X {
		if r.Float64() < rate {
			g.X[i] += sigma * r.NormFloat64()
		}
		if r.Float64() < flip {
			g.Sign[i] = -g.Sign[i]
		}
	}
}

// diversity is the mean over genes of the population's standard deviation
// of the log magnitude [decades].
func diversity(pop []Individual) float64 {
	if len(pop) < 2 {
		return 0
	}
	if len(pop[0].Genome.X) == 0 {
		// Topology genomes: spread of the netlist size [parts].
		mean := 0.0
		for _, p := range pop {
			mean += float64(len(p.Genome.Parts))
		}
		mean /= float64(len(pop))
		v := 0.0
		for _, p := range pop {
			d := float64(len(p.Genome.Parts)) - mean
			v += d * d
		}
		return math.Sqrt(v / float64(len(pop)))
	}
	n := len(pop[0].Genome.X)
	total := 0.0
	for i := 0; i < n; i++ {
		mean := 0.0
		for _, p := range pop {
			mean += p.Genome.X[i]
		}
		mean /= float64(len(pop))
		v := 0.0
		for _, p := range pop {
			d := p.Genome.X[i] - mean
			v += d * d
		}
		total += math.Sqrt(v / float64(len(pop)))
	}
	return total / float64(n)
}
