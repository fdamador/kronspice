// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"fmt"
	"math"
	"sort"

	"github.com/fdamador/kronspice/circuit"
)

// MAP-Elites: J.-B. Mouret, J. Clune, "Illuminating search spaces by
// mapping elites," arXiv:1504.04909 (2015), Alg. 1. The search space is
// mapped onto a grid of behaviour descriptors — any number of dimensions,
// each a measured property of a candidate (its part count, how many
// transistors it holds, one objective's error, …) cut into bins — and the
// archive keeps the best candidate (lowest weighted fitness) found in each
// cell. Parents are drawn uniformly from the occupied cells, so a design
// that is poor on the whole but the only one of its kind (the first
// circuit with a biased transistor, the first that cuts the input–output
// wire) survives as a stepping stone instead of being out-selected, which
// is what an open-ended topology search needs (doc/M8.md).
//
// The archive is the population: State.Population holds the elites in
// cell order, so checkpoints and resume work as for the other algorithms.
// Candidates that fail (a penalty class) are never placed.

// Descriptor is one dimension of the MAP-Elites grid.
type Descriptor struct {
	Name     string
	Min, Max float64 // values outside are placed in the end bins
	Bins     int
	// Of measures an evaluated candidate; ok = false leaves it out of the
	// archive. It must be deterministic and safe for concurrent use.
	Of func(p Problem, ind Individual) (v float64, ok bool)
}

// MapElites configures the MAPElites algorithm.
type MapElites struct {
	Descriptors []Descriptor
}

// DescParts is the number of evolved parts after pruning (topology genomes)
// or genes (value genomes).
func DescParts(min, max float64, bins int) Descriptor {
	return Descriptor{Name: "parts", Min: min, Max: max, Bins: bins,
		Of: func(_ Problem, ind Individual) (float64, bool) { return float64(len(ind.Values)), true }}
}

// DescGenomeSize is the genome's structural size: parts, or program ops
// for developmental genomes.
func DescGenomeSize(min, max float64, bins int) Descriptor {
	return Descriptor{Name: "genome size", Min: min, Max: max, Bins: bins,
		Of: func(_ Problem, ind Individual) (float64, bool) { return float64(ind.Genome.size()), true }}
}

// DescMeshes is the number of independent closed paths (meshes) of the
// whole circuit, m = b − N + S (Kron 1939 p. 75, eqs. 3.1–3.2): the first
// Betti number of its graph. A transistor counts as two branches (from
// the base or gate to the emitter or source, and the main path), a
// MOSFET's bulk tie as none. With MAP-Elites it spreads the archive over
// topologies of increasing loop structure.
func DescMeshes(min, max float64, bins int) Descriptor {
	return Descriptor{Name: "meshes", Min: min, Max: max, Bins: bins,
		Of: func(p Problem, ind Individual) (float64, bool) {
			c, _, err := p.Circuit(ind.Genome)
			if err != nil {
				return 0, false
			}
			return float64(Meshes(c)), true
		}}
}

// Meshes returns m = b − N + S for c's graph.
func Meshes(c *circuit.Circuit) int {
	parent := make([]int, len(c.NodeNames))
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	b, merges := 0, 0
	branch := func(a, z int) {
		b++
		if ra, rz := find(a), find(z); ra != rz {
			parent[ra] = rz
			merges++
		}
	}
	for _, e := range c.Elements {
		branch(e.Nodes[0], e.Nodes[1])
		if len(e.Extra) > 0 && (e.Kind == circuit.BJT || e.Kind == circuit.MOSFET) {
			branch(e.Extra[0], e.Nodes[1])
		}
	}
	// N − S is the number of merges, the edges of a spanning forest.
	return b - merges
}

// DescCount is the number of elements of kind k the candidate adds to its
// template (for example transistors).
func DescCount(k circuit.Kind, min, max float64, bins int) Descriptor {
	return Descriptor{Name: fmt.Sprintf("%v count", k), Min: min, Max: max, Bins: bins,
		Of: func(p Problem, ind Individual) (float64, bool) {
			c, _, err := p.Circuit(ind.Genome)
			if err != nil {
				return 0, false
			}
			n := 0
			for _, e := range c.Elements {
				if e.Kind == k {
					n++
				}
			}
			if p.Template != nil {
				for _, e := range p.Template.Elements {
					if e.Kind == k {
						n--
					}
				}
			}
			return float64(n), true
		}}
}

// DescTerm is objective k's error (Individual.Terms[k]).
func DescTerm(name string, k int, min, max float64, bins int) Descriptor {
	return Descriptor{Name: name, Min: min, Max: max, Bins: bins,
		Of: func(_ Problem, ind Individual) (float64, bool) {
			if k < 0 || k >= len(ind.Terms) {
				return 0, false
			}
			return ind.Terms[k], true
		}}
}

func (m *MapElites) validate(nObjectives int) error {
	if m == nil || len(m.Descriptors) == 0 {
		return fmt.Errorf("%w: MAP-Elites needs at least one descriptor", ErrBadOptions)
	}
	cells := 1
	for _, d := range m.Descriptors {
		if d.Of == nil || d.Bins < 1 || d.Bins > 1000 || !(d.Max > d.Min) || math.IsInf(d.Max, 0) || math.IsInf(d.Min, 0) {
			return fmt.Errorf("%w: descriptor %q: bins %d, range [%g, %g]", ErrBadOptions, d.Name, d.Bins, d.Min, d.Max)
		}
		cells *= d.Bins
		if cells > 10_000_000 {
			return fmt.Errorf("%w: MAP-Elites grid of more than 10⁷ cells", ErrBadOptions)
		}
	}
	return nil
}

// cell returns the flat cell index of an individual, or false when it is
// not placed (failed, or a descriptor declines it).
func (m *MapElites) cell(p Problem, ind Individual) (int, bool) {
	if ind.Failure != "" {
		return 0, false
	}
	idx, stride := 0, 1
	for _, d := range m.Descriptors {
		v, ok := d.Of(p, ind)
		if !ok || math.IsNaN(v) {
			return 0, false
		}
		b := int(math.Floor((v - d.Min) / (d.Max - d.Min) * float64(d.Bins)))
		b = min(max(b, 0), d.Bins-1)
		idx += b * stride
		stride *= d.Bins
	}
	return idx, true
}

// insert places candidates into the archive (cell → elite): a candidate
// takes its cell when the cell is empty or it is strictly better.
// Candidates are inserted in slice order, so ties keep the earlier one.
func (m *MapElites) insert(p Problem, archive map[int]Individual, cands []Individual) {
	for _, ind := range cands {
		k, ok := m.cell(p, ind)
		if !ok {
			continue
		}
		if old, has := archive[k]; !has || ind.Fitness < old.Fitness {
			archive[k] = ind
		}
	}
}

// elites lists the archive in cell order.
func archiveList(archive map[int]Individual) []Individual {
	keys := make([]int, 0, len(archive))
	for k := range archive {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	out := make([]Individual, len(keys))
	for i, k := range keys {
		out[i] = archive[k]
	}
	return out
}

// mapStep inserts a batch into the archive held by pop (the elites of the
// previous step, nil at the start) and returns the new archive in cell
// order. While the archive is still empty (every candidate failed), the
// batch itself is kept so breeding can go on.
func mapStep(p Problem, m *MapElites, pop, batch []Individual) []Individual {
	archive := map[int]Individual{}
	m.insert(p, archive, pop)
	m.insert(p, archive, batch)
	if len(archive) == 0 {
		return batch
	}
	return archiveList(archive)
}
