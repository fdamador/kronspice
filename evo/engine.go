// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"context"
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"runtime"
	"sort"
	"sync"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
)

// ErrBadOptions is returned for invalid engine options.
var ErrBadOptions = errors.New("evo: invalid options")

// Options controls the generational loop. Zero values take the defaults
// noted on each field.
type Options struct {
	Population  int     // default 50
	Generations int     // default 50 (generation 0 is the random start)
	Elite       int     // individuals copied unchanged, default 2
	Tournament  int     // tournament size, default 3
	Crossover   float64 // probability a pair is crossed, default 0.9
	Mutation    float64 // per-gene mutation probability, default 1/genes
	Sigma       float64 // mutation step [decades], default 0.2
	SignFlip    float64 // per-gene sign-flip probability (free-sign genes), default 0.05
	Seed        uint64
	// Algorithm is WeightedSum (default) or NSGA2. Under NSGA2 every
	// objective with a positive weight is a separate objective, Elite is
	// ignored (the combined parent + offspring selection is elitist) and
	// Tournament defaults to 2 (Deb's binary tournament).
	Algorithm Algorithm
	Workers   int // default runtime.NumCPU()
	// Timeout bounds each evaluation (default 10 s), so one bad circuit
	// cannot stall a generation (CLAUDE.md §4). It is wall-clock time, so
	// whether a slow candidate times out depends on the machine and its
	// load (the race detector alone makes evaluations ~10× slower). A
	// timed-out candidate takes PenaltyTimeout and changes the selection,
	// so a run is bit-reproducible across machines only when no candidate
	// timed out: GenStats.TimedOut reports it, and tests that compare an
	// evolved result set a generous Timeout and check it stayed zero.
	Timeout time.Duration
	// Target stops the run once the best fitness is ≤ Target (when
	// HasTarget is set).
	Target    float64
	HasTarget bool
	// Analysis holds the base options for every simulation (tolerances,
	// fallbacks); its Context is set per evaluation.
	Analysis analysis.Options
}

func (o *Options) defaults(genes int) error {
	if o.Population == 0 {
		o.Population = 50
	}
	if o.Generations == 0 {
		o.Generations = 50
	}
	if o.Elite == 0 {
		o.Elite = 2
	}
	if o.Tournament == 0 {
		o.Tournament = 3
		if o.Algorithm == NSGA2 {
			o.Tournament = 2
		}
	}
	if o.Crossover == 0 {
		o.Crossover = 0.9
	}
	if o.Mutation == 0 {
		o.Mutation = 1 / float64(genes)
	}
	if o.Sigma == 0 {
		o.Sigma = 0.2
	}
	if o.SignFlip == 0 {
		o.SignFlip = 0.05
	}
	if o.Workers <= 0 {
		o.Workers = runtime.NumCPU()
	}
	if o.Timeout <= 0 {
		o.Timeout = 10 * time.Second
	}
	switch {
	case o.Population < 2 || o.Population > 1_000_000:
		return fmt.Errorf("%w: population %d", ErrBadOptions, o.Population)
	case o.Generations < 1 || o.Generations > 1_000_000:
		return fmt.Errorf("%w: generations %d", ErrBadOptions, o.Generations)
	case o.Elite < 0 || o.Elite >= o.Population:
		return fmt.Errorf("%w: elite %d", ErrBadOptions, o.Elite)
	case o.Algorithm > NSGA2:
		return fmt.Errorf("%w: algorithm %d", ErrBadOptions, o.Algorithm)
	case o.Tournament < 1:
		return fmt.Errorf("%w: tournament %d", ErrBadOptions, o.Tournament)
	case o.Crossover < 0 || o.Crossover > 1 || o.Mutation < 0 || o.Mutation > 1 || o.SignFlip < 0 || o.SignFlip > 1:
		return fmt.Errorf("%w: rates must be in [0, 1]", ErrBadOptions)
	case !(o.Sigma > 0) || math.IsInf(o.Sigma, 0):
		return fmt.Errorf("%w: sigma %g", ErrBadOptions, o.Sigma)
	}
	return nil
}

// Problem is what is evolved and how it is scored.
//
// Either Template is set and every gene targets one of its elements or
// models, or Build is set and genes are free design parameters: Build turns
// their decoded values into a circuit (for example coil geometry into L, C
// and R via em.Coil). Build must be safe for concurrent use and must not
// panic; an error penalizes the candidate (em.ErrOutsideModelValidity
// counts as an invalid design).
type Problem struct {
	// Topology, when set, evolves the circuit's parts as well as their
	// values: Template is the embryo (fixture, sources, load, analyses) and
	// Space must be empty.
	Topology   *Topology
	Template   *circuit.Circuit
	Build      func(values []float64) (*circuit.Circuit, error)
	Space      Space
	Objectives []Objective
	// Tear, when set, solves the AC analyses of value-genome candidates by
	// Kron's diakoptics (kron.TornAC), with factored pieces cached and
	// shared by every worker: a mutation that changes one piece refactors
	// only that piece (CLAUDE.md §3a; doc/M9b.md).
	Tear *Tear
	// Precheck, when set, runs kron.Precheck on every decoded candidate
	// before any simulation: a candidate whose mesh frame has a
	// zero-impedance loop (a loop of ideal voltage outputs, and with DC of
	// inductors too) is penalized as singular without building a matrix
	// (CLAUDE.md §3b item 1, doc/M9e.md). It never rejects a circuit the
	// simulator could solve.
	Precheck *kron.PrecheckOptions
}

// Tear is a diakoptic cut for Problem.Tear. Torn names the R, C or L
// elements of the template to tear; Cache may be nil (Run then makes one
// for the run) or shared between runs.
//
// The torn solve equals MNA to 1e−9 relative. Candidates it cannot take
// are solved by MNA instead, so tearing never changes which circuits are
// accepted: nonlinear circuits, circuits with idealized active elements
// (controlled sources, transducers, negative resistance — their
// model-supplied power comes from the MNA ledger), and cuts that leave a
// singular or nearly floating piece for this candidate.
type Tear struct {
	Torn  []string
	Cache *kron.TornCache
}

// Circuit decodes a genome into the circuit it stands for and its gene
// (or part) values.
func (p Problem) Circuit(g Genome) (*circuit.Circuit, []float64, error) { return p.circuitFor(g) }

// circuitFor decodes a genome into its circuit and gene values.
func (p Problem) circuitFor(g Genome) (*circuit.Circuit, []float64, error) {
	if p.Topology != nil {
		c, vals := p.Topology.Decode(p.Template, g)
		return c, vals, nil
	}
	if p.Build != nil {
		vals := p.Space.Values(g)
		c, err := p.Build(vals)
		return c, vals, err
	}
	c, vals := p.Space.Decode(p.Template, g)
	return c, vals, nil
}

// Individual is one evaluated candidate.
type Individual struct {
	Genome  Genome    `json:"genome"`
	Values  []float64 `json:"values"`  // decoded gene values (SI)
	Fitness float64   `json:"fitness"` // Σ weight·error, or a penalty; lower is better
	Terms   []float64 `json:"terms"`   // per-objective errors (nil on failure)
	// Failure is the error class for a penalized candidate ("" if none),
	// Error its message.
	Failure string `json:"failure,omitempty"`
	Error   string `json:"error,omitempty"`
	// ModelSupplied is the largest model-supplied power seen while
	// scoring this candidate [W] (negative resistances, controlled sources).
	ModelSupplied float64 `json:"model_supplied"`
	// Rank and Crowding are the NSGA-II Pareto rank (0 = non-dominated) and
	// crowding distance (math.MaxFloat64 for a front's extreme points); both
	// stay zero under WeightedSum.
	Rank     int     `json:"rank,omitempty"`
	Crowding float64 `json:"crowding,omitempty"`
	// Violation is the total weighted constraint violation (NSGA2; 0 when
	// every constraint objective is met).
	Violation float64 `json:"violation,omitempty"`
}

// GenStats summarizes one generation.
type GenStats struct {
	Generation int     `json:"generation"`
	Best       float64 `json:"best"`
	Mean       float64 `json:"mean"`
	Worst      float64 `json:"worst"`
	Diversity  float64 `json:"diversity"` // decades (value genomes) or parts (topology genomes)
	Failed     int     `json:"failed"`
	// TimedOut counts candidates penalized for exceeding Options.Timeout;
	// non-zero means the generation's result depends on machine speed.
	TimedOut     int     `json:"timed_out,omitempty"`
	Evaluations  int     `json:"evaluations"` // cumulative
	BestSupplied float64 `json:"best_model_supplied"`
	BestIndex    int     `json:"best_index"`
	// FrontSize is the number of non-dominated feasible candidates (NSGA2).
	FrontSize int `json:"front_size,omitempty"`
}

// State is the full state after a generation has been evaluated; it is what
// a checkpoint stores and what Observe receives.
type State struct {
	Generation int          `json:"generation"`
	Seed       uint64       `json:"seed"`
	Population []Individual `json:"population"`
	Best       Individual   `json:"best"` // best ever
	History    []GenStats   `json:"history"`
	// Front is the current Pareto front (NSGA2): every non-dominated
	// feasible candidate of the population, in population order.
	Front       []Individual `json:"front,omitempty"`
	Evaluations int          `json:"evaluations"`
	Done        bool         `json:"done"`
}

// Observer is called after each generation is evaluated, in generation
// order, from the goroutine that called Run. Returning an error stops the
// run with that error.
type Observer func(st *State) error

// Run evolves a population. With resume set, it continues from that state
// (a checkpoint) and produces exactly the result an uninterrupted run
// would have.
func Run(ctx context.Context, p Problem, opt Options, resume *State, observe Observer) (*State, error) {
	if (p.Template == nil) == (p.Build == nil) {
		return nil, fmt.Errorf("%w: set exactly one of Template or Build", ErrBadOptions)
	}
	genes := len(p.Space)
	if p.Topology != nil {
		if p.Template == nil || len(p.Space) != 0 {
			return nil, fmt.Errorf("%w: a topology problem needs an embryo Template and no value genes", ErrBadOptions)
		}
		if err := p.Topology.Validate(p.Template); err != nil {
			return nil, err
		}
		genes = 1
	} else if err := p.Space.Validate(p.Template); err != nil {
		return nil, err
	}
	if len(p.Objectives) == 0 {
		return nil, fmt.Errorf("%w: no objectives", ErrBadOptions)
	}
	if p.Tear != nil {
		if p.Topology != nil || p.Template == nil {
			return nil, fmt.Errorf("%w: tearing needs a fixed-topology Template", ErrBadOptions)
		}
		if _, err := kron.Tear(p.Template, kron.Tearing{Torn: p.Tear.Torn}, nil); err != nil && !errors.Is(err, kron.ErrNonlinear) {
			return nil, fmt.Errorf("%w: tearing: %v", ErrBadOptions, err)
		}
		if p.Tear.Cache == nil {
			t := *p.Tear
			t.Cache = kron.NewTornCache(0)
			p.Tear = &t
		}
	}
	for _, ob := range p.Objectives {
		if ob.Eval == nil || !(ob.Weight >= 0) || math.IsInf(ob.Weight, 0) {
			return nil, fmt.Errorf("%w: objective %q", ErrBadOptions, ob.Name)
		}
	}
	if err := opt.defaults(genes); err != nil {
		return nil, err
	}
	if opt.Algorithm == NSGA2 {
		n := 0
		for _, ob := range p.Objectives {
			if ob.Weight > 0 && !ob.Constraint {
				n++
			}
		}
		if n == 0 {
			return nil, fmt.Errorf("%w: NSGA-II needs an objective with a positive weight that is not a constraint", ErrBadOptions)
		}
	}
	ev := &evaluator{p: p, opt: opt}

	var st *State
	if resume != nil {
		if resume.Seed != opt.Seed || len(resume.Population) != opt.Population {
			return nil, fmt.Errorf("%w: checkpoint seed or population size differs from the options", ErrBadOptions)
		}
		for _, ind := range resume.Population {
			if p.Topology == nil && (len(ind.Genome.X) != len(p.Space) || len(ind.Genome.Sign) != len(p.Space)) {
				return nil, fmt.Errorf("%w: checkpoint genome does not match the genes", ErrBadOptions)
			}
		}
		cp := *resume
		st = &cp
		if st.Done {
			return st, nil
		}
	} else {
		r := rngFor(opt.Seed, 0)
		genomes := make([]Genome, opt.Population)
		for i := range genomes {
			if p.Topology != nil {
				genomes[i] = p.Topology.randomGenome(r)
			} else {
				genomes[i] = p.Space.randomGenome(r)
			}
		}
		st = &State{Seed: opt.Seed}
		st.Population = ev.evaluate(ctx, genomes)
		if opt.Algorithm == NSGA2 {
			st.Population = nsgaSurvivors(st.Population, p.Objectives, opt.Population)
		}
		st.Evaluations = len(genomes)
		if err := st.finish(opt, observe); err != nil || st.Done {
			return st, err
		}
	}
	for st.Generation+1 < opt.Generations {
		if err := ctx.Err(); err != nil {
			return st, err
		}
		next := st.Generation + 1
		r := rngFor(opt.Seed, uint64(next))
		if opt.Algorithm == NSGA2 {
			pick := func() int { return crowdedTournament(r, st.Population, opt.Tournament) }
			pop := ev.evaluate(ctx, breed(r, pick, st.Population, opt.Population, p, opt))
			// Parents keep their evaluation: the simulator is deterministic.
			all := append(append([]Individual(nil), st.Population...), pop...)
			st.Population = nsgaSurvivors(all, p.Objectives, opt.Population)
			st.Generation = next
			st.Evaluations += len(pop)
			if err := st.finish(opt, observe); err != nil || st.Done {
				return st, err
			}
			continue
		}
		fit := make([]float64, len(st.Population))
		for i, ind := range st.Population {
			fit[i] = ind.Fitness
		}
		pick := func() int { return tournament(r, fit, opt.Tournament) }
		pop := ev.evaluate(ctx, breed(r, pick, st.Population, opt.Population-opt.Elite, p, opt))
		// Elites keep their evaluation: the simulator is deterministic.
		elite := elites(st.Population, opt.Elite)
		st.Population = append(elite, pop...)
		st.Generation = next
		st.Evaluations += len(pop)
		if err := st.finish(opt, observe); err != nil || st.Done {
			return st, err
		}
	}
	st.Done = true
	if observe != nil {
		if err := observe(st); err != nil {
			return st, err
		}
	}
	return st, nil
}

// finish records statistics and the best individual, calls the observer
// and checks the stopping target.
func (st *State) finish(opt Options, observe Observer) error {
	gs := GenStats{Generation: st.Generation, Best: math.Inf(1), Worst: math.Inf(-1), Evaluations: st.Evaluations}
	sum := 0.0
	for i, ind := range st.Population {
		sum += ind.Fitness
		if ind.Fitness < gs.Best {
			gs.Best, gs.BestIndex = ind.Fitness, i
		}
		gs.Worst = math.Max(gs.Worst, ind.Fitness)
		if ind.Failure != "" {
			gs.Failed++
			if ind.Failure == "timeout" {
				gs.TimedOut++
			}
		}
	}
	gs.Mean = sum / float64(len(st.Population))
	gs.Diversity = diversity(st.Population)
	best := st.Population[gs.BestIndex]
	gs.BestSupplied = best.ModelSupplied
	if opt.Algorithm == NSGA2 {
		st.Front = ParetoFront(st.Population)
		gs.FrontSize = len(st.Front)
	}
	if len(st.History) == 0 || best.Fitness < st.Best.Fitness {
		st.Best = best
	}
	st.History = append(st.History, gs)
	if opt.HasTarget && st.Best.Fitness <= opt.Target {
		st.Done = true
	}
	if observe != nil {
		return observe(st)
	}
	return nil
}

// elites returns the n best individuals (ties: lower index first).
func elites(pop []Individual, n int) []Individual {
	idx := make([]int, len(pop))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool { return pop[idx[a]].Fitness < pop[idx[b]].Fitness })
	out := make([]Individual, n)
	for i := 0; i < n; i++ {
		out[i] = pop[idx[i]]
	}
	return out
}

// breed builds n offspring by selection (pick returns a parent index from
// pop, drawing from r), crossover and mutation, all drawn from
// one random stream in a fixed order.
func breed(r *rand.Rand, pick func() int, pop []Individual, n int, p Problem, opt Options) []Genome {
	s, top := p.Space, p.Topology
	out := make([]Genome, 0, n+1)
	for len(out) < n {
		a := pop[pick()].Genome
		b := pop[pick()].Genome
		c1, c2 := a.Clone(), b.Clone()
		if top != nil {
			if r.Float64() < opt.Crossover {
				c1, c2 = top.crossover(r, a, b)
			}
			// Each child is mutated with probability Mutation (default 1).
			for _, c := range []*Genome{&c1, &c2} {
				if r.Float64() < opt.Mutation {
					top.mutate(r, c, opt.Sigma)
				}
				top.repair(c)
			}
			out = append(out, c1, c2)
			continue
		}
		if r.Float64() < opt.Crossover {
			c1, c2 = crossover(r, a, b)
		}
		mutate(r, &c1, opt.Mutation, opt.Sigma, opt.SignFlip)
		mutate(r, &c2, opt.Mutation, opt.Sigma, opt.SignFlip)
		s.repair(&c1)
		s.repair(&c2)
		out = append(out, c1, c2)
	}
	return out[:n]
}

// evaluator scores genomes on a worker pool.
type evaluator struct {
	p   Problem
	opt Options
}

// evaluate scores genomes in parallel; results are stored by index, so the
// order goroutines finish in never matters.
func (ev *evaluator) evaluate(ctx context.Context, genomes []Genome) []Individual {
	out := make([]Individual, len(genomes))
	jobs := make(chan int)
	var wg sync.WaitGroup
	workers := ev.opt.Workers
	if workers > len(genomes) {
		workers = len(genomes)
	}
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				out[i] = ev.one(ctx, genomes[i])
			}
		}()
	}
	for i := range genomes {
		jobs <- i
	}
	close(jobs)
	wg.Wait()
	return out
}

// one decodes and scores a single genome. It never panics: an unexpected
// panic in the simulator is caught, penalized and reported as a bug.
func (ev *evaluator) one(parent context.Context, g Genome) (ind Individual) {
	ind = Individual{Genome: g}
	defer func() {
		if r := recover(); r != nil {
			ind.Fitness, ind.Failure = PenaltyPanic, "panic"
			ind.Error = fmt.Sprintf("panic while decoding (a bug — please report): %v", r)
			ind.Terms = nil
		}
	}()
	c, vals, err := ev.p.circuitFor(g)
	ind.Values = vals
	if err != nil {
		ind.Fitness, ind.Failure = Penalty(err)
		ind.Error = err.Error()
		return ind
	}
	ctx, cancel := context.WithTimeout(parent, ev.opt.Timeout)
	defer cancel()
	defer func() {
		if r := recover(); r != nil {
			ind.Fitness, ind.Failure = PenaltyPanic, "panic"
			ind.Error = fmt.Sprintf("simulator panic (a bug — please report): %v", r)
			ind.Terms = nil
		}
	}()
	fail := func(err error) Individual {
		ind.Fitness, ind.Failure = Penalty(err)
		ind.Error = err.Error()
		ind.Terms = nil
		return ind
	}
	if err := c.Validate(); err != nil {
		if errors.Is(err, circuit.ErrFloatingNode) {
			// A transient may still be defined (UIC); let the objectives decide.
			if err2 := c.ValidateTransient(); err2 != nil {
				return fail(err)
			}
		} else {
			return fail(err)
		}
	}
	if ev.p.Precheck != nil {
		if _, err := kron.Precheck(c, *ev.p.Precheck); err != nil {
			return fail(err)
		}
	}
	e := &Eval{Ctx: ctx, Circuit: c, Values: vals, Options: ev.opt.Analysis, Parts: -1, tear: ev.p.Tear}
	if ev.p.Topology != nil {
		e.Parts = len(vals)
	}
	defer func() { ind.ModelSupplied = e.ModelSupplied }()
	total := 0.0
	ind.Terms = make([]float64, len(ev.p.Objectives))
	for k, ob := range ev.p.Objectives {
		v, err := ob.Eval(e)
		if err == nil && (math.IsNaN(v) || math.IsInf(v, 0) || v < 0) {
			err = fmt.Errorf("%w: objective %q returned %g", analysis.ErrNoMeasurement, ob.Name, v)
		}
		if err != nil {
			return fail(fmt.Errorf("%s: %w", ob.Name, err))
		}
		ind.Terms[k] = v
		total += ob.Weight * v
	}
	ind.Fitness = total
	return ind
}
