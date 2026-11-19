// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
)

// ErrBadConfig is returned for an invalid experiment configuration.
var ErrBadConfig = errors.New("evo: invalid experiment config")

// Config is a JSON experiment description (cmd/evolve -config). Paths are
// relative to the config file.
type Config struct {
	Name    string       `json:"name"`
	Netlist string       `json:"netlist"`
	Genes   []GeneConfig `json:"genes"`
	// Topology evolves the parts themselves into the netlist (the embryo);
	// genes must then be empty.
	Topology   *TopologyConfig   `json:"topology,omitempty"`
	Objectives []ObjectiveConfig `json:"objectives"`
	GA         GAConfig          `json:"ga"`
	OutDir     string            `json:"out_dir,omitempty"` // default "runs"
}

// GeneConfig is one gene. Sign is "positive" (default), "negative" or
// "either"; Series is "", "E6", "E12", "E24", "E48" or "E96".
type GeneConfig struct {
	Element string  `json:"element,omitempty"`
	Model   string  `json:"model,omitempty"`
	Param   string  `json:"param,omitempty"`
	Min     float64 `json:"min"`
	Max     float64 `json:"max"`
	Sign    string  `json:"sign,omitempty"`
	Series  string  `json:"series,omitempty"`
}

// TopologyConfig maps onto Topology.
type TopologyConfig struct {
	Terminals []string         `json:"terminals"`
	Internal  int              `json:"internal_nodes,omitempty"`
	MaxParts  int              `json:"max_parts,omitempty"`
	InitParts int              `json:"init_parts,omitempty"`
	Parts     []PartKindConfig `json:"parts"`
	// Encoding is "direct" (default) or "developmental" (Koza's
	// circuit-constructing programs); the developmental encoding needs
	// Wires, the embryo's modifiable wires as node pairs.
	Encoding  string      `json:"encoding,omitempty"`
	Wires     [][2]string `json:"wires,omitempty"`
	MaxSize   int         `json:"max_size,omitempty"`
	InitDepth int         `json:"init_depth,omitempty"`
}

// PartKindConfig is one kind of part: "R", "C", "L", "V", "I", "D", "Q" or
// "M", its range and series (PartKind).
type PartKindConfig struct {
	Kind   string  `json:"kind"`
	Min    float64 `json:"min"`
	Max    float64 `json:"max"`
	Series string  `json:"series,omitempty"`
	// KappaA, KappaB: hypothetical ECE scalar terms at the part's two
	// terminals (needs .ece in the netlist; PartKind.KappaA).
	KappaA float64 `json:"kappa_a,omitempty"`
	KappaB float64 `json:"kappa_b,omitempty"`
	// Model names the embryo's .model for D, Q and M parts; Length is a
	// MOSFET part's channel length [m].
	Model  string  `json:"model,omitempty"`
	Length float64 `json:"length,omitempty"`
}

// BandConfig is one band of a mask objective; an omitted bound is open.
type BandConfig struct {
	From  float64  `json:"from"`
	To    float64  `json:"to"`
	MinDB *float64 `json:"min_db,omitempty"`
	MaxDB *float64 `json:"max_db,omitempty"`
}

// ObjectiveConfig is one fitness term. Kind selects the measurement; the
// analysis it needs (.op, .ac or .tran) is taken from the netlist.
//
//	cutoff      node, ref, target_hz, edge ("high" default, or "low")   needs .ac
//	gain_at     node, ref, freq_hz, target_db                          needs .ac
//	op_voltage  node, target                                           needs nothing
//	active_bias v_min, i_min, i_max (a transistor in its active region) needs nothing
//	tran_mean   node, from, to, target                                 needs .tran
//	tran_ripple node, from, to, max                                    needs .tran
//	port_voc    plus, minus, target   (open-circuit voltage)           needs nothing
//	port_isc    plus, minus, target   (short-circuit current)          needs nothing
//	model_supplied  max  (limit on power from negative resistances and
//	                      controlled sources at the operating point)
//	mask        node, ref, bands [{from, to, min_db, max_db}]          needs .ac
//	parts       (number of evolved parts; weight = penalty per part)  topology only
//	gain_max    node, ref, freq_hz, target_db (ceiling)   more gain is better   needs .ac
//	bandwidth_max node, ref, target_hz (ceiling)          more bandwidth is better needs .ac
//	supply_power sources [names]   (power they deliver at the operating point [W])
type ObjectiveConfig struct {
	Kind     string       `json:"kind"`
	Weight   float64      `json:"weight,omitempty"` // default 1
	Node     string       `json:"node,omitempty"`
	Ref      string       `json:"ref,omitempty"`
	Plus     string       `json:"plus,omitempty"`
	Minus    string       `json:"minus,omitempty"`
	TargetHz float64      `json:"target_hz,omitempty"`
	FreqHz   float64      `json:"freq_hz,omitempty"`
	TargetDB float64      `json:"target_db,omitempty"`
	Target   float64      `json:"target,omitempty"`
	From     float64      `json:"from,omitempty"`
	To       float64      `json:"to,omitempty"`
	Max      float64      `json:"max,omitempty"`
	Edge     string       `json:"edge,omitempty"`
	Bands    []BandConfig `json:"bands,omitempty"`
	Sources  []string     `json:"sources,omitempty"`
	// VMin, IMin, IMax set active_bias's region (defaults 1 V, 100 µA,
	// 100 mA).
	VMin float64 `json:"v_min,omitempty"`
	IMin float64 `json:"i_min,omitempty"`
	IMax float64 `json:"i_max,omitempty"`
	// Constraint makes the term a constraint under "nsga2" (see
	// Objective.Constraint): e.g. gain_max with constraint=true is
	// "gain ≥ target_db".
	Constraint bool `json:"constraint,omitempty"`
}

// GAConfig maps onto Options.
type GAConfig struct {
	Population  int      `json:"population,omitempty"`
	Generations int      `json:"generations,omitempty"`
	Elite       int      `json:"elite,omitempty"`
	Tournament  int      `json:"tournament,omitempty"`
	Crossover   float64  `json:"crossover,omitempty"`
	Mutation    float64  `json:"mutation,omitempty"`
	Sigma       float64  `json:"sigma,omitempty"`
	SignFlip    float64  `json:"sign_flip,omitempty"`
	Seed        uint64   `json:"seed"`
	Workers     int      `json:"workers,omitempty"`
	TimeoutMS   int      `json:"timeout_ms,omitempty"`
	Target      *float64 `json:"target_fitness,omitempty"`
	// Algorithm is "weighted" (default: minimize Σ weight·error) or
	// "nsga2" (a Pareto front over the objectives).
	Algorithm string `json:"algorithm,omitempty"`
	// LocalSearch adds the multi-part local search (evo.LocalSearch).
	LocalSearch *LocalSearchConfig `json:"local_search,omitempty"`
	// MapElites gives the grid of algorithm "map_elites".
	MapElites *MapElitesConfig `json:"map_elites,omitempty"`
	// Precheck runs kron.Precheck on every candidate before simulation
	// (Problem.Precheck): loops of ideal voltage outputs — and with dc,
	// loops closed by inductors too — are penalized as singular without
	// building a matrix. It never rejects a circuit the simulator could
	// solve, so results are unchanged; only time is saved.
	Precheck *PrecheckConfig `json:"precheck,omitempty"`
}

// PrecheckConfig maps onto kron.PrecheckOptions.
type PrecheckConfig struct {
	DC bool `json:"dc,omitempty"`
}

// MapElitesConfig is the MAP-Elites grid: one entry per dimension.
type MapElitesConfig struct {
	Descriptors []DescriptorConfig `json:"descriptors"`
}

// DescriptorConfig is one grid dimension. Kind is "parts" (evolved parts
// or genes), "genome_size" (parts or program ops), "count" (elements of
// kind Part added to the template, e.g. "Q") or "objective" (the error of
// objective number Index, 0-based, in config order).
type DescriptorConfig struct {
	Kind  string  `json:"kind"`
	Part  string  `json:"part,omitempty"`
	Index int     `json:"index,omitempty"`
	Min   float64 `json:"min"`
	Max   float64 `json:"max"`
	Bins  int     `json:"bins"`
}

// LocalSearchConfig maps onto LocalSearch; zero fields take its defaults.
type LocalSearchConfig struct {
	Size     int     `json:"size,omitempty"`
	Every    int     `json:"every,omitempty"`
	MaxEvals int     `json:"max_evals,omitempty"`
	Step     float64 `json:"step,omitempty"`
}

// Experiment is a loaded, validated configuration ready to run.
type Experiment struct {
	Config  Config
	Raw     []byte // the config file as read, copied into the run directory
	Problem Problem
	Options Options
	// Control holds the netlist's analysis lines, written into every best
	// netlist so it can be re-run as is in kronspice, ngspice, LTspice or
	// QSpice.
	Control []string
}

// LoadExperiment reads and validates a JSON config and its netlist.
func LoadExperiment(path string) (*Experiment, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var cfg Config
	dec := json.NewDecoder(bytes.NewReader(raw))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&cfg); err != nil {
		return nil, fmt.Errorf("%w: %s: %v", ErrBadConfig, path, err)
	}
	if cfg.Netlist == "" {
		return nil, fmt.Errorf("%w: no netlist", ErrBadConfig)
	}
	np := cfg.Netlist
	if !filepath.IsAbs(np) {
		np = filepath.Join(filepath.Dir(path), np)
	}
	src, err := os.ReadFile(np)
	if err != nil {
		return nil, err
	}
	ex, err := NewExperiment(cfg, string(src))
	if err != nil {
		return nil, err
	}
	ex.Raw = raw
	return ex, nil
}

// NewExperiment builds an experiment from a config and netlist text.
func NewExperiment(cfg Config, netlistText string) (*Experiment, error) {
	d, err := netlist.ParseString(netlistText)
	if err != nil {
		return nil, err
	}
	tpl, err := circuit.FromDeck(d)
	if err != nil {
		return nil, err
	}
	ex := &Experiment{Config: cfg, Problem: Problem{Template: tpl}, Control: controlLines(netlistText)}
	if cfg.Name == "" {
		ex.Config.Name = "experiment"
	}
	for _, g := range cfg.Genes {
		sm, err := ParseSignMode(g.Sign)
		if err != nil {
			return nil, err
		}
		ex.Problem.Space = append(ex.Problem.Space, Gene{Element: g.Element, Model: g.Model, Param: g.Param,
			Min: g.Min, Max: g.Max, Sign: sm, Series: g.Series})
	}
	if tc := cfg.Topology; tc != nil {
		if len(cfg.Genes) != 0 {
			return nil, fmt.Errorf("%w: a topology experiment has no genes", ErrBadConfig)
		}
		enc, err := ParseEncoding(tc.Encoding)
		if err != nil {
			return nil, err
		}
		top := &Topology{Terminals: tc.Terminals, Internal: tc.Internal, MaxParts: tc.MaxParts, InitParts: tc.InitParts,
			Encoding: enc, Wires: tc.Wires, MaxSize: tc.MaxSize, InitDepth: tc.InitDepth}
		for _, pk := range tc.Parts {
			k := strings.ToUpper(strings.TrimSpace(pk.Kind))
			if len(k) != 1 {
				return nil, fmt.Errorf("%w: part kind %q", ErrBadConfig, pk.Kind)
			}
			top.Kinds = append(top.Kinds, PartKind{Kind: circuit.Kind(k[0]), Min: pk.Min, Max: pk.Max, Series: pk.Series,
				KappaA: pk.KappaA, KappaB: pk.KappaB, Model: pk.Model, Length: pk.Length})
		}
		if err := top.Validate(tpl); err != nil {
			return nil, err
		}
		ex.Problem.Topology = top
	} else if err := ex.Problem.Space.Validate(tpl); err != nil {
		return nil, err
	}
	var ac *analysis.ACSpec
	var tran *analysis.TranSpec
	for _, a := range d.Analyses {
		switch a.Kind {
		case "ac":
			if ac == nil {
				ac = &analysis.ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}
			}
		case "tran":
			if tran == nil {
				tran = &analysis.TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC}
			}
		}
	}
	needAC := func(k string) (analysis.ACSpec, error) {
		if ac == nil {
			return analysis.ACSpec{}, fmt.Errorf("%w: objective %s needs an .ac line in the netlist", ErrBadConfig, k)
		}
		return *ac, nil
	}
	needTran := func(k string) (analysis.TranSpec, error) {
		if tran == nil {
			return analysis.TranSpec{}, fmt.Errorf("%w: objective %s needs a .tran line in the netlist", ErrBadConfig, k)
		}
		return *tran, nil
	}
	if len(cfg.Objectives) == 0 {
		return nil, fmt.Errorf("%w: no objectives", ErrBadConfig)
	}
	for _, oc := range cfg.Objectives {
		ref := oc.Ref
		if ref == "" {
			ref = "0"
		}
		minus := oc.Minus
		if minus == "" {
			minus = "0"
		}
		var ob Objective
		switch strings.ToLower(oc.Kind) {
		case "cutoff":
			spec, err := needAC(oc.Kind)
			if err != nil {
				return nil, err
			}
			if !(oc.TargetHz > 0) {
				return nil, fmt.Errorf("%w: cutoff needs target_hz > 0", ErrBadConfig)
			}
			ob = Cutoff(oc.Node, ref, spec, oc.TargetHz, !strings.EqualFold(oc.Edge, "low"))
		case "gain_at":
			spec, err := needAC(oc.Kind)
			if err != nil {
				return nil, err
			}
			ob = GainAt(oc.Node, ref, spec, oc.FreqHz, oc.TargetDB)
		case "op_voltage":
			ob = OPVoltage(oc.Node, oc.Target)
		case "active_bias":
			vmin, imin, imax := oc.VMin, oc.IMin, oc.IMax
			if vmin == 0 {
				vmin = 1
			}
			if imin == 0 {
				imin = 1e-4
			}
			if imax == 0 {
				imax = 0.1
			}
			if !(vmin > 0) || !(imin > 0) || !(imax > imin) {
				return nil, fmt.Errorf("%w: active_bias needs v_min, i_min > 0 and i_max > i_min", ErrBadConfig)
			}
			ob = ActiveBias(vmin, imin, imax)
		case "tran_mean":
			spec, err := needTran(oc.Kind)
			if err != nil {
				return nil, err
			}
			ob = TranAverage(oc.Node, spec, oc.From, oc.To, oc.Target)
		case "tran_ripple":
			spec, err := needTran(oc.Kind)
			if err != nil {
				return nil, err
			}
			if !(oc.Max > 0) {
				return nil, fmt.Errorf("%w: tran_ripple needs max > 0", ErrBadConfig)
			}
			ob = TranRipple(oc.Node, spec, oc.From, oc.To, oc.Max)
		case "port_voc":
			ob = PortVoc(oc.Plus, minus, oc.Target)
		case "port_isc":
			ob = PortIsc(oc.Plus, minus, oc.Target)
		case "model_supplied":
			ob = ModelSuppliedLimit(oc.Max)
		case "mask":
			spec, err := needAC(oc.Kind)
			if err != nil {
				return nil, err
			}
			if len(oc.Bands) == 0 {
				return nil, fmt.Errorf("%w: mask needs bands", ErrBadConfig)
			}
			var bands []Band
			for _, b := range oc.Bands {
				band := Band{From: b.From, To: b.To, MinDB: math.Inf(-1), MaxDB: math.Inf(1)}
				if b.MinDB != nil {
					band.MinDB = *b.MinDB
				}
				if b.MaxDB != nil {
					band.MaxDB = *b.MaxDB
				}
				if !(b.From > 0) || !(b.To >= b.From) || b.From < spec.FStart || b.To > spec.FStop || band.MinDB > band.MaxDB {
					return nil, fmt.Errorf("%w: mask band %+v must lie inside the .ac sweep", ErrBadConfig, b)
				}
				bands = append(bands, band)
			}
			ob = FilterMask(oc.Node, ref, spec, bands)
		case "gain_max":
			spec, err := needAC(oc.Kind)
			if err != nil {
				return nil, err
			}
			ob = GainShortfall(oc.Node, ref, spec, oc.FreqHz, oc.TargetDB)
		case "bandwidth_max":
			spec, err := needAC(oc.Kind)
			if err != nil {
				return nil, err
			}
			if !(oc.TargetHz > 0) {
				return nil, fmt.Errorf("%w: bandwidth_max needs target_hz > 0", ErrBadConfig)
			}
			ob = BandwidthShortfall(oc.Node, ref, spec, oc.TargetHz)
		case "supply_power":
			if len(oc.Sources) == 0 {
				return nil, fmt.Errorf("%w: supply_power needs sources", ErrBadConfig)
			}
			for _, src := range oc.Sources {
				if tpl.Index(src) < 0 {
					return nil, fmt.Errorf("%w: supply_power: no element %q", ErrBadConfig, src)
				}
			}
			ob = SupplyPower(oc.Sources...)
		case "parts":
			if cfg.Topology == nil {
				return nil, fmt.Errorf("%w: parts objective needs a topology", ErrBadConfig)
			}
			ob = PartCount()
		default:
			return nil, fmt.Errorf("%w: unknown objective kind %q", ErrBadConfig, oc.Kind)
		}
		ob.Weight = 1
		if ob.Constraint = oc.Constraint; ob.Constraint {
			ob.Name = "constraint " + ob.Name
		}
		if oc.Weight != 0 {
			ob.Weight = oc.Weight
		}
		ex.Problem.Objectives = append(ex.Problem.Objectives, ob)
	}
	g := cfg.GA
	ex.Options = Options{Population: g.Population, Generations: g.Generations, Elite: g.Elite,
		Tournament: g.Tournament, Crossover: g.Crossover, Mutation: g.Mutation, Sigma: g.Sigma,
		SignFlip: g.SignFlip, Seed: g.Seed, Workers: g.Workers,
		Timeout: time.Duration(g.TimeoutMS) * time.Millisecond}
	if pc := g.Precheck; pc != nil {
		ex.Problem.Precheck = &kron.PrecheckOptions{DC: pc.DC}
	}
	if m := g.MapElites; m != nil {
		me := &MapElites{}
		for _, dc := range m.Descriptors {
			var d Descriptor
			switch dc.Kind {
			case "parts":
				d = DescParts(dc.Min, dc.Max, dc.Bins)
			case "genome_size":
				d = DescGenomeSize(dc.Min, dc.Max, dc.Bins)
			case "count":
				k := strings.ToUpper(strings.TrimSpace(dc.Part))
				if len(k) != 1 {
					return nil, fmt.Errorf("%w: count descriptor part %q", ErrBadConfig, dc.Part)
				}
				d = DescCount(circuit.Kind(k[0]), dc.Min, dc.Max, dc.Bins)
			case "objective":
				if dc.Index < 0 || dc.Index >= len(ex.Problem.Objectives) {
					return nil, fmt.Errorf("%w: objective descriptor index %d", ErrBadConfig, dc.Index)
				}
				d = DescTerm(ex.Problem.Objectives[dc.Index].Name, dc.Index, dc.Min, dc.Max, dc.Bins)
			default:
				return nil, fmt.Errorf("%w: descriptor kind %q (parts, genome_size, count or objective)", ErrBadConfig, dc.Kind)
			}
			me.Descriptors = append(me.Descriptors, d)
		}
		ex.Options.MapElites = me
	}
	if l := g.LocalSearch; l != nil {
		ex.Options.LocalSearch = &LocalSearch{Size: l.Size, Every: l.Every, MaxEvals: l.MaxEvals, Step: l.Step}
	}
	if g.Target != nil {
		ex.Options.Target, ex.Options.HasTarget = *g.Target, true
	}
	if ex.Options.Algorithm, err = ParseAlgorithm(g.Algorithm); err != nil {
		return nil, err
	}
	genes := len(ex.Problem.Space)
	if ex.Problem.Topology != nil {
		genes = 1
	}
	if err := ex.Options.defaults(genes); err != nil {
		return nil, err
	}
	return ex, nil
}

// controlLines returns the netlist's analysis lines (.op, .dc, .ac, .tran),
// continuation lines joined.
func controlLines(src string) []string {
	var out []string
	lines := strings.Split(strings.ReplaceAll(src, "\r\n", "\n"), "\n")
	for i := 1; i < len(lines); i++ {
		l := strings.TrimSpace(lines[i])
		low := strings.ToLower(l)
		for _, k := range []string{".op ", ".dc ", ".ac ", ".tran "} {
			if strings.HasPrefix(low+" ", k) {
				for i+1 < len(lines) && strings.HasPrefix(strings.TrimSpace(lines[i+1]), "+") {
					i++
					l += " " + strings.TrimSpace(strings.TrimPrefix(strings.TrimSpace(lines[i]), "+"))
				}
				out = append(out, l)
				break
			}
		}
	}
	return out
}
