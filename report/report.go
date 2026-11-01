// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package report writes the results sheet that accompanies a netlist and
// its KiCad schematic: the MNA solution and Kron's torn (diakoptic)
// solution of the same circuit, side by side.
//
// The analysis point is the netlist's own: the operating point (for .op,
// .dc and .tran) or one AC frequency (for .ac, the geometric middle of the
// sweep unless Options.Freq is given). The circuit is torn automatically
// (kron.AutoTearSolvable, ground as the common bus); the sheet lists the
// torn pieces, the interconnection system's order, each torn branch (Kron's
// links) with its current and the voltage across it, and the largest
// difference between the torn and the MNA node voltages. Nonlinear
// circuits get the MNA part only, with the reason.
package report

import (
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
)

// Options for a sheet.
type Options struct {
	Source  string  // file name shown in the sheet
	Pieces  int     // pieces asked of AutoTear (default 4)
	Freq    float64 // AC frequency [Hz]; 0 = from the netlist
	MaxRows int     // node rows shown (default 40); the CSV has all
	// Tear is a manual cut set (element names; CLAUDE.md §3a rule 5). Empty:
	// automatic tearing (kron.AutoTearSolvable) into Pieces.
	Tear []string
}

// Sheet is the computed content.
type Sheet struct {
	Title, Source string
	Elements      int
	Nodes         []string
	DC            bool
	Freq          float64
	MNA, Torn     []complex128 // node voltages, aligned with Nodes
	KronNote      string       // why there is no Kron part, if none
	Cut           kron.Tearing
	Manual        bool // Cut was declared, not found by AutoTear
	Stats         kron.TearStats
	Links         []Link
	MaxDiff       float64    // max |torn − MNA| / max |MNA|
	PowerMNA      complex128 // V·I* summed over R, L, C branches, untorn
	PowerTorn     complex128 // the same, torn
}

// Link is one torn branch.
type Link struct {
	Name     string
	From, To string
	I, V     complex128
}

// Compute solves the netlist both ways.
func Compute(src string, opt Options) (*Sheet, error) {
	if opt.Pieces == 0 {
		opt.Pieces = 4
	}
	d, err := netlist.ParseString(src)
	if err != nil {
		return nil, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil, err
	}
	s := &Sheet{Title: d.Title, Source: opt.Source, Elements: len(c.Elements), DC: true, Freq: opt.Freq}
	for _, a := range d.Analyses {
		if a.Kind == "ac" && s.Freq == 0 {
			s.Freq = math.Sqrt(math.Max(a.FStart, 1e-30) * a.FStop)
		}
	}
	if s.Freq > 0 {
		s.DC = false
	}
	s.Nodes = append([]string(nil), c.NodeNames[1:]...)
	// MNA.
	if s.DC {
		r, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			return nil, fmt.Errorf("MNA operating point: %w", err)
		}
		for _, n := range s.Nodes {
			v, _ := r.Voltage(n)
			s.MNA = append(s.MNA, complex(v, 0))
		}
	} else {
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: s.Freq, FStop: s.Freq}, analysis.Options{})
		if err != nil {
			return nil, fmt.Errorf("MNA AC: %w", err)
		}
		for _, n := range s.Nodes {
			v, _ := r.NodeVoltage(n)
			s.MNA = append(s.MNA, v[0])
		}
	}
	// Kron.
	net, err := kron.Build(c)
	if err != nil {
		s.KronNote = fmt.Sprintf("Kron's torn solution is shown for linear circuits; this one is not (%v). The torn Newton solve (kron.JacobianTearing) is checked in kron's tests.", err)
		return s, nil
	}
	w := 2 * math.Pi * s.Freq
	var (
		t   kron.Tearing
		dk  *kron.Diakoptics
		sol *kron.Solution
		st  kron.TearStats
	)
	if len(opt.Tear) > 0 {
		t = kron.Tearing{Torn: opt.Tear}
		s.Manual = true
		if dk, err = kron.Tear(c, t, nil); err == nil {
			sol, st, err = dk.Solve(w, s.DC)
		}
		if err != nil {
			s.KronNote = fmt.Sprintf("the declared cut set does not solve: %v", err)
			return s, nil
		}
	} else if t, dk, sol, st, err = kron.AutoTearSolvable(c, opt.Pieces, w, s.DC); err != nil {
		s.KronNote = fmt.Sprintf("no automatic tearing solved this circuit: %v", err)
		return s, nil
	}
	s.Cut, s.Stats = t, st
	if s.DC {
		r, err := kron.TornOP(c, t, nil)
		if err != nil {
			return nil, err
		}
		for _, n := range s.Nodes {
			v, _ := r.Voltage(n)
			s.Torn = append(s.Torn, complex(v, 0))
		}
	} else {
		r, err := kron.TornAC(c, []float64{s.Freq}, t, nil)
		if err != nil {
			return nil, err
		}
		for _, n := range s.Nodes {
			v, _ := r.NodeVoltage(n)
			s.Torn = append(s.Torn, v[0])
		}
	}
	scale := 0.0
	for _, v := range s.MNA {
		scale = math.Max(scale, cmplx.Abs(v))
	}
	for i := range s.MNA {
		s.MaxDiff = math.Max(s.MaxDiff, cmplx.Abs(s.Torn[i]-s.MNA[i])/math.Max(scale, 1e-300))
	}
	ref, err := net.Solve(kron.NodeForm, w, s.DC)
	if err == nil {
		s.PowerMNA, s.PowerTorn = absorbed(c, net, ref), absorbed(c, dk.Network(), sol)
	}
	for _, name := range t.Torn {
		b, ok := dk.Network().ElementBranch(name)
		if !ok {
			continue
		}
		k := c.Index(name)
		e := c.Elements[k]
		s.Links = append(s.Links, Link{Name: e.Name, From: c.NodeNames[e.Nodes[0]], To: c.NodeNames[e.Nodes[1]], I: sol.I.At(b), V: sol.E.At(b)})
	}
	return s, nil
}

// absorbed is the complex power V·I* taken by the R, L and C branches
// (the power the sources deliver, by Tellegen's theorem).
func absorbed(c *circuit.Circuit, n *kron.Network, sol *kron.Solution) complex128 {
	var p complex128
	for _, e := range c.Elements {
		switch e.Kind {
		case circuit.Resistor, circuit.Inductor, circuit.Capacitor:
		default:
			continue
		}
		if b, ok := n.ElementBranch(e.Name); ok {
			p += sol.E.At(b) * cmplx.Conj(sol.I.At(b))
		}
	}
	return p
}

func cfmt(v complex128, dc bool) string {
	if dc {
		return fmt.Sprintf("%.6g", real(v))
	}
	return fmt.Sprintf("%.6g ∠ %.2f°", cmplx.Abs(v), cmplx.Phase(v)*180/math.Pi)
}

// Markdown writes the sheet.
func (s *Sheet) Markdown(w io.Writer, maxRows int) {
	if maxRows <= 0 {
		maxRows = 40
	}
	fmt.Fprintf(w, "# Results: %s\n\n", s.Title)
	fmt.Fprintf(w, "Netlist `%s` — %d elements, %d nodes. Generated by kronspice (`report` package).\n\n", s.Source, s.Elements, len(s.Nodes))
	if s.DC {
		fmt.Fprintf(w, "**Analysis point:** DC operating point.\n\n")
	} else {
		fmt.Fprintf(w, "**Analysis point:** AC at %.6g Hz (the middle of the netlist's sweep unless chosen).\n\n", s.Freq)
	}
	fmt.Fprintf(w, "## Kron's diakoptic solution\n\n")
	if s.KronNote != "" {
		fmt.Fprintf(w, "%s\n\n", s.KronNote)
	} else {
		how := "automatically"
		if s.Manual {
			how = "along the declared cut set"
		}
		fmt.Fprintf(w, "The circuit is torn %s (ground is the common bus of all pieces, Happ Ch. XII §4) into **%d pieces** of %s unknowns each (node voltages, plus the branch currents of voltage sources and inductors), joined by **%d links** (the torn branches). Each piece is factored on its own; the only coupled solve is the interconnection system of order **%d**, against %d unknowns for the whole circuit.\n\n",
			how, s.Stats.Pieces, sizes(s.Stats.Sizes), len(s.Cut.Torn), s.Stats.Interconnection, sum(s.Stats.Sizes))
		fmt.Fprintf(w, "| Link (torn branch) | From | To | Link current | Voltage across |\n|---|---|---|---|---|\n")
		for k, l := range s.Links {
			if k == maxRows {
				fmt.Fprintf(w, "| … | %d more links | | | |\n", len(s.Links)-maxRows)
				break
			}
			fmt.Fprintf(w, "| %s | %s | %s | %s A | %s V |\n", l.Name, l.From, l.To, cfmt(l.I, s.DC), cfmt(l.V, s.DC))
		}
		fmt.Fprintf(w, "\n**Torn against full (MNA) solve:** largest node-voltage difference %.2g of the largest node voltage", s.MaxDiff)
		if s.PowerMNA != 0 {
			fmt.Fprintf(w, "; power V·I* taken by the R, L, C branches %s (torn) against %s (untorn) — power is invariant under tearing (Roth 1959 §5)", cfmt(s.PowerTorn, s.DC), cfmt(s.PowerMNA, s.DC))
		}
		fmt.Fprintf(w, ".\n\n")
	}
	fmt.Fprintf(w, "## Node voltages\n\n| Node | MNA | Kron (torn) |\n|---|---|---|\n")
	idx := make([]int, len(s.Nodes))
	for i := range idx {
		idx[i] = i
	}
	sort.SliceStable(idx, func(a, b int) bool { return natLess(s.Nodes[idx[a]], s.Nodes[idx[b]]) })
	for k, i := range idx {
		if k == maxRows {
			fmt.Fprintf(w, "| … | %d more nodes in the CSV | |\n", len(idx)-maxRows)
			break
		}
		torn := "—"
		if s.Torn != nil {
			torn = cfmt(s.Torn[i], s.DC)
		}
		fmt.Fprintf(w, "| %s | %s | %s |\n", s.Nodes[i], cfmt(s.MNA[i], s.DC), torn)
	}
}

func sum(z []int) int {
	t := 0
	for _, v := range z {
		t += v
	}
	return t
}

// sizes lists piece sizes, or summarizes them when there are many.
func sizes(z []int) string {
	if len(z) <= 12 {
		return fmt.Sprint(z)
	}
	lo, hi := z[0], z[0]
	for _, v := range z {
		lo, hi = min(lo, v), max(hi, v)
	}
	return fmt.Sprintf("%d to %d", lo, hi)
}

// CSV writes every node: name, MNA real/imag, torn real/imag.
func (s *Sheet) CSV(w io.Writer) {
	fmt.Fprintf(w, "node,mna_re,mna_im,kron_re,kron_im\n")
	for i, n := range s.Nodes {
		t := complex(math.NaN(), math.NaN())
		if s.Torn != nil {
			t = s.Torn[i]
		}
		fmt.Fprintf(w, "%s,%.17g,%.17g,%.17g,%.17g\n", n, real(s.MNA[i]), imag(s.MNA[i]), real(t), imag(t))
	}
}

// natLess orders names with numbers in numeric order (n2 before n10).
func natLess(a, b string) bool {
	for a != "" && b != "" {
		da, db := digits(a), digits(b)
		if da > 0 && db > 0 {
			na, nb := strings.TrimLeft(a[:da], "0"), strings.TrimLeft(b[:db], "0")
			if len(na) != len(nb) {
				return len(na) < len(nb)
			}
			if na != nb {
				return na < nb
			}
			a, b = a[da:], b[db:]
			continue
		}
		if a[0] != b[0] {
			return a[0] < b[0]
		}
		a, b = a[1:], b[1:]
	}
	return len(a) < len(b)
}

func digits(s string) int {
	i := 0
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
	}
	return i
}
