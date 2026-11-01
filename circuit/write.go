// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/fdamador/kronspice/exchange"
)

// Clone returns a deep copy of c, so a genome decoder can change element
// values without touching the template shared by other goroutines.
func (c *Circuit) Clone() *Circuit {
	out := &Circuit{
		Title:     c.Title,
		NodeNames: append([]string(nil), c.NodeNames...),
		Elements:  make([]Element, len(c.Elements)),
	}
	if c.NodeDomains != nil {
		out.NodeDomains = append([]exchange.ChannelID(nil), c.NodeDomains...)
	}
	for i, e := range c.Elements {
		e.Wave.Args = append([]float64(nil), e.Wave.Args...)
		e.Curve = append([]float64(nil), e.Curve...)
		e.Scale.Args = append([]float64(nil), e.Scale.Args...)
		e.Extra = append([]int(nil), e.Extra...)
		if e.Params != nil {
			p := make(map[string]float64, len(e.Params))
			for k, v := range e.Params {
				p[k] = v
			}
			e.Params = p
		}
		out.Elements[i] = e
	}
	if c.Models != nil {
		out.Models = make(map[string]Model, len(c.Models))
		for k, m := range c.Models {
			p := make(map[string]float64, len(m.Params))
			for pk, pv := range m.Params {
				p[pk] = pv
			}
			m.Params = p
			out.Models[k] = m
		}
	}
	return out
}

// Index returns the position of the element with the given name
// (case-insensitive), or −1.
func (c *Circuit) Index(name string) int {
	for i, e := range c.Elements {
		if strings.EqualFold(e.Name, name) {
			return i
		}
	}
	return -1
}

// num formats a value so that it parses back to the same float64 in kronspice
// and reads as a plain number in ngspice, LTspice and QSpice (no
// engineering suffix, so "m" versus "meg" cannot be misread).
func num(v float64) string { return strconv.FormatFloat(v, 'g', -1, 64) }

func joinNums(vs []float64) string {
	parts := make([]string, len(vs))
	for i, v := range vs {
		parts[i] = num(v)
	}
	return strings.Join(parts, " ")
}

// WriteNetlist writes c as a SPICE netlist: title line, .model cards,
// elements, then the given control lines (for example ".ac dec 20 10 100k")
// and .end. The result parses back to an identical circuit and is plain
// SPICE3 syntax that other simulators accept (CLAUDE.md §4, experiment
// outputs) — except for the kronspice extensions of multi-domain circuits:
// .domain lines, Y transducers and EXT= on resistors. A transducer can be
// written for other simulators as two G sources (i1 = g·V2, i2 = −g·V1)
// and .domain dropped; only the ledger's attribution of power is lost.
func (c *Circuit) WriteNetlist(w io.Writer, control ...string) error {
	bw := bufio.NewWriter(w)
	title := strings.ReplaceAll(c.Title, "\n", " ")
	if strings.TrimSpace(title) == "" {
		title = "kronspice circuit"
	}
	fmt.Fprintln(bw, title)

	names := make([]string, 0, len(c.Models))
	for k := range c.Models {
		names = append(names, k)
	}
	sort.Strings(names)
	for _, k := range names {
		m := c.Models[k]
		keys := make([]string, 0, len(m.Params))
		for p := range m.Params {
			keys = append(keys, p)
		}
		sort.Strings(keys)
		parts := make([]string, len(keys))
		for i, p := range keys {
			parts[i] = strings.ToUpper(p) + "=" + num(m.Params[p])
		}
		fmt.Fprintf(bw, ".model %s %s(%s)\n", m.Name, strings.ToUpper(m.Type), strings.Join(parts, " "))
	}

	// Energy domains (kronspice extension): one .domain line per domain.
	byDomain := map[exchange.ChannelID][]string{}
	for i, d := range c.NodeDomains {
		if i > 0 && d != exchange.Electrical {
			byDomain[d] = append(byDomain[d], c.NodeNames[i])
		}
	}
	for _, d := range exchange.All() {
		if nodes := byDomain[d]; len(nodes) > 0 {
			fmt.Fprintf(bw, ".domain %s %s\n", d, strings.Join(nodes, " "))
		}
	}

	// SPICE reads an element's kind from the first letter of its name. A
	// builder name that does not start with its kind's letter ("U" for a
	// voltage source) is written with the letter prefixed ("VU"), made
	// unique, and references to it (controlling sources, coupled
	// inductors) follow.
	spice := map[string]string{}
	taken := map[string]bool{}
	for _, e := range c.Elements {
		taken[strings.ToLower(e.Name)] = true
	}
	for _, e := range c.Elements {
		letter := strings.ToUpper(string(rune(e.Kind)))
		if e.Name != "" && strings.EqualFold(e.Name[:1], letter) {
			continue
		}
		n := letter + e.Name
		for k := 2; taken[strings.ToLower(n)]; k++ {
			n = fmt.Sprintf("%s%s_%d", letter, e.Name, k)
		}
		taken[strings.ToLower(n)] = true
		spice[strings.ToLower(e.Name)] = n
	}
	ref := func(name string) string {
		if n, ok := spice[strings.ToLower(name)]; ok {
			return n
		}
		return name
	}
	for _, e := range c.Elements {
		e.Name = ref(e.Name)
		e.CtrlSource = ref(e.CtrlSource)
		e.Couples = [2]string{ref(e.Couples[0]), ref(e.Couples[1])}
		n1, n2 := c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]
		line := []string{e.Name, n1, n2}
		if e.Kind == Mutual {
			line = []string{e.Name, e.Couples[0], e.Couples[1]}
		}
		if e.Kind == Gyrator {
			line = []string{e.Name, n1, n2, c.NodeNames[e.Extra[0]], c.NodeNames[e.Extra[1]]}
		}
		if e.Kind == BearingMotor {
			line = []string{e.Name, n1, n2, c.NodeNames[e.Extra[0]], c.NodeNames[e.Extra[1]]}
		}
		if e.Kind == MOSFET {
			line = []string{e.Name, n1, c.NodeNames[e.Extra[0]], n2, c.NodeNames[e.Extra[1]]}
		}
		if e.Kind == Machine {
			for _, n := range e.Extra {
				line = append(line, c.NodeNames[n])
			}
		}
		if e.Kind == BJT {
			line = []string{e.Name, n1, c.NodeNames[e.Extra[0]], n2}
			if len(e.Extra) > 1 {
				line = append(line, c.NodeNames[e.Extra[1]])
			}
		}
		switch e.Kind {
		case Resistor, Capacitor, Inductor:
			line = append(line, num(e.Value))
			if e.Kind == Resistor && e.Channel != exchange.Electrical {
				line = append(line, "EXT="+e.Channel.String())
			}
			if e.HasIC {
				line = append(line, "IC="+num(e.IC))
			}
		case NonlinearCapacitor, NonlinearInductor:
			// kronspice extension: C … Q=pwl(v q …) or L … FLUX=pwl(i λ …)
			// (or a value), then optional SCALE=wave(…), EXT= and IC=.
			if len(e.Curve) > 0 {
				key := "Q"
				if e.Kind == NonlinearInductor {
					key = "FLUX"
				}
				line = append(line, key+"=pwl("+joinNums(e.Curve)+")")
			} else {
				line = append(line, num(e.Value))
			}
			if e.Scale.Kind != WaveNone {
				line = append(line, "SCALE="+e.Scale.Kind.String()+"("+joinNums(e.Scale.Args)+")")
			}
			if e.Channel != exchange.Electrical {
				line = append(line, "EXT="+e.Channel.String())
			}
			if e.HasIC {
				line = append(line, "IC="+num(e.IC))
			}
		case VoltageSource, CurrentSource:
			if e.HasDC || e.Wave.Kind == WaveNone {
				line = append(line, "DC", num(e.Value))
			}
			if e.ACMag != 0 || e.ACPhase != 0 {
				line = append(line, "AC", num(e.ACMag), num(e.ACPhase))
			}
			if e.Wave.Kind != WaveNone {
				args := make([]string, len(e.Wave.Args))
				for i, a := range e.Wave.Args {
					args[i] = num(a)
				}
				line = append(line, e.Wave.Kind.String()+"("+strings.Join(args, " ")+")")
			}
		case Behavioral:
			var ctrl string
			if e.BCtrl == 'V' {
				ctrl = "V(" + c.NodeNames[e.CtrlNodes[0]]
				if e.CtrlNodes[1] != 0 {
					ctrl += "," + c.NodeNames[e.CtrlNodes[1]]
				}
				ctrl += ")"
			} else {
				ctrl = "I(" + e.CtrlSource + ")"
			}
			pts := make([]string, len(e.Curve))
			for i, v := range e.Curve {
				pts[i] = num(v)
			}
			line = append(line, string(e.BOut)+"=pwl("+ctrl+", "+strings.Join(pts, ", ")+")")
		case VCVS, VCCS:
			line = append(line, c.NodeNames[e.CtrlNodes[0]], c.NodeNames[e.CtrlNodes[1]], num(e.Value))
		case Mutual, Gyrator:
			line = append(line, num(e.Value))
		case CCCS, CCVS:
			line = append(line, e.CtrlSource, num(e.Value))
		case MOSFET:
			line = append(line, e.Model)
			keys := make([]string, 0, len(e.Params))
			for k := range e.Params {
				keys = append(keys, k)
			}
			sort.Strings(keys)
			for _, k := range keys {
				line = append(line, strings.ToUpper(k)+"="+num(e.Params[k]))
			}
		case BearingMotor:
			line = append(line, num(e.Value))
			keys := make([]string, 0, len(e.Params))
			for k := range e.Params {
				keys = append(keys, k)
			}
			sort.Strings(keys)
			for _, k := range keys {
				line = append(line, strings.ToUpper(k)+"="+num(e.Params[k]))
			}
		case Machine:
			line = append(line, e.Model)
		case Diode, BJT:
			line = append(line, e.Model)
			if e.Area != 0 && e.Area != 1 {
				line = append(line, num(e.Area))
			}
		default:
			return fmt.Errorf("%w: %s", ErrUnknownKind, e.Name)
		}
		fmt.Fprintln(bw, strings.Join(line, " "))
	}
	for _, l := range control {
		fmt.Fprintln(bw, l)
	}
	fmt.Fprintln(bw, ".end")
	return bw.Flush()
}
