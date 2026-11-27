// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package netlist

import (
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// Preprocessing: .include, .param, .global and .subckt/.ends with X
// instances, in the syntax of the ngspice manual (ch. 2.4 ".include",
// 2.6 ".subckt", 2.8 ".param" and "Parameters, functions, expressions").
// The result is a flat list of element and control lines, so every
// element kind parseElement knows works inside a subcircuit.
//
// Scoping is lexical. A subcircuit body sees its own parameters (defaults
// from the .subckt line, values from the X line, .param lines inside the
// body), then the parameters of the scope it was defined in, up to the
// global .param lines. Values on an X line are evaluated where the X line
// is written. A parameter may refer to any other in scope, in any order;
// a circular definition is an error.
//
// Naming follows ngspice. Inside instance X1, element R1 becomes R.X1.R1
// (its kind letter kept first) and internal node n becomes x1.n; ports are
// replaced by the nodes on the X line; ground (0, gnd) and .global nodes
// keep their names. Nested instances add to the path: R.X1.X2.R1, x1.x2.n.
// References to other elements (F, H, W and B controlling sources, K's
// inductors) are renamed the same way, so they point inside the instance.
// A .model inside a .subckt is local: it is renamed x1.name and the
// instance's elements use it.

// ErrSubckt is returned (wrapped with ErrSyntax) for a malformed .subckt
// or X line: unknown subcircuit, wrong node count, missing .ends,
// recursion.
var ErrSubckt = errors.New("netlist: subcircuit error")

// ErrInclude is returned (wrapped with ErrSyntax) when an .include file
// cannot be read, is circular, or no file system was given.
var ErrInclude = errors.New("netlist: include error")

// Limits on preprocessing, so hostile or runaway netlists return an error
// instead of exhausting memory.
const (
	// MaxLines bounds the logical lines of a netlist with its includes.
	MaxLines = 1 << 22
	// MaxIncludeDepth bounds nested .include files.
	MaxIncludeDepth = 16
	// MaxSubcktDepth bounds nested X instances.
	MaxSubcktDepth = 64
	// MaxElements bounds the elements after subcircuit expansion
	// (circuit.MaxElements).
	MaxElements = 1000000
	// MaxInstances bounds the X instances expanded, so empty subcircuits
	// instanced many times over cannot run without end.
	MaxInstances = 1000000
)

// Options control ParseWith.
type Options struct {
	// Open opens an .include file. Nil refuses .include.
	Open func(name string) (io.ReadCloser, error)
	// Dir is the directory relative .include paths of the main netlist
	// resolve against; an included file's own includes resolve against
	// its directory.
	Dir string
}

// ParseFile reads the netlist in path; .include paths resolve against the
// directory of the file that names them, as in ngspice.
func ParseFile(path string) (*Deck, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseWith(f, FileOptions(filepath.Dir(path)))
}

// FileOptions returns Options that read .include files from the operating
// system's file system, relative paths resolving against dir.
func FileOptions(dir string) Options {
	return Options{Open: func(name string) (io.ReadCloser, error) { return os.Open(name) }, Dir: dir}
}

// expandIncludes replaces each .include (or .inc) line with the lines of
// the named file, recursively.
func (o Options) expandIncludes(lines []logicalLine, dir string, depth int, stack []string, total *int) ([]logicalLine, error) {
	var out []logicalLine
	for _, ll := range lines {
		f := strings.Fields(stripComment(ll.text))
		if len(f) == 0 {
			continue
		}
		head := strings.ToLower(f[0])
		if head != ".include" && head != ".inc" {
			out = append(out, ll)
			continue
		}
		where := fmt.Sprintf("%sline %d", filePrefix(ll.file), ll.line)
		name := strings.Trim(strings.TrimSpace(strings.Join(f[1:], " ")), `"'`)
		if name == "" {
			return nil, fmt.Errorf("%w: %w: %s: .include needs a file name", ErrSyntax, ErrInclude, where)
		}
		if o.Open == nil {
			return nil, fmt.Errorf("%w: %w: %s: .include %s: no file system (use ParseFile or ParseWith)", ErrSyntax, ErrInclude, where, name)
		}
		if depth >= MaxIncludeDepth {
			return nil, fmt.Errorf("%w: %w: %s: .include nested deeper than %d", ErrSyntax, ErrInclude, where, MaxIncludeDepth)
		}
		path := name
		if !filepath.IsAbs(path) && dir != "" {
			path = filepath.Join(dir, path)
		}
		path = filepath.Clean(path)
		for _, s := range stack {
			if s == path {
				return nil, fmt.Errorf("%w: %w: %s: .include %s includes itself", ErrSyntax, ErrInclude, where, name)
			}
		}
		rc, err := o.Open(path)
		if err != nil {
			return nil, fmt.Errorf("%w: %w: %s: %v", ErrSyntax, ErrInclude, where, err)
		}
		sub, err := readLogical(rc, false, path)
		rc.Close()
		if err != nil {
			return nil, err
		}
		*total += len(sub)
		if *total > MaxLines {
			return nil, fmt.Errorf("%w: %w: more than %d lines with includes", ErrSyntax, ErrInclude, MaxLines)
		}
		// An included library may end with .end; it ends only that file.
		for i, s := range sub {
			if g := strings.Fields(stripComment(s.text)); len(g) > 0 && strings.EqualFold(g[0], ".end") {
				sub = sub[:i]
				break
			}
		}
		sub, err = o.expandIncludes(sub, filepath.Dir(path), depth+1, append(stack, path), total)
		if err != nil {
			return nil, err
		}
		out = append(out, sub...)
	}
	return out, nil
}

// assign is one name=expression pair.
type assign struct {
	name, expr string
}

// splitWords splits s into plain words and name=value assignments. A value
// is a {braced} or 'quoted' expression (braces and quotes removed) or the
// text up to the next blank; blanks around "=" are allowed. Commas
// between items are blanks.
func splitWords(s string) (words []string, assigns []assign, err error) {
	i := 0
	blank := func(c byte) bool { return c == ' ' || c == '\t' || c == ',' }
	skip := func() {
		for i < len(s) && blank(s[i]) {
			i++
		}
	}
	for {
		skip()
		if i >= len(s) {
			return words, assigns, nil
		}
		j := i
		for j < len(s) && !blank(s[j]) && s[j] != '=' {
			j++
		}
		w := s[i:j]
		i = j
		skip()
		if i >= len(s) || s[i] != '=' {
			if w == "" {
				return nil, nil, fmt.Errorf("%w: %w: unexpected '='", ErrSyntax, ErrParam)
			}
			words = append(words, w)
			continue
		}
		i++ // '='
		skip()
		if w == "" || !validName(w) {
			return nil, nil, fmt.Errorf("%w: %w: %q is not a parameter name", ErrSyntax, ErrParam, w)
		}
		if i >= len(s) {
			return nil, nil, fmt.Errorf("%w: %w: %s= has no value", ErrSyntax, ErrParam, w)
		}
		var v string
		switch s[i] {
		case '{':
			k, err := matchBrace(s, i)
			if err != nil {
				return nil, nil, err
			}
			v, i = s[i+1:k], k+1
		case '\'':
			k := strings.IndexByte(s[i+1:], '\'')
			if k < 0 {
				return nil, nil, fmt.Errorf("%w: %w: unterminated quote in %s=", ErrSyntax, ErrParam, w)
			}
			v, i = s[i+1:i+1+k], i+2+k
		default:
			j := i
			for j < len(s) && !blank(s[j]) {
				j++
			}
			v, i = s[i:j], j
		}
		assigns = append(assigns, assign{name: strings.ToLower(w), expr: v})
	}
}

// matchBrace returns the index of the '}' closing the '{' at s[i].
func matchBrace(s string, i int) (int, error) {
	depth := 0
	for k := i; k < len(s); k++ {
		switch s[k] {
		case '{':
			depth++
		case '}':
			depth--
			if depth == 0 {
				return k, nil
			}
		}
	}
	return 0, fmt.Errorf("%w: %w: unmatched '{'", ErrSyntax, ErrParam)
}

func validName(w string) bool {
	if w == "" || !isIdentStart(w[0]) {
		return false
	}
	for i := 1; i < len(w); i++ {
		if !isIdentChar(w[i]) {
			return false
		}
	}
	return true
}

// scope is one level of parameters, subcircuit definitions and local
// models: the global level, or one subcircuit instance.
type scope struct {
	parent *scope
	exprs  map[string]string  // unevaluated parameters
	vals   map[string]float64 // evaluated (and X-line) values
	busy   map[string]bool    // being evaluated (cycle check)
	defs   map[string]*subckt
	models map[string]string // local model name → renamed model
}

func newScope(parent *scope) *scope {
	return &scope{parent: parent, exprs: map[string]string{}, vals: map[string]float64{},
		busy: map[string]bool{}, defs: map[string]*subckt{}, models: map[string]string{}}
}

// value resolves a parameter through the scope chain, evaluating it in the
// scope that defines it.
func (s *scope) value(name string, depth int) (float64, error) {
	for c := s; c != nil; c = c.parent {
		if v, ok := c.vals[name]; ok {
			return v, nil
		}
		e, ok := c.exprs[name]
		if !ok {
			continue
		}
		if c.busy[name] {
			return 0, fmt.Errorf("%w: %w: parameter %s is defined in terms of itself", ErrSyntax, ErrParam, name)
		}
		c.busy[name] = true
		v, err := evalExpr(e, c, depth)
		delete(c.busy, name)
		if err != nil {
			return 0, fmt.Errorf("parameter %s: %w", name, err)
		}
		c.vals[name] = v
		return v, nil
	}
	return 0, errUndefined
}

func (s *scope) define(a assign) error {
	if _, dup := s.exprs[a.name]; dup {
		return fmt.Errorf("%w: %w: parameter %s defined twice", ErrSyntax, ErrParam, a.name)
	}
	if _, dup := s.vals[a.name]; dup {
		return fmt.Errorf("%w: %w: parameter %s defined twice", ErrSyntax, ErrParam, a.name)
	}
	if _, fn := exprFuncs[a.name]; fn {
		return fmt.Errorf("%w: %w: %s is a built-in name", ErrSyntax, ErrParam, a.name)
	}
	s.exprs[a.name] = a.expr
	return nil
}

// findDef looks a subcircuit up through the scope chain and returns it with
// the scope it was defined in.
func (s *scope) findDef(name string) (*subckt, *scope) {
	for c := s; c != nil; c = c.parent {
		if d, ok := c.defs[name]; ok {
			return d, c
		}
	}
	return nil, nil
}

func (s *scope) model(name string) string {
	for c := s; c != nil; c = c.parent {
		if m, ok := c.models[name]; ok {
			return m
		}
	}
	return name
}

// substitute replaces every {expression} and 'expression' in text by its
// value in scope s.
func substitute(text string, s *scope) (string, error) {
	if !strings.ContainsAny(text, "{'") {
		return text, nil
	}
	var b strings.Builder
	for i := 0; i < len(text); {
		c := text[i]
		var expr string
		switch c {
		case '{':
			k, err := matchBrace(text, i)
			if err != nil {
				return "", err
			}
			expr, i = text[i+1:k], k+1
		case '\'':
			k := strings.IndexByte(text[i+1:], '\'')
			if k < 0 {
				return "", fmt.Errorf("%w: %w: unterminated quote", ErrSyntax, ErrParam)
			}
			expr, i = text[i+1:i+1+k], i+2+k
		default:
			b.WriteByte(c)
			i++
			continue
		}
		v, err := evalExpr(expr, s, 0)
		if err != nil {
			return "", err
		}
		b.WriteString(strconv.FormatFloat(v, 'g', -1, 64))
	}
	return b.String(), nil
}

// subckt is one .subckt definition.
type subckt struct {
	name     string
	ports    []string
	defaults []assign
	body     []logicalLine
	defs     map[string]*subckt // nested definitions
}

func headOf(ll logicalLine) string {
	f := strings.Fields(stripComment(ll.text))
	if len(f) == 0 {
		return ""
	}
	return strings.ToLower(f[0])
}

func where(ll logicalLine) string { return fmt.Sprintf("%sline %d", filePrefix(ll.file), ll.line) }

// splitDefs removes the .subckt … .ends blocks from lines and returns the
// rest and the definitions (each with its own nested definitions split).
func splitDefs(lines []logicalLine) ([]logicalLine, map[string]*subckt, error) {
	var rest []logicalLine
	defs := map[string]*subckt{}
	for i := 0; i < len(lines); i++ {
		ll := lines[i]
		switch headOf(ll) {
		case ".ends":
			return nil, nil, fmt.Errorf("%w: %w: %s: .ends without .subckt", ErrSyntax, ErrSubckt, where(ll))
		case ".subckt":
		default:
			rest = append(rest, ll)
			continue
		}
		words, assigns, err := splitWords(stripComment(ll.text))
		if err != nil {
			return nil, nil, fmt.Errorf("%s: %w", where(ll), err)
		}
		var ports []string
		for _, w := range words[1:] {
			if !strings.EqualFold(w, "params:") {
				ports = append(ports, strings.ToLower(w))
			}
		}
		if len(ports) == 0 {
			return nil, nil, fmt.Errorf("%w: %w: %s: .subckt needs: name node [node …] [params: name=value …]", ErrSyntax, ErrSubckt, where(ll))
		}
		def := &subckt{name: strings.ToLower(ports[0]), ports: ports[1:], defaults: assigns}
		seen := map[string]bool{}
		for _, p := range def.ports {
			if p == "0" || p == "gnd" || seen[p] {
				return nil, nil, fmt.Errorf("%w: %w: %s: .subckt %s: port %q is ground or repeated", ErrSyntax, ErrSubckt, where(ll), def.name, p)
			}
			seen[p] = true
		}
		depth, j := 1, i+1
		for ; j < len(lines); j++ {
			switch headOf(lines[j]) {
			case ".subckt":
				depth++
			case ".ends":
				depth--
			}
			if depth == 0 {
				break
			}
		}
		if depth != 0 {
			return nil, nil, fmt.Errorf("%w: %w: %s: .subckt %s has no .ends", ErrSyntax, ErrSubckt, where(ll), def.name)
		}
		if f := strings.Fields(stripComment(lines[j].text)); len(f) > 1 && !strings.EqualFold(f[1], def.name) {
			return nil, nil, fmt.Errorf("%w: %w: %s: .ends %s closes .subckt %s", ErrSyntax, ErrSubckt, where(lines[j]), f[1], def.name)
		}
		def.body, def.defs, err = splitDefs(lines[i+1 : j])
		if err != nil {
			return nil, nil, err
		}
		if _, dup := defs[def.name]; dup {
			return nil, nil, fmt.Errorf("%w: %w: %s: .subckt %s defined twice", ErrSyntax, ErrSubckt, where(ll), def.name)
		}
		defs[def.name] = def
		i = j
	}
	return rest, defs, nil
}

// instance is one expanded X line.
type instance struct {
	path  string            // "X1", "X1.X2"
	ports map[string]string // port name → actual node (already renamed)
	depth int
}

// expander turns the preprocessed lines into the deck.
type expander struct {
	d       *Deck
	global  *scope
	globals map[string]bool // .global nodes
	insts   int             // X instances expanded so far
}

func newExpander(d *Deck) *expander {
	return &expander{d: d, global: newScope(nil), globals: map[string]bool{}}
}

func (x *expander) run(lines []logicalLine) error {
	for i, ll := range lines { // nothing after .end is read
		if headOf(ll) == ".end" {
			lines = lines[:i]
			break
		}
	}
	rest, defs, err := splitDefs(lines)
	if err != nil {
		return err
	}
	x.global.defs = defs
	if err := x.emit(rest, x.global, nil); err != nil {
		return err
	}
	// Every global parameter must evaluate; record them for callers.
	names := make([]string, 0, len(x.global.exprs))
	for n := range x.global.exprs {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, n := range names {
		v, err := x.global.value(n, 0)
		if err != nil {
			return err
		}
		if x.d.Params == nil {
			x.d.Params = map[string]float64{}
		}
		x.d.Params[n] = v
	}
	return nil
}

// collect takes the .param and .global lines of one level into scope s and
// returns the other lines.
func (x *expander) collect(lines []logicalLine, s *scope) ([]logicalLine, error) {
	var rest []logicalLine
	for _, ll := range lines {
		switch headOf(ll) {
		case ".param":
			words, assigns, err := splitWords(stripComment(ll.text))
			if err == nil && (len(words) != 1 || len(assigns) == 0) {
				err = fmt.Errorf("%w: %w: .param needs: name=value [name=value …]", ErrSyntax, ErrParam)
			}
			for _, a := range assigns {
				if err == nil {
					err = s.define(a)
				}
			}
			if err != nil {
				return nil, fmt.Errorf("%s: %w", where(ll), err)
			}
		case ".global":
			f := strings.Fields(stripComment(ll.text))
			for _, n := range f[1:] {
				x.globals[strings.ToLower(n)] = true
			}
		default:
			rest = append(rest, ll)
		}
	}
	return rest, nil
}

// emit processes one level's lines in scope s (inst nil at the top level).
// Inside an instance the local .model lines go first, so an element may
// name a model defined after it, as at the top level.
func (x *expander) emit(lines []logicalLine, s *scope, inst *instance) error {
	lines, err := x.collect(lines, s)
	if err != nil {
		return err
	}
	if inst != nil {
		var models, rest []logicalLine
		for _, ll := range lines {
			if f := strings.Fields(stripComment(ll.text)); len(f) > 0 && strings.EqualFold(f[0], ".model") {
				models = append(models, ll)
			} else {
				rest = append(rest, ll)
			}
		}
		lines = append(models, rest...)
	}
	for _, ll := range lines {
		if err := x.line(ll, s, inst); err != nil {
			if ll.file != "" && !strings.Contains(err.Error(), ll.file) {
				err = fmt.Errorf("%s: %w", ll.file, err)
			}
			return err
		}
	}
	return nil
}

func (x *expander) line(ll logicalLine, s *scope, inst *instance) error {
	text := stripComment(ll.text)
	fields := strings.Fields(text)
	if len(fields) == 0 {
		return nil
	}
	head := strings.ToLower(fields[0])
	if head != ".title" {
		t, err := substitute(text, s)
		if err != nil {
			return fmt.Errorf("line %d: %w", ll.line, err)
		}
		text = t
		fields = strings.Fields(text)
	}
	if strings.HasPrefix(head, ".") {
		if inst == nil {
			_, err := x.d.control(head, fields, logicalLine{text: text, line: ll.line, file: ll.file})
			return err
		}
		if head != ".model" {
			return fmt.Errorf("%w: %w: line %d: %s is not allowed inside .subckt", ErrSyntax, ErrSubckt, ll.line, fields[0])
		}
		m, err := parseModel(text, ll.line)
		if err != nil {
			return err
		}
		local := strings.ToLower(inst.path) + "." + m.Name
		s.models[m.Name] = local
		m.Name = local
		x.d.Models = append(x.d.Models, m)
		return nil
	}
	if upper(fields[0][0]) == 'X' {
		return x.instantiate(text, ll, s, inst)
	}
	el, err := parseElement(text, ll.line)
	if err != nil {
		return err
	}
	if inst != nil {
		x.rename(&el, s, inst)
	}
	if len(x.d.Elements) >= MaxElements {
		return fmt.Errorf("%w: %w: more than %d elements after subcircuit expansion", ErrSyntax, ErrSubckt, MaxElements)
	}
	x.d.Elements = append(x.d.Elements, el)
	return nil
}

// elementName is the name of element n inside an instance at path.
func elementName(n, path string) string {
	if n == "" {
		return n
	}
	return n[:1] + "." + path + "." + n
}

func (x *expander) node(n string, inst *instance) string {
	if inst == nil || n == "0" || n == "gnd" || x.globals[n] {
		return n
	}
	if a, ok := inst.ports[n]; ok {
		return a
	}
	return strings.ToLower(inst.path) + "." + n
}

// rename moves an element parsed inside a subcircuit into its instance.
// Every field that names a node, an element or a model is listed here; a
// new ElementLine field of that kind must be added.
func (x *expander) rename(el *ElementLine, s *scope, inst *instance) {
	el.Name = elementName(el.Name, inst.path)
	for i, n := range el.Nodes {
		el.Nodes[i] = x.node(n, inst)
	}
	for i, n := range el.CtrlNodes {
		el.CtrlNodes[i] = x.node(n, inst)
	}
	if el.CtrlSource != "" {
		el.CtrlSource = strings.ToLower(elementName(el.CtrlSource, inst.path))
	}
	for i, n := range el.Couples {
		el.Couples[i] = strings.ToLower(elementName(n, inst.path))
	}
	if el.Model != "" {
		el.Model = s.model(el.Model)
	}
}

// instantiate expands "Xname node … subckt [params:] [name=value …]".
func (x *expander) instantiate(text string, ll logicalLine, caller *scope, parent *instance) error {
	words, assigns, err := splitWords(text)
	if err != nil {
		return fmt.Errorf("line %d: %w", ll.line, err)
	}
	var w []string
	for _, v := range words {
		if !strings.EqualFold(v, "params:") {
			w = append(w, v)
		}
	}
	if len(w) < 2 {
		return fmt.Errorf("%w: %w: line %d: X needs: Xname node … subckt [name=value …]", ErrSyntax, ErrSubckt, ll.line)
	}
	xname, subName, nodes := w[0], strings.ToLower(w[len(w)-1]), w[1:len(w)-1]
	def, defScope := caller.findDef(subName)
	if def == nil {
		return fmt.Errorf("%w: %w: line %d: %s: no .subckt %s", ErrSyntax, ErrSubckt, ll.line, xname, subName)
	}
	if len(nodes) != len(def.ports) {
		return fmt.Errorf("%w: %w: line %d: %s has %d nodes, .subckt %s has %d ports", ErrSyntax, ErrSubckt, ll.line, xname, len(nodes), def.name, len(def.ports))
	}
	depth, path := 1, xname
	if parent != nil {
		depth, path = parent.depth+1, parent.path+"."+xname
	}
	if depth > MaxSubcktDepth {
		return fmt.Errorf("%w: %w: line %d: subcircuits nested deeper than %d (recursive .subckt?)", ErrSyntax, ErrSubckt, ll.line, MaxSubcktDepth)
	}
	if x.insts++; x.insts > MaxInstances {
		return fmt.Errorf("%w: %w: line %d: more than %d subcircuit instances", ErrSyntax, ErrSubckt, ll.line, MaxInstances)
	}
	inst := &instance{path: path, ports: make(map[string]string, len(nodes)), depth: depth}
	for i, p := range def.ports {
		inst.ports[p] = x.node(strings.ToLower(nodes[i]), parent)
	}
	s := newScope(defScope)
	s.defs = def.defs
	declared := map[string]bool{}
	for _, a := range def.defaults {
		if err := s.define(a); err != nil {
			return fmt.Errorf("line %d: .subckt %s: %w", ll.line, def.name, err)
		}
		declared[a.name] = true
	}
	for _, a := range assigns {
		if !declared[a.name] {
			return fmt.Errorf("%w: %w: line %d: %s: .subckt %s has no parameter %s", ErrSyntax, ErrParam, ll.line, xname, def.name, a.name)
		}
		if _, dup := s.vals[a.name]; dup {
			return fmt.Errorf("%w: %w: line %d: %s: %s given twice", ErrSyntax, ErrParam, ll.line, xname, a.name)
		}
		v, err := evalExpr(a.expr, caller, 0)
		if err != nil {
			return fmt.Errorf("line %d: %s %s: %w", ll.line, xname, a.name, err)
		}
		s.vals[a.name] = v
	}
	if err := x.emit(def.body, s, inst); err != nil {
		return fmt.Errorf("%s (.subckt %s): %w", path, def.name, err)
	}
	return nil
}
