// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package netlist

import (
	"errors"
	"fmt"
	"math"
	"strings"
)

// ErrParam is returned (wrapped with ErrSyntax) for an undefined, circular
// or non-finite .param, or a malformed expression.
var ErrParam = errors.New("netlist: parameter error")

// maxExprDepth bounds the nesting of an expression and of parameters that
// refer to parameters, so hostile input cannot exhaust the stack.
const maxExprDepth = 200

// resolver looks up a parameter by lower-cased name.
type resolver interface {
	value(name string, depth int) (float64, error)
}

// EvalExpr evaluates an expression with no parameters defined: numbers with
// engineering suffixes, + − * / ^ (or **), comparisons, && || !, the
// conditional c ? a : b, the constant pi and the functions listed at
// exprFuncs.
func EvalExpr(s string) (float64, error) { return evalExpr(s, nil, 0) }

// evalExpr evaluates s, resolving identifiers through r (nil: none).
//
// Grammar (precedence low to high, as in ngspice's .param expressions):
//
//	cond  := or [ "?" cond ":" cond ]
//	or    := and { "||" and }
//	and   := cmp { "&&" cmp }
//	cmp   := sum [ ("=="|"!="|"<"|"<="|">"|">=") sum ]
//	sum   := prod { ("+"|"-") prod }
//	prod  := unary { ("*"|"/"|"%") unary }
//	unary := ("-"|"+"|"!") unary | pow
//	pow   := atom [ ("^"|"**") unary ]          (right-associative)
//	atom  := number | name | name "(" args ")" | "(" cond ")"
//
// A number is a SPICE value with optional suffix ("4.7k", "10meg", "2u").
// Comparisons and logical operators give 1 or 0. The result must be finite.
func evalExpr(s string, r resolver, depth int) (float64, error) {
	if depth > maxExprDepth {
		return 0, fmt.Errorf("%w: %w: expressions nested too deeply", ErrSyntax, ErrParam)
	}
	p := &exprParser{s: s, r: r, depth: depth}
	v, err := p.cond()
	if err != nil {
		return 0, err
	}
	p.space()
	if p.i < len(p.s) {
		return 0, p.fail("unexpected %q", p.s[p.i:])
	}
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, fmt.Errorf("%w: %w: %q is not finite", ErrSyntax, ErrParam, s)
	}
	return v, nil
}

type exprParser struct {
	s     string
	i     int
	r     resolver
	depth int
	nest  int
}

func (p *exprParser) fail(format string, a ...any) error {
	return fmt.Errorf("%w: %w: expression %q: %s", ErrSyntax, ErrParam, p.s, fmt.Sprintf(format, a...))
}

func (p *exprParser) space() {
	for p.i < len(p.s) && (p.s[p.i] == ' ' || p.s[p.i] == '\t') {
		p.i++
	}
}

// eat consumes tok if it comes next.
func (p *exprParser) eat(tok string) bool {
	p.space()
	if strings.HasPrefix(p.s[p.i:], tok) {
		p.i += len(tok)
		return true
	}
	return false
}

func (p *exprParser) enter() error {
	p.nest++
	if p.nest+p.depth > maxExprDepth {
		return p.fail("nested too deeply")
	}
	return nil
}

func b2f(b bool) float64 {
	if b {
		return 1
	}
	return 0
}

func (p *exprParser) cond() (float64, error) {
	if err := p.enter(); err != nil {
		return 0, err
	}
	defer func() { p.nest-- }()
	c, err := p.or()
	if err != nil {
		return 0, err
	}
	if !p.eat("?") {
		return c, nil
	}
	a, err := p.cond()
	if err != nil {
		return 0, err
	}
	if !p.eat(":") {
		return 0, p.fail("missing ':' in c ? a : b")
	}
	b, err := p.cond()
	if err != nil {
		return 0, err
	}
	if c != 0 {
		return a, nil
	}
	return b, nil
}

func (p *exprParser) or() (float64, error) {
	v, err := p.and()
	for err == nil && p.eat("||") {
		var w float64
		w, err = p.and()
		v = b2f(v != 0 || w != 0)
	}
	return v, err
}

func (p *exprParser) and() (float64, error) {
	v, err := p.cmp()
	for err == nil && p.eat("&&") {
		var w float64
		w, err = p.cmp()
		v = b2f(v != 0 && w != 0)
	}
	return v, err
}

func (p *exprParser) cmp() (float64, error) {
	v, err := p.sum()
	if err != nil {
		return 0, err
	}
	for _, op := range []string{"==", "!=", "<=", ">=", "<", ">"} {
		if p.eat(op) {
			w, err := p.sum()
			if err != nil {
				return 0, err
			}
			switch op {
			case "==":
				return b2f(v == w), nil
			case "!=":
				return b2f(v != w), nil
			case "<=":
				return b2f(v <= w), nil
			case ">=":
				return b2f(v >= w), nil
			case "<":
				return b2f(v < w), nil
			default:
				return b2f(v > w), nil
			}
		}
	}
	return v, nil
}

func (p *exprParser) sum() (float64, error) {
	v, err := p.prod()
	for err == nil {
		var w float64
		switch {
		case p.eat("+"):
			w, err = p.prod()
			v += w
		case p.eat("-"):
			w, err = p.prod()
			v -= w
		default:
			return v, nil
		}
	}
	return v, err
}

func (p *exprParser) prod() (float64, error) {
	v, err := p.unary()
	for err == nil {
		p.space()
		if strings.HasPrefix(p.s[p.i:], "**") {
			return v, nil // handled by pow
		}
		var w float64
		switch {
		case p.eat("*"):
			w, err = p.unary()
			v *= w
		case p.eat("/"):
			w, err = p.unary()
			if err == nil && w == 0 {
				return 0, p.fail("division by zero")
			}
			v /= w
		case p.eat("%"):
			w, err = p.unary()
			if err == nil && w == 0 {
				return 0, p.fail("modulo by zero")
			}
			v = math.Mod(v, w)
		default:
			return v, nil
		}
	}
	return v, err
}

func (p *exprParser) unary() (float64, error) {
	if err := p.enter(); err != nil {
		return 0, err
	}
	defer func() { p.nest-- }()
	switch {
	case p.eat("-"):
		v, err := p.unary()
		return -v, err
	case p.eat("+"):
		return p.unary()
	case p.peek("!") && !p.peek("!="):
		p.eat("!")
		v, err := p.unary()
		return b2f(v == 0), err
	}
	return p.pow()
}

func (p *exprParser) peek(tok string) bool {
	p.space()
	return strings.HasPrefix(p.s[p.i:], tok)
}

func (p *exprParser) pow() (float64, error) {
	v, err := p.atom()
	if err != nil {
		return 0, err
	}
	if p.eat("**") || p.eat("^") {
		w, err := p.unary()
		if err != nil {
			return 0, err
		}
		return math.Pow(v, w), nil
	}
	return v, nil
}

func isIdentStart(c byte) bool {
	return c == '_' || (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
}

func isIdentChar(c byte) bool { return isIdentStart(c) || (c >= '0' && c <= '9') }

func (p *exprParser) atom() (float64, error) {
	p.space()
	if p.i >= len(p.s) {
		return 0, p.fail("unexpected end")
	}
	c := p.s[p.i]
	switch {
	case c == '(':
		p.i++
		v, err := p.cond()
		if err != nil {
			return 0, err
		}
		if !p.eat(")") {
			return 0, p.fail("missing ')'")
		}
		return v, nil
	case (c >= '0' && c <= '9') || c == '.':
		n := numericPrefix(p.s[p.i:])
		if n == 0 {
			return 0, p.fail("bad number at %q", p.s[p.i:])
		}
		j := p.i + n
		for j < len(p.s) && isIdentStart(p.s[j]) {
			j++ // engineering suffix and unit letters, as in ParseValue
		}
		v, err := ParseValue(p.s[p.i:j])
		if err != nil {
			return 0, fmt.Errorf("%w: %w: expression %q: %w", ErrSyntax, ErrParam, p.s, err)
		}
		p.i = j
		return v, nil
	case isIdentStart(c):
		j := p.i
		for j < len(p.s) && isIdentChar(p.s[j]) {
			j++
		}
		name := strings.ToLower(p.s[p.i:j])
		p.i = j
		if p.eat("(") {
			return p.call(name)
		}
		if p.r != nil {
			v, err := p.r.value(name, p.depth+p.nest+1)
			if err == nil || !errors.Is(err, errUndefined) {
				return v, err
			}
		}
		if name == "pi" {
			return math.Pi, nil
		}
		return 0, fmt.Errorf("%w: %w: undefined parameter %q", ErrSyntax, ErrParam, name)
	}
	return 0, p.fail("unexpected %q", p.s[p.i:])
}

// errUndefined marks a name no scope defines (internal).
var errUndefined = errors.New("undefined")

// exprFuncs are the functions an expression may call, with their arity
// (ngspice manual, "Parameters, functions, expressions").
var exprFuncs = map[string]struct {
	n int
	f func(a []float64) float64
}{
	"sqrt":  {1, func(a []float64) float64 { return math.Sqrt(a[0]) }},
	"exp":   {1, func(a []float64) float64 { return math.Exp(a[0]) }},
	"ln":    {1, func(a []float64) float64 { return math.Log(a[0]) }},
	"log":   {1, func(a []float64) float64 { return math.Log(a[0]) }}, // natural log, as in ngspice
	"log10": {1, func(a []float64) float64 { return math.Log10(a[0]) }},
	"sin":   {1, func(a []float64) float64 { return math.Sin(a[0]) }},
	"cos":   {1, func(a []float64) float64 { return math.Cos(a[0]) }},
	"tan":   {1, func(a []float64) float64 { return math.Tan(a[0]) }},
	"asin":  {1, func(a []float64) float64 { return math.Asin(a[0]) }},
	"acos":  {1, func(a []float64) float64 { return math.Acos(a[0]) }},
	"atan":  {1, func(a []float64) float64 { return math.Atan(a[0]) }},
	"atan2": {2, func(a []float64) float64 { return math.Atan2(a[0], a[1]) }},
	"sinh":  {1, func(a []float64) float64 { return math.Sinh(a[0]) }},
	"cosh":  {1, func(a []float64) float64 { return math.Cosh(a[0]) }},
	"tanh":  {1, func(a []float64) float64 { return math.Tanh(a[0]) }},
	"abs":   {1, func(a []float64) float64 { return math.Abs(a[0]) }},
	"floor": {1, func(a []float64) float64 { return math.Floor(a[0]) }},
	"ceil":  {1, func(a []float64) float64 { return math.Ceil(a[0]) }},
	"int":   {1, func(a []float64) float64 { return math.Trunc(a[0]) }},
	"nint":  {1, func(a []float64) float64 { return math.RoundToEven(a[0]) }},
	"sgn": {1, func(a []float64) float64 {
		if a[0] == 0 {
			return 0
		}
		return math.Copysign(1, a[0])
	}},
	"min": {2, func(a []float64) float64 { return math.Min(a[0], a[1]) }},
	"max": {2, func(a []float64) float64 { return math.Max(a[0], a[1]) }},
	"pow": {2, func(a []float64) float64 { return math.Pow(a[0], a[1]) }},
	"pwr": {2, func(a []float64) float64 { return math.Copysign(math.Pow(math.Abs(a[0]), a[1]), a[0]) }},
}

func (p *exprParser) call(name string) (float64, error) {
	fn, ok := exprFuncs[name]
	if !ok {
		return 0, p.fail("unknown function %s()", name)
	}
	var args []float64
	if !p.eat(")") {
		for {
			v, err := p.cond()
			if err != nil {
				return 0, err
			}
			args = append(args, v)
			if p.eat(")") {
				break
			}
			if !p.eat(",") {
				return 0, p.fail("missing ',' or ')' in %s()", name)
			}
		}
	}
	if len(args) != fn.n {
		return 0, p.fail("%s() takes %d argument(s), got %d", name, fn.n, len(args))
	}
	return fn.f(args), nil
}
