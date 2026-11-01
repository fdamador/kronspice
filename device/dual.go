// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import "math"

// dual is a forward-mode automatic-differentiation number: a value and its
// partial derivatives with respect to three independent variables (for the
// MOSFET drain current: vgs, vds, vbs). The MOSFET level 2 and 3 currents
// are written once in terms of dual, so their conductances are the exact
// derivatives of the current that is stamped — the SPICE sources carry
// long hand-written derivative chains, several of them approximate
// (for example level 2's dudvds = 0), which a Newton iteration tolerates
// but a small-signal analysis does not. dual is a small value type: no
// allocation.
type dual struct {
	v float64
	d [3]float64
}

// cst is a constant (zero derivatives).
func cst(v float64) dual { return dual{v: v} }

// variable returns the k-th independent variable with value v.
func variable(v float64, k int) dual {
	x := dual{v: v}
	x.d[k] = 1
	return x
}

func (a dual) add(b dual) dual {
	return dual{a.v + b.v, [3]float64{a.d[0] + b.d[0], a.d[1] + b.d[1], a.d[2] + b.d[2]}}
}

func (a dual) sub(b dual) dual {
	return dual{a.v - b.v, [3]float64{a.d[0] - b.d[0], a.d[1] - b.d[1], a.d[2] - b.d[2]}}
}

func (a dual) addc(c float64) dual { a.v += c; return a }

func (a dual) neg() dual { return dual{-a.v, [3]float64{-a.d[0], -a.d[1], -a.d[2]}} }

func (a dual) scale(c float64) dual {
	return dual{a.v * c, [3]float64{a.d[0] * c, a.d[1] * c, a.d[2] * c}}
}

func (a dual) mul(b dual) dual {
	return dual{a.v * b.v, [3]float64{
		a.d[0]*b.v + a.v*b.d[0], a.d[1]*b.v + a.v*b.d[1], a.d[2]*b.v + a.v*b.d[2]}}
}

func (a dual) div(b dual) dual {
	q := a.v / b.v
	return dual{q, [3]float64{
		(a.d[0] - q*b.d[0]) / b.v, (a.d[1] - q*b.d[1]) / b.v, (a.d[2] - q*b.d[2]) / b.v}}
}

// chain applies a scalar function with value f and derivative df at a.v.
func (a dual) chain(f, df float64) dual {
	return dual{f, [3]float64{df * a.d[0], df * a.d[1], df * a.d[2]}}
}

func (a dual) sqrt() dual {
	s := math.Sqrt(a.v)
	return a.chain(s, 0.5/s)
}

func (a dual) exp() dual {
	e := math.Exp(a.v)
	return a.chain(e, e)
}

func (a dual) log() dual { return a.chain(math.Log(a.v), 1/a.v) }

func (a dual) cos() dual { return a.chain(math.Cos(a.v), -math.Sin(a.v)) }

func (a dual) atan() dual { return a.chain(math.Atan(a.v), 1/(1+a.v*a.v)) }

func (a dual) abs() dual {
	if a.v < 0 {
		return a.neg()
	}
	return a
}

// powc is a raised to the constant power c (a > 0).
func (a dual) powc(c float64) dual {
	p := math.Pow(a.v, c)
	return a.chain(p, c*p/a.v)
}

// dmax and dmin pick one operand (ties: the first), derivatives included.
func dmax(a, b dual) dual {
	if b.v > a.v {
		return b
	}
	return a
}

func dmin(a, b dual) dual {
	if b.v < a.v {
		return b
	}
	return a
}
