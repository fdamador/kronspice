// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
)

// Hypothetical spin-connection terms on R, L, C and K (register entry H48;
// circuit.ECEForm gives the field relations and the circuit mapping).
//
// A two-terminal element's body sees the body voltage
//
//	v_b = α_a·φ_a − α_b·φ_b,   α = 1 + s·κ (scalar term, ω·φ),
//
// instead of v = φ_a − φ_b, and an inductor carries the torsion voltage
// r_t·i, r_t = s·c·ω₀·L, in series (ω₀·A term); coupled coils carry
// s·c·ω₀·M times the other coil's current. UFT 456's conductivities add
// −2cκ_Λ0·L to r_t (−2cκ_Λ0·M to the mutual term) and a shunt
// g_e = 2cκ₀·C across a capacitor's body. With every κ, ω₀, κ₀ and κ_Λ0
// zero the stamps are exactly Maxwell's.
//
// Energy: the terms are idealized, so the power they deliver is booked per
// element as model-supplied (ModelSupplier), never as inflow, with the
// label ModelSupplyNote (CLAUDE.md §3c ledger table, §1 hypothesis rules):
//
//	R:  heat v_b·i, model-supplied (v_b − v)·i
//	C:  storage v_b·i_C (charge C·v_b), shunt current g_e·v_b, terminal
//	    current i = i_C + g_e·v_b, model-supplied (v_b − v)·i − g_e·v_b²
//	L:  storage (v_b − r_t·i)·i, model-supplied (v_b − r_t·i − v)·i
//	K:  model-supplied −2·r_m·i₁·i₂ (r_m = s·c·ω₀·M), booked as a
//	    correction to the coils' storage rate; the skew placement
//	    (r_s = s·c·ω_skew·M, +r_s·i₂ on coil 1 and −r_s·i₁ on coil 2)
//	    moves r_s·i₁·i₂ from one coil's storage row to the other's and
//	    delivers nothing (Kron's workless F, Tensors for Circuits eq. 31.5)
//
// so every element satisfies v·i + model-supplied = heat + storage and the
// ledger closes in OP, AC and transient.
type eceTerms struct {
	circuit.ECETerms
	aA, aB float64 // α_a, α_b
	label  string
}

func newECE(t circuit.ECETerms) eceTerms {
	if !t.On {
		return eceTerms{aA: 1, aB: 1}
	}
	return eceTerms{ECETerms: t, aA: t.AlphaA(), aB: t.AlphaB(), label: t.Label()}
}

// ModelSupplyNoter devices label their model-supplied power for reports.
type ModelSupplyNoter interface {
	ModelSupplyNote() string
}

// body is v_b from solution x.
func (b *Base) body(x []float64, e *eceTerms) float64 {
	v := 0.0
	if b.u[0] >= 0 {
		v += e.aA * x[b.u[0]]
	}
	if b.u[1] >= 0 {
		v -= e.aB * x[b.u[1]]
	}
	return v
}

// bodyAC is the phasor v_b from AC solution x.
func (b *Base) bodyAC(x []complex128, e *eceTerms) complex128 {
	var v complex128
	if b.u[0] >= 0 {
		v += complex(e.aA, 0) * x[b.u[0]]
	}
	if b.u[1] >= 0 {
		v -= complex(e.aB, 0) * x[b.u[1]]
	}
	return v
}

// stampGB stamps a current i = g·v_b from node a to node b.
func stampGB(m Stamper, a, b int, g float64, e *eceTerms) {
	m.Add(a, a, g*e.aA)
	m.Add(a, b, -g*e.aB)
	m.Add(b, a, -g*e.aA)
	m.Add(b, b, g*e.aB)
}

// stampYB stamps a current I = y·V_b from node a to node b.
func stampYB(m CStamper, a, b int, y complex128, e *eceTerms) {
	m.Add(a, a, y*complex(e.aA, 0))
	m.Add(a, b, -y*complex(e.aB, 0))
	m.Add(b, a, -y*complex(e.aA, 0))
	m.Add(b, b, y*complex(e.aB, 0))
}

// avgPower is the average power ½·Re(v·i*) of two phasors.
func avgPower(v, i complex128) float64 { return 0.5 * real(v*cmplx.Conj(i)) }

// SetECE attaches resolved spin-connection terms (from circuit.ECEFor).
func (r *Resistor) SetECE(t circuit.ECETerms) { r.ece = newECE(t) }

// SetECE attaches resolved spin-connection terms (from circuit.ECEFor),
// including UFT 456's shunt g_e = 2cκ₀·C across the body.
func (c *Capacitor) SetECE(t circuit.ECETerms) {
	c.ece = newECE(t)
	c.ge = 0
	if t.On {
		c.ge = t.Shunt(c.C)
	}
}

// SetECE attaches resolved spin-connection terms (from circuit.ECEFor). The
// series term is torsion s·c·ω₀·L plus UFT 456's r_h = −2cκ_Λ0·L.
func (l *Inductor) SetECE(t circuit.ECETerms) {
	l.ece = newECE(t)
	l.rt = 0
	if t.On {
		l.rt = t.Torsion(l.L) + t.Homogeneous(l.L)
	}
}

// SetECE attaches resolved spin-connection terms (from circuit.ECEFor):
// mutual torsion s·c·ω₀·M plus UFT 456's −2cκ_Λ0·M.
func (m *Mutual) SetECE(t circuit.ECETerms) {
	m.ece = newECE(t)
	m.rm, m.rs = 0, 0
	if t.On {
		m.rm = t.Torsion(m.m) + t.Homogeneous(m.m)
		m.rs = t.TorsionSkew(m.m)
	}
}

// ModelSupplyNote labels the element's spin-connection power ("" if none).
func (r *Resistor) ModelSupplyNote() string { return r.ece.label }

// ModelSupplyNote labels the element's spin-connection power ("" if none).
func (c *Capacitor) ModelSupplyNote() string { return c.ece.label }

// ModelSupplyNote labels the element's spin-connection power ("" if none).
func (l *Inductor) ModelSupplyNote() string { return l.ece.label }

// ModelSupplyNote labels the element's spin-connection power ("" if none).
func (m *Mutual) ModelSupplyNote() string { return m.ece.label }

// Torsion returns the inductor's series torsion resistance s·c·ω₀·L [Ω].
func (l *Inductor) Torsion() float64 { return l.rt }

// Torsion returns the mutual torsion resistance s·c·ω₀·M [Ω].
func (m *Mutual) Torsion() float64 { return m.rm }

// ModelSuppliedPower is the power the spin-connection terms deliver.
func (c *Capacitor) ModelSuppliedPower(op OperatingPoint) float64 { return op.ECEPower }

// ModelSuppliedPower is the power the spin-connection terms deliver.
func (l *Inductor) ModelSuppliedPower(op OperatingPoint) float64 { return op.ECEPower }

// ModelSuppliedPower is the power the mutual torsion term delivers.
func (m *Mutual) ModelSuppliedPower(op OperatingPoint) float64 { return op.ECEPower }
