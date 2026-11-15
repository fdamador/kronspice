// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Mutual is the mutual inductance M = k·√(L1·L2) between two inductors
// (SPICE K). It adds the flux each inductor's current links through the
// other to their branch equations:
//
//	v1 = L1·di1/dt + M·di2/dt,   v2 = M·di1/dt + L2·di2/dt
//
// so in AC the branch rows gain −jωM off the diagonal, and in a transient
// step the integration formula of each inductor is applied to its total
// flux L1·i1 + M·i2 (backward Euler: −(M/h)·(i2 − i2,prev); trapezoidal:
// −(2M/h)·(i2 − i2,prev), the previous total voltage being already in the
// inductor's own stamp).
//
// Maxwell derivation (CLAUDE.md §3c): the flux of 𝔅 = V.∇𝔄 from one coil
// linking the other — the coupling row of the inductor entry. The element
// stores M·i1·i2 (part of ½LI² for the pair) and has no exterior port: its
// power is already in the two inductors' v·i.
type Mutual struct {
	Base
	k, m           float64
	names          [2]string
	k1, k2         int // the inductors' branch unknowns
	i1Prev, i2Prev float64
	i1, i2         float64 // currents at the most recent Operating

	ece eceTerms // hypothetical spin-connection terms (ece.go)
	rm  float64  // mutual torsion s·c·ω₀·M [Ω]
	rs  float64  // skew mutual torsion s·c·ω_skew·M [Ω]: v1 += rs·i2, v2 −= rs·i1
}

// NewMutual builds a K element; l1, l2 are the coupled inductances [H].
func NewMutual(e circuit.Element, l1, l2 float64) *Mutual {
	m := &Mutual{Base: newBase(e), k: e.Value, m: e.Value * math.Sqrt(l1*l2), names: e.Couples, k1: -1, k2: -1}
	_ = m.SetIntrinsic(exchange.Electrical, exchange.Both,
		"couples two inductors inside the network; stores M·i1·i2")
	return m
}

// M returns the mutual inductance [H].
func (m *Mutual) M() float64 { return m.m }

// Setup allocates nothing: the element acts on the inductors' branches.
func (m *Mutual) Setup(*SetupContext) error { return nil }

// Link finds the two inductors' branch unknowns.
func (m *Mutual) Link(ctx *SetupContext) error {
	for i, n := range m.names {
		b, ok := ctx.BranchOf(n)
		if !ok {
			return fmt.Errorf("%w: %s: inductor %s not found", ErrUnsupported, m.name, n)
		}
		if i == 0 {
			m.k1 = b
		} else {
			m.k2 = b
		}
	}
	return nil
}

// Load adds the coupling terms of a transient step (nothing at DC, where
// inductors are shorts).
func (m *Mutual) Load(ctx *LoadContext) error {
	if m.rs != 0 && ctx.Mode != ModeUICInit {
		// Skew mutual torsion (circuit.ECETerms.W0Skew): −r_s·i2 in coil 1's
		// row and +r_s·i1 in coil 2's, with the trapezoidal history given
		// back as for r_m.
		ctx.M.Add(m.k1, m.k2, -m.rs)
		ctx.M.Add(m.k2, m.k1, m.rs)
		if ctx.Mode == ModeTran && ctx.Method == Trapezoidal {
			ctx.AddRHS(m.k1, m.rs*m.i2Prev)
			ctx.AddRHS(m.k2, -m.rs*m.i1Prev)
		}
	}
	if m.rm != 0 && ctx.Mode != ModeUICInit {
		// Mutual torsion: −r_m·i_other in each coil's branch row, at DC
		// too. Each coil's trapezoidal history −v_prev holds its body
		// voltage less its own torsion, which includes r_m·i_other,prev;
		// that resistive part is not flux, so it is given back here.
		ctx.M.Add(m.k1, m.k2, -m.rm)
		ctx.M.Add(m.k2, m.k1, -m.rm)
		if ctx.Mode == ModeTran && ctx.Method == Trapezoidal {
			ctx.AddRHS(m.k1, m.rm*m.i2Prev)
			ctx.AddRHS(m.k2, m.rm*m.i1Prev)
		}
	}
	switch ctx.Mode {
	case ModeOP, ModeTranOP, ModeUICInit:
		return nil
	}
	f := m.m / ctx.H
	if ctx.Method == Trapezoidal {
		f = 2 * m.m / ctx.H
	}
	ctx.M.Add(m.k1, m.k2, -f)
	ctx.M.Add(m.k2, m.k1, -f)
	ctx.AddRHS(m.k1, -f*m.i2Prev)
	ctx.AddRHS(m.k2, -f*m.i1Prev)
	return nil
}

// Operating records the two currents; the element itself has no voltage,
// current or power of its own.
func (m *Mutual) Operating(x []float64) OperatingPoint {
	m.i1, m.i2 = at(x, m.k1), at(x, m.k2)
	if m.rm != 0 {
		// The coils' storage rates include r_m·i_other·i; the element
		// takes it back out and books it as model-supplied (ece.go).
		p := -2 * m.rm * m.i1 * m.i2
		return OperatingPoint{StoragePower: p, ECEPower: p}
	}
	return OperatingPoint{}
}

// Accept commits the two currents.
func (m *Mutual) Accept(x []float64, _ float64, _ Mode) {
	m.i1Prev, m.i2Prev = at(x, m.k1), at(x, m.k2)
}

// StoredEnergy is M·i1·i2 at the currents of the most recent Operating.
func (m *Mutual) StoredEnergy(OperatingPoint) float64 { return m.m * m.i1 * m.i2 }

// Truncate: the step is limited by the inductors themselves.
func (m *Mutual) Truncate([]float64, *TruncContext) float64 { return math.Inf(1) }

// ExteriorPower: no exterior port.
func (m *Mutual) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// LoadAC stamps −jωM between the two branch rows.
func (m *Mutual) LoadAC(ctx *ACContext) error {
	z := complex(-m.rm, -ctx.Omega*m.m)
	ctx.M.Add(m.k1, m.k2, z-complex(m.rs, 0))
	ctx.M.Add(m.k2, m.k1, z+complex(m.rs, 0))
	return nil
}

// OperatingAC: no phasors of its own; with mutual torsion, the average
// power it delivers, −r_m·Re(I₁·I₂*).
func (m *Mutual) OperatingAC(x []complex128) ACPoint {
	if m.rm == 0 {
		return ACPoint{}
	}
	i1, i2 := x[m.k1], x[m.k2]
	return ACPoint{ECEPower: -2 * m.rm * avgPower(i1, i2)}
}
