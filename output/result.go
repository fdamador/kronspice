// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package output holds simulation results returned to callers, and writers
// for reports and CSV files.
package output

import (
	"math"
	"strings"

	"github.com/fdamador/kronspice/exchange"
)

// ElementResult is one element at the operating point.
type ElementResult struct {
	Name          string
	Kind          string
	Nodes         [2]string
	Voltage       float64 // V(n+) − V(n−) [V]
	Current       float64 // from n+ to n− through the element [A]
	PowerAbsorbed float64 // [W]; negative when supplying
	// Stores is true for elements that store energy inside the circuit
	// (capacitors, inductors, junction charge); StoragePower is the rate.
	Stores       bool
	StoragePower float64 // part of PowerAbsorbed going into storage [W]
	// ModelSupplied is the power an idealized active element (negative
	// resistance) delivers to the circuit [W]; 0 for ordinary elements.
	ModelSupplied float64
	StoredEnergy  float64 // [J]
	Ports         []exchange.ExternalPort
	// ExteriorPower[c] is the power through channel c's port: + into the
	// circuit from the exterior, − out to it. Exactly 0 for unmodeled ports.
	ExteriorPower [exchange.NumChannels]float64
}

// LedgerRow is one element × channel entry of the energy ledger.
type LedgerRow struct {
	Element   string
	Channel   exchange.ChannelID
	Status    exchange.PortStatus
	Intrinsic bool
	Direction exchange.Direction
	PowerW    float64 // + into the circuit from the exterior, − out
	Note      string
}

// SuppliedRow is one idealized active element's model-supplied power.
type SuppliedRow struct {
	Element string
	PowerW  float64
}

// Ledger is the open-system energy ledger (CLAUDE.md §3c, §3.1).
// All totals are rates in watts for an operating point.
type Ledger struct {
	Rows []LedgerRow

	Inflow      float64 // ≥ 0: power entering from the exterior
	Outflow     float64 // ≥ 0: power leaving to the exterior other than heat
	Dissipation float64 // ≥ 0: heat leaving through thermal ports
	StorageRate float64 // d(stored energy)/dt; 0 at an operating point
	// ModelSupplied is power delivered by idealized active elements
	// (negative resistances, later controlled sources). It balances the
	// ledger but is never inflow from the exterior: a real device would take
	// it from a supply the model leaves out.
	ModelSupplied float64
	Supplied      []SuppliedRow // per element
	Hypothetical  float64       // hypothetical-mode power, never counted as real
	// Residual = Inflow + ModelSupplied − Outflow − Dissipation −
	// StorageRate. Must be ≈ 0.
	Residual float64
	// Scale is the sum of |row power| used for the relative residual check.
	Scale float64
	// Roundoff is the residual a correct solution can show from
	// floating-point round-off alone [W], estimated from the solve as
	// 64·ε·Σ_n |x_n|·(Σ_j |A_nj·x_j| + |b_n|): each row's cancellation error
	// times its unknown (Tellegen: residual = Σ x_n·r_n). It matters only
	// when the reported powers are tiny next to the currents that cancel
	// inside the circuit, for example picoamperes through a 1 Ω emitter
	// resistance at 10 V; a residual below it is not evidence of an error.
	Roundoff float64
}

// LedgerAbsFloor is the absolute residual below which a ledger balances
// regardless of its scale [W]. When every power in a circuit is itself at
// round-off level (for example a diode circuit with all sources at 0 V,
// where node voltages are ~1e-30 V), the relative residual is round-off
// divided by round-off and means nothing. 1e-21 W is far below any power the
// engine reports as a result, so this cannot hide a real imbalance.
const LedgerAbsFloor = 1e-21

// RelativeResidual returns |Residual| / Scale, or 0 when nothing flows or
// the residual is below LedgerAbsFloor or the round-off estimate.
func (l *Ledger) RelativeResidual() float64 {
	if l.Scale == 0 || math.Abs(l.Residual) <= math.Max(LedgerAbsFloor, l.Roundoff) {
		return 0
	}
	r := l.Residual / l.Scale
	if r < 0 {
		return -r
	}
	return r
}

// Result is the outcome of one analysis.
type Result struct {
	Title        string
	Analysis     string // "op"
	NodeNames    []string
	NodeVoltages []float64 // aligned with NodeNames; NodeVoltages[0] = 0 (ground)
	Elements     []ElementResult
	Ledger       *Ledger
}

// Voltage returns the voltage of a named node.
func (r *Result) Voltage(node string) (float64, bool) {
	key := strings.ToLower(node)
	if key == "gnd" {
		key = "0"
	}
	for i, n := range r.NodeNames {
		if n == key {
			return r.NodeVoltages[i], true
		}
	}
	return 0, false
}

// Element returns the named element's result.
func (r *Result) Element(name string) (*ElementResult, bool) {
	for i := range r.Elements {
		if strings.EqualFold(r.Elements[i].Name, name) {
			return &r.Elements[i], true
		}
	}
	return nil, false
}
