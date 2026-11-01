// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// ErrLedgerResidual means energy does not balance: inflow minus outflow,
// dissipation and storage change is larger than the tolerance. It signals a
// numerical or modeling bug, never a result (§3c, §1 long-term vision rule 4).
var ErrLedgerResidual = errors.New("em: energy ledger does not balance")

// DefaultTolerance is the relative ledger tolerance for a solved time point.
const DefaultTolerance = 1e-9

// Accumulator sums ledger contributions with one set of classification rules,
// used both for full ledgers and for the per-time-point transient ledger.
type Accumulator struct {
	Inflow, Outflow, Dissipation, StorageRate, ModelSupplied, Hypothetical float64
	Scale                                                                  float64 // Σ|contribution|
}

// AddPort adds the power through one element's port on one channel:
// + into the circuit from the exterior, − out to it.
func (a *Accumulator) AddPort(c exchange.ChannelID, status exchange.PortStatus, pw float64) {
	if pw == 0 {
		return
	}
	if status == exchange.Hypothetical {
		a.Hypothetical += pw
		return
	}
	a.Scale += math.Abs(pw)
	switch {
	case c == exchange.Thermal && pw < 0:
		a.Dissipation += -pw
	case pw > 0:
		a.Inflow += pw
	default:
		a.Outflow += -pw
	}
}

// AddStorage adds the power absorbed by an energy-storing element (dW/dt).
func (a *Accumulator) AddStorage(p float64) {
	a.StorageRate += p
	a.Scale += math.Abs(p)
}

// AddModelSupplied adds power delivered by an idealized active element
// (device.ModelSupplier). It balances the ledger but is kept in its own
// column, never counted as inflow from the exterior.
func (a *Accumulator) AddModelSupplied(p float64) {
	if p == 0 {
		return
	}
	a.ModelSupplied += p
	a.Scale += math.Abs(p)
}

// Residual is inflow + model-supplied − outflow − dissipation − storage rate.
func (a *Accumulator) Residual() float64 {
	return a.Inflow + a.ModelSupplied - a.Outflow - a.Dissipation - a.StorageRate
}

// RelativeResidual is |Residual| / Scale (0 when nothing flows).
func (a *Accumulator) RelativeResidual() float64 {
	if a.Scale == 0 || math.Abs(a.Residual()) <= output.LedgerAbsFloor {
		return 0
	}
	return math.Abs(a.Residual()) / a.Scale
}

// BuildLedger builds the element × channel ledger from a result.
//
// Every element contributes one row per registered channel, including
// unmodeled channels (power exactly 0, status "not modeled"), so reports
// always show which channels were considered. The ledger is the power through
// each element's open path to its channel's exterior junction (§3a.1), plus
// the rate of change of energy stored inside the circuit. Because sources are
// counted as branches, Tellegen's theorem makes the residual zero for a
// correct solution (Tellegen 1952; Poynting's theorem in the lumped limit).
func BuildLedger(r *output.Result) *output.Ledger {
	var acc Accumulator
	l := &output.Ledger{}
	for _, e := range r.Elements {
		for _, p := range e.Ports {
			pw := e.ExteriorPower[p.Channel]
			l.Rows = append(l.Rows, output.LedgerRow{
				Element:   e.Name,
				Channel:   p.Channel,
				Status:    p.Status,
				Intrinsic: p.Intrinsic,
				Direction: p.Direction,
				PowerW:    pw,
				Note:      p.Note,
			})
			acc.AddPort(p.Channel, p.Status, pw)
		}
		if e.Stores {
			acc.AddStorage(e.StoragePower)
		}
		if e.ModelSupplied != 0 {
			l.Supplied = append(l.Supplied, output.SuppliedRow{Element: e.Name, PowerW: e.ModelSupplied})
			acc.AddModelSupplied(e.ModelSupplied)
		}
	}
	l.Inflow, l.Outflow, l.Dissipation = acc.Inflow, acc.Outflow, acc.Dissipation
	l.StorageRate, l.ModelSupplied, l.Hypothetical = acc.StorageRate, acc.ModelSupplied, acc.Hypothetical
	l.Residual, l.Scale = acc.Residual(), acc.Scale
	return l
}

// CheckLedger returns ErrLedgerResidual when the relative residual exceeds tol.
func CheckLedger(l *output.Ledger, tol float64) error {
	if rr := l.RelativeResidual(); rr > tol {
		return fmt.Errorf("%w: residual %.3g W (relative %.3g > %.3g)", ErrLedgerResidual, l.Residual, rr, tol)
	}
	return nil
}
