// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package exchange

import "errors"

// Direction says which way energy may cross a port.
type Direction uint8

const (
	In Direction = iota
	Out
	Both
)

func (d Direction) String() string {
	switch d {
	case In:
		return "in"
	case Out:
		return "out"
	default:
		return "both"
	}
}

// PortStatus says whether a port's exchange is computed.
type PortStatus uint8

const (
	// NotModeled ports contribute exactly zero power (Kron: inactive junction pair).
	NotModeled PortStatus = iota
	// Modeled ports have a cited coupling model, or are computed intrinsically
	// by the solver (see ExternalPort.Intrinsic).
	Modeled
	// Hypothetical ports use a model outside established physics; their power
	// is kept in a separate ledger column and never counted as real inflow.
	Hypothetical
)

func (s PortStatus) String() string {
	switch s {
	case Modeled:
		return "modeled"
	case Hypothetical:
		return "hypothetical"
	default:
		return "not modeled"
	}
}

// ExchangeContext is what a coupling model sees when asked for power.
type ExchangeContext struct {
	Time float64 // simulation time [s]; 0 for operating-point analyses
}

// CouplingModel computes the power exchanged through one port.
type CouplingModel interface {
	// Exchange returns the power entering the element through this port from
	// the exterior (+) or leaving it to the exterior (−), in watts.
	Exchange(ctx *ExchangeContext) (powerW float64, err error)
	// Citation names the published source of the model. It must not be empty.
	Citation() string
}

// ExternalPort is one element's open path to one channel's exterior junction.
type ExternalPort struct {
	Channel   ChannelID
	Direction Direction
	Status    PortStatus
	// Intrinsic is true when the solver itself computes this port's power
	// (for example electrical supply by an independent source, or resistor
	// dissipation to the thermal exterior). Intrinsic ports need no Model.
	Intrinsic bool
	Model     CouplingModel
	// Note is a short human-readable description, e.g. "Joule heating".
	Note string
}

// Exchanger is implemented by every element.
type Exchanger interface {
	// Ports returns one port per registered channel, in channel order.
	Ports() []ExternalPort
}

// Errors returned when configuring ports.
var (
	ErrUnknownChannel = errors.New("exchange: unknown channel")
	ErrNoCitation     = errors.New("exchange: coupling model has no citation")
	ErrNilModel       = errors.New("exchange: modeled port needs a coupling model or must be intrinsic")
)
