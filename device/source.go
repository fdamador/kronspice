// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Independent sources are impressed sources in Kron's sense: they act across
// an open path to the electrical exterior junction, so the power they deliver
// enters the circuit from outside (CLAUDE.md §3a.1). Their electrical port is
// intrinsic; a source that absorbs power (e.g. a battery being charged) sends
// energy out to the exterior through the same port.

// sourceValue holds a source's DC value, optional waveform and AC phasor.
type sourceValue struct {
	dc    float64
	hasDC bool
	spec  circuit.Waveform
	w     wave
	ac    complex128 // AC phasor: magnitude ∠ phase

	override   float64 // DC value forced by a .DC sweep
	overridden bool
}

func newSourceValue(e circuit.Element) sourceValue {
	return sourceValue{dc: e.Value, hasDC: e.HasDC, spec: e.Wave,
		ac: cmplx.Rect(e.ACMag, e.ACPhase*math.Pi/180)}
}

func (s *sourceValue) setup(ctx *SetupContext) {
	s.w = resolveWave(s.spec, ctx.Step, ctx.Stop)
}

// at returns the source value for a solve. A stand-alone operating point uses
// the DC value when one was given; transients use the waveform (SPICE3).
func (s *sourceValue) at(ctx *LoadContext) float64 {
	v := s.raw(ctx)
	if ctx.SourceStepping && ctx.Mode != ModeTran {
		v *= ctx.SrcFact
	}
	return v
}

func (s *sourceValue) raw(ctx *LoadContext) float64 {
	if s.overridden && ctx.Mode == ModeOP {
		return s.override
	}
	if s.spec.Kind == circuit.WaveNone {
		return s.dc
	}
	if ctx.Mode == ModeOP && s.hasDC {
		return s.dc
	}
	t := ctx.Time
	if ctx.Mode != ModeTran {
		t = 0
	}
	return s.w.value(t)
}

func (s *sourceValue) nextBreakpoint(t float64) float64 {
	if s.spec.Kind == circuit.WaveNone {
		return inf
	}
	return s.w.nextBreakpoint(t)
}

// VoltageSource is an independent voltage source with branch current
// unknown, positive when current flows into n+ through the source to n−.
type VoltageSource struct {
	Base
	channel exchange.ChannelID // exterior channel supplying the energy
	src     sourceValue
	branch  int
	value   float64 // value used in the most recent Load
}

// NewVoltageSource builds a voltage source from a circuit element.
func NewVoltageSource(e circuit.Element) *VoltageSource {
	s := &VoltageSource{Base: newBase(e), src: newSourceValue(e), channel: e.Channel}
	s.channel = setSourcePorts(&s.Base, e)
	return s
}

// Setup maps nodes, allocates the branch-current unknown and resolves the waveform.
func (s *VoltageSource) Setup(ctx *SetupContext) error {
	s.setupNodes(ctx)
	s.branch = ctx.AllocBranch()
	ctx.RegisterBranch(s.name, s.branch)
	s.src.setup(ctx)
	return nil
}

// Load stamps the MNA voltage-source pattern (Ho, Ruehli & Brennan 1975).
func (s *VoltageSource) Load(ctx *LoadContext) error {
	a, b, k := s.u[0], s.u[1], s.branch
	s.value = s.src.at(ctx)
	ctx.M.Add(a, k, 1)
	ctx.M.Add(b, k, -1)
	ctx.M.Add(k, a, 1)
	ctx.M.Add(k, b, -1)
	ctx.AddRHS(k, s.value)
	return nil
}

// Operating returns voltage, branch current and absorbed power.
func (s *VoltageSource) Operating(x []float64) OperatingPoint {
	v := s.voltage(x)
	i := x[s.branch]
	return OperatingPoint{Voltage: v, Current: i, PowerAbsorbed: v * i}
}

// NextBreakpoint returns the waveform's next corner after t.
func (s *VoltageSource) NextBreakpoint(t float64) float64 { return s.src.nextBreakpoint(t) }

// ExteriorPower: power delivered by the source enters from the exterior.
func (s *VoltageSource) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == s.channel {
		return -op.PowerAbsorbed
	}
	return 0
}

// CurrentSource is an independent current source. Its current flows from
// n+ through the source to n− (SPICE convention).
type CurrentSource struct {
	Base
	channel exchange.ChannelID // exterior channel supplying the energy
	src     sourceValue
	value   float64
}

// NewCurrentSource builds a current source from a circuit element.
func NewCurrentSource(e circuit.Element) *CurrentSource {
	s := &CurrentSource{Base: newBase(e), src: newSourceValue(e), channel: e.Channel}
	s.channel = setSourcePorts(&s.Base, e)
	return s
}

// Setup maps the source's nodes and resolves the waveform.
func (s *CurrentSource) Setup(ctx *SetupContext) error {
	s.setupNodes(ctx)
	s.src.setup(ctx)
	return nil
}

// Load stamps the current into the right-hand side.
func (s *CurrentSource) Load(ctx *LoadContext) error {
	s.value = s.src.at(ctx)
	ctx.AddRHS(s.u[0], -s.value)
	ctx.AddRHS(s.u[1], s.value)
	return nil
}

// Operating returns voltage, current and absorbed power.
func (s *CurrentSource) Operating(x []float64) OperatingPoint {
	v := s.voltage(x)
	return OperatingPoint{Voltage: v, Current: s.value, PowerAbsorbed: v * s.value}
}

// NextBreakpoint returns the waveform's next corner after t.
func (s *CurrentSource) NextBreakpoint(t float64) float64 { return s.src.nextBreakpoint(t) }

// ExteriorPower: power delivered by the source enters from the exterior.
func (s *CurrentSource) ExteriorPower(c exchange.ChannelID, op OperatingPoint) float64 {
	if c == s.channel {
		return -op.PowerAbsorbed
	}
	return 0
}

// SetDCOverride forces the source's DC value (for .DC sweeps).
func (s *VoltageSource) SetDCOverride(v float64) { s.src.override, s.src.overridden = v, true }

// ClearDCOverride restores the source's own DC value.
func (s *VoltageSource) ClearDCOverride() { s.src.overridden = false }

// LoadAC stamps the source with its AC phasor (zero if it has no AC value).
func (s *VoltageSource) LoadAC(ctx *ACContext) error {
	a, b, k := s.u[0], s.u[1], s.branch
	ctx.M.Add(a, k, 1)
	ctx.M.Add(b, k, -1)
	ctx.M.Add(k, a, 1)
	ctx.M.Add(k, b, -1)
	ctx.AddRHS(k, s.src.ac)
	return nil
}

// OperatingAC returns the source's phasors.
func (s *VoltageSource) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: s.voltageAC(x), Current: x[s.branch]}
}

// SetDCOverride forces the source's DC value (for .DC sweeps).
func (s *CurrentSource) SetDCOverride(v float64) { s.src.override, s.src.overridden = v, true }

// ClearDCOverride restores the source's own DC value.
func (s *CurrentSource) ClearDCOverride() { s.src.overridden = false }

// LoadAC stamps the source's AC phasor current.
func (s *CurrentSource) LoadAC(ctx *ACContext) error {
	ctx.AddRHS(s.u[0], -s.src.ac)
	ctx.AddRHS(s.u[1], s.src.ac)
	return nil
}

// OperatingAC returns the source's phasors.
func (s *CurrentSource) OperatingAC(x []complex128) ACPoint {
	return ACPoint{Voltage: s.voltageAC(x), Current: s.src.ac}
}

// setSourcePorts marks the channel an independent source exchanges energy
// through: the electrical exterior by default, or the physical channel the
// circuit element declares (Seebeck EMF → thermal, galvanic cell →
// chemical). In the latter case the terminals connect inside the network.
func setSourcePorts(b *Base, e circuit.Element) exchange.ChannelID {
	if e.Channel == exchange.Electrical {
		_ = b.SetIntrinsic(exchange.Electrical, exchange.Both, "impressed source: exchanges energy with the exterior")
		return exchange.Electrical
	}
	_ = b.SetIntrinsic(exchange.Electrical, exchange.Both, "terminals connect inside the network")
	note := e.ChannelNote
	if note == "" {
		note = "source energy from the " + e.Channel.String() + " exterior"
	}
	_ = b.SetIntrinsic(e.Channel, exchange.Both, note)
	return e.Channel
}
