// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package exchange defines the external energy channels and ports that every
// circuit element carries (CLAUDE.md §3.1, §3a.1).
//
// Every element is an open system. For each registered channel it has one
// external port: in Kron's terms, an apparent coil (an open mesh) from the
// element to that channel's exterior junction, carrying impressed quantities
// "that originate outside the network" (Kron, Tensor Analysis of Networks,
// p. 76; apparent coils p. 407). A port with no coupling model carries an
// impressed quantity of zero value, as Kron adds to every junction-pair "for
// the sake of completeness" (p. 408); such inactive junction-pairs may be
// ignored when setting up the impedance (p. 429), so unmodeled ports cost
// nothing in the inner loop.
//
// This package imports only the standard library so that every other package
// can depend on it.
package exchange

import "fmt"

// ChannelID identifies one external energy channel.
type ChannelID uint8

// Registered channels. Adding a channel here gives every element a port for it
// automatically; no device code changes are needed.
const (
	Electrical                 ChannelID = iota // circuit terminals; exterior = impressed sources
	Light                                       // photocurrent in, emitted light out
	Thermal                                     // heat exchanged with the surroundings
	Magnetic                                    // coupling to external B(t), radiated magnetic field
	ElectricField                               // coupling to external E(t), emitted field
	Mechanical                                  // shaft power, strain, magnetostriction
	Acoustic                                    // sound and shock in/out
	IonizingRadiation                           // dose and particle strikes
	Chemical                                    // batteries, fuel cells, electrolysis
	Neutrino                                    // energy deposited by neutrinos (established model, ≈0)
	Gravitational                               // frame effects only (time dilation, redshift); no energy
	HypotheticalMagneticCharge                  // dual-pair force/power; hypothetical mode only
	Future1                                     // reserved for channels not yet named
	Future2
	Future3
	Future4

	numChannels // keep last
)

// NumChannels is the number of registered channels.
const NumChannels = int(numChannels)

// Kind says which variables a channel's open path carries.
type Kind uint8

const (
	// AcrossThrough channels have an across and a through variable whose
	// product is power (Kron/Branin network analogies).
	AcrossThrough Kind = iota
	// PowerOnly channels carry only a power balance (light, radiation, neutrinos).
	PowerOnly
	// FrameOnly channels change how signals are observed but carry no energy.
	FrameOnly
)

type channelInfo struct {
	name    string
	kind    Kind
	across  string // across variable and unit, if AcrossThrough
	through string // through variable and unit, if AcrossThrough
}

var channels = [numChannels]channelInfo{
	Electrical:                 {"electrical", AcrossThrough, "voltage [V]", "current [A]"},
	Light:                      {"light", PowerOnly, "", ""},
	Thermal:                    {"thermal", AcrossThrough, "temperature [K]", "heat flow [W]"},
	Magnetic:                   {"magnetic", AcrossThrough, "magnetomotive force [A]", "flux rate [V]"},
	ElectricField:              {"electric-field", AcrossThrough, "voltage [V]", "displacement current [A]"},
	Mechanical:                 {"mechanical", AcrossThrough, "velocity [m/s]", "force [N]"},
	Acoustic:                   {"acoustic", AcrossThrough, "pressure [Pa]", "volume velocity [m³/s]"},
	IonizingRadiation:          {"ionizing-radiation", PowerOnly, "", ""},
	Chemical:                   {"chemical", AcrossThrough, "chemical potential [V]", "reaction current [A]"},
	Neutrino:                   {"neutrino", PowerOnly, "", ""},
	Gravitational:              {"gravitational", FrameOnly, "", ""},
	HypotheticalMagneticCharge: {"hypothetical-magnetic-charge", PowerOnly, "", ""},
	Future1:                    {"future-1", PowerOnly, "", ""},
	Future2:                    {"future-2", PowerOnly, "", ""},
	Future3:                    {"future-3", PowerOnly, "", ""},
	Future4:                    {"future-4", PowerOnly, "", ""},
}

// String returns the channel's short name.
func (c ChannelID) String() string {
	if int(c) >= NumChannels {
		return fmt.Sprintf("channel(%d)", uint8(c))
	}
	return channels[c].name
}

// Kind returns which variables the channel carries.
func (c ChannelID) Kind() Kind {
	if int(c) >= NumChannels {
		return PowerOnly
	}
	return channels[c].kind
}

// Variables returns the across and through variable descriptions for an
// AcrossThrough channel, or empty strings otherwise.
func (c ChannelID) Variables() (across, through string) {
	if int(c) >= NumChannels {
		return "", ""
	}
	return channels[c].across, channels[c].through
}

// Valid reports whether c is a registered channel.
func (c ChannelID) Valid() bool { return int(c) < NumChannels }

// All returns every registered channel in order.
func All() []ChannelID {
	out := make([]ChannelID, NumChannels)
	for i := range out {
		out[i] = ChannelID(i)
	}
	return out
}

// ParseChannel returns the channel with the given short name (as String
// prints it, case-insensitive), and whether one exists.
func ParseChannel(name string) (ChannelID, bool) {
	for i := range channels {
		if equalFold(channels[i].name, name) {
			return ChannelID(i), true
		}
	}
	return 0, false
}

func equalFold(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		x, y := a[i], b[i]
		if 'A' <= x && x <= 'Z' {
			x += 'a' - 'A'
		}
		if 'A' <= y && y <= 'Z' {
			y += 'a' - 'A'
		}
		if x != y {
			return false
		}
	}
	return true
}
