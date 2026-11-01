// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package exchange

import "fmt"

// Base holds the full port set. Every device embeds it (via device.Base), so
// every element exposes every channel from the moment it is created.
type Base struct {
	ports  [numChannels]ExternalPort
	active []ChannelID // ports that have a model (skipped otherwise: zero cost)
	inited bool
}

func (b *Base) init() {
	if b.inited {
		return
	}
	for i := range b.ports {
		b.ports[i] = ExternalPort{Channel: ChannelID(i), Direction: Both, Status: NotModeled}
	}
	b.inited = true
}

// Ports returns a copy of all ports, one per registered channel.
func (b *Base) Ports() []ExternalPort {
	b.init()
	out := make([]ExternalPort, NumChannels)
	copy(out, b.ports[:])
	return out
}

// Port returns the port for one channel.
func (b *Base) Port(c ChannelID) (ExternalPort, error) {
	b.init()
	if !c.Valid() {
		return ExternalPort{}, fmt.Errorf("%w: %d", ErrUnknownChannel, c)
	}
	return b.ports[c], nil
}

// SetIntrinsic marks a port as computed by the solver itself.
func (b *Base) SetIntrinsic(c ChannelID, d Direction, note string) error {
	b.init()
	if !c.Valid() {
		return fmt.Errorf("%w: %d", ErrUnknownChannel, c)
	}
	b.ports[c] = ExternalPort{Channel: c, Direction: d, Status: Modeled, Intrinsic: true, Note: note}
	b.refreshActive()
	return nil
}

// Attach connects a cited coupling model to a port. Hypothetical models must
// be attached with status Hypothetical.
func (b *Base) Attach(c ChannelID, d Direction, m CouplingModel, status PortStatus) error {
	b.init()
	if !c.Valid() {
		return fmt.Errorf("%w: %d", ErrUnknownChannel, c)
	}
	if m == nil {
		return ErrNilModel
	}
	if m.Citation() == "" {
		return ErrNoCitation
	}
	if status == NotModeled {
		status = Modeled
	}
	b.ports[c] = ExternalPort{Channel: c, Direction: d, Status: status, Model: m}
	b.refreshActive()
	return nil
}

// Detach returns a port to "not modeled".
func (b *Base) Detach(c ChannelID) error {
	b.init()
	if !c.Valid() {
		return fmt.Errorf("%w: %d", ErrUnknownChannel, c)
	}
	b.ports[c] = ExternalPort{Channel: c, Direction: Both, Status: NotModeled}
	b.refreshActive()
	return nil
}

// ActiveModels returns the channels whose ports have a coupling model. The
// inner loop visits only these, so unmodeled ports cost nothing.
func (b *Base) ActiveModels() []ChannelID {
	b.init()
	return b.active
}

func (b *Base) refreshActive() {
	b.active = b.active[:0]
	for i := range b.ports {
		if b.ports[i].Model != nil {
			b.active = append(b.active, ChannelID(i))
		}
	}
}
