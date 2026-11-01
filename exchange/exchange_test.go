// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package exchange

import (
	"errors"
	"testing"
)

type fakeModel struct {
	p   float64
	cit string
}

func (f fakeModel) Exchange(*ExchangeContext) (float64, error) { return f.p, nil }
func (f fakeModel) Citation() string                           { return f.cit }

func TestAllChannelsPresent(t *testing.T) {
	var b Base
	ports := b.Ports()
	if len(ports) != NumChannels {
		t.Fatalf("got %d ports, want %d", len(ports), NumChannels)
	}
	for i, p := range ports {
		if p.Channel != ChannelID(i) {
			t.Errorf("port %d has channel %v", i, p.Channel)
		}
		if p.Status != NotModeled || p.Model != nil || p.Intrinsic {
			t.Errorf("new port %v should be not modeled, got %+v", p.Channel, p)
		}
	}
	if len(b.ActiveModels()) != 0 {
		t.Errorf("new base should have no active models")
	}
}

func TestChannelNamesUnique(t *testing.T) {
	seen := map[string]bool{}
	for _, c := range All() {
		n := c.String()
		if n == "" || seen[n] {
			t.Errorf("channel %d has empty or duplicate name %q", c, n)
		}
		seen[n] = true
		if c.Kind() == AcrossThrough {
			a, th := c.Variables()
			if a == "" || th == "" {
				t.Errorf("across/through channel %v missing variables", c)
			}
		}
	}
}

func TestAttachRequiresCitation(t *testing.T) {
	var b Base
	if err := b.Attach(Light, In, fakeModel{p: 1}, Modeled); !errors.Is(err, ErrNoCitation) {
		t.Fatalf("want ErrNoCitation, got %v", err)
	}
	if err := b.Attach(Light, In, nil, Modeled); !errors.Is(err, ErrNilModel) {
		t.Fatalf("want ErrNilModel, got %v", err)
	}
	if err := b.Attach(ChannelID(200), In, fakeModel{cit: "x"}, Modeled); !errors.Is(err, ErrUnknownChannel) {
		t.Fatalf("want ErrUnknownChannel, got %v", err)
	}
}

func TestAttachAndDetach(t *testing.T) {
	var b Base
	if err := b.Attach(Light, In, fakeModel{p: 2, cit: "test model"}, Modeled); err != nil {
		t.Fatal(err)
	}
	if got := b.ActiveModels(); len(got) != 1 || got[0] != Light {
		t.Fatalf("active = %v", got)
	}
	p, _ := b.Port(Light)
	if p.Status != Modeled || p.Model == nil {
		t.Fatalf("light port = %+v", p)
	}
	if err := b.Detach(Light); err != nil {
		t.Fatal(err)
	}
	if len(b.ActiveModels()) != 0 {
		t.Fatal("detach left an active model")
	}
}

func TestPortsIsACopy(t *testing.T) {
	var b Base
	ps := b.Ports()
	ps[0].Status = Modeled
	p, _ := b.Port(0)
	if p.Status != NotModeled {
		t.Fatal("Ports() must return a copy")
	}
}
