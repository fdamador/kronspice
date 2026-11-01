// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"testing"

	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

func element(name string, ch exchange.ChannelID, status exchange.PortStatus, pw float64) output.ElementResult {
	var b exchange.Base
	ports := b.Ports()
	ports[ch].Status = status
	e := output.ElementResult{Name: name, Ports: ports}
	e.ExteriorPower[ch] = pw
	return e
}

func TestLedgerBalances(t *testing.T) {
	r := &output.Result{Elements: []output.ElementResult{
		element("V1", exchange.Electrical, exchange.Modeled, 2),
		element("R1", exchange.Thermal, exchange.Modeled, -1.5),
		element("V2", exchange.Electrical, exchange.Modeled, -0.5),
	}}
	l := BuildLedger(r)
	if l.Inflow != 2 || l.Dissipation != 1.5 || l.Outflow != 0.5 || l.Residual != 0 {
		t.Fatalf("ledger = %+v", l)
	}
	if err := CheckLedger(l, DefaultTolerance); err != nil {
		t.Fatal(err)
	}
	if len(l.Rows) != 3*exchange.NumChannels {
		t.Fatalf("rows = %d", len(l.Rows))
	}
}

func TestLedgerDetectsImbalance(t *testing.T) {
	// Energy appearing from nowhere must be an error, never a result.
	r := &output.Result{Elements: []output.ElementResult{
		element("V1", exchange.Electrical, exchange.Modeled, 1),
		element("R1", exchange.Thermal, exchange.Modeled, -1.1),
	}}
	if err := CheckLedger(BuildLedger(r), DefaultTolerance); !errors.Is(err, ErrLedgerResidual) {
		t.Fatalf("err = %v, want ErrLedgerResidual", err)
	}
}

func TestHypotheticalNotCounted(t *testing.T) {
	r := &output.Result{Elements: []output.ElementResult{
		element("X1", exchange.HypotheticalMagneticCharge, exchange.Hypothetical, 5),
	}}
	l := BuildLedger(r)
	if l.Inflow != 0 || l.Hypothetical != 5 || l.Residual != 0 {
		t.Fatalf("ledger = %+v", l)
	}
}
