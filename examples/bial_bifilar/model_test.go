// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"errors"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
)

func baseX() []float64 {
	// 200 turns, 30 mm coil, 0.4 mm wire at 1.5× pitch, ε 3; 40 mm plates,
	// 0.2 mm gap, dielectric ε 10.
	return []float64{200, 30e-3, 0.4e-3, 1.5, 3, 40e-3, 0.2e-3, 10}
}

// Tesla's comparison at 1,000 turns: the bifilar winding stores
// (N/2)² = 250,000 times the turn-to-turn energy of the ordinary one.
func TestTeslaComparison(t *testing.T) {
	x := baseX()
	x[0], x[1] = 1000, 0.1
	b := designFrom(Dry, DefaultAssumptions(), true, x)
	o := designFrom(Dry, DefaultAssumptions(), false, x)
	if r := b.Coil.SelfCapacitance() / o.Coil.SelfCapacitance(); math.Abs(r-250000) > 1e-6 {
		t.Fatalf("ratio %g", r)
	}
}

// A parallel RLC tank: Resonance finds f0 = 1/(2π√(LC)) and Q = R√(C/L).
func TestResonanceOfKnownTank(t *testing.T) {
	c, err := circuit.NewBuilder("rlc").Resistor("R1", "a", "0", 1e3).Inductor("L1", "a", "0", 1e-3).
		Capacitor("C1", "a", "0", 1e-9).
		Add(circuit.Element{Name: "Ip", Kind: circuit.CurrentSource, HasDC: true, ACMag: 1}, "0", "a").Build()
	if err != nil {
		t.Fatal(err)
	}
	f0, q, zmax, err := Resonance(c, 1, 1e9, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := 1 / (2 * math.Pi * math.Sqrt(1e-12))
	if math.Abs(f0-want) > 1e-9*want || math.Abs(q-1e3*math.Sqrt(1e-6)) > 1e-6 || math.Abs(zmax-1e3) > 1e-6 {
		t.Fatalf("f0 %g (want %g), Q %g (want 1), |Z| %g", f0, want, q, zmax)
	}
}

// Dry: the EMF only charges the plates: stored energy ½·C·EMF², no steady
// current, and the energy comes in through the thermal port.
func TestDrySteadyState(t *testing.T) {
	d := designFrom(Dry, DefaultAssumptions(), true, baseX())
	r, err := Measure(d, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := 0.5 * r.Element.Cplate * r.Element.EMF * r.Element.EMF
	if math.Abs(r.StoredJ-want) > 1e-9*want {
		t.Fatalf("stored %g, want %g", r.StoredJ, want)
	}
	if math.Abs(r.SourcePower) > 1e-30 {
		t.Fatalf("steady source power %g", r.SourcePower)
	}
	c, _ := d.Build(r.Rcoil, true)
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: 1 / (40 * r.F0), Stop: 200 / r.F0, MaxStep: 1 / (40 * r.F0)},
		analysis.TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if wf.Ledger.EnergyIn <= 0 || wf.Ledger.EnergyModelSupplied != 0 || wf.Ledger.MaxRelativeResidual > 1e-9 {
		t.Fatalf("ledger %+v", wf.Ledger)
	}
	op, _ := analysis.OP(c, analysis.Options{})
	if e, _ := op.Element("Vbial"); e.Ports[exchange.Thermal].Status != exchange.Modeled ||
		e.Ports[exchange.Electrical].Note != "terminals connect inside the network" {
		t.Fatalf("Seebeck source ports: %+v", e.Ports[:2])
	}
}

// Wet: a steady current EMF/(R_el + R_ct + R_coil) flows through the coil;
// the stored energy is ½·L·I² plus the double layer's ½·C_dl·V².
func TestWetSteadyState(t *testing.T) {
	d := designFrom(Wet, DefaultAssumptions(), true, baseX())
	r, err := Measure(d, analysis.Options{})
	if err != nil && !errors.Is(err, ErrNoResonance) {
		t.Fatal(err)
	}
	el := r.Element
	rdc := d.Coil.Resistance(0)
	i := el.EMF / (el.Rint + el.Rct + rdc)
	want := 0.5*r.L*i*i + 0.5*el.Cdl*(i*el.Rct)*(i*el.Rct) + 0.5*r.Cself*(i*rdc)*(i*rdc) + 0.5*el.Cgeo*(i*rdc)*(i*rdc)
	if math.Abs(r.StoredJ-want) > 1e-9*want {
		t.Fatalf("stored %g, want %g", r.StoredJ, want)
	}
	if math.Abs(r.SourcePower-el.EMF*i) > 1e-9*el.EMF*i || math.Abs(r.Heat-r.SourcePower) > 1e-9*r.SourcePower {
		t.Fatalf("source %g W, heat %g W, want %g", r.SourcePower, r.Heat, el.EMF*i)
	}
}

func TestOutsideValidityIsRejected(t *testing.T) {
	x := baseX()
	x[0] = 3 // 3 turns of a 30 mm coil: far too short for Wheeler
	_, err := Measure(designFrom(Dry, DefaultAssumptions(), true, x), analysis.Options{})
	if !errors.Is(err, em.ErrOutsideModelValidity) {
		t.Fatalf("err = %v", err)
	}
}

// The whole experiment, small: both variants evolve, fronts are written,
// and the report states that the ledger closes.
func TestQuickExperiment(t *testing.T) {
	dir := t.TempDir()
	if err := run(context.Background(), DefaultAssumptions(), 10, 3, 1, dir, t.Logf); err != nil {
		t.Fatal(err)
	}
	rep, err := os.ReadFile(filepath.Join(dir, "report.md"))
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(rep), "accounted for by the element's declared source") {
		t.Fatalf("ledger did not close:\n%s", rep)
	}
	for _, f := range []string{"front_dry.csv", "front_wet.csv"} {
		if _, err := os.Stat(filepath.Join(dir, f)); err != nil {
			t.Error(err)
		}
	}
}

// Thermopile: n pairs multiply the EMF by n, the stored energy by n², and
// the legs add series resistance, volume and a conductive heat leak.
func TestThermopileScaling(t *testing.T) {
	a := DefaultAssumptions()
	const n, dia, ln = 200, 0.5e-3, 10e-3
	dry, err := designFrom(Dry, a, true, baseX()).ElementModel()
	if err != nil {
		t.Fatal(err)
	}
	d := designFrom(Thermopile, a, true, append(baseX(), n, dia, ln))
	e, err := d.ElementModel()
	if err != nil {
		t.Fatal(err)
	}
	if want := n * 75.5e-6 * a.DeltaT; math.Abs(e.EMF-want) > 1e-12*want {
		t.Fatalf("EMF %g, want %g", e.EMF, want)
	}
	la := math.Pi * dia * dia / 4
	rho := material.Bismuth.Resistivity + material.Aluminium.Resistivity
	if want := dry.Rint + n*rho*ln/la; math.Abs(e.Rint-want) > 1e-12*want {
		t.Fatalf("Rint %g, want %g", e.Rint, want)
	}
	if want := dry.Volume + 2*n*la*ln; math.Abs(e.Volume-want) > 1e-12*want {
		t.Fatalf("volume %g, want %g", e.Volume, want)
	}
	if want := n * (7.97 + 237) * la / ln * a.DeltaT; math.Abs(e.HeatLeak-want) > 1e-12*want {
		t.Fatalf("heat leak %g, want %g", e.HeatLeak, want)
	}
	r, err := Measure(d, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	if want := 0.5 * e.Cplate * e.EMF * e.EMF; math.Abs(r.StoredJ-want) > 1e-9*want {
		t.Fatalf("stored %g, want %g", r.StoredJ, want)
	}
	// Legs shorter than five diameters are outside the 1-D conduction model.
	if _, err := designFrom(Thermopile, a, true, append(baseX(), n, 1e-3, 3e-3)).ElementModel(); !errors.Is(err, em.ErrOutsideModelValidity) {
		t.Fatalf("short legs: %v", err)
	}
}
