// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
)

// The stack's valve: it conducts anode-positive past the electrolysis
// onset, holds anode-negative with the Güntherschulze–Betz leakage reaching
// jForm at the forming voltage, and its film capacitance is ε·A/(k·Vform).
func TestStack(t *testing.T) {
	s := Stack{Ns: 2, Np: 3, Area: 1e-3, Gap: 5e-3, Vform: 400}
	if i := s.Current(1.5); math.Abs(i) > 1e-6 {
		t.Errorf("below the 4 V onset of two cells it conducts %.3g A", i)
	}
	if i := s.Current(10); i <= 0 {
		t.Errorf("forward 10 V gives %.3g A", i)
	}
	jf := -s.Current(-s.Hold()) / (s.Area * float64(s.Np))
	if math.Abs(jf-jForm)/jForm > 0.05 {
		t.Errorf("leakage at the holding voltage %.3g A/m², want %.3g", jf, jForm)
	}
	if r := s.Current(-0.5*s.Hold()) / s.Current(-s.Hold()); r > 1e-3 {
		t.Errorf("leakage at half the holding voltage is %.3g of that at full: not exponential", r)
	}
	if s.Current(-1.2*s.Hold()) >= 20*s.Current(-s.Hold()) {
		t.Error("above the forming voltage the oxide should spark")
	}
	c := epsOx * 8.8541878128e-12 * s.Area * 3 / (kOx * 400) / 2
	if math.Abs(s.C0()-c)/c > 1e-9 {
		t.Errorf("film %.4g F, want %.4g", s.C0(), c)
	}
}

// The photon_electron baseline: silicon diodes pump hundreds of volts but
// the 1 µF tank's circulation heats the coil; aluminium stacks at 100 kHz
// are shorted by their films. Both lose nearly everything.
func TestBaseline(t *testing.T) {
	for _, diode := range []bool{true, false} {
		p := DefaultPump()
		p.Diode = diode
		p.Cycles, p.Window = 300, 60
		r, err := p.Run(analysis.Options{})
		if err != nil {
			t.Fatal(err)
		}
		if r.Out > -100 || r.Eff > 0.01 || r.LossCoil/r.Pin < 0.9 {
			t.Errorf("diode %v: %.0f V, η %.4f, coil %.2f", diode, r.Out, r.Eff, r.LossCoil/r.Pin)
		}
	}
}

// The most efficient evolved design (seed 1): one cell of six 13 mm²
// plates formed at 545 V, driven at 471 Hz. It must pump below −300 V at
// η > 0.8, its losses must account for the rest of the input, and the
// ledger must close.
func TestEvolvedDesign(t *testing.T) {
	p := PumpFrom([]float64{471.3, 50, 2.52e-7, 1, 8.34e-7, 4.13e4, 1, 6, 1.323e-5, 7.91e-3, 544.9})
	r, err := p.Run(analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	if r.Out > -300 || r.Eff < 0.8 || r.Drift > 0.01 || r.RevRatio > 0.9 {
		t.Errorf("%.0f V, η %.3f, drift %.3g, reverse %.2f", r.Out, r.Eff, r.Drift, r.RevRatio)
	}
	sum := (r.Pout + r.LossCoil + r.LossForward + r.LossLeak + r.Rest) / r.Pin
	if math.Abs(sum-1) > 0.01 {
		t.Errorf("load plus losses are %.4f of the input", sum)
	}
	if r.Ledger > 1e-3 {
		t.Errorf("ledger mismatch %.3g", r.Ledger)
	}
	if r.ElPerJ > r.Ideal {
		t.Errorf("%.3g electrons per joule beats the e·|V| bound %.3g", r.ElPerJ, r.Ideal)
	}
}

// A small evolution runs and keeps only designs that meet the constraints.
func TestEvolveSmall(t *testing.T) {
	if testing.Short() {
		t.Skip("evolution")
	}
	ev, err := EvolvePump(8, 2, 1, func(string) {})
	if err != nil {
		t.Fatal(err)
	}
	if len(ev) == 0 {
		t.Fatal("no designs")
	}
	for _, e := range ev {
		if e.R.Drift > 0.01 || e.R.RevRatio > 0.9 || e.R.Eff <= 0 || e.R.Eff >= 1 {
			t.Errorf("design breaks a constraint: η %.3f, drift %.3g, reverse %.2f", e.R.Eff, e.R.Drift, e.R.RevRatio)
		}
	}
}

// Concentric pipes: the parallel tank alternates aluminium faces and holds
// one forming voltage; the series tank puts every aluminium face outside
// its pipe, holds N forming voltages and has the smaller film.
func TestPipesGeometry(t *testing.T) {
	par := Pipes{N: 3, R0: 2e-3, Gap: 3e-3, Wall: wallDef, Len: 0.05, Vform: 500}
	ser := par
	ser.Series = true
	a0, b0 := par.radii(0)
	a1, b1 := par.radii(1)
	if par.alFace(0) != a0 || par.alFace(1) != b1 || ser.alFace(1) != a1 {
		t.Errorf("aluminium faces: parallel %g %g, series %g", par.alFace(0), par.alFace(1), ser.alFace(1))
	}
	_, rw, _ := par.Cell(0)
	if want := rhoWater * math.Log(b0/a0) / (2 * math.Pi * 0.05); math.Abs(rw-want)/want > 1e-12 {
		t.Errorf("annulus resistance %g, want %g", rw, want)
	}
	if par.Hold() != 500 || ser.Hold() != 1500 || ser.Film() >= par.Film()/3 {
		t.Errorf("hold %g / %g, film %g / %g", par.Hold(), ser.Hold(), par.Film(), ser.Film())
	}
}

// A single-cell pipe tank pumps like the best plate stack. A plain
// three-annulus series tank loads its inner annulus most; graded, the
// annuli share equally. Perfect and 1 GΩ seals agree; with pipes merely
// standing in the water (~10 Ω between annuli) the leak shorts the outer
// annuli and the inner one is driven past its forming voltage. Losses
// (seal included) account for the input.
func TestPipesPump(t *testing.T) {
	base := []float64{471.3, 50, 2.52e-7, 1, 8.34e-7, 4.13e4}
	one := PipesFrom(append(append([]float64{}, base...), 1, 1, 1e-3, 7.9e-3, 12.6e-3, 545))
	r, err := one.Run(analysis.Options{})
	if err != nil || r.Eff < 0.8 || r.Out > -300 {
		t.Fatalf("single cell: η %.3f, %.0f V, %v", r.Eff, r.Out, err)
	}
	plain := PipesFrom(append(append([]float64{}, base...), 2, 3, 1e-3, 3e-3, 5e-3, 200))
	rp, err := plain.Run(analysis.Options{})
	if err != nil || rp.CellRev[0] < 2*rp.CellRev[2] {
		t.Errorf("plain series: annuli %.2f, %v", rp.CellRev, err)
	}
	graded := PipesFrom(append(append([]float64{}, base...), 3, 3, 1e-3, 3e-3, 5e-3, 200))
	var res []Result
	for _, rs := range []float64{0, 1e9, 10} {
		p := graded
		q := *p.Pipes
		q.Rseal = rs
		p.Pipes = &q
		r, err := p.Run(analysis.Options{})
		if err != nil && err != ErrNoPump {
			t.Fatal(err)
		}
		sum := (r.Pout + r.LossCoil + r.LossForward + r.LossLeak + r.LossSeal + r.Rest) / r.Pin
		if math.Abs(sum-1) > 1e-9 || math.Abs(r.Rest/r.Pin) > 0.01 {
			t.Errorf("Rseal %g: losses sum %.6f, rest %.3g", rs, sum, r.Rest/r.Pin)
		}
		res = append(res, r)
	}
	c := res[0].CellRev
	if res[0].Eff < 0.8 || math.Abs(c[0]-c[2]) > 0.02*c[0] || res[0].RevRatio > 0.9 {
		t.Errorf("graded, perfect seals: η %.3f, annuli %.2f", res[0].Eff, c)
	}
	if math.Abs(res[1].Eff-res[0].Eff) > 0.01 || res[2].RevRatio < 0.99 || res[2].Eff > res[0].Eff-0.05 {
		t.Errorf("graded: 1 GΩ η %.3f; 10 Ω η %.3f, worst oxide %.2f", res[1].Eff, res[2].Eff, res[2].RevRatio)
	}
}
