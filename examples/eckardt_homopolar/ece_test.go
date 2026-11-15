// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// With no terms, the circuit with the coil's L as a series element is the
// original machine (the loop current is the same everywhere).
func TestECEPlainTwin(t *testing.T) {
	p := Paper()
	a, _, err := p.RunCircuit(p.L, p.Alpha, 0, 3)
	if err != nil {
		t.Fatal(err)
	}
	b, err := p.RunECE(p.ECECases()[0], p.Alpha, 0, 3)
	if err != nil {
		t.Fatal(err)
	}
	n := len(a.T) - 1
	if m := len(b.Tr.T) - 1; m != n || math.Abs(a.W[n]-b.Tr.W[m]) > 1e-3*math.Abs(a.W[n]) || math.Abs(a.I[n]-b.Tr.I[m]) > 1e-3*math.Abs(a.I[n]) {
		t.Errorf("plain twin: ω %.6g vs %.6g, I %.6g vs %.6g", a.W[n], b.Tr.W[len(b.Tr.W)-1], a.I[n], b.Tr.I[len(b.Tr.I)-1])
	}
	if b.EModel != 0 {
		t.Errorf("plain run books %.3g J model-supplied", b.EModel)
	}
}

// Every ECE / ECE2 run closes its books: drive + flywheel + model-supplied
// pays for R, the coil and the reference resistors, and the ECE terms are
// the only model-supplied power.
func TestECEBooks(t *testing.T) {
	p := Paper()
	for _, k := range p.ECECases()[1:] {
		for _, run := range [][3]float64{{p.Alpha, 0, 3}, {0, 40, 0.5}} {
			r, err := p.RunECE(k, run[0], run[1], run[2])
			if err != nil {
				t.Fatalf("%s: %v", k.Name, err)
			}
			if r.EModel == 0 {
				t.Errorf("%s: no model-supplied energy booked", k.Name)
			}
			scale := r.Tr.EDrive + math.Abs(r.Tr.ERot) + math.Abs(r.EModel) + r.EHeatR
			if math.Abs(r.Mismatch) > 1e-5*scale {
				t.Errorf("%s α=%g: ledger mismatch %.3g J of %.3g J", k.Name, run[0], r.Mismatch, scale)
			}
		}
	}
}

// ECE2 torsion is a negative series resistance: the driven machine settles
// near the lowered ω_res, and R heats beyond what drive and flywheel supply.
func TestECE2TorsionLowersResonance(t *testing.T) {
	p := Paper()
	k := p.ECECases()[2]
	r, err := p.RunECE(k, p.Alpha, 0, 3)
	if err != nil {
		t.Fatal(err)
	}
	if w := p.ECEOmegaRes(k); math.Abs(w-14.92) > 0.01 {
		t.Errorf("ECE2 ω_res = %.3f", w)
	}
	if r.EHeatR <= r.Tr.EDrive-r.Tr.ERot || r.EModel <= 0 {
		t.Errorf("heat %.4g J, drive %.4g J, spin-up %.4g J, model %.4g J", r.EHeatR, r.Tr.EDrive, r.Tr.ERot, r.EModel)
	}
}
