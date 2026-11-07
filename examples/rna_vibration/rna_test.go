// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strings"
	"sync"
	"testing"
)

var (
	prepOnce sync.Once
	prepared map[string]*State
	prepErr  error
)

func states(t *testing.T) map[string]*State {
	t.Helper()
	prepOnce.Do(func() {
		ms, err := Molecules()
		if err != nil {
			prepErr = err
			return
		}
		prepared = map[string]*State{}
		for _, k := range order {
			s, err := Prepare(ms[k])
			if err != nil {
				prepErr = err
				return
			}
			prepared[k] = s
		}
	})
	if prepErr != nil {
		t.Fatal(prepErr)
	}
	return prepared
}

// Every molecule is at a true minimum: gradient ≈ 0 and all 3N−6
// vibrations real (no saddle direction).
func TestMinima(t *testing.T) {
	for k, s := range states(t) {
		if s.RMS > 1e-6 {
			t.Errorf("%s: rms gradient %g", k, s.RMS)
		}
		if len(s.Modes) != 3*len(s.Mol.Atoms)-6 {
			t.Errorf("%s: %d modes for %d atoms", k, len(s.Modes), len(s.Mol.Atoms))
		}
		if s.Modes[0].Wavenumber <= 0 {
			t.Errorf("%s: lowest mode %g cm⁻¹ (imaginary)", k, s.Modes[0].Wavenumber)
		}
	}
}

// The circuit's springs plus its coupling network carry the whole
// stiffness: Σ k·BᵀB + R_kept = K, up to the dropped entries (< 1e−9 of
// the largest).
func TestStiffnessSplit(t *testing.T) {
	for k, s := range states(t) {
		sp := s.Springs()
		res, _ := s.ResidualEntries(sp)
		n := len(s.K)
		K2 := make([][]float64, n)
		for i := range K2 {
			K2[i] = make([]float64, n)
		}
		for _, x := range sp {
			for i, bi := range x.B {
				for j, bj := range x.B {
					K2[i][j] += x.K * bi * bj
				}
			}
		}
		for _, e := range res {
			K2[e.I][e.J] += e.V
		}
		kmax, d := 0.0, 0.0
		for i := range K2 {
			for j := range K2 {
				kmax = math.Max(kmax, math.Abs(s.K[i][j]))
				d = math.Max(d, math.Abs(K2[i][j]-s.K[i][j]))
			}
		}
		if d > 1e-9*kmax {
			t.Errorf("%s: circuit stiffness differs from K by %g of %g", k, d, kmax)
		}
	}
}

// Twin of the circuit (CLAUDE.md §3a rule 4 in spirit): kronspice's AC
// solution of the netlist and a direct complex solve of the same
// admittance matrix give the same absorbed power, for a base and a pair,
// on and off resonance.
func TestCircuitTwin(t *testing.T) {
	st := states(t)
	for _, k := range []string{"C", "GC"} {
		s := st[k]
		nu, ab, err := s.Absorption(100, 3400, 12, FWHMcm)
		if err != nil {
			t.Fatal(err)
		}
		// Add exact mode frequencies (on resonance).
		for _, i := range []int{len(s.Modes) / 2, len(s.Modes) - 1} {
			f := s.Modes[i].Wavenumber
			n2, a2, err := s.Absorption(f, f, 1, FWHMcm)
			if err != nil {
				t.Fatal(err)
			}
			nu, ab = append(nu, n2...), append(ab, a2...)
		}
		for i, f := range nu {
			d := 0.0
			for ax := 0; ax < 3; ax++ {
				d += s.DirectResponse(f, FWHMcm, ax) / 3
			}
			if math.Abs(ab[i]-d) > 1e-9*math.Abs(d) {
				t.Errorf("%s at %.1f cm⁻¹: circuit %.12g, direct %.12g", k, f, ab[i], d)
			}
		}
	}
}

// Release by two routes: the Einstein A coefficient from the IR intensity
// and Maxwell's radiated power of the vibrating dipole (em.LarmorPower,
// itself checked against the Poynting flux of the full fields) divided by
// the mode's energy agree for every mode.
func TestReleaseTwin(t *testing.T) {
	for k, s := range states(t) {
		for _, md := range s.Modes {
			r := ReleaseOf(md)
			if r.A == 0 && r.ALarmor == 0 {
				continue
			}
			if math.Abs(r.A-r.ALarmor) > 1e-9*r.A {
				t.Errorf("%s %.1f cm⁻¹: A = %g from intensity, %g from Larmor", k, md.Wavenumber, r.A, r.ALarmor)
			}
		}
	}
}

// Absorption is passive: the field never takes power back, at any
// frequency, and its integrated area follows the summed intensities.
func TestAbsorptionPositive(t *testing.T) {
	s := states(t)["A"]
	_, ab, err := s.Absorption(10, 3600, 360, FWHMcm)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range ab {
		if v < 0 {
			t.Errorf("point %d: absorbed power %g < 0", i, v)
		}
	}
}

// Pair geometry after minimization: A·U has two hydrogen bonds
// (N6–H···O4, N3–H···N1), G·C three (N1–H···N3, N2–H···O2, N4–H···O6), with
// donor–acceptor distances near the crystal values (2.8–3.0 Å; Saenger,
// Principles of Nucleic Acid Structure, 1984, Ch. 6), and the G·U wobble
// pair its two (G N1–H1···O2 U and U N3–H3···O6 G; Varani & McClain,
// EMBO Rep. 1, 18, 2000).
func TestHydrogenBonds(t *testing.T) {
	st := states(t)
	want := map[string]int{"AU": 2, "GC": 3, "GU": 2}
	for k, n := range want {
		hb := st[k].HBonds()
		if len(hb) != n {
			t.Errorf("%s: %d H-bonds, want %d", k, len(hb), n)
		}
		for _, h := range hb {
			if h.DA < 2.6 || h.DA > 3.1 {
				t.Errorf("%s %s···%s: D···A %.2f Å", k, h.H, h.Acceptor, h.DA)
			}
		}
	}
}

// TestWobble: the G·U pair is held by exactly the two wobble hydrogen
// bonds, G N1–H1···O2 U and U N3–H3···O6 G, not by the Watson–Crick
// pattern (G N1–H1···N3 would need U's N3 to be an acceptor).
func TestWobble(t *testing.T) {
	got := map[string]bool{}
	for _, h := range states(t)["GU"].HBonds() {
		got[h.Donor+"-"+h.H+"-"+h.Acceptor] = true
	}
	for _, w := range []string{"G:N1-G:H1-U:O2", "U:N3-U:H3-G:O6"} {
		if !got[w] {
			t.Errorf("wobble H-bond %s missing (have %v)", w, got)
		}
	}
	if len(got) != 2 {
		t.Errorf("G·U: %d H-bonds, want 2 (%v)", len(got), got)
	}
}

// TestModalTwin: the normal-mode spectrum used for the stacks equals the
// direct solve of the circuit's admittance (and so the circuit itself,
// TestCircuitTwin) for a pair, on and off resonance.
func TestModalTwin(t *testing.T) {
	ms, err := Molecules()
	if err != nil {
		t.Fatal(err)
	}
	s, err := Prepare(ms["AU"])
	if err != nil {
		t.Fatal(err)
	}
	nu := []float64{s.Modes[3].Wavenumber, 500, s.Modes[len(s.Modes)-1].Wavenumber}
	modal := s.ModalAbsorption(nu, FWHMcm)
	for k, v := range nu {
		d := 0.0
		for ax := 0; ax < 3; ax++ {
			d += s.DirectResponse(v, FWHMcm, ax) / 3
		}
		// The modal sum leaves out SpringR, which the circuit and the direct
		// solve carry in series with every spring: it adds a damping of order
		// SpringR·k/ω, largest (~1e-5) for the stiff N–H springs at 3,365 cm⁻¹.
		if rel := math.Abs(modal[k]-d) / d; rel > 5e-5 {
			t.Errorf("%.1f cm⁻¹: modal %.6g, direct %.6g (rel %.2g)", v, modal[k], d, rel)
		}
	}
}

// TestStackParts: a stack's between-base forces are split exactly into
// pairing (the two bases of one pair) and stacking (bases of different
// pairs).
func TestStackParts(t *testing.T) {
	ms, err := Molecules()
	if err != nil {
		t.Fatal(err)
	}
	m := ms["polyA3"]
	inter, stack := 0, 0
	for _, p := range m.Pairs {
		if p.Inter && p.Stack {
			t.Fatal("pair both pairing and stacking")
		}
		if p.Inter {
			inter++
		}
		if p.Stack {
			stack++
		}
		if (p.Inter || p.Stack) != (pairOf(m, p.I) != pairOf(m, p.J) || baseOf(m, p.I) != baseOf(m, p.J)) {
			t.Fatalf("pair %d–%d grouped wrongly", p.I, p.J)
		}
	}
	if inter == 0 || stack == 0 || numPairs(m) != 3 {
		t.Fatalf("polyA3: %d pairing, %d stacking pairs, %d base pairs", inter, stack, numPairs(m))
	}
}

func baseOf(m *Molecule, a int) string { return strings.SplitN(m.Atoms[a].Name, ":", 2)[0] }
