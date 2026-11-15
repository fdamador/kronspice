// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"testing"
)

// With a core far better than the air (μr → ∞) every pass carries the
// same flux and the ratio is (N_s·n_s)/(N_p·n_p) exactly.
func TestIdealLimit(t *testing.T) {
	m := Materials[2]
	m.Mur = 1e12
	s := DefaultSpring(m)
	pri, sec := designWindings(s.N)
	for _, ns := range []int{1, 2, 5, 10} {
		q := sec
		q.First, q.Passes = 0, ns
		r, err := s.Ratio(pri, q)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r/IdealRatio(pri, q)-1) > 1e-9 {
			t.Errorf("n_s %d: ratio %.12f, ideal %.12f", ns, r, IdealRatio(pri, q))
		}
	}
}

// The magnetic network is reciprocal (L symmetric) and conserves flux:
// at every node the segment fluxes in and out, the leakage and the return
// leg balance.
func TestNetwork(t *testing.T) {
	for _, m := range Materials {
		s := DefaultSpring(m)
		pri, sec := designWindings(s.N)
		axis := Winding{Name: "axis", Azimuth: -math.Pi / 2, Passes: s.N, Turns: 1}
		L, fs, err := s.Inductance([]Winding{pri, sec, axis})
		if err != nil {
			t.Fatal(err)
		}
		for i := range L {
			for j := range L {
				if d := math.Abs(L[i][j]-L[j][i]) / math.Sqrt(L[i][i]*L[j][j]); d > 1e-10 {
					t.Errorf("%s: L[%d][%d] asymmetry %.2g", m.Name, i, j, d)
				}
			}
		}
		f := fs[1]
		M := s.N * s.S
		for n := 0; n <= M; n++ {
			net := 0.0 // flux leaving node n
			if n < M {
				net += f.Seg[n]
			}
			if n > 0 {
				net -= f.Seg[n-1]
			}
			if n+s.S <= M {
				net += f.Leak[n]
			}
			if n-s.S >= 0 {
				net -= f.Leak[n-s.S]
			}
			if n == M {
				net += f.Return
			}
			if n == 0 {
				net -= f.Return
			}
			if math.Abs(net) > 1e-9*math.Abs(f.Return) {
				t.Fatalf("%s: node %d flux imbalance %.3g of %.3g", m.Name, n, net, f.Return)
			}
		}
	}
}

// The design keeps the flux on the wire: the same dΦ/dt in every pass,
// the ratio 2/N for N core turns, the bore wire N/(N_p·n_p).
func TestDesign(t *testing.T) {
	s := DefaultSpring(Materials[2])
	lk, err := s.Leakage()
	if err != nil {
		t.Fatal(err)
	}
	if lk.Spread > 1e-3 || lk.K < 0.98 {
		t.Errorf("design leakage: spread %.2g, k %.4f", lk.Spread, lk.K)
	}
	for _, n := range []int{2, 10, 40} {
		s.N = n
		pri, sec := designWindings(n)
		r, err := s.Ratio(pri, sec)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(r*float64(n)/2-1) > 1e-4 {
			t.Errorf("N = %d: ratio %.6f, 2/N %.6f", n, r, 2/float64(n))
		}
		axis := Winding{Name: "axis", Azimuth: -math.Pi / 2, Passes: n, Turns: 1}
		ra, err := s.Ratio(pri, axis)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(ra/IdealRatio(pri, axis)-1) > 1e-4 {
			t.Errorf("N = %d: axis ratio %.6g, ideal %.6g", n, ra, IdealRatio(pri, axis))
		}
	}
}

// Every lever that section 1 relies on lowers the leakage.
func TestLevers(t *testing.T) {
	for _, m := range Materials {
		plain := DefaultSpring(m)
		plain.A, plain.P, plain.RetArea = 1e-3, 3e-3, 1
		lp, _ := plain.Leakage()
		ld, _ := DefaultSpring(m).Leakage()
		if !(ld.K > lp.K && ld.Spread < lp.Spread) {
			t.Errorf("%s: design k %.4f spread %.2g vs plain k %.4f spread %.2g", m.Name, ld.K, ld.Spread, lp.K, lp.Spread)
		}
	}
}

// kronspice agrees with the network (up to the declared 1 % air leakage
// and the copper), closes its ledger and never returns more than it takes.
func TestCircuit(t *testing.T) {
	s := DefaultSpring(Materials[2])
	pri, sec := designWindings(s.N)
	cr, err := s.RunCircuit(DefaultRig(), pri, sec)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(cr.RatioAC/0.2-1) > 0.02 {
		t.Errorf("AC ratio %.5f", cr.RatioAC)
	}
	if math.Abs(cr.Mismatch) > 1e-5*cr.EIn {
		t.Errorf("ledger %.3g of %.3g J", cr.Mismatch, cr.EIn)
	}
	if !(cr.Pout < cr.Pin && cr.Pout > 0.9*cr.Pin) {
		t.Errorf("Pin %.4g W, Pout %.4g W", cr.Pin, cr.Pout)
	}
	if !s.EddyOK(DefaultRig().Freq) || cr.Bpk > 0.5*s.Mat.Bsat {
		t.Errorf("design outside its model: eddy %v, B %.3f T", s.EddyOK(DefaultRig().Freq), cr.Bpk)
	}
}
