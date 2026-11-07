// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// Absorption runs kronspice's AC analysis of the circuit driven by a unit
// infrared field along x, y and z in turn (one linear sweep each, from nu0
// to nu1 cm⁻¹ in points steps), and returns, at each wavenumber, the power
// the field puts into the molecule, ½·Re Σ q_a·v_a*, averaged over the
// three directions (an isotropic sample). Arbitrary units: the line areas
// follow the IR intensities, the widths the display FWHM.
func (s *State) Absorption(nu0, nu1 float64, points int, fwhm float64) ([]float64, []float64, error) {
	conv := 2 * math.Pi * wavenumber
	var nu, out []float64
	for d, ax := range []string{"x", "y", "z"} {
		src, _ := s.Netlist(NetlistOptions{FWHM: fwhm, Field: ax})
		deck, err := netlist.ParseString(src)
		if err != nil {
			return nil, nil, err
		}
		c, err := circuit.FromDeck(deck)
		if err != nil {
			return nil, nil, err
		}
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: points, FStart: nu0 / conv, FStop: nu1 / conv}, analysis.Options{})
		if err != nil {
			return nil, nil, fmt.Errorf("%s: %w", s.Mol.Title, err)
		}
		if out == nil {
			out = make([]float64, len(r.Freq))
			nu = make([]float64, len(r.Freq))
			for k, f := range r.Freq {
				nu[k] = f * conv
			}
		}
		for a := range s.Mol.Atoms {
			v, ok := r.NodeVoltage(s.node(3*a + d))
			if !ok {
				return nil, nil, fmt.Errorf("node %s missing", s.node(3*a+d))
			}
			q := complex(s.Mol.Atoms[a].Charge, 0)
			for k := range v {
				out[k] += 0.5 * real(q*cmplx.Conj(v[k])) / 3
			}
		}
	}
	return nu, out, nil
}

// DirectResponse solves (K − ω²M + iωD)·x = F directly (complex Gaussian
// elimination, no circuit) for the same drive as the circuit, and returns
// the absorbed power ½·Re Σ F·v*. It is the independent twin of
// Absorption: the same numbers from the matrices, not from kronspice.
func (s *State) DirectResponse(nu float64, fwhm float64, axis int) float64 {
	n := len(s.K)
	conv := 2 * math.Pi * wavenumber
	w := 2 * math.Pi * nu / conv // circuit rad/s
	gamma := 2 * math.Pi * fwhm / conv
	// Admittance matrix of the circuit, assembled from the same parts as the
	// netlist but with no circuit simulator: masses and damping, each spring
	// y = 1/(iωL + SpringR) between its B-row combination of nodes, and each
	// residual coupling through the integrator x = v/(iω·1 + 1/LeakR).
	a := make([][]complex128, n)
	f := make([]complex128, n)
	for i := 0; i < n; i++ {
		a[i] = make([]complex128, n)
		mass := s.Mol.Atoms[i/3].Mass
		a[i][i] = complex(mass*gamma, w*mass)
		if i%3 == axis {
			f[i] = complex(s.Mol.Atoms[i/3].Charge, 0)
		}
	}
	springs := s.Springs()
	for _, sp := range springs {
		y := 1 / complex(SpringR, w/sp.K)
		for i, bi := range sp.B {
			for j, bj := range sp.B {
				a[i][j] += y * complex(bi*bj, 0)
			}
		}
	}
	res, _ := s.ResidualEntries(springs)
	xi := 1 / complex(1/LeakR, w)
	for _, e := range res {
		a[e.I][e.J] += complex(e.V, 0) * xi
	}
	v := solveComplex(a, append([]complex128(nil), f...)) // the solve overwrites its right-hand side
	p := 0.0
	for i := range v {
		p += 0.5 * real(f[i]*cmplx.Conj(v[i]))
	}
	return p
}

func solveComplex(a [][]complex128, b []complex128) []complex128 {
	n := len(b)
	for k := 0; k < n; k++ {
		p := k
		for i := k + 1; i < n; i++ {
			if cmplx.Abs(a[i][k]) > cmplx.Abs(a[p][k]) {
				p = i
			}
		}
		a[k], a[p] = a[p], a[k]
		b[k], b[p] = b[p], b[k]
		for i := k + 1; i < n; i++ {
			f := a[i][k] / a[k][k]
			for j := k; j < n; j++ {
				a[i][j] -= f * a[k][j]
			}
			b[i] -= f * b[k]
		}
	}
	x := make([]complex128, n)
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for j := i + 1; j < n; j++ {
			s -= a[i][j] * x[j]
		}
		x[i] = s / a[i][i]
	}
	return x
}
