// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// The s states of an electron bound to a charge, in atomic units of the
// reduced mass: u(r) = rψ(r), −½u″ − (Z/r)u = E u, u(0) = u(r_max) = 0.
// On a grid r_i = i·h (i = 1…N) this is Kron's 1945 network (Phys. Rev. 67,
// 39): a ladder whose series inductors between neighbours carry the
// kinetic energy (admittance 1/(2h²)), whose shunt inductors to ground
// carry the potential, and whose unit capacitors carry the energy; the
// boundary nodes are short circuits (Kron's "short circuit at a point of
// zero voltage"). With the potential shifted by V0 so every shunt is a
// real inductor, ω² = E + V0, and the ground state is the network's lowest
// resonance.

// Ladder is the tridiagonal matrix of the network: diagonal D, off-diagonal
// O (the same for every pair).
type Ladder struct {
	H float64
	D []float64
	O float64
}

// NewLadder builds the s-state network for coupling Z, spacing h, edge rmax.
func NewLadder(z, h, rmax float64) Ladder { return NewLadderL(z, 0, h, rmax) }

// NewLadderL builds the network for orbital angular momentum l: the
// centrifugal term l(l+1)/(2r²) joins the potential in the shunt branches
// (u = rψ, −½u″ + [l(l+1)/2r² − Z/r]u = E u). l = 1 gives the p states,
// the two-lobed 2p_z among them.
func NewLadderL(z float64, ang int, h, rmax float64) Ladder {
	n := int(math.Round(rmax/h)) - 1
	l := Ladder{H: h, D: make([]float64, n), O: -1 / (2 * h * h)}
	c := float64(ang*(ang+1)) / 2
	for i := range l.D {
		r := float64(i+1) * h
		l.D[i] = 1/(h*h) - z/r + c/(r*r)
	}
	return l
}

// below counts eigenvalues below x (Sturm sequence of the LDLᵀ pivots).
func (l Ladder) below(x float64) int {
	n, d := 0, 1.0
	for i, di := range l.D {
		if i == 0 {
			d = di - x
		} else {
			if d == 0 {
				d = 1e-300
			}
			d = di - x - l.O*l.O/d
		}
		if d < 0 {
			n++
		}
	}
	return n
}

// Eigen returns the k-th eigenvalue (k = 0: ground) by bisection on the
// Sturm count.
func (l Ladder) Eigen(k int) float64 {
	lo, hi := -1e3, 1e3
	for i := 0; i < 200 && hi-lo > 1e-14*math.Max(1, math.Abs(lo)); i++ {
		m := 0.5 * (lo + hi)
		if l.below(m) > k {
			hi = m
		} else {
			lo = m
		}
	}
	return 0.5 * (lo + hi)
}

// Extrapolated is the 1s and 2s energies [Ha] for coupling Z, Richardson-
// extrapolated from h = 0.02 and 0.01 (the scheme's error is O(h²)), on a
// box of 60/Z.
func Extrapolated(z float64) (e1, e2 float64) {
	rmax := 60 / z
	a, b := NewLadder(z, 0.02, rmax), NewLadder(z, 0.01, rmax)
	e1 = (4*b.Eigen(0) - a.Eigen(0)) / 3
	e2 = (4*b.Eigen(1) - a.Eigen(1)) / 3
	return
}

// TwinResult is the kronspice network's ground state against the matrix.
type TwinResult struct {
	Nodes           int
	Circuit, Matrix float64 // E_1s [Ha]
}

// Network builds the ladder as a circuit for kronspice's AC analysis, with
// a unit AC current source into the node nearest r = 1, and returns it
// with the shift V0. It is the dual of Kron's network (capacitors and
// inductors exchanged), as in examples/kron_schrodinger, because AC starts
// from a DC operating point and the all-inductor network has loops of
// shorts at DC: series capacitors 1/(2h²) between neighbours (and to
// ground at both ends), shunt capacitors V_i + V0 and unit shunt
// inductors. Its admittance times jω is ω²(H + V0) − 1, which vanishes at
// ω² = 1/(E + V0): each state is a resonance, the ground state the highest.
func (l Ladder) Network() (*circuit.Circuit, string, float64, error) {
	n := len(l.D)
	// V0 makes every shunt capacitor positive: the most negative potential
	// is at the first node.
	v0 := 0.0
	for _, d := range l.D {
		v0 = math.Max(v0, -(d - 1/(l.H*l.H)))
	}
	v0++
	b := circuit.NewBuilder("Kron radial network (dual form)")
	node := func(i int) string {
		if i < 0 || i >= n {
			return "0"
		}
		return "n" + strconv.Itoa(i)
	}
	cs := 1 / (2 * l.H * l.H)
	for i := 0; i < n; i++ {
		b.Inductor("LB"+strconv.Itoa(i), node(i), "0", 1)
		b.Capacitor("CV"+strconv.Itoa(i), node(i), "0", (l.D[i]-1/(l.H*l.H))+v0)
		b.Capacitor("CS"+strconv.Itoa(i), node(i), node(i+1), cs) // the last one is the short at r_max
	}
	b.Capacitor("CS0G", node(0), "0", cs) // the short at r = 0
	drive := node(int(1/l.H) - 1)
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "0", drive)
	ckt, err := b.BuildTransient()
	return ckt, drive, v0, err
}

// Twin solves the network in kronspice. The driving-point impedance is
// imaginary: negative above the highest resonance (energies below the
// ground state) and positive just below it, so a scan upward in E from
// −0.7Z² finds the first − → + change and bisection pins it.
func Twin(z, h, rmax float64) (TwinResult, error) { return TwinL(z, 0, h, rmax) }

// TwinL is Twin for orbital angular momentum l (its lowest state: 1s for
// l = 0, 2p for l = 1).
func TwinL(z float64, ang int, h, rmax float64) (TwinResult, error) {
	l := NewLadderL(z, ang, h, rmax)
	ckt, drive, v0, err := l.Network()
	if err != nil {
		return TwinResult{}, err
	}
	im := func(e float64) (float64, error) {
		f := 1 / math.Sqrt(e+v0) / (2 * math.Pi)
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		v, ok := res.NodeVoltage(drive)
		if !ok {
			return 0, fmt.Errorf("no node %s", drive)
		}
		return imag(v[0]) / math.Max(cmplx.Abs(v[0]), 1e-300), nil
	}
	lo := -0.7 * z * z / float64((ang+1)*(ang+1))
	prev, err := im(lo)
	if err != nil {
		return TwinResult{}, err
	}
	if prev >= 0 {
		return TwinResult{}, fmt.Errorf("scan starts above the ground state")
	}
	lo0 := lo
	hi := 1.0
	const steps = 200
	for k := 1; k <= steps; k++ {
		e := lo0 + (-lo0)*float64(k)/steps
		s, err := im(e)
		if err != nil {
			return TwinResult{}, err
		}
		if s > 0 {
			hi = e
			break
		}
		lo = e
	}
	if hi > 0 {
		return TwinResult{}, fmt.Errorf("no resonance found")
	}
	for i := 0; i < 80 && hi-lo > 1e-13; i++ {
		m := 0.5 * (lo + hi)
		s, err := im(m)
		if err != nil {
			return TwinResult{}, err
		}
		if s < 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return TwinResult{Nodes: len(l.D), Circuit: 0.5 * (lo + hi), Matrix: l.Eigen(0)}, nil
}
