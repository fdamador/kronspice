// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/rand/v2"
)

// Hydrogen under pressure: "high pressure, short waves" (Walter Russell's
// octave drawing) read as the confined hydrogen atom of A. Michels,
// J. de Boer, A. Bijl, "Remarks concerning molecular interaction and their
// influence on the polarisability," Physica 4, 981 (1937). In Kron's
// network the outer short circuit (ψ = 0) is a hard wall; moving it in to
// radius r_c compresses the atom. The levels rise, the lines shorten, the
// 2s–2p degeneracy of the free atom splits, and the ground state reaches
// zero binding at r_c ≈ 1.835 a0 (pressure ionization). The pressure the
// atom exerts on the wall is P = −dE/dV, V = (4/3)πr_c³.

// ConfinedLevel is the k-th level of angular momentum l inside a wall at
// rc [a0, Ha; reduced-mass atomic units], Richardson-extrapolated from
// 2000 and 4000 network sections.
func ConfinedLevel(z float64, ang, k int, rc float64) float64 {
	a := NewLadderL(z, ang, rc/2000, rc)
	b := NewLadderL(z, ang, rc/4000, rc)
	return (4*b.Eigen(k) - a.Eigen(k)) / 3
}

// Confined is the atom inside one wall.
type Confined struct {
	Rc            float64 // wall radius [a0]
	E1s, E2s, E2p float64 // [Ha]
	PressureGPa   float64 // on the wall, from the ground state
	LymanNm       float64 // 2p → 1s
}

// AU is the reduced-mass atomic unit of length [m].
var AU = A0 * Recoil

// Compress computes the atom inside a wall of radius rc.
func Compress(rc float64) Confined {
	c := Confined{Rc: rc, E1s: ConfinedLevel(1, 0, 0, rc), E2s: ConfinedLevel(1, 0, 1, rc), E2p: ConfinedLevel(1, 1, 0, rc)}
	d := 1e-3 * rc
	dE := (ConfinedLevel(1, 0, 0, rc+d) - ConfinedLevel(1, 0, 0, rc-d)) / (2 * d) // Ha per a0
	pAU := -dE / (4 * math.Pi * rc * rc)                                          // Ha per a0³
	c.PressureGPa = pAU * HaJ / (AU * AU * AU) / 1e9
	c.LymanNm = nm(c.E2p - c.E1s)
	return c
}

// Walls are the radii of the pressure table and figure.
var Walls = []float64{60, 12, 8, 6, 5, 4, 3.5, 3, 2.5, 2.2, 2, 1.9, 1.835, 1.7, 1.6, 1.5}

// CloudPoints2p samples n points of the 2p_z density r²e^(−r)cos²θ (in a0):
// r ~ r⁴e^(−r) (a sum of five exponentials), cos θ ~ cos²θ (|c| = u^(1/3)).
func CloudPoints2p(n int) [][3]float64 {
	rng := rand.New(rand.NewPCG(2, 64))
	pts := make([][3]float64, n)
	for i := range pts {
		var r float64
		for k := 0; k < 5; k++ {
			r -= math.Log(1 - rng.Float64())
		}
		c := math.Cbrt(rng.Float64())
		if rng.Float64() < 0.5 {
			c = -c
		}
		ph := 2 * math.Pi * rng.Float64()
		s := math.Sqrt(1 - c*c)
		pts[i] = [3]float64{r * s * math.Cos(ph), r * s * math.Sin(ph), r * c}
	}
	return pts
}
