// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Fig. 5 of Barrett 1991 maps pump and signal onto a quaternion and calls
// it the idler: x = x0 + x1·i + x2·j + x3·k with x0 = Es*·e^{−iωs t},
// x1 = Es·e^{−iβs t}, x2 = √Ep·e^{iωp t}, x3 = √Ep·e^{−iβp t}. Two checks:
//
//  1. x is a sum, so each component carries only the signal or the pump
//     frequency; an idler at ωp − ωs needs a product (a nonlinear element),
//     which is what the Fig. 2A varactor run shows.
//  2. β is a propagation constant [rad/m]; β·t (also in the paper's Eq. 24
//     and the patent's Eq. 15) has units rad·s/m, not radians: the phase
//     should be β·z.
//
// FigFiveLines returns, for each component, its largest spectral line at
// fp − fs relative to its own largest line, on the Fig. 2 frequency grid.
func FigFiveLines() [4]float64 {
	p := DefaultFig2
	const m = 20000
	ws, wp := 2*math.Pi*float64(p.MS)*p.F0, 2*math.Pi*float64(p.MP)*p.F0
	comp := [4]func(t float64) float64{
		func(t float64) float64 { return math.Cos(ws * t) },              // Re x0
		func(t float64) float64 { return math.Cos(ws * t) },              // x1 taking β·t as ω·t
		func(t float64) float64 { return math.Sqrt(2) * math.Cos(wp*t) }, // x2
		func(t float64) float64 { return math.Sqrt(2) * math.Cos(wp*t) }, // x3
	}
	var out [4]float64
	T := 1 / p.F0
	for c, f := range comp {
		x := make([]float64, m)
		for j := range x {
			x[j] = f(float64(j) * T / m)
		}
		peak := math.Max(cmplx.Abs(dftReal(x, p.MS)), cmplx.Abs(dftReal(x, p.MP)))
		out[c] = cmplx.Abs(dftReal(x, p.MP-p.MS)) / peak
	}
	return out
}

func printFig5() error {
	l := FigFiveLines()
	fmt.Println("Fig. 5 — pump and signal mapped onto a quaternion 'idler' x = x0 + x1 i + x2 j + x3 k")
	fmt.Printf("  each component's line at fp − fs: %.1e %.1e %.1e %.1e (relative): a sum of the pump and the\n", l[0], l[1], l[2], l[3])
	fmt.Println("  signal holds no idler; one at fp − fs needs a product, i.e. a nonlinear element (Fig. 2A varactor).")
	fmt.Println("  x1 and x3 (and Eq. 24) use β·t: β is a propagation constant [rad/m], so β·t is not a phase;")
	fmt.Println("  the travelling-wave phase is β·z. The Fig. 2A run gives idler ∝ Es*·Ep, not Eq. 24's |Es|²·Ep.")
	fmt.Println()
	return nil
}
