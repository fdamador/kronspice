// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Stage 4: the gap as an acoustic ring in kronspice. Acoustic domain
// (doc/domains.md): pressure across, volume velocity through. The annular
// gap of cross-section A = gap × band length is cut into M sections round
// the circumference: series inertance ρΔx/A, shunt compliance AΔx/(ρc²),
// closed on itself. Its azimuthal modes are exactly
// ω_m = (2/√(LC))·sin(πm/M), which tends to m·c/R (f_m = m·c/2πR) as M grows. The dimple
// ring's pressure pattern of order N travels at the rotor speed U, so it
// drives mode N at the passing frequency, and the two meet when c = U.

// GapRing builds the ring with a 1 m³/s AC volume-velocity source at node 0
// and a leak to ground for the DC path, 10⁴ × the ring's wave impedance
// so it does not load the modes (declared).
func GapRing(r, gap, bandLen, c float64, sections int, lossQ float64) (*circuit.Circuit, float64, float64, error) {
	a := gap * bandLen
	dx := 2 * math.Pi * r / float64(sections)
	// in units of the wave impedance Z = ρc/A (frequencies unchanged:
	// L/Z · C·Z = LC), so the matrix is well scaled
	z := rho * c / a
	l := rho * dx / a / z
	cc := a * dx / (rho * c * c) * z
	b := circuit.NewBuilder("gap ring")
	node := func(k int) string { return fmt.Sprintf("n%d", k%sections) }
	b.Add(circuit.Element{Name: "IDRIVE", Kind: circuit.CurrentSource, ACMag: 1}, "0", node(0))
	b.Resistor("RLEAK", node(0), "0", 1e4)
	w0 := 1 / math.Sqrt(l*cc)
	rs := w0 * l / lossQ // series loss for a declared Q
	for k := 0; k < sections; k++ {
		mid := fmt.Sprintf("m%d", k)
		b.Inductor(fmt.Sprintf("L%d", k), node(k), mid, l)
		b.Resistor(fmt.Sprintf("R%d", k), mid, node(k+1), rs)
		b.Capacitor(fmt.Sprintf("C%d", k), node(k), "0", cc)
	}
	ckt, err := b.Build()
	return ckt, l, cc, err
}

// RingModes runs an AC sweep and returns the peak frequencies of |p| at the
// driven node, with the exact lattice frequencies beside them.
func RingModes(r, gap, bandLen, c float64, sections, modes int) (peaks, exact []float64, err error) {
	ckt, l, cc, err := GapRing(r, gap, bandLen, c, sections, 2000)
	if err != nil {
		return nil, nil, err
	}
	for m := 1; m <= modes; m++ {
		exact = append(exact, 2/math.Sqrt(l*cc)*math.Sin(math.Pi*float64(m)/float64(sections))/(2*math.Pi))
	}
	fmax := exact[len(exact)-1] * 1.08
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 20000, FStart: exact[0] * 0.5, FStop: fmax}, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	v, _ := res.NodeVoltage("n0")
	for k := 1; k+1 < len(v); k++ {
		a, b0, c0 := cmplx.Abs(v[k-1]), cmplx.Abs(v[k]), cmplx.Abs(v[k+1])
		if b0 > a && b0 > c0 {
			// parabolic refinement
			den := a - 2*b0 + c0
			d := 0.0
			if den != 0 {
				d = 0.5 * (a - c0) / den
			}
			df := res.Freq[k+1] - res.Freq[k]
			peaks = append(peaks, res.Freq[k]+d*df)
		}
	}
	return peaks, exact, nil
}

// RingSpectrum returns |p| at the driven node over a linear sweep.
func RingSpectrum(r, gap, bandLen, c float64, sections int, f0, f1 float64, points int) (f, mag []float64, err error) {
	ckt, _, _, err := GapRing(r, gap, bandLen, c, sections, 2000)
	if err != nil {
		return nil, nil, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: points, FStart: f0, FStop: f1}, analysis.Options{})
	if err != nil {
		return nil, nil, err
	}
	v, _ := res.NodeVoltage("n0")
	for k := range v {
		mag = append(mag, cmplx.Abs(v[k]))
	}
	return res.Freq, mag, nil
}
