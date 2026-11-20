// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/em"
)

const (
	mu0   = 4e-7 * math.Pi
	rhoCu = 1.72e-8 // copper at 20 °C [Ω·m]
	cLite = 299792458.0
)

// Helix is a single-layer solenoid as a stack of coaxial circular turns:
// Turns turns of radius R (centre of the wire) at axial pitch Pitch, the
// first at Z0. Every turn of a uniform helix sees its neighbours at whole
// multiples of the pitch, so its self inductance is a sum over offsets.
type Helix struct {
	Name  string
	Turns int
	R     float64 // turn radius to the wire centre [m]
	Wire  float64 // copper diameter [m]
	Pitch float64 // centre-to-centre axial spacing [m]
	Z0    float64 // axial position of the first turn [m]
}

// Length is the winding length [m].
func (h Helix) Length() float64 { return float64(h.Turns) * h.Pitch }

// WireLength is the copper length [m].
func (h Helix) WireLength() float64 { return float64(h.Turns) * 2 * math.Pi * h.R }

// z is the axial position of turn k.
func (h Helix) z(k int) float64 { return h.Z0 + float64(k)*h.Pitch }

// loopSelf is the self inductance of one circular turn of radius r made of
// round wire of radius a, current on the wire surface (skin effect):
// μ₀r(ln(8r/a) − 2) (Rosa & Grover, Bull. Bur. Stand. 8, 1, 1912, eq. 69
// with the internal-inductance term dropped).
func loopSelf(r, a float64) float64 { return mu0 * r * (math.Log(8*r/a) - 2) }

// SelfL is the helix's inductance: N loop self terms plus every pair of
// turns by Maxwell's coaxial-circle formula (em.MutualCircles), which for a
// uniform helix depends only on the turn offset d:
// L = N·L_loop + 2·Σ_{d=1}^{N−1} (N − d)·M(R, R, d·p).
// The same filament sums give every mutual below, so the inductance matrix
// of the three coils is one consistent (positive-definite) Neumann matrix.
func (h Helix) SelfL() float64 {
	n := h.Turns
	L := float64(n) * loopSelf(h.R, h.Wire/2)
	for d := 1; d < n; d++ {
		L += 2 * float64(n-d) * em.MutualCircles(h.R, h.R, float64(d)*h.Pitch)
	}
	return L
}

// Mutual is the mutual inductance between two helices (no shared turns),
// turn by turn.
func Mutual(a, b Helix) float64 {
	m := 0.0
	for i := 0; i < a.Turns; i++ {
		for j := 0; j < b.Turns; j++ {
			m += em.MutualCircles(a.R, b.R, math.Abs(a.z(i)-b.z(j)))
		}
	}
	return m
}

// Resistance is the copper resistance at frequency f: DC value raised by
// the skin effect, R_dc·max(1, d/4δ) (Terman 1943, sec. 2). Proximity
// effect between turns is not included, so the real AC resistance of a
// close-wound coil is higher (typically 1.5–3×); declared limit.
func (h Helix) Resistance(f float64) float64 {
	area := math.Pi * h.Wire * h.Wire / 4
	rdc := rhoCu * h.WireLength() / area
	if f <= 0 {
		return rdc
	}
	delta := math.Sqrt(rhoCu / (math.Pi * f * mu0))
	return rdc * math.Max(1, h.Wire/(4*delta))
}

// SelfC is the coil's self-capacitance by Medhurst's empirical formula
// (Wireless Engineer 24, 35 and 80, 1947), in the fit quoted by Knight
// (2016, §5): C[pF] = D[cm]·(0.1126·l/D + 0.08 + 0.27/√(l/D)), for a
// single-layer solenoid of diameter D and length l away from ground.
func (h Helix) SelfC() float64 {
	D := 2 * h.R
	x := h.Length() / D
	return 1e-12 * D * 100 * (0.1126*x + 0.08 + 0.27/math.Sqrt(x))
}

// Design is one candidate: the three coils and their electrical
// surroundings. Every evolved quantity is a physical choice someone can
// build (turn counts, a spacing, capacitor and load values, the HV drive's
// frequency and phase).
type Design struct {
	N3    int     // output coil turns
	N2    int     // Royer coil turns per strand (bifilar, centre-tapped)
	Gap23 float64 // radial clearance between coil 3's wire and coil 2's [m]
	N1s   int     // HV secondary turns
	N1p   int     // HV primary turns
	C2    float64 // Royer tank capacitor [F]
	C3    float64 // output-coil resonating capacitor [F]
	RL    float64 // load [Ω]
	FHV   float64 // HV generator frequency [Hz]
	PhiHV float64 // HV generator phase at t = 0 [deg]
	// Synchronized mode: the HV amplifier is clocked by the Royer's tank
	// (Sync), with a phase lag SyncPhase [deg, 0–180], a sign SyncSign
	// (±1, so 0–360° overall) and a voltage gain SyncGain.
	Sync      bool
	SyncPhase float64
	SyncSign  int
	SyncGain  float64
	// VgAmp overrides the free-running generator's amplitude when > 0 [V].
	VgAmp float64
	// Mode switches for the comparison runs.
	RoyerOff, HVOff bool
}

// Fixed construction (declared): formers, wires, pitches.
const (
	r3Former = 25e-3  // coil 3 on a Ø50 mm former
	wire3    = 1.0e-3 // coil 3: 1.0 mm enamelled copper (AWG 18)
	pitch3   = 1.1e-3
	wire2    = 2.0e-3 // coil 2: 2.0 mm enamelled copper (AWG 12), two strands side by side
	pitch2   = 2.2e-3 // between neighbouring conductors; one strand's own pitch is 2·pitch2
	r1Former = 37.5e-3
	wire1s   = 0.4e-3 // HV secondary: 0.4 mm enamelled (AWG 26)
	pitch1s  = 0.45e-3
	wire1p   = 2.0e-3 // HV primary: 2.0 mm, 10 mm outside the secondary at its grounded end
	pitch1p  = 3.0e-3
)

// Coils returns the four windings of a design: coil 1's primary and
// secondary, coil 2's two strands (A, B) and coil 3. Coil 2 is centred
// axially over coil 3.
func (d Design) Coils() (c1p, c1s, c2a, c2b, c3 Helix) {
	c3 = Helix{Name: "coil 3", Turns: d.N3, R: r3Former + wire3/2, Wire: wire3, Pitch: pitch3}
	r2 := r3Former + wire3 + d.Gap23 + wire2/2
	l2 := float64(2*d.N2) * pitch2
	z2 := 0.5*(c3.Length()-l2) + 0.5*pitch2
	c2a = Helix{Name: "coil 2A", Turns: d.N2, R: r2, Wire: wire2, Pitch: 2 * pitch2, Z0: z2}
	c2b = Helix{Name: "coil 2B", Turns: d.N2, R: r2, Wire: wire2, Pitch: 2 * pitch2, Z0: z2 + pitch2}
	c1s = Helix{Name: "coil 1 secondary", Turns: d.N1s, R: r1Former + wire1s/2, Wire: wire1s, Pitch: pitch1s}
	c1p = Helix{Name: "coil 1 primary", Turns: d.N1p, R: r1Former + 10e-3, Wire: wire1p, Pitch: pitch1p}
	return
}

// Electrical is the lumped description the circuit is built from.
type Electrical struct {
	L1p, L1s, M1, C1self       float64
	L2a, L2b, M2ab, M2a3, M2b3 float64
	L3                         float64
	R1p, R1s, R2, R3           float64 // at the design frequency
	K1, K2ab, K2a3, K2b3       float64
}

// Electrical computes inductances, couplings and resistances at frequency f.
func (d Design) Electrical(f float64) Electrical {
	c1p, c1s, c2a, c2b, c3 := d.Coils()
	var e Electrical
	e.L1p, e.L1s, e.M1 = c1p.SelfL(), c1s.SelfL(), Mutual(c1p, c1s)
	e.C1self = c1s.SelfC()
	e.L2a, e.L2b, e.M2ab = c2a.SelfL(), c2b.SelfL(), Mutual(c2a, c2b)
	e.L3 = c3.SelfL()
	e.M2a3, e.M2b3 = Mutual(c2a, c3), Mutual(c2b, c3)
	e.R1p, e.R1s, e.R2, e.R3 = c1p.Resistance(f), c1s.Resistance(f), c2a.Resistance(f), c3.Resistance(f)
	e.K1 = e.M1 / math.Sqrt(e.L1p*e.L1s)
	e.K2ab = e.M2ab / math.Sqrt(e.L2a*e.L2b)
	e.K2a3 = e.M2a3 / math.Sqrt(e.L2a*e.L3)
	e.K2b3 = e.M2b3 / math.Sqrt(e.L2b*e.L3)
	return e
}
