// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// AtomField is one hydrogen level's field: the electron's orbital current
// and spin magnetization (by the Legendre series), plus the proton's spin
// dipole at the origin (m_I = +½, along +z; hyperfine coupling left out,
// i.e. the strong-field basis).
type AtomField struct {
	S      State
	Ser    *Series
	Proton float64 // proton moment [μB] along +z
}

// NewAtomField builds the series on r ∈ (0, 80 a0), dr = 0.01 a0, l = 1…5.
func NewAtomField(s State) AtomField {
	return AtomField{S: s, Ser: NewSeries(s.J, []int{1, 2, 3, 4, 5}, 0.01, 8000, 48), Proton: MuPr}
}

// Cyl returns B_ρ, B_z [B0] at (ρ, z).
func (f AtomField) Cyl(rho, z float64) (bp, bz float64) {
	r, th := math.Hypot(rho, z), math.Atan2(rho, z)
	br, bt := f.Ser.Br(r, th), f.Ser.Bt(r, th)
	c, s := math.Cos(th), math.Sin(th)
	bp, bz = br*s+bt*c, br*c-bt*s
	a, b, _ := Dipole(f.Proton, 0, rho, z)
	return bp + a, bz + b
}

// Br is the radial field on a sphere of radius r at polar angle θ.
func (f AtomField) Br(r, th float64) float64 {
	return f.Ser.Br(r, th) + 2*f.Proton*math.Cos(th)/(r*r*r)
}

// Flux is the flux function ρA_φ.
func (f AtomField) Flux(rho, z float64) float64 {
	r, th := math.Hypot(rho, z), math.Atan2(rho, z)
	_, _, p := Dipole(f.Proton, 0, rho, z)
	return f.Ser.Flux(r, th) + p
}

// Dipole and octupole moments, the proton included in the dipole.
func (f AtomField) Moments() (m1, m3 float64) {
	return f.Ser.Moment(1) + f.Proton, f.Ser.Moment(3)
}

// OctupoleRadius is where the exterior field's pattern on a sphere changes
// from four zones (octupole-led) to two (dipole-led). Outside all sources
// B_r = 2m1 P1/r³ + 4m3 P3/r⁵ + … ∝ cosθ[1 − 3u + 5u cos²θ] with
// u = m3/(m1 r²); the two extra zeros exist for u > 1/3 or u < −1/2, so
// r_c² = 3m3/m1 or −2m3/m1. It is a statement about the exterior expansion
// only: where r_c lies inside the electron cloud the interior field decides.
func OctupoleRadius(m1, m3 float64) float64 {
	q := m3 / m1
	if math.Abs(m3) < 1e-6 {
		return 0
	}
	if q > 0 {
		return math.Sqrt(3 * q)
	}
	return math.Sqrt(-2 * q)
}

// Zones counts the N and S bands of B_r on a sphere: sign changes of
// B_r(θ) on (0, π), plus one.
func Zones(br func(th float64) float64, n int) int {
	v := make([]float64, n)
	for i := range v {
		v[i] = br(math.Pi * (float64(i) + 0.5) / float64(n))
	}
	return SignChanges(v) + 1
}

// ZoneScan counts zones on spheres of the given radii.
func ZoneScan(br func(r, th float64) float64, radii []float64, n int) []int {
	out := make([]int, len(radii))
	parallel(len(radii), func(i int) {
		r := radii[i]
		out[i] = Zones(func(th float64) float64 { return br(r, th) }, n)
	})
	return out
}

// LogRadii are n radii spaced evenly in log between lo and hi.
func LogRadii(lo, hi float64, n int) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = lo * math.Pow(hi/lo, float64(i)/float64(n-1))
	}
	return out
}

// Bands summarizes a zone scan as runs: [r_from, r_to, zones].
type Band struct {
	From, To float64
	Zones    int
}

func Bands(radii []float64, z []int) []Band {
	var out []Band
	for i := range z {
		if len(out) > 0 && out[len(out)-1].Zones == z[i] {
			out[len(out)-1].To = radii[i]
			continue
		}
		out = append(out, Band{radii[i], radii[i], z[i]})
	}
	return out
}

// LineScan is B_z along a line parallel to the axis at distance rho —
// the octave bar's bench scan (examples/octave_magnet) applied to the atom.
func LineScan(bz func(rho, z float64) float64, rho, zmax float64, n int) (zs, b []float64) {
	zs, b = make([]float64, n), make([]float64, n)
	parallel(n, func(i int) {
		z := -zmax + 2*zmax*float64(i)/float64(n-1)
		zs[i], b[i] = z, bz(rho, z)
	})
	return
}
