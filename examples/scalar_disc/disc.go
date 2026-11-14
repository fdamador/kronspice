// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/fdamador/kronspice/em"
)

// Disc is the optically gated rotating aluminium disc of Iadicicco &
// Verolino (arXiv:2608.21403): a fibre ring on the disc lights one angular
// sector, which induces a surface charge density in the metal under it;
// the sector turns with the disc. The paper gives no dimensions, speeds or
// charge, so every number here is declared.
type Disc struct {
	RingR   float64 // radius of the lit ring [m]
	RingW   float64 // radial width of the lit ring [m]
	FWHMDeg float64 // angular width of the lit sector, full width at half maximum [deg]
	Q       float64 // charge induced in the sector [C]
	Omega   float64 // angular velocity [rad/s], sign = sense of rotation
	Uniform bool    // light the whole ring (control: axisymmetric source)
	NR, NPh int     // quadrature points across and around the ring
}

// DefaultDisc is the declared bench: a 4 cm lit ring 5 mm wide on a 5 cm
// disc, a 15° sector holding 1 nC, 3000 rpm.
func DefaultDisc() Disc {
	return Disc{RingR: 0.04, RingW: 0.005, FWHMDeg: 15, Q: 1e-9, Omega: 3000 * 2 * math.Pi / 60, NR: 6, NPh: 1440}
}

// Sigma is the sector's Gaussian angular width [rad].
func (d Disc) Sigma() float64 {
	return d.FWHMDeg * math.Pi / 180 / (2 * math.Sqrt(2*math.Ln2))
}

// srcPoint is one quadrature element of the ring.
type srcPoint struct {
	x, y   float64
	area   float64 // [m²]
	sigma  float64 // surface charge density [C/m²]
	source float64 // observer-split source I_G = Ω ∂ρ/∂φ, per area [C/(m² s)]
}

// Sources lays the ring out with the sector centred on φc. The surface
// density is σ_s = Q g(φ − φc)/(w r A_g), with g the Gaussian of width σ
// and A_g its integral, and the source is Iadicicco, Modanese & Verolino's
// weak-rigid-rotation split term (GRG 58, 106, 2026, eqs. 13, 23):
//
//	I_G = Ω ∂_φ ρ = −Ω (φ − φc)/σ² · σ_s
//
// whose peak over φ is |Ω| σ_s,max e^(−½)/σ (their eq. B3) and whose ring
// integral is zero (their eq. 24).
func (d Disc) Sources(phiC float64) []srcPoint {
	sg := d.Sigma()
	dr := d.RingW / float64(d.NR)
	dph := 2 * math.Pi / float64(d.NPh)
	var norm float64 // Σ g dφ
	g := make([]float64, d.NPh)
	u := make([]float64, d.NPh)
	for j := 0; j < d.NPh; j++ {
		phi := (float64(j) + 0.5) * dph
		x := math.Remainder(phi-phiC, 2*math.Pi)
		u[j] = x
		if d.Uniform {
			g[j] = 1
		} else {
			g[j] = math.Exp(-x * x / (2 * sg * sg))
		}
		norm += g[j] * dph
	}
	var out []srcPoint
	for i := 0; i < d.NR; i++ {
		r := d.RingR - d.RingW/2 + (float64(i)+0.5)*dr
		for j := 0; j < d.NPh; j++ {
			phi := (float64(j) + 0.5) * dph
			s := d.Q * g[j] / (d.RingW * r * norm) // ∫σ_s r dr dφ = Q
			var src float64
			if !d.Uniform {
				src = -d.Omega * u[j] / (sg * sg) * s
			}
			out = append(out, srcPoint{x: r * math.Cos(phi), y: r * math.Sin(phi), area: r * dr * dph, sigma: s, source: src})
		}
	}
	return out
}

// soft is the softening length of the point kernels [m] (a quarter of the
// ring's quadrature spacing), so field points on the disc stay finite.
func (d Disc) soft() float64 { return 0.25 * d.RingW / float64(d.NR) }

// SAt is S at x in the quasi-static limit of □S = μ₀I (em/scalar): the
// source changes on the rotation time, about 10⁸ times the light time
// across the disc, so S(x) = (μ₀/4π) Σ I_k A_k / |x − x_k|.
func SAt(src []srcPoint, x [3]float64, soft float64) float64 {
	var s float64
	for _, p := range src {
		dx, dy, dz := x[0]-p.x, x[1]-p.y, x[2]
		s += p.source * p.area / math.Sqrt(dx*dx+dy*dy+dz*dz+soft*soft)
	}
	return em.Mu0 / (4 * math.Pi) * s
}

// EAt is the ordinary Coulomb field of the induced charge at x [V/m] (the
// disc's own image charges are left out: an order-of-magnitude systematic).
func EAt(src []srcPoint, x [3]float64, soft float64) [3]float64 {
	var e [3]float64
	k := 1 / (4 * math.Pi * em.Eps0)
	for _, p := range src {
		dx, dy, dz := x[0]-p.x, x[1]-p.y, x[2]
		r2 := dx*dx + dy*dy + dz*dz + soft*soft
		f := k * p.sigma * p.area / (r2 * math.Sqrt(r2))
		e[0] += f * dx
		e[1] += f * dy
		e[2] += f * dz
	}
	return e
}

// Grid is S sampled on a box around the disc, for the gravity source S².
type Grid struct {
	Lo   [3]float64
	H    float64
	N    [3]int
	S    []float64
	dV   float64
	soft float64
}

// SGrid samples S on a cube of half-width half [m] with spacing h,
// in parallel.
func SGrid(src []srcPoint, half, h, soft float64) *Grid {
	n := int(math.Round(2*half/h)) + 1
	g := &Grid{Lo: [3]float64{-half, -half, -half}, H: h, N: [3]int{n, n, n}, S: make([]float64, n*n*n), dV: h * h * h, soft: h / 2}
	w := runtime.NumCPU()
	var wg sync.WaitGroup
	for p := 0; p < w; p++ {
		wg.Add(1)
		go func(p int) {
			defer wg.Done()
			for k := p; k < n; k += w {
				for j := 0; j < n; j++ {
					for i := 0; i < n; i++ {
						x := [3]float64{g.Lo[0] + float64(i)*h, g.Lo[1] + float64(j)*h, g.Lo[2] + float64(k)*h}
						g.S[i+n*(j+n*k)] = SAt(src, x, soft)
					}
				}
			}
		}(p)
	}
	wg.Wait()
	return g
}

// Coupling is the weak-field gravity coupling of Minotti & Modanese
// (EPJC 86, 593, 2026, eqs. 11–13): □χ = ∂²φ/∂t² + (1/2c²)K_φC_φS²,
// □φ = −K_φS², force per unit mass f = −∇χ. Both constants are their fit
// to the Podkletnov and Poher reports (declared, not measured).
type Coupling struct {
	Kphi float64 // [A² N⁻²]
	Cphi float64 // dimensionless
}

// MMFit is the 2026 paper's fit: K_φ ~ 10⁻¹⁶ A² N⁻², C_φ ≈ 4–9×10²¹ (6×10²¹ taken).
func MMFit() Coupling { return Coupling{Kphi: 1e-16, Cphi: 6e21} }

// Chi is χ at x from the quasi-static source term, the solution of
// −∇²χ = (K_φC_φ/2c²) S² over the grid: χ(x) = (K_φC_φ/8πc²) Σ S² dV/|x − x′|.
// The ∂²φ/∂t² term is smaller by (ΩL/c)², about 10⁻¹⁹, and left out.
func (g *Grid) Chi(x [3]float64, cp Coupling) float64 {
	c2 := em.SpeedOfLight * em.SpeedOfLight
	var s float64
	n := g.N[0]
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				v := g.S[i+n*(j+n*k)]
				if v == 0 {
					continue
				}
				dx := x[0] - (g.Lo[0] + float64(i)*g.H)
				dy := x[1] - (g.Lo[1] + float64(j)*g.H)
				dz := x[2] - (g.Lo[2] + float64(k)*g.H)
				s += v * v / math.Sqrt(dx*dx+dy*dy+dz*dz+g.soft*g.soft)
			}
		}
	}
	return cp.Kphi * cp.Cphi / (8 * math.Pi * c2) * s * g.dV
}

// Force is f = −∇χ at x by central differences of step dx [m/s²].
func (g *Grid) Force(x [3]float64, cp Coupling, dx float64) [3]float64 {
	var f [3]float64
	for a := 0; a < 3; a++ {
		xp, xm := x, x
		xp[a] += dx
		xm[a] -= dx
		f[a] = -(g.Chi(xp, cp) - g.Chi(xm, cp)) / (2 * dx)
	}
	return f
}
