// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// Th. Kaluza, "Zum Unitätsproblem der Physik," Sitzungsber. Preuss. Akad.
// Wiss. Berlin (Math. Phys.) 1921, 966–972 (English: V. T. Toth's revised
// translation). Latin indices 0–4, Greek 1–4; x⁰ is the fifth coordinate.

// Metric5 returns a 5×5 metric at point x (x[0] is x⁰).
type Metric5 func(x []float64) [5][5]float64

// KaluzaMetric is Kaluza's ansatz (eq. 2): the 4-D metric g_κλ framed by
// the electromagnetic potential, g_0κ = 2α·q_κ, g_00 = 2g, with nothing
// depending on x⁰ (the cylinder condition).
func KaluzaMetric(g4 func(x []float64) [4][4]float64, q func(x []float64) [4]float64, g func(x []float64) float64, alpha float64) Metric5 {
	return func(x []float64) [5][5]float64 {
		var m [5][5]float64
		gg, qq := g4(x), q(x)
		for k := 0; k < 4; k++ {
			for l := 0; l < 4; l++ {
				m[k+1][l+1] = gg[k][l]
			}
			m[0][k+1], m[k+1][0] = 2*alpha*qq[k], 2*alpha*qq[k]
		}
		m[0][0] = 2 * g(x)
		return m
	}
}

// derivative returns ∂g_ab/∂x^c by central differences (exact for metrics
// quadratic in x; h relative to the coordinate scale).
func derivative(m Metric5, x []float64, h float64) [5][5][5]float64 {
	var d [5][5][5]float64
	xp := make([]float64, 5)
	for c := 0; c < 5; c++ {
		copy(xp, x)
		xp[c] = x[c] + h
		mp := m(xp)
		xp[c] = x[c] - h
		mm := m(xp)
		for a := 0; a < 5; a++ {
			for b := 0; b < 5; b++ {
				d[a][b][c] = (mp[a][b] - mm[a][b]) / (2 * h)
			}
		}
	}
	return d
}

// KaluzaGamma is Kaluza's three-index symbol in his sign convention
// (eq. 1): 2Γ_ikl = g_ik,l − g_kl,i − g_li,k, i.e. Γ_ikl = −{ik, l}.
func KaluzaGamma(m Metric5, x []float64, h float64) [5][5][5]float64 {
	d := derivative(m, x, h)
	var G [5][5][5]float64
	for i := 0; i < 5; i++ {
		for k := 0; k < 5; k++ {
			for l := 0; l < 5; l++ {
				G[i][k][l] = 0.5 * (d[i][k][l] - d[k][l][i] - d[l][i][k])
			}
		}
	}
	return G
}

// ---------------------------------------------------------------------------
// The geodesic in five dimensions, exactly (beyond Kaluza's approximations).
//
// Modern form of the same ansatz with the corner fixed: coordinates
// (t, x, y, z, w), signature (−, +, +, +, +), units c = G = 1 with
// Gaussian charge,
//
//	ĝ_μν = η_μν + k²·A_μ·A_ν,   ĝ_μw = k·A_μ,   ĝ_ww = 1,   k = 2√G
//
// k is fixed by asking the 5-D curvature to contain Maxwell's Gaussian
// Lagrangian −F²/16π beside Einstein's R/16πG: k²/4 = G. Kaluza's
// 2α = 2·√(κ/2) with κ = 8πG is the same coupling in Heaviside units.
// Nothing depends on w, so p_w = ĝ_wB·u^B is conserved, and the 4-D motion
// is the Lorentz force with
//
//	q/m = k·p_w / √(1 + p_w²)
//
// (Kaluza's eq. 8a, ρ0/μ0 = 2α·u⁰, is its small-p_w form.)

// UniformB is the KK metric for a uniform magnetic field B along z
// (A = (0, −By/2, Bx/2, 0)).
func UniformB(B, k float64) Metric5 {
	return func(x []float64) [5][5]float64 {
		A := [4]float64{0, -B * x[2] / 2, B * x[1] / 2, 0}
		var m [5][5]float64
		eta := [4]float64{-1, 1, 1, 1}
		for a := 0; a < 4; a++ {
			for b := 0; b < 4; b++ {
				m[a][b] = k * k * A[a] * A[b]
			}
			m[a][a] += eta[a]
			m[a][4], m[4][a] = k*A[a], k*A[a]
		}
		m[4][4] = 1
		return m
	}
}

// christoffel2 returns Γ^a_bc for the metric at x.
func christoffel2(m Metric5, x []float64) [5][5][5]float64 {
	d := derivative(m, x, 1e-3)
	g := m(x)
	inv := invert5(g)
	var G [5][5][5]float64
	for a := 0; a < 5; a++ {
		for b := 0; b < 5; b++ {
			for c := 0; c < 5; c++ {
				s := 0.0
				for e := 0; e < 5; e++ {
					s += inv[a][e] * (d[e][b][c] + d[e][c][b] - d[b][c][e])
				}
				G[a][b][c] = 0.5 * s
			}
		}
	}
	return G
}

func invert5(a [5][5]float64) [5][5]float64 {
	var m [5][10]float64
	for i := 0; i < 5; i++ {
		copy(m[i][:5], a[i][:])
		m[i][5+i] = 1
	}
	for c := 0; c < 5; c++ {
		p := c
		for r := c + 1; r < 5; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		m[c], m[p] = m[p], m[c]
		f := m[c][c]
		for k := range m[c] {
			m[c][k] /= f
		}
		for r := 0; r < 5; r++ {
			if r == c || m[r][c] == 0 {
				continue
			}
			f := m[r][c]
			for k := range m[r] {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	var out [5][5]float64
	for i := 0; i < 5; i++ {
		copy(out[i][:], m[i][5:])
	}
	return out
}

// Orbit is a geodesic run in the uniform field.
type Orbit struct {
	Pw               float64 // conserved momentum along the fifth dimension
	PwDrift          float64 // its largest relative change along the run
	NormDrift        float64 // drift of ĝ_AB u^A u^B
	QoverM           float64 // charge-to-mass predicted, k·p_w/√(1+p_w²)
	Radius, RadiusLF float64 // measured orbit radius, Lorentz-force radius γmv/(qB)
	Omega, OmegaLF   float64 // measured angular frequency in coordinate time, qB/(γm)
	X, Y             []float64
}

// RunOrbit integrates the 5-D geodesic (RK4 in the affine parameter) for a
// particle starting at the origin with speed v along x and momentum pw
// along w, in field B, for about two turns.
func RunOrbit(B, k, v, pw float64) Orbit {
	m := UniformB(B, k)
	gam := 1 / math.Sqrt(1-v*v)
	// timelike: ĝ u u = −1; at the origin A = 0, so u_w = u^w = pw and
	// −(u^t)² + (u^x)² + pw² = −1 with u^x = v u^t.
	ut := math.Sqrt((1 + pw*pw) / (1 - v*v))
	x := []float64{0, 0, 0, 0, 0}
	u := []float64{ut, v * ut, 0, 0, pw}
	o := Orbit{Pw: pw}
	o.QoverM = k * pw / math.Sqrt(1+pw*pw)
	o.OmegaLF = o.QoverM * B / gam
	o.RadiusLF = gam * v / (o.QoverM * B)
	acc := func(x, u []float64) []float64 {
		G := christoffel2(m, x)
		a := make([]float64, 5)
		for i := 0; i < 5; i++ {
			for b := 0; b < 5; b++ {
				for c := 0; c < 5; c++ {
					a[i] -= G[i][b][c] * u[b] * u[c]
				}
			}
		}
		return a
	}
	pwOf := func(x, u []float64) float64 {
		g := m(x)
		s := 0.0
		for b := 0; b < 5; b++ {
			s += g[4][b] * u[b]
		}
		return s
	}
	norm := func(x, u []float64) float64 {
		g := m(x)
		s := 0.0
		for a := 0; a < 5; a++ {
			for b := 0; b < 5; b++ {
				s += g[a][b] * u[a] * u[b]
			}
		}
		return s
	}
	period := 2 * math.Pi / o.OmegaLF // coordinate time
	ds := period / ut / 2000
	steps := 4000
	minY, maxY, minX, maxX := 0.0, 0.0, 0.0, 0.0
	angle, prev := 0.0, math.Atan2(u[2], u[1]) // direction of the velocity, unwrapped
	tmpx, tmpu := make([]float64, 5), make([]float64, 5)
	for s := 0; s < steps; s++ {
		k1u := acc(x, u)
		k1x := append([]float64(nil), u...)
		for i := range x {
			tmpx[i], tmpu[i] = x[i]+ds/2*k1x[i], u[i]+ds/2*k1u[i]
		}
		k2u, k2x := acc(tmpx, tmpu), append([]float64(nil), tmpu...)
		for i := range x {
			tmpx[i], tmpu[i] = x[i]+ds/2*k2x[i], u[i]+ds/2*k2u[i]
		}
		k3u, k3x := acc(tmpx, tmpu), append([]float64(nil), tmpu...)
		for i := range x {
			tmpx[i], tmpu[i] = x[i]+ds*k3x[i], u[i]+ds*k3u[i]
		}
		k4u, k4x := acc(tmpx, tmpu), append([]float64(nil), tmpu...)
		for i := range x {
			x[i] += ds / 6 * (k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i])
			u[i] += ds / 6 * (k1u[i] + 2*k2u[i] + 2*k3u[i] + k4u[i])
		}
		o.PwDrift = math.Max(o.PwDrift, math.Abs(pwOf(x, u)-pw)/math.Max(math.Abs(pw), 1e-300))
		o.NormDrift = math.Max(o.NormDrift, math.Abs(norm(x, u)+1))
		a := math.Atan2(u[2], u[1])
		d := a - prev
		d -= 2 * math.Pi * math.Round(d/(2*math.Pi))
		angle += d
		prev = a
		minX, maxX = math.Min(minX, x[1]), math.Max(maxX, x[1])
		minY, maxY = math.Min(minY, x[2]), math.Max(maxY, x[2])
		if s%10 == 0 {
			o.X, o.Y = append(o.X, x[1]), append(o.Y, x[2])
		}
	}
	o.Radius = 0.25 * ((maxX - minX) + (maxY - minY))
	o.Omega = math.Abs(angle) / x[0] // turned angle ÷ coordinate time
	return o
}
