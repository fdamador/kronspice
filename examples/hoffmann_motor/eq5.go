// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/kron"
)

// solve4 solves the 4×4 system a·x = b by Gaussian elimination with
// partial pivoting (a and b are copied).
func solve4(a [][]float64, b []float64) []float64 {
	var m [4][5]float64
	for i := 0; i < 4; i++ {
		copy(m[i][:4], a[i])
		m[i][4] = b[i]
	}
	for c := 0; c < 4; c++ {
		p := c
		for r := c + 1; r < 4; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		m[c], m[p] = m[p], m[c]
		for r := c + 1; r < 4; r++ {
			f := m[r][c] / m[c][c]
			for k := c; k < 5; k++ {
				m[r][k] -= f * m[c][k]
			}
		}
	}
	x := make([]float64, 4)
	for r := 3; r >= 0; r-- {
		v := m[r][4]
		for k := r + 1; k < 4; k++ {
			v -= m[r][k] * x[k]
		}
		x[r] = v / m[r][r]
	}
	return x
}

// Eq5Run integrates Hoffmann's eq. 5′ itself, in the stationary-coil
// frame, with fourth-order Runge–Kutta:
//
//	L̄_λν·ẍ^ν = F̄_λ − R̄_λμ·ẋ^μ − 2·Ω_μλν·ẋ^μ·ẋ^ν
//
// with L̄ from eq. 15 (constant, eq. 29), Ω from eq. 19, R̄ = R (eq. 40:
// equal rotor resistances), F̄ = (supply voltage, 0, 0, −friction·ω). It
// shares no code with kronspice's N element, so it is an independent check
// of the frame-I and frame-II(b) runs.
func (p Params) Eq5Run(omega0, seconds float64) (Run, error) {
	s := p.FrameI()
	f := Frame()
	lb, _, err := kron.NewMetric(s, p.J, f, 0)
	if err != nil {
		return Run{}, err
	}
	r := []float64{p.Rs, p.Rr, p.Rr, p.Friction}
	// state: i1, i2̄, i3̄, ω, θ
	// In this frame Ω does not depend on θ (checked here at eight angles;
	// TestTorsionConstant), so it is computed once.
	om, err := kron.Torsion(s, p.J, f, 0)
	if err != nil {
		return Run{}, err
	}
	if d := torsionSpread(p); d > 1e-12 {
		return Run{}, fmt.Errorf("torsion varies with θ by %g", d)
	}
	deriv := func(t float64, y []float64) ([]float64, error) {
		v := []float64{y[0], y[1], y[2], y[3]}
		rhs := make([]float64, 4)
		for la := 0; la < 4; la++ {
			x := -r[la] * v[la]
			if la == 0 {
				x += p.Vpeak * math.Sin(2*math.Pi*p.F*t)
			}
			for mu := 0; mu < 4; mu++ {
				for nu := 0; nu < 4; nu++ {
					x -= 2 * om[mu][la][nu] * v[mu] * v[nu]
				}
			}
			rhs[la] = x
		}
		a := solve4(lb, rhs)
		return []float64{a[0], a[1], a[2], a[3], y[3]}, nil
	}
	dt := 1 / (p.F * 2000)
	steps := int(math.Round(seconds / dt))
	y := []float64{0, 0, 0, omega0, 0}
	run := Run{Name: "Hoffmann eq. 5′ (own Runge–Kutta)"}
	tmp := make([]float64, 5)
	for k := 0; k < steps; k++ {
		t := float64(k) * dt
		k1, err := deriv(t, y)
		if err != nil {
			return run, err
		}
		for j := range y {
			tmp[j] = y[j] + dt/2*k1[j]
		}
		k2, _ := deriv(t+dt/2, tmp)
		for j := range y {
			tmp[j] = y[j] + dt/2*k2[j]
		}
		k3, _ := deriv(t+dt/2, tmp)
		for j := range y {
			tmp[j] = y[j] + dt*k3[j]
		}
		k4, _ := deriv(t+dt, tmp)
		for j := range y {
			y[j] += dt / 6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
		}
		if (k+1)%20 == 0 {
			run.Time = append(run.Time, t+dt)
			run.Speed = append(run.Speed, y[3])
		}
	}
	run.Final = y[3]
	return run, nil
}

// torsionSpread is the largest change of any Ω component over eight
// rotor angles, relative to the largest component.
func torsionSpread(p Params) float64 {
	s, f := p.FrameI(), Frame()
	ref, _ := kron.Torsion(s, p.J, f, 0)
	top, d := 0.0, 0.0
	for k := 1; k < 8; k++ {
		om, _ := kron.Torsion(s, p.J, f, float64(k)*math.Pi/4+0.1)
		for a := range om {
			for b := range om[a] {
				for c := range om[a][b] {
					top = math.Max(top, math.Abs(ref[a][b][c]))
					d = math.Max(d, math.Abs(om[a][b][c]-ref[a][b][c]))
				}
			}
		}
	}
	return d / top
}
