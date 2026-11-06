// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"runtime"
	"sync"

	"github.com/fdamador/kronspice/em"
)

// LoopField is the exact field of a circular loop (radius a, current I,
// in the plane z = 0) at (ρ, z), in units where μ0 = 4π (so a current in
// μB/a0² gives B in B0) (Smythe, Static and Dynamic Electricity, §7.10).
// It also returns the flux function ρA_φ.
func LoopField(a, I, rho, z float64) (br, bz, psi float64) {
	const mu = 4 * math.Pi
	if rho < 1e-12*a {
		return 0, mu * I * a * a / (2 * math.Pow(a*a+z*z, 1.5)), 0
	}
	s := (a+rho)*(a+rho) + z*z
	k2 := 4 * a * rho / s
	K, E := em.EllipticKE(math.Sqrt(k2))
	d := (a-rho)*(a-rho) + z*z
	pre := mu * I / (2 * math.Pi * math.Sqrt(s))
	bz = pre * (K + (a*a-rho*rho-z*z)/d*E)
	br = pre * z / rho * (-K + (a*a+rho*rho+z*z)/d*E)
	k := math.Sqrt(k2)
	psi = rho * mu * I / (math.Pi * k) * math.Sqrt(a/rho) * ((1-k2/2)*K - E)
	return
}

// Loops is an azimuthal current density written as coaxial loops on a
// (ρ, z) grid of square cells of side h (midpoint rule).
type Loops struct {
	A, Z, I []float64
}

// NewLoops samples J(ρ, z) on ρ ∈ (0, rmax), z ∈ (zlo, zhi).
func NewLoops(J func(rho, z float64) float64, h, rmax, zlo, zhi float64) Loops {
	var l Loops
	nr, nz := int(rmax/h), int((zhi-zlo)/h)
	for i := 0; i < nr; i++ {
		rho := (float64(i) + 0.5) * h
		for k := 0; k < nz; k++ {
			z := zlo + (float64(k)+0.5)*h
			if c := J(rho, z) * h * h; c != 0 {
				l.A, l.Z, l.I = append(l.A, rho), append(l.Z, z), append(l.I, c)
			}
		}
	}
	return l
}

// Field is the summed field and flux function at (ρ, z).
func (l Loops) Field(rho, z float64) (br, bz, psi float64) {
	for i := range l.A {
		a, b, c := LoopField(l.A[i], l.I[i], rho, z-l.Z[i])
		br, bz, psi = br+a, bz+b, psi+c
	}
	return
}

// SphJ turns J(r, θ) into J(ρ, z).
func SphJ(J func(r, th float64) float64) func(rho, z float64) float64 {
	return func(rho, z float64) float64 {
		r := math.Hypot(rho, z)
		return J(r, math.Atan2(rho, z))
	}
}

// Dipole is the field of a point dipole m [μB] along z at z0 on the axis,
// at (ρ, z), in B0, and its flux function ρA_φ.
func Dipole(m, z0, rho, z float64) (br, bz, psi float64) {
	dz := z - z0
	d2 := rho*rho + dz*dz
	d := math.Sqrt(d2)
	d5 := d2 * d2 * d
	br = 3 * m * rho * dz / d5
	bz = m * (3*dz*dz - d2) / d5
	psi = m * rho * rho / (d2 * d)
	return
}

// parallel runs f(0..n−1) on all CPUs; each call writes only its own slots,
// so results do not depend on scheduling.
func parallel(n int, f func(i int)) {
	var wg sync.WaitGroup
	ch := make(chan int)
	for g := 0; g < runtime.NumCPU(); g++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range ch {
				f(i)
			}
		}()
	}
	for i := 0; i < n; i++ {
		ch <- i
	}
	close(ch)
	wg.Wait()
}
