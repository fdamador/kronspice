// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// FlyingDoughnut is the exact single-cycle toroidal pulse of free space
// found by R. W. Hellwarth and P. Nouchi, "Focused one-cycle
// electromagnetic pulses", Phys. Rev. E 54, 889 (1996), in its transverse-
// electric form: E is azimuthal, B lies in the (ρ, z) plane and forms
// closed rings round the axis. It is generated from the scalar solution of
// the wave equation (Ziolkowski's modified power spectrum pulse, the
// paper's eq. 2)
//
//	f = F0 / [ρ² + (q1 + iτ)(q2 − iσ)],   τ = z − ct,  σ = z + ct,
//
// through the magnetic Hertz vector Π = ẑf: E = −μ0 ∂t curl Π,
// H = curl curl Π, so
//
//	E_φ = μ0 ∂²f/∂t∂ρ,  H_ρ = ∂²f/∂ρ∂z,  H_z = ∂²f/∂z² − (1/c²)∂²f/∂t².
//
// Any f that solves the wave equation gives an exact Maxwell field this
// way; the real part is taken. q1 sets the pulse's length (the effective
// wavelength), q2 ≫ q1 its focal depth; it reaches its focus, smallest
// and brightest, at t = 0 and z = 0. Experimentally generated as light
// pulses by A. Zdagkas et al., Nature Photonics 16, 523 (2022).
type FlyingDoughnut struct {
	Q1, Q2 float64 // [m]
	F0     float64 // amplitude [A·m]
}

// Fields returns E [V/m] and B [T] at (x, y, z) and time t.
func (d FlyingDoughnut) Fields(x, y, z, t float64) (E, B [3]float64) {
	c := SpeedOfLight
	rho := math.Hypot(x, y)
	tau, sig := z-c*t, z+c*t
	a := complex(d.Q1, tau)
	b := complex(d.Q2, -sig)
	D := complex(rho*rho, 0) + a*b
	// derivatives of D (see type comment)
	Dr := complex(2*rho, 0)
	Dz := complex(2*z, d.Q2-d.Q1)
	Dt := complex(-2*c*c*t, -c*(d.Q1+d.Q2))
	const Dzz, Dtt = 2, -2 // Dtt in units of c²
	F := complex(d.F0, 0)
	D2, D3 := D*D, D*D*D
	ftr := F * 2 * Dt * Dr / D3                      // ∂²f/∂t∂ρ (D_tρ = 0)
	frz := F * 2 * Dr * Dz / D3                      // ∂²f/∂ρ∂z (D_ρz = 0)
	fzz := F * (2*Dz*Dz/D3 - Dzz/D2)                 // ∂²f/∂z²
	ftt := F * (2*Dt*Dt/D3 - complex(Dtt*c*c, 0)/D2) // ∂²f/∂t²
	Ephi := real(complex(Mu0, 0) * ftr)
	Hr := real(frz)
	Hz := real(fzz - ftt/complex(c*c, 0))
	cosp, sinp := 1.0, 0.0
	if rho > 0 {
		cosp, sinp = x/rho, y/rho
	}
	E = [3]float64{-Ephi * sinp, Ephi * cosp, 0}
	B = [3]float64{Mu0 * Hr * cosp, Mu0 * Hr * sinp, Mu0 * Hz}
	return
}

// Potentials returns A [T·m] with B = curl A, and C [V] with E = curl C:
// A = μ0 curl Π = μ0·φ̂·2F·ρ/D² and C = −μ0 ∂Π/∂t = μ0·ẑ·F·∂D/∂t / D²
// (real parts). Sampling these and taking the discrete curls gives an
// initial FDTD field with zero discrete divergence.
func (d FlyingDoughnut) Potentials(x, y, z, t float64) (A, C [3]float64) {
	c := SpeedOfLight
	rho := math.Hypot(x, y)
	D := complex(rho*rho, 0) + complex(d.Q1, z-c*t)*complex(d.Q2, -(z+c*t))
	Dt := complex(-2*c*c*t, -c*(d.Q1+d.Q2))
	F := complex(d.F0, 0)
	aphi := Mu0 * real(F*complex(2*rho, 0)/(D*D))
	cosp, sinp := 1.0, 0.0
	if rho > 0 {
		cosp, sinp = x/rho, y/rho
	}
	A = [3]float64{-aphi * sinp, aphi * cosp, 0}
	C = [3]float64{0, 0, Mu0 * real(F*Dt/(D*D))}
	return
}

// DoughnutRun is the outcome of RunFlyingDoughnut.
type DoughnutRun struct {
	RelError               float64 // L2 error of E and cB against the exact pulse, over the interior
	DivB                   float64 // max |div B| against max|B|/dx
	Duration               float64 // simulated time [s]
	EnergyStart, EnergyEnd float64 // field energy in the region [J]
	Absorbed               float64 // energy that left through the absorber [J]
}

// RunFlyingDoughnut starts the FDTD solver from the exact pulse at
// t = −T, T = span·Q1/c, runs it to +T through the focus and compares
// with the exact pulse there. perQ1 is the number of cells per Q1;
// halfWidth the interior half-width in units of Q1 (plus an absorber of
// 10 cells).
func RunFlyingDoughnut(d FlyingDoughnut, perQ1, halfWidth, span float64) (DoughnutRun, error) {
	dx := d.Q1 / perQ1
	abs := 10
	n := 2*int(halfWidth*perQ1) + 2*abs + 1
	half := float64(n-1) / 2 * dx
	f, err := NewFDTD(FDTDConfig{N: [3]int{n, n, n}, D: [3]float64{dx, dx, dx}, Origin: [3]float64{-half, -half, -half}, Absorber: abs, Ledger: true})
	if err != nil {
		return DoughnutRun{}, err
	}
	T := span * d.Q1 / SpeedOfLight / 2
	t0 := -T
	f.SetStartTime(t0)
	f.SetEFromPotential(func(x, y, z float64) [3]float64 { _, C := d.Potentials(x, y, z, t0); return C })
	f.SetBFromPotential(func(x, y, z float64) [3]float64 { A, _ := d.Potentials(x, y, z, t0+f.Dt/2); return A })
	inner := half - float64(abs+2)*dx
	in := func(x, y, z float64) bool { return math.Abs(x) < inner && math.Abs(y) < inner && math.Abs(z) < inner }
	run := DoughnutRun{EnergyStart: f.RegionEnergy(in)}
	f.Run(T)
	run.Duration = f.Time() - t0
	run.EnergyEnd = f.RegionEnergy(in)
	run.Absorbed = f.Ledger().Absorbed
	num, den := 0.0, 0.0
	c := SpeedOfLight
	f.each(func(comp, i, j, k, idx int) {
		if x, y, z := f.ePos(comp, i, j, k); in(x, y, z) {
			E, _ := d.Fields(x, y, z, f.Time())
			num += math.Pow(f.E[comp][idx]-E[comp], 2)
			den += E[comp] * E[comp]
		}
		if x, y, z := f.bPos(comp, i, j, k); in(x, y, z) {
			_, B := d.Fields(x, y, z, f.Time()+f.Dt/2)
			num += c * c * math.Pow(f.B[comp][idx]-B[comp], 2)
			den += c * c * B[comp] * B[comp]
		}
	})
	run.RelError = math.Sqrt(num / den)
	w, s := f.DivB()
	run.DivB = w / s
	return run, nil
}
