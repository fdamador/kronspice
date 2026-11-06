// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/kron"
)

// H30 Part 3: the Euler–Heisenberg vacuum as the medium that would supply
// the inward wave.
//
// Parts 1 and 2 showed that a knot kept pulsing needs an inward wave of
// one knot's energy per period, from outside. The one nonlinear vacuum
// with measured support is QED's (Euler–Heisenberg; light-by-light
// scattering seen by ATLAS 2017 and CMS 2019). In it waves interact only
// where they cross, through the invariants F = E² − c²B² and G = E·B; a
// single plane wave has F = G = 0 and never acts on itself. So the
// question is how much of an explosion such a vacuum turns back inward.
// Both measurements run on Kron's plane-wave network with the
// Euler–Heisenberg vacuum as its nonlinear primitive (kron.EHLine), each
// beside Maxwell's linear vacuum, at 0.09 E_s, just inside the weak-field
// expansion (em.EHWeakLimit):
//
//  A. The explosion. A compressed field at rest (E only, the moment of
//     tightest implosion) splits into two halves running apart. While they
//     overlap they cross, the one place this vacuum can act. Measured:
//     the energy of the nonlinear part of the field (EH run less linear
//     run) found between the halves after they separate, i.e. moving back
//     to the centre, as a fraction of the explosion's energy.
//  B. Reflection. A weak probe crosses a strong, sharp-edged counter-
//     propagating field (0.08 E_s; the only thing in this vacuum that can
//     turn a wave round). Measured: the fraction of the probe's energy
//     travelling back with the pump. Compared with the Fresnel bound for
//     an index step of δn = 16a, (δn/2)².
//
// In three dimensions the knot's expanding shell becomes locally a plane
// wave within a few L0, so the vacuum stops acting on it; whatever it can
// turn round is bounded by the nonlinear part of the field energy,
// 6a ≈ 1.3×10⁻⁶ at 0.09 E_s (em.EHEnergyDensity), set beside the one
// knot per period the pulsing knot needs.

// VacuumResult holds Part 3.
type VacuumResult struct {
	A          float64 // a = (α/90π)(E/E_s)² at the peak field
	Returned   float64 // A: nonlinear energy between the separated halves / explosion energy
	ReturnedLn float64 // the same with the linear vacuum (lattice only)
	Altered    float64 // A: energy of the whole nonlinear part / explosion energy
	PumpA      float64 // B: a at the pump's field (0.08 E_s: its sharp edges overshoot on the lattice)
	Reflected  float64 // B: probe energy turned round / probe energy
	Fresnel    float64 // (δn/2)² with δn = 16a at the pump
	Nonlinear  float64 // the nonlinear share of the field energy at the peak, ≈ 6a
}

// cellEnergy is the linear field energy of E (links) and B (faces) in cell
// k per unit area and length, the measure for small difference fields.
func cellEnergy(ex, ey, bx, by float64) float64 {
	return 0.5*em.Eps0*(ex*ex+ey*ey) + 0.5*(bx*bx+by*by)/em.Mu0
}

// RunVacuum runs Part 3.
func RunVacuum() (VacuumResult, error) {
	const (
		dz      = 50e-9
		courant = 0.5
	)
	c := em.SpeedOfLight
	peak := 0.09 * em.SchwingerField
	var r VacuumResult
	r.A = em.EHBirefringenceParameter(peak / c)
	r.Nonlinear = (em.EHEnergyDensity([3]float64{peak}, [3]float64{}) - 0.5*em.Eps0*peak*peak) / (0.5 * em.Eps0 * peak * peak)

	// A. The explosion: E only, a Gaussian of 20 cells at the centre.
	{
		const n, mid, sigma, steps = 2000, 1000.0, 20.0, 2000 // halves 500 cells out
		run := func(linear bool) (*kron.EHLine, error) {
			l, err := kron.NewEHLine(n, dz, courant)
			if err != nil {
				return nil, err
			}
			l.Linear = linear
			if err := l.SetFields(func(z float64) [3]float64 {
				u := (z/dz - mid) / sigma
				return [3]float64{peak * math.Exp(-u*u/2)}
			}, func(float64) [3]float64 { return [3]float64{} }); err != nil {
				return nil, err
			}
			for s := 0; s < steps; s++ {
				if err := l.Step(); err != nil {
					return nil, err
				}
			}
			return l, nil
		}
		nl, err := run(false)
		if err != nil {
			return r, err
		}
		ln, err := run(true)
		if err != nil {
			return r, err
		}
		total := 0.5 * em.Eps0 * peak * peak * sigma * math.Sqrt(math.Pi) // ∫½ε₀E² dz in cells
		// Between the halves: |z − mid| < travelled − 6σ.
		inner := float64(steps)*courant - 6*sigma
		var back, backLn, all float64
		for k := 0; k < n; k++ {
			d := cellEnergy(nl.Ex[k]-ln.Ex[k], nl.Ey[k]-ln.Ey[k], nl.Bx[k]-ln.Bx[k], nl.By[k]-ln.By[k])
			all += d
			if math.Abs(float64(k)-mid) < inner {
				back += d
				backLn += cellEnergy(ln.Ex[k], ln.Ey[k], ln.Bx[k], ln.By[k])
			}
		}
		r.Returned, r.ReturnedLn, r.Altered = back/total, backLn/total, all/total
	}

	// B. Reflection of a probe by a counter-propagating pump with sharp
	// edges (10 cells, half the probe wavelength; sharper ones ring past the
	// weak-field limit on the lattice), the strongest reflector the vacuum offers.
	{
		const (
			n, z0, sigma, lambda = 4000, 300.0, 40.0, 20.0
			zs, lp, ramp         = 600.0, 800.0, 10.0
		)
		probeAmp := 1e-4 * em.SchwingerField
		pumpE := 0.08 * em.SchwingerField
		r.PumpA = em.EHBirefringenceParameter(pumpE / c)
		r.Fresnel = math.Pow(16*r.PumpA/2, 2)
		steps := int((zs + lp + 2*ramp - z0 + 8*sigma) / 2 / courant)
		half := courant / 2
		pump := func(x float64) float64 {
			switch {
			case x < 0 || x > lp+2*ramp:
				return 0
			case x < ramp:
				s := math.Sin(math.Pi * x / (2 * ramp))
				return s * s
			case x > lp+ramp:
				s := math.Sin(math.Pi * (lp + 2*ramp - x) / (2 * ramp))
				return s * s
			}
			return 1
		}
		probe := func(z float64) float64 {
			u := (z - z0) / sigma
			return math.Exp(-u*u/2) * math.Cos(2*math.Pi*(z-z0)/lambda)
		}
		run := func(withPump, withProbe bool) (*kron.EHLine, error) {
			l, err := kron.NewEHLine(n, dz, courant)
			if err != nil {
				return nil, err
			}
			pp, pr := 0.0, 0.0
			if withPump {
				pp = 1
			}
			if withProbe {
				pr = 1
			}
			if err := l.SetFields(func(z float64) [3]float64 {
				zc := z / dz
				return [3]float64{pp*pumpE*pump(zc-zs) + pr*probeAmp*probe(zc)}
			}, func(z float64) [3]float64 {
				zc := z / dz
				// pump along −z (B_y = −E_x/c), probe along +z
				return [3]float64{0, -pp*pumpE*pump(zc+half-zs)/c + pr*probeAmp*probe(zc-half)/c}
			}); err != nil {
				return nil, err
			}
			for s := 0; s < steps; s++ {
				if err := l.Step(); err != nil {
					return nil, err
				}
			}
			return l, nil
		}
		both, err := run(true, true)
		if err != nil {
			return r, err
		}
		pumpOnly, err := run(true, false)
		if err != nil {
			return r, err
		}
		probeOnly, err := run(false, true)
		if err != nil {
			return r, err
		}
		probeEnergy := 0.5 * em.Eps0 * probeAmp * probeAmp * sigma * math.Sqrt(math.Pi) // ½ of ⟨cos²⟩·2 (E and B halves)
		zEnd := z0 + float64(steps)*courant
		var turned float64
		for k := 0; k < n; k++ {
			if math.Abs(float64(k)-zEnd) < 6*sigma {
				continue // the transmitted probe and its phase shift
			}
			turned += cellEnergy(both.Ex[k]-pumpOnly.Ex[k]-probeOnly.Ex[k], 0,
				0, both.By[k]-pumpOnly.By[k]-probeOnly.By[k])
		}
		r.Reflected = turned / probeEnergy
	}
	return r, nil
}

func reportVacuum() {
	fmt.Println("\nPart 3: the Euler–Heisenberg vacuum (QED, measured) as the source of the inward wave")
	fmt.Println("Kron's plane-wave network with the Euler–Heisenberg vacuum as its nonlinear primitive (kron.EHLine), at 0.09 E_s")
	v, err := RunVacuum()
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("  a = (α/90π)(E/E_s)² = %.3g; nonlinear share of the field energy at the peak %.3g (≈ 6a)\n", v.A, v.Nonlinear)
	fmt.Printf("  A. explosion: nonlinear energy back between the separated halves %.2g of the explosion (linear vacuum, lattice only: %.2g);\n", v.Returned, v.ReturnedLn)
	fmt.Printf("     the whole nonlinear part of the field %.2g (phase shifts the halves gave each other while crossing)\n", v.Altered)
	fmt.Printf("  B. a probe crossing a counter-propagating field of 0.08 E_s (a = %.3g): %.2g of it turned round (Fresnel bound (8a)² = %.2g)\n", v.PumpA, v.Reflected, v.Fresnel)
	fmt.Println("  Needed for a knot rebuilt every period: 1 (one knot's energy inward per period).")
	fmt.Println("  So the measured QED vacuum returns none of the explosion that can be told from zero, and could turn")
	fmt.Println("  round at most ~6a ≈ 1e-6 of it even in principle: it is not the source of the inward wave.")
}
