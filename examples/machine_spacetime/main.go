// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command machine_spacetime puts one rotor — rotation (special relativity)
// — into three gravitational fields — curvature (general relativity) — and
// then asks the reverse question: how much does the machine itself curve
// and twist spacetime?
//
// The combination is the equivalence principle: the hub sits in the
// field (Schwarzschild clock rate, moving with its site), and in the hub's
// local inertial frame the rim is ordinary special relativity. Its spin
// axis precesses by Thomas + geodetic + Lense–Thirring (Schiff 1960), and a
// fibre loop on it sees the Sagnac delay relative to the local gyroscopes.
// The model is validated in spacetime/machine_test.go against Gravity
// Probe B (Everitt et al., PRL 106, 221101, 2011).
//
//	go run ./examples/machine_spacetime
package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/spacetime"
)

// Solar GM (IAU 2015 Resolution B3 nominal value).
const gmSun = 1.3271244e20

func main() {
	// The rotor: a 0.5 m T1000G carbon-fibre disk, 2 cm thick, at half its
	// burst speed, carrying a 1 m² fibre loop.
	mat := spacetime.T1000G
	rad := 0.5
	omega := 0.5 * spacetime.BurstSpeed(mat) / rad
	mass := mat.Density * math.Pi * rad * rad * 0.02
	m := spacetime.Machine{Name: "0.5 m T1000G rotor at half burst speed", Radius: rad, Omega: omega, Material: mat, Loop: 1}

	// Sites.
	earthLab := spacetime.EarthSurface(45 * math.Pi / 180)

	a := spacetime.REarthEq + 642e3 // Gravity Probe B's orbit, polar
	v := math.Sqrt(spacetime.GMEarth / a)
	orbit := spacetime.Site{Name: "642 km polar orbit (GP-B), inertially pointed",
		R: spacetime.Vec{a, 0, 0}, V: spacetime.Vec{0, 0, v}, SpinAx: spacetime.Vec{0, 1, 0}}

	// An illustrative neutron star (assumed values, not a measured
	// object): 1.4 solar masses, 12 km radius, 10 ms spin, I = 0.4MR².
	nsM := 1.4 * gmSun / spacetime.G
	nsR := 12e3
	ns := spacetime.Body{Name: "neutron star (illustrative)", GM: 1.4 * gmSun,
		J: spacetime.Vec{0, 0, 0.4 * nsM * nsR * nsR * 2 * math.Pi / 0.01}, R: nsR}
	hover := 1000e3
	g := ns.GM / (hover * hover)
	nsSite := spacetime.Site{Name: "hovering 1000 km above its centre (rocket holds it up)",
		R: spacetime.Vec{hover, 0, 0}, ANG: spacetime.Vec{g, 0, 0}, SpinAx: spacetime.Vec{1, 0, 0}}

	for _, c := range []struct {
		b spacetime.Body
		s spacetime.Site
	}{{spacetime.Earth, earthLab}, {spacetime.Earth, orbit}, {ns, nsSite}} {
		rep, err := spacetime.AnalyzeMachine(c.b, c.s, m)
		if err != nil {
			fmt.Println(c.s.Name, "error:", err)
			continue
		}
		fmt.Printf("\n== %s — %s ==\n", c.b.Name, rep.Site)
		fmt.Printf("  hub clock rate         %.15f  (gravity + site motion)\n", rep.HubRate)
		fmt.Printf("  rim clock rate         %.15f\n", rep.RimRate)
		fmt.Printf("  rim loses vs hub       %.3g s/day   (rim speed %.0f m/s)\n", rep.RimLagPerDay, rep.RimSpeed)
		fmt.Printf("  GM/rc²                 %.3g        v²/c² at rim %.3g\n", rep.PhiOverC2, rep.VOverC2)
		fmt.Printf("  burst margin           %.2f        tidal stress %.3g Pa (margin %.3g)\n", rep.BurstMargin, rep.TidalStress, rep.TidalMargin)
		fmt.Printf("  spin-axis precession   (here) Thomas %.3g, geodetic %.3g, frame dragging %.3g mas/yr\n",
			spacetime.MasPerYear(rep.Thomas.Norm()), spacetime.MasPerYear(rep.Geodetic.Norm()), spacetime.MasPerYear(rep.FrameDrag.Norm()))
		fmt.Printf("  Sagnac (1 m² loop)     %.4g s per pass; local-frame correction %.3g s\n", rep.SagnacPlatform, rep.SagnacFrameCorrection+0)
		if c.s.Name == orbit.Name {
			geo, fd := spacetime.OrbitAveragedPrecession(c.b, a, spacetime.Vec{1, 0, 0}, spacetime.Vec{0, 0, 1}, 360)
			fmt.Printf("  averaged over the orbit geodetic %.1f, frame dragging %.2f mas/yr (frame dragging is 39.1 seen against GP-B's guide star; measured 6601.8 ± 18.3 and 37.2 ± 7.2)\n",
				spacetime.MasPerYear(geo.Norm()), spacetime.MasPerYear(fd.Norm()))
		}
		for _, w := range rep.Warnings {
			fmt.Println("  warning:", w)
		}
	}

	// The reverse question: the machine's own gravity, 1 m from its centre
	// on its axis, beside Earth's at the surface.
	sf := spacetime.MachineSelfField(m, mass, 1)
	earthPhi := spacetime.GMEarth / (spacetime.REarthEq * spacetime.C * spacetime.C)
	earthCurv := spacetime.GMEarth / (spacetime.C * spacetime.C * math.Pow(spacetime.REarthEq, 3))
	_, _, fdE := spacetime.Precession(spacetime.Earth, earthLab.R, spacetime.Vec{}, spacetime.Vec{})
	fmt.Printf("\n== The machine's own spacetime (%.1f kg, J = %.3g kg·m²/s), 1 m away on its axis ==\n", sf.Mass, sf.J)
	fmt.Printf("  clock shift GM/dc²     %.3g   (Earth's at its surface %.3g: ratio %.2g)\n", sf.PhiOverC2, earthPhi, sf.PhiOverC2/earthPhi)
	fmt.Printf("  tidal curvature        %.3g /m² (Earth's %.3g: ratio %.2g)\n", sf.Curvature, earthCurv, sf.Curvature/earthCurv)
	fmt.Printf("  frame dragging         %.3g rad/s = %.3g mas/yr (Earth's at 45° %.3g mas/yr: ratio %.2g)\n",
		sf.FrameDragAxis, spacetime.MasPerYear(sf.FrameDragAxis), spacetime.MasPerYear(fdE.Norm()), sf.FrameDragAxis/fdE.Norm())
	fmt.Println("  Its mass curves spacetime measurably at close range (torsion balances weigh that tidal pull);")
	fmt.Println("  its spin drags frames about 10^-10 of Earth's own frame dragging, far below any instrument.")
}
