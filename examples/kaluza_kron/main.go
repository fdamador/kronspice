// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kaluza_kron works through Th. Kaluza, "On the Unification
// Problem in Physics" (Sitzungsber. Preuss. Akad. Wiss. 1921, 966–972;
// V. T. Toth's revised translation) and its link to Kron's machines.
//
//  1. Kaluza's three-index symbols (eqs. 1–3) checked from his 5-D metric.
//
//  2. The exact 5-D geodesic of a charge in a magnetic field, against the
//     Lorentz force: charge is momentum along the fifth dimension.
//
//  3. Kaluza's electron difficulty (p. 971) in numbers, and Klein's 1926
//     circle that quantizes charge.
//
//  4. The same construction inside a Kron machine: winding charges are
//     cyclic coordinates (Kaluza's "cylinder condition", in Hoffmann's
//     1949 words), their flux linkages are conserved, and eliminating them
//     leaves the rotor's force — solved three ways, kronspice included.
//
//     go run ./examples/kaluza_kron [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()

	fmt.Println("1. Kaluza's equations 1–3")
	a := KaluzaAlpha()
	fmt.Printf("   α = √(κ/2), κ = 8πG/c²: %.4e (Kaluza: 3.06×10⁻¹⁴)\n", a)
	fmt.Printf("   Three-index symbols of his 5-D metric (eq. 1) against his closed forms (eq. 3), sample fields:\n")
	fmt.Printf("   largest difference %.1e with his α, %.1e with α = 0.7 — eq. 3 holds as printed.\n", Eq3Check(a), Eq3Check(0.7))
	fmt.Printf("   His electron e/m = 1.77×10⁷ is e/m in electromagnetic units per gram (today %.4e).\n", Ecgs/Ccgs/MeCgs)

	fmt.Println("\n2. A charge in a magnetic field as free motion in five dimensions (exact geodesic; c = G = 1, Gaussian)")
	fmt.Println("   Metric ĝ = (η + k²AA, kA; kA, 1), k = 2√G, uniform B = 0.05 along z, speed 0.3 c")
	fmt.Printf("   %8s %10s %12s %12s %12s %12s %10s\n", "p_w", "q/m", "radius", "Lorentz", "ω", "Lorentz", "p_w drift")
	var orbits []Orbit
	for _, pw := range []float64{0.1, 0.3, 1, 3, 10} {
		o := RunOrbit(0.05, 2, 0.3, pw)
		orbits = append(orbits, o)
		fmt.Printf("   %8.2f %10.4f %12.5f %12.5f %12.6f %12.6f %10.1e\n", pw, o.QoverM, o.Radius, o.RadiusLF, o.Omega, o.OmegaLF, o.PwDrift)
	}
	fmt.Println("   Nothing depends on w, so the momentum p_w along it is conserved (to round-off), and the 4-D path is a")
	fmt.Println("   Lorentz-force orbit with q/m = k·p_w/√(1 + p_w²): charge is motion along the fifth dimension.")
	fmt.Println("   Kaluza's eq. 8a, ρ0/μ0 = 2α·u⁰, is the small-p_w part of this curve.")

	fmt.Println("\n3. Kaluza's difficulty with the electron (p. 971), and Klein's answer")
	fmt.Printf("   q/m in these units can never exceed k = 2 (√G per unit mass), however fast the motion along w.\n")
	fmt.Printf("   Electron: e/(m√G) = %.3e;  proton: %.3e. Neither fits on a classical 5-D geodesic.\n", ChargeToMass(Ecgs, MeCgs), ChargeToMass(Ecgs, MpCgs))
	fmt.Printf("   Klein (1926): x⁰ is a circle; quantum momentum p_w = nħ/R makes charge come in steps n·e when\n")
	fmt.Printf("   R = 2√G·ħ/(c·e) = 2 l_P/√α = %.3e m. The first charged mode then weighs ħ/(Rc) = %.2e GeV —\n", KleinRadius(), KKMassGeV())
	fmt.Println("   near the Planck mass, not an electron. (In my first reply I gave 1.3×10⁻³³ m: that paired the Heaviside")
	fmt.Println("   coupling √(16πG) with the Gaussian charge, a factor √(4π) too large.)")

	fmt.Println("\n4. The same construction in a Kron machine (stator coil, rotor coil, mutual M·cos θ; windings superconducting)")
	m := DefaultMachine()
	g := m.Geodesic(2, 1e-4)
	r := m.Reduced(2, 1e-4)
	k, mis, err := m.KronSpice(2, 1e-4, 1e-7)
	if err != nil {
		log.Fatal(err)
	}
	p1, p2 := m.Fluxes()
	fmt.Printf("   L1 = L2 = %.1f H, M = %.1f H, J = %.2f kg·m²; start i1 = %.1f A, i2 = %.1f A, θ = %.1f rad, at rest\n", m.L1, m.M, m.J, m.I1, m.I2, m.Theta0)
	fmt.Printf("   conserved flux linkages ψ1 = %.4f Wb, ψ2 = %.4f Wb (Kron geodesic: drift %.0e; energy drift %.0e)\n", p1, p2, g.FluxDrift, g.EnergyDrift)
	fmt.Printf("   %6s %16s %16s %16s\n", "t (s)", "Kron geodesic", "rotor alone", "kronspice")
	dg, dk := 0.0, 0.0
	for _, t := range []float64{0.1, 0.5, 1, 1.5, 2} {
		a, b, c := at(g, t), at(r, t), at(k, t)
		fmt.Printf("   %6.1f %13.6f rad/s %13.6f rad/s %13.6f rad/s\n", t, a, b, c)
	}
	for i := range g.Time {
		dg = math.Max(dg, math.Abs(g.Omega[i]-at(r, g.Time[i])))
		dk = math.Max(dk, math.Abs(g.Omega[i]-at(k, g.Time[i])))
	}
	fmt.Printf("   largest speed difference: rotor alone %.1e rad/s, kronspice %.1e rad/s; kronspice ledger mismatch %.0e\n", dg, dk, mis)
	fmt.Println("   The charges never appear in the inductances (cylinder condition), so their momenta — flux linkages — are")
	fmt.Println("   conserved, like p_w. Eliminating them leaves the rotor in a potential set by the fluxes: the machine's")
	fmt.Println("   version of the force the fifth dimension leaves behind. Kron's machine is a working Kaluza model.")
	fmt.Println("\nNo energy appears from the extra dimension: the 5-D geodesic conserves its norm and the machine its energy.")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "orbits.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(f, "pw,x,y")
	for _, o := range orbits {
		for i := range o.X {
			fmt.Fprintf(f, "%g,%g,%g\n", o.Pw, o.X[i], o.Y[i])
		}
	}
	f.Close()
	h, err := os.Create(filepath.Join(*viz, "qm.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(h, "pw,exact,kaluza_linear")
	for x := 0.0; x <= 5.0001; x += 0.05 {
		fmt.Fprintf(h, "%g,%g,%g\n", x, 2*x/math.Sqrt(1+x*x), 2*x)
	}
	h.Close()
	w, err := os.Create(filepath.Join(*viz, "machine.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(w, "t,geodesic,reduced,kronspice")
	for i := range g.Time {
		fmt.Fprintf(w, "%g,%g,%g,%g\n", g.Time[i], g.Omega[i], at(r, g.Time[i]), at(k, g.Time[i]))
	}
	w.Close()
}

// at returns the trace's speed at the sample nearest t (traces are
// sampled on uniform grids, so a binary search finds it).
func at(tr MTrace, t float64) float64 {
	lo, hi := 0, len(tr.Time)-1
	for hi-lo > 1 {
		mid := (lo + hi) / 2
		if tr.Time[mid] < t {
			lo = mid
		} else {
			hi = mid
		}
	}
	if math.Abs(tr.Time[lo]-t) < math.Abs(tr.Time[hi]-t) {
		return tr.Omega[lo]
	}
	return tr.Omega[hi]
}
