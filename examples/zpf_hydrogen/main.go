// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command zpf_hydrogen puts H. E. Puthoff, "Ground state of hydrogen as a
// zero-point-fluctuation-determined state," Phys. Rev. D 35, 3266 (1987),
// on the electron test bench (hypothesis H24).
//
// Stochastic electrodynamics (SED) treats the vacuum as a real, random
// classical field with energy ½ħω per mode. Puthoff shows that an electron
// on a circular orbit absorbs from it exactly what it radiates when
// m·ω·r² = ħ, Bohr's condition, so the orbit need not collapse.
//
//  1. The paper's numbers (eqs. 16–20) with CODATA values.
//
//  2. The averaged balance around the orbit: does it pull back to a0?
//
//  3. The paper's one-dimensional oscillator in a simulated zero-point field
//     (eqs. 4–16): absorbed = radiated, and energy ½ħω0.
//
//  4. The same oscillator as a kronspice circuit: an LC tank driven by its
//     resistor's zero-point noise (Callen–Welton).
//
//  5. Where SED stands: later full simulations, and the standard answer.
//
//     go run ./examples/zpf_hydrogen [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// Settings used by the report (and lighter ones by the tests).
var (
	ReportOsc  = Osc{Gamma: 0.01, Band: 2, N: 1 << 20, Dt: 0.05, Runs: 16, Seed: 1987}
	ReportTank = Tank{R: 0.01, N: 1 << 16, Dt: 0.25, Step: 0.02, Band: 2, Runs: 8, Seed: 1951}
)

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()

	fmt.Println("1. Puthoff's balance with CODATA values (Phys. Rev. D 35, 3266, 1987)")
	r, w, ev := Orbit()
	fmt.Printf("   m·ω·r² = ħ (eq. 20) with the Coulomb orbit gives r = %.6e m (Bohr radius a0 = %.6e m),\n", r, A0)
	fmt.Printf("   orbit frequency %.4e Hz, energy %.4f eV (Bohr and Schrödinger, infinitely heavy proton: %.4f eV)\n", w/(2*math.Pi), ev, -Ry)
	fmt.Printf("   there: absorbed from the zero-point field (eq. 18) %.6e W, radiated (eq. 19) %.6e W\n", Pabs(w), Prad(r, w))
	fmt.Printf("   Without the zero-point field the orbit would spiral in within %.3g s (classical collapse).\n", Collapse())

	fmt.Println("\n2. Is the balance a stable one? Net power P_abs − P_rad on Coulomb orbits (averaged, as in the paper)")
	for _, f := range []float64{0.5, 0.7, 0.9, 1, 1.1, 1.5, 2} {
		rr := f * A0
		ww := Kepler(rr)
		n := Net(rr)
		fmt.Printf("   r = %.1f a0: absorbed/radiated = %.3f, net %+.3e W\n", f, Pabs(ww)/Prad(rr, ww), n)
	}
	fmt.Printf("   Inside a0 the orbit gains energy, outside it loses: the averaged balance pulls back to a0 with time constant\n")
	fmt.Printf("   %.3g s (≈ %.0f orbits). The ratio absorbed/radiated is ħ/(mωr²), which falls as r^½ on Coulomb orbits.\n",
		Relax(), Relax()*w/(2*math.Pi))

	fmt.Println("\n3. The paper's oscillator in a simulated zero-point field (scaled units m = ω0 = ħ = 1; own Runge–Kutta)")
	o := ReportOsc
	res := o.Run()
	gam := E * E / (6 * math.Pi * Eps0 * Me * C * C * C)
	fmt.Printf("   damping Γω0 = %g (the electron's at the hydrogen orbit frequency is %.1e; larger here so a run spans %.0f relaxation times),\n",
		o.Gamma, gam*w, float64(o.N)*o.Dt*o.Gamma)
	fmt.Printf("   field below %g ω0 (declared cut), %d independent random-phase runs\n", o.Band, o.Runs)
	fmt.Printf("   mean energy %.4f ± %.4f ħω0   (expected for this field %.4f; ½ħω0 = %.4f)\n", res.Energy, res.Spread, res.ExpEnergy, res.HalfHbarW)
	fmt.Printf("   absorbed %.6f, radiated %.6f   (expected %.6f; eq. 16: ħΓω0³/2 = %.6f)\n", res.Absorbed, res.Radiated, res.ExpPower, res.Eq16)
	fmt.Printf("   In steady state the oscillator radiates what it absorbs (difference %.2f %%) and holds ½ħω0 of energy.\n",
		100*math.Abs(res.Absorbed-res.Radiated)/res.Radiated)

	fmt.Println("\n4. The same oscillator as a kronspice circuit: series L = 1 H, C = 1 F, R = 0.01 Ω, driven by R's zero-point noise")
	k := ReportTank
	tr, err := k.Run()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   (one-sided S_V = ħωR/π, Callen–Welton 1951; ħ scaled to 1 J·s — the circuit is linear, energies scale with ħ)\n")
	fmt.Printf("   stored energy %.4f ± %.4f ħω0 (expected %.4f); ledger: noise source %.1f J = resistor %.1f J + stored %.2f J (mismatch %.1g)\n",
		tr.Energy, tr.Spread, tr.ExpEnergy, tr.Supplied, tr.Heat, tr.Stored, tr.Mismatch)
	fmt.Println("   The resistor's zero-point noise and its dissipation balance — the fluctuation–dissipation pair behind Puthoff's")
	fmt.Println("   radiation balance (his footnote 13). Net, the zero-point field supplies nothing: what goes in comes back out.")

	fmt.Println("\n5. Where this stands (H24)")
	fmt.Println("   - Puthoff's result is at the level of Bohr theory, as he says. Quantum mechanics' ground state has zero orbital")
	fmt.Println("     angular momentum and a spread of radii (most probable a0, mean 1.5 a0); a single orbit with L = ħ is not it.")
	fmt.Println("     Extending SED to the Schrödinger ground state was, in his words, \"yet to be achieved\" (footnote 17).")
	fmt.Println("   - Full SED simulations since: Cole & Zou (Phys. Lett. A, 2003) report the Schrödinger ground-state density;")
	fmt.Println("     Nieuwenhuizen & Liska (Phys. Scr. T165, 014006, 2015), with longer runs, find the atom self-ionizes.")
	fmt.Println("   - Energy: the balance is zero net exchange. The zero-point field here is not an energy source.")
	fmt.Printf("   - The energy and size come out right because eq. 20 is Bohr's condition: %.4f eV, the same as Bohr and\n", ev)
	fmt.Println("     Schrödinger for an infinitely heavy proton (measured ionization energy 13.598 eV includes the proton's recoil).")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "balance.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(f, "r_over_a0,absorbed,radiated")
	for x := 0.3; x <= 3.0001; x += 0.02 {
		rr := x * A0
		ww := Kepler(rr)
		fmt.Fprintf(f, "%g,%g,%g\n", x, Pabs(ww), Prad(rr, ww))
	}
	f.Close()
	g, err := os.Create(filepath.Join(*viz, "oscillator.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(g, "t,energy")
	for i, e := range res.EnergyTrace {
		if i%5 == 0 {
			fmt.Fprintf(g, "%g,%g\n", float64(i*20)*o.Dt+10/o.Gamma, e)
		}
	}
	g.Close()
}
