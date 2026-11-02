// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command plasmoid_rings reviews J. Nachamkin, "Force-Free Time-Harmonic
// Plasmoids," Phillips Laboratory report PL-TR-92-3044 (1992), and models
// plasmoids in kronspice as toroidal current rings.
//
//  1. His l = 1 force-free (Beltrami) field: ∇×F = αF, field lines on nested
//     tori around a current ring.
//
//  2. His critical frequency and electron density are the cold-plasma
//     dispersion relation; his Appendix B recomputed.
//
//  3. Does the "resonant size" stop radiation? The plasma sphere's l = 1
//     resonances, found exactly.
//
//  4. The zero-frequency limit — the spheromak, a force-free current ring —
//     decaying through Spitzer resistance, as a kronspice circuit.
//
//  5. A collapsing electron ring: electrons circling a solenoid compressed
//     by a kronspice capacitor-bank discharge.
//
//     go run ./examples/plasmoid_rings [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
)

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()
	e := AppendixB()

	fmt.Println("1. The force-free field (Nachamkin eq. 8, l = 1; the spheromak field)")
	curl, div := CurlCheck(e.Alpha)
	axis := MagneticAxis(e.Alpha)
	_, drift := Trace(e.Alpha, [3]float64{0.75 * axis, 0, 0}, 20000, e.R/2000)
	fmt.Printf("   ∇×F = αF to %.0e, ∇·F = 0 to %.0e; a traced field line keeps its flux ψ to %.0e (it stays on its torus)\n", curl, div, drift)
	fmt.Printf("   the tori nest around a current ring at %.3f of the radius (%.1f mm for R = 10 cm): the plasmoid's \"electron ring\"\n", axis/e.R, axis*1e3)

	fmt.Println("\n2. His critical frequency is the cold-plasma cutoff; Appendix B recomputed (R = 10 cm, k = 2686 /m)")
	fmt.Printf("   charge density ρ = (m/e)ε0c²(k² − α²) = %.2f C/m³ (his −32.63), %.3e electrons (his 8.532×10¹⁷)\n", e.RhoE, e.Electrons)
	fmt.Printf("   that density's plasma frequency satisfies ω² = ωp² + α²c² to %.0e: his plasmoid is a standing wave of a cold\n", e.DispersionMismatch)
	fmt.Println("   plasma whose wavenumber is α, and his critical frequency ω0 = αc is the cutoff of that wave.")
	w126 := 7.920e11
	fmt.Printf("   his field for C = 0.01 at %.0f GHz: E0 = C·R·ω²/|e/m| = %.4g V/m (his 3.566×10⁹)\n", w126/2/math.Pi/1e9, e.C*e.R*w126*w126/EoM)
	fmt.Printf("   but then the electrons' oscillation speed is ω·C·R = %.3g m/s = %.2f c — above the speed of light, so his\n", w126*e.C*e.R, w126*e.C*e.R/C0)
	fmt.Println("   non-relativistic electron model does not hold for his own example; real electrons would be relativistic and")
	fmt.Println("   the field ~1,000 times air's breakdown strength.")

	fmt.Println("\n3. Does the resonant size stop radiation? l = 1 resonances of the plasma sphere (exact, outgoing waves only)")
	sp := Sphere{R: e.R, Wp: e.Wp}
	var modes []Mode
	for _, te := range []bool{true, false} {
		m := sp.FindMode(te)
		modes = append(modes, m)
		kin := cmplx.Sqrt(m.Omega*m.Omega-complex(e.Wp*e.Wp, 0)) / C0 * complex(e.R, 0)
		name := map[bool]string{true: "TE (his M part)", false: "TM (his N part)"}[te]
		fmt.Printf("   %-16s %.4f GHz, needs αR = %.4f; Q = %.2g; energy lifetime %.2g s\n", name, real(m.Omega)/2/math.Pi/1e9, real(kin), m.Q, m.Lifetime)
	}
	fmt.Println("   The plasmoid radiates: with nothing feeding it from outside, its energy leaks away in about a microsecond.")
	fmt.Println("   And the two halves of the force-free field need different sizes (αR = 4.49 and 2.74): an equal mix of them,")
	fmt.Println("   which is what makes E parallel to H, is not one free mode of the sphere. It would need driving or surface")
	fmt.Printf("   currents, the weak point Nachamkin names himself. At 1 MJ the initial radiated power is ~%.0e W.\n", 1e6/modes[0].Lifetime)

	fmt.Println("\n4. The zero-frequency limit: the spheromak, a force-free current ring (1 MJ, R = 10 cm), in kronspice")
	s := Spheromak{R: 0.1, W: 1e6}
	fmt.Printf("   field scale B0 = %.0f T, ring current %.2f MA, magnetic pressure %.2g atm (held only by a wall or outside field)\n",
		s.B0(), s.RingCurrent()/1e6, s.B0()*s.B0()/(2*Mu0)/101325)
	var decays []Decay
	for _, t := range []float64{1, 10, 100} {
		d, err := s.RunDecay(t)
		if err != nil {
			log.Fatal(err)
		}
		decays = append(decays, d)
		fmt.Printf("   T = %3.0f eV: ring L = %.3g H, R = %.3g Ω; energy lifetime %.3g s (kronspice %.3g s); heat %.4g J of %.4g J\n",
			t, d.L, d.Rc, d.Tau/2, d.Fitted, d.Heat, d.Stored0)
	}
	fmt.Printf("   To last one second it would need T = %.0f eV (%.0f million K): fusion-plasma temperatures. Ball lightning\n", s.TempFor(1), s.TempFor(1)*11604.5/1e6)
	fmt.Println("   reported to last seconds cannot be a resistive spheromak in air.")

	fmt.Println("\n5. A collapsing electron ring (electron-ring compressor; bank and coil in kronspice)")
	c := DefaultCompressor()
	r, err := c.Run()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   bank %.0f µF at %.0f kV → coil %.0f µH (%.0f turns/m): field %.3f → %.3f T in a quarter period\n", c.C*1e6, c.V0/1e3, c.L*1e6, c.N, r.B0, r.B1)
	fmt.Printf("   ring radius %.1f → %.2f cm (law R ∝ B^-½: %.2f cm); momentum ×%.4f (law ×%.4f)\n",
		r.R0*100, r.R1*100, r.R0*math.Sqrt(r.B0/r.B1)*100, r.P1/r.P0, math.Sqrt(r.B1/r.B0))
	gEnd := r.Gamma[len(r.Gamma)-1]
	fmt.Printf("   electron energy %.1f → %.1f MeV; work of the induced field %.4e J = energy gain %.4e J per electron\n",
		c.E0/1e6, (gEnd-1)*Me*C0*C0/Qe/1e6, r.WorkE, r.DeltaKE)
	fmt.Printf("   canonical angular momentum (cyclic azimuth) held to %.0e; ring current with 10¹³ electrons %.0f → %.0f A\n", r.Invariant, r.RingCurrent0, r.RingCurrent1)
	fmt.Printf("   the bank (%.3g J) drives the coil; the ring takes %.2g J of it; kronspice ledger %.0e\n", r.BankEnergy, c.Electrons*r.DeltaKE, r.Mismatch)
	fmt.Println("   The ring collapses because the field around it rises, and every joule it gains comes from the bank.")
	fmt.Println("\nNo energy appears from the plasmoid or the void inside the ring: radiation, resistance and the bank account for all of it.")

	if *viz == "" {
		return
	}
	write := func(name, head string, f func(w *os.File)) {
		w, err := os.Create(filepath.Join(*viz, name))
		if err != nil {
			log.Fatal(err)
		}
		defer w.Close()
		fmt.Fprintln(w, head)
		f(w)
	}
	write("fieldlines.csv", "line,x,y,z", func(w *os.File) {
		for i, f := range []float64{0.3, 0.6, 0.85} {
			pts, _ := Trace(e.Alpha, [3]float64{axis + f*(e.R-axis)*0.98, 0, 0}, 6000, e.R/400)
			for _, p := range pts {
				fmt.Fprintf(w, "%d,%g,%g,%g\n", i, p[0], p[1], p[2])
			}
		}
	})
	write("decay.csv", "T_eV,t,energy", func(w *os.File) {
		for _, d := range decays {
			for k := range d.Time {
				fmt.Fprintf(w, "%g,%g,%g\n", d.T, d.Time[k], d.Energy[k])
			}
		}
	})
	write("ring.csv", "t,radius,gamma,B", func(w *os.File) {
		for k := range r.Time {
			fmt.Fprintf(w, "%g,%g,%g,%g\n", r.Time[k], r.Radius[k], r.Gamma[k], r.Bfield[k])
		}
	})
}
