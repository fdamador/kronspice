// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command shoulders_ev tests K. R. Shoulders' "EV" (electrum validum), the
// micrometre bead of about 10¹¹ electrons described in EV: A Tale of
// Discovery (1987), against standard electrostatics, and runs his EV
// source and an EV strike as kronspice circuits with energy ledgers.
//
//  1. The bead's own electrostatics: energy, field, pressure, and how fast
//     a bare cluster flies apart.
//
//  2. Candidate containers: the Casimir (zero-point) pressure, Boyer's
//     sphere, the magnetic pinch at c/10, the field a magnetic trap would
//     need.
//
//  3. Gauss's law against his kilovolt probe measurements.
//
//  4. Charge and energy budgets: his smallest source (1 pF, < 1 µJ) and his
//     larger one (40 pF, 10 kV) against what he reports.
//
//  5. kronspice: the source discharge and an EV strike on a 1 pF anode
//     into a 50 Ω scope.
//
//     go run ./examples/shoulders_ev [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// Sources are Shoulders' two EV sources (2-7): capacitance and voltage.
var Sources = []struct {
	Name string
	C, V float64
}{
	{"smallest (1 pF, 1 µJ)", 1e-12, 1414},
	{"earlier (40 pF, 10 kV)", 40e-12, 10e3},
}

// Claims are the charges Shoulders reports (4-9, 4-63, 8-49).
var Claims = []struct {
	Name string
	N    float64
}{
	{"residual charge of a 3 µm EV on an electrode", 2e10},
	{"one EV bead", 1e11},
	{"emitted by a 10 µm EV over 1 mm of track", 3.5e14},
}

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()

	e := Book()
	fmt.Println("1. The EV as Shoulders describes it: 10¹¹ electrons in a 1 µm sphere moving at c/10")
	fmt.Printf("   charge %.2e C;  density %.2e electrons/cm³ (his 3.5×10¹⁴ in 10 µm gives 6.7×10²³: his 6.6×10²³)\n", e.Charge(), e.Density()/1e6)
	fmt.Printf("   energy to assemble it %.2f J;  surface field %.1e V/m;  outward pressure %.1e Pa (%.0e atm)\n", e.SelfEnergy(), e.SurfaceField(), e.Pressure(), e.Pressure()/101325)
	fmt.Printf("   kinetic energy of all its electrons at c/10: %.1e J — the bead's field holds %.1e times more\n", e.Kinetic(), e.SelfEnergy()/e.Kinetic())
	t2, t10, kev := e.Explode()
	fmt.Printf("   bare, it flies apart: radius ×2 in %.1e s, ×10 in %.1e s; outer electrons leave with %.0f MeV\n", t2, t10, kev/1e6)
	fmt.Printf("   Shoulders' lifetime (1 mm at c/10) is %.1e s: %.0f times longer than the doubling time\n", Lifetime(), Lifetime()/t2)

	fmt.Println("\n2. What could hold it together? (outward pressure to beat:", fmt.Sprintf("%.1e Pa)", e.Pressure()))
	casimir := CasimirPlates(e.D)
	boyer := BoyerSphere(e.R())
	fmt.Printf("   Casimir pull between ideal plates 1 µm apart:  %.1e Pa  (%.0e too weak)\n", casimir, e.Pressure()/casimir)
	fmt.Printf("   zero-point stress on a 1 µm conducting sphere: %.1e Pa, and it pushes OUTWARD (Boyer 1968)\n", boyer)
	fmt.Printf("   magnetic attraction of electrons moving together at c/10: %.0f %% of their repulsion\n", 100*e.MagneticFraction())
	bf := BrillouinField(e.Density())
	fmt.Printf("   a magnetic trap at this density needs %.1e T (Brillouin limit); the indoor record is about 1,200 T\n", bf)
	fmt.Println("   Positive ions would neutralize it, but Shoulders measured fewer than 1 ion per 100,000 electrons.")

	fmt.Println("\n3. Gauss's law and his probe measurements")
	fmt.Printf("   leaving the bare bead an electron gains %.0f MeV; Shoulders sees \"only kilovolts at best\" (9-2).\n", e.EscapeEnergy()/1e6)
	fmt.Println("   His answer: a container hides the core's field. But outside ANY closed surface the potential of a net")
	fmt.Println("   charge q is k·q/r, whatever is inside. For electrons to gain only V volts the container must be at least")
	for _, v := range []float64{1e3, 10e3} {
		fmt.Printf("     %5.0f V:  r ≥ k·q/V = %.1f cm  (q = 10¹¹ e)\n", v, 100*GaussRadius(e.Charge(), v))
	}
	fmt.Printf("   — or the net charge is at most %.0e electrons for a 1 µm bead seen at 1 kV.\n", 1e3*e.R()/Ke/Qe)

	fmt.Println("\n4. Charge and energy budgets")
	for _, s := range Sources {
		q := s.C * s.V
		fmt.Printf("   source %-24s holds %.2e electrons and %.1e J\n", s.Name, q/Qe, 0.5*s.C*s.V*s.V)
	}
	small := Sources[0].C * Sources[0].V / Qe
	big := Sources[1].C * Sources[1].V / Qe
	for _, c := range Claims {
		fmt.Printf("   %-47s %.1e electrons = %6.0f × small source, %8.2f × larger source\n", c.Name, c.N, c.N/small, c.N/big)
	}
	fmt.Printf("   a 1 µm bead of 10¹¹ electrons needs %.2f J to assemble: %.0e × the small source's energy, %.0f × the larger one's\n",
		e.SelfEnergy(), e.SelfEnergy()/(0.5*1e-12*1414*1414), e.SelfEnergy()/(0.5*40e-12*1e8))

	fmt.Println("\n5. kronspice: 1 pF anode into a 50 Ω scope (declared bench)")
	b := DefaultBench()
	src, err := b.Source(1e-12, 1414, 20, 5e-9, 3e-9)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("   a. his smallest source discharging across the gap (20 Ω arc, 5 nH)")
	fmt.Printf("      charge through the scope %.2e electrons;  energy: stored %.3f µJ → scope %.3f µJ, arc heat %.3f µJ;  ledger %.0e\n",
		src.Charge/Qe, 1e6*src.Stored0, 1e6*src.ToScope, 1e6*(src.Heat-src.ToScope), src.Mismatch)
	fmt.Println("   b. an EV strike arriving in 20 ps (an ideal current source brings the charge)")
	var strikes []Shot
	for _, n := range []float64{2e10, 1e11} {
		s, err := b.Strike(n, 20e-12, 1e-9)
		if err != nil {
			log.Fatal(err)
		}
		strikes = append(strikes, s)
		fmt.Printf("      %.0e electrons: peak %6.0f V at the scope; the arriving charge must bring %7.2f µJ (%.0f × the source); ledger %.0e\n",
			n, s.PeakV, 1e6*s.EnergyIn, s.EnergyIn/src.Stored0, s.Mismatch)
	}
	fmt.Println("   Charge is conserved exactly, so everything at the anode is limited by C·V of the source: 8.8×10⁹ electrons here.")
	fmt.Println("   A Faraday cup and this ledger on one shot would settle it.")

	fmt.Println("\nAt the reported size and charge the EV would hold 1,400 to 3 million times its source's energy and fly")
	fmt.Println("apart in femtoseconds; no known container — zero-point, magnetic or ionic — comes close. The bench")
	fmt.Println("effects (pits, rings, X-rays, fast pulses) are real observations; their charge and size need measuring.")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "scope.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(f, "case,t_ns,v")
	put := func(name string, s Shot) {
		for i := range s.Time {
			if i%4 == 0 {
				fmt.Fprintf(f, "%q,%g,%g\n", name, 1e9*s.Time[i], s.Volt[i])
			}
		}
	}
	put("source 1 pF, 1.4 kV", src)
	put("strike 2e10 electrons", strikes[0])
	put("strike 1e11 electrons", strikes[1])
	f.Close()

	g, err := os.Create(filepath.Join(*viz, "pressure.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(g, "name,pa")
	fmt.Fprintf(g, "%q,%g\n", "EV push outward", e.Pressure())
	fmt.Fprintf(g, "%q,%g\n", "magnetic pinch at c/10", e.Pressure()*e.MagneticFraction())
	fmt.Fprintf(g, "%q,%g\n", "1200 T magnet (record)", 1200*1200/(2*4e-7*math.Pi))
	fmt.Fprintf(g, "%q,%g\n", "Casimir, plates 1 µm", casimir)
	g.Close()

	h, err := os.Create(filepath.Join(*viz, "budget.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(h, "name,electrons")
	fmt.Fprintf(h, "%q,%g\n", "small source holds", small)
	fmt.Fprintf(h, "%q,%g\n", "larger source holds", big)
	for _, c := range Claims {
		fmt.Fprintf(h, "%q,%g\n", c.Name, c.N)
	}
	h.Close()
}
