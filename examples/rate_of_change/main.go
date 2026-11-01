// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command rate_of_change shows the two places where higher rates of change
// enter electromagnetism in kronspice:
//
//  1. In the material (standard physics, em/dispersion.go): water's
//     molecules need ≈ 8 ps to turn with a field, so water answers a slow
//     field with ε ≈ 78.6 and a fast one with ε ≈ 5.4 — electrically it
//     "stiffens" when hit fast, like water slapped hard. The polarization
//     equation is solved in time (a field ramped at different speeds) and in
//     frequency, and the same water is built as a plate-cell circuit (C∞ ∥
//     R–C) and solved by kronspice AC; all three agree.
//
//  2. In the field equations (hypothetical, register H19, em/podolsky.go):
//     Bopp–Podolsky electrodynamics, (1 + a²□)□A = μ0J, fourth order in
//     time and space, with one new length a. Its predictions are printed
//     beside standard Maxwell (a = 0).
//
//     go run ./examples/rate_of_change
package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/netlist"
)

// Cell is the plate cell used for the circuit twin: 1 cm² plates, 1 mm apart.
const (
	CellArea = 1e-4
	CellGap  = 1e-3
)

// CircuitTwin solves the water-filled cell in kronspice AC at the given
// frequencies and returns the largest relative difference between the
// source current and jω(ε0A/d)ε(ω).
func CircuitTwin(m em.Dispersive, freqs []float64) (float64, error) {
	worst := 0.0
	for _, f := range freqs {
		src := m.CellNetlist(CellArea, CellGap, fmt.Sprintf(".ac lin 1 %g %g", f, f))
		d, err := netlist.ParseString(src)
		if err != nil {
			return 0, err
		}
		c, err := circuit.FromDeck(d)
		if err != nil {
			return 0, err
		}
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		i, ok := r.ElementCurrent("V1")
		if !ok {
			return 0, fmt.Errorf("no V1 current")
		}
		y := -i[0] // source current n+ → n− through the source; the cell draws the opposite
		want := m.CellAdmittance(CellArea, CellGap, 2*math.Pi*f)
		if d := cmplx.Abs(y-want) / cmplx.Abs(want); d > worst {
			worst = d
		}
	}
	return worst, nil
}

// Slap ramps the field from 0 to 1 over the rise time and returns the
// permittivity water shows at the end of the ramp (D/ε0E) and after a long
// wait (20 τ).
func Slap(m em.Dispersive, rise float64) (atEnd, later float64) {
	tau := m.Debye[0].Tau
	E := func(t float64) float64 {
		if t >= rise {
			return 1
		}
		return t / rise
	}
	// 1000 steps on the ramp, then 1000 more over 20 τ (the ramp's end is
	// a grid point, so the exponential integrator is exact).
	var ts []float64
	for k := 0; k <= 1000; k++ {
		ts = append(ts, rise*float64(k)/1000)
	}
	for k := 1; k <= 1000; k++ {
		ts = append(ts, rise+20*tau*float64(k)/1000)
	}
	d := m.ResponseAt(E, ts)
	return d[1000], d[len(d)-1]
}

func main() {
	w := em.WaterDebye25
	fmt.Println("1. Water: the rate of change lives in the material (standard physics)")
	fmt.Printf("   %s — %s\n", w.Name, w.Source)
	fmt.Printf("   %12s %10s %10s %10s\n", "frequency", "ε′", "ε″", "loss tan")
	for _, f := range []float64{1e6, 1e8, 1e9, 1e10, 19.2e9, 5e10, 1e11} {
		e := w.Eps(2 * math.Pi * f)
		fmt.Printf("   %12.3g %10.2f %10.2f %10.3f\n", f, real(e), -imag(e), w.LossTangent(2*math.Pi*f))
	}
	fmt.Println("   (single Debye: valid below about 100 GHz)")

	fmt.Println("\n   'Slapping' the water with an electric field ramped up at different speeds:")
	for _, rise := range []float64{1e-6, 1e-9, 1e-10, 1e-11, 1e-12} {
		a, b := Slap(w, rise)
		fmt.Printf("   rise %8.3g s: water shows ε = %5.1f at the end of the ramp, %5.1f after it settles\n", rise, a, b)
	}
	fmt.Println("   Slow: the molecules follow (ε ≈ 78.6). Fast: they cannot turn in time and the water answers like a stiff solid (ε → 5.4).")

	worst, err := CircuitTwin(w, []float64{1e6, 1e9, 19.2e9, 1e11})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("\n   Circuit twin: water cell (1 cm², 1 mm) as C∞ ∥ (R–C) solved by kronspice AC vs ε(ω): worst difference %.2g\n", worst)

	fmt.Println("\n2. Bopp–Podolsky electrodynamics (HYPOTHETICAL, H19): higher derivatives in the field equations")
	fmt.Println("   (1 + a²□)□A = μ0J; a = 0 is Maxwell. Point charge: φ = q/(4πε0 r)(1 − e^(−r/a)).")
	const q = 1.602176634e-19
	fmt.Printf("   %10s  %16s  %16s  %s\n", "r / a", "Maxwell φ [V]", "Podolsky φ [V]", "fraction below Coulomb")
	a := 1e-15
	for _, x := range []float64{0, 0.1, 1, 3, 10, 30} {
		r := x * a
		var mx string
		if r == 0 {
			mx = "infinite"
		} else {
			mx = fmt.Sprintf("%.6g", em.PodolskyPotential(q, r, 0))
		}
		fmt.Printf("   %10.1f  %16s  %16.6g  %.3g\n", x, mx, em.PodolskyPotential(q, r, a), em.PodolskyRelativeDeviation(r, a))
	}
	fmt.Println("   (a = 1 fm for illustration; the effect is gone a few a away, so only very short-range experiments can see it)")
	// Electron self-energy equal to its rest energy when a = r_e/2 (Podolsky's motivation).
	const mec2 = 8.1871057769e-14 // J (CODATA 2018)
	aSelf := q * q / (8 * math.Pi * em.Eps0 * mec2)
	fmt.Printf("   electron field energy: Maxwell infinite; Podolsky q²/(8πε0 a) = m_e c² when a = %.4g m\n", aSelf)
	mx, ms := em.PodolskyModes(0, a)
	fmt.Printf("   free waves: the Maxwell wave (ω = ck, here %.3g at k = 0) and a massive mode starting at ω = c/a = %.3g rad/s (%.3g Hz for a = 1 fm)\n", mx, ms, ms/(2*math.Pi))
	fmt.Printf("   a Coulomb-law test good to 1 part in 10⁶ at 1 cm only bounds a < %.3g m; electron–positron scattering (Bufalo et al., PRD 90, 085012, 2014) bounds it far more tightly\n",
		em.PodolskyBound(0.01, 1e-6))
	fmt.Println("   Status: untested hypothesis with the standard prediction (a = 0) beside every result.")
}
