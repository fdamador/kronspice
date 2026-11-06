// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
)

// ringTestPoints are field points outside the magnet: beside the outer face
// (on a pole, between poles, above the mid-plane), over the top face, and in
// the bore.
func ringTestPoints(r Ring) [][3]float64 {
	at := func(rho, phi, z float64) [3]float64 { return [3]float64{rho * math.Cos(phi), rho * math.Sin(phi), z} }
	return [][3]float64{
		at(r.Ro+3e-3, 0, 0), at(r.Ro+3e-3, 0.3, 2e-3), at(r.Ro+1e-3, 0.1, 0),
		at(0.9*r.Ro, 0.11, r.H), at(0.5*r.Ri, 0.38, 1e-3), at(r.Ro+8e-3, 1.0, -6e-3),
	}
}

// TwinError is the largest |B_charge − B_dipole|/|B| over the test points.
func TwinError(r Ring) float64 {
	acc, dip := r.Charges(0.25e-3, 2), r.Dipoles(0.25e-3)
	worst := 0.0
	for _, p := range ringTestPoints(r) {
		a, d := BCharges(acc, p), BDipoles(dip, p)
		worst = math.Max(worst, math.Sqrt(sq(a[0]-d[0])+sq(a[1]-d[1])+sq(a[2]-d[2]))/math.Sqrt(sq(a[0])+sq(a[1])+sq(a[2])))
	}
	return worst
}

func sq(x float64) float64 { return x * x }

// DivCurl is |∇·B| and |∇×B| at p by central differences, each divided by
// the norm of the gradient tensor (both vanish outside the magnet).
func DivCurl(cs []Charge, p [3]float64) (div, curl float64) {
	const d = 2e-6
	var J [3][3]float64 // J[i][j] = ∂B_i/∂x_j
	for j := 0; j < 3; j++ {
		a, b := p, p
		a[j] += d
		b[j] -= d
		ba, bb := BCharges(cs, a), BCharges(cs, b)
		for i := 0; i < 3; i++ {
			J[i][j] = (ba[i] - bb[i]) / (2 * d)
		}
	}
	var n float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			n += sq(J[i][j])
		}
	}
	n = math.Sqrt(n)
	div = math.Abs(J[0][0]+J[1][1]+J[2][2]) / n
	curl = math.Sqrt(sq(J[2][1]-J[1][2])+sq(J[0][2]-J[2][0])+sq(J[1][0]-J[0][1])) / n
	return
}

func degs(v []float64) string {
	s := ""
	for _, x := range v {
		s += fmt.Sprintf(" %.0f°", x*180/math.Pi)
	}
	return s
}

// ringPart prints Part 5.
func ringPart() error {
	rad, tan, ctl := RadialTwelve(), TangentialTwelve(), RadialUniform()
	fmt.Println("\n5. The octave bar bent into a ring: 12 sections of 30°, each magnetized opposite to its neighbours")
	fmt.Printf("   ring %.0f mm across, %.0f mm bore, %.0f mm tall, Br %.2f T (declared); each section uniformly magnetized\n", 2*rad.Ro*1e3, 2*rad.Ri*1e3, rad.H*1e3, Br)
	fmt.Println("   field from the surface charge σ = M·n̂ (outer and inner arcs, cut faces); twin: the volume as point dipoles M·dV")
	for _, r := range []Ring{rad, tan} {
		acc := r.Charges(0.25e-3, 2)
		var q float64
		for _, c := range acc {
			q += math.Abs(c.Q)
		}
		var net float64
		for _, c := range acc {
			net += c.Q
		}
		div, curl := 0.0, 0.0
		for _, p := range ringTestPoints(r)[:3] {
			d, c := DivCurl(acc, p)
			div, curl = math.Max(div, d), math.Max(curl, c)
		}
		fmt.Printf("   %-11s charge vs dipole twin: worst %.1e of |B| at 6 points; net charge %.1e of Σ|Q|; ∇·B %.0e, ∇×B %.0e (of |∇B|)\n",
			r.Short, TwinError(r), math.Abs(net)/q, div, curl)
	}
	fmt.Println("   B_r round the ring, 3 mm outside, mid-plane (the bench scan):")
	for _, r := range []Ring{rad, tan, ctl} {
		cs := r.Charges(0.5e-3, 1)
		phi, br, _, _ := RingScan(cs, r.Ro+3e-3, 0, 720)
		mx := 0.0
		for _, v := range br {
			mx = math.Max(mx, math.Abs(v))
		}
		poles := " none (B_r > 0 all round)"
		if CyclicSignChanges(br) > 0 {
			poles = degs(PeakAngles(phi, br))
		}
		fmt.Printf("     %-58s %2d sign changes, peak |B_r| %.0f mT, N poles at%s\n", r.Name, CyclicSignChanges(br), mx*1e3, poles)
	}
	fmt.Println("   both alternating rings have 12 poles round the outside (6 N, 6 S). Radial: at the section centres;")
	fmt.Println("   along the ring: at the joints, where like poles meet (consequent poles, as in the four-segment bar).")
	fmt.Println("   all outward: the outer face is all N and the bore all S, so B_r beside it never changes sign.")

	p := DefaultPickup(rad)
	res, err := p.Run(rad.Charges(0.5e-3, 1), 360)
	if err != nil {
		return err
	}
	fmt.Printf("   pickup coil (%d turns, %.0f mm, %.1f mm outside; R %.2f Ω, L %.0f µH, load %.1f Ω), radial ring at %.0f rpm, in kronspice:\n",
		p.Turns, 2*p.A*1e3, (p.Rc-rad.Ro)*1e3, p.R, p.L*1e6, p.Load, p.RPM)
	fmt.Printf("     %d electrical cycles per revolution (%.0f Hz); open-circuit EMF %.2f V peak; load %.3f V rms, %.0f mW\n",
		res.Cycles, float64(res.Cycles)*p.RPM/60, res.PeakEMF, res.Vrms, res.Power*1e3)
	fmt.Printf("     ledger over two revolutions: in %.5e J (shaft, through the EMF) = heat %.5e J + Δstored %.1e J (residual %.1e)\n",
		res.Ein, res.Ediss, res.dE, math.Abs(res.Ein-res.Ediss-res.dE)/res.Ein)
	return nil
}
