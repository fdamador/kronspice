// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command hydrogen_poles tests hypothesis H34 (CLAUDE.md §1, user
// hypothesis 2026-11-05): Walter Russell's "four pairs of poles — not one
// pair" (painting "Analysis of the Octave Bar Magnet") read as hydrogen,
// and his two touching spheres as two hydrogen atoms.
//
// Two readings, each with the standard prediction computed beside it:
//
//  1. One atom. Poles come in powers of two: a dipole has one pair, an
//     octupole (2³ poles) four. A state of definite parity has no
//     magnetic quadrupole, and a magnetic 2^k-pole needs 2j ≥ k, so one
//     pair is all the ground state (1s, j = ½) can have, and four pairs
//     first become possible in 2p3⁄2. The electron's field is computed
//     from its orbital current and spin magnetization (Pauli spinors,
//     Clebsch–Gordan coupled), by the exact Legendre-series Green's
//     function; the proton's spin dipole is added. Moments are checked
//     against the Landé factor and the Wigner–Eckart theorem, the field
//     against a direct current-loop sum, and the pole count is read as the
//     bench would: B_r on spheres round the atom, B_z on a line beside it.
//
//  2. Two atoms, H₂: four small magnets, two electron spins and two proton
//     spins. Heitler–London gives the bond and shows which spin state
//     binds; the spin density and the field of the four magnets are
//     computed for each spin arrangement.
//
//     go run ./examples/hydrogen_poles               # about ten seconds
//     go run ./examples/hydrogen_poles -svg viz      # also figures and GIFs
package main

import (
	"flag"
	"fmt"
	"math"
	"path/filepath"
)

func main() {
	svgDir := flag.String("svg", "", "write figures (SVG) and animations (GIF) to this directory")
	flag.Parse()
	fmt.Println("hydrogen_poles (H34): Russell's \"four pairs of poles — not one pair\" read as hydrogen")

	atoms := AtomPart()
	mols := MoleculePart()
	spec, err := SpectrumPart(atoms)
	if err != nil {
		fmt.Println("spectrum:", err)
		return
	}

	if *svgDir != "" {
		if err := AtomFigures(*svgDir, atoms); err != nil {
			fmt.Println("figures:", err)
			return
		}
		if err := MoleculeFigures(*svgDir, mols); err != nil {
			fmt.Println("figures:", err)
			return
		}
		if err := SpectrumFigures(*svgDir, spec); err != nil {
			fmt.Println("figures:", err)
			return
		}
		fmt.Println("\nfigures and animations written to", filepath.Clean(*svgDir))
	}
	fmt.Println("\nStandard prediction: ground-state hydrogen has one pair of poles (a dipole); four pairs (an octupole)")
	fmt.Println("first appear in the 2p3⁄2 level. In H₂ the bond pairs the electron spins, so only the protons' tiny dipoles remain.")
}

// AtomResult is one level's computed field and scans.
type AtomResult struct {
	F           AtomField
	Norm        float64
	M1, M3      float64
	Even        float64 // largest |m2|, |m4| (zero by parity)
	M5          float64
	Rc          float64
	Radii       []float64
	Zones       []int
	LineChanges map[float64]int
}

// AtomPart runs Parts 1–2 and prints them.
func AtomPart() []AtomResult {
	fmt.Println("\n1. One hydrogen atom: the electron's orbital current + spin magnetization (+ the proton's spin dipole)")
	fmt.Println("   Pauli spinors |n l j m⟩ (Clebsch–Gordan), g_s = 2.002319; field from the exact Legendre-series Green's function")
	var out []AtomResult
	for _, s := range States() {
		f := NewAtomField(s)
		r := AtomResult{F: f}
		r.Norm = Norm(s)
		r.M1, r.M3 = f.Moments()
		r.Even = math.Max(math.Abs(f.Ser.Moment(2)), math.Abs(f.Ser.Moment(4)))
		r.M5 = f.Ser.Moment(5)
		r.Rc = OctupoleRadius(r.M1, r.M3)
		out = append(out, r)
	}
	fmt.Printf("   %-16s %8s %12s %12s %14s %10s %10s\n", "level", "∫|ψ|²", "dipole μB", "−g_j m_j", "octupole μB·a0²", "even l", "l = 5")
	for _, r := range out {
		fmt.Printf("   %-16s %8.6f %12.6f %12.6f %14.6f %10.1e %10.1e\n", r.F.S.Name, r.Norm, r.F.Ser.Moment(1), -r.F.S.GJ()*r.F.S.M, r.M3, r.Even, math.Abs(r.M5))
	}
	o32, o12 := out[2].M3, out[3].M3
	fmt.Printf("   Wigner–Eckart: octupole(m = ½)/octupole(m = 3⁄2) = %.6f (theory −3: ⟨j m 3 0|j m⟩ ∝ m(5m² − 3j(j+1) + 1))\n", o12/o32)
	fmt.Println("   split by source (orbital current | spin magnetization), dipole and octupole:")
	for _, r := range out {
		o := NewSeries(r.F.S.JOrb, []int{1, 3}, 0.01, 8000, 48)
		p := NewSeries(r.F.S.JSpin, []int{1, 3}, 0.01, 8000, 48)
		fmt.Printf("     %-16s orbital %+9.6f, %+10.6f | spin %+9.6f, %+10.6f\n", r.F.S.Name, o.Moment(1), o.Moment(3), p.Moment(1), p.Moment(3))
	}
	fmt.Println("   the octupole comes entirely from the electron's spin magnetization (the l = 1 orbital current carries none); j = ½ levels")
	fmt.Println("   (1s, 2p½) have none in either part, as 2j ≥ 3 requires; even-l moments vanish by parity")
	fmt.Printf("   proton spin dipole added: +%.4e μB (m_I = +½, along +z)\n", MuPr)
	tw := TwinAtom(out[3].F)
	fmt.Printf("   twin: series field against a direct sum of %d current loops, 4 points outside the cloud: worst %.1e of |B|\n", tw.N, tw.Err)

	fmt.Println("\n2. The pole count, read as a bench would: N and S zones of B_r on spheres round the atom")
	radii := LogRadii(0.5, 60, 240)
	for i := range out {
		f := out[i].F
		out[i].Radii = radii
		out[i].Zones = ZoneScan(f.Br, radii, 720)
		fmt.Printf("   %-16s", f.S.Name)
		for _, b := range Bands(radii, out[i].Zones) {
			fmt.Printf("  %d zones %.1f–%.1f a0;", b.Zones, b.From, b.To)
		}
		if out[i].Rc > 0 {
			fmt.Printf("  exterior octupole radius r_c = %.2f a0", out[i].Rc)
		}
		fmt.Println()
	}
	fmt.Println("   two zones = one pair (dipole); four zones = the axial octupole's N S N S bands. Outside the cloud the octupole")
	fmt.Println("   leads out to r_c, beyond which the dipole's two zones take over (the octupole falls as 1/r⁵, the dipole as 1/r³).")
	fmt.Println("   for m = 3⁄2, r_c lies inside the electron cloud, where the field is not the exterior expansion: no four zones appear.")
	fmt.Println("   the bar's bench scan (examples/octave_magnet) on the atom: sign changes of B_z on a line beside it, z from −60 to 60 a0")
	for i := range out {
		out[i].LineChanges = map[float64]int{}
		fmt.Printf("   %-16s", out[i].F.S.Name)
		for _, rho := range []float64{2, 5, 10, 20} {
			f := out[i].F
			_, b := LineScan(func(rho, z float64) float64 { _, bz := f.Cyl(rho, z); return bz }, rho, 60, 2401)
			n := SignChanges(b)
			out[i].LineChanges[rho] = n
			fmt.Printf("  ρ = %2.0f a0: %d", rho, n)
		}
		fmt.Println()
	}
	fmt.Println("   a dipole gives 2 sign changes beside it (one pair, as for the bar); more mean more pole sites along the line")
	return out
}

// Norm is ∫|ψ|² dV.
func Norm(s State) float64 {
	x, w := GaussLegendre(48)
	var n float64
	for i := 0; i < 8000; i++ {
		r := (float64(i) + 0.5) * 0.01
		for q := range x {
			n += s.Density(r, math.Acos(x[q])) * w[q] * r * r * 0.01 * 2 * math.Pi
		}
	}
	return n
}

// Twin is a field-twin result.
type Twin struct {
	N   int
	Err float64
}

// TwinAtom compares the series field with a direct loop sum (h = 0.05 a0).
func TwinAtom(f AtomField) Twin {
	lp := NewLoops(SphJ(f.S.J), 0.05, 40, -40, 40)
	worst := 0.0
	for _, p := range [][2]float64{{15, 0.3}, {25, 1.2}, {12, 2.0}, {30, 2.9}} {
		r, th := p[0], p[1]
		bp, bz, _ := lp.Field(r*math.Sin(th), r*math.Cos(th))
		br := bp*math.Sin(th) + bz*math.Cos(th)
		bt := bp*math.Cos(th) - bz*math.Sin(th)
		sr, st := f.Ser.Br(r, th), f.Ser.Bt(r, th)
		worst = math.Max(worst, math.Hypot(br-sr, bt-st)/math.Hypot(sr, st))
	}
	return Twin{len(lp.A), worst}
}

// MolResult is one spin arrangement of H₂.
type MolResult struct {
	F      MolField
	Spin   float64 // ∫ spin density
	Moment float64
	Radii  []float64
	Zones  []int
	Line   int
}

// MoleculePart runs Part 3 and prints it.
func MoleculePart() []MolResult {
	fmt.Println("\n3. Two atoms, H₂ (Heitler–London): four small magnets — two electron spins, two proton spins")
	re, de := HLMinimum()
	fmt.Printf("   singlet (spins opposite) binds: R_e %.4f a0 (%.4f Å), D_e %.4f eV   [measured 1.4011 a0, 4.748 eV]\n", re, re*A0*1e10, de*Hartree)
	fmt.Printf("   triplet (spins parallel) repels at every R: E − 2E_H = %+.3f eV at 1.4 a0, %+.3f eV at 2 a0, %+.4f eV at 4 a0\n",
		(HLEnergy(1.4, -1)+1)*Hartree, (HLEnergy(2, -1)+1)*Hartree, (HLEnergy(4, -1)+1)*Hartree)
	const R = 1.4011
	var out []MolResult
	for _, m := range Molecules(R) {
		f := NewMolField(m, 0.05, 10, 5)
		r := MolResult{F: f, Spin: SpinIntegral(m), Moment: m.Moment()}
		out = append(out, r)
	}
	fmt.Printf("   at R = %.4f a0 (measured bond), spins quantized along the bond:\n", R)
	c1, c2 := GridConvergence(R)
	fmt.Printf("   triplet spin field inside the cloud: corner grid h = 0.05 against h = 0.025, 6 points within 3 a0: worst %.1e of |B| (h = 0.1 against 0.05: %.1e)\n", c1, c2)
	radii := LogRadii(0.2, 20, 120)
	for i := range out {
		f := out[i].F
		out[i].Radii = radii
		br := func(r, th float64) float64 {
			bp, bz, _ := f.Field(r*math.Sin(th), r*math.Cos(th))
			return bp*math.Sin(th) + bz*math.Cos(th)
		}
		if f.M.Moment() == 0 {
			out[i].Zones = make([]int, len(radii))
		} else {
			out[i].Zones = ZoneScan(br, radii, 360)
		}
		bound := "bound"
		if !f.M.Bound {
			bound = "unbound (repulsive)"
		}
		fmt.Printf("   %-52s ∫ spin density %.4f, moment %+.5f μB, %s\n", f.M.Name, out[i].Spin, out[i].Moment, bound)
		if f.M.Moment() == 0 {
			fmt.Println("        no magnetic field at all: every spin expectation vanishes — no poles")
			continue
		}
		fmt.Print("       ")
		for _, b := range Bands(radii, out[i].Zones) {
			fmt.Printf(" %d zones %.2f–%.2f a0;", b.Zones, b.From, b.To)
		}
		_, bl := LineScan(func(rho, z float64) float64 { _, bz, _ := f.Field(rho, z); return bz }, 2, 20, 801)
		out[i].Line = SignChanges(bl)
		fmt.Printf("  B_z beside it (ρ = 2 a0): %d sign changes\n", out[i].Line)
	}
	fmt.Println("   the bound molecule (singlet) has no electron spin density anywhere, so its electrons give no poles; only the")
	fmt.Printf("   protons' dipoles (%.2e μB each, ortho) remain. All four magnets aligned needs the triplet, which does not bind.\n", MuPr)
	return out
}

// SpinIntegral is ∫ spin density dV on the loop grid.
func SpinIntegral(m Molecule) float64 {
	var n float64
	h := 0.02
	for i := 0; i < 600; i++ {
		rho := (float64(i) + 0.5) * h
		for k := 0; k < 1200; k++ {
			z := -12 + (float64(k)+0.5)*h
			n += m.SpinDensity(rho, z) * 2 * math.Pi * rho * h * h
		}
	}
	return n
}

// GridConvergence compares the triplet electron field at six corner points
// inside the cloud for grids h = 0.05 vs 0.025 and h = 0.1 vs 0.05.
func GridConvergence(R float64) (fine, coarse float64) {
	m := Molecules(R)[2]
	pts := [][2]float64{{0.3, 0.2}, {0.5, 0.9}, {1.0, -0.4}, {0.1, 1.6}, {2.0, 1.0}, {0.7, -2.2}}
	field := func(h float64) [][2]float64 {
		l := NewLoops(m.J, h, 10, -10, 10)
		out := make([][2]float64, len(pts))
		parallel(len(pts), func(i int) {
			a, b, _ := l.Field(pts[i][0], pts[i][1])
			out[i] = [2]float64{a, b}
		})
		return out
	}
	f1, f2, f4 := field(0.1), field(0.05), field(0.025)
	for i := range pts {
		n := math.Hypot(f4[i][0], f4[i][1])
		fine = math.Max(fine, math.Hypot(f2[i][0]-f4[i][0], f2[i][1]-f4[i][1])/n)
		coarse = math.Max(coarse, math.Hypot(f1[i][0]-f2[i][0], f1[i][1]-f2[i][1])/n)
	}
	return
}
