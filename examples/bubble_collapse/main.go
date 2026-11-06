// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command bubble_collapse puts the single sonoluminescing bubble on the
// centre node of examples/acoustic_drive's flask and compares its collapse
// with the knot's implosion in examples/knotted_light (H29, H30).
//
// The bubble is an equivalent circuit (CLAUDE.md §1 non-goals allow the
// bubble back in this form): the liquid's inertance M(V) = ρ/4πR in
// parallel with the pulsating sphere's radiation resistance ρc/4πR², the
// wall's viscous resistance, and the argon as a nonlinear compliance with a
// hard core. With M depending on the volume it is a Kron inductance on a
// coordinate (Christoffel term ½M′UU_M); in the incompressible, inviscid
// limit it is the Rayleigh–Plesset equation, and Keller–Miksis runs beside
// it as the twin. Its ledger closes by construction.
//
//	go run ./examples/bubble_collapse          # about 1 min
//	go run ./examples/bubble_collapse -nofield # skip the field-solver part
//
// Parts:
//  1. the bubble over one acoustic cycle, ledger, Keller–Miksis twin;
//  2. the bubble as a load on the drive chain (harmonic balance at the
//     fundamental), with kronspice's ledger;
//  3. reversibility: the knot's implosion is the explosion reversed; the
//     bubble's is only without its two losses;
//  4. the light: a free–free (bremsstrahlung) estimate against the
//     measured flash;
//  5. helicity: does the flash carry the knot's topology? Circular
//     polarization measured on a probe sphere in em/fdtd.go for the
//     hopfion (positive control, against its exact field), a rotating
//     dipole (circular on axis, no net helicity) and a thermal core of
//     random currents.
package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	"github.com/fdamador/kronspice/output"
)

const atm = 101325.0

func main() {
	data := flag.String("data", "", "write the report's data series as JSON to this file and stop")
	rep := flag.String("report", "", "write the interactive report (HTML) to this file and stop")
	nofield := flag.Bool("nofield", false, "skip part 5 (field solver)")
	flag.Parse()
	if *rep != "" {
		if err := writeReport(*rep); err != nil {
			fmt.Fprintln(os.Stderr, "bubble_collapse:", err)
			os.Exit(1)
		}
		return
	}
	if *data != "" {
		if err := writeData(*data); err != nil {
			fmt.Fprintln(os.Stderr, "bubble_collapse:", err)
			os.Exit(1)
		}
		return
	}
	if err := run(!*nofield); err != nil {
		fmt.Fprintln(os.Stderr, "bubble_collapse:", err)
		os.Exit(1)
	}
}

// flashSpectrum is the model flash's spectrum on 250–700 nm in 5 nm steps.
func flashSpectrum(b Bubble, s []Sample) ([]float64, []float64) {
	var ls []float64
	for k := 0; k <= 90; k++ {
		ls = append(ls, 250e-9+float64(k)*5e-9)
	}
	return ls, b.Spectrum(s, ls)
}

func run(field bool) error {
	d0, err := readDrive(0, 1.3*atm)
	if err != nil {
		return err
	}
	b := DefaultBubble(d0.FPeak)
	fmt.Println("bubble_collapse: a sonoluminescing bubble on acoustic_drive's flask, against knotted_light")
	fmt.Printf("argon bubble R0 = %.1f µm in water at 300 K, Pa = %.2f atm at %.2f kHz (the matched flask's peak)\n",
		b.R0*1e6, b.Pa/atm, b.F/1e3)
	fmt.Printf("  ρ = %.3f kg/m³, c = %.2f m/s (IAPWS-95), μ = %.1f µPa·s, σ = %.2f mN/m (IAPWS), γ = 5/3, R0/h = %.2f\n\n",
		b.Rho, b.C, b.Mu*1e6, b.Sigma*1e3, b.Core)

	// Part 1
	cy, samples, err := b.Run(6)
	if err != nil {
		return err
	}
	kmax, kmin, ktmin, kmach, err := b.RunKM(6)
	if err != nil {
		return err
	}
	T := 1 / b.F
	fmt.Println("Part 1: one acoustic cycle (the 6th from rest), equivalent circuit vs Keller–Miksis")
	fmt.Printf("  %-34s %14s %14s\n", "", "circuit", "Keller–Miksis")
	fmt.Printf("  %-34s %14.2f %14.2f\n", "R_max [µm]", cy.Rmax*1e6, kmax*1e6)
	fmt.Printf("  %-34s %14.3f %14.3f\n", "R_min [µm]", cy.Rmin*1e6, kmin*1e6)
	fmt.Printf("  %-34s %14.4f %14.4f\n", "minimum at phase [cycles]", math.Mod(cy.TRmin, T)/T, math.Mod(ktmin, T)/T)
	fmt.Printf("  %-34s %14.3f %14.3f\n", "max wall Mach |Ṙ|/c", cy.MaxMach, kmach)
	fmt.Printf("  peak gas: %.3g Pa, %.0f K (adiabatic); R_max/R_min = %.1f\n", cy.Pmax, cy.Tmax, cy.Rmax/cy.Rmin)
	fmt.Printf("  ledger over the cycle: in from the sound field %s = Δstored %s + viscous %s + radiated sound %s; residual %.1e of in\n",
		output.Eng(cy.In, "J"), output.Eng(cy.DStored, "J"), output.Eng(cy.Viscous, "J"), output.Eng(cy.Radiated, "J"), cy.Residual/cy.In)
	fmt.Printf("  collapse energy P0·V_max = %s; %.0f %% of the radiated sound leaves within ±%.0f ns of the minimum\n",
		output.Eng(cy.CollapseWork, "J"), 100*cy.BurstFrac, BurstWin*1e9)
	uac := b.Pa * b.Pa / (2 * b.Rho * b.C * b.C)
	vmin := 4 * math.Pi / 3 * cy.Rmin * cy.Rmin * cy.Rmin
	fmt.Printf("  energy density: sound field %.2f J/m³ → collapse energy in V_min %.2e J/m³, ×%.1e\n\n",
		uac, cy.CollapseWork/vmin, cy.CollapseWork/vmin/uac)

	// Part 2
	d1, err := readDrive(cy.Y1, 1.3*atm)
	if err != nil {
		return err
	}
	fmt.Println("Part 2: the bubble as a load on the drive chain (fundamental admittance at node p)")
	fmt.Printf("  Y₁ = %.3e %+.3ej m³/(Pa·s): ½·Re(Y₁)·Pa² = %s, cycle work × f = %s\n",
		real(cy.Y1), imag(cy.Y1), output.Eng(0.5*real(cy.Y1)*b.Pa*b.Pa, "W"), output.Eng(cy.In*b.F, "W"))
	fmt.Printf("  %-30s %12s %12s\n", "", "no bubble", "with bubble")
	fmt.Printf("  %-30s %12.3f %12.3f\n", "peak [kHz]", d0.FPeak/1e3, d1.FPeak/1e3)
	fmt.Printf("  %-30s %12.2f %12.2f\n", "pressure per volt [kPa/V]", d0.PerVolt/1e3, d1.PerVolt/1e3)
	fmt.Printf("  %-30s %12.4f %12.4f\n", "generator EMF for 1.3 atm [V]", d0.EMF, d1.EMF)
	fmt.Printf("  %-30s %12s %12s\n", "source power", output.Eng(d0.Source, "W"), output.Eng(d1.Source, "W"))
	fmt.Printf("  %-30s %12s %12s\n", "into the liquid (mode)", output.Eng(d0.Mode, "W"), output.Eng(d1.Mode, "W"))
	fmt.Printf("  %-30s %12s %12s\n", "into the bubble", "—", output.Eng(d1.Bubble, "W"))
	fmt.Printf("  kronspice ledger residual %.1e / %.1e\n\n", d0.Residual, d1.Residual)

	// Part 3
	var tMax float64
	rm := 0.0
	for _, s := range samples {
		if s.T < math.Mod(cy.TRmin, T) && s.R > rm {
			rm, tMax = s.R, s.T
		}
	}
	tStart := 5*T + tMax
	half := cy.TRmin - tStart
	fmt.Println("Part 3: reversibility. From R_max through the collapse and back out (2 × the fall time),")
	fmt.Println("  then the momentum flipped and the drive run backwards: does the bubble return?")
	for _, l := range []Losses{{false, false}, {true, false}, {false, true}, {true, true}} {
		dR, dPi, err := b.Reversal(l, tStart, half, nil)
		if err != nil {
			return err
		}
		fmt.Printf("  viscous %-5v radiation %-5v  |ΔR|/R = %.1e   |ΔΠ|/Π_max = %.1e\n", l.Viscous, l.Radiation, dR, dPi)
	}
	fmt.Println()

	// Part 4
	L := b.FlashEstimate(samples)
	fmt.Println("Part 4: the light (free–free emission from Saha-ionized argon, uniform core, optically thin)")
	fmt.Printf("  energy per flash %s (%.1e of the collapse energy), peak %s, ionization %.1e at the peak\n",
		output.Eng(L.Energy, "J"), L.Energy/cy.CollapseWork, output.Eng(L.Peak, "W"), L.PeakIon)
	fmt.Printf("  flash FWHM %.0f ps; about %.1e photons (energy / kT at the peak)\n", L.FWHM*1e12, L.Photons)
	fmt.Printf("  measured: 1e5–1e7 photons per flash (≈ 1e-12 J), widths 60–250 ps (Gompf et al. 1997)\n")
	ls, spec := flashSpectrum(b, samples)
	mx := 0.0
	for _, v := range spec {
		mx = math.Max(mx, v)
	}
	fmt.Printf("  spectrum (time-integrated energy per nm, relative):")
	for k := 0; k < len(ls); k += 15 {
		fmt.Printf("  %.0f nm %.2f", ls[k]*1e9, spec[k]/mx)
	}
	fmt.Println()
	fmt.Printf("  Planck fit of its shape (250–700 nm): %.0f K; FWHM through a UV band (250–300 nm) %.0f ps, red band (600–700 nm) %.0f ps\n",
		PlanckFit(ls, spec), b.BandWidth(samples, 250e-9, 300e-9)*1e12, b.BandWidth(samples, 600e-9, 700e-9)*1e12)
	fmt.Printf("  measured: a featureless continuum rising into the UV to water's edge near 200 nm, Planck fits of order 1e4 K,\n")
	fmt.Printf("  and the same width in red and UV (Gompf et al. 1997)\n\n")

	if !field {
		return nil
	}
	// Part 5
	fmt.Println("Part 5: helicity of the flash. Degree of circular polarization P_c = ΣV/ΣI over a probe sphere")
	h, err := DefaultHelicity().RunAll()
	if err != nil {
		return err
	}
	printHelicity(h)
	return nil
}
