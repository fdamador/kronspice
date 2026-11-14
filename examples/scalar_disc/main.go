// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command scalar_disc: step 4 of H47 (doc/h47-scalar-gravity.md). The
// optically gated rotating aluminium disc proposed by Iadicicco & Verolino
// (arXiv:2608.21403) as a falsifiable test of the observer-split source of
// the Aharonov–Bohm scalar S (Iadicicco, Modanese & Verolino, GRG 58, 106,
// 2026), with S from em/scalar's field equations in the quasi-static limit
// and its gravitational effect from Minotti & Modanese's weak-field
// coupling (EPJC 86, 593, 2026). Standard prediction beside every number:
// S = 0, no force beyond Newton and Coulomb.
//
//	go run ./examples/scalar_disc
package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
)

const (
	gNewton = 6.67430e-11 // CODATA 2018 [m³ kg⁻¹ s⁻²]
	// declared bench limits
	airBreakdown = 3e6         // [V/m], surface field the induced charge may reach in air
	alUTS        = 310e6       // 6061-T6 tensile strength [Pa]
	alRho        = 2700.0      // [kg/m³]
	alNu         = 0.33        // Poisson's ratio
	discR        = 0.05        // disc radius [m]
	sensorZ      = 0.02        // sensor height above the lit track [m]
	fMinDeclared = 1e-9        // force per mass a good bench can resolve [m/s²] (declared)
	gridHalf     = 0.08        // half-width of the S² box [m]
	gridH        = 0.004       // its spacing [m]
	revSamples   = 72          // sensor samples per revolution
	coulombQTest = 1e-12       // test charge on the sensor for the Coulomb systematic [C]
	sensorMass   = 1e-3        // sensor mass [kg]
	_            = gNewton * 0 // keep the constant documented
)

func main() {
	d := DefaultDisc()
	cp := MMFit()
	soft := d.soft()
	sg := d.Sigma()
	fmt.Println("H47 step 4: optically gated rotating Al disc (Iadicicco & Verolino, arXiv:2608.21403)")
	fmt.Printf("declared: lit ring r = %.0f mm, width %.0f mm, sector FWHM %.0f° (σ = %.4f rad), Q = %.0f nC, %.0f rpm, sensor %.0f mm above the track\n\n",
		d.RingR*1e3, d.RingW*1e3, d.FWHMDeg, sg, d.Q*1e9, d.Omega*60/(2*math.Pi), sensorZ*1e3)

	src := d.Sources(0)
	var smax, imax, net, abs float64
	for _, p := range src {
		smax = math.Max(smax, p.sigma)
		imax = math.Max(imax, math.Abs(p.source))
		net += p.source * p.area
		abs += math.Abs(p.source) * p.area
	}
	fmt.Println("1. The source (IMV eqs. 13, B3, 24)")
	fmt.Printf("   peak σ_s %.4g C/m²; peak |I_G| %.4g C/(m² s) against |Ω|σ_max e^(−½)/σ = %.4g; ring integral %.2g of Σ|I|A\n\n",
		smax, imax, math.Abs(d.Omega)*smax*math.Exp(-0.5)/sg, net/abs)

	// S at a sensor fixed in the lab while the sector turns: by symmetry,
	// the sector at φc = Ωt seen from (r, 0, z) is the sector at 0 seen from
	// angle −Ωt.
	sensor := func(phi float64) [3]float64 {
		return [3]float64{d.RingR * math.Cos(phi), d.RingR * math.Sin(phi), sensorZ}
	}
	var sPeak, ePeak float64
	for k := 0; k < revSamples; k++ {
		x := sensor(-2 * math.Pi * float64(k) / revSamples)
		sPeak = math.Max(sPeak, math.Abs(SAt(src, x, soft)))
		e := EAt(src, x, soft)
		ePeak = math.Max(ePeak, math.Sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]))
	}
	rev := d
	rev.Omega = -d.Omega
	srcRev := rev.Sources(0)
	uni := d
	uni.Uniform = true
	srcUni := uni.Sources(0)
	xTest := sensor(-0.1)
	fmt.Println("2. S at the sensor (□S = μ₀I_G, quasi-static; em/scalar's kernel)")
	fmt.Printf("   peak |S| over a revolution %.4g T (c|S| = %.3g V/m)\n", sPeak, em.SpeedOfLight*sPeak)
	fmt.Printf("   controls: reversed rotation S = %.4g × forward; whole ring lit S = %.2g T; standard prediction S = 0\n\n",
		SAt(srcRev, xTest, soft)/SAt(src, xTest, soft), SAt(srcUni, xTest, soft))

	fmt.Println("3. Gravity from S² (Minotti & Modanese 2026, K_φ = 1e-16 A² N⁻², C_φ = 6e21, their fit)")
	grid := SGrid(src, gridHalf, gridH, soft)
	var fPeak, fMean [3]float64
	for k := 0; k < revSamples; k++ {
		f := grid.Force(sensor(-2*math.Pi*float64(k)/revSamples), cp, gridH/2)
		for a := 0; a < 3; a++ {
			fMean[a] += f[a] / revSamples
			if math.Abs(f[a]) > math.Abs(fPeak[a]) {
				fPeak[a] = f[a]
			}
		}
	}
	gridR := SGrid(srcRev, gridHalf, gridH, soft)
	fF := grid.Force(xTest, cp, gridH/2)
	fR := gridR.Force(xTest, cp, gridH/2)
	big := SGrid(src, 1.5*gridHalf, gridH, soft)
	fBig := big.Force(xTest, cp, gridH/2)
	fmt.Printf("   force per mass at the sensor: peak (%.3g, %.3g, %.3g) m/s², mean over a turn (%.3g, %.3g, %.3g) m/s²\n",
		fPeak[0], fPeak[1], fPeak[2], fMean[0], fMean[1], fMean[2])
	fmt.Printf("   reversed rotation: f_z %.4g of forward (S² is even in Ω); box ×1.5: f_z changes %.2g %%\n",
		fR[2]/fF[2], 100*(fBig[2]/fF[2]-1))
	fmt.Printf("   standard prediction: 0 (Newton's pull of the disc is static and is subtracted by ON/OFF gating)\n\n")

	// limits of the bench
	sector := d.RingW * d.RingR * d.FWHMDeg * math.Pi / 180 * 1.0645 // Gaussian area ≈ FWHM × 1.0645
	qMax := em.Eps0 * airBreakdown * sector
	wMax := math.Sqrt(8 * alUTS / ((3 + alNu) * alRho * discR * discR))
	scale := (qMax / d.Q) * (wMax / d.Omega)
	fMax := math.Abs(fPeak[2]) * scale * scale
	fmt.Println("4. At the bench's limits")
	fmt.Printf("   charge: air breaks down at %.0f MV/m, so the sector holds ≤ %.3g C (%.1f× the declared 1 nC)\n", airBreakdown/1e6, qMax, qMax/d.Q)
	fmt.Printf("   speed: a %.0f mm Al disc bursts at ω = %.3g rad/s (%.0f rpm, %.1f× 3000 rpm)\n", discR*1e3, wMax, wMax*60/(2*math.Pi), wMax/d.Omega)
	fmt.Printf("   S scales as QΩ → peak |S| %.3g T; the S² force as (QΩ)² → %.3g m/s²\n", sPeak*scale, fMax)
	fmt.Printf("   the superconducting junction S of the 2026 fit is ≈ 2e3 T: %.1g× the disc's best\n\n", 2e3/(sPeak*scale))

	// what a null result bounds: f ∝ K_φC_φ with K_φ ∝ C_φ at fixed ω₀, so f ∝ C_φ².
	twoW3 := 8 * math.Pi * gNewton * em.Eps0 * cp.Cphi / (em.SpeedOfLight * em.SpeedOfLight * cp.Kphi)
	cNeed := cp.Cphi * math.Sqrt(fMinDeclared/fMax)
	fmt.Println("5. What the disc can decide")
	fmt.Printf("   the fit implies 2ω₀ + 3 = 8πG₀ε₀C_φ/(c²K_φ) = %.3g, so ω₀ ≈ %.2g; Cassini needs ω > 4e4 unless the matter coupling is suppressed (their App. B)\n", twoW3, (twoW3-3)/2)
	fmt.Printf("   to reach %.0e m/s² at the limits needs C_φ ≈ %.2g (%.2g× the fit): a null at that sensitivity bounds C_φ there\n", fMinDeclared, cNeed, cNeed/cp.Cphi)
	fmt.Printf("   IMV's linear response Λ_opt (response ∝ S, odd in Ω) has no value in either paper: the disc measures it, S = %.3g T per (nC·krpm)\n",
		sPeak/(d.Q*1e9*d.Omega*60/(2*math.Pi)/1000))
	fmt.Println()

	fmt.Println("6. The ordinary systematic beside it")
	fmt.Printf("   Coulomb field of the gated charge at the sensor: peak %.3g V/m; on a %.0f g sensor carrying %.0f pC: %.3g m/s², %.1g× the S² force\n",
		ePeak, sensorMass*1e3, coulombQTest*1e12, coulombQTest*ePeak/sensorMass, coulombQTest*ePeak/sensorMass/math.Abs(fPeak[2]))
	fmt.Println("   it is odd in Q, follows the sector at Ω and does not reverse with rotation; S is odd in Ω and vanishes for the whole ring lit")
}
