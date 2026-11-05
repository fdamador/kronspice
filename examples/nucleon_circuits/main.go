// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command nucleon_circuits models the proton and the neutron as circuits of
// light, hypothesis H31 (CLAUDE.md §1, user hypothesis 2026-11-05): the
// light-ring picture of the electron (H23) and of the neutrino (H28),
// carried to the nucleons, beside the knot picture (H29) and the standard
// accounts.
//
// The light ring. Light of energy mc² going once round a ring of
// circumference λ = h/(mc) has radius R = ħ/(mc), the reduced Compton
// wavelength. A charge q carried round it at c is a current I = qc/(2πR)
// on an area πR², so its magnetic moment is IπR² = qcR/2 = qħ/(2m): the
// Dirac moment, exactly. Each ring is a circuit: an inductor (the ring's
// self-inductance, a thin ring of core radius R/2, μ0R(ln 16 − 7/4)) whose
// current is set by a current source; rings near each other are coupled by
// their mutual inductance (Maxwell's coaxial-circle formula,
// em.MutualCircles). kronspice books every ring's ½LI² and every pair's
// M·I1·I2 on the energy ledger, and the moments follow from I·πR².
//
// What a model has to reproduce (CODATA 2022; PDG 2024):
//
//	                 proton              neutron
//	charge           +e                  0
//	spin             ½                   ½
//	mass             938.272 MeV         939.565 MeV (1.293 MeV heavier)
//	magnetic moment  +2.7928 μN          −1.9130 μN
//	size             0.8408 fm (charge)  ⟨r²⟩ = −0.1155 fm²
//	lifetime         stable (> 10³⁴ y)   878.4 s (n → p e ν̄, 0.782 MeV)
//
// Parts: (1) the proton as one light ring; (2) the neutron as the light
// rings of its decay products, proton + electron + antineutrino (the
// fdamador neutrino ring of H28); (3) the standard constituent-quark
// picture written as three light rings, uud and udd; (4) the knot picture:
// Skyrme's model, where proton and neutron are one knot in two spin states
// (Adkins, Nappi & Witten 1983, quoted); (5) the decay's energy budget and
// lifetime.
package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Constants, SI (CODATA 2022) and particle data (PDG 2024).
const (
	E      = 1.602176634e-19  // C
	Hbar   = 1.054571817e-34  // J s
	C0     = 299792458.0      // m/s
	Eps0   = 8.8541878188e-12 // F/m
	Mu0    = 1.25663706127e-6 // N/A²
	MeV    = 1.602176634e-13  // J
	MpMeV  = 938.27208943     // proton
	MnMeV  = 939.56542194     // neutron
	MeMeV  = 0.51099895069    // electron
	MuP    = 2.79284734463    // proton moment / μN
	MuN    = -1.91304276      // neutron moment / μN
	MuE    = -1.00115965218   // electron moment / μB
	RpFm   = 0.84075          // proton charge radius, fm (CODATA 2022)
	Rn2Fm2 = -0.1155          // neutron mean-square charge radius, fm² (PDG)
	TauN   = 878.4            // neutron lifetime, s (PDG 2024 average)
	Vud    = 0.97373          // CKM element (PDG 2024)
	Lambda = 1.2754           // |g_A/g_V| (PDG 2024)
)

// MuN0 is the nuclear magneton eħ/(2m_p) [J/T]; MuB the Bohr magneton.
var (
	MuN0 = E * Hbar / (2 * MpMeV * MeV / (C0 * C0))
	MuB  = E * Hbar / (2 * MeMeV * MeV / (C0 * C0))
)

// Ring is one light ring: charge Q (in units of e), radius R [m], current
// direction Up (+1) or down (−1) along the axis, at height Z [m] on a
// common axis.
type Ring struct {
	Name string
	Q    float64
	R    float64
	Up   float64
	Z    float64
}

// LightRadius is ħ/(mc) for a rest energy in MeV: the ring that carries
// that energy as one wavelength of light [m].
func LightRadius(mev float64) float64 { return Hbar * C0 / (mev * MeV) }

// Current is the ring's current qc/(2πR) [A], signed by Up.
func (r Ring) Current() float64 { return r.Up * r.Q * E * C0 / (2 * math.Pi * r.R) }

// Moment is I·πR² [J/T].
func (r Ring) Moment() float64 { return r.Current() * math.Pi * r.R * r.R }

// SelfL is the self-inductance of a thin ring of core radius R/2, uniform
// current: μ0R(ln(8R/a) − 7/4) with 8R/a = 16 [H].
func (r Ring) SelfL() float64 { return Mu0 * r.R * (math.Log(16) - 1.75) }

// Mutual is the mutual inductance of two coaxial rings [H].
func Mutual(a, b Ring) float64 { return em.MutualCircles(a.R, b.R, math.Abs(a.Z-b.Z)) }

// Coulomb is the electrostatic energy of the rings [J]: each ring's own
// (thin charged ring of core R/2: q²/(8π²ε0R)·ln(8R/a)) plus each pair's
// exact coaxial-ring energy q1q2/(4πε0)·2K(k)/(π√((a+b)² + d²)),
// k² = 4ab/((a+b)² + d²).
func Coulomb(rings []Ring) float64 {
	var u float64
	for i, a := range rings {
		q := a.Q * E
		u += q * q / (8 * math.Pi * math.Pi * Eps0 * a.R) * math.Log(16)
		for _, b := range rings[i+1:] {
			d := a.Z - b.Z
			s := (a.R+b.R)*(a.R+b.R) + d*d
			k := math.Sqrt(4 * a.R * b.R / s)
			K, _ := em.EllipticKE(k)
			u += a.Q * b.Q * E * E / (4 * math.Pi * Eps0) * 2 * K / (math.Pi * math.Sqrt(s))
		}
	}
	return u
}

// Magnetic is the rings' magnetic energy from the analytic sum
// ½ΣL I² + Σ M I1 I2 [J].
func Magnetic(rings []Ring) float64 {
	var u float64
	for i, a := range rings {
		ia := a.Current()
		u += 0.5 * a.SelfL() * ia * ia
		for _, b := range rings[i+1:] {
			u += Mutual(a, b) * ia * b.Current()
		}
	}
	return u
}

// CircuitRun is the kronspice solution of a set of rings.
type CircuitRun struct {
	Stored   float64 // ledger stored energy at the end [J]
	Analytic float64 // ½ΣLI² + ΣMII [J]
	Residual float64 // ledger: largest relative residual
	Currents []float64
	Moment   float64 // Σ I·πR² from the solved currents [J/T]
}

// Build makes the rings as kronspice inductors, each fed by a current
// source ramped from 0 to the ring's current over 1 ns, coupled pairwise by
// their mutual inductances.
func Build(rings []Ring) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("light rings")
	for i, r := range rings {
		n := fmt.Sprintf("n%d", i)
		b.Inductor("L"+r.Name, n, "0", r.SelfL())
		b.Add(circuit.Element{Name: "I" + r.Name, Kind: circuit.CurrentSource,
			Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-9, r.Current(), 2e-9, r.Current()}}}, "0", n)
	}
	for i, a := range rings {
		for j := i + 1; j < len(rings); j++ {
			c := rings[j]
			k := Mutual(a, c) / math.Sqrt(a.SelfL()*c.SelfL())
			b.Mutual(fmt.Sprintf("K%s%s", a.Name, c.Name), "L"+a.Name, "L"+c.Name, k)
		}
	}
	return b.BuildTransient()
}

// Circuit builds the rings and runs a transient to 2 ns.
func Circuit(rings []Ring) (CircuitRun, error) {
	ckt, err := Build(rings)
	if err != nil {
		return CircuitRun{}, err
	}
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: 1e-11, Stop: 2e-9}, analysis.TranOptions{})
	if err != nil {
		return CircuitRun{}, err
	}
	run := CircuitRun{Stored: w.Ledger.StoredEnd, Analytic: Magnetic(rings), Residual: w.Ledger.MaxRelativeResidual}
	for _, r := range rings {
		i, ok := w.ElementCurrent("L" + r.Name)
		if !ok {
			return run, fmt.Errorf("no current for ring %s", r.Name)
		}
		cur := i[len(i)-1]
		run.Currents = append(run.Currents, cur)
		run.Moment += cur * math.Pi * r.R * r.R
	}
	return run, nil
}

// Config is one spin arrangement of a three-ring nucleon with its weight in
// the spin-½ state (the standard SU(6) quark wavefunction: for the proton
// u↑u↑d↓ with weight 2/3, u↑u↓d↑ and u↓u↑d↑ with 1/6 each).
type Config struct {
	W     float64
	Rings []Ring
}

// QuarkNucleon builds the three-ring arrangements of a nucleon: two quarks
// of charge q2 and one of charge q1 (uud: q2 = 2/3, q1 = −1/3; udd:
// q2 = −1/3, q1 = 2/3), each a light ring of radius R, stacked on the axis
// at heights −d, 0, +d.
func QuarkNucleon(q2, q1, R, d float64) []Config {
	ring := func(name string, q, up, z float64) Ring { return Ring{Name: name, Q: q, R: R, Up: up, Z: z} }
	return []Config{
		{2.0 / 3, []Ring{ring("a", q2, 1, -d), ring("b", q2, 1, 0), ring("c", q1, -1, d)}},
		{1.0 / 6, []Ring{ring("a", q2, 1, -d), ring("b", q2, -1, 0), ring("c", q1, 1, d)}},
		{1.0 / 6, []Ring{ring("a", q2, -1, -d), ring("b", q2, 1, 0), ring("c", q1, 1, d)}},
	}
}

// Weighted is the spin-weighted moment [μN] and electromagnetic energy
// (Coulomb + magnetic, the latter from kronspice) [MeV] of a nucleon, and
// the worst ledger residual.
func Weighted(cfg []Config) (mu, emMeV, resid float64, err error) {
	for _, c := range cfg {
		run, e := Circuit(c.Rings)
		if e != nil {
			return 0, 0, 0, e
		}
		mu += c.W * run.Moment / MuN0
		emMeV += c.W * (run.Stored + Coulomb(c.Rings)) / MeV
		resid = math.Max(resid, run.Residual)
	}
	return
}

// NeutronLifetime is the standard-model lifetime from V_ud and λ = g_A/g_V
// with the radiative corrections in the constant (Czarnecki, Marciano &
// Sirlin, PRL 120, 202002, 2018): τ = 4908.6 s / (|V_ud|²(1 + 3λ²)).
func NeutronLifetime(vud, lambda float64) float64 {
	return 4908.6 / (vud * vud * (1 + 3*lambda*lambda))
}

func main() {
	fmt.Println("nucleon_circuits (H31): the proton and the neutron as circuits of light")
	fmt.Println("A light ring: energy mc² as one wavelength round a ring of radius ħ/mc; charge q carried at c gives")
	fmt.Println("current qc/2πR and moment qħ/2m (the Dirac moment). Rings are kronspice inductors with mutual coupling.")

	// Part 1
	rp := LightRadius(MpMeV)
	p := Ring{Name: "p", Q: 1, R: rp, Up: 1}
	run, err := Circuit([]Ring{p})
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Println("\n1. The proton as one light ring")
	fmt.Printf("   radius ħ/(m_p c) = %.4f fm (measured charge radius %.5f fm: %.1f× larger)\n", rp*1e15, RpFm, RpFm/(rp*1e15))
	fmt.Printf("   current %.4g A; moment %.4f μN (kronspice currents; measured %.4f μN: %.2f× larger)\n", run.Currents[0], run.Moment/MuN0, MuP, MuP/(run.Moment/MuN0))
	fmt.Printf("   magnetic self-energy ½LI² %.4f MeV (ledger %.4f MeV, residual %.1e)\n", run.Analytic/MeV, run.Stored/MeV, run.Residual)
	rq := 2 * MuP * MuN0 / (E * C0)
	fmt.Printf("   a ring carrying +e at c with the measured moment needs radius %.4f fm, i.e. a light ring of %.2f MeV\n", rq*1e15, Hbar*C0/rq/MeV)
	fmt.Printf("   (= m_p c²/%.4f, about a third of the proton: one ring cannot be the whole proton)\n", MuP)

	// Part 2
	fmt.Println("\n2. The neutron as the light rings of its decay products: proton + electron + antineutrino")
	fmt.Printf("   charge: +e − e + 0 = 0 ✓    spin: three spin-½ rings can make ½ ✓ (Rutherford's 1920 p + e could not)\n")
	fmt.Printf("   mass: m_n − (m_p + m_e) = %+.4f MeV: the neutron holds more than its parts, so it can decay ✓\n", MnMeV-MpMeV-MeMeV)
	re := LightRadius(MeMeV)
	muE := MuE * MuB / MuN0
	fmt.Printf("   moment with each ring at its own size (proton as measured, electron ring %.1f fm): %+.1f μN (measured %.4f: %.0f× too large)\n",
		re*1e15, MuP+muE, MuN, (MuP+muE)/MuN)
	need := (MuP - MuN) * MuN0
	reN := 2 * need / (E * C0)
	el := Ring{Name: "e", Q: -1, R: reN, Up: 1}
	pq := Ring{Name: "p", Q: 1, R: rq, Up: 1}
	run2, err := Circuit([]Ring{pq, el})
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Printf("   electron ring that would give the measured moment (beside a proton ring with +%.4f μN): radius %.4f fm\n", MuP, reN*1e15)
	fmt.Printf("     kronspice: currents %+.4g A, %+.4g A; moment %+.4f μN; magnetic energy %.4f MeV (analytic %.4f, residual %.1e)\n",
		run2.Currents[0], run2.Currents[1], run2.Moment/MuN0, run2.Stored/MeV, run2.Analytic/MeV, run2.Residual)
	eRing := Hbar * C0 / reN / MeV
	fmt.Printf("     but light on a ring that small carries %.1f MeV, against %.3f MeV available (m_n − m_p): %.0f× too much\n", eRing, MnMeV-MpMeV, eRing/(MnMeV-MpMeV))
	fmt.Printf("     Coulomb attraction of the two rings: %.3f MeV (and their own Coulomb energies %+.3f MeV): far too weak to bind it\n",
		(Coulomb([]Ring{pq, el})-Coulomb([]Ring{pq})-Coulomb([]Ring{el}))/MeV, (Coulomb([]Ring{pq})+Coulomb([]Ring{el}))/MeV)
	fmt.Println("   the antineutrino ring (H28, an E ring) is neutral: no current, no moment, energy < 0.45 eV — it changes none of this")

	// Part 3
	mq := MpMeV / MuP
	R := LightRadius(mq)
	pc := QuarkNucleon(2.0/3, -1.0/3, R, R)
	nc := QuarkNucleon(-1.0/3, 2.0/3, R, R)
	muP, emP, r1, err := Weighted(pc)
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	muNq, emN, r2, err := Weighted(nc)
	if err != nil {
		fmt.Println("error:", err)
		return
	}
	fmt.Println("\n3. Standard counterpart as circuits: three quark light rings (uud, udd), spin-weighted (SU(6))")
	fmt.Printf("   one free number: the quark ring's energy, fixed by the proton's moment: m_q c² = %.2f MeV, radius %.4f fm\n", mq, R*1e15)
	fmt.Printf("   proton moment  %+.4f μN (fitted)       measured %+.4f\n", muP, MuP)
	fmt.Printf("   neutron moment %+.4f μN (predicted)    measured %+.4f  (%.1f %%)\n", muNq, MuN, 100*(muNq/MuN-1))
	fmt.Printf("   ratio μn/μp    %+.4f (= −2/3)          measured %+.4f\n", muNq/muP, MuN/MuP)
	fmt.Printf("   three rings' light energy %.1f MeV against %.1f (p) / %.1f (n): %.1f %% high; the binding is not in the model\n", 3*mq, MpMeV, MnMeV, 100*(3*mq/MpMeV-1))
	fmt.Printf("   electromagnetic energy (rings spaced %.3f fm, kronspice magnetic + Coulomb): p %.4f MeV, n %.4f MeV\n", R*1e15, emP, emN)
	fmt.Printf("     → electromagnetism alone makes the neutron %.3f MeV lighter than the proton (lattice QED part, Borsanyi et al. 2015: %.2f MeV)\n", emP-emN, -1.00)
	fmt.Print("     the spacing is a declared choice, not fitted; m_n − m_p from electromagnetism at other spacings:")
	for _, f := range []float64{0.5, 2, 4} {
		_, ep, _, e1 := Weighted(QuarkNucleon(2.0/3, -1.0/3, R, f*R))
		_, en, _, e2 := Weighted(QuarkNucleon(-1.0/3, 2.0/3, R, f*R))
		if e1 != nil || e2 != nil {
			fmt.Println(" error")
			return
		}
		fmt.Printf("  %.1fR: %+.3f MeV;", f, en-ep)
	}
	fmt.Println()
	fmt.Printf("     the measured +1.293 MeV needs the d ring heavier than the u ring (lattice: +2.52 MeV from m_d > m_u); ledger residuals %.1e, %.1e\n", r1, r2)

	// Part 4
	fmt.Println("\n4. The knot picture: Skyrme's model, proton and neutron = one knot in two spin states")
	fmt.Println("   (Adkins, Nappi & Witten, Nucl. Phys. B 228, 552, 1983; quoted, two constants fitted to the N and Δ masses)")
	fmt.Println("   quantity                   knot (ANW)   measured")
	fmt.Println("   proton moment (μN)          +1.87        +2.79")
	fmt.Println("   neutron moment (μN)         −1.31        −1.91")
	fmt.Println("   |μp/μn|                      1.43         1.46")
	fmt.Println("   isoscalar radius (fm)        0.59         0.72")
	fmt.Println("   g_A                          0.61         1.23")
	fmt.Println("   the knot gets the neutron's negative moment and the ratio right with no charged parts at all;")
	fmt.Println("   the sizes run 30 % low. Our ether knot (H29) is Faddeev's two-component cousin; the nucleon is Skyrme's")

	// Part 5
	fmt.Println("\n5. Neutron decay, n → p + e + ν̄")
	fmt.Printf("   energy budget: %.4f − %.4f − %.4f = %.4f MeV shared by the electron and the antineutrino\n", MnMeV, MpMeV, MeMeV, MnMeV-MpMeV-MeMeV)
	fmt.Printf("   lifetime: light rings carry no weak coupling, so the light model gives none; standard model\n")
	fmt.Printf("   τ = 4908.6 s/(V_ud²(1 + 3λ²)) = %.1f s (V_ud %.5f, λ %.4f) against measured %.1f s\n", NeutronLifetime(Vud, Lambda), Vud, Lambda, TauN)
}
