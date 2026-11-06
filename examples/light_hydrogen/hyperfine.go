// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// The 21 cm line as a circuit. The electron's spin, spread over the 1s
// cloud, is a magnetization M(r) = μ_e|ψ(r)|² ẑ with |ψ|² = e^(−2r/a)/(πa³).
// A magnetization is a current, J = ∇×M = (2/a)M(r) sinθ φ̂ here, so the
// cloud is a set of coaxial current loops round the proton. Each loop is a
// kronspice inductor driven to its current; the proton is a small ring at
// the centre carrying its moment; their mutual inductances (Maxwell's
// coaxial-circle formula) couple them. The ledger's stored energy minus
// every ring's own ½LI² is the interaction ΣM·I_p·I_k = μ_p·B(0), and
// B(0) = (2/3)μ0 M(0), Fermi's contact field (Fermi 1930; Jackson §5.7).
//
// Quadrature. In the (r, x = cosθ) half-plane a loop's current is
// dI = (2/a)M(r)·r dr dx. With t = 2r/a = −ln(1 − s), M dr ∝ ds, so loops
// at equal steps in s carry the exact radial weights; in x the integrand
// of B(0) is (1 − x²), exact with a 4-point Gauss–Legendre rule.

var glX = []float64{-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526}
var glW = []float64{0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538}

// Loop is one coaxial current loop: radius Rho, height Z, current I.
type Loop struct{ Rho, Z, I float64 }

// SelfL is a thin ring's self-inductance with core radius ρ/16.
func SelfL(rho float64) float64 { return Mu0 * rho * (math.Log(16) - 1.75) }

// SpinLoops is the electron's spin cloud of moment mu [J/T] in a 1s state
// of radius a, as nr × 4 loops.
func SpinLoops(mu, a float64, nr int) []Loop {
	var ls []Loop
	for k := 0; k < nr; k++ {
		s := (float64(k) + 0.5) / float64(nr)
		r := -a / 2 * math.Log(1-s)
		for j, x := range glX {
			ls = append(ls, Loop{Rho: r * math.Sqrt(1-x*x), Z: r * x, I: mu * r * glW[j] / (float64(nr) * math.Pi * a * a * a)})
		}
	}
	return ls
}

// HFResult is the hyperfine circuit's answer.
type HFResult struct {
	Loops       int
	Interaction float64 // ΣM·I_p·I_k from the solved currents [J]
	Analytic    float64 // μ_p·(2/3)μ0 μ_e|ψ(0)|² [J]
	Stored      float64 // ledger stored energy over ½ΣLI² + ΣMII, minus 1
	Residual    float64 // ledger
	Nu          float64 // 4·Interaction/h [Hz]
}

// HyperfineBuild makes the circuit: the proton ring (light-ring radius,
// current set to the model's moment) and the electron's loops, each fed by
// a current source ramped over 1 ns, coupled to the proton ring.
func HyperfineBuild(a Atom, nr int) (*circuit.Circuit, Loop, []Loop, error) {
	rp := ProtonRingRadius()
	p := Loop{Rho: rp, I: a.P.Mu * MuN0 / (math.Pi * rp * rp)}
	loops := SpinLoops(a.E.Moment(), a.BohrRadius(), nr)
	b := circuit.NewBuilder("hydrogen hyperfine loops")
	add := func(name string, l Loop) {
		n := "n" + name
		b.Inductor("L"+name, n, "0", SelfL(l.Rho))
		b.Add(circuit.Element{Name: "I" + name, Kind: circuit.CurrentSource,
			Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-9, l.I, 2e-9, l.I}}}, "0", n)
	}
	add("P", p)
	for k, l := range loops {
		name := fmt.Sprintf("E%d", k)
		add(name, l)
		m := em.MutualCircles(p.Rho, l.Rho, math.Abs(l.Z))
		b.Mutual("K"+name, "LP", "L"+name, m/math.Sqrt(SelfL(p.Rho)*SelfL(l.Rho)))
	}
	ckt, err := b.BuildTransient()
	return ckt, p, loops, err
}

// HyperfineCircuit solves the circuit and returns the interaction energy
// and the 21 cm frequency (spin ½: the splitting is 4 times the zz energy).
func HyperfineCircuit(a Atom, nr int) (HFResult, error) {
	ckt, p, loops, err := HyperfineBuild(a, nr)
	if err != nil {
		return HFResult{}, err
	}
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: 1e-11, Stop: 2e-9}, analysis.TranOptions{})
	if err != nil {
		return HFResult{}, err
	}
	cur := func(name string) (float64, error) {
		i, ok := w.ElementCurrent("L" + name)
		if !ok {
			return 0, fmt.Errorf("no current for %s", name)
		}
		return i[len(i)-1], nil
	}
	ip, err := cur("P")
	if err != nil {
		return HFResult{}, err
	}
	// The proton ring's own ½LI² is about 10¹³ times the interaction, so
	// the interaction is summed from the solved currents and the mutual
	// inductances, not taken as a difference of stored energies; the
	// ledger's total is checked against the same sum.
	self := 0.5 * SelfL(p.Rho) * ip * ip
	var inter float64
	for k, l := range loops {
		ik, err := cur(fmt.Sprintf("E%d", k))
		if err != nil {
			return HFResult{}, err
		}
		self += 0.5 * SelfL(l.Rho) * ik * ik
		inter += em.MutualCircles(p.Rho, l.Rho, math.Abs(l.Z)) * ip * ik
	}
	r := a.BohrRadius()
	res := HFResult{
		Loops:       len(loops),
		Interaction: inter,
		Analytic:    a.P.Mu * MuN0 * 2.0 / 3 * Mu0 * a.E.Moment() / (math.Pi * r * r * r),
		Stored:      w.Ledger.StoredEnd/(self+inter) - 1,
		Residual:    w.Ledger.MaxRelativeResidual,
	}
	res.Nu = 4 * res.Interaction / H
	return res, nil
}

// KnotScale is the ether knot of H29 sized to hold the electron.
type KnotScale struct {
	L, K2, K4 float64 // ℓ [m], κ2 [J/m], κ4 [J·m]
	ExciteMeV float64 // ħω of the squashing tone
}

// Knot fixes the knot by the electron: core radius 1.2ℓ = the ring's
// radius, energy 547.5√(κ2κ4) = m_e c² with κ4 = κ2ℓ² (ether_knot:
// E = 547.5√(κ2κ4), core ≈ 1.2ℓ, squashing ω = 1.416 c/ℓ, Q 7).
func Knot(re float64) KnotScale {
	l := re / 1.2
	k2 := MeMeV * MeV / (547.5 * l)
	return KnotScale{L: l, K2: k2, K4: k2 * l * l, ExciteMeV: Hbar * 1.416 * C0 / l / MeV}
}
