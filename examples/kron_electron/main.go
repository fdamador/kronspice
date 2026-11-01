// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kron_electron puts hypothesis H21 (CLAUDE.md §1) on the electron
// test bench of doc/kelvin-references.md: a closed ring of circulating
// charge written as one of Kron's machines (G. Kron, Equivalent Circuits of
// Electric Machinery, Wiley 1951, Epilogue, reference/), with the ring's
// self-inductance and the circulating mass together as the metric tensor
// (Epilogue eq. 3–4, p. 254: "the metric tensor representing all the self-
// and mutual inductances and the rotor inertia"), an external field as the
// angle-dependent mutual inductance of frame I, and Kron's eq. 11 (p. 259)
// as the formal link: with e = R = 0 the machine equation is the equation
// of motion of a charged particle in an electromagnetic field.
//
// For a ring of charge q and mass m turning at θ̇, the current is
// i = q·θ̇/2π, so the electrical and the mechanical coordinate are one
// axis and the kinetic energy ½·m·a²·θ̇² = ½·(m·a²·(2π/q)²)·i² is a magnetic
// energy on that axis: the inertia enters the metric as an inductance
// m·a²·(2π/q)², added to the ring's own self-inductance L. All of Kron's
// machine equations then apply unchanged (N element, device/machine.go).
//
//  1. Twin test (standard). Kelvin's 1899 charged ring (tube +e, mass w′;
//     fluid −e, mass w; coil switched on below it), already solved charge
//     by charge with the Lorentz force in examples/kelvin_ring, is written
//     as a three-winding machine (coil, tube, fluid) with constant metric
//     and solved by the circuit engine. Kron's frame must give Kelvin's
//     moments of momentum ∓eΦ/2π, the fluid at Larmor's rate eB̄/2w and a
//     closed ledger, to the bench's own tolerance.
//
//  2. Torque from the metric (standard). The ring tilts about a diameter
//     in a steady field B: frame I with M(θ) = B·π·a²·cos θ between a
//     field winding and the ring. Kron's torque ½·iᵀ·(∂L/∂θ)·i is
//     −μ·B·sin θ (μ = i·π·a²), and the period of the small oscillation is
//     checked against 2π·√(J/μB).
//
//  3. The electron as the ring (hypothesis H21). With the charge e fixed
//     and the total angular momentum of the metric set to ħ/2, the model
//     gives the moment, g-factor, speed and energies of rings of several
//     radii, with the electron's mass circulating or with inertia from the
//     field alone, and solves for the massless ring that has the measured
//     moment. Each is printed beside the Dirac/QED value and the classical
//     obstacle it meets (g = 1, v > c, size bound).
//
//  4. The field of the ring (standard), drawn with -png: field lines of the
//     ring current from Biot–Savart (em.Loop), in the meridional plane and
//     in an oblique view, for comparison with pictures of toroidal
//     electron models.
//
//     go run ./examples/kron_electron [-png viz/ring_field.png]
package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/gif"
	"image/png"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// Measured constants, SI (CODATA 2018 recommended values, NIST). The size
// bound is the LEP contact-interaction limit on electron substructure,
// below 10⁻¹⁸ m (D. Bourilkov, Phys. Rev. D 64, 071701, 2001).
const (
	QE         = 1.602176634e-19  // elementary charge [C]
	ME         = 9.1093837015e-31 // electron mass [kg]
	HBAR       = 1.054571817e-34  // reduced Planck constant [J s]
	C0         = 299792458.0      // speed of light [m/s]
	MuB        = 9.2740100783e-24 // Bohr magneton [J/T]
	GMeasured  = 2.00231930436256 // electron g-factor (magnitude)
	MuMeasured = 1.00115965218128 * MuB
	LambdaC    = 3.8615926796e-13 // reduced Compton wavelength ħ/(m_e c) [m]
	RClassical = 2.8179403262e-15 // classical electron radius [m]
	SizeBound  = 1e-18            // [m]
)

// FineStructure is α = μ₀ c e²/(4πħ) from the constants above.
var FineStructure = em.Mu0 * C0 * QE * QE / (4 * math.Pi * HBAR)

// GammaElectronHz is the electron's measured spin resonance at 1 T,
// |γ_e|/2π = 28.024 951 GHz/T (CODATA 2018).
const GammaElectronHz = 28.0249514242e9

// ThinRingInductance is the self-inductance of a circular ring of radius R
// made of round wire of radius r (R ≫ r), L = μ₀R(ln(8R/r) − 7/4)
// (Maxwell, Treatise Art. 706; uniform current over the wire section).
func ThinRingInductance(R, r float64) float64 { return em.Mu0 * R * (math.Log(8*R/r) - 1.75) }

// ThinRingCapacitance is the capacitance of a thin conducting torus of
// ring radius R and tube radius r to a distant return, 4π²ε₀R/ln(8R/r)
// (leading term; W. R. Smythe, Static and Dynamic Electricity, §5.17).
func ThinRingCapacitance(R, r float64) float64 {
	return 4 * math.Pi * math.Pi * em.Eps0 * R / math.Log(8*R/r)
}

// ---------- 1. Kelvin's ring as a Kron machine ----------

// Ring is Kelvin's ring: radius A, tube charge +E and mass MTube, fluid
// charge −E and mass MFluid (as in examples/kelvin_ring).
type Ring struct{ A, E, MTube, MFluid float64 }

// Coil is a single-turn-equivalent coil of radius R carrying NI
// ampere-turns, coaxial with the ring, Z below it.
type Coil struct{ R, Z, NI float64 }

// InertialInductance is the inductance that a mass m circulating at
// radius a with charge q presents on its current axis, m·a²·(2π/q)².
func InertialInductance(m, a, q float64) float64 {
	return m * a * a * (2 * math.Pi / q) * (2 * math.Pi / q)
}

// KelvinSpec is Kron's metric tensor over the three current axes (coil,
// tube, fluid): the coil's self-inductance, Maxwell's mutual inductance of
// coaxial circles between coil and ring (em.MutualCircles), and on the
// charge axes the inertia of tube and fluid as inductances. With
// selfInduction the ring's own inductance is added (self and tube–fluid
// mutual, the two currents sharing one ring); Kelvin's 1899 analysis
// leaves it out, and so does the Lorentz-force bench.
func KelvinSpec(rg Ring, c Coil, selfInduction bool) (*circuit.MachineSpec, float64) {
	s := circuit.NewMachineSpec(3)
	m := em.MutualCircles(c.R, rg.A, c.Z)
	s.L0[0][0] = ThinRingInductance(c.R, c.R/100) // declared wire radius R/100
	s.L0[0][1], s.L0[1][0] = m, m
	s.L0[0][2], s.L0[2][0] = m, m
	s.L0[1][1] = InertialInductance(rg.MTube, rg.A, rg.E)
	s.L0[2][2] = InertialInductance(rg.MFluid, rg.A, rg.E)
	if selfInduction {
		l := ThinRingInductance(rg.A, rg.A/100)
		s.L0[1][1] += l
		s.L0[2][2] += l
		s.L0[1][2], s.L0[2][1] = l, l
	}
	return s, m
}

// KelvinCircuit builds the circuit: the coil winding driven by a current
// source ramping linearly to NI over T and holding; tube and fluid windings
// closed through 0 V sources (frictionless: no resistance), so their flux
// linkages are conserved; the mechanical port idle (the metric does not
// depend on an angle here).
func KelvinCircuit(rg Ring, c Coil, T float64, selfInduction bool) (*circuit.Circuit, error) {
	s, _ := KelvinSpec(rg, c, selfInduction)
	ramp := circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, T, c.NI, 2 * T, c.NI}}
	return circuit.NewBuilder("Kelvin's charged ring as a Kron machine").
		Domain(exchange.Mechanical, "w").
		MachineModel("ring", s).
		CurrentWave("IC", "0", "c", ramp).
		Resistor("RC", "c", "0", 1e9). // fixes the source node's potential at UIC (10⁻¹¹ A against 20 kA)
		Machine("N1", "ring", [][2]string{{"c", "0"}, {"t", "0"}, {"f", "0"}}, "w", "0").
		VoltageSource("VT", "t", "0", 0).
		VoltageSource("VF", "f", "0", 0).
		Resistor("RW", "w", "0", 1).
		BuildTransient()
}

// KelvinResult is the outcome of the switch-on, read from the machine's
// currents.
type KelvinResult struct {
	Flux               float64 // M·NI through the ring [Wb]
	LTube, LFluid      float64 // moments of momentum J·θ̇ [kg m²/s]
	LKelvin            float64 // eΦ/2π
	OmegaFluid, Larmor float64 // [rad/s]
	KE                 float64 // ½Jθ̇² of tube and fluid [J]
	WRing              float64 // source work against the ring's back-EMF, ∫M·I_c·d(i_t + i_f) [J]
	WInteraction       float64 // M·I_c·(i_t + i_f) at the end [J]
	RingResidual       float64 // (WRing − KE − WInteraction)/KE, the bench's ledger
	Stored             float64 // ½iᵀLi of the machine at the end, mostly the coil's own field [J]
	EnergyIn           float64 // work by the coil's source [J]
	Mismatch           float64 // engine ledger integration mismatch [J]
	Points             int
	// Time series for the animation: coil current and the angles turned
	// by tube and fluid (∫θ̇ dt, trapezoidal).
	Time, CoilCurrent, ThetaTube, ThetaFluid []float64
}

func last(x []float64) float64 { return x[len(x)-1] }

// integrate returns ∫y dt over the points by the trapezoidal rule,
// starting at y0.
func integrate(t, y []float64, y0 float64) []float64 {
	out := make([]float64, len(t))
	out[0] = y0
	for k := 1; k < len(t); k++ {
		out[k] = out[k-1] + (y[k]+y[k-1])/2*(t[k]-t[k-1])
	}
	return out
}

// RunKelvin solves the circuit and reads Kelvin's quantities from it.
func RunKelvin(rg Ring, c Coil, T float64, selfInduction bool) (KelvinResult, error) {
	cir, err := KelvinCircuit(rg, c, T, selfInduction)
	if err != nil {
		return KelvinResult{}, err
	}
	wf, err := analysis.Tran(cir, analysis.TranSpec{Step: T / 200, Stop: 2 * T, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-6}})
	if err != nil {
		return KelvinResult{}, err
	}
	return kelvinFromWaveforms(wf, rg, c), nil
}

func kelvinFromWaveforms(wf *output.Waveforms, rg Ring, c Coil) KelvinResult {
	_, m := KelvinSpec(rg, c, false)
	// A 0 V source between the winding's + node and ground carries the
	// winding current back, so I(V) = −i.
	vt, _ := wf.ElementCurrent("VT")
	vf, _ := wf.ElementCurrent("VF")
	iT, iF := -last(vt), -last(vf)
	thT := 2 * math.Pi / rg.E * iT  // tube: +e
	thF := -2 * math.Pi / rg.E * iF // fluid: −e
	jt, jf := rg.MTube*rg.A*rg.A, rg.MFluid*rg.A*rg.A
	r := KelvinResult{Flux: m * c.NI, Points: len(wf.Time)}
	r.LTube, r.LFluid = jt*thT, jf*thF
	r.LKelvin = rg.E * r.Flux / (2 * math.Pi)
	r.OmegaFluid = thF
	r.Larmor = rg.E * (r.Flux / (math.Pi * rg.A * rg.A)) / (2 * rg.MFluid)
	r.KE = 0.5*jt*thT*thT + 0.5*jf*thF*thF
	// The bench's ledger: the coil source's work against the ring's
	// back-EMF, M·I_c·d(i_t + i_f), by the trapezoidal rule over the points.
	ic, _ := wf.ElementCurrent("IC")
	for k := 1; k < len(wf.Time); k++ {
		dRing := -(vt[k] + vf[k]) + (vt[k-1] + vf[k-1])
		r.WRing += m * (ic[k] + ic[k-1]) / 2 * dRing
	}
	r.WInteraction = m * last(ic) * (iT + iF)
	r.RingResidual = (r.WRing - r.KE - r.WInteraction) / r.KE
	r.Time, r.CoilCurrent = wf.Time, ic
	omT := make([]float64, len(vt))
	omF := make([]float64, len(vf))
	for k := range vt {
		omT[k] = -2 * math.Pi / rg.E * vt[k]
		omF[k] = 2 * math.Pi / rg.E * vf[k]
	}
	r.ThetaTube, r.ThetaFluid = integrate(wf.Time, omT, 0), integrate(wf.Time, omF, 0)
	r.Stored = wf.Ledger.StoredEnd
	r.EnergyIn = wf.Ledger.EnergyIn - wf.Ledger.EnergyOut
	r.Mismatch = wf.Ledger.IntegrationMismatch
	return r
}

// ---------- 2. Torque from the metric: the ring tilting in a field ----------

// Tilt is a ring of charge E and mass Mass at radius A, spinning at
// OmegaSpin (its current E·ω/2π), free to tilt about a diameter in a steady
// field B, released from rest at angle Theta0.
type Tilt struct{ A, E, Mass, OmegaSpin, B, Theta0 float64 }

// Current, Moment and Inertia of the tilting ring.
func (t Tilt) Current() float64 { return t.E * t.OmegaSpin / (2 * math.Pi) }
func (t Tilt) Moment() float64  { return t.Current() * math.Pi * t.A * t.A }
func (t Tilt) Inertia() float64 { return 0.5 * t.Mass * t.A * t.A } // ring about a diameter

// TiltSpec is frame I: a field winding (1 A standing for B) and the ring
// winding, coupled by M(θ) = B·π·a²·cos θ; the spin inertia on the ring's
// current axis keeps its current nearly constant (flux conservation).
func TiltSpec(t Tilt) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.Theta0 = t.Theta0
	s.L0[0][0] = 1 // field winding self-inductance [H], declared
	s.L0[1][1] = InertialInductance(t.Mass, t.A, t.E) + ThinRingInductance(t.A, t.A/100)
	phi := t.B * math.Pi * t.A * t.A
	h := s.Harmonic(1)
	h.Cos[0][1], h.Cos[1][0] = phi, phi
	s.IC[0], s.IC[1] = 1, t.Current()
	return s
}

// TiltCircuit: field winding on a 1 A source, ring winding closed through
// 0 V, the tilt inertia as a capacitor on the mechanical port.
func TiltCircuit(t Tilt) (*circuit.Circuit, error) {
	return circuit.NewBuilder("ring tilting in a field: Kron frame I").
		Domain(exchange.Mechanical, "w").
		MachineModel("ring", TiltSpec(t)).
		CurrentSource("IB", "0", "b", 1).
		Resistor("RB", "b", "0", 1e9).
		Machine("N1", "ring", [][2]string{{"b", "0"}, {"r", "0"}}, "w", "0").
		VoltageSource("VR", "r", "0", 0).
		CapacitorIC("CJ", "w", "0", t.Inertia(), 0).
		BuildTransient()
}

// TiltResult: the measured period against the pendulum value.
type TiltResult struct {
	Period, Analytic float64 // [s]; analytic includes the amplitude correction
	CurrentDrift     float64 // relative change of the ring current over the run
	Mismatch         float64 // ledger [J]
	Points           int
	// Time series for the animation: tilt angle (∫ω dt from Theta0) and ω.
	Time, Theta, Omega []float64
}

// RunTilt simulates three periods and times the zero crossings of ω.
func RunTilt(t Tilt) (TiltResult, error) {
	cir, err := TiltCircuit(t)
	if err != nil {
		return TiltResult{}, err
	}
	w0 := math.Sqrt(t.Moment() * t.B / t.Inertia())
	t0 := 2 * math.Pi / w0
	// Exact pendulum period 4K(k)/ω₀, k = sin(θ₀/2).
	kk, _ := em.EllipticKE(math.Sin(t.Theta0 / 2))
	an := 4 * kk / w0
	wf, err := analysis.Tran(cir, analysis.TranSpec{Step: an / 400, Stop: 3.2 * an, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-6}})
	if err != nil {
		return TiltResult{}, err
	}
	om, _ := wf.NodeVoltage("w")
	var cross []float64
	for k := 1; k < len(om); k++ {
		if wf.Time[k] > t0/8 && om[k-1] != 0 && (om[k-1] < 0) != (om[k] < 0) {
			f := om[k-1] / (om[k-1] - om[k])
			cross = append(cross, wf.Time[k-1]+f*(wf.Time[k]-wf.Time[k-1]))
		}
	}
	res := TiltResult{Analytic: an, Mismatch: wf.Ledger.IntegrationMismatch, Points: len(wf.Time)}
	if n := len(cross); n >= 2 {
		res.Period = 2 * (cross[n-1] - cross[0]) / float64(n-1)
	}
	ir, _ := wf.ElementCurrent("VR")
	res.CurrentDrift = (last(ir) - ir[0]) / ir[0]
	res.Time, res.Omega = wf.Time, om
	res.Theta = integrate(wf.Time, om, t.Theta0)
	return res, nil
}

// Harmonics of the tilt: the ratios of the 3rd and 5th harmonic
// amplitudes of ω(t) to the fundamental, measured by Fourier sums over
// the last two periods of the run, and the exact pendulum values from the
// Jacobi series θ(t) = Σ_{n odd} (8/n)·qⁿᐟ²/(1 + qⁿ)·cos(nωt), nome
// q = exp(−πK′/K), modulus k = sin(θ₀/2) (so the velocity ratios are
// n·A_n/A_1).
type Harmonics struct {
	Fundamental float64    // [rad/s]
	Measured    [2]float64 // |ω₃|/|ω₁|, |ω₅|/|ω₁|
	Exact       [2]float64
	Even        float64 // |ω₂|/|ω₁|, zero by symmetry
}

// TiltHarmonics runs the tilt and analyses its spectrum.
func TiltHarmonics(t Tilt) (Harmonics, error) {
	r, err := RunTilt(t)
	if err != nil {
		return Harmonics{}, err
	}
	w1 := 2 * math.Pi / r.Period
	// Fourier sums over exactly two periods ending at the last point,
	// Simpson's rule on the interpolated record, so no spectral leakage.
	const n = 4000
	ta := last(r.Time) - 2*r.Period
	amp := func(m int) float64 {
		var re, im float64
		for j := 0; j <= n; j++ {
			tt := ta + 2*r.Period*float64(j)/n
			w := 2.0
			if j == 0 || j == n {
				w = 1
			} else if j%2 == 1 {
				w = 4
			}
			om := sample(r.Time, r.Omega, tt)
			re += w * om * math.Cos(float64(m)*w1*tt)
			im += w * om * math.Sin(float64(m)*w1*tt)
		}
		return math.Hypot(re, im)
	}
	a1 := amp(1)
	h := Harmonics{Fundamental: w1, Measured: [2]float64{amp(3) / a1, amp(5) / a1}, Even: amp(2) / a1}
	h.Exact = PendulumHarmonics(t.Theta0)
	return h, nil
}

// PendulumHarmonics returns the exact velocity-harmonic ratios
// |ω₃/ω₁|, |ω₅/ω₁| of a pendulum released from rest at θ₀, from the Jacobi
// series above; at small θ₀ they tend to θ₀²/64 and θ₀⁴/3072·(5/3)…
func PendulumHarmonics(theta0 float64) [2]float64 {
	k := math.Sin(theta0 / 2)
	kk, _ := em.EllipticKE(k)
	kp, _ := em.EllipticKE(math.Sqrt(1 - k*k))
	q := math.Exp(-math.Pi * kp / kk)
	coef := func(n int) float64 { return math.Pow(q, float64(n)/2) / (1 + math.Pow(q, float64(n))) / float64(n) }
	return [2]float64{3 * coef(3) / coef(1), 5 * coef(5) / coef(1)}
}

// ---------- 3. The electron as the ring (H21) ----------

// RingModel is a ring of charge QE: radius A, wire (core) radius Core, and
// circulating mass Mass (ME, or 0 for inertia from the field alone).
type RingModel struct {
	Name          string
	A, Core, Mass float64
}

// RingNumbers are the model's values at a given total angular momentum.
type RingNumbers struct {
	L             float64 // ring self-inductance [H]
	Omega, V      float64 // angular velocity [rad/s], rim speed [m/s]
	I, Mu         float64 // current [A], magnetic moment [J/T]
	S, G          float64 // angular momentum [J s], g-factor μ = g·(e/2m_e)·S
	KE, WMag, WEs float64 // kinetic, magnetic and electrostatic energy [J]
	Precession    float64 // μB/S per tesla [rad/s/T] (Larmor eB/2m_e for g = 1)
}

// AtSpin imposes the total angular momentum S of the metric,
// (Mass·A² + L·e²/4π²)·ω = S, and computes the rest.
func (r RingModel) AtSpin(S float64) RingNumbers {
	l := ThinRingInductance(r.A, r.Core)
	inertia := r.Mass*r.A*r.A + l*QE*QE/(4*math.Pi*math.Pi)
	w := S / inertia
	i := QE * w / (2 * math.Pi)
	mu := i * math.Pi * r.A * r.A
	return RingNumbers{
		L: l, Omega: w, V: w * r.A, I: i, Mu: mu, S: S,
		G:          mu * 2 * ME / (QE * S),
		KE:         0.5 * r.Mass * r.A * r.A * w * w,
		WMag:       0.5 * l * i * i,
		WEs:        QE * QE / (2 * ThinRingCapacitance(r.A, r.Core)),
		Precession: mu / S,
	}
}

// Frequencies are the natural frequencies of a ring at angular momentum
// S in a field B [Hz]:
//   - Rotation: the charge's turning rate ω/2π, i.e. the current's
//     fundamental;
//   - Precession: the spin's precession in B, (μ/S)·B/2π = g·eB/4πm, which
//     is the magnetic-resonance frequency (examples/kron_nmr);
//   - Libration: the tilt oscillation of part 2, √(μB/J_tilt)/2π with
//     J_tilt = ½·m·a² (rings with mass only; it is what a gyroscope shows
//     as nutation once it spins);
//   - SelfResonance: the ring as an LC circuit, 1/(2π√(LC)) with its
//     thin-ring inductance and capacitance;
//   - Circumferential: a wave once round the ring at c, c/2πa, the
//     fundamental of the ring as a transmission line.
type Frequencies struct {
	Rotation, Precession, Libration, SelfResonance, Circumferential float64
}

// Frequencies evaluates them at angular momentum S in field B.
func (r RingModel) Frequencies(S, B float64) Frequencies {
	n := r.AtSpin(S)
	f := Frequencies{
		Rotation:        n.Omega / (2 * math.Pi),
		Precession:      n.Precession * B / (2 * math.Pi),
		SelfResonance:   1 / (2 * math.Pi * math.Sqrt(n.L*ThinRingCapacitance(r.A, r.Core))),
		Circumferential: C0 / (2 * math.Pi * r.A),
	}
	if r.Mass > 0 {
		f.Libration = math.Sqrt(n.Mu*B/(0.5*r.Mass*r.A*r.A)) / (2 * math.Pi)
	}
	return f
}

// Reference frequencies [Hz]: the Compton frequency m_e c²/h, which equals
// c/2πλ̄_C, and the electron's measured spin resonance at 1 T.
const (
	HPlanck   = 6.62607015e-34
	ComptonHz = ME * C0 * C0 / HPlanck
)

// FitMassless finds the ring radius at which a massless ring (inertia
// from its own field) carries the measured moment at S = ħ/2, for a given
// core-to-radius ratio, by bisection on a²/L(a).
func FitMassless(coreRatio float64) RingModel {
	target := MuMeasured
	f := func(a float64) float64 {
		return RingModel{A: a, Core: a * coreRatio}.AtSpin(HBAR/2).Mu - target
	}
	lo, hi := 1e-20, 1e-9
	for k := 0; k < 200; k++ {
		mid := math.Sqrt(lo * hi)
		if f(mid) < 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	a := math.Sqrt(lo * hi)
	return RingModel{Name: "massless ring fitted to the measured moment", A: a, Core: a * coreRatio}
}

// Scenarios are the rings reported in part 3.
func Scenarios() []RingModel {
	const core = 0.01
	return []RingModel{
		{"Compton ring (a = ħ/m_e c), electron mass circulating", LambdaC, LambdaC * core, ME},
		{"half-Compton ring, electron mass circulating (v → c)", LambdaC / 2, LambdaC / 2 * core, ME},
		{"Compton ring, massless (inertia from the field)", LambdaC, LambdaC * core, 0},
		{"classical-radius ring, massless", RClassical, RClassical * core, 0},
		FitMassless(core),
	}
}

// ---------- 4. The field of the ring ----------

// FieldLine traces one field line of a unit ring current of radius a in
// the meridional plane from (rho0, z0), returning (ρ, z) points. It stops
// when the line closes on its start or after maxSteps.
func FieldLine(a, rho0, z0 float64, maxSteps int) [][2]float64 {
	loop := em.Loop{P: em.Circle([3]float64{0, 0, 0}, [3]float64{0, 0, 1}, a, 360), I: 1}
	dir := func(p [2]float64) [2]float64 {
		b := loop.B([3]float64{p[0], 0, p[1]})
		n := math.Hypot(b[0], b[2])
		if n == 0 {
			return [2]float64{0, 1}
		}
		return [2]float64{b[0] / n, b[2] / n}
	}
	h := 0.02 * a
	p := [2]float64{rho0, z0}
	pts := [][2]float64{p}
	for k := 0; k < maxSteps; k++ {
		k1 := dir(p)
		k2 := dir([2]float64{p[0] + h/2*k1[0], p[1] + h/2*k1[1]})
		k3 := dir([2]float64{p[0] + h/2*k2[0], p[1] + h/2*k2[1]})
		k4 := dir([2]float64{p[0] + h*k3[0], p[1] + h*k3[1]})
		p = [2]float64{p[0] + h/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]), p[1] + h/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])}
		pts = append(pts, p)
		if k > 50 && math.Hypot(p[0]-rho0, p[1]-z0) < 1.5*h {
			pts = append(pts, [2]float64{rho0, z0})
			break
		}
	}
	return pts
}

// canvas is a minimal line drawer on an RGBA image.
type canvas struct {
	img  *image.RGBA
	w, h int
}

func newCanvas(w, h int) *canvas {
	c := &canvas{img: image.NewRGBA(image.Rect(0, 0, w, h)), w: w, h: h}
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			c.img.Set(x, y, color.White)
		}
	}
	return c
}

// line draws from (x0, y0) to (x1, y1), keeping to the clip rectangle
// [xa, xb) × [0, h).
func (c *canvas) line(x0, y0, x1, y1 float64, col color.Color, xa, xb float64) {
	n := int(math.Max(math.Abs(x1-x0), math.Abs(y1-y0))) + 1
	for k := 0; k <= n; k++ {
		f := float64(k) / float64(n)
		x, y := x0+f*(x1-x0), y0+f*(y1-y0)
		if x >= xa && x < xb && y >= 0 && int(y) < c.h {
			c.img.Set(int(x), int(y), col)
		}
	}
}

func (c *canvas) dot(x, y, r float64, col color.Color) {
	for dy := -r; dy <= r; dy++ {
		for dx := -r; dx <= r; dx++ {
			if dx*dx+dy*dy <= r*r {
				c.img.Set(int(x+dx), int(y+dy), col)
			}
		}
	}
}

// DrawRingField writes a two-panel picture: left, field lines of the ring
// current in the meridional plane (the ring's wire as two dots, the axis
// vertical); right, the same lines at twelve azimuths in an oblique view.
func DrawRingField(path string) error {
	const W, H = 1400, 700
	cv := newCanvas(W, H)
	blue := color.RGBA{30, 60, 200, 255}
	grey := color.RGBA{150, 150, 150, 255}
	red := color.RGBA{200, 40, 40, 255}
	a := 1.0
	var lines [][][2]float64
	for _, rho := range []float64{0.3, 0.45, 0.6, 0.72, 0.82, 0.9} {
		lines = append(lines, FieldLine(a, rho*a, 0, 40000))
	}
	// Left panel: meridional plane, 3a half-width, clipped to its half.
	const half = 0.5 * W
	sc := 0.5 * H / 3.0
	cx, cy := 0.25*W, 0.5*H
	cv.line(cx, 20, cx, H-20, grey, 0, half)
	cv.line(cx-2.9*sc, cy, cx+2.9*sc, cy, grey, 0, half)
	for _, ln := range lines {
		for k := 1; k < len(ln); k++ {
			for _, s := range []float64{1, -1} { // both sides of the axis
				cv.line(cx+s*ln[k-1][0]*sc, cy-ln[k-1][1]*sc, cx+s*ln[k][0]*sc, cy-ln[k][1]*sc, blue, 0, half)
			}
		}
	}
	cv.dot(cx+a*sc, cy, 5, red)
	cv.dot(cx-a*sc, cy, 5, red)
	// Right panel: oblique projection, ring axis vertical, the three
	// middle lines at eight azimuths.
	cx2 := 0.75 * W
	proj := func(rho, phi, z float64) (float64, float64) {
		x, y := rho*math.Cos(phi), rho*math.Sin(phi)
		return cx2 + (x-0.45*y)*sc, cy - (z-0.35*y)*sc
	}
	cv.line(cx2, 20, cx2, H-20, grey, half, W)
	for k := 0; k < 8; k++ {
		phi := 2 * math.Pi * float64(k) / 8
		for _, ln := range lines[0:4] {
			for j := 1; j < len(ln); j++ {
				x0, y0 := proj(ln[j-1][0], phi, ln[j-1][1])
				x1, y1 := proj(ln[j][0], phi, ln[j][1])
				cv.line(x0, y0, x1, y1, blue, half, W)
			}
		}
	}
	for k := 0; k < 360; k++ {
		x0, y0 := proj(a, 2*math.Pi*float64(k)/360, 0)
		x1, y1 := proj(a, 2*math.Pi*float64(k+1)/360, 0)
		cv.line(x0, y0, x1, y1, red, half, W)
		cv.line(x0, y0+1, x1, y1+1, red, half, W)
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	return png.Encode(f, cv.img)
}

// ---------- 5. Animation ----------

// sample returns y at time tt by linear interpolation on the points.
func sample(t, y []float64, tt float64) float64 {
	k := sort.SearchFloat64s(t, tt)
	if k <= 0 {
		return y[0]
	}
	if k >= len(t) {
		return last(y)
	}
	f := (tt - t[k-1]) / (t[k] - t[k-1])
	return y[k-1] + f*(y[k]-y[k-1])
}

func (c *canvas) arrow(x0, y0, x1, y1 float64, col color.Color, xa, xb float64) {
	c.line(x0, y0, x1, y1, col, xa, xb)
	ang := math.Atan2(y1-y0, x1-x0)
	for _, d := range []float64{2.6, -2.6} {
		c.line(x1, y1, x1-9*math.Cos(ang+d), y1-9*math.Sin(ang+d), col, xa, xb)
	}
}

func (c *canvas) circle(cx, cy, r float64, col color.Color, xa, xb float64) {
	for k := 0; k < 360; k++ {
		a0, a1 := 2*math.Pi*float64(k)/360, 2*math.Pi*float64(k+1)/360
		c.line(cx+r*math.Cos(a0), cy+r*math.Sin(a0), cx+r*math.Cos(a1), cy+r*math.Sin(a1), col, xa, xb)
	}
}

// DrawMotionGIF writes an animation of the two transients, three panels
// progressing together over their own runs: Kelvin's ring from above as
// the coil switches on (tube red and fluid blue, turning opposite ways,
// angles scaled so the fluid ends at half a turn; the coil current and
// the two moments of momentum as bars), the ring tilting in the field
// seen edge-on (angle scaled by tiltScale; moment arrow red, field grey)
// and its dipole field lines tilting with it. Each panel has a progress
// bar along its bottom edge.
func DrawMotionGIF(path string, k KelvinResult, t TiltResult, frames int, tiltScale float64) error {
	const W, H = 1200, 400
	const pw = W / 3.0
	blue := color.RGBA{30, 60, 200, 255}
	red := color.RGBA{200, 40, 40, 255}
	grey := color.RGBA{150, 150, 150, 255}
	black := color.RGBA{0, 0, 0, 255}
	pal := color.Palette{color.White, blue, red, grey, black, color.RGBA{60, 160, 60, 255}}
	var lines [][][2]float64
	for _, rho := range []float64{0.3, 0.45, 0.6, 0.75, 0.9} {
		lines = append(lines, FieldLine(1, rho, 0, 40000))
	}
	angScale := math.Pi / math.Abs(last(k.ThetaFluid))
	lMax := 0.0
	for i := range k.Time {
		lMax = math.Max(lMax, math.Abs(k.ThetaFluid[i]))
	}
	g := &gif.GIF{LoopCount: 0}
	for f := 0; f < frames; f++ {
		prog := float64(f) / float64(frames-1)
		cv := newCanvas(W, H)
		// Panel 1: Kelvin's ring from above.
		tk := prog * last(k.Time)
		cx, cy, r := pw/2, 0.47*H, 0.33*H
		cv.circle(cx, cy, r, grey, 0, pw)
		thT := angScale * sample(k.Time, k.ThetaTube, tk)
		thF := angScale * sample(k.Time, k.ThetaFluid, tk)
		cv.dot(cx+r*math.Cos(thT), cy-r*math.Sin(thT), 7, red)
		cv.dot(cx+0.8*r*math.Cos(thF), cy-0.8*r*math.Sin(thF), 7, blue)
		cv.line(cx, cy, cx+r*math.Cos(thT), cy-r*math.Sin(thT), red, 0, pw)
		cv.line(cx, cy, cx+0.8*r*math.Cos(thF), cy-0.8*r*math.Sin(thF), blue, 0, pw)
		// bars: coil current (green), moment of momentum of tube (red,
		// leftwards) and fluid (blue, rightwards), from the centre line.
		ic := sample(k.Time, k.CoilCurrent, tk) / last(k.CoilCurrent)
		for y := 0.86 * H; y < 0.89*H; y++ {
			cv.line(20, y, 20+ic*(pw-40), y, pal[5], 0, pw)
		}
		lt := math.Abs(sample(k.Time, k.ThetaTube, tk)) / lMax * (pw/2 - 20)
		lf := math.Abs(sample(k.Time, k.ThetaFluid, tk)) / lMax * (pw/2 - 20)
		for y := 0.91 * H; y < 0.94*H; y++ {
			cv.line(pw/2, y, pw/2-lt, y, red, 0, pw)
			cv.line(pw/2, y, pw/2+lf, y, blue, 0, pw)
		}
		// Panel 2: the tilting ring edge-on.
		tt := prog * last(t.Time)
		th := tiltScale * sample(t.Time, t.Theta, tt)
		cx2 := 1.5 * pw
		for x := cx2 - 0.4*pw; x <= cx2+0.4*pw; x += 0.1 * pw {
			cv.arrow(x, 0.85*H, x, 0.12*H, grey, pw, 2*pw)
		}
		ex, ey := r*math.Cos(th), -r*math.Sin(th)
		for d := -2.0; d <= 2; d++ {
			cv.line(cx2-ex+d*math.Sin(th), cy-ey+d*math.Cos(th), cx2+ex+d*math.Sin(th), cy+ey+d*math.Cos(th), black, pw, 2*pw)
		}
		cv.dot(cx2-ex, cy-ey, 5, red)
		cv.dot(cx2+ex, cy+ey, 5, red)
		cv.arrow(cx2, cy, cx2-0.6*r*math.Sin(th), cy-0.6*r*math.Cos(th), red, pw, 2*pw)
		// Panel 3: the field lines, tilted with the ring.
		cx3 := 2.5 * pw
		sc := 0.5 * H / 3.0
		cs, sn := math.Cos(th), math.Sin(th)
		rot := func(x, z float64) (float64, float64) { return cx3 + (x*cs-z*sn)*sc, cy - (x*sn+z*cs)*sc }
		for _, ln := range lines {
			for j := 1; j < len(ln); j++ {
				for _, s := range []float64{1, -1} {
					x0, y0 := rot(s*ln[j-1][0], ln[j-1][1])
					x1, y1 := rot(s*ln[j][0], ln[j][1])
					cv.line(x0, y0, x1, y1, blue, 2*pw, W)
				}
			}
		}
		x0, y0 := rot(1, 0)
		x1, y1 := rot(-1, 0)
		cv.dot(x0, y0, 5, red)
		cv.dot(x1, y1, 5, red)
		// Progress bars and panel borders.
		for p := 0.0; p < 3; p++ {
			cv.line(p*pw, H-3, p*pw+prog*pw, H-3, black, 0, W)
			cv.line(p*pw, H-2, p*pw+prog*pw, H-2, black, 0, W)
			cv.line(p*pw, 0, p*pw, H, grey, 0, W)
		}
		fr := image.NewPaletted(image.Rect(0, 0, W, H), pal)
		draw.Draw(fr, fr.Bounds(), cv.img, image.Point{}, draw.Src)
		g.Image = append(g.Image, fr)
		g.Delay = append(g.Delay, 6)
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	return gif.EncodeAll(fh, g)
}

func main() {
	pngPath := flag.String("png", "", "write the ring's field-line picture to this PNG")
	gifPath := flag.String("gif", "", "write the animation of the two transients to this GIF")
	flag.Parse()

	fmt.Println("H21: a ring of circulating charge as one of Kron's machines (Equivalent Circuits of Electric Machinery, 1951, Epilogue)")
	fmt.Println()
	fmt.Println("1. Twin test: Kelvin's 1899 ring (a = 5 cm, ±1 µC, tube and fluid 10 g; coil r = 10 cm, 5 cm below, 20 kA-turns, 1 s ramp)")
	rg := Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.01}
	c := Coil{R: 0.1, Z: 0.05, NI: 2e4}
	s, m := KelvinSpec(rg, c, false)
	fmt.Printf("   metric tensor [H]: coil %.4g, coil–ring mutual %.4g, tube %.4g, fluid %.4g (inertia as inductance; ring's own L %.3g, dropped as in Kelvin)\n",
		s.L0[0][0], m, s.L0[1][1], s.L0[2][2], ThinRingInductance(rg.A, rg.A/100))
	var kelvin KelvinResult
	for _, si := range []bool{false, true} {
		r, err := RunKelvin(rg, c, 1, si)
		if err != nil {
			log.Fatal(err)
		}
		if !si {
			kelvin = r
		}
		tag := "Kelvin's metric"
		if si {
			tag = "with the ring's self-inductance"
		}
		fmt.Printf("   %-32s flux %.6g Wb; moment of momentum tube %.6g, fluid %.6g kg m²/s (eΦ/2π = %.6g)\n", tag, r.Flux, r.LTube, r.LFluid, r.LKelvin)
		fmt.Printf("   %-32s fluid %.6g rad/s = Larmor eB̄/2w %.6g; source work against the ring %.6g J = kinetic %.6g J + mutual-field energy %.6g J (residual %.1e)\n",
			"", r.OmegaFluid, r.Larmor, r.WRing, r.KE, r.WInteraction, r.RingResidual)
		fmt.Printf("   %-32s engine ledger: source %.6g J, machine's ½iᵀLi %.6g J (the coil's own field), mismatch %.1e J, %d points\n",
			"", r.EnergyIn, r.Stored, r.Mismatch, r.Points)
	}
	fmt.Println("   examples/kelvin_ring (Lorentz force on 24 charges, Biot–Savart flux 0.000698736 Wb): tube −1.11207e-10, fluid 1.11207e-10 kg m²/s, 4.44829e-06 rad/s, kinetic 4.94683e-16 J")
	fmt.Println()

	fmt.Println("2. Torque from the metric: the ring (a = 5 cm, 1 µC, 10 g, spinning 1000 rad/s) tilts about a diameter in 1 T, released at 0.05 rad")
	tl := Tilt{A: 0.05, E: 1e-6, Mass: 0.01, OmegaSpin: 1000, B: 1, Theta0: 0.05}
	tr, err := RunTilt(tl)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   current %.4g A, moment %.4g A m², torque −μB sin θ from ½iᵀ(∂L/∂θ)i; period %.5g s vs the exact pendulum period 4K(k)/ω₀ = %.5g s (%.1e); ring current drift %.1e; ledger mismatch %.1e J\n",
		tl.Current(), tl.Moment(), tr.Period, tr.Analytic, tr.Period/tr.Analytic-1, tr.CurrentDrift, tr.Mismatch)
	fmt.Println("   harmonics of the tilt (nonlinear pendulum, odd only): |ω₃/ω₁|, |ω₅/ω₁| from the run against the Jacobi series")
	for _, th := range []float64{0.05, 1, 2.5} {
		h, err := TiltHarmonics(Tilt{A: tl.A, E: tl.E, Mass: tl.Mass, OmegaSpin: tl.OmegaSpin, B: tl.B, Theta0: th})
		if err != nil {
			log.Fatal(err)
		}
		fmt.Printf("   θ₀ = %.2f rad: fundamental %.4g Hz; 3rd %.4g (exact %.4g), 5th %.3g (exact %.3g), 2nd %.1e\n",
			th, h.Fundamental/(2*math.Pi), h.Measured[0], h.Exact[0], h.Measured[1], h.Exact[1], h.Even)
	}
	fmt.Println()

	fmt.Println("3. The electron as the ring: charge e, total angular momentum of the metric S = ħ/2 (measured: g = 2.00232, μ = 1.00116 μ_B, no size above 1e-18 m)")
	fmt.Printf("   %-54s %9s %8s %8s %8s %10s %10s %10s %10s\n", "ring", "a [m]", "v/c", "μ/μ_B", "g", "KE/mc²", "½LI²/mc²", "W_es/mc²", "Ω/(eB/2m)")
	for _, sc := range Scenarios() {
		n := sc.AtSpin(HBAR / 2)
		mc2 := ME * C0 * C0
		fmt.Printf("   %-54s %9.3g %8.3g %8.4g %8.4g %10.3g %10.3g %10.3g %10.4g\n",
			sc.Name, sc.A, n.V/C0, n.Mu/MuB, n.G, n.KE/mc2, n.WMag/mc2, n.WEs/mc2, n.Precession/(QE/(2*ME)))
	}
	fmt.Printf("   A massless thin ring with S = ħ/2 moves at v/c = π/(2α(ln(8a/r) − 7/4)) = %.4g whatever its radius: the speed is set by the fine-structure constant.\n",
		math.Pi/(2*FineStructure*(math.Log(8/0.01)-1.75)))
	fmt.Println("   Standard counterpart: Dirac + QED, g = 2.00232, Ω/(eB/2m) = g. A ring whose mass circulates with its charge has μ/S ≤ e/2m_e (g ≤ 1):")
	fmt.Println("   its spin precession is Larmor's eB/2m, half the electron's. Only a massless ring reaches the measured moment at S = ħ/2, and it does so")
	fmt.Println("   at a few femtometres with v ≫ c; at v ≤ c (a ≥ ħ/2m_e c) the massless ring's moment is far below μ_B. All radii lie above the 1e-18 m bound.")
	fmt.Println()

	fmt.Println("3b. Natural frequencies of the rings at S = ħ/2 in B = 1 T [Hz]: rotation of the charge, spin precession (= magnetic resonance),")
	fmt.Println("    tilt libration (rings with mass), electrical self-resonance 1/2π√(LC), and a wave once round the ring at c")
	fmt.Printf("   %-54s %10s %10s %10s %10s %10s\n", "ring", "rotation", "precession", "libration", "LC", "c/2πa")
	for _, sc := range Scenarios() {
		f := sc.Frequencies(HBAR/2, 1)
		lib := "—"
		if f.Libration > 0 {
			lib = fmt.Sprintf("%10.3g", f.Libration)
		}
		fmt.Printf("   %-54s %10.3g %10.3g %10s %10.3g %10.3g\n", sc.Name, f.Rotation, f.Precession, lib, f.SelfResonance, f.Circumferential)
	}
	fmt.Printf("   %-54s %10s %10.4g %10s %10s %10.4g\n", "measured / reference", "", GammaElectronHz, "", "", ComptonHz)
	fmt.Println("   The Compton ring's c/2πa is the Compton frequency m_e c²/h by identity (a = ħ/m_e c); its charge turns at about half of it (v ≈ c/2).")
	fmt.Println("   The precession row is the ESR frequency each ring would show; only the fitted massless ring reaches the measured 28.02 GHz/T, by construction.")
	kl := ThinRingInductance(tl.A, tl.A/100)
	fmt.Printf("   Kelvin's 5 cm ring of part 2 (1 µC, 10 g, 1000 rad/s, 1 T): rotation %.3g Hz, precession (μ/S)B/2π = e·B/4πm = %.3g Hz, libration %.4g Hz, LC %.3g MHz\n",
		tl.OmegaSpin/(2*math.Pi), tl.E*tl.B/(4*math.Pi*tl.Mass), math.Sqrt(tl.Moment()*tl.B/tl.Inertia())/(2*math.Pi),
		1/(2*math.Pi*math.Sqrt(kl*ThinRingCapacitance(tl.A, tl.A/100)))/1e6)
	fmt.Println()

	if *pngPath != "" {
		if err := DrawRingField(*pngPath); err != nil {
			log.Fatal(err)
		}
		fmt.Printf("4. Field lines of the ring current written to %s (meridional plane; oblique view at twelve azimuths)\n", *pngPath)
	} else {
		fmt.Println("4. Field lines of the ring current: run with -png viz/ring_field.png")
	}
	if *gifPath != "" {
		if err := DrawMotionGIF(*gifPath, kelvin, tr, 90, 8); err != nil {
			log.Fatal(err)
		}
		fmt.Printf("5. Animation written to %s: Kelvin's switch-on (angles scaled to a half turn), the tilt (angle ×8) and its field\n", *gifPath)
	} else {
		fmt.Println("5. Animation: run with -gif viz/ring_motion.gif")
	}
}
