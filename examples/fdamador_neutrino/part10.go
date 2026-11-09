// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 10 (H45): is the weak interaction the ring's own electromagnetic
// field acting on the electron's field (magnetic, or electric)?
//
// The ring's electromagnetic vertex is computed from its own source
// distribution (the polarization current J = ∂P/∂t round a thin ring of
// radius a, and the bound charge −∇·P that the current's azimuthal
// variation leaves), state by state (m = 0, ±1): net charge, charge radius
// ∫r²ρ dV/e, electric dipole, magnetic dipole and the toroidal (anapole)
// dipole T = (1/10)∫[(r·J)r − 2r²J] dV (V. M. Dubovik, V. V. Tugushev,
// Phys. Rep. 187, 145, 1990). Their couplings to an electron are then the
// standard ones (C. Giunti, A. Studenikin, Rev. Mod. Phys. 87, 531, 2015):
// the dipole term dσ/dT = (πα²/m_e²)(μ_eff/μB)²(1/T − 1/E_ν) with
// μ_eff² = μ² + d² (P. Vogel, J. Engel, Phys. Rev. D 39, 3378, 1989),
// multiplied here by the ring's own form factor; the charge radius as the
// shift g_V → g_V + (2/3)M_W²⟨r²⟩sin²θ_W; beside the Standard Model ν–e
// elastic cross-section.

import (
	"fmt"
	"math"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Particle-physics constants (PDG 2024, CODATA 2018), MeV units.
const (
	alphaFS  = 7.2973525693e-3
	meMeV    = 0.51099895
	hcMeVcm  = 1.973269804e-11               // ħc [MeV·cm]
	hcMeVfm  = 197.3269804                   // ħc [MeV·fm]
	gFMeV    = 1.1663787e-11                 // G_F/(ħc)³ [MeV⁻²]
	sin2W    = 0.23122                       // sin²θ_W (MS-bar at M_Z)
	mWMeV    = 80369.2                       // M_W [MeV]
	muBecm   = hcMeVfm / (2 * meMeV) * 1e-13 // μB as an EDM, e·ħ/(2m_e c) [e·cm]
	bDeut    = 2.224566                      // deuteron binding energy [MeV]
	mNucMeV  = 938.918                       // mean nucleon mass [MeV]
	rhoTrip  = 1.75                          // triplet effective range [fm]
	dipoleK  = 2.4951e-25                    // πα²(ħc)²/m_e² [cm²·MeV]: dσ/dT = dipoleK·f²(1/T − 1/E)
	boundBx  = 2.8e-11                       // Borexino 2017 [μB]
	boundXe  = 6.4e-12                       // XENONnT 2022 [μB]
	sigNCref = 1e-42                         // SM ν d → ν n p at 10 MeV, order of magnitude [cm²]
)

// Recorded Part 8d field-solver values (run_output.txt): |p|c with one
// quantum stored in the turning states, at their own frequency.
var p8dGrids = []struct {
	Label    string
	F, PcMuB float64
	Ratio    float64 // f(±1)/f(0)
}{{"r/2", 4.4396e9, 2.89e10, 2.1132}, {"r/3", 4.4354e9, 7.35e10, 2.1561}}

// Src is one source element: position R, current element I·dl [A·m] and
// charge Q [C].
type Src struct {
	R, J [3]float64
	Q    float64
}

// Mom is the multipole content of a source distribution.
type Mom struct {
	Q    float64    // net charge [C]
	P    [3]float64 // electric dipole Σ Q r [C·m]
	M    [3]float64 // magnetic dipole ½ Σ r × J [A·m²]
	T    [3]float64 // toroidal dipole (1/10) Σ [(r·J) r − 2 r² J] [A·m³]
	QR2  float64    // Σ Q r² [C·m²]; ⟨r²⟩ = QR2/e
	Scal float64    // Σ |Q| + Σ |J| scale for round-off judgements
}

// Moments sums the moments of discrete sources.
func Moments(src []Src) Mom {
	var m Mom
	for _, s := range src {
		r, j := s.R, s.J
		rj := r[0]*j[0] + r[1]*j[1] + r[2]*j[2]
		r2 := r[0]*r[0] + r[1]*r[1] + r[2]*r[2]
		m.Q += s.Q
		m.QR2 += s.Q * r2
		for k := 0; k < 3; k++ {
			m.P[k] += s.Q * r[k]
			m.T[k] += 0.1 * (rj*r[k] - 2*r2*j[k])
		}
		m.M[0] += 0.5 * (r[1]*j[2] - r[2]*j[1])
		m.M[1] += 0.5 * (r[2]*j[0] - r[0]*j[2])
		m.M[2] += 0.5 * (r[0]*j[1] - r[1]*j[0])
		m.Scal += math.Abs(s.Q) + math.Sqrt(j[0]*j[0]+j[1]*j[1]+j[2]*j[2])
	}
	return m
}

func norm3(v [3]float64) float64 { return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) }

// RingSources is the thin ring of radius a in state n (0, ±1) at phase
// θ = ωt: current I(φ) = I0·cos(nφ − θ) round the ring and the line charge
// that continuity (∂λ/∂t = −(1/a)∂I/∂φ) leaves, λ = (n·I0/(aω))·cos(nφ − θ).
func RingSources(a float64, n int, i0, w, theta float64, np int) []Src {
	out := make([]Src, np)
	dphi := 2 * math.Pi / float64(np)
	for k := 0; k < np; k++ {
		phi := (float64(k) + 0.5) * dphi
		s, c := math.Sincos(phi)
		cur := i0 * math.Cos(float64(n)*phi-theta)
		lam := float64(n) * i0 / (a * w) * math.Cos(float64(n)*phi-theta)
		out[k] = Src{R: [3]float64{a * c, a * s, 0}, J: [3]float64{-s * cur * a * dphi, c * cur * a * dphi, 0}, Q: lam * a * dphi}
	}
	return out
}

// ToroidalSolenoid is n turns of radius b round a torus of major radius a,
// each in a meridian plane, carrying current i; np segments per turn.
// Its toroidal dipole is exactly ½·n·i·π·b²·a along the axis.
func ToroidalSolenoid(a, b, i float64, n, np int) []Src {
	var out []Src
	dpsi := 2 * math.Pi / float64(np)
	for t := 0; t < n; t++ {
		ph := 2 * math.Pi * (float64(t) + 0.5) / float64(n)
		sp, cp := math.Sincos(ph)
		for k := 0; k < np; k++ {
			ps := (float64(k) + 0.5) * dpsi
			s, c := math.Sincos(ps)
			rho := a + b*c
			// tangent of the turn: d/dψ (ρ, z) = (−b sin ψ, b cos ψ)
			dr, dz := -b*s*dpsi, b*c*dpsi
			out = append(out, Src{R: [3]float64{rho * cp, rho * sp, b * s}, J: [3]float64{i * dr * cp, i * dr * sp, i * dz}})
		}
	}
	return out
}

// StateMom is one state's moments over a period (amplitudes and averages).
type StateMom struct {
	N                int
	A, Omega, I0     float64
	QMax, R2Max      float64 // largest |net charge| [C] and |Σ Q r²| [C·m²] over the period
	PAmp, MAmp, TAmp float64 // largest |p| [C·m], |m| [A·m²], |T| [A·m³]
	PAvg, MAvg, TAvg float64 // |time average| of each vector
	MuEff2           float64 // ⟨μ²⟩ + ⟨(pc)²⟩ in μB² (time average)
	ToroidalOverP    float64 // ω|T|/(c²|p|), = (ka)²/5 for the thin ring
	Scale            float64
}

// ComputeState samples the period at 64 phases.
func ComputeState(a float64, n int, i0, w float64) StateMom {
	s := StateMom{N: n, A: a, Omega: w, I0: i0}
	const nph = 64
	var pav, mav, tav [3]float64
	var mu2, p2 float64
	for k := 0; k < nph; k++ {
		th := 2 * math.Pi * float64(k) / nph
		m := Moments(RingSources(a, n, i0, w, th, 720))
		s.Scale = math.Max(s.Scale, m.Scal)
		s.QMax = math.Max(s.QMax, math.Abs(m.Q))
		s.R2Max = math.Max(s.R2Max, math.Abs(m.QR2))
		s.PAmp = math.Max(s.PAmp, norm3(m.P))
		s.MAmp = math.Max(s.MAmp, norm3(m.M))
		s.TAmp = math.Max(s.TAmp, norm3(m.T))
		for j := 0; j < 3; j++ {
			pav[j] += m.P[j] / nph
			mav[j] += m.M[j] / nph
			tav[j] += m.T[j] / nph
		}
		mu2 += norm3(m.M) * norm3(m.M) / nph
		p2 += norm3(m.P) * norm3(m.P) / nph
	}
	s.PAvg, s.MAvg, s.TAvg = norm3(pav), norm3(mav), norm3(tav)
	s.MuEff2 = (mu2 + p2*C0*C0) / (MuB * MuB)
	if s.PAmp > 0 {
		s.ToroidalOverP = w * s.TAmp / (C0 * C0 * s.PAmp)
	}
	return s
}

// BoundChargeK1 is ∫0^{2π} cos ψ dψ / √(4 sin²(ψ/2) + (r/a)²), the thin
// ring's n = 1 potential coefficient (reduced kernel): a line charge
// λ0·cos φ has potential λ0·cos φ·K1/(4πε0) and energy λ0²·a·K1/(8ε0).
func BoundChargeK1(aOverR float64) float64 {
	eps := 1 / aOverR
	const n = 400000
	h := math.Pi / n
	var s float64
	for k := 0; k < n; k++ {
		p := (float64(k) + 0.5) * h
		sh := math.Sin(p / 2)
		s += math.Cos(p) / math.Sqrt(4*sh*sh+eps*eps) * h
	}
	return 2 * s
}

// I0Quantum0 is the m = 0 current amplitude with one quantum ħω stored
// (½LI0² = ħω, Part 1's L).
func I0Quantum0(a, r, w float64) float64 {
	l := Mu0 * a * (math.Log(8*a/r) - 1.75)
	return math.Sqrt(2 * Hbar * w / l)
}

// I0Bound1 is the m = ±1 current amplitude at which the bound charges'
// vacuum field holds the whole electric half of one quantum, ħω/2. Since
// ε0∇·E = −∇·P, that field is the longitudinal part of E, and
// ½ε0∫|E_long|² ≤ ½ε0∫|E|² ≤ U_E = ħω/2; so this is an upper estimate of
// the turning states' electric dipole (exact up to the thin-ring kernel).
func I0Bound1(a, r, w float64) float64 {
	k1 := BoundChargeK1(a / r)
	lam0 := math.Sqrt(4 * Eps0 * Hbar * w / (a * k1))
	return lam0 * a * w
}

// FormFactor2 is the orientation average of the thin ring's form factor
// squared at x = q·a, weighted as a dipole: (3/4)∫ sin³θ [2J1(x sinθ)/(x sinθ)]² dθ.
func FormFactor2(x float64) float64 {
	if x < 1e-6 {
		return 1
	}
	n := 400 + int(40*x)
	if n%2 == 1 {
		n++
	}
	h := math.Pi / float64(n)
	f := func(t float64) float64 {
		s := math.Sin(t)
		y := x * s
		var g float64 = 1
		if y > 1e-8 {
			g = 2 * math.J1(y) / y
		}
		return s * s * s * g * g
	}
	sum := f(0) + f(math.Pi)
	for k := 1; k < n; k++ {
		w := 2.0
		if k%2 == 1 {
			w = 4
		}
		sum += w * f(float64(k)*h)
	}
	return 0.75 * sum * h / 3
}

// ffTable interpolates FormFactor2 on a log grid (for the deuteron
// double integrals).
type ffTable struct {
	lx0, dlx float64
	v        []float64
}

func newFFTable() *ffTable {
	t := &ffTable{lx0: math.Log(1e-4), dlx: math.Log(1e4/1e-4) / 1600}
	for k := 0; k <= 1600; k++ {
		t.v = append(t.v, math.Log(FormFactor2(math.Exp(t.lx0+float64(k)*t.dlx))))
	}
	return t
}

func (t *ffTable) at(x float64) float64 {
	if x <= 1e-4 {
		return 1
	}
	u := (math.Log(x) - t.lx0) / t.dlx
	k := int(u)
	if k >= len(t.v)-1 {
		// beyond the table the average falls as x⁻³
		return math.Exp(t.v[len(t.v)-1]) * math.Pow(1e4/x, 3)
	}
	f := u - float64(k)
	return math.Exp(t.v[k]*(1-f) + t.v[k+1]*f)
}

// qOfT is the momentum transfer to an electron at rest given kinetic
// energy T, √(T² + 2m_eT) [MeV].
func qOfT(t float64) float64 { return math.Sqrt(t*t + 2*meMeV*t) }

// TMax is the largest electron recoil, 2E²/(m_e + 2E) [MeV].
func TMax(e float64) float64 { return 2 * e * e / (meMeV + 2*e) }

// Flavour of the Standard Model ν–e elastic channel.
type Flavour int

const (
	NuE  Flavour = iota // charged and neutral current
	NuMu                // neutral current only (also ν_τ)
)

// WeakCouplings returns g_V, g_A (with the charged current for ν_e), and
// the charge-radius shift of g_V.
func WeakCouplings(f Flavour, r2cm2 float64) (gv, ga float64) {
	gv, ga = 2*sin2W-0.5, -0.5
	if f == NuE {
		gv, ga = gv+1, ga+1
	}
	return gv + DeltaGV(r2cm2), ga
}

// DeltaGV is the charge-radius shift (2/3)M_W²⟨r²⟩sin²θ_W =
// √2πα⟨r²⟩/(3G_F) (Vogel & Engel 1989; Giunti & Studenikin 2015), ⟨r²⟩ in cm².
func DeltaGV(r2cm2 float64) float64 {
	return 2.0 / 3 * mWMeV * mWMeV * r2cm2 / (hcMeVcm * hcMeVcm) * sin2W
}

// WeakDSigmaDT is the SM ν–e elastic cross-section [cm²/MeV]:
// (G_F²m_e/2π)[(g_V ± g_A)² + (g_V ∓ g_A)²(1 − T/E)² − (g_A² − g_V²)m_eT/E²]
// (upper signs ν, lower ν̄; 't Hooft 1971, Vogel & Engel 1989 eq. 1).
func WeakDSigmaDT(f Flavour, anti bool, e, t, r2 float64) float64 {
	if t <= 0 || t > TMax(e)*(1+1e-12) {
		return 0
	}
	gv, ga := WeakCouplings(f, r2)
	if anti {
		ga = -ga
	}
	pre := gFMeV * gFMeV * meMeV / (2 * math.Pi) * hcMeVcm * hcMeVcm
	y := 1 - t/e
	return pre * ((gv+ga)*(gv+ga) + (gv-ga)*(gv-ga)*y*y - (ga*ga-gv*gv)*meMeV*t/(e*e))
}

// DipoleDSigmaDT is the dipole term (Vogel & Engel 1989 eq. 2) with
// μ_eff² = f2 μB² and the ring's form factor at q·a (a = 0: point) [cm²/MeV].
func DipoleDSigmaDT(f2, e, t, a float64) float64 {
	if t <= 0 || t > TMax(e)*(1+1e-12) {
		return 0
	}
	ff := 1.0
	if a > 0 {
		ff = FormFactor2(qOfT(t) * a * 1e15 / hcMeVfm)
	}
	return dipoleK * f2 * (1/t - 1/e) * ff
}

// integrateT integrates g(T) over [t0, TMax(e)] in ln T (Simpson).
func integrateT(g func(t float64) float64, e, t0 float64) float64 {
	t1 := TMax(e)
	if t0 >= t1 {
		return 0
	}
	const n = 800
	h := math.Log(t1/t0) / n
	var s float64
	for k := 0; k <= n; k++ {
		t := t0 * math.Exp(float64(k)*h)
		if k == n {
			t = t1 // exp may round past T_max, where the integrands are cut to 0
		}
		w := 2.0
		if k == 0 || k == n {
			w = 1
		} else if k%2 == 1 {
			w = 4
		}
		s += w * g(t) * t
	}
	return s * h / 3
}

// SigmaWeak is the total SM cross-section above recoil t0 [cm²].
func SigmaWeak(f Flavour, anti bool, e, t0 float64) float64 {
	return integrateT(func(t float64) float64 { return WeakDSigmaDT(f, anti, e, t, 0) }, e, t0)
}

// SigmaDipole is the dipole cross-section above recoil t0 [cm²].
func SigmaDipole(f2, e, t0, a float64) float64 {
	return integrateT(func(t float64) float64 { return DipoleDSigmaDT(f2, e, t, a) }, e, t0)
}

// SigmaE1 is the deuteron's E1 photodisintegration cross-section in the
// zero-range (Bethe–Peierls) form with the effective-range factor
// 1/(1 − κρ_t) (H. A. Bethe, C. Longmire, Phys. Rev. 77, 647, 1950):
// (8π/3)α(ħc)²/(Mc²)·√B(ω − B)^{3/2}/ω³/(1 − κρ_t) [cm²].
func SigmaE1(w float64) float64 {
	if w <= bDeut {
		return 0
	}
	kappa := math.Sqrt(mNucMeV*bDeut) / hcMeVfm
	fm2 := 8 * math.Pi / 3 * alphaFS * hcMeVfm * hcMeVfm / mNucMeV * math.Sqrt(bDeut) * math.Pow(w-bDeut, 1.5) / (w * w * w) / (1 - kappa*rhoTrip)
	return fm2 * 1e-26
}

// SigmaDeutEM is the deuteron break-up ν d → ν n p by the dipole moment
// (μ_eff² = f2 μB²) through the Coulomb (longitudinal) channel in the dipole
// approximation, R_L(q, ω) = q²σ_E1(ω)/(4π²αω) (Siegert): dσ/dω =
// κ²σ_E1(ω)/(4π²ω)·∫dq² F²(qa), κ² = f2·πα/m_e², q from ω to 2E − ω.
// Order of magnitude: the transverse (real-photon) channel adds a
// comparable amount [cm²].
func SigmaDeutEM(f2, e, a float64, ff *ffTable) float64 {
	if e <= bDeut {
		return 0
	}
	kap2 := f2 * math.Pi * alphaFS / (meMeV * meMeV)
	const nw, nq = 400, 400
	hw := (e - bDeut) / nw
	var s float64
	for i := 0; i <= nw; i++ {
		w := bDeut + float64(i)*hw
		if w <= bDeut || w >= e {
			continue
		}
		q0, q1 := w, 2*e-w
		var g float64
		if a == 0 {
			g = q1*q1 - q0*q0
		} else {
			hq := math.Log(q1/q0) / nq
			for k := 0; k <= nq; k++ {
				q := q0 * math.Exp(float64(k)*hq)
				wt := 2.0
				if k == 0 || k == nq {
					wt = 1
				} else if k%2 == 1 {
					wt = 4
				}
				g += wt * 2 * q * q * ff.at(q*a*1e15/hcMeVfm)
			}
			g *= hq / 3
		}
		wt := 2.0
		if i == 0 || i == nw {
			wt = 1
		} else if i%2 == 1 {
			wt = 4
		}
		s += wt * kap2 / (4 * math.Pi * math.Pi * w) * SigmaE1(w) * g
	}
	return s * hw / 3
}

// ProbeTwin is the kronspice twin: the bench ring as one loop (C charged,
// L, R_rad + R_d) and the electron as a small coaxial probe loop whose
// moment is held fixed by a DC current source, coupled by a K element.
type ProbeTwin struct {
	M, K, Ie                  float64
	F0, QWith, QWithout       float64
	MaxDiff                   float64 // max |I_ring(with) − I_ring(without)| / max |I_ring|
	WorkSource, WorkPredicted float64 // ∫V·I_e dt and M·I_e·ΔI_ring from the first computed step [J]
	PeakInteraction, URing0   float64 // max |M·I_ring·I_e| and the ring's start energy [J]
	MeanPowerOverPeak         float64 // |⟨V·I_e⟩| over the last 10 periods / peak |V·I_e|
	LedgerRel                 float64
	Time, IRing, PElectron    []float64
}

// RunProbeTwin runs the twin for `periods` at T/200, with and without the
// probe current.
func RunProbeTwin(e ERing, periods float64) (ProbeTwin, error) {
	var tw ProbeTwin
	b := e.A / 10
	lp := Mu0 * b * (math.Log(80) - 2)
	tw.M = em.MutualCircles(e.A, b, e.A)
	tw.K = tw.M / math.Sqrt(e.L()*lp)
	tw.Ie = 0.01
	w := e.Omega0()
	tw.F0 = w / (2 * math.Pi)
	T := 2 * math.Pi / w
	run := func(ie float64) (*circuitRun, error) {
		bb := circuit.NewBuilder("E-ring and electron probe loop (fdamador_neutrino Part 10)")
		bb.CapacitorIC("C", "a", "b", e.C(), 1)
		bb.InductorIC("L", "b", "c", e.L(), 0)
		bb.Resistor("R", "c", "a", e.RRad()+e.RDiel())
		bb.Resistor("RG", "a", "0", 1e-3)
		bb.CurrentSource("IE", "0", "n", ie)
		bb.InductorIC("LP", "n", "0", lp, ie)
		// a shunt so the UIC start is not a current source in series with a
		// pinned inductor current; it carries V/1e8 ≈ 1e-10 of I_e
		bb.Resistor("RP", "n", "0", 1e8)
		bb.Mutual("K", "L", "LP", tw.K)
		ckt, err := bb.BuildTransient()
		if err != nil {
			return nil, err
		}
		h := T / 200
		wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: h, Stop: periods * T, MaxStep: h, UIC: true}, analysis.TranOptions{MaxSteps: 5_000_000})
		if err != nil {
			return nil, err
		}
		i, _ := wf.ElementCurrent("L")
		v, _ := wf.NodeVoltage("n")
		return &circuitRun{t: wf.Time, i: i, v: v, led: wf.Ledger.MaxRelativeResidual}, nil
	}
	with, err := run(tw.Ie)
	if err != nil {
		return tw, err
	}
	without, err := run(0)
	if err != nil {
		return tw, err
	}
	qOf := func(r *circuitRun) float64 {
		var xs, ys []float64
		for k := 1; k+1 < len(r.i); k++ {
			a := math.Abs(r.i[k])
			if r.t[k] > T && a > math.Abs(r.i[k-1]) && a >= math.Abs(r.i[k+1]) {
				xs, ys = append(xs, r.t[k]), append(ys, math.Log(a))
			}
		}
		return -w / (2 * fitSlope(xs, ys))
	}
	tw.QWith, tw.QWithout = qOf(with), qOf(without)
	var imax float64
	for k := range with.i {
		imax = math.Max(imax, math.Abs(with.i[k]))
	}
	n := len(with.i)
	if len(without.i) < n {
		n = len(without.i)
	}
	for k := 0; k < n; k++ {
		tw.MaxDiff = math.Max(tw.MaxDiff, math.Abs(with.i[k]-without.i[k])/imax)
	}
	var pk, sum, dt float64
	t0 := (periods - 10) * T
	// from the first computed point on: the imposed UIC point at t = 0
	// carries no node voltage, so the half-step before it is not booked
	for k := 2; k < len(with.t); k++ {
		p0, p1 := with.v[k-1]*tw.Ie, with.v[k]*tw.Ie
		h := with.t[k] - with.t[k-1]
		tw.WorkSource += 0.5 * (p0 + p1) * h
		pk = math.Max(pk, math.Abs(p1))
		if with.t[k-1] >= t0 {
			sum += 0.5 * (p0 + p1) * h
			dt += h
		}
		tw.PeakInteraction = math.Max(tw.PeakInteraction, math.Abs(tw.M*with.i[k]*tw.Ie))
	}
	tw.MeanPowerOverPeak = math.Abs(sum/dt) / pk
	tw.WorkPredicted = tw.M * tw.Ie * (with.i[len(with.i)-1] - with.i[1])
	tw.URing0 = 0.5 * e.C()
	tw.LedgerRel = math.Max(with.led, without.led)
	tw.Time, tw.IRing = with.t, with.i
	for k := range with.v {
		tw.PElectron = append(tw.PElectron, with.v[k]*tw.Ie)
	}
	return tw, nil
}

type circuitRun struct {
	t, i, v []float64
	led     float64
}

// Part10 holds what the report computed (also for viz10.go).
type Part10 struct {
	Bench, Nu            [3]StateMom // m = 0, +1, −1 at the bench scale and at Part 2's scale
	K1Bench, K1Nu        float64
	RoundQ, RoundR2      float64   // largest |q|/(|p|/a) and |Σλr²|/(|p|a) of the turning states
	PFieldFrac           []float64 // Part 8d field |p| / the bound-charge upper estimate, per grid
	ANu, RNu, WNu        float64
	SigWeak              map[string]float64
	SigRing              map[string]float64
	EffXe, EffBx         [2]float64 // effective moment (with form factor) of m = 0 and ±1 at T = 10 keV, 300 keV [μB]
	AMaxXe               float64    // ring radius meeting XENONnT (point dipole, one quantum of 0.1 eV) [m]
	DGVNaive             float64
	SNOPoint1            float64 // σ(ν d → ν n p) at 10 MeV for a point dipole of 1 μB [cm²]
	SNORing              [2]float64
	SNOMuMax             float64 // point moment for which the EM break-up is 0.1 of the SM NC cross-section [μB]
	EDep                 []float64
	SigEWeak, SigERing   [][]float64
	SlopeWeak, SlopeRing float64
	AntiRatioWeak        float64
	Twin                 ProbeTwin
	ContactRing          float64 // ring field energy at the electron, one radius along the axis [eV]
	ContactWeak          float64 // √2G_F(ħc)³ over one ring's volume [eV]
	Ponder               float64 // time-averaged ponderomotive energy of a free electron at ρ = 2a in the plane [eV]
	Seconds              float64
}

func runPart10(viz string) error {
	t0 := time.Now()
	var p Part10
	e := DefaultERing()
	bench := DefaultBench()
	bres := bench.Run()
	p.ANu = bres.RadiusNeeded
	p.RNu = p.ANu / bench.AoverR
	p.WNu = bres.Omega0
	fmt.Println("\n10. The weak interaction as the ring's own field (H45, user hypothesis 2026-11-08): \"the weak interaction is the")
	fmt.Println("    ring's electromagnetic field acting on the electron's field\" — first as magnetic, then as electric. In the model the")
	fmt.Println("    ring's only sources are the polarization current J = ∂P/∂t round the ring and the bound charge −∇·P that an")
	fmt.Println("    azimuthally varying current leaves (∂λ/∂t = −(1/a)∂I/∂φ). Its electromagnetic vertex is computed from them, state by")
	fmt.Println("    state, and coupled to an electron by the standard formulas (Giunti & Studenikin, Rev. Mod. Phys. 87, 531, 2015),")
	fmt.Println("    beside the Standard Model weak cross-section.")

	// 10a moments
	fmt.Println("10a. Moments of the thin ring (720 elements, 64 phases per period). Formulas: q = Σλ dl, ⟨r²⟩ = Σ λ r² dl / e,")
	fmt.Println("    p = Σ λ r dl, μ = ½Σ r × I dl, toroidal T = (1/10)Σ[(r·I dl)r − 2r²I dl] (Dubovik & Tugushev, Phys. Rep. 187, 145,")
	fmt.Println("    1990). Normalization, one quantum ħω per state: m = 0 by ½LI0² = ħω (Part 1's L, Part 2's moment); m = ±1 by the")
	fmt.Println("    bound charges' vacuum field holding the whole electric half, ħω/2 (ε0∇·E = −∇·P makes that field the longitudinal")
	fmt.Println("    part of E, so ½ε0∫E_long² ≤ U_E: an upper estimate of p; thin-ring n = 1 kernel K1).")
	p.K1Bench = BoundChargeK1(e.A / e.R)
	p.K1Nu = BoundChargeK1(bench.AoverR)
	w0 := e.Omega0()
	states := []int{0, 1, -1}
	for k, n := range states {
		var iB, iN, wB float64
		if n == 0 {
			wB = w0
			iB = I0Quantum0(e.A, e.R, wB)
			iN = I0Quantum0(p.ANu, p.RNu, p.WNu)
		} else {
			wB = w0 * p8dGrids[1].Ratio
			iB = I0Bound1(e.A, e.R, wB)
			iN = I0Bound1(p.ANu, p.RNu, p.WNu)
		}
		p.Bench[k] = ComputeState(e.A, n, iB, wB)
		p.Nu[k] = ComputeState(p.ANu, n, iN, p.WNu)
	}
	for _, g := range p8dGrids {
		pmax := ComputeState(e.A, 1, I0Bound1(e.A, e.R, 2*math.Pi*g.F), 2*math.Pi*g.F).PAmp
		p.PFieldFrac = append(p.PFieldFrac, g.PcMuB*MuB/C0/pmax)
	}
	fmt.Printf("    K1 = %.4f (bench, a/r = 4; thin-ring estimate 2(ln(8a/r) − 2) = %.4f). Amplitudes with one quantum:\n", p.K1Bench, 2*(math.Log(32)-2))
	fmt.Println("    scale  state  |q|max/e     ⟨r²⟩max [cm²]  |μ| [μB]     |pc| [μB]    |p| [e·cm]   |T| [A·m³]   ω|T|/(c²|p|) (ka)²/5  |⟨dipole⟩|/amp")
	for sc, set := range [][3]StateMom{p.Bench, p.Nu} {
		name := []string{"bench", "Part 2"}[sc]
		for k := range set {
			s := set[k]
			ka := s.Omega * s.A / C0
			avg := s.MAvg / s.MAmp // the state's own dipole: magnetic for m = 0, electric for ±1
			if s.N != 0 {
				avg = s.PAvg / s.PAmp
			}
			fmt.Printf("    %-6s %+d     %9.2e    %10.2e     %10.4g   %10.4g   %10.4g   %10.3e   %10.4g  %10.4g   %.1e\n", name, s.N, s.QMax/E,
				s.R2Max/E*1e4, s.MAmp/MuB, s.PAmp*C0/MuB, s.PAmp/E*100, s.TAmp, s.ToroidalOverP, ka*ka/5, avg)
		}
	}
	var rq, rr float64
	for _, set := range [][3]StateMom{p.Bench, p.Nu} {
		for _, s := range set[1:] {
			rq = math.Max(rq, s.QMax/(s.PAmp/s.A))
			rr = math.Max(rr, s.R2Max/(s.PAmp*s.A))
		}
	}
	p.RoundQ, p.RoundR2 = rq, rr
	fmt.Printf("    Round-off: |q| ≤ %.1e of |p|/a and |Σλr²dl| ≤ %.1e of |p|·a in the turning states (m = 0 has λ = 0 identically).\n", rq, rr)
	fmt.Printf("    Field-solver check (Part 8d, recorded): the turning states' |p| is %.3f (r/2) and %.3f (r/3) of the upper estimate\n", p.PFieldFrac[0], p.PFieldFrac[1])
	fmt.Println("    at the bench (below it, as it must be, and not converged between the grids). Charge and ⟨r²⟩ are zero to round-off in")
	fmt.Println("    every state and at every instant: m = 0 has no charge at all, m = ±1 only a dipole (l = 1) pattern, which r² (l = 0)")
	fmt.Println("    does not see. So the ring has no charge radius and no anapole: the only couplings are the dipoles — magnetic for")
	fmt.Println("    m = 0, electric for m = ±1 — and they oscillate (turn) at ω, averaging to zero, as Part 8d found. The toroidal moment")
	fmt.Println("    of the turning states is the (ka)²/5 correction of their electric dipole (computed = formula), 1e-13 at Part 2's ring;")
	fmt.Println("    for m = 0 it is zero (the |T| printed is round-off of I0·a³). The dipoles' time averages are round-off of 64 phases.")

	// 10b cross-sections
	fmt.Println("10b. Electron scattering. Weak (SM): dσ/dT = (G_F²m_e/2π)[(g_V+g_A)² + (g_V−g_A)²(1−T/E)² − (g_A²−g_V²)m_eT/E²], ν_e with the")
	fmt.Println("    charged current (g_V = ½ + 2sin²θ_W, g_A = ½), ν_μ,τ neutral only (g_V = −½ + 2sin²θ_W, g_A = −½), ν̄: g_A → −g_A.")
	fmt.Println("    Ring: the dipole term (Vogel & Engel, PRD 39, 3378, 1989) with μ_eff² = ⟨μ²⟩ + ⟨(pc)²⟩ averaged over the period (the")
	fmt.Println("    oscillation shifts the recoil by ±ħω = ±0.1 eV, negligible) times the ring's own form factor ⟨F²(qa)⟩, q² = T² + 2m_eT:")
	fmt.Println("    the ring (a = 1460 fm) is larger than the probe's ħ/q above T ≈ 10 keV, so the point-dipole formula overstates it.")
	p.SigWeak, p.SigRing = map[string]float64{}, map[string]float64{}
	ring := []float64{p.Nu[0].MuEff2, p.Nu[1].MuEff2}
	fmt.Println("    E_ν      T >     σ weak ν_e     σ weak ν_μ     σ ring m=0 (point)   σ ring m=0 (ring F)   σ ring m=±1 (ring F)   ring/weak(ν_e)")
	for _, ee := range []float64{0.862, 10} {
		for _, tt := range []float64{1e-3, 0.25} {
			if tt >= TMax(ee) {
				continue
			}
			swe, swm := SigmaWeak(NuE, false, ee, tt), SigmaWeak(NuMu, false, ee, tt)
			s0p := SigmaDipole(ring[0], ee, tt, 0)
			s0 := SigmaDipole(ring[0], ee, tt, p.ANu)
			s1 := SigmaDipole(ring[1], ee, tt, p.ANu)
			key := fmt.Sprintf("%g_%g", ee, tt)
			p.SigWeak[key], p.SigRing[key] = swe, s0
			fmt.Printf("    %-6g  %4.0f keV  %.3e      %.3e      %.3e            %.3e             %.3e              %.2e – %.2e\n",
				ee, tt*1e3, swe, swm, s0p, s0, s1, s0/swe, s1/swe)
		}
	}
	fmt.Println("    Shapes: the weak spectrum is flat in T (to (1 − T/E)²), the dipole one rises as 1/T toward low recoil — the measured")
	fmt.Println("    low-energy spectra (Borexino, XENONnT) follow the flat weak shape.")
	p.DGVNaive = DeltaGV(math.Pow(p.ANu*100, 2))
	fmt.Printf("    Charge radius: the ring's ⟨r²⟩ is 0, so g_V is unshifted. For scale: had the ring a charge radius ⟨r²⟩ = a² = %.2g cm²,\n", math.Pow(p.ANu*100, 2))
	fmt.Printf("    the shift (2/3)M_W²⟨r²⟩sin²θ_W would be %.2g (g_V itself is O(1)); the SM value ⟨r²⟩ ≈ 4e-33 cm² shifts it by %.3f.\n", p.DGVNaive, DeltaGV(4e-33))

	// 10c bounds
	fmt.Println("10c. Against the bounds (Part 2 ring: a = 1.46e-12 m, one quantum of 0.1 eV per state):")
	for k, x := range []struct {
		lab string
		t   float64
	}{{"XENONnT 2022 (Aprile et al., PRL 129, 161805), recoil ~10 keV, μ < 6.4e-12 μB", 0.01}, {"Borexino 2017 (Agostini et al., PRD 96, 091103), recoil ~300 keV, μ < 2.8e-11 μB", 0.3}} {
		ff := FormFactor2(qOfT(x.t) * p.ANu * 1e15 / hcMeVfm)
		e0, e1 := math.Sqrt(ring[0]*ff), math.Sqrt(ring[1]*ff)
		bnd := []float64{boundXe, boundBx}[k]
		if k == 0 {
			p.EffXe = [2]float64{e0, e1}
		} else {
			p.EffBx = [2]float64{e0, e1}
		}
		fmt.Printf("    %s:\n      form factor %.3g; effective moment m = 0 %.3g μB (%.1e× the bound), m = ±1 ≤ %.3g μB (%.1e×)\n", x.lab, ff, e0, e0/bnd, e1, e1/bnd)
	}
	p.AMaxXe = p.ANu * math.Pow(boundXe/math.Sqrt(ring[0]), 2.0/3)
	fmt.Printf("    At fixed ħω the moment per quantum grows as a^{3/2} (μ = π√(2ħωa³/(μ0Λ))): to meet XENONnT the ring would need a ≤ %.2g m,\n", p.AMaxXe)
	fmt.Printf("    %.1e× below Part 2's ring, where it would radiate faster still (Q_rad ∝ (ka)⁻³ needs εr ∝ a⁻²).\n", p.ANu/p.AMaxXe)
	fmt.Println("    Charge: model 0 (exact; ∇·D = 0), bound |q| ≲ 1e-21 e (neutrality of matter; Giunti & Studenikin 2015). Charge radius:")
	fmt.Println("    model 0, SM ≈ 4e-33 cm² (Bernabéu, Papavassiliou & Vidal, PRL 89, 101802, 2002), TEXONO −2.1e-32 < ⟨r²⟩ < 3.3e-32 cm²")
	fmt.Println("    (Deniz et al., PRD 81, 072001, 2010). The model meets the charge and radius bounds trivially and fails the dipole")
	fmt.Println("    bounds by ~1e9–1e10; static (time-averaged) moments are 0, as Part 8d found.")

	// 10d SNO
	ff := newFFTable()
	p.SNOPoint1 = SigmaDeutEM(1, 10, 0, ff)
	p.SNORing = [2]float64{SigmaDeutEM(ring[0], 10, p.ANu, ff), SigmaDeutEM(ring[1], 10, p.ANu, ff)}
	p.SNOMuMax = math.Sqrt(0.1 * sigNCref / p.SNOPoint1)
	fmt.Println("10d. SNO (deuterium). Charged current ν_e d → e⁻ p p: an electromagnetic field conserves every particle's charge, so the")
	fmt.Println("    ring cannot turn ν into e⁻ or n into p: model CC rate 0 exactly. SNO measured CC at about a third of NC (P_ee ≈ 0.317;")
	fmt.Println("    total ⁸B flux 5.25e6 cm⁻²s⁻¹ from NC; Aharmim et al., PRC 88, 025501, 2013). Neutral current ν d → ν n p (≥ 2.224 MeV;")
	fmt.Println("    SM: axial Gamow–Teller spin flip): the ring's dipole can break the deuteron through the photon (E1, as in")
	fmt.Println("    photodisintegration). Estimate, Coulomb channel in the dipole approximation, σ_E1 Bethe–Peierls with effective range")
	fmt.Printf("    (peak %.2f mb at %.2f MeV): point dipole of 1 μB at E = 10 MeV σ = %.2g cm² (order; the transverse channel adds as much).\n",
		SigmaE1(2*bDeut)*1e27, 2*bDeut, p.SNOPoint1)
	fmt.Printf("    Ring (with its form factor, q ≥ ω ≥ 2.2 MeV ⇒ qa ≥ 16): m = 0 %.2g cm², m = ±1 ≤ %.2g cm², against the SM NC cross-section\n", p.SNORing[0], p.SNORing[1])
	fmt.Printf("    of order 1e-42 cm² at 10 MeV (Nakamura et al., Nucl. Phys. A 707, 561, 2002): %.0e× too many break-ups. SNO's NC agreeing\n", p.SNORing[0]/sigNCref)
	fmt.Printf("    with the SM to ~10%% allows a point moment ≲ %.1g μB, far weaker than XENONnT.\n", p.SNOMuMax)

	// 10e parity and energy
	fmt.Println("10e. Parity and energy dependence.")
	sNuE := SigmaWeak(NuE, false, 10, 1e-3)
	sAnE := SigmaWeak(NuE, true, 10, 1e-3)
	p.AntiRatioWeak = sAnE / sNuE
	fmt.Printf("    σ(ν̄_e e)/σ(ν_e e) at 10 MeV: weak %.3f; ring 1 (the dipole terms are parity- and C-even: ν and ν̄ of either helicity\n", p.AntiRatioWeak)
	fmt.Println("    scatter alike). The weak interaction acts on left-handed ν only (asymmetry 1; C. S. Wu et al., Phys. Rev. 105, 1413, 1957,")
	fmt.Println("    for parity violation); the ring's electromagnetism gives asymmetry 0, and a dipole flips helicity, so it cannot interfere")
	fmt.Println("    with the weak amplitude.")
	p.EDep = []float64{0.3, 0.862, 3, 10, 30}
	fmt.Println("    E_ν [MeV]   σ weak ν_e (T > 1 keV)   σ/E          σ ring m = 0 (ring F, T > 1 keV)")
	for _, ee := range p.EDep {
		sw := SigmaWeak(NuE, false, ee, 1e-3)
		sr := SigmaDipole(ring[0], ee, 1e-3, p.ANu)
		p.SigEWeak = append(p.SigEWeak, []float64{ee, sw})
		p.SigERing = append(p.SigERing, []float64{ee, sr})
		fmt.Printf("    %-9g   %.4e               %.4e   %.4e\n", ee, sw, sw/ee, sr)
	}
	n := len(p.EDep)
	p.SlopeWeak = math.Log(p.SigEWeak[n-1][1]/p.SigEWeak[n-2][1]) / math.Log(p.EDep[n-1]/p.EDep[n-2])
	p.SlopeRing = math.Log(p.SigERing[n-1][1]/p.SigERing[n-2][1]) / math.Log(p.EDep[n-1]/p.EDep[n-2])
	fmt.Printf("    Slope d ln σ/d ln E over 10–30 MeV: weak %.3f (σ ∝ E, measured from solar to accelerator energies), ring %.3f (the\n", p.SlopeWeak, p.SlopeRing)
	fmt.Println("    dipole cross-section grows only as ln E for a point, and the ring's form factor flattens even that).")

	// 10f twin
	tw, err := RunProbeTwin(e, 40)
	if err != nil {
		return err
	}
	p.Twin = tw
	fmt.Println("10f. kronspice twin: the bench ring as one loop (C charged to 1 V, L, R_rad + R_d) and the electron as a small coaxial loop")
	fmt.Printf("    (radius a/10, one ring radius along the axis, M = %.4g H by em.MutualCircles, k = %.3g) whose moment is held by a DC\n", tw.M, tw.K)
	fmt.Printf("    current source (%.0f mA), 40 periods at T/200. Ring Q with the electron %.5g, without %.5g; ring current changed by\n", tw.Ie*1e3, tw.QWith, tw.QWithout)
	fmt.Printf("    %.1e of its peak (a fixed moment induces nothing back). Work done on the electron's source ∫V·I_e dt = %.3e J against\n", tw.MaxDiff, tw.WorkSource)
	fmt.Printf("    M·I_e·ΔI_ring = %.3e J: the interaction energy M·I_ring·I_e (peak %.2e J, %.1e of the ring's 1 V start energy) passes\n", tw.WorkPredicted, tw.PeakInteraction, tw.PeakInteraction/tw.URing0)
	fmt.Printf("    back and forth and its time average is %.1e of the peak power over the last 10 periods. Ledger %.1e. So the ring's\n", tw.MeanPowerOverPeak, tw.LedgerRel)
	fmt.Println("    field gives an electron no static potential; it acts only through the oscillating (sideband) terms of 10b.")
	p.ContactWeak = nuPair(bres).ContactV
	i0 := I0Quantum0(p.ANu, p.RNu, p.WNu)
	bAxis := Mu0 * i0 * p.ANu * p.ANu / (2 * math.Pow(2*p.ANu*p.ANu, 1.5))
	p.ContactRing = MuB * bAxis / E
	aphi := Mu0 * i0 * math.Pi * p.ANu * p.ANu / (4 * math.Pi * 4 * p.ANu * p.ANu)
	p.Ponder = E * E * aphi * aphi / (4 * 9.1093837015e-31) / E
	fmt.Printf("    At Part 2's scale: the m = 0 ring's field on a μB moment one ring radius along the axis is ±%.3g eV, oscillating at ω0\n", p.ContactRing)
	fmt.Printf("    (time average 0), against the weak contact √2G_F(ħc)³ spread over the ring's volume, %.2g eV (Part 9e). The one\n", p.ContactWeak)
	fmt.Printf("    time-averaged term, the ponderomotive energy e²A²/4m_e of a free electron in the ring's induced E at ρ = 2a in its plane\n")
	fmt.Printf("    (dipole A_φ), is %.2g eV and repulsive, a second-order effect falling as ρ⁻⁴.\n", p.Ponder)

	fmt.Println("    Status H45: bounded. The ring's own field gives no charge, no charge radius and no anapole (exact, by its symmetry);")
	fmt.Printf("    its couplings to an electron are dipoles (magnetic for m = 0, electric for m = ±1) of ~0.05–0.08 μB per quantum, which\n")
	fmt.Printf("    oscillate (no static potential, kronspice twin) and scatter electrons with a 1/T spectrum %.0e× the weak rate at\n", p.SigRing["0.862_0.25"]/p.SigWeak["0.862_0.25"])
	fmt.Printf("    862 keV (T > 250 keV), %.0e× the XENONnT bound in moment, equal for ν and ν̄ (weak: %.2f), flat in E (weak ∝ E),\n", p.EffXe[0]/boundXe, p.AntiRatioWeak)
	fmt.Println("    and with no charged current. Standard prediction: W and Z exchange, V−A (parity-violating), σ ∝ G_F²E, moments")
	fmt.Println("    ~3e-20 μB, ⟨r²⟩ ≈ 4e-33 cm² — every column measured agrees with it.")
	p.Seconds = time.Since(t0).Seconds()
	fmt.Printf("   (Part 10: %s)\n", time.Duration(p.Seconds*1e9).Round(time.Second))
	if viz != "" {
		if err := writeViz10(viz, p); err != nil {
			return err
		}
		fmt.Println("\nviz data for Part 10 written to", viz, "— run python3 viz/render_part10.py")
	}
	return nil
}
