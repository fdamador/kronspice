// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command toroidal_coherence puts the bench tests proposed in J. Konstapel,
// "Toroidal Coherence Engineering: From Hall's Photon Theory to Practical
// Propulsion Systems" (Leiden, 1 Nov 2025, 17 pp.; register entry H35)
// beside the standard prediction for each.
//
// The paper claims that heated fiber-optic toroidal coils, phase-locked to
// the Earth's electromagnetic field, invert a scalar "Sn field" (identified
// with the scalar part of Maxwell's quaternion equations) and reduce the
// inertial mass of the rig. Its Tier 1 bench (§6.1) has three measurements:
//
//  1. radiated field with harmonics at (3n ± 1) × the fundamental, called a
//     "topological signature" of toroidal photon organisation;
//  2. a weight change of 0.001 % to 0.1 % of the rig on a precision balance;
//  3. discrete infrared lines (mode-locking) in the fiber's emission.
//
// The paper gives no equation with a value for its coupling σ, so the
// hypothesis itself cannot be put into the simulator (CLAUDE.md §3f: a
// hypothetical model needs stated equations). What is computed here is
// the standard prediction for each measurement and the size of the
// ordinary effects a bench must remove before reading a result:
//
//	Part 1  three coils 120° apart with an ordinary asymmetric loss, in
//	        kronspice transient with the ledger: the combined field carries
//	        exactly the 3n ± 1 harmonics, a single coil carries all of them;
//	Part 2  the closed fiber loop: its free spectral range and the
//	        coherence length of 1000 K thermal light;
//	Part 3  the balance budget for the Tier 1 rig (declared values);
//	Part 4  the paper's own figures recomputed;
//	Part 5  Hall's "FTIR lensing" of galaxy images: total reflection at
//	        hydrogen clouds against the observed lens angles.
package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Physical constants (CODATA 2018 exact or recommended values).
const (
	c0      = 299792458.0    // speed of light [m/s]
	gn      = 9.80665        // standard gravity [m/s²]
	sigmaSB = 5.670374419e-8 // Stefan–Boltzmann [W/(m²K⁴)]
	wienB   = 2.897771955e-3 // Wien displacement constant [m·K]
	rhoAir0 = 1.204          // dry air at 293 K, 101.325 kPa [kg/m³]
	kgToMg  = 1e6            // kilogram to milligram
	tAmb    = 293.15         // laboratory temperature [K]
)

// ---------------------------------------------------------------------------
// Part 1 — the (3n ± 1) harmonic signature.

// Coils is an array of N identical coils driven by a balanced N-phase
// sinusoidal source (phase of coil p: −360°·p/N). Each coil is L in series
// with Rs, closed through a load that is either Rp alone (linear) or Rp in
// parallel with a diode (an ordinary rectifying contact: the asymmetric
// nonlinearity any hot joint, arc or semiconductor junction can give).
type Coils struct {
	N            int
	F, Vpk       float64 // drive frequency [Hz] and amplitude [V]
	L, Rs, Rp    float64 // coil inductance [H], series and load resistance [Ω]
	Diode        bool
	DiodeIS, DiN float64 // SPICE3 diode saturation current [A] and emission coefficient
}

// DefaultCoils is a 10 kHz bench (inside the paper's 1–100 kHz range, §4.3).
var DefaultCoils = Coils{N: 3, F: 10e3, Vpk: 10, L: 1e-3, Rs: 2, Rp: 20, Diode: true, DiodeIS: 1e-12, DiN: 1.5}

// Build returns the circuit.
func (k Coils) Build() (*circuit.Circuit, error) {
	b := circuit.NewBuilder(fmt.Sprintf("toroidal_coherence: %d coils", k.N))
	if k.Diode {
		b.Model("dmod", "d", map[string]float64{"is": k.DiodeIS, "n": k.DiN})
	}
	for p := 0; p < k.N; p++ {
		ph := -360 * float64(p) / float64(k.N)
		src, mid, load := fmt.Sprintf("p%d", p), fmt.Sprintf("q%d", p), fmt.Sprintf("m%d", p)
		b.VoltageWave(fmt.Sprintf("V%d", p), src, "0",
			circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, k.Vpk, k.F, 0, 0, ph}})
		b.Inductor(fmt.Sprintf("L%d", p), src, mid, k.L)
		b.Resistor(fmt.Sprintf("RS%d", p), mid, load, k.Rs)
		b.Resistor(fmt.Sprintf("RP%d", p), load, "0", k.Rp)
		if k.Diode {
			b.Diode(fmt.Sprintf("D%d", p), load, "0", "dmod")
		}
	}
	return b.Build()
}

// Signature is the harmonic content of one run. Single[n] is |I_n|/|I_1| of
// coil 0's current. Fwd[n] and Bwd[n] are the forward- and backward-rotating
// parts of the array's field (the space vector (2/N)Σ i_p e^{j2πp/N}, which
// is the field on the axis of N coils set at 360°/N), relative to the
// forward fundamental. Index 0 is unused.
type Signature struct {
	Single, Fwd, Bwd []float64
	// Ledger over the run [J] and its worst per-point relative residual.
	EnergyIn, EnergyOut, Dissipated, StoredStart, StoredEnd, Mismatch, Residual float64
}

// Run simulates warm-up periods, then analyses the last `window` periods
// for harmonics 1..nh.
func (k Coils) Run(nh int) (Signature, error) {
	var s Signature
	c, err := k.Build()
	if err != nil {
		return s, err
	}
	const warm, window, per = 70, 10, 400 // periods, periods, points per period
	T := 1 / k.F
	stop := float64(warm+window) * T
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: T / per, Stop: stop, MaxStep: T / per}, analysis.TranOptions{MaxSteps: 200000})
	if err != nil {
		return s, err
	}
	l := wf.Ledger
	s.EnergyIn, s.EnergyOut, s.Dissipated = l.EnergyIn, l.EnergyOut, l.EnergyDissipated
	s.StoredStart, s.StoredEnd, s.Mismatch, s.Residual = l.StoredStart, l.StoredEnd, l.IntegrationMismatch, l.MaxRelativeResidual

	// Resample each coil current on a uniform grid over the window.
	m := window * per
	t0 := stop - float64(window)*T
	cur := make([][]float64, k.N)
	for p := range cur {
		i, ok := wf.ElementCurrent(fmt.Sprintf("L%d", p))
		if !ok {
			return s, fmt.Errorf("no current for L%d", p)
		}
		cur[p] = resample(wf.Time, i, t0, T/per, m)
	}
	sv := make([]complex128, m)
	for p := range cur {
		w := cmplx.Exp(complex(0, 2*math.Pi*float64(p)/float64(k.N)))
		for j := range sv {
			sv[j] += complex(2/float64(k.N)*cur[p][j], 0) * w
		}
	}
	single := make([]complex128, m)
	for j := range single {
		single[j] = complex(cur[0][j], 0)
	}
	s.Single, s.Fwd, s.Bwd = make([]float64, nh+1), make([]float64, nh+1), make([]float64, nh+1)
	a1 := cmplx.Abs(dft(single, window, 1))
	f1 := cmplx.Abs(dft(sv, window, 1))
	for n := 1; n <= nh; n++ {
		// A real signal's n-th harmonic splits equally between ±n.
		s.Single[n] = cmplx.Abs(dft(single, window, n)) / a1
		s.Fwd[n] = cmplx.Abs(dft(sv, window, n)) / f1
		s.Bwd[n] = cmplx.Abs(dft(sv, window, -n)) / f1
	}
	return s, nil
}

// resample interpolates (t, y) linearly at t0 + j·dt, j = 0..m−1.
func resample(t, y []float64, t0, dt float64, m int) []float64 {
	out := make([]float64, m)
	k := 0
	for j := range out {
		tj := t0 + float64(j)*dt
		for k < len(t)-2 && t[k+1] < tj {
			k++
		}
		u := (tj - t[k]) / (t[k+1] - t[k])
		out[j] = y[k] + u*(y[k+1]-y[k])
	}
	return out
}

// dft returns the coefficient of e^{+j2π·n·f·t} over a window of whole
// periods (bin n·periods of a length-len(x) DFT).
func dft(x []complex128, periods, n int) complex128 {
	var acc complex128
	m := len(x)
	for j, v := range x {
		acc += v * cmplx.Exp(complex(0, -2*math.Pi*float64(n*periods*j)/float64(m)))
	}
	return acc / complex(float64(m), 0)
}

// ---------------------------------------------------------------------------
// Part 2 — the closed fiber loop (Measurement 3).

// Loop is the Tier 1 fiber coil as an optical ring.
type Loop struct {
	Length float64 // circumference [m] (paper §6.1: 10 m)
	Ng     float64 // group index of silica (declared 1.46; Malitson 1965 gives n ≈ 1.42–1.45 over 1–3 µm)
	T      float64 // coil temperature [K] (paper §6.1: 1000 K)
}

// FSR is the spacing of the ring's longitudinal modes c/(n_g·ℓ) [Hz]: the
// discrete lines every closed loop shows (examples/photon_electron, ring
// lines: AC peaks at n·v/ℓ to 3e-6).
func (l Loop) FSR() float64 { return c0 / (l.Ng * l.Length) }

// WienPeak is the wavelength of peak thermal emission [m].
func (l Loop) WienPeak() float64 { return wienB / l.T }

// CoherenceLength of blackbody light, ≈ λ_peak (its bandwidth is of order
// its centre wavelength; Mandel & Wolf 1995 §13.1) [m].
func (l Loop) CoherenceLength() float64 { return l.WienPeak() }

// ---------------------------------------------------------------------------
// Part 3 — the balance budget (Measurement 2).

// Rig is the Tier 1 bench. The paper gives the heater power (100 kW RF), the
// coil temperature (1000 K) and circumference (10 m); everything else is
// declared here.
type Rig struct {
	Mass       float64 // coil, housing and heater applicator on the balance [kg] (declared)
	Power      float64 // heater power [W] (paper)
	T          float64 // coil temperature [K] (paper)
	Length     float64 // coil circumference [m] (paper)
	Dia        float64 // outer diameter of the insulated coil bundle [m] (declared)
	Emissivity float64 // coil surface (declared; silica ≈ 0.8 in the thermal IR)
	Cp         float64 // mean specific heat 293–1000 K [J/(kg·K)] (fused silica ≈ 1000)
	LeadAmps   float64 // DC feed current crossing onto the balance [A] (declared)
	BEarth     float64 // horizontal geomagnetic field [T]
}

// DefaultRig is the declared Tier 1 bench.
var DefaultRig = Rig{Mass: 50, Power: 100e3, T: 1000, Length: 10, Dia: 0.05, Emissivity: 0.8, Cp: 1000, LeadAmps: 250, BEarth: 50e-6}

// Budget lists forces on the balance as equivalent mass [kg] (F/g).
type Budget struct {
	ClaimLow, ClaimHigh float64 // paper: 0.001 % and 0.1 % of Mass (sign: lighter)
	HeatMass            float64 // standard prediction: +E/c², heavier
	RadiatedPower       float64 // σεA(T⁴ − T0⁴) [W]
	RecoilMax           float64 // P/c, all emission one way
	LeadPerMetre        float64 // I·B per metre of feed lead
	AirBuoyancy         float64 // ρ0·V·g·(1 − T0/T): buoyancy lost when the air round the coil heats to T
}

// Budget computes the balance budget.
func (r Rig) Budget() Budget {
	area := math.Pi * r.Dia * r.Length
	vol := math.Pi * r.Dia * r.Dia / 4 * r.Length
	return Budget{
		ClaimLow:      1e-5 * r.Mass,
		ClaimHigh:     1e-3 * r.Mass,
		HeatMass:      r.Mass * r.Cp * (r.T - tAmb) / (c0 * c0),
		RadiatedPower: r.Emissivity * sigmaSB * area * (math.Pow(r.T, 4) - math.Pow(tAmb, 4)),
		RecoilMax:     r.Power / c0 / gn,
		LeadPerMetre:  r.LeadAmps * r.BEarth / gn,
		AirBuoyancy:   rhoAir0 * vol * (1 - tAmb/r.T),
	}
}

// ---------------------------------------------------------------------------
// Part 4 — the paper's own figures recomputed.

// Accel returns F/m [m/s²].
func Accel(force, mass float64) float64 { return force / mass }

// ResidualMassFraction is the fraction of the rest mass a vehicle of mass m
// may keep to reach acceleration a under force F: F/(m·a).
func ResidualMassFraction(force, mass, a float64) float64 { return force / (mass * a) }

// ---------------------------------------------------------------------------
// Part 5 — Hall's "FTIR lensing": galaxy images bent by (frustrated) total
// internal reflection at hydrogen clouds instead of by gravity.

// loschmidt is the number density of an ideal gas at 273.15 K, 101.325 kPa [m⁻³].
const loschmidt = 2.6867811e25

// H2Refractivity returns n − 1 of molecular hydrogen at 273.15 K, 101.325 kPa
// for vacuum wavelength lambda [m] (Peck & Huang, J. Opt. Soc. Am. 67, 1550,
// 1977: n − 1 = 0.0148956/(180.7 − σ²) + 0.0049037/(92 − σ²), σ in µm⁻¹).
func H2Refractivity(lambda float64) float64 {
	s2 := math.Pow(1e-6/lambda, 2)
	return 0.0148956/(180.7-s2) + 0.0049037/(92-s2)
}

// TIRMaxDeflection is the largest deflection [rad] that total reflection at
// a sharp boundary between hydrogen of number density nH2 [m⁻³] and vacuum
// can give: only rays within the critical grazing angle √(2(n−1)) of the
// surface reflect, and a reflection turns them by twice that.
func TIRMaxDeflection(nH2, lambda float64) float64 {
	dn := H2Refractivity(lambda) * nH2 / loschmidt
	return 2 * math.Sqrt(2*dn)
}

// radToArcsec converts radians to arcseconds.
const radToArcsec = 180 / math.Pi * 3600

// ---------------------------------------------------------------------------

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "toroidal_coherence:", err)
		os.Exit(1)
	}
}

func run() error {
	const nh = 9
	fmt.Println("H35 — Toroidal coherence propulsion (Konstapel 2025): Tier 1 bench against the standard prediction")
	fmt.Println()
	fmt.Println("Part 1 — the (3n ± 1) harmonic signature (paper §6.1, Measurement 1)")
	k := DefaultCoils
	fmt.Printf("  %d coils at 120°, %.0f kHz, %.0f V, L %.0f mH, Rs %.0f Ω, load %.0f Ω ∥ diode (IS %.0e A, N %.1f)\n",
		k.N, k.F/1e3, k.Vpk, k.L*1e3, k.Rs, k.Rp, k.DiodeIS, k.DiN)
	s, err := k.Run(nh)
	if err != nil {
		return err
	}
	lin := k
	lin.Diode = false
	sl, err := lin.Run(nh)
	if err != nil {
		return err
	}
	fmt.Println("    n   one coil |I_n|/|I_1|   array field fwd   array field bwd   linear control (one coil)")
	for n := 1; n <= nh; n++ {
		tag := ""
		if n%3 == 0 {
			tag = "  ← multiple of 3"
		}
		fmt.Printf("  %3d   %18.3e   %15.3e   %15.3e   %15.3e%s\n", n, s.Single[n], s.Fwd[n], s.Bwd[n], sl.Single[n], tag)
	}
	fmt.Printf("  ledger (diode run): in %.4g J − returned %.4g J = dissipated %.4g J + Δstored %.4g J; mismatch %.2g J, worst point residual %.2g\n",
		s.EnergyIn, s.EnergyOut, s.Dissipated, s.StoredEnd-s.StoredStart, s.Mismatch, s.Residual)
	fmt.Println("  standard prediction: any balanced three-coil array with any ordinary nonlinearity radiates only 3n ± 1")
	fmt.Println("  (3n+1 rotating forward, 3n−1 backward; multiples of 3 cancel by symmetry); one coil radiates every n.")
	fmt.Println("  So 3n ± 1 in Tier 2 cannot distinguish the hypothesis; in Tier 1 (one coil) the paper does not say whether")
	fmt.Println("  multiples of 3 are absent, which is the only form that would be unexpected.")
	fmt.Println()

	fmt.Println("Part 2 — the closed fiber loop (Measurement 3)")
	lp := Loop{Length: 10, Ng: 1.46, T: 1000}
	fmt.Printf("  %.0f m loop, n_g %.2f: longitudinal modes every %.2f MHz (c/n_g ℓ) — discrete lines from any closed loop\n",
		lp.Length, lp.Ng, lp.FSR()/1e6)
	fmt.Printf("  thermal light at %.0f K peaks at %.2f µm with coherence length ≈ %.1f µm = %.1e of the loop:\n",
		lp.T, lp.WienPeak()*1e6, lp.CoherenceLength()*1e6, lp.CoherenceLength()/lp.Length)
	fmt.Println("  the loop filters thermal light into its modes but cannot phase-lock it; lines spaced by the FSR are the")
	fmt.Println("  standard prediction, so the test needs a stated line spacing or linewidth that differs from c/(n_g ℓ).")
	fmt.Println()

	fmt.Println("Part 3 — the balance budget (Measurement 2); forces as equivalent mass")
	r := DefaultRig
	b := r.Budget()
	fmt.Printf("  rig (declared): %.0f kg, %.0f kW heater, coil %.0f K, %.0f m × %.0f mm, ε %.1f, feed %.0f A, B %.0f µT\n",
		r.Mass, r.Power/1e3, r.T, r.Length, r.Dia*1e3, r.Emissivity, r.LeadAmps, r.BEarth*1e6)
	fmt.Printf("  claimed weight loss (0.001 %%–0.1 %%)            %10.4g – %.4g mg\n", b.ClaimLow*kgToMg, b.ClaimHigh*kgToMg)
	fmt.Printf("  standard prediction: heat adds E/c² (heavier)  %+10.4g mg\n", b.HeatMass*kgToMg)
	fmt.Printf("  radiation recoil, all emission one way (P/c)  %10.4g mg\n", b.RecoilMax*kgToMg)
	fmt.Printf("  feed lead in the geomagnetic field, per metre %10.4g mg\n", b.LeadPerMetre*kgToMg)
	fmt.Printf("  air buoyancy lost round the hot coil (in air) %10.4g mg\n", b.AirBuoyancy*kgToMg)
	fmt.Printf("  (radiated power at %.0f K: %.0f kW of the %.0f kW heater — the heater sustains the temperature)\n",
		r.T, b.RadiatedPower/1e3, r.Power/1e3)
	fmt.Println("  ordinary effects are as large as the claim: run in vacuum, cancel lead forces (twisted, reversed feed),")
	fmt.Println("  compare with a dummy coil at the same power, and fix the 'coherence protocol' before the run.")
	fmt.Println()

	fmt.Println("Part 4 — the paper's figures recomputed")
	fmt.Printf("  1 kN on 1000 t: %.3g m/s² (paper §4.3: 1 m/s²; ×1000 high, and its 5 and 100 m/s² rows with it)\n", Accel(1e3, 1e6))
	fmt.Printf("  6000 g on 1000 t with 1 kN needs a residual mass fraction %.2g (paper: '>99.99 %%' reduction, i.e. 1e-4)\n",
		ResidualMassFraction(1e3, 1e6, 6000*gn))
	fmt.Println("  coil length: 3–10 km (§4.1), 3,000–5,000 km per coil (§4.3), ~2,000 km in total (§6.3)")
	fmt.Println("  aluminium cladding (melts 933 K) run at 1,500–2,000 K; fused silica softens near 1,940 K (paper: melts 2,400 K)")
	fmt.Println("  lock to the 7.8 Hz Schumann mode within 0.1 %: that mode has Q ≈ 5 (Sentman 1995), a linewidth of ~20 %")
	fmt.Println("  E_eff = (1 − f)γm0c² for v > c: γ = 1/√(1 − v²/c²) is not real above c whatever the mass")
	fmt.Println()
	fmt.Println("Part 5 — Hall's FTIR lensing: total reflection at hydrogen clouds (H₂, Peck & Huang 1977)")
	fmt.Printf("  H₂ at STP, 589 nm: n − 1 = %.4g\n", H2Refractivity(589e-9))
	fmt.Println("  density [cm⁻³]          n − 1     max TIR deflection   ratio 400 nm / 1 µm")
	for _, d := range []float64{1, 1e4, 1e6} {
		dn := H2Refractivity(589e-9) * d * 1e6 / loschmidt
		fmt.Printf("  %-14.0e  %14.3g   %10.3g arcsec   %.4f\n", d, dn, TIRMaxDeflection(d*1e6, 589e-9)*radToArcsec,
			TIRMaxDeflection(d*1e6, 400e-9)/TIRMaxDeflection(d*1e6, 1e-6))
	}
	fmt.Println("  observed strong-lensing arcs and multiple images are separated by ~1–30 arcsec; even a sharp-edged")
	fmt.Println("  dense core (10⁶ cm⁻³) turns light by at most ~1e-3 arcsec, and real cloud edges are graded over")
	fmt.Println("  parsecs, not the sub-wavelength step TIR needs (the FTIR gap must be ≲ λ). Reflection is chromatic")
	fmt.Println("  (blue turned ~2 % more than near-IR); gravitational lensing is achromatic. Test: image positions of a")
	fmt.Println("  multiply imaged quasar across bands — positions that agree across bands rule the reflection reading out.")
	fmt.Println()
	fmt.Println("Register (CLAUDE.md §1): the paper's toroidal electron is H23 (bounded: charge 0.910 e, anomaly 17 % low, not")
	fmt.Println("self-confined under linear Maxwell); inertial-mass reduction of a craft is H18 (thrust ≤ P/c); the quaternion")
	fmt.Println("scalar part adds no physics (§3c). MICROSCOPE bounds composition-dependent free fall at ~1e-15 (Touboul et al. 2022).")
	return nil
}
