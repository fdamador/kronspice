// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command light_hydrogen builds the hydrogen atom from the light models of
// the project, hypothesis H32 (CLAUDE.md §1, user hypothesis 2026-11-05):
// the electron as the toroidal photon of Williamson & van der Mark (H23,
// examples/photon_electron), the proton as a circuit of light (H31,
// examples/nucleon_circuits: one light ring, three quark rings, or
// Skyrme's knot), and the knot of H29 (examples/knotted_light,
// examples/ether_knot) as what would hold the electron's light together.
// Each prediction is set beside the measured hydrogen and beside the
// standard account (Dirac electron, QED, measured proton).
//
// The constituents, as the three examples give them:
//
//	electron (H23)  one Compton wavelength of light on a double loop of
//	                radius λ_C/4π = 193 fm; charge q = √(3ε0ħc)/2π = 0.910 e
//	                (eq. 5), g = 2(1 + α′/2π), α′ = (q/e)²α (eqs. 8, 12)
//	proton (H31)    a light ring of radius ħ/(m_p c) = 0.210 fm carrying +e
//	                at c (moment 1 μN); or three quark rings fitted to the
//	                proton's moment (2.793 μN, rms radius 0.678 fm); or
//	                Skyrme's B = 1 knot (1.880 μN, computed in skyrme/)
//	binding (H29)   the ether knot: energy 547.5√(κ2κ4), core radius 1.2ℓ,
//	                squashing tone ω = 1.416 c/ℓ (ether_knot)
//
// Parts: (1) scales — both rings are far smaller than the atom, so at the
// atom's scale they are point charges and moments; (2) the bound state as
// Kron's 1945 Schrödinger network (radial ladder, network.go), solved by
// Sturm counting and twinned in kronspice AC; (3) the spectrum, neutrality
// and the proton-size shift of 1S–2S for each model; (4) the 21 cm
// hyperfine line, with the contact field of the electron's spin cloud
// built as coaxial current loops round the proton ring in kronspice
// (hyperfine.go); (5) the knot that would hold the electron, its
// parameters fixed by the electron's mass and size, and its first
// vibration as a predicted excited electron.
//
// Part 6 (spectrum.go) reads the spectrum from the frequency response of
// the whole network: an AC sweep shows every level as a pole, the poles
// are pinned and extrapolated, and their differences are set beside the
// measured Lyman and Balmer lines and 1S–2S.
//
//	go run ./examples/light_hydrogen                 # about 40 s
//	go run ./examples/light_hydrogen -svg viz        # also the report's figures and the 3-D scene
package main

import (
	"flag"
	"fmt"
	"math"
	"path/filepath"
)

// Constants, SI (CODATA 2022) and particle data (PDG 2024).
const (
	E       = 1.602176634e-19   // C (exact)
	H       = 6.62607015e-34    // J s (exact)
	Hbar    = H / (2 * math.Pi) // J s
	C0      = 299792458.0       // m/s (exact)
	Eps0    = 8.8541878188e-12  // F/m
	Mu0     = 1.25663706127e-6  // N/A²
	AlphaFS = 7.2973525643e-3
	RInf    = 10973731.568157  // m⁻¹, Rydberg constant
	MeV     = 1.602176634e-13  // J
	MpMeV   = 938.27208943     // proton rest energy
	MeMeV   = 0.51099895069    // electron rest energy
	MuPmeas = 2.79284734463    // proton moment / μN
	GeMeas  = 2.00231930436092 // |g| of the free electron
	RpMeas  = 0.84075e-15      // proton charge radius, m (CODATA 2022)
)

// Measured hydrogen (the targets).
const (
	HFMeasHz   = 1420405751.768      // 21 cm line (Hellwig et al., IEEE Trans. Instrum. Meas. 19, 200, 1970)
	OneSTwoSHz = 2466061413187.018e3 // 1S–2S (Parthey et al., PRL 107, 203001, 2011), ±10 Hz
	LymanANm   = 121.5668            // Lyman-α, vacuum, fine-structure-weighted (NIST ASD)
	IonEV      = 13.598434599702     // ionization energy of H (NIST ASD)
	Neutrality = 1e-21               // |q_p + q_e|/e bound (Bressi et al., PRA 83, 052101, 2011)
)

// Derived constants.
var (
	MeKg   = MeMeV * MeV / (C0 * C0)
	MpKg   = MpMeV * MeV / (C0 * C0)
	MuB    = E * Hbar / (2 * MeKg) // Bohr magneton [J/T]
	MuN0   = E * Hbar / (2 * MpKg) // nuclear magneton [J/T]
	A0     = Hbar / (MeKg * C0 * AlphaFS)
	RyJ    = RInf * H * C0 // Rydberg energy [J]
	Recoil = 1 + MeKg/MpKg // m_e/μ for hydrogen
)

// Electron is a model electron: charge Q [e], g-factor, size R [m].
type Electron struct {
	Name string
	Q, G float64
	R    float64
}

// Moment is the electron's spin moment (g/2)(q/e)μ_B [J/T].
func (e Electron) Moment() float64 { return e.G / 2 * e.Q * MuB }

// WvdM is the toroidal-photon electron of H23 (Williamson & van der Mark
// 1997, eqs. 5, 8, 12; examples/photon_electron Part 1).
func WvdM() Electron {
	q := math.Sqrt(3*Eps0*Hbar*C0) / (2 * math.Pi) / E
	ap := q * q * AlphaFS
	return Electron{Name: "toroidal photon (H23)", Q: q, G: 2 * (1 + ap/(2*math.Pi)), R: Hbar / (2 * MeKg * C0)}
}

// WvdMUnit is the same electron with its charge set to e, as the paper's
// eq. 13 compares (g = 2(1 + α/2π)).
func WvdMUnit() Electron {
	return Electron{Name: "toroidal photon, q set to e", Q: 1, G: 2 * (1 + AlphaFS/(2*math.Pi)), R: Hbar / (2 * MeKg * C0)}
}

// Measured is the electron as measured (point-like to ~10⁻¹⁸ m).
func Measured() Electron { return Electron{Name: "measured (Dirac + QED)", Q: 1, G: GeMeas} }

// Proton is a model proton: charge Q [e], moment Mu [μN], rms charge
// radius R [m] (0 when the model gives none).
type Proton struct {
	Name  string
	Q, Mu float64
	R     float64
}

// ProtonRingRadius is ħ/(m_p c), the light-ring proton of H31.
func ProtonRingRadius() float64 { return Hbar / (MpKg * C0) }

// QuarkRingRadius is the quark ring of H31: light of energy m_p c²/μ_p on a
// ring of radius ħ/(m_q c).
func QuarkRingRadius() float64 { return ProtonRingRadius() * MuPmeas }

// Protons are the H31 protons and the measured one.
func Protons() []Proton {
	return []Proton{
		{Name: "light ring (H31)", Q: 1, Mu: 1, R: ProtonRingRadius()},
		// uud rings of radius R at z = −R, 0, +R: ⟨r²⟩ = Σq(R² + z²) = (4/3)R².
		{Name: "quark rings (H31)", Q: 1, Mu: MuPmeas, R: QuarkRingRadius() * math.Sqrt(4.0/3)},
		{Name: "Skyrme knot (H31)", Q: 1, Mu: 1.880},
		{Name: "measured", Q: 1, Mu: MuPmeas, R: RpMeas},
	}
}

// Atom is one model hydrogen atom.
type Atom struct {
	E Electron
	P Proton
}

// Z is the coupling q_e·q_p/e² that replaces 1 in every Coulomb quantity.
func (a Atom) Z() float64 { return a.E.Q * a.P.Q }

// Charge is the atom's net charge [e].
func (a Atom) Charge() float64 { return a.P.Q - a.E.Q }

// Level is the Bohr level n with reduced mass [J] (negative: bound).
func (a Atom) Level(n int) float64 { z := a.Z(); return -z * z * RyJ / Recoil / float64(n*n) }

// LymanAlpha is the 1→2 wavelength [nm] of the Bohr levels.
func (a Atom) LymanAlpha() float64 { return H * C0 / (a.Level(2) - a.Level(1)) * 1e9 }

// SizeShift is the finite-size shift of the nS level [J]:
// (2/3n³)(Zα)⁴ m_r c² (r/ƛ_r)², ƛ_r = ħ/(m_r c) (Eides, Grotch &
// Shelyuto, Phys. Rep. 342, 63, 2001, eq. 2.4 with the rms radius).
func (a Atom) SizeShift(n int) float64 {
	mr := MeKg / Recoil
	lam := Hbar / (mr * C0)
	za := a.Z() * AlphaFS
	return 2.0 / (3 * float64(n*n*n)) * math.Pow(za, 4) * mr * C0 * C0 * (a.P.R / lam) * (a.P.R / lam)
}

// BohrRadius is the model atom's radius a0·Recoil/Z [m].
func (a Atom) BohrRadius() float64 { return A0 * Recoil / a.Z() }

// Hyperfine is the Fermi contact splitting of the ground state [Hz]:
// ΔE = (8/3)μ0 μ_e μ_p |ψ(0)|², |ψ(0)|² = 1/(πa³) (Fermi, Z. Phys. 60,
// 320, 1930; Bethe & Salpeter §22).
func (a Atom) Hyperfine() float64 {
	r := a.BohrRadius()
	return 8.0 / 3 * Mu0 * a.E.Moment() * a.P.Mu * MuN0 / (math.Pi * r * r * r) / H
}

func main() {
	svgDir := flag.String("svg", "", "write the report's figures (SVG) and the 3-D scene (JSON) to this directory")
	flag.Parse()
	fmt.Println("light_hydrogen (H32): hydrogen from the light models — electron H23, proton H31, binding H29")

	wv := WvdM()
	pr := Protons()
	ring, quark := pr[0], pr[1]

	// Part 1
	fmt.Println("\n1. Scales: the constituents against the atom")
	fmt.Printf("   electron ring (H23) radius λ_C/4π      %10.1f fm\n", wv.R*1e15)
	fmt.Printf("   proton light ring (H31) radius ħ/m_p c %10.4f fm   (measured charge radius %.5f fm)\n", ring.R*1e15, RpMeas*1e15)
	fmt.Printf("   quark rings (H31), rms charge radius   %10.4f fm\n", quark.R*1e15)
	fmt.Printf("   Bohr radius a0                          %10.0f fm\n", A0*1e15)
	fmt.Printf("   electron ring / a0 = %.4f, proton ring / a0 = %.1e: at the atom's scale both are point charges\n", wv.R/A0, ring.R/A0)
	fmt.Println("   and point moments, so the atom's equation is the Coulomb problem with charges q_e and q_p (Part 2)")

	// Part 2
	fmt.Println("\n2. The bound state as Kron's 1945 Schrödinger network (radial ladder for s states, network.go)")
	fmt.Println("   series inductors = kinetic energy, shunt inductors = −Z/r, unit capacitors = energy; ψ = node voltage")
	for _, z := range []float64{1, wv.Q} {
		e1, e2 := Extrapolated(z)
		fmt.Printf("   Z = %.6f: E_1s = %.9f Ha, E_2s = %.9f Ha (h → 0); exact −Z²/2 = %.9f, −Z²/8 = %.9f\n", z, e1, e2, -z*z/2, -z*z/8)
	}
	tw, err := Twin(1, 0.1, 30)
	if err != nil {
		fmt.Println("   kronspice twin error:", err)
		return
	}
	fmt.Printf("   kronspice twin (h = 0.1 a0, %d nodes, AC resonance by bisection): E_1s = %.10f Ha, network matrix %.10f (diff %.1e)\n",
		tw.Nodes, tw.Circuit, tw.Matrix, tw.Circuit-tw.Matrix)
	fmt.Println("   the network is linear in the potential, so the light models change only Z = q_e q_p / e²")

	// Part 3
	fmt.Println("\n3. Spectrum, neutrality and proton size, model by model")
	atoms := []Atom{{wv, ring}, {wv, quark}, {WvdMUnit(), ring}, {WvdMUnit(), quark}, {Measured(), pr[3]}}
	fmt.Printf("   %-44s %10s %10s %10s %14s\n", "electron + proton", "charge (e)", "1S (eV)", "Lyα (nm)", "1S–2S size")
	ref := Atom{Measured(), pr[3]}
	for _, a := range atoms {
		size := "—"
		if a.P.R > 0 {
			d := (a.SizeShift(1) - a.SizeShift(2)) - (ref.SizeShift(1) - ref.SizeShift(2))
			size = fmt.Sprintf("%+.3f MHz", d/H/1e6)
		}
		fmt.Printf("   %-44s %+10.4f %10.4f %10.3f %14s\n", a.E.Name+" + "+a.P.Name, a.Charge(), -a.Level(1)/E, a.LymanAlpha(), size)
	}
	fmt.Printf("   measured: charge |q| < %.0e e (Bressi et al. 2011); ionization %.4f eV; Lyman-α %.4f nm (NIST)\n", Neutrality, IonEV, LymanANm)
	fmt.Printf("   1S–2S measured to 10 Hz (Parthey et al. 2011): the size column is the shift from the measured-radius value\n")
	fmt.Printf("   standard (Bohr, reduced mass) 1S–2S %.6f THz against measured %.6f THz: the remaining %.1f GHz is Dirac + QED (Lamb)\n",
		(ref.Level(2)-ref.Level(1))/H/1e12, OneSTwoSHz/1e12, (OneSTwoSHz-(ref.Level(2)-ref.Level(1))/H)/1e9)

	// Part 4
	fmt.Println("\n4. The 21 cm line: the electron's spin cloud as coaxial current loops round the proton ring (hyperfine.go)")
	fmt.Printf("   %-44s %9s %9s %12s %10s %9s\n", "electron + proton", "μ_e (μB)", "μ_p (μN)", "kronspice zz", "ν (MHz)", "/measured")
	for _, a := range []Atom{{wv, ring}, {wv, quark}, {wv, pr[2]}, {WvdMUnit(), ring}, {WvdMUnit(), quark}, {WvdMUnit(), pr[2]}, {Measured(), pr[3]}} {
		hf, err := HyperfineCircuit(a, 8)
		if err != nil {
			fmt.Println("   error:", err)
			return
		}
		fmt.Printf("   %-44s %9.6f %9.4f %9.3e J %10.3f %9.5f\n", a.E.Name+" + "+a.P.Name, a.E.Moment()/MuB, a.P.Mu, hf.Interaction, hf.Nu/1e6, hf.Nu/HFMeasHz)
	}
	hf, _ := HyperfineCircuit(Atom{Measured(), pr[3]}, 8)
	fmt.Printf("   circuit: %d electron loops + the proton ring, K couplings to the ring. Interaction ΣM·I_p·I_k from the\n", hf.Loops)
	fmt.Printf("   solved currents against μ_p·B(0) = (2/3)μ0 μ_p μ_e|ψ(0)|²: diff %.1e (the proton ring's finite size)\n", hf.Interaction/hf.Analytic-1)
	fmt.Printf("   ledger stored energy against ½ΣLI² + ΣMII: %.1e; ledger residual %.1e\n", hf.Stored, hf.Residual)
	fmt.Println("   the circuit gives the aligned/anti-aligned (zz) energy ±U; spin ½ (S·S, not S_zS_z) makes the splitting 4U")
	fmt.Printf("   measured %.6f MHz; the standard Fermi value is %.4f%% off (Breit, recoil and proton-structure terms)\n",
		HFMeasHz/1e6, 100*(hf.Nu/HFMeasHz-1))

	// Part 5
	fmt.Println("\n5. What would hold the electron's light: the ether knot of H29, fixed by the electron's mass and size")
	k := Knot(wv.R)
	fmt.Printf("   core radius 1.2ℓ = %.1f fm → ℓ = %.1f fm; energy 547.5√(κ2κ4) = m_e c² → κ2 = %.4g J/m, κ4 = %.4g J·m\n", wv.R*1e15, k.L*1e15, k.K2, k.K4)
	fmt.Printf("   its squashing tone ω = 1.416 c/ℓ: ħω = %.3f MeV = %.2f m_e c² — an excited electron at %.3f MeV\n", k.ExciteMeV, k.ExciteMeV/MeMeV, MeMeV+k.ExciteMeV)
	fmt.Println("   no excited electron has been seen: collider searches exclude excited electrons to the TeV scale (PDG 2024),")
	fmt.Println("   and g − 2 limits the electron's size to about 10⁻¹⁸ m, against the knot's 10⁻¹³ m")
	fmt.Println("   the knot also carries no charge and no spin ½; those still come from H23's assumption, not from the knot")

	// Part 6
	fmt.Println("\n6. The spectrum from the frequency response of the whole network (spectrum.go)")
	hyd, err := RunSpectrum(1, 8000)
	if err != nil {
		fmt.Println("   error:", err)
		return
	}
	light, err := RunSpectrum(wv.Q, 8000)
	if err != nil {
		fmt.Println("   error:", err)
		return
	}
	for _, sp := range []Spectrum{hyd, light} {
		fmt.Printf("   Z = %.6f: AC sweep of %d frequencies over %d nodes; poles pinned, h = 0.1 and 0.05 a0, extrapolated h → 0\n", sp.Z, sp.SweepPts, sp.Nodes)
		fmt.Printf("     %3s %15s %15s %15s %15s %11s\n", "n", "pole h=0.1", "matrix h=0.1", "pole h→0", "exact −Z²/2n²", "rel. diff")
		for k := 0; k < NLevels; k++ {
			n := float64(k + 1)
			ex := -sp.Z * sp.Z / (2 * n * n)
			fmt.Printf("     %3d %15.10f %15.10f %15.10f %15.10f %11.1e\n", k+1, sp.Circuit[k], sp.Matrix[k], sp.Extrap[k], ex, sp.Extrap[k]/ex-1)
		}
	}
	fmt.Printf("   %-15s %12s %14s %14s %16s\n", "line", "measured nm", "q = e (poles)", "q = 0.910 e", "q = e − measured")
	ll := light.Lines()
	for i, l := range hyd.Lines() {
		fmt.Printf("   %-15s %12.2f %14.4f %14.4f %+13.4f nm\n", l.Name, l.MeasuredNm, l.Circuit, ll[i].Circuit, l.Circuit-l.MeasuredNm)
	}
	fmt.Printf("   1S–2S from the poles: q = e %.6f THz, q = 0.910 e %.6f THz; measured %.6f THz (q = e: %.1f ppm low,\n",
		hyd.OneSTwoS()/1e12, light.OneSTwoS()/1e12, OneSTwoSHz/1e12, 1e6*(1-hyd.OneSTwoS()/OneSTwoSHz))
	fmt.Println("   of which ~9 ppm is fine structure + Lamb shift, the rest the extrapolation)")
	fmt.Println("   The frequency response carries the spectrum: each pole is a level, and the poles' differences are the lines.")
	fmt.Println("   With q = e they fall on measured hydrogen to the 0.01 nm listed; with the toroidal-photon charge every line is 1/Z² = 1.207× long.")

	// Part 7
	fmt.Println("\n7. Hydrogen under pressure: the network's outer short circuit as a wall moved in (pressure.go)")
	fmt.Println("   Russell's \"high pressure, short waves\" read as the confined atom (Michels, de Boer & Bijl, Physica 4, 981, 1937)")
	fmt.Printf("   %8s %11s %11s %11s %12s %12s\n", "wall (a0)", "E_1s (Ha)", "E_2s (Ha)", "E_2p (Ha)", "Lyman α (nm)", "pressure")
	for _, rc := range Walls {
		c := Compress(rc)
		p := fmt.Sprintf("%.3g GPa", c.PressureGPa)
		if rc >= 12 {
			p = "≈ 0"
		}
		fmt.Printf("   %8.3f %+11.6f %+11.6f %+11.6f %12.2f %12s\n", c.Rc, c.E1s, c.E2s, c.E2p, c.LymanNm, p)
	}
	fmt.Println("   compression raises every level and shortens every line; the free atom's 2s = 2p splits (2p stays lower);")
	fmt.Println("   at r_c = 2 a0 the ground state is exactly −1/8 Ha, and at 1.835 a0 it reaches zero: the atom is ionized by pressure.")
	fmt.Println("   For scale: the pressure at the centre of the Earth is about 364 GPa (PREM), close to the r_c = 2 a0 row.")

	// Part 8
	fmt.Println("\n8. The two-lobed state: 2p_z, two spheres meeting at a nodal plane (the network with l = 1)")
	a2, b2 := NewLadderL(1, 1, 0.02, 60), NewLadderL(1, 1, 0.01, 60)
	e2p := (4*b2.Eigen(0) - a2.Eigen(0)) / 3
	tw2, err := TwinL(1, 1, 0.1, 40)
	if err != nil {
		fmt.Println("   error:", err)
		return
	}
	fmt.Printf("   network, centrifugal term l(l+1)/2r² in the shunts: E_2p = %.9f Ha (h → 0), exact −1/8; equal to E_2s in the free atom\n", e2p)
	fmt.Printf("   kronspice twin (h = 0.1 a0, %d nodes): %.12f Ha against the matrix %.12f (diff %.1e)\n", tw2.Nodes, tw2.Circuit, tw2.Matrix, tw2.Circuit-tw2.Matrix)
	fmt.Println("   density r²e^(−r)cos²θ: two lobes on the axis, zero on the central plane; 1s → 2p is Lyman α, the line Fig. 2 matches")

	if *svgDir != "" {
		if err := WriteFigures(*svgDir, hyd, light); err != nil {
			fmt.Println("figures:", err)
			return
		}
		if err := WriteScene(filepath.Join(*svgDir, "scene3d.json")); err != nil {
			fmt.Println("scene:", err)
			return
		}
		fmt.Println("\nfigures written to", *svgDir)
	}

	fmt.Println("\nSummary: the light atom is a Coulomb atom whose charge is q_e = 0.910 e. That leaves it charged (0.090 e,")
	fmt.Println("excluded at 10⁻²¹ e), shifts every line by Z² = 0.829, and lowers the 21 cm line by Z³·(q_e/e) ≈ 0.69 times the proton factor.")
	fmt.Println("With q set to e, the levels are Bohr's; the light-ring proton then misses 21 cm by its moment (1 μN against 2.793)")
	fmt.Println("and 1S–2S by its size; the quark-ring proton, fitted to the moment, gets 21 cm and misses 1S–2S because it is smaller")
	fmt.Println("(0.678 fm against 0.841): 0.34 MHz, against a measurement good to 10 Hz.")
}
