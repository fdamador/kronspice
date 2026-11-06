// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Part 4: the spectral resonances.
//
// The atom's levels are resonances of Kron's 1945 network of the radial
// Schrödinger equation (Phys. Rev. 67, 39; built as in
// examples/light_hydrogen, network.go): u = rψ on a grid r_i = i·h, series
// branches for the kinetic energy, shunt branches for the potential, the
// energy as ω². kronspice solves it in its dual form (capacitors and
// inductors exchanged, so AC starts from a clean DC point) and each level is
// where the driving-point reactance changes sign. For the p states the shunt
// potential carries the spin–orbit term (g_s − 1)(α²/2)⟨L·S⟩/r³ (the
// Thomas-corrected Pauli term, Bethe & Salpeter §12), so the networks for
// j = ½ and j = 3⁄2 resonate at different frequencies and their difference
// is the fine-structure splitting. H₂'s vibration is the resonance of a
// two-mass LC circuit (mass ↔ capacitor, spring ↔ inductor) with the spring
// constant from the Heitler–London curve.

// Physical constants (CODATA 2018).
const (
	Alpha   = 7.2973525693e-3
	HaHz    = 6.579683920502e15 // hartree in Hz
	HaCm    = 219474.6313632    // hartree in cm⁻¹
	CLight  = 299792458.0
	MpMe    = 1836.15267343
	MdMe    = 3670.48296788
	MuN     = 5.0507837461e-27 // J/T
	GP      = 5.5856946893     // proton g-factor
	HPlanck = 6.62607015e-34
	MuRedH  = 1 / (1 + 1/MpMe) // hydrogen reduced mass / m_e
)

// Ladder is the radial network as a tridiagonal matrix.
type Ladder struct {
	H float64
	D []float64
	O float64
}

// NewLadderV builds the network for the radial potential V(r) (including
// the centrifugal term), spacing h, edge rmax (short circuits at both ends).
func NewLadderV(v func(r float64) float64, h, rmax float64) Ladder {
	n := int(math.Round(rmax/h)) - 1
	l := Ladder{H: h, D: make([]float64, n), O: -1 / (2 * h * h)}
	for i := range l.D {
		l.D[i] = 1/(h*h) + v(float64(i+1)*h)
	}
	return l
}

// Potential is the radial potential for orbital l and, for l > 0, total j:
// −1/r + l(l+1)/2r² + (g_s − 1)(α²/2)⟨L·S⟩/r³ with
// ⟨L·S⟩ = [j(j+1) − l(l+1) − ¾]/2. so = false leaves the spin–orbit out.
func Potential(l int, j float64, so bool) func(r float64) float64 {
	L := float64(l)
	ls := (j*(j+1) - L*(L+1) - 0.75) / 2
	return func(r float64) float64 {
		v := -1/r + L*(L+1)/(2*r*r)
		if so && l > 0 {
			v += (GS - 1) * Alpha * Alpha / 2 * ls / (r * r * r)
		}
		return v
	}
}

func (l Ladder) below(x float64) int {
	n, d := 0, 1.0
	for i, di := range l.D {
		if i == 0 {
			d = di - x
		} else {
			if d == 0 {
				d = 1e-300
			}
			d = di - x - l.O*l.O/d
		}
		if d < 0 {
			n++
		}
	}
	return n
}

// Eigen is the k-th eigenvalue by Sturm bisection (the matrix twin).
func (l Ladder) Eigen(k int) float64 {
	lo, hi := -10.0, 10.0
	for i := 0; i < 200 && hi-lo > 1e-15; i++ {
		m := 0.5 * (lo + hi)
		if l.below(m) > k {
			hi = m
		} else {
			lo = m
		}
	}
	return 0.5 * (lo + hi)
}

// network is the dual-form circuit: series capacitors 1/(2h²) between
// neighbours (and to ground at both ends), shunt capacitors V_i + V0, unit
// shunt inductors, a unit AC current into the node nearest r = rd. Its
// admittance times jω is ω²(H + V0) − 1: a level E is the resonance at
// ω² = 1/(E + V0).
func (l Ladder) network(rd float64) (*circuit.Circuit, string, float64, error) {
	n := len(l.D)
	v0 := 0.0
	for _, d := range l.D {
		v0 = math.Max(v0, -(d - 1/(l.H*l.H)))
	}
	v0++
	b := circuit.NewBuilder("Kron radial network (dual form)")
	node := func(i int) string {
		if i < 0 || i >= n {
			return "0"
		}
		return "n" + strconv.Itoa(i)
	}
	cs := 1 / (2 * l.H * l.H)
	for i := 0; i < n; i++ {
		b.Inductor("LB"+strconv.Itoa(i), node(i), "0", 1)
		b.Capacitor("CV"+strconv.Itoa(i), node(i), "0", (l.D[i]-1/(l.H*l.H))+v0)
		b.Capacitor("CS"+strconv.Itoa(i), node(i), node(i+1), cs)
	}
	b.Capacitor("CS0G", node(0), "0", cs)
	drive := node(int(rd/l.H) - 1)
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "0", drive)
	ckt, err := b.BuildTransient()
	return ckt, drive, v0, err
}

// Resonance finds the network's lowest level E in kronspice: a scan up in E
// from lo for the reactance's − → + change, then bisection. Returns the
// level [Ha] and the number of AC solves.
func (l Ladder) Resonance(lo, rd float64) (float64, int, error) {
	ckt, drive, v0, err := l.network(rd)
	if err != nil {
		return 0, 0, err
	}
	solves := 0
	im := func(e float64) (float64, error) {
		solves++
		f := 1 / math.Sqrt(e+v0) / (2 * math.Pi)
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		v, ok := res.NodeVoltage(drive)
		if !ok {
			return 0, fmt.Errorf("no node %s", drive)
		}
		return imag(v[0]) / math.Max(cmplx.Abs(v[0]), 1e-300), nil
	}
	s, err := im(lo)
	if err != nil {
		return 0, solves, err
	}
	if s >= 0 {
		return 0, solves, fmt.Errorf("scan starts above the level")
	}
	hi := 0.0
	const steps = 120
	found := false
	for k := 1; k <= steps; k++ {
		e := lo * (1 - float64(k)/steps)
		s, err := im(e)
		if err != nil {
			return 0, solves, err
		}
		if s > 0 {
			hi, found = e, true
			break
		}
		lo = e
	}
	if !found {
		return 0, solves, fmt.Errorf("no resonance found")
	}
	for i := 0; i < 90 && hi-lo > 1e-15; i++ {
		m := 0.5 * (lo + hi)
		s, err := im(m)
		if err != nil {
			return 0, solves, err
		}
		if s < 0 {
			lo = m
		} else {
			hi = m
		}
	}
	return 0.5 * (lo + hi), solves, nil
}

// Level is one hydrogen level found by kronspice.
type Level struct {
	Name   string
	N, L   int
	J      float64
	SO     bool
	Coarse float64 // E at h
	Fine   float64 // E at h/2
	Rich   float64 // Richardson (4·fine − coarse)/3, O(h²) removed
	Matrix float64 // Sturm twin at h/2
	Solves int
}

// FindLevel resonates the network at h and h/2 (Richardson) for the lowest
// level of orbital l (the circular state n = l + 1).
func FindLevel(name string, l int, j float64, so bool, h, rmax float64) (Level, error) {
	lv := Level{Name: name, N: l + 1, L: l, J: j, SO: so}
	n := float64(l + 1)
	lo := -0.7 / (n * n)
	rd := n * n // drive near the orbit's radius
	var es [2]float64
	for k, hh := range []float64{h, h / 2} {
		e, s, err := NewLadderV(Potential(l, j, so), hh, rmax).Resonance(lo, rd)
		if err != nil {
			return lv, fmt.Errorf("%s: %v", name, err)
		}
		es[k] = e
		lv.Solves += s
	}
	lv.Coarse, lv.Fine = es[0], es[1]
	lv.Rich = (4*es[1] - es[0]) / 3
	lv.Matrix = NewLadderV(Potential(l, j, so), h/2, rmax).Eigen(0)
	return lv, nil
}

// Line is a computed spectral line against its measured value.
type Line struct {
	Name, Where  string
	Hz, Measured float64 // frequency [Hz]
	How          string  // "kronspice resonance" or "model"
}

// VacNm is the vacuum wavelength in nm.
func VacNm(hz float64) float64 { return CLight / hz * 1e9 }

// H2Vib is the vibration of a diatomic as a kronspice LC resonance:
// two masses m1, m2 (capacitors, in m_e) joined by the spring k (inductor
// 1/k, in Eh/a0²), driven by an AC current (a force) at mass 1 and read at
// mass 2. The transfer response has one pole, the internal mode
// ω = √(k/μ); the driving-point response would also cross zero at the
// antiresonance √(k/m2), where mass 1 stands still.
// Returns ω_e in cm⁻¹ and the number of AC solves.
func H2Vib(k, m1, m2 float64) (float64, int, error) {
	b := circuit.NewBuilder("diatomic vibration: masses as capacitors, spring as inductor")
	b.Capacitor("CM1", "a", "0", m1)
	b.Capacitor("CM2", "b", "0", m2)
	b.Inductor("LK", "a", "b", 1/k)
	// DC paths for the operating point: leaks 10⁶ times the spring's
	// reactance, whose damping (Q ≈ 10⁶) moves the resonance by ~10⁻¹², below the bisection
	// tolerance
	b.Resistor("RA", "a", "0", 1e6*math.Sqrt(m1/k))
	b.Resistor("RB", "b", "0", 1e6*math.Sqrt(m2/k))
	b.Add(circuit.Element{Name: "IF", Kind: circuit.CurrentSource, ACMag: 1}, "0", "a")
	ckt, err := b.BuildTransient()
	if err != nil {
		return 0, 0, err
	}
	solves := 0
	im := func(w float64) (float64, error) {
		solves++
		f := w / (2 * math.Pi)
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return 0, err
		}
		v, _ := res.NodeVoltage("b")
		return imag(v[0]), nil
	}
	// the transfer reactance changes sign only at the pole; scan up from
	// 0.2 of the guess
	mu := m1 * m2 / (m1 + m2)
	guess := math.Sqrt(k / mu)
	lo, hi := 0.2*guess, 0.0
	s0, err := im(lo)
	if err != nil {
		return 0, solves, err
	}
	for w := lo; w < 5*guess; w *= 1.01 {
		s, err := im(w)
		if err != nil {
			return 0, solves, err
		}
		if s*s0 < 0 {
			hi = w
			break
		}
		lo = w
	}
	if hi == 0 {
		return 0, solves, fmt.Errorf("no resonance")
	}
	for i := 0; i < 80 && hi-lo > 1e-11*hi; i++ {
		m := 0.5 * (lo + hi)
		s, err := im(m)
		if err != nil {
			return 0, solves, err
		}
		if s*s0 > 0 {
			lo = m
		} else {
			hi = m
		}
	}
	w := 0.5 * (lo + hi)
	return w * HaCm, solves, nil
}

// HLCurvature is d²E/dR² of the Heitler–London singlet at its minimum
// [Eh/a0²].
func HLCurvature() (re, k float64) {
	re, _ = HLMinimum()
	const d = 1e-3
	k = (HLEnergy(re+d, 1) - 2*HLEnergy(re, 1) + HLEnergy(re-d, 1)) / (d * d)
	return
}

// RotB is the rotational constant 1/(2μR²) in cm⁻¹ (μ in m_e, R in a0).
func RotB(mu, r float64) float64 { return HaCm / (2 * mu * r * r) }

// HyperfineHz is the 1s hyperfine (21 cm) frequency from the Fermi contact
// interaction, ΔE = (2μ0/3) g_s μB g_p μN |ψ(0)|², with |ψ(0)|² from the
// model's 1s spinor and the reduced-mass factor (μ/m)³ (Bethe & Salpeter
// §22).
func HyperfineHz(s State) float64 {
	u, d := s.Spinor(1e-9, 0.3)
	psi0 := (u*u + d*d) / (A0 * A0 * A0)
	return 2 * Mu0 / 3 * GS * MuB * GP * MuN * psi0 * math.Pow(MuRedH, 3) / HPlanck
}

// ESRHz is the electron spin resonance frequency per tesla, |μ_z|/(ħ m_j)
// from the model's dipole moment (μB/h × g).
func ESRHz(f AtomField) float64 {
	g := math.Abs(f.Ser.Moment(1)) / f.S.M
	return g * MuB / HPlanck
}

// Spectrum is Part 4's results.
type Spectrum struct {
	Levels     []Level
	FS, FSNoSO float64 // 2p3⁄2 − 2p½ [Hz], and with the spin–orbit off
	FSMatrix   float64
	Lines      []Line
	K, Re      float64
	Vib        [3]float64 // ω_e H₂, HD, D₂ [cm⁻¹]
	VibSolves  int
}

// Measured values used for comparison.
var (
	measLya = CLight / 121.567e-9 // NIST ASD, vacuum, fine-structure weighted
	measHa  = CLight / 656.461e-9
	measPaa = CLight / 1875.61e-9
	measBra = CLight / 4052.27e-9
	measFS  = 10969.1e6                             // 2P3⁄2 − 2P½ (e.g. Hagley & Pipkin, PRL 72, 1172, 1994)
	meas21  = 1420.405751768e6                      // Hellwig et al., IEEE Trans. IM 19, 200 (1970)
	measESR = 2.0022838 * MuB / HPlanck             // bound g_J(1s) = g_s(1 − 17.709×10⁻⁶) (Tiedemann & Robinson 1977)
	measWe  = [3]float64{4401.21, 3813.15, 3115.50} // ω_e, Huber & Herzberg (1979)
	measBe  = 60.853                                // B_e of H₂, Huber & Herzberg (1979)
	measRe  = 1.4011
)

// SpectrumPart runs Part 4 and prints it.
func SpectrumPart(as []AtomResult) (Spectrum, error) {
	var sp Spectrum
	fmt.Println("\n4. Spectral resonances found by kronspice, against measurement")
	fmt.Println("   Kron's 1945 radial network in kronspice (dual form, AC resonance by bisection), h = 0.1 and 0.05 a0, Richardson;")
	fmt.Println("   p networks carry the spin–orbit term (g_s − 1)(α²/2)⟨L·S⟩/r³; energies scaled by the reduced mass μ/m = 0.999456")
	const h, rmax = 0.1, 100
	specs := []struct {
		name string
		l    int
		j    float64
		so   bool
	}{{"1s", 0, 0.5, false}, {"2p (no spin–orbit)", 1, 0.5, false}, {"2p½", 1, 0.5, true}, {"2p3⁄2", 1, 1.5, true}, {"3d", 2, 2.5, false}, {"4f", 3, 3.5, false}, {"5g", 4, 4.5, false}}
	out := make([]Level, len(specs))
	errs := make([]error, len(specs))
	parallel(len(specs), func(i int) {
		s := specs[i]
		out[i], errs[i] = FindLevel(s.name, s.l, s.j, s.so, h, rmax)
	})
	for _, e := range errs {
		if e != nil {
			return sp, e
		}
	}
	sp.Levels = out
	fmt.Printf("   %-20s %16s %16s %16s %10s %8s\n", "level", "kronspice h", "kronspice h/2", "Richardson", "exact", "AC solves")
	for _, l := range out {
		exact := -0.5 / float64(l.N*l.N)
		fmt.Printf("   %-20s %16.10f %16.10f %16.10f %10.6f %8d\n", l.Name, l.Coarse, l.Fine, l.Rich, exact, l.Solves)
	}
	worst := 0.0
	for _, l := range out {
		worst = math.Max(worst, math.Abs(l.Fine-l.Matrix))
	}
	fmt.Printf("   twin: kronspice resonance against the network matrix (Sturm), h/2: worst %.1e Ha\n", worst)
	E := func(i int) float64 { return out[i].Rich * MuRedH * HaHz }
	sp.FS = (out[3].Rich - out[2].Rich) * MuRedH * HaHz
	sp.FSNoSO = 0
	sp.FSMatrix = (out[3].Matrix - out[2].Matrix) * MuRedH * HaHz
	sp.Lines = []Line{
		{"Lyman α (2p → 1s)", "solar and interstellar ultraviolet; the Lyman-α forest of quasar spectra", E(1) - E(0), measLya, "kronspice resonance"},
		{"Balmer α, Hα (3d → 2p)", "the red glow of emission nebulae and the solar chromosphere", E(4) - E(1), measHa, "kronspice resonance"},
		{"Paschen α (4f → 3d)", "infrared hydrogen line seen through dust in star-forming regions", E(5) - E(4), measPaa, "kronspice resonance"},
		{"Brackett α (5g → 4f)", "infrared line of ionized gas round young stars", E(6) - E(5), measBra, "kronspice resonance"},
		{"fine structure 2p3⁄2 − 2p½", "the microwave interval of the Lamb–Retherford experiments", sp.FS, measFS, "kronspice resonance"},
	}
	hf := HyperfineHz(as[0].F.S)
	esr := ESRHz(as[0].F)
	sp.Lines = append(sp.Lines,
		Line{"hyperfine 1s, 21 cm (F = 1 → 0)", "radio astronomy of galactic hydrogen; the hydrogen-maser clock", hf, meas21, "model (Fermi contact, ψ(0) of Part 1)"},
		Line{"electron spin resonance of H, per tesla", "ESR of hydrogen atoms in solids and gases", esr, measESR, "model (dipole of Part 1)"})

	re, k := HLCurvature()
	sp.K, sp.Re = k, re
	masses := [3][2]float64{{MpMe, MpMe}, {MpMe, MdMe}, {MdMe, MdMe}}
	for i, m := range masses {
		w, n, err := H2Vib(k, m[0], m[1])
		if err != nil {
			return sp, err
		}
		sp.Vib[i] = w
		sp.VibSolves += n
	}
	fsAt := func(a, b float64) float64 { return (a - b) * MuRedH * HaHz / 1e6 }
	fmt.Printf("   fine structure 2p3⁄2 − 2p½ from the two networks' resonances: h %.1f MHz, h/2 %.1f MHz, Richardson %.1f MHz\n",
		fsAt(out[3].Coarse, out[2].Coarse), fsAt(out[3].Fine, out[2].Fine), sp.FS/1e6)
	fmt.Println("   (with the spin–orbit term off, 2p½ and 2p3⁄2 are one network and one resonance)")
	fmt.Printf("   %-36s %16s %16s %12s  %s\n", "line", "found", "measured", "difference", "how")
	for _, l := range sp.Lines {
		fmt.Printf("   %-36s %16s %16s %+11.4f%%  %s\n", l.Name, fmtFreq(l.Hz), fmtFreq(l.Measured), 100*(l.Hz/l.Measured-1), l.How)
	}
	fmt.Printf("   H₂ vibration as a kronspice LC circuit (masses = capacitors, Heitler–London spring k = %.4f Eh/a0² at %.4f a0 = inductor 1/k):\n", k, re)
	names := []string{"H₂", "HD", "D₂"}
	for i := range names {
		fmt.Printf("     %-3s ω_e %8.1f cm⁻¹ against measured %8.2f (%+.1f%%)\n", names[i], sp.Vib[i], measWe[i], 100*(sp.Vib[i]/measWe[i]-1))
	}
	fmt.Printf("     isotope ratios, which do not depend on the spring: HD/H₂ %.4f (measured %.4f), D₂/H₂ %.4f (measured %.4f)\n",
		sp.Vib[1]/sp.Vib[0], measWe[1]/measWe[0], sp.Vib[2]/sp.Vib[0], measWe[2]/measWe[0])
	fmt.Printf("     rotation B = 1/(2μR²): %.2f cm⁻¹ at the Heitler–London R_e, %.2f at the measured R_e, against measured B_e %.3f\n",
		RotB(MpMe/2, re), RotB(MpMe/2, measRe), measBe)
	return sp, nil
}

func fmtFreq(hz float64) string {
	switch {
	case hz > 1e13:
		return fmt.Sprintf("%.4f nm", VacNm(hz))
	case hz > 5e9:
		return fmt.Sprintf("%.2f MHz", hz/1e6)
	default:
		return fmt.Sprintf("%.4f MHz", hz/1e6)
	}
}

// Response is Im(1/V) of a network's driving point at the energies es
// (h/2 grid): it crosses zero, linearly, at each resonance.
func (l Ladder) Response(es []float64, rd float64) ([]float64, error) {
	ckt, drive, v0, err := l.network(rd)
	if err != nil {
		return nil, err
	}
	out := make([]float64, len(es))
	for i, e := range es {
		f := 1 / math.Sqrt(e+v0) / (2 * math.Pi)
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return nil, err
		}
		v, _ := res.NodeVoltage(drive)
		out[i] = imag(1 / v[0])
	}
	return out, nil
}

// VibResponse is |V| at mass 2 of the diatomic circuit over wavenumbers
// ws [cm⁻¹], with a light damping (Q = 200) so the peak has a width.
func VibResponse(k, m1, m2 float64, ws []float64) ([]float64, error) {
	b := circuit.NewBuilder("diatomic vibration, damped for display")
	b.Capacitor("CM1", "a", "0", m1)
	b.Capacitor("CM2", "b", "0", m2)
	b.Inductor("LK", "a", "x", 1/k)
	mu := m1 * m2 / (m1 + m2)
	w0 := math.Sqrt(k / mu)
	b.Resistor("RK", "x", "b", w0/k/200) // series loss: Q = ω0 L/R = 200
	b.Resistor("RA", "a", "0", 1e6*math.Sqrt(m1/k))
	b.Resistor("RB", "b", "0", 1e6*math.Sqrt(m2/k))
	b.Add(circuit.Element{Name: "IF", Kind: circuit.CurrentSource, ACMag: 1}, "0", "a")
	ckt, err := b.BuildTransient()
	if err != nil {
		return nil, err
	}
	out := make([]float64, len(ws))
	for i, w := range ws {
		f := w / HaCm / (2 * math.Pi)
		res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return nil, err
		}
		v, _ := res.NodeVoltage("b")
		out[i] = cmplx.Abs(v[0])
	}
	return out, nil
}
