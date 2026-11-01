// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command pais_cavity puts the annular resonant cavity of U.S. Patent
// 10,144,532 B2 (S. C. Pais, "Craft using an inertial mass reduction
// device", assignee U.S. Navy, filed 28 Apr 2016, issued 4 Dec 2018) through
// a Kron circuit analysis and the energy ledger, beside the ring towers
// (examples/ring_towers) and the double hull (examples/double_hull).
//
// What the patent states (claims 1–4; description): an inner and an outer
// resonant-cavity wall forming an annular cavity, filled with a noble gas
// (xenon), driven by microwave emitters (300 MHz–300 GHz, Q 10⁴–10⁹); the
// outer wall electrically charged and "vibrating in an accelerated mode"
// (10⁹–10¹⁸ Hz) and optionally rotated (example: 2 m radius, 30,000 rpm,
// σ = 50,000 C/m²), to "create a local polarized vacuum" and reduce the
// craft's inertial mass; sections of the shroud driven separately to steer.
// Its energy-flux figure is S = f_G·(σ²/ε0)·[R_r·ω + R_v·ν + v·R], giving
// ≈ 10²⁴ W/m² for the rotation example and ≈ 10³³ W/m² at σ = 1 C/m²,
// 10¹¹ Hz. No wall gap, height, voltage or drive power is given, so the
// geometry below is declared, not the patent's.
//
// Hypothesis H18 (CLAUDE.md §1): inertial mass reduction from a polarized
// vacuum. Standard prediction: none — a closed cavity's fields exert no net
// time-averaged force on the craft (momentum conservation), mass is
// unchanged, and the only way to push is to emit momentum, at most P/c. The
// gyroscope weight reduction the patent cites (H. Hayasaka, S. Takeuchi,
// PRL 63, 2701, 1989) was not reproduced (J. E. Faller et al., PRL 64,
// 825, 1990; T. J. Quinn, A. Picard, Nature 343, 732, 1990; J. M. Nitschke,
// P. A. Wilmarth, PRL 64, 2115, 1990).
//
// Circuit: the annulus is a closed parallel-plate TEM line around the
// craft (as in double_hull; Paul §4.2): L′ = μ0·d/h, C′ = ε0·h/d,
// R′ = 2·Rs/h with Rs the walls' surface resistance; N sections close the
// ring, the inner wall is the common bus. Emitters are current
// injections. Kron analysis: Happ's two-level diakoptics with a single
// common bus (kron.TwoLevel; Happ 1971 Ch. XII §4, Fig. 12.12): the ring is
// torn into S sectors (each a radially hinged subdivision) joined by S link
// branches, so the interconnection system has order S instead of the
// 2N unknowns of MNA. The same circuit through kronspice's MNA AC analysis
// is the twin.
//
//	go run ./examples/pais_cavity
package main

import (
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/spacetime"
)

// Physical constants (CODATA 2018).
const (
	mu0  = 1.25663706212e-6
	eps0 = 8.8541878128e-12
	c0   = 299792458.0
	qe   = 1.602176634e-19
	me   = 9.1093837015e-31
	hbar = 1.054571817e-34
	kB   = 1.380649e-23
	hP   = 6.62607015e-34
)

// Geometry: declared (the patent gives only the 2 m example radius).
type Geometry struct {
	R float64 // mean ring radius [m] (patent example disc radius 2 m)
	D float64 // wall gap [m]
	H float64 // wall height [m]
	N int     // ring sections
	S int     // Kron sectors (subdivisions); N must be a multiple of S
	// Rs is the wall surface resistance [Ω]; 0 means copper at the drive
	// frequency (ρ = 1.724e−8 Ω·m, 100 % IACS).
	Rs   float64
	Mode int // ring mode number n (wavelengths around the ring)
}

// Default: 2 m radius, 5 cm gap, 0.5 m walls, mode 16 (≈ 382 MHz, inside
// the patent's 300 MHz–300 GHz band), 640 sections (40 per wavelength).
var Default = Geometry{R: 2, D: 0.05, H: 0.5, N: 640, S: 16, Mode: 16}

// Line returns the per-section L, C and the circumference.
func (g Geometry) line() (dl, l, cc float64) {
	dl = 2 * math.Pi * g.R / float64(g.N)
	return dl, mu0 * g.D / g.H * dl, eps0 * g.H / g.D * dl
}

// Freq is the exact resonance of mode n of the discrete ring,
// ω = (2/√(LC))·sin(πn/N) (lossless lumped ladder).
func (g Geometry) Freq() float64 {
	_, l, c := g.line()
	return 2 / math.Sqrt(l*c) * math.Sin(math.Pi*float64(g.Mode)/float64(g.N)) / (2 * math.Pi)
}

// SurfaceR is the wall surface resistance at f.
func (g Geometry) SurfaceR(f float64) float64 {
	if g.Rs > 0 {
		return g.Rs
	}
	return math.Sqrt(math.Pi * f * mu0 * 1.724e-8)
}

// Emitter is a current injection at a ring node.
type Emitter struct {
	Node  int
	Amp   float64 // peak current [A] before power scaling
	Phase float64 // [deg]
}

// Result of one drive case.
type Result struct {
	F, Q         float64
	V            []complex128 // node voltages (peak phasors)
	I            []complex128 // series currents, node k → k+1
	Pin, Ploss   float64      // emitter power, wall loss [W]
	W            float64      // time-averaged stored energy [J]
	Quadrant     [4]float64   // fraction of W in the quadrants centred on 0°, 90°, 180°, 270° (node 0 at 0°)
	Vmax, Emax   float64      // peak gap voltage, field
	SigmaRF      float64      // peak RF surface charge on the walls ε0·Emax [C/m²]
	Interconnect int          // order of the Kron interconnection system
	Scale        float64      // factor applied to the emitter currents
}

// Solve drives the ring at f with the emitters by Happ's two-level
// diakoptics. With p > 0 every phasor is scaled so the emitters deliver p
// [W]; with p ≤ 0 the emitter currents are used as given.
func Solve(g Geometry, f float64, em []Emitter, p float64) (*Result, error) {
	if g.N%g.S != 0 {
		return nil, fmt.Errorf("N = %d is not a multiple of S = %d", g.N, g.S)
	}
	dl, l, cs := g.line()
	w := 2 * math.Pi * f
	z := complex(2*g.SurfaceR(f)/g.H*dl, w*l)
	ysh := complex(0, w*cs)
	m := g.N / g.S
	y := make([][][]complex128, g.S)
	for k := range y {
		yk := make([][]complex128, m)
		for i := range yk {
			yk[i] = make([]complex128, m)
			yk[i][i] = ysh
		}
		for i := 0; i+1 < m; i++ { // series branches inside the sector
			a := 1 / z
			yk[i][i] += a
			yk[i+1][i+1] += a
			yk[i][i+1] -= a
			yk[i+1][i] -= a
		}
		y[k] = yk
	}
	links := make([]kron.Link, g.S)
	for k := range links { // last node of sector k → first of k+1
		links[k] = kron.Link{From: kron.BusRef{Sub: k, Bus: m - 1}, To: kron.BusRef{Sub: (k + 1) % g.S, Bus: 0}, Z: z}
	}
	tl, err := kron.NewTwoLevel(y, links)
	if err != nil {
		return nil, err
	}
	inj := make([][]complex128, g.S)
	for k := range inj {
		inj[k] = make([]complex128, m)
	}
	for _, e := range em {
		n := ((e.Node % g.N) + g.N) % g.N
		inj[n/m][n%m] += cmplx.Rect(e.Amp, e.Phase*math.Pi/180)
	}
	s, err := tl.Solve(inj)
	if err != nil {
		return nil, err
	}
	r := &Result{F: f, V: make([]complex128, g.N), I: make([]complex128, g.N), Interconnect: g.S}
	for k := range s.V {
		copy(r.V[k*m:], s.V[k])
	}
	for k := 0; k < g.N; k++ {
		if (k+1)%m == 0 {
			r.I[k] = s.Ic[k/m] // link current
		} else {
			r.I[k] = (r.V[k] - r.V[k+1]) / z
		}
	}
	var pin float64
	for _, e := range em {
		n := ((e.Node % g.N) + g.N) % g.N
		pin += 0.5 * real(r.V[n]*cmplx.Conj(cmplx.Rect(e.Amp, e.Phase*math.Pi/180)))
	}
	scale := 1.0
	if p > 0 {
		scale = math.Sqrt(p / pin) // linear: scale every phasor to input power p
	}
	for k := range r.V {
		r.V[k] *= complex(scale, 0)
		r.I[k] *= complex(scale, 0)
	}
	r.Pin = pin * scale * scale
	r.Scale = scale
	for k := 0; k < g.N; k++ {
		r.Ploss += 0.5 * real(z) * real(r.I[k]*cmplx.Conj(r.I[k]))
		wk := 0.25*cs*real(r.V[k]*cmplx.Conj(r.V[k])) + 0.25*l*real(r.I[k]*cmplx.Conj(r.I[k]))
		r.W += wk
		r.Quadrant[(k+g.N/8)%g.N*4/g.N] += wk
		r.Vmax = math.Max(r.Vmax, cmplx.Abs(r.V[k]))
	}
	for i := range r.Quadrant {
		r.Quadrant[i] /= r.W
	}
	r.Q = w * r.W / r.Ploss
	r.Emax = r.Vmax / g.D
	r.SigmaRF = eps0 * r.Emax
	return r, nil
}

// Netlist builds the same ring for kronspice's MNA AC analysis (the twin).
func Netlist(g Geometry, f float64, em []Emitter) string {
	dl, l, cs := g.line()
	rsec := 2 * g.SurfaceR(f) / g.H * dl
	var b strings.Builder
	b.WriteString("pais annular cavity ring\n")
	for k := 0; k < g.N; k++ {
		fmt.Fprintf(&b, "R%d n%d m%d %.17g\nL%d m%d n%d %.17g\nC%d n%d 0 %.17g\n", k, k, k, rsec, k, k, (k+1)%g.N, l, k, k, cs)
	}
	for i, e := range em {
		fmt.Fprintf(&b, "IE%d 0 n%d AC %.17g %.17g\n", i, ((e.Node%g.N)+g.N)%g.N, e.Amp, e.Phase)
	}
	// The operating point that precedes AC needs a DC path to ground: a
	// 1 MH choke, |Y| ≈ 4e−16 S at the drive frequency, against a section
	// admittance of ~1e−2 S.
	b.WriteString("LDC n0 0 1e6\n.end\n")
	return b.String()
}

// TwinMNA solves Netlist at f and returns the node voltages (unscaled).
func TwinMNA(g Geometry, f float64, em []Emitter) ([]complex128, error) {
	d, err := netlist.ParseString(Netlist(g, f, em))
	if err != nil {
		return nil, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil, err
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		return nil, err
	}
	v := make([]complex128, g.N)
	for k := range v {
		s, ok := r.NodeVoltage(fmt.Sprintf("n%d", k))
		if !ok {
			return nil, fmt.Errorf("node n%d missing", k)
		}
		v[k] = s[0]
	}
	return v, nil
}

// residual returns max_k |Y·V − I|_k / max_k (|Y|·|V|)_k for the ring's
// nodal equations, the backward error of a solution.
func residual(g Geometry, f float64, em []Emitter, v []complex128, scale float64) float64 {
	dl, l, cs := g.line()
	w := 2 * math.Pi * f
	z := complex(2*g.SurfaceR(f)/g.H*dl, w*l)
	ysh := complex(0, w*cs)
	inj := make([]complex128, g.N)
	for _, e := range em {
		inj[e.Node] += cmplx.Rect(e.Amp*scale, e.Phase*math.Pi/180)
	}
	var worst, den float64
	for k := 0; k < g.N; k++ {
		a, b := v[(k+g.N-1)%g.N], v[(k+1)%g.N]
		r := ysh*v[k] + (v[k]-a)/z + (v[k]-b)/z - inj[k]
		worst = math.Max(worst, cmplx.Abs(r))
		den = math.Max(den, cmplx.Abs(ysh*v[k])+2*cmplx.Abs(v[k])/cmplx.Abs(z)+cmplx.Abs(a)/cmplx.Abs(z)+cmplx.Abs(b)/cmplx.Abs(z))
	}
	return worst / den
}

// Limits are the physical ceilings that the patent's figures meet.
type Limits struct {
	Sigma        float64 // the surface charge examined [C/m²]
	E            float64 // its field σ/ε0 [V/m]
	Pressure     float64 // electrostatic pressure on the wall σ²/(2ε0) [Pa]
	CopperUTS    float64 // annealed copper tensile strength [Pa]
	FNCurrent    float64 // Fowler–Nordheim emission at E, φ = 4.5 eV [A/m²]
	DrainTime    float64 // σ / J_FN [s]: time to lose the charge by field emission
	Schwinger    float64 // m²c³/(eħ) [V/m]
	Kilpatrick   float64 // Kilpatrick RF limit at the drive frequency [V/m]
	PowerNeeded  float64 // RF drive power to reach σ on the cavity walls at this Q [W]
	SpeedForFlux float64 // R_v·ν needed for the patent's 10³³ W/m² at σ = 1 C/m² [m/s]
}

// fowlerNordheim is the elementary FN current density (no image-charge
// correction): J = (a·E²/φ)·exp(−b·φ^{3/2}/E), a = 1.541434e−6 A·eV·V⁻²,
// b = 6.830890e9 eV^{−3/2}·V/m (R. G. Forbes, J. Vac. Sci. Technol. B 26,
// 788, 2008, eq. 1).
func fowlerNordheim(e, phi float64) float64 {
	return 1.541434e-6 * e * e / phi * math.Exp(-6.830890e9*math.Pow(phi, 1.5)/e)
}

// kilpatrick solves f[MHz] = 1.64·E²·exp(−8.5/E), E in MV/m (T. J.
// Boyd's 1982 fit of Kilpatrick 1957; as in double_hull).
func kilpatrick(f float64) float64 {
	fm := f / 1e6
	lo, hi := 0.01, 1000.0
	for i := 0; i < 100; i++ {
		e := 0.5 * (lo + hi)
		if 1.64*e*e*math.Exp(-8.5/e) < fm {
			lo = e
		} else {
			hi = e
		}
	}
	return 0.5 * (lo + hi) * 1e6
}

// LimitsFor evaluates surface charge sigma against the ceilings, using a
// solved case r for the drive-power scaling (P ∝ E² in a linear cavity).
func LimitsFor(sigma float64, r *Result) Limits {
	e := sigma / eps0
	j := fowlerNordheim(e, 4.5)
	l := Limits{Sigma: sigma, E: e, Pressure: sigma * sigma / (2 * eps0), CopperUTS: 210e6,
		FNCurrent: j, DrainTime: sigma / j, Schwinger: me * me * c0 * c0 * c0 / (qe * hbar),
		Kilpatrick: kilpatrick(r.F), SpeedForFlux: 1e33 * eps0}
	l.PowerNeeded = r.Pin * (e / r.Emax) * (e / r.Emax)
	return l
}

// PatentFlux is the patent's S_max = (σ²/ε0)·(R·ω) for rotation (f_G = 1).
func PatentFlux(sigma, r, omega float64) float64 { return sigma * sigma / eps0 * r * omega }

// Debye frequencies ν_D = k_B·θ_D/h: the highest lattice vibration a solid
// supports (θ_D: Cu 343 K, Al 428 K; C. Kittel, Introduction to Solid State
// Physics, 8th ed., Table 5.1).
func Debye(theta float64) float64 { return kB * theta / hP }

func main() {
	log.SetFlags(0)
	g := Default
	f := g.Freq()
	fmt.Printf("US 10,144,532 B2 annular cavity — Kron (Happ two-level) analysis\n")
	fmt.Printf("geometry (declared): R %.2f m, gap %.0f mm, height %.2f m, %d sections, %d Kron sectors; copper walls\n",
		g.R, g.D*1e3, g.H, g.N, g.S)
	fmt.Printf("mode n = %d at %.3f MHz (λ = %.3f m); Z0 = %.1f Ω\n\n", g.Mode, f/1e6, c0/f, 376.730313*g.D/g.H)

	const P = 1e6 // 1 MW total drive, declared
	q := g.N / 4
	four := []Emitter{{0, 1, 0}, {q, 1, 0}, {2 * q, 1, 0}, {3 * q, 1, 0}}
	show := func(name string, r *Result, n int) {
		fmt.Printf("%s\n  f %.4f MHz  Q %.4g  in %.4g W  stored %.4g J  Vmax %.4g kV  Emax %.4g MV/m  RF σ %.3g µC/m²\n",
			name, r.F/1e6, r.Q, r.Pin, r.W, r.Vmax/1e3, r.Emax/1e6, r.SigmaRF*1e6)
		fmt.Printf("  energy by quadrant %.3f %.3f %.3f %.3f   ledger residual (in − wall loss)/in %.1e\n",
			r.Quadrant[0], r.Quadrant[1], r.Quadrant[2], r.Quadrant[3], (r.Pin-r.Ploss)/r.Pin)
		fmt.Printf("  Kron interconnection order %d (MNA unknowns %d)\n\n", r.Interconnect, 2*n)
	}
	base, err := Solve(g, f, four, P)
	if err != nil {
		log.Fatal(err)
	}
	show("A. 4 emitters at 0/90/180/270°, in phase (the ring-tower synchronization)", base, g.N)
	// Steering: the patent drives shroud sections separately. Same emitter
	// current as one of A's four, top emitter only.
	amp := base.Scale
	top, err := Solve(g, f, []Emitter{{0, amp, 0}}, 0)
	if err != nil {
		log.Fatal(err)
	}
	show("B. top emitter only (0°), same current, on resonance", top, g.N)
	fOff := (g.Freq() + Geometry{R: g.R, D: g.D, H: g.H, N: g.N, S: g.S, Mode: g.Mode + 1}.Freq()) / 2
	off, err := Solve(g, fOff, []Emitter{{0, amp, 0}}, 0)
	if err != nil {
		log.Fatal(err)
	}
	show("C. top emitter only, same current, between modes n and n+1 (steering attempt)", off, g.N)
	// Lossier walls localize more: decay length 1/α = 2·Z0/R′.
	lossy := g
	lossy.Rs = 2000 * g.SurfaceR(f)
	lz, err := Solve(lossy, f, []Emitter{{0, amp, 0}}, 0)
	if err != nil {
		log.Fatal(err)
	}
	show("D. top emitter only, walls 2000× copper's surface resistance (a lossy, localized ring)", lz, g.N)
	// Superconducting walls: all loss is in R′, so Q ∝ 1/Rs exactly. The
	// direct solve becomes ill-conditioned (condition number grows with Q),
	// so the niobium case is scaled from A.
	rsNb := 20e-9
	qNb := base.Q * g.SurfaceR(f) / rsNb
	fmt.Printf("E. niobium walls, Rs = %.0f nΩ (declared; SRF cavities at 2 K): Q %.3g; the same %.3g J stored needs %.3g W\n\n",
		rsNb*1e9, qNb, base.W, P*rsNb/g.SurfaceR(f))
	// Twin: MNA AC on the same netlist. Off resonance the solutions agree
	// to round-off; on the Q ≈ 1.5e4 resonance both are checked by their
	// backward error (the system is ill-conditioned by about Q·n).
	for _, fx := range []float64{fOff, f} {
		r, err := Solve(g, fx, four, 0)
		if err != nil {
			log.Fatal(err)
		}
		vm, err := TwinMNA(g, fx, four)
		if err != nil {
			log.Fatal(err)
		}
		var d, vmax float64
		for i := range vm {
			vmax = math.Max(vmax, cmplx.Abs(vm[i]))
		}
		for i := range vm {
			d = math.Max(d, cmplx.Abs(vm[i]-r.V[i])/vmax)
		}
		fmt.Printf("twin at %.4f MHz: Kron vs MNA %.2g relative; backward error Kron %.2g, MNA %.2g\n",
			fx/1e6, d, residual(g, fx, four, r.V, 1), residual(g, fx, four, vm, 1))
	}
	fmt.Println()
	fmt.Println("the patent's figures against physical ceilings (copper case, 1 MW scaling):")
	for _, s := range []float64{1, 50000} {
		l := LimitsFor(s, base)
		fmt.Printf("  σ = %g C/m²: E = %.3g V/m (Kilpatrick %.3g, Schwinger %.3g); wall pressure %.3g Pa (copper UTS %.3g);\n",
			s, l.E, l.Kilpatrick, l.Schwinger, l.Pressure, l.CopperUTS)
		fmt.Printf("      Fowler–Nordheim emission %.3g A/m² drains it in %.2g s; drive power to reach it in this cavity %.3g W\n",
			l.FNCurrent, l.DrainTime, l.PowerNeeded)
	}
	omega := 30000 * 2 * math.Pi / 60
	fmt.Printf("  rotation example: S = σ²/ε0·Rω = %.3g W/m² (patent ≈ 10²⁴) — magnitude of E×H at a moving charged surface, divergence-free: no energy leaves\n",
		PatentFlux(50000, 2, omega))
	fmt.Printf("  vibration example: 10³³ W/m² at σ = 1 C/m² needs R_v·ν = %.3g m/s = %.2g c\n", 1e33*eps0, 1e33*eps0/c0)
	rot, err := spacetime.AnalyzeRotor(omega, 2, spacetime.T1000G)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  30,000 rpm at 2 m: rim %.0f m/s, (v/c)² %.2g, burst margin %.2f against a T1000G fibre hoop (>1: flies apart)\n",
		rot.RimSpeed, rot.Beta2, rot.BurstMargin)
	fmt.Printf("  wall vibration 10⁹–10¹⁸ Hz: lattice vibrations end at the Debye frequency, Cu %.2g Hz, Al %.2g Hz\n",
		Debye(343), Debye(428))
	fmt.Printf("  thrust: a closed cavity exerts zero net time-averaged force on itself; emitting all %.0f W as a beam gives at most P/c = %.3g N\n",
		P, P/c0)
}
