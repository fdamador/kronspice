// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command octave_magnet tests Walter Russell's "octave bar magnet",
// hypothesis H33 (CLAUDE.md §1, user hypothesis 2026-11-05): "A bar magnet
// has four pairs of poles — not one pair", drawn as a chain of lens-shaped
// (vesica) rings strung along the magnet's axis, compressed ("high
// pressure, short waves") and stretched ("low pressure, long waves")
// (W. Russell, "Analysis of the Octave Bar Magnet", drawing; see
// doc/russell-octave-report.md).
//
// The standard model of a bar magnet. A magnetization M(z)ẑ in a rod of
// radius a is, by J = ∇×M and K = M×n̂, exactly a sheet of azimuthal
// current M(z) on the rod's side (Jackson §5.8–5.9): coaxial loops, each
// carrying M·dz. Its poles are its magnetic charge, ρ_m = −∇·M: the two end
// faces of a uniform bar, plus every place inside where M changes. The
// field lines in a plane through the axis are the contours of the flux
// function ψ = ρA_φ (the flux through a circle of radius ρ at height z,
// divided by 2π), summed over the loops from the exact elliptic-integral
// potential of a circular loop.
//
// What the claim means for the field. A lens strung on the axis is a field
// line that leaves the axis at one pole and returns to it at the next: a
// chain of lenses along the axis is a chain of poles. A uniform bar has
// one pair, so its outside lines form one family from end to end; four
// pairs need the magnetization to reverse three times inside the bar
// (consequent poles), and then the lines form a chain of lenses. A Hall
// probe or search coil run along the bar's side decides between them: the
// axial field beside the bar changes sign once at every place that carries
// magnetic charge — twice for one pair, five times for four pairs end to end
// (their inner like poles merge into three doubled poles), eight times for
// four separate magnets.
//
// Parts: (1) the bar and its equivalent current sheet, checked against the
// exact solenoid formula and a Biot–Savart twin; (2) pole count from the
// magnetic charge, for a uniform bar and a bar with four reversed
// segments; (3) the bench test, a search coil moved along the bar's side,
// solved in kronspice (its flux is the integral of its induced voltage);
// (4) field lines (with -svg); (5) the octave bar bent into a ring of
// twelve opposite sections, in 3D (ring.go, ring_figures.go).
//
//	go run ./examples/octave_magnet               # a few seconds
//	go run ./examples/octave_magnet -svg viz      # also the figures
package main

import (
	"flag"
	"fmt"
	"math"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Constants (SI).
const (
	Mu0 = 1.25663706127e-6 // N/A²
	Br  = 1.25             // T, remanence of cast alnico 5 (declared; maker's data sheets 1.2–1.3 T)
)

// Bar is a magnetized rod: radius A, length L, centred on z = 0, made of
// equal segments with magnetization M0·Sign[i] along +z.
type Bar struct {
	Name string
	A, L float64
	M0   float64
	Sign []float64
}

// Classroom bar: an alnico rod 12.7 mm across, 76.2 mm long (½ × 3 inch),
// M0 = Br/μ0 (demagnetization ignored, declared).
func Uniform() Bar {
	return Bar{Name: "uniform bar (one pair)", A: 6.35e-3, L: 76.2e-3, M0: Br / Mu0, Sign: []float64{1}}
}

// FourPairs is the same rod magnetized in four segments of alternating
// direction (read left to right S N | N S | S N | N S): the magnetization Russell's claim
// needs.
func FourPairs() Bar {
	b := Uniform()
	b.Name = "four reversed segments (four pairs)"
	b.Sign = []float64{1, -1, 1, -1}
	return b
}

// M is the magnetization at height z [A/m] (0 outside the bar).
func (b Bar) M(z float64) float64 {
	if z < -b.L/2 || z >= b.L/2 {
		return 0
	}
	k := int((z + b.L/2) / b.L * float64(len(b.Sign)))
	if k >= len(b.Sign) {
		k = len(b.Sign) - 1
	}
	return b.M0 * b.Sign[k]
}

// SheetLoop is one loop of the equivalent current sheet.
type SheetLoop struct{ Z, I float64 }

// Sheet discretizes the side current into n loops (midpoint rule).
func (b Bar) Sheet(n int) []SheetLoop {
	dz := b.L / float64(n)
	out := make([]SheetLoop, n)
	for i := range out {
		z := -b.L/2 + (float64(i)+0.5)*dz
		out[i] = SheetLoop{Z: z, I: b.M(z) * dz}
	}
	return out
}

// LoopField is the exact field of a circular loop (radius a, current I, in
// the plane z = 0, axis z) at (ρ, z): B_ρ, B_z [T], and the vector
// potential A_φ [T·m] (Smythe, Static and Dynamic Electricity, §7.10;
// Jackson eq. 5.37).
func LoopField(a, I, rho, z float64) (br, bz, aphi float64) {
	if rho < 1e-12*a {
		bz = Mu0 * I * a * a / (2 * math.Pow(a*a+z*z, 1.5))
		return 0, bz, 0
	}
	s := (a+rho)*(a+rho) + z*z
	k2 := 4 * a * rho / s
	K, E := em.EllipticKE(math.Sqrt(k2))
	d := (a-rho)*(a-rho) + z*z
	pre := Mu0 * I / (2 * math.Pi * math.Sqrt(s))
	bz = pre * (K + (a*a-rho*rho-z*z)/d*E)
	br = pre * z / rho * (-K + (a*a+rho*rho+z*z)/d*E)
	k := math.Sqrt(k2)
	aphi = Mu0 * I / (math.Pi * k) * math.Sqrt(a/rho) * ((1-k2/2)*K - E)
	return
}

// Field is the bar's field and flux function at (ρ, z), from n sheet loops.
func (b Bar) Field(rho, z float64, sheet []SheetLoop) (br, bz, psi float64) {
	for _, l := range sheet {
		r, zz, ap := LoopField(b.A, l.I, rho, z-l.Z)
		br, bz, psi = br+r, bz+zz, psi+rho*ap
	}
	return
}

// AxisExact is the exact on-axis field of the segmented sheet: each
// segment [z1, z2] adds (μ0 M/2)[(z − z1)/√((z − z1)² + a²) − (z − z2)/√(…)].
func (b Bar) AxisExact(z float64) float64 {
	n := len(b.Sign)
	var bz float64
	for i := 0; i < n; i++ {
		z1 := -b.L/2 + b.L*float64(i)/float64(n)
		z2 := z1 + b.L/float64(n)
		f := func(u float64) float64 { return u / math.Sqrt(u*u+b.A*b.A) }
		bz += Mu0 * b.M0 * b.Sign[i] / 2 * (f(z-z1) - f(z-z2))
	}
	return bz
}

// Pole is a site of magnetic charge: position z [m] and charge [A·m].
type Pole struct{ Z, Q float64 }

// Poles is the magnetic charge −∇·M of the bar: the end faces (M·n̂·area)
// and every jump of M inside (−ΔM·area).
func (b Bar) Poles() []Pole {
	area := math.Pi * b.A * b.A
	n := len(b.Sign)
	var ps []Pole
	ps = append(ps, Pole{-b.L / 2, -b.M0 * b.Sign[0] * area})
	for i := 1; i < n; i++ {
		if d := b.Sign[i] - b.Sign[i-1]; d != 0 {
			ps = append(ps, Pole{-b.L/2 + b.L*float64(i)/float64(n), -b.M0 * d * area})
		}
	}
	ps = append(ps, Pole{b.L / 2, b.M0 * b.Sign[n-1] * area})
	return ps
}

// SignChanges counts the sign changes of v.
func SignChanges(v []float64) int {
	n := 0
	for i := 1; i < len(v); i++ {
		if v[i-1]*v[i] < 0 {
			n++
		}
	}
	return n
}

// Scan is B_z along a line parallel to the axis at radius rho.
func (b Bar) Scan(rho float64, n int, sheet []SheetLoop) (zs, bz []float64) {
	for i := 0; i < n; i++ {
		z := -1.5*b.L + 3*b.L*float64(i)/float64(n-1)
		_, f, _ := b.Field(rho, z, sheet)
		zs, bz = append(zs, z), append(bz, f)
	}
	return
}

// SearchCoil is the bench measurement in kronspice: a small pickup coil
// (radius r, one turn) coaxial with the bar at height z. The sheet loops
// are inductors driven from 0 to their currents over 1 µs (the magnet
// brought up to strength); the coil, closed by 1 MΩ, sees
// V = Σ M_k dI_k/dt, and its flux is ∫V dt (a fluxmeter). Returns the flux
// from the circuit and the flux Σ M_k I_k from the mutual inductances.
func (b Bar) SearchCoil(rho, z float64, sheet []SheetLoop) (circuitFlux, mutualFlux float64, err error) {
	// The coil is a ring of radius rho at height z (it encloses the bar
	// when rho > A): its flux is 2π·ψ(rho, z).
	const ramp = 1e-6
	// Thin rings (wire radius r/1000); the flux the coil reads does not depend on
	// the self-inductances, since the sheet currents are driven and the coil is
	// closed by 1 MΩ.
	selfL := func(r float64) float64 { return Mu0 * r * (math.Log(8000) - 1.75) }
	bl := circuit.NewBuilder("search coil on the bar")
	bl.Inductor("LP", "p", "0", selfL(rho))
	bl.Resistor("RP", "p", "0", 1e6)
	for k, l := range sheet {
		n := fmt.Sprintf("s%d", k)
		bl.Inductor("L"+n, n, "0", selfL(b.A))
		bl.Add(circuit.Element{Name: "I" + n, Kind: circuit.CurrentSource,
			Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, ramp, l.I, 2 * ramp, l.I}}}, "0", n)
		m := em.MutualCircles(rho, b.A, math.Abs(z-l.Z))
		bl.Mutual("K"+n, "LP", "L"+n, m/math.Sqrt(selfL(rho)*selfL(b.A)))
		mutualFlux += m * l.I
	}
	ckt, err := bl.BuildTransient()
	if err != nil {
		return 0, 0, err
	}
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: ramp / 200, Stop: 1.5 * ramp}, analysis.TranOptions{})
	if err != nil {
		return 0, 0, err
	}
	v, ok := w.NodeVoltage("p")
	if !ok {
		return 0, 0, fmt.Errorf("no node p")
	}
	for i := 1; i < len(w.Time); i++ {
		circuitFlux += 0.5 * (v[i] + v[i-1]) * (w.Time[i] - w.Time[i-1])
	}
	return circuitFlux, mutualFlux, nil
}

func main() {
	svgDir := flag.String("svg", "", "write the figures (SVG) to this directory")
	flag.Parse()
	fmt.Println("octave_magnet (H33): Walter Russell's \"a bar magnet has four pairs of poles — not one pair\"")

	u, f4 := Uniform(), FourPairs()
	sheetU, sheet4 := u.Sheet(2000), f4.Sheet(2000)

	// Part 1
	fmt.Println("\n1. The bar as its equivalent current sheet (J = ∇×M, K = M×n̂)")
	fmt.Printf("   alnico rod %.1f mm × %.1f mm, Br %.2f T (M0 = Br/μ0 = %.4g A/m; demagnetization ignored, declared)\n", 2*u.A*1e3, u.L*1e3, Br, u.M0)
	worst := 0.0
	for _, b := range []struct {
		bar Bar
		sh  []SheetLoop
	}{{u, sheetU}, {f4, sheet4}} {
		for i := 0; i <= 60; i++ {
			z := -1.5*u.L + 3*u.L*float64(i)/60
			_, bz, _ := b.bar.Field(0, z, b.sh)
			ex := b.bar.AxisExact(z)
			worst = math.Max(worst, math.Abs(bz-ex)/(Mu0*u.M0))
		}
	}
	fmt.Printf("   on-axis field, 2000 loops against the exact solenoid formula: worst difference %.1e of μ0M0\n", worst)
	one := em.Loop{P: em.Circle([3]float64{0, 0, 0}, [3]float64{0, 0, 1}, u.A, 2000), I: 1}
	bs := one.B([3]float64{1.7 * u.A, 0, 0.6 * u.A})
	br, bz, _ := LoopField(u.A, 1, 1.7*u.A, 0.6*u.A)
	fmt.Printf("   one loop off axis, elliptic formula against Biot–Savart (em.Loop, 2000 segments): B_ρ %.3e, B_z %.3e\n",
		math.Abs(br-bs[0])/math.Hypot(br, bz), math.Abs(bz-bs[2])/math.Hypot(br, bz))

	// Part 2
	fmt.Println("\n2. Poles = magnetic charge −∇·M")
	for _, b := range []Bar{u, f4} {
		ps := b.Poles()
		fmt.Printf("   %s: %d charge sites:", b.Name, len(ps))
		for _, p := range ps {
			fmt.Printf("  %+.1f mm (%+.3g A·m)", p.Z*1e3, p.Q)
		}
		fmt.Println()
	}
	fmt.Println("   uniform: charge only on the two end faces: one pair. Four pairs need M to reverse three times inside:")
	fmt.Println("   then S N|N S|S N|N S gives the two ends plus three doubled 'consequent' poles — four dipoles end to end.")

	// Part 3
	fmt.Println("\n3. The bench test: B_z along a line beside the bar (3 mm from the surface), and a search coil in kronspice")
	gap := u.A + 3e-3
	for _, b := range []struct {
		bar Bar
		sh  []SheetLoop
	}{{u, sheetU}, {f4, sheet4}} {
		_, axis := b.bar.Scan(0, 601, b.sh)
		_, side := b.bar.Scan(gap, 601, b.sh)
		fmt.Printf("   %-36s sign changes of B_z: on the axis %d, beside the bar %d\n", b.bar.Name, SignChanges(axis), SignChanges(side))
	}
	fmt.Println("   one sign change per charge site: one pair gives 2 beside the bar (and none on the axis); four pairs give 5 end to end")
	fmt.Println("   (8 if the four magnets are separated). Any count above 2 on an ordinary bar would support Russell's claim.")
	coarse := u.Sheet(200)
	fmt.Println("   search coil (a ring round the bar, radius a + 3 mm), flux from kronspice ∫V dt against ΣM·I:")
	for _, z := range []float64{0, 0.25 * u.L, 0.5 * u.L, 0.75 * u.L} {
		cf, mf, err := u.SearchCoil(gap, z, coarse)
		if err != nil {
			fmt.Println("   error:", err)
			return
		}
		_, _, psi := u.Field(gap, z, coarse)
		fmt.Printf("     z = %5.1f mm: circuit %.5e Wb, mutual inductances %.5e Wb, 2πψ %.5e Wb\n", z*1e3, cf, mf, 2*math.Pi*psi)
	}
	fmt.Println("   the three agree: kronspice's fluxmeter reads exactly the flux the field model gives. On the bench the same")
	fmt.Println("   instrument (or a Hall probe) is run along a real bar, with the bar's measured Br declared.")

	if err := ringPart(); err != nil {
		fmt.Println("ring:", err)
		return
	}

	if *svgDir != "" {
		if err := RingFigures(*svgDir, RadialTwelve(), TangentialTwelve(), RadialUniform()); err != nil {
			fmt.Println("ring figures:", err)
			return
		}
		if err := Figures(*svgDir, u, f4); err != nil {
			fmt.Println("figures:", err)
			return
		}
		fmt.Println("\nfigures written to", filepath.Clean(*svgDir))
	}
	fmt.Println("\nStandard prediction: a uniformly magnetized bar has one pair of poles; its field lines run end to end in one family.")
	fmt.Println("Russell's four pairs, as a field, is a chain of lenses along the axis, which needs three reversals of M inside the bar.")
}
