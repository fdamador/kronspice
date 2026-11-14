// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command ether_horizon tests hypothesis H49 (beside examples/ether_cavitation,
// the ether that tears): a nonlinear ether can open an
// empty void and a horizon (a "black" region), and such a horizon at the
// centre of the photon electron (H23) acts on its light ring as a resistor
// — the 377 Ω of the black-hole membrane.
//
// Part 1, the void. The nonlinear ether is taken as a superfluid whose
// order parameter can vanish (Gross–Pitaevskii). A straight vortex of
// winding n is solved on the radial profile equation
//
//	f'' + f'/ρ − n²f/ρ² + f − f³ = 0,   f(0) = 0, f(∞) = 1,   ρ = r/ξ,
//
// by Newton's method on a finite-difference grid; the near-axis slope, the
// core size and the line energy are compared with the literature.
//
// Part 2, the horizon. Where the ether flows faster than its own wave
// speed, waves cannot get back out (Unruh 1981; Visser 1998). The vortex of
// Part 1 has an ergoregion (supersonic swirl) but no horizon; a horizon
// needs a radial inflow into a sink — the draining bathtub, whose horizon,
// ergosurface and Hawking temperature are closed-form. If the ether's
// waves are light (c_ether = c), the horizon radiates; the lifetime this
// gives an electron-sized horizon is compared with the age of the universe
// and the Borexino electron-stability bound.
//
// Part 3, the coupling. A horizon inside the photon ring absorbs the
// ring's light. Two estimates: the membrane paradigm (the horizon a sheet
// of resistance Z0 = 376.73 Ω, Thorne, Price & Macdonald 1986) and the
// exact low-frequency photon absorption of a Schwarzschild hole (Page
// 1976). For general relativity the horizon is fixed by the mass; the
// electron's own numbers (Schwarzschild radius, naked Kerr–Newman
// parameters, Hawking lifetime, Carter's g = 2) are computed. For the ether
// horizon the radius is free (the "fat horizon" knob), and the window
// between too absorbing (small Q) and too hot (Hawking) is found.
//
// Part 4, the circuit. The crossed two-wire ring line of
// examples/photon_electron (its 720° double loop) with a horizon resistor
// across every section: kronspice AC gives the loaded Q and the mode
// frequency against the ladder's complex dispersion, the power in the
// horizon resistors is booked as its own column and the ledger closes, and
// the half-integer (spin-½) modes are checked to survive.
//
// Part 5 (sink.go) gives the inflow somewhere to go: a localized loss −iγψ
// in the Gross–Pitaevskii ether, solved for its steady state; the loss at
// which a horizon first forms, the flux and the mass ledger. Part 6 books
// the sink's columns and asks who pays the Hawking radiation. Part 7
// (far.go) is the sink electron from afar: the Lagally force between sinks
// (gravity's sign, not Coulomb's), the flux fitted to gravity against the
// ether density, and a kronspice resistor-lattice twin of the potential
// flow (KCL = continuity). Part 8 (fountain.go, spin.go, part8.go) lets the
// ether circulate instead — in round the equator, out at the poles — as an
// exact steady Euler flow with a line vortex on the axis, finds the
// streamline torus that winds 2:1 like the photon ring, and builds the
// spinning core as a battery in kronspice (membrane circuit, extended Smith
// chart, superradiance band, spin reservoir).
//
//	go run ./examples/ether_horizon [-viz DIR]
//
// With -viz, the steady sink profiles (n = 0, r_s = 1 ξ, at γ_c and at
// γ = 10) are also written to DIR/sink_profile.csv for viz/plot_sink3d.py.
package main

import (
	"encoding/csv"
	"errors"
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Constants, SI (CODATA 2018).
const (
	Hbar   = 1.054571817e-34  // J s
	C0     = 299792458.0      // m/s
	G      = 6.67430e-11      // m³ kg⁻¹ s⁻²
	Me     = 9.1093837015e-31 // kg
	E      = 1.602176634e-19  // C
	Eps0   = 8.8541878128e-12 // F/m
	Mu0    = 1.25663706212e-6 // H/m
	KB     = 1.380649e-23     // J/K
	LPl    = 1.616255e-35     // m, Planck length
	Year   = 3.15576e7        // s, Julian year
	AgeU   = 13.787e9 * Year  // s, age of the universe (Planck 2018)
	TauBX  = 6.6e28 * Year    // s, Borexino e → νγ lifetime bound (PRL 115, 231802, 2015)
	GDirac = 2.0
	GE     = 2.00231930436256 // free-electron |g|
)

// Z0 is the impedance of free space, √(μ0/ε0) [Ω].
func Z0() float64 { return math.Sqrt(Mu0 / Eps0) }

// Compton is h/(m_e c) [m]; RingRadius is the photon_electron loop radius λ/4π.
func Compton() float64    { return 2 * math.Pi * Hbar / (Me * C0) }
func RingRadius() float64 { return Compton() / (4 * math.Pi) }

// OmegaC is the electron's Compton angular frequency m c²/ħ [rad/s].
func OmegaC() float64 { return Me * C0 * C0 / Hbar }

// ---------------------------------------------------------------- Part 1

// ErrNoConvergence is returned when the vortex profile's Newton iteration fails.
var ErrNoConvergence = errors.New("ether_horizon: vortex profile did not converge")

// Vortex is the radial profile f(ρ) of a straight vortex of winding N in a
// Gross–Pitaevskii ether, ρ in healing lengths ξ = ħ/√(2mμ).
type Vortex struct {
	N   int
	Rho []float64
	F   []float64
}

// SolveVortex solves the profile on [0, R] with M intervals.
func SolveVortex(n int, R float64, M int) (*Vortex, error) {
	h := R / float64(M)
	rho := make([]float64, M+1)
	f := make([]float64, M+1)
	nn := float64(n * n)
	for i := range rho {
		rho[i] = float64(i) * h
		x := rho[i] * rho[i]
		f[i] = math.Pow(x/(x+2), float64(n)/2)
	}
	f[M] = math.Sqrt(math.Max(0, 1-nn/(R*R)))
	a, b, c, d := make([]float64, M+1), make([]float64, M+1), make([]float64, M+1), make([]float64, M+1)
	for it := 0; it < 200; it++ {
		maxd := 0.0
		for i := 1; i < M; i++ {
			r := rho[i]
			res := (f[i+1]-2*f[i]+f[i-1])/(h*h) + (f[i+1]-f[i-1])/(2*h*r) - nn*f[i]/(r*r) + f[i] - f[i]*f[i]*f[i]
			a[i] = 1/(h*h) - 1/(2*h*r)
			c[i] = 1/(h*h) + 1/(2*h*r)
			b[i] = -2/(h*h) - nn/(r*r) + 1 - 3*f[i]*f[i]
			d[i] = -res
		}
		// Thomas algorithm on interior points 1..M-1 (δf = 0 at both ends).
		cp, dp := make([]float64, M+1), make([]float64, M+1)
		for i := 1; i < M; i++ {
			den := b[i]
			if i > 1 {
				den -= a[i] * cp[i-1]
			}
			cp[i] = c[i] / den
			num := d[i]
			if i > 1 {
				num -= a[i] * dp[i-1]
			}
			dp[i] = num / den
		}
		dx := make([]float64, M+1)
		for i := M - 1; i >= 1; i-- {
			dx[i] = dp[i] - cp[i]*dx[i+1]
			f[i] += dx[i]
			maxd = math.Max(maxd, math.Abs(dx[i]))
		}
		if maxd < 1e-13 {
			return &Vortex{N: n, Rho: rho, F: f}, nil
		}
	}
	return nil, ErrNoConvergence
}

// AxisCoefficient is c in f ≈ c·ρⁿ near the axis (for n = 1 the series
// f = p1ρ − (p1/8)ρ³ + … is used on the first grid point).
func (v *Vortex) AxisCoefficient() float64 {
	r := v.Rho[1]
	if v.N == 1 {
		return v.F[1] / (r - r*r*r/8)
	}
	return v.F[1] / math.Pow(r, float64(v.N))
}

// Crossing returns the first ρ where g(ρ, f) changes sign (linear interpolation).
func (v *Vortex) Crossing(g func(rho, f float64) float64) float64 {
	for i := 1; i+1 < len(v.Rho); i++ {
		g0, g1 := g(v.Rho[i], v.F[i]), g(v.Rho[i+1], v.F[i+1])
		if g0 == 0 {
			return v.Rho[i]
		}
		if g0*g1 < 0 {
			return v.Rho[i] + (v.Rho[i+1]-v.Rho[i])*g0/(g0-g1)
		}
	}
	return math.NaN()
}

// CoreRadius is where the density f² reaches half its bulk value.
func (v *Vortex) CoreRadius() float64 {
	return v.Crossing(func(_, f float64) float64 { return f*f - 0.5 })
}

// ErgoRadius is where the swirl v = nħ/(m r) equals the local wave speed
// c∞·f: in healing lengths, √2·n/ρ = f.
func (v *Vortex) ErgoRadius() float64 {
	n := float64(v.N)
	return v.Crossing(func(r, f float64) float64 { return f - math.Sqrt2*n/r })
}

// LineEnergy is the energy per length out to ρ = R in units of πħ²n0/m:
// ∫[f'² + n²f²/ρ² + ½(1 − f²)²] ρ dρ. Far out it grows as n²·ln R.
func (v *Vortex) LineEnergy() float64 {
	nn := float64(v.N * v.N)
	s := 0.0
	for i := 0; i+1 < len(v.Rho); i++ {
		h := v.Rho[i+1] - v.Rho[i]
		rm := 0.5 * (v.Rho[i] + v.Rho[i+1])
		fm := 0.5 * (v.F[i] + v.F[i+1])
		fp := (v.F[i+1] - v.F[i]) / h
		s += (fp*fp + nn*fm*fm/(rm*rm) + 0.5*(1-fm*fm)*(1-fm*fm)) * rm * h
	}
	return s
}

// ---------------------------------------------------------------- Part 2

// Bathtub is the draining-bathtub flow v = (−A r̂ + B φ̂)/r in an ether of
// wave speed C (Visser, Class. Quantum Grav. 15, 1767, 1998).
type Bathtub struct{ A, B, C float64 }

// Horizon is the radius where the inflow equals the wave speed, A/c.
func (b Bathtub) Horizon() float64 { return b.A / b.C }

// Ergosurface is where the whole flow speed equals c, √(A² + B²)/c.
func (b Bathtub) Ergosurface() float64 { return math.Hypot(b.A, b.B) / b.C }

// HawkingT is the horizon temperature ħκ/(2π k_B) with κ = d(c − |v_r|)/dr = c/r_h [K].
func HawkingT(rh, c float64) float64 { return Hbar * c / (2 * math.Pi * KB * rh) }

// OneChannelPower is the thermal power of one 1-D channel, π(k_B T)²/(6ħ)
// [W] — the least a horizon of temperature T can emit (a 2-D or 3-D
// horizon has more channels and emits more).
func OneChannelPower(T float64) float64 { return math.Pi * (KB * T) * (KB * T) / (6 * Hbar) }

// HawkingLifetime is m_e c² over the one-channel power: an upper estimate
// of how long an ether horizon of radius rh (c_ether = c) takes to radiate
// away the electron's rest energy [s]. It grows as rh².
func HawkingLifetime(rh float64) float64 { return Me * C0 * C0 / OneChannelPower(HawkingT(rh, C0)) }

// HawkingRadiusFor is the smallest horizon radius whose HawkingLifetime is tau.
func HawkingRadiusFor(tau float64) float64 { return math.Sqrt(tau / HawkingLifetime(1)) }

// ---------------------------------------------------------------- Part 3

// GR numbers of an electron-mass hole.
type GRElectron struct {
	Rs, RG, SpinA, ChargeR float64 // Schwarzschild radius, GM/c², a = J/(Mc), r_Q = √(G e²/(4πε0 c⁴)) [m]
	HawkingLifetime        float64 // 5120πG²M³/(ħc⁴) [s]
	HawkingT               float64 // ħc³/(8πGMk_B) [K]
}

// ElectronHole computes them.
func ElectronHole() GRElectron {
	rg := G * Me / (C0 * C0)
	return GRElectron{
		Rs: 2 * rg, RG: rg,
		SpinA:           Hbar / (2 * Me * C0),
		ChargeR:         math.Sqrt(G * E * E / (4 * math.Pi * Eps0 * math.Pow(C0, 4))),
		HawkingLifetime: 5120 * math.Pi * G * G * Me * Me * Me / (Hbar * math.Pow(C0, 4)),
		HawkingT:        Hbar * math.Pow(C0, 3) / (8 * math.Pi * G * Me * KB),
	}
}

// FieldVolume is the volume the photon_electron model puts the ring's
// field in: a sphere of diameter λ/2 (Williamson & van der Mark, eq. 5).
func FieldVolume() float64 { r := Compton() / 4; return 4 * math.Pi / 3 * r * r * r }

// MembraneRate is the fraction of the ring's energy a horizon of radius rh
// absorbs per second in the membrane paradigm: a sphere of sheet
// resistance Z0 in the ring's central field takes ∫E_t²/Z0 dA =
// (8π/3)rh²·c·u, with u = m c²/V the ring's mean energy density [1/s].
func MembraneRate(rh float64) float64 { return 8 * math.Pi / 3 * rh * rh * C0 / FieldVolume() }

// PageRate is the same with the exact low-frequency photon absorption
// cross-section of a Schwarzschild hole of horizon rh,
// σ = (64π/3)(r_h/2)⁴(ω/c)² (Page 1976), at the ring's wave number 2π/λ.
func PageRate(rh float64) float64 {
	k := 2 * math.Pi / Compton()
	sigma := 64 * math.Pi / 3 * math.Pow(rh/2, 4) * k * k
	return sigma * C0 / FieldVolume()
}

// HorizonQ is the ring's Q against horizon absorption alone, ω_C/rate.
func HorizonQ(rate float64) float64 { return OmegaC() / rate }

// MembraneRadiusFor is the largest horizon (membrane estimate) that keeps
// the ring's absorption time above tau.
func MembraneRadiusFor(tau float64) float64 { return math.Sqrt(1 / (tau * MembraneRate(1))) }

// PageRadiusFor is the same with Page's cross-section.
func PageRadiusFor(tau float64) float64 { return math.Pow(1/(tau*PageRate(1)), 0.25) }

// ---------------------------------------------------------------- Part 4

// Ring is the two-wire ring line of examples/photon_electron (copied, so
// the two examples stay independent): N sections, series L and R per
// wire, Cd between the wires and Cg from each wire to ground, crossed at
// one point for the 720° double loop. RH is the horizon resistor across
// each section (0: none).
type Ring struct {
	N         int
	L, Cd, Cg float64
	R, RH     float64
	Crossed   bool
}

// DefaultRing is photon_electron's: 200 sections, 1 µH, 1 nF, 0.1 nF, 0.01 Ω, crossed.
func DefaultRing() Ring { return Ring{N: 200, L: 1e-6, Cd: 1e-9, Cg: 0.1e-9, R: 0.01, Crossed: true} }

// Cdiff is the differential capacitance per section, Cd + Cg/2.
func (r Ring) Cdiff() float64 { return r.Cd + r.Cg/2 }

// Mode returns the complex mode m of the lossy ladder from its dispersion
// Z_s·Y = −4 sin²(θ/2), Z_s = 2(R + sL), Y = sC + 1/RH, θ = 2π(m + ½)/N
// (crossed) or 2πm/N: the damped frequency [Hz] and Q = ω_d/(2α).
func (r Ring) Mode(m int) (f, q float64) {
	k := float64(m)
	if r.Crossed {
		k += 0.5
	}
	K := 4 * math.Pow(math.Sin(math.Pi*k/float64(r.N)), 2)
	g := 0.0
	if r.RH > 0 {
		g = 1 / r.RH
	}
	C := r.Cdiff()
	a2, a1, a0 := 2*r.L*C, 2*r.R*C+2*r.L*g, 2*r.R*g+K
	s := (-complex(a1, 0) + cmplx.Sqrt(complex(a1*a1-4*a2*a0, 0))) / complex(2*a2, 0)
	if imag(s) < 0 {
		s = cmplx.Conj(s)
	}
	return imag(s) / (2 * math.Pi), imag(s) / (2 * -real(s))
}

// HorizonRFor is the horizon resistor that alone gives the mode a Q of qh:
// Q_h = ω·Cdiff·RH.
func (r Ring) HorizonRFor(m int, qh float64) float64 {
	f, _ := Ring{N: r.N, L: r.L, Cd: r.Cd, Cg: r.Cg, Crossed: r.Crossed}.Mode(m)
	return qh / (2 * math.Pi * f * r.Cdiff())
}

// Build makes the ring with a 1 A AC current source across (a0, b0).
func (r Ring) Build() (*circuit.Circuit, error) {
	b := circuit.NewBuilder("crossed ring with horizon resistors")
	a := func(k int) string { return "a" + strconv.Itoa(k%r.N) }
	bb := func(k int) string { return "b" + strconv.Itoa(k%r.N) }
	for k := 0; k < r.N; k++ {
		ks := strconv.Itoa(k)
		na, nb := a(k+1), bb(k+1)
		if k == r.N-1 && r.Crossed {
			na, nb = bb(0), a(0)
		}
		b.Inductor("LA"+ks, a(k), "ma"+ks, r.L)
		b.Resistor("RA"+ks, "ma"+ks, na, r.R)
		b.Inductor("LB"+ks, bb(k), "mb"+ks, r.L)
		b.Resistor("RB"+ks, "mb"+ks, nb, r.R)
		b.Capacitor("CD"+ks, a(k), bb(k), r.Cd)
		b.Capacitor("CGA"+ks, a(k), "0", r.Cg)
		b.Capacitor("CGB"+ks, bb(k), "0", r.Cg)
		if r.RH > 0 {
			b.Resistor("RH"+ks, a(k), bb(k), r.RH) // the horizon, one share per section
		}
	}
	b.Resistor("RBA", a(0), "0", 1e6) // DC bleeders, as in photon_electron
	b.Resistor("RBB", bb(0), "0", 1e6)
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "b0", "a0")
	return b.Build()
}

// Measured is one ring's mode seen by kronspice AC.
type Measured struct {
	F, Q              float64 // peak frequency [Hz], loaded Q from the half-power width
	PIn, PWire, PHor  float64 // at the peak: power in, in the wire resistors, in the horizon resistors [W]
	PBleed            float64 // in the two 1 MΩ bleeders [W]
	MaxLedgerResidual float64 // kronspice's own AC ledger
}

// Measure sweeps ±span around the analytic mode m and reads the peak of
// the differential impedance, its half-power width and the power split.
func (r Ring) Measure(m int, span float64, points int) (Measured, error) {
	c, err := r.Build()
	if err != nil {
		return Measured{}, err
	}
	f0, _ := r.Mode(m)
	// The grid is set off-centre so that no point falls on f0 by construction.
	res, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: points, FStart: f0 * (1 - span), FStop: f0 * (1 + 1.07*span)}, analysis.Options{})
	if err != nil {
		return Measured{}, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	z := make([]float64, len(res.Freq))
	k := 0
	for i := range res.Freq {
		z[i] = cmplx.Abs(va[i] - vb[i])
		if z[i] > z[k] {
			k = i
		}
	}
	half := z[k] / math.Sqrt2
	edge := func(dir int) float64 {
		for i := k; i+dir >= 0 && i+dir < len(z); i += dir {
			if z[i+dir] < half {
				return res.Freq[i] + (res.Freq[i+dir]-res.Freq[i])*(z[i]-half)/(z[i]-z[i+dir])
			}
		}
		return math.NaN()
	}
	// Parabolic refinement of the peak on its three points.
	fpk := res.Freq[k]
	if k > 0 && k+1 < len(z) {
		y0, y1, y2 := z[k-1], z[k], z[k+1]
		fpk += 0.5 * (y0 - y2) / (y0 - 2*y1 + y2) * (res.Freq[1] - res.Freq[0])
	}
	out := Measured{F: fpk, Q: res.Freq[k] / (edge(1) - edge(-1)), MaxLedgerResidual: res.Ledger.MaxRelativeResidual}
	out.PIn = 0.5 * real(va[k]-vb[k]) // 1 A into a0, out of b0
	for s := 0; s < r.N; s++ {
		ss := strconv.Itoa(s)
		for _, w := range []string{"LA", "LB"} {
			i, _ := res.ElementCurrent(w + ss)
			out.PWire += 0.5 * r.R * math.Pow(cmplx.Abs(i[k]), 2)
		}
		if r.RH > 0 {
			x, _ := res.NodeVoltage("a" + ss)
			y, _ := res.NodeVoltage("b" + ss)
			out.PHor += 0.5 * math.Pow(cmplx.Abs(x[k]-y[k]), 2) / r.RH
		}
	}
	out.PBleed = 0.5 * (math.Pow(cmplx.Abs(va[k]), 2) + math.Pow(cmplx.Abs(vb[k]), 2)) / 1e6
	return out, nil
}

// SignAfterTrip reports, from the ladder's wave number, whether a wave
// comes back inverted after one trip: the crossed ring's half-integer
// modes (θ·N = π(2m + 1)) do, with or without the horizon resistors (they
// only add loss, they do not change the boundary condition).
func (r Ring) SignAfterTrip(m int) float64 {
	k := float64(m)
	if r.Crossed {
		k += 0.5
	}
	return math.Cos(2 * math.Pi * k)
}

// ---------------------------------------------------------------- main

// WriteSinkProfiles writes the n = 0, r_s = 1 ξ sink at γ_c and at γ = 10
// as columns r_xi, rho_c, mach_c, rho_10, mach_10 (r in ξ, ρ in ρ∞).
func WriteSinkProfiles(path string) error {
	rs := 1 / math.Sqrt2
	_, sc, err := CriticalGamma(0, rs, 60, 0.002)
	if err != nil {
		return err
	}
	s10, err := SolveSink(0, 10, rs, 60, 0.002)
	if err != nil {
		return err
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	w := csv.NewWriter(f)
	w.Write([]string{"r_xi", "rho_c", "mach_c", "rho_10", "mach_10"})
	n := len(sc.R)
	if len(s10.R) < n {
		n = len(s10.R)
	}
	for i := 0; i < n; i += 10 {
		row := []float64{sc.R[i] * math.Sqrt2, sc.A[i] * sc.A[i], math.Abs(sc.V[i]) / sc.A[i], s10.A[i] * s10.A[i], math.Abs(s10.V[i]) / s10.A[i]}
		rec := make([]string, len(row))
		for j, v := range row {
			rec[j] = strconv.FormatFloat(v, 'g', 8, 64)
		}
		w.Write(rec)
	}
	w.Flush()
	if err := w.Error(); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

func main() {
	log.SetFlags(0)
	viz := flag.String("viz", "", "also write the sink profiles to DIR/sink_profile.csv")
	flag.Parse()
	if *viz != "" {
		if err := WriteSinkProfiles(filepath.Join(*viz, "sink_profile.csv")); err != nil {
			log.Fatal(err)
		}
	}
	fmt.Println("H49: void and horizon in a nonlinear ether, tested on the photon electron")
	fmt.Println()

	fmt.Println("Part 1. The void: a Gross–Pitaevskii vortex (ρ in healing lengths ξ)")
	fmt.Println("  n   axis: f ≈ c·ρⁿ   f(0)   core (f² = ½)   ergoregion (swirl = local wave speed)   line energy − n²ln R")
	for n := 1; n <= 3; n++ {
		v, err := SolveVortex(n, 40, 8000)
		if err != nil {
			log.Fatal(err)
		}
		e := v.LineEnergy() - float64(n*n)*math.Log(40)
		fmt.Printf("  %d   c = %.6f       %.0f      %.3f ξ          %.3f ξ                                %.4f\n",
			n, v.AxisCoefficient(), v.F[0], v.CoreRadius(), v.ErgoRadius(), e)
	}
	fmt.Println("  n = 1: c against 0.5827811878 (Rorai, Sreenivasan & Fisher, arXiv:1212.6389);")
	fmt.Printf("  line energy π(ħ²n0/m)·ln(%.3f R/ξ) (Fetter, RMP 81, 647, 2009: ln(1.46 R/ξ)).\n", math.Exp(func() float64 {
		v, _ := SolveVortex(1, 40, 8000)
		return v.LineEnergy() - math.Log(40)
	}()))
	fmt.Println("  The axis is empty for every winding n ≥ 1, and only there: the void is a line (a thin")
	fmt.Println("  torus for a vortex ring) one healing length thick, not a sphere. Below one quantum of")
	fmt.Println("  circulation (n = 0) there is no vortex and no void.")
	fmt.Println()

	fmt.Println("Part 2. The horizon")
	fmt.Println("  The vortex alone: supersonic inside the ergoregion, but the flow is all swirl,")
	fmt.Println("  so waves can still leave radially — an ergoregion, no horizon.")
	bt := Bathtub{A: RingRadius() * C0, B: 2 * RingRadius() * C0, C: C0}
	fmt.Printf("  Draining bathtub with c_ether = c, inflow set for r_h = ring radius %.4g m:\n", RingRadius())
	fmt.Printf("    horizon %.4g m, ergosurface %.4g m (B = 2A), T_H = %.4g K (k_B T = %.4g keV)\n",
		bt.Horizon(), bt.Ergosurface(), HawkingT(bt.Horizon(), C0), KB*HawkingT(bt.Horizon(), C0)/E/1e3)
	fmt.Printf("    one-channel power %.4g W → electron's rest energy gone in %.3g s\n",
		OneChannelPower(HawkingT(bt.Horizon(), C0)), HawkingLifetime(bt.Horizon()))
	fmt.Printf("    stable for the age of the universe needs r_h ≥ %.3g m; for Borexino's bound, ≥ %.3g m\n",
		HawkingRadiusFor(AgeU), HawkingRadiusFor(TauBX))
	fmt.Println("    The inflow needs a sink: 2πA·ρ0 of ether per unit length per second leaves the")
	fmt.Println("    picture there — a ledger column the model must name.")
	fmt.Println()

	fmt.Println("Part 3. The coupling: a horizon inside the photon ring")
	h := ElectronHole()
	fmt.Printf("  GR, electron mass: r_s = %.4g m (%.3g Planck lengths), ring radius %.4g m\n", h.Rs, h.Rs/LPl, RingRadius())
	fmt.Printf("    spin a/M = %.3g, charge r_Q/M = %.3g: both ≫ 1, Kerr–Newman has no horizon (naked)\n", h.SpinA/h.RG, h.ChargeR/h.RG)
	fmt.Printf("    if it had one: Hawking T = %.3g K, evaporates in %.3g s\n", h.HawkingT, h.HawkingLifetime)
	fmt.Printf("    Carter (1968): Kerr–Newman g = %.0f exactly; Dirac 2; measured %.14g (anomaly is QED)\n", GDirac, GE)
	fmt.Printf("    absorption Q of the ring at r_s: membrane %.3g, Page %.3g\n", HorizonQ(MembraneRate(h.Rs)), HorizonQ(PageRate(h.Rs)))
	fmt.Println("  Ether horizon, radius free (the fat-horizon knob), c_ether = c:")
	fmt.Println("  (Page's column: if it absorbed like a GR hole of that radius; low-frequency formula, k·r_h ≪ 1)")
	fmt.Println("    r_h [m]        Q_h membrane   Q_h Page      absorbs ring in [s]   Hawking lifetime [s]")
	for _, rh := range []float64{1e-35, 1e-32, 1e-28, 1e-24, 1e-20, 1e-16, RingRadius()} {
		fmt.Printf("    %-13.3g  %-13.3g  %-13.3g %-21.3g %.3g\n", rh, HorizonQ(MembraneRate(rh)), HorizonQ(PageRate(rh)), 1/MembraneRate(rh), HawkingLifetime(rh))
	}
	fmt.Printf("  Ring survives the age of the universe: r_h ≤ %.3g m (membrane), ≤ %.3g m (Page)\n", MembraneRadiusFor(AgeU), PageRadiusFor(AgeU))
	fmt.Printf("  Horizon survives its own Hawking emission as long: r_h ≥ %.3g m\n", HawkingRadiusFor(AgeU))
	fmt.Printf("  Window: none — the two bounds miss by %.1f (membrane) and %.1f (Page) orders of magnitude.\n",
		math.Log10(HawkingRadiusFor(AgeU)/MembraneRadiusFor(AgeU)), math.Log10(HawkingRadiusFor(AgeU)/PageRadiusFor(AgeU)))
	fmt.Println()

	fmt.Println("Part 4. The circuit: photon_electron's crossed ring with a horizon resistor per section")
	fmt.Println("  Q_h set     RH [Ω]       f peak [Hz]   f ladder [Hz]   Q AC     Q ladder   P_in [W]     wires        horizon      residual")
	for _, qh := range []float64{0, 1000, 100, 30} {
		r := DefaultRing()
		if qh > 0 {
			r.RH = r.HorizonRFor(0, qh)
		}
		fl, ql := r.Mode(0)
		m, err := r.Measure(0, 0.25, 6001)
		if err != nil {
			log.Fatal(err)
		}
		qs := "∞"
		if qh > 0 {
			qs = fmt.Sprintf("%.0f", qh)
		}
		fmt.Printf("  %-10s  %-11.4g  %-12.2f  %-14.2f  %-7.2f  %-9.2f  %-11.5g  %-11.5g  %-11.5g  %.2g\n",
			qs, r.RH, m.F, fl, m.Q, ql, m.PIn, m.PWire, m.PHor, (m.PIn-m.PWire-m.PHor-m.PBleed)/m.PIn)
	}
	r := DefaultRing()
	fmt.Printf("  sign after one trip, lowest mode: %+.0f with and without the horizon (spin-½ boundary kept)\n", r.SignAfterTrip(0))
	f0, _ := r.Mode(0)
	r.RH = r.HorizonRFor(0, 30)
	f30, _ := r.Mode(0)
	fmt.Printf("  ladder frequency shift at Q_h = 30: %.3g, second order in the losses (≈ 1/(Q_wire·Q_h)/8);\n", f30/f0-1)
	fmt.Println("  at the physical Q_h > 10³⁸ (Part 3) the shift, and with it any change of g, is below 10⁻³⁸")
	fmt.Println()
	part5()
	part6()
	part7()
	part8(*viz)
}

// SinkRadii are the sink sizes scanned, in healing lengths ξ.
var SinkRadii = []float64{0.5, 1, 2, 5}

func part5() {
	fmt.Println("Part 5. The sink: a localized loss −iγψ in the Gross–Pitaevskii ether, steady inflow")
	fmt.Println("  (units: radii in ξ, γ in μ/ħ, flux in ρ∞ħ/m per unit length; bulk wave speed 1)")
	fmt.Println("  n  r_s [ξ]  γ_c (Mach 1 at the edge)  flux at γ_c  loss at γ_c   γ = 10: flux   outer horizon [ξ]  Mach at edge")
	for _, n := range []int{0, 1} {
		for _, rx := range SinkRadii {
			rs := rx / math.Sqrt2 // code units: ξ = 1/√2
			gc, sc, err := CriticalGamma(n, rs, 60, 0.002)
			if err != nil {
				log.Fatal(err)
			}
			s10, err := SolveSink(n, 10, rs, 60, 0.002)
			if err != nil {
				log.Fatal(err)
			}
			fmt.Printf("  %d  %-7.1f  %-24.4f  %-11.4f  %-11.4f  %-13.4f  %-17.3f  %.2f\n",
				n, rx, gc, -sc.FluxFar, sc.LossRate, -s10.FluxFar, s10.OuterHorizon()*math.Sqrt2, s10.MachAt(rs))
		}
	}
	fmt.Println("  Above γ_c the inflow passes the local wave speed outside the sink, in loss-free ether:")
	fmt.Println("  a true acoustic horizon. Past the flux maximum a stronger sink takes less (the density")
	fmt.Println("  is pushed out of it, the quantum Zeno effect of localized loss). A vortex (n = 1) feeds")
	fmt.Println("  a small sink poorly: its axis is already empty.")
	fmt.Println()
}

func part6() {
	fmt.Println("Part 6. Where the ether goes: the sink's columns")
	s, err := SolveSink(0, 10, 1/math.Sqrt2, 60, 0.002)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("  mass (n = 0, r_s = 1 ξ, γ = 10): inflow %.6f = loss %.6f (closes to %.1e)\n",
		-s.FluxFar, s.LossRate, math.Abs(s.LossRate+s.FluxFar)/s.LossRate)
	fmt.Printf("  energy: each particle brings the Bernoulli constant μ = 1, so the sink column takes %.6f per unit time\n", s.EnergyColumn())
	fmt.Println("  The column is ether destroyed (or turned into something the model does not hold yet).")
	fmt.Println("  Who pays the Hawking radiation: in a sink flow the flow does, as the pump does in the")
	fmt.Println("  laboratory analogues — so Part 2's lifetime (paid from the electron's rest energy) is not")
	fmt.Println("  the bound. The radiation still leaves, and matter at rest does not glow:")
	fmt.Printf("    Earth's heat flow %.0f TW over its %.3g electrons: ≤ %.2g W per electron\n", EarthHeat/1e12, EarthElectrons, HeatBoundPerElectron())
	fmt.Printf("    an electron-sized horizon (%.3g m) emits %.3g W; a horizon must be ≥ %.3g m to stay under the bound\n",
		RingRadius(), HawkingPowerPaidByFlow(RingRadius()), HeatHorizonBound())
	fmt.Println("    Horizons that large round every electron would overlap through all matter: no ether")
	fmt.Println("    wave (light, if c_ether = c) could leave a planet. So a sink electron must stay subsonic:")
	fmt.Println("    its core larger than the radius where the inflow would reach c — no horizon.")
	fmt.Println()
}

func part7() {
	fmt.Println("Part 7. A sink electron seen from afar")
	fmt.Println("  Lagally: two sinks attract with ρe·Q1·Q2/(4π d²) — inverse square, like attracts like:")
	fmt.Printf("  gravity's sign, the opposite of Coulomb's (two electrons repel %.2g times harder).\n", CoulombOverGravity())
	fmt.Println("  Fitted to gravity, ρe·Q²/4π = G m_e², the flux depends only on the ether density:")
	fmt.Println("    ρe [kg/m³]            Q [m³/s]     swallowed [kg/s]  growth [1/yr]  energy column [W]  r_h if point [m]")
	for _, rho := range []float64{RhoForGrowth(GMSunRate), RhoCrit, 1e3, 1e17} {
		e := SinkElectron{RhoE: rho}
		fmt.Printf("    %-20.3g  %-11.3g  %-16.3g  %-13.3g  %-17.3g  %.3g\n", rho, e.Q(), e.MassRate(), e.GrowthRate()*Year, e.EnergyColumn(), e.Horizon())
	}
	fmt.Printf("  If the swallowed ether stayed, every gravitating body would grow at √(4πGρe), the Sun too;\n")
	fmt.Printf("  ephemerides hold d(GM_sun)/dt/GM_sun to ~1e-13 per year, so ρe ≤ %.2g kg/m³ (≈ %.1e of critical).\n",
		RhoForGrowth(GMSunRate), RhoForGrowth(GMSunRate)/RhoCrit)
	fmt.Println("  Otherwise the ether is destroyed and the energy column above leaves the books — and if it")
	fmt.Printf("  turned into heat instead, the heat bound would need ρe ≤ %.2g kg/m³.\n", RhoForHeat())
	fmt.Println("  At every density listed the point-sink horizon is below 10⁻¹⁴ m, inside the ring radius")
	fmt.Println("  (1.93×10⁻¹³ m): a sink whose core is the photon ring is subsonic, with no horizon.")
	fmt.Println("  Antimatter: a positron is a sink too (it gravitates); ALPHA-g saw antihydrogen fall at")
	fmt.Println("  0.75 ± 0.13 ± 0.16 g (Nature 621, 716, 2023) — consistent, not a test of sinks.")
	g := &SinkGrid{N: 21, Sinks: [][3]int{{10, 10, 10}}, Q: []float64{1}}
	if err := g.Solve(); err != nil {
		log.Fatal(err)
	}
	fmt.Println("  kronspice twin: 21³ resistor lattice (potential flow: φ = voltage, flow = current,")
	fmt.Println("  continuity = KCL), faces at rest, one sink drawing Q = 1:")
	fmt.Printf("    inflow through the faces %.15f (= Q), source power %.12g = resistor power %.12g, ledger %.1g\n", g.Bound, g.PIn, g.PRes, g.Ledger)
	fmt.Println("    inflow speed × 4πd² (1 in the continuum):")
	for _, d := range []int{2, 3, 4, 5, 6} {
		fmt.Printf("      d = %d cells: %.4f\n", d, -g.Speed(10+d, 10, 10)*4*math.Pi*float64(d*d))
	}
	fmt.Println("    (lattice discreteness near the sink, the faces far from it). A second sink at d takes")
	fmt.Println("    ρe·Q2·v1 towards the first: the attraction above, from the lattice's own v1.")
}
