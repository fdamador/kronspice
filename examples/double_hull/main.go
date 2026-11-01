// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command double_hull models the ring-tower experiment (examples/
// ring_towers) inside a craft's double wall: the gap between an inner and
// an outer conducting hull is a closed parallel-plate waveguide around the
// craft, driven at one of its ring resonances by radiators in the gap, with
// a DC bias between the hulls added to the RF.
//
// Circuit (telegrapher's equations for the parallel-plate TEM mode; Paul,
// Analysis of Multiconductor Transmission Lines, §4.2): per unit length
// around the hull, L′ = μ0·d/w, C′ = ε0·w/d, Z0 = η0·d/w, with d the gap and
// w the hull height; wall loss R′ = 2·Rs/w, Rs = √(π·f·μ0·ρ) the surface
// resistance of each wall (skin effect), so the attenuation is
// α = Rs/(η0·d). N sections close the ring. Steady state at the drive
// frequency comes from the AC analysis (exact for this linear circuit),
// the bias from the operating point; the gap voltage is their sum.
//
// Resonant voltage: with P_in fed in at mode n, the circulating power is
// P_in·Q/(2πn) and the traveling-wave gap voltage √(Z0·P_in·Q/(π·n)); a
// standing wave doubles it at the antinodes.
//
// DC bias: a supply charges the outer hull to V_dc through a high
// resistance (no RF loading). The gap then swings between V_dc ± V_rf, so
// the peak rises to V_dc + V_rf; the bias stores ½·C_gap·V_dc² once and
// draws no steady power. Measured use: a DC bias across a coaxial gap
// suppresses multipacting — 2.5 kV in the LEP2 352 MHz couplers ("the
// multipacting events being completely suppressed for all operating
// conditions"; Kindermann, Haebel, Stirbet, Veshcherevich, CERN, EPAC 1996)
// and 2 kV in the KEKB 509 MHz couplers (Kijima et al., EPAC 2000).
//
// Multipactor check (vacuum): two-surface resonance between parallel walls
// (Hatch & Williams 1954; Vaughan, IEEE Trans. Electron Devices 35, 1172,
// 1988), by tracking electrons launched from either wall at every RF phase
// with 2 eV emission energy through the combined RF + DC field, impact after
// impact, and multiplying the secondary yield at each impact. Yield: Furman's
// fit δ = δmax·s·x/(s − 1 + x^s), x = E/Emax, with the measured copper values
// (V. Baglin, I. Collins, B. Henrist, N. Hilleret, G. Vorlaufer, CERN LHC
// Project Report 472, 2002): as received δmax 2.03, Emax 262 eV, s 1.39;
// conditioned δmax 1.13, Emax 318 eV, s 1.35. A starting phase whose mean
// yield per impact stays above 1 is a growing discharge. Space charge and the
// spread of emission energies are left out, so the check marks where a
// discharge can start, not how large it grows.
//
//	go run ./examples/double_hull
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/material"
)

// Physical constants.
const (
	Mu0  = 4e-7 * math.Pi
	Eps0 = 8.8541878128e-12
	C0   = 299792458.0
	Eta0 = Mu0 * C0
	Qe   = 1.602176634e-19
	Me   = 9.1093837015e-31
)

// Hull is the double wall and its drive.
type Hull struct {
	R, H, Gap float64 // inner-hull radius, hull height (plate width), gap [m]
	Rho       float64 // wall resistivity [Ω·m]
	Mode      int     // wavelengths around the hull
	N         int     // ring sections
	Pin       float64 // RF power fed in [W]
	Vdc       float64 // DC bias, outer hull to inner hull [V]
	Standing  bool    // standing wave (else traveling)
}

// Example is a 5 m radius, 10 m tall cylindrical craft with a 10 cm copper
// gap, mode 4, 1 kW.
func Example() Hull {
	return Hull{R: 5, H: 10, Gap: 0.1, Rho: material.Copper.Resistivity, Mode: 4, N: 96, Pin: 1000}
}

// Circumference along the middle of the gap [m].
func (h Hull) Circumference() float64 { return 2 * math.Pi * (h.R + h.Gap/2) }

// Z0 of the parallel-plate gap [Ω].
func (h Hull) Z0() float64 { return Eta0 * h.Gap / h.H }

func (h Hull) perSection() (l, c, dx float64) {
	dx = h.Circumference() / float64(h.N)
	return Mu0 * h.Gap / h.H * dx, Eps0 * h.H / h.Gap * dx, dx
}

// Freq is the ring resonance of mode n on the N-section ladder
// (lattice dispersion; → n·c/circumference for many sections).
func (h Hull) Freq() float64 {
	l, c, _ := h.perSection()
	return 2 / math.Sqrt(l*c) * math.Sin(math.Pi*float64(h.Mode)/float64(h.N)) / (2 * math.Pi)
}

// Rs is the wall surface resistance at f [Ω per square].
func (h Hull) Rs(f float64) float64 { return math.Sqrt(math.Pi * f * Mu0 * h.Rho) }

// Q of the ring from wall loss: β/(2α), α = Rs/(η0·d).
func (h Hull) Q() float64 {
	f := h.Freq()
	return 2 * math.Pi * f / C0 / (2 * h.Rs(f) / (Eta0 * h.Gap))
}

// PredictedVrf is the traveling-wave gap amplitude √(Z0·P·Q/(π·n)).
func (h Hull) PredictedVrf() float64 {
	return math.Sqrt(h.Z0() * h.Pin * h.Q() / (math.Pi * float64(h.Mode)))
}

func (h Hull) node(k int) string { return fmt.Sprintf("n%d", ((k%h.N)+h.N)%h.N) }

// circuit builds the ring with four radiators of amplitude amps (A, peak)
// and the bias feed.
func (h Hull) circuit(amps float64) (*circuit.Circuit, error) {
	l, c, _ := h.perSection()
	f := h.Freq()
	rsec := 2 * h.Rs(f) / h.H * h.Circumference() / float64(h.N)
	b := circuit.NewBuilder("double hull")
	for k := 0; k < h.N; k++ {
		mid := fmt.Sprintf("m%d", k)
		b.Inductor(fmt.Sprintf("L%d", k), h.node(k), mid, l)
		b.Resistor(fmt.Sprintf("R%d", k), mid, h.node(k+1), rsec)
		b.Capacitor(fmt.Sprintf("C%d", k), h.node(k), "0", c)
	}
	q := h.N / (4 * h.Mode)
	half := h.N / 2
	phB := -360.0 * float64(h.Mode) * float64(half) / float64(h.N)
	b2 := phB - 90
	if h.Standing {
		b2 = phB + 90
	}
	for i, p := range []struct {
		at int
		ph float64
	}{{0, 0}, {q, -90}, {half, phB}, {half + q, b2}} {
		b.Add(circuit.Element{Name: []string{"IA1", "IA2", "IB1", "IB2"}[i], Kind: circuit.CurrentSource,
			HasDC: true, ACMag: amps, ACPhase: math.Mod(p.ph, 360)}, "0", h.node(p.at))
	}
	// Bias: outer hull (the ring's conductor) to inner hull (ground) through
	// RFeed, far above Z0; at RF the supply is a short, so the feed is a
	// small load V_rf²/(2·RFeed) on the ring, booked in the ledger.
	b.VoltageSource("VDC", "bias", "0", h.Vdc).Resistor("RB", "bias", h.node(0), RFeed)
	return b.Build()
}

// RFeed is the bias feed resistance [Ω]; it charges the gap capacitance
// (≈ 1.4 nF for the example hull) in RFeed·C ≈ 14 ms. (Much higher values make the conductance
// fall below the solver's pivot threshold.)
const RFeed = 1e7

// Result is a solved hull.
type Result struct {
	F, Q, Z0       float64
	Vrf            []float64 // RF gap amplitude per node [V]
	VrfMax, VrfMin float64
	Vdc            float64 // bias on the outer hull from the operating point [V]
	PeakGap        float64 // V_dc + max V_rf [V]
	Pin, PLoss     float64 // radiators' power, wall loss [W]
	PFeed          float64 // RF power lost in the bias feed resistor [W]
	Stored         float64 // RF energy in the gap [J]
	EDC            float64 // bias energy ½·C_gap·V_dc² [J]
	Predicted      float64
}

// Solve runs the operating point (bias) and the AC steady state, scaling
// the radiators so they deliver h.Pin.
func (h Hull) Solve() (*Result, error) {
	run := func(amps float64) (*circuit.Circuit, *Result, error) {
		ckt, err := h.circuit(amps)
		if err != nil {
			return nil, nil, err
		}
		f := h.Freq()
		ac, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return nil, nil, err
		}
		r := &Result{F: f, Q: h.Q(), Z0: h.Z0(), Vrf: make([]float64, h.N), VrfMin: math.Inf(1), Predicted: h.PredictedVrf()}
		l, c, _ := h.perSection()
		rsec := 2 * h.Rs(f) / h.H * h.Circumference() / float64(h.N)
		for k := 0; k < h.N; k++ {
			v, _ := ac.NodeVoltage(h.node(k))
			il, _ := ac.ElementCurrent(fmt.Sprintf("L%d", k))
			r.Vrf[k] = cmplx.Abs(v[0])
			r.VrfMax, r.VrfMin = math.Max(r.VrfMax, r.Vrf[k]), math.Min(r.VrfMin, r.Vrf[k])
			r.PLoss += 0.5 * rsec * cmplx.Abs(il[0]) * cmplx.Abs(il[0])
			// Time-averaged stored energy: ¼C|V|² + ¼L|I|².
			r.Stored += 0.25*c*r.Vrf[k]*r.Vrf[k] + 0.25*l*cmplx.Abs(il[0])*cmplx.Abs(il[0])
		}
		v0, _ := ac.NodeVoltage(h.node(0))
		r.PFeed = 0.5 * cmplx.Abs(v0[0]) * cmplx.Abs(v0[0]) / RFeed
		for _, n := range []string{"IA1", "IA2", "IB1", "IB2"} {
			is, _ := ac.ElementCurrent(n)
			var at string
			for _, e := range ckt.Elements {
				if e.Name == n {
					at = ckt.NodeNames[e.Nodes[1]]
				}
			}
			v, _ := ac.NodeVoltage(at)
			r.Pin += 0.5 * real(v[0]*cmplx.Conj(is[0]))
		}
		return ckt, r, nil
	}
	_, r0, err := run(1e-3)
	if err != nil {
		return nil, err
	}
	if r0.Pin <= 0 {
		return nil, fmt.Errorf("radiators deliver no power (%g W)", r0.Pin)
	}
	ckt, r, err := run(1e-3 * math.Sqrt(h.Pin/r0.Pin))
	if err != nil {
		return nil, err
	}
	op, err := analysis.OP(ckt, analysis.Options{})
	if err != nil {
		return nil, err
	}
	r.Vdc = op.NodeVoltages[ckt.NodeIndex(h.node(0))]
	r.PeakGap = math.Abs(r.Vdc) + r.VrfMax
	cg := Eps0 * h.H * h.Circumference() / h.Gap
	r.EDC = 0.5 * cg * r.Vdc * r.Vdc
	return r, nil
}

// SEY, the copper yields and Multipactor live in em/multipactor.go.
type SEY = em.SEY

var (
	CopperAsReceived  = em.CopperAsReceived
	CopperConditioned = em.CopperConditioned
)

// Multipactor is em.Multipactor.
func Multipactor(d, f, vrf, vdc float64, sey SEY) (growth, eImpact float64) {
	return em.Multipactor(d, f, vrf, vdc, sey)
}

// RingUp is the worst growth the gap meets on the way to amplitude vrf:
// while the ring builds up (and at every point of a standing wave, from
// node to antinode) the gap sees every RF amplitude from 0 to vrf. It
// scans `steps` amplitudes evenly from vrf/steps to vrf at bias vdc and
// returns the worst growth and the amplitude where it occurs. The
// resonant bands are a few hundred volts wide at 10 cm and 38 MHz, so
// keep the step below about 25 V there.
func RingUp(d, f, vrf, vdc float64, sey SEY, steps int) (worst, at float64) {
	for i := 1; i <= steps; i++ {
		v := vrf * float64(i) / float64(steps)
		if g, _ := Multipactor(d, f, v, vdc, sey); g > worst {
			worst, at = g, v
		}
	}
	return worst, at
}

// ringSteps gives a ring-up scan step of about 20 V (at least 20 points).
func ringSteps(vrf float64) int { return max(20, int(math.Ceil(vrf/20))) }

// BiasGrid is the set of bias voltages MinBias tries [V].
var BiasGrid = []float64{0, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 7500, 10000}

// MinBias returns the smallest bias on BiasGrid that keeps the whole
// ring-up below growth 1, or −1 if none does. (A grid, not a bisection:
// suppression is not guaranteed to be monotonic in bias.)
func MinBias(d, f, vrf float64, sey SEY) float64 {
	for _, b := range BiasGrid {
		if g, _ := RingUp(d, f, vrf, b, sey, ringSteps(vrf)); g <= 1 {
			return b
		}
	}
	return -1
}

// KilpatrickField is the RF surface field [V/m] at which vacuum sparking
// becomes likely at frequency f: f[MHz] = 1.64·E²·exp(−8.5/E), E in MV/m
// (W. D. Kilpatrick, Rev. Sci. Instrum. 28, 824, 1957, in the form of
// T. J. Boyd, LANL report LA-UR-82-1166, 1982). A conservative RF limit;
// applied here to the peak (DC + RF) gap field.
func KilpatrickField(f float64) float64 {
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

// Problem evolves the gap, mode, RF power and bias for a target peak gap
// voltage, with at least minRF of it as RF amplitude, using the least RF
// power, keeping the whole ring-up out of the
// multipactor zones of as-received copper and the peak field below
// Kilpatrick's limit. The hull radius and height are fixed (the craft).
func Problem(target, minRF float64, base Hull) evo.Problem {
	hull := func(v []float64) Hull {
		h := base
		h.Gap, h.Mode, h.Pin, h.Vdc = v[0], int(math.Round(v[1])), v[2], v[3]
		h.N = 24 * h.Mode
		return h
	}
	solve := func(e *evo.Eval) (*Result, error) {
		r, err := e.Memo("solve", func() (any, error) { return hull(e.Values).Solve() })
		if err != nil {
			return nil, err
		}
		return r.(*Result), nil
	}
	return evo.Problem{
		Space: evo.Space{
			{Name: "gap", Min: 0.02, Max: 0.5},
			{Name: "mode", Min: 1, Max: 8},
			{Name: "P_in", Min: 10, Max: 1e5},
			{Name: "V_dc", Min: 10, Max: 2e4},
		},
		Build: func(v []float64) (*circuit.Circuit, error) { return hull(v).circuit(1e-3) },
		Objectives: []evo.Objective{
			{Name: "peak off target", Weight: 10, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := solve(e)
				if err != nil {
					return 0, err
				}
				return math.Abs(r.PeakGap-target) / target, nil
			}},
			{Name: "RF amplitude short of minimum", Weight: 10, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := solve(e)
				if err != nil {
					return 0, err
				}
				return math.Max(0, (minRF-r.VrfMax)/minRF), nil
			}},
			{Name: "RF power [decades]", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				return math.Log10(e.Values[2]), nil
			}},
			{Name: "multipactor growth over 1", Weight: 100, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := solve(e)
				if err != nil {
					return 0, err
				}
				g, _ := RingUp(e.Values[0], r.F, r.VrfMax, r.Vdc, CopperAsReceived, ringSteps(r.VrfMax))
				return math.Max(0, g-1), nil
			}},
			{Name: "field over Kilpatrick", Weight: 100, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				r, err := solve(e)
				if err != nil {
					return 0, err
				}
				return math.Max(0, r.PeakGap/e.Values[0]/KilpatrickField(r.F)-1), nil
			}},
		},
	}
}

func main() {
	target := flag.Float64("target", 20000, "target peak gap voltage for the GA [V]")
	pop := flag.Int("pop", 24, "GA population")
	gens := flag.Int("gens", 15, "GA generations")
	minRF := flag.Float64("minrf", 5000, "least RF amplitude in the GA design [V]")
	flag.Parse()
	h := Example()
	fmt.Printf("Double hull: radius %.1f m, height %.1f m, gap %.2f m, copper walls, mode %d, %d sections\n", h.R, h.H, h.Gap, h.Mode, h.N)
	for _, st := range []bool{false, true} {
		for _, vdc := range []float64{0, 2500} {
			hh := h
			hh.Standing, hh.Vdc = st, vdc
			r, err := hh.Solve()
			if err != nil {
				log.Fatal(err)
			}
			name := "traveling"
			if st {
				name = "standing"
			}
			fmt.Printf("\n%s, bias %.0f V: f = %.3f MHz, Q = %.0f, Z0 = %.2f Ω\n", name, vdc, r.F/1e6, r.Q, r.Z0)
			fmt.Printf("  radiators %.1f W = wall loss %.1f W + bias feed %.3f W; RF energy in gap %.1f mJ; bias energy %.3f mJ (stored once)\n",
				r.Pin, r.PLoss, r.PFeed, r.Stored*1e3, r.EDC*1e3)
			fmt.Printf("  RF gap amplitude %.0f–%.0f V (formula, traveling: %.0f V); bias %.0f V; peak gap voltage %.0f V\n",
				r.VrfMin, r.VrfMax, r.Predicted, r.Vdc, r.PeakGap)
		}
	}

	fmt.Printf("\nMultipactor in the vacuum gap (d = %.2f m, f = %.2f MHz; f·d = %.2f MHz·m)\n", h.Gap, h.Freq()/1e6, h.Freq()*h.Gap/1e6)
	fmt.Printf("  mean secondary yield per impact on the worst phase (> 1 = discharge grows)\n")
	fmt.Printf("  %8s", "V_rf")
	biases := []float64{0, 500, 2500, 5000}
	for _, b := range biases {
		fmt.Printf("  %16s", fmt.Sprintf("bias %.0f V", b))
	}
	fmt.Println()
	for _, vrf := range []float64{100, 300, 1000, 1700, 3400, 10000} {
		fmt.Printf("  %6.0f V", vrf)
		for _, b := range biases {
			g, e := Multipactor(h.Gap, h.Freq(), vrf, b, CopperAsReceived)
			mark := " "
			if g > 1 {
				mark = "*"
			}
			fmt.Printf("  %6.2f @%5.0f eV%s", g, e, mark)
		}
		fmt.Println()
	}
	fmt.Println("  (* = susceptible, as-received copper; Baglin et al. 2002)")

	vrf := h.PredictedVrf()
	fmt.Printf("\nRing-up to %.0f V (traveling) and %.0f V (standing antinode): every amplitude below is passed through\n", vrf, vrf*math.Sqrt2)
	for _, sey := range []SEY{CopperAsReceived, CopperConditioned} {
		for _, v := range []float64{vrf, vrf * math.Sqrt2} {
			g0, a0 := RingUp(h.Gap, h.Freq(), v, 0, sey, ringSteps(v))
			fmt.Printf("  %-20s to %5.0f V, no bias: worst growth %.2f at %.0f V; least bias that clears it: %.0f V\n",
				sey.Name, v, g0, a0, MinBias(h.Gap, h.Freq(), v, sey))
		}
	}
	fmt.Printf("  Kilpatrick RF limit at %.1f MHz: %.1f MV/m (the 2.5 kV-biased standing wave peaks at %.3f MV/m)\n",
		h.Freq()/1e6, KilpatrickField(h.Freq())/1e6, (2500+vrf*math.Sqrt2)/h.Gap/1e6)

	fmt.Printf("\nEvolving gap, mode, RF power and bias for a %.0f V peak (≥ %.0f V of it RF) with the least RF power, multipactor-free ring-up, under Kilpatrick\n", *target, *minRF)
	s, err := evo.Run(context.Background(), Problem(*target, *minRF, h), evo.Options{Population: *pop, Generations: *gens, Seed: 1, Timeout: 120e9}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	v := s.Best.Values
	best := h
	best.Gap, best.Mode, best.Pin, best.Vdc = v[0], int(math.Round(v[1])), v[2], v[3]
	best.N = 24 * best.Mode
	r, err := best.Solve()
	if err != nil {
		log.Fatal(err)
	}
	g, at := RingUp(best.Gap, r.F, r.VrfMax, r.Vdc, CopperAsReceived, ringSteps(r.VrfMax))
	fmt.Printf("  best of %d evaluations: gap %.3f m, mode %d (%.2f MHz, Q %.0f), RF %.0f W, bias %.0f V\n",
		s.Evaluations, best.Gap, best.Mode, r.F/1e6, r.Q, r.Pin, r.Vdc)
	fmt.Printf("  RF %.0f V + bias %.0f V = peak %.0f V; ring-up worst growth %.2f (at %.0f V); field %.3f MV/m vs Kilpatrick %.1f\n",
		r.VrfMax, r.Vdc, r.PeakGap, g, at, r.PeakGap/best.Gap/1e6, KilpatrickField(r.F)/1e6)
	fmt.Printf("  ledger: radiators %.1f W = wall loss %.1f W + bias feed %.3f W; the bias supply delivers C·V² = %.2f J once: %.2f J stored, %.2f J lost in the feed resistor; no steady DC current (no gap leakage modeled)\n",
		r.Pin, r.PLoss, r.PFeed, 2*r.EDC, r.EDC, r.EDC)
}
