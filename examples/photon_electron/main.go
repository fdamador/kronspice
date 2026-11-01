// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command photon_electron puts the toroidal-photon electron of J. G.
// Williamson and M. B. van der Mark (hypothesis H23) on the electron test
// bench, and builds its central topological idea as a circuit. Source:
// J. G. Williamson, M. B. van der Mark, "Is the electron a photon with
// toroidal topology?", Annales de la Fondation Louis de Broglie 22(2),
// 133–160 (1997), https://fondationlouisdebroglie.org/AFLB-222/MARK.TEX2.pdf.
//
// The model: one photon of wavelength λ confined to a closed path that
// goes twice round a torus (the paper assumes the confinement; it adds no
// field equation). With λ the Compton wavelength h/(m_e c) it states:
//
//	q  = (1/2π)·√(3ε0ħc)                 (eq. 5, from the field in a sphere of diameter λ/2)
//	α′ = q²/(4πε0ħc) = (q/e)²α           (eq. 8)
//	g  = 2(1 + α′/2π)                    (eq. 12; eq. 13 compares QED's 2(1 + α/2π))
//	spin ½ from the double loop: the field returns to itself only after 720°
//
// Part 1 compares each prediction with the measured value (CODATA 2018).
//
// Part 2 is the circuit: a two-wire ring transmission line (an LC ladder
// closed on itself, Kron's closed path; as in examples/ring_towers), once
// plain and once with its two conductors crossed at one point. On the
// plain ring a wave returns to itself after one trip: its differential
// modes have a whole number of wavelengths round the ring. On the crossed
// ring the differential voltage changes sign at the crossing, so a wave
// comes back inverted after one trip and only matches itself after two —
// the 720° property of the paper's double loop. Its modes have half-integer
// numbers of wavelengths round the ring; the lowest has a wavelength of
// two circumferences. Both rings are solved in kronspice: an AC sweep finds
// the modes (twin of the ladder's analytic dispersion), and a transient
// sends a pulse round and shows its sign after each trip, with the energy
// ledger.
//
// Part 3 (field.go) tests the paper's central assumption, that the light
// holds itself in that loop, with the 3-D field solver (em/fdtd.go,
// validated on the exact "flying doughnut" toroidal pulse): one wavelength
// of light is set on the double loop in empty space and followed.
//
//	go run ./examples/photon_electron
package main

import (
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

// Measured constants, SI (CODATA 2018, NIST).
const (
	E       = 1.602176634e-19  // C (exact)
	Hbar    = 1.054571817e-34  // J s
	Eps0    = 8.8541878128e-12 // F/m
	C0      = 299792458.0      // m/s (exact)
	Me      = 9.1093837015e-31 // kg
	AlphaFS = 7.2973525693e-3
	GE      = 2.00231930436256 // |g| of the free electron
	RClass  = 2.8179403262e-15 // m, classical electron radius
)

// Compton is the electron's Compton wavelength h/(m_e c) [m].
func Compton() float64 { return 2 * math.Pi * Hbar / (Me * C0) }

// Charge is the model's charge, eq. 5 [C].
func Charge() float64 { return math.Sqrt(3*Eps0*Hbar*C0) / (2 * math.Pi) }

// AlphaPrime is the model's coupling, eq. 8.
func AlphaPrime() float64 { q := Charge() / E; return q * q * AlphaFS }

// GFactor is the model's g, eq. 12.
func GFactor() float64 { return 2 * (1 + AlphaPrime()/(2*math.Pi)) }

// Ring is a two-wire ring line of N sections: series inductance L in each
// wire per section, capacitance Cd between the wires and Cg from each wire
// to ground at each junction, a small series resistance R per wire.
type Ring struct {
	N         int
	L, Cd, Cg float64
	R         float64
	Crossed   bool
}

// DefaultRing: 200 sections; L = 1 µH, Cd = 1 nF, Cg = 0.1 nF, R = 0.01 Ω.
func DefaultRing(crossed bool) Ring {
	return Ring{N: 200, L: 1e-6, Cd: 1e-9, Cg: 0.1e-9, R: 0.01, Crossed: crossed}
}

// Delay is the time a differential wave takes round the ring once (long
// wavelengths): N·√(2L·(Cd + Cg/2)).
func (r Ring) Delay() float64 { return float64(r.N) * math.Sqrt(2*r.L*(r.Cd+r.Cg/2)) }

// ModeFrequency is the differential mode m of the ladder [Hz]: whole m for
// the plain ring, m + ½ for the crossed one (the antiperiodic boundary),
// ω = 2/√(2L(Cd + Cg/2))·sin(πk/N), k = m or m + ½.
func (r Ring) ModeFrequency(m int) float64 {
	k := float64(m)
	if r.Crossed {
		k += 0.5
	}
	return 2 / math.Sqrt(2*r.L*(r.Cd+r.Cg/2)) * math.Sin(math.Pi*k/float64(r.N)) / (2 * math.Pi)
}

// Build makes the ring with a source element between a0 and b0 (an AC
// current or a current pulse). The last section closes the ring onto
// (a0, b0), or onto (b0, a0) when crossed.
func (r Ring) Build(src circuit.Element) (*circuit.Circuit, error) {
	b := circuit.NewBuilder(map[bool]string{false: "plain two-wire ring", true: "crossed two-wire ring"}[r.Crossed])
	a := func(k int) string { return "a" + strconv.Itoa(k%r.N) }
	bb := func(k int) string { return "b" + strconv.Itoa(k%r.N) }
	for k := 0; k < r.N; k++ {
		na, nb := a(k+1), bb(k+1)
		if k == r.N-1 && r.Crossed {
			na, nb = bb(0), a(0)
		}
		b.Inductor("LA"+strconv.Itoa(k), a(k), "ma"+strconv.Itoa(k), r.L)
		b.Resistor("RA"+strconv.Itoa(k), "ma"+strconv.Itoa(k), na, r.R)
		b.Inductor("LB"+strconv.Itoa(k), bb(k), "mb"+strconv.Itoa(k), r.L)
		b.Resistor("RB"+strconv.Itoa(k), "mb"+strconv.Itoa(k), nb, r.R)
		b.Capacitor("CD"+strconv.Itoa(k), a(k), bb(k), r.Cd)
		b.Capacitor("CGA"+strconv.Itoa(k), a(k), "0", r.Cg)
		b.Capacitor("CGB"+strconv.Itoa(k), bb(k), "0", r.Cg)
	}
	// The capacitors give the wires no DC path to ground; a 1 MΩ bleeder on
	// each wire fixes the operating point and does nothing at the ring's
	// frequencies (it is 10⁴ times the line impedance).
	b.Resistor("RBA", a(0), "0", 1e6)
	b.Resistor("RBB", bb(0), "0", 1e6)
	src.Name = "IS"
	b.Add(src, "b0", "a0") // current into a0, out of b0
	return b.BuildTransient()
}

// Modes sweeps the differential impedance seen at (a0, b0) and returns
// the frequencies of its first n peaks.
func (r Ring) Modes(n int, fStop float64) ([]float64, error) {
	ckt, err := r.Build(circuit.Element{Kind: circuit.CurrentSource, ACMag: 1})
	if err != nil {
		return nil, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 20001, FStart: fStop / 20000, FStop: fStop}, analysis.Options{})
	if err != nil {
		return nil, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	z := func(i int) float64 { return cmplx.Abs(va[i] - vb[i]) }
	var peaks []float64
	for i := 1; i+1 < len(res.Freq) && len(peaks) < n; i++ {
		if z(i) > z(i-1) && z(i) >= z(i+1) {
			// parabolic refinement on the three points
			y0, y1, y2 := z(i-1), z(i), z(i+1)
			d := 0.5 * (y0 - y2) / (y0 - 2*y1 + y2)
			peaks = append(peaks, res.Freq[i]+d*(res.Freq[1]-res.Freq[0]))
		}
	}
	return peaks, nil
}

// Trip is the pulse seen back at the start after a number of trips.
type Trip struct {
	N        int
	Time     float64 // arrival time [s]
	Vdiff    float64 // differential voltage at (a0, b0) at the arrival [V]
	Relative float64 // as a fraction of the launched pulse
}

// PulseRun is a pulse sent round the ring.
type PulseRun struct {
	Trips              []Trip
	Launch             float64 // peak differential voltage at launch [V]
	EIn, EDiss, EStore float64
	Mismatch           float64
	Time, V            []float64 // differential voltage at (a0, b0)
}

// Pulse injects a short differential current pulse at (a0, b0) and
// follows it for `trips` round trips. A pulse launched at one point runs
// both ways round the ring; the two halves meet back at the start after
// each full trip.
func (r Ring) Pulse(trips int) (PulseRun, error) {
	td := r.Delay()
	w := td / 40 // pulse length: a twentieth of the ring
	src := circuit.Element{Kind: circuit.CurrentSource, Wave: circuit.Waveform{Kind: circuit.WavePulse,
		Args: []float64{0, 0.1, 0, w / 4, w / 4, w / 2, 0}}}
	ckt, err := r.Build(src)
	if err != nil {
		return PulseRun{}, err
	}
	stop := (float64(trips) + 0.3) * td
	step := td / 2000
	res, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return PulseRun{}, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	run := PulseRun{EIn: res.Ledger.EnergyIn, EDiss: res.Ledger.EnergyDissipated,
		EStore: res.Ledger.StoredEnd - res.Ledger.StoredStart, Mismatch: res.Ledger.IntegrationMismatch, Time: res.Time}
	for i := range res.Time {
		run.V = append(run.V, va[i]-vb[i])
	}
	// Launch peak in the first quarter trip; each return within ±¼ trip of n·td.
	peakIn := func(t0, t1 float64) (float64, float64) {
		best, at := 0.0, 0.0
		for i, t := range res.Time {
			if t >= t0 && t <= t1 && math.Abs(run.V[i]) > math.Abs(best) {
				best, at = run.V[i], t
			}
		}
		return best, at
	}
	run.Launch, _ = peakIn(0, 0.25*td)
	for n := 1; n <= trips; n++ {
		v, t := peakIn((float64(n)-0.25)*td, (float64(n)+0.25)*td)
		run.Trips = append(run.Trips, Trip{N: n, Time: t, Vdiff: v, Relative: v / run.Launch})
	}
	return run, nil
}

// writeTrace writes time [µs] and differential voltage [V] for viz/render.py.
func writeTrace(dir, name string, p PulseRun) error {
	f, err := os.Create(filepath.Join(dir, name+"_trace.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "time_us,vdiff")
	for i := range p.Time {
		fmt.Fprintf(f, "%.6g,%.6g\n", p.Time[i]*1e6, p.V[i])
	}
	return nil
}

func main() {
	viz := flag.String("viz", "", "directory for the round-trip traces (viz/render.py)")
	flag.Parse()
	fmt.Println("Part 1 — the toroidal-photon electron (Williamson & van der Mark 1997) against measurement (CODATA 2018)")
	q, g := Charge(), GFactor()
	aModel, aMeas, aQED1 := (g-2)/2, (GE-2)/2, AlphaFS/(2*math.Pi)
	fmt.Printf("   %-36s %20s %22s\n", "quantity", "model", "measured")
	fmt.Printf("   %-36s %17.6g e %20s\n", "charge (eq. 5)", q/E, "1 e (exact)")
	fmt.Printf("   %-36s %20.6g %22.6g\n", "coupling α′ (eq. 8) vs α", AlphaPrime(), AlphaFS)
	fmt.Printf("   %-36s %20.9f %22.11f\n", "magnetic anomaly (g − 2)/2", aModel, aMeas)
	fmt.Printf("   %-36s %20.9f %22s\n", "  with q set to e (= α/2π, QED 1st)", aQED1, fmt.Sprintf("off by %.2g", (aQED1-aMeas)/aMeas))
	fmt.Printf("   %-36s %20s %22s\n", "spin", "1/2 (double loop)", "1/2")
	fmt.Printf("   %-36s %18.4g m %20s\n", "size: Compton wavelength λ", Compton(), "point-like to ~1e-18 m")
	fmt.Printf("   %-36s %18.4g m %20s\n", "  loop radius λ/4π", Compton()/(4*math.Pi), "")
	fmt.Printf("   model anomaly is %.1f %% below the measured one; the charge is %.1f %% below e.\n", 100*(aMeas-aModel)/aMeas, 100*(1-q/E))
	fmt.Println("   Standard counterpart: the Dirac electron with QED corrections, which gives the measured values.")

	fmt.Println("\nPart 2 — a closed path as a circuit: plain and crossed two-wire ring lines")
	for _, crossed := range []bool{false, true} {
		r := DefaultRing(crossed)
		name := map[bool]string{false: "plain ring", true: "crossed ring (conductors swapped once)"}[crossed]
		fmt.Printf("\n   == %s: %d sections, one trip %.3g µs ==\n", name, r.N, r.Delay()*1e6)
		fStop := 1.05 * r.ModeFrequency(4)
		got, err := r.Modes(4, fStop)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Printf("   %6s %16s %16s %10s %24s\n", "mode", "AC peak [Hz]", "ladder theory [Hz]", "diff", "wavelengths round ring")
		for m, f := range got {
			mm := m + 1
			if crossed {
				mm = m
			}
			k := float64(mm)
			if crossed {
				k += 0.5
			}
			th := r.ModeFrequency(mm)
			fmt.Printf("   %6d %16.6g %16.6g %9.2e %24.1f\n", m+1, f, th, (f-th)/th, k)
		}
		p, err := r.Pulse(4)
		if err != nil {
			log.Fatal(err)
		}
		if *viz != "" {
			if err := writeTrace(*viz, map[bool]string{false: "plain", true: "crossed"}[crossed], p); err != nil {
				log.Fatal(err)
			}
		}
		fmt.Printf("   pulse launched at %.3g V (differential); back at the start after each trip:\n", p.Launch)
		for _, t := range p.Trips {
			fmt.Printf("     trip %d at %6.3f µs: %+8.3g V (%+.2f of launch)\n", t.N, t.Time*1e6, t.Vdiff, t.Relative)
		}
		fmt.Printf("   ledger: in %.4g J = heat %.4g J + stored %.4g J (+ time-stepping error %.2g J)\n", p.EIn, p.EDiss, p.EStore, p.Mismatch)
	}
	fmt.Println("\nOn the plain ring the pulse comes back the same way up after every trip; on the crossed ring it comes back")
	fmt.Println("inverted after one trip and upright after two — the field returns to itself only after a 720° double loop.")

	fmt.Println("\nPart 3 — light set on the double loop, left alone in empty space (em/fdtd.go, linear Maxwell)")
	k := DefaultKnot()
	k.Periods = 2
	r, err := k.Run()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   ring radius R = λ/4π = %.4g m (λ = %.3g m), path length %.3f λ, grid %.3g m (R/%.0f), region ±%.0fR\n",
		r.R, k.Lambda, r.PathLength/k.Lambda, r.Dx, k.CellsR, k.HalfR)
	fmt.Printf("   %10s %22s %20s %18s\n", "time (T)", "energy within 2R", "in the region", "left the region")
	for _, s := range r.Samples {
		fmt.Printf("   %10.3g %21.4g%% %19.4g%% %17.4g%%\n", s.T, 100*s.Inside, 100*s.InRegion, 100*s.Absorbed)
	}
	fmt.Printf("   ledger residual %.2g of the starting energy; div B %.2g (round-off)\n", r.Residual, r.DivB)
	fmt.Println("   Standard prediction (linear Maxwell): nothing holds light in a loop, so it leaves at the speed of light.")
	fmt.Println("   The model needs 100 % to stay for ever; it would need a confining mechanism the paper does not give.")
	if *viz != "" {
		if err := writeKnot(*viz, r); err != nil {
			log.Fatal(err)
		}
	}
	if err := part4(*viz); err != nil {
		log.Fatal(err)
	}
	if err := part5(*viz); err != nil {
		log.Fatal(err)
	}
	if err := part6(*viz); err != nil {
		log.Fatal(err)
	}
}

// writeKnot writes the energy-versus-time samples for viz/render.py.
func writeKnot(dir string, r KnotRun) error {
	f, err := os.Create(filepath.Join(dir, "knot_energy.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, "periods,inside_2R,in_region,absorbed")
	for _, s := range r.Samples {
		fmt.Fprintf(f, "%g,%g,%g,%g\n", s.T, s.Inside, s.InRegion, s.Absorbed)
	}
	return nil
}
