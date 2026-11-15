// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command ece_rlc runs one series RLC circuit three ways: plain (Maxwell),
// with the ECE spin-connection terms (.ece form=ece) and with the ECE2 terms
// (.ece form=ece2). These are kronspice's hypothetical H48 option
// (circuit/ece.go, device/ece.go; CLAUDE.md §1 register entry H48 and §3c).
//
// The circuit: 1 V source, R = 1 Ω, L = 10 mH, C = 84.6346 µF
// (f0 = 173.000 Hz, Q0 = ω0L/R = 10.87). Two terms are switched on:
//
//   - torsion: a series resistance s·c·ω₀·L in the coil, s = +1 (ECE) or
//     −2 (ECE2). Here ω₀ = 8.339×10⁻⁸ m⁻¹, so c·ω₀·L = 0.25 Ω;
//   - scalar: the capacitor's top terminal enters its body as (1 + s·κ)·φ,
//     here with κ = 0.05.
//
// Each case is checked against closed forms, Q = ω0L/(R + s·c·ω₀·L) and
// f_r = f0/√(1 + s·κ), and against a twin built from ordinary elements.
// The ledger is checked in AC and in a 0.3 s transient. Power the terms
// deliver or absorb is booked as model-supplied, never as inflow.
package main

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// Circuit values.
const (
	R  = 1.0
	L  = 10e-3
	C  = 84.6346e-6
	c0 = 299792458.0
	// W0 is the torsion spin connection: c·W0·L = 0.25 Ω.
	W0 = 0.25 / (c0 * L)
	// Kappa is the scalar term on the capacitor's top terminal.
	Kappa = 0.05
)

// F0 is the plain circuit's resonance [Hz].
var F0 = 1 / (2 * math.Pi * math.Sqrt(L*C))

// Case is one netlist.
type Case struct {
	Name   string
	Form   string  // "", "ece" or "ece2"
	S      float64 // 0 (plain), +1 (ECE), −2 (ECE2)
	W0     float64 // torsion [1/m]
	KappaA float64 // scalar term on the capacitor's top terminal
}

// Cases are the runs, in print order.
var Cases = []Case{
	{Name: "plain (Maxwell)"},
	{Name: "ECE torsion", Form: "ece", S: 1, W0: W0},
	{Name: "ECE2 torsion", Form: "ece2", S: -2, W0: W0},
	{Name: "ECE scalar", Form: "ece", S: 1, KappaA: Kappa},
	{Name: "ECE2 scalar", Form: "ece2", S: -2, KappaA: Kappa},
}

// Netlist is a case's kronspice netlist.
func (k Case) Netlist() string {
	var b strings.Builder
	fmt.Fprintf(&b, "series RLC, %s\n", k.Name)
	if k.Form != "" {
		fmt.Fprintf(&b, ".ece form=%s w0=%.17g\n", k.Form, k.W0)
	}
	b.WriteString("V1 in 0 DC 0 AC 1 SIN(0 1 173)\nR1 in a 1\nL1 a b 10m\n")
	if k.KappaA != 0 {
		fmt.Fprintf(&b, "C1 b 0 84.6346u kappa_a=%g\n", k.KappaA)
	} else {
		b.WriteString("C1 b 0 84.6346u\n")
	}
	b.WriteString(".end\n")
	return b.String()
}

// Twin is the same case from ordinary elements: the torsion as a series
// resistor (negative in ECE2), the scalar term as a voltage-controlled
// source (E) in the capacitor's lead.
func (k Case) Twin() string {
	rt := k.S * c0 * k.W0 * L
	var b strings.Builder
	fmt.Fprintf(&b, "series RLC twin, %s\nV1 in 0 DC 0 AC 1 SIN(0 1 173)\nR1 in a %.17g\nL1 a b 10m\n", k.Name, R+rt)
	if k.KappaA != 0 {
		// node x = (1 + s·κ)·φ_b drives the capacitor body
		fmt.Fprintf(&b, "E1 x b b 0 %.17g\nC1 x 0 84.6346u\n", k.S*k.KappaA)
	} else {
		b.WriteString("C1 b 0 84.6346u\n")
	}
	b.WriteString(".end\n")
	return b.String()
}

// Result is one case's numbers.
type Result struct {
	Case
	FPeak, Q, IPeak float64 // AC: resonance [Hz], Q, peak |I| [A]
	FWant, QWant    float64 // closed forms
	TwinDiff        float64 // max |I − I_twin| / max |I_twin| over the sweep
	ACResidual      float64 // worst relative AC ledger residual
	// At the AC point nearest f_r [W]:
	PSource, PHeat, PModel float64
	// Transient, 0.3 s of the 1 V sine at 173 Hz [J]:
	EIn, EOut, EHeat, EModel, EStored, Mismatch float64
}

func build(text string) (*circuit.Circuit, error) {
	d, err := netlist.ParseString(text)
	if err != nil {
		return nil, err
	}
	return circuit.FromDeck(d)
}

func sweep(text string) ([]float64, []complex128, *output.ACResult, error) {
	ck, err := build(text)
	if err != nil {
		return nil, nil, nil, err
	}
	r, err := analysis.AC(ck, analysis.ACSpec{Variation: "lin", Points: 4001, FStart: 120, FStop: 230}, analysis.Options{})
	if err != nil {
		return nil, nil, nil, err
	}
	i, ok := r.ElementCurrent("V1")
	if !ok {
		return nil, nil, nil, fmt.Errorf("no current for V1")
	}
	return r.Freq, i, r, nil
}

// Run computes one case.
func Run(k Case) (Result, error) {
	res := Result{Case: k, FWant: F0, QWant: 2 * math.Pi * F0 * L / R}
	if k.Form != "" {
		res.QWant = 2 * math.Pi * F0 * L / (R + k.S*c0*k.W0*L)
		res.FWant = F0 / math.Sqrt(1+k.S*k.KappaA)
		if k.KappaA != 0 {
			// The scalar term scales the capacitance: C(1 + sκ), so Q = ω_rL/R.
			res.QWant = 2 * math.Pi * res.FWant * L / R
		}
	}
	f, i, r, err := sweep(k.Netlist())
	if err != nil {
		return res, fmt.Errorf("%s: %w", k.Name, err)
	}
	_, it, _, err := sweep(k.Twin())
	if err != nil {
		return res, fmt.Errorf("%s twin: %w", k.Name, err)
	}
	m := make([]float64, len(i))
	mx, d := 0.0, 0.0
	for n := range i {
		m[n] = cmplx.Abs(i[n])
		mx = math.Max(mx, cmplx.Abs(it[n]))
		d = math.Max(d, cmplx.Abs(i[n]-it[n]))
	}
	res.TwinDiff = d / mx
	var kp int
	res.FPeak, res.Q, kp = peak(f, m)
	res.IPeak = m[kp]
	res.ACResidual = r.Ledger.MaxRelativeResidual
	res.PSource, res.PHeat, res.PModel = r.Ledger.Inflow[kp], r.Ledger.Dissipation[kp], r.Ledger.ModelSupplied[kp]

	ck, err := build(k.Netlist())
	if err != nil {
		return res, err
	}
	wf, err := analysis.Tran(ck, analysis.TranSpec{Step: 20e-6, Stop: 0.3, MaxStep: 20e-6}, analysis.TranOptions{})
	if err != nil {
		return res, fmt.Errorf("%s transient: %w", k.Name, err)
	}
	l := wf.Ledger
	res.EIn, res.EOut, res.EHeat, res.EModel = l.EnergyIn, l.EnergyOut, l.EnergyDissipated, l.EnergyModelSupplied
	res.EStored, res.Mismatch = l.StoredEnd-l.StoredStart, l.IntegrationMismatch
	return res, nil
}

// peak returns the frequency of the largest magnitude, Q from the
// interpolated −3 dB points, and the index of the peak.
func peak(f, m []float64) (fp, q float64, k int) {
	for n := range m {
		if m[n] > m[k] {
			k = n
		}
	}
	h := m[k] / math.Sqrt2
	i, j := k, k
	for i > 0 && m[i] >= h {
		i--
	}
	for j < len(m)-1 && m[j] >= h {
		j++
	}
	fl := f[i] + (h-m[i])*(f[i+1]-f[i])/(m[i+1]-m[i])
	fh := f[j-1] + (h-m[j-1])*(f[j]-f[j-1])/(m[j]-m[j-1])
	return f[k], f[k] / (fh - fl), k
}

func main() {
	data := flag.String("data", "", "also write the report's CSV files to this directory")
	flag.Parse()
	if *data != "" {
		if err := WriteData(*data); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
	}
	fmt.Printf("Series RLC: R %g Ω, L %g mH, C %g µF, f0 %.3f Hz, Q0 %.4f\n", R, L*1e3, C*1e6, F0, 2*math.Pi*F0*L/R)
	fmt.Printf("ω₀ = %.4e 1/m (c·ω₀·L = %.2f Ω), κ = %g on the capacitor's top terminal\n\n", W0, c0*W0*L, Kappa)
	fmt.Println("AC sweep 120–230 Hz (4001 points), 1 V; power at the peak")
	fmt.Println("  case               f_r (Hz)  closed form   Q        closed form  peak |I| (A)  vs twin   source (W)  heat (W)  model-supplied (W)  residual")
	var rs []Result
	for _, k := range Cases {
		r, err := Run(k)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		rs = append(rs, r)
		fmt.Printf("  %-17s  %8.3f  %8.3f     %7.4f  %7.4f      %.4f        %.0e    %.4f      %.4f    %+.4f             %.0e\n",
			r.Name, r.FPeak, r.FWant, r.Q, r.QWant, r.IPeak, r.TwinDiff, r.PSource, r.PHeat, r.PModel, r.ACResidual)
	}
	fmt.Println("\nTransient 0–0.3 s, 1 V sine at 173 Hz (trapezoidal, 20 µs steps); energies [J]")
	fmt.Println("  case               source in  back to source  heat       model-supplied  Δstored     mismatch")
	for _, r := range rs {
		fmt.Printf("  %-17s  %8.5f   %8.5f        %8.5f   %+9.5f       %+.2e   %.1e\n", r.Name, r.EIn, r.EOut, r.EHeat, r.EModel, r.EStored, r.Mismatch)
	}
}
