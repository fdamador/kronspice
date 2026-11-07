// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Layer 7: a nonlinear ionosphere in the Earth–ionosphere cavity. Two
// transmitters drive the cavity ladder (layer 3, fitted to mode 1); the
// ionosphere's loss, a shunt conductance G = 2α·C at every band, is made to
// depend on the field by electron heating:
//
//	G(E) = G0 / (1 + (E/Ep)²)
//
// A strong field heats the electrons, their collision frequency rises and
// the conductivity σ = n·e²/(m·ν) falls (A. V. Gurevich, Nonlinear Phenomena
// in the Ionosphere, Springer, 1978, ch. 2: in the collision-dominated
// steady state Te/T = 1 + (E/Ep)², with ν ∝ Te). This is the mechanism of
// the Luxembourg effect, a strong station modulating the absorption of
// another (B. D. H. Tellegen, Nature 131, 840, 1933), and of HF heating
// experiments (HAARP, EISCAT).
//
// The plasma field Ep = √(3·k·T·m·δ·(ω² + ν²))/e is about 0.02 V/m for
// T = 200 K, δ = 10⁻³ (fraction of energy lost per collision), ν = 10⁶ s⁻¹
// (night-time D region near 70 km) — declared. The field at a band is
// taken as V/h_e, the cavity's vertical field, an UPPER BOUND on the field
// that reaches the plasma (the conducting ionosphere screens most of it),
// so the heating voltage is Vh = Ep·h_e.
//
// The two transmitters are current sources at the source band, at the
// measured first two Schumann frequencies, 7.8 and 14.3 Hz. A symmetric
// nonlinearity gives odd-order products; the third-order "four-wave" ones
// are 2f1 − f2 = 1.3 Hz and 2f2 − f1 = 20.8 Hz — the measured third
// Schumann frequency.

package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// NLIonosphere declares the run.
type NLIonosphere struct {
	Ep      float64 // plasma heating field [V/m]
	F0      float64 // frequency grid [Hz]: both drives are multiples of it
	M1, M2  int     // drives at M1·F0, M2·F0
	Drive   float64 // each transmitter's current [A] (amplitude)
	Station int     // receiving band (0 = the source band)
	Warm    float64 // warm-up [s]
	Dt      float64 // step [s]
	Linear  bool    // keep G constant (the standard linear cavity)
	Bands   int     // ladder bands (layer 3 uses 181)
}

// DefaultNLIonosphere: 7.8 and 14.3 Hz on a 0.1 Hz grid, receiver 30° from
// the source, drive chosen so the source band reaches about Vh.
func DefaultNLIonosphere() NLIonosphere {
	return NLIonosphere{Ep: 0.02, F0: 0.1, M1: 78, M2: 143, Drive: 0, Station: 30, Warm: 3, Dt: 1e-3, Bands: 181}
}

// cavityNL is the lossy cavity of layer 3 (mode 1 fitted, Q = 5) at the
// declared band count.
func (p NLIonosphere) cavity() Cavity {
	c := FitLossy(MeasuredSchumann[0], 5)
	c.Bands = p.Bands
	return c
}

// Vh is the heating voltage Ep·h_e.
func (p NLIonosphere) Vh() float64 { return p.Ep * p.cavity().He }

// Circuit builds the ladder with the two transmitters.
func (p NLIonosphere) Circuit() (*circuit.Circuit, error) {
	c := p.cavity()
	C, L := c.Elements()
	vh := p.Vh()
	b := circuit.NewBuilder("Earth-ionosphere cavity, nonlinear ionosphere")
	nd := func(k int) string { return fmt.Sprintf("b%d", k) }
	for k := range C {
		b.Capacitor(fmt.Sprintf("C%d", k), nd(k), "0", C[k])
		g := 2 * c.Alpha * C[k]
		if p.Linear {
			b.Resistor(fmt.Sprintf("G%d", k), nd(k), "0", 1/g)
		} else {
			// i = G0·v/(1 + (v/Vh)²) on ±4 Vh, 801 points.
			var cv []float64
			for j := -400; j <= 400; j++ {
				v := 4 * vh * float64(j) / 400
				cv = append(cv, v, g*v/(1+(v/vh)*(v/vh)))
			}
			b.NonlinearVC(fmt.Sprintf("G%d", k), nd(k), "0", cv)
		}
		if k < len(L) {
			b.Inductor(fmt.Sprintf("L%d", k), nd(k), nd(k+1), L[k])
		}
	}
	for i, m := range []int{p.M1, p.M2} {
		f := float64(m) * p.F0
		b.Add(circuit.Element{Name: fmt.Sprintf("ITX%d", i+1), Kind: circuit.CurrentSource,
			Wave: circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, p.Drive, f, 0, 0, 0}}}, "0", nd(0))
	}
	return b.Build()
}

// NLRun is one run: the window's time, source-band and station voltages,
// their spectra on the grid, and the ledger.
type NLRun struct {
	T, Source, Station []float64
	Lines              map[int]complex128 // station spectrum: bin k (k·F0 Hz) → coefficient
	SrcLines           map[int]complex128
	SourcePeak         float64 // max |V| at the source band [V]
	Ledger             float64 // relative integration mismatch
	DissFrac           float64 // energy dissipated / energy in
}

// Run simulates warm-up plus one window 1/F0.
func (p NLIonosphere) Run() (NLRun, error) {
	var r NLRun
	c, err := p.Circuit()
	if err != nil {
		return r, err
	}
	T := 1 / p.F0
	stop := p.Warm + T
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: p.Dt, Stop: stop, MaxStep: p.Dt}, analysis.TranOptions{MaxSteps: 2_000_000})
	if err != nil {
		return r, err
	}
	l := wf.Ledger
	r.Ledger = math.Abs(l.IntegrationMismatch) / math.Max(l.EnergyIn, 1e-300)
	r.DissFrac = l.EnergyDissipated / math.Max(l.EnergyIn, 1e-300)
	vs, _ := wf.NodeVoltage("b0")
	vr, _ := wf.NodeVoltage(fmt.Sprintf("b%d", p.Station))
	m := int(math.Round(T / p.Dt))
	r.Source = nlResample(wf.Time, vs, p.Warm, p.Dt, m)
	r.Station = nlResample(wf.Time, vr, p.Warm, p.Dt, m)
	r.T = make([]float64, m)
	for j := range r.T {
		r.T[j] = p.Warm + float64(j)*p.Dt
	}
	for _, v := range r.Source {
		r.SourcePeak = math.Max(r.SourcePeak, math.Abs(v))
	}
	r.Lines, r.SrcLines = map[int]complex128{}, map[int]complex128{}
	kmax := int(45 / p.F0)
	for k := 1; k <= kmax; k++ {
		r.Lines[k] = nlHann(r.Station, k)
		r.SrcLines[k] = nlHann(r.Source, k)
	}
	return r, nil
}

func nlResample(t, y []float64, t0, dt float64, m int) []float64 {
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

// nlHann is the Hann-windowed (unit gain) coefficient of e^{+j2πkt/T}.
func nlHann(x []float64, k int) complex128 {
	var acc complex128
	m := len(x)
	for j, v := range x {
		w := 1 - math.Cos(2*math.Pi*float64(j)/float64(m))
		acc += complex(v*w, 0) * cmplx.Exp(complex(0, -2*math.Pi*float64(k*j)/float64(m)))
	}
	return acc / complex(float64(m), 0)
}

// Products are the third-order (four-wave) lines and the drives.
func (p NLIonosphere) Products() map[string]int {
	return map[string]int{"f1": p.M1, "f2": p.M2, "2f1-f2": 2*p.M1 - p.M2, "2f2-f1": 2*p.M2 - p.M1,
		"2f1+f2": 2*p.M1 + p.M2, "3f1": 3 * p.M1}
}

// absBin returns |coefficient| at bin k, folding negative bins.
func absBin(lines map[int]complex128, k int) float64 {
	if k < 0 {
		k = -k
	}
	return 2 * cmplx.Abs(lines[k]) // amplitude of the real sine
}

// DriveFor finds the transmitter current at which the source band peaks at
// frac·Vh (linear cavity, both drives on).
func (p NLIonosphere) DriveFor(frac float64) (float64, error) {
	q := p
	q.Linear, q.Drive = true, 1
	r, err := q.Run()
	if err != nil {
		return 0, err
	}
	return frac * p.Vh() / r.SourcePeak, nil
}

func runNonlinearIonosphere(dataDir string) error {
	p := DefaultNLIonosphere()
	fmt.Println("Layer 7 — a nonlinear ionosphere: two transmitters in the cavity, loss G(E) = G0/(1 + (E/Ep)²) (electron heating)")
	c := p.cavity()
	fmt.Printf("  cavity of layer 3 fitted to mode 1 (7.8 Hz, Q 5); Ep = %.2f V/m (declared), h_e = %.1f km → Vh = %.0f V\n",
		p.Ep, c.He/1e3, p.Vh())
	fmt.Printf("  transmitters at the source band: %.1f and %.1f Hz; station %d° away\n", float64(p.M1)*p.F0, float64(p.M2)*p.F0, p.Station)
	type row struct {
		name string
		lin  NLRun
		nl   NLRun
	}
	var strong row
	fmt.Println("  source peak / Vh   ledger(lin, nl)   station lines [V]: f1, f2 | 2f1−f2 (1.3 Hz), 2f2−f1 (20.8 Hz), 2f1+f2, 3f1  (nonlinear; linear in [])")
	for _, frac := range []float64{0.01, 0.1, 0.5} {
		d, err := p.DriveFor(frac)
		if err != nil {
			return err
		}
		q := p
		q.Drive, q.Linear = d, true
		lin, err := q.Run()
		if err != nil {
			return err
		}
		q.Linear = false
		nl, err := q.Run()
		if err != nil {
			return err
		}
		pr := p.Products()
		fmt.Printf("  %5.2f (%.3g A)   %.0e %.0e   %.3g, %.3g | %.3g [%.0e], %.3g [%.0e], %.3g, %.3g\n",
			nl.SourcePeak/p.Vh(), d, lin.Ledger, nl.Ledger,
			absBin(nl.Lines, pr["f1"]), absBin(nl.Lines, pr["f2"]),
			absBin(nl.Lines, pr["2f1-f2"]), absBin(lin.Lines, pr["2f1-f2"]),
			absBin(nl.Lines, pr["2f2-f1"]), absBin(lin.Lines, pr["2f2-f1"]),
			absBin(nl.Lines, pr["2f1+f2"]), absBin(nl.Lines, pr["3f1"]))
		if frac == 0.5 {
			strong = row{"strong", lin, nl}
		}
	}
	pr := p.Products()
	r3 := absBin(strong.nl.Lines, pr["2f2-f1"]) / absBin(strong.nl.Lines, pr["f2"])
	fmt.Printf("  at 0.5 Vh the 20.8 Hz product is %.1f %% of the 14.3 Hz line at the station; the 7.8 Hz line changes by %+.1f %%\n",
		100*r3, 100*(absBin(strong.nl.Lines, pr["f1"])/absBin(strong.lin.Lines, pr["f1"])-1))
	fmt.Println("  (cross- and self-modulation through the shared loss; its sign follows from taking the loss ∝ σ(E),")
	fmt.Println("  an assumption of this shunt model). Linear cavity: the voltages simply add; no new line.")
	fmt.Println("  Scale: Vh corresponds to a vertical field of Ep ≈ 20 mV/m over the whole cavity height; natural")
	fmt.Println("  Schumann fields are far smaller, and real ELF transmitters (antennas kilometres long, radiating watts)")
	fmt.Println("  stay in the linear regime: at ELF the cavity is linear in practice; heating needs HF power (HAARP).")
	fmt.Println()
	if dataDir != "" {
		if err := os.MkdirAll(dataDir, 0o755); err != nil {
			return err
		}
		f, err := os.Create(filepath.Join(dataDir, "ionosphere_nl.csv"))
		if err != nil {
			return err
		}
		w := csv.NewWriter(f)
		_ = w.Write([]string{"t", "station_linear", "station_nonlinear", "source_linear", "source_nonlinear"})
		for j := range strong.nl.T {
			_ = w.Write([]string{ff(strong.nl.T[j]), ff(strong.lin.Station[j]), ff(strong.nl.Station[j]), ff(strong.lin.Source[j]), ff(strong.nl.Source[j])})
		}
		w.Flush()
		f.Close()
		g, err := os.Create(filepath.Join(dataDir, "ionosphere_nl_spectrum.csv"))
		if err != nil {
			return err
		}
		w = csv.NewWriter(g)
		_ = w.Write([]string{"f", "station_linear", "station_nonlinear"})
		for k := 1; k <= int(45/p.F0); k++ {
			_ = w.Write([]string{ff(float64(k) * p.F0), ff(absBin(strong.lin.Lines, k)), ff(absBin(strong.nl.Lines, k))})
		}
		w.Flush()
		g.Close()
	}
	return nil
}

func ff(x float64) string { return strconv.FormatFloat(x, 'g', 8, 64) }
