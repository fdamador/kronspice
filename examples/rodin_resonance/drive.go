// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/output"
)

// Sense of the counter-wound pair: which end of B is driven.
const (
	Toroidal = +1 // through-the-hole parts of A and B add: toroidal, confined B
	Dipole   = -1 // round-the-axis parts add: a loop-like field through the hole (the paper's series connection)
)

// SourceOhms is the drive amplifier's output resistance (declared).
const SourceOhms = 1.0

// Drive is one lumped two-winding run: winding A from a sine source, B
// from a source of phase PhaseB (90° = cosine), each through SourceOhms,
// optionally through a series capacitor tuned to the winding's own L.
type Drive struct {
	Name   string
	F      float64
	VA, VB float64 // amplitudes [V] (0 = winding open)
	PhaseB float64 // degrees
	Sense  int
	Tuned  bool
	// TuneLA, TuneLB set the inductance each capacitor is tuned to (0 =
	// the winding's own L), e.g. L + M to tune a mode of the coupled pair.
	TuneLA, TuneLB float64
	Cycles         int
	StepsPerCyc    int
}

// DriveRun is the result.
type DriveRun struct {
	D                 Drive
	IA, IB            complex128 // steady-state path-sense current phasors (peak) [A]
	PA, PB            float64    // average power delivered by each source over the last cycle [W]
	Copper, SrcR      float64    // average heat in the winding wire and in the source resistances [W]
	Exchange          float64    // average power passed from A to B through M [W]
	ExchangePhasor    float64    // ω·M·I_A,rms·I_B,rms·sin φ [W]
	WmMin, WmMax      float64    // magnetic energy over the last cycle [J]
	WmMean            float64
	EIn, EDiss, EStor float64 // ledger over the whole run [J]
	Mismatch          float64
	Time, Wm          []float64 // magnetic energy vs time
	IAt, IBt          []float64
}

// Build makes the circuit.
func (m *Model) Build(d Drive) (*circuit.Circuit, error) {
	b := circuit.NewBuilder(d.Name)
	r := WireResistance(m.LenA, d.F)
	rb := WireResistance(m.LenB, d.F)
	if d.VA != 0 {
		b.Resistor("RSA", sine(b, "VA", "sa", d.VA, d.F, 0), "xa", SourceOhms)
		y := "xa"
		if d.Tuned {
			b.Capacitor("CA", "xa", "ya", 1/(math.Pow(2*math.Pi*d.F, 2)*orL(d.TuneLA, m.LA)))
			y = "ya"
		}
		b.Inductor("LA", y, "ma", m.LA)
		b.Resistor("RWA", "ma", "0", r)
	}
	if d.VB != 0 {
		b.Resistor("RSB", sine(b, "VB", "sb", d.VB, d.F, d.PhaseB), "xb", SourceOhms)
		y := "xb"
		if d.Tuned {
			b.Capacitor("CB", "xb", "yb", 1/(math.Pow(2*math.Pi*d.F, 2)*orL(d.TuneLB, m.LB)))
			y = "yb"
		}
		if d.Sense == Toroidal {
			b.Inductor("LB", y, "mb", m.LB)
			b.Resistor("RWB", "mb", "0", rb)
		} else { // B driven from its far end: path-sense current = −(circuit current)
			b.Resistor("RWB", y, "mb", rb)
			b.Inductor("LB", "0", "mb", m.LB)
		}
	}
	if d.VA != 0 && d.VB != 0 {
		b.Mutual("KAB", "LA", "LB", m.K())
	}
	return b.BuildTransient()
}

func orL(x, def float64) float64 {
	if x > 0 {
		return x
	}
	return def
}

// sine adds V·sin(2πft + phase°) from ground to node lo and returns the
// node to connect to: when the sinusoid is not zero at t = 0 (the cosine
// starts at its peak), a series EXP source cancels that value and fades
// over two cycles, so the run starts from rest; both elements count as
// the source (named name and name+"X").
func sine(b *circuit.Builder, name, lo string, v, f, phaseDeg float64) string {
	b.VoltageWave(name, lo, "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, v, f, 0, 0, phaseDeg}})
	v0 := v * math.Sin(phaseDeg*math.Pi/180)
	if math.Abs(v0) < 1e-12*math.Abs(v) {
		return lo
	}
	top := lo + "x"
	b.VoltageWave(name+"X", top, lo, circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{-v0, 0, 0, 2 / f, 1e3, 1}})
	return top
}

// TranOpts: absolute tolerances at the scale of these circuits.
var TranOpts = analysis.TranOptions{Options: analysis.Options{AbsTol: 1e-9, VnTol: 1e-6}}

// Run simulates a drive and extracts the steady state from the last cycle.
func (m *Model) Run(d Drive) (*DriveRun, error) {
	ckt, err := m.Build(d)
	if err != nil {
		return nil, err
	}
	T := 1 / d.F
	stop := float64(d.Cycles) * T
	step := T / float64(d.StepsPerCyc)
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, TranOpts)
	if err != nil {
		return nil, err
	}
	r := &DriveRun{D: d, Time: w.Time}
	get := func(n string) []float64 {
		v, ok := w.ElementCurrent(n)
		if !ok {
			return make([]float64, len(w.Time))
		}
		return v
	}
	ia, ib := get("LA"), get("LB")
	r.IAt, r.IBt = ia, ib
	r.Wm = make([]float64, len(w.Time))
	for k := range w.Time {
		r.Wm[k] = 0.5*m.LA*ia[k]*ia[k] + 0.5*m.LB*ib[k]*ib[k] + m.M*ia[k]*ib[k]
	}
	t0 := stop - T
	r.IA = phasor(w.Time, ia, d.F, t0)
	r.IB = phasor(w.Time, ib, d.F, t0)
	avg := func(name string) float64 { return cycleAvg(w, name, t0, stop) }
	r.PA, r.PB = -avg("VA")-avg("VAX"), -avg("VB")-avg("VBX")
	r.Copper = avg("RWA") + avg("RWB")
	r.SrcR = avg("RSA") + avg("RSB")
	// Exchange: EMF induced in B by A's current, times B's current.
	if d.VA != 0 && d.VB != 0 {
		n := 2000
		s := 0.0
		for k := 0; k < n; k++ {
			t := t0 + (float64(k)+0.5)*T/float64(n)
			h := T / float64(n) / 4
			dia := (output.Interp(w.Time, ia, t+h) - output.Interp(w.Time, ia, t-h)) / (2 * h)
			s += m.M * dia * output.Interp(w.Time, ib, t)
		}
		r.Exchange = s / float64(n)
		phi := cmplx.Phase(r.IB) - cmplx.Phase(r.IA) // P_A→B = ⟨i_B·M·di_A/dt⟩ = ω·M·I_A·I_B·sin(φ_B − φ_A) (rms)
		r.ExchangePhasor = 2 * math.Pi * d.F * m.M * cmplx.Abs(r.IA) / math.Sqrt2 * cmplx.Abs(r.IB) / math.Sqrt2 * math.Sin(phi)
	}
	r.WmMin, r.WmMax = math.Inf(1), math.Inf(-1)
	cnt := 0
	for k, t := range w.Time {
		if t < t0 {
			continue
		}
		r.WmMin = math.Min(r.WmMin, r.Wm[k])
		r.WmMax = math.Max(r.WmMax, r.Wm[k])
	}
	nn := 400
	for k := 0; k < nn; k++ {
		t := t0 + (float64(k)+0.5)*T/float64(nn)
		a, b := output.Interp(w.Time, ia, t), output.Interp(w.Time, ib, t)
		r.WmMean += 0.5*m.LA*a*a + 0.5*m.LB*b*b + m.M*a*b
		cnt++
	}
	r.WmMean /= float64(cnt)
	l := w.Ledger
	r.EIn, r.EDiss, r.EStor, r.Mismatch = l.EnergyIn, l.EnergyDissipated, l.StoredEnd-l.StoredStart, l.IntegrationMismatch
	return r, nil
}

// phasor is the complex amplitude (peak) of y at frequency f over one
// period from t0, y ≈ Re(Y·e^{jωt}) + c0 + c1·(t − t0): the offset and
// slope absorb what is left of the start-up transient (a slowly decaying
// L/R current), by least squares on 1000 points.
func phasor(t, y []float64, f, t0 float64) complex128 {
	n := 1000
	T := 1 / f
	var ata [4][4]float64
	var atb [4]float64
	for k := 0; k < n; k++ {
		tt := t0 + (float64(k)+0.5)*T/float64(n)
		s, c := math.Sincos(2 * math.Pi * f * tt)
		row := [4]float64{c, -s, 1, (tt - t0) / T}
		v := output.Interp(t, y, tt)
		for i := 0; i < 4; i++ {
			atb[i] += row[i] * v
			for j := 0; j < 4; j++ {
				ata[i][j] += row[i] * row[j]
			}
		}
	}
	x := solve4(ata, atb)
	return complex(x[0], x[1])
}

func solve4(a [4][4]float64, b [4]float64) [4]float64 {
	for c := 0; c < 4; c++ {
		p := c
		for r := c + 1; r < 4; r++ {
			if math.Abs(a[r][c]) > math.Abs(a[p][c]) {
				p = r
			}
		}
		a[c], a[p] = a[p], a[c]
		b[c], b[p] = b[p], b[c]
		for r := c + 1; r < 4; r++ {
			f := a[r][c] / a[c][c]
			for j := c; j < 4; j++ {
				a[r][j] -= f * a[c][j]
			}
			b[r] -= f * b[c]
		}
	}
	var x [4]float64
	for r := 3; r >= 0; r-- {
		s := b[r]
		for j := r + 1; j < 4; j++ {
			s -= a[r][j] * x[j]
		}
		x[r] = s / a[r][r]
	}
	return x
}

// cycleAvg is the mean absorbed power of an element over [t0, t1].
func cycleAvg(w *output.Waveforms, name string, t0, t1 float64) float64 {
	for e, n := range w.ElementNames {
		if n != name && !equalFold(n, name) {
			continue
		}
		// trapezoids clipped to [t0, t1], P linear within each step
		s := 0.0
		for k := 1; k < len(w.Time); k++ {
			ta, tb := w.Time[k-1], w.Time[k]
			if tb <= t0 || ta >= t1 {
				continue
			}
			pa, pb := w.P[e][k-1], w.P[e][k]
			ca, cb := math.Max(ta, t0), math.Min(tb, t1)
			p := func(t float64) float64 { return pa + (pb-pa)*(t-ta)/(tb-ta) }
			s += 0.5 * (p(ca) + p(cb)) * (cb - ca)
		}
		return s / (t1 - t0)
	}
	return 0
}

func equalFold(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		x, y := a[i], b[i]
		if 'A' <= x && x <= 'Z' {
			x += 32
		}
		if 'A' <= y && y <= 'Z' {
			y += 32
		}
		if x != y {
			return false
		}
	}
	return true
}

// FieldStats describes the field of a drive at a point or on a plane.
type FieldStats struct {
	Peak, RMS, Ellipticity float64 // [T], [T], minor/major
}

// pointStats returns the ellipse of B(t) = Re((I_A·b_A + I_B·b_B)e^{jωt}).
func pointStats(ia, ib complex128, ba, bb P3) FieldStats {
	var c [3]complex128
	for i := 0; i < 3; i++ {
		c[i] = ia*complex(ba[i], 0) + ib*complex(bb[i], 0)
	}
	mag2, cc := 0.0, complex(0, 0)
	for i := 0; i < 3; i++ {
		mag2 += real(c[i] * cmplx.Conj(c[i]))
		cc += c[i] * c[i]
	}
	a2 := 0.5 * (mag2 + cmplx.Abs(cc))
	b2 := math.Max(0, 0.5*(mag2-cmplx.Abs(cc)))
	e := 0.0
	if a2 > 0 {
		e = math.Sqrt(b2 / a2)
	}
	return FieldStats{Peak: math.Sqrt(a2), RMS: math.Sqrt(mag2 / 2), Ellipticity: e}
}

// Plane is a grid of points with the per-amp fields of A and B.
type Plane struct {
	Name   string
	Pts    []P3
	BA, BB []P3
	Near   []bool // within 3 mm of a bundle path (excluded from statistics)
}

// NewPlane samples the plane "xy" (z = 0) or "xz" (y = 0) on an n × n grid
// over ±half.
func (m *Model) NewPlane(name string, n int, half float64) *Plane {
	pl := &Plane{Name: name}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			u := -half + 2*half*float64(i)/float64(n-1)
			v := -half + 2*half*float64(j)/float64(n-1)
			p := P3{u, v, 0}
			if name == "xz" {
				p = P3{u, 0, v}
			}
			pl.Pts = append(pl.Pts, p)
		}
	}
	pl.BA = make([]P3, len(pl.Pts))
	pl.BB = make([]P3, len(pl.Pts))
	pl.Near = make([]bool, len(pl.Pts))
	parallelFor(len(pl.Pts), func(k int) {
		p := pl.Pts[k]
		pl.BA[k] = m.FieldA(p)
		pl.BB[k] = m.FieldB(p)
		pl.Near[k] = distToPath(m.P.A, p) < 3e-3 || distToPath(m.P.B, p) < 3e-3
	})
	return pl
}

// PlaneStats returns the largest peak |B| on the plane (away from the
// wires), its position, the plane-averaged ⟨B²⟩/2μ0 and the mean
// ellipticity weighted by ⟨B²⟩.
func (pl *Plane) Stats(ia, ib complex128) (peak float64, at P3, energy, ell float64) {
	wsum, n := 0.0, 0
	for k := range pl.Pts {
		if pl.Near[k] {
			continue
		}
		s := pointStats(ia, ib, pl.BA[k], pl.BB[k])
		if s.Peak > peak {
			peak, at = s.Peak, pl.Pts[k]
		}
		w := s.RMS * s.RMS
		energy += w
		ell += w * s.Ellipticity
		wsum += w
		n++
	}
	return peak, at, energy / float64(n) / (2 * em.Mu0), ell / wsum
}

// RingTravel tests whether the field pattern turns round the axis. On the
// circle ρ = R0, z = 0 it takes the 12-fold Fourier component of |B|²
// (the poles' symmetry allows only multiples of 12) at nt instants of one
// cycle and unwraps its phase: a pattern that travels round the axis
// advances it steadily (360° per 30° of travel); a standing pattern does
// not. It returns the net advance over the cycle [degrees of pattern
// phase], the angular travel it means (advance/12), and the ratio of the
// smallest to the largest 12-fold amplitude during the cycle.
func (m *Model) RingTravel(ia, ib complex128, nphi, nt int) (advance, travel, ampRatio float64) {
	ba := make([]P3, nphi)
	bb := make([]P3, nphi)
	parallelFor(nphi, func(k int) {
		ph := 2 * math.Pi * float64(k) / float64(nphi)
		p := P3{TorusR0 * math.Cos(ph), TorusR0 * math.Sin(ph), 0}
		ba[k], bb[k] = m.FieldA(p), m.FieldB(p)
	})
	prev, acc := 0.0, 0.0
	lo, hi := math.Inf(1), 0.0
	for it := 0; it <= nt; it++ {
		wt := 2 * math.Pi * float64(it) / float64(nt)
		ea := real(ia * cmplx.Exp(complex(0, wt)))
		eb := real(ib * cmplx.Exp(complex(0, wt)))
		var c complex128
		for k := 0; k < nphi; k++ {
			b := add(scl(ba[k], ea), scl(bb[k], eb))
			c += complex(dot(b, b), 0) * cmplx.Exp(complex(0, -12*2*math.Pi*float64(k)/float64(nphi)))
		}
		a := cmplx.Abs(c)
		lo, hi = math.Min(lo, a), math.Max(hi, a)
		ph := cmplx.Phase(c)
		if it > 0 {
			d := ph - prev
			for d > math.Pi {
				d -= 2 * math.Pi
			}
			for d < -math.Pi {
				d += 2 * math.Pi
			}
			acc += d
		}
		prev = ph
	}
	advance = acc * 180 / math.Pi
	return advance, advance / 12, lo / hi
}

func parallelFor(n int, f func(k int)) {
	ch := make(chan int, 256)
	done := make(chan struct{})
	workers := numWorkers()
	for w := 0; w < workers; w++ {
		go func() {
			for k := range ch {
				f(k)
			}
			done <- struct{}{}
		}()
	}
	for k := 0; k < n; k++ {
		ch <- k
	}
	close(ch)
	for w := 0; w < workers; w++ {
		<-done
	}
}

func (r *DriveRun) String() string {
	return fmt.Sprintf("|IA| %.4g A, |IB| %.4g A, PA %.4g W, PB %.4g W", cmplx.Abs(r.IA), cmplx.Abs(r.IB), r.PA, r.PB)
}
