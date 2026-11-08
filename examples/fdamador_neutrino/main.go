// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command fdamador_neutrino tests hypothesis H28 (CLAUDE.md §1): the
// field of the toroidal-photon electron (examples/photon_electron, H23)
// turned inside out, so that the electric field is the ring and the
// magnetic field circles round it. The user's proposal (2026-11-04) is that
// this inverted structure is the neutrino.
//
// # What the inversion is
//
// In the charged-ring and photon pictures the ring carries a current and B
// wraps round the tube. Here E runs along the ring (E_φ, toroidal) and B
// wraps round it (poloidal). Three things follow from Maxwell's equations
// alone, before any model:
//
//  1. No charge. A closed loop of E has ∇·E = 0 everywhere, so by Gauss's
//     law the object encloses no charge. It is neutral, as a neutrino is.
//  2. It cannot be static. ∮E·dl ≠ 0 round the ring needs −dΦ_B/dt
//     through the hole (Faraday). The poloidal B threads the hole, so the
//     ring sustains itself only as an oscillation: E_φ is a ring of
//     displacement current ε·∂E_φ/∂t (Ampère–Maxwell), whose B threads the
//     hole, whose changing flux drives E_φ. That is a resonance at ω0.
//  3. Its magnetic moment oscillates. The displacement-current ring is a
//     magnetic dipole m(t) = I_d(t)·πa², with time average zero. A
//     time-varying magnetic dipole radiates (magnetic-dipole radiation), so
//     the object loses its energy at a rate set by its Q.
//
// Standard physics already contains this field pattern: it is the TE01δ
// mode of a dielectric ring resonator (E azimuthal in the ring, H poloidal
// through the hole), used in microwave filters, where the dielectric is
// what keeps the field in place (D. Kajfez, P. Guillon, Dielectric
// Resonators, Artech House 1986, ch. 1); bent into a travelling wave it is
// the TM-type guided wave with E along the path and B circling it.
//
// # Part 1: the E-ring resonator in kronspice
//
// The ring as a lumped series loop of N sections: each section a series
// capacitor (the displacement path along the ring, total
// C = ε0·εr·πr²/(2πa) = ε0·εr·r²/(2a)), the ring's self-inductance (thin
// ring with uniform current, L = μ0·a·(ln(8a/r) − 7/4), Maxwell, Treatise
// Art. 706, as in examples/kron_electron), and a series resistance for the
// radiated power of a small magnetic dipole, R_rad = (η0·π/6)·(ka)⁴
// (Balanis, Antenna Theory, 3rd ed., eq. 5-24), plus dielectric loss
// tan δ·ω0·L. Then Q = ω0·L/(R_rad + R_d), and with no dielectric loss
//
//	Q_rad = 6Λ/(π·(ka)³),  ka = √2·(a/r)/√(εr·Λ),  Λ = ln(8a/r) − 7/4,
//
// so Q_rad ∝ εr^(3/2): the permittivity is the confinement. The loop is
// solved in kronspice: AC (resonance and bandwidth) and a ring-down from a
// charged section (decay rate, ledger, net charge, moment).
//
// # Part 2: the neutrino test bench
//
// In vacuum (εr = 1) the ring is electrically large at its own resonance
// (ka ≈ 4 for a/r = 4), so the lumped model fails and the field leaves
// within half a period, as the field solve of examples/photon_electron
// Part 3 shows for the uninverted loop. The bench then asks what
// confinement a neutrino would need: SN1987A neutrinos (≈ 10 MeV) crossed
// 168,000 light-years; at a declared illustrative mass of 0.1 eV
// (below KATRIN's m < 0.45 eV) that is a proper time τ0 = t/γ, so the ring
// must hold its energy for ω0·τ0 radians. The required Q gives the
// required ka, and so the effective permittivity and the size.
//
// # Part 3: the inverted travelling wave as a circuit
//
// The ring line of examples/photon_electron with E and B exchanged is its
// circuit dual (Kron, Branin): series capacitors (E along the path) and
// shunt inductors (B round it), a high-pass backward-wave ladder. Once
// plain and once with its conductors crossed, as in photon_electron.
// Dispersion of the differential mode: ω = 1/(2√(C·Lg)·sin(πk/N)), k
// whole (plain) or k + ½ (crossed). The crossed ring keeps the 720° (spin ½)
// property in the inverted form.
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"sort"
	"strconv"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// Constants (CODATA 2018, SI 2019 exact where so).
const (
	C0    = 299792458.0
	Mu0   = 1.25663706212e-6
	Eps0  = 8.8541878128e-12
	Hbar  = 1.054571817e-34
	E     = 1.602176634e-19
	MuB   = 9.2740100783e-24
	Eta0  = Mu0 * C0 // vacuum impedance [Ω]
	Year  = 365.25 * 86400
	LYear = C0 * Year
)

// ERing is the E-ring: a torus of major radius A and tube radius R whose
// electric field runs round the ring inside a material of relative
// permittivity EpsR and loss tangent TanD, split into N lumped sections.
type ERing struct {
	A, R, EpsR, TanD float64
	N                int
}

// DefaultERing is a bench-scale ceramic ring: a = 10 mm, r = 2.5 mm,
// εr = 100, tan δ = 1e-4 (typical microwave ceramics), 24 sections.
func DefaultERing() ERing { return ERing{A: 10e-3, R: 2.5e-3, EpsR: 100, TanD: 1e-4, N: 24} }

// Lambda is ln(8a/r) − 7/4.
func (e ERing) Lambda() float64 { return math.Log(8*e.A/e.R) - 1.75 }

// L is the ring's self-inductance for a current spread uniformly over the
// tube, μ0·a·(ln(8a/r) − 7/4) [H].
func (e ERing) L() float64 { return Mu0 * e.A * e.Lambda() }

// C is the capacitance along the ring, ε0·εr·πr²/(2πa) [F]: the
// displacement path of length 2πa and cross-section πr².
func (e ERing) C() float64 { return Eps0 * e.EpsR * e.R * e.R / (2 * e.A) }

// Omega0 is the resonance 1/√(LC) [rad/s].
func (e ERing) Omega0() float64 { return 1 / math.Sqrt(e.L()*e.C()) }

// KA is the ring's electrical size k·a in the surrounding vacuum.
func (e ERing) KA() float64 { return e.Omega0() * e.A / C0 }

// RRad is the radiation resistance of a small magnetic dipole loop,
// (η0·π/6)·(ka)⁴ [Ω] (Balanis eq. 5-24, 20π²(C/λ)⁴ with η0 = 120π).
func (e ERing) RRad() float64 { ka := e.KA(); return Eta0 * math.Pi / 6 * ka * ka * ka * ka }

// RDiel is the dielectric loss as a series resistance at resonance,
// tan δ/(ω0·C) = tan δ·ω0·L [Ω].
func (e ERing) RDiel() float64 { return e.TanD * e.Omega0() * e.L() }

// Q is ω0·L/(R_rad + R_d).
func (e ERing) Q() float64 { return e.Omega0() * e.L() / (e.RRad() + e.RDiel()) }

// QRad is the radiation-only Q, 6Λ/(π·(ka)³).
func (e ERing) QRad() float64 { ka := e.KA(); return 6 * e.Lambda() / (math.Pi * ka * ka * ka) }

// MomentPerAmp is the magnetic moment per ampere of ring current, πa² [A·m²/A].
func (e ERing) MomentPerAmp() float64 { return math.Pi * e.A * e.A }

// Build makes the loop: section k is CS_k (N·C) from n_k to m_k, LS_k
// (L/N) from m_k to p_k, RS_k ((R_rad + R_d)/N) from p_k to n_{k+1}. The
// source element sits in section 0 between p_0 and q_0 (RS_0 then runs from
// q_0). Each section's capacitor isolates it at DC, and at the start of a
// ring-down (inductors open, uncharged capacitors open) each n_k too, so
// every m_k and n_k has a 1 GΩ bleeder to ground; at the ring's impedance
// (≈ 300 Ω) they are 10⁶ times too large to matter. If v0 ≠ 0, CS_0
// starts charged to v0.
func (e ERing) Build(src circuit.Element, v0 float64) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("E-ring resonator (fdamador_neutrino)")
	s := strconv.Itoa
	n := func(k int) string { return "n" + s(k%e.N) }
	cs, ls, rs := e.C()*float64(e.N), e.L()/float64(e.N), (e.RRad()+e.RDiel())/float64(e.N)
	for k := 0; k < e.N; k++ {
		if k == 0 && v0 != 0 {
			b.CapacitorIC("CS0", n(0), "m0", cs, v0)
		} else {
			b.Capacitor("CS"+s(k), n(k), "m"+s(k), cs)
		}
		b.Inductor("LS"+s(k), "m"+s(k), "p"+s(k), ls)
		from := "p" + s(k)
		if k == 0 {
			from = "q0"
		}
		b.Resistor("RS"+s(k), from, n(k+1), rs)
		b.Resistor("RB"+s(k), "m"+s(k), "0", 1e9)
		b.Resistor("RN"+s(k), n(k), "0", 1e9)
	}
	src.Name = "VS"
	b.Add(src, "p0", "q0")
	return b.BuildTransient()
}

// ACRun is the loop driven by 1 V in series, swept round its resonance.
type ACRun struct {
	F0, Q       float64 // peak of |I| and Q = f0/bandwidth(−3 dB)
	MaxResidual float64
	Freq, IMag  []float64
}

// Resonance sweeps ±(6/Q)·f0 round the analytic resonance.
func (e ERing) Resonance() (ACRun, error) {
	ckt, err := e.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, 0)
	if err != nil {
		return ACRun{}, err
	}
	f0 := e.Omega0() / (2 * math.Pi)
	span := 6 * f0 / e.Q()
	if span > 0.9*f0 {
		span = 0.9 * f0
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 6001, FStart: f0 - span, FStop: f0 + span}, analysis.Options{})
	if err != nil {
		return ACRun{}, err
	}
	i, _ := res.ElementCurrent("LS0")
	run := ACRun{MaxResidual: res.Ledger.MaxRelativeResidual, Freq: res.Freq}
	best := 0
	for k := range i {
		run.IMag = append(run.IMag, cmplx.Abs(i[k]))
		if run.IMag[k] > run.IMag[best] {
			best = k
		}
	}
	run.F0 = res.Freq[best]
	half := run.IMag[best] / math.Sqrt2
	cross := func(dir int) float64 {
		for k := best; k > 0 && k < len(i)-1; k += dir {
			if run.IMag[k+dir] < half {
				// linear interpolation between k and k+dir
				fr := (run.IMag[k] - half) / (run.IMag[k] - run.IMag[k+dir])
				return res.Freq[k] + fr*(res.Freq[k+dir]-res.Freq[k])
			}
		}
		return math.NaN()
	}
	run.Q = run.F0 / (cross(+1) - cross(-1))
	return run, nil
}

// RingDown is the loop kicked by a short voltage pulse and left to ring.
type RingDown struct {
	Q              float64 // from the decay of the current peaks after the kick
	EnergyIn       float64 // delivered by the kick [J]
	Dissipated     float64 // radiated + dielectric heat [J]
	StoredEnd      float64
	Mismatch       float64 // time-stepping error [J]
	MaxResidual    float64
	PeakMoment     float64 // largest |m| = |I|·πa² [A·m²]
	MeanMoment     float64 // time average of m over the last whole cycles [A·m²]
	PeakMomentLate float64 // largest |m| over the same cycles [A·m²]
	GroundCharge   float64 // |∫ current to ground| / (C_s·1 V): net charge lost
	Cycles         int
	Time, I        []float64
}

// RingDown kicks the loop with a 1 V pulse a quarter period long in the
// series source, then follows it for `cycles` periods with the source at
// 0 V (a short, so the loop is closed).
func (e ERing) RingDown(cycles int) (RingDown, error) {
	const v0 = 1.0
	T := 2 * math.Pi / e.Omega0()
	kick := circuit.Element{Kind: circuit.VoltageSource, Wave: circuit.Waveform{Kind: circuit.WavePulse,
		Args: []float64{0, v0, 0, T / 40, T / 40, T / 4, 0}}}
	ckt, err := e.Build(kick, 0)
	if err != nil {
		return RingDown{}, err
	}
	step := T / 200
	res, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: float64(cycles) * T, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return RingDown{}, err
	}
	i, _ := res.ElementCurrent("LS0")
	l := res.Ledger
	rd := RingDown{EnergyIn: l.EnergyIn, Dissipated: l.EnergyDissipated, StoredEnd: l.StoredEnd, Mismatch: l.IntegrationMismatch,
		MaxResidual: l.MaxRelativeResidual, Cycles: cycles, Time: res.Time, I: i}
	// Decay: fit ln|peak| against time over the local maxima of |i|.
	var ts, ys []float64
	for k := 1; k+1 < len(i); k++ {
		a := math.Abs(i[k])
		if res.Time[k] > T && a > math.Abs(i[k-1]) && a >= math.Abs(i[k+1]) && a > 0 {
			ts, ys = append(ts, res.Time[k]), append(ys, math.Log(a))
		}
	}
	slope := fitSlope(ts, ys)
	rd.Q = -e.Omega0() / (2 * slope)
	mpa := e.MomentPerAmp()
	t0 := float64(cycles-10) * T
	var sum, dt float64
	for k := 1; k < len(i); k++ {
		rd.PeakMoment = math.Max(rd.PeakMoment, math.Abs(i[k])*mpa)
		if res.Time[k-1] >= t0 {
			h := res.Time[k] - res.Time[k-1]
			sum += 0.5 * (i[k] + i[k-1]) * mpa * h
			dt += h
			rd.PeakMomentLate = math.Max(rd.PeakMomentLate, math.Abs(i[k])*mpa)
		}
	}
	rd.MeanMoment = sum / dt
	// Net charge: whatever left the ring went through the bleeders.
	var q float64
	for k := 0; k < e.N; k++ {
		for _, name := range []string{"RB", "RN"} {
			ib, _ := res.ElementCurrent(name + strconv.Itoa(k))
			for j := 1; j < len(ib); j++ {
				q += 0.5 * (ib[j] + ib[j-1]) * (res.Time[j] - res.Time[j-1])
			}
		}
	}
	rd.GroundCharge = math.Abs(q) / (e.C() * float64(e.N) * v0)
	return rd, nil
}

func fitSlope(x, y []float64) float64 {
	n := float64(len(x))
	var sx, sy, sxx, sxy float64
	for k := range x {
		sx += x[k]
		sy += y[k]
		sxx += x[k] * x[k]
		sxy += x[k] * y[k]
	}
	return (n*sxy - sx*sy) / (n*sxx - sx*sx)
}

// Sizes checks the ring against the quasi-static limit with em: each
// section's length 2πa/N and the ring's diameter, at f0, in vacuum outside.
func (e ERing) Sizes() (*em.SizeReport, error) {
	ckt, err := e.Build(circuit.Element{Kind: circuit.VoltageSource}, 0)
	if err != nil {
		return nil, err
	}
	lengths := map[string]float64{}
	for k := 0; k < e.N; k++ {
		lengths["LS"+strconv.Itoa(k)] = 2 * math.Pi * e.A / float64(e.N)
	}
	return em.CheckQuasiStatic(ckt, lengths, 2*e.A, e.Omega0()/(2*math.Pi), em.SizeOptions{WarnOnly: true})
}

// Bench is the neutrino test: a neutrino of rest mass MassEV and energy
// EnergyEV that travels for TravelS (lab time) must not lose its energy.
type Bench struct {
	MassEV, EnergyEV, TravelS, AoverR float64
}

// DefaultBench: 0.1 eV (declared, illustrative; KATRIN bounds m < 0.45 eV),
// 10 MeV (typical SN1987A neutrino), 168,000 years (SN1987A in the LMC,
// 51.4 kpc), ring proportions a/r = 4 as in Part 1.
func DefaultBench() Bench {
	return Bench{MassEV: 0.1, EnergyEV: 10e6, TravelS: 168000 * Year, AoverR: 4}
}

// BenchResult is what the E-ring needs to be that neutrino.
type BenchResult struct {
	Omega0, Gamma, ProperTime float64
	QNeeded                   float64
	VacuumKA, VacuumQ         float64 // εr = 1 (lumped model outside its validity)
	VacuumLifetime            float64 // Q/ω0 [s]
	KANeeded, EpsRNeeded      float64
	RadiusNeeded              float64 // a [m]
	MomentAmplitude           float64 // |m| when the stored energy is ħω0, in μB
}

// Run evaluates the bench.
func (b Bench) Run() BenchResult {
	lam := math.Log(8*b.AoverR) - 1.75
	var r BenchResult
	r.Omega0 = b.MassEV * E / Hbar
	r.Gamma = b.EnergyEV / b.MassEV
	r.ProperTime = b.TravelS / r.Gamma
	r.QNeeded = r.Omega0 * r.ProperTime
	qOf := func(ka float64) float64 { return 6 * lam / (math.Pi * ka * ka * ka) }
	r.VacuumKA = math.Sqrt2 * b.AoverR / math.Sqrt(lam)
	r.VacuumQ = qOf(r.VacuumKA)
	r.VacuumLifetime = r.VacuumQ / r.Omega0
	r.KANeeded = math.Cbrt(6 * lam / (math.Pi * r.QNeeded))
	r.EpsRNeeded = 2 * b.AoverR * b.AoverR / (lam * r.KANeeded * r.KANeeded)
	r.RadiusNeeded = r.KANeeded * C0 / r.Omega0
	// One quantum ħω0 in the ring: ½·L·I0² = ħω0, L = μ0·a·Λ.
	l := Mu0 * r.RadiusNeeded * lam
	i0 := math.Sqrt(2 * Hbar * r.Omega0 / l)
	r.MomentAmplitude = i0 * math.Pi * r.RadiusNeeded * r.RadiusNeeded / MuB
	return r
}

// DualRing is the ring line of examples/photon_electron with E and B
// exchanged: per section a series capacitor C in each wire and an
// inductor Lg (with series loss Rg) from each wire to ground.
type DualRing struct {
	N       int
	C, Lg   float64
	Rg      float64
	Crossed bool
}

// DefaultDualRing: 100 sections, C = 1 nF, Lg = 1 µH, Rg = 0.1 Ω.
func DefaultDualRing(crossed bool) DualRing {
	return DualRing{N: 100, C: 1e-9, Lg: 1e-6, Rg: 0.1, Crossed: crossed}
}

// ModeFrequency is differential mode m [Hz], counted from the top:
// ω = 1/(2√(C·Lg)·sin(πk/N)) with k = m (plain) or m − ½ (crossed), m ≥ 1.
// The highest mode has the longest wavelength: a backward-wave line.
func (d DualRing) ModeFrequency(m int) float64 {
	k := float64(m)
	if d.Crossed {
		k -= 0.5
	}
	return 1 / (2 * math.Sqrt(d.C*d.Lg) * math.Sin(math.Pi*k/float64(d.N))) / (2 * math.Pi)
}

// Build makes the dual ring with an AC current source at (a0, b0).
func (d DualRing) Build() (*circuit.Circuit, error) {
	b := circuit.NewBuilder(map[bool]string{false: "plain dual (CL) ring", true: "crossed dual (CL) ring"}[d.Crossed])
	s := strconv.Itoa
	a := func(k int) string { return "a" + s(k%d.N) }
	bb := func(k int) string { return "b" + s(k%d.N) }
	for k := 0; k < d.N; k++ {
		na, nb := a(k+1), bb(k+1)
		if k == d.N-1 && d.Crossed {
			na, nb = bb(0), a(0)
		}
		b.Capacitor("CA"+s(k), a(k), na, d.C)
		b.Capacitor("CB"+s(k), bb(k), nb, d.C)
		b.Inductor("LGA"+s(k), a(k), "ga"+s(k), d.Lg)
		b.Resistor("RGA"+s(k), "ga"+s(k), "0", d.Rg)
		b.Inductor("LGB"+s(k), bb(k), "gb"+s(k), d.Lg)
		b.Resistor("RGB"+s(k), "gb"+s(k), "0", d.Rg)
	}
	b.Add(circuit.Element{Name: "IS", Kind: circuit.CurrentSource, ACMag: 1}, "b0", "a0")
	return b.BuildTransient()
}

// Modes returns the top n peaks of the differential impedance at (a0, b0),
// highest first, and the sweep's worst ledger residual.
func (d DualRing) Modes(n int, fStart, fStop float64) ([]float64, float64, error) {
	ckt, err := d.Build()
	if err != nil {
		return nil, 0, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 19001, FStart: fStart, FStop: fStop}, analysis.Options{})
	if err != nil {
		return nil, 0, err
	}
	va, _ := res.NodeVoltage("a0")
	vb, _ := res.NodeVoltage("b0")
	z := func(i int) float64 { return cmplx.Abs(va[i] - vb[i]) }
	var peaks []float64
	for i := 1; i+1 < len(res.Freq); i++ {
		if z(i) > z(i-1) && z(i) >= z(i+1) {
			y0, y1, y2 := z(i-1), z(i), z(i+1)
			dd := 0.5 * (y0 - y2) / (y0 - 2*y1 + y2)
			peaks = append(peaks, res.Freq[i]+dd*(res.Freq[1]-res.Freq[0]))
		}
	}
	sort.Sort(sort.Reverse(sort.Float64Slice(peaks)))
	if len(peaks) > n {
		peaks = peaks[:n]
	}
	return peaks, res.Ledger.MaxRelativeResidual, nil
}

func main() {
	viz := flag.String("viz", "", "also write the data for viz/render.py (chart and sound) and viz/render_part56.py to this directory")
	field := flag.Bool("field", false, "also run Part 5b, the field-solver implosion (about 10 minutes for the default grids r/2, r/3)")
	fieldCells := flag.String("fieldcells", "2,3", "Part 5b grids, in cells per tube radius")
	only := flag.String("only", "", "run only this part (\"7\": Part 7, with 7b when -field is set; \"8\": Part 8, with 8d when -field is set; \"9\": Part 9, with the 9c field runs when -field is set)")
	spinCells := flag.String("spincells", "2", "Part 7b grids, in cells per tube radius")
	tripletCells := flag.String("tripletcells", "2,3", "Part 8d grids, in cells per tube radius")
	cpAll := flag.Bool("cpall", false, "Part 7a: search every defect section pair for the CP phase (about 5 minutes) instead of Part 6's pair")
	pairCells := flag.String("paircells", "2,3", "Part 9c grids, in cells per tube radius (the first grid runs every configuration, the others the single ring and d = a)")
	pairCSV := flag.String("pairfieldcsv", "", "Part 9c: with -field, read the field rows from this p9_field.csv (written by an earlier -field -viz run) instead of running them")
	flag.Parse()
	if *only == "9" {
		if err := runPart9(*field, parseCells(*pairCells), *viz, *pairCSV); err != nil {
			fmt.Println("   Part 9 error:", err)
		}
		return
	}
	if *only == "8" {
		if err := runPart8(*field, parseCells(*tripletCells), *viz); err != nil {
			fmt.Println("   Part 8 error:", err)
		}
		return
	}
	if *only == "7" {
		if err := runPart7(*field, *spinCells, *cpAll, *viz); err != nil {
			fmt.Println("   Part 7 error:", err)
		}
		return
	}
	if *viz != "" {
		if err := writeViz(*viz); err != nil {
			log.Fatal(err)
		}
		fmt.Println("viz data written to", *viz, "— run python3 viz/render.py for the chart and the sound")
	}
	fmt.Println("fdamador_neutrino (H28): the toroidal field inverted — E is the ring, B circles it")
	fmt.Println("Hypothesis: this structure is the neutrino. Standard prediction beside every result.")

	// Part 1.
	e := DefaultERing()
	fmt.Println("\n1. The E-ring resonator (a = 10 mm, r = 2.5 mm, εr = 100, tan δ = 1e-4)")
	fmt.Printf("   L = %.4g H, C = %.4g F, f0 = %.5g GHz, ka = %.3f\n", e.L(), e.C(), e.Omega0()/2/math.Pi/1e9, e.KA())
	fmt.Printf("   R_rad = %.4g Ω, R_d = %.4g Ω, Q = %.4g (radiation alone %.4g)\n", e.RRad(), e.RDiel(), e.Q(), e.QRad())
	ac, err := e.Resonance()
	if err != nil {
		fmt.Println("   AC error:", err)
		return
	}
	fmt.Printf("   kronspice AC: f0 = %.5g GHz, Q = %.4g, ledger residual %.1e\n", ac.F0/1e9, ac.Q, ac.MaxResidual)
	rd, err := e.RingDown(40)
	if err != nil {
		fmt.Println("   ring-down error:", err)
		return
	}
	fmt.Printf("   ring-down (40 periods): Q = %.4g; kick %.4g J = radiated + heat %.4g J + stored %.3g J, mismatch %.1e J, residual %.1e\n",
		rd.Q, rd.EnergyIn, rd.Dissipated, rd.StoredEnd, rd.Mismatch, rd.MaxResidual)
	fmt.Printf("   net charge: 0 by construction; lost to ground %.1e of one section's charge\n", rd.GroundCharge)
	fmt.Printf("   magnetic moment: peak %.3g A·m², over the last 10 periods mean %.2g of their peak (it oscillates)\n",
		rd.PeakMoment, rd.MeanMoment/rd.PeakMomentLate)
	if sz, err := e.Sizes(); err == nil {
		var dia, sec float64
		for _, en := range sz.Entries {
			if en.Name == "circuit" {
				dia = en.Ratio
			} else {
				sec = math.Max(sec, en.Ratio)
			}
		}
		fmt.Printf("   quasi-static check at f0 (em.CheckQuasiStatic): diameter %.3f λ, section %.4f λ\n", dia, sec)
	}
	fmt.Println("   Q against εr (formula; kronspice AC):")
	for _, er := range []float64{30, 100, 300, 1000} {
		x := e
		x.EpsR = er
		a, err := x.Resonance()
		if err != nil {
			fmt.Println("   ", err)
			continue
		}
		fmt.Printf("     εr %5g: f0 %.4g GHz, ka %.3f, Q_rad %.4g, Q %.4g; AC Q %.4g\n", er, x.Omega0()/2/math.Pi/1e9, x.KA(), x.QRad(), x.Q(), a.Q)
	}

	// Part 2.
	b := DefaultBench().Run()
	fmt.Println("\n2. Neutrino test bench (mass 0.1 eV declared, 10 MeV, SN1987A 168,000 years)")
	fmt.Printf("   ω0 = mc²/ħ = %.4g rad/s, γ = %.3g, proper time %.4g s, Q needed %.3g\n", b.Omega0, b.Gamma, b.ProperTime, b.QNeeded)
	fmt.Printf("   vacuum ring (εr = 1, a/r = 4): ka = %.2f (lumped model invalid), Q ≈ %.2g, lifetime ≈ %.2g s\n", b.VacuumKA, b.VacuumQ, b.VacuumLifetime)
	fmt.Printf("   to survive: ka = %.3g, effective εr = %.3g, ring radius %.3g m\n", b.KANeeded, b.EpsRNeeded, b.RadiusNeeded)
	fmt.Printf("   oscillating moment with one quantum ħω0 stored: %.3g μB (time average 0)\n", b.MomentAmplitude)
	fmt.Println("   Measured: charge consistent with 0; static moment < 2.8e-11 μB (Borexino 2017); m < 0.45 eV (KATRIN 2025).")
	fmt.Println("   Standard prediction: Dirac or Majorana neutrino, weak interaction only, moment ~3e-20 μB·(m/0.1 eV).")

	// Part 3.
	fmt.Println("\n3. The inverted ring line: series C (E along the path), shunt L (B round it)")
	for _, crossed := range []bool{false, true} {
		d := DefaultDualRing(crossed)
		peaks, resid, err := d.Modes(4, 10e6, 200e6)
		if err != nil {
			fmt.Println("  ", err)
			continue
		}
		name := map[bool]string{false: "plain  ", true: "crossed"}[crossed]
		fmt.Printf("   %s ring (ledger residual %.1e):\n", name, resid)
		for m := 1; m <= len(peaks); m++ {
			k := float64(m)
			if crossed {
				k -= 0.5
			}
			fmt.Printf("     mode %d: AC %.6g MHz, theory %.6g MHz, %g wavelengths round the ring\n", m, peaks[m-1]/1e6, d.ModeFrequency(m)/1e6, k)
		}
	}

	// Part 4.
	fmt.Println("\n4. Inflow: the zero-point field of the radiation resistance (Callen–Welton), booked on the magnetic channel")
	fmt.Println("   Scaled loop: L = 1 H, C = 1 F (ω0 = 1), R_rad = 1/Q, ħ → 1 J·s; energies in ħ, so ½ħω0 = 0.5")
	z := DefaultZPFRing()
	zr, err := z.Run()
	if err != nil {
		fmt.Println("   ", err)
		return
	}
	fmt.Printf("   Q = %.4g (the bench ring's radiation Q), %d realizations:\n", z.Q, z.Runs)
	fmt.Printf("     stored energy %.4f ± %.4f ħω0 (expected %.4f for the same noise modes)\n", zr.Energy, zr.Spread, zr.ExpEnergy)
	fmt.Printf("     ledger, magnetic channel: zero-point in %.2f J = radiated out %.2f J + stored %.2f J, time-stepping %.1e of the outflow\n",
		zr.ZPIn, zr.RadOut, zr.Stored, zr.Mismatch)
	fmt.Printf("     inflow / outflow = %.4f (the difference is the energy the ring filled up with)\n", zr.InflowOverOutflow)
	fmt.Printf("     energy fluctuation σ/⟨E⟩ = %.2f (≈ 1: random, not a fixed quantum); mean moment %.1e of its rms\n", zr.Fluctuation, zr.MeanMomentRatio)
	fmt.Println("   Boundaries change the radiation rate (Purcell factor F: R_rad → F·R_rad, the same mode-density change")
	fmt.Println("   that gives the Casimir force); the noise changes with it, so the held energy should not:")
	for _, f := range []float64{4, 0.2} {
		x := z
		x.Q = z.Q / f
		x.Runs = 2
		xr, err := x.Run()
		if err != nil {
			fmt.Println("   ", err)
			continue
		}
		fmt.Printf("     F = %-4g (Q = %6.4g): stored %.4f ± %.4f ħω0 (expected %.4f); in/out %.4f; lifetime ×%g\n",
			f, x.Q, xr.Energy, xr.Spread, xr.ExpEnergy, xr.InflowOverOutflow, 1/f)
	}
	fmt.Printf("   For the neutrino ring: the zero-point field holds it at ½ħω0 = %.3g eV (half the declared mass-energy) whatever F is,\n",
		0.5*DefaultBench().MassEV)
	fmt.Printf("   but that state is noise. To keep a coherent ring for SN1987A's proper time needs F ≈ %.1e (vacuum Q %.2g → %.2g);\n",
		b.VacuumQ/b.QNeeded, b.VacuumQ, b.QNeeded)
	fmt.Println("   the classic measurement of inhibited emission found ~20× (Hulet, Hilfer & Kleppner, PRL 55, 2137, 1985), and a")
	fmt.Println("   free neutrino has no walls. Standard prediction: the vacuum is the ground state; net exchange with it is zero.")

	// Parts 5 and 6.
	p5a, err := reportPart5a()
	if err != nil {
		fmt.Println("   Part 5a error:", err)
		return
	}
	var p5b []FieldResult
	if *field {
		var cells []float64
		for _, f := range strings.Split(*fieldCells, ",") {
			c, err := strconv.ParseFloat(strings.TrimSpace(f), 64)
			if err != nil || c <= 0 {
				fmt.Println("   bad -fieldcells:", f)
				return
			}
			cells = append(cells, c)
		}
		p5b, err = reportPart5b(cells)
		if err != nil {
			fmt.Println("   Part 5b error:", err)
			return
		}
	} else {
		fmt.Println("\n5b. Field solver: run with -field (about 10 minutes); see run_output.txt for the recorded run.")
	}
	p6, err := reportPart6()
	if err != nil {
		fmt.Println("   Part 6 error:", err)
		return
	}
	if *viz != "" {
		if err := writeViz56(*viz, p5a, p5b, p6); err != nil {
			log.Fatal(err)
		}
		fmt.Println("\nviz data for Parts 5 and 6 written to", *viz, "— run python3 viz/render_part56.py")
	}
	if err := runPart7(*field, *spinCells, *cpAll, *viz); err != nil {
		fmt.Println("   Part 7 error:", err)
	}
	if err := runPart8(*field, parseCells(*tripletCells), *viz); err != nil {
		fmt.Println("   Part 8 error:", err)
	}
	if err := runPart9(*field, parseCells(*pairCells), *viz, *pairCSV); err != nil {
		fmt.Println("   Part 9 error:", err)
	}
}

// parseCells reads a comma-separated list of grid resolutions (bad
// entries are skipped).
func parseCells(s string) []float64 {
	var out []float64
	for _, f := range strings.Split(s, ",") {
		if c, err := strconv.ParseFloat(strings.TrimSpace(f), 64); err == nil && c > 0 {
			out = append(out, c)
		}
	}
	return out
}

// runPart7 is Part 7 (H40): 7a always, 7b with -field, then 7c.
func runPart7(field bool, spinCells string, cpAll bool, viz string) error {
	p7a, err := reportPart7a(cpAll)
	if err != nil {
		return err
	}
	var p7b []Part7b
	if field {
		var cells []float64
		for _, f := range strings.Split(spinCells, ",") {
			c, err := strconv.ParseFloat(strings.TrimSpace(f), 64)
			if err != nil || c <= 0 {
				return fmt.Errorf("bad -spincells: %q", f)
			}
			cells = append(cells, c)
		}
		p7b, err = reportPart7b(cells)
		if err != nil {
			return err
		}
	} else {
		fmt.Println("7b. Field solver: run with -field; see run_output.txt for the recorded run.")
	}
	reportMass7(p7a, p7b)
	if viz != "" {
		if err := writeViz7(viz, p7a, p7b); err != nil {
			return err
		}
		fmt.Println("\nviz data for Part 7 written to", viz, "— run python3 viz/render_part7.py")
	}
	return nil
}
