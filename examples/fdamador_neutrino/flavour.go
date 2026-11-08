// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 6: the user's flavour hypothesis, H39 (CLAUDE.md §1, 2026-11-07):
// "for the flavors of the neutrino, it seems to me to be three different
// types of standing wave oscillations on the wave: constructive,
// destructive and standing wave. Check for that and see what the
// implosion/explosion would do to this coherent unit."
//
// The uniform E-ring of Part 1 (TE01δ, azimuthal order 0) has no
// travelling/standing distinction: it is one state. The three wave types
// live on the ring's azimuthal harmonics n ≥ 1, i.e. on the inverted ring
// lines of Part 3 (series C, shunt L). On a ring, harmonic n is a
// degenerate pair, z_j = A·e^{−iκj} + B·e^{+iκj} (κ = 2πn/N; crossed
// ring κ = π(2m−1)/N), with z_j = u_j + iωL_g·g_j built from the node
// voltage u_j and the shunt current g_j, so that A is the wave running one
// way (CW) and B the other (CCW). One fixed definition is used in every
// test, relative to a reference point c on the ring (the production /
// detection site):
//
//	constructive (C): travelling, B = 0 — every point reaches the full
//	                  amplitude, nothing cancels anywhere;
//	standing (S):     |A| = |B|, antinode of u at c;
//	destructive (D):  |A| = |B|, node of u at c — the two waves cancel at c.
//
// A state that cancels everywhere does not exist: Σ_j|z_j|² = N(|A|²+|B|²)
// exactly (the two waves are orthogonal), so only zero cancels
// everywhere, and "destructive" is a standing wave with its node at the
// chosen point. On a perfect ring S and D are the same state turned by a
// quarter wavelength; and C = (S + iD)/√2, so the three types span only a
// two-dimensional space.
//
// Test 1 (persistence), Test 2 (flavour oscillation: a defect lifts the
// degeneracy; the band-edge three-state cluster for 3×3 mixing) and Test 3
// (implosion and explosion of each type through a port) follow.

import (
	"math"
	"math/cmplx"
	"runtime"
	"sort"
	"strconv"
	"sync"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// FRing is an inverted (CL) two-wire ring line for Part 6: per section a
// series capacitor C·(1 + δ_j) in each wire, a shunt inductor Lg (with
// series loss Rg) from every node to ground, and optional ports.
type FRing struct {
	N       int
	C, Lg   float64
	Cg      float64 // stray capacitance from every node to ground
	Rg      float64
	Crossed bool
	Def     map[int]float64 // section → δ (both wires)
	Ports   []FPort
	// Gyro, when non-zero, puts a gyrator of transconductance Gyro [S]
	// across every section, from node pair j to node pair j+1: the ring's
	// rotation in its own frame (Part 7, spin.go). Zero in Parts 1–6.
	Gyro float64
	// NodeC, when set, adds capacitance from node pair j to ground (both
	// wires, beside Cg): a shunt defect (Part 8, triplet.go). Nil in
	// Parts 1–7.
	NodeC map[int]float64
}

// FPort is a differential port at node Node: a source of 2·v_inc behind R
// between a_Node and b_Node (an incident wave on a line of impedance R).
// Drive nil and PWL nil: a passive (matched) port.
type FPort struct {
	Node  int
	R     float64
	Drive *Drive    // SIN segments
	PWL   []float64 // or a PWL v_inc(t), as pairs (time, v_inc)
}

// Comp is one harmonic's content: z_j = A·e^{−iκj} + B·e^{iκj} at ω.
type Comp struct {
	Kappa, Omega float64
	A, B         complex128
}

// DefaultFRing: the element values of Part 3 (C = 1 nF, Lg = 1 µH) on 24
// sections (Part 3 uses 100; 24 keeps the beat runs short), with a shunt
// loss of 1 mΩ (Q ≈ 1.2×10⁵ at the first harmonic) and 1 pF of stray
// capacitance from each node to ground. The stray capacitance is what lets
// the state be imposed exactly at t = 0 (UIC): without a capacitive path
// to ground the nodes' common level is undetermined at the imposed point.
// It lowers ω by 1.5 % at the first harmonic, uniformly, and is in every
// formula below.
func DefaultFRing(crossed bool) FRing {
	return FRing{N: 24, C: 1e-9, Lg: 1e-6, Cg: 1e-12, Rg: 1e-3, Crossed: crossed}
}

// Kappa is the phase per section of harmonic m (m ≥ 1).
func (r FRing) Kappa(m int) float64 {
	if r.Crossed {
		return math.Pi * float64(2*m-1) / float64(r.N)
	}
	return 2 * math.Pi * float64(m) / float64(r.N)
}

// OmegaK is the unperturbed ring's ω at phase κ per section.
func (r FRing) OmegaK(k float64) float64 {
	return 1 / math.Sqrt(r.Lg*(4*r.C*sq(math.Sin(k/2))+r.Cg))
}

// Type is a wave type at harmonic κ relative to reference point c (in
// sections; node j sits at j), unit |A|² + |B|².
func (r FRing) Type(name string, kappa, c float64) Comp {
	cp := Comp{Kappa: kappa, Omega: r.OmegaK(kappa)}
	e := cmplx.Exp(complex(0, kappa*c))
	s2 := complex(math.Sqrt2/2, 0)
	switch name {
	case "C":
		cp.A = 1
	case "CCW":
		cp.B = 1
	case "S":
		cp.A, cp.B = s2*e, s2*cmplx.Conj(e)
	case "D":
		cp.A, cp.B = complex(0, 1)*s2*e, complex(0, -1)*s2*cmplx.Conj(e)
	}
	return cp
}

// Reversed is the time reverse of a state: z → z̄, i.e. A → B̄, B → Ā.
func Reversed(cs []Comp) []Comp {
	out := make([]Comp, len(cs))
	for i, c := range cs {
		out[i] = Comp{Kappa: c.Kappa, Omega: c.Omega, A: cmplx.Conj(c.B), B: cmplx.Conj(c.A)}
	}
	return out
}

// state returns u_j and g_j for a superposition.
func (r FRing) state(cs []Comp) (u, g []float64) {
	u, g = make([]float64, r.N), make([]float64, r.N)
	for _, c := range cs {
		for j := 0; j < r.N; j++ {
			z := c.A*cmplx.Exp(complex(0, -c.Kappa*float64(j))) + c.B*cmplx.Exp(complex(0, c.Kappa*float64(j)))
			u[j] += real(z)
			g[j] += imag(z) / (c.Omega * r.Lg)
		}
	}
	return
}

// Build makes the ring with initial state cs (UIC) and its ports.
func (r FRing) Build(cs []Comp) (*circuit.Circuit, error) {
	b := circuit.NewBuilder(map[bool]string{false: "plain CL ring (flavour test)", true: "crossed CL ring (flavour test)"}[r.Crossed])
	s := strconv.Itoa
	u, g := r.state(cs)
	va := func(j int) float64 { return u[j] / 2 }
	vb := func(j int) float64 { return -u[j] / 2 }
	for j := 0; j < r.N; j++ {
		c := r.C * (1 + r.Def[j])
		na, nb := "a"+s((j+1)%r.N), "b"+s((j+1)%r.N)
		vna, vnb := va((j+1)%r.N), vb((j+1)%r.N)
		if j == r.N-1 && r.Crossed {
			na, nb = "b0", "a0"
			vna, vnb = vb(0), va(0)
		}
		b.CapacitorIC("CA"+s(j), "a"+s(j), na, c, va(j)-vna)
		b.CapacitorIC("CB"+s(j), "b"+s(j), nb, c, vb(j)-vnb)
		b.InductorIC("LGA"+s(j), "a"+s(j), "ga"+s(j), r.Lg, g[j]/2)
		b.Resistor("RGA"+s(j), "ga"+s(j), "0", r.Rg)
		b.InductorIC("LGB"+s(j), "b"+s(j), "gb"+s(j), r.Lg, -g[j]/2)
		b.Resistor("RGB"+s(j), "gb"+s(j), "0", r.Rg)
		b.CapacitorIC("CGA"+s(j), "a"+s(j), "0", r.Cg+r.NodeC[j], va(j))
		b.CapacitorIC("CGB"+s(j), "b"+s(j), "0", r.Cg+r.NodeC[j], vb(j))
		if r.Gyro != 0 {
			b.Gyrator("YR"+s(j), "a"+s(j), "b"+s(j), na, nb, r.Gyro)
		}
	}
	for _, p := range r.Ports {
		x := "x" + s(p.Node)
		b.ResistorTo("RP"+s(p.Node), x, "a"+s(p.Node), p.R, exchange.Magnetic)
		switch {
		case p.Drive != nil:
			addSources(b, p.Drive.Segs, x, "b"+s(p.Node))
		case p.PWL != nil:
			w := make([]float64, len(p.PWL))
			for k := range p.PWL {
				w[k] = p.PWL[k]
				if k%2 == 1 {
					w[k] *= 2
				}
			}
			b.Add(circuit.Element{Name: "VP" + x, Kind: circuit.VoltageSource, Channel: exchange.Magnetic,
				ChannelNote: "incident wave on the ring's port (declared)", Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: w}}, x, "b"+s(p.Node))
		default:
			addSources(b, nil, x, "b"+s(p.Node))
		}
	}
	return b.BuildTransient()
}

// FRun is a transient of the ring with its readout.
type FRun struct {
	Time        []float64
	U, G        [][]float64 // [t][j]
	Stored      []float64
	StoredStart float64
	StoredEnd   float64
	Dissipated  float64
	SourceWork  float64
	Mismatch    float64
	LedgerRel   float64 // |work + stored_start − dissipated − stored_end| / max(work, stored_start)
	PortV       [][]float64
	PortI       [][]float64
	PortVs      [][]float64
}

// Run simulates for `stop` seconds at `per` steps per period of ω.
func (r FRing) Run(cs []Comp, omega, stop float64, per int) (FRun, error) {
	ckt, err := r.Build(cs)
	if err != nil {
		return FRun{}, err
	}
	h := 2 * math.Pi / omega / float64(per)
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: h, Stop: stop, MaxStep: h, UIC: true}, analysis.TranOptions{MaxSteps: 5_000_000})
	if err != nil {
		return FRun{}, err
	}
	return r.read(wf)
}

func (r FRing) read(wf *output.Waveforms) (FRun, error) {
	s := strconv.Itoa
	fr := FRun{Time: wf.Time, Stored: wf.Ledger.StoredEnergy, StoredStart: wf.Ledger.StoredStart, StoredEnd: wf.Ledger.StoredEnd,
		Mismatch: wf.Ledger.IntegrationMismatch}
	nt := len(wf.Time)
	fr.U, fr.G = make([][]float64, nt), make([][]float64, nt)
	for t := range fr.U {
		fr.U[t], fr.G[t] = make([]float64, r.N), make([]float64, r.N)
	}
	for j := 0; j < r.N; j++ {
		a, _ := wf.NodeVoltage("a" + s(j))
		bb, _ := wf.NodeVoltage("b" + s(j))
		ia, _ := wf.ElementCurrent("LGA" + s(j))
		ib, _ := wf.ElementCurrent("LGB" + s(j))
		for t := 0; t < nt; t++ {
			fr.U[t][j] = a[t] - bb[t]
			fr.G[t][j] = ia[t] - ib[t]
		}
	}
	for _, n := range wf.ElementNames {
		if n[0] == 'R' {
			fr.Dissipated += elementEnergy(wf.Time, wf.ElementNames, wf.P, n)
		}
	}
	fr.SourceWork = sourceWork(wf.Time, wf.ElementNames, wf.P)
	for _, p := range r.Ports {
		x, _ := wf.NodeVoltage("x" + s(p.Node))
		a, _ := wf.NodeVoltage("a" + s(p.Node))
		bb, _ := wf.NodeVoltage("b" + s(p.Node))
		v, i, vs := make([]float64, nt), make([]float64, nt), make([]float64, nt)
		for t := 0; t < nt; t++ {
			v[t] = a[t] - bb[t]
			i[t] = (x[t] - a[t]) / p.R
			vs[t] = x[t] - bb[t]
		}
		fr.PortV, fr.PortI, fr.PortVs = append(fr.PortV, v), append(fr.PortI, i), append(fr.PortVs, vs)
	}
	den := math.Max(fr.SourceWork, fr.StoredStart)
	fr.LedgerRel = math.Abs(fr.SourceWork+fr.StoredStart-fr.Dissipated-fr.StoredEnd) / den
	return fr, nil
}

// Amp projects the state at sample t onto harmonic κ (at ω): A, B.
func (r FRing) Amp(fr FRun, t int, kappa, omega float64) (A, B complex128) {
	for j := 0; j < r.N; j++ {
		z := complex(fr.U[t][j], omega*r.Lg*fr.G[t][j])
		A += z * cmplx.Exp(complex(0, kappa*float64(j)))
		B += z * cmplx.Exp(complex(0, -kappa*float64(j)))
	}
	n := complex(float64(r.N), 0)
	return A / n, B / n
}

// Fractions of a (A, B) state in the four reference states at c.
type Fractions struct{ CW, CCW, S, D float64 }

// Frac returns the energy fractions of (A, B) in CW, CCW, S(c), D(c).
func Frac(A, B complex128, kappa, c float64) Fractions {
	e := cmplx.Exp(complex(0, -kappa*c))
	al := A*e + B*cmplx.Conj(e)
	be := complex(0, -1) * (A*e - B*cmplx.Conj(e))
	s := sq(cmplx.Abs(A)) + sq(cmplx.Abs(B))
	return Fractions{CW: sq(cmplx.Abs(A)) / s, CCW: sq(cmplx.Abs(B)) / s, S: sq(cmplx.Abs(al)) / (2 * s), D: sq(cmplx.Abs(be)) / (2 * s)}
}

func sq(x float64) float64 { return x * x }

// purity is the share of Σ|z_j|² carried by harmonic κ at sample t.
func (r FRing) purity(fr FRun, t int, kappa, omega float64) float64 {
	A, B := r.Amp(fr, t, kappa, omega)
	var tot float64
	for j := 0; j < r.N; j++ {
		tot += sq(fr.U[t][j]) + sq(omega*r.Lg*fr.G[t][j])
	}
	return float64(r.N) * (sq(cmplx.Abs(A)) + sq(cmplx.Abs(B))) / tot
}

// ---------------------------------------------------------------- Test 1

// Persist is one type left alone on the perfect ring.
type Persist struct {
	Ring, Type             string
	Periods                float64
	Start, MinOwn, MaxOwn  float64 // own-type fraction at the start, min and max over the run
	MaxOther               float64 // largest fraction in the type it should not have (CCW for C, D for S, S for D)
	MinPurity              float64 // smallest share in the harmonic
	Rotation               float64 // for C: phase advance of A per period / (2π) — the pattern turns
	EnergyEnd              float64 // ring energy at the end / start
	LedgerRel, MismatchRel float64
	Snap                   [][]float64 // u_j at t = 0, T/8, T/4, 3T/8, T/2
}

// Persistence runs C, S, D (reference c = ½, the midpoint of section 0)
// on the plain and crossed rings, harmonic 1, for `periods` periods.
func Persistence(periods float64) ([]Persist, error) {
	var out []Persist
	for _, crossed := range []bool{false, true} {
		r := DefaultFRing(crossed)
		k := r.Kappa(1)
		w := r.OmegaK(k)
		for _, ty := range []string{"C", "S", "D"} {
			fr, err := r.Run([]Comp{r.Type(ty, k, 0.5)}, w, periods*2*math.Pi/w, 60)
			if err != nil {
				return out, err
			}
			p := Persist{Ring: map[bool]string{false: "plain", true: "crossed"}[crossed], Type: ty, Periods: periods, MinOwn: 1, MinPurity: 1}
			var a0 complex128
			for t := range fr.Time {
				A, B := r.Amp(fr, t, k, w)
				f := Frac(A, B, k, 0.5)
				own, other := f.CW, f.CCW
				switch ty {
				case "S":
					own, other = f.S, f.D
				case "D":
					own, other = f.D, f.S
				}
				if t == 0 {
					p.Start = own
					a0 = A
				}
				p.MinOwn, p.MaxOwn = math.Min(p.MinOwn, own), math.Max(p.MaxOwn, own)
				p.MaxOther = math.Max(p.MaxOther, other)
				p.MinPurity = math.Min(p.MinPurity, r.purity(fr, t, k, w))
				if q := len(p.Snap); q <= 4 && fr.Time[t] >= float64(q)*math.Pi/(4*w)-1e-15 {
					p.Snap = append(p.Snap, append([]float64(nil), fr.U[t]...))
				}
				_ = a0
			}
			if ty == "C" {
				// count turns of A: unwrap its phase.
				var ph, prev float64
				for t := range fr.Time {
					A, _ := r.Amp(fr, t, k, w)
					x := cmplx.Phase(A)
					if t > 0 {
						d := x - prev
						for d > math.Pi {
							d -= 2 * math.Pi
						}
						for d < -math.Pi {
							d += 2 * math.Pi
						}
						ph += d
					}
					prev = x
				}
				p.Rotation = ph / (2 * math.Pi) / periods
			}
			p.EnergyEnd = fr.StoredEnd / fr.StoredStart
			p.LedgerRel = fr.LedgerRel
			p.MismatchRel = math.Abs(fr.Mismatch) / fr.StoredStart
			out = append(out, p)
		}
	}
	return out, nil
}

// ---------------------------------------------------------------- eigen

// eigSym is the cyclic Jacobi method for a symmetric matrix: eigenvalues
// and orthonormal eigenvectors (columns of v).
func eigSym(a [][]float64) (vals []float64, v [][]float64) {
	n := len(a)
	m := make([][]float64, n)
	v = make([][]float64, n)
	for i := range m {
		m[i] = append([]float64(nil), a[i]...)
		v[i] = make([]float64, n)
		v[i][i] = 1
	}
	for sweep := 0; sweep < 100; sweep++ {
		var off float64
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				off += m[i][j] * m[i][j]
			}
		}
		if off < 1e-30 {
			break
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if math.Abs(m[p][q]) < 1e-300 {
					continue
				}
				th := (m[q][q] - m[p][p]) / (2 * m[p][q])
				t := 1 / (math.Abs(th) + math.Sqrt(th*th+1))
				if th < 0 {
					t = -t
				}
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					mkp, mkq := m[k][p], m[k][q]
					m[k][p], m[k][q] = c*mkp-s*mkq, s*mkp+c*mkq
				}
				for k := 0; k < n; k++ {
					mpk, mqk := m[p][k], m[q][k]
					m[p][k], m[q][k] = c*mpk-s*mqk, s*mpk+c*mqk
				}
				for k := 0; k < n; k++ {
					vkp, vkq := v[k][p], v[k][q]
					v[k][p], v[k][q] = c*vkp-s*vkq, s*vkp+c*vkq
				}
			}
		}
	}
	vals = make([]float64, n)
	for i := range vals {
		vals[i] = m[i][i]
	}
	return vals, v
}

// Modes returns the ring's exact differential modes: ω (ascending) and
// node-voltage eigenvectors (columns), from ω²·Cmat·u = u/Lg with Cmat the
// capacitance-weighted ring Laplacian (antiperiodic for the crossed ring)
// plus the stray capacitance on the diagonal.
func (r FRing) Modes() (omega []float64, vec [][]float64) {
	n := r.N
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	for j := 0; j < n; j++ {
		c := r.C * (1 + r.Def[j])
		k := (j + 1) % n
		sign := -1.0
		if r.Crossed && j == n-1 {
			sign = 1
		}
		m[j][j] += c + r.Cg + r.NodeC[j]
		m[k][k] += c
		m[j][k] += sign * c
		m[k][j] += sign * c
	}
	vals, v := eigSym(m)
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	sort.Slice(idx, func(a, b int) bool { return vals[idx[a]] > vals[idx[b]] }) // largest μ = lowest ω
	vec = make([][]float64, n)
	for i := range vec {
		vec[i] = make([]float64, n)
	}
	for col, i := range idx {
		omega = append(omega, 1/math.Sqrt(r.Lg*vals[i]))
		for j := 0; j < n; j++ {
			vec[j][col] = v[j][i]
		}
	}
	return
}

// pairNear returns the two exact ω nearest w.
func (r FRing) pairNear(w float64) (lo, hi float64) {
	om, _ := r.Modes()
	sort.Slice(om, func(a, b int) bool { return math.Abs(om[a]-w) < math.Abs(om[b]-w) })
	return math.Min(om[0], om[1]), math.Max(om[0], om[1])
}

// ---------------------------------------------------------------- Test 2

// Beat is a two-state flavour oscillation run.
type Beat struct {
	Ring           string
	Delta          float64 // defect δ on section 0
	Start          string  // "C" (CW) or the standing state at offset s
	Offset         float64 // production offset s from the defect centre, sections
	ThetaPred      float64 // predicted mixing angle, degrees (45 for C; κ·s for standing)
	DfExact        float64 // exact |f_a − f_b| from the eigenproblem [Hz]
	DfFirst        float64 // first order: δ·f/(N(1+δ)) · 4C·sin²(κ/2)/(4C·sin²(κ/2) + Cg) [Hz]
	DfFit          float64 // fitted from the kronspice run [Hz]
	AmpFit         float64 // fitted sin²2θ
	ThetaFit       float64 // degrees
	F0             float64 // unperturbed f [Hz]
	MaxResid       float64 // largest |P − fit|
	LedgerRel      float64
	MismatchRel    float64
	T, P           []float64 // conversion probability against time (s), thinned
	Periods        float64
	FitPeriods     float64
	EnergyEnd      float64
	ProductionNote string
}

// FlavourBeat runs the two-state oscillation for defect δ: start in the CW
// travelling wave (start = "C") or in the standing wave with its
// antinode s sections from the defect's centre (start = "S"); P is the
// probability of the conjugate state (CCW, or the standing wave turned by
// a quarter wavelength).
func FlavourBeat(crossed bool, delta float64, start string, s float64, beats float64) (Beat, error) {
	r := DefaultFRing(crossed)
	r.Def = map[int]float64{0: delta}
	k := r.Kappa(1)
	w := r.OmegaK(k)
	lo, hi := r.pairNear(w)
	bt := Beat{Ring: map[bool]string{false: "plain", true: "crossed"}[crossed], Delta: delta, Start: start, Offset: s, F0: w / (2 * math.Pi)}
	bt.DfExact = (hi - lo) / (2 * math.Pi)
	sk := 4 * r.C * sq(math.Sin(k/2))
	bt.DfFirst = delta / (1 + delta) / float64(r.N) * bt.F0 * sk / (sk + r.Cg)
	c := 0.5 + s
	var init Comp
	if start == "C" {
		init = r.Type("C", k, c)
		bt.ThetaPred = 45
	} else {
		init = r.Type("S", k, c)
		bt.ThetaPred = k * s * 180 / math.Pi
	}
	stop := beats / bt.DfExact
	fr, err := r.Run([]Comp{init}, w, stop, 60)
	if err != nil {
		return bt, err
	}
	bt.Periods = stop * bt.F0
	var ts, ps []float64
	for t := range fr.Time {
		A, B := r.Amp(fr, t, k, w)
		f := Frac(A, B, k, c)
		p := f.CCW
		if start != "C" {
			p = f.D
		}
		ts, ps = append(ts, fr.Time[t]), append(ps, p)
	}
	bt.DfFit, bt.AmpFit, bt.MaxResid = fitBeat(ts, ps, bt.DfExact)
	bt.ThetaFit = 0.5 * math.Asin(math.Sqrt(math.Min(1, bt.AmpFit))) * 180 / math.Pi
	for t := 0; t < len(ts); t += 20 {
		bt.T, bt.P = append(bt.T, ts[t]), append(bt.P, ps[t])
	}
	bt.LedgerRel = fr.LedgerRel
	bt.MismatchRel = math.Abs(fr.Mismatch) / fr.StoredStart
	bt.EnergyEnd = fr.StoredEnd / fr.StoredStart
	return bt, nil
}

// fitBeat fits P = a·sin²(πΔf·t): scan Δf round guess, a by linear least
// squares, then refine by golden section.
func fitBeat(t, p []float64, guess float64) (df, a, maxRes float64) {
	cost := func(f float64) (float64, float64) {
		var sgg, spg float64
		for k := range t {
			g := sq(math.Sin(math.Pi * f * t[k]))
			sgg += g * g
			spg += p[k] * g
		}
		a := spg / sgg
		var c float64
		for k := range t {
			d := p[k] - a*sq(math.Sin(math.Pi*f*t[k]))
			c += d * d
		}
		return c, a
	}
	best, bf := math.Inf(1), guess
	for x := 0.5; x <= 1.5; x += 0.002 {
		if c, _ := cost(guess * x); c < best {
			best, bf = c, guess*x
		}
	}
	lo, hi := bf-0.002*guess, bf+0.002*guess
	g := (math.Sqrt(5) - 1) / 2
	for it := 0; it < 80; it++ {
		m1, m2 := hi-g*(hi-lo), lo+g*(hi-lo)
		c1, _ := cost(m1)
		c2, _ := cost(m2)
		if c1 < c2 {
			hi = m2
		} else {
			lo = m1
		}
	}
	df = (lo + hi) / 2
	_, a = cost(df)
	for k := range t {
		maxRes = math.Max(maxRes, math.Abs(p[k]-a*sq(math.Sin(math.Pi*df*t[k]))))
	}
	return df, a, maxRes
}

// ---------------------------------------------------------------- three states

// Cluster is the band-edge three-state system of the plain ring with even
// N: the π-mode (κ = π, one state, a standing wave that cannot travel:
// neighbouring nodes in antiphase) and the pair at κ = π − 2π/N, nearly
// degenerate with it (the CL ladder's dispersion is flat at the band
// edge). Flavour states (fixed, relative to the production point c = ½):
// f1 = standing cos at κ11 (antinode at c), f2 = standing sin (node at c),
// f3 = the π-mode. Two capacitor defects mix them.
type Cluster struct {
	D1, D2     float64 // δ on section P1 and on section P2
	P1, P2     int
	Omega      [3]float64 // the three exact mode ω, ascending
	U          [3][3]float64
	Leak       float64 // largest share of a mass state outside the three flavour states
	Assign     [3]int  // flavour index used as e, μ, τ
	Mass       [3]int  // mode index used as 1, 2, 3
	Th12, Th23 float64 // degrees (PDG convention from |U|)
	Th13       float64
	Ratio      float64 // |ω3 − ω1| / |ω2 − ω1|
	Cost       float64
}

func (r FRing) clusterFlavours() [3][]float64 {
	n := r.N
	k := r.Kappa(n/2 - 1)
	var f [3][]float64
	for a := range f {
		f[a] = make([]float64, n)
	}
	for j := 0; j < n; j++ {
		x := float64(j) - 0.5
		f[0][j] = math.Cos(k * x)
		f[1][j] = math.Sin(k * x)
		f[2][j] = math.Cos(math.Pi * float64(j))
	}
	for a := range f {
		var s float64
		for _, x := range f[a] {
			s += x * x
		}
		for j := range f[a] {
			f[a][j] /= math.Sqrt(s)
		}
	}
	return f
}

// ClusterMix computes the exact modes for the defects and the 3×3 mixing.
func ClusterMix(d1, d2 float64, p1, p2 int) Cluster {
	r := DefaultFRing(false)
	r.Def = map[int]float64{p1: d1}
	r.Def[p2] += d2
	om, vec := r.Modes()
	f := r.clusterFlavours()
	n := r.N
	type cand struct {
		i int
		w float64
	}
	var cs []cand
	for i := 0; i < n; i++ {
		var s float64
		for a := 0; a < 3; a++ {
			var d float64
			for j := 0; j < n; j++ {
				d += f[a][j] * vec[j][i]
			}
			s += d * d
		}
		cs = append(cs, cand{i, s})
	}
	sort.Slice(cs, func(a, b int) bool { return cs[a].w > cs[b].w })
	ids := []int{cs[0].i, cs[1].i, cs[2].i}
	sort.Slice(ids, func(a, b int) bool { return om[ids[a]] < om[ids[b]] })
	cl := Cluster{D1: d1, D2: d2, P1: p1, P2: p2}
	for m, i := range ids {
		cl.Omega[m] = om[i]
		var s float64
		for a := 0; a < 3; a++ {
			var d float64
			for j := 0; j < n; j++ {
				d += f[a][j] * vec[j][i]
			}
			cl.U[a][m] = d
			s += d * d
		}
		cl.Leak = math.Max(cl.Leak, 1-s)
	}
	return cl
}

// angles fills the PDG angles for flavour assignment perm (e, μ, τ) and
// returns the cost against the measured values.
func (cl *Cluster) angles(perm [3]int, target PMNS) float64 {
	// mass labels: the closest pair is (1, 2), the other 3.
	w := cl.Omega
	pairs := [][2]int{{0, 1}, {1, 2}, {0, 2}}
	best := 0
	for q := range pairs {
		if math.Abs(w[pairs[q][0]]-w[pairs[q][1]]) < math.Abs(w[pairs[best][0]]-w[pairs[best][1]]) {
			best = q
		}
	}
	i1, i2 := pairs[best][0], pairs[best][1]
	i3 := 3 - i1 - i2
	col := func(i int) [3]float64 {
		v := [3]float64{cl.U[perm[0]][i], cl.U[perm[1]][i], cl.U[perm[2]][i]}
		n := math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
		return [3]float64{v[0] / n, v[1] / n, v[2] / n}
	}
	if math.Abs(col(i1)[0]) < math.Abs(col(i2)[0]) { // ν1 has the larger e content
		i1, i2 = i2, i1
	}
	c1, c2, c3 := col(i1), col(i2), col(i3)
	cl.Assign, cl.Mass = perm, [3]int{i1, i2, i3}
	cl.Th13 = math.Asin(math.Min(1, math.Abs(c3[0]))) * 180 / math.Pi
	cl.Th12 = math.Atan2(math.Abs(c2[0]), math.Abs(c1[0])) * 180 / math.Pi
	cl.Th23 = math.Atan2(math.Abs(c3[1]), math.Abs(c3[2])) * 180 / math.Pi
	cl.Ratio = math.Abs(w[i3]-w[i1]) / math.Abs(w[i2]-w[i1])
	cl.Cost = sq((cl.Th12-target.Th12)/target.S12) + sq((cl.Th23-target.Th23)/target.S23) + sq((cl.Th13-target.Th13)/target.S13) +
		sq(math.Log(cl.Ratio/target.Ratio)/target.SRatio)
	return cl.Cost
}

// PMNS holds the measured mixing (NuFIT 6.0, normal ordering, with
// Super-Kamiokande atmospheric data; Esteban et al., JHEP 12 (2024) 216,
// arXiv:2410.05380, Table 1): angles and 1σ widths in degrees, and the
// splitting ratio Δm²31/Δm²21 = 2.513e-3/7.49e-5 with its relative width.
type PMNS struct {
	Th12, S12, Th23, S23, Th13, S13 float64
	Ratio, SRatio                   float64
	DeltaCP                         float64 // degrees (best fit; consistent with CP conservation at 1σ in NO)
}

// NuFIT6 is the measured target.
func NuFIT6() PMNS {
	return PMNS{Th12: 33.68, S12: 0.72, Th23: 43.3, S23: 0.9, Th13: 8.56, S13: 0.11,
		Ratio: 2.513e-3 / 7.49e-5, SRatio: math.Hypot(0.02/2.513, 0.19/7.49), DeltaCP: 212}
}

// FitPMNS searches two capacitor defects (sizes |δ| ≤ 0.5, sections P1 <
// P2 anywhere on the ring) and the assignment of the three ring states to
// e, μ, τ for the closest match to the measured angles and splitting
// ratio (cost: squared deviations in units of the measured 1σ). Each
// section pair is searched by Nelder–Mead from nine starts, in parallel;
// results are kept by pair index, so the outcome does not depend on
// scheduling.
func FitPMNS() Cluster {
	target := NuFIT6()
	perms := [][3]int{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}
	eval := func(d1, d2 float64, p1, p2 int) Cluster {
		if math.Abs(d1) > 0.5 || math.Abs(d2) > 0.5 {
			return Cluster{Cost: math.Inf(1)}
		}
		cl := ClusterMix(d1, d2, p1, p2)
		bc := Cluster{Cost: math.Inf(1)}
		for _, pm := range perms {
			c := cl
			c.angles(pm, target)
			if c.Cost < bc.Cost {
				bc = c
			}
		}
		return bc
	}
	n := DefaultFRing(false).N
	var pairs [][2]int
	for p1 := 0; p1 < n/2; p1++ { // the reflection j → −j maps the rest onto these
		for p2 := p1 + 1; p2 < n; p2++ {
			pairs = append(pairs, [2]int{p1, p2})
		}
	}
	results := make([]Cluster, len(pairs))
	sem := make(chan struct{}, runtime.NumCPU())
	var wg sync.WaitGroup
	for q, pp := range pairs {
		wg.Add(1)
		sem <- struct{}{}
		go func(q, p1, p2 int) {
			defer wg.Done()
			defer func() { <-sem }()
			best := Cluster{Cost: math.Inf(1)}
			for _, s1 := range []float64{-0.2, 0.05, 0.2} {
				for _, s2 := range []float64{-0.2, 0.05, 0.2} {
					x := [3][2]float64{{s1, s2}, {s1 * 1.3, s2}, {s1, s2 * 1.3}}
					var fx [3]float64
					for i := range x {
						fx[i] = eval(x[i][0], x[i][1], p1, p2).Cost
					}
					for it := 0; it < 150; it++ {
						for i := 0; i < 3; i++ {
							for j := i + 1; j < 3; j++ {
								if fx[j] < fx[i] {
									fx[i], fx[j] = fx[j], fx[i]
									x[i], x[j] = x[j], x[i]
								}
							}
						}
						c := [2]float64{(x[0][0] + x[1][0]) / 2, (x[0][1] + x[1][1]) / 2}
						xr := [2]float64{2*c[0] - x[2][0], 2*c[1] - x[2][1]}
						fr := eval(xr[0], xr[1], p1, p2).Cost
						switch {
						case fr < fx[0]:
							xe := [2]float64{3*c[0] - 2*x[2][0], 3*c[1] - 2*x[2][1]}
							if fe := eval(xe[0], xe[1], p1, p2).Cost; fe < fr {
								x[2], fx[2] = xe, fe
							} else {
								x[2], fx[2] = xr, fr
							}
						case fr < fx[1]:
							x[2], fx[2] = xr, fr
						default:
							xc := [2]float64{(c[0] + x[2][0]) / 2, (c[1] + x[2][1]) / 2}
							if fc := eval(xc[0], xc[1], p1, p2).Cost; fc < fx[2] {
								x[2], fx[2] = xc, fc
							} else {
								for i := 1; i < 3; i++ {
									x[i] = [2]float64{(x[0][0] + x[i][0]) / 2, (x[0][1] + x[i][1]) / 2}
									fx[i] = eval(x[i][0], x[i][1], p1, p2).Cost
								}
							}
						}
					}
					if c := eval(x[0][0], x[0][1], p1, p2); c.Cost < best.Cost {
						best = c
					}
				}
			}
			results[q] = best
		}(q, pp[0], pp[1])
	}
	wg.Wait()
	best := Cluster{Cost: math.Inf(1)}
	for _, c := range results {
		if c.Cost < best.Cost {
			best = c
		}
	}
	return best
}

// ThreeRun is a kronspice run of the cluster started in one flavour.
type ThreeRun struct {
	From      int          // flavour index started in (0, 1, 2 = ring states f1, f2, f3)
	T         []float64    // seconds, thinned
	P         [3][]float64 // measured energy in f1, f2, f3 / the ring's energy (all harmonics)
	Pred      [3][]float64 // 3×3 oscillation formula with the fitted U and ω
	PredFull  [3][]float64 // all N modes of the ring, read out the same way
	TAll      []float64    // every time step (for the T-symmetry comparison)
	PAll      [3][]float64
	MaxDev    float64 // largest |measured − 3×3 formula|
	MaxDevAll float64 // largest |measured − all-mode prediction|
	LedgerRel float64
	Periods   float64
	Reverse   float64 // max |P(a→b)(t) − P(b→a)(t)| between this run and the reverse run (T symmetry: no CP phase)
}

// ThreeFlavour runs the cluster cl in kronspice from flavour `from`, for
// `slow` periods of the slow beat, and compares it with the 3×3 formula
// (U column-normalized) and with the ring's full mode expansion read out
// by the same projection.
func ThreeFlavour(cl Cluster, from int, slow float64) (ThreeRun, error) {
	r := DefaultFRing(false)
	r.Def = map[int]float64{cl.P1: cl.D1}
	r.Def[cl.P2] += cl.D2
	n := r.N
	k11 := r.Kappa(n/2 - 1)
	w11, w12 := r.OmegaK(k11), r.OmegaK(math.Pi)
	init := map[int]Comp{0: r.Type("S", k11, 0.5), 1: r.Type("D", k11, 0.5), 2: {Kappa: math.Pi, Omega: w12, A: 1}}[from]
	// slow beat: smallest |Δω|
	dw := math.Inf(1)
	for i := 0; i < 3; i++ {
		for j := i + 1; j < 3; j++ {
			dw = math.Min(dw, math.Abs(cl.Omega[i]-cl.Omega[j]))
		}
	}
	stop := slow * 2 * math.Pi / dw
	fr, err := r.Run([]Comp{init}, w11, stop, 120)
	if err != nil {
		return ThreeRun{}, err
	}
	tr := ThreeRun{From: from, LedgerRel: fr.LedgerRel, Periods: stop * w11 / (2 * math.Pi)}
	// Energy per |amplitude|² of harmonic κ is ∝ 1/ω² ∝ sin²(κ/2) + Cg/4C.
	wk := func(k float64) float64 { return sq(math.Sin(k/2)) + r.Cg/(4*r.C) }
	energies := func(fr FRun, t int) [3]float64 {
		A, B := r.Amp(fr, t, k11, w11)
		P, _ := r.Amp(fr, t, math.Pi, w12)
		fs := Frac(A, B, k11, 0.5)
		s11 := sq(cmplx.Abs(A)) + sq(cmplx.Abs(B))
		return [3]float64{wk(k11) * fs.S * s11, wk(k11) * fs.D * s11, wk(math.Pi) * sq(cmplx.Abs(P))}
	}
	// The whole ring's energy in the same units (every harmonic), so the
	// fractions do not drift with the 1 mΩ loss.
	total := func(fr FRun, t int) float64 {
		var s float64
		for m := 1; m <= n/2; m++ {
			k := r.Kappa(m)
			A, B := r.Amp(fr, t, k, r.OmegaK(k))
			if m == n/2 {
				s += wk(k) * sq(cmplx.Abs(A))
			} else {
				s += wk(k) * (sq(cmplx.Abs(A)) + sq(cmplx.Abs(B)))
			}
		}
		return s
	}
	// Column-normalized U for the 3×3 formula.
	var U [3][3]float64
	for i := 0; i < 3; i++ {
		var s float64
		for a := 0; a < 3; a++ {
			s += sq(cl.U[a][i])
		}
		for a := 0; a < 3; a++ {
			U[a][i] = cl.U[a][i] / math.Sqrt(s)
		}
	}
	// Full expansion: u(0) = the initial pattern, g(0) = 0, so each mode
	// goes as cos(ω_i t) in u and sin(ω_i t)/(ω_i·Lg) in g. The
	// trapezoidal rule at a fixed step h runs mode ω at (2/h)·atan(ωh/2);
	// over a thousand periods that shifts the beat phases by ~0.1 rad, so
	// the expansion uses the same warped frequencies (amplitudes unchanged).
	om, vec := r.Modes()
	hstep := 2 * math.Pi / w11 / 120
	for i := range om {
		om[i] = 2 / hstep * math.Atan(om[i]*hstep/2)
	}
	u0, _ := r.state([]Comp{init})
	c := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			c[i] += vec[j][i] * u0[j]
		}
	}
	one := FRun{U: [][]float64{make([]float64, n)}, G: [][]float64{make([]float64, n)}}
	for t := range fr.Time {
		tt := fr.Time[t]
		e := energies(fr, t)
		for j := 0; j < n; j++ {
			var u, g float64
			for i := 0; i < n; i++ {
				u += c[i] * vec[j][i] * math.Cos(om[i]*tt)
				g += c[i] * vec[j][i] * math.Sin(om[i]*tt) / (om[i] * r.Lg)
			}
			one.U[0][j], one.G[0][j] = u, g
		}
		ef := energies(one, 0)
		tm, tf := total(fr, t), total(one, 0)
		var pred [3]float64
		for b := 0; b < 3; b++ {
			var amp complex128
			for i := 0; i < 3; i++ {
				amp += complex(U[from][i]*U[b][i], 0) * cmplx.Exp(complex(0, cl.Omega[i]*tt))
			}
			pred[b] = sq(cmplx.Abs(amp))
			tr.MaxDev = math.Max(tr.MaxDev, math.Abs(e[b]/tm-pred[b]))
			tr.MaxDevAll = math.Max(tr.MaxDevAll, math.Abs(e[b]/tm-ef[b]/tf))
		}
		tr.TAll = append(tr.TAll, tt)
		for b := 0; b < 3; b++ {
			tr.PAll[b] = append(tr.PAll[b], e[b]/tm)
		}
		if t%20 == 0 {
			tr.T = append(tr.T, tt)
			for b := 0; b < 3; b++ {
				tr.P[b] = append(tr.P[b], e[b]/tm)
				tr.Pred[b] = append(tr.Pred[b], pred[b])
				tr.PredFull[b] = append(tr.PredFull[b], ef[b]/tf)
			}
		}
	}
	return tr, nil
}

// ---------------------------------------------------------------- Test 3

// PortR is the port resistance that gives a standing wave with its
// antinode at the port the external Q `q`: Q = ω·W/P with
// W = V²N(C·sin²(κ/2) + Cg/4)/2 and P = V²/(2R) for a differential
// antinode amplitude V.
func (r FRing) PortR(kappa, q float64) float64 {
	return q / (r.OmegaK(kappa) * float64(r.N) * (r.C*sq(math.Sin(kappa/2)) + r.Cg/4))
}

// Overlap is |⟨x|y⟩|²/(⟨x|x⟩⟨y|y⟩) for states in (A, B) of one harmonic.
func Overlap(xa, xb, ya, yb complex128) float64 {
	ip := cmplx.Conj(xa)*ya + cmplx.Conj(xb)*yb
	return sq(cmplx.Abs(ip)) / ((sq(cmplx.Abs(xa)) + sq(cmplx.Abs(xb))) * (sq(cmplx.Abs(ya)) + sq(cmplx.Abs(yb))))
}

// Explosion is a type left to radiate through passive ports.
type Explosion struct {
	Case          string
	Remaining     float64   // ring energy at the end / start
	End           Fractions // at the end, relative to c
	T, E, PConj   []float64 // time (periods), energy / start, conjugate-type fraction (thinned)
	LedgerRel     float64
	Periods       float64
	DarkPredicted string
}

// Explode runs the free explosion of cs through passive ports at `nodes`.
func (r FRing) Explode(name string, cs []Comp, nodes []int, R, periods, c float64, conj func(Fractions) float64) (Explosion, FRun, error) {
	k := cs[0].Kappa
	w := cs[0].Omega
	x := r
	x.Ports = nil
	for _, n := range nodes {
		x.Ports = append(x.Ports, FPort{Node: n, R: R})
	}
	stop := periods * 2 * math.Pi / w
	fr, err := x.Run(cs, w, stop, 60)
	if err != nil {
		return Explosion{}, fr, err
	}
	ex := Explosion{Case: name, Remaining: fr.StoredEnd / fr.StoredStart, LedgerRel: fr.LedgerRel, Periods: periods}
	A, B := x.Amp(fr, len(fr.Time)-1, k, w)
	ex.End = Frac(A, B, k, c)
	for t := 0; t < len(fr.Time); t += 10 {
		A, B := x.Amp(fr, t, k, w)
		ex.T = append(ex.T, fr.Time[t]*w/(2*math.Pi))
		ex.E = append(ex.E, fr.Stored[t]/fr.StoredStart)
		ex.PConj = append(ex.PConj, conj(Frac(A, B, k, c)))
	}
	return ex, fr, nil
}

// Rebuild is an implosion driven through the ports, then the explosion.
type Rebuild struct {
	Case         string
	EInc         float64 // incident energy [J]
	ERefDrive    float64 // outgoing while driven [J]
	Absorbed     float64 // 1 − ERefDrive/EInc
	StoredPeak   float64 // ring energy when the drive stops / EInc
	Fidelity     float64 // overlap of the ring state at the peak with the target
	PeakFrac     Fractions
	EOutAfter    float64 // outgoing after the drive / EInc
	RecRemaining float64 // energy the recording left in the ring / its start (what the reversal cannot return)
	LedgerRel    float64
	Mismatch     float64 // time-stepping mismatch / EInc
	T, E, PConj  []float64
	Periods      float64
}

// RebuildByReversal records the explosion of the time-reversed target
// through passive ports, plays each port's outgoing wave back reversed in
// time as its incident wave (PWL, 20 points a period), and follows the
// ring through the implosion and `after` periods of explosion, at 400
// steps a period: every PWL corner is a breakpoint where the transient
// takes a backward-Euler step, and at 60 steps a period those steps lost
// 3 % of the energy to time-stepping error (3e-4 at 400).
func (r FRing) RebuildByReversal(name string, target []Comp, nodes []int, R, recPeriods, after, c float64, conj func(Fractions) float64) (Rebuild, error) {
	k := target[0].Kappa
	w := target[0].Omega
	T := 2 * math.Pi / w
	trec := recPeriods * T
	_, rec, err := r.Explode(name, Reversed(target), nodes, R, recPeriods, c, conj)
	if err != nil {
		return Rebuild{}, err
	}
	x := r
	x.Ports = nil
	for p, n := range nodes {
		var pwl []float64
		nPts := int(recPeriods * 20)
		for q := 0; q <= nPts; q++ {
			tp := float64(q) * trec / float64(nPts)
			pwl = append(pwl, tp, interp(rec.Time, rec.PortV[p], trec-tp))
		}
		pwl = append(pwl, trec*(1+1e-9), 0)
		x.Ports = append(x.Ports, FPort{Node: n, R: R, PWL: pwl})
	}
	fr, err := x.Run(nil, w, trec+after*T, 400)
	if err != nil {
		return Rebuild{}, err
	}
	rb := Rebuild{Case: name, RecRemaining: rec.StoredEnd / rec.StoredStart, LedgerRel: fr.LedgerRel, Periods: recPeriods + after}
	x.fillRebuild(&rb, fr, target, trec, k, w, c, conj)
	return rb, nil
}

// RebuildByGrowth drives one port with the single-tone CPA growth at ω
// and rate γ (a SIN drive), for comparison.
func (r FRing) RebuildByGrowth(name string, target []Comp, node int, R, gamma, recPeriods, after, c float64, conj func(Fractions) float64) (Rebuild, error) {
	k := target[0].Kappa
	w := target[0].Omega
	T := 2 * math.Pi / w
	f := math.Sqrt(w*w-gamma*gamma) / (2 * math.Pi)
	K := math.Round(recPeriods * T * f)
	tp := K / f
	d := Drive{Start: 0, Stop: tp, Segs: piece(0, tp, math.Exp(-gamma*tp), f, -gamma)}
	x := r
	x.Ports = []FPort{{Node: node, R: R, Drive: &d}}
	fr, err := x.Run(nil, w, tp+after*T, 400)
	if err != nil {
		return Rebuild{}, err
	}
	rb := Rebuild{Case: name, LedgerRel: fr.LedgerRel, Periods: recPeriods + after}
	x.fillRebuild(&rb, fr, target, tp, k, w, c, conj)
	return rb, nil
}

func (r FRing) fillRebuild(rb *Rebuild, fr FRun, target []Comp, tstop, k, w, c float64, conj func(Fractions) float64) {
	var pin, pref float64
	ip := -1
	for t := range fr.Time {
		var in, rf float64
		for p := range r.Ports {
			vin := fr.PortVs[p][t] / 2
			vr := vin - r.Ports[p].R*fr.PortI[p][t]
			in += vin * vin / r.Ports[p].R
			rf += vr * vr / r.Ports[p].R
		}
		if t > 0 {
			h := fr.Time[t] - fr.Time[t-1]
			rb.EInc += 0.5 * (in + pin) * h
			if fr.Time[t] <= tstop*(1+1e-12) {
				rb.ERefDrive += 0.5 * (rf + pref) * h
			} else {
				rb.EOutAfter += 0.5 * (rf + pref) * h
			}
		}
		pin, pref = in, rf
		if fr.Time[t] <= tstop*(1+1e-12) {
			ip = t
		}
	}
	A, B := r.Amp(fr, ip, k, w)
	rb.Fidelity = Overlap(target[0].A, target[0].B, A, B)
	rb.PeakFrac = Frac(A, B, k, c)
	rb.StoredPeak = fr.Stored[ip] / rb.EInc
	rb.EOutAfter /= rb.EInc
	rb.Absorbed = 1 - rb.ERefDrive/rb.EInc
	rb.Mismatch = math.Abs(fr.Mismatch) / rb.EInc
	for t := 0; t < len(fr.Time); t += 10 {
		A, B := r.Amp(fr, t, k, w)
		rb.T = append(rb.T, fr.Time[t]*w/(2*math.Pi))
		rb.E = append(rb.E, fr.Stored[t]/rb.EInc)
		rb.PConj = append(rb.PConj, conj(Frac(A, B, k, c)))
	}
}

// interp is linear interpolation of y(x) at x0 (x ascending).
func interp(x, y []float64, x0 float64) float64 {
	if x0 <= x[0] {
		return y[0]
	}
	if x0 >= x[len(x)-1] {
		return y[len(y)-1]
	}
	k := sort.SearchFloat64s(x, x0)
	f := (x0 - x[k-1]) / (x[k] - x[k-1])
	return y[k-1] + f*(y[k]-y[k-1])
}
