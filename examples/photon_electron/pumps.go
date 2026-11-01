// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 4: three electron pumps that work, each with its energy ledger.
//
// The electron's own frequency in the toroidal-photon model (the Compton
// frequency m_e c²/h = 1.24×10²⁰ Hz) is 10⁸ times beyond any circuit, so a
// pump cannot be tuned to it. These pumps move electrons by resonance or
// by a periodic drive at frequencies circuits do reach:
//
//  1. Single-electron pump (Pothier et al., Europhys. Lett. 17, 249,
//     1992): two metal islands between three tunnel junctions, two gates
//     driven 90° apart; Coulomb blockade lets exactly one electron through
//     per cycle, I = e·f. Solved with the orthodox theory of single-
//     electron tunnelling (Averin & Likharev 1986; Grabert & Devoret,
//     Single Charge Tunneling, Plenum 1992, Ch. 2): a master equation over
//     the islands' electron numbers with rates Γ = −ΔF/(e²R(1 − e^{ΔF/kT})).
//  2. Multipactor pump (Farnsworth, US 2,091,439; §3i): electrons bounce
//     across a vacuum gap in step with the RF field of a resonator and
//     multiply at each impact (measured copper yield, em.CopperAsReceived);
//     a grid passes part of them to a collector. The gap is the resonator's
//     capacitor, the electron sheet loads it (Kishek & Lau, Phys. Rev. Lett.
//     80, 193, 1998), and space charge caps the sheet.
//  3. Ground charge pump (§3j): a resonant LC step-up drives a negative
//     Cockcroft–Walton ladder (Proc. R. Soc. A 136, 619, 1932) that pulls
//     electrons from ground into a storage bank, solved by kronspice.
//
// In every case the ground or the metal supplies the electrons and a
// declared source supplies the energy; the ledgers show it.

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

const (
	kB = 1.380649e-23 // J/K (exact)
	mE = Me
)

// ---------------------------------------------------------------------------
// 1. Single-electron pump

// SEPump is a two-island, three-junction pump with two gates.
type SEPump struct {
	C  [3]float64 // junction capacitances [F]: lead L–island 1, 1–2, 2–lead R
	R  [3]float64 // tunnel resistances [Ω]
	Cg [2]float64 // gate capacitances [F]
	T  float64    // temperature [K]
	V  float64    // bias: left lead at +V/2, right lead at −V/2 [V]
	// Gate charges ν_i = Cg_i·Vg_i/e: centre, amplitude, phase lag of gate 2.
	Centre [2]float64
	Amp    float64
	Phase  float64
}

// DefaultSEPump: aluminium-oxide junctions as in Pothier 1992 (C ≈ 0.2 fF,
// R ≈ 200 kΩ, well above h/e² = 25.8 kΩ), 20 mK, centred on the triple
// point found by TriplePoint.
func DefaultSEPump() SEPump {
	p := SEPump{C: [3]float64{200e-18, 200e-18, 200e-18}, R: [3]float64{200e3, 200e3, 200e3},
		Cg: [2]float64{50e-18, 50e-18}, T: 0.02, Amp: 0.25, Phase: math.Pi / 2}
	p.Centre = p.TriplePoint()
	return p
}

const nMax = 2 // electron numbers per island: −nMax … nMax

func (p SEPump) cmat() (a, b, c float64) { // island capacitance matrix [[a, −b], [−b, c]]
	return p.C[0] + p.C[1] + p.Cg[0], p.C[1], p.C[1] + p.C[2] + p.Cg[1]
}

// potentials of the islands for electron numbers n and gate charges ν.
func (p SEPump) phi(n [2]int, nu [2]float64) (float64, float64) {
	a, b, c := p.cmat()
	det := a*c - b*b
	q1 := p.C[0]*p.V/2 + nu[0]*E - float64(n[0])*E
	q2 := -p.C[2]*p.V/2 + nu[1]*E - float64(n[1])*E
	return (c*q1 + b*q2) / det, (b*q1 + a*q2) / det
}

// energy: electrostatic energy of the islands' charge configuration
// (for finding the triple point only).
func (p SEPump) energy(n [2]int, nu [2]float64) float64 {
	a, b, c := p.cmat()
	det := a*c - b*b
	q1 := (nu[0] - float64(n[0])) * E
	q2 := (nu[1] - float64(n[1])) * E
	return 0.5 * (c*q1*q1 + 2*b*q1*q2 + a*q2*q2) / det
}

// TriplePoint is the gate-charge point where the states (0,0), (1,0) and
// (0,1) have equal energy (at zero bias); the gates circle it.
func (p SEPump) TriplePoint() [2]float64 {
	best, at := math.Inf(1), [2]float64{}
	for i := 0; i <= 400; i++ {
		for j := 0; j <= 400; j++ {
			nu := [2]float64{float64(i) / 400, float64(j) / 400}
			e0, e1, e2 := p.energy([2]int{0, 0}, nu), p.energy([2]int{1, 0}, nu), p.energy([2]int{0, 1}, nu)
			if d := math.Abs(e0-e1) + math.Abs(e0-e2); d < best {
				best, at = d, nu
			}
		}
	}
	return at
}

// rate of one tunnelling event with free-energy change dF through R.
func (p SEPump) rate(dF, R float64) float64 {
	x := dF / (kB * p.T)
	if math.Abs(x) < 1e-9 {
		return kB * p.T / (E * E * R)
	}
	if x > 700 {
		return 0
	}
	return dF / (E * E * R * math.Expm1(x))
}

type seEvent struct {
	to   [2]int  // state after
	dF   float64 // free-energy change [J]
	rate float64 // [1/s]
	j    int     // junction 0, 1, 2
	dir  int     // +1: electron moves left → right through the junction
}

// events from state n at gate charges nu.
func (p SEPump) events(n [2]int, nu [2]float64) []seEvent {
	a, b, c := p.cmat()
	det := a*c - b*b
	inv := [2][2]float64{{c / det, b / det}, {b / det, a / det}}
	f1, f2 := p.phi(n, nu)
	pot := []float64{p.V / 2, f1, f2, -p.V / 2} // conductors: L, 1, 2, R
	// inverse-capacitance entry for conductors (leads have none)
	ci := func(x, y int) float64 {
		if x < 1 || x > 2 || y < 1 || y > 2 {
			return 0
		}
		return inv[x-1][y-1]
	}
	var ev []seEvent
	move := func(from, to, j, dir int) {
		nn := n
		if from >= 1 && from <= 2 {
			nn[from-1]--
		}
		if to >= 1 && to <= 2 {
			nn[to-1]++
		}
		if nn[0] < -nMax || nn[0] > nMax || nn[1] < -nMax || nn[1] > nMax {
			return
		}
		// electron (−e) from conductor `from` to `to`
		dF := -E*(pot[to]-pot[from]) + E*E/2*(ci(from, from)+ci(to, to)-2*ci(from, to))
		ev = append(ev, seEvent{to: nn, dF: dF, rate: p.rate(dF, p.R[j]), j: j, dir: dir})
	}
	for j := 0; j < 3; j++ { // junction j joins conductors j and j+1
		move(j, j+1, j, +1)
		move(j+1, j, j, -1)
	}
	return ev
}

// SEResult is one pumping run.
type SEResult struct {
	F        float64 // drive frequency [Hz]
	PerCycle float64 // electrons moved left → right per cycle (ideal 1)
	Current  float64 // [A]
	Heat     float64 // per cycle [J]
	GateWork float64 // per cycle [J]
	BiasWork float64 // per cycle [J]
	Residual float64 // ledger residual per cycle / heat
}

// Pump runs the master equation for `cycles` drive cycles at frequency f
// (after `settle` cycles), counting electrons through the middle junction.
func (p SEPump) Pump(f float64, settle, cycles int) SEResult {
	nS := 2*nMax + 1
	idx := func(n [2]int) int { return (n[0]+nMax)*nS + n[1] + nMax }
	P := make([]float64, nS*nS)
	P[idx([2]int{0, 0})] = 1
	nu := func(t float64) [2]float64 {
		w := 2 * math.Pi * f * t
		return [2]float64{p.Centre[0] + p.Amp*math.Cos(w), p.Centre[1] + p.Amp*math.Cos(w-p.Phase)}
	}
	// step size from the fastest rate met on the drive cycle
	gMax := 0.0
	for k := 0; k < 16; k++ {
		v := nu(float64(k) / 16 / f)
		for a := -nMax; a <= nMax; a++ {
			for b := -nMax; b <= nMax; b++ {
				for _, e := range p.events([2]int{a, b}, v) {
					if math.Abs(float64(a))+math.Abs(float64(b)) <= 2 {
						gMax = math.Max(gMax, e.rate)
					}
				}
			}
		}
	}
	period := 1 / f
	steps := int(math.Max(400, math.Ceil(period*gMax*2)))
	dt := period / float64(steps)
	deriv := func(t float64, P []float64, out []float64, flow *[3]float64, heat *float64) {
		for i := range out {
			out[i] = 0
		}
		v := nu(t)
		for a := -nMax; a <= nMax; a++ {
			for b := -nMax; b <= nMax; b++ {
				n := [2]int{a, b}
				pn := P[idx(n)]
				if pn == 0 {
					continue
				}
				for _, e := range p.events(n, v) {
					r := pn * e.rate
					out[idx(n)] -= r
					out[idx(e.to)] += r
					if flow != nil {
						flow[e.j] += float64(e.dir) * r
						*heat += -e.dF * r
					}
				}
			}
		}
	}
	// gate and lead charges for the source work: q = C(V_src − φ) on each
	// capacitor plate at a source, minus e per electron the lead gave away.
	plates := func(t float64, P []float64) (qL, qR, qg1, qg2 float64) {
		v := nu(t)
		for a := -nMax; a <= nMax; a++ {
			for b := -nMax; b <= nMax; b++ {
				n := [2]int{a, b}
				pn := P[idx(n)]
				if pn == 0 {
					continue
				}
				f1, f2 := p.phi(n, v)
				qL += pn * p.C[0] * (p.V/2 - f1)
				qR += pn * p.C[2] * (-p.V/2 - f2)
				qg1 += pn * p.Cg[0] * (v[0]*E/p.Cg[0] - f1)
				qg2 += pn * p.Cg[1] * (v[1]*E/p.Cg[1] - f2)
			}
		}
		return
	}
	stored := func(t float64, P []float64) float64 {
		v := nu(t)
		u := 0.0
		for a := -nMax; a <= nMax; a++ {
			for b := -nMax; b <= nMax; b++ {
				n := [2]int{a, b}
				pn := P[idx(n)]
				if pn == 0 {
					continue
				}
				f1, f2 := p.phi(n, v)
				vg1, vg2 := v[0]*E/p.Cg[0], v[1]*E/p.Cg[1]
				u += pn * 0.5 * (p.C[0]*sq(p.V/2-f1) + p.C[1]*sq(f1-f2) + p.C[2]*sq(f2+p.V/2) +
					p.Cg[0]*sq(vg1-f1) + p.Cg[1]*sq(vg2-f2))
			}
		}
		return u
	}
	k1, k2, k3, k4 := make([]float64, len(P)), make([]float64, len(P)), make([]float64, len(P)), make([]float64, len(P))
	tmp := make([]float64, len(P))
	t := 0.0
	var res SEResult
	res.F = f
	var nL, nR float64 // electrons the leads have given to the islands
	var flowMid, heat, wGate, wBias float64
	qL0, qR0, qg10, qg20 := plates(0, P)
	u0 := 0.0
	for c := 0; c < settle+cycles; c++ {
		if c == settle {
			flowMid, heat, wGate, wBias = 0, 0, 0, 0
			u0 = stored(t, P)
		}
		for s := 0; s < steps; s++ {
			var fl [4][3]float64
			var hh [4]float64
			deriv(t, P, k1, &fl[0], &hh[0])
			for i := range P {
				tmp[i] = P[i] + dt/2*k1[i]
			}
			deriv(t+dt/2, tmp, k2, &fl[1], &hh[1])
			for i := range P {
				tmp[i] = P[i] + dt/2*k2[i]
			}
			deriv(t+dt/2, tmp, k3, &fl[2], &hh[2])
			for i := range P {
				tmp[i] = P[i] + dt*k3[i]
			}
			deriv(t+dt, tmp, k4, &fl[3], &hh[3])
			for i := range P {
				P[i] += dt / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
			}
			var flow [3]float64
			for j := 0; j < 3; j++ {
				flow[j] = dt / 6 * (fl[0][j] + 2*fl[1][j] + 2*fl[2][j] + fl[3][j])
			}
			h := dt / 6 * (hh[0] + 2*hh[1] + 2*hh[2] + hh[3])
			t += dt
			// electrons from the left lead into island 1 = flow through J1 (+);
			// from the right lead into island 2 = −flow through J3.
			nL += flow[0]
			nR -= flow[2]
			qL, qR, qg1, qg2 := plates(t, P)
			qL -= E * nL
			qR -= E * nR
			v := nu(t - dt/2)
			vg1, vg2 := v[0]*E/p.Cg[0], v[1]*E/p.Cg[1]
			wBias += p.V/2*(qL-qL0) - p.V/2*(qR-qR0)
			wGate += vg1*(qg1-qg10) + vg2*(qg2-qg20)
			qL0, qR0, qg10, qg20 = qL, qR, qg1, qg2
			flowMid += flow[1]
			heat += h
		}
	}
	n := float64(cycles)
	res.PerCycle = flowMid / n
	res.Current = -E * flowMid / n * f // electrons left → right = conventional current right → left
	res.Heat, res.GateWork, res.BiasWork = heat/n, wGate/n, wBias/n
	res.Residual = (wGate + wBias - heat - (stored(t, P) - u0)) / heat
	return res
}

// ---------------------------------------------------------------------------
// 2. Multipactor pump

// Multipactor is a parallel-plate vacuum gap that is the capacitor of a
// parallel RLC resonator driven by an RF current source.
type Multipactor struct {
	Gap, Area float64 // [m], [m²]
	F, Q      float64 // resonance [Hz], unloaded Q
	V0        float64 // unloaded RF gap amplitude [V] (sets the drive)
	Grid      float64 // fraction of electrons reaching plate B passed to the collector
	Seed      float64 // starting electrons
	// Background: one stray electron (field emission, cosmic rays) appears
	// at plate A every SeedEvery periods while the gap is empty, at RF
	// phases spread by the golden ratio; 0 = none.
	SeedEvery int
	// Rext is the coupling port's resistance seen across the gap (the
	// generator side, through a circulator); 0 = none. Part of R is then
	// power returned out of the port rather than heat in the cavity.
	Rext float64
	// Steps per RF period (RK4); default 800.
	Steps   int
	Phase   float64 // RF phase at which the seed leaves plate A [rad]
	Surface em.SEY
}

// DefaultMultipactor: 1 mm copper gap (measured as-received yield), 1 cm²,
// 1 GHz, Q = 1000, a 10 % grid. The gap amplitude and the seed's launch
// phase are those that grow fastest (Tune); em.Multipactor puts the
// first-order zone of this gap near 100 V.
func DefaultMultipactor() Multipactor {
	m := Multipactor{Gap: 1e-3, Area: 1e-4, F: 1e9, Q: 1000, Grid: 0.1, Seed: 1000, SeedEvery: 20, Surface: em.CopperAsReceived}
	m.Tune(80, 160)
	return m
}

// Tune sets V0 (in 5 V steps between lo and hi) and Phase (5° steps) to
// the pair whose sheet grows most in 12 periods.
func (m *Multipactor) Tune(lo, hi float64) {
	best := 0.0
	bv, bp := m.V0, m.Phase
	for v := lo; v <= hi; v += 5 {
		for ph := 0.0; ph < math.Pi; ph += math.Pi / 36 {
			t := *m
			t.V0, t.Phase = v, ph
			r := t.Run(12, true)
			if n := r.Samples[len(r.Samples)-1].Electrons; n > best {
				best, bv, bp = n, v, ph
			}
		}
	}
	m.V0, m.Phase = bv, bp
}

// Tank values: C = ε0·A/d, L = 1/(ω²C), R = Q/(ωC), drive amplitude V0/R.
func (m Multipactor) Tank() (L, C, R, Is float64) {
	C = em.Eps0 * m.Area / m.Gap
	w := 2 * math.Pi * m.F
	L = 1 / (w * w * C)
	R = m.Q / (w * C)
	return L, C, R, m.V0 / R
}

// MPSample is the state once per RF period.
type MPSample struct {
	Period    int
	Electrons float64
	GapV      float64 // RF amplitude over that period [V]
}

// MPResult is a multipactor run.
type MPResult struct {
	Samples    []MPSample
	Collected  float64 // electrons collected
	Current    float64 // average collector current over the last half of the run [A]
	SourceWork float64 // [J]
	TankHeat   float64 // in R [J]
	WallHeat   float64 // electron impacts on the plates [J]
	Collector  float64 // electron energy delivered to the collector [J]
	Stored     float64 // change of resonator + electron kinetic energy [J]
	PortOut    float64 // part of TankHeat that went back out of the coupling port [J]
	PeakGapV   float64 // largest gap voltage met [V]
	StepsUsed  int     // RK4 steps per period actually used
	Strays     float64 // stray electrons added, and their 2 eV each [J] below
	StrayWork  float64
	Residual   float64 // ledger residual / source work
}

// Run follows the resonator and the electron sheet for `periods` RF
// periods, starting from the unloaded steady state with Seed electrons
// leaving plate A as the field turns to push them across. With electrons
// = false the tank runs alone (twin test against kronspice).
func (m Multipactor) Run(periods int, electrons bool) MPResult {
	L, C, R, Is := m.Tank()
	w := 2 * math.Pi * m.F
	d := m.Gap
	v0 := math.Sqrt(2 * 2 * E / mE) // 2 eV secondaries
	// steady state (source Is sin(ωt) in phase with V at resonance): V = V0 sin(ωt), IL = −V0 cos(ωt)/(ωL), started at phase m.Phase
	t := m.Phase / w
	V, IL := m.V0*math.Sin(m.Phase), -m.V0*math.Cos(m.Phase)/(w*L)
	x, vel, N := 0.0, v0, 0.0
	if electrons {
		N = m.Seed
	}
	steps := m.Steps
	if steps == 0 {
		steps = 800
	}
	// at least 40 steps per gap transit at the fastest possible electron
	// (all of V0 gained), so a coarse grid cannot invent impacts
	vmax := math.Sqrt(2*E*m.V0/mE) + v0
	if need := int(math.Ceil(40 * vmax / (m.F * d))); need > steps {
		steps = need
	}
	dt := 1 / m.F / float64(steps)
	stored := func(V, IL float64) float64 { return 0.5*C*V*V + 0.5*L*IL*IL }
	ke := func(N, v float64) float64 { return 0.5 * mE * v * v * N }
	var res MPResult
	res.StepsUsed = steps
	u0 := stored(V, IL) + ke(N, vel)
	type st struct{ V, IL, x, v float64 }
	der := func(t float64, s st) st {
		ib := 0.0
		if N > 0 {
			ib = E * N * s.v / d // Ramo: current the moving sheet draws from the gap
		}
		a := E * s.V / (mE * d) // V = plate B minus plate A: pushes electrons toward B
		return st{V: (Is*math.Sin(w*t) - s.V/R - s.IL - ib) / C, IL: s.V / L, x: s.v, v: a}
	}
	add := func(s st, k st, h float64) st { return st{s.V + h*k.V, s.IL + h*k.IL, s.x + h*k.x, s.v + h*k.v} }
	var collectedLate float64
	seedStep := int(m.Phase/(2*math.Pi)*float64(steps)) % steps
	for p := 0; p < periods; p++ {
		amp := 0.0
		for s := 0; s < steps; s++ {
			// strays arrive at phases spread by the golden ratio (not all
			// at the best phase)
			if electrons && N == 0 && m.SeedEvery > 0 && p%m.SeedEvery == 0 && s == (seedStep+int(math.Mod(res.Strays*0.6180339887, 1)*float64(steps)))%steps {
				N, x, vel = 1, 0, v0
				res.Strays++
				res.StrayWork += 0.5 * mE * v0 * v0
			}
			y := st{V, IL, x, vel}
			k1 := der(t, y)
			k2 := der(t+dt/2, add(y, k1, dt/2))
			k3 := der(t+dt/2, add(y, k2, dt/2))
			k4 := der(t+dt, add(y, k3, dt))
			n := st{
				y.V + dt/6*(k1.V+2*k2.V+2*k3.V+k4.V), y.IL + dt/6*(k1.IL+2*k2.IL+2*k3.IL+k4.IL),
				y.x + dt/6*(k1.x+2*k2.x+2*k3.x+k4.x), y.v + dt/6*(k1.v+2*k2.v+2*k3.v+k4.v),
			}
			// energy flows over the step (Simpson, matching RK4)
			pv := func(tt, v float64) float64 { return Is * math.Sin(w*tt) * v }
			vm := 0.5 * (y.V + n.V)
			res.SourceWork += dt / 6 * (pv(t, y.V) + 4*pv(t+dt/2, vm) + pv(t+dt, n.V))
			res.TankHeat += dt / 6 * (y.V*y.V + 4*vm*vm + n.V*n.V) / R
			if m.Rext > 0 {
				res.PortOut += dt / 6 * (y.V*y.V + 4*vm*vm + n.V*n.V) / m.Rext
			}
			V, IL, x, vel = n.V, n.IL, n.x, n.v
			t += dt
			amp = math.Max(amp, math.Abs(V))
			if N > 0 && (x <= 0 || x >= d) {
				eImp := 0.5 * mE * vel * vel // [J]
				arrive := N
				if x >= d {
					c := m.Grid * N
					res.Collected += c
					res.Collector += c * eImp
					if p >= periods/2 {
						collectedLate += c
					}
					arrive -= c
				}
				newN := m.Surface.Yield(eImp/E) * arrive
				// space charge: a sheet of charge σ leaves the wall only while
				// the RF field exceeds its own field σ/ε0
				if cap := em.Eps0 * m.Area * math.Abs(V) / (d * E); newN > cap {
					newN = cap
				}
				res.WallHeat += arrive*eImp - newN*0.5*mE*v0*v0
				if newN < 1 { // fewer than one electron: the discharge is out
					res.WallHeat += newN * 0.5 * mE * v0 * v0
					newN = 0
				}
				N = newN
				if x <= 0 {
					x, vel = 0, v0
				} else {
					x, vel = d, -v0
				}
			}
		}
		res.Samples = append(res.Samples, MPSample{Period: p + 1, Electrons: N, GapV: amp})
		res.PeakGapV = math.Max(res.PeakGapV, amp)
	}
	late := float64(periods-periods/2) / m.F
	res.Current = E * collectedLate / late
	res.Stored = stored(V, IL) + ke(N, vel) - u0
	res.Residual = (res.SourceWork + res.StrayWork - res.TankHeat - res.WallHeat - res.Collector - res.Stored) / res.SourceWork
	return res
}

// TankTwin runs the bare resonator in kronspice (UIC from the same
// steady state) and returns the largest difference in gap voltage from
// Run(…, false) over `periods`, relative to V0.
func (m Multipactor) TankTwin(periods int) (float64, error) {
	m.Phase = 0
	L, C, R, Is := m.Tank()
	w := 2 * math.Pi * m.F
	b := circuit.NewBuilder("multipactor resonator, bare")
	b.CapacitorIC("CGAP", "g", "0", C, 0)
	b.InductorIC("LT", "g", "0", L, -m.V0/(w*L))
	b.Resistor("RT", "g", "0", R)
	b.CurrentWave("IS", "0", "g", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, Is, m.F}})
	ckt, err := b.BuildTransient()
	if err != nil {
		return 0, err
	}
	stop := float64(periods) / m.F
	step := 1 / m.F / 400
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return 0, err
	}
	vg, _ := wv.NodeVoltage("g")
	// own solver, sampled at the same times
	// the same integrator as Run, recording V(t)
	steps := 800
	dt := 1 / m.F / float64(steps)
	V, IL, t := 0.0, -m.V0/(w*L), 0.0
	var tt, vv []float64
	for s := 0; s < periods*steps; s++ {
		f := func(t, V, IL float64) (float64, float64) {
			return (Is*math.Sin(w*t) - V/R - IL) / C, V / L
		}
		a1, b1 := f(t, V, IL)
		a2, b2 := f(t+dt/2, V+dt/2*a1, IL+dt/2*b1)
		a3, b3 := f(t+dt/2, V+dt/2*a2, IL+dt/2*b2)
		a4, b4 := f(t+dt, V+dt*a3, IL+dt*b3)
		V += dt / 6 * (a1 + 2*a2 + 2*a3 + a4)
		IL += dt / 6 * (b1 + 2*b2 + 2*b3 + b4)
		t += dt
		tt, vv = append(tt, t), append(vv, V)
	}
	worst := 0.0
	for i, ts := range wv.Time {
		k := int(ts/dt+0.5) - 1
		if k < 0 || k >= len(vv) {
			continue
		}
		worst = math.Max(worst, math.Abs(vg[i]-vv[k]))
	}
	return worst / m.V0, nil
}

// ---------------------------------------------------------------------------
// 3. Ground charge pump

// GroundPump: a 10 V, 100 kHz source drives a series resonator (L, Ct,
// unloaded Q = 20) whose capacitor voltage feeds an n-stage negative
// Cockcroft–Walton ladder; the output charges a storage bank with a 10 MΩ
// load, all returned to earth (node 0). Earth lends the electrons that
// end up on the bank's negative plate, and the bank's earthed plate gives
// the same number back, so the earth's net charge does not change (§3j,
// energy versus charge); the energy comes only from the 10 V source.
type GroundPump struct {
	Vs, F, Q     float64
	Ct           float64 // tank capacitor [F]
	Stages       int
	Cstage       float64 // ladder capacitors [F]
	Cbank, Rload float64
	Cycles       int
}

// DefaultGroundPump as described above: 3 stages of 100 nF, a 100 nF
// bank, 600 cycles.
func DefaultGroundPump() GroundPump {
	return GroundPump{Vs: 10, F: 100e3, Q: 20, Ct: 1e-6, Stages: 3, Cstage: 100e-9, Cbank: 100e-9, Rload: 10e6, Cycles: 600}
}

// GPResult is a ground-pump run.
type GPResult struct {
	TankPeak, Bank, Ideal  float64 // V: tank peak, final bank voltage, 2·n·tank peak
	BankCharge, LoadCharge float64 // on the bank's negative plate, and passed through the load [C]
	BankEnergy             float64 // ½·Cbank·V² [J]
	EnergyIn, Heat, Stored float64 // [J]
	Mismatch               float64 // time-stepping error [J]
	Time, BankV            []float64
}

// Build the circuit.
func (g GroundPump) Build() (*circuit.Circuit, error) {
	w := 2 * math.Pi * g.F
	L := 1 / (w * w * g.Ct)
	Rs := math.Sqrt(L/g.Ct) / g.Q
	b := circuit.NewBuilder("resonant-driven ground charge pump")
	b.Model("dhv", "D", map[string]float64{"is": 1e-9, "n": 1.8, "rs": 0.05})
	b.VoltageWave("VS", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, g.Vs, g.F}})
	b.Resistor("RS", "s", "m", Rs)
	b.Inductor("LT", "m", "ac", L)
	b.Capacitor("CT", "ac", "0", g.Ct)
	// negative Cockcroft–Walton: column of "push" caps from ac, column of
	// "hold" caps from ground; diodes alternate down the ladder.
	prevPush, prevHold := "ac", "0"
	for k := 1; k <= g.Stages; k++ {
		push, hold := fmt.Sprintf("p%d", k), fmt.Sprintf("h%d", k)
		b.Capacitor(fmt.Sprintf("CP%d", k), prevPush, push, g.Cstage)
		b.Capacitor(fmt.Sprintf("CH%d", k), prevHold, hold, g.Cstage)
		b.Diode(fmt.Sprintf("DA%d", k), push, prevHold, "dhv") // electrons from the hold column into push
		b.Diode(fmt.Sprintf("DB%d", k), hold, push, "dhv")
		prevPush, prevHold = push, hold
	}
	b.Capacitor("CBANK", prevHold, "0", g.Cbank)
	b.Resistor("RLOAD", prevHold, "0", g.Rload)
	return b.BuildTransient()
}

// Run simulates the pump.
func (g GroundPump) Run() (GPResult, error) {
	ckt, err := g.Build()
	if err != nil {
		return GPResult{}, err
	}
	stop := float64(g.Cycles) / g.F
	step := 1 / g.F / 200
	wv, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, analysis.TranOptions{})
	if err != nil {
		return GPResult{}, err
	}
	ac, _ := wv.NodeVoltage("ac")
	out, _ := wv.NodeVoltage(fmt.Sprintf("h%d", g.Stages))
	r := GPResult{Time: wv.Time, BankV: out}
	n := len(wv.Time)
	for i := n * 3 / 4; i < n; i++ {
		r.TankPeak = math.Max(r.TankPeak, math.Abs(ac[i]))
	}
	r.Bank = out[n-1]
	r.Ideal = -2 * float64(g.Stages) * r.TankPeak
	for i := 1; i < n; i++ {
		dt := wv.Time[i] - wv.Time[i-1]
		r.LoadCharge += 0.5 * (out[i] + out[i-1]) / g.Rload * dt
	}
	r.BankCharge = g.Cbank * r.Bank
	r.BankEnergy = 0.5 * g.Cbank * r.Bank * r.Bank
	l := wv.Ledger
	r.EnergyIn, r.Heat, r.Stored, r.Mismatch = l.EnergyIn, l.EnergyDissipated, l.StoredEnd-l.StoredStart, l.IntegrationMismatch
	return r, nil
}
