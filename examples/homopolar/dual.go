// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// The dual motor: two of the copper disks above, each on its own motor,
// cross-excited as in Rikitake's machine (examples/earth_machine) — each
// disk's current runs through the field coil of the OTHER unit — with a
// feedback coil on each unit carrying its own current (the feedback inductor
// of Evans, Amador & Eckardt, paper107 §14.5.2.3), side by side or stacked
// on one axis. No magnets and no separate excitation: the machine must
// excite itself.
//
// Every number comes from the geometry:
//
//   - a disk's EMF is ω·Φ/2π, Φ the flux through the annulus r_i … r_o, so
//     the coupling per ampere of coil current is
//     M = (N/2π)·⟨M_circ(r_c, r_o, u) − M_circ(r_c, r_i, u)⟩ over the turns'
//     axial positions u (em.MutualCircles, the disk in the coil's mid-plane);
//   - a coil's inductance is Nagaoka's (em.Nagaoka); its resistance
//     ρ·2πr_c·N²/A_w for a copper window A_w;
//   - side by side, the feedback coil has K·N turns on the same former in a
//     window K·A_w: coupling K·M, inductance K²·L, resistance |K|·R_coil,
//     and mutual inductance −k·K·L with the cross coil on that disk (k
//     declared);
//   - stacked, disk 2 sits inside coil 1 and disk 1 inside coil 2 one gap
//     above: each coil also threads its OWN disk through its stray field
//     (MFar, the same circle formula a gap away) and the two coils couple
//     (MCoils): a feedback coil built into the geometry, without copper;
//   - a shaft's inertia is the disk's ½·m·r_o² plus a declared motor rotor.
//
// The machine element books ω·G·i as a drop along each winding (see
// examples/earth_machine/main.go): a winding that runs against its loop
// turns the sign. Side by side, winding 0 is the disk path, winding 1 the
// cross coil (the other loop, wound against it), winding 2 the feedback coil
// (this loop, along it): G[0][1] = +M, G[0][2] = −K·M. Stacked, the coils are
// inductors with a Mutual element and the machine elements carry 1 nH sense
// windings along their loops: G[0][1] = −M, G[0][2] = −MFar.

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
)

// Dual is the cross-excited two-motor machine.
type Dual struct {
	Disk    Design  // disk geometry and shaft contacts (Defaults)
	Rc, Lc  float64 // coil mean radius and axial length [m]
	N       float64 // cross-coil turns
	Aw      float64 // copper area of the cross coil's winding window [m²]
	K       float64 // side by side: feedback-coil turns ratio (> 0 aids, < 0 opposes)
	Kc      float64 // coupling coefficient of the two coils on one former (declared)
	JMotor  float64 // motor rotor inertia per shaft (declared) [kg·m²]
	Torque  [2]float64
	Speed   *[2]float64 // when set, each shaft is held at this speed (VFD) [rad/s]
	I0, W0  [2]float64  // initial currents [A] and speeds [rad/s]
	Windage float64     // viscous loss per shaft [N·m·s] (declared)
	Stacked bool        // disks on one axis, a gap apart (else side by side)
	Gap     float64     // distance between the disk planes when stacked [m]
	// RLoad is a load in each loop [Ω]. CParallel is a capacitor across the
	// load [F] (0: none), with its series resistance ESR [Ω]: the load and
	// the capacitor make a parallel tank fed by the coils' inductance
	// (resonance.go). DC still reaches the load; the ESR is outside the DC
	// path.
	RLoad, CParallel, ESR float64
	// LPrim, when set, puts the primary of a transformer (winding resistance
	// RPrim) across the load too, so the output is R_L ∥ C ∥ L_p; the
	// secondary LSec (coupling KT) feeds an AC load RAC. The primary is a
	// short at DC: it carries the DC, and R_L sees only AC.
	LPrim, RPrim, LSec, KT, RAC float64
	// CBank > 0 rectifies the secondary instead: a diode bridge into a
	// capacitor bank CBank with a DC load RDC across it (engine only; the
	// bridge is nonlinear, so RunRK4 does not model it).
	CBank, RDC, BankV0 float64 // BankV0: the bank's starting voltage
	// Output (stacked, with a tank) chooses how the machines meet the
	// load: a tank per loop (OutPerLoop), one tank shared across both
	// (OutShared) or one tank between the two coils (OutSeries). The field
	// values above are always one loop's; a shared or series tank is built
	// from them (see OutputMode).
	Output OutputMode
	// TorqueAmp [N·m] and TorqueFreq [Hz] modulate both motor torques,
	// G(t) = G + TorqueAmp·sin(2π·TorqueFreq·t), so the shaft speed
	// oscillates; SpeedAmp (relative) and SpeedFreq [Hz] modulate held
	// speeds, ω(t) = ω·(1 + SpeedAmp·sin(2π·SpeedFreq·t)).
	TorqueAmp, TorqueFreq float64
	SpeedAmp, SpeedFreq   float64
	SpeedPhase2           float64 // phase of shaft 2's speed oscillation [degrees] (180: antiphase)
	// IronGain multiplies every coupling and inductance of the coils, as an
	// iron magnetic circuit would below saturation (declared; 0 or 1: air).
	IronGain float64
}

// OutputMode is how the loops feed the output tank. In OutShared and
// OutSeries the two shafts are the common (ground) node (a homopolar disk's
// EMF is between axle and rim), and the transformer secondary has its common
// end tied to it, the bridge's DC side floating.
type OutputMode int

const (
	// OutPerLoop: each loop has its own R_L ∥ C ∥ L_p across its terminals.
	OutPerLoop OutputMode = iota
	// OutShared: each rim feeds R_int and its coil to one node, and one tank
	// sits between that node and the shafts: the two per-loop tanks in
	// parallel (impedances ÷ 2, capacitances × 2). The loops' sum passes
	// through the tank, their difference circulates rim to rim.
	OutShared
	// OutSeries: one current, shaft → disk 1 → R_int → coil 1 → tank →
	// coil 2 → R_int → disk 2 → shaft; the tank sits between the coils'
	// output ends: the two per-loop tanks in series (impedances × 2,
	// capacitances ÷ 2). The disks' EMFs add round the loop only if shaft 2
	// counter-rotates (Speed[1] or Torque[1] negative) and coil 2 is wound to
	// aid coil 1, as built here. With one current there is no difference
	// mode: a series self-excited dynamo with the per-loop threshold.
	OutSeries
)

func (m Dual) gain() float64 {
	if m.IronGain > 0 {
		return m.IronGain
	}
	return 1
}

// DefaultDual: the reference disks in single-layer copper coils of four
// turns, a 6 cm × 2 cm window (declared), and 20 N·m motors.
func DefaultDual() Dual {
	d := Defaults()
	d.RShaft = 5e-6 // declared: liquid-metal shaft contacts, as in high-current homopolar machines
	return Dual{Disk: d, Rc: 0.18, Lc: 0.06, N: 4, Aw: 0.06 * 0.02, K: 0, Kc: 0.95, JMotor: 0.02,
		Torque: [2]float64{20, 20}, I0: [2]float64{100, 50}, W0: [2]float64{900, 600}, Windage: 1e-4, Gap: 0.08}
}

// shifted is the annulus coupling of a coil whose mid-plane is z away from
// the disk: the disk's EMF per rad/s per ampere [V·s].
func (m Dual) shifted(z float64) float64 {
	const n = 64
	var s float64
	for k := 0; k < n; k++ {
		u := z + m.Lc*((float64(k)+0.5)/n-0.5)
		s += em.MutualCircles(m.Rc, m.Disk.Ro, u) - em.MutualCircles(m.Rc, m.Disk.Ri, u)
	}
	return m.gain() * m.N / (2 * math.Pi) * s / n
}

// Coupling is M, a coil's coupling to the disk in its own mid-plane.
func (m Dual) Coupling() float64 { return m.shifted(0) }

// MFar is a coil's coupling to the disk one gap away (stacked).
func (m Dual) MFar() float64 { return m.shifted(m.Gap) }

// MCoils is the mutual inductance of the two coaxial coils a gap apart.
func (m Dual) MCoils() float64 {
	const n = 32
	var s float64
	for a := 0; a < n; a++ {
		for b := 0; b < n; b++ {
			s += em.MutualCircles(m.Rc, m.Rc, m.Gap+m.Lc*float64(a-b)/n)
		}
	}
	return m.gain() * m.N * m.N * s / (n * n)
}

// Coil returns the cross coil's inductance and resistance.
func (m Dual) Coil() (L, R float64) {
	L = m.gain() * em.Nagaoka(m.Rc, m.Lc, m.N)
	R = material.Copper.Resistivity * 2 * math.Pi * m.Rc * m.N * m.N / m.Aw
	return
}

// LoopR is one loop's resistance: disk, two shaft contacts, cross coil and,
// side by side, the feedback coil.
func (m Dual) LoopR() float64 {
	_, rc := m.Coil()
	if m.Stacked {
		return m.Disk.DiskResistance() + 2*m.Disk.RShaft + rc + m.RLoad
	}
	return m.Disk.DiskResistance() + 2*m.Disk.RShaft + rc*(1+math.Abs(m.K)) + m.RLoad
}

// loop returns the loop equations' constants: self and mutual inductance of
// the two loops, the self-coupling (feedback) as a fraction of M, and R.
func (m Dual) loop() (l11, l12, kself, r float64) {
	L, _ := m.Coil()
	if m.Stacked {
		return L, m.MCoils(), m.MFar() / m.Coupling(), m.LoopR()
	}
	return L * (1 + m.K*m.K), 2 * m.Kc * m.K * L, m.K, m.LoopR()
}

// Inertia is one shaft's moment of inertia.
func (m Dual) Inertia() float64 {
	mass := material.Copper.Density * math.Pi * (m.Disk.Ro*m.Disk.Ro - m.Disk.Ri*m.Disk.Ri) * m.Disk.Thickness
	return 0.5*mass*m.Disk.Ro*m.Disk.Ro + m.JMotor
}

// Threshold is the common shaft speed above which the machine excites
// itself: (−R + (1 + K)·M·ω)·I = 0 for the symmetric mode.
func (m Dual) Threshold() float64 {
	_, _, ks, r := m.loop()
	return r / ((1 + ks) * m.Coupling())
}

// Mu is Rikitake's dimensionless resistance μ = R·√(J/(G·M·L)): the machine
// reverses near μ ≈ 1 and settles to a steady field when μ is large.
func (m Dual) Mu() float64 {
	l11, _, _, r := m.loop()
	return r * math.Sqrt(m.Inertia()/(m.Torque[0]*m.Coupling()*l11))
}

// Circuit builds the machine in its layout.
func (m Dual) Circuit() (*circuit.Circuit, error) {
	if m.Stacked {
		return m.stackedCircuit()
	}
	return m.sideCircuit()
}

// sideCircuit: three windings per machine element.
//
//	loop i: a_i →[N_i disk]→ b_i →R_i→ c_i →[N_i feedback coil]→ d_i →[N_j cross coil i]→ a_i
func (m Dual) sideCircuit() (*circuit.Circuit, error) {
	M := m.Coupling()
	L, _ := m.Coil()
	R := m.LoopR()
	spec := func(iOwn, iOther float64) *circuit.MachineSpec {
		s := circuit.NewMachineSpec(3)
		s.L0[0][0] = 1e-9 // disk path
		s.L0[1][1] = L    // cross coil (the other loop)
		s.L0[2][2] = m.K * m.K * L
		if m.K == 0 {
			s.L0[2][2] = 1e-9 // no feedback coil: a short, inert winding
		}
		s.L0[1][2] = -m.Kc * m.K * L // same former, wound for the field signs below
		s.L0[2][1] = s.L0[1][2]
		s.G[0][1] = M
		s.G[0][2] = -m.K * M
		s.IC[0], s.IC[1], s.IC[2] = iOwn, -iOther, iOwn
		return s
	}
	b := circuit.NewBuilder("dual-motor cross-excited homopolar machine").Domain(exchange.Mechanical, "w1", "w2")
	b.MachineModel("u1", spec(m.I0[0], m.I0[1])).MachineModel("u2", spec(m.I0[1], m.I0[0]))
	b.Machine("N1", "u1", [][2]string{{"a1", "b1"}, {"a2", "d2"}, {"c1", "d1"}}, "w1", "0")
	b.Machine("N2", "u2", [][2]string{{"a2", "b2"}, {"a1", "d1"}, {"c2", "d2"}}, "w2", "0")
	b.Resistor("R1", "b1", "c1", R).Resistor("R2", "b2", "c2", R)
	for _, n := range []string{"a1", "b1", "c1", "d1", "a2", "b2", "c2", "d2"} {
		b.Resistor("RREF"+n, n, "0", 1e6)
	}
	m.shafts(b)
	return b.Build()
}

// stackedCircuit: the coils as inductors LC1, LC2 with their mutual, the
// machine elements carrying only 1 nH sense windings for the speed EMFs.
//
//	loop i: a_i →[N_i disk]→ b_i →R_i→ c_i →[N_i sense of coil i]→ d_i →[N_j sense of coil i]→ e_i →LC_i→ a_i
func (m Dual) stackedCircuit() (*circuit.Circuit, error) {
	M, MF := m.Coupling(), m.MFar()
	L, _ := m.Coil()
	R := m.LoopR()
	spec := func(iOwn, iOther float64) *circuit.MachineSpec {
		s := circuit.NewMachineSpec(3)
		s.L0[0][0], s.L0[1][1], s.L0[2][2] = 1e-9, 1e-9, 1e-9
		// All windings run along their loops: w1 senses the OTHER loop's
		// coil (the main field on this disk), w2 this loop's coil (its stray
		// field on this disk).
		s.G[0][1] = -M
		s.G[0][2] = -MF
		s.IC[0], s.IC[1], s.IC[2] = iOwn, iOther, iOwn
		return s
	}
	b := circuit.NewBuilder("dual-motor homopolar machine, disks stacked").Domain(exchange.Mechanical, "w1", "w2")
	if m.CBank > 0 {
		// declared: a power rectifier diode (saturation current, emission coefficient, series resistance)
		b.Model("drect", "D", map[string]float64{"is": 1e-9, "n": 1.5, "rs": 2e-3})
	}
	if m.series() {
		// loop: 0 (shafts) →[N1 disk]→ b1 →R1→ c1 →[senses]→ e1 →LC1→ g1 →tank→ g2 ←LC2← e2 ←[senses]← c2 ←R2← b2 ←[N2 disk]← 0
		// Loop 2 is run backwards (its current is −I) and coil 2 is wound
		// the other way, so the coefficients on coil 2's sense windings and
		// the coils' mutual change sign: the coils' fields aid.
		i := m.I0[0]
		s1, s2 := spec(i, -i), spec(-i, i)
		s1.G[0][1] = M  // coil 2's main field on disk 1
		s2.G[0][2] = MF // coil 2's stray field on disk 2
		b.MachineModel("u1", s1).MachineModel("u2", s2)
		b.Machine("N1", "u1", [][2]string{{"0", "b1"}, {"d2", "e2"}, {"c1", "d1"}}, "w1", "0")
		b.Machine("N2", "u2", [][2]string{{"0", "b2"}, {"d1", "e1"}, {"c2", "d2"}}, "w2", "0")
		b.Resistor("R1", "b1", "c1", R-m.RLoad).Resistor("R2", "b2", "c2", R-m.RLoad)
		// LC2 entered from g2 (coil 2 reversed), carrying +I in the loop's direction: the coils' mutual aids
		b.InductorIC("LC1", "e1", "g1", L, i).InductorIC("LC2", "g2", "e2", L, i)
		b.Mutual("KC", "LC1", "LC2", m.MCoils()/L)
		m.addTank(b, "", "g1", "g2", 0.5, i, true)
		for _, n := range []string{"b1", "c1", "d1", "e1", "b2", "c2", "d2", "e2", "g2"} {
			b.Resistor("RREF"+n, n, "0", 1e3)
		}
		m.shafts(b)
		return b.Build()
	}
	b.MachineModel("u1", spec(m.I0[0], m.I0[1])).MachineModel("u2", spec(m.I0[1], m.I0[0]))
	if m.shared() {
		// loop i: 0 (shafts, common) →[N_i disk]→ b_i (rim) →R_i→ c_i →[senses]→ e_i →LC_i→ g →tank→ 0
		b.Machine("N1", "u1", [][2]string{{"0", "b1"}, {"d2", "e2"}, {"c1", "d1"}}, "w1", "0")
		b.Machine("N2", "u2", [][2]string{{"0", "b2"}, {"d1", "e1"}, {"c2", "d2"}}, "w2", "0")
		b.Resistor("R1", "b1", "c1", R-m.RLoad).Resistor("R2", "b2", "c2", R-m.RLoad)
		b.InductorIC("LC1", "e1", "g", L, m.I0[0]).InductorIC("LC2", "e2", "g", L, m.I0[1])
		b.Mutual("KC", "LC1", "LC2", m.MCoils()/L)
		m.addTank(b, "", "g", "0", 2, m.I0[0]+m.I0[1], true)
		for _, n := range []string{"b1", "c1", "d1", "e1", "b2", "c2", "d2", "e2"} {
			b.Resistor("RREF"+n, n, "0", 1e3)
		}
		m.shafts(b)
		return b.Build()
	}
	b.Machine("N1", "u1", [][2]string{{"a1", "b1"}, {"d2", "e2"}, {"c1", "d1"}}, "w1", "0")
	b.Machine("N2", "u2", [][2]string{{"a2", "b2"}, {"d1", "e1"}, {"c2", "d2"}}, "w2", "0")
	if (m.CParallel > 0 || m.LPrim > 0) && m.RLoad > 0 {
		// The load gets its own element with the capacitor (and the
		// transformer primary) across it:
		// b_i →R_i (internal)→ g_i →(RL_i ∥ [CP_i + ESR] ∥ [RPRIM_i + LP_i])→ c_i.
		for k, n := range []string{"1", "2"} {
			b.Resistor("R"+n, "b"+n, "g"+n, R-m.RLoad)
			m.addTank(b, n, "g"+n, "c"+n, 1, m.I0[k], false)
		}
	} else {
		b.Resistor("R1", "b1", "c1", R).Resistor("R2", "b2", "c2", R)
	}
	b.InductorIC("LC1", "e1", "a1", L, m.I0[0]).InductorIC("LC2", "e2", "a2", L, m.I0[1])
	b.Mutual("KC", "LC1", "LC2", m.MCoils()/L)
	// The loops float; reference them. 1 kΩ against the loop's 10⁻⁴ Ω leaks
	// a few parts in 10⁸, and keeps the UIC start solvable when a capacitor
	// is pinned by its large start conductance.
	for _, n := range []string{"a1", "b1", "c1", "d1", "e1", "a2", "b2", "c2", "d2", "e2"} {
		b.Resistor("RREF"+n, n, "0", 1e3)
	}
	m.shafts(b)
	return b.Build()
}

// held returns both shafts held at speed w; with a series tank shaft 2
// counter-rotates.
func (m Dual) held(w float64) *[2]float64 {
	if m.series() {
		return &[2]float64{w, -w}
	}
	return &[2]float64{w, w}
}

func (m Dual) hasTank() bool { return m.Stacked && (m.CParallel > 0 || m.LPrim > 0) && m.RLoad > 0 }

// shared reports whether both loops feed one tank in parallel, series
// whether they feed one tank between their coils.
func (m Dual) shared() bool { return m.Output == OutShared && m.hasTank() }
func (m Dual) series() bool { return m.Output == OutSeries && m.hasTank() }

// TankScale is how the built tank relates to the per-loop values: its
// impedances are z × theirs (capacitances ÷ z), and it carries i × one
// loop's current.
func (m Dual) TankScale() (z, i float64) {
	switch {
	case m.shared():
		return 0.5, 2
	case m.series():
		return 2, 1
	}
	return 1, 1
}

// addTank puts R_L ∥ [C + ESR] ∥ [R_p + L_p] between nodes g and c, with
// the transformer's secondary into R_ac or a bridge, bank and DC load. s
// scales the per-loop values (1: one loop's tank; 2: one tank for both
// loops in parallel, 0.5 in series: impedances ÷ s and capacitances × s);
// i0 is its starting DC current. secCommon ties the secondary's common end to
// the shafts' common (ground); the bridge's DC side then floats.
func (m Dual) addTank(b *circuit.Builder, n, g, c string, s, i0 float64, secCommon bool) {
	rl := m.RLoad / s
	b.Resistor("RL"+n, g, c, rl)
	v0 := rl * i0 // start the capacitor at the output's DC voltage
	if m.LPrim > 0 {
		rp := m.RPrim / s
		v0 = rp * i0
		if rp <= 0 {
			rp = 1e-9
		}
		b.Resistor("RPRIM"+n, g, "p"+n, rp).InductorIC("LP"+n, "p"+n, c, m.LPrim/s, i0)
		tn := "t" + n
		if secCommon {
			tn = "0"
		}
		if m.CBank > 0 {
			// secondary s – t into a bridge: DA s→dp, DB t→dp, DC dn→s, DD dn→t; bank and load dp–dn
			// declared: secondary winding 0.05 Ω; RC snubber (10 Ω + 100 µF) across the bridge input (per loop)
			b.InductorIC("LS"+n, "s"+n, "sw"+n, m.LSec/s, 0).Resistor("RSEC"+n, "sw"+n, tn, 0.05/s)
			b.Resistor("RSN"+n, "s"+n, "x"+n, 10/s).Capacitor("CSN"+n, "x"+n, tn, 100e-6*s)
			b.Resistor("RREFsw"+n, "sw"+n, "0", 1e3).Resistor("RREFx"+n, "x"+n, "0", 1e3)
			b.Diode("DA"+n, "s"+n, "dp"+n, "drect").Diode("DB"+n, tn, "dp"+n, "drect")
			b.Diode("DC"+n, "dn"+n, "s"+n, "drect").Diode("DD"+n, "dn"+n, tn, "drect")
			b.CapacitorIC("CB"+n, "dp"+n, "dn"+n, m.CBank*s, m.BankV0).Resistor("RDC"+n, "dp"+n, "dn"+n, m.RDC/s)
			b.Resistor("RREFdn"+n, "dn"+n, "0", 1e3).Resistor("RREFs"+n, "s"+n, "0", 1e3)
			if !secCommon {
				b.Resistor("RREFt"+n, tn, "0", 1e3)
			}
		} else {
			b.InductorIC("LS"+n, "s"+n, "0", m.LSec/s, 0).Resistor("RAC"+n, "s"+n, "0", m.RAC/s)
		}
		b.Mutual("KT"+n, "LP"+n, "LS"+n, m.KT)
		b.Resistor("RREFp"+n, "p"+n, "0", 1e3)
	}
	if m.CParallel > 0 {
		if m.ESR > 0 {
			b.CapacitorIC("CP"+n, g, "h"+n, m.CParallel*s, v0).Resistor("RESR"+n, "h"+n, c, m.ESR/s)
			b.Resistor("RREFh"+n, "h"+n, "0", 1e3)
		} else {
			b.CapacitorIC("CP"+n, g, c, m.CParallel*s, v0)
		}
	}
	if g != "0" {
		b.Resistor("RREFg"+n, g, "0", 1e3)
	}
}

// shafts adds the two shafts: inertia, windage and a motor torque each, or
// a speed drive (VFD) holding each at a set speed.
func (m Dual) shafts(b *circuit.Builder) {
	for k, w := range []string{"w1", "w2"} {
		n := fmt.Sprint(k + 1)
		if m.Speed != nil {
			e := circuit.Element{Name: "VFD" + n, Kind: circuit.VoltageSource, Value: m.Speed[k], HasDC: true, Channel: exchange.Mechanical}
			if m.SpeedAmp != 0 {
				ph := 0.0
				if k == 1 {
					ph = m.SpeedPhase2
				}
				e.Wave = circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{m.Speed[k], m.SpeedAmp * m.Speed[k], m.SpeedFreq, 0, 0, ph}}
			}
			b.Add(e, w, "0")
			continue
		}
		b.CapacitorIC("J"+n, w, "0", m.Inertia(), m.W0[k])
		b.Resistor("RW"+n, w, "0", 1/m.Windage)
		e := circuit.Element{Name: "MOTOR" + n, Kind: circuit.CurrentSource, Value: m.Torque[k], HasDC: true, Channel: exchange.Mechanical}
		if m.TorqueAmp != 0 {
			e.Wave = circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{m.Torque[k], m.TorqueAmp, m.TorqueFreq, 0, 0, 0}}
		}
		b.Add(e, "0", w)
	}
}

// DualTrace is a run.
type DualTrace struct {
	T, I1, I2, W1, W2 []float64
	V1, IP1, IS1      []float64 // loop 1's output voltage, primary and secondary currents (RK4, with an output tank)
	Mismatch, In      float64
}

// Run integrates the machine in the engine.
func (m Dual) Run(tEnd, dt float64) (DualTrace, error) {
	c, err := m.Circuit()
	if err != nil {
		return DualTrace{}, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: tEnd, MaxStep: dt, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return DualTrace{}, err
	}
	var tr DualTrace
	tr.T = wf.Time
	i1, _ := wf.ElementCurrent("R1")
	i2, _ := wf.ElementCurrent("R2")
	// R_i runs b_i → c_i, along loop i: its current is the loop current
	// (checked against RunRK4, including the sign the initial currents set).
	tr.I1, tr.I2 = i1, i2
	tr.W1, _ = wf.NodeVoltage("w1")
	tr.W2, _ = wf.NodeVoltage("w2")
	tr.Mismatch, tr.In = wf.Ledger.IntegrationMismatch, wf.Ledger.EnergyIn
	return tr, nil
}

// RunRK4 integrates the same loop and shaft equations directly (twin):
//
//	l11·İ1 + l12·İ2 = −R_int·I1 + M·ω1·(I2 + kself·I1) − v1,  v1 the output voltage,
//	J·ω̇1 = G1(t) − b·ω1 − M·I1·(I2 + kself·I1)
//
// and the same with 1 ↔ 2; held speeds replace the shaft equations. The
// output across each loop's terminals is R_L alone (v = R_L·I, lumped), or
// R_L ∥ C (C·v̇_c = I − v/R_L through the ESR), or R_L ∥ C ∥ L_p with L_p the
// primary of a transformer (winding R_p) whose secondary feeds R_ac:
//
//	L_p·i̇_p + M_T·i̇_s = v − R_p·i_p,   L_s2·i̇_s + M_T·i̇_p = −R_ac·i_s,
//	I = v/R_L + i_c + i_p,   i_c = (v − v_c)/ESR,   C·v̇_c = i_c.
func (m Dual) RunRK4(tEnd, h float64) DualTrace {
	M := m.Coupling()
	l11, l12, ks, R := m.loop()
	J := m.Inertia()
	det := l11*l11 - l12*l12
	speed := func(k int, t float64) float64 {
		ph := 0.0
		if k == 1 {
			ph = m.SpeedPhase2 * math.Pi / 180
		}
		return m.Speed[k] * (1 + m.SpeedAmp*math.Sin(2*math.Pi*m.SpeedFreq*t+ph))
	}
	torque := func(k int, t float64) float64 {
		return m.Torque[k] + m.TorqueAmp*math.Sin(2*math.Pi*m.TorqueFreq*t)
	}
	xfmr := m.LPrim > 0 && m.RLoad > 0
	tank := (m.CParallel > 0 || xfmr) && m.RLoad > 0
	sh, se := m.shared(), m.series()
	z, _ := m.TankScale() // a shared tank is the loops' tanks in parallel, a series one in series
	rl, cp, esr := m.RLoad*z, m.CParallel/z, m.ESR*z
	lp, ls, rp, rac := m.LPrim*z, m.LSec*z, m.RPrim*z, m.RAC*z
	mt := m.KT * math.Sqrt(lp*ls)
	dt := lp*ls - mt*mt
	// out returns the output voltage and the capacitor-branch current from the
	// current fed into the tank, the capacitor voltage and the primary current.
	out := func(i, vc, ip float64) (v, ic float64) {
		switch {
		case cp == 0:
			v = rl * (i - ip)
			return v, 0
		case esr > 0:
			v = (i - ip + vc/esr) / (1/rl + 1/esr)
		default:
			v = vc
		}
		return v, i - ip - v/rl
	}
	// tankD sets the tank's derivatives (states 4+j, 6+j, 8+j) and returns v.
	tankD := func(s, d *[10]float64, j int, i float64) float64 {
		v, ic := out(i, s[4+j], s[6+j])
		if cp > 0 {
			d[4+j] = ic / cp
		}
		if xfmr {
			a := v - rp*s[6+j]
			bs := -rac * s[8+j]
			d[6+j] = (ls*a - mt*bs) / dt
			d[8+j] = (lp*bs - mt*a) / dt
		}
		return v
	}
	// state: I1, I2, ω1, ω2, v_c1, v_c2, i_p1, i_p2, i_s1, i_s2
	f := func(t float64, s [10]float64) [10]float64 {
		i1, i2, w1, w2 := s[0], s[1], s[2], s[3]
		if m.Speed != nil {
			w1, w2 = speed(0, t), speed(1, t)
		}
		e1 := -R*i1 + M*w1*(i2+ks*i1)
		e2 := -R*i2 + M*w2*(i1+ks*i2)
		var d [10]float64
		if se {
			// one current I = I1 = −I2; both coils' fields carry +I (coil 2
			// wound to aid), disk 2's EMF counts backwards round the loop:
			// 2·(l11 + l12)·İ = M·(1+ks)·(ω1 − ω2)·I − 2·R_int·I − v
			v := tankD(&s, &d, 0, i1)
			d[0] = (M*(1+ks)*(w1-w2)*i1 - 2*(R-m.RLoad)*i1 - v) / (2 * (l11 + l12))
			d[1] = -d[0]
			if m.Speed == nil {
				d[2] = (torque(0, t) - m.Windage*w1 - M*(1+ks)*i1*i1) / J
				d[3] = (torque(1, t) - m.Windage*w2 + M*(1+ks)*i1*i1) / J
			}
			return d
		}
		if tank {
			// R above includes the load; take it out and use the output voltage.
			if sh {
				v := tankD(&s, &d, 0, i1+i2)
				e1 += m.RLoad*i1 - v
				e2 += m.RLoad*i2 - v
			} else {
				for j, e := range []*float64{&e1, &e2} {
					*e += m.RLoad*s[j] - tankD(&s, &d, j, s[j])
				}
			}
		}
		d[0] = (l11*e1 - l12*e2) / det
		d[1] = (l11*e2 - l12*e1) / det
		if m.Speed == nil {
			d[2] = (torque(0, t) - m.Windage*w1 - M*i1*(i2+ks*i1)) / J
			d[3] = (torque(1, t) - m.Windage*w2 - M*i2*(i1+ks*i2)) / J
		}
		return d
	}
	var s [10]float64
	s[0], s[1], s[2], s[3] = m.I0[0], m.I0[1], m.W0[0], m.W0[1]
	if se {
		s[1] = -m.I0[0]
	}
	switch {
	case se && xfmr:
		s[6] = m.I0[0]
		s[4] = rp * s[6]
	case se:
		s[4] = rl * m.I0[0]
	case sh && xfmr: // one tank: the sum of the loop currents through its primary
		s[6] = m.I0[0] + m.I0[1]
		s[4] = rp * s[6]
	case sh:
		s[4] = rl * (m.I0[0] + m.I0[1])
	case xfmr: // DC through the primary, the capacitor at its DC voltage
		s[6], s[7] = m.I0[0], m.I0[1]
		s[4], s[5] = m.RPrim*m.I0[0], m.RPrim*m.I0[1]
	default:
		s[4], s[5] = m.RLoad*m.I0[0], m.RLoad*m.I0[1]
	}
	var tr DualTrace
	n := int(math.Round(tEnd / h))
	for k := 0; k <= n; k++ {
		t := float64(k) * h
		w1, w2 := s[2], s[3]
		if m.Speed != nil {
			w1, w2 = speed(0, t), speed(1, t)
		}
		tr.T, tr.I1, tr.I2, tr.W1, tr.W2 = append(tr.T, t), append(tr.I1, s[0]), append(tr.I2, s[1]), append(tr.W1, w1), append(tr.W2, w2)
		if tank {
			i := s[0]
			if sh {
				i += s[1]
			}
			v, _ := out(i, s[4], s[6])
			tr.V1, tr.IP1, tr.IS1 = append(tr.V1, v), append(tr.IP1, s[6]), append(tr.IS1, s[8])
		}
		k1 := f(t, s)
		var a, bb, c [10]float64
		for j := range s {
			a[j] = s[j] + h/2*k1[j]
		}
		k2 := f(t+h/2, a)
		for j := range s {
			bb[j] = s[j] + h/2*k2[j]
		}
		k3 := f(t+h/2, bb)
		for j := range s {
			c[j] = s[j] + h*k3[j]
		}
		k4 := f(t+h, c)
		for j := range s {
			s[j] += h / 6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
		}
	}
	return tr
}

func reversals(t, x []float64, tSkip, floor float64) []float64 {
	var out []float64
	last := 0.0
	for k := range t {
		if t[k] < tSkip || math.Abs(x[k]) < floor {
			continue
		}
		s := math.Copysign(1, x[k])
		if last != 0 && s != last {
			out = append(out, t[k])
		}
		last = s
	}
	return out
}

// burst is the rim speed of a thin copper hoop at its tensile strength,
// √(σ/ρ) (declared σ = 220 MPa, annealed copper): a conservative limit.
func burst() float64 { return math.Sqrt(220e6 / material.Copper.Density) }

// growth measures a speed-held run's late growth rate of the field.
func growth(tr DualTrace) float64 {
	n := len(tr.T) - 1
	return math.Log(math.Hypot(tr.I1[n], tr.I2[n])/math.Hypot(tr.I1[n/2], tr.I2[n/2])) / (tr.T[n] - tr.T[n/2])
}

func runDual() error {
	m := DefaultDual()
	M := m.Coupling()
	L, rc := m.Coil()
	fmt.Println("Dual motor, cross-excited (the earth_machine construction in copper), side by side and stacked")
	fmt.Printf("  disks as above (r_o %.2f m, %.0f mm Cu, liquid-metal shaft contacts %.0f µΩ, declared); coils %.0f turns,\n",
		m.Disk.Ro, m.Disk.Thickness*1e3, m.Disk.RShaft*1e6, m.N)
	fmt.Printf("  r_c %.2f m, %.0f mm long, %.0f mm² Cu window: M = %.4g V·s/A (flux through the disk annulus), L %.3g µH, R %.3g µΩ;\n",
		m.Rc, m.Lc*1e3, m.Aw*1e6, M, L*1e6, rc*1e6)
	fmt.Printf("  shaft J %.4f kg·m²; copper hoop burst speed %.0f m/s (declared 220 MPa)\n", m.Inertia(), burst())

	fmt.Println("  side by side, a feedback coil of K·N turns on each unit: threshold ω* = R/((1+K)·M)")
	for _, k := range []float64{-0.5, 0, 0.5, 1} {
		q := m
		q.K = k
		w := q.Threshold()
		fmt.Printf("    K = %+.1f: loop R %.1f µΩ, ω* = %.0f rad/s (%.0f rpm), rim %.0f m/s\n", k, q.LoopR()*1e6, w, w*60/(2*math.Pi), w*m.Disk.Ro)
	}
	fmt.Println("    a feedback coil adds coupling and its own copper in proportion: ω* barely moves (opposing raises it).")

	s := m
	s.Stacked = true
	ws := s.Threshold()
	fmt.Printf("  stacked on one axis, %.0f mm between disk planes: stray coupling MFar = %.3g V·s/A (%.3f of M), coil mutual %.3g µH (k = %.3f)\n",
		s.Gap*1e3, s.MFar(), s.MFar()/M, s.MCoils()*1e6, s.MCoils()/L)
	fmt.Printf("    threshold ω* = R/(M + MFar) = %.0f rad/s (%.0f rpm), rim %.0f m/s: a feedback coil of K = %.3f with no copper\n",
		ws, ws*60/(2*math.Pi), ws*m.Disk.Ro, s.MFar()/M)

	// Speed-held twins of the thresholds.
	for _, c := range []struct {
		name string
		q    Dual
	}{{"side by side", m}, {"stacked", s}} {
		for _, x := range []float64{0.9, 1.1} {
			q := c.q
			w := x * q.Threshold()
			q.Speed = &[2]float64{w, w}
			q.I0 = [2]float64{100, 100} // the symmetric mode alone
			tr, err := q.Run(0.05, 1e-5)
			if err != nil {
				return err
			}
			l11, l12, ks, r := q.loop()
			fmt.Printf("    %s, both shafts held at %.1f ω*: field grows at %+.3f 1/s; exact (−R + (1+K)·M·ω)/(l11 + l12) = %+.3f 1/s\n",
				c.name, x, growth(tr), (-r+(1+ks)*M*w)/(l11+l12))
		}
	}

	fmt.Println("  dual motor, a torque on each shaft, from 100 A / 50 A at 900 / 600 rad/s, 0.5 s skipped:")
	fmt.Println("    layout        K      torque [N·m]  μ      reversals  mean interval [s]  RMS I1 [kA]  mean ω1 [rpm]  rim mean / peak [m/s]  motors [kW]  ledger/drive")
	for _, c := range []struct {
		name    string
		stacked bool
		k, tq   float64
		tEnd    float64
	}{
		{"side by side", false, 0, 20, 6}, {"side by side", false, 0.5, 20, 6},
		{"stacked", true, 0, 20, 6},
		{"side by side", false, 0, 1800, 12}, {"side by side", false, 0.5, 1800, 12},
		{"stacked", true, 0, 1800, 12},
	} {
		q := m
		q.K, q.Torque, q.Stacked = c.k, [2]float64{c.tq, c.tq}, c.stacked
		tr, err := q.Run(c.tEnd, 5e-4)
		if err != nil {
			return fmt.Errorf("%s K %g: %w", c.name, c.k, err)
		}
		rev := reversals(tr.T, tr.I1, 0.5, 50)
		mean := math.NaN()
		if len(rev) > 1 {
			mean = (rev[len(rev)-1] - rev[0]) / float64(len(rev)-1)
		}
		var s2, wsum, pin, wmax float64
		cnt := 0
		for j := range tr.T {
			if tr.T[j] < 0.5 {
				continue
			}
			s2 += tr.I1[j] * tr.I1[j]
			wsum += tr.W1[j]
			wmax = math.Max(wmax, math.Max(math.Abs(tr.W1[j]), math.Abs(tr.W2[j])))
			pin += c.tq * (tr.W1[j] + tr.W2[j])
			cnt++
		}
		rim := wsum / float64(cnt) * m.Disk.Ro
		flag := ""
		if rim > burst() {
			flag = " BURST"
		}
		fmt.Printf("    %-12s  %+.1f   %6.0f        %5.2f   %5d      %8.3f           %6.1f       %7.0f        %4.0f / %4.0f%-6s  %8.0f     %.0e\n",
			c.name, c.k, c.tq, q.Mu(), len(rev), mean, math.Sqrt(s2/float64(cnt))/1e3, wsum/float64(cnt)*60/(2*math.Pi), rim, wmax*m.Disk.Ro, flag, pin/float64(cnt)/1e3, tr.Mismatch/tr.In)
	}
	fmt.Println("  rim: running mean / peak (the peak is the start-up overshoot, before the field takes the load: a real")
	fmt.Println("  machine starts on its speed drives). BURST: the running rim passes the copper hoop limit.")
	fmt.Println("  at 20 N·m the copper machine has μ ≈ 10: it excites itself and settles to a steady field. Reversals,")
	fmt.Println("  the earth_machine behaviour, need μ near 1, about 90× the torque and megawatts at the shafts.")
	fmt.Println("  Every watt of field and heat is paid by the two motors (ledger); stacking and feedback coils move the")
	fmt.Println("  threshold and the dynamics, they add no input.")
	fmt.Println()
	return nil
}
