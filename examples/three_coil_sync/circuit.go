// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// Declared electronics (not evolved).
const (
	vcc       = 24.0   // Royer supply [V]
	vGate     = 12.0   // gate-drive rail [V]
	rSupply   = 0.03   // supply + choke ESR [Ω]
	lChoke    = 100e-6 // Royer feed choke [H]
	rGate     = 470.0  // gate pull-up resistors [Ω]
	rGateSkew = 1.02   // RG2 = 1.02·RG1 (part tolerance): which side conducts first is fixed
	cGate     = 3e-9   // MOSFET gate capacitance (IRFP260N C_iss ≈ 4 nF) [F]
	cOss      = 1e-9   // MOSFET output capacitance C_oss, drain to source (IRFP260N ≈ 0.5–1 nF) [F]
	esrC2     = 5e-3   // tank capacitor ESR [Ω]
	vgHV      = 48.0   // HV generator amplitude [V] (a 48 V full bridge or audio amplifier)
	rgHV      = 1.0    // HV generator source resistance [Ω]

	tStop    = 2.5e-3 // simulated time [s]
	tMeasure = 2.0e-3 // steady-state window starts here [s]
	maxStep  = 0.2e-6
	softUp   = 20e-6 // supply and gate-rail ramp at switch-on [s]
)

// RoyerFreq0 is the unloaded tank frequency 1/(2π√((L2a + L2b + 2M)·C2)).
func (d Design) RoyerFreq0() float64 {
	e := d.Electrical(50e3)
	lt := e.L2a + e.L2b + 2*e.M2ab
	return 1 / (2 * math.Pi * math.Sqrt(lt*d.C2))
}

// Build turns a design into the kronspice circuit. Nodes:
//
//	coil 1 (HV):    hvg —Rg— hvp —R1p— p1 —L1p— 0;  hv —L1s— s1 —R1s— 0, C1self hv–0
//	coil 3 (out):   hv —L3— a3 —R3— b3 —C3— out —RL— 0
//	coil 2 (Royer): vcc —Rs— v1 —Lchoke— ct; c1 —L2a— a2 —R2a— ct —L2b— b2 —R2b— c2;
//	                C2 c1–c2; M1 c1/g1, M2 c2/g2; gates from 12 V via 470 Ω,
//	                pulled down through D1 (g1→c2) and D2 (g2→c1)
//
// Couplings: L1p–L1s (coil 1 is its own air-core transformer, set apart),
// L2a–L2b (bifilar strands), L2a–L3 and L2b–L3 (coil 2 wound over coil 3).
func (d Design) Build() (*circuit.Circuit, Electrical, error) {
	f := d.FHV
	if f <= 0 {
		f = 50e3
	}
	e := d.Electrical(f)
	b := circuit.NewBuilder("three coils: HV coil 1 wired to resonant coil 3, Royer coil 2 wound over coil 3")

	// Coil 1: HV generator into an air-core step-up transformer.
	amp := vgHV
	if d.VgAmp > 0 {
		amp = d.VgAmp
	}
	if d.HVOff {
		amp = 0
	}
	if d.Sync {
		// Slaved HV drive: the tank voltage c1 − c2, buffered to ground (Es),
		// through a first-order all-pass phase shifter (Rph, Cph and the
		// half divider Rh1, Rh2), drives the HV amplifier Eg with gain G:
		// v(hvg) = G·(2·v(lag) − v(s)) = G·v(s)·(1 − jωRC)/(1 + jωRC).
		// Eg stands for a class-D power amplifier clocked by a pickup on
		// coil 2; its power is the amplifier's supply, booked as
		// model-supplied and counted as input (PInHV).
		gain := d.SyncGain * float64(d.SyncSign)
		if d.HVOff {
			gain = 0
		}
		b.VCVS("Es", "s", "0", "c1", "c2", 1)
		b.Resistor("Rsb", "s", "0", 1e6)
		b.Resistor("Rph", "s", "lag", d.syncR())
		b.Capacitor("Cph", "lag", "0", cPh)
		b.Resistor("Rh1", "s", "half", 1e5)
		b.Resistor("Rh2", "half", "0", 1e5)
		b.VCVS("Eg", "hvg", "0", "lag", "half", 2*gain)
	} else {
		b.VoltageWave("Vg", "hvg", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, amp, d.FHV, 0, 0, d.PhiHV}})
	}
	b.Resistor("Rg", "hvg", "hvp", rgHV)
	b.Resistor("R1p", "hvp", "p1", e.R1p)
	b.Inductor("L1p", "p1", "0", e.L1p)
	b.Inductor("L1s", "hv", "s1", e.L1s)
	b.Resistor("R1s", "s1", "0", e.R1s)
	b.Capacitor("C1self", "hv", "0", e.C1self)
	b.Mutual("K1", "L1p", "L1s", e.K1)

	// Coil 3: the resonant output loop, wired to coil 1's HV terminal.
	b.Inductor("L3", "hv", "a3", e.L3)
	b.Resistor("R3", "a3", "b3", e.R3)
	b.Capacitor("C3", "b3", "out", d.C3)
	if d.Regen {
		regenLoad(b, d)
	} else {
		b.Resistor("RL", "out", "0", d.RL)
	}

	// Coil 2: resonant (ZVS, Mazzilli) Royer oscillator, tank coil over coil 3.
	supply := vcc
	gate := vGate
	if d.RoyerOff {
		supply, gate = 0, 0
	}
	b.VoltageWave("Vcc", "vcc", "0", circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, supply, 0, softUp, softUp, 1, 2}})
	b.VoltageWave("Vgate", "v12", "0", circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, gate, 0, softUp, softUp, 1, 2}})
	b.Resistor("Rs", "vcc", "v1", rSupply)
	b.Inductor("Lchoke", "v1", "ct", lChoke)
	b.Inductor("L2a", "c1", "a2", e.L2a)
	b.Resistor("R2a", "a2", "ct", e.R2)
	b.Inductor("L2b", "ct", "b2", e.L2b)
	b.Resistor("R2b", "b2", "c2", e.R2)
	b.Resistor("RC2", "c1", "t2", esrC2)
	b.Capacitor("C2", "t2", "c2", d.C2)
	b.Model("NPWR", "nmos", map[string]float64{"level": 1, "vto": 4, "kp": 30, "lambda": 0.001, "rd": 0.02, "rs": 0.005})
	b.Model("DFAST", "d", map[string]float64{"is": 1e-9, "n": 1.5, "rs": 0.1})
	b.MOSFET("M1", "c1", "g1", "0", "0", "NPWR", nil)
	b.MOSFET("M2", "c2", "g2", "0", "0", "NPWR", nil)
	b.Resistor("RG1", "v12", "g1", rGate)
	b.Resistor("RG2", "v12", "g2", rGate*rGateSkew)
	b.Capacitor("CG1", "g1", "0", cGate)
	b.Capacitor("CG2", "g2", "0", cGate)
	b.Capacitor("COSS1", "c1", "0", cOss)
	b.Capacitor("COSS2", "c2", "0", cOss)
	b.Diode("D1", "g1", "c2", "DFAST")
	b.Diode("D2", "g2", "c1", "DFAST")
	b.Mutual("K2ab", "L2a", "L2b", e.K2ab)
	b.Mutual("K2a3", "L2a", "L3", e.K2a3)
	b.Mutual("K2b3", "L2b", "L3", e.K2b3)

	c, err := b.Build()
	return c, e, err
}

// cPh is the phase shifter's capacitor [F].
const cPh = 1e-9

// syncR is the all-pass resistor giving a phase lag of SyncPhase degrees at
// the unloaded tank frequency: lag = 2·atan(ωRC).
func (d Design) syncR() float64 {
	w := 2 * math.Pi * d.RoyerFreq0()
	ph := math.Max(0.5, math.Min(179.5, d.SyncPhase))
	return math.Tan(ph*math.Pi/360) / (w * cPh)
}

// Regenerative load (declared parts): the load resistor is replaced by a
// step-down transformer (ferrite, k = 0.995, primary L_p, ratio n:1,
// winding resistances), a full bridge of low-drop rectifiers (Schottky or
// synchronous MOSFETs, about 0.35 V at 20 A, 2 nF junction capacitance), an
// RC snubber across the secondary and a 0 V ammeter into the
// 24 V supply rail. The energy that would have heated the load returns to
// the supply, which then delivers only the circuit's losses.
const (
	regenLp  = 2e-3  // transformer primary (magnetizing) inductance [H]
	regenK   = 0.995 // transformer coupling
	regenRp  = 0.05  // primary winding resistance [Ω]
	regenRs  = 0.005 // secondary winding resistance [Ω]
	regenRsn = 2.0   // snubber resistor [Ω]
	regenCsn = 47e-9 // snubber capacitor [F]
)

func regenLoad(b *circuit.Builder, d Design) {
	n := d.RegenN
	b.Inductor("Lr1", "out", "r1", regenLp)
	b.Resistor("Rr1", "r1", "0", regenRp)
	b.Inductor("Lr2", "x1", "r2", regenLp/(n*n))
	b.Resistor("Rr2", "r2", "x2", regenRs)
	b.Mutual("Kr", "Lr1", "Lr2", regenK)
	b.Model("DSYN", "d", map[string]float64{"is": 1e-5, "n": 1.05, "rs": 0.005, "cjo": 2e-9})
	b.Diode("DB1", "x1", "rb", "DSYN")
	b.Diode("DB2", "x2", "rb", "DSYN")
	b.Diode("DB3", "0", "x1", "DSYN")
	b.Diode("DB4", "0", "x2", "DSYN")
	// RC snubber across the secondary (damps the leakage-inductance ring
	// at each diode commutation; its loss is counted).
	b.Resistor("Rsn", "x1", "sn", regenRsn)
	b.Capacitor("Csn", "sn", "x2", regenCsn)
	// The secondary floats between the bridge diodes; a 1 MΩ leak to
	// ground gives its nodes a DC path.
	b.Resistor("Rleak", "x2", "0", 1e6)
	b.VoltageSource("Vret", "rb", "vcc", 0)
}
