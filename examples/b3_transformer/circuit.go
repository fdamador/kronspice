// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Rig is the electrical set-up round the core.
type Rig struct {
	Freq    float64 // drive frequency [Hz]
	Vpk     float64 // drive amplitude [V]
	Rsrc    float64 // source resistance [Ω]
	AirLeak float64 // declared winding air leakage, fraction of each self inductance not linked to the core
	CuDia   float64 // copper wire diameter [m]
	Load    float64 // fixed load [Ω]; 0 = the output winding's reactance at Freq
}

// DefaultRig: 1 kHz, 20 V peak (0.25 T in the core), 0.5 Ω source, 1 % winding air leakage,
// Ø0.5 mm copper.
func DefaultRig() Rig { return Rig{Freq: 1e3, Vpk: 20, Rsrc: 0.5, AirLeak: 0.01, CuDia: 0.5e-3} }

const rhoCu = 1.72e-8

// WindingResistance is the copper resistance of a winding wrapped round
// the wall: each turn is a rectangle round the enclosed passes (wire
// diameter plus 0.5 mm clearance each side, height (n−1)·p + the same).
// The axis wire is one straight run up the bore closed outside.
func (s Spring) WindingResistance(w Winding, cu float64) float64 {
	area := math.Pi * cu * cu / 4
	if w.Name == "axis" {
		return rhoCu * 2 * (float64(s.N)*s.P + 2*s.R) / area
	}
	side := 2*s.A + 1e-3
	h := float64(w.Passes-1)*s.P + side
	return rhoCu * float64(w.Turns) * 2 * (side + h) / area
}

// CircuitRun is one kronspice run of the transformer.
type CircuitRun struct {
	RatioAC   float64 // |v_out/v_in| open circuit, AC analysis
	RatioAxis float64 // |v_axis/v_in| open circuit
	Rload     float64
	VinRMS    float64
	VoutRMS   float64 // loaded
	Pin, Pout float64 // mean over the last 5 cycles [W]
	EIn, EDis float64 // transient ledger [J]
	EOut      float64 // energy returned to the source (reactive part of each cycle) [J]
	ELoad     float64 // energy into the load over the run [J]
	ESrcR     float64 // heat in the source resistance [J]
	EStore    float64
	Mismatch  float64
	Bpk       float64 // peak core flux density at the input winding's passes [T]
	K         [3]float64
}

// build makes the circuit: source → Rsrc → input winding; output winding
// into rload; the axis wire into 1 MΩ. Windings carry the network
// inductances, the declared air leakage and their copper resistance.
func (s Spring) build(r Rig, in, out, axis Winding, rload float64, ac bool) (*circuit.Circuit, [3]float64, [][]float64, error) {
	ws := []Winding{in, out, axis}
	L, _, err := s.Inductance(ws)
	var ks [3]float64
	if err != nil {
		return nil, ks, nil, err
	}
	for i := range L {
		L[i][i] *= 1 + r.AirLeak
	}
	b := circuit.NewBuilder("B(3) transformer: N-turn spring core")
	src := circuit.Element{Name: "VSRC", Kind: circuit.VoltageSource}
	if ac {
		src.ACMag = 1
	} else {
		src.Wave = circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, r.Vpk, r.Freq}}
	}
	b.Add(src, "s", "0")
	b.Resistor("RSRC", "s", "a", r.Rsrc)
	names := []string{"LIN", "LOUT", "LAX"}
	top := []string{"a1", "o1", "x1"}
	bot := []string{"a", "o", "x"}
	for i, w := range ws {
		b.Resistor("RW"+names[i][1:], bot[i], top[i], s.WindingResistance(w, r.CuDia))
		b.Inductor(names[i], top[i], "0", L[i][i])
	}
	b.Resistor("RLOAD", "o", "0", rload)
	b.Resistor("RAX", "x", "0", 1e6)
	pairs := [][2]int{{0, 1}, {0, 2}, {1, 2}}
	for q, p := range pairs {
		k := L[p[0]][p[1]] / math.Sqrt(L[p[0]][p[0]]*L[p[1]][p[1]])
		if math.Abs(k) >= 1 {
			return nil, ks, nil, fmt.Errorf("coupling %s–%s = %.9f", names[p[0]], names[p[1]], k)
		}
		ks[q] = k
		b.Mutual(fmt.Sprintf("K%d", q), names[p[0]], names[p[1]], k)
	}
	c, err := b.BuildTransient()
	return c, ks, L, err
}

// RunCircuit drives `in` and loads `out`, in kronspice: an AC run open
// circuit for the ratios, then a transient with the load for the powers and
// the ledger.
func (s Spring) RunCircuit(r Rig, in, out Winding) (CircuitRun, error) {
	axis := Winding{Name: "axis", Azimuth: -math.Pi / 2, First: 0, Passes: s.N, Turns: 1}
	var cr CircuitRun
	c, ks, L, err := s.build(r, in, out, axis, 1e9, true)
	if err != nil {
		return cr, err
	}
	cr.K = ks
	res, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: r.Freq, FStop: r.Freq}, analysis.Options{})
	if err != nil {
		return cr, err
	}
	va, _ := res.NodeVoltage("a")
	vo, _ := res.NodeVoltage("o")
	vx, _ := res.NodeVoltage("x")
	cr.RatioAC = cmplx.Abs(vo[0]) / cmplx.Abs(va[0])
	cr.RatioAxis = cmplx.Abs(vx[0]) / cmplx.Abs(va[0])

	// Load: the output winding's reactance at the drive frequency.
	cr.Rload = 2 * math.Pi * r.Freq * L[1][1]
	if r.Load > 0 {
		cr.Rload = r.Load
	}
	c, _, _, err = s.build(r, in, out, axis, cr.Rload, false)
	if err != nil {
		return cr, err
	}
	per := 1 / r.Freq
	stop := 20 * per
	w, err := analysis.Tran(c, analysis.TranSpec{Step: per / 200, Stop: stop, MaxStep: per / 200}, analysis.TranOptions{})
	if err != nil {
		return cr, err
	}
	vs, _ := w.NodeVoltage("s")
	vin, _ := w.NodeVoltage("a")
	vout, _ := w.NodeVoltage("o")
	is, _ := w.ElementCurrent("VSRC")
	if is == nil {
		is, _ = w.ElementCurrent("V1")
	}
	t0 := stop - 5*per
	var sin2, sout2, pin, pout, tt float64
	for i := 1; i < len(w.Time); i++ {
		if w.Time[i] <= t0 {
			continue
		}
		dt := w.Time[i] - w.Time[i-1]
		tt += dt
		sin2 += vin[i] * vin[i] * dt
		sout2 += vout[i] * vout[i] * dt
		pout += vout[i] * vout[i] / cr.Rload * dt
		if is != nil {
			pin += -vs[i] * is[i] * dt
		}
	}
	cr.VinRMS, cr.VoutRMS = math.Sqrt(sin2/tt), math.Sqrt(sout2/tt)
	cr.Pin, cr.Pout = pin/tt, pout/tt
	for i := 1; i < len(w.Time); i++ { // whole-run split of the dissipation
		dt := w.Time[i] - w.Time[i-1]
		cr.ELoad += vout[i] * vout[i] / cr.Rload * dt
		cr.ESrcR += (vs[i] - vin[i]) * (vs[i] - vin[i]) / r.Rsrc * dt
	}
	cr.EIn, cr.EOut, cr.EDis = w.Ledger.EnergyIn, w.Ledger.EnergyOut, w.Ledger.EnergyDissipated
	cr.EStore = w.Ledger.StoredEnd - w.Ledger.StoredStart
	cr.Mismatch = w.Ledger.IntegrationMismatch
	// Peak flux density: dΦ/dt per pass = v_in/(N_in·n_in), Φ = that / ω.
	cr.Bpk = math.Sqrt2 * cr.VinRMS / (2 * math.Pi * r.Freq * float64(in.Turns*in.Passes)) / s.Area()
	return cr, nil
}
