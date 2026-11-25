// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"strings"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

// Gas discharges are outside the engine's scope (CLAUDE.md §1, non-goals),
// so the gas enters the circuits the way that section allows: as a cited
// breakdown threshold (the Paschen voltage, gas.go) and as a declared
// time-varying element. A spark gap is a switch with memory; kronspice has
// no switch element yet, so each circuit is run in its two phases, before
// and after the gaps fire, with the state handed over through IC= values.

// Marx is an N-stage Marx generator whose spark gaps are filled with one
// gas. Classic topology (Marx 1923; Bluhm, Pulsed Power Systems, Springer
// 2006, §4.1): each stage capacitor C sits between a low node L_k and a
// high node H_k; charging resistors join H_k to H_{k+1} and L_k to L_{k+1};
// gap k joins H_k to L_{k+1}. Charged in parallel to Vc, then the gaps
// fire and the stages stand in series: about N·Vc at the output.
type Marx struct {
	Gas          Gas
	P, Gap       float64 // gap fill pressure [Pa] and gap length [m]
	Stages       int
	C, Rch, Rs   float64 // stage capacitance [F], charging resistors [Ω], supply resistor [Ω]
	Rarc, Lstray float64 // arc resistance per gap [Ω] (declared), loop inductance [H] (declared)
	Rload        float64 // load [Ω]
	Margin       float64 // charge to Margin·V_breakdown, then trigger the first gap
	ECE          string  // "" plain, or a .ece line for the hypothetical H48 terms
}

// DefaultMarx: four stages, 10 nF, 0.2 mm gaps at 500 psi, charged to
// 90 % of breakdown (a triggered Marx; self-breaking fires with the gap's
// shot-to-shot jitter), 0.5 Ω per arc and 1 µH loop (declared), 50 Ω load.
func DefaultMarx(g Gas) Marx {
	return Marx{Gas: g, P: 500 * psi, Gap: 0.2e-3, Stages: 4, C: 10e-9, Rch: 100e3, Rs: 1e6,
		Rarc: 0.5, Lstray: 1e-6, Rload: 50, Margin: 0.9}
}

// Vb is the gap's breakdown voltage [V]; Vc the charging voltage [V].
func (m Marx) Vb() float64 { return m.Gas.Breakdown(m.P, m.Gap) }

// Vc is the stage charging voltage.
func (m Marx) Vc() float64 { return m.Margin * m.Vb() }

// ChargeNetlist is the charging phase: the supply Vc through Rs into H1,
// gaps open (their leakage, a few pA at most, is left out); the supply
// steps on at t = 0, so the capacitors start empty.
func (m Marx) ChargeNetlist() string {
	var b strings.Builder
	fmt.Fprintf(&b, "Marx generator charging, %d stages, %s gaps\n", m.Stages, m.Gas.Formula)
	fmt.Fprintf(&b, "VS s 0 PULSE(0 %.6g 0 1u 1u 1 2)\nRS s h1 %.6g\n", m.Vc(), m.Rs)
	for k := 1; k <= m.Stages; k++ {
		fmt.Fprintf(&b, "C%d h%d %s %.6g\n", k, k, low(k), m.C)
		if k < m.Stages {
			fmt.Fprintf(&b, "RH%d h%d h%d %.6g\nRL%d %s %s %.6g\n", k, k, k+1, m.Rch, k, low(k), low(k+1), m.Rch)
		}
	}
	fmt.Fprintf(&b, ".tran %.6g %.6g\n.end\n", m.ChargeTime()/400, m.ChargeTime())
	return b.String()
}

// ChargeTime is 8 time constants of the slowest stage (the last one,
// through Rs and the ladder).
func (m Marx) ChargeTime() float64 {
	return 8 * (m.Rs + 2*float64(m.Stages-1)*m.Rch) * float64(m.Stages) * m.C
}

func low(k int) string {
	if k == 1 {
		return "0"
	}
	return fmt.Sprintf("l%d", k)
}

// FireNetlist is the discharge: every stage at Vc (IC=), each gap an arc
// of Rarc, the last high node into the load through the loop inductance.
func (m Marx) FireNetlist() string {
	var b strings.Builder
	fmt.Fprintf(&b, "Marx generator erected, %d stages, %s gaps\n", m.Stages, m.Gas.Formula)
	if m.ECE != "" {
		b.WriteString(m.ECE + "\n")
	}
	for k := 1; k <= m.Stages; k++ {
		fmt.Fprintf(&b, "C%d h%d %s %.6g IC=%.6g\n", k, k, low(k), m.C, m.Vc())
		if k < m.Stages {
			fmt.Fprintf(&b, "RH%d h%d h%d %.6g\nRL%d %s %s %.6g\n", k, k, k+1, m.Rch, k, low(k), low(k+1), m.Rch)
			fmt.Fprintf(&b, "RG%d h%d %s %.6g\n", k, k, low(k+1), m.Rarc)
		}
	}
	fmt.Fprintf(&b, "RG%d h%d x %.6g\nLS x out %.6g\nRLOAD out 0 %.6g\n", m.Stages, m.Stages, m.Rarc, m.Lstray, m.Rload)
	fmt.Fprintf(&b, ".tran 1n %.6g UIC\n.end\n", m.FireTime())
	return b.String()
}

// FireTime covers the pulse: 10 decay times of the erected bank in the load.
func (m Marx) FireTime() float64 {
	return 10 * (m.Rload + float64(m.Stages)*m.Rarc) * m.C / float64(m.Stages)
}

// MarxResult is one gas's Marx shot.
type MarxResult struct {
	Marx
	Vb, Vc         float64 // breakdown and charging voltage per stage [V]
	StageV         []float64
	ChargeIn       float64 // energy from the supply while charging [J]
	ChargeHeat     float64 // in the supply and charging resistors [J]
	Stored         float64 // ½·N·C·Vc² [J]
	Peak           float64 // output peak [V]
	Gain           float64 // Peak / Vc
	Rise           float64 // 10–90 % rise time of the output [s]
	ELoad, EGaps   float64 // energy into the load and the arcs [J]
	EModel         float64 // model-supplied (ECE terms) [J]
	EHeat          float64 // all heat while firing (load, arcs, charging resistors) [J]
	StoredStart    float64 // stored energy when the gaps fire, by the ledger [J]
	ChargeResidual float64 // worst relative point residuals of the ledger
	FireResidual   float64
	ChargeMismatch float64 // integration mismatch / energy in, charging
	Mismatch       float64 // integration mismatch / stored, firing
}

func run(text string, stop float64) (*output.Waveforms, error) {
	d, err := netlist.ParseString(text)
	if err != nil {
		return nil, err
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		return nil, err
	}
	var spec analysis.TranSpec
	for _, a := range d.Analyses {
		if a.Kind == "tran" {
			spec = analysis.TranSpec{Step: a.Step, Stop: a.Stop, UIC: a.UIC}
		}
	}
	if spec.Stop == 0 {
		spec.Stop = stop
	}
	spec.MaxStep = spec.Step
	ctx, cancel := context.WithTimeout(context.Background(), time.Minute)
	defer cancel()
	return analysis.Tran(c, spec, analysis.TranOptions{Options: analysis.Options{Context: ctx}})
}

// energyIn integrates an element's delivered power (−absorbed) [J].
func energy(w *output.Waveforms, name string, sign float64) float64 {
	for e, n := range w.ElementNames {
		if strings.EqualFold(n, name) {
			s := 0.0
			for k := 1; k < len(w.Time); k++ {
				s += 0.5 * (w.P[e][k] + w.P[e][k-1]) * (w.Time[k] - w.Time[k-1])
			}
			return sign * s
		}
	}
	return math.NaN()
}

// Run charges and fires the generator.
func (m Marx) Run() (MarxResult, error) {
	r := MarxResult{Marx: m, Vb: m.Vb(), Vc: m.Vc()}
	ch, err := run(m.ChargeNetlist(), 0)
	if err != nil {
		return r, fmt.Errorf("%s charge: %w", m.Gas.Formula, err)
	}
	for k := 1; k <= m.Stages; k++ {
		h, _ := ch.NodeVoltage(fmt.Sprintf("h%d", k))
		l, _ := ch.NodeVoltage(low(k))
		r.StageV = append(r.StageV, h[len(h)-1]-l[len(l)-1])
	}
	r.ChargeIn = ch.Ledger.EnergyIn
	r.ChargeHeat = ch.Ledger.EnergyDissipated
	r.ChargeResidual = ch.Ledger.MaxRelativeResidual
	r.ChargeMismatch = ch.Ledger.IntegrationMismatch / ch.Ledger.EnergyIn
	r.Stored = 0.5 * float64(m.Stages) * m.C * r.Vc * r.Vc
	fw, err := run(m.FireNetlist(), 0)
	if err != nil {
		return r, fmt.Errorf("%s fire: %w", m.Gas.Formula, err)
	}
	out, _ := fw.NodeVoltage("out")
	ipk := 0
	for k, v := range out {
		if v > r.Peak {
			r.Peak, ipk = v, k
		}
	}
	r.Gain = r.Peak / r.Vc
	t10, t90 := math.NaN(), math.NaN()
	for k := 0; k <= ipk; k++ {
		if math.IsNaN(t10) && out[k] >= 0.1*r.Peak {
			t10 = fw.Time[k]
		}
		if math.IsNaN(t90) && out[k] >= 0.9*r.Peak {
			t90 = fw.Time[k]
		}
	}
	r.Rise = t90 - t10
	r.ELoad = energy(fw, "RLOAD", 1)
	for k := 1; k <= m.Stages; k++ {
		r.EGaps += energy(fw, fmt.Sprintf("RG%d", k), 1)
	}
	r.EModel = fw.Ledger.EnergyModelSupplied
	r.EHeat = fw.Ledger.EnergyDissipated
	r.StoredStart = fw.Ledger.StoredStart
	r.FireResidual = fw.Ledger.MaxRelativeResidual
	r.Mismatch = fw.Ledger.IntegrationMismatch / fw.Ledger.StoredStart
	return r, nil
}

// Pump is a charge pump built on a gas-gap varactor (Fig. "pump" in the
// README): a DC source charges the varactor through D1 while its
// capacitance is high, the capacitance then falls with its charge held, so
// its voltage rises, and D2 hands the charge to the output capacitor.
// Ideal output = r·Vin for a capacitance ratio r (charge Q = r·C0·Vin,
// then V = Q/C0), less the diode drops and the load's draw.
//
// The varactor's capacitance changes because a plasma column, lit by a
// separate RF supply (not modelled), fills a fraction f of the plate gap
// and acts as a conductor: C_on/C_off = 1/(1 − f) (parallel plates with a
// conducting slab; the sheath is neglected, declared). The work that
// changes C is booked on the varactor's electric-field port (EXT=electric-field):
// it is what the plasma supply pays, and the ledger shows it.
type Pump struct {
	Gas         Gas
	P, Gap      float64 // varactor fill and plate gap [m]
	Vin         float64 // DC input [V]
	C0          float64 // varactor capacitance with no plasma [F]
	Fill        float64 // plasma fraction of the gap when lit (declared)
	F           float64 // pump frequency [Hz]
	Cout, Rload float64
	Cycles      int
}

// DefaultPump: 100 V in, 10 nF varactor with plasma filling 2/3 of the gap
// (ratio 3), 10 kHz, 100 nF output into 1 MΩ, 1 mm plates at 500 psi.
func DefaultPump(g Gas) Pump {
	return Pump{Gas: g, P: 500 * psi, Gap: 1e-3, Vin: 100, C0: 10e-9, Fill: 2.0 / 3, F: 10e3,
		Cout: 100e-9, Rload: 1e6, Cycles: 200}
}

// Ratio is C_on/C_off.
func (p Pump) Ratio() float64 { return 1 / (1 - p.Fill) }

// Netlist is the pump. The varactor's SCALE starts at the ratio (plasma
// on), falls to 1 for the second half of each period.
func (p Pump) Netlist() string {
	per := 1 / p.F
	edge := 0.1 * per // edges of 10 % of the period: the time step resolves dC/dt (ledger 6e-4; 2 % edges at 1/200 of the period left 13 %)
	var b strings.Builder
	fmt.Fprintf(&b, "gas-gap varactor charge pump, %s at %.0f psi, ratio %.3g\n", p.Gas.Formula, p.P/psi, p.Ratio())
	fmt.Fprintf(&b, "VIN in 0 DC %.6g\n", p.Vin)
	b.WriteString("D1 in v DSI\n")
	fmt.Fprintf(&b, "CV v 0 %.6g SCALE=PULSE(%.6g 1 %.6g %.6g %.6g %.6g %.6g) EXT=electric-field\n",
		p.C0, p.Ratio(), per/2-edge, edge, edge, per/2-edge, per)
	b.WriteString("D2 v out DSI\n")
	fmt.Fprintf(&b, "COUT out 0 %.6g\nRLOAD out 0 %.6g\n", p.Cout, p.Rload)
	b.WriteString(".model DSI D(IS=1e-14 N=1 RS=0.1)\n")
	fmt.Fprintf(&b, ".tran %.6g %.6g\n.end\n", per/2000, float64(p.Cycles)*per)
	return b.String()
}

// PumpResult is one pump run.
type PumpResult struct {
	Pump
	Ideal, Vout  float64 // r·Vin and the last cycle's mean output [V]
	Vhold        float64 // what the varactor's gap holds with the plasma off [V]
	PIn, PPlasma float64 // mean power from VIN and through the varactor's port, last 10 % [W]
	PLoad        float64
	Residual     float64
	Mismatch     float64 // integration mismatch / energy in
	NeutralRatio float64 // C ratio from compressing the neutral gas 1 atm → P instead
}

// Run simulates the pump.
func (p Pump) Run() (PumpResult, error) {
	r := PumpResult{Pump: p, Ideal: p.Ratio() * p.Vin, Vhold: p.Gas.Breakdown(p.P, p.Gap)}
	r.NeutralRatio = (1 + p.Gas.EpsMinus1*p.P/atm) / (1 + p.Gas.EpsMinus1)
	w, err := run(p.Netlist(), 0)
	if err != nil {
		return r, err
	}
	out, _ := w.NodeVoltage("out")
	tEnd := w.Time[len(w.Time)-1]
	t0 := 0.9 * tEnd
	var T, sv float64
	idx := map[string]int{}
	for e, n := range w.ElementNames {
		idx[strings.ToLower(n)] = e
	}
	var pin, pload, pcv float64
	for k := 1; k < len(w.Time); k++ {
		if w.Time[k] <= t0 {
			continue
		}
		dt := w.Time[k] - w.Time[k-1]
		T += dt
		sv += 0.5 * (out[k] + out[k-1]) * dt
		avg := func(e int) float64 { return 0.5 * (w.P[e][k] + w.P[e][k-1]) * dt }
		pin -= avg(idx["vin"])
		pload += avg(idx["rload"])
		pcv += avg(idx["cv"])
	}
	r.Vout = sv / T
	r.PIn, r.PLoad = pin/T, pload/T
	// The varactor absorbs electrical power pcv; over whole cycles its
	// stored energy repeats, so −pcv is what its port (the plasma supply)
	// delivers into the circuit.
	r.PPlasma = -pcv / T
	r.Residual = w.Ledger.MaxRelativeResidual
	r.Mismatch = w.Ledger.IntegrationMismatch / math.Max(w.Ledger.EnergyIn, 1e-30)
	return r, nil
}

// MarxECE returns the argon Marx with the hypothetical ECE2 terms (H48),
// using the declared values of examples/ece_rlc: ω₀ such that c·ω₀·10 mH =
// 0.25 Ω (8.339×10⁻⁸ m⁻¹) for the torsion run, κ = 0.05 for the scalar
// run. Standard prediction: the plain run.
func MarxECE(g Gas) []Marx {
	w0 := 0.25 / (circuit.SpeedOfLight * 10e-3)
	plain := DefaultMarx(g)
	tor, sca := plain, plain
	tor.ECE = fmt.Sprintf(".ece form=ece2 w0=%.10g", w0)
	sca.ECE = ".ece form=ece2 kappa=0.05"
	return []Marx{plain, tor, sca}
}
