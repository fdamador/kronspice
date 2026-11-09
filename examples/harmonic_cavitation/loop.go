// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Part 10: the closed steam loop (user's layout, 2026-11-08). Water enters
// each drum at its outer top, is heated by the drums' cavitation and shaft
// work, the steam drives the Tesla turbine on the same shaft, is condensed
// at the outlet (the condenser heats the home) and pumped back to the drums.
// The motor supplies what the turbine does not.
//
// Two ways to raise the steam are compared:
//   A. boiling in the drums (the drums are the boiler at p_s);
//   B. the drums heat pressurized, subcooled water (p_d above its vapour
//      pressure), which is throttled into a flash tank at p_f; the flash
//      steam drives the turbine and the hot liquid goes back to the drums.
// Cavitation collapses only in subcooled water: the collapse energy is
// (p − p_v)·V, and steam that boils off and leaves the drum has not
// collapsed. So in A the latent heat (≈ 90 % of the heat) is added with no
// collapse, and the water's own p − p_v falls to zero as it nears boiling.

// SatPressure is water's saturation pressure [Pa] at T [K], IAPWS-IF97
// region 4, eq. 30 (Wagner et al., J. Eng. Gas Turbines Power 122, 150,
// 2000). Verification: 300 K → 3536.58941 Pa (tested).
func SatPressure(T float64) float64 {
	n := [...]float64{0, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
		0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4,
		0.40511340542057e6, -0.23855557567849, 0.65017534844798e3}
	th := T + n[9]/(T-n[10])
	a := th*th + n[1]*th + n[2]
	b := n[3]*th*th + n[4]*th + n[5]
	c := n[6]*th*th + n[7]*th + n[8]
	x := 2 * c / (-b + math.Sqrt(b*b-4*a*c))
	return x * x * x * x * 1e6
}

// Sat is one saturated state from the steam tables (IAPWS-95 values as
// tabulated in Çengel & Boles, Thermodynamics, Table A-5): p [Pa], T [°C],
// h [kJ/kg], s [kJ/kg/K]. To be checked against an IF97 implementation.
type Sat struct {
	P, T, Hf, Hfg, Hg, Sf, Sfg, Sg float64
}

var (
	sat20k  = Sat{20e3, 60.06, 251.42, 2357.5, 2608.9, 0.8320, 7.0752, 7.9072}
	sat101k = Sat{101.325e3, 99.97, 419.06, 2256.4, 2675.5, 1.3069, 6.0476, 7.3545}
	sat150k = Sat{150e3, 111.35, 467.13, 2226.0, 2693.1, 1.4337, 5.7894, 7.2231}
	sat200k = Sat{200e3, 120.21, 504.71, 2201.6, 2706.3, 1.5302, 5.5968, 7.1270}
)

// TurbineWork returns the work of saturated steam at s expanding to the
// condenser c [kJ/kg], at isentropic efficiency eta, and the exit enthalpy.
func TurbineWork(s, c Sat, eta float64) (w, hOut float64) {
	x := (s.Sg - c.Sf) / c.Sfg
	h2s := c.Hf + x*c.Hfg
	w = eta * (s.Hg - h2s)
	return w, s.Hg - w
}

// LoopResult is one cycle per kilogram of steam, scaled to a shaft power.
type LoopResult struct {
	Name                    string
	QIn, WTurbine, WPump    float64 // per kg of steam [kJ/kg]
	QCondenser              float64 // [kJ/kg], to the home
	EtaSensible, EtaLoop    float64 // collapse per joule of shaft heat: on the sensible heating, and over all heat
	FracSensible            float64 // share of the heat that is sensible (subcooled) heating
	SteamKgH                float64 // at the shaft power
	PShaft, PTurbine        float64 // [W]
	PMotorShaft, PElec      float64 // [W]
	PHome, PRoom, PCollapse float64 // [W]
	CavPerElecJ             float64
	RPM                     float64 // speed chosen for the heating path
	Residual                float64 // books: electricity − (home + room) [W]
	TwinResidual            float64 // kronspice shaft circuit ledger, relative
}

// pathEta averages η over a heating path from tIn to tOut [°C] at drum
// pressure pd, with the machine's geometry, choosing the one shaft speed
// that maximizes the average (the drum runs at one speed for the whole
// path). Equal heat per degree, so the average is over temperature.
func pathEta(p Params, m Machine, pd, tIn, tOut float64) (eta, rpm float64) {
	m.POut = pd
	const pts = 16
	for r := 500.0; r <= 15000; r += 100 {
		m.RPM = r
		sum := 0.0
		for k := 0; k < pts; k++ {
			t := tIn + (tOut-tIn)*(float64(k)+0.5)/pts
			q := p
			q.PVap = SatPressure(t + 273.15)
			if e, err := q.Evaluate(m); err == nil {
				sum += e.Eta
			}
		}
		if a := sum / pts; a > eta {
			eta, rpm = a, r
		}
	}
	return
}

// LoopA: the drums boil the condensate at ps; turbine to the condenser c.
func LoopA(p Params, m Machine, ps, c Sat, etaT, pShaft float64) LoopResult {
	r := LoopResult{Name: fmt.Sprintf("A: boil in the drums at %.2f bar", ps.P/1e5)}
	wp := 0.00102 * (ps.P - c.P) / 1e3 / 0.7 // pump: v·Δp / η_pump [kJ/kg]
	r.WPump = wp
	r.QIn = ps.Hg - c.Hf - wp*0.7 // shaft heat per kg of steam (the pump's useful work also enters)
	r.WTurbine, _ = TurbineWork(ps, c, etaT)
	qs := ps.Hf - c.Hf
	r.FracSensible = qs / r.QIn
	r.EtaSensible, r.RPM = pathEta(p, m, ps.P, c.T, ps.T)
	r.EtaLoop = r.EtaSensible * r.FracSensible // latent heat: steam that leaves has not collapsed
	return finish(r, c, pShaft)
}

// LoopB: the drums heat subcooled water at pd from the flash tank's liquid
// to tTop; it flashes at pf; the flash steam drives the turbine to c; the
// condensate returns to the flash tank.
func LoopB(p Params, m Machine, pd float64, top, pf, c Sat, etaT, pShaft float64) LoopResult {
	r := LoopResult{Name: fmt.Sprintf("B: subcooled drums at %.1f bar, flash at %.2f bar", pd/1e5, pf.P/1e5)}
	x := (top.Hf - pf.Hf) / pf.Hfg // flash fraction per kg circulated
	circ := 1 / x                  // kg circulated through the drums per kg of steam
	wp := (circ*0.00104*(pd-pf.P) + 0.00102*(pf.P-c.P)) / 1e3 / 0.7
	r.WPump = wp
	r.QIn = pf.Hg - c.Hf - wp*0.7 // per kg of steam, same end states as a boiler at pf
	r.WTurbine, _ = TurbineWork(pf, c, etaT)
	r.FracSensible = 1 // every joule goes into subcooled water
	// the drum water runs from the flash tank's liquid (pf.T) to top.T
	r.EtaSensible, r.RPM = pathEta(p, m, pd, pf.T, top.T)
	r.EtaLoop = r.EtaSensible
	return finish(r, c, pShaft)
}

// finish scales a cycle to a shaft power (the drums' shaft work, all of it
// heat in the water) and books it: the turbine gives back w_T per kg, the
// motor the rest; the condenser takes h_out − h_f to the home.
func finish(r LoopResult, c Sat, pShaft float64) LoopResult {
	const etaMotor = 0.90
	mdot := pShaft / (r.QIn * 1e3) // kg/s of steam
	r.SteamKgH = mdot * 3600
	r.PShaft = pShaft
	r.PTurbine = mdot * r.WTurbine * 1e3
	r.PMotorShaft = pShaft - r.PTurbine
	pPump := mdot * r.WPump * 1e3
	r.PElec = r.PMotorShaft/etaMotor + pPump
	r.PRoom = r.PMotorShaft/etaMotor - r.PMotorShaft + pPump*0.3 // motor losses; the pump's losses outside the water
	r.QCondenser = r.QIn + r.WPump*0.7 - r.WTurbine
	r.PHome = mdot * r.QCondenser * 1e3
	r.PCollapse = r.EtaLoop * pShaft
	r.CavPerElecJ = r.PCollapse / r.PElec
	r.Residual = r.PElec - (r.PHome + r.PRoom)
	return r
}

// LoopCircuit is the loop's shaft in kronspice (mobility analogy): the
// motor as the speed source ω, the turbine as a torque (current) source
// into the shaft node, the drums' loads as conductances (collapse, the rest
// of the drum heat). The operating point must show the motor giving
// P_shaft − P_turbine.
func LoopCircuit(r LoopResult, rpm float64) (motor, resid float64, err error) {
	w := omega(rpm)
	b := circuit.NewBuilder("steam loop shaft")
	b.VoltageSource("VMOTOR", "w", "0", w)
	b.CurrentSource("ITURBINE", "0", "w", r.PTurbine/w)
	if r.PCollapse > 0 {
		b.Resistor("RCOLLAPSE", "w", "0", w*w/r.PCollapse)
	}
	b.Resistor("RDRUMHEAT", "w", "0", w*w/(r.PShaft-r.PCollapse))
	c, err := b.Build()
	if err != nil {
		return 0, 0, err
	}
	res, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return 0, 0, err
	}
	src, _ := res.Element("VMOTOR")
	return -src.PowerAbsorbed, res.Ledger.RelativeResidual(), nil
}
