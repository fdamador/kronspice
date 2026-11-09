// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 11: the heater sized for a home at −10 °F outside, 24/7 (user
// request, 2026-11-08). The closed steam loop of part 10, variant B
// (subcooled drums at 3 bar, flash, Tesla turbine, condenser to the home),
// with the shaft speed matched to the demand (user's choice): a variable-
// frequency drive sets the rpm so that the shaft work equals the heat the
// house needs at each outdoor temperature.
//
// Declared inputs (replace with the real house):
//   - indoor 70 °F, design outdoor −10 °F (ΔT 80 °F = 44.4 K);
//   - 3,000 ft² (278.7 m², user) at three envelope levels, rule-of-thumb design
//     loads 15 / 25 / 40 BTU/h per ft² at that ΔT (tight new / average /
//     older leaky). A room-by-room ACCA Manual J calculation replaces them;
//   - load ∝ (T_in − T_out) (UA constant; solar and internal gains ignored,
//     which is conservative);
//   - condenser U = 2,500 W/m²·K (condensing steam to water, typical
//     1,500–4,000);
//   - emitters: radiant floor (EN 1264 basic characteristic q = 8.92·
//     (θ_F − θ_i)^1.1 W/m², floor surface ≤ 29 °C) at 45/35 °C water, with
//     the condenser at 0.2 bar (60 °C); or panel radiators (EN 442, rated at
//     ΔT 50 K, exponent 1.3) at 70/50 °C water, with the condenser at 1.01 bar
//     (100 °C) and the flash at 1.5 bar;
//   - cold-climate heat pump for comparison: COP 1.5 at −10 °F (declared;
//     the NEEP cold-climate specification asks COP ≥ 1.75 at 5 °F).
//
// The condenser carries the whole load; the motor's and pump's losses
// (machine in the basement) warm the house besides, as margin.

import (
	"fmt"
	"math"
)

const (
	tInF       = 70.0
	tDesignF   = -10.0
	houseFt2   = 3000.0 // user, 2026-11-08
	ft2ToM2    = 0.09290304
	btuhToW    = 0.29307107
	uCondenser = 2500.0 // W/m²K, declared
	copHPCold  = 1.5    // declared, −10 °F
	floorMaxC  = 29.0   // EN 1264 floor surface limit, occupied zones
)

func fToC(f float64) float64 { return (f - 32) * 5 / 9 }

// Emitter is the house side of the condenser.
type Emitter struct {
	Name       string
	Cond       Sat     // condenser state
	Flash, Top Sat     // flash tank and drum outlet states
	TSup, TRet float64 // hydronic supply / return [°C]
	Floor      bool    // radiant floor (else radiators)
}

var emitters = []Emitter{
	{"radiant floor 45/35 °C, condenser 0.2 bar (60 °C)", sat20k, sat101k, sat200k, 45, 35, true},
	{"panel radiators 70/50 °C, condenser 1.01 bar (100 °C)", sat101k, sat150k, sat200k, 70, 50, false},
}

// Sizing is one house × emitter.
type Sizing struct {
	LoadW                float64
	Loop                 LoopResult // at the design shaft power
	Rows                 int
	RPM, RPMFull         float64 // matched speed; speed of the built rows at the loop's best speed
	DrumLen              float64 // band length per drum [m]
	EtaPath              float64 // collapse per shaft J along the heating path at the matched rpm
	CondArea, LMTD       float64
	FlowLmin             float64
	FloorM2, RadRatedW   float64
	KWhDay, KWhDayResist float64
	KWhDayHP             float64
}

// shaftAt is the drums' shaft power with rows rows at rpm (drum pressure pd).
func shaftAt(p Params, m Machine, rows int, rpm, pd float64) float64 {
	mm := withRows(m, rows)
	mm.RPM, mm.POut = rpm, pd
	e, err := p.Evaluate(mm)
	if err != nil {
		return math.NaN()
	}
	return e.PShaft
}

// withRows copies m with rows identical rows.
func withRows(m Machine, rows int) Machine {
	mm := m
	mm.Rows = rows
	mm.Rotor = make([][]Cavity, rows)
	for i := range mm.Rotor {
		mm.Rotor[i] = m.Rotor[0]
	}
	return mm
}

// etaAt averages η along the heating path tIn → tOut at a fixed rpm.
func etaAt(p Params, m Machine, rpm, pd, tIn, tOut float64) float64 {
	m.RPM, m.POut = rpm, pd
	const pts = 16
	var s float64
	for k := 0; k < pts; k++ {
		q := p
		q.PVap = SatPressure(tIn + (tOut-tIn)*(float64(k)+0.5)/pts + 273.15)
		if e, err := q.Evaluate(m); err == nil {
			s += e.Eta
		}
	}
	return s / pts
}

// matchRPM finds the rpm at which rows rows absorb pShaft (shaft power
// rises monotonically with speed).
func matchRPM(p Params, m Machine, rows int, pShaft, pd float64) float64 {
	lo, hi := 100.0, 20000.0
	for it := 0; it < 60; it++ {
		mid := (lo + hi) / 2
		if shaftAt(p, m, rows, mid, pd) < pShaft {
			lo = mid
		} else {
			hi = mid
		}
	}
	return (lo + hi) / 2
}

// SizeHouse sizes the heater for a design load with emitter em.
func SizeHouse(p Params, best Scored, load float64, em Emitter) Sizing {
	const pd = 3e5
	s := Sizing{LoadW: load}
	// shaft power per watt to the condenser (finish is linear in shaft power)
	unit := LoopB(p, best.M, pd, em.Top, em.Flash, em.Cond, 0.30, 1000)
	pShaft := 1000 * load / unit.PHome
	s.Loop = LoopB(p, best.M, pd, em.Top, em.Flash, em.Cond, 0.30, pShaft)
	// rows: as many as the loop's best speed needs, so the design day runs
	// at (just under) that speed and colder-than-design days have no reserve
	// but milder days slow down.
	perRow := shaftAt(p, best.M, 1, unit.RPM, pd)
	s.Rows = int(math.Ceil(pShaft / perRow))
	s.RPMFull = unit.RPM
	s.RPM = matchRPM(p, best.M, s.Rows, pShaft, pd)
	s.DrumLen = float64(s.Rows) * best.M.RowLen
	s.EtaPath = etaAt(p, withRows(best.M, s.Rows), s.RPM, pd, em.Flash.T, em.Top.T)
	// condenser: steam at Cond.T to hydronic water TRet → TSup
	d1, d2 := em.Cond.T-em.TRet, em.Cond.T-em.TSup
	s.LMTD = (d1 - d2) / math.Log(d1/d2)
	s.CondArea = load / (uCondenser * s.LMTD)
	s.FlowLmin = load / (4186 * (em.TSup - em.TRet)) * 60 // ≈ kg/s → L/min
	tIn := fToC(tInF)
	if em.Floor {
		s.FloorM2 = load / (8.92 * math.Pow(floorMaxC-tIn, 1.1))
	} else {
		s.RadRatedW = load / math.Pow(((em.TSup+em.TRet)/2-tIn)/50, 1.3)
	}
	s.KWhDay = s.Loop.PElec * 24 / 1e3
	s.KWhDayResist = load * 24 / 1e3
	s.KWhDayHP = load / copHPCold * 24 / 1e3
	return s
}

func part11(p Params, best Scored) {
	dT := fToC(tInF) - fToC(tDesignF)
	fmt.Printf("11. Home heating at %.0f °F outside, 24/7 (indoor %.0f °F, ΔT %.1f K); loop B, drum speed matched to the demand by a VFD\n",
		tDesignF, tInF, dT)
	fmt.Printf("   declared: %.0f ft² (%.0f m²), loads 15/25/40 BTU/h·ft² (replace with Manual J), condenser U %.0f W/m²K,\n",
		houseFt2, houseFt2*ft2ToM2, uCondenser)
	fmt.Printf("   heat-pump COP %.1f at −10 °F; heater drums of part 4 (R %.1f mm, Ø %.1f mm), rows added to carry the load\n",
		copHPCold, best.M.R*1e3, best.M.Rotor[0][0].W*1e3)
	houses := []struct {
		name string
		bpf  float64
	}{{"tight new", 15}, {"average", 25}, {"older, leaky", 40}}
	for _, em := range emitters {
		fmt.Println("   " + em.Name)
		fmt.Println("     house          load kW  BTU/h   shaft kW  rows  drum band  rpm (best)  η path  steam kg/h  electric kW  kWh/day  resist. kWh/day  HP kWh/day  condenser m²  flow L/min  emitters")
		for _, h := range houses {
			load := h.bpf * houseFt2 * btuhToW
			s := SizeHouse(p, best, load, em)
			emit := fmt.Sprintf("radiators rated %.1f kW at ΔT50", s.RadRatedW/1e3)
			if em.Floor {
				emit = fmt.Sprintf("floor %.0f m² at 29 °C", s.FloorM2)
				if s.FloorM2 > houseFt2*ft2ToM2 {
					emit += " (more than the house has)"
				}
			}
			fmt.Printf("     %-13s %7.2f  %6.0f  %8.2f  %4d  %6.0f mm  %4.0f (%4.0f)  %6.3f  %10.1f  %11.2f  %7.0f  %15.0f  %10.0f  %12.2f  %10.1f  %s\n",
				h.name, load/1e3, load/btuhToW, s.Loop.PShaft/1e3, s.Rows, s.DrumLen*1e3, s.RPM, s.RPMFull, s.EtaPath,
				s.Loop.SteamKgH, s.Loop.PElec/1e3, s.KWhDay, s.KWhDayResist, s.KWhDayHP, s.CondArea, s.FlowLmin, emit)
			fmt.Printf("       books: electricity %.0f W = home (condenser) %.0f W + basement (motor, pump) %.0f W; residual %.1e W\n",
				s.Loop.PElec, s.Loop.PHome, s.Loop.PRoom, s.Loop.Residual)
		}
	}
	// the season: speed follows the load
	em := emitters[1]
	load := 25 * houseFt2 * btuhToW
	s := SizeHouse(p, best, load, em)
	unit := LoopB(p, best.M, 3e5, em.Top, em.Flash, em.Cond, 0.30, 1000)
	fmt.Printf("   through the winter (average house, radiators, %d rows): load ∝ (70 °F − outside)\n", s.Rows)
	fmt.Println("     outside °F  load kW  shaft kW   rpm   U m/s   σ      η path   cavitating")
	m := withRows(best.M, s.Rows)
	for _, tf := range []float64{-20, -10, 0, 10, 20, 30, 40, 50, 60} {
		l := load * (tInF - tf) / (tInF - tDesignF)
		ps := 1000 * l / unit.PHome
		rpm := matchRPM(p, best.M, s.Rows, ps, 3e5)
		mm := m
		mm.RPM, mm.POut = rpm, 3e5
		e, err := p.Evaluate(mm)
		eta := etaAt(p, m, rpm, 3e5, em.Flash.T, em.Top.T)
		cav := "yes"
		if err != nil || eta < 0.01 {
			cav = "no: below inception speed, cycle on/off at the best speed instead"
		}
		if rpm > s.RPMFull*1.0001 {
			cav += " (beyond the design: over speed)"
		}
		fmt.Printf("     %10.0f  %7.2f  %8.2f  %5.0f  %5.1f  %5.2f  %7.3f   %s\n", tf, l/1e3, ps/1e3, rpm, e.U, e.Sigma, eta, cav)
	}
	// hole size, drum radius and the turbine's speed
	fmt.Println("   dimple size and drum radius against the best speed (one row, loop B heating path 111 → 120 °C at 3 bar):")
	for _, R := range []float64{best.M.R, 0.040} {
		for _, w := range []float64{0.003, 0.006, 0.010, 0.025} {
			mm := DimpleMachine(R, 1, w, 0.5, 3.33, 1.20, best.M.Gap, 2000, 3e5)
			eta, rpm := pathEta(p, mm, 3e5, sat150k.T, sat200k.T)
			mm.RPM, mm.POut = rpm, 3e5
			e, _ := p.Evaluate(mm)
			fmt.Printf("     R %4.1f mm, Ø %4.1f mm: best %5.0f rpm, rim %4.1f m/s, η %.3f, shaft per row %5.0f W\n", R*1e3, w*1e3, rpm, e.U, eta, e.PShaft)
		}
	}
	fmt.Println("     the best speed is set by the rim speed (σ = (p − p_v)/½ρU²), not by the hole size; a smaller drum turns faster")
	fmt.Println("   the Tesla turbine's speed: tip ≈ half the isentropic steam-jet speed (Rice 1965: best near the jet's tangential speed)")
	for _, em := range emitters {
		dh, _ := TurbineWork(em.Flash, em.Cond, 1)
		v := math.Sqrt(2 * dh * 1e3)
		fmt.Printf("     %.2f → %.2f bar: Δh_s %3.0f kJ/kg, jet %3.0f m/s; tip v/2 at R %.1f mm needs %5.0f rpm (drums best ≈ %4.0f)\n",
			em.Flash.P/1e5, em.Cond.P/1e5, dh, v, best.M.R*1e3, v/2/best.M.R*60/(2*math.Pi), SizeHouse(p, best, 25*houseFt2*btuhToW, em).RPMFull)
	}
	fmt.Println("     so the turbine on the drums' shaft runs at ~1/10 of its speed; on its own fast shaft (generator back to the motor),")
	fmt.Println("     the most it can return is its isentropic work (radiators, average house):")
	for _, em := range emitters {
		for _, et := range []float64{0.30, 0.70, 1.0} {
			load := 25 * houseFt2 * btuhToW
			u := LoopB(p, best.M, 3e5, em.Top, em.Flash, em.Cond, et, 1000)
			r := LoopB(p, best.M, 3e5, em.Top, em.Flash, em.Cond, et, 1000*load/u.PHome)
			fmt.Printf("     %-55s turbine η %.2f: turbine %5.0f W, electric %6.0f W for %6.0f W home heat (home/electric %.3f)\n",
				em.Name, et, r.PTurbine, r.PElec, r.PHome, r.PHome/r.PElec)
		}
	}
	fmt.Println("   all heat ≤ the electricity (COP 1 with the machine indoors, as a resistance heater); a cold-climate heat pump")
	fmt.Println("   uses about 1/1.5 of it at −10 °F and less in milder weather. The cavitation does not add heat; it changes where it forms.")
}
