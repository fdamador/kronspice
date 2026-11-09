// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Stage 6: the energy books from the wall socket to the radiator.
//
// Two drives for the same shaft power P_s:
//   - direct: a motor on the shaft (η_motor declared 0.90);
//   - turbine: a pump (η_pump declared 0.70, motor included) drives water
//     through the Tesla turbine (η_turbine declared 0.30; reported Tesla
//     turbine efficiencies in water are low — to be verified, QUESTIONS.md).
// Every watt the pump and turbine lose in the water also heats the water;
// the motor losses heat the room. Heat to the water is booked; cavitation
// per electrical joule is η·P_s/P_e. Steam from 15 °C water at 1 atm:
// c_p 4186 J/kg/K × 85 K + h_fg 2.257 MJ/kg (IAPWS-IF97).

// Drive is one way of turning the shaft.
type Drive struct {
	Name                   string
	PElec, PShaft          float64 // [W]
	PWaterHeat, PRoomHeat  float64 // [W]
	PCollapse, CavPerJoule float64
	SteamKgH, WaterLmin45K float64
	HeatPumpW              float64 // heat a COP 3.5 heat pump gives on the same electricity
}

// Drives books both drives for an evaluated machine.
func Drives(e Eval) []Drive {
	const etaMotor, etaPump, etaTurb = 0.90, 0.70, 0.30
	mk := func(name string, pe, water, room float64) Drive {
		return Drive{Name: name, PElec: pe, PShaft: e.PShaft, PWaterHeat: water, PRoomHeat: room,
			PCollapse: e.PCollapse, CavPerJoule: e.PCollapse / pe,
			SteamKgH: water / (4186*85 + 2.257e6) * 3600, WaterLmin45K: water / (4186 * 45) * 60,
			HeatPumpW: 3.5 * pe}
	}
	pe1 := e.PShaft / etaMotor
	direct := mk("motor on the shaft", pe1, e.PShaft, pe1-e.PShaft)
	ph := e.PShaft / etaTurb // hydraulic power into the turbine
	pe2 := ph / etaPump
	// pump hydraulic losses heat the water; its motor's ~10 % heats the room
	room := 0.10 * pe2
	turb := mk("pump + Tesla turbine", pe2, pe2-room, room)
	return []Drive{direct, turb}
}
