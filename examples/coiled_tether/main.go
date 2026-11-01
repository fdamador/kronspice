// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command coiled_tether takes a 5-mile (8,046.72 m) copper tether in low
// Earth orbit, coils it up, and asks whether "multiplication between
// windings" makes it an energy source. Three windings of the same wire are
// compared:
//
//   - straight: the ordinary electrodynamic tether (examples/tether);
//   - coiled: an ordinary single-layer coil, 2 m across, 1,281 turns;
//   - bifilar: the same coil wound as Tesla's US 512,340 (two strands side
//     by side, the end of the first joined to the start of the second), so
//     each turn lies next to the turn half a coil away (TurnVoltage = V/2).
//
// Energy comes only from the declared exterior sources the user named:
//
//   - the Sun: a solar array (total solar irradiance 1,360.8 W/m², G. Kopp
//     & J. Lean, Geophys. Res. Lett. 38, L01706, 2011; array area and
//     efficiency declared below) feeding a drive coupled into the winding
//     loop at the coil's resonance — the high-Q LC case;
//   - Earth's magnetic field: the motional EMF ∫(v×B)·dl along the actual
//     wire path (em.MotionalEMF), and a coil spinning in the field;
//   - Earth's electric field: in the craft's frame the same physics appears
//     as the electric field E′ = v×B (Lorentz transformation, em/frames.go),
//     so it is the motional EMF again, not a second source.
//
// The resonant cases are simulated in kronspice (.TRAN) with the energy
// ledger on. Between adjacent turns the run reports the voltage, the field,
// the electric energy stored there and a vacuum multipactor check
// (em.Multipactor: electron multiplication by secondary emission, Furman
// copper yield from Baglin et al. 2002). Multiplication of voltage and of
// electrons is real; the ledger shows where every joule comes from.
//
// Declared, not looked up: v = 7,700 m/s and B = 30 µT (order of magnitude
// in LEO, as in examples/tether; a real orbit takes B from IGRF); array
// 2 m² at 30 % efficiency; bare turns 1 mm apart in vacuum.
//
//	go run ./examples/coiled_tether
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
	"github.com/fdamador/kronspice/spacetime"
)

// Declared inputs.
const (
	FiveMiles   = 5 * 1609.344 // m (international mile)
	OrbitSpeed  = 7700.0       // m/s, relative to the co-rotating field
	FieldB      = 30e-6        // T
	TSI         = 1360.8       // W/m² (Kopp & Lean 2011)
	ArrayArea   = 2.0          // m²
	ArrayEff    = 0.30         // declared
	CoilDiam    = 2.0          // m
	WireDiam    = 2e-3         // m
	Pitch       = 3e-3         // m: 1 mm bare vacuum gap between turns
	InsulEps    = 1.0          // vacuum between bare turns
	SpinPerSec  = 1.0          // rev/s, for the spinning-coil case
	DesignField = 10e6         // V/m: declared limit for the turn gap, to be replaced by the chosen insulation's data-sheet value
)

// Coils returns the ordinary and bifilar windings of the five-mile wire.
func Coils() (ord, bif em.Coil) {
	n := int(math.Round(FiveMiles / (math.Pi * CoilDiam)))
	ord = em.Coil{Turns: n, Diameter: CoilDiam, WireDiameter: WireDiam, Pitch: Pitch, InsulationEps: InsulEps, Conductor: material.Copper}
	bif = ord
	bif.Bifilar = true
	return
}

// SolarPower available to the drive [W].
func SolarPower() float64 { return TSI * ArrayArea * ArrayEff }

// helix returns the wire path of a coil: axis along z, turns about it,
// starting at z = 0. For the bifilar coil the two strands are interleaved
// and the return lead from the end of strand 1 to the start of strand 2
// runs back along the coil.
func helix(c em.Coil) [][3]float64 {
	const perTurn = 36
	r := c.Diameter / 2
	var p [][3]float64
	strand := func(turns int, z0, pitch float64) {
		for k := 0; k <= turns*perTurn; k++ {
			th := 2 * math.Pi * float64(k) / perTurn
			p = append(p, [3]float64{r * math.Cos(th), r * math.Sin(th), z0 + pitch*th/(2*math.Pi)})
		}
	}
	if !c.Bifilar {
		strand(c.Turns, 0, c.Pitch)
		return p
	}
	h := c.Turns / 2
	strand(h, 0, 2*c.Pitch)
	end := p[len(p)-1]
	for k := 1; k <= 20; k++ { // return lead back to z = Pitch along the outside
		f := float64(k) / 20
		p = append(p, [3]float64{end[0], end[1], end[2] + (c.Pitch-end[2])*f})
	}
	strand(c.Turns-h, c.Pitch, 2*c.Pitch)
	return p
}

// MotionalEMFs from Earth's field for each winding, translating at orbital
// speed along x through B along y (v×B along z, the coil axis — the
// orientation that gives a coil its largest EMF).
func MotionalEMFs() (straight, coiled, bifilar, folded float64) {
	vel := func([3]float64) [3]float64 { return [3]float64{OrbitSpeed, 0, 0} }
	fld := func([3]float64) [3]float64 { return [3]float64{0, FieldB, 0} }
	line := func(z0, z1 float64) [][3]float64 {
		var p [][3]float64
		for k := 0; k <= 100; k++ {
			p = append(p, [3]float64{0, 0, z0 + (z1-z0)*float64(k)/100})
		}
		return p
	}
	straight = em.MotionalEMF(line(0, FiveMiles), vel, fld)
	ord, bif := Coils()
	coiled = em.MotionalEMF(helix(ord), vel, fld)
	bifilar = em.MotionalEMF(helix(bif), vel, fld)
	// A straight bifilar tether: out 2.5 miles along one strand and back
	// along the other, the two ends joined at the far end.
	out := line(0, FiveMiles/2)
	back := line(FiveMiles/2, 0)
	for i := range back {
		back[i][0] = 0.01 // the second strand 1 cm away
	}
	folded = em.MotionalEMF(append(out, back...), vel, fld)
	return
}

// Resonant is one solar-driven resonant run.
type Resonant struct {
	Name                     string
	L, C, R, F, Q            float64
	VDrive                   float64 // drive amplitude [V]
	VCoil, ICoil             float64 // peak at steady state
	VTurn, ETurn             float64 // adjacent-turn voltage [V] and field [V/m]
	WTurns, WTotal           float64 // electric energy between turns, total stored [J]
	PIn, PHeat, PLoad        float64 // mean over the last cycles [W]
	EIn, EDiss, EStored, Mis float64 // whole-run ledger [J]
	WireOverLambda           float64
	Multipactor, EImpact     float64
	Note                     string
}

// Loop is the resonant circuit: the solar-fed drive in series in the
// winding loop, wire resistance, coil inductance, C across the coil (and
// an optional load in series below it).
func Loop(name string, rw, l, c, rLoad, vDrive, f float64) (*circuit.Circuit, error) {
	b := circuit.NewBuilder(name)
	b.Add(circuit.Element{Name: "VSUN", Kind: circuit.VoltageSource,
		Wave:    circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, vDrive, f}},
		Channel: exchange.Light, ChannelNote: "solar array → drive inverter"}, "d", "0")
	b.Resistor("RWIRE", "d", "m", rw).Inductor("LCOIL", "m", "top", l)
	if rLoad > 0 {
		b.Capacitor("CWIND", "top", "x", c).Resistor("RLOAD", "x", "0", rLoad)
	} else {
		b.Capacitor("CWIND", "top", "0", c)
	}
	return b.BuildTransient()
}

// RunResonant drives a coil at its resonance from the solar array: the
// drive sits in series in the winding loop (as a primary loop would couple
// it), with C across the coil terminals (the bifilar coil's own
// inter-winding capacitance, or a tank capacitor for the ordinary coil).
// rLoad > 0 adds a load in series in the loop that takes power out.
func RunResonant(name string, c em.Coil, cTotal, rLoad float64) (Resonant, error) {
	l := c.Inductance()
	f := 1 / (2 * math.Pi * math.Sqrt(l*cTotal))
	rw := c.Resistance(f)
	r := Resonant{Name: name, L: l, C: cTotal, R: rw, F: f}
	r.Q = 2 * math.Pi * f * l / (rw + rLoad)
	// At resonance the loop draws V²/(2R_total): set the amplitude so the
	// array's whole output goes in.
	r.VDrive = math.Sqrt(2 * (rw + rLoad) * SolarPower())
	ckt, err := Loop(name, rw, l, cTotal, rLoad, r.VDrive, f)
	if err != nil {
		return r, err
	}
	tau := 2 * l / (rw + rLoad)
	stop := 7 * tau
	per := 1 / f
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: per / 200, Stop: stop, MaxStep: per / 200}, analysis.TranOptions{})
	if err != nil {
		return r, err
	}
	vt, _ := w.NodeVoltage("top")
	vx, okx := w.NodeVoltage("x")
	il, _ := w.ElementCurrent("LCOIL")
	// Mean powers and peaks over the last 20 cycles. (200 points per
	// cycle keep the trapezoidal rule's frequency warping, (ωh)²/12, well
	// inside the resonance's bandwidth 1/Q.)
	t0 := stop - 20*per
	var span float64
	pl := []float64(nil)
	for i, n := range w.ElementNames {
		if n == "RLOAD" {
			pl = w.P[i]
		}
	}
	for k := 1; k < len(w.Time); k++ {
		if w.Time[k] < t0 {
			continue
		}
		dt := w.Time[k] - w.Time[k-1]
		span += dt
		r.PIn += 0.5 * (w.Ledger.Inflow[k] + w.Ledger.Inflow[k-1]) * dt
		r.PHeat += 0.5 * (w.Ledger.Dissipation[k] + w.Ledger.Dissipation[k-1]) * dt
		if pl != nil {
			r.PLoad += 0.5 * (pl[k] + pl[k-1]) * dt
		}
		v := vt[k]
		if okx {
			v -= vx[k]
		}
		r.VCoil = math.Max(r.VCoil, math.Abs(v))
		r.ICoil = math.Max(r.ICoil, math.Abs(il[k]))
		r.WTotal = math.Max(r.WTotal, w.Ledger.StoredEnergy[k])
	}
	r.PIn /= span
	r.PHeat /= span
	r.PLoad /= span
	if pl != nil {
		r.PHeat -= r.PLoad // dissipation includes the load
	}
	r.EIn, r.EDiss, r.EStored, r.Mis = w.Ledger.EnergyIn, w.Ledger.EnergyDissipated, w.Ledger.StoredEnd-w.Ledger.StoredStart, w.Ledger.IntegrationMismatch
	r.VTurn = c.TurnVoltage(r.VCoil)
	r.ETurn = r.VTurn / (c.Pitch - c.WireDiameter)
	r.WTurns = c.TurnToTurnEnergy(r.VCoil)
	r.WireOverLambda = c.WireLength() * f / spacetime.C
	// While the coil rings up the gap sees every amplitude from 0 to VTurn:
	// keep the worst.
	for i := 1; i <= 100; i++ {
		if g, e := em.Multipactor(c.Pitch-c.WireDiameter, f, r.VTurn*float64(i)/100, 0, em.CopperAsReceived); g > r.Multipactor {
			r.Multipactor, r.EImpact = g, e
		}
	}
	return r, nil
}

func main() {
	ord, bif := Coils()
	if err := ord.Validate(); err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Five-mile copper tether, %.1f mm wire, %.0f m long, %.1f Ω\n", WireDiam*1e3, ord.WireLength(), ord.Resistance(0))
	fmt.Printf("Coiled: %d turns, %.0f m across, %.2f m long, %.1f mm bare gap between turns (vacuum)\n\n", ord.Turns, CoilDiam, ord.Length(), (Pitch-WireDiam)*1e3)

	// 1. Earth's magnetic field (and the electric field E' = v×B the craft sees).
	s, c, b, fo := MotionalEMFs()
	fmt.Println("1. Earth's magnetic field: motional EMF ∫(v×B)·dl along the wire (orbit 7.7 km/s, 30 µT)")
	fmt.Printf("   craft-frame electric field E' = v×B = %.3f V/m (the same coupling seen as an electric field)\n", OrbitSpeed*FieldB)
	fmt.Printf("   straight 5-mile tether      %10.3f V\n", s)
	fmt.Printf("   coiled (ordinary)           %10.3f V   = E' × the coil's %.2f m length: the turns cancel\n", c, ord.Length())
	fmt.Printf("   coiled (Tesla bifilar)      %10.3f V   same: both strands run the same way through the field\n", b)
	fmt.Printf("   straight bifilar (out/back) %10.3g V   out and back cancel\n", fo)
	grad := float64(ord.Turns) * math.Pi * CoilDiam * CoilDiam / 4 * OrbitSpeed * 3 * FieldB / 6.78e6
	fmt.Printf("   field gradient across the orbit adds at most N·A·v·3B/r ≈ %.2g V (a dipole field)\n", grad)
	fmt.Printf("   All of it is quasi-DC: the field seen by the craft changes over an orbit (~90 min), far below any LC resonance.\n")
	omega := 2 * math.Pi * SpinPerSec
	fmt.Printf("   coil spinning at %.0f rev/s in the field: AC EMF %.3f V peak at %.0f Hz (energy from the spin, which slows)\n\n",
		SpinPerSec, float64(ord.Turns)*FieldB*math.Pi*CoilDiam*CoilDiam/4*omega, SpinPerSec)

	// 2. The Sun: solar-driven resonance, ordinary and bifilar.
	cb := bif.SelfCapacitance()
	fmt.Printf("2. The Sun: solar array %.0f m² × %.0f %% × %.1f W/m² = %.0f W into a drive at resonance (high-Q LC)\n",
		ArrayArea, ArrayEff*100, TSI, SolarPower())
	fmt.Printf("   bifilar inter-winding capacitance %.3g F; ordinary coil's own %.3g F, so it gets a tank capacitor of the same value\n", cb, ord.SelfCapacitance())
	var runs []Resonant
	for _, k := range []struct {
		name string
		c    em.Coil
		load float64
	}{{"ordinary coil + tank capacitor", ord, 0}, {"Tesla bifilar coil", bif, 0}, {"Tesla bifilar coil, matched load", bif, 0}} {
		load := k.load
		if k.name == "Tesla bifilar coil, matched load" {
			load = k.c.Resistance(1 / (2 * math.Pi * math.Sqrt(k.c.Inductance()*cb)))
		}
		r, err := RunResonant(k.name, k.c, cb, load)
		if err != nil {
			log.Fatal(k.name, ": ", err)
		}
		runs = append(runs, r)
	}
	for _, r := range runs {
		fmt.Printf("\n   == %s ==\n", r.Name)
		fmt.Printf("   L %.3f H, C %.3g F, f %.1f Hz, wire R %.1f Ω, Q %.1f; wire = %.3f λ (lumped model valid below 0.1)\n", r.L, r.C, r.F, r.R, r.Q, r.WireOverLambda)
		fmt.Printf("   drive %.1f V → coil %.0f V peak: voltage multiplied ×%.0f; current %.2f A\n", r.VDrive, r.VCoil, r.VCoil/r.VDrive, r.ICoil)
		fmt.Printf("   between adjacent turns: %.1f V, %.3g V/m across the gap", r.VTurn, r.ETurn)
		if r.ETurn > DesignField {
			fmt.Printf("  — above the %.0f MV/m declared gap limit: needs insulation or wider spacing", DesignField/1e6)
		}
		fmt.Printf("\n   electric energy between the turns %.4g J of %.3f J stored in all\n", r.WTurns, r.WTotal)
		fmt.Printf("   multipactor between turns (vacuum, copper), worst over the ring-up: growth per impact %.3f (%s)\n", r.Multipactor, map[bool]string{true: "electrons multiply", false: "no avalanche at this frequency"}[r.Multipactor > 1])
		fmt.Printf("   steady state: Sun in %.1f W = heat %.1f W", r.PIn, r.PHeat)
		if r.PLoad > 0 {
			fmt.Printf(" + load %.1f W", r.PLoad)
		}
		fmt.Printf("\n   whole run: in %.2f J = heat %.2f J + stored %.3f J (+ integration error %.2g J)\n", r.EIn, r.EDiss, r.EStored, r.Mis)
	}

	// What it would take for Earth's field to drive the resonance.
	f := runs[1].F
	rim := math.Pi * CoilDiam * f
	fmt.Printf("\n3. Driving the resonance from Earth's field instead would need the coil to spin at %.0f rev/s:\n", f)
	fmt.Printf("   rim speed %.0f m/s against %.0f m/s, the burst speed of T1000G carbon fibre (the strongest rotor material in spacetime/) — it would fly apart.\n",
		rim, spacetime.BurstSpeed(spacetime.T1000G))

	// Where multipactor would happen in the same 1 mm gap: GHz drive, a
	// few hundred volts (scan of amplitudes at each frequency).
	fmt.Println("\n4. Where electrons would multiply in the same 1 mm vacuum gap (worst amplitude up to 2 kV):")
	for _, ff := range []float64{566, 1e6, 1e8, 1e9, 3e9} {
		var g, at float64
		for i := 1; i <= 100; i++ {
			v := 2000 * float64(i) / 100
			if gi, _ := em.Multipactor(Pitch-WireDiam, ff, v, 0, em.CopperAsReceived); gi > g {
				g, at = gi, v
			}
		}
		fmt.Printf("   %8.3g Hz: worst growth %.3f at %.0f V", ff, g, at)
		if g > 1 {
			fmt.Printf("  — avalanche; its electrons take their energy from the drive (and at this frequency the 8 km wire is many wavelengths: not a coil resonance)")
		}
		fmt.Println()
	}
	fmt.Println("\nMultiplication is real for voltage (Q) and for electrons (multipactor); the ledger traces every joule to the Sun's array.")
}
