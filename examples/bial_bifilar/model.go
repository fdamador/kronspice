// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
)

// Variant selects how the bismuth–aluminium element is built.
type Variant int

const (
	// Dry: Bi and Al plates with a solid dielectric. The plates' leads are
	// a Bi–Al thermocouple: with the plates ΔT hotter than the coil
	// connections, a Seebeck EMF (S_Al − S_Bi)·ΔT acts in series with the
	// plate capacitance. Energy enters through the thermal port.
	Dry Variant = iota
	// Wet: Bi and Al plates with an electrolyte between them: a galvanic
	// cell (Al anode, Bi cathode) behind the electrolyte resistance and a
	// Randles electrode model, with the plates' geometric capacitance in
	// parallel. Energy enters through the chemical port.
	Wet
	// Thermopile: the dry element charged by n Bi–Al junction pairs in
	// series (legs of bismuth and aluminium wire running from the hot plates
	// to the cold coil connections): EMF n·(S_Al − S_Bi)·ΔT behind the legs'
	// resistance. Energy enters through the thermal port; the legs also
	// conduct heat from hot to cold, which is reported.
	Thermopile
)

func (v Variant) String() string {
	switch v {
	case Wet:
		return "wet (galvanic)"
	case Thermopile:
		return "thermopile (Seebeck, n pairs)"
	}
	return "dry (Seebeck)"
}

// seebeck reports whether the element is driven by the Seebeck effect
// (dry or thermopile): plates in series with a thermoelectric EMF.
func (v Variant) seebeck() bool { return v != Wet }

// short is the variant's one-word name, used in file names.
func (v Variant) short() string { return strings.Fields(v.String())[0] }

// Fixed assumptions of the experiment, each stated so it can be replaced by
// a measured value.
type Assumptions struct {
	DeltaT         float64 // dry: plates above the coil connections [K]
	PlateThickness float64 // each plate [m]
	EMFFactor      float64 // wet: fraction of the standard EMF reached (1 = thermodynamic upper bound)
	ElectrolyteRho float64 // wet: electrolyte resistivity [Ω·m]
	DoubleLayer    float64 // wet: double-layer capacitance per electrode area [F/m²]
	ChargeTransfer float64 // wet: charge-transfer resistance × area per electrode [Ω·m²]
	WaterEps       float64 // wet: relative permittivity of the electrolyte
	Line           material.SpectralLine
}

// DefaultAssumptions are typical values, not measurements of a built device.
func DefaultAssumptions() Assumptions {
	return Assumptions{
		DeltaT:         100,
		PlateThickness: 1e-3,
		EMFFactor:      1,
		// Sea water (S = 35, 20 °C) conducts about 4.8 S/m (UNESCO 1981,
		// practical salinity scale): ρ ≈ 0.21 Ω·m.
		ElectrolyteRho: 0.21,
		// 10–40 µF/cm² is typical for metal electrodes in aqueous
		// electrolytes (Bard & Faulkner, Electrochemical Methods, ch. 13):
		// 20 µF/cm² = 0.2 F/m².
		DoubleLayer: 0.2,
		// 10 Ω·cm² (= 1e-3 Ω·m²): an assumed charge-transfer resistance;
		// it depends strongly on the surface and must be measured.
		ChargeTransfer: 1e-3,
		WaterEps:       80,
		Line:           material.LineCopper,
	}
}

// Design is one candidate: coil geometry plus the Bi–Al element.
type Design struct {
	Variant       Variant
	Coil          em.Coil
	PlateSide     float64 // square plates [m]
	Gap           float64 // dielectric or electrolyte thickness [m]
	DielectricEps float64 // dry and thermopile
	// Thermopile only: junction pairs and the size of each Bi and Al leg.
	Pairs                  int
	LegDiameter, LegLength float64
	A                      Assumptions
}

// geneNames lists the experiment's genes for a variant, in order.
func geneNames(v Variant) []string {
	n := []string{"turns", "coil_diameter_m", "wire_diameter_m", "pitch_ratio", "insulation_eps",
		"plate_side_m", "gap_m", "dielectric_eps"}
	if v == Thermopile {
		n = append(n, "pairs", "leg_diameter_m", "leg_length_m")
	}
	return n
}

// designFrom maps decoded gene values onto a design.
func designFrom(v Variant, a Assumptions, bifilar bool, x []float64) Design {
	d := Design{
		Variant: v,
		Coil: em.Coil{Turns: int(x[0]), Diameter: x[1], WireDiameter: x[2], Pitch: x[2] * x[3],
			InsulationEps: x[4], Bifilar: bifilar, Conductor: material.Copper},
		PlateSide: x[5], Gap: x[6], DielectricEps: x[7], A: a,
	}
	if v == Thermopile && len(x) >= 11 {
		d.Pairs, d.LegDiameter, d.LegLength = int(x[8]), x[9], x[10]
	}
	return d
}

// Element is the lumped model of the Bi–Al element.
type Element struct {
	EMF     float64            // open-circuit EMF [V]
	Channel exchange.ChannelID // where its energy comes from
	Note    string
	Rint    float64 // series resistance: plate leads (dry) or electrolyte (wet) [Ω]
	Cplate  float64 // dry: plate capacitance in series with the EMF [F]
	Cgeo    float64 // wet: geometric plate capacitance across the terminals [F]
	Cdl     float64 // wet: two double layers in series [F]
	Rct     float64 // wet: two charge-transfer resistances in series [Ω]
	// CapacityC is the charge the aluminium plate can deliver (wet) [C].
	CapacityC float64
	Volume    float64 // [m³]
	// HeatLeak is the heat the thermopile legs conduct from the hot plates
	// to the cold connections at ΔT [W]. It is not electrical: it is the
	// heat a source must supply just to hold ΔT, reported for scale.
	HeatLeak float64
}

// ElementModel computes the Bi–Al element's lumped parameters.
func (d Design) ElementModel() (Element, error) {
	a := d.PlateSide * d.PlateSide
	t := d.A.PlateThickness
	var e Element
	e.Volume = a * (d.Gap + 2*t)
	switch d.Variant {
	case Dry, Thermopile:
		c, err := em.ParallelPlate(a, d.Gap, d.DielectricEps)
		if err != nil {
			return e, err
		}
		e.Cplate = c
		e.EMF = material.ThermocoupleEMF(material.Bismuth, material.Aluminium, d.A.DeltaT)
		e.Channel = exchange.Thermal
		e.Note = fmt.Sprintf("Bi–Al Seebeck EMF, (S_Al − S_Bi)·ΔT, ΔT = %g K: heat from the hot plates", d.A.DeltaT)
		// Current crosses each plate edge to edge (length ≈ side, cross
		// section side × t): R = ρ/t per plate.
		e.Rint = (material.Bismuth.Resistivity + material.Aluminium.Resistivity) / t
		if d.Variant == Thermopile {
			if d.Pairs < 1 || !(d.LegDiameter > 0) || !(d.LegLength > 0) {
				return e, fmt.Errorf("%w: thermopile needs pairs ≥ 1 and positive leg size", em.ErrOutsideModelValidity)
			}
			if d.LegLength < 5*d.LegDiameter {
				return e, fmt.Errorf("%w: thermopile legs must be at least 5 diameters long (1-D conduction)", em.ErrOutsideModelValidity)
			}
			n := float64(d.Pairs)
			la := math.Pi * d.LegDiameter * d.LegDiameter / 4
			e.EMF *= n
			e.Rint += n * (material.Bismuth.Resistivity + material.Aluminium.Resistivity) * d.LegLength / la
			e.Volume += n * 2 * la * d.LegLength
			e.HeatLeak = n * (material.Bismuth.ThermalConductivity + material.Aluminium.ThermalConductivity) *
				la / d.LegLength * d.A.DeltaT
			e.Note = fmt.Sprintf("thermopile, %d Bi–Al pairs × (S_Al − S_Bi)·ΔT, ΔT = %g K: heat from the hot plates",
				d.Pairs, d.A.DeltaT)
		}
	case Wet:
		cg, err := em.ParallelPlate(a, d.Gap, d.A.WaterEps)
		if err != nil {
			return e, err
		}
		e.Cgeo = cg
		e.EMF = d.A.EMFFactor * material.CellEMF(material.Aluminium, material.Bismuth)
		e.Channel = exchange.Chemical
		e.Note = "Al|electrolyte|Bi galvanic cell (Al anode): chemical energy of the aluminium plate"
		e.Rint = d.A.ElectrolyteRho * d.Gap / a
		e.Cdl = d.A.DoubleLayer * a / 2
		e.Rct = 2 * d.A.ChargeTransfer / a
		e.CapacityC = material.ChargeCapacity(material.Aluminium, a*t*material.Aluminium.Density)
	}
	return e, nil
}

// Build returns the circuit: the coil (series R and L, self-capacitance
// across its terminals) in parallel with the Bi–Al element, between node
// "a" and ground, plus a 1 A AC probe current into "a" (DC 0) to measure the
// tank impedance. rCoil is the coil resistance to use (it depends on the
// frequency through the skin effect). With step set, the EMF switches on at
// t = 0 (an idealized step; a real thermal or chemical source turns on
// over milliseconds or longer).
func (d Design) Build(rCoil float64, step bool) (*circuit.Circuit, error) {
	if err := d.Coil.Validate(); err != nil {
		return nil, err
	}
	el, err := d.ElementModel()
	if err != nil {
		return nil, err
	}
	b := circuit.NewBuilder(fmt.Sprintf("Bi–Al %s element with %s coil", d.Variant, coilKind(d.Coil.Bifilar)))
	b.Resistor("Rcoil", "a", "c", rCoil).
		Inductor("Lcoil", "c", "0", d.Coil.Inductance()).
		Capacitor("Cself", "a", "0", d.Coil.SelfCapacitance())
	src := circuit.Element{Name: "Vbial", Kind: circuit.VoltageSource, Value: el.EMF, HasDC: true,
		Channel: el.Channel, ChannelNote: el.Note}
	if step {
		src.HasDC = false
		src.Value = 0
		src.Wave = circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, 1e-12, el.EMF}}
	}
	switch d.Variant {
	case Dry, Thermopile:
		b.Resistor("Rint", "a", "n1", el.Rint).Add(src, "n1", "n2").Capacitor("Cplate", "n2", "0", el.Cplate)
	case Wet:
		b.Resistor("Relec", "a", "n1", el.Rint).Add(src, "n1", "n2").
			Resistor("Rct", "n2", "0", el.Rct).Capacitor("Cdl", "n2", "0", el.Cdl).
			Capacitor("Cgeo", "a", "0", el.Cgeo)
	}
	b.Add(circuit.Element{Name: "Iprobe", Kind: circuit.CurrentSource, HasDC: true, ACMag: 1}, "0", "a")
	if step {
		return b.BuildTransient()
	}
	return b.Build()
}

func coilKind(bif bool) string {
	if bif {
		return "bifilar"
	}
	return "ordinary"
}

// ErrNoResonance means the tank has no parallel resonance (its impedance
// phase never crosses from inductive to capacitive) in the searched range.
var ErrNoResonance = fmt.Errorf("%w: no parallel resonance", analysis.ErrNoMeasurement)

// impedance is Z(f) of the tank at node a (1 A probe).
func impedance(c *circuit.Circuit, f float64, o analysis.Options) (complex128, error) {
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, o)
	if err != nil {
		return 0, err
	}
	v, _ := r.NodeVoltage("a")
	return v[0], nil
}

// Resonance finds the tank's parallel resonance — where the impedance
// phase falls through zero, from inductive to capacitive — and its quality
// factor from the phase slope, Q = (f0/2)·|dφ/df| (exact for a parallel RLC
// and equal there to f0/bandwidth; defined also when Q is too low for −3 dB
// points to exist). A coarse sweep from fLo to fHi brackets the crossing
// with the largest |Z|; bisection then fixes f0 to 1e-12 relative.
func Resonance(c *circuit.Circuit, fLo, fHi float64, o analysis.Options) (f0, q, zmax float64, err error) {
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 40, FStart: fLo, FStop: fHi}, o)
	if err != nil {
		return 0, 0, 0, err
	}
	v, _ := r.NodeVoltage("a")
	best := -1
	for k := 1; k < len(v); k++ {
		if cmplx.Phase(v[k-1]) > 0 && cmplx.Phase(v[k]) <= 0 {
			if best < 0 || cmplx.Abs(v[k]) > cmplx.Abs(v[best]) {
				best = k
			}
		}
	}
	if best < 0 {
		return 0, 0, 0, ErrNoResonance
	}
	lo, hi := r.Freq[best-1], r.Freq[best]
	for i := 0; i < 80 && hi/lo-1 > 1e-13; i++ {
		mid := math.Sqrt(lo * hi)
		z, err := impedance(c, mid, o)
		if err != nil {
			return 0, 0, 0, err
		}
		if cmplx.Phase(z) > 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	f0 = math.Sqrt(lo * hi)
	h := 1e-5 * f0
	zp, err1 := impedance(c, f0+h, o)
	zm, err2 := impedance(c, f0-h, o)
	z0, err3 := impedance(c, f0, o)
	if err := errors.Join(err1, err2, err3); err != nil {
		return 0, 0, 0, err
	}
	slope := (cmplx.Phase(zp) - cmplx.Phase(zm)) / (2 * h)
	return f0, f0 / 2 * math.Abs(slope), cmplx.Abs(z0), nil
}

// Result is everything measured for one design.
type Result struct {
	Design          Design
	Element         Element
	L, Cself, Rcoil float64
	F0, Q, Zmax     float64
	StoredJ         float64 // stored energy at the steady state with the EMF on [J]
	Volume          float64 // coil + element [m³]
	EnergyDensity   float64 // [J/m³]
	SubN            float64 // nearest whole-number subharmonic of the spectral line
	SubOffset       float64 // (f0 − f_line/n)/f0
	SourcePower     float64 // power from the Bi–Al source at the steady state [W]
	Heat            float64 // dissipation at the steady state [W]
}

// Measure builds the design, computes its stored energy at the steady
// state, then finds its resonance (with the coil resistance taken at that
// frequency). When the tank has no resonance the steady-state fields are
// still filled in and ErrNoResonance is returned.
func Measure(d Design, o analysis.Options) (Result, error) {
	res := Result{Design: d, L: d.Coil.Inductance(), Cself: d.Coil.SelfCapacitance()}
	el, err := d.ElementModel()
	if err != nil {
		return res, err
	}
	res.Element = el
	res.Volume = d.Coil.Volume() + el.Volume
	// Steady state with the EMF on: the DC operating point.
	c, err := d.Build(d.Coil.Resistance(0), false)
	if err != nil {
		return res, err
	}
	op, err := analysis.OP(c, o)
	if err != nil {
		return res, err
	}
	for _, e := range op.Elements {
		res.StoredJ += e.StoredEnergy
		if e.Name == "Vbial" {
			res.SourcePower = -e.PowerAbsorbed
		}
	}
	res.Heat = op.Ledger.Dissipation
	res.EnergyDensity = res.StoredJ / res.Volume
	// Resonance with the DC resistance first, then with R at f0.
	rc := d.Coil.Resistance(0)
	var f0 float64
	for pass := 0; pass < 2; pass++ {
		c, err := d.Build(rc, false)
		if err != nil {
			return res, err
		}
		f0, res.Q, res.Zmax, err = Resonance(c, 1, 1e10, o)
		if err != nil {
			return res, err
		}
		rc = d.Coil.Resistance(f0)
	}
	res.F0, res.Rcoil = f0, rc
	// Lumped (quasi-static) validity, CLAUDE.md §3c: the wire must be short
	// against the wavelength at the frequencies that matter (λ/10).
	if lambda := material.SpeedOfLight / f0; d.Coil.WireLength() > 0.1*lambda {
		return res, fmt.Errorf("%w: wire length %.3g m exceeds λ/10 = %.3g m at the resonance %.4g Hz (not quasi-static)",
			em.ErrOutsideModelValidity, d.Coil.WireLength(), 0.1*lambda, f0)
	}
	res.SubN, res.SubOffset = d.A.Line.Subharmonic(f0)
	return res, nil
}
