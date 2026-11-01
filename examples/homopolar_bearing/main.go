// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command homopolar_bearing designs Tesla's two-disk homopolar generator
// (U.S. Patent 406,968) with its load current taken through the shaft ball
// bearings, so that each bearing is also a Huber ball-bearing motor
// (kronspice P element). Circuit simulation only: every part is a lumped
// element, the shafts are mechanical nodes (across = ω, through = torque).
//
// The goal: deliver energy to an external resistor or capacitor with the
// least drive torque per watt. The trade-off: the bearing bore sets the
// axle radius r_i, so bigger bearings remove magnetic area from each disk
// (EMF = ½Bω(r_o² − r_i²)); but bigger, more heavily loaded bearings carry
// more current before they are damaged, and have more friction.
//
// Measured inputs (no coefficient is assumed; unmeasured ones are swept):
//
//   - Bearing voltage under DC: almost constant, 0.73 V average at 4 A
//     (range 0.5–1.5 V reported), from repeated breakdown of the oxide and
//     lubricant film (O. Safdarzadeh, A. Farahi, A. Binder et al.,
//     Lubricants 12(5), 148, 2024, DOI 10.3390/lubricants12050148; greased
//     axial bearings 51208, d = 40 mm).
//   - Damage limit: apparent current density 0.1–0.3 A/mm² of Hertzian
//     contact area; A_Hz = 1.0 mm² at 400 N and 2.5 mm² at 1600 N for that
//     bearing (same paper). Other sizes and loads scale by Hertz contact
//     theory, A ∝ (F·R)^(2/3), with the ball radius R proportional to the
//     bore within a bearing series (K. L. Johnson, Contact Mechanics, 1985,
//     §4.2; T. A. Harris, Rolling Bearing Analysis).
//   - Friction torque M = ½·f·F·d, f = 0.0015 for single-row ball bearings
//     (SKF rolling-bearing catalogue simplified friction model; RoyMech
//     "Rolling Bearing Friction").
//   - Rim-speed limit of a solid copper disk: hoop stress (3+ν)/8·ρv² at the
//     centre (S. Timoshenko, Strength of Materials II, §44) ≤ yield of
//     annealed copper ≈ 69 MPa (ASM Handbook vol. 2, C10100 annealed),
//     ν = 0.34, with a declared safety factor.
//
// Unmeasured, therefore swept and reported as a function: the Huber
// coefficient κψ (examples/bearing_motor: no open-access table) in EM and
// thermal modes, the belt's rolling-contact resistance, and a shaft
// collector's resistance for the insulated-bearing comparison.
//
//	go run ./examples/homopolar_bearing [-p 20] [-gens 25] [-pop 30]
package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/material"
)

// Measured bearing data (see the package comment).
const (
	BearingVolt    = 0.73   // V, DC bearing voltage, Lubricants 12, 148 (2024)
	BearingVoltLo  = 0.5    // V, low end of reported range
	BearingVoltHi  = 1.5    // V, high end
	JLimit         = 0.3e6  // A/m², upper damage threshold (0.1–0.3 A/mm²)
	JLimitLow      = 0.1e6  // A/m², lower threshold
	AHzRef         = 1.0e-6 // m², at FRef on a dRef bore
	FRef           = 400    // N
	DRef           = 0.040  // m, 51208 bore
	FrictionCoeff  = 0.0015 // single-row ball bearing
	CopperYield    = 69e6   // Pa, annealed C10100
	CopperPoisson  = 0.34
	SafetyFactor   = 1.5 // declared design margin on yield
	gravity        = 9.80665
	transitionAmps = 1e-3 // width of the zero crossing of the bearing's v(i)
)

// Collector says how the load current reaches each shaft.
type Collector int

const (
	ThroughBearings   Collector = iota // current through the ball bearings (Huber motors)
	InsulatedBearings                  // bearings insulated; a shaft collector of declared resistance
)

// Design is one machine.
type Design struct {
	Ro, Thickness, B, Omega float64 // disk outer radius [m], thickness [m], field [T], speed [rad/s]
	Bore                    float64 // bearing bore = axle diameter [m]; r_i = Bore/2
	NBear                   int     // bearings per shaft, in parallel in the current path
	Preload                 float64 // load per bearing [N] (at least the rotor weight share)
	RL                      float64 // external load resistor [Ω]
	Collector               Collector
	RCollector              float64 // shaft collector resistance when insulated [Ω] (declared, swept)
	RBelt                   float64 // rim–belt–rim rolling contact [Ω] (declared, swept)
	VBearing                float64 // measured DC bearing voltage [V]
	RShare                  float64 // bearing series resistance [Ω] (declared, small: the measured voltage is almost constant; it also shares current between parallel bearings)
	KappaPsi                float64 // Huber κ·ψ(ω) per bearing [N·m/A²] (unmeasured: swept)
	Thermal                 bool    // Huber mode: thermal (η) instead of EM (back-EMF)
	Eta                     float64 // thermal-mode efficiency ceiling (Marinov 1989: 0.10)
}

// Defaults is a reference unit: 15 cm copper disks, 0.8 T, 40 mm bearings.
func Defaults() Design {
	return Design{Ro: 0.15, Thickness: 5e-3, B: 0.8, Omega: 600, Bore: 0.040, NBear: 2, Preload: 400,
		RL: 2, Collector: ThroughBearings, RCollector: 50e-6, RBelt: 100e-6,
		VBearing: BearingVolt, RShare: 1e-3, Eta: 0.10}
}

// Ri is the axle radius, where the disk's active annulus begins.
func (d Design) Ri() float64 { return d.Bore / 2 }

// K is one disk's EMF per unit speed, ½B(r_o² − r_i²) [V·s].
func (d Design) K() float64 { return 0.5 * d.B * (d.Ro*d.Ro - d.Ri()*d.Ri()) }

// DiskResistance is the radial spreading resistance ln(r_o/r_i)/(2πσt).
func (d Design) DiskResistance() float64 {
	return math.Log(d.Ro/d.Ri()) * material.Copper.Resistivity / (2 * math.Pi * d.Thickness)
}

// BearingLoad is the load on each bearing: the preload, but at least the
// disk's weight shared over the shaft's bearings [N].
func (d Design) BearingLoad() float64 {
	m := material.Copper.Density * math.Pi * (d.Ro*d.Ro - d.Ri()*d.Ri()) * d.Thickness
	return math.Max(d.Preload, m*gravity/float64(d.NBear))
}

// HertzArea is one bearing's apparent contact area [m²], scaled from the
// measured 51208 point by Hertz theory.
func (d Design) HertzArea() float64 {
	return AHzRef * math.Pow(d.BearingLoad()/FRef*d.Bore/DRef, 2.0/3)
}

// CurrentLimit is the largest load current the bearing bank carries at the
// damage threshold j [A/m²].
func (d Design) CurrentLimit(j float64) float64 {
	return j * d.HertzArea() * float64(d.NBear)
}

// FrictionTorque per shaft [N·m].
func (d Design) FrictionTorque() float64 {
	return float64(d.NBear) * 0.5 * FrictionCoeff * d.BearingLoad() * d.Bore
}

// RimSpeedLimit is the largest rim speed the copper disk carries [m/s].
func RimSpeedLimit() float64 {
	return math.Sqrt(8 * CopperYield / SafetyFactor / ((3 + CopperPoisson) * material.Copper.Density))
}

// machine returns the machine spec of one disk (as examples/homopolar).
func (d Design) machine(polarity float64) *circuit.MachineSpec {
	s := circuit.NewMachineSpec(2)
	s.L0[0][0] = 1e-9
	s.L0[1][1] = 1e-3
	s.G[0][1] = polarity * d.K()
	return s
}

// bearingBank adds the shaft's bearings between node a (frame terminal)
// and node b (axle): P elements with the measured contact drop, driving
// shaft node w with their Huber torque.
func (d Design) bearingBank(bl *circuit.Builder, tag, a, b, w string) {
	if d.Collector == InsulatedBearings {
		bl.Resistor("RC"+tag, a, b, d.RCollector)
		return
	}
	w0 := d.Omega // κ is scaled so κ·tanh(ω/ω0) = KappaPsi at the running speed
	kappa := d.KappaPsi / math.Tanh(1)
	for k := 0; k < d.NBear; k++ {
		bl.BearingMotorVF(fmt.Sprintf("P%s%d", tag, k), a, b, w, "0", d.RShare, d.VBearing, kappa, w0, d.Thermal, d.thermalEta())
	}
}

// HuberTorque is the Huber torque per shaft at load current i [N·m],
// the element's law at the running speed.
func (d Design) HuberTorque(i float64) float64 {
	if d.Collector != ThroughBearings || d.KappaPsi == 0 {
		return 0
	}
	ib := math.Abs(i) / float64(d.NBear)
	a := d.KappaPsi * ib * ib
	if d.Thermal {
		h := d.VBearing*math.Tanh(ib/1e-3)*ib + d.RShare*ib*ib
		a /= 1 + a*d.Omega/(d.Eta*h)
	}
	return float64(d.NBear) * a
}

func (d Design) thermalEta() float64 {
	if d.Thermal {
		return d.Eta
	}
	return 0
}

// Circuit builds the two-disk machine with load resistor RL across the
// terminals t1, t2. Shafts w1, w2 are driven at Omega by VW1, VW2 (their
// currents are the drive torques); bearing friction is a mechanical
// resistance ω/M_f on each shaft (exact at the driven speed).
func (d Design) Circuit() (*circuit.Circuit, error) { return d.build(false, 0, 0) }

// build makes the circuit; with capacitor it puts C in series with RL and
// ramps the field over tRamp (0 = step).
func (d Design) build(capacitor bool, c, tRamp float64) (*circuit.Circuit, error) {
	bl := circuit.NewBuilder("homopolar with ball-bearing current path").Domain(exchange.Mechanical, "w1", "w2")
	bl.MachineModel("dA", d.machine(1)).MachineModel("dB", d.machine(-1))
	rd := d.DiskResistance()
	field := func(name, n string) {
		if capacitor {
			// From zero field (the DC operating point has C uncharged) to
			// 1 A over tRamp; a "step" is a 1 ms ramp.
			bl.Add(circuit.Element{Name: name, Kind: circuit.CurrentSource,
				Wave: circuit.Waveform{Kind: circuit.WavePWL, Args: []float64{0, 0, math.Max(tRamp, 1e-3), 1}}}, "0", n)
		} else {
			bl.CurrentSource(name, "0", n, 1)
		}
		bl.Resistor("R"+name, n, "0", 1e3)
	}
	field("IFA", "fa")
	field("IFB", "fb")
	bl.VoltageSource("VW1", "w1", "0", d.Omega).VoltageSource("VW2", "w2", "0", d.Omega)
	if mf := d.FrictionTorque(); mf > 0 {
		bl.Resistor("RF1", "w1", "0", d.Omega/mf).Resistor("RF2", "w2", "0", d.Omega/mf)
	}
	d.bearingBank(bl, "1", "t1", "ax1", "w1")
	bl.Resistor("RD1", "ax1", "c1", rd)
	bl.Machine("NA", "dA", [][2]string{{"rim1", "c1"}, {"fa", "0"}}, "w1", "0")
	bl.Resistor("RBELT", "rim1", "rim2", d.RBelt)
	bl.Machine("NB", "dB", [][2]string{{"rim2", "c2"}, {"fb", "0"}}, "w2", "0")
	bl.Resistor("RD2", "c2", "ax2", rd)
	d.bearingBank(bl, "2", "ax2", "t2", "w2")
	if capacitor {
		bl.Resistor("RL", "t1", "tc", d.RL).Capacitor("CL", "tc", "t2", c)
	} else {
		bl.Resistor("RL", "t1", "t2", d.RL)
	}
	// The armature floats; a megohm references it.
	bl.Resistor("RREF", "t2", "0", 1e6)
	return bl.Build()
}

// Result is a solved operating point.
type Result struct {
	EMF          float64 // both disks, open circuit [V]
	ILoad, PLoad float64
	PShaft       float64    // net mechanical input from the drives [W]
	Torque       [2]float64 // drive torques [N·m]
	Friction     float64    // friction torque per shaft [N·m]
	Huber        float64    // Huber torque per shaft [N·m] (aiding)
	PBearingHeat float64    // heat in the bearings' contacts [W]
	JRatio       float64    // bearing current density / JLimit
	TorquePerW   float64    // (T1 + T2) / P_load [N·m/W]
	Efficiency   float64    // P_load / P_shaft
	Residual     float64
}

// Solve finds the operating point.
func (d Design) Solve() (*Result, error) {
	c, err := d.Circuit()
	if err != nil {
		return nil, err
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return nil, err
	}
	res := &Result{EMF: 2 * d.K() * d.Omega, Friction: d.FrictionTorque(), Residual: r.Ledger.RelativeResidual()}
	rl, _ := r.Element("RL")
	res.ILoad, res.PLoad = rl.Current, rl.PowerAbsorbed
	for k, n := range []string{"VW1", "VW2"} {
		el, _ := r.Element(n)
		res.Torque[k] = -el.Current
		res.PShaft += -el.Current * d.Omega
	}
	if d.Collector == ThroughBearings {
		for _, e := range r.Elements {
			if e.Kind == "bearing motor" {
				res.PBearingHeat += e.PowerAbsorbed
			}
		}
		res.Huber = d.HuberTorque(res.ILoad)
		res.JRatio = math.Abs(res.ILoad) / d.CurrentLimit(JLimit)
	}
	if res.PLoad > 0 {
		res.TorquePerW = (res.Torque[0] + res.Torque[1]) / res.PLoad
	}
	if res.PShaft > 0 {
		res.Efficiency = res.PLoad / res.PShaft
	}
	return res, nil
}

// Charge charges capacitor c through RL from rest, with the field stepped
// (tRamp = 0: a 1 ms ramp) or ramped over tRamp, and returns the energy stored in C and
// the drives' energy over tEnd.
func (d Design) Charge(c, tRamp, tEnd float64) (stored, shaft, mismatch float64, err error) {
	ckt, err := d.build(true, c, tRamp)
	if err != nil {
		return 0, 0, 0, err
	}
	wf, err := analysis.Tran(ckt, analysis.TranSpec{Step: tEnd / 400, Stop: tEnd}, analysis.TranOptions{})
	if err != nil {
		return 0, 0, 0, err
	}
	vc1, _ := wf.NodeVoltage("tc")
	vc2, _ := wf.NodeVoltage("t2")
	n := len(wf.Time) - 1
	v := vc1[n] - vc2[n]
	stored = 0.5 * c * v * v
	for _, name := range []string{"VW1", "VW2"} {
		i, _ := wf.ElementCurrent(name)
		for k := 1; k <= n; k++ {
			shaft += -0.5 * (i[k] + i[k-1]) * d.Omega * (wf.Time[k] - wf.Time[k-1])
		}
	}
	return stored, shaft, wf.Ledger.IntegrationMismatch, nil
}

// Problem evolves the unit for the least drive torque per load watt at an
// output of at least pTarget watts into the load resistor, within the
// measured bearing current limit and the disk's rim-speed limit.
func Problem(pTarget float64, base Design) evo.Problem {
	design := func(v []float64) Design {
		d := base
		d.Bore, d.Preload, d.Thickness, d.Omega, d.RL = v[0], v[1], v[2], v[3], v[4]
		d.NBear = int(math.Round(v[5]))
		return d
	}
	return evo.Problem{
		Space: evo.Space{
			{Name: "bore", Min: 0.010, Max: 0.12},
			{Name: "preload", Min: 50, Max: 5000},
			{Name: "thickness", Min: 1e-3, Max: 30e-3},
			{Name: "omega", Min: 50, Max: RimSpeedLimit() / base.Ro},
			{Name: "R_load", Min: 1e-3, Max: 100},
			{Name: "bearings", Min: 1, Max: 6},
		},
		Build: func(v []float64) (*circuit.Circuit, error) { return design(v).Circuit() },
		Objectives: []evo.Objective{
			{Name: "drive torque per load watt", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
				res, err := design(e.Values).Solve()
				if err != nil {
					return 0, err
				}
				if res.PLoad <= 0 {
					return 1e3, nil
				}
				return res.TorquePerW * 1e3, nil // mN·m per W
			}},
			{Name: "power shortfall", Weight: 100, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				res, err := design(e.Values).Solve()
				if err != nil {
					return 0, err
				}
				return math.Max(0, (pTarget-res.PLoad)/pTarget), nil
			}},
			{Name: "bearing current over limit", Weight: 100, Constraint: true, Eval: func(e *evo.Eval) (float64, error) {
				res, err := design(e.Values).Solve()
				if err != nil {
					return 0, err
				}
				return math.Max(0, res.JRatio-1), nil
			}},
		},
	}
}

func report(name string, d Design) *Result {
	res, err := d.Solve()
	if err != nil {
		log.Fatalf("%s: %v", name, err)
	}
	fmt.Printf("%s\n", name)
	fmt.Printf("  bore %.0f mm × %d/shaft, load/bearing %.0f N, disk r_o %.3f m r_i %.3f m t %.1f mm, ω %.0f rad/s (rim %.0f m/s), R_L %.3g Ω\n",
		d.Bore*1e3, d.NBear, d.BearingLoad(), d.Ro, d.Ri(), d.Thickness*1e3, d.Omega, d.Omega*d.Ro, d.RL)
	fmt.Printf("  EMF %.3f V; load %.3f A, %.2f W; shaft %.2f W; efficiency %.3f; drive torque %.4f + %.4f N·m = %.3g N·m/W\n",
		res.EMF, res.ILoad, res.PLoad, res.PShaft, res.Efficiency, res.Torque[0], res.Torque[1], res.TorquePerW)
	if d.Collector == ThroughBearings {
		fmt.Printf("  bearings: contact heat %.2f W, friction %.4f N·m/shaft, Huber %.4f N·m/shaft; current %.2f× the 0.3 A/mm² limit (%.2f× the 0.1 A/mm² one)\n",
			res.PBearingHeat, res.Friction, res.Huber, res.JRatio, res.JRatio*JLimit/JLimitLow)
	}
	fmt.Printf("  ledger residual %.1e\n", res.Residual)
	return res
}

func main() {
	gens := flag.Int("gens", 25, "generations")
	pop := flag.Int("pop", 30, "population")
	ptarget := flag.Float64("p", 20, "required load power [W]")
	flag.Parse()

	fmt.Printf("Tesla two-disk homopolar, load current through the shaft ball bearings (circuit simulation)\n")
	fmt.Printf("copper rim-speed limit %.0f m/s (yield %.0f MPa / %.1f)\n\n", RimSpeedLimit(), CopperYield/1e6, SafetyFactor)
	d := Defaults()
	report("Reference unit", d)

	fmt.Println("\nTrade-off: bearing bore (axle radius) against disk magnetic area, at the reference speed and load")
	fmt.Printf("  %6s %8s %9s %9s %10s %10s %10s\n", "bore", "area %", "EMF V", "I_max A", "P_load W", "N·m/W", "J/limit")
	for _, bore := range []float64{0.010, 0.020, 0.040, 0.060, 0.080, 0.100, 0.120} {
		x := d
		x.Bore = bore
		x.RL = 0
		// Match the load to deliver the bearing bank's current limit.
		ilim := x.CurrentLimit(JLimit)
		x.RL = math.Max(1e-3, (2*x.K()*x.Omega-2*x.VBearing)/ilim)
		r, err := x.Solve()
		if err != nil {
			log.Fatal(err)
		}
		area := (x.Ro*x.Ro - x.Ri()*x.Ri()) / (x.Ro * x.Ro) * 100
		fmt.Printf("  %4.0f mm %8.1f %9.3f %9.3f %10.2f %10.4f %10.2f\n", bore*1e3, area, r.EMF, ilim, r.PLoad, r.TorquePerW, r.JRatio)
	}

	fmt.Println("\nHuber torque (κψ unmeasured → swept), reference unit:")
	fmt.Printf("  %10s %8s %12s %12s %12s %10s\n", "κψ N·m/A²", "mode", "Huber N·m", "drive N·m", "P_load W", "N·m/W")
	for _, kp := range []float64{0, 1e-3, 1e-2, 1e-1} {
		for _, th := range []bool{false, true} {
			x := d
			x.KappaPsi, x.Thermal = kp, th
			r, err := x.Solve()
			if err != nil {
				log.Fatal(err)
			}
			mode := "EM"
			if th {
				mode = "thermal"
			}
			fmt.Printf("  %10.0e %8s %12.4f %12.4f %12.3f %10.4f\n", kp, mode, r.Huber, r.Torque[0]+r.Torque[1], r.PLoad, r.TorquePerW)
		}
	}

	fmt.Println("\nMeasured bearing-voltage range (0.5–1.5 V) and the insulated-bearing alternative:")
	for _, v := range []float64{BearingVoltLo, BearingVolt, BearingVoltHi} {
		x := d
		x.VBearing = v
		r, _ := x.Solve()
		fmt.Printf("  through bearings, %.2f V each: P_load %.2f W, efficiency %.3f, %.4f N·m/W\n", v, r.PLoad, r.Efficiency, r.TorquePerW)
	}
	for _, rc := range []float64{10e-6, 100e-6, 1e-3} {
		x := d
		x.Collector, x.RCollector = InsulatedBearings, rc
		r, _ := x.Solve()
		fmt.Printf("  insulated bearings, collector %.0e Ω: P_load %.2f W, efficiency %.3f, %.4f N·m/W\n", rc, r.PLoad, r.Efficiency, r.TorquePerW)
	}

	fmt.Println("\nCapacitor load: 1 F through R_L, field stepped (1 ms) or ramped, then 5·R_L·C to settle")
	fmt.Printf("  %8s %10s %10s %12s %14s %12s\n", "ramp s", "stored J", "shaft J", "friction J", "conversion %", "mismatch J")
	for _, ramp := range []float64{0, 10 * d.RL, 50 * d.RL} {
		tEnd := ramp + 5*d.RL*1
		st, sh, mm, err := d.Charge(1, ramp, tEnd)
		if err != nil {
			log.Fatal(err)
		}
		fr := 2 * d.FrictionTorque() * d.Omega * tEnd
		fmt.Printf("  %8.1f %10.2f %10.2f %12.2f %14.1f %12.1e\n", ramp, st, sh, fr, 100*st/(sh-fr), mm)
	}
	vc := 2*d.K()*d.Omega - 2*d.VBearing
	fmt.Printf("  the bearing drops cost 2·V_b·Q whatever the ramp: at best ½CV²/(½CV² + 2V_bCV) = %.1f %% (V = %.2f V)\n",
		100*vc/(vc+4*d.VBearing), vc)

	fmt.Printf("\nEvolving for the least drive torque per watt at ≥ %.0f W, within the bearing current limit and rim speed\n", *ptarget)
	s, err := evo.Run(context.Background(), Problem(*ptarget, d), evo.Options{Population: *pop, Generations: *gens, Seed: 1}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	v := s.Best.Values
	best := d
	best.Bore, best.Preload, best.Thickness, best.Omega, best.RL, best.NBear = v[0], v[1], v[2], v[3], v[4], int(math.Round(v[5]))
	rb := report(fmt.Sprintf("Evolved, current through the bearings (%d evaluations)", s.Evaluations), best)

	ins := d
	ins.Collector = InsulatedBearings
	s2, err := evo.Run(context.Background(), Problem(*ptarget, ins), evo.Options{Population: *pop, Generations: *gens, Seed: 1}, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	v = s2.Best.Values
	bi := ins
	bi.Bore, bi.Preload, bi.Thickness, bi.Omega, bi.RL, bi.NBear = v[0], v[1], v[2], v[3], v[4], int(math.Round(v[5]))
	ri := report(fmt.Sprintf("Evolved, bearings insulated, shaft collector %.0e Ω declared (%d evaluations)", ins.RCollector, s2.Evaluations), bi)
	fmt.Printf("\nDrive torque per load watt: %.4f N·m/W through the bearings, %.4f N·m/W insulated (%.1f× less)\n",
		rb.TorquePerW, ri.TorquePerW, rb.TorquePerW/ri.TorquePerW)
}
