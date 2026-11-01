// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command bearing_motor adds the ball-bearing (Huber-effect) motor to the
// homopolar and SPG analyses, using measurements only.
//
// A current through a rotating ball bearing drives the shaft with a torque
// that scales as I², aids the rotation and is zero at rest (device.Bearing,
// kronspice P element). In a homopolar machine whose output current passes
// through its bearings, that torque helps the drive motor — so the drive
// power measured under load would understate what the electrical output
// took from the shaft. This program asks whether that could account for
// the Space Power Generator's published figures
// (examples/space_power_generator), using:
//
//   - the measurement registry below: every number is a measured value,
//     an estimate, or a claim, with its source, and only measured values
//     are used in a calculation;
//
//   - the element's coefficients κ and ω0, which must be fitted to a
//     measured torque–current–speed table (FitTorque) or back-EMF table
//     (FitBackEMF). No such table is open access; the program says so and
//     uses no invented coefficients. The SPG check needs none: it derives
//     the bearing torque the SPG figures would require and tests it against
//     the measured efficiency and the measured back-EMF.
//
//     go run ./examples/bearing_motor
package main

import (
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Status says what kind of number a record is.
type Status int

const (
	Measured Status = iota // reported from an experiment
	Estimate               // observer's estimate, no instrument named
	Claim                  // asserted, no data shown
)

func (s Status) String() string { return [...]string{"measured", "estimate", "claim"}[s] }

// Record is one published number or finding about ball-bearing motors.
type Record struct {
	Source   string
	Quantity string
	Lo, Hi   float64 // range (Lo = Hi for a single value); NaN if qualitative
	Unit     string
	Status   Status
	Note     string
}

var nan = math.NaN()

// Registry is the public record as read (open-access text only; the
// paywalled papers' tables are listed in Needed).
var Registry = []Record{
	{"P. G. Moyssides, P. Hatzikonstantinou, IEEE Trans. Magn. 26(4), 1274–1281 (1990), DOI 10.1109/20.54012 (abstract)",
		"steady currents studied", 43.5, 70.15, "A", Measured, "self-rotates at high currents, with small efficiency; at low currents needs a driving motor"},
	{"Moyssides & Hatzikonstantinou 1990 (abstract)",
		"back-EMF", nan, nan, "", Measured, "depends on the angular velocity of shaft and inner race → EM mode"},
	{"Moyssides & Hatzikonstantinou 1990 (abstract)",
		"generator action", nan, nan, "", Measured, "none: driven at low current it does not act as a generator"},
	{"Moyssides (1992), two discs, as cited by Lauterbach, Soong & Abbott, Proc. SPIE 4236 (2001)",
		"self-rotation range", 70, 140, "A", Measured, "secondary citation; with a driving motor 18–140 A"},
	{"J. L. Choo, W. L. Soong, D. Abbott, Proc. SPIE 5649, 700–720 (2005), as reported by K. T. McDonald, Ball-Bearing Motor (Princeton, 2020)",
		"torque scaling", 2, 2, "power of I", Measured, "T ∝ I²"},
	{"S. Marinov, \"The intriguing ball-bearing motor\", Electronics & Wireless World, April 1989",
		"mechanical output / electrical input", 0.10, 0.10, "", Measured, "calorimeter; underlying data not published"},
	{"Marinov 1989",
		"heat with rotation vs without", 1.08, 1.08, "ratio", Claim, "read by Marinov as energy from nothing; no data shown; H-register item, not used here"},
	{"Marinov 1989",
		"current needed", 25, math.Inf(1), "A", Claim, "\"5 A\" and \"25 A\" not enough (in F. Ogden's account, same issue)"},
	{"F. Ogden, \"Great balls of fire!\", Electronics & Wireless World, April 1989",
		"speed", 1500, 1500, "r/min", Estimate, "\"at least 1500 rev/min\""},
	{"R. A. Milroy, J. Appl. Mech. 34, 525 (1967), as cited by Wikipedia",
		"speed", 1000, 1000, "r/min", Measured, "1959 implementation"},
	{"M. Harrison, electricstuff.co.uk (archived 2006)",
		"start from rest", nan, nan, "", Measured, "usually needs a hand spin; sense set by the spin"},
	{"K. T. McDonald, Ball-Bearing Motor (2020), fn. 6 and March 2016 experiments",
		"start-up mechanism", nan, nan, "", Measured, "thermal effects start it; Lorentz force cannot start it from rest but can contribute once running"},
	{"Univ. Adelaide project 2019s1-142 (projects wiki)",
		"currents run", 58.5, 62.8, "A", Measured, "torque increased with current; torque–speed relation inconsistent"},
	{"Univ. Adelaide ball bearing motor mystery 2022 (wiki)",
		"bearing temperature", 80, 80, "°C", Measured, "maximum, 30 A supply"},
}

// Needed lists the papers whose measured tables would fix κ and ω0.
var Needed = []string{
	"H. Gruenberg, \"The ball bearing as a motor\", Am. J. Phys. 46, 1213–1219 (1978), DOI 10.1119/1.11455",
	"P. G. Moyssides, P. Hatzikonstantinou, IEEE Trans. Magn. 26, 1274–1281 (1990): current, back-EMF and speed",
	"D. B. Watson, M. R. Williams, C. S. Crimp, \"Ball-bearing motors\", IEE Proc. A 140(4) (1993): speed vs current",
	"P. G. Moyssides, P. Hatzikonstantinou, \"Ball bearing motors\", IEEE Trans. Magn. 33, 4566–4569 (1997)",
	"J. L. Choo, W. L. Soong, D. Abbott, Proc. SPIE 5649, 700–720 (2005): acceleration torque vs speed",
	"Y. Shen et al., Proc. SPIE 3891, 178–183 (1999), DOI 10.1117/12.364436",
	"D. B. Watson, S. M. Patel, N. P. Sedcole, IEE Proc. Sci. Meas. Technol. 146(2), 83 (1999)",
	"K. M. Polivanov, N. Tatarinova, A. V. Netushil, Elektrichestvo no. 8, 72–76 (1973)",
	"J. Huber, Elektrotechnik und Maschinenbau 76(8), 169–174 (1959)",
}

// TorquePoint is one measured (I, ω, T) point.
type TorquePoint struct{ I, Omega, T float64 }

// EMFPoint is one measured (I, ω, back-EMF) point.
type EMFPoint struct{ I, Omega, E float64 }

// Dataset is a measured table with its source.
type Dataset struct {
	Source string
	Torque []TorquePoint
	EMF    []EMFPoint
}

// Datasets holds the measured tables transcribed from the papers in Needed.
// Empty until a paper's table is entered from its page.
var Datasets []Dataset

// fit1 finds ω0 on a log grid (refined by golden section) and the
// least-squares κ for model y_k = κ·a_k(ω0).
func fit1(n int, a func(k int, w0 float64) float64, y func(k int) float64) (kappa, w0, rms float64) {
	cost := func(w0 float64) (float64, float64) {
		var sa, sy float64
		for k := 0; k < n; k++ {
			sa += a(k, w0) * a(k, w0)
			sy += a(k, w0) * y(k)
		}
		if sa == 0 {
			return 0, math.Inf(1)
		}
		kap := sy / sa
		var r float64
		for k := 0; k < n; k++ {
			d := kap*a(k, w0) - y(k)
			r += d * d
		}
		return kap, r
	}
	best, bw := math.Inf(1), 0.0
	for e := -3.0; e <= 5; e += 0.05 {
		if _, r := cost(math.Pow(10, e)); r < best {
			best, bw = r, math.Pow(10, e)
		}
	}
	lo, hi := math.Log(bw)-0.06*math.Ln10, math.Log(bw)+0.06*math.Ln10
	for i := 0; i < 80; i++ {
		m1, m2 := lo+(hi-lo)*0.382, lo+(hi-lo)*0.618
		_, r1 := cost(math.Exp(m1))
		_, r2 := cost(math.Exp(m2))
		if r1 < r2 {
			hi = m2
		} else {
			lo = m1
		}
	}
	w0 = math.Exp(0.5 * (lo + hi))
	kappa, r := cost(w0)
	return kappa, w0, math.Sqrt(r / float64(n))
}

// FitTorque fits T = κ·I²·tanh(ω/ω0) (EM-mode torque law) to measured points.
func FitTorque(p []TorquePoint) (kappa, w0, rms float64, err error) {
	if len(p) < 3 {
		return 0, 0, 0, fmt.Errorf("need ≥ 3 measured points, have %d", len(p))
	}
	k, w, r := fit1(len(p), func(k int, w0 float64) float64 { return p[k].I * p[k].I * math.Tanh(p[k].Omega/w0) },
		func(k int) float64 { return p[k].T })
	return k, w, r, nil
}

// FitBackEMF fits E = κ·tanh(ω/ω0)·ω·I to measured back-EMF points. In EM
// mode it must give the same κ, ω0 as FitTorque on the same bearing
// (energy: E·I = ω·T); a measured back-EMF with no matching torque, or the
// reverse, rules the EM mode out for that bearing.
func FitBackEMF(p []EMFPoint) (kappa, w0, rms float64, err error) {
	if len(p) < 3 {
		return 0, 0, 0, fmt.Errorf("need ≥ 3 measured points, have %d", len(p))
	}
	k, w, r := fit1(len(p), func(k int, w0 float64) float64 { return p[k].I * p[k].Omega * math.Tanh(p[k].Omega/w0) },
		func(k int) float64 { return p[k].E })
	return k, w, r, nil
}

// SPGCheck asks what bearing torque the SPG figures would need.
type SPGCheck struct {
	Omega, I, E0        float64 // measured SPG speed [rad/s], current [A], no-load EMF [V]
	TClaimed, TMeasured float64 // k·I implied by the claim; DM shaft-torque increase [N·m]
	TBearing            float64 // aiding torque needed to bridge them [N·m]
	PBearing            float64 // its shaft power [W]
	KappaPsi            float64 // κ·ψ(ω) it needs if all of I crosses the bearings [N·m/A²]
	BackEMF             float64 // EM mode: voltage the bearings drop for it [V]
	PLoadEM             float64 // EM mode: power left for the rest of the loop [W]
	PElecThermal        float64 // thermal mode at the measured η: electrical input to the bearings [W]
	VThermal            float64 // and the voltage across them [V]
	EtaMeasured         float64
}

// CheckSPG uses the SPG test (Tewari, UPSM pp. 159–160) and the measured
// bearing-motor efficiency.
func CheckSPG(eta float64) SPGCheck {
	const (
		rpm, w0, wl, etaDM = 2700, 2943, 8620, 0.83 // UPSM p. 159–160
		e0, amps           = 3.41, 3438.0
	)
	var c SPGCheck
	c.Omega, c.I, c.E0, c.EtaMeasured = rpm*2*math.Pi/60, amps, e0, eta
	c.TClaimed = e0 * amps / c.Omega
	c.TMeasured = etaDM * (wl - w0) / c.Omega
	c.TBearing = c.TClaimed - c.TMeasured
	c.PBearing = c.TBearing * c.Omega
	c.KappaPsi = c.TBearing / (amps * amps)
	c.BackEMF = c.PBearing / amps
	c.PLoadEM = (e0 - c.BackEMF) * amps
	c.PElecThermal = c.PBearing / eta
	c.VThermal = c.PElecThermal / amps
	return c
}

// SimulateSPGLoop runs the EM-mode case through the engine: the no-load EMF
// as a source (the claim's own premise), the bearing in the loop at the
// SPG's speed with the κψ CheckSPG derived, and the rest of the loop sized
// so 3438 A flows. Rb is the bearing's contact resistance, which no source
// gives, so it is swept: the result must not depend on it.
func SimulateSPGLoop(c SPGCheck, rb float64) (pLoad, pShaft, residual float64, err error) {
	w0 := c.Omega // any ω0: κ is scaled so κ·tanh(ω/ω0) = KappaPsi
	kappa := c.KappaPsi / math.Tanh(c.Omega/w0)
	rRest := c.E0/c.I - rb - c.KappaPsi*c.Omega
	if rRest <= 0 {
		return 0, 0, 0, fmt.Errorf("Rb = %g Ω leaves no loop resistance", rb)
	}
	ckt, err := circuit.NewBuilder("SPG loop with bearing motor").Domain(exchange.Mechanical, "w").
		VoltageSource("E0", "e", "0", c.E0).
		BearingMotor("PB", "e", "x", "w", "0", rb, kappa, w0, false, 0).
		Resistor("RLOOP", "x", "0", rRest).
		VoltageSource("VW", "w", "0", c.Omega).Build()
	if err != nil {
		return 0, 0, 0, err
	}
	r, err := analysis.OP(ckt, analysis.Options{})
	if err != nil {
		return 0, 0, 0, err
	}
	vx := r.NodeVoltages[ckt.NodeIndex("x")]
	pLoad = vx * vx / rRest
	for _, e := range r.Elements {
		if e.Name == "VW" {
			pShaft = e.PowerAbsorbed // shaft power the speed source absorbs = bearing torque × ω
		}
	}
	return pLoad, pShaft, r.Ledger.RelativeResidual(), nil
}

func main() {
	fmt.Println("Ball-bearing (Huber) motor — measurement registry")
	for _, r := range Registry {
		v := ""
		switch {
		case math.IsNaN(r.Lo):
		case r.Lo == r.Hi:
			v = fmt.Sprintf("%g %s", r.Lo, r.Unit)
		case math.IsInf(r.Hi, 1):
			v = fmt.Sprintf("> %g %s", r.Lo, r.Unit)
		default:
			v = fmt.Sprintf("%g–%g %s", r.Lo, r.Hi, r.Unit)
		}
		fmt.Printf("  [%-8s] %-38s %-18s %s\n             %s\n", r.Status, r.Quantity, v, r.Note, r.Source)
	}
	fmt.Println()
	if len(Datasets) == 0 {
		fmt.Println("κ, ω0: no measured torque or back-EMF table is open access. Needed (any one fixes κ and ω0):")
		for _, n := range Needed {
			fmt.Println("  -", n)
		}
		fmt.Println("  The element refuses to run without fitted KAPPA and W0; no coefficients are invented here.")
	}
	for _, d := range Datasets {
		if k, w, r, err := FitTorque(d.Torque); err == nil {
			fmt.Printf("%s: torque fit κ = %.4g N·m/A², ω0 = %.4g rad/s, rms %.3g N·m\n", d.Source, k, w, r)
		}
		if k, w, r, err := FitBackEMF(d.EMF); err == nil {
			fmt.Printf("%s: back-EMF fit κ = %.4g N·m/A², ω0 = %.4g rad/s, rms %.3g V\n", d.Source, k, w, r)
		}
	}

	eta := 0.0
	for _, r := range Registry {
		if r.Status == Measured && r.Quantity == "mechanical output / electrical input" {
			eta = r.Lo
		}
	}
	c := CheckSPG(eta)
	fmt.Println()
	fmt.Println("Could bearing torque explain the SPG figures? (all inputs measured)")
	fmt.Printf("  SPG at %.1f rad/s, %.0f A: claimed E0·I needs %.1f N·m on the shaft; drive motor supplied %.1f N·m\n", c.Omega, c.I, c.TClaimed, c.TMeasured)
	fmt.Printf("  gap %.1f N·m (%.2f kW) would have to come from bearing torque aiding the rotation\n", c.TBearing, c.PBearing/1e3)
	fmt.Printf("  → κ·ψ(ω) ≥ %.3g N·m/A² if the whole %.0f A crosses the bearings\n", c.KappaPsi, c.I)
	fmt.Println()
	fmt.Println("  EM mode (measured back-EMF, Moyssides & Hatzikonstantinou 1990):")
	fmt.Printf("    the bearings drop %.2f V of the %.2f V EMF; the rest of the loop gets %.2f kW, not %.2f kW\n", c.BackEMF, c.E0, c.PLoadEM/1e3, c.E0*c.I/1e3)
	fmt.Printf("    (= the %.2f V loaded-EMF bound of the SPG ledger: the bearing transfers shaft work, it adds none)\n", c.E0-c.BackEMF)
	fmt.Printf("  thermal mode at the measured η = %.2f (Marinov 1989):\n", eta)
	fmt.Printf("    the bearings would need %.0f kW electrical, %.1f V across them — more than the %.2f V the machine makes: excluded\n", c.PElecThermal/1e3, c.VThermal, c.E0)
	fmt.Println()
	fmt.Println("  Engine run (EM mode, bearing contact resistance unknown → swept):")
	for _, rb := range []float64{1e-6, 1e-5, 1e-4} {
		pl, ps, res, err := SimulateSPGLoop(c, rb)
		if err != nil {
			fmt.Printf("    Rb = %g Ω: %v\n", rb, err)
			continue
		}
		fmt.Printf("    Rb = %-6g Ω: load %.2f kW, bearing shaft work %.2f kW, E0·I %.2f kW, ledger residual %.1e\n", rb, pl/1e3, ps/1e3, c.E0*c.I/1e3, res)
	}
	fmt.Println()
	fmt.Println(strings.TrimSpace(`
  Decisive bench tests: (1) insulate the SPG bearings (ceramic balls or insulated housings) and repeat the loaded
  run — if bearing torque mattered, the drive-motor power under load rises toward k·I·ω; (2) measure the voltage
  across each bearing under load; (3) a shaft torque transducer between drive motor and rotor.`))
}
