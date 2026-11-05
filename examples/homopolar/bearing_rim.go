// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Ball bearings at the disks' rims as the current collectors of the series
// dual-motor machine, and an adjustable positive feedback that sets the
// current build-up.
//
// Each disk's rim current passes through a ball-bearing bank (the kronspice
// P element, device/bearing.go) on its way to R_int and the field coil; the
// bank's inner race turns with the disk, so its Huber torque acts on that
// disk's shaft. Measured inputs, as in examples/homopolar_bearing:
//
//   - DC bearing voltage: almost constant, 0.73 V at 4 A, 0.5–1.5 V
//     (O. Safdarzadeh, A. Farahi, A. Binder et al., Lubricants 12(5), 148,
//     2024; greased 51208, d = 40 mm).
//   - Damage limit: 0.1–0.3 A/mm² of Hertzian contact area; A = 1.0 mm² at
//     400 N on that bearing (same paper), scaled by Hertz theory,
//     A ∝ (F·d)^(2/3) (K. L. Johnson, Contact Mechanics, §4.2). A rim
//     bearing of the disk's 300 mm diameter is 7.5× the measured bore: the
//     scaling is an extrapolation, labelled as such.
//   - The Huber coefficient κψ has never been measured: it is swept, in EM
//     mode (its shaft work taken from the loop as a back-EMF) and in thermal
//     mode (at most η = 10 % of the contact heat, declared as in
//     examples/homopolar_bearing).
//
// The adjustable feedback: K_v·M of own-disk coupling added to the stacked
// self-coupling K_S (a feedback winding moved or tapped). The bearing drop
// is a constant voltage, while the disks' EMF is proportional to the current,
// so at any fixed K_v the loop has one equilibrium current
//
//	I* = 2·V_b / (2·M_s·ω − R_loop)
//
// and it is unstable: above I* the current grows, below it collapses. A
// regulator that trims K_v holds the current at a set point, typically at
// the bearings' damage limit; the iron's saturation would also stop a runaway,
// but the model is linear, so the field is checked against 1.6 T instead.

package main

import (
	"fmt"
	"math"
)

// Measured and declared bearing inputs (see the file comment).
const (
	rimJLimit = 0.3e6  // A/m², upper damage threshold
	rimAHzRef = 1.0e-6 // m² at rimFRef on a rimDRef bore
	rimFRef   = 400.0  // N
	rimDRef   = 0.040  // m
	rimVolt   = 0.73   // V, DC bearing voltage
	rimEta    = 0.10   // declared thermal-mode efficiency ceiling
)

// RimBearing is one disk's rim collector: Races ball bearings in parallel
// around the rim, each preloaded to Preload.
type RimBearing struct {
	Races    int
	Bore     float64 // m: the rim's diameter
	Preload  float64 // N per race
	VF       float64 // V: measured DC bearing voltage
	R        float64 // Ω per race (declared; shares current between races)
	KappaPsi float64 // N·m/A² per race at the running speed (unmeasured: swept)
	W0       float64 // rad/s: ψ = tanh(ω/ω0); set to the running speed
	Thermal  bool
	Eta      float64
}

// defaultRim: four 300 mm races per rim at 5 kN each (declared).
func defaultRim(w0 float64) *RimBearing {
	return &RimBearing{Races: 4, Bore: 0.30, Preload: 5e3, VF: rimVolt, R: 1e-3, W0: w0, Eta: rimEta}
}

// HertzArea is one race's apparent contact area, scaled from the measured
// 51208 point (an extrapolation for bores far from 40 mm).
func (r *RimBearing) HertzArea() float64 {
	return rimAHzRef * math.Pow(r.Preload/rimFRef*r.Bore/rimDRef, 2.0/3)
}

// ILimit is the bank's current at the 0.3 A/mm² damage threshold.
func (r *RimBearing) ILimit() float64 { return rimJLimit * r.HertzArea() * float64(r.Races) }

// bank returns the bank's Huber κ and resistance: Races in parallel share
// the current, so the torque Races·κ·(I/Races)² = (κ/Races)·I².
func (r *RimBearing) bank() (kappa, res float64) {
	k := r.KappaPsi / math.Tanh(1)
	return k / float64(r.Races), r.R / float64(r.Races)
}

// drop is one bank's voltage and Huber torque at current i and speed w
// (the P element's equations, circuit/bearing.go).
func (r *RimBearing) drop(i, w float64) (v, t float64) {
	k, res := r.bank()
	psi := math.Tanh(w / r.W0)
	h := r.VF*math.Abs(i) + res*i*i // contact heat
	v = r.VF*math.Tanh(i/1e-3) + res*i
	a := k * i * i * psi
	if r.Thermal {
		if h > 0 {
			t = a / (1 + a*w/(r.Eta*h))
		}
		return
	}
	return v + k*psi*w*i, a
}

// loop returns both banks' drop round the series loop (current i along it)
// and the Huber torques on shafts 1 and 2; a nil bank contributes nothing.
func (r *RimBearing) loop(i, w1, w2 float64) (v, t1, t2 float64) {
	if r == nil {
		return 0, 0, 0
	}
	v1, t1 := r.drop(i, w1)
	v2, t2 := r.drop(i, w2) // bank 2 is crossed rim to shaft side backwards; its drop is odd in i, so it adds
	return v1 + v2, t1, t2
}

// Regulator trims K_v to hold the loop current at Target, a PI law on the
// error e = 1 − |I|/Target: K_v = x + Kp·e within [Kmin, Kmax], with
// dx/dt = e/Tau (held while K_v sits at a limit). The proportional part is
// required: at fixed K_v the loop is unstable (rate 2·V_b/(L·I)), and an
// integral law alone cannot hold an unstable loop. Kp must exceed
// 2·V_b/(c·I), c·I the disks' EMF per unit of (1 + K).
type Regulator struct {
	Target, Kmin, Kmax, Kp, Tau float64
}

// k is the coupling the regulator sets from the current and its integral x.
func (g *Regulator) k(i, x float64) float64 {
	return math.Max(g.Kmin, math.Min(g.Kmax, x+g.Kp*(1-math.Abs(i)/g.Target)))
}

func (g *Regulator) rate(i, x float64) float64 {
	e := 1 - math.Abs(i)/g.Target
	if k := x + g.Kp*e; (k >= g.Kmax && e > 0) || (k <= g.Kmin && e < 0) {
		return 0
	}
	return e / g.Tau
}

// fieldB is the mean axial field on a disk at loop current i, from
// EMF = M_s·ω·I = ½·B·ω·(r_o² − r_i²).
func (m Dual) fieldB(i float64) float64 {
	_, _, ks, _ := m.loop()
	return 2 * (1 + ks) * m.Coupling() * math.Abs(i) / (m.Disk.Ro*m.Disk.Ro - m.Disk.Ri*m.Disk.Ri)
}

// equilibrium is the series loop's one current at fixed K_v (held speeds
// ±ω, no sine): I* = 2·V_b / (2·M_s·ω − R_DC), R_DC the loop's DC resistance.
func (m Dual) equilibrium(w float64) float64 {
	_, _, ks, r := m.loop()
	z, _ := m.TankScale()
	rdc := 2*(r-m.RLoad) + z*m.RPrim // the primary is the tank's DC path when present
	if m.LPrim == 0 {
		rdc = 2 * r // the load is in the DC path
	}
	_, rb := m.Rim.bank()
	rdc += 2 * rb
	return 2 * m.Rim.VF / (2*(1+ks)*m.Coupling()*w - rdc)
}

// superMachine is the bearing-matched design: superconducting field coils of
// N turns chosen so the disk field is 1.6 T at the set current iSet, and a
// series load R_L chosen so that iSet is the equilibrium at K_v = 0.
func superMachine(rim *RimBearing, iSet, w, gain float64) Dual {
	m := DefaultDual()
	m.Stacked, m.Windage, m.IronGain, m.SuperField = true, 5e-3, gain, true
	m.Output, m.Rim = OutSeries, rim
	m.RLoad = 1 // placeholder for loop()
	m.N *= bSatIron / m.fieldB(iSet)
	_, _, ks, r := m.loop()
	_, rb := rim.bank()
	rint := r - m.RLoad
	rl := 2*(1+ks)*m.Coupling()*w - 2*rint - 2*rb - 2*rim.VF/iSet // total series load
	m.RLoad = rl / 2                                              // the series tank doubles the per-loop value
	m.Speed = m.held(w)
	return m
}

// energy integrates a series RK4 trace's books from t0: drive work (the
// disks' generated power), load, bearings, the rest of the loop, and the
// change of stored magnetic energy; it returns them and the residual.
func (m Dual) energy(tr DualTrace, t0 float64) (gen, load, bear, other, dW, resid float64) {
	l11, l12, ks0, r := m.loop()
	M := m.Coupling()
	lt := 2 * (l11 + l12)
	z, _ := m.TankScale()
	k0 := -1
	for k := 1; k < len(tr.T); k++ {
		if tr.T[k] < t0 {
			continue
		}
		if k0 < 0 {
			k0 = k - 1
		}
		dt := tr.T[k] - tr.T[k-1]
		p := func(j int) (g, l, b, o float64) {
			i, ks := tr.I1[j], ks0
			if tr.K != nil {
				ks = ks0 - m.FeedbackK + tr.K[j]
			}
			vb, _, _ := m.Rim.loop(i, tr.W1[j], tr.W2[j])
			return M * (1 + ks) * (tr.W1[j] - tr.W2[j]) * i * i, z * m.RLoad * i * i, vb * i, 2 * (r - m.RLoad) * i * i
		}
		g0, l0, b0, o0 := p(k - 1)
		g1, l1, b1, o1 := p(k)
		gen += dt / 2 * (g0 + g1)
		load += dt / 2 * (l0 + l1)
		bear += dt / 2 * (b0 + b1)
		other += dt / 2 * (o0 + o1)
	}
	n := len(tr.T) - 1
	dW = 0.5 * lt * (tr.I1[n]*tr.I1[n] - tr.I1[k0]*tr.I1[k0])
	resid = gen - load - bear - other - dW
	return
}

func runRimBearing() error {
	w0 := rpmBase * 2 * math.Pi / 60
	fmt.Println("  I. ball bearings at the rims as the current collectors, and an adjustable positive feedback")
	rim := defaultRim(w0)
	ilim := rim.ILimit()
	fmt.Printf("     rim bank per disk: %d × %.0f mm races at %.0f N (declared), V_b = %.2f V (measured, Safdarzadeh 2024)\n",
		rim.Races, rim.Bore*1e3, rim.Preload, rim.VF)
	fmt.Printf("     Hertz area %.1f mm² per race (scaled from the measured 40 mm bearing: a 7.5× extrapolation);\n", rim.HertzArea()*1e6)
	fmt.Printf("     damage limit %.1f A per bank at 0.3 A/mm² (%.1f A at 0.1 A/mm²); speed factor d·n = %.0f mm·rpm at 1768 rpm\n",
		ilim, ilim/3, rim.Bore*1e3*rpmBase)

	// 1. The present copper machine (part G design) with rim bearings.
	m, _ := outputMachine(10e-6, 3, 10, 0.98)
	m.Rim = rim
	fmt.Println("     1. the part G machine (copper coils, 4 turns) with its rims on bearings:")
	for _, x := range []float64{1.0, 1.5, 2.0, 2.5} {
		q := m
		q.Speed = q.held(x * w0)
		if is := q.equilibrium(x * w0); is > 0 {
			fmt.Printf("        at %.1f× the onset speed the loop holds only above I* = %.1f kA\n", x, is/1e3)
		} else {
			fmt.Printf("        at %.1f× the onset speed: no equilibrium (the field cannot build at any current)\n", x)
		}
	}
	fmt.Printf("        at the 10 kA design current: bearings %.0f kW of heat, %.0f× their damage limit\n",
		2*rim.VF*10e3/1e3, 10e3/ilim)
	{
		q := m
		w := 2.5 * w0
		q.Speed = q.held(w)
		is := q.equilibrium(w)
		for _, f := range []float64{0.9, 1.1} {
			q.I0 = [2]float64{f * is, f * is}
			ct, err := q.Run(0.5, 2e-4)
			if err != nil {
				return err
			}
			rk := q.RunRK4(0.5, 1e-4)
			n := len(ct.T) - 1
			fmt.Printf("        engine at 2.5× onset (I* = %.1f kA) from %.2f I*: I(0.5 s)/I(0) = %.4f (RK4 %.4f), ledger mismatch %.0e\n",
				is/1e3, f, ct.I1[n]/ct.I1[0], rk.I1[len(rk.I1)-1]/rk.I1[0], ct.Mismatch/ct.In)
		}
	}
	_, rc := DefaultDual().Coil()
	l0 := rpmMachine(1)
	l0.IronGain = m.IronGain
	isat := bSatIron / l0.fieldB(1)
	fmt.Printf("        copper field: 1.6 T needs %.0f kA·turns whatever the turn count; in this window that is %.0f kW of\n",
		4*isat/1e3, isat*isat*rc/1e3)
	fmt.Printf("        coil heat per coil, against %.0f W generated per disk at the bearing limit: no copper-coil design\n",
		0.5*bSatIron*w0*(m.Disk.Ro*m.Disk.Ro-m.Disk.Ri*m.Disk.Ri)*ilim)
	fmt.Println("        can excite itself at a current a bearing carries (the field costs more than the disks make).")

	// 2. The bearing-matched design: superconducting field coils.
	iSet := 0.8 * ilim
	s := superMachine(rim, iSet, w0, m.IronGain)
	l11, l12, ks, rs := s.loop()
	fmt.Println("     2. bearing-matched: superconducting field coils (declared; cryocooler power not modeled), set current")
	fmt.Printf("        %.1f A = 0.8 × the damage limit: %.0f turns per coil for 1.6 T, loop L = %.0f H, series load %.3f Ω\n",
		iSet, s.N, 2*(l11+l12), 2*s.RLoad)
	fmt.Printf("        EMF %.2f V per disk, bearings %.2f V each; equilibrium at K_v = 0: %.2f A (set %.2f A)\n",
		(1+ks)*s.Coupling()*w0*iSet, rim.VF, s.equilibrium(w0), iSet)
	fmt.Println("        fixed K_v = 0, started 5 % off the equilibrium (engine, 600 s; RK4 beside it):")
	for _, f := range []float64{0.95, 1.05} {
		q := s
		q.I0 = [2]float64{f * iSet, f * iSet}
		ct, err := q.Run(600, 0.5)
		if err != nil {
			return err
		}
		rk := q.RunRK4(600, 0.05)
		n := len(ct.T) - 1
		fmt.Printf("          from %.2f I*: I(600 s) = %.3f A (RK4 %.3f A), ledger mismatch %.0e of the drive work\n",
			f, ct.I1[n], rk.I1[len(rk.I1)-1], ct.Mismatch/ct.In)
	}
	fmt.Println("          the equilibrium is unstable: the current runs away from it either way.")

	// The regulator: start pulse, build-up, hold.
	reg := &Regulator{Target: iSet, Kmin: -0.5, Kmax: 1.0, Kp: 2, Tau: 100}
	q := s
	q.Reg = reg
	q.FeedbackK = reg.Kmax
	iStart := q.equilibrium(w0) // the least start current with K_v at its top
	q.FeedbackK = 0
	q.I0 = [2]float64{1.5 * iStart, 1.5 * iStart}
	q.FeedbackK = 0 // the integral starts at 0; the proportional part opens the feedback fully at the start
	tr := q.RunRK4(3000, 0.05)
	n := len(tr.T) - 1
	t95 := math.NaN()
	for k := range tr.T {
		if tr.I1[k] >= 0.95*iSet {
			t95 = tr.T[k]
			break
		}
	}
	gen, load, bear, other, dW, res := q.energy(tr, 0)
	fmt.Printf("        variable feedback K_v ∈ [%.1f, %.1f], PI regulator (K_p = %.0f, τ = %.0f s): start pulse %.2f A (1.5× the least,\n",
		reg.Kmin, reg.Kmax, reg.Kp, reg.Tau, q.I0[0])
	fmt.Printf("        %.2f A, at full feedback); 95 %% of the set current after %.0f s; at 3000 s I = %.3f A, K_v = %+.3f,\n",
		iStart, t95, tr.I1[n], tr.K[n])
	q.FeedbackK = tr.K[n]
	fmt.Printf("        field %.2f T; load %.1f W, bearings %.1f W, the rest of the loop %.3f W (drive %.1f W)\n",
		q.fieldB(tr.I1[n]), 2*s.RLoad*tr.I1[n]*tr.I1[n], 2*rim.VF*tr.I1[n], 2*(rs-s.RLoad)*tr.I1[n]*tr.I1[n], 2*(1+ks-s.FeedbackK+tr.K[n])*s.Coupling()*w0*tr.I1[n]*tr.I1[n])
	fmt.Printf("        books over the run: drive %.0f J = load %.0f + bearings %.0f + loop %.2g + stored %.0f J (residual %.0e)\n",
		gen, load, bear, other, dW, res/gen)

	// Huber torque swept at the regulated point.
	fmt.Println("        Huber torque at the set current, κψ per race swept (unmeasured), regulated steady state at 3000 s:")
	fmt.Println("          mode      κψ [N·m/A²]   K_v      I [A]    load [W]   Huber torque per shaft [N·m]   drive torque [N·m]")
	for _, th := range []bool{false, true} {
		for _, kp := range []float64{0, 1e-4, 1e-3, 1e-2} {
			r2 := *rim
			r2.KappaPsi, r2.Thermal = kp, th
			p := s
			p.Rim, p.Reg, p.FeedbackK = &r2, reg, 0
			p.I0 = [2]float64{iSet, iSet}
			tq := p.RunRK4(3000, 0.05)
			j := len(tq.T) - 1
			i := tq.I1[j]
			_, t1, _ := r2.loop(i, w0, -w0)
			_, _, ksr, _ := p.loop()
			te := (1 + ksr - p.FeedbackK + tq.K[j]) * p.Coupling() * i * i
			mode := "EM"
			if th {
				mode = "thermal"
			}
			fmt.Printf("          %-8s  %-11.0e   %+.3f   %6.3f   %7.1f    %8.4f                       %8.4f\n",
				mode, kp, tq.K[j], i, 2*p.RLoad*i*i, t1, te-t1)
		}
	}
	{
		r2 := *rim
		r2.KappaPsi = 1e-3
		p := s
		p.Rim = &r2
		p.I0 = [2]float64{iSet, iSet}
		ct, err := p.Run(100, 0.5)
		if err != nil {
			return err
		}
		fmt.Printf("        engine with the P elements at κψ = 1e-3 (EM), fixed K_v = 0, 100 s: ledger mismatch %.0e of the drive work\n",
			ct.Mismatch/ct.In)
	}
	fmt.Println("        EM mode: at the regulated current the drives pay the same torque whatever κψ is; the Huber torque on")
	fmt.Println("        each shaft equals the extra disk torque the feedback must add to cover its back-EMF. Past what the")
	fmt.Println("        feedback can make up (K_v at its top) the current collapses. Thermal mode returns at most 10 % of the")
	fmt.Println("        bearing heat as shaft work. Either way every joule comes from the two drives.")
	fmt.Println()
	return nil
}
