// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Layer 2 refinements: the global circuit with
//
//   - a two-scale conductivity profile (a boundary layer under the free
//     troposphere, continuous at the top of the boundary layer);
//   - land and ocean as separate columns in parallel between the ground and
//     the ionosphere;
//   - storm generators with their own resistances: Wilson's generator, a
//     charge-separating current between the cloud's lower and upper charge
//     centres, with a path up to the ionosphere and a path down to the
//     ground (C. T. R. Wilson 1920; the equivalent circuit of Rycroft,
//     Israelsson & Price 2000, Fig. 2).
//
// The parameters marked "declared" are illustrative, not measured; the
// three global anchors of StandardAtmosphere still fix the calibration.

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Column is one fair-weather column of the global circuit.
type Column struct {
	Name  string
	Frac  float64               // fraction of the Earth's surface
	Sigma func(float64) float64 // conductivity σ(z) [S/m]
}

// Storms is the generator side: Wilson's generator current and its two
// paths (declared).
type Storms struct {
	Is          float64 // charge-separation current, all storms together [A]
	RUp, RDown  float64 // cloud top → ionosphere, cloud base → ground (all storms in parallel) [Ω]
	RCloud      float64 // local return between the charge centres: intracloud lightning, conduction [Ω]
	Distributed bool    // false: ideal source straight to the ionosphere (RUp = RDown = ignored)
}

// columnShells discretizes one column with the same height grid as Shells,
// integrating ∫dz/σ numerically (Simpson, 64 panels per shell).
func (a Atmosphere) columnShells(col Column) []Shell {
	const stretch = 8.0
	z := func(k int) float64 {
		return a.Top * (math.Exp(stretch*float64(k)/float64(a.Layers)) - 1) / (math.Exp(stretch) - 1)
	}
	out := make([]Shell, a.Layers)
	for k := range out {
		z0, z1 := z(k), z(k+1)
		const m = 64
		hstep := (z1 - z0) / m
		var integ float64
		for i := 0; i <= m; i++ {
			w := 2.0
			if i == 0 || i == m {
				w = 1
			} else if i%2 == 1 {
				w = 4
			}
			integ += w / col.Sigma(z0+float64(i)*hstep)
		}
		integ *= hstep / 3
		r0, r1 := rEarth+z0, rEarth+z1
		rm := 0.5 * (r0 + r1)
		area := col.Frac * 4 * math.Pi * rm * rm
		out[k] = Shell{Z0: z0, Z1: z1, R: integ / area, C: col.Frac * 4 * math.Pi * eps0 * r0 * r1 / (r1 - r0)}
	}
	return out
}

// GlobalCircuit builds the columns in parallel between the ground and
// "ion", and the storms.
func (a Atmosphere) GlobalCircuit(cols []Column, st Storms) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("Global circuit, columns and storms")
	for _, col := range cols {
		for k, s := range a.columnShells(col) {
			lo, hi := col.Name+fmt.Sprint(k), col.Name+fmt.Sprint(k+1)
			if k == 0 {
				lo = "0"
			}
			if k == a.Layers-1 {
				hi = "ion"
			}
			b.Resistor(fmt.Sprintf("R%s%d", col.Name, k), hi, lo, s.R)
			b.Capacitor(fmt.Sprintf("C%s%d", col.Name, k), hi, lo, s.C)
		}
	}
	src := circuit.Element{Name: "ISTORM", Kind: circuit.CurrentSource, Value: st.Is, HasDC: true, Channel: exchange.Thermal}
	if !st.Distributed {
		b.Add(src, "0", "ion")
	} else {
		// lower charge centre "cb" → (generator) → upper charge centre "ct"
		b.Add(src, "cb", "ct")
		b.Resistor("RUP", "ct", "ion", st.RUp)
		b.Resistor("RDOWN", "cb", "0", st.RDown)
		b.Resistor("RCLOUD", "ct", "cb", st.RCloud)
	}
	return b.Build()
}

// GlobalState is the operating point of a refined circuit.
type GlobalState struct {
	VI      float64
	IFair   []float64 // fair-weather current per column [A]
	E0      []float64 // surface field per column [V/m]
	IUp     float64   // generator current reaching the ionosphere [A]
	PGen    float64   // generator power [W]
	PStorm  float64   // dissipated inside the storm paths [W]
	PFair   float64   // dissipated in fair weather [W]
	ZHalf   []float64 // height below which half of each column's space charge lies
	Ledger  float64
	columns []Column
}

// SolveGlobal runs the operating point.
func (a Atmosphere) SolveGlobal(cols []Column, st Storms) (GlobalState, error) {
	g := GlobalState{columns: cols}
	c, err := a.GlobalCircuit(cols, st)
	if err != nil {
		return g, err
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return g, err
	}
	g.VI, _ = op.Voltage("ion")
	for _, col := range cols {
		sh := a.columnShells(col)
		v := make([]float64, len(sh))
		for k := range sh {
			e, _ := op.Element(fmt.Sprintf("R%s%d", col.Name, k))
			v[k] = e.Current * sh[k].R
		}
		i := v[0] / sh[0].R
		g.IFair = append(g.IFair, math.Abs(i))
		g.E0 = append(g.E0, math.Abs(v[0])/(sh[0].Z1-sh[0].Z0))
		var tot, cum float64
		n := len(sh) - 1
		tot = math.Abs(sh[0].C*v[0]) - math.Abs(sh[n].C*v[n])
		zh := 0.0
		for k := 1; k <= n; k++ {
			cum += math.Abs(sh[k-1].C*v[k-1]) - math.Abs(sh[k].C*v[k])
			if cum >= tot/2 {
				zh = sh[k].Z0
				break
			}
		}
		g.ZHalf = append(g.ZHalf, zh)
		for k := range sh {
			e, _ := op.Element(fmt.Sprintf("R%s%d", col.Name, k))
			g.PFair += e.PowerAbsorbed
		}
	}
	if st.Distributed {
		up, _ := op.Element("RUP")
		dn, _ := op.Element("RDOWN")
		cl, _ := op.Element("RCLOUD")
		g.IUp = math.Abs(up.Current)
		g.PStorm = up.PowerAbsorbed + dn.PowerAbsorbed + cl.PowerAbsorbed
	} else {
		g.IUp = st.Is
	}
	gen, _ := op.Element("ISTORM")
	g.PGen = -gen.PowerAbsorbed
	if op.Ledger != nil && op.Ledger.Scale > 0 {
		g.Ledger = op.Ledger.Residual / op.Ledger.Scale
	}
	return g, nil
}

// TwoScale returns a continuous two-scale profile: σ0·exp(z/hBL) below zBL,
// then exp((z − zBL)/hFree) above, with hBL solved so that the column has
// the anchors' V_I at the anchors' J (σ0 = J/E0 as before).
func (a Atmosphere) TwoScale(zBL, hFree float64) (func(float64) float64, float64) {
	s0 := a.J / a.E0
	prof := func(hBL float64) func(float64) float64 {
		return func(z float64) float64 {
			if z < zBL {
				return s0 * math.Exp(z/hBL)
			}
			return s0 * math.Exp(zBL/hBL) * math.Exp((z-zBL)/hFree)
		}
	}
	colV := func(hBL float64) float64 { // flat column: J·∫dz/σ
		s := prof(hBL)
		const m = 20000
		var integ float64
		dz := a.Top / m
		for i := 0; i < m; i++ {
			integ += dz / s((float64(i)+0.5)*dz)
		}
		return a.J * integ
	}
	lo, hi := 50.0, 1e5 // bisection on hBL: V decreases as hBL falls
	for it := 0; it < 100; it++ {
		mid := math.Sqrt(lo * hi)
		if colV(mid) > a.VI {
			hi = mid
		} else {
			lo = mid
		}
	}
	h := math.Sqrt(lo * hi)
	return prof(h), h
}

func runAtmosphereRefined() error {
	a := StandardAtmosphere()
	s0, H := a.Calibrate()
	single := []Column{{"g", 1, func(z float64) float64 { return s0 * math.Exp(z/H) }}}
	fmt.Println("Layer 2 refinements")
	base, err := a.SolveGlobal(single, Storms{Is: a.GlobalCurrent()})
	if err != nil {
		return err
	}
	fmt.Printf("  one exponential column (as above): V_I %.2f kV, half the space charge below %.2f km\n", base.VI/1e3, base.ZHalf[0]/1e3)

	const zBL, hFree = 1500.0, 5000.0 // declared: boundary layer 1.5 km, free-troposphere scale height 5 km
	two, hBL := a.TwoScale(zBL, hFree)
	ts, err := a.SolveGlobal([]Column{{"g", 1, two}}, Storms{Is: a.GlobalCurrent()})
	if err != nil {
		return err
	}
	fmt.Printf("  two-scale column (declared: boundary layer %.1f km, free-air scale height %.0f km; solved H_BL = %.2f km):\n", zBL/1e3, hFree/1e3, hBL/1e3)
	fmt.Printf("    V_I %.2f kV, E0 %.1f V/m, half the space charge below %.2f km; τ = ε0/σ0 = %.0f s unchanged —\n",
		ts.VI/1e3, ts.E0[0], ts.ZHalf[0]/1e3, eps0/s0)
	fmt.Println("    Q = ε0·E0·A and I = J·A make R·C_eff = ε0/σ0 for any profile.")

	// Land and ocean columns (declared surface conductivities, common profile shape).
	const landFrac = 0.29
	sl, so := 1.0e-14, 2.2e-14 // declared: aerosol-laden continental air conducts less than maritime air
	mk := func(sg float64) func(float64) float64 {
		return func(z float64) float64 { return sg * math.Exp(z/H) }
	}
	cols := []Column{{"land", landFrac, mk(sl)}, {"ocean", 1 - landFrac, mk(so)}}
	lo, err := a.SolveGlobal(cols, Storms{Is: a.GlobalCurrent()})
	if err != nil {
		return err
	}
	fmt.Printf("  land (%.0f %%, σ0 %.2g S/m) and ocean (σ0 %.2g S/m) columns in parallel, same storm current:\n", 100*landFrac, sl, so)
	fmt.Printf("    V_I %.1f kV; land: %.0f A, E0 %.0f V/m; ocean: %.0f A, E0 %.0f V/m (the same V_I over both: J ∝ σ, E0 equal)\n",
		lo.VI/1e3, lo.IFair[0], lo.E0[0], lo.IFair[1], lo.E0[1])

	// Wilson generators with their own resistances.
	fmt.Println("  storm generators with an upward and a downward path (declared, all storms in parallel):")
	fmt.Println("    R_up [Ω]  R_down [Ω]  R_cloud [Ω]  I_s needed [A]  reaching ionosphere  P_gen [MW]  in storms [MW]  fair weather [MW]  ledger")
	for _, r := range [][3]float64{{50, 100, 1000}, {100, 100, 500}, {100, 500, 250}} {
		// Solve for I_s that keeps the anchors' global current: linear, two runs.
		st := Storms{Is: 1, RUp: r[0], RDown: r[1], RCloud: r[2], Distributed: true}
		g1, err := a.SolveGlobal(single, st)
		if err != nil {
			return err
		}
		st.Is = a.GlobalCurrent() / g1.IUp
		g, err := a.SolveGlobal(single, st)
		if err != nil {
			return err
		}
		fmt.Printf("    %-8.0f  %-10.0f  %-11.0f  %8.0f        %5.1f %%              %6.0f      %6.0f          %6.0f       %.0e\n",
			r[0], r[1], r[2], st.Is, 100*g.IUp/st.Is, g.PGen/1e6, g.PStorm/1e6, g.PFair/1e6, g.Ledger)
		if math.Abs(g.PGen-g.PStorm-g.PFair) > 1e-5*g.PGen { // engine round-off ≈ 5e-7 of Σ|P|
			return fmt.Errorf("storm books: P_gen %g ≠ %g + %g", g.PGen, g.PStorm, g.PFair)
		}
	}
	fmt.Println("  the generator must lift more current, and do more work, than the fair-weather circuit carries;")
	fmt.Println("  the rest returns locally (lightning to ground, point discharge, precipitation currents).")
	fmt.Println()
	return nil
}
