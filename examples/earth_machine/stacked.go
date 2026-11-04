// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// The stacked design: the two disks of Rikitake's machine on one axis, one
// above the other, each inside the coil fed by the other disk, instead of
// side by side. In the core's terms: two eddies on one axis, as convection
// columns aligned with the rotation axis stack, rather than two eddies side
// by side.
//
// Stacking changes the machine in two ways, both fixed by the geometry (the
// copper machine of examples/homopolar/dual.go, scaled to the disk radius:
// r_i/r_o = 1/15, coil radius 1.2 r_o, coil length 0.4 r_o):
//
//   - each coil's stray field threads its OWN disk too, so every loop feeds
//     itself as well as the other: a built-in feedback coil of strength
//     KS = M_far/M, with no copper of its own;
//   - the two coils share flux: the loops have a mutual inductance KM·L.
//
// Both fall off with the gap between the disks. The question asked here is
// the one the side-by-side study answered for a feedback coil (task 4c, where
// any aiding K ≥ 0.125 locked the field to one polarity): how far apart must
// stacked eddies sit for the field to reverse again?

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
)

// Stacked geometry, in units of the disk radius (examples/homopolar).
const (
	stRi, stRc, stLc = 1.0 / 15, 1.2, 0.4
	stN              = 4 // turns: the ratios KS and KM do not depend on it
	// StackedGapHomopolar is the copper machine's gap, 8 cm over a 15 cm disk.
	StackedGapHomopolar = 0.08 / 0.15
)

// annulus is a coil's coupling to a disk whose plane is z from the coil's
// mid-plane: (N/2π)·⟨M_circ(r_c, r_o, u) − M_circ(r_c, r_i, u)⟩ over the turns.
func annulus(z float64) float64 {
	const n = 64
	var s float64
	for k := 0; k < n; k++ {
		u := z + stLc*((float64(k)+0.5)/n-0.5)
		s += em.MutualCircles(stRc, 1, u) - em.MutualCircles(stRc, stRi, u)
	}
	return stN / (2 * math.Pi) * s / n
}

// StackedCouplings returns KS = M_far/M (a coil on the disk one gap away)
// and KM = M_coils/L (the two coils a gap apart) for a gap in disk radii.
func StackedCouplings(gap float64) (ks, km float64) {
	ks = annulus(gap) / annulus(0)
	const n = 32
	var s float64
	for a := 0; a < n; a++ {
		for b := 0; b < n; b++ {
			s += em.MutualCircles(stRc, stRc, gap+stLc*float64(a-b)/n)
		}
	}
	km = stN * stN * s / (n * n) / em.Nagaoka(stRc, stLc, stN)
	return
}

// StackedRikitake is Rikitake's machine stacked with the given gap.
func StackedRikitake(gap float64) Params {
	p := Rikitake()
	p.Stacked = true
	p.KS, p.KM = StackedCouplings(gap)
	return p
}

// stackedCircuit: the coils are inductors LC1, LC2 with a Mutual element
// (KM), and the machine elements carry 1 nH sense windings, all along their
// loops, for the speed EMFs: G[0][1] = −M (the other loop's coil, the main
// field on this disk) and G[0][2] = −(K + KS)·M (this loop's own coil, its
// stray field). The element books ω·G·i as a drop along a winding.
//
//	loop i: a_i →[N_i disk]→ b_i →R_i→ c_i →[N_i own-coil sense]→ d_i →[N_j coil-i sense]→ e_i →LC_i→ a_i
func (p Params) stackedCircuit() (*circuit.Circuit, error) {
	spec := func(iOwn, iOther float64) *circuit.MachineSpec {
		s := circuit.NewMachineSpec(3)
		s.L0[0][0], s.L0[1][1], s.L0[2][2] = 1e-9, 1e-9, 1e-9
		s.G[0][1] = -p.M
		s.G[0][2] = -(p.K + p.KS) * p.M
		s.IC[0], s.IC[1], s.IC[2] = iOwn, iOther, iOwn
		return s
	}
	b := circuit.NewBuilder("Rikitake machine, disks stacked").Domain(exchange.Mechanical, "w1", "w2")
	b.MachineModel("disk1", spec(p.I1, p.I2)).MachineModel("disk2", spec(p.I2, p.I1))
	b.Machine("N1", "disk1", [][2]string{{"a1", "b1"}, {"d2", "e2"}, {"c1", "d1"}}, "w1", "0")
	b.Machine("N2", "disk2", [][2]string{{"a2", "b2"}, {"d1", "e1"}, {"c2", "d2"}}, "w2", "0")
	b.Resistor("R1", "b1", "c1", p.R).Resistor("R2", "b2", "c2", p.R)
	b.InductorIC("LC1", "e1", "a1", p.L, p.I1).InductorIC("LC2", "e2", "a2", p.L, p.I2)
	if p.KM != 0 {
		b.Mutual("KC", "LC1", "LC2", p.KM)
	}
	for _, n := range []string{"a1", "b1", "c1", "d1", "e1", "a2", "b2", "c2", "d2", "e2"} {
		b.Resistor("RREF"+n, n, "0", 1e7)
	}
	p.shafts(b)
	return b.Build()
}

// StackedRow is one gap of the study.
type StackedRow struct {
	Gap, KS, KM      float64
	Reversals        int
	Mean, CV, RMS    float64
	MeanW1, MeanW2   float64
	ThresholdSpeedUp float64 // side-by-side threshold speed / stacked one = 1 + KS
}

// StackedStudy runs RK4 statistics (t = 20 000, μ = 1, A = 1) at each gap.
func StackedStudy(gaps []float64) []StackedRow {
	var out []StackedRow
	for _, g := range gaps {
		p := StackedRikitake(g)
		r := p.StatRun(20000, 0.002, 200, 1e-2)
		m, cv, _ := meanCV(r.Intervals)
		tr := p.RunRK4(2000, 0.01)
		var w1, w2 float64
		h := len(tr.W1) / 2
		for k := h; k < len(tr.W1); k++ {
			w1 += tr.W1[k]
			w2 += tr.W2[k]
		}
		n := float64(len(tr.W1) - h)
		out = append(out, StackedRow{Gap: g, KS: p.KS, KM: p.KM, Reversals: len(r.Rev), Mean: m, CV: cv, RMS: r.RMS,
			MeanW1: w1 / n, MeanW2: w2 / n, ThresholdSpeedUp: 1 + p.KS})
	}
	return out
}

func runStacked() error {
	fmt.Println("Stacked design — the two disks (eddies) on one axis, each inside the other's coil")
	ks, km := StackedCouplings(StackedGapHomopolar)
	fmt.Printf("  geometry of examples/homopolar scaled to the disk radius (coil 1.2 r_o, 0.4 r_o long); at its gap %.2f r_o:\n", StackedGapHomopolar)
	fmt.Printf("  each coil's stray field on its own disk KS = %.3f of M (homopolar: 0.519), coil–coil coupling KM = %.3f (0.387)\n", ks, km)
	fmt.Println("  so a stacked loop excites itself at (1 + KS)·M·ω: a feedback coil built into the geometry.")

	// Twin: circuit (inductors + Mutual, sense windings) against the RK4 equations.
	p := StackedRikitake(StackedGapHomopolar)
	tr, mis, err := p.RunCircuit(10, 0.01)
	if err != nil {
		return err
	}
	ref := p.RunRK4(10, 0.01)
	var dmax float64
	for k := 1; k < len(tr.T); k++ {
		dmax = math.Max(dmax, math.Abs(tr.I1[k]-ref.I1[int(math.Round(tr.T[k]/0.01))]))
	}
	fmt.Printf("  circuit vs RK4 at that gap, t ≤ 10: |ΔI1| ≤ %.2g; engine ledger mismatch %.1g J\n", dmax, mis)

	fmt.Println("  μ = 1, A = 1, t = 20 000 (RK4), against the gap:")
	fmt.Println("    gap [r_o]   KS      KM      reversals  mean interval  CV     RMS I1   mean ω1, ω2")
	base := Rikitake().StatRun(20000, 0.002, 200, 1e-2)
	bm, bcv, _ := meanCV(base.Intervals)
	fmt.Printf("    side by side  0       0       %5d      %8.2f       %.2f   %.3f\n", len(base.Rev), bm, bcv, base.RMS)
	for _, r := range StackedStudy([]float64{StackedGapHomopolar, 1, 1.5, 2, 3, 4, 6}) {
		fmt.Printf("    %-9.2f   %.4f  %.4f  %5d      %8.2f       %.2f   %.3f    %.2f, %.2f\n",
			r.Gap, r.KS, r.KM, r.Reversals, r.Mean, r.CV, r.RMS, r.MeanW1, r.MeanW2)
	}
	fmt.Println("  no gap reverses: Rikitake's chaos rests on a conserved speed difference A = ω1 − ω2 (both disks feel")
	fmt.Println("  the same back torque M·I1·I2). Any self-coupling adds KS·M·Ii², different on each disk, so A drains away:")
	p3 := StackedRikitake(6)
	ta := p3.RunRK4(2000, 0.5)
	fmt.Printf("    gap 6 r_o (KS = %.4f): A = %.2f at t = 0, %.2f at 50, %.2f at 200, %.3f at 1000\n",
		p3.KS, ta.W1[0]-ta.W2[0], ta.W1[100]-ta.W2[100], ta.W1[400]-ta.W2[400], ta.W1[2000]-ta.W2[2000])
	fmt.Println("  and with A gone the field locks to one polarity, a steady dynamo, however weak the stray coupling.")
	fmt.Println("  Holding A with unequal drive torques, G1 = 1 + δ, G2 = 1 − δ:")
	fmt.Println("    gap [r_o]   δ      reversals  mean interval  CV     A held at")
	for _, g := range []float64{StackedGapHomopolar, 1.5} {
		for _, d := range []float64{0.1, 0.3, 0.5} {
			q := StackedRikitake(g)
			q.G1, q.G2 = 1+d, 1-d
			r := q.StatRun(20000, 0.002, 200, 1e-2)
			m, cv, _ := meanCV(r.Intervals)
			t5 := q.RunRK4(5000, 0.5)
			k := len(t5.T) - 1
			fmt.Printf("    %-9.2f   %.1f    %5d      %8.2f       %.2f   %.3f\n", g, d, len(r.Rev), m, cv, t5.W1[k]-t5.W2[k])
		}
	}
	fmt.Println("  only close stacking with a strong imbalance reverses, and then like a clock (interval ≈ the 3.9")
	fmt.Println("  oscillation period, CV 0.13), not chaotically. Stacked eddies make a stable dynamo; irregular")
	fmt.Println("  reversals need the two loops to couple only through each other's field, as side by side.")
	fmt.Println("  The ledger is the same as side by side: the drive pays for field and heat, stacking adds no input.")
	fmt.Println()
	return nil
}

// stackedRunEngine is used by the tests: the circuit at a gap, with its
// engine ledger.
func stackedRunEngine(gap, tEnd, dt float64) (Trace, float64, error) {
	p := StackedRikitake(gap)
	return p.RunCircuit(tEnd, dt)
}

var _ = analysis.Tran // the engine runs through RunCircuit
