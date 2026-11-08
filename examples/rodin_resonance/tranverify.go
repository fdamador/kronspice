// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// The kronspice transient of a net (lumped or reduced sectional, at most
// 64 inductive branches including the groups): every inductive branch is
// a winding of one N element (kronspice's coupled-inductor element with a
// full inductance matrix, mutuals of either sign; its mechanical port is
// held still, so it is a fixed coupled inductor), with its wire resistance
// in series; node capacitances as capacitors; the sine (A) and cosine (B)
// sources through SourceOhms; the groups through their tuning capacitors
// into Schottky bridges and one capacitor bank.

// DiodeModel is a small Schottky rectifier (declared, 1N5819-class).
var DiodeModel = map[string]float64{"is": 3e-6, "n": 1.05, "rs": 0.03, "cjo": 100e-12}

// LumpedNet is the two-winding lumped network (no capacitance), valid well
// below the windings' self-resonance.
func (m *Model) LumpedNet(sense int) *Net {
	n := &Net{Top: DrivenPair, NS: 1, Nodes: 2, UseB: true}
	n.A = []int{0, -1}
	if sense == Toroidal {
		n.B = []int{1, -1}
	} else {
		n.B = []int{-1, 1}
	}
	n.L = [][]float64{{m.LA, m.M}, {m.M, m.LB}}
	n.Len = []float64{m.LenA, m.LenB}
	n.R = []float64{WireResistance(m.LenA, 0), WireResistance(m.LenB, 0)}
	n.Branches = [][2]int{{n.A[0], n.A[1]}, {n.B[0], n.B[1]}}
	n.CNode = newMat(2)
	return n
}

// lumpedSectionMutual: in a lumped net the two branches are the whole
// windings, so a group's mutual is the pass mutual summed (× Nc).
func (n *Net) groupMSec(m *Model, pass []float64) []float64 {
	if n.NS == 1 {
		a, b := 0.0, 0.0
		for q := 0; q < m.P.Sec; q++ {
			a += pass[q]
			b += pass[m.P.Sec+q]
		}
		return []float64{a * float64(m.P.Nc), b * float64(m.P.Nc)}
	}
	return n.SectionMutual(m, pass)
}

// TranSpec of a verification run.
type TranRun struct {
	Time, VBank, EBank []float64
	ESrc               []float64 // energy delivered by each source over the run [J]
	EBankEnd           float64
	ECopper, ESrcR     float64
	EDiode, EBleed     float64
	ECoilsEnd          float64 // energy left in inductors and other capacitors at the end [J]
	EIn, Mismatch      float64
	MaxRel             float64
	IA, IB             []float64
	ISrc               [2]complex128 // source current phasors over the last cycle
}

// TranVerify builds and runs the net with the layout's groups (nil groups:
// drive only, with a 1 Ω-free resistive check).
func (m *Model) TranVerify(n *Net, l Layout, withGroups bool, cycles, steps int, bank float64) (*TranRun, error) {
	f := l.F
	n.SetFreq(f)
	var gs []Group
	var mg [][]float64
	if withGroups {
		th := l.Thetas()
		lg := l.Shape.SelfL()
		gs = make([]Group, l.NG)
		for k := range gs {
			gs[k] = Group{L: lg, R: l.Shape.Resistance(f), Ctune: 1 / math.Pow(2*math.Pi*f, 2) / lg,
				MSec: n.groupMSec(m, m.PassMutual(l.Shape, th[k]))}
		}
		mg = newMat(l.NG)
		for i := range mg {
			mg[i][i] = lg
			for j := i + 1; j < l.NG; j++ {
				v := GroupMutual(l.Shape, th[i], th[j])
				mg[i][j], mg[j][i] = v, v
			}
		}
	}
	nb := len(n.L)
	nw := nb + len(gs)
	if nw > circuit.MaxWindings {
		return nil, fmt.Errorf("%d windings exceed the N element's %d", nw, circuit.MaxWindings)
	}
	spec := circuit.NewMachineSpec(nw)
	for i := 0; i < nb; i++ {
		for j := 0; j < nb; j++ {
			spec.L0[i][j] = n.L[i][j]
		}
		for k := range gs {
			spec.L0[i][nb+k] = gs[k].MSec[i]
			spec.L0[nb+k][i] = gs[k].MSec[i]
		}
	}
	for i := range gs {
		for j := range gs {
			spec.L0[nb+i][nb+j] = mg[i][j]
		}
	}
	b := circuit.NewBuilder("rodin pair with output groups").Domain(exchange.Mechanical, "w")
	b.MachineModel("coils", spec)
	b.Model("dsch", "D", DiodeModel)
	var wind [][2]string
	for g := 0; g < nb; g++ {
		wind = append(wind, [2]string{NodeName(n.Branches[g][0]), fmt.Sprintf("m%d", g)})
		b.Resistor(fmt.Sprintf("RW%d", g), fmt.Sprintf("m%d", g), NodeName(n.Branches[g][1]), n.R[g])
	}
	for k := range gs {
		wind = append(wind, [2]string{fmt.Sprintf("g%da", k), fmt.Sprintf("g%dr", k)})
		b.Resistor(fmt.Sprintf("RG%d", k), fmt.Sprintf("g%dr", k), fmt.Sprintf("g%dt", k), gs[k].R)
		b.Capacitor(fmt.Sprintf("CT%d", k), fmt.Sprintf("g%dt", k), fmt.Sprintf("g%db", k), gs[k].Ctune)
	}
	b.Machine("NC", "coils", wind, "w", "0")
	b.Resistor("RMECH", "w", "0", 1)
	// node capacitances
	c := 0
	for i := 0; i < n.Nodes; i++ {
		row := 0.0
		for j := 0; j < n.Nodes; j++ {
			row += n.CNode[i][j]
			if j > i && n.CNode[i][j] < 0 {
				b.Capacitor(fmt.Sprintf("CN%d", c), NodeName(i), NodeName(j), -n.CNode[i][j])
				c++
			}
		}
		if row > 1e-18 {
			b.Capacitor(fmt.Sprintf("CN%d", c), NodeName(i), "0", row)
			c++
		}
	}
	// sources: A = sin, B = cos
	src := n.DriveSources(l.Sense, Part3Volts)
	for i, s := range src {
		nm := []string{"A", "B"}[i]
		lo := "s" + strings.ToLower(nm)
		b.VoltageWave("V"+nm, lo, "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, s.V, f, 0, 0, s.PhaseDeg + 90}})
		// soft start: a series EXP source cancels the sinusoid's value at
		// t = 0 (the cosine starts at its peak) and fades out over two
		// cycles, so the run starts from rest; both belong to the source
		v0 := s.V * math.Sin((s.PhaseDeg+90)*math.Pi/180)
		top := lo
		if math.Abs(v0) > 1e-12 {
			top = lo + "x"
			b.VoltageWave("V"+nm+"X", top, lo, circuit.Waveform{Kind: circuit.WaveExp, Args: []float64{-v0, 0, 0, 2 / f, 1e3, 1}})
		}
		b.Resistor("RS"+nm, top, NodeName(s.Node), s.R)
	}
	// groups → bridges → bank
	if len(gs) > 0 {
		var inputs [][2]string
		switch l.Comb {
		case Separate:
			for k := range gs {
				inputs = append(inputs, [2]string{fmt.Sprintf("g%da", k), fmt.Sprintf("g%db", k)})
			}
		case Series:
			for k := 0; k+1 < len(gs); k++ {
				b.Resistor(fmt.Sprintf("RJ%d", k), fmt.Sprintf("g%db", k), fmt.Sprintf("g%da", k+1), 1e-6)
			}
			inputs = append(inputs, [2]string{"g0a", fmt.Sprintf("g%db", len(gs)-1)})
		case Parallel:
			for k := range gs {
				b.Resistor(fmt.Sprintf("RJA%d", k), fmt.Sprintf("g%da", k), "pa", 1e-6)
				b.Resistor(fmt.Sprintf("RJB%d", k), fmt.Sprintf("g%db", k), "pb", 1e-6)
			}
			inputs = append(inputs, [2]string{"pa", "pb"})
		}
		for k, in := range inputs {
			b.Diode(fmt.Sprintf("D%d_1", k), in[0], "bk", "dsch")
			b.Diode(fmt.Sprintf("D%d_2", k), in[1], "bk", "dsch")
			b.Diode(fmt.Sprintf("D%d_3", k), "0", in[0], "dsch")
			b.Diode(fmt.Sprintf("D%d_4", k), "0", in[1], "dsch")
		}
		b.Capacitor("CBANK", "bk", "0", bank)
		// DC reference for the floating group loops: a choke from each
		// terminal to ground, reactance 10⁴ × the larger of R_L and the
		// coil's resistance at f (stores, does not dissipate)
		lch := 1e4 * math.Max(l.RLoad, gs[0].R) / (2 * math.Pi * f)
		for k := range gs {
			b.Inductor(fmt.Sprintf("LCHA%d", k), fmt.Sprintf("g%da", k), "0", lch)
			b.Inductor(fmt.Sprintf("LCHB%d", k), fmt.Sprintf("g%db", k), "0", lch)
		}
	}
	ckt, err := b.BuildTransient()
	if err != nil {
		return nil, err
	}
	T := 1 / f
	stop := float64(cycles) * T
	step := T / float64(steps)
	w, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step},
		analysis.TranOptions{Options: analysis.Options{AbsTol: 1e-9, VnTol: 1e-6}})
	if err != nil {
		return nil, err
	}
	r := &TranRun{Time: w.Time}
	if vb, ok := w.NodeVoltage("bk"); ok {
		r.VBank = vb
		for _, v := range vb {
			r.EBank = append(r.EBank, 0.5*bank*v*v)
		}
		r.EBankEnd = r.EBank[len(r.EBank)-1]
	}
	energy := func(e int) float64 {
		s := 0.0
		for k := 1; k < len(w.Time); k++ {
			s += 0.5 * (w.P[e][k] + w.P[e][k-1]) * (w.Time[k] - w.Time[k-1])
		}
		return s
	}
	r.ESrc = make([]float64, 2)
	for e, name := range w.ElementNames {
		u := strings.ToUpper(name)
		switch {
		case u == "VA" || u == "VAX":
			r.ESrc[0] -= energy(e)
		case u == "VB" || u == "VBX":
			r.ESrc[1] -= energy(e)
		case strings.HasPrefix(u, "RW") || strings.HasPrefix(u, "RG"):
			r.ECopper += energy(e)
		case strings.HasPrefix(u, "RS"):
			r.ESrcR += energy(e)
		case strings.HasPrefix(u, "D"):
			r.EDiode += energy(e)
		case strings.HasPrefix(u, "RJ") || u == "RMECH":
			r.EBleed += energy(e)
		}
	}
	l2 := w.Ledger
	r.EIn, r.Mismatch, r.MaxRel = l2.EnergyIn, l2.IntegrationMismatch, l2.MaxRelativeResidual
	r.ECoilsEnd = l2.StoredEnd - r.EBankEnd
	ia, _ := w.ElementCurrent("RSA")
	ib, _ := w.ElementCurrent("RSB")
	r.IA, r.IB = ia, ib
	r.ISrc[0] = phasor(w.Time, ia, f, stop-T)
	r.ISrc[1] = phasor(w.Time, ib, f, stop-T)
	return r, nil
}
