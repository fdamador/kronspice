// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// VIPoint is one frequency of a coil's AC sweep, per volt of drive on the
// primary's source (amplitudes).
type VIPoint struct {
	Freq    float64
	VEnds   float64 // |V(top) − V(base)| [V]
	VTurn   float64 // largest |V| between physically neighbouring turns (interpolated in each section) [V]
	IWire   float64 // largest wire current in any section [A]
	IPrim   float64 // primary current [A]
	PSource float64 // mean power from the source [W]
}

// VISweep sweeps the coil's AC response from fLo to fHi (points per
// decade) with the drive through the primary, as Resonance does, and
// returns voltage, current and power at every frequency per volt of drive.
// The wire resistance is taken at 1 MHz throughout (as Resonance's coarse
// sweep), so Q is not exact away from that frequency.
func (c *Coil) VISweep(fLo, fHi float64, perDecade int) ([]VIPoint, error) {
	ckt, err := c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, c.SectionResistance(1e6))
	if err != nil {
		return nil, err
	}
	res, err := analysis.AC(ckt, analysis.ACSpec{Variation: "dec", Points: perDecade, FStart: fLo, FStop: fHi}, analysis.Options{})
	if err != nil {
		return nil, err
	}
	ns := NumSections()
	nodes := make([][]complex128, c.NNodes)
	for k := range nodes {
		nodes[k], _ = res.NodeVoltage("n" + strconv.Itoa(k))
	}
	wire := make([][]complex128, ns)
	for s := range wire {
		wire[s], _ = res.ElementCurrent(fmt.Sprintf("R%d", s))
	}
	prim, _ := res.ElementCurrent("RS")
	src, _ := res.NodeVoltage("s")
	at := make([]int, c.Geom.N()) // wire turn at each physical position
	for t, p := range c.Pos {
		at[p] = t
	}
	out := make([]VIPoint, len(res.Freq))
	for i, f := range res.Freq {
		p := VIPoint{Freq: f, VEnds: cabs(nodes[c.Top][i] - nodes[0][i]), IPrim: cabs(prim[i])}
		for s := 0; s < ns; s++ {
			p.IWire = math.Max(p.IWire, cabs(wire[s][i]))
		}
		turn := func(t int) complex128 {
			from, to, fr := c.turnFrac(t)
			return nodes[from][i]*complex(1-fr, 0) + nodes[to][i]*complex(fr, 0)
		}
		for q := 0; q+1 < len(at); q++ {
			p.VTurn = math.Max(p.VTurn, cabs(turn(at[q])-turn(at[q+1])))
		}
		// Mean power from a 1 V amplitude source: ½·Re(V·I*), I into the primary branch.
		v, ii := src[i], prim[i]
		p.PSource = 0.5 * (real(v)*real(ii) + imag(v)*imag(ii))
		out[i] = p
	}
	return out, nil
}

// WriteVISweep sweeps the five solenoid windings (ordinary, Tesla
// bifilar, opposite-current bifilar in series, and the two parallel
// bifilars) from 1 kHz to the sections' λ/10 limit and writes
// viz/vi_sweep.csv with one column group per coil.
func WriteVISweep(dir string) error {
	coils := []*Coil{NewCoil(Solenoid(), false), NewCoil(Solenoid(), true), NewCounterCoil(Solenoid()),
		NewParallelCoil(Solenoid(), false), NewParallelCoil(Solenoid(), true)}
	fHi := coils[0].SectionLimit()
	var sweeps [][]VIPoint
	for _, c := range coils {
		s, err := c.VISweep(1e3, fHi, 400)
		if err != nil {
			return fmt.Errorf("%s: %w", c.Name, err)
		}
		sweeps = append(sweeps, s)
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, "vi_sweep.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	head := []string{"freq"}
	for _, c := range coils {
		k := c.Kind()
		head = append(head, k+"_v_ends", k+"_v_turn", k+"_i_wire", k+"_i_prim", k+"_p_source")
	}
	w.Write(head)
	for i := range sweeps[0] {
		rec := []string{strconv.FormatFloat(sweeps[0][i].Freq, 'g', 7, 64)}
		for _, s := range sweeps {
			p := s[i]
			for _, x := range []float64{p.VEnds, p.VTurn, p.IWire, p.IPrim, p.PSource} {
				rec = append(rec, strconv.FormatFloat(x, 'g', 7, 64))
			}
		}
		w.Write(rec)
	}
	w.Flush()
	if err := w.Error(); err != nil {
		return err
	}
	fmt.Println("V–I sweep per volt of drive (1 kHz to the λ/10 limit), at each coil's peak:")
	for k, c := range coils {
		best := sweeps[k][0]
		for _, p := range sweeps[k] {
			if p.VEnds > best.VEnds {
				best = p
			}
		}
		fmt.Printf("   %-56s f %.4g Hz: coil ends %.4g V, wire %.4g A, V/I %.4g Ω, ω·L = sqrt(L/C) %.4g Ω\n",
			c.Name, best.Freq, best.VEnds, best.IWire, best.VEnds/best.IWire, math.Sqrt(c.TotalL()*sq(2*math.Pi*best.Freq)*c.TotalL()))
	}
	return nil
}

func sq(x float64) float64 { return x * x }
