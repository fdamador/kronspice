// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// twinF1 is the Tesla bifilar solenoid's open-ends resonance; twinCycles the
// pulse periods per transient run.
const (
	twinF1     = 41.61e3
	twinCycles = 151
)

var twinFreqs = []float64{1e4, twinF1, 1e5}

// TestNgspiceTwinShort writes the Tesla bifilar solenoid with its ends open
// and joined through 1 mΩ as SPICE netlists (AC at three frequencies, and
// the pulse transient at f1) for ngspice, and prints kronspice's values of
// the same quantities: per volt of AC drive the coil-end voltage, the
// current through the join, the mid-coil wire current and the primary
// current; in the transient their peaks over the last 10 % of the run.
// Runs only with KRONSPICE_NGSPICE_TWIN=<dir>; compare with
// viz/ngspice_twin.py (docker image kronspice-tools runs ngspice).
func TestNgspiceTwinShort(t *testing.T) {
	dir := os.Getenv("KRONSPICE_NGSPICE_TWIN")
	if dir == "" {
		t.Skip("set KRONSPICE_NGSPICE_TWIN=<dir> to write the ngspice twin netlists")
	}
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	for _, short := range []float64{0, 1e-3} {
		c := NewCoil(Solenoid(), true)
		c.Short = short
		tag := map[bool]string{false: "open", true: "joined"}[short > 0]
		top := "n" + strconv.Itoa(c.Top)
		mid := NumSections()/2 - 1 // section 9: m9 → n10
		// AC.
		rAC := c.SectionResistance(1e4)
		ckt, err := c.Build(circuit.Element{Kind: circuit.VoltageSource, ACMag: 1}, rAC)
		if err != nil {
			t.Fatal(err)
		}
		var ctl []string
		ctl = append(ctl, ".control", "set wr_singlescale", "set wr_vecnames")
		for i, f := range twinFreqs {
			ctl = append(ctl, fmt.Sprintf("ac lin 1 %g %g", f, f),
				fmt.Sprintf("wrdata %s_ac%d.txt vr(%s) vi(%s) vr(n0) vi(n0) vr(m%d) vi(m%d) vr(n%d) vi(n%d) vr(s) vi(s) vr(p) vi(p)",
					tag, i, top, top, mid, mid, mid+1, mid+1))
		}
		ctl = append(ctl, "quit", ".endc")
		if err := writeNet(filepath.Join(dir, tag+"_ac.cir"), ckt, ctl); err != nil {
			t.Fatal(err)
		}
		for i, f := range twinFreqs {
			r, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
			if err != nil {
				t.Fatal(err)
			}
			v := func(n string) complex128 { x, _ := r.NodeVoltage(n); return x[0] }
			ends := cabs(v(top) - v("n0"))
			fmt.Fprintf(&out, "kronspice %s ac%d f=%g ends=%.6g join=%.6g wire=%.6g prim=%.6g\n", tag, i, f,
				ends, joinCurrent(ends, short), cabs(v(fmt.Sprintf("m%d", mid))-v(fmt.Sprintf("n%d", mid+1)))/rAC,
				cabs(v("s")-v("p"))/SourceOhms)
		}
		// Transient.
		rTr := c.SectionResistance(twinF1)
		per := 1 / twinF1
		src := circuit.Element{Kind: circuit.VoltageSource, Wave: circuit.Waveform{Kind: circuit.WavePulse,
			Args: []float64{0, PulseVolts, 0, PulseEdge / twinF1, PulseEdge / twinF1, 0.1 / twinF1, per}}}
		ckt, err = c.Build(src, rTr)
		if err != nil {
			t.Fatal(err)
		}
		stop := twinCycles / twinF1
		step := per / 50
		ctl = []string{
			fmt.Sprintf(".options abstol=1e-6 vntol=1e-3 chgtol=1e-12 reltol=1e-3 method=trap"),
			".control", "set wr_singlescale", "set wr_vecnames",
			fmt.Sprintf("tran %g %g 0 %g", step, stop, step),
			fmt.Sprintf("wrdata %s_tran.txt v(%s) v(n0) v(m%d) v(n%d) v(s) v(p)", tag, top, mid, mid+1),
			"quit", ".endc"}
		if err := writeNet(filepath.Join(dir, tag+"_tran.cir"), ckt, ctl); err != nil {
			t.Fatal(err)
		}
		w, err := analysis.Tran(ckt, analysis.TranSpec{Step: step, Stop: stop, MaxStep: step}, TranOpts)
		if err != nil {
			t.Fatal(err)
		}
		nv := func(n string) []float64 { x, _ := w.NodeVoltage(n); return x }
		vt, v0, vm, vn, vs, vp := nv(top), nv("n0"), nv(fmt.Sprintf("m%d", mid)), nv(fmt.Sprintf("n%d", mid+1)), nv("s"), nv("p")
		var ends, wire, prim float64
		for i, tt := range w.Time {
			if tt < 0.9*stop {
				continue
			}
			ends = math.Max(ends, math.Abs(vt[i]-v0[i]))
			wire = math.Max(wire, math.Abs(vm[i]-vn[i])/rTr)
			prim = math.Max(prim, math.Abs(vs[i]-vp[i])/SourceOhms)
		}
		fmt.Fprintf(&out, "kronspice %s tran ends=%.6g join=%.6g wire=%.6g prim=%.6g\n", tag, ends, joinCurrent(ends, short), wire, prim)
		fmt.Fprintf(&out, "rsection %s ac=%.9g tran=%.9g\n", tag, rAC, rTr)
	}
	if err := os.WriteFile(filepath.Join(dir, "kronspice.txt"), out.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	t.Log("\n" + out.String())
}

func joinCurrent(v, short float64) float64 {
	if short == 0 {
		return 0
	}
	return v / short
}

func writeNet(path string, c *circuit.Circuit, control []string) error {
	var b bytes.Buffer
	if err := c.WriteNetlist(&b, control...); err != nil {
		return err
	}
	return os.WriteFile(path, b.Bytes(), 0o644)
}
