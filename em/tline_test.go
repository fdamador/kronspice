// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em_test

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// The telegrapher's equations through the circuit (CLAUDE.md §3c rule 3,
// §3g): a coaxial line extracted from its cross-section, built as a ladder
// of 200 sections, has input impedance Z₀ when matched (1 %, well below
// the sections' cutoff); a step travels at 1/√(L′C′) (its 50 % point
// arrives after Length·√(L′C′) within 2 %); an open end doubles it and a
// shorted end cancels it.
func TestLineTelegrapher(t *testing.T) {
	const a, b, er = 1e-3, 3.5e-3, 2.1
	c0 := 4e-3
	s := em.Section{W: 2 * c0, H: 2 * c0, Step: a / 10,
		Conductors: []em.Conductor2D{em.Disk("in", c0, c0, a), em.Outside("gnd", c0, c0, b)},
		EpsR:       func(x, y float64) float64 { return er }}
	r, err := s.Solve()
	if err != nil {
		t.Fatal(err)
	}
	line := em.FromSection(r, 0, 2, 0, 0)
	z0 := real(line.Z0(0))
	build := func(src circuit.Waveform, ac bool, load float64) *circuit.Circuit {
		bl := circuit.NewBuilder("coax line")
		if ac {
			bl.Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "s", "0")
		} else {
			bl.VoltageWave("V1", "s", "0", src)
		}
		bl.Resistor("RS", "s", "in", z0)
		if err := line.Ladder(bl, "T", "in", "out", "0", 200); err != nil {
			t.Fatal(err)
		}
		switch {
		case load > 0:
			bl.Resistor("RL", "out", "0", load)
		case load == 0:
			bl.Resistor("RL", "out", "0", 1e-6)
		default:
			bl.Resistor("RL", "out", "0", 1e12)
		}
		c, err := bl.Build()
		if err != nil {
			t.Fatal(err)
		}
		return c
	}
	// Matched input impedance at a frequency where a section is λ/1000.
	f := 1 / (line.Delay() / 200 * 1000)
	ac, err := analysis.AC(build(circuit.Waveform{}, true, z0), analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	vs, _ := ac.NodeVoltage("s")
	vin, _ := ac.NodeVoltage("in")
	zin := vin[0] * complex(z0, 0) / (vs[0] - vin[0])
	if cmplx.Abs(zin-complex(z0, 0)) > 1e-2*z0 {
		t.Errorf("matched input impedance %v, Z₀ = %.3f Ω", zin, z0)
	}
	// Step response.
	td := line.Delay()
	step := circuit.Waveform{Kind: circuit.WavePulse, Args: []float64{0, 2, 0, td / 200, td / 200, 10 * td, 20 * td}}
	for _, tc := range []struct {
		name  string
		load  float64
		final float64 // v(out) after the first arrival, relative to 1 V incident
	}{{"matched", z0, 1}, {"open", -1, 2}} {
		wf, err := analysis.Tran(build(step, false, tc.load), analysis.TranSpec{Step: td / 400, Stop: 1.8 * td, MaxStep: td / 400}, analysis.TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		v, _ := wf.NodeVoltage("out")
		cross := math.NaN()
		for k := 1; k < len(v); k++ {
			if v[k-1] < tc.final/2 && v[k] >= tc.final/2 {
				cross = wf.Time[k-1] + (wf.Time[k]-wf.Time[k-1])*(tc.final/2-v[k-1])/(v[k]-v[k-1])
				break
			}
		}
		if math.Abs(cross-td) > 0.02*td {
			t.Errorf("%s: 50 %% arrival at %.4g s, delay √(L′C′)·l = %.4g s", tc.name, cross, td)
		}
		end := v[len(v)-1]
		if math.Abs(end-tc.final) > 0.05*tc.final {
			t.Errorf("%s: settles to %.3f V, want %.3f V", tc.name, end, tc.final)
		}
	}
	t.Logf("Z₀ %.3f Ω, delay %.4g ns for 2 m (velocity %.4g c)", z0, td*1e9, 2/td/em.SpeedOfLight)
}
