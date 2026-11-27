// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"bytes"
	"math"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
	"github.com/fdamador/kronspice/output"
)

func deckCircuit(t *testing.T, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func runTran(t *testing.T, c *circuit.Circuit, step, stop, maxStep float64) *output.Waveforms {
	t.Helper()
	w, err := Tran(c, TranSpec{Step: step, Stop: stop, MaxStep: maxStep}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if r := w.Ledger.MaxRelativeResidual; r > 1e-9 {
		t.Errorf("ledger residual %.3g", r)
	}
	return w
}

// sampleAt interpolates y(t) linearly.
func sampleAt(ts, y []float64, t float64) float64 {
	for k := 1; k < len(ts); k++ {
		if ts[k] >= t {
			f := (t - ts[k-1]) / (ts[k] - ts[k-1])
			return y[k-1] + f*(y[k]-y[k-1])
		}
	}
	return y[len(y)-1]
}

// Under a constant current every law relaxes to its steady state with the
// exact solution of its equation: Mayr g = q + (g0 − q)e^(−t/τ), q = I²/P0,
// so v → P0/I; Cassie g² = q + (g0² − q)e^(−2t/τ), q = I²/U0², so v → U0;
// Habedank v → U0 + P0/I.
func TestArcConstantCurrent(t *testing.T) {
	const amps, tau, p0, u0, g0 = 10.0, 10e-6, 1e3, 50.0, 0.05
	cases := []struct {
		name, line string
		exact      func(t float64) float64 // voltage
	}{
		{"mayr", "A1 a 0 MAYR TAU=10u P0=1k G0=0.05", func(t float64) float64 {
			q := amps * amps / p0
			return amps / (q + (g0-q)*math.Exp(-t/tau))
		}},
		{"cassie", "A1 a 0 CASSIE TAU=10u U0=50 G0=0.05", func(t float64) float64 {
			q := amps * amps / (u0 * u0)
			return amps / math.Sqrt(q+(g0*g0-q)*math.Exp(-2*t/tau))
		}},
		{"habedank", "A1 a 0 HABEDANK TAUC=10u U0=50 TAUM=10u P0=1k G0=0.05", nil},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			c := deckCircuit(t, "arc\nI1 0 a DC 10\n"+tc.line+"\n.end\n")
			w := runTran(t, c, 1e-6, 200e-6, 0.2e-6)
			v, _ := w.NodeVoltage("a")
			if math.Abs(v[0]-amps/g0) > 1e-9*amps/g0 {
				t.Errorf("operating point %g V, want %g (frozen G0)", v[0], amps/g0)
			}
			if tc.exact != nil {
				for _, tt := range []float64{5e-6, 20e-6, 60e-6} {
					got, want := sampleAt(w.Time, v, tt), tc.exact(tt)
					if math.Abs(got-want) > 2e-3*want {
						t.Errorf("v(%g) = %g, want %g", tt, got, want)
					}
				}
			}
			want := map[string]float64{"mayr": p0 / amps, "cassie": u0, "habedank": u0 + p0/amps}[tc.name]
			if got := v[len(v)-1]; math.Abs(got-want) > 1e-4*want {
				t.Errorf("steady voltage %g, want %g", got, want)
			}
		})
	}
}

// A Mayr arc on a sinusoidal current (τ ≪ period): near the current peak
// it sits on its static characteristic v = P0/|i|, and its voltage peaks
// just before current zero (the extinction peak), the published current-zero
// behavior of the Mayr arc (Mayr 1943; Habedank 1993, Fig. 2). The arc is
// passive at every point.
func TestArcMayrCurrentZero(t *testing.T) {
	const ipk, f, tau, p0 = 100.0, 50.0, 20e-6, 1e4
	c := deckCircuit(t, "arc\nI1 0 a SIN(0 100 50)\nA1 a 0 MAYR TAU=20u P0=10k G0=1\n.end\n")
	w := runTran(t, c, 10e-6, 20e-3, 20e-6)
	v, _ := w.NodeVoltage("a")
	ia := w.I[1]
	tPeak := 0.25 / f
	if got, want := sampleAt(w.Time, v, tPeak), p0/ipk; math.Abs(got-want) > 0.01*want {
		t.Errorf("v at current peak %g, want static %g", got, want)
	}
	// Largest |v| in the first half cycle, and the current zero at T/2.
	vmax, tmax := 0.0, 0.0
	for k, tt := range w.Time {
		if tt > 0.5/f {
			break
		}
		if math.Abs(v[k]) > vmax {
			vmax, tmax = math.Abs(v[k]), tt
		}
	}
	if d := 0.5/f - tmax; d < 0 || d > 20*tau {
		t.Errorf("extinction peak at %g s, %g s before current zero (want within 20τ)", tmax, d)
	}
	if vmax < 5*p0/ipk {
		t.Errorf("extinction peak %g V not above the burning voltage %g V", vmax, p0/ipk)
	}
	for k := range w.Time {
		if v[k]*ia[k] < -1e-9 {
			t.Fatalf("arc delivers power at t = %g", w.Time[k])
		}
	}
}

// A spark-gap relaxation oscillator: a supply charges C through R until the
// gap breaks down at VBD, the Mayr arc discharges C through the lead
// inductance and extinguishes, the gap recovers and C charges again. (Without
// the inductance a Mayr arc on a charged capacitor runs away, g growing at
// v·i/(τ·P0), and the step falls to the femtosecond limit: that returns
// ErrTimestepTooSmall, as an evolved circuit of that kind must.) Every joule is booked: the energy from
// the supply equals the heat in R and the gap plus the energy left in C.
func TestArcSparkGapRelaxation(t *testing.T) {
	src := `spark gap
V1 in 0 PWL(0 0 1u 10k)
R1 in c 1meg
C1 c 0 1n
L1 c g 1u
A1 g 0 MAYR TAU=1u P0=10k VBD=5k GON=10m
.end
`
	c := deckCircuit(t, src)
	w, err := Tran(c, TranSpec{Step: 1e-6, Stop: 3e-3, MaxStep: 2e-6}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if r := w.Ledger.MaxRelativeResidual; r > 1e-9 {
		t.Errorf("ledger residual %.3g", r)
	}
	v, _ := w.NodeVoltage("c")
	// Count breakdowns as falls through VBD/2 and check the voltage never
	// stays far above VBD.
	falls, above := 0, false
	for k := range v {
		if v[k] > 5.2e3 {
			t.Fatalf("gap voltage %g V at %g s, above VBD", v[k], w.Time[k])
		}
		if v[k] > 2.5e3 {
			above = true
		} else if above {
			falls++
			above = false
		}
	}
	// Charging from near 0 V to 5 kV of a 10 kV supply takes RC·ln 2.
	if want := int(math.Floor(3e-3 / (1e-3 * math.Ln2))); falls < want-1 || falls > want+1 {
		t.Errorf("%d breakdowns in 3 ms, want about %d", falls, want)
	}
	l := w.Ledger
	if m := math.Abs(l.IntegrationMismatch) / l.EnergyIn; m > 5e-3 {
		t.Errorf("integrated energy mismatch %.3g of the input", m)
	}
}

func TestArcNetlistRoundTrip(t *testing.T) {
	c := deckCircuit(t, "arc\nV1 a 0 1\nR1 a b 1\nA1 b 0 habedank tauc=1u u0=40 taum=2u p0=3k vbd=1k\n.end\n")
	var buf bytes.Buffer
	if err := c.WriteNetlist(&buf); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(buf.String(), "A1 b 0 HABEDANK P0=3000 TAUC=1e-06 TAUM=2e-06 U0=40 VBD=1000") {
		t.Errorf("written line:\n%s", buf.String())
	}
	c2 := deckCircuit(t, buf.String())
	if got := circuit.ArcParamsOf(c2.Elements[2]); got.P0 != 3000 || got.Vbd != 1000 || got.Model != "habedank" {
		t.Errorf("re-read %+v", got)
	}
}

func TestArcValidation(t *testing.T) {
	for _, line := range []string{
		"A1 a 0 MAYR TAU=1u",           // no P0
		"A1 a 0 CASSIE TAU=1u P0=1",    // P0 not a Cassie parameter
		"A1 a 0 SPARK TAU=1u P0=1",     // unknown model
		"A1 a 0 MAYR TAU=-1u P0=1",     // negative τ
		"A1 a 0 MAYR TAU=1u P0=1 G0=0", // G0 below GOFF
	} {
		d, err := netlist.ParseString("t\nV1 a 0 1\n" + line + "\n.end\n")
		if err != nil {
			continue
		}
		if _, err := circuit.FromDeck(d); err == nil {
			t.Errorf("%q accepted", line)
		}
	}
}
