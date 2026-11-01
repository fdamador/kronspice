// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

func relDiff(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(math.Max(math.Abs(a), math.Abs(b)), 1e-300)
}

// Synthetic band-limited waveforms with random harmonics: the components
// are recovered, and the GA power equals the IEEE 1459 quantities computed
// conventionally (CLAUDE.md §3b rule 1: GA power vs IEEE-1459, 1e−9).
func TestGAPowerMatchesIEEE1459(t *testing.T) {
	r := rand.New(rand.NewSource(31))
	const f0, H, periods = 50.0, 9, 2
	for trial := 0; trial < 50; trial++ {
		var vc, ic [H][2]float64
		for k := 0; k < H; k++ {
			sc := 1 / float64(k+1)
			vc[k] = [2]float64{100 * sc * r.NormFloat64(), 100 * sc * r.NormFloat64()}
			ic[k] = [2]float64{10 * sc * r.NormFloat64(), 10 * sc * r.NormFloat64()}
		}
		vdc, idc := r.NormFloat64(), r.NormFloat64()
		// Uniform samples that land exactly on HarmonicsOf's grid.
		n := 4096
		T := periods / f0
		tm := make([]float64, n+1)
		v := make([]float64, n+1)
		i := make([]float64, n+1)
		for s := 0; s <= n; s++ {
			tt := T * float64(s) / float64(n)
			tm[s] = tt
			v[s], i[s] = vdc, idc
			for k := 0; k < H; k++ {
				sn, cs := math.Sincos(2 * math.Pi * float64(k+1) * f0 * tt)
				v[s] += math.Sqrt2 * (vc[k][0]*cs + vc[k][1]*sn)
				i[s] += math.Sqrt2 * (ic[k][0]*cs + ic[k][1]*sn)
			}
		}
		h, err := HarmonicsOf(tm, v, i, f0, periods, 16)
		if err != nil {
			t.Fatal(err)
		}
		for k := 0; k < H; k++ {
			for c := 0; c < 2; c++ {
				if math.Abs(h.V[k][c]-vc[k][c]) > 1e-11*100 || math.Abs(h.I[k][c]-ic[k][c]) > 1e-11*10 {
					t.Fatalf("harmonic %d component %d: %g, %g; want %g, %g", k+1, c, h.V[k][c], h.I[k][c], vc[k][c], ic[k][c])
				}
			}
		}
		g, err := GeometricPower(h)
		if err != nil {
			t.Fatal(err)
		}
		q, err := IEEE1459Of(tm, v, i, f0, periods)
		if err != nil {
			t.Fatal(err)
		}
		for _, c := range []struct {
			name      string
			got, want float64
		}{
			{"P", g.P, q.P}, {"S", g.S, q.S}, {"S1", g.S1, q.S1}, {"DI", g.DI, q.DI},
			{"DV", g.DV, q.DV}, {"SH", g.SH, q.SH}, {"P1", g.PerHarmonic[0].P, q.P1}, {"Q1", g.PerHarmonic[0].Q, q.Q1},
			{"S² = P²+Q²+D²", g.P*g.P + g.Q*g.Q + g.D*g.D, g.S * g.S},
			{"S² = S1²+DI²+DV²+SH²", q.S1*q.S1 + q.DI*q.DI + q.DV*q.DV + q.SH*q.SH, q.S * q.S},
			{"PDC", g.PDC, vdc * idc},
			{"Σ pair² = D²", sumSq(g.Distortion), g.D * g.D},
		} {
			if relDiff(c.got, c.want) > 1e-9 {
				t.Fatalf("%s: GA %.12g, IEEE 1459 %.12g", c.name, c.got, c.want)
			}
		}
	}
}

// A sinusoidal source into R–L, from the engine's transient: no
// distortion, P = I²R, Q = I²ωL (positive, inductive), and P, Q equal the
// engine's own AC analysis of the same circuit.
func TestGAPowerRL(t *testing.T) {
	const f0, R, L = 60.0, 10.0, 20e-3
	c, err := circuit.NewBuilder("RL").
		VoltageWave("V1", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 100, f0}}).
		Resistor("R1", "s", "a", R).Inductor("L1", "a", "0", L).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 1 / (f0 * 2000), Stop: 20 / f0, MaxStep: 1 / (f0 * 2000)}, TranOptions{Options: Options{RelTol: 1e-7}})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("s")
	iv, _ := wf.ElementCurrent("V1")
	i := make([]float64, len(iv))
	for k := range iv {
		i[k] = -iv[k] // current delivered into the load
	}
	h, err := HarmonicsOf(wf.Time, v, i, f0, 2, 8)
	if err != nil {
		t.Fatal(err)
	}
	g, _ := GeometricPower(h)
	w := 2 * math.Pi * f0
	irms := 100 / math.Sqrt2 / math.Hypot(R, w*L)
	if relDiff(g.P, irms*irms*R) > 1e-4 || relDiff(g.PerHarmonic[0].Q, irms*irms*w*L) > 1e-4 || g.D > 1e-3*g.S {
		t.Errorf("P %g (want %g), Q1 %g (want %g), D %g", g.P, irms*irms*R, g.PerHarmonic[0].Q, irms*irms*w*L, g.D)
	}
}

// A single-phase rectifier with a capacitor load: the engine's waveforms
// give GA and IEEE 1459 quantities that agree (to the harmonics above the
// 32nd, which the GA product leaves out), and the GA names the pair of
// harmonics behind the largest distortion power: fundamental voltage with
// third-harmonic current.
func TestGAPowerRectifier(t *testing.T) {
	const f0 = 60.0
	c, err := circuit.NewBuilder("rectifier").
		VoltageWave("V1", "s", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, 170, f0}}).
		Resistor("RS", "s", "a", 0.5).
		Model("dm", "D", map[string]float64{"is": 1e-12, "n": 1}).
		Diode("D1", "a", "p", "dm").Diode("D2", "0", "p", "dm").Diode("D3", "m", "a", "dm").Diode("D4", "m", "0", "dm").
		Capacitor("CL", "p", "m", 470e-6).Resistor("RL", "p", "m", 100).Resistor("RM", "m", "0", 1e6).Build()
	if err != nil {
		t.Fatal(err)
	}
	wf, err := Tran(c, TranSpec{Step: 1 / (f0 * 1000), Stop: 30 / f0, MaxStep: 1 / (f0 * 1000)}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := wf.NodeVoltage("s")
	iv, _ := wf.ElementCurrent("V1")
	i := make([]float64, len(iv))
	for k := range iv {
		i[k] = -iv[k]
	}
	h, err := HarmonicsOf(wf.Time, v, i, f0, 2, 32)
	if err != nil {
		t.Fatal(err)
	}
	g, _ := GeometricPower(h)
	q, _ := IEEE1459Of(wf.Time, v, i, f0, 2)
	if relDiff(g.P, q.P) > 1e-3 || relDiff(g.S, q.S) > 1e-2 || relDiff(g.S1, q.S1) > 1e-9 || relDiff(g.DI, q.DI) > 2e-2 {
		t.Errorf("GA P %g S %g S1 %g DI %g; IEEE 1459 P %g S %g S1 %g DI %g", g.P, g.S, g.S1, g.DI, q.P, q.S, q.S1, q.DI)
	}
	top := g.Distortion[0]
	if top.H != 1 || top.K != 3 {
		t.Errorf("largest distortion pair is harmonics %d–%d, want 1–3", top.H, top.K)
	}
	if !(q.THDI > 0.5) || q.PF > 0.9 {
		t.Errorf("rectifier THD_I %g, PF %g: expected a strongly distorted current", q.THDI, q.PF)
	}
	t.Logf("P %.1f W, S %.1f VA, PF %.3f, THD_I %.2f; GA: Q %.1f var, D %.1f VA, largest pair V%d×I%d %.1f VA",
		q.P, q.S, q.PF, q.THDI, g.Q, g.D, top.H, top.K, top.Value)
}

func sumSq(p []HarmonicPair) float64 {
	s := 0.0
	for _, x := range p {
		s += x.Value * x.Value
	}
	return s
}
