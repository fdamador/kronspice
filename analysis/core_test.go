// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// Jiles & Atherton 1986 (J. Magn. Magn. Mater. 61, 48, Fig. 4 parameters
// for a ferromagnet): Ms = 1.6e6 A/m, a = 1100 A/m, k = 400 A/m, c = 0.2,
// α = 1.6e-3. Core: Ae = 1 cm², le = 10 cm, N = 100, so H = 1000·i.
const jaCore = "U1 a 0 MS=1.6meg A=1100 K=400 C=0.2 ALPHA=1.6m AREA=1e-4 LEN=0.1 N1=100"

const (
	jaMs, jaA, jaK, jaC, jaAlpha = 1.6e6, 1100.0, 400.0, 0.2, 1.6e-3
	jaArea, jaLen, jaN           = 1e-4, 0.1, 100.0
	mu0Test                      = 1.25663706212e-6
)

// jaRef integrates the Jiles–Atherton equations along a given H path with
// classical RK4 in many small steps: a reference written apart from the
// device (Jiles & Atherton 1986 eqs. 1–8; the M, Mirr form of Jiles,
// Thoelke & Devine 1992; Lederer et al. 1999 zero-slope rule).
type jaRef struct{ h, m, mirr float64 }

func (r *jaRef) slopes(h, m, mirr, delta float64) (float64, float64) {
	x := (h + jaAlpha*m) / jaA
	var l, dl float64
	if math.Abs(x) < 1e-3 {
		l, dl = x/3, 1.0/3
	} else {
		l, dl = 1/math.Tanh(x)-1/x, 1/(x*x)-1/math.Pow(math.Sinh(x), 2)
	}
	man, dman := jaMs*l, jaMs/jaA*dl
	dirr := 0.0
	if delta*(man-mirr) > 0 {
		dirr = (man - mirr) / (delta * jaK)
	}
	chi := (1-jaC)*dirr + jaC*dman
	den := 1 - jaAlpha*chi
	return chi / den, dirr / den
}

func (r *jaRef) to(h1 float64) {
	delta := math.Copysign(1, h1-r.h)
	n := int(math.Ceil(math.Abs(h1-r.h)/(jaA/200))) + 1
	dh := (h1 - r.h) / float64(n)
	for s := 0; s < n; s++ {
		h, m, mi := r.h, r.m, r.mirr
		a1, b1 := r.slopes(h, m, mi, delta)
		a2, b2 := r.slopes(h+dh/2, m+dh/2*a1, mi+dh/2*b1, delta)
		a3, b3 := r.slopes(h+dh/2, m+dh/2*a2, mi+dh/2*b2, delta)
		a4, b4 := r.slopes(h+dh, m+dh*a3, mi+dh*b3, delta)
		r.m += dh / 6 * (a1 + 2*a2 + 2*a3 + a4)
		r.mirr += dh / 6 * (b1 + 2*b2 + 2*b3 + b4)
		r.h += dh
	}
}

// A sinusoidal winding current drives the core round its major loop. The
// flux density from the winding voltage, B = ∫v dt/(N·Ae), follows the
// reference integration; the loop is symmetric with coercivity and
// remanence; and the heat booked by the ledger over the second cycle equals
// the loop area ∮H dB times the core volume.
func TestCoreHysteresisLoop(t *testing.T) {
	const ipk, f = 10.0, 50.0
	c := deckCircuit(t, "core\nI1 0 a SIN(0 10 50)\n"+jaCore+"\n.end\n")
	w := runTran(t, c, 20e-6, 2/f, 20e-6)
	v, _ := w.NodeVoltage("a")
	ref := &jaRef{}
	b := 0.0
	worst, bmax := 0.0, 0.0
	hs, bs := make([]float64, len(w.Time)), make([]float64, len(w.Time))
	for k, tt := range w.Time {
		if k > 0 {
			b += 0.5 * (v[k] + v[k-1]) * (tt - w.Time[k-1]) / (jaN * jaArea)
		}
		h := jaN * ipk * math.Sin(2*math.Pi*f*tt) / jaLen
		ref.to(h)
		bref := mu0Test * (h + ref.m)
		worst = math.Max(worst, math.Abs(b-bref))
		bmax = math.Max(bmax, math.Abs(bref))
		hs[k], bs[k] = h, b
	}
	if worst > 0.01*bmax {
		t.Errorf("B differs from the reference by %.3g T (peak %.3g T)", worst, bmax)
	}
	// Coercivity and remanence on the second cycle.
	var hc, br []float64
	area := 0.0
	for k := 1; k < len(hs); k++ {
		if w.Time[k-1] < 1/f {
			continue
		}
		if (bs[k-1] < 0) != (bs[k] < 0) {
			hc = append(hc, hs[k-1]-bs[k-1]*(hs[k]-hs[k-1])/(bs[k]-bs[k-1]))
		}
		if (hs[k-1] < 0) != (hs[k] < 0) {
			br = append(br, bs[k-1]-hs[k-1]*(bs[k]-bs[k-1])/(hs[k]-hs[k-1]))
		}
		area += 0.5 * (hs[k] + hs[k-1]) * (bs[k] - bs[k-1])
	}
	if len(hc) != 2 || len(br) < 1 {
		t.Fatalf("crossings: Hc %v, Br %v", hc, br)
	}
	if hc[0]*hc[1] >= 0 || math.Abs(math.Abs(hc[0])-math.Abs(hc[1])) > 0.03*math.Abs(hc[0]) || math.Abs(hc[0]) < jaK/4 {
		t.Errorf("coercivities %v A/m", hc)
	}
	if math.Abs(br[0]) < 0.1 {
		t.Errorf("remanences %v T", br)
	}
	// Heat on the second cycle, from the ledger.
	l := w.Ledger
	heat := 0.0
	for k := 1; k < len(w.Time); k++ {
		if w.Time[k-1] >= 1/f {
			heat += 0.5 * (l.Dissipation[k] + l.Dissipation[k-1]) * (w.Time[k] - w.Time[k-1])
		}
	}
	vol := jaArea * jaLen
	if want := area * vol; !(want > 0) || math.Abs(heat-want) > 0.02*want {
		t.Errorf("heat per cycle %.4g J, loop area × volume %.4g J", heat, want)
	}
	t.Logf("Hc %.0f A/m, Br %.3f T, Bpk %.3f T, loss %.3g J/cycle (%.3g J/m³)", hc[0], br[0], bmax, heat, area)
}

// At the demagnetized operating point the winding is an inductor of
// N²·Ae·μ0·(1 + χ0)/le with the initial susceptibility
// χ0 = (c·Ms/3a)/(1 − α·c·Ms/3a).
func TestCoreSmallSignalInductance(t *testing.T) {
	c := deckCircuit(t, "core\nI1 0 a DC 0 AC 1\n"+jaCore+"\n.end\n")
	r, err := AC(c, ACSpec{Variation: "lin", Points: 1, FStart: 1e3, FStop: 1e3}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := r.NodeVoltage("a")
	x := jaC * jaMs / (3 * jaA)
	l := jaN * jaN * jaArea * mu0Test * (1 + x/(1-jaAlpha*x)) / jaLen
	if got, want := cmplx.Abs(v[0]), 2*math.Pi*1e3*l; math.Abs(got-want) > 1e-9*want {
		t.Errorf("|Z| = %g Ω, want ωL = %g Ω", got, want)
	}
}

// Two windings share one flux: an open secondary carries N2/N1 of the
// primary's induced voltage, through saturation and hysteresis alike.
func TestCoreTransformerRatio(t *testing.T) {
	c := deckCircuit(t, `xfmr
V1 s 0 SIN(0 20 1k)
R1 s p 1
U1 p 0 q 0 MS=1.6meg A=1100 K=400 C=0.2 ALPHA=1.6m AREA=1e-4 LEN=0.1 N1=100 N2=25
RL q 0 1meg
.end
`)
	w := runTran(t, c, 2e-6, 3e-3, 2e-6)
	vp, _ := w.NodeVoltage("p")
	vq, _ := w.NodeVoltage("q")
	worst, peak := 0.0, 0.0
	for k := range vp {
		worst = math.Max(worst, math.Abs(vq[k]-0.25*vp[k]))
		peak = math.Max(peak, math.Abs(vp[k]))
	}
	if worst > 1e-3*peak {
		t.Errorf("secondary differs from N2/N1·v1 by %.3g V (primary peak %.3g V)", worst, peak)
	}
}

func TestCoreValidation(t *testing.T) {
	for _, line := range []string{
		"U1 a 0 A=1100 K=400 AREA=1e-4 LEN=0.1 N1=100",                     // no MS
		"U1 a 0 MS=1.6meg A=1100 K=400 AREA=1e-4 LEN=0.1",                  // no N1
		"U1 a 0 b 0 MS=1.6meg A=1100 K=400 AREA=1e-4 LEN=0.1 N1=100",       // no N2
		"U1 a 0 MS=1.6meg A=1100 K=400 ALPHA=0.01 AREA=1e-4 LEN=0.1 N1=10", // α·Ms/3a > 1
		"U1 a 0 MS=1.6meg A=1100 K=400 C=1 AREA=1e-4 LEN=0.1 N1=10",        // c = 1
		"U1 a 0 MS=1.6meg A=1100 K=400 AREA=1e-4 LEN=0.1 N1=10 N3=4",       // stray turns
	} {
		d, err := netlist.ParseString("t\nV1 a 0 1\nR9 b 0 1\n" + line + "\n.end\n")
		if err != nil {
			continue
		}
		if _, err := circuit.FromDeck(d); err == nil {
			t.Errorf("%q accepted", line)
		}
	}
}
