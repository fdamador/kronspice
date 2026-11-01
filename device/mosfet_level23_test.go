// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
)

// Level 2 and 3 model variants exercising every branch: plain, weak
// inversion (NFS), velocity saturation (VMAX), short and narrow channel
// (XJ, DELTA), given LAMBDA, mobility degradation (UCRIT/UEXP, THETA) and
// static feedback (ETA).
var mosLevelModels = []struct {
	name   string
	params map[string]float64
}{
	{"l2-plain", map[string]float64{"level": 2, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6}},
	{"l2-full", map[string]float64{"level": 2, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6,
		"xj": 0.3e-6, "delta": 1.2, "ucrit": 5e4, "uexp": 0.15, "nfs": 1e11, "neff": 2}},
	{"l2-vmax", map[string]float64{"level": 2, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6,
		"xj": 0.3e-6, "vmax": 5e4, "neff": 1.5, "ucrit": 5e4, "uexp": 0.1}},
	{"l2-lambda", map[string]float64{"level": 2, "vto": 0.8, "kp": 5e-5, "gamma": 0.5, "phi": 0.65, "lambda": 0.02}},
	{"l3-plain", map[string]float64{"level": 3, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6}},
	{"l3-full", map[string]float64{"level": 3, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6,
		"xj": 0.3e-6, "delta": 0.8, "theta": 0.08, "eta": 0.6, "kappa": 0.5, "nfs": 5e11, "wd": 0.1e-6}},
	{"l3-vmax", map[string]float64{"level": 3, "tox": 25e-9, "nsub": 2e16, "uo": 600, "ld": 0.1e-6,
		"xj": 0.3e-6, "theta": 0.05, "eta": 0.3, "vmax": 1.5e5, "kappa": 0.3}},
}

var mosLevelInst = map[string]float64{"l": 2e-6, "w": 10e-6}

// gm, gds and gmbs are the derivatives of the drain current (central
// differences), away from region boundaries, for every variant.
func TestMOSLevel23Derivatives(t *testing.T) {
	rng := rand.New(rand.NewSource(5))
	for _, mv := range mosLevelModels {
		q, _ := testMOS(t, "nmos", mv.params, mosLevelInst)
		checked := 0
		for trial := 0; trial < 400; trial++ {
			v := [3]float64{-0.5 + 5*rng.Float64(), 5 * rng.Float64(), -4 + 4.4*rng.Float64()}
			if trial%10 == 0 {
				v[1] = 1e-3 * rng.Float64() // near vds = 0
			}
			id, gm, gds, gmbs, _, _ := q.ids(v[0], v[1], v[2])
			if math.IsNaN(id) || math.IsInf(id, 0) {
				t.Fatalf("%s at %v: current %g", mv.name, v, id)
			}
			for k, g := range [3]float64{gm, gds, gmbs} {
				h := 1e-6
				vp, vm := v, v
				vp[k] += h
				vm[k] -= h
				ip, _, _, _, _, _ := q.ids(vp[0], vp[1], vp[2])
				im, _, _, _, _, _ := q.ids(vm[0], vm[1], vm[2])
				// One-sided slopes that disagree mark a region boundary
				// (or the vds ≤ 1e−10 special case) inside ±h.
				fwd, bwd := (ip-id)/h, (id-im)/h
				if vm[1] <= 1e-10 || math.Abs(fwd-bwd) > 1e-3*math.Max(math.Abs(fwd), math.Abs(bwd))+1e-12 {
					continue
				}
				num := (ip - im) / (2 * h)
				// The level-2 velocity-saturated vdsat comes from a closed-form
				// quartic root whose round-off shows in tiny currents: 1e−4·|id|.
				if math.Abs(num-g) > 1e-5*math.Abs(num)+1e-4*math.Abs(id)+1e-12 {
					t.Errorf("%s at vgs %.4f vds %.4f vbs %.4f: d id/d v[%d] = %g, model %g", mv.name, v[0], v[1], v[2], k, num, g)
				}
				checked++
			}
		}
		if checked < 600 {
			t.Errorf("%s: only %d derivatives checked", mv.name, checked)
		}
	}
}

// Expected level-2 and level-3 behavior: cutoff below threshold, current
// rising with vgs and vds, a weak-inversion tail with NFS, velocity
// saturation lowering the current, and the body effect raising the
// threshold.
func TestMOSLevel23Behavior(t *testing.T) {
	for _, mv := range mosLevelModels {
		t.Run(mv.name, func(t *testing.T) {
			q, _ := testMOS(t, "nmos", mv.params, mosLevelInst)
			id := func(vgs, vds, vbs float64) float64 { c, _, _, _, _, _ := q.ids(vgs, vds, vbs); return c }
			_, _, _, _, von, _ := q.ids(0, 1, 0)
			if !(von > 0.2 && von < 2) {
				t.Fatalf("von %g", von)
			}
			if _, _, _, _, vonB, _ := q.ids(0, 1, -2); !(vonB > von) {
				t.Errorf("body effect: von %g at vbs = −2 V, %g at 0", vonB, von)
			}
			prev := 0.0
			for _, vgs := range []float64{von + 0.2, von + 0.5, von + 1, von + 2} {
				c := id(vgs, 3, 0)
				if !(c > prev) {
					t.Errorf("id(vgs = %.2f) = %g not above %g", vgs, c, prev)
				}
				prev = c
			}
			prev = 0
			for _, vds := range []float64{0.05, 0.2, 0.5, 1, 2, 4} {
				c := id(von+1.5, vds, 0)
				if !(c > prev) {
					t.Errorf("id(vds = %.2f) = %g not above %g", vds, c, prev)
				}
				prev = c
			}
			below := id(von-0.2, 2, 0)
			if mv.params["nfs"] != 0 {
				if !(below > 0 && below < id(von+0.2, 2, 0)*0.1) {
					t.Errorf("weak inversion %g", below)
				}
			} else if below != 0 && !(mv.params["level"] == 2 && mv.params["vmax"] > 0) {
				// Level 2 without NFS cuts off at vbin, not von; between them
				// the Baum–Beneking vdsat (VMAX) leaves a small current, as
				// in SPICE. Otherwise vdsat = 0 there and the current is 0.
				t.Errorf("current %g below threshold without NFS", below)
			}
			if mv.params["vmax"] > 0 {
				p := map[string]float64{}
				for k, v := range mv.params {
					if k != "vmax" {
						p[k] = v
					}
				}
				q0, _ := testMOS(t, "nmos", p, mosLevelInst)
				c0, _, _, _, _, _ := q0.ids(von+2, 4, 0)
				if c := id(von+2, 4, 0); !(c < c0) {
					t.Errorf("velocity saturation: %g with VMAX, %g without", c, c0)
				}
			}
			t.Logf("von %.3f V, id(von+1.5, 2) = %s", von, fmt.Sprintf("%.4g A", id(von+1.5, 2, 0)))
		})
	}
}
