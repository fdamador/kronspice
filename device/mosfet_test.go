// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// Terminal layout: nodes 1 = d, 2 = s, 3 = g, 4 = b; unknowns 0..3, then
// the internal d', s'.
func testMOS(t testing.TB, typ string, params map[string]float64, inst map[string]float64) (*MOSFET, int) {
	m := circuit.Model{Name: "m", Type: typ, Params: params}
	q, err := NewMOSFET(circuit.Element{Name: "M1", Kind: circuit.MOSFET, Nodes: [2]int{1, 2}, Extra: []int{3, 4},
		Model: "m", Params: inst}, m)
	if err != nil {
		t.Fatal(err)
	}
	ctx := NewSetupContext(5)
	if err := q.Setup(ctx); err != nil {
		t.Fatal(err)
	}
	return q, ctx.Size()
}

var mosModel = map[string]float64{"tox": 20e-9, "nsub": 3e16, "uo": 550, "ld": 0.1e-6, "lambda": 0.03,
	"rd": 15, "rs": 10, "cgso": 2e-10, "cgdo": 2e-10, "cgbo": 1e-10, "cj": 4e-4, "mj": 0.45, "cjsw": 3e-10,
	"mjsw": 0.3, "js": 1e-4}
var mosInst = map[string]float64{"l": 1e-6, "w": 20e-6, "ad": 1e-10, "as": 1e-10, "pd": 30e-6, "ps": 30e-6}

// mosResidual loads q linearized exactly at x and returns F(x) = A·x − b
// and A.
func mosResidual(q *MOSFET, n int, x []float64, mode Mode) ([]float64, []float64) {
	vgs, vds, vbs, vgb := q.voltages(x)
	q.acc.vgs, q.acc.vds, q.acc.vbs = vgs, vds, vbs // no limiting
	m := &denseStub{n: n, a: make([]float64, n*n)}
	rhs := make([]float64, n)
	_ = q.Load(&LoadContext{M: m, RHS: rhs, X: x, Mode: mode, H: 1e-10, Method: Trapezoidal, Gmin: 1e-12})
	_ = vgb
	f := make([]float64, n)
	for i := 0; i < n; i++ {
		s := -rhs[i]
		for j := 0; j < n; j++ {
			s += m.a[i*n+j] * x[j]
		}
		f[i] = s
	}
	return f, m.a
}

// The stamps are the Jacobian of the device's currents (and, in a
// transient step, its charge currents) in every region, for NMOS and PMOS,
// with drain and source in either order. Meyer's model ignores ∂C/∂v by
// definition, so the check is made where the gate voltages equal their
// previous values (the ∂C/∂v·Δv term vanishes there).
func TestMOSJacobian(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	withJunctions := func(m map[string]float64) map[string]float64 {
		out := map[string]float64{"rd": 15, "rs": 10, "cgso": 2e-10, "cgdo": 2e-10, "cgbo": 1e-10, "cj": 4e-4,
			"mj": 0.45, "cjsw": 3e-10, "mjsw": 0.3, "js": 1e-4}
		for k, v := range m {
			out[k] = v
		}
		return out
	}
	models := []map[string]float64{mosModel, withJunctions(mosLevelModels[1].params),
		withJunctions(mosLevelModels[2].params), withJunctions(mosLevelModels[5].params)}
	for c := 0; c < 2*len(models); c++ {
		typ := []string{"nmos", "pmos"}[c%2]
		q, n := testMOS(t, typ, models[c/2], mosInst)
		sgn := q.typ
		for trial := 0; trial < 60; trial++ {
			vgs := -1 + 5*rng.Float64()
			vds := -3 + 6*rng.Float64()
			vbs := -3 + 3.3*rng.Float64()
			// For a forward-biased bulk (normal-mode vbs > 0) SPICE keeps
			// gmbs = gm·γ/(2·sarg), not the exact derivative of its
			// linearly continued threshold; TestMOSForwardBodyGmbs covers it.
			vbsN := vbs // bulk–source voltage in the normal-mode frame
			if vds < 0 {
				vbsN = vbs - vds
			}
			if vbsN > 0 {
				continue
			}
			x := make([]float64, n)
			x[q.sp] = 0.01 * rng.NormFloat64()
			x[q.dp] = x[q.sp] + sgn*vds
			x[q.g] = x[q.sp] + sgn*vgs
			x[q.b] = x[q.sp] + sgn*vbs
			x[q.u[0]] = x[q.dp] + 0.01*rng.NormFloat64()
			x[q.u[1]] = x[q.sp] + 0.01*rng.NormFloat64()
			vg, vd, vb, vgb := q.voltages(x)
			q.vPrev = [5]float64{vg, vg - vd, vgb, vb, vb - vd}
			q.qPrev = [5]float64{1e-15, -2e-15, 1e-16, 3e-15, -1e-15}
			q.iqPrev = [5]float64{1e-6, -2e-6, 1e-7, 0, 1e-6}
			for _, mode := range []Mode{ModeOP, ModeTran} {
				_, a := mosResidual(q, n, x, mode)
				for j := 0; j < n; j++ {
					h := 1e-7
					xp, xm := append([]float64(nil), x...), append([]float64(nil), x...)
					xp[j] += h
					xm[j] -= h
					fp, _ := mosResidual(q, n, xp, mode)
					fm, _ := mosResidual(q, n, xm, mode)
					for i := 0; i < n; i++ {
						rowMax := 0.0
						for k := 0; k < n; k++ {
							rowMax = math.Max(rowMax, math.Abs(a[i*n+k]))
						}
						num := (fp[i] - fm[i]) / (2 * h)
						// Region boundaries (cutoff, saturation) make the
						// derivative one-sided within h of them.
						if d := math.Abs(num - a[i*n+j]); d > 1e-5*math.Abs(num)+1e-7*rowMax+1e-12 {
							_, _, _, _, von, vdsat := q.ids(vgs, math.Abs(vds), vbs)
							if math.Abs(vgs-von) < 1e-5 || math.Abs(math.Abs(vds)-vdsat) < 1e-5 || math.Abs(vds) < 1e-5 {
								continue
							}
							t.Fatalf("%s mode %d vgs %.3f vds %.3f vbs %.3f: dF[%d]/dx[%d] = %g, stamp %g",
								typ, mode, vgs, vds, vbs, i, j, num, a[i*n+j])
						}
					}
				}
			}
		}
	}
}

// PMOS is the mirror of NMOS: negate every voltage and current.
func TestMOSPolarityMirror(t *testing.T) {
	qn, n := testMOS(t, "nmos", mosModel, mosInst)
	qp, _ := testMOS(t, "pmos", map[string]float64{"tox": 20e-9, "vto": -qn.p.VTO, "gamma": qn.p.GAMMA,
		"phi": qn.p.PHI, "kp": qn.p.KP, "ld": 0.1e-6, "lambda": 0.03, "rd": 15, "rs": 10, "js": 1e-4}, mosInst)
	for _, x := range [][]float64{{3, 0.1, 2.5, -1, 2.9, 0.12}, {-1, 1, 3, 0, -0.9, 0.95}} {
		fn, _ := mosResidual(qn, n, x, ModeOP)
		xm := make([]float64, len(x))
		for i := range x {
			xm[i] = -x[i]
		}
		fp, _ := mosResidual(qp, n, xm, ModeOP)
		for i := range fn {
			if math.Abs(fn[i]+fp[i]) > 1e-15+1e-12*math.Abs(fn[i]) {
				t.Fatalf("x %v unknown %d: nmos %g, pmos %g", x, i, fn[i], fp[i])
			}
		}
	}
}

// Level 1 regions: cutoff, square law in saturation, triode, and the body
// effect on the threshold.
func TestMOSLevel1Regions(t *testing.T) {
	q, _ := testMOS(t, "nmos", map[string]float64{"vto": 0.8, "kp": 50e-6, "lambda": 0.02, "gamma": 0.4, "phi": 0.65},
		map[string]float64{"l": 2e-6, "w": 20e-6})
	beta := 50e-6 * 20e-6 / 2e-6
	cases := []struct{ vgs, vds, vbs, want float64 }{
		{0.5, 3, 0, 0},
		{2, 3, 0, beta / 2 * 1.2 * 1.2 * (1 + 0.02*3)},
		{2, 0.5, 0, beta * 0.5 * (1.2 - 0.25) * (1 + 0.02*0.5)},
	}
	for _, c := range cases {
		id, _, _, _, _, _ := q.ids(c.vgs, c.vds, c.vbs)
		if math.Abs(id-c.want) > 1e-15+1e-12*c.want {
			t.Errorf("%+v: id %g, want %g", c, id, c.want)
		}
	}
	_, _, _, _, von, _ := q.ids(2, 3, -2)
	if want := 0.8 + 0.4*(math.Sqrt(2.65)-math.Sqrt(0.65)); math.Abs(von-want) > 1e-12 {
		t.Errorf("threshold at vbs = −2 V: %g, want %g", von, want)
	}
}

// Process parameters: KP = UO·Cox, φ = 2·vt·ln(NSUB/ni), γ = √(2·ε_Si·q·NSUB)/Cox
// (Vladimirescu & Liu 1980, Table 2).
func TestMOSProcessParameters(t *testing.T) {
	q, _ := testMOS(t, "nmos", map[string]float64{"tox": 20e-9, "nsub": 3e16, "uo": 550}, nil)
	cox := 3.9 * 8.854214871e-12 / 20e-9
	vt := BoltzmanK * TNomK / ChargeQ
	checks := map[string][2]float64{
		"KP":    {q.p.KP, 550e-4 * cox},
		"PHI":   {q.p.PHI, 2 * vt * math.Log(3e22/1.45e16)},
		"GAMMA": {q.p.GAMMA, math.Sqrt(2*11.7*8.854214871e-12*ChargeQ*3e22) / cox},
	}
	for name, v := range checks {
		if math.Abs(v[0]-v[1]) > 1e-12*math.Abs(v[1]) {
			t.Errorf("%s = %g, want %g", name, v[0], v[1])
		}
	}
	// An n-channel device on a 3e16 cm⁻³ substrate with an n+ poly gate
	// has a positive threshold of order 0.5–1 V.
	if q.p.VTO < 0.3 || q.p.VTO > 1.2 {
		t.Errorf("VTO = %g V", q.p.VTO)
	}
}

// With the bulk forward-biased (vbs > 0), SPICE continues √(φ − vbs) by its
// tangent, sarg = √φ − vbs/(2√φ), and keeps gmbs = gm·γ/(2·sarg). That is
// SPICE's (and ngspice's) small-signal body transconductance there, kept
// for compatibility; the exact derivative of the continued threshold would
// be gm·γ/(2√φ).
func TestMOSForwardBodyGmbs(t *testing.T) {
	q, _ := testMOS(t, "nmos", map[string]float64{"vto": 0.8, "kp": 50e-6, "gamma": 0.4, "phi": 0.64}, nil)
	_, gm, _, gmbs, _, _ := q.ids(2.5, 3, 0.2)
	sarg := 0.8 - 0.2/1.6
	if want := gm * 0.4 / (2 * sarg); math.Abs(gmbs-want) > 1e-15*want {
		t.Fatalf("gmbs %g, want SPICE's %g", gmbs, want)
	}
}

func TestFetlimLimvds(t *testing.T) {
	// fetlim: from off, a large rise is cut at vto + 0.5; within the middle
	// region a rise is cut at vto + 4; going down from strong inversion is
	// cut at vto + 2.
	if v := fetlim(10, 0, 1); v != 1.5 {
		t.Errorf("fetlim off→on: %g", v)
	}
	if v := fetlim(10, 2, 1); v != 5 {
		t.Errorf("fetlim middle rise: %g", v)
	}
	if v := fetlim(0, 6, 1); v != 3 {
		t.Errorf("fetlim strong inversion fall: %g", v)
	}
	if v := limvds(20, 1); v != 4 {
		t.Errorf("limvds low rise: %g", v)
	}
	if v := limvds(30, 5); v != 17 {
		t.Errorf("limvds high rise: %g", v)
	}
	if v := limvds(-5, 1); v != -0.5 {
		t.Errorf("limvds fall: %g", v)
	}
}

// Meyer capacitances are continuous across their region boundaries and
// sum to Cox in strong inversion at vds = 0.
func TestMeyerContinuity(t *testing.T) {
	cox, phi, von := 1e-12, 0.7, 0.8
	f := func(vgs, vds float64) [3]float64 {
		vdsat := math.Max(vgs-von, 0)
		a, b, c := meyer(vgs, vgs-vds, von, vdsat, phi, cox)
		return [3]float64{a, b, c}
	}
	for _, edge := range []float64{von - phi, von - phi/2, von} {
		lo, hi := f(edge-1e-9, 0.5), f(edge+1e-9, 0.5)
		for k := 0; k < 3; k++ {
			if math.Abs(lo[k]-hi[k]) > 1e-20 {
				t.Errorf("discontinuity at vgs = %g: %v vs %v", edge, lo, hi)
			}
		}
	}
	if c := f(3, 0); math.Abs(c[0]+c[1]-cox) > 1e-24 {
		t.Errorf("cgs + cgd at vds = 0: %g, want %g", c[0]+c[1], cox)
	}
	if c := f(3, 5); math.Abs(c[0]-2*cox/3) > 1e-24 || c[1] != 0 {
		t.Errorf("saturation: %v, want cgs = 2/3 Cox, cgd = 0", c)
	}
}

func BenchmarkMOSLoad(b *testing.B) {
	for _, bc := range []struct {
		name  string
		model map[string]float64
	}{{"level1", mosModel}, {"level2", mosLevelModels[2].params}, {"level3", mosLevelModels[6].params}} {
		b.Run(bc.name, func(b *testing.B) {
			q, n := testMOS(b, "nmos", bc.model, mosInst)
			x := make([]float64, n)
			copy(x, []float64{3, 0.1, 2.5, -1, 2.9, 0.12})
			m := &denseStub{n: n, a: make([]float64, n*n)}
			ctx := &LoadContext{M: m, RHS: make([]float64, n), X: x, Mode: ModeTran, H: 1e-10, Method: Trapezoidal, Gmin: 1e-12}
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				_ = q.Load(ctx)
			}
		})
	}
}
