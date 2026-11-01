// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/circuit"
)

// Terminal layout used by these tests: nodes 1 = c, 2 = e, 3 = b, 4 = s, so
// unknowns 0..3, then the internal c', b', e'.
func testBJT(t testing.TB, typ string, params map[string]float64) (*BJT, int) {
	m := circuit.Model{Name: "q", Type: typ, Params: params}
	q := NewBJT(circuit.Element{Name: "Q1", Kind: circuit.BJT, Nodes: [2]int{1, 2}, Extra: []int{3, 4}, Model: "q"}, m)
	ctx := NewSetupContext(5)
	if err := q.Setup(ctx); err != nil {
		t.Fatal(err)
	}
	return q, ctx.Size()
}

var fullModel = map[string]float64{"is": 3e-16, "bf": 120, "nf": 1.02, "vaf": 70, "ikf": 20e-3, "ise": 4e-14,
	"ne": 1.6, "br": 3, "nr": 1.01, "var": 30, "ikr": 5e-3, "isc": 2e-14, "nc": 1.8, "rb": 80, "re": 1.5,
	"rc": 7, "cje": 5e-12, "mje": 0.35, "tf": 0.3e-9, "xtf": 2, "vtf": 3, "itf": 30e-3, "cjc": 3e-12,
	"mjc": 0.4, "xcjc": 0.7, "tr": 5e-9, "cjs": 2e-12, "mjs": 0.3}

// residual loads q linearized at x and returns F(x) = A·x − b, the current
// leaving each unknown's node into the device, and A.
func residual(q *BJT, n int, x []float64, mode Mode) ([]float64, []float64) {
	vbe, vbc, _, _ := q.voltages(x)
	q.vbeAcc, q.vbcAcc = vbe, vbc // no limiting: linearize exactly at x
	m := &denseStub{n: n, a: make([]float64, n*n)}
	rhs := make([]float64, n)
	_ = q.Load(&LoadContext{M: m, RHS: rhs, X: x, Mode: mode, H: 1e-10, Method: Trapezoidal, Gmin: 1e-12})
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

// The analytic stamps are the Jacobian of the device's own currents (and,
// in a transient step, of its trapezoidal charge currents), checked by
// central differences in every operating region, for NPN and PNP.
func TestBJTJacobian(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	for _, typ := range []string{"npn", "pnp"} {
		q, n := testBJT(t, typ, fullModel)
		sgn := 1.0
		if typ == "pnp" {
			sgn = -1
		}
		q.qPrev = [4]float64{1e-13, -2e-13, -1e-13, 3e-13}
		q.iqPrev = [4]float64{1e-4, -2e-4, 1e-5, -1e-5}
		for trial := 0; trial < 40; trial++ {
			// Terminal and internal voltages around each region.
			vbe := -1 + 1.8*rng.Float64()
			vbc := -6 + 6.8*rng.Float64()
			x := make([]float64, n)
			x[q.ep] = rng.NormFloat64() * 0.01
			x[q.bp] = x[q.ep] + sgn*vbe
			x[q.cp] = x[q.bp] - sgn*vbc
			x[q.u[0]] = x[q.cp] + 0.01*rng.NormFloat64()
			x[q.u[1]] = x[q.ep] + 0.01*rng.NormFloat64()
			x[q.b] = x[q.bp] + 0.01*rng.NormFloat64()
			x[q.s] = -3 * sgn
			for _, mode := range []Mode{ModeOP, ModeTran} {
				_, a := residual(q, n, x, mode)
				for j := 0; j < n; j++ {
					h := 1e-7
					xp, xm := append([]float64(nil), x...), append([]float64(nil), x...)
					xp[j] += h
					xm[j] -= h
					fp, _ := residual(q, n, xp, mode)
					fm, _ := residual(q, n, xm, mode)
					for i := 0; i < n; i++ {
						// Finite-difference noise scales with the row's
						// largest entry (round-off in F over 2h).
						rowMax := 0.0
						for k := 0; k < n; k++ {
							rowMax = math.Max(rowMax, math.Abs(a[i*n+k]))
						}
						num := (fp[i] - fm[i]) / (2 * h)
						if d := math.Abs(num - a[i*n+j]); d > 1e-5*math.Abs(num)+1e-8*rowMax {
							t.Fatalf("%s mode %d vbe %.3f vbc %.3f: dF[%d]/dx[%d] = %g, stamp %g",
								typ, mode, vbe, vbc, i, j, num, a[i*n+j])
						}
					}
				}
			}
		}
	}
}

// A PNP is the mirror image of an NPN: negate every voltage and every
// current.
func TestBJTPolarityMirror(t *testing.T) {
	// The mirror of a vertical NPN is a vertical PNP (ngspice makes a PNP
	// lateral by default: substrate junction to the base).
	pm := map[string]float64{"subs": 1}
	for k, v := range fullModel {
		pm[k] = v
	}
	qn, n := testBJT(t, "npn", fullModel)
	qp, _ := testBJT(t, "pnp", pm)
	x := []float64{3, 0.1, 0.8, -2, 2.95, 0.79, 0.12}
	fn, _ := residual(qn, n, x, ModeOP)
	xm := make([]float64, n)
	for i := range x {
		xm[i] = -x[i]
	}
	fp, _ := residual(qp, n, xm, ModeOP)
	for i := range fn {
		if math.Abs(fn[i]+fp[i]) > 1e-15+1e-12*math.Abs(fn[i]) {
			t.Fatalf("unknown %d: npn %g, pnp %g", i, fn[i], fp[i])
		}
	}
}

// With every second-order effect off, the model is Ebers–Moll (transport
// form): ic = IS(e^(vbe/vt) − e^(vbc/vt)) − IS/BR(e^(vbc/vt) − 1), plus gmin
// (for junction voltages above −3·vt; below, SPICE's cubic reverse tail
// replaces the exponential, as in the diode).
func TestBJTEbersMollLimit(t *testing.T) {
	q, _ := testBJT(t, "npn", map[string]float64{"is": 1e-15, "bf": 80, "br": 2})
	q.gmin = 1e-12
	vt := q.vt
	for _, v := range [][2]float64{{0.65, -4}, {0.7, 0.3}, {-0.05, -0.07}, {0.5, 0.6}} {
		d := q.evalDC(v[0], v[1])
		ef, er := math.Exp(v[0]/vt)-1, math.Exp(v[1]/vt)-1
		ic := 1e-15*(ef-er) - 1e-15/2*er - 1e-12*v[1]
		ib := 1e-15/80*ef + 1e-15/2*er + 1e-12*(v[0]+v[1])
		if got := d.it - d.ibc; math.Abs(got-ic) > 1e-12*math.Abs(ic)+1e-24 {
			t.Errorf("vbe %g vbc %g: ic %g, want %g", v[0], v[1], got, ic)
		}
		if got := d.ibe + d.ibc; math.Abs(got-ib) > 1e-12*math.Abs(ib)+1e-24 {
			t.Errorf("vbe %g vbc %g: ib %g, want %g", v[0], v[1], got, ib)
		}
	}
}

// Capacitances (including the forward transit charge's dependence on vbc)
// are the derivatives of the charges.
func TestBJTChargeDerivatives(t *testing.T) {
	q, _ := testBJT(t, "npn", fullModel)
	q.gmin = 1e-12
	for _, v := range [][4]float64{{0.7, -2, -2.1, -3}, {0.3, 0.5, 0.45, 1}, {-0.5, -5, -5, -1}, {0.75, 0.2, 0.3, 0.2}} {
		ch := func(a, b, c, e float64) bjtCharge {
			d := q.evalDC(a, b)
			return q.charges(a, b, c, e, &d)
		}
		c0 := ch(v[0], v[1], v[2], v[3])
		h := 1e-7
		for k := 0; k < 4; k++ {
			vp, vm := v, v
			vp[k] += h
			vm[k] -= h
			num := (ch(vp[0], vp[1], vp[2], vp[3]).q[k] - ch(vm[0], vm[1], vm[2], vm[3]).q[k]) / (2 * h)
			if math.Abs(num-c0.c[k]) > 1e-5*math.Abs(num)+1e-20 {
				t.Errorf("v %v: dq%d/dv = %g, c = %g", v, k, num, c0.c[k])
			}
		}
		vp, vm := v, v
		vp[1] += h
		vm[1] -= h
		num := (ch(vp[0], vp[1], vp[2], vp[3]).q[qBE] - ch(vm[0], vm[1], vm[2], vm[3]).q[qBE]) / (2 * h)
		if math.Abs(num-c0.dqbeDvbc) > 1e-5*math.Abs(num)+1e-20 {
			t.Errorf("v %v: dqbe/dvbc = %g, model %g", v, num, c0.dqbeDvbc)
		}
	}
}

// For charges that depend only on their own voltage, StoredEnergy is the
// exact ∫₀^v u·C(u) du (checked here on the emitter depletion charge below
// and above FC·VJ).
func TestBJTStoredEnergyDepletion(t *testing.T) {
	q, _ := testBJT(t, "npn", map[string]float64{"cje": 4e-12, "vje": 0.7, "mje": 0.4})
	for _, v := range []float64{-2, 0.2, 0.6} {
		q.vSol = [4]float64{v, 0, 0, 0}
		got := q.StoredEnergy(OperatingPoint{})
		// Midpoint rule on u·C(u), 2·10⁵ panels.
		want, n := 0.0, 200000
		for i := 0; i < n; i++ {
			u := v * (float64(i) + 0.5) / float64(n)
			_, c := depletion(u, 4e-12, 0.7, 0.4, 0.5, q.f1e, q.f2e, q.f3e)
			want += u * c * v / float64(n)
		}
		if math.Abs(got-want) > 1e-8*math.Abs(want) {
			t.Errorf("v = %g: stored %g, want %g", v, got, want)
		}
	}
}

func BenchmarkBJTLoad(b *testing.B) {
	q, n := testBJT(b, "npn", fullModel)
	x := []float64{3, 0.1, 0.8, -2, 2.95, 0.79, 0.12}
	m := &denseStub{n: n, a: make([]float64, n*n)}
	ctx := &LoadContext{M: m, RHS: make([]float64, n), X: x, Mode: ModeTran, H: 1e-10, Method: Trapezoidal, Gmin: 1e-12}
	q.hasCharge = true
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = q.Load(ctx)
	}
}
