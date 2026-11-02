// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/kron"
)

// The unbalanced currents of the 1945 Fig. 15 follow from the potentials
// of Fig. 14 by eq. (1) on the lattice as built: this fixes the topology,
// the signs, y_c < 0 and the drive independently of any solve.
func TestFig15FromFig14(t *testing.T) {
	worst, at := 0.0, ""
	for j := 0; j < 5; j++ {
		for i := 0; i < 5; i++ {
			if Fig15Exception(j, i) {
				continue
			}
			d := math.Abs(UnbalancedCurrent(kron.Kron1945Potentials, j, i) - kron.Kron1945Residuals[j][i])
			if d > worst {
				worst, at = d, fmt.Sprintf("row %d col %d", j, i)
			}
		}
	}
	t.Logf("worst %.3f A at %s over 23 junctions", worst, at)
	if worst > 0.35 {
		t.Errorf("Fig. 15 not reproduced: %.3f A at %s", worst, at)
	}
	if d := UnbalancedCurrent(kron.Kron1945Potentials, 2, 3) + kron.Kron1945Residuals[2][3]; math.Abs(d) > 0.05 {
		t.Errorf("row 2 col 3 should be the printed value with its sign flipped: off by %.3f", d)
	}
}

// The given columns of Tables I–VI are the divided differences of
// kron.CurveFit on the printed z, and the printed z are the functions
// (the logarithmic ones rounded to four decimals).
func TestGivenColumns(t *testing.T) {
	p := kron.Kron1962Polyhedron()
	if w := worstAgainst(Given(p), 0); w > 3e-6 { // the worst is his eighth figure on the second differences
		t.Errorf("worst %.2g", w)
	}
	for f, fn := range Functions {
		for k, pt := range p.Points {
			if d := math.Abs(fn(pt[0], pt[1]) - Tables[f].F[k][0]); d > 2e-3 {
				t.Errorf("function %d at point %d: off by %.3g", f, k, d)
			}
		}
	}
}

// Frame 1 is the whole-number residual direction: every entry of the
// printed column follows, to the table's last figures.
func TestFrame1(t *testing.T) {
	p := kron.Kron1962Polyhedron()
	ys := Given(p)
	est := make([][13]float64, len(ys))
	for f, y := range ys {
		est[f], _ = p.Fit(kron.Kron1962Frame1, y)
	}
	if w := worstAgainst(est, 2); w > 1e-5 {
		t.Errorf("worst %.2g", w)
	}
}

// The mesh circuit, solved by MNA and by Kron's mesh form, gives the
// estimates as node voltages and dissipates the error-square: Kron p. 67
// on the engine's own ledger.
func TestMeshCircuit(t *testing.T) {
	p := kron.Kron1962Polyhedron()
	ys := Given(p)
	c, err := p.MeshCircuit(p.Residual(kron.Kron1962Frame1), ys)
	if err != nil {
		t.Fatal(err)
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	sum := 0.0
	for f, y := range ys {
		h, e := p.Fit(kron.Kron1962Frame1, y)
		sum += e
		for k := range h {
			v, ok := r.Voltage(fmt.Sprintf("u%d_%d", f, k))
			if !ok || math.Abs(v-h[k]) > 1e-9*(1+math.Abs(h[k])) {
				t.Errorf("function %d row %d: %g, want %g", f, k, v, h[k])
			}
		}
	}
	if d := math.Abs(r.Ledger.Dissipation-sum) / sum; d > 1e-9 {
		t.Errorf("dissipation %g, error-square %g", r.Ledger.Dissipation, sum)
	}
	nw, err := kron.Build(c)
	if err != nil {
		t.Fatal(err)
	}
	sol, err := nw.Solve(kron.MeshForm, 0, true)
	if err != nil {
		t.Fatal(err)
	}
	for f := range ys {
		h, _ := p.Fit(kron.Kron1962Frame1, ys[f])
		for k := range h {
			b, _ := nw.ElementBranch(fmt.Sprintf("R%d_%d", f, k))
			if d := math.Abs(math.Abs(real(sol.E.At(b))) - math.Abs(ys[f][k]-h[k])); d > 1e-9*(1+math.Abs(h[k])) {
				t.Errorf("Kron mesh form: function %d row %d voltage off by %g", f, k, d)
			}
		}
	}
}

// With the two misprints corrected the LS+linear column is one linear map
// of the four z values: the six tables span four dimensions, so two exact
// relations hold among any such column. As printed they do not.
func TestLinearColumnWithMisprints(t *testing.T) {
	corrected := Tables
	corrected[Misprints[0].Table].F[Misprints[0].Row][2] = Misprints[0].Value
	corrected[Misprints[1].Table].D[Misprints[1].Row-4][1] = Misprints[1].Value
	var zm [6][4]float64
	for f := range corrected {
		for k := 0; k < 4; k++ {
			zm[f][k] = corrected[f].F[k][0]
		}
	}
	null := leftNull(zm)
	if len(null) != 2 {
		t.Fatalf("left null space dimension %d", len(null))
	}
	check := func(tb [6]Table) float64 {
		w := 0.0
		for _, n := range null {
			for row := 0; row < 10; row++ {
				s := 0.0
				for f := range tb {
					s += n[f] * Column(tb[f], 1)[row]
				}
				w = math.Max(w, math.Abs(s))
			}
		}
		return w
	}
	if w := check(corrected); w > 1e-4 {
		t.Errorf("corrected column not linear: %.3g", w)
	}
	if w := check(Tables); w < 0.1 {
		t.Errorf("as printed the column should not be linear: %.3g", w)
	}
}

// leftNull returns an orthonormal basis of the left null space of m.
func leftNull(m [6][4]float64) [][6]float64 {
	var orth, basis [][6]float64
	add := func(v [6]float64) bool {
		for pass := 0; pass < 2; pass++ {
			for _, o := range orth {
				d := 0.0
				for i := range v {
					d += v[i] * o[i]
				}
				for i := range v {
					v[i] -= d * o[i]
				}
			}
		}
		n := 0.0
		for i := range v {
			n += v[i] * v[i]
		}
		n = math.Sqrt(n)
		if n < 1e-9 {
			return false
		}
		for i := range v {
			v[i] /= n
		}
		orth = append(orth, v)
		return true
	}
	for k := 0; k < 4; k++ {
		var c [6]float64
		for f := 0; f < 6; f++ {
			c[f] = m[f][k]
		}
		add(c)
	}
	for k := 0; k < 6; k++ {
		var e [6]float64
		e[k] = 1
		if add(e) {
			basis = append(basis, orth[len(orth)-1])
		}
	}
	return basis
}

// The GA, seeing only the printed column through the circuit, recovers
// frame 1 exactly.
func TestGARecoversFrame1(t *testing.T) {
	if testing.Short() {
		t.Skip("GA run")
	}
	p := kron.Kron1962Polyhedron()
	pr, err := FrameProblem(p, Given(p), 12, MatchObjective(2))
	if err != nil {
		t.Fatal(err)
	}
	st, err := Search(pr, 120, 150, 1, 0.2)
	if err != nil {
		t.Fatal(err)
	}
	if got := frameFromGenes(st.Best.Values); got != kron.Kron1962Frame1 {
		t.Errorf("best frame %v (fitness %.3g), want %v", got, st.Best.Fitness, kron.Kron1962Frame1)
	}
}

// The closed-form frame meets Kron's objective exactly: the circuit
// dissipates nothing for all six functions, the printed oscillatory column
// follows to its printed accuracy, and arbitrary numbers fit too (Kron
// p. 69: "the significance of this is not yet understood").
func TestOscillatoryFrame(t *testing.T) {
	p := kron.Kron1962Polyhedron()
	ys := Given(p)
	for _, m := range [][3]int{{0, 0, 0}, {2, 2, 2}, {-4, 6, 10}} {
		f, err := p.OscillatoryFrame(m)
		if err != nil {
			t.Fatal(err)
		}
		est := make([][13]float64, len(ys))
		for k, y := range ys {
			var e float64
			est[k], e = p.Fit(f, y)
			if e > 1e-18 {
				t.Errorf("frame %v function %d: error-square %.3g", f, k, e)
			}
		}
		if w := worstAgainst(est, 4); w > 3e-5 {
			t.Errorf("frame %v against the printed oscillatory column: %.3g", f, w)
		}
		if _, e := p.Fit(f, p.Characteristics([4]float64{3.7, -41.2, 18.9, 0.06})); e > 1e-18 {
			t.Errorf("frame %v arbitrary numbers: error-square %.3g", f, e)
		}
	}
	f, err := p.OscillatoryFrame([3]int{2, 2, 2})
	if err != nil {
		t.Fatal(err)
	}
	c, err := p.MeshCircuit(p.Residual(f), ys)
	if err != nil {
		t.Fatal(err)
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	if r.Ledger.Dissipation > 1e-15 {
		t.Errorf("circuit dissipation %.3g W", r.Ledger.Dissipation)
	}
	if _, err := p.OscillatoryFrame([3]int{1, 2, 2}); err == nil {
		t.Error("odd plane weight accepted")
	}
}

// The GA on Kron's objective improves on frame 1 by orders of magnitude
// but does not reach the exact lattice point; this records that, so a
// later search that does will show up here.
func TestGAMinimumPower(t *testing.T) {
	if testing.Short() {
		t.Skip("GA run")
	}
	p := kron.Kron1962Polyhedron()
	pr, err := FrameProblem(p, Given(p), 20, PowerObjective())
	if err != nil {
		t.Fatal(err)
	}
	st, err := Search(pr, 120, 150, 2, 0.04)
	if err != nil {
		t.Fatal(err)
	}
	if st.Best.Fitness > 1 {
		t.Errorf("best power %.3g W, frame 1 has 90 W", st.Best.Fitness)
	}
	t.Logf("best power %.3g W, frame %v", st.Best.Fitness, frameFromGenes(st.Best.Values))
}
