// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kron_polyhedron works through Kron's field network and his 1962
// polyhedron (register entry H20; doc/kron-figure-catalogue.md):
//
//  1. The 1945 cavity (J. Appl. Phys. 16, 172, Figs. 14–15): the lattice
//     solved exactly against the potentials he read on the network
//     analyzer, and his published unbalanced currents recomputed from them.
//
//  2. The waveguide of Balasubramanian's thesis (Figs. 15–16): matching
//     impedance against Δz/λ.
//
//  3. The curve-fitting tables (J. Math. Anal. Appl. 5, 46, 1962, Tables
//     I–VI): the first planar frame recovered as whole numbers, and the fit
//     solved as a circuit in Kron's mesh form, the error-square as its
//     dissipation.
//
//  4. The genetic algorithm, on the same circuit: rediscovering frame 1
//     from the printed column alone, then Kron's stated objective, minimum
//     power input, which it approaches but does not reach.
//
//  5. The oscillatory frame in closed form, zero power for every function.
//
//     go run ./examples/kron_polyhedron
package main

import (
	"context"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"time"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/evo"
	"github.com/fdamador/kronspice/kron"
)

// UnbalancedCurrent is eq. (1) of the 1945 paper at junction (row, col) of
// Fig. 14: i_u = Σ(e_k − e₀)y_L − e₀y_c + I, with y_L = 1.25, y_c = −0.77
// (his "acapacitance, a negative number") and the impressed current I.
func UnbalancedCurrent(e [5][5]float64, row, col int) float64 {
	const yl, yc = 1.25, -0.77
	sum := 0.0
	for _, d := range [][2]int{{-1, 0}, {1, 0}, {0, -1}, {0, 1}} {
		r, c := row+d[0], col+d[1]
		if r >= 0 && r < 5 && c >= 0 && c < 5 {
			sum += (e[r][c] - e[row][col]) * yl
		}
	}
	i := 0.0
	if row == 0 {
		i = 50 * math.Cos(math.Pi*float64(col)/4)
	}
	return sum - e[row][col]*yc + i
}

// Fig15Exception names the two printed entries of Fig. 15 that do not
// follow from Fig. 14: row 2 col 3, printed 2.15 where eq. (1) gives −2.16
// and the figure's own e_u is −0.51 (a dropped minus sign); and row 1 col 0,
// printed 2.12 (e_u 0.5) where the printed potentials give −5.94, which
// would need a potential near 22.95 for the 25.66 printed there.
func Fig15Exception(row, col int) bool { return (row == 1 && col == 0) || (row == 2 && col == 3) }

// Column gathers a table column as 13 characteristics; col is 0 given,
// 1 LS+L (first 10 rows only), 2 frame 1, 3 frame 2, 4 oscillatory.
func Column(t Table, col int) [13]float64 {
	var y [13]float64
	fc := [5]int{0, 2, 3, 4, 5}[col]
	for k := 0; k < 4; k++ {
		y[k] = t.F[k][fc]
	}
	for k := 0; k < 6; k++ {
		y[4+k] = t.D[k][col]
	}
	if col == 0 {
		for k := 0; k < 3; k++ {
			y[10+k] = t.S[k][0]
		}
	} else if col >= 2 {
		for k := 0; k < 3; k++ {
			y[10+k] = t.S[k][col-1]
		}
	}
	return y
}

// Given returns the characteristics of the six functions from the z
// values as printed in the tables: the fits start from them, and the
// logarithmic functions are rounded there to four decimals.
func Given(p *kron.CurveFit) [][13]float64 {
	out := make([][13]float64, len(Tables))
	for f := range Tables {
		var z [4]float64
		for k := range z {
			z[k] = Tables[f].F[k][0]
		}
		out[f] = p.Characteristics(z)
	}
	return out
}

// worstAgainst is the largest |estimate − printed| over all functions and
// rows of a table column.
func worstAgainst(est [][13]float64, col int) float64 {
	w := 0.0
	for f := range est {
		y := Column(Tables[f], col)
		for k := range y {
			w = math.Max(w, math.Abs(est[f][k]-y[k]))
		}
	}
	return w
}

// frameFromGenes maps the nine signed integer genes onto a frame.
func frameFromGenes(v []float64) kron.PolyFrame {
	var f kron.PolyFrame
	for k := 0; k < 6; k++ {
		f.N[k] = int(math.Round(v[k]))
	}
	for k := 0; k < 3; k++ {
		f.M[k] = int(math.Round(v[6+k]))
	}
	return f
}

// FrameProblem is the search over integer frames, every candidate solved
// as the mesh circuit by the circuit engine. Magnitudes run 1..maxN on
// the branches and 1..2·maxN on the planes.
func FrameProblem(p *kron.CurveFit, ys [][13]float64, maxN float64, objective evo.Objective) (evo.Problem, error) {
	var space evo.Space
	for k := 0; k < 6; k++ {
		space = append(space, evo.Gene{Name: fmt.Sprintf("N%d", k), Min: 1, Max: maxN, Sign: evo.Either, Integer: true})
	}
	for k := 0; k < 3; k++ {
		space = append(space, evo.Gene{Name: fmt.Sprintf("M%d", k), Min: 1, Max: 2 * maxN, Sign: evo.Either, Integer: true})
	}
	build := func(v []float64) (*circuit.Circuit, error) {
		return p.MeshCircuit(p.Residual(frameFromGenes(v)), ys)
	}
	_, err := build([]float64{1, 1, 1, 1, 1, 1, 1, 1, 1})
	if err != nil {
		return evo.Problem{}, err
	}
	return evo.Problem{Build: build, Space: space, Objectives: []evo.Objective{objective}}, nil
}

// PowerObjective is the objective Kron states (p. 67): the network's
// total dissipation, here the error-square summed over the functions.
func PowerObjective() evo.Objective {
	return evo.Objective{Name: "power input", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		return r.Ledger.Dissipation, nil
	}}
}

// MatchObjective is the mean squared distance of the circuit's node
// voltages from a printed column.
func MatchObjective(col int) evo.Objective {
	return evo.Objective{Name: "printed column", Weight: 1, Eval: func(e *evo.Eval) (float64, error) {
		r, err := e.OP()
		if err != nil {
			return 0, err
		}
		s, n := 0.0, 0
		for f := range Tables {
			y := Column(Tables[f], col)
			for k := range y {
				v, ok := r.Voltage(fmt.Sprintf("u%d_%d", f, k))
				if !ok {
					return 0, fmt.Errorf("no node u%d_%d", f, k)
				}
				s += (v - y[k]) * (v - y[k])
				n++
			}
		}
		return s / float64(n), nil
	}}
}

// Search runs the GA. sigma is the mutation step in decades: on whole
// numbers near 10 a step of 0.03 is about one unit, while the default 0.2
// jumps several and rarely lands on an exact lattice point.
func Search(p evo.Problem, pop, gens int, seed uint64, sigma float64) (*evo.State, error) {
	return evo.Run(context.Background(), p, evo.Options{Population: pop, Generations: gens, Seed: seed, Sigma: sigma, Mutation: 0.3, Timeout: 2 * time.Minute}, nil, nil)
}

func main() {
	fmt.Println("1. Kron 1945 Fig. 14: the 5 × 5 cavity lattice, ω = 1, L = 0.8, C = 1/1.3")
	c, err := kron.Kron1945Cavity(1e9)
	if err != nil {
		log.Fatal(err)
	}
	nw, err := kron.Build(c)
	if err != nil {
		log.Fatal(err)
	}
	sol, err := nw.Solve(kron.NodeForm, 1, false)
	if err != nil {
		log.Fatal(err)
	}
	var exact [5][5]float64
	worst, top := 0.0, 0.0
	for i := 0; i < 5; i++ {
		for j := 0; j < 5; j++ {
			b, _ := nw.ElementBranch(fmt.Sprintf("C_%d_%d", i, j))
			exact[j][i] = imag(sol.E.At(b))
		}
	}
	sign := 1.0
	if exact[0][0]*kron.Kron1945Potentials[0][0] < 0 {
		sign = -1
	}
	for j := 0; j < 5; j++ {
		for i := 0; i < 5; i++ {
			top = math.Max(top, math.Abs(exact[j][i]))
			worst = math.Max(worst, math.Abs(sign*exact[j][i]-kron.Kron1945Potentials[j][i]))
		}
	}
	fmt.Printf("   exact H_y against his analyzer values: worst %.2f of %.1f (%.1f %%)\n", worst, top, 100*worst/top)
	wr := 0.0
	for j := 0; j < 5; j++ {
		for i := 0; i < 5; i++ {
			if !Fig15Exception(j, i) {
				wr = math.Max(wr, math.Abs(UnbalancedCurrent(kron.Kron1945Potentials, j, i)-kron.Kron1945Residuals[j][i]))
			}
		}
	}
	fmt.Printf("   Fig. 15 recomputed from his potentials by his eq. (1): worst %.2f A over 23 junctions; two printed entries disagree (row 1 col 0, row 2 col 3)\n\n", wr)

	fmt.Println("2. Balasubramanian's waveguide: 11 cm guide, 2 GHz, five cells across")
	mu0 := 1.25663706212e-6
	eps0 := 1 / (mu0 * 299792458.0 * 299792458.0)
	w := kron.Waveguide{A: 0.11, Dy: 0.01, N: 5, Sections: 24, Freq: 2e9, Eps: eps0, Mu: mu0}
	lambda := 299792458.0 / w.Freq
	w0 := w
	w0.Dz = 1e-5
	for _, dz := range []float64{0.0075, 0.011, 0.0165, 0.022} {
		w.Dz = dz
		cw, _, err := w.Circuit(w.LatticeImpedance())
		if err != nil {
			log.Fatal(err)
		}
		n2, err := kron.Build(cw)
		if err != nil {
			log.Fatal(err)
		}
		s2, err := n2.Solve(kron.NodeForm, 2*math.Pi*w.Freq, false)
		if err != nil {
			log.Fatal(err)
		}
		lo, hi := math.Inf(1), 0.0
		for j := 3; j <= w.Sections-3; j++ {
			b, _ := n2.ElementBranch(fmt.Sprintf("C_2_%d", j))
			m := cmplx.Abs(s2.E.At(b))
			lo, hi = math.Min(lo, m), math.Max(hi, m)
		}
		fmt.Printf("   Δz/λ = %4.1f %%: longitudinal error %+.2f %%, ripple at that match %.1e\n", 100*dz/lambda, 100*(w.LatticeImpedance()/w0.LatticeImpedance()-1), (hi-lo)/(hi+lo))
	}
	fmt.Println("   thesis Fig. 15, read from the curve: about 1.0, 1.6, 3.4, 6.2 %")
	fmt.Println()

	fmt.Println("3. Kron 1962 Tables I–VI")
	p := kron.Kron1962Polyhedron()
	ys := Given(p)
	fmt.Printf("   given column, from the printed z: worst %.1e\n", worstAgainst(ys, 0))
	est := make([][13]float64, len(ys))
	sum := 0.0
	for f, y := range ys {
		var e float64
		est[f], e = p.Fit(kron.Kron1962Frame1, y)
		sum += e
	}
	fmt.Printf("   frame 1 from the integers N = %v, M = %v: worst %.1e against the printed column\n", kron.Kron1962Frame1.N, kron.Kron1962Frame1.M, worstAgainst(est, 2))
	wo := 0.0
	for f := range Tables {
		a, b := Column(Tables[f], 4), Column(Tables[f], 0)
		for k := range a {
			wo = math.Max(wo, math.Abs(a[k]-b[k]))
		}
	}
	fmt.Printf("   printed oscillatory column against the given column: worst %.1e\n", wo)
	cf, err := p.MeshCircuit(p.Residual(kron.Kron1962Frame1), ys)
	if err != nil {
		log.Fatal(err)
	}
	r, err := analysis.OP(cf, analysis.Options{})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   mesh circuit: dissipation %.6f W, summed error-square %.6f; ledger residual %.1e W\n\n", r.Ledger.Dissipation, sum, r.Ledger.Residual)

	fmt.Println("4. The genetic algorithm: integer genes, every candidate solved by the circuit engine")
	pm, err := FrameProblem(p, ys, 12, MatchObjective(2))
	if err != nil {
		log.Fatal(err)
	}
	st, err := Search(pm, 120, 150, 1, 0.2)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   against the printed frame-1 column: best %+v, fitness %.2e (%d evaluations)\n", frameFromGenes(st.Best.Values), st.Best.Fitness, st.Evaluations)
	pp, err := FrameProblem(p, ys, 20, PowerObjective())
	if err != nil {
		log.Fatal(err)
	}
	st, err = Search(pp, 120, 150, 2, 0.04)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   minimum power input by the GA: best %+v, power %.2e W against frame 1's %.1f W; it does not reach zero (see README)\n", frameFromGenes(st.Best.Values), st.Best.Fitness, sum)
	fmt.Println()
	fmt.Println("5. Zero power in closed form: w orthogonal to the range of the characteristics")
	best, err := p.OscillatoryFrame([3]int{2, 2, 2})
	if err != nil {
		log.Fatal(err)
	}
	co, err := p.MeshCircuit(p.Residual(best), ys)
	if err != nil {
		log.Fatal(err)
	}
	ro, err := analysis.OP(co, analysis.Options{})
	if err != nil {
		log.Fatal(err)
	}
	for f, y := range ys {
		est[f], _ = p.Fit(best, y)
	}
	fmt.Printf("   frame %+v: circuit dissipation %.1e W; against the printed oscillatory column worst %.1e\n", best, ro.Ledger.Dissipation, worstAgainst(est, 4))
	rng := rand.New(rand.NewPCG(7, 11))
	wrand := 0.0
	for n := 0; n < 1000; n++ {
		var z [4]float64
		for k := range z {
			z[k] = 100 * (rng.Float64() - 0.5)
		}
		y := p.Characteristics(z)
		h, _ := p.Fit(best, y)
		for k := range y {
			wrand = math.Max(wrand, math.Abs(h[k]-y[k]))
		}
	}
	fmt.Printf("   1000 sets of random numbers at the four points: worst fit error %.1e\n", wrand)
}
