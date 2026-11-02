// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
)

// CurveFit is the worked polyhedron of G. Kron, "Multidimensional
// curve-fitting with self-organizing automata," J. Math. Anal. Appl. 5,
// 46–69, 1962, Fig. 3 (p. 56): four points on a plane, six branches (the
// 1-network), three planes (the 2-network) and the least-square model (the
// −1 network), giving 4 + 6 + 3 = 13 characteristics, 12 regression
// coefficients and one residual degree of freedom (p. 57). Branch order
// and direction and plane order are those of the rows of Tables I–VI,
// found from the "given" columns, which they reproduce to 1e-6.
// See doc/kron-figure-catalogue.md §3.
type CurveFit struct {
	Points   [4][2]float64
	Branches [6][2]int // from → to
	// Planes lists each plane's branch opposite the interior vertex 1 and
	// the plane's area.
	Planes [3]struct {
		Opposite int
		Area     float64
	}
}

// Kron1962Polyhedron returns the polyhedron of Fig. 3. Vertex 1 (index 0)
// is interior; the three planes 1-2-4, 1-2-3, 1-3-4 each have area 4.5.
func Kron1962Polyhedron() *CurveFit {
	p := &CurveFit{
		Points:   [4][2]float64{{1, 1}, {1, -2}, {-2, 3}, {4, 2}},
		Branches: [6][2]int{{1, 0}, {1, 2}, {1, 3}, {0, 2}, {0, 3}, {2, 3}},
	}
	for k, b := range []int{2, 1, 5} {
		p.Planes[k].Opposite = b
		p.Planes[k].Area = 4.5
	}
	return p
}

// Length is the length of branch b.
func (p *CurveFit) Length(b int) float64 {
	u, v := p.Points[p.Branches[b][0]], p.Points[p.Branches[b][1]]
	return math.Hypot(v[0]-u[0], v[1]-u[1])
}

// Characteristics is the vector y of step 2 of the algorithm (pp. 54–55):
// the four function values, the six first divided differences (Δz/L along
// each branch) and the three second divided differences. The paper does
// not write the last formula; the tables fix it: S = −(Δz/L)·(A/L) on the
// plane's branch opposite the interior vertex. All 13 are linear in the
// four values z, so y lies in a 4-dimensional subspace of R¹³.
func (p *CurveFit) Characteristics(z [4]float64) [13]float64 {
	var y [13]float64
	copy(y[:4], z[:])
	dd := func(b int) float64 { return (z[p.Branches[b][1]] - z[p.Branches[b][0]]) / p.Length(b) }
	for b := 0; b < 6; b++ {
		y[4+b] = dd(b)
	}
	for k, pl := range p.Planes {
		y[10+k] = -dd(pl.Opposite) * pl.Area / p.Length(pl.Opposite)
	}
	return y
}

// PolyFrame is a reference frame of the polyhedron, given by the left null
// vector of its 13 × 12 model, the residual direction w. With one degree
// of freedom the least-square fit is ŷ = y − w·(wᵀy)/(wᵀw), whichever 12
// columns span the model. The function part of w is the hyperplane's own
// residual direction (−3, 1, 1, 1); the rest are whole numbers: w_b =
// N[b]/L_b on the branches and w_p = M[p]/(2A) on the planes.
type PolyFrame struct {
	N [6]int
	M [3]int
}

// Kron1962Frame1 is the frame of the first "two nonoscillatory polyhedra"
// column of Tables I–VI. It is recovered from that column (the paper
// states no frame): the column's error map has rank one to the table's
// rounding, and its direction is these integers.
var Kron1962Frame1 = PolyFrame{N: [6]int{4, -1, -3, -4, -4, 2}, M: [3]int{2, -10, -4}}

// OscillatoryFrame returns a frame whose error vanishes for every set of
// values at the four points: wᵀy = 0 on the whole 4-dimensional range of
// Characteristics, so ŷ = y and the power input is zero whatever the
// function, which is Kron's objective (minimum stored energy, p. 67) met
// exactly. The inner branches close the hyperplane's residual
// z₂ + z₃ + z₄ − 3z₁ with weights N = ±L² (a Kirchhoff sum at the interior
// vertex), and each outer branch is balanced against its plane by
// N = M/2. M must be even; M = 0 is allowed, and the family is why the
// printed oscillatory column cannot say which member Kron used.
func (p *CurveFit) OscillatoryFrame(m [3]int) (PolyFrame, error) {
	f := PolyFrame{M: m}
	for b := 0; b < 6; b++ {
		l2 := int(math.Round(p.Length(b) * p.Length(b)))
		switch {
		case p.Branches[b][1] == 0:
			f.N[b] = l2
		case p.Branches[b][0] == 0:
			f.N[b] = -l2
		}
	}
	for k, pl := range p.Planes {
		if m[k]%2 != 0 {
			return f, fmt.Errorf("%w: plane weight %d is odd", ErrPolyhedron, m[k])
		}
		f.N[pl.Opposite] = m[k] / 2
	}
	return f, nil
}

// Residual returns w for the frame.
func (p *CurveFit) Residual(f PolyFrame) [13]float64 {
	w := [13]float64{-3, 1, 1, 1}
	for b := 0; b < 6; b++ {
		w[4+b] = float64(f.N[b]) / p.Length(b)
	}
	for k, pl := range p.Planes {
		w[10+k] = float64(f.M[k]) / (2 * pl.Area)
	}
	return w
}

// Fit is the fit in the frame: the estimate ŷ and the error-square
// (y − ŷ)ᵀ(y − ŷ), which Kron identifies with the network's power input
// (p. 67).
func (p *CurveFit) Fit(f PolyFrame, y [13]float64) (yhat [13]float64, errSq float64) {
	w := p.Residual(f)
	var wy, ww float64
	for k := range w {
		wy += w[k] * y[k]
		ww += w[k] * w[k]
	}
	s := wy / ww
	for k := range w {
		yhat[k] = y[k] - s*w[k]
	}
	return yhat, s * s * ww
}

// MeshCircuit is the fit as one circuit per function, in Kron's mesh form
// with C = w (i = C·i′, Z′ = Cₜ·Z·C, e′ = Cₜ·e; Kron 1939 Ch. IV): each of
// the 13 primitive branches is an impressed voltage y_k and 1 Ω, joined to
// the single mesh through an ideal transformer of ratio w_k (a VCVS in the
// mesh and a CCCS on the branch, the "L = M = ∞" transformers of Kron 1962
// Fig. 1). The mesh current is i′ = wᵀy/wᵀw, the branch current and the
// resistor's voltage are the error w_k·i′, and the node voltage "u<f>_<k>"
// is the estimate ŷ_k. The resistors dissipate exactly the error-square.
func (p *CurveFit) MeshCircuit(w [13]float64, ys [][13]float64) (*circuit.Circuit, error) {
	b := circuit.NewBuilder("Kron 1962 polyhedron: one mesh per function, C = w")
	for f, y := range ys {
		mesh := "0"
		for k := 0; k < 13; k++ {
			u, v := fmt.Sprintf("u%d_%d", f, k), fmt.Sprintf("v%d_%d", f, k)
			b.VoltageSource(fmt.Sprintf("VY%d_%d", f, k), v, "0", y[k])
			b.Resistor(fmt.Sprintf("R%d_%d", f, k), v, u, 1)
			b.CCCS(fmt.Sprintf("FT%d_%d", f, k), u, "0", fmt.Sprintf("VM%d", f), w[k])
			next := fmt.Sprintf("m%d_%d", f, k)
			b.VCVS(fmt.Sprintf("ET%d_%d", f, k), next, mesh, u, "0", w[k])
			mesh = next
		}
		b.VoltageSource(fmt.Sprintf("VM%d", f), mesh, "0", 0)
	}
	return b.Build()
}
