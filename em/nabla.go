// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/hyper/quat"
)

// ErrGrid is a grid too small or two fields on different grids.
var ErrGrid = errors.New("em: grid mismatch or too small")

// Grid is a uniform Cartesian grid: N points along x, y, z with spacing H,
// the first point at Origin. Point (i, j, k) is at Origin + (i, j, k)·H.
type Grid struct {
	N      [3]int
	H      [3]float64
	Origin [3]float64
}

// stencil is the half-width of the derivative stencil.
const stencil = 2

// Check verifies the grid has an interior (at least 5 points per axis) and
// positive spacing.
func (g Grid) Check() error {
	for a := 0; a < 3; a++ {
		if g.N[a] < 2*stencil+1 || !(g.H[a] > 0) {
			return fmt.Errorf("%w: %d points, spacing %g on axis %d", ErrGrid, g.N[a], g.H[a], a)
		}
	}
	return nil
}

// Len is the number of points.
func (g Grid) Len() int { return g.N[0] * g.N[1] * g.N[2] }

// Index of point (i, j, k), x fastest.
func (g Grid) Index(i, j, k int) int { return i + g.N[0]*(j+g.N[1]*k) }

// IJK is the inverse of Index.
func (g Grid) IJK(n int) (int, int, int) {
	i := n % g.N[0]
	n /= g.N[0]
	return i, n % g.N[1], n / g.N[1]
}

// Point returns the coordinates of point (i, j, k).
func (g Grid) Point(i, j, k int) [3]float64 {
	return [3]float64{g.Origin[0] + float64(i)*g.H[0], g.Origin[1] + float64(j)*g.H[1], g.Origin[2] + float64(k)*g.H[2]}
}

// Interior reports whether the derivative stencil fits around point n.
func (g Grid) Interior(n int) bool {
	i, j, k := g.IJK(n)
	idx := [3]int{i, j, k}
	for a := 0; a < 3; a++ {
		if idx[a] < stencil || idx[a] >= g.N[a]-stencil {
			return false
		}
	}
	return true
}

// QField is a quaternion-valued field on a grid (a vector field when every
// scalar part is 0).
type QField struct {
	G Grid
	Q []quat.Q
}

// SampleQ samples f at every grid point.
func SampleQ(g Grid, f func(x, y, z float64) quat.Q) QField {
	out := QField{G: g, Q: make([]quat.Q, g.Len())}
	for n := range out.Q {
		p := g.Point(g.IJK(n))
		out.Q[n] = f(p[0], p[1], p[2])
	}
	return out
}

// SampleS samples a scalar field.
func SampleS(g Grid, f func(x, y, z float64) float64) []float64 {
	out := make([]float64, g.Len())
	for n := range out {
		p := g.Point(g.IJK(n))
		out[n] = f(p[0], p[1], p[2])
	}
	return out
}

// d4 is the fourth-order central difference from the four neighbours.
func d4(m2, m1, p1, p2, h float64) float64 { return (m2 - 8*m1 + 8*p1 - p2) / (12 * h) }

// stride is the index step along axis a.
func (g Grid) stride(a int) int {
	switch a {
	case 0:
		return 1
	case 1:
		return g.N[0]
	}
	return g.N[0] * g.N[1]
}

var nanQ = quat.Q{W: math.NaN(), X: math.NaN(), Y: math.NaN(), Z: math.NaN()}

// Partial returns ∂f/∂x_a (a = 0, 1, 2) at interior points, NaN at the
// edges.
func Partial(f QField, a int) QField {
	g := f.G
	out := QField{G: g, Q: make([]quat.Q, g.Len())}
	s := g.stride(a)
	h := g.H[a]
	for n := range out.Q {
		if !g.Interior(n) {
			out.Q[n] = nanQ
			continue
		}
		m2, m1, p1, p2 := f.Q[n-2*s], f.Q[n-s], f.Q[n+s], f.Q[n+2*s]
		out.Q[n] = quat.Q{
			W: d4(m2.W, m1.W, p1.W, p2.W, h),
			X: d4(m2.X, m1.X, p1.X, p2.X, h),
			Y: d4(m2.Y, m1.Y, p1.Y, p2.Y, h),
			Z: d4(m2.Z, m1.Z, p1.Z, p2.Z, h),
		}
	}
	return out
}

// Nabla is Hamilton's ∇f = i·∂f/∂x + j·∂f/∂y + k·∂f/∂z (quaternion
// products, ∇ on the left). For a vector field its scalar part is −div f
// and its vector part curl f; for a scalar field (W only) its vector part
// is grad f.
func Nabla(f QField) QField {
	dx, dy, dz := Partial(f, 0), Partial(f, 1), Partial(f, 2)
	i, j, k := quat.Vec(1, 0, 0), quat.Vec(0, 1, 0), quat.Vec(0, 0, 1)
	out := QField{G: f.G, Q: make([]quat.Q, f.G.Len())}
	for n := range out.Q {
		out.Q[n] = i.Mul(dx.Q[n]).Add(j.Mul(dy.Q[n])).Add(k.Mul(dz.Q[n]))
	}
	return out
}

// TimeDerivative is ∂f/∂t at the middle time from snapshots at t − 2δ,
// t − δ, t + δ, t + 2δ (fourth order).
func TimeDerivative(m2, m1, p1, p2 QField, dt float64) (QField, error) {
	for _, f := range []QField{m1, p1, p2} {
		if f.G != m2.G || len(f.Q) != len(m2.Q) {
			return QField{}, ErrGrid
		}
	}
	out := QField{G: m2.G, Q: make([]quat.Q, len(m2.Q))}
	for n := range out.Q {
		a, b, c, d := m2.Q[n], m1.Q[n], p1.Q[n], p2.Q[n]
		out.Q[n] = quat.Q{W: d4(a.W, b.W, c.W, d.W, dt), X: d4(a.X, b.X, c.X, d.X, dt), Y: d4(a.Y, b.Y, c.Y, d.Y, dt), Z: d4(a.Z, b.Z, c.Z, d.Z, dt)}
	}
	return out, nil
}

// Scalars wraps a scalar field as a quaternion field (W part).
func Scalars(g Grid, s []float64) QField {
	out := QField{G: g, Q: make([]quat.Q, len(s))}
	for n, v := range s {
		out.Q[n] = quat.Q{W: v}
	}
	return out
}
