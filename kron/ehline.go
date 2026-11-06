// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
)

// Kron's field network with the Euler–Heisenberg vacuum as a nonlinear,
// state-dependent primitive network, for fields uniform across a beam.
//
// Reduction. For fields that do not vary across x and y (plane waves
// along z, in a box with periodic sides one cell wide), Kron's
// three-dimensional cell (MaxwellNetwork) keeps only its transverse
// branches: per cell along z an x link and a y link (capacitors, voltage
// E·Δz, charge D·A) and the faces normal to x and to y (inductors, flux
// B·A′, current H·Δ). They form two LC ladders, one per polarization, the
// reduction the coaxial TEM pulse already used (TestMaxwellNetworkCoaxPulse).
//
// Nonlinear primitive. In the Euler–Heisenberg vacuum (em.EHConstitutive)
// each capacitor's charge depends on its own voltage, on the other
// polarization's and on the neighbouring faces' fluxes, and each face's
// current on its flux and on the links beside it:
//
//	q = A·D(E, B̄),   I = Δ·H(B, Ē),
//
// B̄ and Ē averaged to the element from its two neighbours (the Yee
// staggering). The connection tensor, the ladder's wiring, stays fixed:
// Kron's network with a state-dependent primitive tensor (Kron 1951
// Epilogue; kron/nonlinear.go), integrated in time by the leapfrog of the
// lattice field (gauge.go, em/fdtd.go):
//
//	H^{n+½} = H(B^{n+½}, Ē^n)                       (faces, explicit)
//	D^{n+1} = D^n − Δt·∂H/∂z                       (Ampère: KCL at each link)
//	E^{n+1} from D(E^{n+1}, B̄^{n+½}) = D^{n+1}     (Newton per cell)
//	B^{n+3/2} = B^{n+½} − Δt·∂E^{n+1}/∂z           (Faraday: KVL round each face)
//
// in the 1-D curl: (∇×F)_x = −∂F_y/∂z, (∇×F)_y = ∂F_x/∂z. The line is a
// ring (periodic in z), so there are no ends to reflect from.
//
// With Linear set the elements are ε₀ and μ₀: the plain lattice, the
// standard (Maxwell) prediction beside every nonlinear run.

// ErrEHLine reports an invalid set-up or a cell whose field left the
// weak-field range during a run.
var ErrEHLine = errors.New("kron: Euler–Heisenberg line")

// EHLine is a ring of N cells of length Dz.
type EHLine struct {
	N      int
	Dz, Dt float64
	Linear bool // ε₀, μ₀ only (Maxwell's vacuum)

	// E at the links, z = k·Dz [V/m]; D its displacement [C/m²]; B at the
	// faces, z = (k + ½)·Dz, half a step ahead of E in time [T].
	Ex, Ey, Dx, Dy, Bx, By []float64
	Steps                  int
	hx, hy                 []float64
}

// NewEHLine makes a ring with all fields zero; courant is c·Δt/Δz (< 1).
func NewEHLine(n int, dz, courant float64) (*EHLine, error) {
	if n < 4 || !(dz > 0) || !(courant > 0 && courant < 1) {
		return nil, fmt.Errorf("%w: %d cells of %g m, Courant %g", ErrEHLine, n, dz, courant)
	}
	l := &EHLine{N: n, Dz: dz, Dt: courant * dz / em.SpeedOfLight}
	for _, p := range []*[]float64{&l.Ex, &l.Ey, &l.Dx, &l.Dy, &l.Bx, &l.By, &l.hx, &l.hy} {
		*p = make([]float64, n)
	}
	return l, nil
}

// Time is the time of E [s] (B is half a step later).
func (l *EHLine) Time() float64 { return float64(l.Steps) * l.Dt }

func (l *EHLine) next(k int) int { return (k + 1) % l.N }
func (l *EHLine) prev(k int) int { return (k + l.N - 1) % l.N }

// constitutive returns D and H for the line's medium.
func (l *EHLine) constitutive(e, b [3]float64) (d, h [3]float64) {
	if l.Linear {
		for i := 0; i < 3; i++ {
			d[i], h[i] = em.Eps0*e[i], b[i]/em.Mu0
		}
		return d, h
	}
	return em.EHConstitutive(e, b)
}

// bAtLink is B̄ at link k: the mean of the faces either side.
func (l *EHLine) bAtLink(k int) [3]float64 {
	p := l.prev(k)
	return [3]float64{(l.Bx[p] + l.Bx[k]) / 2, (l.By[p] + l.By[k]) / 2, 0}
}

// eAtFace is Ē at face k (between links k and k+1).
func (l *EHLine) eAtFace(k int) [3]float64 {
	q := l.next(k)
	return [3]float64{(l.Ex[k] + l.Ex[q]) / 2, (l.Ey[k] + l.Ey[q]) / 2, 0}
}

// SetFields sets E (at the links, time Time()) and B (at the faces, half a
// step later) from functions of z, and D from them.
func (l *EHLine) SetFields(e, b func(z float64) [3]float64) error {
	for k := 0; k < l.N; k++ {
		ek := e(float64(k) * l.Dz)
		bk := b((float64(k) + 0.5) * l.Dz)
		l.Ex[k], l.Ey[k], l.Bx[k], l.By[k] = ek[0], ek[1], bk[0], bk[1]
	}
	for k := 0; k < l.N; k++ {
		d, _ := l.constitutive([3]float64{l.Ex[k], l.Ey[k]}, l.bAtLink(k))
		l.Dx[k], l.Dy[k] = d[0], d[1]
	}
	return l.checkWeak()
}

func (l *EHLine) checkWeak() error {
	if l.Linear {
		return nil
	}
	for k := 0; k < l.N; k++ {
		if r := em.EHStrength([3]float64{l.Ex[k], l.Ey[k]}, [3]float64{l.Bx[k], l.By[k]}); r > em.EHWeakLimit {
			return fmt.Errorf("%w: cell %d at %.3g E_s: %w", ErrEHLine, k, r, em.ErrEHStrongField)
		}
	}
	return nil
}

// Step advances one time step.
func (l *EHLine) Step() error {
	// Face currents H^{n+½} from B^{n+½} and Ē^n.
	for k := 0; k < l.N; k++ {
		_, h := l.constitutive(l.eAtFace(k), [3]float64{l.Bx[k], l.By[k]})
		l.hx[k], l.hy[k] = h[0], h[1]
	}
	// Ampère at each link: Ḋ = ∇×H.
	r := l.Dt / l.Dz
	for k := 0; k < l.N; k++ {
		p := l.prev(k)
		l.Dx[k] -= r * (l.hy[k] - l.hy[p])
		l.Dy[k] += r * (l.hx[k] - l.hx[p])
	}
	// E from D, cell by cell (Newton on the 2 × 2 small-signal tensor).
	for k := 0; k < l.N; k++ {
		if err := l.solveE(k); err != nil {
			return err
		}
	}
	// Faraday round each face: Ḃ = −∇×E.
	for k := 0; k < l.N; k++ {
		q := l.next(k)
		l.Bx[k] += r * (l.Ey[q] - l.Ey[k])
		l.By[k] -= r * (l.Ex[q] - l.Ex[k])
	}
	l.Steps++
	return nil
}

// solveE finds E at link k with D(E, B̄) = (Dx, Dy), starting from the
// previous E.
func (l *EHLine) solveE(k int) error {
	if l.Linear {
		l.Ex[k], l.Ey[k] = l.Dx[k]/em.Eps0, l.Dy[k]/em.Eps0
		return nil
	}
	b := l.bAtLink(k)
	e := [3]float64{l.Ex[k], l.Ey[k]}
	// Converged at 1e-14 of the field (a few round-offs of ε₀E), or below
	// 1e-30 E_s where there is no field to speak of.
	tol := 1e-14*math.Max(math.Hypot(l.Dx[k], l.Dy[k])/em.Eps0, math.Hypot(e[0], e[1])) + 1e-30*em.SchwingerField
	for it := 0; it < 20; it++ {
		d, _ := em.EHConstitutive(e, b)
		rx, ry := l.Dx[k]-d[0], l.Dy[k]-d[1]
		s, err := em.EHLinearize(e, b)
		if err != nil {
			return fmt.Errorf("%w: cell %d at step %d: %w", ErrEHLine, k, l.Steps, err)
		}
		a, bb, c, dd := s.Eps[0][0], s.Eps[0][1], s.Eps[1][0], s.Eps[1][1]
		det := a*dd - bb*c
		ux, uy := (dd*rx-bb*ry)/det, (a*ry-c*rx)/det
		e[0] += ux
		e[1] += uy
		if math.Hypot(ux, uy) <= tol {
			l.Ex[k], l.Ey[k] = e[0], e[1]
			return nil
		}
	}
	return fmt.Errorf("%w: Newton for E at cell %d, step %d did not converge", ErrEHLine, k, l.Steps)
}

// Energy is the field energy per unit transverse area [J/m²]: Σ (E·D −
// L)·Δz at the links with B̄ there, the Hamiltonian density of the medium
// (em.EHEnergyDensity; ½ε₀E² + ½B²/μ₀ when Linear). B̄ is half a step ahead
// of E, so this is the energy to O(Δt), enough to watch for drift.
func (l *EHLine) Energy() float64 {
	w := 0.0
	for k := 0; k < l.N; k++ {
		e := [3]float64{l.Ex[k], l.Ey[k]}
		b := l.bAtLink(k)
		if l.Linear {
			w += 0.5*em.Eps0*(e[0]*e[0]+e[1]*e[1]) + 0.5*(b[0]*b[0]+b[1]*b[1])/em.Mu0
		} else {
			w += em.EHEnergyDensity(e, b)
		}
	}
	return w * l.Dz
}
