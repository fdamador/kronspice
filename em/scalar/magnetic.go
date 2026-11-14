// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package scalar

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
)

// The magnetic sector of the region (hypothetical, H2 + H47; equations in
// the package comment). Layout on the Yee grid of em.FDTD:
//
//   - S_m and ρ_m on the cell centres, indexed by the cell's lower corner
//     (i, j, k), 0 ≤ i ≤ N₀−2 and so on, where div B lives;
//   - S_m at the whole steps with E, ρ_m at the half steps with B;
//   - grad S_m on the B faces, with S_m = 0 outside the box (all faces,
//     the wall-normal ones too), so Σ_faces B·grad S_m = −Σ_cells S_m div B
//     exactly and the energy theorem closes step by step;
//   - magnetic currents on B faces, entering Faraday as em.FDTD's
//     FaceSource does (M = μ₀J_m, applied with the time of E^{n+1}, work
//     −dt B^{n+½}·(J_m^n + J_m^{n+1})/2 dV with J_m^0 = 0), moving
//     magnetic charge between the two cells they separate, so that it is
//     conserved.

// MagMode selects the magnetic sector.
type MagMode int

const (
	// MagNone: no magnetic charge (standard). Magnetic sources are refused.
	MagNone MagMode = iota
	// MagDual: magnetic charge and current with S_m held at 0, the
	// dual-symmetric Maxwell equations of em/duality.go.
	MagDual
	// MagAB: magnetic charge with S_m physical and evolved.
	MagAB
)

func (m MagMode) String() string { return [...]string{"none", "dual", "AB"}[m] }

// ErrMagneticDisabled is returned when a magnetic source or S_m is set on a
// region without magnetic charge enabled (hypothetical mode, as
// em.ErrHypotheticalDisabled).
var ErrMagneticDisabled = errors.New("scalar: magnetic charge needs Config.Magnetic (hypothetical mode)")

// MagFaceSource is an impressed magnetic current through one B face:
// component Comp of face (I, J, K), total current [A·m/s]; J_m = Current /
// face area. It moves magnetic charge from the cell below the face (along
// Comp) to the cell above, so it is conserved and does not source S_m.
type MagFaceSource struct {
	Comp    int
	I, J, K int
	Current func(t float64) float64
}

// MagCellSource is magnetic charge appearing (or vanishing) in one cell at
// Rate [A·m/s] with no current to bring it: the source of S_m,
// □S_m = (μ₀/c)(∂ρ_m/∂t + div J_m).
type MagCellSource struct {
	I, J, K int
	Rate    func(t float64) float64
}

type magState struct {
	fsrc   []MagFaceSource
	fPrev  []float64 // each face source's density J_m in the previous step
	csrc   []MagCellSource
	defect float64 // relative magnetic Gauss defect of the last step
}

// cell reports whether (i, j, k) is a cell of the box.
func (r *Region) cell(i, j, k int) bool {
	return i >= 0 && j >= 0 && k >= 0 && i <= r.N[0]-2 && j <= r.N[1]-2 && k <= r.N[2]-2
}

func (r *Region) smAt(i, j, k int) float64 {
	if !r.cell(i, j, k) {
		return 0
	}
	return r.Sm[r.index(i, j, k)]
}

// gradSm is component c of grad S_m on B face (i, j, k): between the cell
// (i, j, k) and the cell one step below it along c.
func (r *Region) gradSm(c, i, j, k int) float64 {
	if r.Magnetic != MagAB {
		return 0
	}
	ijk := [3]int{i, j, k}
	ijk[c]--
	return (r.smAt(i, j, k) - r.smAt(ijk[0], ijk[1], ijk[2])) / r.D[c]
}

// divB is the discrete divergence of B in cell n.
func (r *Region) divB(n int) float64 {
	sy, sz := r.N[0], r.N[0]*r.N[1]
	return (r.B[0][n+1]-r.B[0][n])/r.D[0] + (r.B[1][n+sy]-r.B[1][n])/r.D[1] + (r.B[2][n+sz]-r.B[2][n])/r.D[2]
}

func (r *Region) eachCell(fn func(i, j, k, n int)) {
	for k := 0; k <= r.N[2]-2; k++ {
		for j := 0; j <= r.N[1]-2; j++ {
			for i := 0; i <= r.N[0]-2; i++ {
				fn(i, j, k, r.index(i, j, k))
			}
		}
	}
}

// energySm is Σ (S_m^n)²/(2μ₀) dV [J].
func (r *Region) energySm() float64 {
	if r.Magnetic != MagAB {
		return 0
	}
	var w float64
	r.eachCell(func(i, j, k, n int) { w += r.Sm[n] * r.Sm[n] / (2 * em.Mu0) })
	return w * r.dV()
}

// stepSm advances S_m by one step from B^{n+½} and ρ_m^{n+½}, books the work
// cρ_mS̄_m, records the magnetic Gauss defect, and returns the new S_m
// energy [J].
func (r *Region) stepSm() float64 {
	if r.Magnetic == MagNone {
		return 0
	}
	c, dt, dV := em.SpeedOfLight, r.Dt, r.dV()
	var work, worst, scale, w float64
	r.eachCell(func(i, j, k, n int) {
		g := em.Mu0*r.RhoM[n] - r.divB(n)
		scale = math.Max(scale, math.Max(math.Abs(em.Mu0*r.RhoM[n]), math.Abs(r.divB(n))))
		if r.Magnetic == MagAB {
			old := r.Sm[n]
			nw := old + c*dt*g
			work += r.RhoM[n] * 0.5 * (old + nw)
			r.Sm[n] = nw
			w += nw * nw / (2 * em.Mu0)
			g -= (nw - old) / (c * dt) // zero by construction: the AB Gauss law
		}
		worst = math.Max(worst, math.Abs(g))
	})
	r.mag.defect = 0
	if scale > 0 {
		r.mag.defect = worst / scale
	}
	r.led.MagScalarWork += dt * c * work * dV
	return w * dV
}

// applyMagSources adds −dt μ₀J_m^{n+1} on the source faces (B^{n+½} is in
// bPrev after the Faraday sweep), books their work, moves and creates
// magnetic charge (ρ_m to n + 3/2), and returns the change it makes to
// Σ B^{n+½}·B^{n+3/2}/(2μ₀) (per unit volume).
func (r *Region) applyMagSources() float64 {
	if r.Magnetic == MagNone {
		return 0
	}
	dt, dV := r.Dt, r.dV()
	tE := r.Time() + dt
	var dw float64
	for q, s := range r.mag.fsrc {
		n := r.index(s.I, s.J, s.K)
		area := 1.0
		for a := 0; a < 3; a++ {
			if a != s.Comp {
				area *= r.D[a]
			}
		}
		jm := s.Current(tE) / area
		bh := r.bPrev[s.Comp][n]
		r.B[s.Comp][n] -= dt * em.Mu0 * jm
		dw -= bh * dt * jm / 2
		r.led.MagSourceWork -= dt * bh * 0.5 * (r.mag.fPrev[q] + jm) * dV
		r.mag.fPrev[q] = jm
		// the magnetic charge crosses the face, from the cell below to the cell above
		ijk := [3]int{s.I, s.J, s.K}
		ijk[s.Comp]--
		dq := dt * jm / r.D[s.Comp]
		r.RhoM[r.index(ijk[0], ijk[1], ijk[2])] -= dq
		r.RhoM[n] += dq
	}
	for _, s := range r.mag.csrc {
		r.RhoM[r.index(s.I, s.J, s.K)] += dt * s.Rate(tE) / dV
	}
	return dw
}

// AddMagSource adds a magnetic current on an interior B face.
func (r *Region) AddMagSource(s MagFaceSource) error {
	if r.Magnetic == MagNone {
		return ErrMagneticDisabled
	}
	if s.Comp < 0 || s.Comp > 2 || s.Current == nil {
		return fmt.Errorf("%w: magnetic source component %d", ErrConfig, s.Comp)
	}
	below := [3]int{s.I, s.J, s.K}
	below[s.Comp]--
	if !r.cell(s.I, s.J, s.K) || !r.cell(below[0], below[1], below[2]) {
		return fmt.Errorf("%w: magnetic source face (%d, %d, %d) needs a cell on each side", ErrConfig, s.I, s.J, s.K)
	}
	r.mag.fsrc = append(r.mag.fsrc, s)
	r.mag.fPrev = append(r.mag.fPrev, 0)
	return nil
}

// AddMagCellSource adds created magnetic charge in a cell.
func (r *Region) AddMagCellSource(s MagCellSource) error {
	if r.Magnetic == MagNone {
		return ErrMagneticDisabled
	}
	if !r.cell(s.I, s.J, s.K) || s.Rate == nil {
		return fmt.Errorf("%w: magnetic cell source not in a cell", ErrConfig)
	}
	r.mag.csrc = append(r.mag.csrc, s)
	return nil
}

// SetSm samples S_m [T] at the cell centres (time Time()). MagAB only.
func (r *Region) SetSm(s func(x, y, z float64) float64) error {
	if r.Magnetic != MagAB {
		return ErrMagneticDisabled
	}
	r.eachCell(func(i, j, k, n int) {
		x, y, z := r.pos([3]float64{0.5, 0.5, 0.5}, i, j, k)
		r.Sm[n] = s(x, y, z)
	})
	return nil
}

// SmAt returns S_m in cell (i, j, k) at Time().
func (r *Region) SmAt(i, j, k int) float64 { return r.Sm[r.index(i, j, k)] }

// MagGaussDefect is the last step's max over cells of
// |μ₀ρ_m − div B − (1/c)∂S_m/∂t|, relative to the largest |μ₀ρ_m| or
// |div B|. In MagAB it is round-off by construction; in MagDual created
// magnetic charge makes it O(1), the dual of GaussDefect in standard mode.
func (r *Region) MagGaussDefect() float64 { return r.mag.defect }
