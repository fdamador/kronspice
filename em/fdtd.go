// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
	"runtime"
	"sync"
)

// Local full-wave field solver: Yee's finite-difference time-domain scheme
// (K. S. Yee, IEEE Trans. Antennas Propag. 14, 302, 1966; A. Taflove,
// S. C. Hagness, Computational Electrodynamics, 3rd ed., Artech 2005,
// Ch. 3), CLAUDE.md §3c rule 4 and §3g step 3.
//
// The region is a box of Nx×Ny×Nz grid nodes with spacing Dx, Dy, Dz. E
// lives on cell edges and B on cell faces, staggered half a cell, and B is
// half a time step ahead of E:
//
//	Ex (i+½, j, k)   Ey (i, j+½, k)   Ez (i, j, k+½)       at t = n·dt
//	Bx (i, j+½, k+½) By (i+½, j, k+½) Bz (i+½, j+½, k)     at t = (n+½)·dt
//
// Every update is one of the two curl equations of the Treatise in SI
// (doc.go):  ∂B/∂t = −curl E − (σ*/μ₀)B,  ε ∂E/∂t = curl H − σE − J.
// Because the discrete divergence of a discrete curl is exactly zero, div B
// stays at its initial value (zero) to round-off, and ε div E changes only
// by the charge the sources carry (discrete Gauss law). Both are checked by
// DivB and NodeCharge.
//
// Boundaries: the outer faces are perfect conductors (tangential E = 0).
// An optional graded matched absorbing layer lines every face: electric
// loss rate σ/ε and magnetic loss rate σ*/μ₀ are equal, so the layer has
// the impedance of vacuum, and the rate grows as (depth/L)³ so a wave
// entering it is absorbed with little reflection (Taflove & Hagness §7.2,
// lossy-layer form). Energy absorbed there is energy that left the region.
//
// Energy ledger. With W^n = Σ ½εE^n·E^n dV + Σ B^{n−½}·B^{n+½}/(2μ₀) dV,
// the scheme satisfies exactly, step by step,
//
//	W^{n+1} − W^n = −dt Σ J^{n+½}·Ē dV − dt Σ σ Ē² dV
//	                − dt Σ (σ*/μ₀²)/4 · B^{n+½}·(B^{n−½} + 2B^{n+½} + B^{n+3/2}) dV
//
// with Ē = (E^n + E^{n+1})/2: source work in, heat in materials, and energy
// absorbed by the layer. It follows from the two update equations and the
// discrete identity Σ E·curl_h H = Σ H·curl_e E (the curls are adjoint with
// conducting walls), the discrete form of Poynting's theorem. The residual
// of this ledger is round-off; a larger one is a bug.

// FDTDConfig sets up a field region.
type FDTDConfig struct {
	N      [3]int     // grid nodes per axis (≥ 3)
	D      [3]float64 // node spacing [m]
	Origin [3]float64 // position of node (0, 0, 0) [m]
	// Courant number S: dt = S / (c·√(1/Dx² + 1/Dy² + 1/Dz²)); default 0.99,
	// stable for S < 1 in vacuum.
	Courant float64
	// Absorber is the thickness in cells of the absorbing layer on every
	// face; 0 gives a closed box with conducting walls.
	Absorber int
	// Reflection is the absorbing layer's design reflection at normal
	// incidence; default 1e-6.
	Reflection float64
	// Material gives relative permittivity and conductivity [S/m] at a
	// point; nil is vacuum. Permeability is μ₀ everywhere.
	Material func(x, y, z float64) (epsR, sigma float64)
	// Ledger keeps the energy ledger (needs one more copy of B).
	Ledger bool
}

// FieldLedger is the energy account of a field region [J].
type FieldLedger struct {
	W0, W        float64 // field energy at the start and now
	SourceWork   float64 // work done on the field by the impressed currents
	MaterialLoss float64 // heat in conducting materials
	Absorbed     float64 // energy absorbed by the boundary layer (left the region)
}

// Residual is W − W0 − SourceWork + MaterialLoss + Absorbed (zero up to round-off).
func (l FieldLedger) Residual() float64 {
	return l.W - l.W0 - l.SourceWork + l.MaterialLoss + l.Absorbed
}

// EdgeSource is an impressed current I(t) [A] along one E edge: Comp 0, 1,
// 2 for an x, y or z edge starting at node (I, J, K).
type EdgeSource struct {
	Comp    int
	I, J, K int
	Current func(t float64) float64
}

type eCoef struct {
	ca, cb     float64
	eps        float64 // ε [F/m]
	rMat, rAbs float64 // loss rates σ/ε of material and absorber [1/s]
}

type bCoef struct{ da, db, rAbs float64 }

// FDTD is a field region being stepped in time.
type FDTD struct {
	N      [3]int
	D      [3]float64
	Origin [3]float64
	Dt     float64
	Steps  int // completed steps; E is at Time(), B half a step later

	E, B   [3][]float64
	bPrev  [3][]float64 // B half a step behind B (ledger)
	eIdx   [3][]uint8
	bIdx   [3][]uint8
	eTab   []eCoef
	bTab   []bCoef
	src    []EdgeSource
	ledger bool
	led    FieldLedger
	absorb int
	start  bool
	t0     float64
}

var (
	// ErrFDTDConfig reports an invalid field-region set-up.
	ErrFDTDConfig = errors.New("em: invalid FDTD configuration")
	// ErrFDTDSource reports a source on an edge the solver cannot drive.
	ErrFDTDSource = errors.New("em: FDTD source not on a free edge")
)

// NewFDTD builds a field region from cfg with all fields zero.
func NewFDTD(cfg FDTDConfig) (*FDTD, error) {
	for a := 0; a < 3; a++ {
		if cfg.N[a] < 3 || !(cfg.D[a] > 0) {
			return nil, fmt.Errorf("%w: axis %d has %d nodes, spacing %g", ErrFDTDConfig, a, cfg.N[a], cfg.D[a])
		}
		if 2*cfg.Absorber >= cfg.N[a]-2 {
			return nil, fmt.Errorf("%w: absorber of %d cells leaves no interior on axis %d", ErrFDTDConfig, cfg.Absorber, a)
		}
	}
	n := cfg.N[0] * cfg.N[1] * cfg.N[2]
	if n > 200_000_000 {
		return nil, fmt.Errorf("%w: %d cells is too many", ErrFDTDConfig, n)
	}
	s := cfg.Courant
	if s == 0 {
		s = 0.99
	}
	if !(s > 0 && s < 1) {
		return nil, fmt.Errorf("%w: Courant number %g must be in (0, 1)", ErrFDTDConfig, s)
	}
	inv := 0.0
	for a := 0; a < 3; a++ {
		inv += 1 / (cfg.D[a] * cfg.D[a])
	}
	f := &FDTD{N: cfg.N, D: cfg.D, Origin: cfg.Origin, Dt: s / (SpeedOfLight * math.Sqrt(inv)), ledger: cfg.Ledger, absorb: cfg.Absorber}
	for c := 0; c < 3; c++ {
		f.E[c] = make([]float64, n)
		f.B[c] = make([]float64, n)
		f.eIdx[c] = make([]uint8, n)
		f.bIdx[c] = make([]uint8, n)
		if cfg.Ledger {
			f.bPrev[c] = make([]float64, n)
		}
	}
	// Absorber loss-rate profile: rate(d) = rmax·(d/L)³, d the depth into
	// the layer in cells; normal-incidence reflection exp(−2∫rate/c dx)
	// = exp(−2·rmax·L·dx/(4c)) sets rmax.
	L := float64(cfg.Absorber)
	rmax := 0.0
	if cfg.Absorber > 0 {
		R := cfg.Reflection
		if R == 0 {
			R = 1e-6
		}
		if !(R > 0 && R < 1) {
			return nil, fmt.Errorf("%w: reflection %g", ErrFDTDConfig, R)
		}
		dmin := math.Min(cfg.D[0], math.Min(cfg.D[1], cfg.D[2]))
		rmax = -4 * SpeedOfLight * math.Log(R) / (2 * L * dmin)
	}
	depth := func(p [3]float64) float64 { // p in index units
		d := 0.0
		for a := 0; a < 3; a++ {
			hi := float64(f.N[a]-1) - L
			d = math.Max(d, math.Max(L-p[a], p[a]-hi))
		}
		return d
	}
	eKey := map[[3]float64]uint8{}
	bKey := map[float64]uint8{}
	addE := func(epsR, sigma, rAbs float64) (uint8, error) {
		k := [3]float64{epsR, sigma, rAbs}
		if id, ok := eKey[k]; ok {
			return id, nil
		}
		if len(f.eTab) == 256 {
			return 0, fmt.Errorf("%w: more than 256 distinct material/absorber cells", ErrFDTDConfig)
		}
		if !(epsR >= 1) || sigma < 0 {
			return 0, fmt.Errorf("%w: material εr %g, σ %g", ErrFDTDConfig, epsR, sigma)
		}
		eps := Eps0 * epsR
		r := sigma/eps + rAbs
		x := r * f.Dt / 2
		f.eTab = append(f.eTab, eCoef{ca: (1 - x) / (1 + x), cb: f.Dt / eps / (1 + x), eps: eps, rMat: sigma / eps, rAbs: rAbs})
		eKey[k] = uint8(len(f.eTab) - 1)
		return eKey[k], nil
	}
	addB := func(rAbs float64) uint8 {
		if id, ok := bKey[rAbs]; ok {
			return id
		}
		x := rAbs * f.Dt / 2
		f.bTab = append(f.bTab, bCoef{da: (1 - x) / (1 + x), db: f.Dt / (1 + x), rAbs: rAbs})
		bKey[rAbs] = uint8(len(f.bTab) - 1)
		return bKey[rAbs]
	}
	// Quantize absorber depth to sixteenths of a cell so the tables stay small.
	rate := func(p [3]float64) float64 {
		if cfg.Absorber == 0 {
			return 0
		}
		d := math.Round(depth(p)*16) / 16
		if d <= 0 {
			return 0
		}
		return rmax * math.Pow(d/L, 3)
	}
	addB(0)
	if _, err := addE(1, 0, 0); err != nil {
		return nil, err
	}
	for c := 0; c < 3; c++ {
		for k := 0; k < f.N[2]; k++ {
			for j := 0; j < f.N[1]; j++ {
				for i := 0; i < f.N[0]; i++ {
					n := f.index(i, j, k)
					pe := [3]float64{float64(i), float64(j), float64(k)}
					pe[c] += 0.5
					epsR, sigma := 1.0, 0.0
					if cfg.Material != nil {
						x := f.Origin[0] + pe[0]*f.D[0]
						y := f.Origin[1] + pe[1]*f.D[1]
						z := f.Origin[2] + pe[2]*f.D[2]
						epsR, sigma = cfg.Material(x, y, z)
					}
					id, err := addE(epsR, sigma, rate(pe))
					if err != nil {
						return nil, err
					}
					f.eIdx[c][n] = id
					pb := [3]float64{float64(i) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
					pb[c] -= 0.5
					if len(f.bTab) == 256 {
						return nil, fmt.Errorf("%w: absorber profile needs more than 256 levels", ErrFDTDConfig)
					}
					f.bIdx[c][n] = addB(rate(pb))
				}
			}
		}
	}
	return f, nil
}

func (f *FDTD) index(i, j, k int) int { return i + f.N[0]*(j+f.N[1]*k) }

// Time is the time of E [s]; B is at Time() + Dt/2.
func (f *FDTD) Time() float64 { return f.t0 + float64(f.Steps)*f.Dt }

// SetStartTime sets the time of the initial E field (default 0).
func (f *FDTD) SetStartTime(t float64) { f.t0 = t }

// Ledger returns the energy account (zero unless the config asked for it).
func (f *FDTD) Ledger() FieldLedger { return f.led }

// cellVolume is dV [m³].
func (f *FDTD) cellVolume() float64 { return f.D[0] * f.D[1] * f.D[2] }

// ePos and bPos are the positions of field components [m].
func (f *FDTD) ePos(c, i, j, k int) (x, y, z float64) {
	p := [3]float64{float64(i), float64(j), float64(k)}
	p[c] += 0.5
	return f.Origin[0] + p[0]*f.D[0], f.Origin[1] + p[1]*f.D[1], f.Origin[2] + p[2]*f.D[2]
}

func (f *FDTD) bPos(c, i, j, k int) (x, y, z float64) {
	p := [3]float64{float64(i) + 0.5, float64(j) + 0.5, float64(k) + 0.5}
	p[c] -= 0.5
	return f.Origin[0] + p[0]*f.D[0], f.Origin[1] + p[1]*f.D[1], f.Origin[2] + p[2]*f.D[2]
}

// eFree reports whether E component c at node (i, j, k) is updated (not
// tangential to a conducting wall and inside the array).
func (f *FDTD) eFree(c, i, j, k int) bool {
	ijk := [3]int{i, j, k}
	for a := 0; a < 3; a++ {
		if a == c {
			if ijk[a] < 0 || ijk[a] > f.N[a]-2 {
				return false
			}
		} else if ijk[a] < 1 || ijk[a] > f.N[a]-2 {
			return false
		}
	}
	return true
}

// SetE samples E [V/m] at the edge positions (time Time()); edges on the
// conducting walls stay zero.
func (f *FDTD) SetE(e func(x, y, z float64) [3]float64) {
	f.start = false
	f.each(func(c, i, j, k, n int) {
		if !f.eFree(c, i, j, k) {
			f.E[c][n] = 0
			return
		}
		x, y, z := f.ePos(c, i, j, k)
		f.E[c][n] = e(x, y, z)[c]
	})
}

// SetB samples B [T] at the face positions (time Time() + Dt/2).
func (f *FDTD) SetB(b func(x, y, z float64) [3]float64) {
	f.start = false
	f.each(func(c, i, j, k, n int) {
		x, y, z := f.bPos(c, i, j, k)
		f.B[c][n] = b(x, y, z)[c]
	})
}

// SetBFromPotential sets B = curl A with A [T·m] sampled on the edges and
// the discrete curl, so the discrete div B is exactly zero.
func (f *FDTD) SetBFromPotential(A func(x, y, z float64) [3]float64) {
	var a [3][]float64
	for c := 0; c < 3; c++ {
		a[c] = make([]float64, len(f.E[c]))
	}
	f.each(func(c, i, j, k, n int) {
		x, y, z := f.ePos(c, i, j, k)
		a[c][n] = A(x, y, z)[c]
	})
	sy, sz := f.N[0], f.N[0]*f.N[1]
	f.each(func(c, i, j, k, n int) {
		ok := [3]bool{i < f.N[0]-1, j < f.N[1]-1, k < f.N[2]-1}
		ok[c] = true
		if !(ok[0] && ok[1] && ok[2]) {
			f.B[c][n] = 0
			return
		}
		switch c {
		case 0:
			f.B[0][n] = (a[2][n+sy]-a[2][n])/f.D[1] - (a[1][n+sz]-a[1][n])/f.D[2]
		case 1:
			f.B[1][n] = (a[0][n+sz]-a[0][n])/f.D[2] - (a[2][n+1]-a[2][n])/f.D[0]
		default:
			f.B[2][n] = (a[1][n+1]-a[1][n])/f.D[0] - (a[0][n+sy]-a[0][n])/f.D[1]
		}
	})
	f.start = false
}

// SetEFromPotential sets E = curl_h C with C [V] sampled on the faces and
// the discrete dual curl, so the discrete div E is exactly zero (vacuum).
func (f *FDTD) SetEFromPotential(C func(x, y, z float64) [3]float64) {
	var h [3][]float64
	for c := 0; c < 3; c++ {
		h[c] = make([]float64, len(f.B[c]))
	}
	f.each(func(c, i, j, k, n int) {
		x, y, z := f.bPos(c, i, j, k)
		h[c][n] = C(x, y, z)[c]
	})
	f.each(func(c, i, j, k, n int) {
		if f.eFree(c, i, j, k) {
			f.E[c][n] = f.curlH(h, c, n, 1)
		} else {
			f.E[c][n] = 0
		}
	})
	f.start = false
}

// AddSource adds an impressed edge current.
func (f *FDTD) AddSource(s EdgeSource) error {
	if s.Comp < 0 || s.Comp > 2 || s.Current == nil || !f.eFree(s.Comp, s.I, s.J, s.K) {
		return fmt.Errorf("%w: component %d at (%d, %d, %d)", ErrFDTDSource, s.Comp, s.I, s.J, s.K)
	}
	f.src = append(f.src, s)
	return nil
}

func (f *FDTD) each(fn func(c, i, j, k, n int)) {
	for c := 0; c < 3; c++ {
		for k := 0; k < f.N[2]; k++ {
			for j := 0; j < f.N[1]; j++ {
				for i := 0; i < f.N[0]; i++ {
					fn(c, i, j, k, f.index(i, j, k))
				}
			}
		}
	}
}

// curlH is component c of curl_h(h) at E edge n, times scale.
func (f *FDTD) curlH(h [3][]float64, c, n int, scale float64) float64 {
	sy, sz := f.N[0], f.N[0]*f.N[1]
	dx, dy, dz := f.D[0], f.D[1], f.D[2]
	switch c {
	case 0:
		return scale * ((h[2][n]-h[2][n-sy])/dy - (h[1][n]-h[1][n-sz])/dz)
	case 1:
		return scale * ((h[0][n]-h[0][n-sz])/dz - (h[2][n]-h[2][n-1])/dx)
	default:
		return scale * ((h[1][n]-h[1][n-1])/dx - (h[0][n]-h[0][n-sy])/dy)
	}
}

// parallel runs fn over k-slabs; part collects per-slab sums in order so
// results do not depend on scheduling.
func (f *FDTD) parallel(nk int, fn func(k0, k1 int, part *[4]float64)) [4]float64 {
	w := runtime.NumCPU()
	if w > nk {
		w = nk
	}
	if w < 1 {
		w = 1
	}
	parts := make([][4]float64, w)
	var wg sync.WaitGroup
	for p := 0; p < w; p++ {
		k0, k1 := p*nk/w, (p+1)*nk/w
		wg.Add(1)
		go func(p, k0, k1 int) {
			defer wg.Done()
			fn(k0, k1, &parts[p])
		}(p, k0, k1)
	}
	wg.Wait()
	var s [4]float64
	for _, q := range parts {
		for i := range s {
			s[i] += q[i]
		}
	}
	return s
}

// faraday applies B ← da·B − db·scale·curl_e(e) to every face; with prev
// (ledger) it also returns the magnetic absorber loss and the magnetic
// energy B_old·B_new/(2μ₀) (both ×dV), and moves B_old into prev.
func (f *FDTD) faraday(e [3][]float64, scale float64, b [3][]float64, prev *[3][]float64) (loss, wm float64) {
	nx, ny, nz := f.N[0], f.N[1], f.N[2]
	sy, sz := nx, nx*ny
	dx, dy, dz := f.D[0], f.D[1], f.D[2]
	tab := f.bTab
	sums := f.parallel(nz, func(k0, k1 int, part *[4]float64) {
		for c := 0; c < 3; c++ {
			bc, idx := b[c], f.bIdx[c]
			var pc []float64
			if prev != nil {
				pc = (*prev)[c]
			}
			imax, jmax, kmax := nx-1, ny-1, nz-1
			switch c {
			case 0:
				imax = nx
			case 1:
				jmax = ny
			case 2:
				kmax = nz
			}
			for k := k0; k < k1 && k < kmax; k++ {
				for j := 0; j < jmax; j++ {
					n := f.index(0, j, k)
					for i := 0; i < imax; i, n = i+1, n+1 {
						var curl float64
						switch c {
						case 0:
							curl = (e[2][n+sy]-e[2][n])/dy - (e[1][n+sz]-e[1][n])/dz
						case 1:
							curl = (e[0][n+sz]-e[0][n])/dz - (e[2][n+1]-e[2][n])/dx
						default:
							curl = (e[1][n+1]-e[1][n])/dx - (e[0][n+sy]-e[0][n])/dy
						}
						t := &tab[idx[n]]
						old := bc[n]
						nw := t.da*old - t.db*scale*curl
						bc[n] = nw
						if pc != nil {
							if t.rAbs != 0 {
								part[0] += t.rAbs * old * (pc[n] + 2*old + nw)
							}
							part[1] += old * nw
							pc[n] = old
						}
					}
				}
			}
		}
	})
	return sums[0], sums[1]
}

// ampere applies E ← ca·E + cb·curl_h(B/μ₀) on every free edge; with the
// ledger it returns Σ σĒ² (material, absorber) and Σ ½εE² (all ×1).
func (f *FDTD) ampere() (lossMat, lossAbs, we float64) {
	nx, ny, nz := f.N[0], f.N[1], f.N[2]
	sy, sz := nx, nx*ny
	dx, dy, dz := f.D[0], f.D[1], f.D[2]
	b := f.B
	tab := f.eTab
	led := f.ledger
	inv := 1 / Mu0
	sums := f.parallel(nz, func(k0, k1 int, part *[4]float64) {
		for c := 0; c < 3; c++ {
			ec, idx := f.E[c], f.eIdx[c]
			i0, i1, j0, j1, kk0, kk1 := 1, nx-1, 1, ny-1, 1, nz-1
			switch c {
			case 0:
				i0, i1 = 0, nx-1
			case 1:
				j0, j1 = 0, ny-1
			default:
				kk0, kk1 = 0, nz-1
			}
			for k := max(k0, kk0); k < min(k1, kk1); k++ {
				for j := j0; j < j1; j++ {
					n := f.index(i0, j, k)
					for i := i0; i < i1; i, n = i+1, n+1 {
						var curl float64
						switch c {
						case 0:
							curl = (b[2][n]-b[2][n-sy])/dy - (b[1][n]-b[1][n-sz])/dz
						case 1:
							curl = (b[0][n]-b[0][n-sz])/dz - (b[2][n]-b[2][n-1])/dx
						default:
							curl = (b[1][n]-b[1][n-1])/dx - (b[0][n]-b[0][n-sy])/dy
						}
						t := &tab[idx[n]]
						old := ec[n]
						nw := t.ca*old + t.cb*inv*curl
						ec[n] = nw
						if led {
							if t.rMat != 0 || t.rAbs != 0 {
								m := 0.5 * (old + nw)
								part[0] += t.rMat * t.eps * m * m
								part[1] += t.rAbs * t.eps * m * m
							}
							part[2] += 0.5 * t.eps * nw * nw
						}
					}
				}
			}
		}
	})
	return sums[0], sums[1], sums[2]
}

// begin prepares the ledger: B^{−½} by running the first B update
// backwards with the initial E, and W^0.
func (f *FDTD) begin() {
	f.start = true
	if !f.ledger {
		return
	}
	// B^{½} = da·B^{−½} − db·curl E^0  ⇒  B^{−½} = (B^{½} + db·curl E^0)/da.
	var tmp [3][]float64
	for c := 0; c < 3; c++ {
		tmp[c] = make([]float64, len(f.B[c]))
	}
	f.faraday(f.E, -1, tmp, nil) // tmp = +db·curl E^0 (with da = 1 on tmp's zeros: tmp = db·curl)
	for c := 0; c < 3; c++ {
		for n := range tmp[c] {
			t := f.bTab[f.bIdx[c][n]]
			f.bPrev[c][n] = (f.B[c][n] + tmp[c][n]) / t.da
		}
	}
	dV := f.cellVolume()
	we, wm := 0.0, 0.0
	f.each(func(c, i, j, k, n int) {
		t := f.eTab[f.eIdx[c][n]]
		we += 0.5 * t.eps * f.E[c][n] * f.E[c][n]
		wm += f.bPrev[c][n] * f.B[c][n] / (2 * Mu0)
	})
	f.led = FieldLedger{W0: (we + wm) * dV, W: (we + wm) * dV}
}

// Step advances the region by one time step: E from Time() to Time()+Dt,
// then B to half a step later.
func (f *FDTD) Step() {
	if !f.start {
		f.begin()
	}
	dt, dV := f.Dt, f.cellVolume()
	tHalf := f.Time() + dt/2
	type was struct {
		c, n int
		old  float64
		j    float64
	}
	var srcs []was
	for _, s := range f.src {
		n := f.index(s.I, s.J, s.K)
		area := 1.0
		for a := 0; a < 3; a++ {
			if a != s.Comp {
				area *= f.D[a]
			}
		}
		srcs = append(srcs, was{s.Comp, n, f.E[s.Comp][n], s.Current(tHalf) / area})
	}
	lm, la, we := f.ampere()
	for _, s := range srcs {
		t := f.eTab[f.eIdx[s.c][s.n]]
		raw := f.E[s.c][s.n]
		nw := raw - t.cb*s.j
		f.E[s.c][s.n] = nw
		if f.ledger {
			mr, mn := 0.5*(s.old+raw), 0.5*(s.old+nw)
			lm += t.rMat * t.eps * (mn*mn - mr*mr)
			la += t.rAbs * t.eps * (mn*mn - mr*mr)
			we += 0.5 * t.eps * (nw*nw - raw*raw)
			f.led.SourceWork -= dt * s.j * mn * dV
		}
	}
	var prev *[3][]float64
	if f.ledger {
		prev = &f.bPrev
	}
	lb, wm := f.faraday(f.E, 1, f.B, prev)
	f.Steps++
	if f.ledger {
		f.led.MaterialLoss += dt * lm * dV
		f.led.Absorbed += dt*la*dV + dt*lb/(4*Mu0)*dV
		f.led.W = (we + wm/(2*Mu0)) * dV
	}
}

// Run steps until Time() ≥ t.
func (f *FDTD) Run(t float64) {
	for f.Time() < t-1e-9*f.Dt {
		f.Step()
	}
}

// DivB returns the largest |div B| over the cells inside the absorbing
// layer's inner edge, and the scale max|B|/min(D) it is to be compared
// with. (Inside the graded layer the loss rate varies from cell to cell,
// so div(rate·B) is not rate·div B and div B is not kept there; the layer
// only absorbs.)
func (f *FDTD) DivB() (worst, scale float64) {
	nx, ny, nz := f.N[0], f.N[1], f.N[2]
	sy, sz := nx, nx*ny
	bmax := 0.0
	for c := 0; c < 3; c++ {
		for _, v := range f.B[c] {
			bmax = math.Max(bmax, math.Abs(v))
		}
	}
	a := f.absorb
	for k := a; k < nz-1-a; k++ {
		for j := a; j < ny-1-a; j++ {
			for i := a; i < nx-1-a; i++ {
				n := f.index(i, j, k)
				d := (f.B[0][n+1]-f.B[0][n])/f.D[0] + (f.B[1][n+sy]-f.B[1][n])/f.D[1] + (f.B[2][n+sz]-f.B[2][n])/f.D[2]
				worst = math.Max(worst, math.Abs(d))
			}
		}
	}
	return worst, bmax / math.Min(f.D[0], math.Min(f.D[1], f.D[2]))
}

// NodeCharge returns the charge [C] the field implies at interior node
// (i, j, k): ∮εE·dA over the dual cell (discrete Gauss law).
func (f *FDTD) NodeCharge(i, j, k int) float64 {
	n := f.index(i, j, k)
	sy, sz := f.N[0], f.N[0]*f.N[1]
	eps := func(c, m int) float64 { return f.eTab[f.eIdx[c][m]].eps }
	dx, dy, dz := f.D[0], f.D[1], f.D[2]
	return (eps(0, n)*f.E[0][n]-eps(0, n-1)*f.E[0][n-1])*dy*dz +
		(eps(1, n)*f.E[1][n]-eps(1, n-sy)*f.E[1][n-sy])*dx*dz +
		(eps(2, n)*f.E[2][n]-eps(2, n-sz)*f.E[2][n-sz])*dx*dy
}

// FieldAt interpolates E and B (each from its own staggered grid,
// trilinearly) at a point [m]. B is half a step ahead of E.
func (f *FDTD) FieldAt(x, y, z float64) (E, B [3]float64) {
	interp := func(arr []float64, off [3]float64) float64 {
		p := [3]float64{(x - f.Origin[0]) / f.D[0], (y - f.Origin[1]) / f.D[1], (z - f.Origin[2]) / f.D[2]}
		var i0 [3]int
		var w [3]float64
		for a := 0; a < 3; a++ {
			q := p[a] - off[a]
			fl := math.Floor(q)
			i0[a] = int(fl)
			w[a] = q - fl
			if i0[a] < 0 {
				i0[a], w[a] = 0, 0
			}
			if i0[a] > f.N[a]-2 {
				i0[a], w[a] = f.N[a]-2, 1
			}
		}
		s := 0.0
		for dk := 0; dk < 2; dk++ {
			for dj := 0; dj < 2; dj++ {
				for di := 0; di < 2; di++ {
					wt := (1 - w[0] + float64(di)*(2*w[0]-1)) * (1 - w[1] + float64(dj)*(2*w[1]-1)) * (1 - w[2] + float64(dk)*(2*w[2]-1))
					s += wt * arr[f.index(i0[0]+di, i0[1]+dj, i0[2]+dk)]
				}
			}
		}
		return s
	}
	for c := 0; c < 3; c++ {
		oe := [3]float64{}
		oe[c] = 0.5
		ob := [3]float64{0.5, 0.5, 0.5}
		ob[c] = 0
		E[c] = interp(f.E[c], oe)
		B[c] = interp(f.B[c], ob)
	}
	return
}

// RegionEnergy is the field energy [J] inside a region, with the same
// discrete energy as the ledger (B^{n−½}·B^{n+½} when the ledger is kept,
// B² otherwise).
func (f *FDTD) RegionEnergy(in func(x, y, z float64) bool) float64 {
	if !f.start {
		f.begin()
	}
	w := 0.0
	f.each(func(c, i, j, k, n int) {
		if x, y, z := f.ePos(c, i, j, k); in(x, y, z) {
			t := f.eTab[f.eIdx[c][n]]
			w += 0.5 * t.eps * f.E[c][n] * f.E[c][n]
		}
		if x, y, z := f.bPos(c, i, j, k); in(x, y, z) {
			if f.ledger {
				w += f.bPrev[c][n] * f.B[c][n] / (2 * Mu0)
			} else {
				w += f.B[c][n] * f.B[c][n] / (2 * Mu0)
			}
		}
	})
	return w * f.cellVolume()
}
