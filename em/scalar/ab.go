// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package scalar holds the tier-2 extended electrodynamics of CLAUDE.md §3k:
// theories in which the scalar part of Maxwell's quaternion equations, the
// convergence S.∇𝔄 that Maxwell set to zero by choice (Treatise Vol. II,
// Art. 619), is kept as a physical field. Hypothetical mode only (H8, H47):
// every result is reported beside the standard prediction S = 0.
//
// This file is Aharonov–Bohm (Ohmura–Stueckelberg) electrodynamics in flat
// space, the electromagnetic sector of Minotti & Modanese (arXiv:2408.05230,
// eq. 1 with the β terms off): the action carries −(ε₀c/2)(∂_μA^μ)² beside
// −(ε₀c/4)F_μνF^μν, so S = ∂_μA^μ = ∇·A + (1/c²)∂φ/∂t is physical and
//
//	∂_μF^μν = μ₀j^ν − ∂^νS        (their eq. 7, flat space, β = 0)
//	□S      = μ₀ ∂_ν j^ν           (their eq. 21)
//
// In E, B and S (SI, signature +−−−) these are
//
//	∂B/∂t = −curl E                         div B = 0
//	∂E/∂t = c²(curl B − μ₀J − grad S)       div E = ρ/ε₀ − ∂S/∂t
//
// The second Gauss law, read as ∂S/∂t = ρ/ε₀ − div E, is the update for S,
// and with ∂E/∂t it gives □S = μ₀(∂ρ/∂t + div J): S has a source only where
// local charge conservation fails (the "extra current" I = ∂ρ/∂t + div J).
// With I = 0 and S = 0 initially, S stays 0 and the scheme is Maxwell's
// (Yee's, em.FDTD) exactly.
//
// Energy theorem (derived here; neither paper states it). From the three
// update equations, as for Poynting's theorem,
//
//	∂u/∂t + div P = −J·E + c²ρS,
//	u = ½ε₀E² + B²/(2μ₀) + S²/(2μ₀),   P = (E×B + S E)/μ₀
//
// so S carries energy (equal to an E-wave's when E = cS, as Minotti and
// Modanese state) and does work c²ρS on charge. On the Yee grid with S on
// the nodes, half a step ahead of E like B, and the walls perfect
// conductors with S = 0 on them (the Dirichlet condition, which gives the
// TM cutoffs of Pinchuk et al., arXiv 2610.28916), the scheme satisfies
// exactly, step by step,
//
//	W^{n+1} − W^n = −dt Σ J^{n+½}·Ē dV + dt c² Σ ρ̄ S^{n+½} dV,
//	W^n = Σ ½ε₀(E^n)² + B^{n−½}·B^{n+½}/(2μ₀) + S^{n−½}S^{n+½}/(2μ₀),
//
// with Ē, ρ̄ the means of the two whole steps, because the discrete grad on
// edges and div on nodes are adjoint (Σ E·grad S = −Σ S div E with S = 0 on
// the walls), as curl_e and curl_h are. The c²ρS term is booked in its own
// column (ScalarWork), the hypothetical column of CLAUDE.md §3c.
//
// # The magnetic scalar (magnetic.go)
//
// With magnetic charge allowed (Config.Magnetic, the hypothetical mode of
// em/duality.go, H2), the same construction applied to the dual potential
// ℨ gives a second scalar S_m = div ℨ + (1/c²)∂ψ_m/∂t [T], obtained from
// the electric one by the duality E → cB, cB → −E:
//
//	div B  = μ₀ρ_m − (1/c)∂S_m/∂t            curl E = −∂B/∂t − μ₀J_m − c grad S_m
//	□S_m   = (μ₀/c)(∂ρ_m/∂t + div J_m)
//
// with ρ_m [A·m/m³] and J_m [A·m/(m²·s)] as in em/duality.go. Its energy
// theorem joins the electric one:
//
//	∂u/∂t + div P = −J·E − J_m·B + c²ρS + cρ_mS_m,
//	u = ½ε₀E² + B²/(2μ₀) + (S² + S_m²)/(2μ₀),  P = (E×B + S E + c S_m B)/μ₀.
//
// S_m lives at the cell centres, where the Yee grid's div B is, at the
// whole steps with E; ρ_m at the cell centres at the half steps with B.
// The cells' discrete grad (on the B faces, with S_m = 0 outside the box)
// and div are adjoint, so the theorem holds exactly step by step, with
// the magnetic-source work in MagSourceWork and cρ_mS_m in
// MagScalarWork. Neither source paper (Minotti & Modanese) has magnetic
// charge; this sector is the duality image of theirs, derived here.
package scalar

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
)

// Mode selects the electrodynamics.
type Mode int

const (
	// AB is Aharonov–Bohm electrodynamics: S physical and evolved.
	AB Mode = iota
	// Standard is Maxwell's: S held at zero. Charge that appears without a
	// current (an extra current) then breaks Gauss's law, which GaussDefect
	// shows: the Lorenz condition dropped with nothing added to the action
	// leaves no field to carry it.
	Standard
)

func (m Mode) String() string { return [...]string{"AB", "standard"}[m] }

// Config sets up a field region.
type Config struct {
	N       [3]int     // grid nodes per axis (≥ 3)
	D       [3]float64 // node spacing [m]
	Origin  [3]float64 // position of node (0, 0, 0) [m]
	Courant float64    // dt = Courant / (c·√Σ1/D²); default 0.99, as em.FDTD
	Mode    Mode
	// Magnetic enables magnetic charge and current (hypothetical, H2):
	// MagNone (default: no magnetic sources accepted), MagDual (magnetic
	// charge with S_m held at 0: em/duality.go's equations) or MagAB (S_m
	// physical and evolved).
	Magnetic MagMode
	// Spin enables ECE2's spin connection in vacuum (hypothetical, H48;
	// spin.go): the geometric current (2/μ₀)(κ₀E/c + κ×B). Nil = Maxwell.
	Spin *SpinConnection
}

// EdgeSource is an impressed conduction current on one E edge, as
// em.EdgeSource: component Comp from node (I, J, K), total current [A].
// The charge it moves is conserved: ρ changes by −dt·div J at the two end
// nodes.
type EdgeSource struct {
	Comp    int
	I, J, K int
	Current func(t float64) float64
}

// NodeSource is an extra current at one node: charge appearing (or
// vanishing) there at Rate [C/s] with no current to bring it, I = ∂ρ/∂t +
// div J ≠ 0. It is the source of S (□S = μ₀I) and the hypothesis' input.
type NodeSource struct {
	I, J, K int
	Rate    func(t float64) float64
}

// Ledger is the energy account [J]. Residual is round-off.
type Ledger struct {
	W0, W      float64 // field energy at the start and now
	WS         float64 // S part of W (S^{n−½}S^{n+½}/2μ₀)
	SourceWork float64 // −dt Σ J·Ē dV: work of the edge sources on the field
	ScalarWork float64 // dt c² Σ ρ̄ S dV: hypothetical column (zero when S = 0)
	// Magnetic sector (hypothetical; zero unless Config.Magnetic is set).
	WSm           float64 // S_m part of W ((S_m^n)²/2μ₀)
	MagSourceWork float64 // −dt Σ B^{n+½}·(J_m^n + J_m^{n+1})/2 dV: work of magnetic currents
	MagScalarWork float64 // dt c Σ ρ_m S̄_m dV: work of S_m on magnetic charge
	// ECE2 spin connection (hypothetical; zero unless Config.Spin is set).
	SpinWork float64 // −dt Σ J_ω·Ē dV: work of the spin-connection current on the field
}

// Residual is W − W0 − SourceWork − ScalarWork − MagSourceWork −
// MagScalarWork − SpinWork.
func (l Ledger) Residual() float64 {
	return l.W - l.W0 - l.SourceWork - l.ScalarWork - l.MagSourceWork - l.MagScalarWork - l.SpinWork
}

// Region is a box of vacuum with perfectly conducting walls.
type Region struct {
	N      [3]int
	D      [3]float64
	Origin [3]float64
	Dt     float64
	Steps  int
	Mode   Mode

	E, B   [3][]float64 // E at Time(), B at Time() + dt/2 (em.FDTD layout)
	S, Rho []float64    // S at Time() + dt/2, ρ at Time() on the nodes

	Magnetic MagMode
	Sm, RhoM []float64 // S_m at Time(), ρ_m at Time() + dt/2 on the cells (index of the cell's lower corner)
	mag      magState
	spin     *spinState

	bPrev [3][]float64
	sPrev []float64
	src   []EdgeSource
	nsrc  []NodeSource
	led   Ledger
	start bool
	t0    float64
}

// ErrConfig reports an invalid set-up.
var ErrConfig = errors.New("scalar: invalid configuration")

// New builds a region with all fields zero.
func New(cfg Config) (*Region, error) {
	inv := 0.0
	for a := 0; a < 3; a++ {
		if cfg.N[a] < 3 || !(cfg.D[a] > 0) {
			return nil, fmt.Errorf("%w: axis %d needs ≥ 3 nodes and spacing > 0", ErrConfig, a)
		}
		inv += 1 / (cfg.D[a] * cfg.D[a])
	}
	s := cfg.Courant
	if s == 0 {
		s = 0.99
	}
	if !(s > 0 && s < 1) {
		return nil, fmt.Errorf("%w: Courant number %g must be in (0, 1)", ErrConfig, s)
	}
	if cfg.Mode != AB && cfg.Mode != Standard {
		return nil, fmt.Errorf("%w: unknown mode %d", ErrConfig, cfg.Mode)
	}
	if cfg.Magnetic > MagAB {
		return nil, fmt.Errorf("%w: unknown magnetic mode %d", ErrConfig, cfg.Magnetic)
	}
	r := &Region{N: cfg.N, D: cfg.D, Origin: cfg.Origin, Mode: cfg.Mode, Magnetic: cfg.Magnetic,
		Dt: s / (em.SpeedOfLight * math.Sqrt(inv))}
	n := cfg.N[0] * cfg.N[1] * cfg.N[2]
	for c := 0; c < 3; c++ {
		r.E[c] = make([]float64, n)
		r.B[c] = make([]float64, n)
		r.bPrev[c] = make([]float64, n)
	}
	r.S = make([]float64, n)
	r.Rho = make([]float64, n)
	r.sPrev = make([]float64, n)
	r.Sm = make([]float64, n)
	r.RhoM = make([]float64, n)
	if cfg.Spin != nil {
		r.setupSpin(cfg.Spin)
	}
	return r, nil
}

func (r *Region) index(i, j, k int) int { return i + r.N[0]*(j+r.N[1]*k) }

// Time is the time of E and ρ.
func (r *Region) Time() float64 { return r.t0 + float64(r.Steps)*r.Dt }

// Ledger returns the energy account.
func (r *Region) Ledger() Ledger { return r.led }

func (r *Region) dV() float64 { return r.D[0] * r.D[1] * r.D[2] }

// eFree is em.FDTD's rule: E edges tangential to a wall stay zero.
func (r *Region) eFree(c, i, j, k int) bool {
	ijk := [3]int{i, j, k}
	for a := 0; a < 3; a++ {
		if a == c {
			if ijk[a] < 0 || ijk[a] > r.N[a]-2 {
				return false
			}
		} else if ijk[a] < 1 || ijk[a] > r.N[a]-2 {
			return false
		}
	}
	return true
}

// interior reports a node off the walls (where S is evolved).
func (r *Region) interior(i, j, k int) bool {
	return i > 0 && j > 0 && k > 0 && i < r.N[0]-1 && j < r.N[1]-1 && k < r.N[2]-1
}

func (r *Region) pos(off [3]float64, i, j, k int) (x, y, z float64) {
	return r.Origin[0] + (float64(i)+off[0])*r.D[0], r.Origin[1] + (float64(j)+off[1])*r.D[1],
		r.Origin[2] + (float64(k)+off[2])*r.D[2]
}

// SetE samples E [V/m] at the edges (time Time()); wall edges stay zero.
func (r *Region) SetE(e func(x, y, z float64) [3]float64) {
	r.each(func(i, j, k, n int) {
		for c := 0; c < 3; c++ {
			if r.eFree(c, i, j, k) {
				off := [3]float64{}
				off[c] = 0.5
				x, y, z := r.pos(off, i, j, k)
				r.E[c][n] = e(x, y, z)[c]
			}
		}
	})
}

// SetS samples S [T] at the interior nodes (time Time() + dt/2).
func (r *Region) SetS(s func(x, y, z float64) float64) {
	r.each(func(i, j, k, n int) {
		if r.interior(i, j, k) {
			x, y, z := r.pos([3]float64{}, i, j, k)
			r.S[n] = s(x, y, z)
		}
	})
}

// AddSource adds an edge current.
func (r *Region) AddSource(s EdgeSource) error {
	if s.Comp < 0 || s.Comp > 2 || !r.eFree(s.Comp, s.I, s.J, s.K) || s.Current == nil {
		return fmt.Errorf("%w: edge source not on a free edge", ErrConfig)
	}
	r.src = append(r.src, s)
	return nil
}

// AddNodeSource adds an extra current at an interior node.
func (r *Region) AddNodeSource(s NodeSource) error {
	if !r.interior(s.I, s.J, s.K) || s.Rate == nil {
		return fmt.Errorf("%w: node source not at an interior node", ErrConfig)
	}
	r.nsrc = append(r.nsrc, s)
	return nil
}

func (r *Region) each(fn func(i, j, k, n int)) {
	for k := 0; k < r.N[2]; k++ {
		for j := 0; j < r.N[1]; j++ {
			for i := 0; i < r.N[0]; i++ {
				fn(i, j, k, r.index(i, j, k))
			}
		}
	}
}

// divE is the discrete divergence of E at node n (interior).
func (r *Region) divE(n int) float64 {
	sy, sz := r.N[0], r.N[0]*r.N[1]
	return (r.E[0][n]-r.E[0][n-1])/r.D[0] + (r.E[1][n]-r.E[1][n-sy])/r.D[1] + (r.E[2][n]-r.E[2][n-sz])/r.D[2]
}

// curlE is component c of curl_e(E) at B face n (em.FDTD's stencil).
func (r *Region) curlE(c, n int) float64 {
	sy, sz := r.N[0], r.N[0]*r.N[1]
	e, dx, dy, dz := r.E, r.D[0], r.D[1], r.D[2]
	switch c {
	case 0:
		return (e[2][n+sy]-e[2][n])/dy - (e[1][n+sz]-e[1][n])/dz
	case 1:
		return (e[0][n+sz]-e[0][n])/dz - (e[2][n+1]-e[2][n])/dx
	default:
		return (e[1][n+1]-e[1][n])/dx - (e[0][n+sy]-e[0][n])/dy
	}
}

// faceRange is the index range of B component c (em.FDTD's).
func (r *Region) faceRange(c int) (imax, jmax, kmax int) {
	imax, jmax, kmax = r.N[0]-1, r.N[1]-1, r.N[2]-1
	switch c {
	case 0:
		imax = r.N[0]
	case 1:
		jmax = r.N[1]
	default:
		kmax = r.N[2]
	}
	return
}

func (r *Region) begin() {
	r.start = true
	dt := r.Dt
	// B^{−½} = B^{½} + dt curl E^0;  S^{−½} = S^{½} − dt(ρ⁰/ε₀ − div E⁰).
	for c := 0; c < 3; c++ {
		imax, jmax, kmax := r.faceRange(c)
		for k := 0; k < kmax; k++ {
			for j := 0; j < jmax; j++ {
				for i := 0; i < imax; i++ {
					n := r.index(i, j, k)
					r.bPrev[c][n] = r.B[c][n] + dt*(r.curlE(c, n)+em.SpeedOfLight*r.gradSm(c, i, j, k))
				}
			}
		}
	}
	if r.Mode == AB {
		r.each(func(i, j, k, n int) {
			if r.interior(i, j, k) {
				r.sPrev[n] = r.S[n] - dt*(r.Rho[n]/em.Eps0-r.divE(n))
			}
		})
	}
	w, ws := r.energy()
	wsm := r.energySm()
	r.led = Ledger{W0: w + wsm, W: w + wsm, WS: ws, WSm: wsm}
}

// energy is W^n and its S part from the stored copies.
func (r *Region) energy() (w, ws float64) {
	var we, wm float64
	r.each(func(i, j, k, n int) {
		for c := 0; c < 3; c++ {
			we += 0.5 * em.Eps0 * r.E[c][n] * r.E[c][n]
			wm += r.bPrev[c][n] * r.B[c][n] / (2 * em.Mu0)
		}
		ws += r.sPrev[n] * r.S[n] / (2 * em.Mu0)
	})
	dV := r.dV()
	return (we + wm + ws) * dV, ws * dV
}

// Step advances one time step: E and ρ to Time()+dt, then S and B half a
// step beyond.
func (r *Region) Step() {
	if !r.start {
		r.begin()
	}
	dt, dV := r.Dt, r.dV()
	tHalf := r.Time() + dt/2
	sy, sz := r.N[0], r.N[0]*r.N[1]
	step := [3]int{1, sy, sz}
	cb, inv := dt/em.Eps0, 1/em.Mu0

	// Ampère: E^{n+1} = E^n + (dt/ε₀)(curl_h B/μ₀ − J − grad S/μ₀).
	eOld := [3][]float64{}
	var we float64
	for c := 0; c < 3; c++ {
		eOld[c] = append([]float64(nil), r.E[c]...)
	}
	r.each(func(i, j, k, n int) {
		for c := 0; c < 3; c++ {
			if !r.eFree(c, i, j, k) {
				continue
			}
			b := r.B
			var curl float64
			switch c {
			case 0:
				curl = inv * ((b[2][n]-b[2][n-sy])/r.D[1] - (b[1][n]-b[1][n-sz])/r.D[2])
			case 1:
				curl = inv * ((b[0][n]-b[0][n-sz])/r.D[2] - (b[2][n]-b[2][n-1])/r.D[0])
			default:
				curl = inv * ((b[1][n]-b[1][n-1])/r.D[0] - (b[0][n]-b[0][n-sy])/r.D[1])
			}
			e := r.E[c][n] + cb*curl
			if r.Mode == AB {
				e -= cb * inv * (r.S[n+step[c]] - r.S[n]) / r.D[c]
			}
			r.E[c][n] = e
		}
	})
	rhoOld := append([]float64(nil), r.Rho...)
	srcJ := make([]float64, 0, len(r.src))
	for _, s := range r.src {
		n := r.index(s.I, s.J, s.K)
		area := 1.0
		for a := 0; a < 3; a++ {
			if a != s.Comp {
				area *= r.D[a]
			}
		}
		j := s.Current(tHalf) / area
		r.E[s.Comp][n] -= cb * j
		srcJ = append(srcJ, j)
		// conserved: the charge leaves node n and arrives at its neighbour.
		q := dt * j / r.D[s.Comp]
		r.Rho[n] -= q
		r.Rho[n+step[s.Comp]] += q
	}
	for _, s := range r.nsrc {
		r.Rho[r.index(s.I, s.J, s.K)] += dt * s.Rate(tHalf) / dV
	}
	// Spin-connection current (spin.go), then the work of every current
	// against the same final Ē, so the energy theorem stays exact.
	var jw [3][]float64
	if r.spin != nil {
		jw = r.applySpin(eOld, tHalf)
		var w float64
		for c := 0; c < 3; c++ {
			for n, j := range jw[c] {
				if j != 0 {
					w += j * 0.5 * (eOld[c][n] + r.E[c][n])
				}
			}
		}
		r.led.SpinWork -= dt * w * dV
	}
	for i, s := range r.src {
		n := r.index(s.I, s.J, s.K)
		r.led.SourceWork -= dt * srcJ[i] * 0.5 * (eOld[s.Comp][n] + r.E[s.Comp][n]) * dV
	}
	for c := 0; c < 3; c++ {
		for _, e := range r.E[c] {
			we += 0.5 * em.Eps0 * e * e
		}
	}

	// S^{n+3/2} = S^{n+½} + dt(ρ^{n+1}/ε₀ − div E^{n+1}); work c²ρ̄S^{n+½}.
	var ws float64
	if r.Mode == AB {
		c2 := em.SpeedOfLight * em.SpeedOfLight
		var work float64
		r.each(func(i, j, k, n int) {
			if !r.interior(i, j, k) {
				return
			}
			old := r.S[n]
			work += 0.5 * (rhoOld[n] + r.Rho[n]) * old
			nw := old + dt*(r.Rho[n]/em.Eps0-r.divE(n))
			r.sPrev[n] = old
			r.S[n] = nw
			ws += old * nw / (2 * em.Mu0)
		})
		r.led.ScalarWork += dt * c2 * work * dV
	}

	// Magnetic scalar: S_m^{n+1} = S_m^n + c dt(μ₀ρ_m^{n+½} − div B^{n+½}).
	wsm := r.stepSm()

	// Faraday: B^{n+3/2} = B^{n+½} − dt(curl_e E^{n+1} + c grad S_m^{n+1})
	// − dt μ₀J_m^{n+1} on the magnetic-source faces (magnetic.go).
	csm := 0.0
	if r.Magnetic == MagAB {
		csm = em.SpeedOfLight
	}
	var wm float64
	for c := 0; c < 3; c++ {
		imax, jmax, kmax := r.faceRange(c)
		for k := 0; k < kmax; k++ {
			for j := 0; j < jmax; j++ {
				for i := 0; i < imax; i++ {
					n := r.index(i, j, k)
					old := r.B[c][n]
					nw := old - dt*r.curlE(c, n)
					if csm != 0 {
						nw -= dt * csm * r.gradSm(c, i, j, k)
					}
					r.bPrev[c][n] = old
					r.B[c][n] = nw
					wm += old * nw / (2 * em.Mu0)
				}
			}
		}
	}
	wm += r.applyMagSources()
	r.Steps++
	r.led.W = (we+wm+ws)*dV + wsm
	r.led.WS = ws * dV
	r.led.WSm = wsm
}

// Run steps until Time() ≥ t.
func (r *Region) Run(t float64) {
	for r.Time() < t-1e-9*r.Dt {
		r.Step()
	}
}

// SAt returns S at node (i, j, k), time Time() + dt/2.
func (r *Region) SAt(i, j, k int) float64 { return r.S[r.index(i, j, k)] }

// GaussDefect is max over interior nodes of |ρ/ε₀ − div E − ∂S/∂t| at the
// time of E, relative to the largest |ρ/ε₀| or |div E| there. In AB mode it
// is round-off by construction; in standard mode an extra current makes it
// O(1): Maxwell's equations cannot hold charge that arrives without a
// current.
func (r *Region) GaussDefect() float64 {
	var worst, scale float64
	r.each(func(i, j, k, n int) {
		if !r.interior(i, j, k) {
			return
		}
		g := r.Rho[n]/em.Eps0 - r.divE(n)
		if r.Mode == AB && r.start {
			g -= (r.S[n] - r.sPrev[n]) / r.Dt
		}
		worst = math.Max(worst, math.Abs(g))
		scale = math.Max(scale, math.Max(math.Abs(r.Rho[n]/em.Eps0), math.Abs(r.divE(n))))
	})
	if scale == 0 {
		return worst
	}
	return worst / scale
}
