// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"
)

// Machine is a set of n windings whose inductance matrix depends on a
// mechanical coordinate θ, with a mechanical port (kronspice N element):
// Kron's rotating electric machine (G. Kron, Equivalent Circuits of
// Electric Machinery, Wiley 1951, Epilogue pp. 253–259).
//
// Nodes = winding 1 (+, −); Extra = windings 2…n (+, −) then the
// mechanical port (+, −), whose across variable is the angular velocity
// ω = dθ/dt [rad/s] and whose through variable is torque [N·m]. Model
// names a ".model … machine" card (MachineSpec).
//
// Winding k obeys
//
//	v_k = d/dt (Σ_m L_km(θ)·i_m) + ω·Σ_m G_km·i_m
//
// and the element applies the torque
//
//	T = ½·iᵀ·(∂L/∂θ)·i + iᵀ·G·i
//
// to its mechanical port. With G = 0 this is Kron's holonomic frame I
// (windings fixed to the iron, e = Ri + p(Li), Epilogue eq. 2); with L
// constant and G ≠ 0 it is his frame II(b) (axes not fixed to the iron,
// e = Ri + L·pi + G·i·pθ, eq. 9), as produced by a rotating-frame
// transformation of a frame-I machine (kron.RotatingFrame). Either way the
// power into the windings is the change of ½·iᵀ·L·i plus ω·T, so the
// element is lossless and its energy balance closes exactly. A coordinate
// driven at a prescribed speed makes the same element a parametric
// (time-varying) inductance, the pump power entering at the port.
const Machine Kind = 'N'

// MaxWindings bounds the windings of one machine.
const MaxWindings = 64

// MaxHarmonic bounds the harmonic order of a machine's inductances.
const MaxHarmonic = 64

// MachineSpec is the parameter set of a machine:
//
//	L(θ) = L0 + Σ_h ( C_h·cos(k_h·θ) + S_h·sin(k_h·θ) )
//
// with every matrix symmetric, a constant speed matrix G (not necessarily
// symmetric), and the angle θ0 at t = 0.
//
// In a netlist: ".model name machine (n=… theta0=… l_i_j=… lc_k_i_j=…
// ls_k_i_j=… g_i_j=… ic_k=…)", windings counted from 1; l_, lc_ and ls_
// entries set both (i, j) and (j, i); ic_k is winding k's current at the
// start of a UIC transient.
type MachineSpec struct {
	N         int
	Theta0    float64
	L0        [][]float64
	Harmonics []MachineHarmonic
	G         [][]float64
	// IC are the winding currents at the start of a UIC transient
	// (default 0, like an inductor without IC=).
	IC []float64
}

// MachineHarmonic is one harmonic of the inductance matrix.
type MachineHarmonic struct {
	K        int
	Cos, Sin [][]float64
}

func squareZero(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	return m
}

// NewMachineSpec returns an empty spec of n windings.
func NewMachineSpec(n int) *MachineSpec {
	return &MachineSpec{N: n, L0: squareZero(n), G: squareZero(n), IC: make([]float64, n)}
}

// Harmonic returns the harmonic of order k, adding it if needed.
func (s *MachineSpec) Harmonic(k int) *MachineHarmonic {
	for i := range s.Harmonics {
		if s.Harmonics[i].K == k {
			return &s.Harmonics[i]
		}
	}
	s.Harmonics = append(s.Harmonics, MachineHarmonic{K: k, Cos: squareZero(s.N), Sin: squareZero(s.N)})
	sort.Slice(s.Harmonics, func(a, b int) bool { return s.Harmonics[a].K < s.Harmonics[b].K })
	return s.Harmonic(k)
}

// Inductance writes L(θ), ∂L/∂θ and ∂²L/∂θ² into l, d1, d2 (any may be
// nil).
func (s *MachineSpec) Inductance(theta float64, l, d1, d2 [][]float64) {
	n := s.N
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if l != nil {
				l[i][j] = s.L0[i][j]
			}
			if d1 != nil {
				d1[i][j] = 0
			}
			if d2 != nil {
				d2[i][j] = 0
			}
		}
	}
	for _, h := range s.Harmonics {
		k := float64(h.K)
		c, sn := math.Cos(k*theta), math.Sin(k*theta)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				a, b := h.Cos[i][j], h.Sin[i][j]
				if a == 0 && b == 0 {
					continue
				}
				if l != nil {
					l[i][j] += a*c + b*sn
				}
				if d1 != nil {
					d1[i][j] += k * (-a*sn + b*c)
				}
				if d2 != nil {
					d2[i][j] += -k * k * (a*c + b*sn)
				}
			}
		}
	}
}

// Params returns the spec as .model parameters (the inverse of
// ParseMachine).
func (s *MachineSpec) Params() map[string]float64 {
	p := map[string]float64{"n": float64(s.N)}
	if s.Theta0 != 0 {
		p["theta0"] = s.Theta0
	}
	for k, v := range s.IC {
		if v != 0 {
			p[fmt.Sprintf("ic_%d", k+1)] = v
		}
	}
	for i := 0; i < s.N; i++ {
		for j := i; j < s.N; j++ {
			if v := s.L0[i][j]; v != 0 {
				p[fmt.Sprintf("l_%d_%d", i+1, j+1)] = v
			}
			for _, h := range s.Harmonics {
				if v := h.Cos[i][j]; v != 0 {
					p[fmt.Sprintf("lc_%d_%d_%d", h.K, i+1, j+1)] = v
				}
				if v := h.Sin[i][j]; v != 0 {
					p[fmt.Sprintf("ls_%d_%d_%d", h.K, i+1, j+1)] = v
				}
			}
		}
		for j := 0; j < s.N; j++ {
			if v := s.G[i][j]; v != 0 {
				p[fmt.Sprintf("g_%d_%d", i+1, j+1)] = v
			}
		}
	}
	return p
}

// ParseMachine reads a machine .model card and checks it: windings 1 to
// MaxWindings, harmonic orders 1 to MaxHarmonic, finite values, and L(θ)
// positive definite at 4·(highest order)+8 angles over a turn (the stored
// energy ½·iᵀ·L·i must be positive for every current).
func ParseMachine(m Model) (*MachineSpec, error) {
	bad := func(format string, a ...interface{}) error {
		return fmt.Errorf("%w: machine %s: %s", ErrBadModel, m.Name, fmt.Sprintf(format, a...))
	}
	nv, ok := m.Params["n"]
	if !ok || nv != math.Trunc(nv) || nv < 1 || nv > MaxWindings {
		return nil, bad("N must be a whole number of windings, 1 to %d", MaxWindings)
	}
	s := NewMachineSpec(int(nv))
	idx := func(t string) (int, error) {
		v, err := strconv.Atoi(t)
		if err != nil || v < 1 || v > s.N {
			return 0, bad("winding %q out of range 1…%d", t, s.N)
		}
		return v - 1, nil
	}
	keys := make([]string, 0, len(m.Params))
	for k := range m.Params {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		v := m.Params[k]
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, bad("%s is not finite", k)
		}
		f := strings.Split(k, "_")
		switch {
		case k == "n":
		case k == "theta0":
			s.Theta0 = v
		case f[0] == "ic" && len(f) == 2:
			i, err := idx(f[1])
			if err != nil {
				return nil, err
			}
			s.IC[i] = v
		case (f[0] == "l" || f[0] == "g") && len(f) == 3:
			i, err := idx(f[1])
			if err != nil {
				return nil, err
			}
			j, err := idx(f[2])
			if err != nil {
				return nil, err
			}
			if f[0] == "g" {
				s.G[i][j] = v
			} else {
				s.L0[i][j], s.L0[j][i] = v, v
			}
		case (f[0] == "lc" || f[0] == "ls") && len(f) == 4:
			h, err := strconv.Atoi(f[1])
			if err != nil || h < 1 || h > MaxHarmonic {
				return nil, bad("%s: harmonic order must be 1 to %d", k, MaxHarmonic)
			}
			i, err := idx(f[2])
			if err != nil {
				return nil, err
			}
			j, err := idx(f[3])
			if err != nil {
				return nil, err
			}
			hm := s.Harmonic(h)
			if f[0] == "lc" {
				hm.Cos[i][j], hm.Cos[j][i] = v, v
			} else {
				hm.Sin[i][j], hm.Sin[j][i] = v, v
			}
		default:
			return nil, bad("unknown parameter %q (n, theta0, l_i_j, lc_k_i_j, ls_k_i_j, g_i_j, ic_k)", k)
		}
	}
	if err := s.CheckPositive(); err != nil {
		return nil, bad("%v", err)
	}
	return s, nil
}

// CheckPositive verifies that L(θ) is positive definite at 4·K+8 angles
// over one turn (K the highest harmonic order), by Cholesky factorization.
func (s *MachineSpec) CheckPositive() error {
	kmax := 0
	for _, h := range s.Harmonics {
		if h.K > kmax {
			kmax = h.K
		}
	}
	steps := 4*kmax + 8
	l := squareZero(s.N)
	for a := 0; a < steps; a++ {
		th := 2 * math.Pi * float64(a) / float64(steps)
		s.Inductance(th, l, nil, nil)
		if !cholesky(l) {
			return fmt.Errorf("inductance matrix is not positive definite at θ = %.4g rad", th)
		}
	}
	return nil
}

// cholesky reports whether symmetric a is positive definite (a is
// overwritten).
func cholesky(a [][]float64) bool {
	n := len(a)
	for j := 0; j < n; j++ {
		d := a[j][j]
		for k := 0; k < j; k++ {
			d -= a[j][k] * a[j][k]
		}
		if !(d > 0) {
			return false
		}
		d = math.Sqrt(d)
		a[j][j] = d
		for i := j + 1; i < n; i++ {
			s := a[i][j]
			for k := 0; k < j; k++ {
				s -= a[i][k] * a[j][k]
			}
			a[i][j] = s / d
		}
	}
	return true
}

// MachineWindings returns a machine element's winding node pairs and its
// mechanical port.
func MachineWindings(e Element) (windings [][2]int, port [2]int) {
	windings = append(windings, e.Nodes)
	n := len(e.Extra)
	for i := 0; i+1 < n-2; i += 2 {
		windings = append(windings, [2]int{e.Extra[i], e.Extra[i+1]})
	}
	if n >= 2 {
		port = [2]int{e.Extra[n-2], e.Extra[n-1]}
	}
	return windings, port
}

// MachineFor returns a machine element's spec.
func (c *Circuit) MachineFor(e Element) (*MachineSpec, error) {
	m, ok := c.ModelFor(e)
	if !ok {
		return nil, fmt.Errorf("%w: %s uses undefined model %q", ErrBadModel, e.Name, e.Model)
	}
	if m.Type != "machine" {
		return nil, fmt.Errorf("%w: %s needs a MACHINE model, %q is type %s", ErrBadModel, e.Name, e.Model, strings.ToUpper(m.Type))
	}
	return ParseMachine(m)
}

// Machine adds a machine element: windings as consecutive (+, −) node
// pairs, then the mechanical port (+, −), and a machine model defined
// with Model or MachineModel.
func (b *Builder) Machine(name, model string, windings [][2]string, portP, portN string) *Builder {
	if len(windings) == 0 {
		windings = [][2]string{{"0", "0"}} // rejected by Validate
	}
	var extra []string
	for _, w := range windings[1:] {
		extra = append(extra, w[0], w[1])
	}
	extra = append(extra, portP, portN)
	return b.AddMulti(Element{Name: name, Kind: Machine, Model: strings.ToLower(model)}, windings[0][0], windings[0][1], extra...)
}

// MachineModel defines a machine model from a spec.
func (b *Builder) MachineModel(name string, s *MachineSpec) *Builder {
	return b.Model(name, "machine", s.Params())
}

func (c *Circuit) validateMachine(e Element) error {
	s, err := c.MachineFor(e)
	if err != nil {
		return err
	}
	if len(e.Extra) != 2*s.N {
		return fmt.Errorf("%w: %s: %d windings need %d nodes, got %d", ErrBadValue, e.Name, s.N, 2*s.N+2, len(e.Extra)+2)
	}
	w, port := MachineWindings(e)
	for k, p := range w {
		if p[0] == p[1] {
			return fmt.Errorf("%w: %s winding %d", ErrSelfConnected, e.Name, k+1)
		}
	}
	if port[0] == port[1] {
		return fmt.Errorf("%w: %s mechanical port", ErrSelfConnected, e.Name)
	}
	return nil
}
