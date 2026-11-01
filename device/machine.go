// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Machine is the kronspice N element: n windings whose inductance matrix
// L(θ) depends on a mechanical coordinate θ, a speed matrix G, and a
// mechanical port of angular velocity ω = dθ/dt (across) and torque
// (through). See circuit.Machine for the equations.
//
// Unknowns: one branch current per winding (current n+ → n− through the
// winding) and the angle θ. Equations, at a transient step of length h
// (trapezoidal: c = 2/h and the previous dψ/dt carried; backward Euler:
// c = 1/h):
//
//	winding k:  V_k = c·(ψ_k − ψ_k,prev) − hist_k + ω·(G·i)_k,  ψ = L(θ)·i
//	angle:      θ = θ_prev + (h/2)·(ω + ω_prev)      (backward Euler: + h·ω)
//	port:       the element drives torque T = ½·iᵀ·L′(θ)·i + iᵀ·G·i into
//	            the port's + node (current −T from + to − through it)
//
// The products L(θ)·i, ω·G·i and the torque make the element nonlinear; it
// is linearized at every Newton iteration. At an operating point the
// currents are steady (di/dt = 0) at the angle θ0:
// V_k = ω·((L′(θ0) + G)·i)_k. At the start of a UIC transient the winding
// currents are their initial values (spec IC, default 0), as for an
// inductor.
//
// Maxwell derivation (CLAUDE.md §3c, inductor row): the flux 𝔅 = V.∇𝔄 of
// every winding linking every other, through a geometry that moves with
// θ; the force on the moving iron and copper is the Lorentz force on the
// currents, whose work is the ω·T term. Kron's dynamical equations of
// rotating machinery (G. Kron, Equivalent Circuits of Electric Machinery,
// Wiley 1951, Epilogue eqs. 2, 7 and 9, pp. 254–257).
//
// Driven at a prescribed speed the element is a parametric inductance.
// Kouřil and Vrba note that in a circuit with time-varying elements the
// energy delivered to the loads and turned into heat comes partly from
// the source that changes the parameter (Non-linear and Parametric
// Circuits, 1988, §4.4.1 p. 129); here that source is the mechanical port,
// so the pump's share is on the ledger. Its powers obey the Manley–Rowe
// relations (ibid. eqs. 4.95, 4.101; analysis.TestManleyRowePumpedInductor).
type Machine struct {
	Base
	spec  *circuit.MachineSpec
	n     int
	extra []int
	w     [][2]int // winding unknowns
	p     [2]int   // port unknowns
	br    []int    // winding branch unknowns
	th    int      // θ unknown

	// work space (allocated once)
	l, d1, d2          [][]float64
	i, li, d1i, gi, ti []float64

	// accepted state
	thPrev, omPrev float64
	psiPrev, hist  []float64 // ψ and dψ/dt at the previous accepted point
	hpsi           []history
	hth            history
	// linearization of the torque at the last Load (Operating reports the
	// port current the linear solve used, so the ledger's KCL closes)
	linOK          bool
	linT, linTh0   float64
	linDTh         float64
	linI0, linDTdi []float64
	// most recent Operating
	iOp     []float64
	thOp    float64
	omOp    float64
	acOmega float64
}

// NewMachine builds a machine element.
func NewMachine(e circuit.Element, s *circuit.MachineSpec) *Machine {
	n := s.N
	m := &Machine{Base: newBase(e), spec: s, n: n, extra: append([]int(nil), e.Extra...)}
	sq := func() [][]float64 {
		a := make([][]float64, n)
		for i := range a {
			a[i] = make([]float64, n)
		}
		return a
	}
	m.l, m.d1, m.d2 = sq(), sq(), sq()
	m.i, m.li, m.d1i, m.gi, m.ti = make([]float64, n), make([]float64, n), make([]float64, n), make([]float64, n), make([]float64, n)
	m.psiPrev, m.hist, m.iOp = make([]float64, n), make([]float64, n), make([]float64, n)
	m.hpsi = make([]history, n)
	m.linI0, m.linDTdi = make([]float64, n), make([]float64, n)
	m.thPrev, m.thOp = s.Theta0, s.Theta0
	_ = m.SetIntrinsic(exchange.Electrical, exchange.Both,
		"windings connect inside the network; stores ½iᵀL(θ)i")
	_ = m.SetIntrinsic(exchange.Mechanical, exchange.Both,
		"electromagnetic torque on the shaft: energy converted inside the network")
	return m
}

// Setup maps the windings and the port and allocates n + 1 unknowns.
func (m *Machine) Setup(ctx *SetupContext) error {
	m.setupNodes(ctx)
	e := circuit.Element{Nodes: m.nodes, Extra: m.extra}
	w, port := circuit.MachineWindings(e)
	if len(w) != m.n || len(m.extra) != 2*m.n {
		return fmt.Errorf("%w: %s: %d windings need %d nodes", ErrUnsupported, m.name, m.n, 2*m.n+2)
	}
	m.w = make([][2]int, len(w))
	for k, p := range w {
		m.w[k] = [2]int{ctx.Unknown(p[0]), ctx.Unknown(p[1])}
	}
	m.p = [2]int{ctx.Unknown(port[0]), ctx.Unknown(port[1])}
	m.br = make([]int, m.n)
	for k := range m.br {
		m.br[k] = ctx.AllocBranch()
	}
	m.th = ctx.AllocBranch()
	return nil
}

// omega is the port's angular velocity at x.
func (m *Machine) omega(x []float64) float64 { return at(x, m.p[0]) - at(x, m.p[1]) }

// eval fills L, L′, L″, i and the products at x.
func (m *Machine) eval(x []float64, theta float64) {
	m.spec.Inductance(theta, m.l, m.d1, m.d2)
	for k := 0; k < m.n; k++ {
		m.i[k] = at(x, m.br[k])
	}
	g := m.spec.G
	for a := 0; a < m.n; a++ {
		var s, s1, sg, st float64
		for b := 0; b < m.n; b++ {
			s += m.l[a][b] * m.i[b]
			s1 += m.d1[a][b] * m.i[b]
			sg += g[a][b] * m.i[b]
			st += (g[a][b] + g[b][a]) * m.i[b]
		}
		m.li[a], m.d1i[a], m.gi[a], m.ti[a] = s, s1, sg, st
	}
}

// torque is T = ½·iᵀ·L′·i + iᵀ·G·i for the current eval.
func (m *Machine) torque() float64 {
	t := 0.0
	for a := 0; a < m.n; a++ {
		t += m.i[a] * (0.5*m.d1i[a] + m.gi[a])
	}
	return t
}

// Load linearizes and stamps the element.
func (m *Machine) Load(ctx *LoadContext) error {
	x := ctx.X
	om := m.omega(x)
	var theta float64
	steady := false
	switch ctx.Mode {
	case ModeOP, ModeTranOP, ModeUICInit:
		theta, steady = m.spec.Theta0, true
	default:
		theta = at(x, m.th)
	}
	m.eval(x, theta)
	pp, pn := m.p[0], m.p[1]
	uic := ctx.Mode == ModeUICInit
	// Winding KCL: i_k leaves w+ and enters w−.
	for k := 0; k < m.n; k++ {
		ctx.M.Add(m.w[k][0], m.br[k], 1)
		ctx.M.Add(m.w[k][1], m.br[k], -1)
	}
	// Angle row.
	switch {
	case steady:
		ctx.M.Add(m.th, m.th, 1)
		ctx.AddRHS(m.th, m.spec.Theta0)
	case ctx.Method == Trapezoidal:
		ctx.M.Add(m.th, m.th, 1)
		ctx.M.Add(m.th, pp, -ctx.H/2)
		ctx.M.Add(m.th, pn, ctx.H/2)
		ctx.AddRHS(m.th, m.thPrev+ctx.H/2*m.omPrev)
	default:
		ctx.M.Add(m.th, m.th, 1)
		ctx.M.Add(m.th, pp, -ctx.H)
		ctx.M.Add(m.th, pn, ctx.H)
		ctx.AddRHS(m.th, m.thPrev)
	}
	// Winding rows, linearized: F_k = V_k − c·ψ_k − ω·(G i)_k − ω·s·(L′ i)_k + r_k = 0,
	// where s = 1 at a steady point (V = ω·(L′+G)·i) and c = 0 there.
	c, hist := 0.0, 0.0
	if !steady {
		c = 1 / ctx.H
		if ctx.Method == Trapezoidal {
			c = 2 / ctx.H
		}
	}
	g := m.spec.G
	for k := 0; k < m.n; k++ {
		r := m.br[k]
		if uic {
			// Winding currents pinned to their initial values (as an
			// inductor at UIC); θ = θ0.
			ctx.M.Add(r, r, 1)
			ctx.AddRHS(r, m.spec.IC[k])
			continue
		}
		ctx.M.Add(r, m.w[k][0], 1)
		ctx.M.Add(r, m.w[k][1], -1)
		// F(x0) and its partials; stamp ∂F and RHS = ∂F·x0 − F(x0).
		var f, dx float64 // F without the node terms, and Σ ∂F·x0 over non-node unknowns
		if steady {
			f = -om * (m.d1i[k] + m.gi[k])
			for b := 0; b < m.n; b++ {
				d := -om * (m.d1[k][b] + g[k][b])
				ctx.M.Add(r, m.br[b], d)
				dx += d * m.i[b]
			}
			dom := -(m.d1i[k] + m.gi[k])
			ctx.M.Add(r, pp, dom)
			ctx.M.Add(r, pn, -dom)
			dx += dom * om
		} else {
			if ctx.Method == Trapezoidal {
				hist = m.hist[k]
			}
			f = -c*(m.li[k]-m.psiPrev[k]) + hist - om*m.gi[k]
			for b := 0; b < m.n; b++ {
				d := -c*m.l[k][b] - om*g[k][b]
				ctx.M.Add(r, m.br[b], d)
				dx += d * m.i[b]
			}
			dth := -c * m.d1i[k]
			ctx.M.Add(r, m.th, dth)
			dx += dth * theta
			dom := -m.gi[k]
			ctx.M.Add(r, pp, dom)
			ctx.M.Add(r, pn, -dom)
			dx += dom * om
		}
		ctx.AddRHS(r, dx-f)
	}
	// Port: current −T leaves p+ (torque T drives the + node).
	t := m.torque()
	m.linOK, m.linT, m.linTh0, m.linDTh = true, t, theta, 0
	copy(m.linI0, m.i)
	var dx float64
	for b := 0; b < m.n; b++ {
		d := -(m.d1i[b] + m.ti[b])
		m.linDTdi[b] = -d
		ctx.M.Add(pp, m.br[b], d)
		ctx.M.Add(pn, m.br[b], -d)
		dx += d * m.i[b]
	}
	if !steady {
		var tt float64
		for a := 0; a < m.n; a++ {
			for b := 0; b < m.n; b++ {
				tt += m.i[a] * m.d2[a][b] * m.i[b]
			}
		}
		d := -0.5 * tt
		m.linDTh = -d
		ctx.M.Add(pp, m.th, d)
		ctx.M.Add(pn, m.th, -d)
		dx += d * theta
	}
	ctx.AddRHS(pp, dx-(-t))
	ctx.AddRHS(pn, -(dx - (-t)))
	return nil
}

// Limited: the machine has no junction limiting.
func (m *Machine) Limited() bool { return false }

// ConvTest: convergence of the unknowns is enough (the element's currents
// are unknowns themselves).
func (m *Machine) ConvTest([]float64, float64, float64) bool { return true }

// Operating returns winding 1's voltage and current and the element's
// total absorbed power Σ v·i − ω·T, which all goes into stored energy.
func (m *Machine) Operating(x []float64) OperatingPoint {
	m.omOp = m.omega(x)
	m.thOp = at(x, m.th)
	m.eval(x, m.thOp)
	copy(m.iOp, m.i)
	p := -m.omOp * m.portTorque(x)
	for k := 0; k < m.n; k++ {
		p += (at(x, m.w[k][0]) - at(x, m.w[k][1])) * m.i[k]
	}
	return OperatingPoint{Voltage: at(x, m.w[0][0]) - at(x, m.w[0][1]), Current: m.i[0], PowerAbsorbed: p, StoragePower: p}
}

// portTorque is the torque the last linear solve used: the linearization
// of the last Load evaluated at x (the exact torque before any Load).
func (m *Machine) portTorque(x []float64) float64 {
	if !m.linOK {
		return m.torque()
	}
	t := m.linT
	for b := 0; b < m.n; b++ {
		t += m.linDTdi[b] * (at(x, m.br[b]) - m.linI0[b])
	}
	if m.linDTh != 0 {
		t += m.linDTh * (at(x, m.th) - m.linTh0)
	}
	return t
}

// Torque returns the electromagnetic torque at solution x [N·m].
func (m *Machine) Torque(x []float64) float64 {
	m.eval(x, at(x, m.th))
	return m.torque()
}

// Angle returns θ at solution x [rad].
func (m *Machine) Angle(x []float64) float64 { return at(x, m.th) }

// Accept commits ψ, dψ/dt, θ and ω. dψ/dt is taken from the solved
// winding equation, V_k − ω·(G·i)_k, so it is the derivative the
// integration formula produced, whatever the step and method.
func (m *Machine) Accept(x []float64, t float64, mode Mode) {
	om := m.omega(x)
	theta := at(x, m.th)
	if mode == ModeUICInit {
		theta = m.spec.Theta0
	}
	m.eval(x, theta)
	reset := mode == ModeOP || mode == ModeTranOP || mode == ModeUICInit
	for k := 0; k < m.n; k++ {
		m.hist[k] = at(x, m.w[k][0]) - at(x, m.w[k][1]) - om*m.gi[k]
		m.psiPrev[k] = m.li[k]
		if reset {
			m.hpsi[k].reset()
		}
		m.hpsi[k].push(t, m.li[k])
	}
	if reset {
		m.hth.reset()
	}
	m.thPrev, m.omPrev = theta, om
	m.hth.push(t, theta)
}

// StoredEnergy is ½·iᵀ·L(θ)·i at the most recent Operating.
func (m *Machine) StoredEnergy(OperatingPoint) float64 {
	m.spec.Inductance(m.thOp, m.l, nil, nil)
	e := 0.0
	for a := 0; a < m.n; a++ {
		for b := 0; b < m.n; b++ {
			e += m.iOp[a] * m.l[a][b] * m.iOp[b]
		}
	}
	return 0.5 * e
}

// Truncate limits the step by the local truncation error of every
// winding's flux and of the angle.
func (m *Machine) Truncate(x []float64, tc *TruncContext) float64 {
	theta := at(x, m.th)
	m.eval(x, theta)
	om := m.omega(x)
	step := math.Inf(1)
	for k := 0; k < m.n; k++ {
		dpsi := at(x, m.w[k][0]) - at(x, m.w[k][1]) - om*m.gi[k]
		step = math.Min(step, lteStep(&m.hpsi[k], m.li[k], dpsi, m.hist[k], tc.VnTol, tc))
	}
	return math.Min(step, lteStep(&m.hth, theta, om, m.omPrev, tc.VnTol, tc))
}

// ExteriorPower: the machine exchanges nothing with the exterior; its
// conversion between windings and shaft happens inside the network.
func (m *Machine) ExteriorPower(exchange.ChannelID, OperatingPoint) float64 { return 0 }

// LoadAC stamps the small-signal model about the operating point (θ0, ω0,
// I0): δV = jΩ·(L·δi + L′·I0·δθ) + ω0·G·δi + (G·I0)·δω, jΩ·δθ = δω, and
// δT = ((L′ + G + Gᵀ)·I0)ᵀ·δi + ½·I0ᵀ·L″·I0·δθ. It is exact when the
// operating point is stationary (ω0 = 0, or L′·I0 = 0 as in frame II).
func (m *Machine) LoadAC(ctx *ACContext) error {
	jw := complex(0, ctx.Omega)
	m.acOmega = ctx.Omega
	th, om := m.thOp, m.omOp
	i0 := m.iOp
	m.spec.Inductance(th, m.l, m.d1, m.d2)
	g := m.spec.G
	pp, pn := m.p[0], m.p[1]
	for k := 0; k < m.n; k++ {
		r := m.br[k]
		ctx.M.Add(m.w[k][0], r, 1)
		ctx.M.Add(m.w[k][1], r, -1)
		ctx.M.Add(r, m.w[k][0], 1)
		ctx.M.Add(r, m.w[k][1], -1)
		var li0, gi0 float64
		for b := 0; b < m.n; b++ {
			ctx.M.Add(r, m.br[b], -jw*complex(m.l[k][b], 0)-complex(om*g[k][b], 0))
			li0 += m.d1[k][b] * i0[b]
			gi0 += g[k][b] * i0[b]
		}
		ctx.M.Add(r, m.th, -jw*complex(li0, 0))
		ctx.M.Add(r, pp, complex(-gi0, 0))
		ctx.M.Add(r, pn, complex(gi0, 0))
	}
	ctx.M.Add(m.th, m.th, jw)
	ctx.M.Add(m.th, pp, -1)
	ctx.M.Add(m.th, pn, 1)
	var tt float64
	for b := 0; b < m.n; b++ {
		var d float64
		for a := 0; a < m.n; a++ {
			d += (m.d1[b][a] + g[a][b] + g[b][a]) * i0[a]
			tt += i0[a] * m.d2[a][b] * i0[b]
		}
		ctx.M.Add(pp, m.br[b], complex(-d, 0))
		ctx.M.Add(pn, m.br[b], complex(d, 0))
	}
	ctx.M.Add(pp, m.th, complex(-0.5*tt, 0))
	ctx.M.Add(pn, m.th, complex(0.5*tt, 0))
	return nil
}

// OperatingAC returns winding 1's phasors.
func (m *Machine) OperatingAC(x []complex128) ACPoint {
	v := cat(x, m.w[0][0]) - cat(x, m.w[0][1])
	return ACPoint{Voltage: v, Current: x[m.br[0]]}
}

// ACPower is ½·Σ V·I* over the windings and the port (torque current −δT).
func (m *Machine) ACPower(x []complex128) complex128 {
	var s complex128
	for k := 0; k < m.n; k++ {
		s += (cat(x, m.w[k][0]) - cat(x, m.w[k][1])) * cmplx.Conj(x[m.br[k]])
	}
	// Port current: −δT, from the linearization at the operating point.
	m.spec.Inductance(m.thOp, nil, m.d1, m.d2)
	var dT complex128
	var tt float64
	g := m.spec.G
	for b := 0; b < m.n; b++ {
		var d float64
		for a := 0; a < m.n; a++ {
			d += (m.d1[b][a] + g[a][b] + g[b][a]) * m.iOp[a]
			tt += m.iOp[a] * m.d2[a][b] * m.iOp[b]
		}
		dT += complex(d, 0) * x[m.br[b]]
	}
	dT += complex(0.5*tt, 0) * x[m.th]
	s += (cat(x, m.p[0]) - cat(x, m.p[1])) * cmplx.Conj(-dT)
	return 0.5 * s
}

func cat(x []complex128, i int) complex128 {
	if i < 0 {
		return 0
	}
	return x[i]
}
