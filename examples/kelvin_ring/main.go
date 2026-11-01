// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kelvin_ring simulates Lord Kelvin's 1899 charged-ring molecule
// ("Magnetism and Molecular Rotation", Math. Phys. Papers Vol. V, Art. 175,
// pp. 354–358; reference/mathematicalphys05kelvuoft.pdf; notes in
// doc/kelvin-references.md) with standard electrodynamics only.
//
// The ring: a rigid non-conducting tube of radius a carrying charge +e
// (mass w′), filled with an incompressible fluid carrying −e (mass w) that
// can slide round the tube without friction.
//
//  1. Induction (Kelvin §§1–3). A coil's field is switched on. Every charge
//     feels the Lorentz force q(E + v×B) of the coil's changing field,
//     evaluated point by point with em.LorentzVector; the field comes from
//     Biot–Savart (em.Loop) and its flux is checked against Maxwell's
//     elliptic-integral mutual inductance (em.MutualCircles). Tube and fluid
//     turn in opposite senses with moments of momentum ∓eΦ/2π (Kelvin's
//     corrected eAM/2π). The energy and angular-momentum ledgers close;
//     the fluid turns at Larmor's rate eB/2w; switching the field off
//     returns everything to rest.
//
//  2. Zeeman test (Kelvin §§5–7). The ring is held by isotropic angular
//     springs (Kelvin's spring connexion to the ether) in a steady field and
//     made to vibrate; the spectrum of its magnetic moment is computed for
//     rings in every orientation, and compared with Lorentz's bound charge
//     (the standard classical model of the normal Zeeman effect). Three
//     models are compared:
//     - Lorentz's bound charge: m r̈ = −k r + q ṙ×B;
//     - Kelvin's gyrostat as he wrote it: fixed internal spin
//     S = (ea²/2) n·B along the ring axis n;
//     - the full ring: every charge under the Lorentz force, the fluid free
//     to slip, nothing linearised (rigid-body rotation with hyper/quat).
//
//     go run ./examples/kelvin_ring
package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand/v2"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/hyper/quat"
)

// ---------- small vector helpers ----------

type v3 = [3]float64

func add(a, b v3) v3           { return v3{a[0] + b[0], a[1] + b[1], a[2] + b[2]} }
func sub(a, b v3) v3           { return v3{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func scale(a v3, s float64) v3 { return v3{a[0] * s, a[1] * s, a[2] * s} }
func dot(a, b v3) float64      { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
func cross(a, b v3) v3 {
	return v3{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}
func norm(a v3) float64 { return math.Sqrt(dot(a, a)) }

func rot(r quat.Q, v v3) v3 {
	x, y, z := r.Rotate(v[0], v[1], v[2])
	return v3{x, y, z}
}

// solve3 solves m·x = b by Cramer's rule.
func solve3(m [3][3]float64, b v3) v3 {
	det := func(a [3][3]float64) float64 {
		return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
	}
	d := det(m)
	var x v3
	for c := 0; c < 3; c++ {
		mc := m
		for r := 0; r < 3; r++ {
			mc[r][c] = b[r]
		}
		x[c] = det(mc) / d
	}
	return x
}

// ---------- part 1: induction ----------

// Ring is Kelvin's ring: radius A, tube charge +E and mass MTube, fluid
// charge −E and mass MFluid.
type Ring struct{ A, E, MTube, MFluid float64 }

// Coil is a single-turn-equivalent coil of radius R carrying NI
// ampere-turns, coaxial with the ring at axial distance Z below it.
type Coil struct{ R, Z, NI float64 }

// Induction is the outcome of switching the coil on (hold) and off again.
type Induction struct {
	Flux, FluxMaxwell   float64 // flux through the ring at full current, Biot–Savart and elliptic [Wb]
	LTube, LFluid       float64 // moments of momentum after switch-on [kg m²/s]
	LKelvin             float64 // eΦ/2π
	OmegaFluid, Larmor  float64 // fluid angular velocity, eB̄/2w with B̄ = Φ/πa² [rad/s]
	Moment              float64 // magnetic moment of the ring currents [A m²]
	KE, WSource         float64 // kinetic energy gained, work by the coil's source [J]
	WInteraction        float64 // M·I_coil·I_ring stored in the mutual field [J]
	LedgerResidual      float64 // (W_source − KE − W_int)/KE
	AngularImpulseNet   float64 // net angular impulse on the ring (tube + fluid) [kg m²/s]
	AxialForce          float64 // force on the ring along the axis at full field [N] (+: away from coil)
	LTubeEnd, LFluidEnd float64 // after switch-off
	WSourceEnd          float64 // net work by the source over on + off [J]
}

const nCharges = 24

// Run switches the coil current on with a smooth ramp of duration T,
// holds it for T and ramps it off over T, integrating tube and fluid with
// RK4 under the Lorentz force on nCharges points each.
func (rg Ring) Run(c Coil, T float64, steps int) Induction {
	loop := em.Loop{P: em.Circle(v3{0, 0, -c.Z}, v3{0, 0, 1}, c.R, 720), I: 1}
	// Field per ampere-turn at the ring's points (axisymmetric, so one
	// azimuth suffices; rotate for the others).
	b1 := loop.B(v3{rg.A, 0, 0})
	// Flux per ampere-turn by Gauss–Legendre quadrature of B_z over the
	// ring's disk.
	var flux1 float64
	const nq = 64
	for k := 0; k < nq; k++ {
		r := rg.A * (float64(k) + 0.5) / nq
		flux1 += loop.B(v3{r, 0, 0})[2] * 2 * math.Pi * r * rg.A / nq
	}
	res := Induction{Flux: flux1 * c.NI, FluxMaxwell: em.MutualCircles(c.R, rg.A, c.Z) * c.NI}

	// Current profile: smooth ramp on [0,T], hold [T,2T], ramp off [2T,3T].
	s := func(x float64) float64 { return x - math.Sin(2*math.Pi*x)/(2*math.Pi) }
	ds := func(x float64) float64 { return 1 - math.Cos(2*math.Pi*x) }
	cur := func(t float64) (i, di float64) {
		switch {
		case t < T:
			return c.NI * s(t/T), c.NI * ds(t/T) / T
		case t < 2*T:
			return c.NI, 0
		case t < 3*T:
			return c.NI * (1 - s((t-2*T)/T)), -c.NI * ds((t-2*T)/T) / T
		}
		return 0, 0
	}
	// Torque and axial force on a rotor of charge q turning at ω.
	rotor := func(q, w, i, di float64) (tz, fz float64) {
		ephi := -di * flux1 / (2 * math.Pi * rg.A) // Faraday: ∮E·dl = −dΦ/dt
		for k := 0; k < nCharges; k++ {
			ph := 2 * math.Pi * float64(k) / nCharges
			cs, sn := math.Cos(ph), math.Sin(ph)
			p := v3{rg.A * cs, rg.A * sn, 0}
			v := v3{-w * rg.A * sn, w * rg.A * cs, 0}
			e := v3{-ephi * sn, ephi * cs, 0}
			b := v3{b1[0] * cs * i, b1[0] * sn * i, b1[2] * i}
			f, _ := em.LorentzVector(q/nCharges, v, e, b)
			tz += cross(p, f)[2]
			fz += f[2]
		}
		return tz, fz
	}
	it, iff := rg.MTube*rg.A*rg.A, rg.MFluid*rg.A*rg.A
	ringCurrent := func(wt, wf float64) float64 { return (rg.E*wt - rg.E*wf) / (2 * math.Pi) }
	deriv := func(t, wt, wf float64) (float64, float64) {
		i, di := cur(t)
		tt, _ := rotor(rg.E, wt, i, di)
		tf, _ := rotor(-rg.E, wf, i, di)
		return tt / it, tf / iff
	}
	var wt, wf, wsrc, impulse float64
	dt := T / float64(steps)
	for n := 0; n < 3*steps; n++ {
		t := n * 1
		_ = t
		t0 := float64(n) * dt
		k1t, k1f := deriv(t0, wt, wf)
		k2t, k2f := deriv(t0+dt/2, wt+dt/2*k1t, wf+dt/2*k1f)
		k3t, k3f := deriv(t0+dt/2, wt+dt/2*k2t, wf+dt/2*k2f)
		k4t, k4f := deriv(t0+dt, wt+dt*k3t, wf+dt*k3f)
		nwt := wt + dt/6*(k1t+2*k2t+2*k3t+k4t)
		nwf := wf + dt/6*(k1f+2*k2f+2*k3f+k4f)
		// Source work against the ring's back-EMF, M·I_c·dI_r/dt
		// (trapezoid), and the net angular impulse on tube + fluid.
		ia, _ := cur(t0)
		ib, _ := cur(t0 + dt)
		dIr := ringCurrent(nwt, nwf) - ringCurrent(wt, wf)
		wsrc += flux1 * (ia + ib) / 2 * dIr
		impulse += it*(nwt-wt) + iff*(nwf-wf)
		wt, wf = nwt, nwf
		if n == 2*steps-1 { // end of the hold
			res.LTube, res.LFluid = it*wt, iff*wf
			res.OmegaFluid = wf
			res.KE = 0.5*it*wt*wt + 0.5*iff*wf*wf
			res.WSource = wsrc
			ir := ringCurrent(wt, wf)
			res.WInteraction = flux1 * c.NI * ir
			res.LedgerResidual = (res.WSource - res.KE - res.WInteraction) / res.KE
			res.Moment = ir * math.Pi * rg.A * rg.A
			res.AngularImpulseNet = impulse
			_, ft := rotor(rg.E, wt, c.NI, 0)
			_, ff := rotor(-rg.E, wf, c.NI, 0)
			res.AxialForce = ft + ff
		}
	}
	res.LKelvin = rg.E * res.Flux / (2 * math.Pi)
	res.Larmor = rg.E * (res.Flux / (math.Pi * rg.A * rg.A)) / (2 * rg.MFluid)
	res.LTubeEnd, res.LFluidEnd = it*wt, iff*wf
	res.WSourceEnd = wsrc
	return res
}

// ---------- part 2: Zeeman test ----------

// Model is one of the three molecule models of the Zeeman test.
type Model int

const (
	Lorentz  Model = iota // bound charge on an isotropic spring
	Gyrostat              // Kelvin's gyrostat as written (fixed spin along n)
	FullRing              // Kelvin's ring under the full Lorentz force
)

func (m Model) String() string {
	return [...]string{"Lorentz bound charge", "Kelvin gyrostat (as written)", "Kelvin ring (full Lorentz force)"}[m]
}

// Zeeman holds the dimensionless set-up: a = w = w′ = e = 1, spring
// constants giving ω₀ = 1 for every model, field B along z with Larmor
// frequency ωL = eB/2w.
type Zeeman struct {
	OmegaL     float64 // eB/2w
	T, Dt      float64 // record length and step
	Stride     int     // samples recorded every Stride steps
	Kick       float64 // initial angular (or linear) velocity
	Orients    int     // orientations, uniform in cos θ on (0, 1)
	Seed       uint64
	FMin, FMax float64 // analysis band
	NFreq      int
}

// DefaultZeeman is the set-up used by the example and its test.
func DefaultZeeman() Zeeman {
	return Zeeman{OmegaL: 0.05, T: 2400, Dt: 0.04, Stride: 5, Kick: 0.01, Orients: 12, Seed: 1899,
		FMin: 0.9, FMax: 1.1, NFreq: 801}
}

// series integrates one molecule with axis n (cos θ = n_z) and returns,
// sampled every Stride·Dt, its mode coordinate (displacement r for Lorentz,
// rotation vector θ for the rings), which shows every mode frequency, and
// its radiating dipole (electric for Lorentz, the ring's own magnetic
// moment for the rings), which sets the polarisation an observer sees.
func (z Zeeman) series(m Model, n v3, kick v3) (mode, rad []v3) {
	B := v3{0, 0, 2 * z.OmegaL} // e = w = 1 ⇒ B = 2ωL
	steps := int(z.T / z.Dt)
	out := make([]v3, 0, steps/z.Stride+1)
	var modes []v3
	switch m {
	case Lorentz:
		// m r̈ = −k r + q ṙ×B, m = k = q = 1.
		r, v := v3{}, kick
		acc := func(r, v v3) v3 {
			f, _ := em.LorentzVector(1, v, v3{}, B)
			return add(scale(r, -1), f)
		}
		for s := 0; s < steps; s++ {
			if s%z.Stride == 0 {
				out = append(out, r)
				modes = append(modes, r)
			}
			a1 := acc(r, v)
			r2, v2 := add(r, scale(v, z.Dt/2)), add(v, scale(a1, z.Dt/2))
			a2 := acc(r2, v2)
			r3, v3_ := add(r, scale(v2, z.Dt/2)), add(v, scale(a2, z.Dt/2))
			a3 := acc(r3, v3_)
			r4, v4 := add(r, scale(v3_, z.Dt)), add(v, scale(a3, z.Dt))
			a4 := acc(r4, v4)
			r = add(r, scale(add(add(v, scale(add(v2, v3_), 2)), v4), z.Dt/6))
			v = add(v, scale(add(add(a1, scale(add(a2, a3), 2)), a4), z.Dt/6))
		}
	case Gyrostat:
		// I θ̈ = −κθ + S n×θ̇ (Kelvin §5, fixed spin S along n), I = κ = 1,
		// S = (ea²/2) n·B. Radiator: the fluid's moment −(e/2w)S·n carried
		// by the tilting axis, plus the tube's spin about n.
		S := 0.5 * dot(n, B)
		th, om := v3{}, kick
		acc := func(th, om v3) v3 { return add(scale(th, -1), scale(cross(n, om), S)) }
		mom := func(th, om v3) v3 {
			nt := add(n, cross(th, n))
			return add(scale(n, 0.5*dot(om, n)), scale(nt, -0.5*S))
		}
		for s := 0; s < steps; s++ {
			if s%z.Stride == 0 {
				out = append(out, mom(th, om))
				modes = append(modes, th)
			}
			a1 := acc(th, om)
			t2, o2 := add(th, scale(om, z.Dt/2)), add(om, scale(a1, z.Dt/2))
			a2 := acc(t2, o2)
			t3, o3 := add(th, scale(o2, z.Dt/2)), add(om, scale(a2, z.Dt/2))
			a3 := acc(t3, o3)
			t4, o4 := add(th, scale(o3, z.Dt)), add(om, scale(a3, z.Dt))
			a4 := acc(t4, o4)
			th = add(th, scale(add(add(om, scale(add(o2, o3), 2)), o4), z.Dt/6))
			om = add(om, scale(add(add(a1, scale(add(a2, a3), 2)), a4), z.Dt/6))
		}
	case FullRing:
		modes, out = z.fullRing(n, kick, B, steps)
	}
	return modes, out
}

// ringState is the full ring: orientation R (body z = ring axis), total
// angular momentum H (world), fluid spin u about the axis.
type ringState struct {
	R quat.Q
	H v3
	U float64
}

// fullRing integrates the tube (a = w′ = 1, inertia ½,½,1) and the fluid
// (w = 1, slipping freely along the tube) under the Lorentz force on
// eight charges each, held by isotropic angular springs κ = 1 about the
// starting orientation. Transverse inertia of tube + fluid = 1 and axial
// inertia of the tube = 1, so ω₀ = 1 for every mode. The fluid starts at
// the spin the switch-on gave it, u = e(n·B)/2w (part 1); the tube at rest.
func (z Zeeman) fullRing(n v3, kick v3, B v3, steps int) (modes, rad []v3) {
	const nq = 8
	// Orientation taking body z to n.
	zh := v3{0, 0, 1}
	ax := cross(zh, n)
	R0 := quat.Scalar(1)
	if s := norm(ax); s > 1e-12 {
		R0, _ = quat.Rotor(ax[0], ax[1], ax[2], math.Atan2(s, dot(zh, n)))
	}
	R0inv := R0.Conj()
	omegaOf := func(st ringState) (v3, v3) {
		nn := rot(st.R, zh)
		m := st.R.Matrix()
		// World inertia: tube R·diag(½,½,1)·Rᵀ + fluid transverse ½(1 − nnᵀ).
		var M [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				t := 0.5*m[i][0]*m[j][0] + 0.5*m[i][1]*m[j][1] + 1*m[i][2]*m[j][2]
				f := -0.5 * nn[i] * nn[j]
				if i == j {
					f += 0.5
				}
				M[i][j] = t + f
			}
		}
		return solve3(M, sub(st.H, scale(nn, st.U))), nn
	}
	type deriv struct {
		dR quat.Q
		dH v3
		dU float64
	}
	f := func(st ringState) (deriv, v3, v3) {
		w, nn := omegaOf(st)
		wf := add(w, scale(nn, st.U-dot(w, nn))) // fluid angular velocity
		var tauT, tauF, mu v3
		for k := 0; k < nq; k++ {
			ph := 2 * math.Pi * float64(k) / nq
			p := rot(st.R, v3{math.Cos(ph), math.Sin(ph), 0})
			vt, vf := cross(w, p), cross(wf, p)
			ft, _ := em.LorentzVector(1.0/nq, vt, v3{}, B)
			ff, _ := em.LorentzVector(-1.0/nq, vf, v3{}, B)
			tauT = add(tauT, cross(p, ft))
			tauF = add(tauF, cross(p, ff))
			mu = add(mu, add(scale(cross(p, vt), 0.5/nq), scale(cross(p, vf), -0.5/nq)))
		}
		// Spring torque −κθ, θ the rotation vector of R·R0⁻¹.
		lg, _ := st.R.Mul(R0inv).Log()
		th := v3{2 * lg.X, 2 * lg.Y, 2 * lg.Z}
		return deriv{
			dR: quat.Vec(w[0], w[1], w[2]).Mul(st.R).Scale(0.5),
			dH: add(scale(th, -1), add(tauT, tauF)),
			dU: dot(nn, tauF), // fluid axial inertia w·a² = 1
		}, mu, th
	}
	step := func(st ringState, d deriv, h float64) ringState {
		return ringState{R: st.R.Add(d.dR.Scale(h)), H: add(st.H, scale(d.dH, h)), U: st.U + h*d.dU}
	}
	u0 := dot(n, B) / 2
	st := ringState{R: R0, U: u0}
	// Initial angular velocity = kick: H = I·ω + u n, total inertia
	// isotropic (= 1).
	st.H = add(kick, scale(n, u0))
	for s := 0; s < steps; s++ {
		d1, mu, th := f(st)
		if s%z.Stride == 0 {
			rad = append(rad, mu)
			modes = append(modes, th)
		}
		d2, _, _ := f(step(st, d1, z.Dt/2))
		d3, _, _ := f(step(st, d2, z.Dt/2))
		d4, _, _ := f(step(st, d3, z.Dt))
		st = ringState{
			R: st.R.Add(d1.dR.Add(d2.dR.Scale(2)).Add(d3.dR.Scale(2)).Add(d4.dR).Scale(z.Dt / 6)),
			H: add(st.H, scale(add(add(d1.dH, scale(add(d2.dH, d3.dH), 2)), d4.dH), z.Dt/6)),
			U: st.U + z.Dt/6*(d1.dU+2*d2.dU+2*d3.dU+d4.dU),
		}
		st.R, _ = st.R.Unit()
	}
	return modes, rad
}

// Spectrum is an ensemble spectrum: power viewed along B in the two
// circular senses (counter-clockwise and clockwise about +B), and viewed
// across B (along x) polarised parallel and perpendicular to B.
type Spectrum struct {
	Freq               []float64
	CCW, CW, Par, Perp []float64
}

func (s Spectrum) total(i int) float64 { return s.CCW[i] + s.CW[i] + s.Par[i] }

// dft accumulates the Hann-windowed spectrum of x into sp.
func (z Zeeman) dft(x []v3, sp *Spectrum) {
	n := len(x)
	dt := z.Dt * float64(z.Stride)
	var mean v3
	for _, v := range x {
		mean = add(mean, v)
	}
	mean = scale(mean, 1/float64(n))
	for k, w := range sp.Freq {
		var cx, cy, cz complex128
		rot := cmplx.Exp(complex(0, -w*dt))
		ph := complex(1, 0)
		for j, v := range x {
			hann := 0.5 - 0.5*math.Cos(2*math.Pi*float64(j)/float64(n-1))
			d := sub(v, mean)
			cx += complex(hann*d[0], 0) * ph
			cy += complex(hann*d[1], 0) * ph
			cz += complex(hann*d[2], 0) * ph
			ph *= rot
			if j%256 == 255 { // renormalise the rotating phasor
				ph /= complex(cmplx.Abs(ph), 0)
			}
		}
		// e^{+iωt} motion (counter-clockwise about +z) shows in x + iy.
		ccw := (cx + 1i*cy) / math.Sqrt2
		cw := (cx - 1i*cy) / math.Sqrt2
		// Radiated power ∝ ω⁴|dipole|²; ω ≈ 1 across the band, kept anyway.
		w4 := w * w * w * w
		sp.CCW[k] += w4 * real(ccw*cmplx.Conj(ccw))
		sp.CW[k] += w4 * real(cw*cmplx.Conj(cw))
		sp.Par[k] += w4 * real(cz*cmplx.Conj(cz))
		sp.Perp[k] += w4 * real(cy*cmplx.Conj(cy))
	}
}

// Line is a spectral peak.
type Line struct{ Freq, Power float64 }

// peaks returns local maxima of p above frac of the largest, with
// parabolic interpolation.
func peaks(freq, p []float64, frac float64) []Line {
	mx := 0.0
	for _, v := range p {
		mx = math.Max(mx, v)
	}
	var out []Line
	for i := 1; i+1 < len(p); i++ {
		if p[i] > p[i-1] && p[i] >= p[i+1] && p[i] > frac*mx {
			d := 0.5 * (p[i-1] - p[i+1]) / (p[i-1] - 2*p[i] + p[i+1])
			out = append(out, Line{freq[i] + d*(freq[1]-freq[0]), p[i]})
		}
	}
	return out
}

// GyroVector is the gyroscopic vector of the full ring held at axis n in
// a uniform field B (a = 1 units scaled by a²): the linearised equation of
// the tilt θ is I θ̈ = −κθ + θ̇ × G with
//
//	G = (e a²/2) (B − 2 (n·B) n),
//
// the reflection of B in the ring's plane, times e a²/2. Derivation: the
// fluid's canonical angular momentum about n, w a² u − e Φ(n)/2π, is
// conserved (zero from switch-on), so its spin follows the flux through
// the tilted ring, u = e(n·B)/2w; eliminating it (Routh) leaves the term
// (e a²/2)(n·B)(ω·n) linear in ω, whose curl in the rotation vector — with
// ω = θ̇ + ½ θ×θ̇ — is G. Kelvin's gyrostat keeps the spin fixed at its
// switch-on value and so has G = −(e a²/2)(n·B) n instead: the term
// B − 3(n·B)n from the flux following the tilt is what it leaves out.
// |G| = e a² |B|/2 for every orientation, so the splitting
// |G|/2I = e B/(2(w + w′)) does not depend on how the ring lies.
func GyroVector(e, a float64, n, B v3) v3 {
	return scale(sub(B, scale(n, 2*dot(n, B))), e*a*a/2)
}

// MotionPolarisationFull is the ensemble circular polarisation of the full
// ring's side lines seen along B, predicted from G: each mode circles Ĝ,
// whose component along B̂ is g = 1 − 2cos²θ; a circle about Ĝ seen along
// B̂ has intensity (1+g²)/2 and Stokes V = g, so over cos θ uniform on
// (0, 1) the ensemble value is ⟨g⟩/⟨(1+g²)/2⟩ = (1/3)/(11/15) = 5/11.
const MotionPolarisationFull = 5.0 / 11

// Orientation is the outcome for one ring orientation.
type Orientation struct {
	CosTheta float64
	Lines    []Line  // mode frequencies (peaks of the mode-coordinate spectrum)
	Split    float64 // largest distance of a line from the unshifted centre
}

// Polarisation of a spectrum's side bands and centre.
type Polarisation struct {
	UpperV, LowerV float64 // (CCW − CW)/(CCW + CW) viewed along B
	CentrePar      float64 // fraction of the centre polarised ∥ B, viewed across B
	SidePerp       float64 // fraction of the side bands polarised ⊥ B, viewed across B
}

// Result of the Zeeman test for one model.
type Result struct {
	Model              Model
	Mode, Rad          Spectrum // ensembles of the mode coordinate and the radiating dipole
	PerOrient          []Orientation
	Resolution         float64 // 2π/T
	SplitMin, SplitMax float64 // over orientations
	UpperWidth         float64 // FWHM of the upper side feature of the ensemble mode spectrum
	SharpFrac          float64 // side power within ±2·resolution of the main side pair
	ModePol, RadPol    Polarisation
}

// Run computes one model over the ensemble of orientations.
func (z Zeeman) Run(m Model) Result {
	rng := rand.New(rand.NewPCG(z.Seed, uint64(m)))
	freq := make([]float64, z.NFreq)
	for i := range freq {
		freq[i] = z.FMin + (z.FMax-z.FMin)*float64(i)/float64(z.NFreq-1)
	}
	newSp := func() Spectrum {
		return Spectrum{Freq: freq, CCW: make([]float64, z.NFreq), CW: make([]float64, z.NFreq),
			Par: make([]float64, z.NFreq), Perp: make([]float64, z.NFreq)}
	}
	accum := func(dst *Spectrum, src Spectrum) {
		for i := range freq {
			dst.CCW[i] += src.CCW[i]
			dst.CW[i] += src.CW[i]
			dst.Par[i] += src.Par[i]
			dst.Perp[i] += src.Perp[i]
		}
	}
	res := Result{Model: m, Mode: newSp(), Rad: newSp(), Resolution: 2 * math.Pi / z.T, SplitMin: math.Inf(1)}
	for k := 0; k < z.Orients; k++ {
		c := (float64(k) + 0.5) / float64(z.Orients)
		phi := 2 * math.Pi * rng.Float64()
		s := math.Sqrt(1 - c*c)
		n := v3{s * math.Cos(phi), s * math.Sin(phi), c}
		kick := v3{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		kick = scale(kick, z.Kick/norm(kick))
		mode, rad := z.series(m, n, kick)
		ms, rs := newSp(), newSp()
		z.dft(mode, &ms)
		z.dft(rad, &rs)
		accum(&res.Mode, ms)
		accum(&res.Rad, rs)
		tot := make([]float64, z.NFreq)
		for i := range tot {
			tot[i] = ms.total(i)
		}
		o := Orientation{CosTheta: c, Lines: peaks(freq, tot, 0.02)}
		for _, l := range o.Lines { // largest distance of a line from the centre
			o.Split = math.Max(o.Split, math.Abs(l.Freq-centre))
		}
		res.SplitMin = math.Min(res.SplitMin, o.Split)
		res.SplitMax = math.Max(res.SplitMax, o.Split)
		res.PerOrient = append(res.PerOrient, o)
	}
	res.UpperWidth, res.SharpFrac = lineShape(res.Mode, res.Resolution)
	res.ModePol = polarisation(res.Mode, res.Resolution)
	res.RadPol = polarisation(res.Rad, res.Resolution)
	return res
}

// centre of the band: the unshifted line (√(1+δ²) ≈ 1 to second order).
const centre = 1.0

// lineShape measures the upper side feature of an ensemble spectrum: its
// full width at half maximum, and the fraction of all side-band power
// (|ω − 1| > 2·resolution) lying within ±2·resolution of the strongest
// upper line or its mirror below.
func lineShape(sp Spectrum, res float64) (width, sharp float64) {
	n := len(sp.Freq)
	guard := 2 * res
	upMax, up := 0.0, -1
	var side float64
	for i, w := range sp.Freq {
		if math.Abs(w-centre) <= guard {
			continue
		}
		side += sp.total(i)
		if w > centre && sp.total(i) > upMax {
			upMax, up = sp.total(i), i
		}
	}
	if up < 0 || side == 0 {
		return 0, 0
	}
	lo, hi := up, up
	for lo > 0 && sp.total(lo-1) > upMax/2 && sp.Freq[lo-1] > centre+guard {
		lo--
	}
	for hi < n-1 && sp.total(hi+1) > upMax/2 {
		hi++
	}
	width = sp.Freq[hi] - sp.Freq[lo] + (sp.Freq[1] - sp.Freq[0])
	wu := sp.Freq[up]
	wl := 2*centre - wu
	var near float64
	for i, w := range sp.Freq {
		if math.Abs(w-centre) <= guard {
			continue
		}
		if math.Abs(w-wu) <= guard || math.Abs(w-wl) <= guard {
			near += sp.total(i)
		}
	}
	return width, near / side
}

func polarisation(sp Spectrum, res float64) Polarisation {
	guard := 2 * res
	var upCCW, upCW, loCCW, loCW, sPerp, sAcross, cPar, cAcross float64
	for i, w := range sp.Freq {
		across := sp.Par[i] + sp.Perp[i]
		switch {
		case w > centre+guard:
			upCCW += sp.CCW[i]
			upCW += sp.CW[i]
			sPerp += sp.Perp[i]
			sAcross += across
		case w < centre-guard:
			loCCW += sp.CCW[i]
			loCW += sp.CW[i]
			sPerp += sp.Perp[i]
			sAcross += across
		default:
			cPar += sp.Par[i]
			cAcross += across
		}
	}
	return Polarisation{
		UpperV: (upCCW - upCW) / (upCCW + upCW), LowerV: (loCCW - loCW) / (loCCW + loCW),
		CentrePar: cPar / cAcross, SidePerp: sPerp / sAcross,
	}
}

func main() {
	fmt.Println("Kelvin 1899, \"Magnetism and Molecular Rotation\" (Math. Phys. Papers V, Art. 175)")
	fmt.Println()
	fmt.Println("1. Induction: ring a = 5 cm, ±1 µC, tube and fluid 10 g each; coil r = 10 cm, 5 cm below, 20 kA-turns, 1 s ramps")
	rg := Ring{A: 0.05, E: 1e-6, MTube: 0.01, MFluid: 0.01}
	r := rg.Run(Coil{R: 0.1, Z: 0.05, NI: 2e4}, 1, 4000)
	fmt.Printf("   flux through ring %.6g Wb (Biot–Savart) vs %.6g Wb (Maxwell, elliptic)\n", r.Flux, r.FluxMaxwell)
	fmt.Printf("   moment of momentum: tube %.6g, fluid %.6g kg m²/s; Kelvin's corrected eAM/2π = %.6g\n", r.LTube, r.LFluid, r.LKelvin)
	fmt.Printf("   net angular impulse on the ring %.2g (the coil gets none back: the two halves cancel)\n", r.AngularImpulseNet)
	fmt.Printf("   fluid turns at %.6g rad/s = Larmor eB̄/2w %.6g rad/s (B̄ = Φ/πa²)\n", r.OmegaFluid, r.Larmor)
	fmt.Printf("   ring moment %.3g A m² (diamagnetic: opposes the coil field), axial force %.3g N away from the coil\n", r.Moment, r.AxialForce)
	fmt.Printf("   ledger: source work %.6g J = kinetic %.6g J + mutual-field energy %.6g J (residual %.1e)\n", r.WSource, r.KE, r.WInteraction, r.LedgerResidual)
	fmt.Printf("   after switch-off: tube %.2g, fluid %.2g kg m²/s; net source work %.2g J\n", r.LTubeEnd, r.LFluidEnd, r.WSourceEnd)
	fmt.Printf("   the same ring at atomic scale (electron charge and mass, 1 T): fluid at %.4g rad/s = Larmor's eB/2m\n", 1.602176634e-19/(2*9.1093837e-31))
	fmt.Println()

	z := DefaultZeeman()
	fmt.Printf("2. Zeeman test: ω₀ = 1, Larmor ωL = %.3g, %d orientations, resolution 2π/T = %.4f\n", z.OmegaL, z.Orients, 2*math.Pi/z.T)
	for _, m := range []Model{Lorentz, Gyrostat, FullRing} {
		res := z.Run(m)
		fmt.Printf("   %s\n      splitting by cos θ:", m)
		for _, o := range res.PerOrient {
			fmt.Printf(" %.2f→%.4f", o.CosTheta, o.Split)
		}
		fmt.Printf("\n      ensemble: splitting %.4f…%.4f, upper side FWHM %.4f, %.0f %% of side power in the main pair\n",
			res.SplitMin, res.SplitMax, res.UpperWidth, 100*res.SharpFrac)
		for _, p := range []struct {
			name string
			pol  Polarisation
		}{{"motion", res.ModePol}, {"radiation", res.RadPol}} {
			fmt.Printf("      %-9s along B: circular upper %+.2f, lower %+.2f; across B: centre %.0f %% ∥B, sides %.0f %% ⊥B\n",
				p.name, p.pol.UpperV, p.pol.LowerV, 100*p.pol.CentrePar, 100*p.pol.SidePerp)
		}
	}
}

// Resolution is the frequency resolution 2π/T of the record.
func (z Zeeman) Resolution() float64 { return 2 * math.Pi / z.T }
