// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package spacetime

import (
	"fmt"
	"math"
)

// Rotation and curvature together: a spinning machine (rotor, flywheel,
// homopolar disk, ring gyroscope) sitting in, or moving through, the
// gravitational field of a spinning body.
//
// How the two are combined (and why it is consistent):
//
//   - Clocks. The machine's hub follows a worldline in the body's
//     spacetime; its proper-time rate dτ_hub/dt comes from the metric
//     (gravity.go: Schwarzschild, plus the hub's own motion). A point on
//     the rim moves at u = ωρ relative to the hub's local inertial frame,
//     so over a turn dτ_rim/dτ_hub = √(1 − u²/c²) exactly (special
//     relativity in the local frame — the equivalence principle). The
//     product is the rim's rate seen from far away. To first order this is
//     the weak-field metric's 1 − U/c² − |v_hub + u|²/2c² averaged over a
//     turn, the twin checked in the tests.
//   - Gyroscopic directions. A spinning rotor's axis is a gyroscope; in
//     curved spacetime it turns slowly relative to the distant stars:
//     Ω = (1/2c²)(a_ng × v) + (3GM/2c²r³)(r × v) + (G/c²r³)[3(J·r̂)r̂ − J]
//     — Thomas precession from the non-gravitational acceleration that
//     holds it (support, thrust), geodetic (de Sitter) precession from
//     moving through the curved field, and frame dragging (Lense–Thirring)
//     from the body's spin J (L. I. Schiff, Phys. Rev. Lett. 4, 215, 1960;
//     J. D. Jackson, Classical Electrodynamics, 3rd ed., eq. 11.119 for
//     the Thomas term). Validated against Gravity Probe B.
//   - Light on the machine. A loop of fibre or a ring laser on the rotor
//     measures its rotation relative to the local gyroscopes (the local
//     inertial frame), so the Sagnac delay is 4A·n̂·(ω_vs_stars − Ω)/c²:
//     the machine's spin plus the body's rotation of the platform, minus
//     the precession of the local frame (the GINGER ring-laser idea; A. Di
//     Virgilio et al., C. R. Physique 15, 866, 2014, states the rotation
//     a ring laser on Earth measures as Earth's minus the geodetic and
//     Lense–Thirring terms).
//   - Fields on the machine. A probe riding the rim sees the fields of the
//     local co-moving inertial frame (em/frames.go boost at ω×ρ; Mashhoon's
//     locality hypothesis), in a medium of index n(r) for the gravity
//     (Plebanski 1960, optics.go). The circuit itself is solved in its own
//     proper time with ordinary Maxwell; frequencies are converted to the
//     distant observer at the end by the clock rate.
//   - The windings. Kron's curved configuration space (kron/machine.go)
//     stays as it is: it is the machine's own geometry, not spacetime.

// Vec is a 3-vector [SI].
type Vec [3]float64

func (a Vec) add(b Vec) Vec       { return Vec{a[0] + b[0], a[1] + b[1], a[2] + b[2]} }
func (a Vec) sub(b Vec) Vec       { return Vec{a[0] - b[0], a[1] - b[1], a[2] - b[2]} }
func (a Vec) scale(s float64) Vec { return Vec{a[0] * s, a[1] * s, a[2] * s} }
func (a Vec) dot(b Vec) float64   { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
func (a Vec) norm() float64       { return math.Sqrt(a.dot(a)) }

// Norm is the length of a.
func (a Vec) Norm() float64 { return a.norm() }

func (a Vec) cross(b Vec) Vec {
	return Vec{a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]}
}

// Body is a spinning mass at the origin: GM [m³/s²] and spin angular
// momentum J [kg·m²/s].
type Body struct {
	Name string
	GM   float64
	J    Vec
	R    float64 // surface radius [m], for checks
}

// Earth with J = (I/MR²)·M·R_mean²·ω along +z: I/MR² = 0.3308, M =
// 5.9722×10²⁴ kg, R_mean = 6371.000 km (NASA Earth Fact Sheet), ω from
// IERS 2010.
var Earth = Body{
	Name: "Earth",
	GM:   GMEarth,
	J:    Vec{0, 0, 0.3308 * 5.9722e24 * 6371000.0 * 6371000.0 * OmegaEarth},
	R:    REarthEq,
}

// Precession returns the three precession rates [rad/s] of a gyroscope
// (a rotor's spin axis) at position r with velocity v relative to the
// body's rest frame, held by a non-gravitational acceleration aNG (zero in
// free fall).
func Precession(b Body, r, v, aNG Vec) (thomas, geodetic, frameDrag Vec) {
	rr := r.norm()
	c2 := C * C
	thomas = aNG.cross(v).scale(1 / (2 * c2))
	geodetic = r.cross(v).scale(1.5 * b.GM / (c2 * rr * rr * rr))
	n := r.scale(1 / rr)
	frameDrag = n.scale(3 * b.J.dot(n)).sub(b.J).scale(G / (c2 * rr * rr * rr))
	return
}

// OrbitAveragedPrecession averages the geodetic and frame-dragging rates
// over a circular orbit of radius a whose plane contains the unit vectors
// e1, e2 (free fall: no Thomas term), over the orbit phase
// (midpoint rule, exact for these trigonometric integrands once
// steps ≥ 3).
func OrbitAveragedPrecession(b Body, a float64, e1, e2 Vec, steps int) (geodetic, frameDrag Vec) {
	v := math.Sqrt(b.GM / a)
	for k := 0; k < steps; k++ {
		th := 2 * math.Pi * (float64(k) + 0.5) / float64(steps)
		r := e1.scale(a * math.Cos(th)).add(e2.scale(a * math.Sin(th)))
		vv := e1.scale(-v * math.Sin(th)).add(e2.scale(v * math.Cos(th)))
		_, g, f := Precession(b, r, vv, Vec{})
		geodetic = geodetic.add(g.scale(1 / float64(steps)))
		frameDrag = frameDrag.add(f.scale(1 / float64(steps)))
	}
	return
}

// MasPerYear converts rad/s to milliarcseconds per year (Julian year).
func MasPerYear(w float64) float64 { return w * 365.25 * 86400 * 180 / math.Pi * 3600 * 1000 }

// RimRateInField is dτ/dt (distant-observer time) of a point at radius ρ
// on a rotor turning at ω, whose hub sits at radius r from the body and
// moves tangentially at vHub (0 if hovering/static; for a body-fixed lab,
// the body's surface speed): the hub's Schwarzschild rate times the
// special-relativistic rim factor √(1 − (ωρ)²/c²), averaged over a turn.
func RimRateInField(b Body, r, vHub, omega, rho float64) (float64, error) {
	hub, err := MovingRate(b.GM, r, vHub)
	if err != nil {
		return 0, err
	}
	rim, err := RimRate(omega, rho)
	if err != nil {
		return 0, err
	}
	return hub * rim, nil
}

// RimRateWeak is the same rate from the weak-field metric to first order,
// 1 − GM/(rc²) − ⟨|v_hub + u|²⟩/(2c²) with |u| = ωρ, averaged over a turn
// (the cross term averages to zero): the twin of RimRateInField.
func RimRateWeak(b Body, r, vHub, omega, rho float64) float64 {
	u := omega * rho
	return 1 - b.GM/(r*C*C) - (vHub*vHub+u*u)/(2*C*C)
}

// Site is a machine's place: position and velocity relative to the body,
// and the non-gravitational acceleration holding it there (a lab floor
// pushes up; free fall is zero).
type Site struct {
	Name   string
	R, V   Vec
	ANG    Vec
	SpinAx Vec // unit normal of the machine's rotor plane (its spin axis)
	// Turn is the rotation of the lab platform itself relative to the
	// distant stars [rad/s]: Earth's spin for a lab on the ground, zero
	// for an inertially pointed spacecraft, the orbital rate r×v/r² for
	// one that keeps facing the planet.
	Turn Vec
}

// EarthSurface returns a lab at rest on the rotating Earth at latitude
// lat [rad] on the meridian x–z, with its rotor axis along the local
// vertical (spherical Earth of equatorial radius; the floor supplies
// g plus the centripetal term).
func EarthSurface(lat float64) Site {
	r := Vec{REarthEq * math.Cos(lat), 0, REarthEq * math.Sin(lat)}
	w := Vec{0, 0, OmegaEarth}
	v := w.cross(r)
	up := r.scale(1 / r.norm())
	g := up.scale(GMEarth / (REarthEq * REarthEq)) // the floor's push, cancelling gravity
	cent := w.cross(v)                             // centripetal acceleration of the floor
	return Site{Name: fmt.Sprintf("Earth surface, latitude %.1f°", lat*180/math.Pi), R: r, V: v, ANG: g.add(cent), SpinAx: up, Turn: w}
}

// Machine is a rotor: radius, angular speed and material.
type Machine struct {
	Name     string
	Radius   float64 // m
	Omega    float64 // rad/s
	Material Material
	Loop     float64 // area of a fibre loop / ring laser riding the rotor [m²] (0: none)
}

// MachineReport gathers what rotation and curvature do to a machine.
type MachineReport struct {
	Site, Machine string
	// Clocks.
	HubRate, RimRate float64 // dτ/dt relative to a distant observer
	RimLagPerDay     float64 // s/day the rim clock loses against the hub
	GravShift        float64 // fractional shift of the machine's signals received far away
	// Rotor and field.
	RimSpeed, VOverC2, PhiOverC2 float64
	BurstMargin                  float64 // burst speed / rim speed
	TidalAcc, TidalStress        float64 // tidal acceleration across the rotor [m/s²], hoop-stress equivalent ρ·a·R [Pa]
	TidalMargin                  float64 // material strength / tidal stress
	// Gyroscopic axis (rad/s) and Sagnac.
	Thomas, Geodetic, FrameDrag Vec
	SagnacPlatform              float64 // delay from the platform's rotation (body and machine) [s]
	SagnacFrameCorrection       float64 // change from the local frame's precession [s]
	WeakField                   bool    // U/c² < 1e-3: first-order formulas accurate
	Warnings                    []string
}

// AnalyzeMachine combines the rotor, its site and the body.
func AnalyzeMachine(b Body, s Site, m Machine) (MachineReport, error) {
	rep := MachineReport{Site: s.Name, Machine: m.Name}
	r := s.R.norm()
	vHub := s.V.norm()
	var err error
	if rep.HubRate, err = MovingRate(b.GM, r, vHub); err != nil {
		return rep, err
	}
	if rep.RimRate, err = RimRateInField(b, r, vHub, m.Omega, m.Radius); err != nil {
		return rep, err
	}
	rep.RimLagPerDay = (1 - rep.RimRate/rep.HubRate) * 86400
	rep.GravShift = rep.HubRate - 1
	rep.RimSpeed = m.Omega * m.Radius
	rep.VOverC2 = rep.RimSpeed * rep.RimSpeed / (C * C)
	rep.PhiOverC2 = b.GM / (r * C * C)
	rep.WeakField = rep.PhiOverC2 < 1e-3
	if !rep.WeakField {
		rep.Warnings = append(rep.Warnings, fmt.Sprintf("strong field (GM/rc² = %.3g): clock rates use the exact Schwarzschild form; precession and Sagnac terms are first-order and only indicative", rep.PhiOverC2))
	}
	if m.Material.Density > 0 {
		rep.BurstMargin = BurstSpeed(m.Material) / rep.RimSpeed
		if rep.BurstMargin < 1 {
			rep.Warnings = append(rep.Warnings, fmt.Sprintf("rim speed %.0f m/s exceeds the burst speed %.0f m/s", rep.RimSpeed, BurstSpeed(m.Material)))
		}
		rep.TidalAcc = Tidal(b.GM, r, 2*m.Radius)
		rep.TidalStress = m.Material.Density * rep.TidalAcc * m.Radius
		rep.TidalMargin = m.Material.Strength / rep.TidalStress
		if rep.TidalMargin < 1 {
			rep.Warnings = append(rep.Warnings, fmt.Sprintf("tidal stress %.3g Pa exceeds the material strength", rep.TidalStress))
		}
	}
	rep.Thomas, rep.Geodetic, rep.FrameDrag = Precession(b, s.R, s.V, s.ANG)
	if m.Loop > 0 {
		// Platform rotation relative to the stars: the machine's spin about
		// its axis plus the lab platform's own turning.
		wPlat := s.SpinAx.scale(m.Omega).add(s.Turn)
		prec := rep.Thomas.add(rep.Geodetic).add(rep.FrameDrag)
		rep.SagnacPlatform = SagnacDelay(m.Loop, wPlat.dot(s.SpinAx))
		rep.SagnacFrameCorrection = -SagnacDelay(m.Loop, prec.dot(s.SpinAx))
	}
	return rep, nil
}

// SelfField is the machine's own gravity: what its mass and spin do to
// spacetime at distance d on its spin axis, from the far-field
// (weak-field, multipole) metric of a spinning mass — the same
// Lense–Thirring term as Precession. For d comparable to the rotor
// radius it is an order-of-magnitude estimate.
type SelfField struct {
	Mass, J       float64 // kg, kg·m²/s (uniform disk: ½MR²ω)
	PhiOverC2     float64 // GM/(dc²): fractional clock shift its mass causes at d
	FrameDragAxis float64 // 2GJ/(c²d³) [rad/s]: how fast it drags a gyroscope on its axis
	Curvature     float64 // tidal curvature GM/(c²d³) [1/m²], the Riemann component a ruler feels
}

// MachineSelfField returns SelfField for a uniform disk rotor of the
// machine's radius and speed with mass mass, at distance d.
func MachineSelfField(m Machine, mass, d float64) SelfField {
	j := 0.5 * mass * m.Radius * m.Radius * m.Omega
	return SelfField{Mass: mass, J: j,
		PhiOverC2:     G * mass / (d * C * C),
		FrameDragAxis: 2 * G * j / (C * C * d * d * d),
		Curvature:     G * mass / (C * C * d * d * d)}
}
