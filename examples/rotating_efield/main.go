// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command rotating_efield studies the "inverse homopolar" question
// (CLAUDE.md §3f, milestone M9n): can a rotating electric field make a
// magnetic field? It separates the two physical questions:
//
//  1. Real physics. A magnetic field comes from moving charge (Rowland's
//     rotating charged disk, 1876) and from a changing electric field (the
//     displacement current of the Ampère–Maxwell law). Four cases:
//     a. Rowland disk: B on the axis of a spinning charged disk.
//     b. A rotating uniform E field (as between quadrature-driven plate
//     pairs): its displacement current gives B = (μ₀ε₀/2)·∂E/∂t × r.
//     c. A charged ring at rest vs the same ring spinning: the E field is
//     identical and does not change in time, yet B is zero in one case
//     and not in the other — B follows the moving charge, not the
//     "rotation" of a symmetric E pattern (the mirror of the Faraday
//     disk paradox).
//     d. A rotating pair of charges: its E pattern really rotates, and
//     the B of the moving charges matches the Ampère–Maxwell law with
//     the displacement current of that rotating E (checked with the M9h
//     operators).
//
//  2. Hypothetical physics (opt-in): the force a magnetic charge would
//     feel in these fields — the dual pair of Joot's Theorem 3.16 —
//     labelled hypothetical and kept out of the real momentum ledger.
//
//     go run ./examples/rotating_efield
package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/hyper/quat"
)

// RowlandB is B_z on the axis of a disk of surface charge sigma and radius
// R spinning at w, at height z: by summing current rings (Biot–Savart)
// and analytically (μ₀σω/2)[(R² + 2z²)/√(R² + z²) − 2|z|].
func RowlandB(sigma, R, w, z float64) (numeric, analytic float64) {
	const rings = 400
	for k := 0; k < rings; k++ {
		r := R * (float64(k) + 0.5) / rings
		l := em.Loop{P: em.Circle([3]float64{}, [3]float64{0, 0, 1}, r, 720), I: sigma * w * r * R / rings}
		numeric += l.B([3]float64{0, 0, z})[2]
	}
	analytic = em.Mu0 * sigma * w / 2 * ((R*R+2*z*z)/math.Sqrt(R*R+z*z) - 2*math.Abs(z))
	return numeric, analytic
}

// RotatingUniformE returns the B of a uniform E field of amplitude e0
// rotating at w, at distance r from the axis, when the displacement
// current closes symmetrically: |B| = (μ₀ε₀/2)·e0·w·r. It also returns the
// Ampère–Maxwell and Faraday residuals of that field pair on a small grid,
// relative to their scales (the Faraday one measures the neglected
// (ωr/c)² correction).
func RotatingUniformE(e0, w, r float64) (b, ampere, faraday float64) {
	b = em.Mu0 * em.Eps0 / 2 * e0 * w * r
	g := em.Grid{N: [3]int{9, 9, 9}, H: [3]float64{r / 4, r / 4, r / 4}, Origin: [3]float64{-r, -r, -r}}
	field := func(t float64) (em.QField, em.QField) {
		s, c := math.Sincos(w * t)
		edot := [3]float64{-e0 * w * s, e0 * w * c, 0}
		E := em.SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(e0*c, e0*s, 0) })
		B := em.SampleQ(g, func(x, y, z float64) quat.Q {
			// (μ₀ε₀/2)·Ė × r
			k := em.Mu0 * em.Eps0 / 2
			return quat.Vec(k*(edot[1]*z-edot[2]*y), k*(edot[2]*x-edot[0]*z), k*(edot[0]*y-edot[1]*x))
		})
		return E, B
	}
	t0, dt := 0.2/w, 0.01/w
	E, B := field(t0)
	snap := func(k float64) (em.QField, em.QField) { return field(t0 + k*dt) }
	em2, bm2 := snap(-2)
	em1, bm1 := snap(-1)
	ep1, bp1 := snap(1)
	ep2, bp2 := snap(2)
	dE, _ := em.TimeDerivative(em2, em1, ep1, ep2, dt)
	dB, _ := em.TimeDerivative(bm2, bm1, bp1, bp2, dt)
	res := em.MaxwellQuaternion(E, B, dE, dB, nil, nil)
	for n := 0; n < g.Len(); n++ {
		if !g.Interior(n) {
			continue
		}
		ampere = math.Max(ampere, res.R2.Q[n].Abs())
		faraday = math.Max(faraday, res.R1.Q[n].Abs())
	}
	return b, ampere / (em.Mu0 * em.Eps0 * e0 * w), faraday / (e0 / r)
}

// RingCases returns B at the centre of a ring of charge q, radius a, at
// rest and spinning at w: 0 and μ₀I/(2a) with I = qω/2π, though E is the
// same static pattern in both.
func RingCases(q, a, w float64) (atRest, spinning float64) {
	l := em.Loop{P: em.Circle([3]float64{}, [3]float64{0, 0, 1}, a, 720), I: q * w / (2 * math.Pi)}
	return 0, l.B([3]float64{})[2]
}

// RotatingPair returns, for charges +q and −q on a rod of half-length a
// spinning at w, |B| on the axis at height a/2 from the moving charges
// (low-velocity Biot–Savart, μ₀q v×r̂/(4πr²); at the centre itself the two
// charges' contributions cancel, their currents circulating in opposite
// senses), and the relative
// Ampère–Maxwell residual of the fields on a grid around the centre, away
// from the charges, using the displacement current of the rotating E.
func RotatingPair(q, a, w float64) (bCentre, ampere float64) {
	pos := func(t float64) ([2][3]float64, [2][3]float64) {
		s, c := math.Sincos(w * t)
		p := [3]float64{a * c, a * s, 0}
		v := [3]float64{-a * w * s, a * w * c, 0}
		return [2][3]float64{p, {-p[0], -p[1], 0}}, [2][3]float64{v, {-v[0], -v[1], 0}}
	}
	charges := [2]float64{q, -q}
	fields := func(t float64, x [3]float64) ([3]float64, [3]float64) {
		ps, vs := pos(t)
		var e, b [3]float64
		for k := 0; k < 2; k++ {
			r := [3]float64{x[0] - ps[k][0], x[1] - ps[k][1], x[2] - ps[k][2]}
			r3 := math.Pow(r[0]*r[0]+r[1]*r[1]+r[2]*r[2], 1.5)
			ek := charges[k] / (4 * math.Pi * em.Eps0 * r3)
			bk := em.Mu0 * charges[k] / (4 * math.Pi * r3)
			v := vs[k]
			vxr := [3]float64{v[1]*r[2] - v[2]*r[1], v[2]*r[0] - v[0]*r[2], v[0]*r[1] - v[1]*r[0]}
			for i := 0; i < 3; i++ {
				e[i] += ek * r[i]
				b[i] += bk * vxr[i]
			}
		}
		return e, b
	}
	_, bc := fields(0, [3]float64{0, 0, a / 2})
	bCentre = math.Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2])
	h := a / 20
	g := em.Grid{N: [3]int{9, 9, 9}, H: [3]float64{h, h, h}, Origin: [3]float64{-4 * h, -4 * h, -4 * h}}
	at := func(t float64) (em.QField, em.QField) {
		E := em.SampleQ(g, func(x, y, z float64) quat.Q {
			e, _ := fields(t, [3]float64{x, y, z})
			return quat.Vec(e[0], e[1], e[2])
		})
		B := em.SampleQ(g, func(x, y, z float64) quat.Q {
			_, b := fields(t, [3]float64{x, y, z})
			return quat.Vec(b[0], b[1], b[2])
		})
		return E, B
	}
	dt := 1e-3 / w
	E, B := at(0)
	e2, b2 := at(-2 * dt)
	e1, b1 := at(-dt)
	e3, b3 := at(dt)
	e4, b4 := at(2 * dt)
	dE, _ := em.TimeDerivative(e2, e1, e3, e4, dt)
	dB, _ := em.TimeDerivative(b2, b1, b3, b4, dt)
	res := em.MaxwellQuaternion(E, B, dE, dB, nil, nil)
	scale := 0.0
	for n := 0; n < g.Len(); n++ {
		if g.Interior(n) {
			ampere = math.Max(ampere, quat.Vec(res.R2.Q[n].X, res.R2.Q[n].Y, res.R2.Q[n].Z).Abs())
			scale = math.Max(scale, dE.Q[n].Abs()/(em.SpeedOfLight*em.SpeedOfLight))
		}
	}
	return bCentre, ampere / scale
}

func main() {
	fmt.Println("1a. Rowland disk (σ = 1e−5 C/m², R = 10 cm, 20 000 rpm):")
	w := 20000 * 2 * math.Pi / 60
	for _, z := range []float64{0, 0.05} {
		n, a := RowlandB(1e-5, 0.1, w, z)
		fmt.Printf("    B_z(z = %.2f m) = %.4g T (rings) vs %.4g T (analytic)\n", z, n, a)
	}
	fmt.Println("1b. Uniform E of 1 MV/m rotating at 1 MHz:")
	b, amp, far := RotatingUniformE(1e6, 2*math.Pi*1e6, 0.05)
	fmt.Printf("    |B| at r = 5 cm: %.3g T; Ampère–Maxwell residual %.2g, Faraday residual %.2g (the neglected (ωr/c)² term)\n", b, amp, far)
	rest, spin := RingCases(1e-6, 0.05, 2*math.Pi*1000)
	fmt.Printf("1c. Ring of 1 µC, 5 cm: same static E pattern; B at centre %.3g T at rest, %.3g T spinning at 1000 rev/s\n", rest, spin)
	bc, ar := RotatingPair(1e-9, 0.05, 2*math.Pi*1000)
	fmt.Printf("1d. ±1 nC on a 10 cm rod at 1000 rev/s: |B| on the axis 2.5 cm up %.3g T; Ampère–Maxwell residual %.2g of the displacement term (stencil truncation)\n", bc, ar)

	fmt.Println("2. HYPOTHETICAL — a magnetic charge q_m = 1 A·m (none has ever been observed):")
	h := em.Hypothetical{MagneticCharge: true}
	_, bRow := RowlandB(1e-5, 0.1, w, 0)
	f, _, _ := em.DyonForce(h, 0, 1, [3]float64{}, [3]float64{}, [3]float64{0, 0, bRow})
	var l em.MomentumLedger
	l.AddExterior("monopole", "Rowland disk field", f, true)
	fmt.Printf("    at rest at the Rowland disk's centre: force %.3g N (hypothetical column; real net force %v)\n", f[2], l.Net())
	v := [3]float64{0, 100, 0}
	f, _, _ = em.DyonForce(h, 0, 1, v, [3]float64{1e6, 0, 0}, [3]float64{})
	fmt.Printf("    moving at 100 m/s through 1 MV/m: force %.3g N (q_m·v×E/c²)\n", math.Sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]))
}
