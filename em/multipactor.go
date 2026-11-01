// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import "math"

// Multipactor: two-surface resonant secondary-electron avalanche between
// parallel walls in vacuum under an RF field (Hatch & Williams 1954;
// Vaughan, IEEE Trans. Electron Devices 35, 1172, 1988). Electrons are
// tracked from either wall at every RF phase with 2 eV emission energy
// through the combined RF + DC field, impact after impact, multiplying the
// secondary yield at each impact. Yield: Furman's fit δ = δmax·s·x/(s − 1 +
// x^s), x = E/Emax. Space charge and the spread of emission energies are
// left out, so the check marks where a discharge can start, not how large
// it grows. (Moved here from examples/double_hull.)

const (
	qe = 1.602176634e-19  // C (CODATA 2018, exact)
	me = 9.1093837015e-31 // kg (CODATA 2018)
)

// SEY is a secondary-emission yield curve (Furman fit).
type SEY struct {
	Name              string
	DeltaMax, Emax, S float64 // —, eV, —
}

// Measured copper (Baglin et al., LHC Project Report 472, 2002, Table 1).
var (
	CopperAsReceived  = SEY{"copper, as received", 2.03, 262, 1.39}
	CopperConditioned = SEY{"copper, conditioned", 1.13, 318, 1.35}
)

// Yield at impact energy e [eV].
func (s SEY) Yield(e float64) float64 {
	if e <= 0 {
		return 0
	}
	x := e / s.Emax
	return s.DeltaMax * s.S * x / (s.S - 1 + math.Pow(x, s.S))
}

// Multipactor reports the largest mean secondary yield per impact over all
// launch phases, for a gap d, frequency f, RF amplitude vrf and bias vdc:
// above 1 a discharge grows. It returns that yield and the mean impact
// energy on the worst phase.
func Multipactor(d, f, vrf, vdc float64, sey SEY) (growth, eImpact float64) {
	w := 2 * math.Pi * f
	k := qe / (me * d) // acceleration per volt across the gap
	v0 := math.Sqrt(2 * 2 * qe / me)
	// fly launches from wall 0 (x = 0, moving +x) or wall 1 (x = d, moving
	// −x) at time t0; returns impact time, wall and energy [eV], or ok =
	// false if the electron is still in flight after 30 periods.
	fly := func(t0 float64, wall int) (t1 float64, w1 int, e float64, ok bool) {
		x, v := 0.0, v0
		if wall == 1 {
			x, v = d, -v0
		}
		dt := 1 / f / 400
		t := t0
		for step := 0; step < 400*30; step++ {
			a := k * (vdc + vrf*math.Sin(w*t))
			// velocity Verlet
			xn := x + v*dt + 0.5*a*dt*dt
			an := k * (vdc + vrf*math.Sin(w*(t+dt)))
			vn := v + 0.5*(a+an)*dt
			if xn <= 0 || xn >= d {
				// interpolate the crossing
				var frac float64
				if xn <= 0 {
					frac = x / (x - xn)
					w1 = 0
				} else {
					frac = (d - x) / (xn - x)
					w1 = 1
				}
				vi := v + (vn-v)*frac
				return t + frac*dt, w1, 0.5 * me * vi * vi / qe, true
			}
			x, v, t = xn, vn, t+dt
		}
		return 0, 0, 0, false
	}
	growth = 0
	for wall := 0; wall < 2; wall++ {
		for ph := 0; ph < 72; ph++ {
			t, wl := float64(ph)/72/f, wall
			var sum, esum float64
			n := 0
			for imp := 0; imp < 20; imp++ {
				t1, w1, e, ok := fly(t, wl)
				if !ok {
					sum = math.Inf(-1)
					break
				}
				y := sey.Yield(e)
				if y <= 0 {
					sum = math.Inf(-1)
					break
				}
				sum += math.Log(y)
				esum += e
				n++
				t, wl = t1, w1
			}
			if n == 0 || math.IsInf(sum, -1) {
				continue
			}
			if g := math.Exp(sum / float64(n)); g > growth {
				growth, eImpact = g, esum/float64(n)
			}
		}
	}
	return growth, eImpact
}
