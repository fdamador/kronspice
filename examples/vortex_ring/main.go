// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command vortex_ring computes Kelvin's vortex atom — a closed vortex ring
// in a perfect liquid ("On Vortex Atoms", Proc. R. Soc. Edinb. 6, 1867;
// Mathematical and Physical Papers Vol. IV, 1910, Art. 1, pp. 1–12) —
// through the analogy Kelvin draws in that paper: the velocity of the
// liquid round a vortex ring is the magnetic field of the "corresponding
// ring electromagnet", and the ring's impulse ("effective momentum") equals
// that electromagnet's magnetic moment (pp. 9–10). Circulation Γ plays the
// current I, velocity u plays B/μ₀, kinetic energy ½ρ∫u² plays the
// magnetic energy ½LI²/μ₀·μ₀ — so the engine's inductance and Biot–Savart
// code (em.MutualCircles, em.Loop) computes the vortex ring directly.
//
//  1. Kelvin's ring speed. The energy of a ring of radius R with a core of
//     radius a and uniform vorticity is ½ρΓ²·L/μ₀, L the self-inductance of
//     a ring conductor with uniform current. L is computed here by averaging
//     Maxwell's mutual inductance of coaxial circles over the core cross
//     section; the speed is U = dE/dP at fixed circulation and core volume
//     (P = ρΓπR², the impulse) and is compared with Kelvin's
//     U = Γ/4πR·(ln 8R/a − ¼) (Kelvin 1867, Vol. IV Art. 3 p. 67, as quoted
//     by Saffman, Vortex Dynamics §10.2; to be checked against the page).
//
//  2. Leapfrogging. Two coaxial rings (Helmholtz 1858; Tait's smoke rings,
//     Kelvin pp. 2, 11–12) pass through each other in turn. Each ring moves
//     with Kelvin's self speed plus the velocity of the other ring's flow
//     (Biot–Savart, em.Loop); impulse and energy must stay constant.
//
//  3. Where the analogy stops. The flow and the magnetic field share their
//     shape, not their dynamics: vorticity is carried along by the flow it
//     makes, so the ring travels; a current loop's field does not carry the
//     loop, so the "ring electromagnet" stays put.
//
//     go run ./examples/vortex_ring
package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/em"
)

type v3 = [3]float64

// gaussLegendre returns n nodes and weights on (−1, 1).
func gaussLegendre(n int) (x, w []float64) {
	x, w = make([]float64, n), make([]float64, n)
	for i := 0; i < (n+1)/2; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var dp float64
		for it := 0; it < 100; it++ {
			p0, p1 := 1.0, z
			for k := 2; k <= n; k++ {
				p0, p1 = p1, ((2*float64(k)-1)*z*p1-(float64(k)-1)*p0)/float64(k)
			}
			dp = float64(n) * (z*p1 - p0) / (z*z - 1)
			dz := p1 / dp
			z -= dz
			if math.Abs(dz) < 1e-15 {
				break
			}
		}
		x[i], x[n-1-i] = -z, z
		w[i] = 2 / ((1 - z*z) * dp * dp)
		w[n-1-i] = w[i]
	}
	return x, w
}

// diskPoints returns quadrature points (ρ, φ) and area weights for a disk
// of radius a: Gauss–Legendre in ρ², midpoint in φ (offset rotates the
// angular grid, so two sets never share a point).
func diskPoints(a float64, nr, nphi int, offset float64) (r, z, w []float64) {
	x, wx := gaussLegendre(nr)
	for i := range x {
		s := (x[i] + 1) / 2 // ρ²/a² uniform in area
		rho := a * math.Sqrt(s)
		for j := 0; j < nphi; j++ {
			ph := 2 * math.Pi * (float64(j) + offset) / float64(nphi)
			r = append(r, rho*math.Cos(ph))
			z = append(z, rho*math.Sin(ph))
			w = append(w, wx[i]/2/float64(nphi)) // fraction of the disk's area
		}
	}
	return r, z, w
}

// SelfInductance of a circular ring of radius R whose conductor has a
// round section of radius a carrying uniform current: the mutual
// inductance of coaxial circles (Maxwell, Treatise Art. 701) averaged over
// pairs of filaments in the section.
func SelfInductance(R, a float64, nr, nphi int) float64 {
	r1, z1, w1 := diskPoints(a, nr, nphi, 0)
	r2, z2, w2 := diskPoints(a, nr, nphi, 0.5)
	var L float64
	for i := range r1 {
		for j := range r2 {
			L += w1[i] * w2[j] * em.MutualCircles(R+r1[i], R+r2[j], z1[i]-z2[j])
		}
	}
	return L
}

// ThinRingInductance is μ₀R(ln 8R/a − 7/4), the uniform-current ring in
// the thin limit (internal inductance μ₀R/4 plus the external part).
func ThinRingInductance(R, a float64) float64 { return em.Mu0 * R * (math.Log(8*R/a) - 1.75) }

// KelvinSpeed is Γ/4πR·(ln 8R/a − ¼).
func KelvinSpeed(gamma, R, a float64) float64 {
	return gamma / (4 * math.Pi * R) * (math.Log(8*R/a) - 0.25)
}

// Speed is the ring speed U = dE/dP at fixed Γ and core volume
// (a²R constant), from the numerical self-inductance: E = ½ρΓ²L/μ₀,
// P = ρΓπR² (ρ cancels).
func Speed(gamma, R, a float64, nr, nphi int) float64 {
	h := 1e-3 * R
	E := func(r float64) float64 {
		return 0.5 * gamma * gamma * SelfInductance(r, a*math.Sqrt(R/r), nr, nphi) / em.Mu0
	}
	P := func(r float64) float64 { return gamma * math.Pi * r * r }
	return (E(R+h) - E(R-h)) / (P(R+h) - P(R-h))
}

// Velocity of the liquid at p from a vortex filament loop of circulation
// Γ: Biot–Savart, i.e. B/μ₀ of the loop carrying current Γ.
func Velocity(loop em.Loop, p v3) v3 {
	b := loop.B(p)
	return v3{b[0] / em.Mu0, b[1] / em.Mu0, b[2] / em.Mu0}
}

// Ring is one coaxial ring of a leapfrog run.
type Ring struct{ R, Z, A float64 } // radius, axial position, core radius

// Leapfrog is the outcome of a two-ring run.
type Leapfrog struct {
	Passes             int     // times one ring passed through the other
	Period             float64 // mean time between passes
	ImpulseDrift       float64 // max |ΔΣR²|/ΣR²
	EnergyDrift        float64 // max |ΔH|/|H|
	Advance, MeanSpeed float64 // distance travelled by the pair, mean speed
	SingleSpeed        float64 // speed of one ring alone (Kelvin)
	RMin, RMax         float64 // range of ring radii
}

// RunLeapfrog integrates two coaxial rings of equal circulation Γ with
// RK4. Each ring keeps its core volume (a²R constant), moves at Kelvin's
// self speed and is carried by the other ring's flow at its own centre
// line. H = Σ ½Γ²R(ln 8R/a − 7/4) + Γ²M₁₂/μ₀ (ρ = 1) is the energy.
func RunLeapfrog(gamma float64, r0 [2]Ring, T, dt float64, sides int) Leapfrog {
	vol := [2]float64{r0[0].A * r0[0].A * r0[0].R, r0[1].A * r0[1].A * r0[1].R}
	core := func(k int, R float64) float64 { return math.Sqrt(vol[k] / R) }
	type st [4]float64 // R0, Z0, R1, Z1
	der := func(s st) st {
		var d st
		for k := 0; k < 2; k++ {
			j := 1 - k
			Rk, Zk, Rj, Zj := s[2*k], s[2*k+1], s[2*j], s[2*j+1]
			loop := em.Loop{P: em.Circle(v3{0, 0, Zj}, v3{0, 0, 1}, Rj, sides), I: gamma}
			u := Velocity(loop, v3{Rk, 0, Zk})
			d[2*k] = u[0]
			d[2*k+1] = KelvinSpeed(gamma, Rk, core(k, Rk)) + u[2]
		}
		return d
	}
	H := func(s st) float64 {
		h := gamma * gamma * em.MutualCircles(s[0], s[2], s[1]-s[3]) / em.Mu0
		for k := 0; k < 2; k++ {
			R := s[2*k]
			h += 0.5 * gamma * gamma * R * (math.Log(8*R/core(k, R)) - 1.75)
		}
		return h
	}
	s := st{r0[0].R, r0[0].Z, r0[1].R, r0[1].Z}
	h0, p0 := H(s), s[0]*s[0]+s[2]*s[2]
	res := Leapfrog{RMin: math.Inf(1), SingleSpeed: KelvinSpeed(gamma, r0[0].R, r0[0].A)}
	ahead := s[1] > s[3]
	var passT []float64
	z0 := (s[1] + s[3]) / 2
	n := int(T / dt)
	for i := 0; i < n; i++ {
		k1 := der(s)
		var t2, t3, t4 st
		for m := range s {
			t2[m] = s[m] + dt/2*k1[m]
		}
		k2 := der(t2)
		for m := range s {
			t3[m] = s[m] + dt/2*k2[m]
		}
		k3 := der(t3)
		for m := range s {
			t4[m] = s[m] + dt*k3[m]
		}
		k4 := der(t4)
		for m := range s {
			s[m] += dt / 6 * (k1[m] + 2*k2[m] + 2*k3[m] + k4[m])
		}
		res.ImpulseDrift = math.Max(res.ImpulseDrift, math.Abs(s[0]*s[0]+s[2]*s[2]-p0)/p0)
		res.EnergyDrift = math.Max(res.EnergyDrift, math.Abs(H(s)-h0)/math.Abs(h0))
		res.RMin = math.Min(res.RMin, math.Min(s[0], s[2]))
		res.RMax = math.Max(res.RMax, math.Max(s[0], s[2]))
		if now := s[1] > s[3]; now != ahead {
			ahead = now
			res.Passes++
			passT = append(passT, float64(i+1)*dt)
		}
	}
	if len(passT) > 1 {
		res.Period = (passT[len(passT)-1] - passT[0]) / float64(len(passT)-1)
	}
	res.Advance = (s[1]+s[3])/2 - z0
	res.MeanSpeed = res.Advance / T
	return res
}

func main() {
	fmt.Println("Kelvin, \"On Vortex Atoms\" (1867), Math. Phys. Papers IV Art. 1 — through his ring-electromagnet analogy")
	fmt.Println()
	fmt.Println("1. Ring speed from the engine's inductance (Γ = 1, R = 1)")
	for _, a := range []float64{0.2, 0.1, 0.05} {
		L := SelfInductance(1, a, 24, 48)
		u := Speed(1, 1, a, 24, 48)
		fmt.Printf("   a = %.2f: L/μ₀R = %.5f (thin-ring ln 8R/a − 7/4 = %.5f); U = %.5f vs Kelvin %.5f (%+.2f %%)\n",
			a, L/em.Mu0, math.Log(8/a)-1.75, u, KelvinSpeed(1, 1, a), 100*(u/KelvinSpeed(1, 1, a)-1))
	}
	loop := em.Loop{P: em.Circle(v3{}, v3{0, 0, 1}, 1, 720), I: 1}
	fmt.Printf("   flow at the ring's centre %.5f = Γ/2R; impulse P = ρΓπR² ↔ magnetic moment IπR² = %.5f (Kelvin p. 10)\n",
		Velocity(loop, v3{})[2], math.Pi)
	fmt.Println()
	fmt.Println("2. Leapfrogging: two rings Γ = 1, R = 1, core a = 0.1, 0.6 apart")
	lf := RunLeapfrog(1, [2]Ring{{R: 1, Z: 0, A: 0.1}, {R: 1, Z: 0.6, A: 0.1}}, 60, 0.02, 360)
	fmt.Printf("   %d passes in 60 time units, period %.2f; radii swing %.3f…%.3f\n", lf.Passes, lf.Period, lf.RMin, lf.RMax)
	fmt.Printf("   pair advances at %.4f (single ring %.4f); impulse drift %.1e, energy drift %.1e\n", lf.MeanSpeed, lf.SingleSpeed, lf.ImpulseDrift, lf.EnergyDrift)
	fmt.Println()
	fmt.Println("3. Same shape, different dynamics: the ring travels at U because the liquid's vorticity is carried by")
	fmt.Println("   the flow it makes; a current loop with the identical (magnetic) field is not moved by that field.")
}
