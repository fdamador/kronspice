// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"math"
	"math/rand"
	"testing"

	"github.com/fdamador/kronspice/hyper/quat"
)

func rvec(r *rand.Rand, s float64) [3]float64 {
	return [3]float64{s * r.NormFloat64(), s * r.NormFloat64(), s * r.NormFloat64()}
}

func near3(a, b [3]float64, tol float64) bool {
	return norm3(sub3(a, b)) <= tol*math.Max(1e-300, math.Max(norm3(a), norm3(b)))
}

// Grade selection ⟨F·Q†⟩ equals the vector forms for both pairs (CLAUDE.md
// §3f tests; Joot Theorem 3.16): q(E + v×B), qE·v; q_m(B − v×E/c²), q_m B·v.
func TestLorentzGradeSelection(t *testing.T) {
	r := rand.New(rand.NewSource(61))
	c := SpeedOfLight
	for n := 0; n < 2000; n++ {
		e, b, v := rvec(r, 1e3), rvec(r, 1e-3), rvec(r, 0.1*c)
		q := r.NormFloat64() * 1e-6
		fg, pg := LorentzGA(FieldMV(e, b), ElectricQ(q, v))
		fv, pv := LorentzVector(q, v, e, b)
		if !near3(fg, fv, 1e-12) || math.Abs(pg-pv) > 1e-12*math.Abs(q)*norm3(e)*norm3(v) {
			t.Fatalf("electric: GA %v %g, vector %v %g", fg, pg, fv, pv)
		}
		qm := r.NormFloat64()
		fg, pg = LorentzGA(FieldMV(e, b), MagneticQ(qm, v))
		fv, pv = MagneticLorentzVector(qm, v, e, b)
		if !near3(fg, fv, 1e-12) || math.Abs(pg-pv) > 1e-12*math.Abs(qm)*norm3(b)*norm3(v) {
			t.Fatalf("magnetic: GA %v %g, vector %v %g", fg, pg, fv, pv)
		}
	}
}

// Electric–magnetic duality: the force on a dyon is invariant under the
// duality rotation of fields and charges (Jackson eq. 6.151), and the
// rotation by π/2 maps the electric pair onto the magnetic pair.
func TestDuality(t *testing.T) {
	r := rand.New(rand.NewSource(62))
	h := Hypothetical{MagneticCharge: true}
	c := SpeedOfLight
	for n := 0; n < 500; n++ {
		e, b, v := rvec(r, 1e3), rvec(r, 1e-5), rvec(r, 0.05*c)
		q, qm := r.NormFloat64()*1e-6, r.NormFloat64()*1e2
		f, p, err := DyonForce(h, q, qm, v, e, b)
		if err != nil {
			t.Fatal(err)
		}
		xi := 2 * math.Pi * r.Float64()
		e2, b2, q2, qm2 := DualRotate(xi, e, b, q, qm)
		f2, p2, _ := DyonForce(h, q2, qm2, v, e2, b2)
		if !near3(f, f2, 1e-10) || math.Abs(p-p2) > 1e-10*(math.Abs(p)+1e-30) {
			t.Fatalf("ξ = %g: force %v → %v, power %g → %g", xi, f, f2, p, p2)
		}
	}
	// π/2: an electric charge in (E, B) feels what the dual magnetic charge
	// feels in the dual fields.
	e, b, v := [3]float64{1e3, 0, 0}, [3]float64{0, 2e-6, 0}, [3]float64{0, 0, 1e5}
	fe, _ := LorentzVector(1e-6, v, e, b)
	e2, b2, q2, qm2 := DualRotate(-math.Pi/2, e, b, 1e-6, 0)
	if math.Abs(q2) > 1e-20 {
		t.Fatalf("rotated electric charge %g", q2)
	}
	fm, _ := MagneticLorentzVector(qm2, v, e2, b2)
	if !near3(fe, fm, 1e-12) {
		t.Errorf("electric pair %v, dual magnetic pair %v", fe, fm)
	}
	if _, _, err := DyonForce(Hypothetical{}, 1, 1, v, e, b); !errors.Is(err, ErrHypotheticalDisabled) {
		t.Error("magnetic charge accepted outside hypothetical mode")
	}
}

// The homopolar EMF from the Lorentz force on the charges of a disk
// spinning in an axial field: integrating (v×B)·dr from axle to rim, with
// v = ω×r, gives ½Bω(r_o² − r_i²) (Tesla US 406,968; CLAUDE.md §3f).
func TestHomopolarEMF(t *testing.T) {
	const bz, w, ri, ro = 0.8, 300.0, 0.01, 0.15
	var path [][3]float64
	for k := 0; k <= 50; k++ {
		r := ri + (ro-ri)*float64(k)/50
		path = append(path, [3]float64{r * math.Cos(0.3), r * math.Sin(0.3), 0})
	}
	vel := func(p [3]float64) [3]float64 { return cross([3]float64{0, 0, w}, p) }
	field := func([3]float64) [3]float64 { return [3]float64{0, 0, bz} }
	got := MotionalEMF(path, vel, field)
	want := 0.5 * bz * w * (ro*ro - ri*ri)
	if math.Abs(got-want) > 1e-12*want {
		t.Errorf("EMF %.12g V, want ½Bω(r_o² − r_i²) = %.12g V", got, want)
	}
}

// Dual-symmetric Maxwell: inside a ball of uniform magnetic charge
// density, B = μ₀ρ_m·r/3 satisfies div B = μ₀ρ_m exactly; ordinary Maxwell
// flags the same field (S.R2 = −μ₀ρ_m ≠ 0), and magnetic sources outside
// hypothetical mode are refused.
func TestDualMaxwell(t *testing.T) {
	g := cubeGrid(7, 1e-3)
	const rhoM = 1e3
	B := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(x, y, z).Scale(Mu0 * rhoM / 3) })
	zero := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Q{} })
	rm := SampleS(g, func(x, y, z float64) float64 { return rhoM })
	d, err := MaxwellDual(Hypothetical{MagneticCharge: true}, zero, B, zero, zero, nil, nil, rm, nil)
	if err != nil || !d.Hypothetical {
		t.Fatal(err)
	}
	if m := maxInterior(g, func(n int) float64 { return d.R2.Q[n].Abs() }); m > 1e-9*Mu0*rhoM {
		t.Errorf("dual residual %g", m)
	}
	std := MaxwellQuaternion(zero, B, zero, zero, nil, nil)
	if m := maxInterior(g, func(n int) float64 { return math.Abs(std.R2.Q[n].W + Mu0*rhoM) }); m > 1e-9*Mu0*rhoM {
		t.Errorf("standard Maxwell should show −div B = −μ₀ρ_m: off by %g", m)
	}
	if _, err := MaxwellDual(Hypothetical{}, zero, B, zero, zero, nil, nil, rm, nil); !errors.Is(err, ErrHypotheticalDisabled) {
		t.Error("magnetic charge density accepted outside hypothetical mode")
	}
}

// A changing E field makes B (Ampère–Maxwell): between the plates of a
// charging circular capacitor, E = Ė·t ẑ and B = μ₀ε₀Ė·(−y, x, 0)/2, which
// satisfy all four equations with no conduction current (M9h operators).
// And Rowland's rotating charged disk: B on its axis from summing the
// current rings equals (μ₀σω/2)[(R² + 2z²)/√(R² + z²) − 2|z|].
func TestFieldsFromChangingE(t *testing.T) {
	g := cubeGrid(9, 1e-3)
	const edot = 1e12 // V/(m·s)
	at := func(tt float64) (QField, QField) {
		return SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(0, 0, edot*tt) }),
			SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(-y, x, 0).Scale(Mu0 * Eps0 * edot / 2) })
	}
	E, B := at(1e-6)
	dE := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Vec(0, 0, edot) })
	zero := SampleQ(g, func(x, y, z float64) quat.Q { return quat.Q{} })
	r := MaxwellQuaternion(E, B, dE, zero, nil, nil)
	scale := Mu0 * Eps0 * edot
	if m := maxInterior(g, func(n int) float64 { return r.R2.Q[n].Abs() }); m > 1e-9*scale {
		t.Errorf("charging capacitor: Ampère–Maxwell residual %g of %g", m, scale)
	}
	if m := maxInterior(g, func(n int) float64 { return r.R1.Q[n].Abs() }); m > 1e-9 {
		t.Errorf("charging capacitor: Faraday/Gauss residual %g", m)
	}
	// Rowland disk: surface charge σ, radius R, angular velocity ω.
	const sigma, R, w = 1e-5, 0.1, 2000.0
	const rings = 400
	for _, z := range []float64{0, 0.02, 0.1} {
		var bz float64
		for k := 0; k < rings; k++ {
			rr := R * (float64(k) + 0.5) / rings
			dI := sigma * w * rr * R / rings
			l := Loop{P: Circle([3]float64{}, [3]float64{0, 0, 1}, rr, 720), I: dI}
			bz += l.B([3]float64{0, 0, z})[2]
		}
		want := Mu0 * sigma * w / 2 * ((R*R+2*z*z)/math.Sqrt(R*R+z*z) - 2*math.Abs(z))
		if math.Abs(bz-want) > 1e-4*want {
			t.Errorf("Rowland disk B_z(%g) = %g T, want %g T", z, bz, want)
		}
	}
}

// Momentum ledger: two coaxial current loops push each other with equal
// and opposite forces (Newton's third law for closed circuits), the force
// at large spacing tends to the dipole–dipole value 3μ₀m₁m₂/(2πz⁴), and a
// force booked from an exterior source appears in Net and ByExterior.
func TestMomentumLedger(t *testing.T) {
	a := Loop{P: Circle([3]float64{}, [3]float64{0, 0, 1}, 0.05, 360), I: 10}
	for _, z := range []float64{0.03, 2.0} {
		b := Loop{P: Circle([3]float64{0.004, -0.002, z}, [3]float64{0, 0.1, 1}, 0.04, 360), I: -7}
		var l MomentumLedger
		l.AddInternal("A", "B", a.ForceIn(b.B))
		l.AddInternal("B", "A", b.ForceIn(a.B))
		if err := l.Check(1e-6); err != nil {
			t.Errorf("z = %g: %v", z, err)
		}
		if z == 2.0 {
			// Dipole–dipole (coaxial, aligned up to the small tilt).
			b2 := Loop{P: Circle([3]float64{0, 0, z}, [3]float64{0, 0, 1}, 0.04, 360), I: -7}
			fz := a.ForceIn(b2.B)[2]
			m1, m2 := 10*math.Pi*0.05*0.05, -7*math.Pi*0.04*0.04
			// Force on loop 1 from a coaxial dipole at +z: +3μ₀m₁m₂/(2πz⁴)·ẑ
			// (parallel moments attract; these are antiparallel and repel).
			want := 3 * Mu0 * m1 * m2 / (2 * math.Pi * math.Pow(z, 4))
			if math.Abs(fz-want) > 1e-2*math.Abs(want) {
				t.Errorf("far-field force %g N, dipole value %g N", fz, want)
			}
		}
	}
	var l MomentumLedger
	l.AddInternal("A", "B", [3]float64{1, 0, 0})
	l.AddInternal("B", "A", [3]float64{-1, 0, 0})
	l.AddExterior("tether", "geomagnetic field", [3]float64{0, 0.2, 0}, false)
	l.AddExterior("tether", "monopole (hypothetical)", [3]float64{5, 5, 5}, true)
	if n := l.Net(); !near3(n, [3]float64{0, 0.2, 0}, 1e-15) {
		t.Errorf("net %v", n)
	}
	if h := l.Hypothetical(); h != [3]float64{5, 5, 5} {
		t.Errorf("hypothetical %v", h)
	}
	if ex := l.ByExterior(); len(ex) != 1 || ex["geomagnetic field"] != [3]float64{0, 0.2, 0} {
		t.Errorf("by exterior %v", ex)
	}
	var bad MomentumLedger
	bad.AddInternal("A", "B", [3]float64{1, 0, 0})
	if bad.Check(1e-9) == nil {
		t.Error("an uncancelled internal force passed")
	}
}
