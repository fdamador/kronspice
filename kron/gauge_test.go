// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"
	"math/rand/v2"
	"testing"
)

// The cell complex of a Kron network is its connection tensors: the
// plaquettes' signed incidence is C, Curl∘Grad = 0 is Aᵀ·C = 0, and
// Div∘Curlᵀ = 0 is Kirchhoff's current law holding identically for mesh
// currents (Kron 1939 p. 75, eqs 3.1–3.2).
func TestLatticeIsKronsConnection(t *testing.T) {
	r := rand.New(rand.NewPCG(1, 2))
	for trial := 0; trial < 50; trial++ {
		nw, err := Build(parse(t, randomCircuit(r, 3+r.IntN(6))))
		if err != nil {
			t.Fatal(err)
		}
		lat, err := LatticeFromNetwork(nw)
		if err != nil {
			t.Fatal(err)
		}
		cnt := nw.Counts()
		if len(lat.Plaquettes) != cnt.Meshes || lat.Nodes != cnt.Junctions || lat.Links() != cnt.Branches {
			t.Fatalf("trial %d: lattice %d/%d/%d, Kron b/N/m = %d/%d/%d", trial,
				lat.Links(), lat.Nodes, len(lat.Plaquettes), cnt.Branches, cnt.Junctions, cnt.Meshes)
		}
		c := nw.g.loops()
		for p, path := range lat.Plaquettes {
			inc := make([]float64, lat.Links())
			for _, s := range path {
				inc[s.Link] += float64(s.Dir)
			}
			for k := range inc {
				if inc[k] != real(c.At(k, p)) {
					t.Fatalf("trial %d: plaquette %d link %d incidence %g, C %g", trial, p, k, inc[k], real(c.At(k, p)))
				}
			}
		}
		ip := make([]float64, len(lat.Plaquettes))
		for i := range ip {
			ip[i] = r.NormFloat64()
		}
		for n, v := range lat.Div(lat.CurlT(ip)) {
			if math.Abs(v) > 1e-12 {
				t.Fatalf("trial %d: KCL residual %g at node %d", trial, v, n)
			}
		}
	}
}

// A gauge transformation a → a + Grad χ leaves every holonomy, every
// Wilson loop and the Gauss generator unchanged; a flux of one flux
// quantum turns a Cooper pair's Wilson loop through exactly 2π.
func TestU1GaugeInvariance(t *testing.T) {
	r := rand.New(rand.NewPCG(3, 4))
	nw, err := Build(parse(t, randomCircuit(r, 6)))
	if err != nil {
		t.Fatal(err)
	}
	lat, err := LatticeFromNetwork(nw)
	if err != nil {
		t.Fatal(err)
	}
	a := make([]float64, lat.Links())
	for k := range a {
		a[k] = r.NormFloat64() * FluxQuantum
	}
	chi := make([]float64, lat.Nodes)
	for n := range chi {
		chi[n] = r.NormFloat64() * FluxQuantum
	}
	b := lat.Gauge(a, chi)
	ha, hb := lat.Curl(a), lat.Curl(b)
	for p := range ha {
		if d := math.Abs(ha[p]-hb[p]) / FluxQuantum; d > 1e-12 {
			t.Errorf("plaquette %d: holonomy changed by %g Φ₀ under a gauge transformation", p, d)
		}
		wa, wb := lat.WilsonLoop(a, p, 2*ElementaryCharge), lat.WilsonLoop(b, p, 2*ElementaryCharge)
		if cmplx.Abs(wa-wb) > 1e-9 {
			t.Errorf("plaquette %d: Wilson loop %v → %v", p, wa, wb)
		}
	}
	// One flux quantum through plaquette 0 and nothing else: set a on its
	// first link only (the link belongs to that fundamental loop alone).
	one := make([]float64, lat.Links())
	s := lat.Plaquettes[0][0]
	one[s.Link] = float64(s.Dir) * FluxQuantum
	if w := lat.WilsonLoop(one, 0, 2*ElementaryCharge); cmplx.Abs(w-1) > 1e-9 {
		t.Errorf("Wilson loop at one flux quantum = %v, want 1", w)
	}
	if w := lat.WilsonLoop(one, 0, ElementaryCharge); cmplx.Abs(w+1) > 1e-9 {
		t.Errorf("single-charge Wilson loop at Φ₀ = %v, want −1", w)
	}
}

// Yee's cavity on the gauge lattice: the TM110 mode of a closed box rings
// at the Yee scheme's own frequency, the exact discrete energy is
// constant to round-off, and the Gauss generator stays zero at every
// interior node (the wall nodes carry the induced surface charge).
func TestU1YeeCavity(t *testing.T) {
	const c0, mu0 = 299792458.0, 1.25663706212e-6
	eps0 := 1 / (mu0 * c0 * c0)
	n := [3]int{21, 17, 6}
	dx := 0.01
	g, err := NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		t.Fatal(err)
	}
	cap, ind := g.Constitutive(eps0, mu0)
	dt := 0.99 / (c0 * math.Sqrt(3/(dx*dx)))
	f, err := NewU1Field(g.Lattice, cap, ind, dt)
	if err != nil {
		t.Fatal(err)
	}
	a, b := float64(n[0]-1)*dx, float64(n[1]-1)*dx
	v := make([]float64, g.Links())
	for k := range v {
		if g.Axis(k) == 2 {
			p := g.Position(k)
			v[k] = math.Sin(math.Pi*p[0]/a) * math.Sin(math.Pi*p[1]/b) * dx
		}
	}
	f.SetVoltage(v)
	we0, wm0 := f.Energy()
	if wm0 != 0 || we0 <= 0 {
		t.Fatalf("initial energy electric %g magnetic %g", we0, wm0)
	}
	interior := func(node int) bool {
		i, j, k := node%n[0], (node/n[0])%n[1], node/(n[0]*n[1])
		return i > 0 && i < n[0]-1 && j > 0 && j < n[1]-1 && k > 0 && k < n[2]-1
	}
	probe := g.Link(2, n[0]/2, n[1]/2, 2)
	prev := f.Voltage()[probe]
	var cross []float64
	worstW, worstG := 0.0, 0.0
	for f.Steps < 4000 {
		f.Step()
		vv := f.Voltage()[probe]
		if (prev > 0) != (vv > 0) {
			cross = append(cross, f.Time()-f.Dt*vv/(vv-prev))
		}
		prev = vv
		we, wm := f.Energy()
		worstW = math.Max(worstW, math.Abs(we+wm-we0)/we0)
		if f.Steps%500 == 0 {
			qmax := 0.0
			for _, q := range f.Q {
				qmax = math.Max(qmax, math.Abs(q))
			}
			for node, gg := range f.Gauss(nil) {
				if interior(node) {
					worstG = math.Max(worstG, math.Abs(gg)/qmax)
				}
			}
		}
	}
	period := 2 * (cross[len(cross)-1] - cross[0]) / float64(len(cross)-1)
	w := 2 * math.Pi / period
	s := c0 * dt * math.Sqrt(math.Pow(2/dx*math.Sin(math.Pi*dx/(2*a)), 2)+math.Pow(2/dx*math.Sin(math.Pi*dx/(2*b)), 2)) / 2
	wYee := 2 / dt * math.Asin(s)
	wCont := c0 * math.Pi * math.Hypot(1/a, 1/b)
	if r := math.Abs(w-wYee) / wYee; r > 1e-5 {
		t.Errorf("ω %.8g, Yee %.8g (%.2g)", w, wYee, r)
	}
	if r := math.Abs(w-wCont) / wCont; r > 2e-3 {
		t.Errorf("ω %.8g, continuum %.8g (%.2g)", w, wCont, r)
	}
	if worstW > 1e-12 {
		t.Errorf("energy drift %.3g", worstW)
	}
	if worstG > 1e-12 {
		t.Errorf("Gauss generator %.3g of max Q at an interior node", worstG)
	}
	t.Logf("TM110 on the gauge lattice: ω/ω_Yee − 1 = %.2g, ω/ω_cont − 1 = %.2g, energy drift %.2g, Gauss %.2g",
		w/wYee-1, w/wCont-1, worstW, worstG)
}

// The SU(3) algebra and field: generators normalized, exponentials in the
// group, Wilson-loop traces gauge invariant, the Gauss generator
// covariant, and the composition defect zero exactly when the field is
// abelian (Cartan).
func TestSU3(t *testing.T) {
	for a := 0; a < 8; a++ {
		for b := 0; b < 8; b++ {
			want := 0.0
			if a == b {
				want = 2
			}
			if tr := GellMann[a].Mul(GellMann[b]).Trace(); cmplx.Abs(tr-complex(want, 0)) > 1e-15 {
				t.Errorf("Tr(λ%d λ%d) = %v, want %g", a+1, b+1, tr, want)
			}
		}
	}
	// exp(i θ λ₃/2) = diag(e^{iθ/2}, e^{−iθ/2}, 1).
	theta := 2.3
	u := Exp(Algebra([8]float64{0, 0, theta}))
	want := SU3{{cmplx.Exp(complex(0, theta/2)), 0, 0}, {0, cmplx.Exp(complex(0, -theta/2)), 0}, {0, 0, 1}}
	if d := u.Sub(want).Norm(); d > 1e-14 {
		t.Errorf("exp(iθλ₃/2) off by %g", d)
	}
	r := rand.New(rand.NewPCG(5, 6))
	for i := 0; i < 100; i++ {
		if d := RandomSU3(r, 3, false).UnitaryDefect(); d > 1e-12 {
			t.Errorf("random element: unitary defect %g", d)
		}
	}

	nw, err := Build(parse(t, "two meshes\nV1 1 0 1\nR1 1 2 10\nL1 2 0 1m\nL2 2 3 2m\nR2 3 0 5\n"))
	if err != nil {
		t.Fatal(err)
	}
	lat, err := LatticeFromNetwork(nw)
	if err != nil {
		t.Fatal(err)
	}
	if len(lat.Plaquettes) != 2 {
		t.Fatalf("%d plaquettes, want 2", len(lat.Plaquettes))
	}
	f := NewSU3Field(lat)
	f.Randomize(r, 1, false)
	e := make([]SU3, lat.Links())
	for k := range e {
		e[k] = RandomAlgebra(r, 1, false)
	}
	rho := make([]SU3, lat.Nodes)
	for n := range rho {
		rho[n] = RandomAlgebra(r, 1, false)
	}
	tr0 := []complex128{f.WilsonLoop(0).Trace(), f.WilsonLoop(1).Trace()}
	g0 := f.Gauss(e, rho)
	gauge := make([]SU3, lat.Nodes)
	for n := range gauge {
		gauge[n] = RandomSU3(r, 1, false)
	}
	f.Gauge(gauge)
	for k := range e {
		gk := gauge[lat.From[k]]
		e[k] = gk.Mul(e[k]).Mul(gk.Dagger())
	}
	for n := range rho {
		rho[n] = gauge[n].Mul(rho[n]).Mul(gauge[n].Dagger())
	}
	for p, tr := range tr0 {
		if d := cmplx.Abs(f.WilsonLoop(p).Trace() - tr); d > 1e-12 {
			t.Errorf("Tr W_%d changed by %g under a gauge transformation", p, d)
		}
	}
	for n, g1 := range f.Gauss(e, rho) {
		cov := gauge[n].Mul(g0[n]).Mul(gauge[n].Dagger())
		if d := g1.Sub(cov).Norm(); d > 1e-12 {
			t.Errorf("node %d: Gauss generator is not covariant (defect %g)", n, d)
		}
		if n > 0 && g1.Sub(g0[n]).Norm() < 1e-6 {
			t.Errorf("node %d: Gauss generator is invariant, which only an abelian field gives", n)
		}
	}

	// Composition at the node the two meshes share (node 2, where R1, L1
	// and L2 meet).
	node := -1
	for n := 0; n < lat.Nodes && node < 0; n++ {
		if _, ok1 := lat.Rebase(lat.Plaquettes[0], n); ok1 {
			if _, ok2 := lat.Rebase(lat.Plaquettes[1], n); ok2 {
				node = n
			}
		}
	}
	if node < 0 {
		t.Fatal("the two meshes share no node")
	}
	_, _, defect, ok := f.Compose(0, 1, node)
	if !ok || defect < 0.05 {
		t.Errorf("SU(3) composition defect %g (ok %v), expected order 1", defect, ok)
	}
	f.Randomize(r, 1, true)
	_, _, abelian, _ := f.Compose(0, 1, node)
	if abelian > 1e-14 {
		t.Errorf("Cartan (abelian) composition defect %g, want 0", abelian)
	}
	// In the abelian field the holonomy is additive: the eigenphases of
	// W_p are the signed sums of the link angles.
	angles := make([]float64, lat.Links())
	for k, u := range f.U {
		angles[k] = cmplx.Phase(u[0][0])
	}
	for p := range lat.Plaquettes {
		w := f.WilsonLoop(p)
		sum := lat.Curl(angles)[p]
		if d := cmplx.Abs(w[0][0] - cmplx.Exp(complex(0, sum))); d > 1e-12 {
			t.Errorf("plaquette %d: abelian holonomy %v, sum of link angles gives %v", p, w[0][0], cmplx.Exp(complex(0, sum)))
		}
	}
}
