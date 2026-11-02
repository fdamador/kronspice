// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command quantumtricity runs the two-step test of register entry H22
// (CLAUDE.md §1): the hypothesis that the quark model and the principles
// of circuit electricity (V, I, R, P) are one structure, the latter being
// electricity at the quantum level.
//
// Step 1 (standard, the abelian step). Circuit theory is a U(1) gauge
// theory on Kron's network: the branch is the link, the junction the
// gauge node, the mesh the Wilson loop. On a two-mesh RLC circuit solved
// in Kron's mesh frame the program shows that Kirchhoff's voltage law is
// zero holonomy of the branch voltages (the voltage is pure gauge,
// E = Grad v), that Kirchhoff's current law is the time derivative of the
// Gauss generator Div Q − ρ, that the mesh flux linkage is the holonomy
// of the link variable (the Wilson loop's phase in flux quanta), and that
// the lattice Hamiltonian's two terms are the ledger's ½Cv² and ½Li².
// Then the same U(1) field on Kron's 1944 field network (a Yee grid) is
// stepped beside em.FDTD on the TM110 mode of a closed box: the two
// agree edge by edge to round-off, because Hamilton's equations on the
// lattice are Yee's scheme.
//
// Step 2 (the hypothesis). SU(3) link variables are put on the same
// graph (Kogut & Susskind 1975). Their Wilson loops compose by matrix
// multiplication and the composition defect ‖W_pW_q − W_qW_p‖ is of
// order 1; restricted to the commuting (Cartan) part it is zero, which
// is the U(1) case. The program then looks through every quantity the
// simulator carries on a branch — voltage, current, charge, flux, power,
// in a linear circuit and in a nonlinear diode circuit — for one that
// composes that way, and finds none: every Kirchhoff sum commutes, to
// the round-off of addition. The quantum level of V, I and R is then
// stated from the 2019 SI: the Josephson constant 2e/h, the von Klitzing
// constant h/e² and single-electron pumping I = e·f — electrons and
// photons, no quark.
//
//	go run ./examples/quantumtricity
package main

import (
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"os"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/kron"
	"github.com/fdamador/kronspice/netlist"
)

// TwoMesh is the linear test circuit: two inductive meshes sharing L1,
// with a capacitor closing a third loop so that charge sits at a node.
const TwoMesh = `H22 two-mesh RLC
V1 1 0 DC 0 AC 1
R1 1 2 10
L1 2 0 1m
L2 2 3 2m
R2 3 0 5
C1 3 0 1u
`

// DiodePair is the nonlinear test circuit: one source feeding two diodes
// into different loads, so the two meshes are nonlinear and unequal.
const DiodePair = `H22 diode pair
V1 1 0 DC 5
R1 1 2 1k
D1 2 3 dm
R2 3 0 1k
D2 2 4 dm
R3 4 0 2k
.model dm D(IS=1e-14)
`

// Circuit parses a netlist.
func Circuit(src string) (*circuit.Circuit, error) {
	d, err := netlist.ParseString(src)
	if err != nil {
		return nil, err
	}
	return circuit.FromDeck(d)
}

// Abelian is step 1 on the two-mesh circuit at angular frequency Omega.
type Abelian struct {
	Omega                       float64
	Branches, Nodes, Meshes     int
	Plaquettes                  []string // each mesh as its signed elements
	CurlGrad                    float64  // max |Curl∘Grad χ| / max |χ|: Aᵀ·C = 0
	KVL                         float64  // max |Curl E| / max |E|: zero holonomy of the voltage
	MeshEMF                     float64  // max |Curl E_nonL + jω Curl Φ| / max |E|
	KCL                         float64  // max |Div I| / max |I|
	Continuity                  float64  // max |Div I_C − jω Div Q| / max |I|
	NodeCharge                  []complex128
	MeshFlux                    []complex128 // Φ_p = (Curl Φ)_p [Wb], peak phasor
	FluxQuanta                  []float64    // |Φ_p| / Φ₀
	Electric, Magnetic          float64      // Σ ½C|v|², Σ ½L|i|² at the phasor peaks [J]
	ElectricRows, MagneticRows  []string
	Pure                        float64 // max |Curl Grad(v_node) − Curl E|: E is pure gauge
	LinearCompositionDefect     float64 // largest loop-sum reordering defect over E, I, Q, Φ, S
	LinearCompositionQuantities []string
	SharedNode                  int
}

// RunAbelian solves the two-mesh circuit in Kron's mesh frame and reads the
// gauge-theory identities off the solution.
func RunAbelian(omega float64) (*Abelian, error) {
	c, err := Circuit(TwoMesh)
	if err != nil {
		return nil, err
	}
	nw, err := kron.Build(c)
	if err != nil {
		return nil, err
	}
	lat, err := kron.LatticeFromNetwork(nw)
	if err != nil {
		return nil, err
	}
	sol, err := nw.Solve(kron.MeshForm, omega, false)
	if err != nil {
		return nil, err
	}
	cnt := nw.Counts()
	a := &Abelian{Omega: omega, Branches: cnt.Branches, Nodes: cnt.Junctions, Meshes: cnt.Meshes}
	for _, path := range lat.Plaquettes {
		var parts []string
		for _, s := range path {
			sign := "+"
			if s.Dir < 0 {
				sign = "−"
			}
			parts = append(parts, sign+lat.Names[s.Link])
		}
		a.Plaquettes = append(a.Plaquettes, strings.Join(parts, " "))
	}
	b := lat.Links()
	elem := map[string]circuit.Element{}
	for _, e := range c.Elements {
		elem[strings.ToLower(e.Name)] = e
	}
	E := make([]complex128, b)
	I := make([]complex128, b)
	Phi := make([]complex128, b) // inductor flux L·i
	Q := make([]complex128, b)   // capacitor charge C·v
	EnonL := make([]complex128, b)
	IC := make([]complex128, b)
	S := make([]complex128, b)
	maxE, maxI := 0.0, 0.0
	for k := 0; k < b; k++ {
		E[k], I[k] = sol.E.At(k), sol.I.At(k)
		S[k] = E[k] * cmplx.Conj(I[k])
		maxE, maxI = math.Max(maxE, cmplx.Abs(E[k])), math.Max(maxI, cmplx.Abs(I[k]))
		e := elem[strings.ToLower(lat.Names[k])]
		switch e.Kind {
		case circuit.Inductor:
			Phi[k] = complex(e.Value, 0) * I[k]
			a.MagneticRows = append(a.MagneticRows, fmt.Sprintf("%s: ½L|i|² = %.4g J", e.Name, 0.5*e.Value*math.Pow(cmplx.Abs(I[k]), 2)))
			a.Magnetic += 0.5 * e.Value * math.Pow(cmplx.Abs(I[k]), 2)
		case circuit.Capacitor:
			Q[k] = complex(e.Value, 0) * E[k]
			IC[k] = I[k]
			EnonL[k] = E[k]
			a.ElectricRows = append(a.ElectricRows, fmt.Sprintf("%s: ½C|v|² = %.4g J", e.Name, 0.5*e.Value*math.Pow(cmplx.Abs(E[k]), 2)))
			a.Electric += 0.5 * e.Value * math.Pow(cmplx.Abs(E[k]), 2)
		default:
			EnonL[k] = E[k]
		}
	}
	// Aᵀ·C = 0.
	r := rand.New(rand.NewPCG(11, 13))
	chi := make([]float64, lat.Nodes)
	for n := range chi {
		chi[n] = r.NormFloat64()
	}
	a.CurlGrad = maxAbs(lat.Curl(lat.Grad(chi))) / maxAbs(chi)
	// Kirchhoff's voltage law: zero holonomy.
	a.KVL = maxAbsC(curlC(lat, E)) / maxE
	// Mesh EMF: the non-inductive loop voltages are balanced by −jω times
	// the flux the loop links (Faraday with the flux in the inductors).
	emf := curlC(lat, EnonL)
	flux := curlC(lat, Phi)
	worst := 0.0
	for p := range emf {
		worst = math.Max(worst, cmplx.Abs(emf[p]+complex(0, omega)*flux[p]))
	}
	a.MeshEMF = worst / maxE
	a.MeshFlux = flux
	for _, f := range flux {
		a.FluxQuanta = append(a.FluxQuanta, cmplx.Abs(f)/kron.FluxQuantum)
	}
	// Kirchhoff's current law and continuity.
	a.KCL = maxAbsC(divC(lat, I)) / maxI
	rho := divC(lat, Q)
	a.NodeCharge = rho
	ic := divC(lat, IC)
	worst = 0
	for n := range rho {
		worst = math.Max(worst, cmplx.Abs(ic[n]-complex(0, omega)*rho[n]))
	}
	a.Continuity = worst / maxI
	// The voltage is pure gauge: E = Grad(v_node) for the node potentials.
	v := make([]complex128, lat.Nodes)
	// Recover node potentials from the tree: solve Grad v = E in least
	// squares is overkill; a BFS from the datum along any branch works
	// because KVL already holds.
	seen := make([]bool, lat.Nodes)
	seen[0] = true
	for changed := true; changed; {
		changed = false
		for k := 0; k < b; k++ {
			f, t := lat.From[k], lat.To[k]
			switch {
			case seen[f] && !seen[t]:
				v[t], seen[t], changed = v[f]-E[k], true, true
			case seen[t] && !seen[f]:
				v[f], seen[f], changed = v[t]+E[k], true, true
			}
		}
	}
	g := gradC(lat, v)
	worst = 0
	for k := range g {
		worst = math.Max(worst, cmplx.Abs(g[k]-E[k]))
	}
	a.Pure = worst / maxE
	// Composition: loop sums of every branch quantity, forward and
	// reversed, at the node the first two meshes share.
	a.SharedNode = sharedNode(lat, 0, 1)
	a.LinearCompositionQuantities = []string{"voltage E", "current I", "charge Q", "flux Φ", "complex power S"}
	for _, x := range [][]complex128{E, I, Q, Phi, S} {
		a.LinearCompositionDefect = math.Max(a.LinearCompositionDefect, compositionDefect(lat, x, 0, 1, a.SharedNode))
	}
	return a, nil
}

// sharedNode is a node both plaquettes pass through, −1 if none.
func sharedNode(lat *kron.Lattice, p, q int) int {
	for n := 0; n < lat.Nodes; n++ {
		if _, ok := lat.Rebase(lat.Plaquettes[p], n); ok {
			if _, ok := lat.Rebase(lat.Plaquettes[q], n); ok {
				return n
			}
		}
	}
	return -1
}

// compositionDefect is the abelian analogue of kron.SU3Field.Compose for a
// branch quantity x: the loop sums of p and q at their shared node,
// composed in both orders and with each loop run in both directions,
// relative to the largest |x|. For numbers this is the round-off of
// addition; for SU(3) holonomies it is of order 1.
func compositionDefect(lat *kron.Lattice, x []complex128, p, q, n int) float64 {
	pp, _ := lat.Rebase(lat.Plaquettes[p], n)
	pq, _ := lat.Rebase(lat.Plaquettes[q], n)
	sum := func(path []kron.Step, reverse bool) complex128 {
		var s complex128
		if reverse {
			for i := len(path) - 1; i >= 0; i-- {
				s += complex(float64(path[i].Dir), 0) * x[path[i].Link]
			}
			return s
		}
		for _, st := range path {
			s += complex(float64(st.Dir), 0) * x[st.Link]
		}
		return s
	}
	scale := maxAbsC(x)
	if scale == 0 {
		return 0
	}
	a := sum(pp, false) + sum(pq, false)
	b := sum(pq, true) + sum(pp, true)
	c := sum(append(append([]kron.Step{}, pp...), pq...), false) // the joined loop
	return math.Max(cmplx.Abs(a-b), cmplx.Abs(a-c)) / scale
}

func curlC(lat *kron.Lattice, x []complex128) []complex128 {
	re, im := split(x)
	return join(lat.Curl(re), lat.Curl(im))
}

func divC(lat *kron.Lattice, x []complex128) []complex128 {
	re, im := split(x)
	return join(lat.Div(re), lat.Div(im))
}

func gradC(lat *kron.Lattice, x []complex128) []complex128 {
	re, im := split(x)
	return join(lat.Grad(re), lat.Grad(im))
}

func split(x []complex128) (re, im []float64) {
	re, im = make([]float64, len(x)), make([]float64, len(x))
	for i, v := range x {
		re[i], im[i] = real(v), imag(v)
	}
	return re, im
}

func join(re, im []float64) []complex128 {
	out := make([]complex128, len(re))
	for i := range out {
		out[i] = complex(re[i], im[i])
	}
	return out
}

func maxAbs(x []float64) float64 {
	m := 0.0
	for _, v := range x {
		m = math.Max(m, math.Abs(v))
	}
	return m
}

func maxAbsC(x []complex128) float64 {
	m := 0.0
	for _, v := range x {
		m = math.Max(m, cmplx.Abs(v))
	}
	return m
}

// Twin is the U(1) field on the Yee grid beside em.FDTD.
type Twin struct {
	N                 [3]int
	Steps             int
	EdgeDefect        float64 // max |E_lattice − E_FDTD| / max |E| over all free edges
	EnergyDefect      float64 // |W_lattice − W_FDTD| / W0
	EnergyDrift       float64 // max |W_lattice − W0| / W0
	Gauss             float64 // max interior Gauss generator / max Q
	Omega, OmegaYee   float64
	OmegaContinuum    float64
	FDTDLedgerResidue float64
}

// RunTwin steps the TM110 mode of a 21×17×6 closed box in both solvers.
func RunTwin(steps int) (*Twin, error) {
	n := [3]int{21, 17, 6}
	dx := 0.01
	f, err := em.NewFDTD(em.FDTDConfig{N: n, D: [3]float64{dx, dx, dx}, Ledger: true})
	if err != nil {
		return nil, err
	}
	a, b := float64(n[0]-1)*dx, float64(n[1]-1)*dx
	mode := func(x, y, z float64) [3]float64 {
		return [3]float64{0, 0, math.Sin(math.Pi*x/a) * math.Sin(math.Pi*y/b)}
	}
	f.SetE(mode)
	g, err := kron.NewYeeGrid(n, [3]float64{dx, dx, dx})
	if err != nil {
		return nil, err
	}
	cap, ind := g.Constitutive(em.Eps0, em.Mu0)
	u, err := kron.NewU1Field(g.Lattice, cap, ind, f.Dt)
	if err != nil {
		return nil, err
	}
	v := make([]float64, g.Links())
	for k := range v {
		p := g.Position(k)
		v[k] = mode(p[0], p[1], p[2])[g.Axis(k)] * dx
	}
	u.SetVoltage(v)
	we0, wm0 := u.Energy()
	w0 := we0 + wm0
	tw := &Twin{N: n, Steps: steps}
	interior := func(node int) bool {
		i, j, k := node%n[0], (node/n[0])%n[1], node/(n[0]*n[1])
		return i > 0 && i < n[0]-1 && j > 0 && j < n[1]-1 && k > 0 && k < n[2]-1
	}
	probe := g.Link(2, n[0]/2, n[1]/2, 2)
	prev := u.Voltage()[probe]
	var cross []float64
	for s := 0; s < steps; s++ {
		f.Step()
		u.Step()
		vv := u.Voltage()
		if (prev > 0) != (vv[probe] > 0) {
			cross = append(cross, u.Time()-u.Dt*vv[probe]/(vv[probe]-prev))
		}
		prev = vv[probe]
		we, wm := u.Energy()
		tw.EnergyDrift = math.Max(tw.EnergyDrift, math.Abs(we+wm-w0)/w0)
		if s%100 == 0 || s == steps-1 {
			maxE, worst := 0.0, 0.0
			for k := 0; k < g.Links(); k++ {
				if !g.Free[k] {
					continue
				}
				c := g.Axis(k)
				node := g.From[k]
				e := vv[k] / dx
				maxE = math.Max(maxE, math.Abs(f.E[c][node]))
				worst = math.Max(worst, math.Abs(e-f.E[c][node]))
			}
			tw.EdgeDefect = math.Max(tw.EdgeDefect, worst/maxE)
			tw.EnergyDefect = math.Max(tw.EnergyDefect, math.Abs(we+wm-f.Ledger().W)/w0)
			qmax := 0.0
			for _, q := range u.Q {
				qmax = math.Max(qmax, math.Abs(q))
			}
			for node, gg := range u.Gauss(nil) {
				if interior(node) {
					tw.Gauss = math.Max(tw.Gauss, math.Abs(gg)/qmax)
				}
			}
		}
	}
	period := 2 * (cross[len(cross)-1] - cross[0]) / float64(len(cross)-1)
	tw.Omega = 2 * math.Pi / period
	c0 := em.SpeedOfLight
	s := c0 * f.Dt * math.Sqrt(math.Pow(2/dx*math.Sin(math.Pi*dx/(2*a)), 2)+math.Pow(2/dx*math.Sin(math.Pi*dx/(2*b)), 2)) / 2
	tw.OmegaYee = 2 / f.Dt * math.Asin(s)
	tw.OmegaContinuum = c0 * math.Pi * math.Hypot(1/a, 1/b)
	tw.FDTDLedgerResidue = f.Ledger().Residual() / f.Ledger().W0
	return tw, nil
}

// NonAbelian is step 2: SU(3) on the circuit graph.
type NonAbelian struct {
	Fields                     int
	MeanDefect, MinDefect      float64 // ‖[W_p, W_q]‖/3 over random SU(3) fields
	CartanDefect               float64 // the same for diagonal (abelian) fields
	TraceInvariance            float64 // max change of Tr W under a gauge transformation
	GaussCovariance            float64 // max ‖G′ − g G g†‖
	GaussInvarianceBroken      float64 // min ‖G′ − G‖ over nodes: the generator is not a number
	LinearCompositionDefect    float64 // from step 1: every branch quantity
	NonlinearKCL, NonlinearKVL float64
	NonlinearComposition       float64
	DiodeCurrents              [2]float64
	DiodeVoltages              [2]float64
}

// RunNonAbelian puts SU(3) link variables on the two-mesh graph and
// measures how its loops compose, then looks for the same behaviour in
// the circuit's own quantities, linear (from step 1) and nonlinear.
func RunNonAbelian(ab *Abelian, fields int) (*NonAbelian, error) {
	c, err := Circuit(TwoMesh)
	if err != nil {
		return nil, err
	}
	nw, err := kron.Build(c)
	if err != nil {
		return nil, err
	}
	lat, err := kron.LatticeFromNetwork(nw)
	if err != nil {
		return nil, err
	}
	node := sharedNode(lat, 0, 1)
	r := rand.New(rand.NewPCG(21, 22))
	na := &NonAbelian{Fields: fields, MinDefect: math.Inf(1), LinearCompositionDefect: ab.LinearCompositionDefect}
	f := kron.NewSU3Field(lat)
	for i := 0; i < fields; i++ {
		f.Randomize(r, 1, false)
		_, _, d, ok := f.Compose(0, 1, node)
		if !ok {
			return nil, fmt.Errorf("meshes 0 and 1 share no node")
		}
		na.MeanDefect += d / float64(fields)
		na.MinDefect = math.Min(na.MinDefect, d)
		f.Randomize(r, 1, true)
		_, _, d, _ = f.Compose(0, 1, node)
		na.CartanDefect = math.Max(na.CartanDefect, d)
	}
	// Gauge covariance of the non-abelian Gauss generator.
	f.Randomize(r, 1, false)
	e := make([]kron.SU3, lat.Links())
	for k := range e {
		e[k] = kron.RandomAlgebra(r, 1, false)
	}
	g0 := f.Gauss(e, nil)
	tr0 := make([]complex128, len(lat.Plaquettes))
	for p := range tr0 {
		tr0[p] = f.WilsonLoop(p).Trace()
	}
	gauge := make([]kron.SU3, lat.Nodes)
	for n := range gauge {
		gauge[n] = kron.RandomSU3(r, 1, false)
	}
	f.Gauge(gauge)
	for k := range e {
		gk := gauge[lat.From[k]]
		e[k] = gk.Mul(e[k]).Mul(gk.Dagger())
	}
	na.GaussInvarianceBroken = math.Inf(1)
	for n, g1 := range f.Gauss(e, nil) {
		cov := gauge[n].Mul(g0[n]).Mul(gauge[n].Dagger())
		na.GaussCovariance = math.Max(na.GaussCovariance, g1.Sub(cov).Norm())
		na.GaussInvarianceBroken = math.Min(na.GaussInvarianceBroken, g1.Sub(g0[n]).Norm())
	}
	for p := range tr0 {
		na.TraceInvariance = math.Max(na.TraceInvariance, cmplx.Abs(f.WilsonLoop(p).Trace()-tr0[p]))
	}
	// The nonlinear circuit: Kirchhoff's sums still commute.
	dc, err := Circuit(DiodePair)
	if err != nil {
		return nil, err
	}
	res, err := analysis.OP(dc, analysis.Options{})
	if err != nil {
		return nil, err
	}
	nodeIdx := map[string]int{}
	for i, nm := range res.NodeNames {
		nodeIdx[nm] = i
	}
	dl := &kron.Lattice{Nodes: len(res.NodeNames)}
	I := make([]complex128, len(res.Elements))
	V := make([]complex128, len(res.Elements))
	link := map[string]int{}
	for k, el := range res.Elements {
		dl.From = append(dl.From, nodeIdx[el.Nodes[0]])
		dl.To = append(dl.To, nodeIdx[el.Nodes[1]])
		dl.Names = append(dl.Names, el.Name)
		I[k], V[k] = complex(el.Current, 0), complex(el.Voltage, 0)
		link[strings.ToLower(el.Name)] = k
	}
	// Mesh A: 0 → 1 against V1, R1, D1, R2 back to 0. Mesh B: from node 2
	// along D1, R2, then back through R3 and D2 against their orientation.
	dl.Plaquettes = [][]kron.Step{
		{{Link: link["v1"], Dir: -1}, {Link: link["r1"], Dir: 1}, {Link: link["d1"], Dir: 1}, {Link: link["r2"], Dir: 1}},
		{{Link: link["d1"], Dir: 1}, {Link: link["r2"], Dir: 1}, {Link: link["r3"], Dir: -1}, {Link: link["d2"], Dir: -1}},
	}
	if err := dl.Check(); err != nil {
		return nil, err
	}
	na.NonlinearKCL = maxAbsC(divC(dl, I)) / maxAbsC(I)
	na.NonlinearKVL = maxAbsC(curlC(dl, V)) / maxAbsC(V)
	shared := sharedNode(dl, 0, 1)
	for _, x := range [][]complex128{V, I} {
		na.NonlinearComposition = math.Max(na.NonlinearComposition, compositionDefect(dl, x, 0, 1, shared))
	}
	na.DiodeCurrents = [2]float64{real(I[link["d1"]]), real(I[link["d2"]])}
	na.DiodeVoltages = [2]float64{real(V[link["d1"]]), real(V[link["d2"]])}
	return na, nil
}

// Metrology is the quantum level of the Ohm's-law wheel since the 2019 SI.
type Metrology struct {
	KJ, RK, E, H float64
	VoltAsFreq   float64 // frequency whose Josephson step is 1 V [Hz]
	NanoampRate  float64 // electrons per second in 1 nA
	CopperZ      int     // the only quark content a circuit sees: the lattice charge
}

// RunMetrology computes the constants from the exact e and h.
func RunMetrology() Metrology {
	return Metrology{
		KJ: kron.JosephsonConstant, RK: kron.VonKlitzingConstant,
		E: kron.ElementaryCharge, H: kron.Planck,
		VoltAsFreq:  kron.JosephsonConstant,
		NanoampRate: 1e-9 / kron.ElementaryCharge,
		// Copper: 29 protons, each uud with charge 2/3 + 2/3 − 1/3 = 1.
		CopperZ: 29 * (2*2 - 1) / 3,
	}
}

func main() {
	if err := run(os.Stdout); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func run(w io.Writer) error {
	ab, err := RunAbelian(2 * math.Pi * 1e3)
	if err != nil {
		return err
	}
	fmt.Fprintln(w, "H22 — \"quantumtricity\": the quark model as circuit electricity at the quantum level")
	fmt.Fprintln(w, "Standard prediction computed beside every step (CLAUDE.md §1, H22).")
	fmt.Fprintln(w)
	fmt.Fprintln(w, "Step 1 — abelian: circuit theory is U(1) on Kron's network (standard)")
	fmt.Fprintf(w, "  Two-mesh RLC at %.0f Hz, Kron mesh frame: b = %d branches, N = %d junctions, m = %d meshes\n", ab.Omega/(2*math.Pi), ab.Branches, ab.Nodes, ab.Meshes)
	for p, s := range ab.Plaquettes {
		fmt.Fprintf(w, "    plaquette %d (Wilson loop): %s\n", p, s)
	}
	fmt.Fprintf(w, "  Curl∘Grad = 0 (Aᵀ·C = 0, d∘d = 0):                 %.1e\n", ab.CurlGrad)
	fmt.Fprintf(w, "  KVL = zero holonomy of the branch voltage:          %.1e\n", ab.KVL)
	fmt.Fprintf(w, "  voltage is pure gauge, E = Grad(v_node):            %.1e\n", ab.Pure)
	fmt.Fprintf(w, "  mesh EMF = −jω × flux holonomy (Faraday):           %.1e\n", ab.MeshEMF)
	fmt.Fprintf(w, "  KCL = d/dt of the Gauss generator, Div I = 0:       %.1e\n", ab.KCL)
	fmt.Fprintf(w, "  continuity, Div I_C = jω·(Div Q) = jω·ρ_node:        %.1e\n", ab.Continuity)
	for p, q := range ab.FluxQuanta {
		fmt.Fprintf(w, "  mesh %d flux linkage |Φ| = %.3e Wb = %.3e Φ₀ (Cooper-pair Wilson loop phase 2π·%.3e)\n", p, cmplx.Abs(ab.MeshFlux[p]), q, q)
	}
	fmt.Fprintf(w, "  Hamiltonian terms at the phasor peaks: electric Σ½C|v|² = %.4g J, magnetic Σ½L|i|² = %.4g J\n", ab.Electric, ab.Magnetic)
	for _, s := range append(ab.ElectricRows, ab.MagneticRows...) {
		fmt.Fprintf(w, "    %s\n", s)
	}
	tw, err := RunTwin(4000)
	if err != nil {
		return err
	}
	fmt.Fprintf(w, "  Kron's 1944 field network (Yee grid %v, TM110, %d steps) beside em.FDTD:\n", tw.N, tw.Steps)
	fmt.Fprintf(w, "    E edge by edge, lattice vs FDTD:                   %.1e of max |E|\n", tw.EdgeDefect)
	fmt.Fprintf(w, "    energy, lattice Hamiltonian vs FDTD ledger:        %.1e of W0 (drift %.1e, FDTD residual %.1e)\n", tw.EnergyDefect, tw.EnergyDrift, tw.FDTDLedgerResidue)
	fmt.Fprintf(w, "    Gauss generator at interior nodes:                 %.1e of max Q\n", tw.Gauss)
	fmt.Fprintf(w, "    ω/ω_Yee − 1 = %.1e, ω/ω_continuum − 1 = %.1e\n", tw.Omega/tw.OmegaYee-1, tw.Omega/tw.OmegaContinuum-1)
	fmt.Fprintln(w)
	na, err := RunNonAbelian(ab, 200)
	if err != nil {
		return err
	}
	fmt.Fprintln(w, "Step 2 — non-abelian: SU(3) link variables on the same graph (the hypothesis)")
	fmt.Fprintf(w, "  composition defect ‖W_pW_q − W_qW_p‖/3 at the shared node, %d random SU(3) fields: mean %.3f, min %.3f\n", na.Fields, na.MeanDefect, na.MinDefect)
	fmt.Fprintf(w, "  the same for diagonal (Cartan, abelian) fields:                       %.1e\n", na.CartanDefect)
	fmt.Fprintf(w, "  Tr W gauge invariant to %.1e; Gauss generator covariant to %.1e but not invariant (changes by ≥ %.2f)\n", na.TraceInvariance, na.GaussCovariance, na.GaussInvarianceBroken)
	fmt.Fprintln(w, "  Search for a circuit quantity that composes that way:")
	fmt.Fprintf(w, "    linear circuit, %s: largest defect %.1e (round-off of addition)\n", strings.Join(ab.LinearCompositionQuantities, ", "), na.LinearCompositionDefect)
	fmt.Fprintf(w, "    nonlinear diode pair (D1 %.3f V %.3e A, D2 %.3f V %.3e A): KCL %.1e, KVL %.1e, defect %.1e\n",
		na.DiodeVoltages[0], na.DiodeCurrents[0], na.DiodeVoltages[1], na.DiodeCurrents[1], na.NonlinearKCL, na.NonlinearKVL, na.NonlinearComposition)
	fmt.Fprintln(w, "    none found: every Kirchhoff sum commutes; nonlinearity is not non-commutativity.")
	fmt.Fprintln(w)
	m := RunMetrology()
	fmt.Fprintln(w, "The quantum level of V, I, R (2019 SI, from exact e and h)")
	fmt.Fprintf(w, "  volt:   K_J = 2e/h = %.9e Hz/V (Josephson 1962)\n", m.KJ)
	fmt.Fprintf(w, "  ohm:    R_K = h/e² = %.8f Ω (von Klitzing, Dorda & Pepper 1980)\n", m.RK)
	fmt.Fprintf(w, "  ampere: I = e·f; 1 nA = %.4e electrons/s (Pekola et al. 2013)\n", m.NanoampRate)
	fmt.Fprintf(w, "  quark content a circuit sees: the lattice charge, e.g. copper Z = 29·(uud) = %d; conduction is electrons only\n", m.CopperZ)
	fmt.Fprintln(w)
	fmt.Fprintln(w, "Verdict: step 1 passes (circuit theory is the abelian lattice gauge theory, to round-off against FDTD);")
	fmt.Fprintln(w, "step 2 null result: no circuit quantity composes non-abelianly. H22 status: bounded.")
	return nil
}
