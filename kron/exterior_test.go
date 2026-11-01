// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"math"
	"math/cmplx"
	"math/rand/v2"
	"path/filepath"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// linearCircuits returns every linear golden circuit plus random ones.
func linearCircuits(t *testing.T, random int) map[string]*circuit.Circuit {
	out := map[string]*circuit.Circuit{}
	files, _ := filepath.Glob(filepath.Join("..", "testdata", "netlists", "*.cir"))
	for _, f := range files {
		c, _ := load(t, f)
		if _, err := Build(c); err == nil {
			out[filepath.Base(f)] = c
		}
	}
	r := rand.New(rand.NewPCG(7, 8))
	for i := 0; i < random; i++ {
		out["random"+string(rune('A'+i%26))+string(rune('a'+i/26))] = parse(t, randomCircuit(r, 3+r.IntN(8)))
	}
	return out
}

// The §3a.1 ledger twin: the power through every open path to every
// exterior junction, from Kron's mesh and node-pair solutions, equals the
// MNA ledger row by row (element × channel) to 1e−9 of the ledger's scale;
// inactive paths carry exactly zero; the Kron ledger balances; and the net
// exterior input equals Kron's source power e·i − E·I computed in the
// orthogonal frame.
func TestExteriorLedgerTwin(t *testing.T) {
	circuits := linearCircuits(t, 120)
	worst := 0.0
	rows := 0
	for name, c := range circuits {
		mna, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			continue
		}
		n, err := Build(c)
		if err != nil {
			t.Fatal(err)
		}
		for _, form := range []Form{MeshForm, NodeForm} {
			ext, err := ExteriorDC(n, form)
			if err != nil {
				t.Fatalf("%s/%s: %v", name, form, err)
			}
			ml, kl := mna.Ledger, ext.Ledger
			if len(ml.Rows) != len(kl.Rows) || len(ext.Paths) != len(kl.Rows) {
				t.Fatalf("%s/%s: %d MNA rows, %d Kron rows, %d paths", name, form, len(ml.Rows), len(kl.Rows), len(ext.Paths))
			}
			scale := math.Max(ml.Scale, 1e-300)
			for i, mr := range ml.Rows {
				kr := kl.Rows[i]
				if mr.Element != kr.Element || mr.Channel != kr.Channel || mr.Status != kr.Status {
					t.Fatalf("%s/%s row %d: %+v against %+v", name, form, i, kr, mr)
				}
				d := math.Abs(mr.PowerW-kr.PowerW) / scale
				worst = math.Max(worst, d)
				if d > twinTol {
					t.Errorf("%s/%s: %s %v: Kron %g W, MNA %g W", name, form, kr.Element, kr.Channel, kr.PowerW, mr.PowerW)
				}
				if p := ext.Paths[i]; !p.Active && p.PowerW != 0 {
					t.Errorf("%s/%s: inactive path %s %v carries %g W", name, form, p.Element, p.Channel, p.PowerW)
				}
				rows++
			}
			for _, pair := range [][2]float64{{ml.Inflow, kl.Inflow}, {ml.Outflow, kl.Outflow},
				{ml.Dissipation, kl.Dissipation}, {ml.ModelSupplied, kl.ModelSupplied}} {
				if math.Abs(pair[0]-pair[1]) > twinTol*scale {
					t.Errorf("%s/%s: totals %v", name, form, pair)
				}
			}
			if rr := math.Abs(kl.Residual) / scale; rr > 1e-9 && math.Abs(kl.Residual) > 1e-21 {
				t.Errorf("%s/%s: Kron ledger residual %g (relative %g)", name, form, kl.Residual, rr)
			}
			iv := ext.Invariants
			for _, pr := range [][2]complex128{{iv.InputPrimitive, iv.InputFrame}, {iv.OutputPrimitive, iv.OutputFrame}} {
				if cmplx.Abs(pr[0]-pr[1]) > 1e-12*scale {
					t.Errorf("%s/%s: invariant %v in the primitive frame, %v in the orthogonal frame", name, form, pr[0], pr[1])
				}
			}
			if d := math.Abs(real(iv.SourcePower()) - (kl.Inflow - kl.Outflow)); d > twinTol*scale {
				t.Errorf("%s/%s: e·i − E·I = %g, net exterior input %g", name, form, real(iv.SourcePower()), kl.Inflow-kl.Outflow)
			}
		}
	}
	t.Logf("%d circuits, %d element × channel rows per form; worst row difference %.3g of the ledger scale", len(circuits), rows/2, worst)
}

// The AC average-power ledger from Kron's solutions equals the MNA AC
// ledger's totals at every frequency.
func TestExteriorACTwin(t *testing.T) {
	circuits := linearCircuits(t, 60)
	for name, c := range circuits {
		mna, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 2, FStart: 10, FStop: 1e7}, analysis.Options{})
		if err != nil {
			continue
		}
		n, _ := Build(c)
		for _, form := range []Form{MeshForm, NodeForm} {
			for k, f := range mna.Freq {
				ext, err := ExteriorAC(n, form, f)
				if err != nil {
					t.Fatalf("%s/%s at %g Hz: %v", name, form, f, err)
				}
				ml, kl := mna.Ledger, ext.Ledger
				scale := math.Max(kl.Scale, 1e-300)
				for _, pr := range [][2]float64{{ml.Inflow[k], kl.Inflow}, {ml.Outflow[k], kl.Outflow},
					{ml.Dissipation[k], kl.Dissipation}, {ml.ModelSupplied[k], kl.ModelSupplied}} {
					if math.Abs(pr[0]-pr[1]) > twinTol*scale {
						t.Fatalf("%s/%s at %g Hz: MNA %g W, Kron %g W (scale %g)", name, form, f, pr[0], pr[1], scale)
					}
				}
			}
		}
	}
}

// The orthogonal frame is square and non-singular, counts m = b − N + S
// closed and p = N − S open paths also for several sub-networks, and its
// two families are orthogonal: currents that satisfy KCL have no
// open-path coordinates, voltages that satisfy KVL no closed-path ones.
func TestOrthogonalFrame(t *testing.T) {
	r := rand.New(rand.NewPCG(9, 10))
	for trial := 0; trial < 200; trial++ {
		n := 2 + r.IntN(10)
		b := r.IntN(20)
		from, to := make([]int, b), make([]int, b)
		for k := range from {
			from[k] = r.IntN(n)
			to[k] = r.IntN(n)
			for to[k] == from[k] {
				to[k] = r.IntN(n)
			}
		}
		nw := &Network{g: newGraph(n, from, to)}
		f, err := nw.Frame()
		if err != nil {
			t.Fatal(err)
		}
		_, m := f.Closed.Dims()
		_, p := f.Open.Dims()
		S := nw.g.subNetworks()
		if m != b-n+S || p != n-S {
			t.Fatalf("trial %d: m = %d, p = %d for b = %d, N = %d, S = %d", trial, m, p, b, n, S)
		}
		PC, _ := TMul(f.Open, f.Closed)
		for _, v := range PC.a {
			if v != 0 {
				t.Fatalf("trial %d: Pᵀ·C ≠ 0", trial)
			}
		}
		// A KCL current (a random mesh combination) and a KVL voltage (a
		// random potential difference).
		ip := NewVec[Mesh](m)
		for i := range ip.v {
			ip.v[i] = complex(r.NormFloat64(), r.NormFloat64())
		}
		vp := NewVec[Node](p)
		for i := range vp.v {
			vp.v[i] = complex(r.NormFloat64(), 0)
		}
		i, _ := MulVec(f.Closed, ip)
		e, _ := MulVec(f.Open, vp)
		if b == 0 {
			continue
		}
		ic, io, err := f.Coordinates(i)
		if err != nil {
			t.Fatalf("trial %d: frame singular: %v", trial, err)
		}
		ec, eo, _ := f.Coordinates(e)
		scale := math.Max(maxAbs(i.v), maxAbs(e.v)) + 1
		if maxAbs(io.v) > 1e-12*scale || maxAbs(ec.v) > 1e-12*scale {
			t.Fatalf("trial %d: KCL current has open part %g, KVL voltage closed part %g", trial, maxAbs(io.v), maxAbs(ec.v))
		}
		for k := range ic.v {
			if cmplx.Abs(ic.v[k]-ip.v[k]) > 1e-12*scale {
				t.Fatalf("trial %d: closed coordinates not recovered", trial)
			}
		}
		for k := range eo.v {
			if cmplx.Abs(eo.v[k]-vp.v[k]) > 1e-12*scale {
				t.Fatalf("trial %d: open coordinates not recovered", trial)
			}
		}
	}
}

// A TAN chord z = jωM between two inductors is exactly the MNA K element;
// a chord with a propagation delay matches the analytic two-loop
// solution; chords cannot attach to current-source branches.
func TestChords(t *testing.T) {
	src := "two loops\nV1 a 0 DC 0 AC 1\nR1 a p 50\nL1 p 0 10u\nL2 s 0 40u\nR2 s 0 200\n"
	k, L1, L2 := 0.7, 10e-6, 40e-6
	M := k * math.Sqrt(L1*L2)
	withK := parse(t, src+"K1 L1 L2 0.7\n.end\n")
	plain := parse(t, src+".end\n")
	freqs := []float64{1e3, 1e5, 1e6, 3e6}
	mna, err := analysis.AC(withK, analysis.ACSpec{Variation: "lin", Points: 1, FStart: 1e6, FStop: 1e6}, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	for _, form := range []Form{MeshForm, NodeForm} {
		n, err := Build(plain)
		if err != nil {
			t.Fatal(err)
		}
		if err := n.AddChord("L1", "L2", func(w float64) complex128 { return complex(0, w*M) }, true); err != nil {
			t.Fatal(err)
		}
		sol, err := n.Solve(form, 2*math.Pi*1e6, false)
		if err != nil {
			t.Fatal(err)
		}
		v, _ := n.voltages(sol)
		for node := 1; node < len(mna.V); node++ {
			if d := cmplx.Abs(v[node] - mna.V[node][0]); d > twinTol*cmplx.Abs(mna.V[1][0]) {
				t.Errorf("%s: chord V[%s] = %v, MNA with K %v", form, mna.NodeNames[node], v[node], mna.V[node][0])
			}
		}
		// Delayed chord: the secondary's emf is z·i1 with z = jωM·e^(−jωτ).
		const tau = 50e-9
		nd, _ := Build(plain)
		_ = nd.AddChord("L2", "L1", func(w float64) complex128 { return complex(0, w*M) * cmplx.Exp(complex(0, -w*tau)) }, false)
		for _, f := range freqs {
			w := 2 * math.Pi * f
			sol, err := nd.Solve(form, w, false)
			if err != nil {
				t.Fatal(err)
			}
			vv, _ := nd.voltages(sol)
			// Primary is unaffected by the one-way chord; the secondary
			// loop gives (jωL2 ∥ R2) driven by the emf: v_s = emf·R2/(R2 + jωL2).
			i1 := 1 / complex(50, w*L1)
			emf := complex(0, w*M) * cmplx.Exp(complex(0, -w*tau)) * i1
			want := emf * 200 / complex(200, w*L2)
			si := 0
			for i, nn := range plain.NodeNames {
				if nn == "s" {
					si = i
				}
			}
			if d := cmplx.Abs(vv[si] - want); d > 1e-12*cmplx.Abs(want) {
				t.Errorf("%s at %g Hz: delayed chord v(s) = %v, analytic %v", form, f, vv[si], want)
			}
		}
	}
	n, _ := Build(parse(t, "i\nI1 0 a 1m\nR1 a 0 1k\nR2 a 0 2k\n.end\n"))
	if err := n.AddChord("I1", "R1", func(float64) complex128 { return 1 }, false); !errors.Is(err, ErrUnsupported) {
		t.Errorf("chord on a current source: %v", err)
	}
}

// Radiation-resistance edge: a small loop (inductor plus wire resistance)
// fed by a current source radiates ½·R_rad·|I|² with R_rad from Balanis;
// the edge's power leaves through the loop's light/radiation open path,
// the ledger balances, and R_rad agrees with the circumference form
// 20π²·(C/λ)⁴ for a circular loop.
func TestRadiationEdge(t *testing.T) {
	a := 0.05 // loop radius [m]
	area, circ := math.Pi*a*a, 2*math.Pi*a
	rr := SmallLoopRadiation(area, 1)
	const c0 = 299792458.0
	for _, f := range []float64{1e6, 1e7, 3e7} {
		lam := c0 / f
		want := 20 * math.Pi * math.Pi * math.Pow(circ/lam, 4)
		if got := rr(2 * math.Pi * f); math.Abs(got-want) > 1e-12*want {
			t.Fatalf("R_rad at %g Hz: %g Ω, circumference form %g Ω", f, got, want)
		}
	}
	c := parse(t, "loop\nI1 0 a DC 0 AC 0.1\nL1 a w 1u\nRw w 0 0.2\nC1 a 0 100p\n.end\n")
	for _, form := range []Form{MeshForm, NodeForm} {
		n, _ := Build(c)
		if err := n.AddRadiation("L1", rr); err != nil {
			t.Fatal(err)
		}
		f := 1.5e7
		ext, err := ExteriorAC(n, form, f)
		if err != nil {
			t.Fatal(err)
		}
		sol, _ := n.Solve(form, 2*math.Pi*f, false)
		iL := sol.I.At(n.elemBranch[c.Index("L1")])
		want := 0.5 * rr(2*math.Pi*f) * cmplx.Abs(iL) * cmplx.Abs(iL)
		var radiated float64
		for _, p := range ext.Paths {
			if p.Element == "L1" && p.Channel == exchange.Light {
				radiated = -p.PowerW
				if !p.Active {
					t.Errorf("%s: radiation path inactive", form)
				}
			}
		}
		if math.Abs(radiated-want) > 1e-12*want || want == 0 {
			t.Errorf("%s: radiated %g W, ½·R_rad·|I|² = %g W", form, radiated, want)
		}
		l := ext.Ledger
		if math.Abs(l.Residual) > 1e-12*l.Scale || math.Abs(l.Outflow-want) > 1e-12*want {
			t.Errorf("%s: ledger inflow %g, dissipation %g, outflow %g, residual %g", form, l.Inflow, l.Dissipation, l.Outflow, l.Residual)
		}
		t.Logf("%s: %.4g W radiated of %.4g W in (R_rad = %.3g Ω)", form, radiated, l.Inflow, rr(2*math.Pi*f))
		if _, err := n.Solve(form, 0, true); err != nil {
			t.Errorf("%s: DC with a radiation edge: %v", form, err)
		}
	}
}

// Every element has one open path to every channel's exterior junction;
// the electrical junction receives the sources' power and the thermal
// junction the resistors' heat, and the junction totals sum to the
// ledger's net exchange.
func TestExteriorJunctions(t *testing.T) {
	c := parse(t, "div\nV1 a 0 10\nR1 a b 1k\nR2 b 0 4k\nC1 b 0 1n\n.end\n")
	n, _ := Build(c)
	ext, err := ExteriorDC(n, MeshForm)
	if err != nil {
		t.Fatal(err)
	}
	total := 0.0
	for ch := exchange.ChannelID(0); int(ch) < exchange.NumChannels; ch++ {
		paths, net := ext.Junction(ch)
		if len(paths) != len(c.Elements) {
			t.Fatalf("channel %v: %d open paths for %d elements", ch, len(paths), len(c.Elements))
		}
		total += net
		switch ch {
		case exchange.Electrical:
			if math.Abs(net-0.02) > 1e-15 {
				t.Errorf("electrical junction %g W, want 20 mW in", net)
			}
		case exchange.Thermal:
			if math.Abs(net+0.02) > 1e-15 {
				t.Errorf("thermal junction %g W, want 20 mW out", net)
			}
		default:
			if net != 0 {
				t.Errorf("channel %v carries %g W", ch, net)
			}
		}
	}
	if math.Abs(total) > 1e-15 {
		t.Errorf("junction totals sum to %g W", total)
	}
}
