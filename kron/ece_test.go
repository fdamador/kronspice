// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"errors"
	"fmt"
	"math"
	"math/rand/v2"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Tests of the hypothetical ECE / ECE2 spin-connection terms in Kron's
// frame (ece.go; circuit/ece.go; register entry H48). Every result is a
// twin of the MNA engine, which stamps the same terms from the device
// models (CLAUDE.md §3a rule 4: 1e-9 relative).

// eceRandom is randomCircuit with spin-connection terms on every R, C, L
// and K: torsion on the coils and couplings, scalar terms on every
// non-ground terminal, defaults from .ece and per-element overrides.
func eceRandom(r *rand.Rand, form string) string {
	src := randomCircuit(r, 3+r.IntN(8))
	lines := strings.Split(src, "\n")
	w0 := 1e-6 * (0.2 + r.Float64())
	head := fmt.Sprintf(".ece form=%s w0=%g kappa=%g", form, w0, 0.01+0.05*r.Float64())
	if r.IntN(2) == 0 { // UFT 456's terms on every capacitor (k0) and coil and K (kl0)
		head += fmt.Sprintf(" k0=%g kl0=%g", 1e-3*(r.Float64()-0.3), 1e-6*(r.Float64()-0.5))
	}
	out := []string{lines[0], head}
	for _, l := range lines[1:] {
		f := strings.Fields(l)
		if len(f) >= 4 && (l[0] == 'R' || l[0] == 'C' || l[0] == 'L') && r.IntN(3) == 0 {
			l += fmt.Sprintf(" kappa_a=%g kappa_b=%g", 0.1*(r.Float64()-0.5), 0.1*(r.Float64()-0.5))
		}
		if len(f) >= 4 && l[0] == 'L' && r.IntN(2) == 0 {
			l += fmt.Sprintf(" w0=%g", 2*w0*r.Float64())
		}
		out = append(out, l)
	}
	return strings.Join(out, "\n")
}

// Mesh and node forms equal MNA with the terms on, OP and AC, both forms
// of the field relation.
func TestECETwinRandom(t *testing.T) {
	r := rand.New(rand.NewPCG(7, 11))
	checked, worst := 0, 0.0
	for trial := 0; trial < 200; trial++ {
		form := []string{"ece", "ece2"}[trial%2]
		src := eceRandom(r, form)
		c := parse(t, src)
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			continue // a random controlled-source loop with no solution
		}
		var mnaI []float64
		for _, e := range op.Elements {
			mnaI = append(mnaI, e.Current)
		}
		for _, f := range []Form{MeshForm, NodeForm} {
			k, err := OP(c, f)
			if err != nil {
				t.Fatalf("trial %d %s OP: %v\n%s", trial, f, err, src)
			}
			var kI []float64
			for _, e := range k.Elements {
				kI = append(kI, e.Current)
			}
			w1, _ := worstRel(real2c(k.NodeVoltages), real2c(op.NodeVoltages))
			w2, at := worstRel(real2c(kI), real2c(mnaI))
			worst = math.Max(worst, math.Max(w1, w2))
			if w1 > twinTol || w2 > twinTol {
				t.Fatalf("trial %d %s OP: voltages %.3g, currents %.3g (%s)\n%s", trial, f, w1, w2, op.Elements[at].Name, src)
			}
			for _, fr := range []float64{10, 1e4, 1e7} {
				m, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: fr, FStop: fr}, analysis.Options{})
				if err != nil {
					continue
				}
				ka, err := AC(c, m.Freq, f)
				if err != nil {
					t.Fatalf("trial %d %s AC %g: %v\n%s", trial, f, fr, err, src)
				}
				var got, want []complex128
				for n := 1; n < len(m.V); n++ {
					got, want = append(got, ka.V[n][0]), append(want, m.V[n][0])
				}
				w, _ := worstRel(got, want)
				worst = math.Max(worst, w)
				if w > twinTol {
					t.Fatalf("trial %d %s AC %g Hz: voltages %.3g\n%s", trial, f, fr, w, src)
				}
			}
		}
		checked++
	}
	if checked < 150 {
		t.Fatalf("only %d of 200 random ECE circuits solved", checked)
	}
	t.Logf("%d random ECE/ECE2 circuits, mesh and node forms, OP and AC: worst relative difference to MNA %.3g", checked, worst)
}

// The ledger built in Kron's frame (open paths, model-supplied power of the
// terms) equals the MNA ledger row by row, DC and AC, and both close.
func TestECELedgerTwin(t *testing.T) {
	r := rand.New(rand.NewPCG(3, 5))
	checked := 0
	for trial := 0; trial < 120; trial++ {
		src := eceRandom(r, []string{"ece", "ece2"}[trial%2])
		c := parse(t, src)
		mna, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			continue
		}
		mac, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 1, FStart: 10, FStop: 1e6}, analysis.Options{})
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
				t.Fatalf("trial %d %s: %v", trial, form, err)
			}
			ml, kl := mna.Ledger, ext.Ledger
			scale := math.Max(ml.Scale, 1e-300)
			for i, mr := range ml.Rows {
				if d := math.Abs(mr.PowerW-kl.Rows[i].PowerW) / scale; d > twinTol {
					t.Fatalf("trial %d %s: row %s %v: Kron %g W, MNA %g W\n%s", trial, form, mr.Element, mr.Channel, kl.Rows[i].PowerW, mr.PowerW, src)
				}
			}
			for _, p := range [][2]float64{{ml.ModelSupplied, kl.ModelSupplied}, {ml.StorageRate, kl.StorageRate},
				{ml.Dissipation, kl.Dissipation}, {ml.Inflow, kl.Inflow}} {
				if math.Abs(p[0]-p[1]) > twinTol*scale {
					t.Fatalf("trial %d %s: totals MNA %g, Kron %g\n%s", trial, form, p[0], p[1], src)
				}
			}
			if math.Abs(kl.Residual) > 1e-9*scale && math.Abs(kl.Residual) > 1e-21 {
				t.Fatalf("trial %d %s: Kron DC ledger residual %g of %g", trial, form, kl.Residual, scale)
			}
			for k, f := range mac.Freq {
				ea, err := ExteriorAC(n, form, f)
				if err != nil {
					t.Fatalf("trial %d %s AC: %v", trial, form, err)
				}
				al := ea.Ledger
				sc := math.Max(al.Scale, 1e-300)
				for _, p := range [][2]float64{{mac.Ledger.ModelSupplied[k], al.ModelSupplied},
					{mac.Ledger.Dissipation[k], al.Dissipation}, {mac.Ledger.Inflow[k], al.Inflow}} {
					if math.Abs(p[0]-p[1]) > twinTol*sc {
						t.Fatalf("trial %d %s AC %g Hz: MNA %g, Kron %g\n%s", trial, form, f, p[0], p[1], src)
					}
				}
			}
		}
		if mna.Ledger.ModelSupplied == 0 {
			t.Fatalf("trial %d: the terms booked no model-supplied power\n%s", trial, src)
		}
		checked++
	}
	if checked < 80 {
		t.Fatalf("only %d circuits checked", checked)
	}
}

const eceTorn = `ece torn
.ece form=ece2 w0=5e-7 kappa=0.03
V1 in 0 DC 2 AC 1 SIN(2 1 1k)
R1 in a 10 kappa_a=0.02 kappa_b=-0.01
L1 a b 10m
C1 b 0 10u
R3 b 0 50
L3 a d 2m
R4 d 0 40
K1 L1 L3 0.4 w0=8e-7
LX b c 1m kappa=0
RX b c 30 kappa=0
L2 c 0 5m w0=1e-6
R2 c 0 20
C2 c 0 4u
.end
`

// Diakoptics with the terms: a coil with torsion is torn as the admittance
// 1/(jωL + r_t), also at DC; tearing an element with a scalar term is
// refused with the reason; the torn Newton Jacobian tears the device stamps
// and needs neither.
func TestECETearing(t *testing.T) {
	c := parse(t, eceTorn)
	ref, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	// LX and RX in parallel join the two halves; tear both.
	d, err := Tear(c, Tearing{Torn: []string{"LX", "RX"}}, nil)
	if err != nil {
		t.Fatal(err)
	}
	for _, dc := range []bool{true, false} {
		w := 2 * math.Pi * 1e3
		if dc {
			w = 0
		}
		sol, st, err := d.Solve(w, dc)
		if err != nil {
			t.Fatalf("dc=%v: %v", dc, err)
		}
		if st.Pieces < 2 {
			t.Errorf("dc=%v: %d pieces", dc, st.Pieces)
		}
		un, err := d.n.Solve(NodeForm, w, dc)
		if err != nil {
			t.Fatal(err)
		}
		for k := 0; k < len(d.n.branches); k++ {
			a, b := sol.E.At(k), un.E.At(k)
			if d := absC(a - b); d > twinTol*math.Max(absC(b), 1) {
				t.Fatalf("dc=%v branch %d: torn %v, untorn %v", dc, k, a, b)
			}
		}
		if dc {
			r, err := d.n.opResult(sol)
			if err != nil {
				t.Fatal(err)
			}
			if w, _ := worstRel(real2c(r.NodeVoltages), real2c(ref.NodeVoltages)); w > twinTol {
				t.Fatalf("torn DC vs MNA: %.3g", w)
			}
		}
	}
	if _, err := Tear(c, Tearing{Torn: []string{"R1"}}, nil); !errors.Is(err, ErrUnsupported) || !strings.Contains(err.Error(), "scalar term") {
		t.Fatalf("tearing R1 (scalar term): %v", err)
	}
	// The torn Jacobian (it cannot cut an inductor: its coupling runs through
	// a branch current): the same circuit with the coil bridge LX replaced by
	// a capacitor, torn at the bridge, OP and transient.
	c = parse(t, strings.Replace(eceTorn, "LX b c 1m kappa=0", "CX b c 2u kappa=0", 1))
	if ref, err = analysis.OP(c, analysis.Options{}); err != nil {
		t.Fatal(err)
	}
	jt, err := NewJacobianTearing(c, Tearing{Torn: []string{"RX", "CX"}})
	if err != nil {
		t.Fatal(err)
	}
	tor, err := analysis.OP(c, analysis.Options{Solver: jt})
	if err != nil {
		t.Fatal(err)
	}
	if w, _ := worstRel(real2c(tor.NodeVoltages), real2c(ref.NodeVoltages)); w > twinTol {
		t.Fatalf("torn Jacobian OP vs MNA: %.3g", w)
	}
	ts := analysis.TranSpec{Step: 1e-5, Stop: 2e-3, MaxStep: 2e-6}
	wu, err := analysis.Tran(c, ts, analysis.TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	wt, err := analysis.Tran(c, ts, analysis.TranOptions{Options: analysis.Options{Solver: jt}})
	if err != nil {
		t.Fatal(err)
	}
	if len(wu.Time) != len(wt.Time) {
		t.Fatalf("transient points %d vs %d", len(wu.Time), len(wt.Time))
	}
	worst, peak := 0.0, 0.0
	for n := 1; n < len(wu.V); n++ {
		for k := range wu.Time {
			worst = math.Max(worst, math.Abs(wu.V[n][k]-wt.V[n][k]))
			peak = math.Max(peak, math.Abs(wu.V[n][k]))
		}
	}
	if worst > 1e-9*peak {
		t.Fatalf("torn Jacobian transient differs by %.3g of %.3g", worst, peak)
	}
	if s := jt.Stats(); s.Fallbacks != 0 || s.MaxPieces < 2 {
		t.Errorf("torn Jacobian stats %+v", s)
	}
}

// Without .ece, or with every term zero, the primitive network is the
// Maxwell one: no control coils, no torsion, the same solution bit for bit.
func TestECEZeroIsPlainKron(t *testing.T) {
	zero := strings.NewReplacer(".ece form=ece2 w0=5e-7 kappa=0.03", ".ece form=ece2 w0=0 kappa=0",
		" kappa_a=0.02 kappa_b=-0.01", "", " w0=1e-6", "", " w0=8e-7", "").Replace(eceTorn)
	bare := strings.NewReplacer(".ece form=ece2 w0=5e-7 kappa=0.03\n", "",
		" kappa_a=0.02 kappa_b=-0.01", "", " w0=1e-6", "", " w0=8e-7", "", " kappa=0", "").Replace(eceTorn)
	nz, err := Build(parse(t, zero))
	if err != nil {
		t.Fatal(err)
	}
	nb, err := Build(parse(t, bare))
	if err != nil {
		t.Fatal(err)
	}
	if len(nz.branches) != len(nb.branches) || len(nz.eceK) != 0 {
		t.Fatalf("zero terms: %d branches (plain %d), %d mutual torsion entries", len(nz.branches), len(nb.branches), len(nz.eceK))
	}
	for _, f := range []Form{MeshForm, NodeForm} {
		a, err := nz.Solve(f, 2*math.Pi*1e3, false)
		if err != nil {
			t.Fatal(err)
		}
		b, err := nb.Solve(f, 2*math.Pi*1e3, false)
		if err != nil {
			t.Fatal(err)
		}
		for k := range nb.branches {
			if a.E.At(k) != b.E.At(k) || a.I.At(k) != b.I.At(k) {
				t.Fatalf("%s branch %d differs", f, k)
			}
		}
	}
}

// Precheck: a loop of inductors and a voltage source is a zero-impedance
// DC loop, unless torsion gives the coil a resistance.
func TestECEPrecheck(t *testing.T) {
	loop := "loop\nV1 a 0 1\nL1 a 0 1m\nR1 a 0 1\n.end\n"
	if _, err := Precheck(parse(t, loop), PrecheckOptions{DC: true}); !errors.Is(err, ErrSingular) {
		t.Fatalf("plain loop: %v", err)
	}
	ece := strings.Replace(loop, "loop\n", "loop\n.ece form=ece2 w0=1e-6\n", 1)
	if _, err := Precheck(parse(t, ece), PrecheckOptions{DC: true}); err != nil {
		t.Fatalf("loop with torsion: %v", err)
	}
	_ = circuit.ECE2
}

func absC(z complex128) float64 { return math.Hypot(real(z), imag(z)) }

// UFT 456's shunt makes a capacitor conduct at DC: the linear torn solve
// cuts it with admittance g_e (DC) and g_e + jωC (AC), equal to the untorn
// Kron solution and to MNA; kl0 on the coils rides along in r_t.
func TestECEUFT456Torn(t *testing.T) {
	src := strings.NewReplacer(".ece form=ece2 w0=5e-7 kappa=0.03", ".ece form=ece2 w0=5e-7 kappa=0.03 kl0=2e-7",
		"LX b c 1m kappa=0", "CX b c 2u kappa=0 k0=0.05").Replace(eceTorn)
	c := parse(t, src)
	ref, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		t.Fatal(err)
	}
	d, err := Tear(c, Tearing{Torn: []string{"CX", "RX"}}, nil)
	if err != nil {
		t.Fatal(err)
	}
	for _, dc := range []bool{true, false} {
		w := 2 * math.Pi * 1e3
		if dc {
			w = 0
		}
		sol, st, err := d.Solve(w, dc)
		if err != nil {
			t.Fatalf("dc=%v: %v", dc, err)
		}
		if st.Pieces < 2 {
			t.Errorf("dc=%v: %d pieces", dc, st.Pieces)
		}
		un, err := d.n.Solve(MeshForm, w, dc)
		if err != nil {
			t.Fatal(err)
		}
		for k := range d.n.branches {
			if dd := absC(sol.E.At(k) - un.E.At(k)); dd > twinTol*math.Max(absC(un.E.At(k)), 1) {
				t.Fatalf("dc=%v branch %d: torn %v, untorn %v", dc, k, sol.E.At(k), un.E.At(k))
			}
		}
		if dc {
			r, err := d.n.opResult(sol)
			if err != nil {
				t.Fatal(err)
			}
			if wr, _ := worstRel(real2c(r.NodeVoltages), real2c(ref.NodeVoltages)); wr > twinTol {
				t.Fatalf("torn DC vs MNA: %.3g", wr)
			}
		}
	}
}

// The skew mutual torsion (K's w0skew, circuit.ECETerms.W0Skew) is a
// non-reciprocal pair of chords in Kron's frame: mesh and node forms equal
// MNA (OP and AC), the ledgers agree row by row, and the skew term alone
// delivers no power (Kron's workless F, Tensors for Circuits eq. 31.5).
func TestECESkewTwinRandom(t *testing.T) {
	r := rand.New(rand.NewPCG(13, 17))
	checked := 0
	for trial := 0; trial < 150; trial++ {
		form := []string{"ece", "ece2"}[trial%2]
		src := randomCircuit(r, 3+r.IntN(8))
		if !strings.Contains(src, "\nK1 ") {
			continue
		}
		skew := 1e-6 * (r.Float64() - 0.5)
		src = strings.Replace(src, "\nK1 ", "\n.ece form="+form+"\nK1 ", 1)
		lines := strings.Split(src, "\n")
		for i, l := range lines {
			if strings.HasPrefix(l, "K1 ") {
				lines[i] = l + fmt.Sprintf(" w0skew=%g", skew)
			}
		}
		src = strings.Join(lines, "\n")
		c := parse(t, src)
		op, err := analysis.OP(c, analysis.Options{})
		if err != nil {
			continue
		}
		// Random circuits carry their own model-supplied elements (negative
		// R, controlled sources); the skew term must add nothing to them.
		for _, row := range op.Ledger.Supplied {
			if row.Element == "K1" && row.PowerW != 0 {
				t.Fatalf("trial %d: K1 booked %g W model-supplied\n%s", trial, row.PowerW, src)
			}
		}
		n, err := Build(c)
		if err != nil {
			t.Fatal(err)
		}
		for _, f := range []Form{MeshForm, NodeForm} {
			k, err := OP(c, f)
			if err != nil {
				t.Fatalf("trial %d %s OP: %v\n%s", trial, f, err, src)
			}
			if w, _ := worstRel(real2c(k.NodeVoltages), real2c(op.NodeVoltages)); w > twinTol {
				t.Fatalf("trial %d %s OP: voltages %.3g\n%s", trial, f, w, src)
			}
			ext, err := ExteriorDC(n, f)
			if err != nil {
				t.Fatalf("trial %d %s: %v", trial, f, err)
			}
			scale := math.Max(op.Ledger.Scale, 1e-300)
			for i, mr := range op.Ledger.Rows {
				if d := math.Abs(mr.PowerW-ext.Ledger.Rows[i].PowerW) / scale; d > twinTol {
					t.Fatalf("trial %d %s: row %s: Kron %g W, MNA %g W\n%s", trial, f, mr.Element, ext.Ledger.Rows[i].PowerW, mr.PowerW, src)
				}
			}
			for _, fr := range []float64{10, 1e4, 1e7} {
				m, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: fr, FStop: fr}, analysis.Options{})
				if err != nil {
					continue
				}
				ka, err := AC(c, m.Freq, f)
				if err != nil {
					t.Fatalf("trial %d %s AC %g: %v\n%s", trial, f, fr, err, src)
				}
				var got, want []complex128
				for nn := 1; nn < len(m.V); nn++ {
					got, want = append(got, ka.V[nn][0]), append(want, m.V[nn][0])
				}
				if w, _ := worstRel(got, want); w > twinTol {
					t.Fatalf("trial %d %s AC %g Hz: voltages %.3g\n%s", trial, f, fr, w, src)
				}
			}
		}
		checked++
	}
	if checked < 40 {
		t.Fatalf("only %d random skew circuits checked", checked)
	}
	t.Logf("%d random circuits with skew mutual torsion: Kron mesh and node forms equal MNA, no model-supplied power", checked)
}
