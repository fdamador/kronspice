// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// Tests of the hypothetical ECE / ECE2 spin-connection terms on R, L, C and
// K (circuit/ece.go, device/ece.go; register entry H48).

const cLight = circuit.SpeedOfLight

// eceDeck is a coupled RLC network with every kind of term: scalar κ on
// every non-ground terminal (defaults and per-element overrides), torsion
// on both coils (default and override ω₀) and mutual torsion on K1.
func eceDeck(form string) string {
	return `ece network
.ece form=` + form + ` w0=5e-7 kappa=0.03
V1 in 0 DC 2 AC 1 SIN(2 1 1k)
R1 in a 10 kappa_a=0.02 kappa_b=-0.01
L1 a b 10m
C1 b 0 10u
R3 b 0 50
L2 c 0 5m w0=1e-6
R2 c 0 20
K1 L1 L2 0.5 w0=8e-7
.end
`
}

// eceTwin is the same network built from ordinary elements only, as the
// AIAS netlist expander does: a series VCVS for each terminal's scalar EMF
// (the body sees (1 + s·κ)·φ), a series resistor s·c·ω₀·L for torsion, and
// a CCVS s·c·ω₀·M driven by the other coil's current for mutual torsion.
func eceTwin(s float64) string {
	g := func(x float64) string { return fmt.Sprintf("%.17g", x) }
	rt1 := s * cLight * 5e-7 * 10e-3
	rt2 := s * cLight * 1e-6 * 5e-3
	rm := s * cLight * 8e-7 * 0.5 * math.Sqrt(10e-3*5e-3)
	return `ece twin
V1 in 0 DC 2 AC 1 SIN(2 1 1k)
Ea_R1 in x1 in 0 ` + g(-s*0.02) + `
R1 x1 x2 10
Eb_R1 x2 a a 0 ` + g(s*-0.01) + `
Ea_L1 a y1 a 0 ` + g(-s*0.03) + `
Vs1 y1 y2 0
L1 y2 y3 10m
Rt1 y3 y4 ` + g(rt1) + `
Ht1 y4 y5 Vs2 ` + g(rm) + `
Eb_L1 y5 b b 0 ` + g(s*0.03) + `
Ea_C1 b z1 b 0 ` + g(-s*0.03) + `
C1 z1 0 10u
Ea_R3 b w1 b 0 ` + g(-s*0.03) + `
R3 w1 0 50
Ea_L2 c u1 c 0 ` + g(-s*0.03) + `
Vs2 u1 u2 0
L2 u2 u3 5m
Rt2 u3 u4 ` + g(rt2) + `
Ht2 u4 0 Vs1 ` + g(rm) + `
Ea_R2 c v1 c 0 ` + g(-s*0.03) + `
R2 v1 0 20
K1 L1 L2 0.5
.end
`
}

func eceCircuit(t *testing.T, src string) *circuit.Circuit {
	t.Helper()
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func eceRel(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(math.Max(math.Abs(a), math.Abs(b)), 1e-12)
}

// interp is y(t) by linear interpolation on (ts, ys).
func interp(ts, ys []float64, t float64) float64 {
	k := 1
	for k < len(ts)-1 && ts[k] < t {
		k++
	}
	f := (t - ts[k-1]) / (ts[k] - ts[k-1])
	return ys[k-1] + f*(ys[k]-ys[k-1])
}

// The engine's terms equal the same terms built from ordinary elements, in
// OP, AC and transient, for both forms; and the ledger closes with the
// terms' power in the model-supplied column, labelled hypothetical.
func TestECETwin(t *testing.T) {
	for _, tc := range []struct {
		form string
		s    float64
	}{{"ece", 1}, {"ece2", -2}} {
		t.Run(tc.form, func(t *testing.T) {
			ce, tw := eceCircuit(t, eceDeck(tc.form)), eceCircuit(t, eceTwin(tc.s))

			// Operating point.
			re, err := OP(ce, Options{})
			if err != nil {
				t.Fatal(err)
			}
			rt, err := OP(tw, Options{})
			if err != nil {
				t.Fatal(err)
			}
			for _, n := range []string{"in", "a", "b", "c"} {
				ve, _ := re.Voltage(n)
				vt, _ := rt.Voltage(n)
				if d := eceRel(ve, vt); d > 1e-9 {
					t.Errorf("OP V(%s) = %.12g, twin %.12g (rel %.2g)", n, ve, vt, d)
				}
			}
			for _, n := range []string{"L1", "L2", "R1"} {
				ee, _ := re.Element(n)
				et, _ := rt.Element(n)
				if d := eceRel(ee.Current, et.Current); d > 1e-9 {
					t.Errorf("OP I(%s) = %.12g, twin %.12g", n, ee.Current, et.Current)
				}
			}
			if rr := re.Ledger.RelativeResidual(); rr > 1e-9 {
				t.Errorf("OP ledger residual %.3g", rr)
			}
			if re.Ledger.ModelSupplied == 0 {
				t.Error("OP: the ECE terms booked no model-supplied power")
			}
			label := strings.ToUpper(tc.form) + " spin connection (hypothetical, H48)"
			for _, row := range re.Ledger.Supplied {
				if !strings.HasPrefix(row.Note, label) {
					t.Errorf("OP supplied row %s note %q, want prefix %q", row.Element, row.Note, label)
				}
			}

			// AC.
			spec := ACSpec{Variation: "dec", Points: 5, FStart: 10, FStop: 1e5}
			ae, err := AC(ce, spec, Options{})
			if err != nil {
				t.Fatal(err)
			}
			at, err := AC(tw, spec, Options{})
			if err != nil {
				t.Fatal(err)
			}
			for _, n := range []string{"a", "b", "c"} {
				xe, _ := ae.NodeVoltage(n)
				xt, _ := at.NodeVoltage(n)
				for k := range xe {
					if d := cmplx.Abs(xe[k]-xt[k]) / math.Max(cmplx.Abs(xt[k]), 1e-12); d > 1e-9 {
						t.Errorf("AC V(%s) at %g Hz = %v, twin %v", n, ae.Freq[k], xe[k], xt[k])
					}
				}
			}
			for k := range ae.Freq {
				l := ae.Ledger
				scale := l.Inflow[k] + l.Outflow[k] + l.Dissipation[k] + math.Abs(l.ModelSupplied[k])
				if math.Abs(l.Residual[k]) > 1e-9*scale {
					t.Errorf("AC ledger at %g Hz: residual %.3g of %.3g", ae.Freq[k], l.Residual[k], scale)
				}
			}

			// Transient.
			ts := TranSpec{Step: 1e-5, Stop: 3e-3, MaxStep: 2e-6}
			we, err := Tran(ce, ts, TranOptions{})
			if err != nil {
				t.Fatal(err)
			}
			wt, err := Tran(tw, ts, TranOptions{})
			if err != nil {
				t.Fatal(err)
			}
			for _, n := range []string{"b", "c"} {
				ye, _ := we.NodeVoltage(n)
				yt, _ := wt.NodeVoltage(n)
				peak := 0.0
				for _, v := range yt {
					peak = math.Max(peak, math.Abs(v))
				}
				worst := 0.0
				for k, tt := range we.Time {
					worst = math.Max(worst, math.Abs(ye[k]-interp(wt.Time, yt, tt)))
				}
				if worst > 2e-4*peak {
					t.Errorf("TRAN V(%s): max difference %.3g of peak %.3g", n, worst, peak)
				}
			}
			if we.Ledger.MaxRelativeResidual > 1e-7 {
				t.Errorf("TRAN ledger max relative residual %.3g", we.Ledger.MaxRelativeResidual)
			}
			if we.Ledger.EnergyModelSupplied == 0 {
				t.Error("TRAN: no model-supplied energy")
			}
			in := math.Abs(we.Ledger.EnergyIn) + math.Abs(we.Ledger.EnergyDissipated)
			if m := math.Abs(we.Ledger.IntegrationMismatch); m > 1e-3*in {
				t.Errorf("TRAN integration mismatch %.3g J of %.3g J", m, in)
			}
		})
	}
}

// With the terms declared but all zero, or with no declaration, the
// engine is Maxwell's bit for bit.
func TestECEZeroIsMaxwell(t *testing.T) {
	plain := strings.Replace(eceDeck("ece2"), ".ece form=ece2 w0=5e-7 kappa=0.03\n", "", 1)
	for _, k := range []string{" kappa_a=0.02 kappa_b=-0.01", " w0=1e-6", " w0=8e-7"} {
		plain = strings.Replace(plain, k, "", 1)
	}
	zero := strings.Replace(plain, "ece network\n", "ece network\n.ece form=ece2 w0=0 kappa=0\n", 1)
	cp, cz := eceCircuit(t, plain), eceCircuit(t, zero)
	rp, err := OP(cp, Options{})
	if err != nil {
		t.Fatal(err)
	}
	rz, err := OP(cz, Options{})
	if err != nil {
		t.Fatal(err)
	}
	for i := range rp.NodeVoltages {
		if rp.NodeVoltages[i] != rz.NodeVoltages[i] {
			t.Errorf("OP node %d: %v vs %v", i, rp.NodeVoltages[i], rz.NodeVoltages[i])
		}
	}
	ts := TranSpec{Step: 1e-5, Stop: 1e-3}
	wp, err := Tran(cp, ts, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	wz, err := Tran(cz, ts, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if len(wp.Time) != len(wz.Time) {
		t.Fatalf("TRAN: %d vs %d points", len(wp.Time), len(wz.Time))
	}
	for n := 1; n < len(wp.V); n++ {
		for k := range wp.Time {
			if wp.V[n][k] != wz.V[n][k] {
				t.Fatalf("TRAN node %d at %g: %v vs %v", n, wp.Time[k], wp.V[n][k], wz.V[n][k])
			}
		}
	}
	if len(rz.Ledger.Supplied) != 0 {
		t.Errorf("zero terms booked model-supplied power: %+v", rz.Ledger.Supplied)
	}
}

// Torsion is a series resistance s·c·ω₀·L: at a series RLC's resonance the
// current is 1/(R + r_t). The same ω₀ gives r_t = +cω₀L under ECE and
// −2cω₀L under ECE2 (Q falls under ECE, rises under ECE2 for ω₀ > 0).
func TestECETorsionSeriesRLC(t *testing.T) {
	const R, L = 1.0, 10e-3
	w0 := 0.5 / (2 * cLight * L) // ECE2: r_t = −0.5 Ω
	C := 1 / (math.Pow(2*math.Pi*173, 2) * L)
	f0 := 1 / (2 * math.Pi * math.Sqrt(L*C))
	for _, tc := range []struct {
		form string
		rt   float64
	}{{"ece", 0.25}, {"ece2", -0.5}} {
		src := fmt.Sprintf("rlc\n.ece form=%s w0=%.17g\nV1 in 0 AC 1\nR1 in a 1\nL1 a b 10m\nC1 b 0 %.17g\n.end\n", tc.form, w0, C)
		c := eceCircuit(t, src)
		r, err := AC(c, ACSpec{Variation: "lin", Points: 1, FStart: f0, FStop: f0}, Options{})
		if err != nil {
			t.Fatal(err)
		}
		var i complex128
		for k, n := range r.ElementNames {
			if n == "L1" {
				i = r.I[k][0]
			}
		}
		want := 1 / (R + tc.rt)
		if d := math.Abs(cmplx.Abs(i) - want); d > 1e-9*want {
			t.Errorf("%s: |I(f0)| = %.12g A, want 1/(R + r_t) = %.12g", tc.form, cmplx.Abs(i), want)
		}
		// Average power: the source pays R + r_t; the terms supply −r_t·|I|²/2.
		ms := r.Ledger.ModelSupplied[0]
		if d := math.Abs(ms - (-tc.rt * want * want / 2)); d > 1e-9 {
			t.Errorf("%s: model-supplied %.12g W, want %.12g", tc.form, ms, -tc.rt*want*want/2)
		}
	}
}

// Declarations and parameters are checked, and survive a netlist round trip.
func TestECENetlist(t *testing.T) {
	bad := []struct{ name, src string }{
		{"param without .ece", "t\nV1 a 0 1\nR1 a 0 1 kappa=0.1\n.end\n"},
		{"unknown form", "t\n.ece form=ece3\nV1 a 0 1\nR1 a 0 1\n.end\n"},
		{"kappa on K", "t\n.ece form=ece\nV1 a 0 1\nL1 a 0 1m\nL2 b 0 1m\nR2 b 0 1\nK1 L1 L2 0.5 kappa=0.1\n.end\n"},
		{"nonlinear C", "t\n.ece form=ece\nV1 a 0 1\nR1 a b 1\nC1 b 0 Q=pwl(-1 -1n 1 1n) kappa=0.1\n.end\n"},
	}
	for _, b := range bad {
		d, err := netlist.ParseString(b.src)
		if err == nil {
			_, err = circuit.FromDeck(d)
		}
		if err == nil {
			t.Errorf("%s: accepted", b.name)
		} else if !errors.Is(err, circuit.ErrBadValue) && !errors.Is(err, netlist.ErrSyntax) {
			t.Errorf("%s: error %v", b.name, err)
		}
	}
	c := eceCircuit(t, eceDeck("ece2"))
	var sb strings.Builder
	if err := c.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	c2 := eceCircuit(t, sb.String())
	if c2.ECE == nil || *c2.ECE != *c.ECE {
		t.Fatalf("round trip .ece: %+v, want %+v\n%s", c2.ECE, c.ECE, sb.String())
	}
	for i, e := range c.Elements {
		if c2.ECEFor(c2.Elements[i]) != c.ECEFor(e) {
			t.Errorf("round trip %s: %+v, want %+v", e.Name, c2.ECEFor(c2.Elements[i]), c.ECEFor(e))
		}
	}
}

// UFT 456's conductivities (k0 on C, kl0 on L and K) against the same
// network built from ordinary elements: a resistor 1/g_e across C1
// (g_e = 2cκ₀C), a series resistor r_h = −2cκ_Λ0·L in each coil (negative:
// the paper's gain) and a CCVS −2cκ_Λ0·M for the mutual term. OP, AC and
// transient agree; the engine books every watt of the terms as
// model-supplied, the OP total in closed form.
func TestECEUFT456Twin(t *testing.T) {
	const (
		L1, L2, C1, kc = 10e-3, 5e-3, 10e-6, 0.5
		ge             = 1e-3 // S
		rh1            = -0.5 // Ω, on L1
	)
	k0 := ge / (2 * cLight * C1)
	kl0 := -rh1 / (2 * cLight * L1)
	rh2 := -2 * cLight * kl0 * L2
	rm := -2 * cLight * kl0 * kc * math.Sqrt(L1*L2)
	g := func(x float64) string { return fmt.Sprintf("%.17g", x) }
	body := `V1 in 0 DC 2 AC 1 SIN(2 1 1k)
R1 in a 10
R3 b 0 50
R2 c 0 20
`
	deck := "uft456\n.ece form=ece2 k0=" + g(k0) + " kl0=" + g(kl0) + "\n" + body +
		"L1 a b 10m\nC1 b 0 10u\nL2 c 0 5m\nK1 L1 L2 0.5\n.end\n"
	twin := "uft456 twin\n" + body + `Vs1 a y1 0
L1 y1 y2 10m
Rh1 y2 y3 ` + g(rh1) + `
Hh1 y3 b Vs2 ` + g(rm) + `
C1 b 0 10u
Rg1 b 0 ` + g(1/ge) + `
Vs2 c u1 0
L2 u1 u2 5m
Rh2 u2 u3 ` + g(rh2) + `
Hh2 u3 0 Vs1 ` + g(rm) + `
K1 L1 L2 0.5
.end
`
	ce, tw := eceCircuit(t, deck), eceCircuit(t, twin)

	re, err := OP(ce, Options{})
	if err != nil {
		t.Fatal(err)
	}
	rt, err := OP(tw, Options{})
	if err != nil {
		t.Fatal(err)
	}
	for _, n := range []string{"in", "a", "b", "c"} {
		ve, _ := re.Voltage(n)
		vt, _ := rt.Voltage(n)
		if d := eceRel(ve, vt); d > 1e-9 {
			t.Errorf("OP V(%s) = %.12g, twin %.12g", n, ve, vt)
		}
	}
	e1, _ := re.Element("L1")
	e2, _ := re.Element("L2")
	vb, _ := re.Voltage("b")
	want := -rh1*e1.Current*e1.Current - rh2*e2.Current*e2.Current - 2*rm*e1.Current*e2.Current - ge*vb*vb
	if d := eceRel(re.Ledger.ModelSupplied, want); d > 1e-9 {
		t.Errorf("OP model-supplied %.12g W, want −r_h·I² − 2r_m·I₁I₂ − g_e·V² = %.12g", re.Ledger.ModelSupplied, want)
	}
	if rr := re.Ledger.RelativeResidual(); rr > 1e-9 {
		t.Errorf("OP ledger residual %.3g", rr)
	}
	for _, row := range re.Ledger.Supplied {
		if !strings.Contains(row.Note, "(UFT 456)") {
			t.Errorf("supplied row %s note %q lacks the UFT 456 label", row.Element, row.Note)
		}
	}

	spec := ACSpec{Variation: "dec", Points: 5, FStart: 10, FStop: 1e5}
	ae, err := AC(ce, spec, Options{})
	if err != nil {
		t.Fatal(err)
	}
	at, err := AC(tw, spec, Options{})
	if err != nil {
		t.Fatal(err)
	}
	for _, n := range []string{"a", "b", "c"} {
		xe, _ := ae.NodeVoltage(n)
		xt, _ := at.NodeVoltage(n)
		for k := range xe {
			if d := cmplx.Abs(xe[k]-xt[k]) / math.Max(cmplx.Abs(xt[k]), 1e-12); d > 1e-9 {
				t.Errorf("AC V(%s) at %g Hz = %v, twin %v", n, ae.Freq[k], xe[k], xt[k])
			}
		}
	}
	for k := range ae.Freq {
		l := ae.Ledger
		scale := l.Inflow[k] + l.Outflow[k] + l.Dissipation[k] + math.Abs(l.ModelSupplied[k])
		if math.Abs(l.Residual[k]) > 1e-9*scale {
			t.Errorf("AC ledger at %g Hz: residual %.3g of %.3g", ae.Freq[k], l.Residual[k], scale)
		}
	}

	ts := TranSpec{Step: 1e-5, Stop: 3e-3, MaxStep: 2e-6}
	we, err := Tran(ce, ts, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	wt, err := Tran(tw, ts, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	for _, n := range []string{"b", "c"} {
		ye, _ := we.NodeVoltage(n)
		yt, _ := wt.NodeVoltage(n)
		peak, worst := 0.0, 0.0
		for _, v := range yt {
			peak = math.Max(peak, math.Abs(v))
		}
		for k, tt := range we.Time {
			worst = math.Max(worst, math.Abs(ye[k]-interp(wt.Time, yt, tt)))
		}
		if worst > 2e-4*peak {
			t.Errorf("TRAN V(%s): max difference %.3g of peak %.3g", n, worst, peak)
		}
	}
	if we.Ledger.MaxRelativeResidual > 1e-7 {
		t.Errorf("TRAN ledger max relative residual %.3g", we.Ledger.MaxRelativeResidual)
	}
	in := math.Abs(we.Ledger.EnergyIn) + math.Abs(we.Ledger.EnergyDissipated)
	if m := math.Abs(we.Ledger.IntegrationMismatch); m > 1e-3*in {
		t.Errorf("TRAN integration mismatch %.3g J of %.3g J", m, in)
	}
	t.Logf("OP model-supplied %.6g W; TRAN: source %.4g J, model-supplied %.4g J, mismatch %.2g J",
		re.Ledger.ModelSupplied, we.Ledger.EnergyIn, we.Ledger.EnergyModelSupplied, we.Ledger.IntegrationMismatch)

	for _, b := range []struct{ name, src string }{
		{"k0 on L", "t\n.ece form=ece2\nV1 a 0 1\nL1 a 0 1m k0=1\n.end\n"},
		{"kl0 on R", "t\n.ece form=ece2\nV1 a 0 1\nR1 a 0 1 kl0=1\n.end\n"},
		{"k0 without .ece", "t\nV1 a 0 1\nR1 a b 1\nC1 b 0 1u k0=1\n.end\n"},
	} {
		d, err := netlist.ParseString(b.src)
		if err == nil {
			_, err = circuit.FromDeck(d)
		}
		if err == nil {
			t.Errorf("%s: accepted", b.name)
		}
	}
	var sb strings.Builder
	if err := ce.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	c2 := eceCircuit(t, sb.String())
	if c2.ECE == nil || *c2.ECE != *ce.ECE {
		t.Errorf("round trip .ece: %+v, want %+v", c2.ECE, ce.ECE)
	}
}

// A K element takes W0= and KL0= in either order.
func TestECEMutualParams(t *testing.T) {
	for _, line := range []string{"K1 L1 L2 0.5 w0=1e-7 kl0=2e-8", "K1 L1 L2 0.5 KL0=2e-8 W0=1e-7"} {
		c := eceCircuit(t, "t\n.ece form=ece2\nV1 a 0 AC 1\nR1 a b 1\nL1 b 0 1m\nL2 c 0 1m\nR2 c 0 1\n"+line+"\n.end\n")
		for _, e := range c.Elements {
			if e.Name == "K1" {
				tm := c.ECEFor(e)
				if tm.W0 != 1e-7 || tm.KL0 != 2e-8 {
					t.Errorf("%q: w0 %g, kl0 %g", line, tm.W0, tm.KL0)
				}
			}
		}
	}
	if _, err := netlist.ParseString("t\n.ece form=ece2\nL1 a 0 1m\nL2 b 0 1m\nK1 L1 L2 0.5 kl0=1 kl0=2\n.end\n"); err == nil {
		t.Error("K1 with kl0 twice accepted")
	}
}
