// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/output"
)

// A series RLC port against its closed form Z = R + jωL + 1/(jωC) (C with
// a 1 GΩ leak, its DC path), for a
// positive, a zero and a negative R. The negative one lies outside the
// unit circle, |Γ| > 1, and the power it returns is booked by the ledger
// as model-supplied.
func TestSmithSeriesRLC(t *testing.T) {
	const L, C, rleak, z0 = 1e-6, 1e-9, 1e9, 50.0
	spec := ACSpec{Variation: "dec", Points: 20, FStart: 1e5, FStop: 1e8}
	for _, R := range []float64{20, 0, -20} {
		src := "rlc\nR1 a b 20\nL1 b c 1u\nC1 c 0 1n\nRleak c 0 1e9\n.end\n"
		if R == 0 {
			src = "rlc\nL1 a c 1u\nC1 c 0 1n\nRleak c 0 1e9\n.end\n"
		}
		if R < 0 {
			src = "rlc\nR1 a b -20\nL1 b c 1u\nC1 c 0 1n\nRleak c 0 1e9\n.end\n"
		}
		c := mustDeck(t, src)
		n := len(c.Elements)
		r, err := Smith(c, "a", "0", z0, spec, Options{})
		if err != nil {
			t.Fatalf("R = %g: %v", R, err)
		}
		if len(c.Elements) != n {
			t.Fatalf("Smith changed the circuit")
		}
		for _, p := range r.Points {
			w := 2 * math.Pi * p.Freq
			zc := 1 / complex(1/rleak, w*C) // C ∥ Rleak, the DC path to ground
			want := complex(R, w*L) + zc
			if cmplx.Abs(p.Z-want) > 1e-9*cmplx.Abs(want) {
				t.Fatalf("R = %g, f = %g: Z = %v, want %v", R, p.Freq, p.Z, want)
			}
			g := (want - z0) / (want + z0)
			if cmplx.Abs(p.Gamma-g) > 1e-9*math.Max(1, cmplx.Abs(g)) {
				t.Fatalf("R = %g: Γ = %v, want %v", R, p.Gamma, g)
			}
			mag := cmplx.Abs(p.Gamma)
			switch {
			case R > 0 && (p.Region != output.SmithPassive || mag >= 1):
				t.Fatalf("R = %g: region %v |Γ| %g", R, p.Region, mag)
			case R == 0 && (p.Region == output.SmithActive || math.Abs(mag-1) > 1e-6):
				t.Fatalf("R = 0: region %v |Γ| %g", p.Region, mag)
			case R < 0:
				if p.Region != output.SmithActive || mag <= 1 {
					t.Fatalf("R = %g: region %v |Γ| %g", R, p.Region, mag)
				}
				// 1 A peak through −20 Ω returns ½·20·1² = 10 W, paid by R1
				// (less what the leak dissipates).
				if !near(p.PortPower, 0.5*real(want), 1e-9) || !near(p.ModelSupplied, 10, 1e-9) {
					t.Fatalf("port power %g, model-supplied %g, want −10, 10", p.PortPower, p.ModelSupplied)
				}
			}
		}
		if R < 0 {
			b, n := r.ActiveBands(), len(r.Points)
			if len(b) != 1 || b[0][0] != r.Points[0].Freq || b[0][1] != r.Points[n-1].Freq {
				t.Fatalf("active bands %v", b)
			}
		}
	}
}

// The circuit's own AC sources are set to zero: a voltage source becomes
// a short, a current source an open, so the port sees R1 alone.
func TestSmithZeroesSources(t *testing.T) {
	c := mustDeck(t, "src\nV1 x 0 DC 1 AC 5\nR1 a x 75\nI1 0 a AC 2\n.end\n")
	r, err := Smith(c, "a", "0", 50, ACSpec{Variation: "lin", Points: 3, FStart: 1e3, FStop: 1e6}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range r.Points {
		if cmplx.Abs(p.Z-75) > 1e-9 || cmplx.Abs(p.Gamma-0.2) > 1e-12 {
			t.Fatalf("Z = %v Γ = %v, want 75 Ω, 0.2", p.Z, p.Gamma)
		}
	}
	if e := c.Elements[c.Index("V1")]; e.ACMag != 5 {
		t.Fatalf("original circuit changed: V1 AC %g", e.ACMag)
	}
}

// A negative conductance made from a VCCS: G1 draws −10 mS × V(a) out of
// the port, in parallel with 200 Ω, so Z = 1/(1/200 − 1/100) = −200 Ω,
// active, with Γ = (−200 − 50)/(−200 + 50) = 5/3.
func TestSmithNegativeConductance(t *testing.T) {
	c := mustDeck(t, "g\nG1 a 0 0 a 0.01\nR1 a 0 200\n.end\n")
	r, err := Smith(c, "a", "0", 50, ACSpec{Variation: "lin", Points: 2, FStart: 1e3, FStop: 2e3}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	p := r.Points[0]
	if cmplx.Abs(p.Z+200) > 1e-9 || cmplx.Abs(p.Gamma-5.0/3) > 1e-12 || p.Region != output.SmithActive {
		t.Fatalf("Z = %v Γ = %v region %v", p.Z, p.Gamma, p.Region)
	}
}

func TestSmithErrors(t *testing.T) {
	c := mustDeck(t, "e\nR1 a 0 50\n.end\n")
	spec := ACSpec{Variation: "lin", Points: 2, FStart: 1, FStop: 2}
	for _, z0 := range []float64{0, -50, math.Inf(1), math.NaN()} {
		if _, err := Smith(c, "a", "0", z0, spec, Options{}); !errors.Is(err, ErrBadReference) {
			t.Fatalf("z0 = %g: err %v", z0, err)
		}
	}
	if _, err := Smith(c, "nope", "0", 50, spec, Options{}); !errors.Is(err, ErrUnknownNode) {
		t.Fatalf("unknown node: %v", err)
	}
	if _, err := Smith(c, "a", "a", 50, spec, Options{}); !errors.Is(err, ErrUnknownNode) {
		t.Fatalf("same node: %v", err)
	}
	// A port matched to Z0 sits at the centre.
	r, err := Smith(c, "a", "gnd", 50, spec, Options{})
	if err != nil || cmplx.Abs(r.Points[0].Gamma) > 1e-12 {
		t.Fatalf("matched port: %v %v", r, err)
	}
}
