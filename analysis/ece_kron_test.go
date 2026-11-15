// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/netlist"
)

// The hypothetical ECE / ECE2 terms (H48) read against Kron, Tensors for
// Circuits (Dover 1959), chapters 31 and 32:
//
//   - p. 245: an equation has an equivalent circuit only if it is a tensor
//     equation, so it must not depend on the reference chosen. The datum
//     tests show which terms pass and which fail that criterion.
//   - eq. 31.5: in the compound equation of motion p = r·ẋ + a·pẋ + F·ẋ
//     the skew F does no work (ẋ·F·ẋ = 0). eqs. 31.9–31.11: a speed matrix
//     G that is independent of the metric (the commutator machine) is
//     genuine torsion, and its power goes through the shaft.
//   - eq. 31.14: the affine connection is not a tensor. A frame change adds
//     the inhomogeneous term Cᵀ·L·dC/dt, whose symmetric part is
//     ½·d(CᵀLC)/dt.

// liftDeck is a coupled RLC network whose only connection to node 0 is a
// source V0 that lifts every potential by V0 (a shift of the datum).
func liftDeck(ece string, v0 float64, r1 string) string {
	return fmt.Sprintf(`lift
%s
V0 g 0 DC %.17g
V1 in g DC 2 AC 1
%s
L1 a b 10m
C1 b g 10u
R3 b g 50
L2 c g 5m
R2 c g 20
K1 L1 L2 0.5
.end
`, ece, v0, r1)
}

// Item 1: shift of the datum. Maxwell, torsion, UFT 456's k0 and kl0 and a
// UNIFORM κ do not depend on it. A κ that differs between an element's two
// leads does: under φ → φ + φ0 its body voltage α_a·φ_a − α_b·φ_b gains
// s·(κ_a − κ_b)·φ0, an EMF set by the absolute potential. That is not a
// tensor equation under a change of reference (Kron, p. 245), so the
// result depends on where the datum is. The test records the dependence
// (linear in φ0, equal to the predicted EMF) instead of asserting
// invariance.
func TestECEDatumShift(t *testing.T) {
	const plainR1 = "R1 in a 10"
	cases := []struct {
		name, ece, r1 string
		invariant     bool
	}{
		{"maxwell", "", plainR1, true},
		{"torsion", ".ece form=ece2 w0=5e-7", plainR1, true},
		{"uft456", ".ece form=ece2 k0=2e-4 kl0=3e-7", plainR1, true},
		{"uniform-kappa", ".ece form=ece2 kappa=0.03", plainR1, true},
		{"lead-kappa", ".ece form=ece2", "R1 in a 10 kappa_a=0.02 kappa_b=-0.01", false},
		{"lead-kappa-ece", ".ece form=ece", "R1 in a 10 kappa_a=0.02 kappa_b=-0.01", false},
	}
	branch := []string{"R1", "L1", "L2", "R3"}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			cur := func(v0 float64) []float64 {
				r, err := OP(eceCircuit(t, liftDeck(tc.ece, v0, tc.r1)), Options{})
				if err != nil {
					t.Fatal(err)
				}
				if rr := r.Ledger.RelativeResidual(); rr > 1e-9 {
					t.Fatalf("V0 = %g: ledger residual %.3g", v0, rr)
				}
				var out []float64
				for _, n := range branch {
					e, _ := r.Element(n)
					out = append(out, e.Current)
				}
				return out
			}
			i0, i5, i50 := cur(0), cur(5), cur(50)
			worst := 0.0
			for k := range i0 {
				worst = math.Max(worst, math.Abs(i50[k]-i0[k])/math.Max(math.Abs(i0[k]), 1e-12))
			}
			if tc.invariant {
				if worst > 1e-12 {
					t.Errorf("datum shift by 50 V changed the currents by %.3g (relative)", worst)
				}
				return
			}
			if worst < 1e-3 {
				t.Fatalf("lead κ: the datum shift changed nothing (%.3g)", worst)
			}
			for k := range i0 {
				d5, d50 := i5[k]-i0[k], i50[k]-i0[k]
				if math.Abs(d50-10*d5) > 1e-9*math.Abs(d50) {
					t.Errorf("%s: ΔI not linear in φ0: %g at 50 V, %g at 5 V", branch[k], d50, d5)
				}
			}
			// Twin of the lifted circuit: V0 = 0 with the predicted EMF
			// s·(κ_a − κ_b)·φ0 placed at R1's b lead as an ordinary source
			// (φ_x = φ_a + e, so −α_b·e = s·Δκ·φ0).
			s := 1.0
			if strings.Contains(tc.ece, "ece2") {
				s = -2
			}
			const v0, ka, kb = 50.0, 0.02, -0.01
			e := -s * (ka - kb) * v0 / (1 + s*kb)
			twin := strings.Replace(liftDeck(tc.ece, 0, "R1 in x 10 kappa_a=0.02 kappa_b=-0.01"),
				"L1 a b 10m", fmt.Sprintf("VE x a DC %.17g\nL1 a b 10m", e), 1)
			rt, err := OP(eceCircuit(t, twin), Options{})
			if err != nil {
				t.Fatal(err)
			}
			for k, n := range branch {
				el, _ := rt.Element(n)
				if d := eceRel(el.Current, i50[k]); d > 1e-9 {
					t.Errorf("%s: lifted %.12g A, twin with EMF s·Δκ·φ0 %.12g A", n, i50[k], el.Current)
				}
			}
			t.Logf("form %s: lifting the datum by 50 V changes the currents by up to %.3g (relative), = an EMF %.4g V at R1's lead", tc.ece, worst, -e*(1+s*kb))
		})
	}
}

// Item 1, corollary: a uniform κ (the same on every lead, 0 at ground) is a
// pure scaling. Every passive's body voltage is α·v, α = 1 + s·κ, so the
// network equals the Maxwell network with R/α, L/α and C·α (and torsion
// r_t/α). That is why it passes the datum test.
func TestECEUniformKappaIsScaling(t *testing.T) {
	const kappa, w0 = 0.03, 5e-7
	for _, tc := range []struct {
		form string
		s    float64
	}{{"ece", 1}, {"ece2", -2}} {
		t.Run(tc.form, func(t *testing.T) {
			a := 1 + tc.s*kappa
			g := func(x float64) string { return fmt.Sprintf("%.17g", x) }
			rt1 := tc.s * cLight * w0 * 10e-3
			ece := fmt.Sprintf(`u
.ece form=%s kappa=%g w0=%g
V1 in 0 DC 2 AC 1
R1 in a 10
L1 a b 10m
C1 b 0 10u
R3 b 0 50
.end
`, tc.form, kappa, w0)
			plain := `u
V1 in 0 DC 2 AC 1
R1 in a ` + g(10/a) + `
L1 a y ` + g(10e-3/a) + `
RT y b ` + g(rt1/a) + `
C1 b 0 ` + g(10e-6*a) + `
R3 b 0 ` + g(50/a) + `
.end
`
			ce, cp := eceCircuit(t, ece), eceCircuit(t, plain)
			re, err := OP(ce, Options{})
			if err != nil {
				t.Fatal(err)
			}
			rp, err := OP(cp, Options{})
			if err != nil {
				t.Fatal(err)
			}
			for _, n := range []string{"a", "b"} {
				ve, _ := re.Voltage(n)
				vp, _ := rp.Voltage(n)
				if d := eceRel(ve, vp); d > 1e-9 {
					t.Errorf("OP V(%s) = %.12g, scaled Maxwell %.12g", n, ve, vp)
				}
			}
			spec := ACSpec{Variation: "dec", Points: 5, FStart: 10, FStop: 1e5}
			ae, err := AC(ce, spec, Options{})
			if err != nil {
				t.Fatal(err)
			}
			ap, err := AC(cp, spec, Options{})
			if err != nil {
				t.Fatal(err)
			}
			xe, _ := ae.NodeVoltage("b")
			xp, _ := ap.NodeVoltage("b")
			for k := range xe {
				if d := cmplx.Abs(xe[k]-xp[k]) / cmplx.Abs(xp[k]); d > 1e-9 {
					t.Errorf("AC V(b) at %g Hz: %v, scaled Maxwell %v", ae.Freq[k], xe[k], xp[k])
				}
			}
		})
	}
}

// Item 2: the ECE torsion term r_t = s·c·ω₀·L against its standard
// counterpart, the series commutator machine (Tensors for Circuits
// eqs. 31.9–31.11): a winding of inductance L whose speed matrix G is
// independent of L, turned at a held speed ω. The winding equation is
// v = L·di/dt + ω·G·i, which is the torsion relation when ω·G = r_t, so the
// currents agree. The books differ. ECE books r_t·i² as model-supplied
// power from nowhere. The machine passes the same power through its shaft
// to the speed source, which pays it when r_t < 0 (ECE2 with ω₀ > 0: the
// self-exciting series dynamo, H15).
func TestECECommutatorTwin(t *testing.T) {
	const (
		R, L  = 10.0, 10e-3
		omega = 100.0 // rad/s held by the speed source
	)
	sin := circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{1, 1, 1e3}}
	for _, tc := range []struct {
		form string
		s    float64
	}{{"ece", 1}, {"ece2", -2}} {
		t.Run(tc.form, func(t *testing.T) {
			w0 := 3 / (cLight * L) // |r_t| = 3 Ω (ECE) or 6 Ω (ECE2)
			rt := tc.s * cLight * w0 * L
			ce := eceCircuit(t, fmt.Sprintf("ece coil\n.ece form=%s w0=%.17g\nV1 in 0 SIN(1 1 1k)\nR1 in a %g\nL1 a 0 %g\n.end\n", tc.form, w0, R, L))
			spec := circuit.NewMachineSpec(1)
			spec.L0[0][0] = L
			spec.G[0][0] = rt / omega
			cm, err := circuit.NewBuilder("commutator").Domain(exchange.Mechanical, "w").
				MachineModel("cm", spec).
				VoltageWave("V1", "in", "0", sin).Resistor("R1", "in", "a", R).
				Machine("N1", "cm", [][2]string{{"a", "0"}}, "w", "0").
				VoltageSource("VW", "w", "0", omega).BuildTransient()
			if err != nil {
				t.Fatal(err)
			}

			// Operating point: i = v/(R + r_t) in both.
			re, err := OP(ce, Options{})
			if err != nil {
				t.Fatal(err)
			}
			rm, err := OP(cm, Options{})
			if err != nil {
				t.Fatal(err)
			}
			ie, _ := re.Element("L1")
			im, _ := rm.Element("N1")
			want := 1 / (R + rt)
			if eceRel(ie.Current, want) > 1e-9 || eceRel(im.Current, want) > 1e-9 {
				t.Errorf("OP: ECE %.12g A, machine %.12g A, want 1/(R + r_t) = %.12g", ie.Current, im.Current, want)
			}
			shaft := 0.0
			for _, row := range rm.Ledger.Rows {
				if row.Element == "VW" {
					shaft += row.PowerW
				}
			}
			if d := math.Abs(shaft - re.Ledger.ModelSupplied); d > 1e-9*math.Abs(shaft) {
				t.Errorf("OP: ECE model-supplied %.12g W, machine shaft source %.12g W", re.Ledger.ModelSupplied, shaft)
			}
			if rm.Ledger.ModelSupplied != 0 {
				t.Errorf("OP: the machine booked model-supplied power %g", rm.Ledger.ModelSupplied)
			}
			for _, r := range []*struct{ rr float64 }{{re.Ledger.RelativeResidual()}, {rm.Ledger.RelativeResidual()}} {
				if r.rr > 1e-9 {
					t.Errorf("OP ledger residual %.3g", r.rr)
				}
			}

			// Transient: same currents; the shaft source's energy equals the
			// energy ECE books as model-supplied.
			ts := TranSpec{Step: 1e-5, Stop: 3e-3, MaxStep: 2e-6}
			we, err := Tran(ce, ts, TranOptions{})
			if err != nil {
				t.Fatal(err)
			}
			wm, err := Tran(cm, ts, TranOptions{})
			if err != nil {
				t.Fatal(err)
			}
			ye, _ := we.ElementCurrent("L1")
			ym, _ := wm.ElementCurrent("N1")
			peak, worst := 0.0, 0.0
			for k, tt := range we.Time {
				peak = math.Max(peak, math.Abs(ye[k]))
				worst = math.Max(worst, math.Abs(ye[k]-interp(wm.Time, ym, tt)))
			}
			if worst > 1e-4*peak {
				t.Errorf("TRAN I: max difference %.3g of peak %.3g", worst, peak)
			}
			eShaft := (wm.Ledger.EnergyIn - wm.Ledger.EnergyOut) - (we.Ledger.EnergyIn - we.Ledger.EnergyOut) // V1 is common
			eECE := we.Ledger.EnergyModelSupplied
			if d := math.Abs(eShaft - eECE); d > 1e-3*math.Abs(eECE) {
				t.Errorf("TRAN: ECE model-supplied %.6g J, machine shaft %.6g J", eECE, eShaft)
			}
			if wm.Ledger.EnergyModelSupplied != 0 {
				t.Errorf("TRAN: the machine booked model-supplied energy %g J", wm.Ledger.EnergyModelSupplied)
			}
			if wm.Ledger.MaxRelativeResidual > 1e-7 || we.Ledger.MaxRelativeResidual > 1e-7 {
				t.Errorf("TRAN ledger residuals %.3g (ECE), %.3g (machine)", we.Ledger.MaxRelativeResidual, wm.Ledger.MaxRelativeResidual)
			}
			t.Logf("%s: r_t = %g Ω; over 3 ms ECE books %.4g J model-supplied, the commutator machine passes %.4g J through its shaft", tc.form, rt, eECE, eShaft)
		})
	}
}

// skewDeck: two equal coils, each discharging into its own resistor from
// an initial current, joined only by a K element whose coupling k is
// negligible (k = 1e-9, so M·di/dt is 1e-9 of L·di/dt), carrying the
// mutual torsion in symmetric ("w0") or skew ("w0skew") placement.
func skewDeck(form, key string, w float64) string {
	return fmt.Sprintf(`skew
.ece form=%s
L1 a 0 10m IC=1
R1 a 0 2
L2 b 0 10m IC=0
R2 b 0 2
K1 L1 L2 1e-9 %s=%.17g
.end
`, form, key, w)
}

// Items 3 and 4: placement of the mutual torsion.
//
// The skew placement (v1 += r_s·i2, v2 −= r_s·i1, Kron's F) does no work:
// the coils' energy ½L(i1² + i2²) decays exactly as without it,
// e^(−2Rt/L), and the ledger books no model-supplied energy. Its effect is
// that the current vector turns at Ω = r_s/L. That is the plain network
// seen from a reference frame rotating at Ω, the frame-generated term
// L·Cᵀ·dC/dt (eq. 31.14) being skew, so this part of a declared spin
// connection can be transformed away.
//
// The symmetric placement (v1 += r_m·i2, v2 += r_m·i1, the engine's w0 on
// K) splits the decay rates of i1 ± i2 to (R ± r_m)/L. It changes |i|,
// which no rotation does, and books model-supplied energy. No frame with a
// constant transformed inductance can generate or remove it.
func TestECETorsionPlacement(t *testing.T) {
	const (
		L, R = 10e-3, 2.0
		k    = 1e-9
	)
	M := k * L
	omegaRot := 2 * math.Pi * 50 // the frame's rate Ω = r_s/L
	rs := omegaRot * L
	stop := 20e-3 // one turn
	ts := TranSpec{Step: 1e-5, Stop: stop, MaxStep: 1e-5, UIC: true}
	run := func(form, key string, r float64) (tt, i1, i2 []float64, en, ems float64) {
		s := 1.0
		if form == "ece2" {
			s = -2
		}
		w, err := Tran(eceCircuit(t, skewDeck(form, key, r/(s*cLight*M))), ts, TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		i1, _ = w.ElementCurrent("L1")
		i2, _ = w.ElementCurrent("L2")
		if w.Ledger.MaxRelativeResidual > 1e-7 {
			t.Errorf("%s %s: ledger residual %.3g", form, key, w.Ledger.MaxRelativeResidual)
		}
		return w.Time, i1, i2, w.Ledger.EnergyDissipated, w.Ledger.EnergyModelSupplied
	}
	for _, form := range []string{"ece", "ece2"} {
		tt, i1, i2, heat, ems := run(form, "w0skew", rs)
		worstMag, worstAng := 0.0, 0.0
		for n, x := range tt {
			env := math.Exp(-R * x / L)
			worstMag = math.Max(worstMag, math.Abs(math.Hypot(i1[n], i2[n])-env))
			worstAng = math.Max(worstAng, math.Hypot(i1[n]-env*math.Cos(omegaRot*x), i2[n]-env*math.Sin(omegaRot*x)))
		}
		if worstMag > 1e-4 || worstAng > 2e-3 {
			t.Errorf("%s skew: |i| off the plain decay by %.3g, off the rotating frame by %.3g", form, worstMag, worstAng)
		}
		if math.Abs(ems) > 1e-9*heat {
			t.Errorf("%s skew: model-supplied %.3g J (heat %.3g J); the skew term must do no work", form, ems, heat)
		}
		w0 := 0.5 * L * (1 - math.Exp(-2*R*stop/L)) // energy the plain network turns into heat
		if d := math.Abs(heat - w0); d > 1e-4*w0 {
			t.Errorf("%s skew: heat %.6g J, plain network %.6g J", form, heat, w0)
		}
		t.Logf("%s skew r_s = %.4g Ω: |i| equals the plain decay to %.2g, the current turns at r_s/L to %.2g, model-supplied %.2g J", form, rs, worstMag, worstAng, ems)

		// Symmetric placement of the same size: modes i1 ± i2 decay at
		// (R ± r_m)/L, |i| leaves the plain envelope, energy is supplied.
		tt, i1, i2, _, ems = run(form, "w0", rs)
		worstSym, dev := 0.0, 0.0
		for n, x := range tt {
			sum := math.Exp(-(R + rs) * x / L) // i1 + i2 from (1, 0)
			dif := math.Exp(-(R - rs) * x / L) // i1 − i2
			worstSym = math.Max(worstSym, math.Max(math.Abs(i1[n]-(sum+dif)/2), math.Abs(i2[n]-(sum-dif)/2)))
			dev = math.Max(dev, math.Abs(math.Hypot(i1[n], i2[n])-math.Exp(-R*x/L)))
		}
		if worstSym > 2e-3 {
			t.Errorf("%s symmetric: modes off (R ± r_m)/L by %.3g", form, worstSym)
		}
		if dev < 0.1 || ems == 0 {
			t.Errorf("%s symmetric: |i| deviation %.3g, model-supplied %g J: expected a change no frame gives", form, dev, ems)
		}
	}
}

// Item 4, the matrix form of eq. 31.14 behind TestECETorsionPlacement:
// for a frame C(t) the frame-generated term S = Cᵀ·L·dC/dt has symmetric
// part ½·d(CᵀLC)/dt. So when the transformed inductance stays constant
// (here: C orthogonal and L = ℓ·1), S is skew. A frame change then
// generates or removes only skew (workless) terms, never a symmetric
// resistance such as the ECE torsion r_t = s·c·ω₀·L.
func TestFrameTermSymmetricPart(t *testing.T) {
	// L: a fixed symmetric positive-definite 3×3; C(t): a rotation about a
	// turning axis times a time-varying scale on one axis.
	Lm := [3][3]float64{{4, 1, 0.5}, {1, 3, 0.2}, {0.5, 0.2, 2}}
	rot := func(t float64) [3][3]float64 {
		a := 0.7 * t
		ax := [3]float64{math.Cos(t), math.Sin(t), 0.5}
		n := math.Sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2])
		for i := range ax {
			ax[i] /= n
		}
		c, s := math.Cos(a), math.Sin(a)
		var r [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				d := 0.0
				if i == j {
					d = 1
				}
				r[i][j] = c*d + (1-c)*ax[i]*ax[j]
			}
		}
		r[0][1] -= s * ax[2]
		r[1][0] += s * ax[2]
		r[0][2] += s * ax[1]
		r[2][0] -= s * ax[1]
		r[1][2] -= s * ax[0]
		r[2][1] += s * ax[0]
		return r
	}
	general := func(t float64) [3][3]float64 {
		r := rot(t)
		for i := 0; i < 3; i++ {
			r[i][0] *= 1 + 0.3*math.Sin(2*t)
		}
		return r
	}
	mul := func(a, b [3][3]float64) (c [3][3]float64) {
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				for k := 0; k < 3; k++ {
					c[i][j] += a[i][k] * b[k][j]
				}
			}
		}
		return
	}
	tr := func(a [3][3]float64) (b [3][3]float64) {
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				b[i][j] = a[j][i]
			}
		}
		return
	}
	const h = 1e-5
	deriv := func(f func(float64) [3][3]float64, t float64) (d [3][3]float64) {
		p, m := f(t+h), f(t-h)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				d[i][j] = (p[i][j] - m[i][j]) / (2 * h)
			}
		}
		return
	}
	for _, t0 := range []float64{0.1, 0.9, 2.3} {
		// General frame: sym(S) = ½·d(CᵀLC)/dt.
		C, dC := general(t0), deriv(general, t0)
		S := mul(tr(C), mul(Lm, dC))
		dLp := deriv(func(t float64) [3][3]float64 { c := general(t); return mul(tr(c), mul(Lm, c)) }, t0)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				if d := math.Abs(S[i][j] + S[j][i] - dLp[i][j]); d > 1e-8 {
					t.Errorf("t = %g: S + Sᵀ − d(CᵀLC)/dt [%d][%d] = %.3g", t0, i, j, d)
				}
			}
		}
		// Orthogonal frame, L = 2·1: S is skew.
		C, dC = rot(t0), deriv(rot, t0)
		S = mul(tr(C), dC)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				if d := math.Abs(2 * (S[i][j] + S[j][i])); d > 1e-8 {
					t.Errorf("t = %g: orthogonal frame, symmetric part [%d][%d] = %.3g", t0, i, j, d)
				}
			}
		}
	}
}

// w0skew is declared on K only, and survives a netlist round trip.
func TestECESkewNetlist(t *testing.T) {
	c := eceCircuit(t, skewDeck("ece2", "w0skew", 0.25))
	var sb strings.Builder
	if err := c.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(sb.String(), "W0SKEW=0.25") {
		t.Errorf("round trip lost W0SKEW:\n%s", sb.String())
	}
	bad := strings.Replace(skewDeck("ece2", "w0skew", 0.25), "L2 b 0 10m IC=0", "L2 b 0 10m IC=0 w0skew=1", 1)
	d, err := netlist.ParseString(bad)
	if err == nil {
		if _, err = circuit.FromDeck(d); err == nil {
			t.Error("w0skew on an inductor was accepted")
		}
	}
}
