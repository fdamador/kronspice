// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Machine checks from L. V. Bewley, Tensor Analysis of Electric Circuits
// and Machines, Ronald Press 1961 ("Bewley 1961").

// Problem 15‑3: a salient-pole rotor with slotting has
// L = ½(L_d + L_q) + ½(L_d − L_q)·cos 2θ + ΔL·cos sθ. With a DC current I
// at standstill the electric torque is ½·I²·∂L/∂θ (Eq. 15‑5) =
// −½·I²·[(L_d − L_q)·sin 2θ + s·ΔL·sin sθ], the reluctance torque; through
// a damper of mobility Rd the shaft settles at ω = Rd·T, and the winding
// shows the speed voltage ω·I·∂L/∂θ (Eq. 15‑2).
func TestBewleyReluctanceTorque(t *testing.T) {
	const (
		ld, lq, dl = 12e-3, 7e-3, 0.4e-3
		sl         = 6
		I          = 3.0
		th0        = 0.4
		rd         = 2.0
	)
	spec := circuit.NewMachineSpec(1)
	spec.L0[0][0] = (ld + lq) / 2
	spec.Harmonic(2).Cos[0][0] = (ld - lq) / 2
	spec.Harmonic(sl).Cos[0][0] = dl
	spec.Theta0 = th0
	c, err := circuit.NewBuilder("reluctance").Domain(exchange.Mechanical, "w").
		MachineModel("mm", spec).
		CurrentSource("I1", "0", "a", I).Resistor("RA", "a", "0", 1e6).
		Machine("N1", "mm", [][2]string{{"a", "0"}}, "w", "0").
		Resistor("RD", "w", "0", rd).Build()
	if err != nil {
		t.Fatal(err)
	}
	r, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	dL := -(ld-lq)*math.Sin(2*th0) - sl*dl*math.Sin(sl*th0)
	T := 0.5 * I * I * dL
	w := r.NodeVoltages[c.NodeIndex("w")]
	if math.Abs(w-rd*T) > 1e-9*math.Abs(rd*T) {
		t.Errorf("ω = %g, want Rd·T = %g", w, rd*T)
	}
	va := r.NodeVoltages[c.NodeIndex("a")]
	if want := w * I * dL; math.Abs(va-want) > 1e-6*math.Abs(want) {
		t.Errorf("speed voltage %g, want %g", va, want)
	}
	if rr := r.Ledger.RelativeResidual(); rr > 1e-9 {
		t.Errorf("ledger residual %g", rr)
	}
}

// Eq. 19‑53, 19‑54, hunting: a rotor between DC-fed stator and rotor
// windings (mutual M·cos θ) has the synchronizing torque coefficient
// T_S = M·I₁·I₂ at θ = 0; with inertia J and damping T_D (a damper of
// mobility 1/T_D on the shaft) a small swing obeys
// J·p²Δθ + T_D·pΔθ + T_S·Δθ = 0, so the rotor hunts at
// ω_h = √(T_S/J − (T_D/2J)²) and the swing decays at T_D/2J. The
// transient ω(t) from a small initial angle matches that solution, and
// the ledger closes.
func TestBewleyHunting(t *testing.T) {
	const (
		l1, l2, m = 10e-3, 20e-3, 5e-3
		i1, i2    = 2.0, 3.0
		th0       = 0.01
		J         = 1e-4
		alpha     = 2.0 // T_D/2J [1/s]
	)
	ts := m * i1 * i2
	td := 2 * J * alpha
	wh := math.Sqrt(ts/J - alpha*alpha)
	c, err := circuit.NewBuilder("hunting").Domain(exchange.Mechanical, "w").
		MachineModel("mm", withIC(rotatingMutual(l1, l2, m, th0), i1, i2)).
		CurrentSource("I1", "0", "a", i1).CurrentSource("I2", "0", "b", i2).
		Resistor("RA", "a", "0", 1e9).Resistor("RB", "b", "0", 1e9).
		Machine("N1", "mm", [][2]string{{"a", "0"}, {"b", "0"}}, "w", "0").
		CapacitorIC("CJ", "w", "0", J, 0).
		Resistor("RD", "w", "0", 1/td).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	period := 2 * math.Pi / wh
	wf, err := Tran(c, TranSpec{Step: period / 400, Stop: 4 * period, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	w, _ := wf.NodeVoltage("w")
	// Linear solution with θ(0) = θ₀, θ̇(0) = 0.
	amp := th0 * (alpha*alpha + wh*wh) / wh
	worst := 0.0
	var crossings []float64
	for i, tt := range wf.Time {
		ref := -amp * math.Exp(-alpha*tt) * math.Sin(wh*tt)
		worst = math.Max(worst, math.Abs(w[i]-ref))
		if i > 0 && tt > period/4 && (w[i-1] < 0) != (w[i] < 0) {
			f := w[i-1] / (w[i-1] - w[i])
			crossings = append(crossings, wf.Time[i-1]+f*(tt-wf.Time[i-1]))
		}
	}
	if worst > 1e-2*amp {
		t.Errorf("ω differs from Eq. 19‑54's solution by %.3g of its amplitude", worst/amp)
	}
	if len(crossings) < 5 {
		t.Fatalf("only %d zero crossings", len(crossings))
	}
	half := (crossings[len(crossings)-1] - crossings[0]) / float64(len(crossings)-1)
	if got := math.Pi / half; math.Abs(got-wh) > 2e-3*wh {
		t.Errorf("hunting frequency %.5g rad/s, Eq. 19‑54 gives %.5g", got, wh)
	}
	if l := wf.Ledger; l.MaxRelativeResidual > 1e-9 {
		t.Errorf("ledger residual %g", l.MaxRelativeResidual)
	}
	t.Logf("ω_h = %.5g rad/s (√(T_S/J) = %.5g), decay %.3g 1/s; worst error %.2g of amplitude",
		math.Pi/half, math.Sqrt(ts/J), alpha, worst/amp)
}

// Ch. 21, Eqs. 21‑7 and 21‑11: a symmetrical two-phase induction motor
// (stator self reactance X_s, rotor X_r, mutual X_m, resistances r_s, r)
// at slip s takes per phase the rotor current of the equivalent circuit
// Z₁ = r_s + jX_s + X_m²/(r/s + jX_r), and Bewley's Eq. 21‑11
//
//	T₁₁ = e²·s·r·X_m² / ([r_s·r − s(X_sX_r − X_m²)]² + (r·X_s + s·r_s·X_r)²)
//
// equals |I_r|²·r/s, the air-gap power of one phase per unit of phase
// voltage amplitude squared. With peak phasors each phase contributes
// ½·T₁₁, so the two phases give a shaft torque T = T₁₁/ω (two-pole). The
// motor is built from the general machine element (windings a, b on the
// stator; A, B on the rotor with mutual X_m·cos θ, rotated 90°), fed a
// balanced two-phase supply, its shaft held at (1 − s)·ω by a velocity
// source; once the transient has died away the torque is constant and
// matches the formula.
func TestBewleyInductionTorque(t *testing.T) {
	const (
		f                    = 50.0
		rs, xs, r, xr, xm, s = 0.5, 10.4, 0.4, 10.6, 10.0, 0.05
		E                    = 100.0
	)
	w := 2 * math.Pi * f
	spec := circuit.NewMachineSpec(4) // a, b, A, B
	spec.L0[0][0], spec.L0[1][1] = xs/w, xs/w
	spec.L0[2][2], spec.L0[3][3] = xr/w, xr/w
	h := spec.Harmonic(1)
	M := xm / w
	// L_aA = L_bB = M cos θ; L_aB = −M sin θ; L_bA = M sin θ.
	h.Cos[0][2], h.Cos[2][0], h.Cos[1][3], h.Cos[3][1] = M, M, M, M
	h.Sin[0][3], h.Sin[3][0] = -M, -M
	h.Sin[1][2], h.Sin[2][1] = M, M
	c, err := circuit.NewBuilder("induction").Domain(exchange.Mechanical, "w").
		MachineModel("im", spec).
		VoltageWave("VA", "pa", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, E, f, 0, 0, 90}}).
		VoltageWave("VB", "pb", "0", circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{0, E, f, 0, 0, 0}}).
		Resistor("RSA", "pa", "sa", rs).Resistor("RSB", "pb", "sb", rs).
		Resistor("RRA", "ra", "0", r).Resistor("RRB", "rb", "0", r).
		Machine("N1", "im", [][2]string{{"sa", "0"}, {"sb", "0"}, {"ra", "0"}, {"rb", "0"}}, "w", "0").
		VoltageSource("VW", "w", "0", (1-s)*w).BuildTransient()
	if err != nil {
		t.Fatal(err)
	}
	stop := 0.6
	wf, err := Tran(c, TranSpec{Step: 2e-5, Stop: stop, MaxStep: 2e-5, UIC: true}, TranOptions{})
	if err != nil {
		t.Fatal(err)
	}
	ivw, _ := wf.ElementCurrent("VW")
	var sum, lo, hi float64
	n := 0
	lo, hi = math.Inf(1), math.Inf(-1)
	for i, tt := range wf.Time {
		if tt < stop-0.1 {
			continue
		}
		sum += ivw[i]
		lo, hi = math.Min(lo, ivw[i]), math.Max(hi, ivw[i])
		n++
	}
	T := sum / float64(n)
	t11 := E * E * s * r * xm * xm / (math.Pow(rs*r-s*(xs*xr-xm*xm), 2) + math.Pow(r*xs+s*rs*xr, 2))
	want := t11 / w
	if math.Abs(T-want) > 1e-4*want { // motoring: positive torque into the shaft
		t.Errorf("torque %.6g N·m, Eq. 21‑11 gives %.6g", T, want)
	}
	if ripple := (hi - lo) / want; ripple > 1e-2 {
		t.Errorf("torque ripple %.3g of mean in steady state", ripple)
	}
	t.Logf("torque %.6g N·m (velocity-source current), Eq. 21‑11: %.6g; ripple %.2g", T, want, (hi-lo)/want)
}
