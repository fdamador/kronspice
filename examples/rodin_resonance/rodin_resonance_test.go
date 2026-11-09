// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"math/cmplx"
	"testing"

	"github.com/fdamador/kronspice/em"
)

func rodinModel(t *testing.T) *Model {
	t.Helper()
	return NewModel(NewPair(Rodin, Circuits))
}

// The {12/5} path closes after 12 passes, goes 5 times round the axis and
// its centre field matches the paper's Table II at 1.60 and 1.25 A within
// the spread of the paper's own theory, FEM and measurement (±6 %).
func TestCentreFieldTableII(t *testing.T) {
	m := rodinModel(t)
	if d := dist(m.P.A[0], m.P.A[len(m.P.A)-1]); d > 1e-12 {
		t.Fatalf("path not closed: %g m", d)
	}
	// winding number round the axis: unwrap φ along one circuit
	turns := 0.0
	for i := 0; i+1 < len(m.P.A); i++ {
		a, b := m.P.A[i], m.P.A[i+1]
		d := math.Atan2(b[1], b[0]) - math.Atan2(a[1], a[0])
		d = math.Remainder(d, 2*math.Pi)
		turns += d / (2 * math.Pi)
	}
	if math.Abs(turns-5) > 1e-9 {
		t.Errorf("circuit goes %.6f times round the axis, want 5", turns)
	}
	bc := sub(m.FieldA(P3{}), m.FieldB(P3{}))[2]
	for _, r := range TableII[:2] {
		b := bc * r[0] * 1e3
		if math.Abs(b/r[3]-1) > 0.06 {
			t.Errorf("I = %.2f A: %.3f mT against measured %.2f", r[0], b, r[3])
		}
	}
	if !testing.Short() {
		best, _ := FitCircuits(10, 14)
		if best != Circuits {
			t.Errorf("least-squares fit gives %d circuits, code uses %d", best, Circuits)
		}
	}
}

// The segment method reproduces Kirchhoff's ring formula; the flux of the
// Biot–Savart field through a circle equals Neumann's mutual inductance.
func TestInductanceChecks(t *testing.T) {
	ring := em.Circle([3]float64{}, [3]float64{0, 0, 1}, 0.05, 200)
	pr := make([]P3, len(ring))
	copy(pr, ring)
	if e := LoopInductance(pr, 3e-4)/CircleRing(0.05, 3e-4) - 1; math.Abs(e) > 1e-3 {
		t.Errorf("ring inductance off by %.2e", e)
	}
	m := rodinModel(t)
	n, f := m.FluxCheck(P3{}, P3{0, 0, 1}, 0.01)
	if math.Abs(f/n-1) > 1e-3 {
		t.Errorf("axial flux check %.3e", f/n-1)
	}
	mab, mba := m.Reciprocity()
	if math.Abs(mab-mba) > 1e-9*math.Abs(mab) {
		t.Errorf("reciprocity %g vs %g", mab, mba)
	}
	if k := m.K(); k <= 0 || k >= 1 {
		t.Errorf("coupling %g", k)
	}
}

// Field energy outside a tube equals ½·L_ext·I² for the (2, 1) circuit.
func TestEnergyCheck(t *testing.T) {
	h := 1.5e-3
	if testing.Short() {
		h = 3e-3
	}
	q := NewPair(Knot21, 1)
	c, f := EnergyCheck(q.A, 2.5e-3, h)
	if math.Abs(f/c-1) > 0.03 {
		t.Errorf("field energy %g J against ½LI² %g J", f, c)
	}
}

// Quadrature drive: per-source books close, the power passed through M
// equals ω·M·I_A·I_B·sin φ, the ledger closes, and the 12-fold pattern
// stands still round the axis.
func TestQuadratureDrive(t *testing.T) {
	m := rodinModel(t)
	f := 1e3
	v := 1.25 * cmplx.Abs(complex(SourceOhms+m.RDC, 2*math.Pi*f*m.LA))
	for _, sense := range []int{Toroidal, Dipole} {
		r, err := m.Run(Drive{Name: "q", F: f, VA: v, VB: v, PhaseB: 90, Sense: sense, Cycles: 12, StepsPerCyc: 100})
		if err != nil {
			t.Fatal(err)
		}
		if d := r.PA + r.PB - r.Copper - r.SrcR; math.Abs(d) > 1e-3*(r.PA+r.PB) {
			t.Errorf("sense %d: books open by %g W", sense, d)
		}
		if math.Abs(r.Exchange-r.ExchangePhasor) > 0.02*math.Abs(r.ExchangePhasor) {
			t.Errorf("sense %d: exchange %g W, phasor formula %g W", sense, r.Exchange, r.ExchangePhasor)
		}
		if math.Abs(r.Mismatch) > 1e-5*r.EIn {
			t.Errorf("sense %d: ledger mismatch %g of %g J", sense, r.Mismatch, r.EIn)
		}
		adv, _, _ := m.RingTravel(r.IA, r.IB, 120, 48)
		if math.Abs(adv) > 1e-6 {
			t.Errorf("sense %d: pattern advanced %g°", sense, adv)
		}
	}
}

// The tuned winding builds up to Q = ω0·L/R (from the energy rise time).
func TestResonantBuildUp(t *testing.T) {
	m := rodinModel(t)
	f0 := 1e5
	r0 := SourceOhms + WireResistance(m.LenA, f0)
	r, err := m.Run(Drive{Name: "t", F: f0, VA: 1.25 * r0, Tuned: true, Cycles: 150, StepsPerCyc: 60})
	if err != nil {
		t.Fatal(err)
	}
	q := math.Pi * f0 * riseTime(r, f0)
	want := 2 * math.Pi * f0 * m.LA / r0
	if math.Abs(q/want-1) > 0.1 {
		t.Errorf("Q from build-up %.1f, ω0L/R %.1f", q, want)
	}
	if a := cmplx.Abs(r.IA); math.Abs(a-1.25) > 0.02 {
		t.Errorf("tuned current %.4f A, want 1.25", a)
	}
}

// The complex AC solver of the sectional net agrees with a kronspice
// transient of the same (reduced) net, N element and all.
func TestACTwin(t *testing.T) {
	m := NewModel(NewPairSec(Rodin, Circuits, 2))
	f := 1e5
	n := m.NewNet(DrivenPair, Dipole, f)
	n.SetFreq(f)
	r, err := n.SolveAC(f, n.DriveSources(Dipole, 3.4), nil, nil, Separate, 1)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := m.TranVerify(n, Layout{F: f, Sense: Dipole}, false, 40, 100, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 2; i++ {
		if e := cmplx.Abs(tr.ISrc[i]-r.ISrc[i]) / cmplx.Abs(r.ISrc[i]); e > 3e-3 {
			t.Errorf("source %d: transient vs AC %.2e", i, e)
		}
	}
	if d := r.PSrc[0] + r.PSrc[1] - r.PCopper - r.PSrcR; math.Abs(d) > 1e-9*(math.Abs(r.PSrc[0])+math.Abs(r.PSrc[1])) {
		t.Errorf("AC books open by %g W", d)
	}
}

// The crossed (series counter-wound) ring has 720° modes, for the Rodin
// path and for a (2, 1) winding alike.
func TestCrossedRing720(t *testing.T) {
	for _, kind := range []Kind{Rodin, Knot21} {
		nc := Circuits
		if testing.Short() {
			nc = 3
		}
		m := NewModel(NewPair(kind, nc))
		n := m.NewNet(CrossedRing, 0, 1e6)
		f, v, err := n.Modes()
		if err != nil {
			t.Fatal(err)
		}
		i := 0
		for f[i] < 1e3 {
			i++
		}
		if p := n.TripParity(v, i); p > -0.9 {
			t.Errorf("%v: first crossed-ring mode parity %.3f, want ≈ −1", kind, p)
		}
	}
}

// A layout with groups: the AC books close, and the kronspice transient
// with bridges and bank closes its books (sources = bank + losses +
// stored, to the integration mismatch).
func TestGroupsLedger(t *testing.T) {
	m := NewModel(NewPairSec(Rodin, Circuits, 2))
	p3 := NewPart3(m)
	l := Layout{F: 2e4, Sense: Toroidal, NG: 2, Theta0: 0.1, DTh: math.Pi,
		Shape: GroupShape{Rho: TorusR0, Z: 0, Beta: 0, Radius: m.P.OuterShell() + 4e-3, Turns: 60}, Comb: Separate, RLoad: 50}
	if c := m.Feasible(l); c < 0 {
		t.Fatalf("test layout infeasible by %g m", -c)
	}
	pw, r, err := p3.Score(l)
	if err != nil {
		t.Fatal(err)
	}
	ps := r.PSrc[0] + r.PSrc[1]
	if d := ps - pw - r.PCopper - r.PSrcR; math.Abs(d) > 1e-9*ps {
		t.Errorf("AC books open by %g W", d)
	}
	cyc := 120
	if testing.Short() {
		cyc = 40
	}
	tr, err := m.TranVerify(p3.net(l.Sense), l, true, cyc, 40, 1e-6)
	if err != nil {
		t.Fatal(err)
	}
	es := tr.ESrc[0] + tr.ESrc[1]
	books := es - tr.EBankEnd - tr.ECopper - tr.ESrcR - tr.EDiode - tr.EBleed - tr.ECoilsEnd
	if math.Abs(books-tr.Mismatch) > 1e-6*tr.EIn || math.Abs(tr.Mismatch) > 1e-3*tr.EIn {
		t.Errorf("transient books %g J, ledger mismatch %g J of %g J", books, tr.Mismatch, tr.EIn)
	}
	if tr.EBankEnd <= 0 || tr.EBankEnd > es {
		t.Errorf("bank %g J against sources %g J", tr.EBankEnd, es)
	}
}

// The toroidal secondary (Part 3t): its inductance tables reproduce a
// direct computation, the current-sheet toroid and an explicit helix; the
// AC books close; 12 alternating sections in series cancel by symmetry;
// the kronspice transient with bridges and bank closes its books.
func TestToroidalSecondary(t *testing.T) {
	kit, err := NewTorKit(NewPair(Rodin, Circuits))
	if err != nil {
		t.Fatal(err)
	}
	c := kit.Checks(rodinModel(t), 114, 8, 20)
	if c.TableErr > 1e-3 {
		t.Errorf("table interpolation error %.2e", c.TableErr)
	}
	if e := c.TorRings/c.TorSheet - 1; math.Abs(e) > 5e-3 {
		t.Errorf("uniform toroid %.4g H vs sheet %.4g H", c.TorRings, c.TorSheet)
	}
	if math.Abs(c.RingMA/c.HelMA-1) > 1e-2 || math.Abs(c.RingMB/c.HelMB-1) > 1e-2 || math.Abs(c.RingL/c.HelL-1) > 5e-2 {
		t.Errorf("rings vs helix: M_A %g/%g, M_B %g/%g, L %g/%g", c.RingMA, c.HelMA, c.RingMB, c.HelMB, c.RingL, c.HelL)
	}
	m := NewModel(NewPairSec(Rodin, Circuits, 2))
	tp := NewTorPart(kit, m)
	l := TorLayout{F: 4e3, Sense: Toroidal, Theta0: 0.1, Arc: 25 * math.Pi / 180, Turns: 12, Comb: Series, RLoad: 10}
	pw, r, err := tp.Score(l, 0)
	if err != nil {
		t.Fatal(err)
	}
	ps := r.PSrc[0] + r.PSrc[1]
	if d := ps - pw - r.PCopper - r.PSrcR; math.Abs(d) > 1e-9*ps || pw <= 0 || pw >= ps {
		t.Errorf("AC books open by %g W (load %g of %g W)", d, pw, ps)
	}
	alt := l
	alt.Alt = true
	if pa, _, err := tp.Score(alt, 0); err != nil || pa > 1e-12*pw {
		t.Errorf("alternating series sections give %g W against %g W (%v)", pa, pw, err)
	}
	if testing.Short() {
		return
	}
	n := tp.net(l.Sense)
	n.SetFreq(l.F)
	gs, mg, err := kit.Groups(m, n, l)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := TranGroups(n, l.F, l.Sense, 0, gs, mg, l.Comb, l.RLoad, 30, 30, 1e-4)
	if err != nil {
		t.Fatal(err)
	}
	es := tr.ESrc[0] + tr.ESrc[1]
	books := es - tr.EBankEnd - tr.ECopper - tr.ESrcR - tr.EDiode - tr.EBleed - tr.ECoilsEnd
	if math.Abs(books-tr.Mismatch) > 1e-6*tr.EIn || math.Abs(tr.Mismatch) > 2e-3*tr.EIn {
		t.Errorf("transient books %g J, ledger mismatch %g J of %g J", books, tr.Mismatch, tr.EIn)
	}
	if tr.EBankEnd <= 0 || tr.EBankEnd > es {
		t.Errorf("bank %g J against sources %g J", tr.EBankEnd, es)
	}
}
