// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"
	"strconv"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/fluid"
	"github.com/fdamador/kronspice/netlist"
)

// A mass–spring–damper in the mechanical domain (mobility analogy: mass =
// capacitance m, spring = inductance 1/k, damper = resistance 1/b; SPICE
// names start with the element letter: CM, LK, RB) driven
// by a 1 N force: the velocity is F/(jωm + b + k/(jω)), peaking at F/b at
// ω0 = √(k/m); the damper's heat leaves through the thermal port and the
// force source draws on the mechanical exterior.
func TestMassSpringDamper(t *testing.T) {
	m, k, b := 0.02, 2e4, 0.5 // kg, N/m, N·s/m
	c, err := circuit.NewBuilder("msd").Domain(exchange.Mechanical, "v").
		Add(circuit.Element{Name: "IF", Kind: circuit.CurrentSource, ACMag: 1}, "0", "v").
		Capacitor("CM", "v", "0", m).Inductor("LK", "v", "0", 1/k).Resistor("RB", "v", "0", 1/b).Build()
	if err != nil {
		t.Fatal(err)
	}
	w0 := math.Sqrt(k / m)
	f0 := w0 / (2 * math.Pi)
	r, err := AC(c, ACSpec{Variation: "lin", Points: 3, FStart: f0 / 2, FStop: 1.5 * f0}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	v, _ := r.NodeVoltage("v")
	for i, f := range r.Freq {
		w := 2 * math.Pi * f
		want := 1 / complex(b, w*m-k/w)
		if d := cmplx.Abs(v[i] - want); d > 1e-12*cmplx.Abs(want) {
			t.Errorf("%g Hz: v = %v m/s, want %v", f, v[i], want)
		}
	}
	if math.Abs(cmplx.Abs(v[1])-1/b) > 1e-12 {
		t.Errorf("velocity at resonance %g m/s, want F/b = %g", cmplx.Abs(v[1]), 1/b)
	}
	// Average power: F·v/2 in, b·|v|²/2 as heat.
	pin := 0.5 * real(v[1])
	if math.Abs(r.Ledger.Inflow[1]-pin) > 1e-12*pin || math.Abs(r.Ledger.Dissipation[1]-pin) > 1e-12*pin {
		t.Errorf("ledger in %g, heat %g, want %g W", r.Ledger.Inflow[1], r.Ledger.Dissipation[1], pin)
	}
}

// A piezoelectric transducer: clamped capacitance C0 across a gyrator of
// force factor φ onto a mechanical resonator. Its electrical admittance
// is the Butterworth–Van Dyke circuit, Y = jωC0 + 1/(R1 + jωL1 + 1/(jωC1))
// with L1 = m/φ², C1 = φ²/k, R1 = b/φ² (Beranek & Mellow 2012, Ch. 3;
// IEEE Std 176-1987); the power taken from the electrical source equals the
// power dissipated in the mechanical damper, and the transducer itself
// absorbs nothing.
func TestPiezoButterworthVanDyke(t *testing.T) {
	C0, phi, m, k, b := 2e-9, 0.8, 5e-3, 1.2e8, 4.0
	src := `piezo
.domain mechanical v
V1 e 0 DC 0 AC 1
C0 e 0 2n
Y1 e 0 v 0 0.8
CM v 0 5m
LK v 0 {1/k}
RB v 0 {1/b}
.end
`
	src = strings.NewReplacer("{1/k}", strconv.FormatFloat(1/k, 'g', -1, 64), "{1/b}", strconv.FormatFloat(1/b, 'g', -1, 64)).Replace(src)
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	if c.NodeDomain(c.NodeIndex("v")) != exchange.Mechanical {
		t.Fatal("v is not mechanical")
	}
	f0 := math.Sqrt(k/m) / (2 * math.Pi)
	r, err := AC(c, ACSpec{Variation: "lin", Points: 5, FStart: 0.9 * f0, FStop: 1.1 * f0}, Options{})
	if err != nil {
		t.Fatal(err)
	}
	iv, _ := r.ElementCurrent("V1")
	for i, f := range r.Freq {
		w := 2 * math.Pi * f
		L1, C1, R1 := m/(phi*phi), phi*phi/k, b/(phi*phi)
		Y := complex(0, w*C0) + 1/complex(R1, w*L1-1/(w*C1))
		if d := cmplx.Abs(-iv[i] - Y); d > 1e-12*cmplx.Abs(Y) {
			t.Errorf("%g Hz: admittance %v, BVD %v", f, -iv[i], Y)
		}
		if l := r.Ledger; math.Abs(l.Inflow[i]-l.Dissipation[i]) > 1e-12*l.Inflow[i] || l.Inflow[i] <= 0 {
			t.Errorf("%g Hz: electrical in %g W, mechanical heat %g W", f, l.Inflow[i], l.Dissipation[i])
		}
	}
	op, err := OP(c, Options{})
	if err != nil {
		t.Fatal(err)
	}
	y, _ := op.Element("Y1")
	if y.PowerAbsorbed != 0 {
		t.Errorf("transducer absorbs %g W at DC", y.PowerAbsorbed)
	}
}

// Electrical → mechanical → acoustic: a piezo disc radiating into water
// through a piston (gyrator of area S) loaded by the piston's radiation
// network. The acoustic radiation resistance's power leaves through the
// acoustic exterior and equals the analytic radiated power ½·Re(Z)·|U|²;
// the ledger balances; H₂O and D₂O loads differ as their ρc does.
func TestPiezoRadiatesIntoWater(t *testing.T) {
	a := 0.01
	for _, med := range []fluid.Medium{fluid.Water, fluid.HeavyWater} {
		Ma, Rpar, fmax, err := fluid.PistonRadiation(a, med)
		if err != nil {
			t.Fatal(err)
		}
		S := math.Pi * a * a
		c, err := circuit.NewBuilder("piston").Domain(exchange.Mechanical, "v").Domain(exchange.Acoustic, "p").
			Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 10}, "e", "0").
			Capacitor("C0", "e", "0", 2e-9).Gyrator("Y1", "e", "0", "v", "0", 0.8).
			Capacitor("CM", "v", "0", 2e-3).Inductor("LK", "v", "0", 1/5e7).Resistor("RB", "v", "0", 1/2.0).
			Gyrator("Y2", "v", "0", "p", "0", S).
			Inductor("LMA", "p", "0", Ma).ResistorTo("RA", "p", "0", Rpar, exchange.Acoustic).Build()
		if err != nil {
			t.Fatal(err)
		}
		f := math.Min(math.Sqrt(5e7/2e-3)/(2*math.Pi), 0.4*fmax)
		r, err := AC(c, ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, Options{})
		if err != nil {
			t.Fatal(err)
		}
		p, _ := r.NodeVoltage("p")
		ua, _ := r.ElementCurrent("RA")
		um, _ := r.ElementCurrent("LMA")
		U := ua[0] + um[0]
		z := fluid.PistonImpedance(a, med, f)
		radiated := 0.5 * real(p[0]*cmplx.Conj(ua[0]))
		analytic := 0.5 * real(z) * cmplx.Abs(U) * cmplx.Abs(U)
		ka := 2 * math.Pi * f * a / med.SoundSpeed
		if math.Abs(radiated-analytic) > 2*ka*ka*analytic {
			t.Errorf("%s: radiated %g W, ½·Re(Z)·|U|² = %g W (ka = %.3f)", med.Name, radiated, analytic, ka)
		}
		l := r.Ledger
		if math.Abs(l.Outflow[0]-radiated) > 1e-12*radiated || math.Abs(l.Residual[0]) > 1e-12*l.Inflow[0] {
			t.Errorf("%s: ledger in %g, heat %g, acoustic out %g, residual %g", med.Name, l.Inflow[0], l.Dissipation[0], l.Outflow[0], l.Residual[0])
		}
		t.Logf("%s at %.1f kHz (ka = %.3f): %.3g W in, %.3g W heat, %.3g W radiated", med.Name, f/1e3, ka,
			l.Inflow[0], l.Dissipation[0], l.Outflow[0])
	}
}

// Elements that join two domains (other than a transducer) and
// semiconductors outside the electrical domain are rejected.
func TestDomainValidation(t *testing.T) {
	ok := func() *circuit.Builder {
		return circuit.NewBuilder("d").Domain(exchange.Mechanical, "v").VoltageSource("V1", "e", "0", 1).
			Resistor("R1", "e", "0", 1).Gyrator("Y1", "e", "0", "v", "0", 1).Inductor("LK", "v", "0", 1)
	}
	if _, err := ok().Build(); err != nil {
		t.Fatal(err)
	}
	for _, b := range []*circuit.Builder{
		ok().Resistor("R2", "e", "v", 1),
		ok().Model("dm", "d", nil).Diode("D1", "v", "0", "dm"),
		ok().Gyrator("Y2", "e", "v", "v", "0", 1),
		ok().ResistorTo("R3", "v", "0", 1, exchange.Neutrino),
		ok().Domain(exchange.Light, "q").Resistor("R4", "q", "0", 1),
	} {
		if _, err := b.Build(); !errors.Is(err, circuit.ErrBadValue) {
			t.Errorf("want ErrBadValue, got %v", err)
		}
	}
}

// Netlists with domains, transducers and EXT= round-trip.
func TestDomainNetlistRoundTrip(t *testing.T) {
	src := "rt\n.domain mechanical v\n.domain acoustic p\nV1 e 0 DC 1 AC 1\nC0 e 0 1n\nY1 e 0 v 0 0.5\nLK v 0 1u\nCM v 0 1m\nY2 v 0 p 0 3e-4\nRA p 0 1e6 EXT=acoustic\nLA p 0 1\n.end\n"
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := circuit.FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	d2, err := netlist.ParseString(sb.String())
	if err != nil {
		t.Fatalf("%v\n%s", err, sb.String())
	}
	c2, err := circuit.FromDeck(d2)
	if err != nil {
		t.Fatalf("%v\n%s", err, sb.String())
	}
	for i, n := range c.NodeNames {
		if c2.NodeDomain(c2.NodeIndex(n)) != c.NodeDomain(i) {
			t.Errorf("node %s domain lost:\n%s", n, sb.String())
		}
	}
	if e := c2.Elements[c2.Index("RA")]; e.Channel != exchange.Acoustic {
		t.Errorf("RA channel %v", e.Channel)
	}
	if e := c2.Elements[c2.Index("Y2")]; e.Kind != circuit.Gyrator || e.Value != 3e-4 {
		t.Errorf("Y2 = %+v", e)
	}
}
