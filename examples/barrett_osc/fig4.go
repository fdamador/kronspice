// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"strings"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Fig. 4 of Barrett 1991: coil networks said to carry the algebra of
// quaternions (A, "Q: αβ ≠ βα"), Cayley numbers (B, "H: αβ ≠ βα,
// (αβ)γ ≠ γ(αβ)") and the next Cayley–Dickson algebra (C, "J: … αx + β = 0
// but x not unique"). Read from the drawing (U = top of a loop, N = nodes
// down the stem, open circles = terminals):
//
//	A      γ, α both U–N1 (a closed loop on N1); β N1–T
//	B(i)   loop γ, α on N1; β, δ both N1–N2; ε N2–T
//	B(ii)  loop γ, α on N1; β N1–N; δ N–T1; N–T2
//	C(i)   loop γ, α on N1; β, δ N1–N2; ε, φ N2–N3; ψ N3–T
//	C(ii)  loop γ, α on N1; β N1–N; δ N–T1; ε T1–T3 (open); N–T2
//
// Coils have 1 mΩ each (Q ≈ 3·10⁵ at 5 MHz; ideal coils in a closed loop
// make the DC operating point singular). Values are declared in Fig4Params. A, B(i) and C(i)
// have one terminal; B(ii) and C(ii) have two (with G1 as reference they
// are two-ports, port 1 at T2, port 2 at T1). Every node has a small
// capacitance to G1 (declared); an open terminal has a larger one.
type Fig4Params struct {
	Alpha, Beta, Gam, Delta, Eps, Phi, Psi float64 // [H]
	Cnode, Cterm                           float64 // node and terminal capacitance to G1 [F]
	ESR                                    float64 // series resistance of each coil [Ω]
	F                                      float64 // test frequency [Hz]
	R0                                     float64 // reference resistance for S [Ω]
}

// DefaultFig4: 10…70 µH, 1 mΩ each, 1 pF per node, 10 pF per open
// terminal, 5 MHz.
var DefaultFig4 = Fig4Params{Alpha: 10e-6, Beta: 20e-6, Gam: 30e-6, Delta: 40e-6, Eps: 50e-6, Phi: 60e-6, Psi: 70e-6,
	Cnode: 1e-12, Cterm: 10e-12, ESR: 1e-3, F: 5e6, R0: 50}

// net4 adds one Fig. 4 network with element names prefixed by pre. in and
// out name the two terminals of B(ii)/C(ii) (T2, T1); one-terminal networks
// use in only. Every internal node gets Cnode to ground; open terminals
// Cterm.
func (p Fig4Params) net4(b *circuit.Builder, v, pre, in, out string) {
	n := func(s string) string { return pre + s }
	k := 0
	coil := func(name, a, c string, h float64) { // with its resistance
		k++
		m := fmt.Sprintf("%sr%d", pre, k)
		b.Inductor(n(name), a, m, h)
		b.Resistor(n("R"+name), m, c, p.ESR)
	}
	cn := func(s string) {
		if p.Cnode > 0 {
			b.Capacitor(n("C"+s), n(s), "0", p.Cnode)
		}
	}
	loop := func() { // γ and α from U to N1
		coil("Lgam", n("U"), n("N1"), p.Gam)
		coil("Lalp", n("U"), n("N1"), p.Alpha)
		cn("U")
		cn("N1")
	}
	switch v {
	case "A":
		loop()
		coil("Lbet", n("N1"), in, p.Beta)
	case "B(i)":
		loop()
		coil("Lbet", n("N1"), n("N2"), p.Beta)
		coil("Ldel", n("N1"), n("N2"), p.Delta)
		cn("N2")
		coil("Leps", n("N2"), in, p.Eps)
	case "C(i)":
		loop()
		coil("Lbet", n("N1"), n("N2"), p.Beta)
		coil("Ldel", n("N1"), n("N2"), p.Delta)
		cn("N2")
		coil("Leps", n("N2"), n("N3"), p.Eps)
		coil("Lphi", n("N2"), n("N3"), p.Phi)
		cn("N3")
		coil("Lpsi", n("N3"), in, p.Psi)
	case "B(ii)", "C(ii)":
		loop()
		coil("Lbet", n("N1"), in, p.Beta) // in = N = T2
		coil("Ldel", in, out, p.Delta)    // out = T1
		if v == "C(ii)" {
			coil("Leps", out, n("T3"), p.Eps)
			b.Capacitor(n("CT3"), n("T3"), "0", p.Cterm)
		}
	}
}

// OneTerminal drives the terminal of A, B(i) or C(i) against G1 and returns
// the terminal admittance relative to the terminal's own capacitance alone,
// Y/(jωCterm), and the largest coil current relative to the terminal current.
func (p Fig4Params) OneTerminal(v string) (yRel complex128, coilMax float64, err error) {
	b := circuit.NewBuilder("Fig. 4" + v)
	b.Add(circuit.Element{Name: "VD", Kind: circuit.VoltageSource, ACMag: 1}, "t", "0")
	b.Capacitor("CT", "t", "0", p.Cterm)
	p.net4(b, v, "x", "t", "")
	c, err := b.Build()
	if err != nil {
		return 0, 0, err
	}
	r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: p.F, FStop: p.F}, analysis.Options{})
	if err != nil {
		return 0, 0, err
	}
	iv, _ := r.ElementCurrent("VD")
	y := -iv[0] // current delivered into the terminal for 1 V
	w := 2 * math.Pi * p.F
	yRel = y / complex(0, w*p.Cterm)
	for k, name := range r.ElementNames {
		if strings.HasPrefix(name, "xL") {
			coilMax = math.Max(coilMax, cmplx.Abs(r.I[k][0])/cmplx.Abs(y))
		}
	}
	return yRel, coilMax, nil
}

// ABCD is a two-port's chain matrix [[A, B], [C, D]]: V1 = A·V2 + B·I2,
// I1 = C·V2 + D·I2, with I2 leaving port 2.
type ABCD [2][2]complex128

// Mul is the cascade a then b.
func (a ABCD) Mul(b ABCD) ABCD {
	var m ABCD
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			m[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j]
		}
	}
	return m
}

// Det is A·D − B·C (1 for a reciprocal network).
func (a ABCD) Det() complex128 { return a[0][0]*a[1][1] - a[0][1]*a[1][0] }

// Diff is the largest relative element difference.
func (a ABCD) Diff(b ABCD) float64 {
	d, s := 0.0, 0.0
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			d = math.Max(d, cmplx.Abs(a[i][j]-b[i][j]))
			s = math.Max(s, cmplx.Abs(a[i][j]))
		}
	}
	return d / s
}

// S returns the scattering matrix for reference resistance r0.
func (a ABCD) S(r0 float64) [2][2]complex128 {
	A, B, C, D := a[0][0], a[0][1], a[1][0], a[1][1]
	z := complex(r0, 0)
	den := A + B/z + C*z + D
	return [2][2]complex128{
		{(A + B/z - C*z - D) / den, 2 * (A*D - B*C) / den},
		{2 / den, (-A + B/z - C*z + D) / den},
	}
}

// Unitarity is max |(S^H S − I)_ij| (0 for a lossless two-port).
func Unitarity(s [2][2]complex128) float64 {
	m := 0.0
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			var acc complex128
			for k := 0; k < 2; k++ {
				acc += cmplx.Conj(s[k][i]) * s[k][j]
			}
			if i == j {
				acc--
			}
			m = math.Max(m, cmplx.Abs(acc))
		}
	}
	return m
}

// Section is one two-port network with its own values.
type Section struct {
	Kind string // "B(ii)" or "C(ii)"
	P    Fig4Params
}

// Chain measures, with kronspice AC, the ABCD matrix of the cascade of
// sections in order (port 2 of each to port 1 of the next): A, C with
// port 2 open, B, D with port 2 shorted.
func Chain(f float64, secs ...Section) (ABCD, error) {
	var m ABCD
	for _, short := range []bool{false, true} {
		b := circuit.NewBuilder("Fig. 4 cascade")
		b.Add(circuit.Element{Name: "V1", Kind: circuit.VoltageSource, ACMag: 1}, "p0", "0")
		for k, s := range secs {
			s.P.net4(b, s.Kind, fmt.Sprintf("s%d", k), fmt.Sprintf("p%d", k), fmt.Sprintf("p%d", k+1))
		}
		last := fmt.Sprintf("p%d", len(secs))
		if short {
			b.Add(circuit.Element{Name: "V2", Kind: circuit.VoltageSource}, last, "0")
		}
		c, err := b.Build()
		if err != nil {
			return m, err
		}
		r, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 1, FStart: f, FStop: f}, analysis.Options{})
		if err != nil {
			return m, err
		}
		i1v, _ := r.ElementCurrent("V1")
		i1 := -i1v[0]
		if !short {
			v2, _ := r.NodeVoltage(last)
			m[0][0], m[1][0] = 1/v2[0], i1/v2[0]
		} else {
			i2v, _ := r.ElementCurrent("V2")
			i2 := i2v[0] // into V2's + terminal = leaving port 2
			m[0][1], m[1][1] = 1/i2, i1/i2
		}
	}
	return m, nil
}

// Fig4Sections are three different two-ports from the drawing: X = B(ii),
// Y = C(ii), Z = B(ii) with its coil values reversed.
func Fig4Sections() (x, y, z Section) {
	p := DefaultFig4
	q := p
	q.Alpha, q.Beta, q.Gam, q.Delta = p.Delta, p.Gam, p.Beta, p.Alpha
	return Section{"B(ii)", p}, Section{"C(ii)", p}, Section{"B(ii)", q}
}

// Fig4Result collects the two-port checks.
type Fig4Result struct {
	CascadeVsProduct     float64 // kronspice cascade XY vs ABCD(X)·ABCD(Y)
	Commutator           float64 // |XY − YX| relative
	Associator           float64 // kronspice XYZ vs (XY)Z and X(YZ), largest
	DetErr, UnitarityErr float64 // largest |det − 1| and S unitarity error over X, Y, Z, XY
}

// TwoPorts runs the two-port checks.
func TwoPorts() (Fig4Result, error) {
	var r Fig4Result
	x, y, z := Fig4Sections()
	f := DefaultFig4.F
	X, err := Chain(f, x)
	if err != nil {
		return r, err
	}
	Y, err := Chain(f, y)
	if err != nil {
		return r, err
	}
	Z, err := Chain(f, z)
	if err != nil {
		return r, err
	}
	XY, err := Chain(f, x, y)
	if err != nil {
		return r, err
	}
	YX, err := Chain(f, y, x)
	if err != nil {
		return r, err
	}
	XYZ, err := Chain(f, x, y, z)
	if err != nil {
		return r, err
	}
	r.CascadeVsProduct = math.Max(XY.Diff(X.Mul(Y)), YX.Diff(Y.Mul(X)))
	r.Commutator = XY.Diff(YX)
	r.Associator = math.Max(XYZ.Diff(X.Mul(Y).Mul(Z)), XYZ.Diff(X.Mul(Y.Mul(Z))))
	for _, m := range []ABCD{X, Y, Z, XY} {
		r.DetErr = math.Max(r.DetErr, cmplx.Abs(m.Det()-1))
		r.UnitarityErr = math.Max(r.UnitarityErr, Unitarity(m.S(DefaultFig4.R0)))
	}
	return r, nil
}

func printFig4() error {
	p := DefaultFig4
	fmt.Println("Fig. 4 — the coil networks A (quaternion), B (Cayley), C (sedenion-like)")
	fmt.Printf("  α…ψ = %.0f…%.0f µH, %.0f pF per node, %.0f pF per open terminal, %.0f MHz (declared)\n",
		p.Alpha*1e6, p.Psi*1e6, p.Cnode*1e12, p.Cterm*1e12, p.F/1e6)
	fmt.Println("  one-terminal networks, driven at their terminal against G1:")
	for _, v := range []string{"A", "B(i)", "C(i)"} {
		y, cm, err := p.OneTerminal(v)
		if err != nil {
			return err
		}
		fmt.Printf("    %-5s  Y/(jωCterm) = %.6f%+.2ei   largest coil current %.2e of the terminal current\n", v, real(y), imag(y), cm)
	}
	q := p
	q.Cnode = 0
	for _, v := range []string{"A", "B(i)", "C(i)"} {
		_, cm, err := q.OneTerminal(v)
		if err != nil {
			return err
		}
		if cm > 1e-12 {
			return fmt.Errorf("Fig. 4%s: coil current %.2g without node capacitance", v, cm)
		}
	}
	fmt.Println("    with no capacitance at the inner nodes every coil current is 0: a network with one terminal")
	fmt.Println("    is only its terminal's capacitance; the coils, and so the α, β, γ 'products', act only through")
	fmt.Println("    stray capacitance (or as antennas in an outside field).")
	r, err := TwoPorts()
	if err != nil {
		return err
	}
	fmt.Println("  two-port networks B(ii) = X, C(ii) = Y, B(ii) with values reversed = Z (ports T2 → T1):")
	fmt.Printf("    cascade XY in kronspice vs ABCD(X)·ABCD(Y)        %.1e\n", r.CascadeVsProduct)
	fmt.Printf("    XY vs YX (order matters: 'αβ ≠ βα')               %.3f\n", r.Commutator)
	fmt.Printf("    XYZ vs (XY)Z and X(YZ) (associative)              %.1e\n", r.Associator)
	fmt.Printf("    |det ABCD − 1| (reciprocal: invertible, no zero divisors) %.1e\n", r.DetErr)
	fmt.Printf("    S†S − I (lossless: S unitary, an element of U(2))  %.1e\n", r.UnitarityErr)
	fmt.Println("  order dependence is real and is ordinary two-port theory: cascades multiply as 2×2 matrices,")
	fmt.Println("  which do not commute. They do associate, and a reciprocal network's matrix has determinant 1,")
	fmt.Println("  so cascades can carry neither the Cayley numbers' non-associativity nor the zero divisors of")
	fmt.Println("  the box J. (Box H's '(αβ)γ ≠ γ(αβ)' is a commutator, not an associator.)")
	fmt.Println()
	return nil
}
