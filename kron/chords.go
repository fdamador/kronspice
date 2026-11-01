// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"fmt"
	"math"
	"strings"
)

// chord is an interaction term of the Tensorial Analysis of Networks (TAN):
// an off-diagonal entry z_ab(ω) of the primitive impedance tensor, "replaced
// by functions of any complexity" (P. Durand, O. Maurice, A. Reineix,
// "Generalized interaction principle implemented in the Kron's method,"
// Proc. WCE 2013; O. Maurice et al., arXiv:1412.0187). The voltage of
// branch a gains z·i_b; a reciprocal chord also gives branch b the term
// z·i_a. Mutual inductance is the classic chord, z = jωM; a coupling with
// a propagation delay τ is z = Z·e^(−jωτ), TAN's affine (delayed) source
// term seen in the frequency domain.
type chord struct {
	a, b       int
	z          func(omega float64) complex128
	reciprocal bool
}

// radiationEdge is a radiation-resistance edge: energy radiated away by an
// element, booked as power through that element's open path to the
// light/radiation exterior (CLAUDE.md §3a.1 item 5). It is a primitive
// branch of its own, in series with the element, whose resistance depends
// on frequency.
type radiationEdge struct {
	branch int // the edge's own branch
	elem   int // the radiating element
	r      func(omega float64) float64
}

// chordable reports whether a branch's relation has an impedance form a
// chord can add to (a branch whose current is impressed or controlled has
// no voltage equation).
func chordable(k branchKind) bool {
	switch k {
	case bI, bVCCS, bCCCS, bControl, bGyr1, bGyr2:
		return false
	}
	return true
}

// AddChord adds a TAN chord from element b's current into element a's
// voltage: E_a gains z(ω)·i_b (and E_b gains z(ω)·i_a when reciprocal).
// At DC the chord is z(0). Branches whose current is impressed (current
// sources, F, G) cannot carry a chord.
func (n *Network) AddChord(a, b string, z func(omega float64) complex128, reciprocal bool) error {
	ka, err := n.branchOf(a)
	if err != nil {
		return err
	}
	kb, err := n.branchOf(b)
	if err != nil {
		return err
	}
	if !chordable(n.branches[ka].kind) || (reciprocal && !chordable(n.branches[kb].kind)) {
		return fmt.Errorf("%w: chord between %s and %s: a current-source branch has no voltage equation", ErrUnsupported, a, b)
	}
	n.chords = append(n.chords, chord{a: ka, b: kb, z: z, reciprocal: reciprocal})
	return nil
}

// AddRadiation gives element name a radiation-resistance edge R_rad(ω) in
// series: a new primitive branch from a new internal junction, so the
// element's own relation is untouched and its radiated power is separate
// from its heat. The power through the edge leaves through the element's
// open path to the light/radiation exterior.
func (n *Network) AddRadiation(name string, r func(omega float64) float64) error {
	k, err := n.branchOf(name)
	if err != nil {
		return err
	}
	b := n.branches[k]
	switch b.kind {
	case bR, bC, bL:
	default:
		return fmt.Errorf("%w: radiation edge on %s: only R, C and L elements radiate here", ErrUnsupported, name)
	}
	x := n.nodes // the new junction between the element and its edge
	n.nodes++
	n.branches[k].to = x
	n.radiation = append(n.radiation, radiationEdge{branch: len(n.branches), elem: b.elem, r: r})
	n.branches = append(n.branches, branch{kind: bRadiation, elem: -1, from: x, to: b.to, ctrl: -1})
	n.rebuildGraph()
	return nil
}

// branchOf finds an element's branch by name.
func (n *Network) branchOf(name string) (int, error) {
	for ei, e := range n.c.Elements {
		if strings.EqualFold(e.Name, name) {
			if k := n.elemBranch[ei]; k >= 0 {
				return k, nil
			}
			return 0, fmt.Errorf("%w: %s is not a branch", ErrUnsupported, name)
		}
	}
	return 0, fmt.Errorf("%w: no element %q", ErrUnsupported, name)
}

func (n *Network) radiationR(k int, omega float64) float64 {
	for _, e := range n.radiation {
		if e.branch == k {
			return e.r(omega)
		}
	}
	return 0
}

// applyChords adds every chord to the branch relations. A relation
// ce·E_a + … = s in impedance form reads E_a = (s − …)/ce; adding z·i_b
// to E_a adds −ce·z·i_b to the relation. At DC a capacitor has ce = 0 (it
// is open) and a chord on it vanishes, as it should.
func (n *Network) applyChords(rel []relation, omega float64) {
	add := func(a, b int, z complex128) {
		ce, _ := rel[a].own(a)
		if ce == 0 || z == 0 {
			return
		}
		rel[a].i = append(rel[a].i, term{b, -ce * z})
	}
	for _, c := range n.chords {
		z := c.z(omega)
		add(c.a, c.b, z)
		if c.reciprocal {
			add(c.b, c.a, z)
		}
	}
}

// SmallLoopRadiation is the radiation resistance of an electrically small
// loop of n turns enclosing area A [m²]: R_rad = 320·π⁴·n²·A²/λ⁴ [Ω]
// (C. A. Balanis, Antenna Theory, 3rd ed., Wiley 2005, eqs 5-24 and 5-24a),
// with λ = 2πc/ω. Valid while the loop's circumference is below about
// λ/10 (same reference, §5.2).
func SmallLoopRadiation(area float64, turns int) func(omega float64) float64 {
	const c0 = 299792458.0
	return func(omega float64) float64 {
		if omega == 0 {
			return 0
		}
		lam := 2 * math.Pi * c0 / omega
		nA := float64(turns) * area
		return 320 * math.Pow(math.Pi, 4) * nA * nA / math.Pow(lam, 4)
	}
}
