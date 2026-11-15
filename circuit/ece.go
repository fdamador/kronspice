// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
	"strings"
)

// ECEForm selects the field–potential relation of the hypothetical
// spin-connection terms on R, L, C and K elements (register entry H48).
//
// Maxwell (form ECENone):  E = −∇φ − ∂A/∂t
// ECE  (UFT 455 fn. 4; Eckardt, ECE UFT Vol. 1, eq. 4.211):
//
//	E = −∇φ − ∂A/∂t − cω₀A + ωφ
//
// ECE2 (Vol. 1, eq. 6.189):
//
//	E = −∇φ − ∂A/∂t + 2(cω₀A − ωφ)
//
// so the ECE2 spin-connection terms are −2 × the ECE ones (Factor). In a
// circuit (a proposed mapping, not the papers'): integrating the ω₀ term
// along a coil's wire gives a torsion voltage s·cω₀·L·i, a resistance
// s·cω₀·L in series with every inductance, and s·cω₀·M between coupled
// coils; integrating the ω term along a lead from the terminal into the
// body gives an EMF s·κ·φ(terminal), κ = ∫ω·dl, so the body sees
// (1 + s·κ)·φ at each terminal. s = +1 (ECE) or −2 (ECE2). Ideal R and C
// have no series inductance, so they carry only the scalar term.
//
// UFT 456 (Eckardt & Goel 2024, eqs. 28–43) adds the field-built currents
// of ECE2 as conductivities, written with its wave numbers (not scaled by
// s): J = σ₀E with σ₀/ε₀ = ω_e = 2cκ₀ in a dielectric, and V = η₀H with
// η₀/μ₀ = ω_p = 2cκ_Λ0 in a magnetic path. By Kron's field-network
// analogy (ε₀ ↔ C, μ₀ ↔ L) they are a shunt conductance g_e = 2cκ₀·C
// across a capacitor's body (parameter k0) and a series resistance
// r_h = −2cκ_Λ0·L in a coil, −2cκ_Λ0·M between coupled coils (kl0). The
// sign is the paper's: κ_Λ0 > 0 makes r_h NEGATIVE, a gain (its "+η₀H²",
// eq. 57). Both are booked as model-supplied, as the torsion term is.
type ECEForm uint8

const (
	ECENone ECEForm = iota // Maxwell: no spin-connection terms (default)
	ECE                    // Einstein–Cartan–Evans field relation
	ECE2                   // ECE2 field relation
)

func (f ECEForm) String() string {
	switch f {
	case ECE:
		return "ECE"
	case ECE2:
		return "ECE2"
	}
	return "none"
}

// Factor is s, the multiplier of the spin-connection terms relative to the
// ECE form: +1 for ECE, −2 for ECE2, 0 for none.
func (f ECEForm) Factor() float64 {
	switch f {
	case ECE:
		return 1
	case ECE2:
		return -2
	}
	return 0
}

// ParseECEForm reads "ece", "ece2" or "none".
func ParseECEForm(s string) (ECEForm, bool) {
	switch strings.ToLower(s) {
	case "ece":
		return ECE, true
	case "ece2":
		return ECE2, true
	case "none", "maxwell":
		return ECENone, true
	}
	return ECENone, false
}

// ECESpec switches the hypothetical spin-connection terms on for a circuit
// and gives their defaults. Element parameters override them per element:
// w0 [1/m] on R, L, C and K; kappa, kappa_a, kappa_b [dimensionless] on R,
// L and C; k0 [1/m] on C; kl0 [1/m] on L and K. Netlist:
// ".ece form=ece2 w0=… kappa=… k0=… kl0=…". Hypothetical (H48): off unless
// declared, labelled in every report, and the power the terms deliver is
// booked per element as model-supplied, never as inflow.
type ECESpec struct {
	Form  ECEForm
	W0    float64 // scalar spin connection ω₀ [1/m]
	Kappa float64 // lead integral κ = ∫ω·dl from a terminal into the body
	K0    float64 // UFT 456 κ₀ [1/m]: shunt g_e = 2cκ₀·C on capacitors
	KL0   float64 // UFT 456 κ_Λ0 [1/m]: series r_h = −2cκ_Λ0·L on coils (and ·M on K)
}

// ECEParamKeys are the element parameters the spin-connection terms read.
var ECEParamKeys = []string{"w0", "kappa", "kappa_a", "kappa_b", "k0", "kl0", "w0skew"}

// ECETerms are one element's spin-connection terms, resolved.
type ECETerms struct {
	On     bool
	Form   ECEForm
	S      float64 // Form.Factor()
	W0     float64 // ω₀ [1/m]
	KappaA float64 // κ at Nodes[0] (0 on ground)
	KappaB float64 // κ at Nodes[1] (0 on ground)
	K0     float64 // UFT 456 κ₀ [1/m] (capacitors only)
	KL0    float64 // UFT 456 κ_Λ0 [1/m] (inductors and K only)
	// W0Skew [1/m] (K only) places mutual torsion antisymmetrically between
	// the two coils the K element declares: v1 += r_s·i2, v2 −= r_s·i1,
	// r_s = s·c·ω_skew·M. This is the placement of Kron's skew field tensor
	// F in the compound equation of motion p = r·ẋ + a·pẋ + F·ẋ (Tensors
	// for Circuits, 1959, ch. 31, eq. 31.5): ẋ·F·ẋ = 0, so it does no work.
	// It is also the form a rotating reference frame generates (L·Cᵀ·dC/dt
	// is skew for orthogonal C, eq. 31.14's inhomogeneous term), so it can
	// be transformed away; the symmetric W0 placement cannot.
	W0Skew float64
}

// AlphaA and AlphaB are the factors the body sees on each terminal's
// potential, 1 + s·κ.
func (t ECETerms) AlphaA() float64 { return 1 + t.S*t.KappaA }

// AlphaB is 1 + s·κ_b.
func (t ECETerms) AlphaB() float64 { return 1 + t.S*t.KappaB }

// Torsion is the series resistance s·c·ω₀·L of an inductance L [Ω].
func (t ECETerms) Torsion(l float64) float64 { return t.S * SpeedOfLight * t.W0 * l }

// Shunt is UFT 456's conductance g_e = 2cκ₀·C across a capacitance C [S]
// (σ₀/ε₀ = 2cκ₀, eqs. 33, 36).
func (t ECETerms) Shunt(c float64) float64 { return 2 * SpeedOfLight * t.K0 * c }

// TorsionSkew is the antisymmetric mutual torsion s·c·ω_skew·M of a
// mutual inductance M [Ω] (W0Skew).
func (t ECETerms) TorsionSkew(m float64) float64 { return t.S * SpeedOfLight * t.W0Skew * m }

// Homogeneous is UFT 456's series resistance r_h = −2cκ_Λ0·L of an
// inductance (or mutual inductance) L [Ω] (η₀/μ₀ = 2cκ_Λ0, eqs. 39, 42,
// with V = η₀H a magnetic conductivity of the opposite sign to a loss).
func (t ECETerms) Homogeneous(l float64) float64 { return -2 * SpeedOfLight * t.KL0 * l }

// Label names the terms for reports.
func (t ECETerms) Label() string {
	var parts []string
	if t.W0 != 0 {
		parts = append(parts, fmt.Sprintf("w0=%g 1/m", t.W0))
	}
	if t.KappaA != 0 {
		parts = append(parts, fmt.Sprintf("kappa_a=%g", t.KappaA))
	}
	if t.KappaB != 0 {
		parts = append(parts, fmt.Sprintf("kappa_b=%g", t.KappaB))
	}
	if t.K0 != 0 {
		parts = append(parts, fmt.Sprintf("k0=%g 1/m (UFT 456)", t.K0))
	}
	if t.KL0 != 0 {
		parts = append(parts, fmt.Sprintf("kl0=%g 1/m (UFT 456)", t.KL0))
	}
	if t.W0Skew != 0 {
		parts = append(parts, fmt.Sprintf("w0skew=%g 1/m (skew, workless)", t.W0Skew))
	}
	return fmt.Sprintf("%s spin connection (hypothetical, H48): %s", t.Form, strings.Join(parts, " "))
}

// SpeedOfLight is c [m/s].
const SpeedOfLight = 299792458.0

// ECEFor resolves element e's spin-connection terms. On is false when the
// circuit declares none, or when the element's terms are all zero (it then
// behaves exactly as without them).
func (c *Circuit) ECEFor(e Element) ECETerms {
	if c.ECE == nil || c.ECE.Form == ECENone {
		return ECETerms{}
	}
	switch e.Kind {
	case Resistor, Capacitor, Inductor, Mutual:
	default:
		return ECETerms{}
	}
	p := func(k string, def float64) float64 {
		if v, ok := e.Params[k]; ok {
			return v
		}
		return def
	}
	t := ECETerms{Form: c.ECE.Form, S: c.ECE.Form.Factor(), W0: p("w0", c.ECE.W0)}
	if e.Kind != Mutual {
		k := p("kappa", c.ECE.Kappa)
		t.KappaA, t.KappaB = p("kappa_a", k), p("kappa_b", k)
		if e.Nodes[0] == 0 {
			t.KappaA = 0
		}
		if e.Nodes[1] == 0 {
			t.KappaB = 0
		}
	}
	if e.Kind != Inductor && e.Kind != Mutual {
		t.W0 = 0 // no series inductance: no torsion
	}
	if e.Kind == Capacitor {
		t.K0 = p("k0", c.ECE.K0) // a dielectric: electric conductivity
	}
	if e.Kind == Inductor || e.Kind == Mutual {
		t.KL0 = p("kl0", c.ECE.KL0) // a magnetic path: magnetic conductivity
	}
	if e.Kind == Mutual {
		t.W0Skew = p("w0skew", 0) // per K element only: the K declares the partner coil
	}
	t.On = t.W0 != 0 || t.KappaA != 0 || t.KappaB != 0 || t.K0 != 0 || t.KL0 != 0 || t.W0Skew != 0
	return t
}

// validateECE checks the spin-connection declaration and parameters.
func (c *Circuit) validateECE() error {
	if c.ECE != nil {
		if c.ECE.Form > ECE2 {
			return fmt.Errorf("%w: .ece: unknown form", ErrBadValue)
		}
		if !finite(c.ECE.W0) || !finite(c.ECE.Kappa) || !finite(c.ECE.K0) || !finite(c.ECE.KL0) {
			return fmt.Errorf("%w: .ece: w0, kappa, k0 and kl0 must be finite", ErrBadValue)
		}
	}
	for _, e := range c.Elements {
		switch e.Kind {
		case Resistor, Capacitor, Inductor, Mutual, NonlinearCapacitor, NonlinearInductor:
		default:
			continue // other elements' parameters are their own (a P element's w0 is its ω0)
		}
		for _, k := range ECEParamKeys {
			v, ok := e.Params[k]
			if !ok {
				continue
			}
			switch {
			case c.ECE == nil:
				return fmt.Errorf("%w: %s: %s= needs a .ece declaration", ErrBadValue, e.Name, k)
			case e.Kind != Resistor && e.Kind != Capacitor && e.Kind != Inductor && e.Kind != Mutual:
				return fmt.Errorf("%w: %s: %s= is only for R, C, L and K", ErrBadValue, e.Name, k)
			case e.Kind == Mutual && k != "w0" && k != "kl0" && k != "w0skew":
				return fmt.Errorf("%w: %s: a K element takes only w0=, kl0= and w0skew=", ErrBadValue, e.Name)
			case k == "w0skew" && e.Kind != Mutual:
				return fmt.Errorf("%w: %s: w0skew= (skew mutual torsion) is only for K", ErrBadValue, e.Name)
			case k == "k0" && e.Kind != Capacitor:
				return fmt.Errorf("%w: %s: k0= (UFT 456 shunt) is only for C", ErrBadValue, e.Name)
			case k == "kl0" && e.Kind != Inductor && e.Kind != Mutual:
				return fmt.Errorf("%w: %s: kl0= (UFT 456 series term) is only for L and K", ErrBadValue, e.Name)
			case !finite(v):
				return fmt.Errorf("%w: %s: %s must be finite", ErrBadValue, e.Name, k)
			}
		}
	}
	return nil
}

func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

// ECE declares the hypothetical spin-connection terms for the circuit
// (see ECESpec). Per-element values go in Element.Params.
func (b *Builder) ECE(s ECESpec) *Builder {
	sp := s
	b.ece = &sp
	return b
}
