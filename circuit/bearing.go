// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"fmt"
	"math"
)

// BearingMotor is the ball-bearing (Huber-effect) motor, kronspice P:
//
//	Pname a b m+ m- R KAPPA=κ W0=ω0 [VF=v] [THERMAL=1 ETAMAX=η]
//
// A current I through a rotating ball bearing (a → b, electrical) produces a
// torque on the shaft (port m+ m−, mechanical: across = angular velocity ω,
// through = torque) that aids the rotation and scales as I², with no torque
// at rest (J. Huber, Elektrotechnik und Maschinenbau 76, 169–174, 1959;
// R. A. Milroy, J. Appl. Mech. 34, 525, 1967; H. Gruenberg, Am. J. Phys.
// 46, 1213, 1978; J. L. Choo, W. L. Soong, D. Abbott, Proc. SPIE 5649,
// 700, 2005, torque ∝ I² per K. T. McDonald, "Ball-Bearing Motor",
// Princeton, 2020):
//
//	T = κ·I²·ψ(ω),   ψ(ω) = tanh(ω/ω0)                  (EM mode)
//	V = VF·tanh(I/1 mA) + R·I + κ·ψ(ω)·ω·I               back-EMF
//
// VF is the bearing's measured, almost constant DC contact voltage (0.5–1.5
// V for greased bearings; O. Safdarzadeh, A. Farahi, A. Binder et al.,
// Lubricants 12(5), 148, 2024), optional (0 = purely resistive). So
// V·I = H + ω·T exactly, H = VF·|I| + R·I² the contact heat: the shaft work comes from the electrical
// circuit, as a back-EMF (P. G. Moyssides, P. Hatzikonstantinou, IEEE
// Trans. Magn. 26, 1274, 1990, measured a speed-dependent back-EMF and no
// generator action). In THERMAL=1 mode (S. Marinov, Electronics & Wireless
// World, April 1989; McDonald 2020: thermal start-up) there is no back-EMF,
// and the shaft work is taken from the contact heat H with the measured
// efficiency ceiling η:
//
//	T = a/(1 + a·ω/(η·H)),  a = κ·I²·tanh(ω/ω0),   ω·T < η·H
//
// The coefficients carry no defaults: KAPPA [N·m/A²] and W0 [rad/s] (and
// ETAMAX in thermal mode) must come from a fit to a cited measurement
// (examples/bearing_motor). Value is the bearing's resistance R > 0 [Ω].
const BearingMotor Kind = 'P'

// BearingParams returns κ, ω0, thermal mode and η of a bearing element.
func BearingParams(e Element) (kappa, w0 float64, thermal bool, eta float64) {
	return e.Params["kappa"], e.Params["w0"], e.Params["thermal"] == 1, e.Params["etamax"]
}

// BearingMotor adds a ball-bearing motor between electrical nodes a, b and
// the mechanical port mp, mn. kappa, w0 (and eta in thermal mode) come from
// a measurement fit; see Kind BearingMotor. BearingMotorVF also sets the
// measured contact drop.
func (b *Builder) BearingMotor(name, a, c, mp, mn string, r, kappa, w0 float64, thermal bool, eta float64) *Builder {
	return b.BearingMotorVF(name, a, c, mp, mn, r, 0, kappa, w0, thermal, eta)
}

// BearingMotorVF is BearingMotor with a measured contact drop vf [V].
func (b *Builder) BearingMotorVF(name, a, c, mp, mn string, r, vf, kappa, w0 float64, thermal bool, eta float64) *Builder {
	p := map[string]float64{"kappa": kappa, "w0": w0}
	if vf != 0 {
		p["vf"] = vf
	}
	if thermal {
		p["thermal"], p["etamax"] = 1, eta
	}
	return b.AddMulti(Element{Name: name, Kind: BearingMotor, Value: r, Params: p}, a, c, mp, mn)
}

func validateBearing(e Element) error {
	if len(e.Extra) != 2 || e.Extra[0] == e.Extra[1] {
		return fmt.Errorf("%w: %s needs two distinct mechanical nodes", ErrBadValue, e.Name)
	}
	if !(e.Value > 0) {
		return fmt.Errorf("%w: %s bearing resistance must be > 0, got %g", ErrBadValue, e.Name, e.Value)
	}
	for k := range e.Params {
		switch k {
		case "kappa", "w0", "thermal", "etamax", "vf":
		default:
			return fmt.Errorf("%w: %s: unknown parameter %q", ErrBadValue, e.Name, k)
		}
	}
	kappa, w0, thermal, eta := BearingParams(e)
	if vf := e.Params["vf"]; !(vf >= 0) || math.IsInf(vf, 0) {
		return fmt.Errorf("%w: %s VF must be a finite measured drop ≥ 0", ErrBadValue, e.Name)
	}
	for _, v := range []float64{kappa, w0, eta} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return fmt.Errorf("%w: %s has a non-finite parameter", ErrBadValue, e.Name)
		}
	}
	// No defaults: an unset κ or ω0 is an error, not an assumption.
	_, hasK := e.Params["kappa"]
	if !hasK || !(kappa >= 0) || !(w0 > 0) {
		return fmt.Errorf("%w: %s needs measured KAPPA ≥ 0 and W0 > 0", ErrBadValue, e.Name)
	}
	if t := e.Params["thermal"]; t != 0 && t != 1 {
		return fmt.Errorf("%w: %s THERMAL must be 0 or 1", ErrBadValue, e.Name)
	}
	if thermal && !(eta > 0 && eta < 1) {
		return fmt.Errorf("%w: %s thermal mode needs a measured ETAMAX in (0, 1)", ErrBadValue, e.Name)
	}
	if !thermal && eta != 0 {
		return fmt.Errorf("%w: %s ETAMAX applies to thermal mode only", ErrBadValue, e.Name)
	}
	return nil
}
