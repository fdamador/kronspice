// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math/cmplx"

	"github.com/fdamador/kronspice/device"
)

// Hypothetical ECE / ECE2 spin-connection terms in Kron's frame (register
// entry H48; circuit/ece.go gives the field relations and the mapping).
//
// In the primitive network each R, C or L carrying the terms keeps its own
// branch, and its relation sees the body voltage
//
//	v_b = E + s·κ_a·φ_a − s·κ_b·φ_b,
//
// where φ_a, φ_b are the terminals' potentials: junction-pair voltages to
// the datum, carried by control coils (open, i = 0) from each terminal to
// the datum, as a VCVS carries its controlling pair (Kron 1939 Ch. XVI:
// open paths, "active" junction pairs). The relation is non-reciprocal,
// like the gyrator's. An inductor's relation adds the torsion resistance
// r_t = s·cω₀L (v_b = (jωL + r_t)·i); a K element's chord becomes
// z = jωM + s·cω₀M, so mutual torsion is a chord term, as Durand et al.'s
// TAN chords are (CLAUDE.md §3a.1). UFT 456's terms (kronspice-0f's
// device-side definitions) enter the same way, as Kron's 1944 analogy reads
// them: k0 a shunt conductance g_e = 2cκ₀·C across a capacitor's body
// (i = (jωC + g_e)·v_b), kl0 a series resistance −2cκ_Λ0·L added to r_t
// and −2cκ_Λ0·M to the K chord.
//
// The ledger books the terms' power as the device models do (device/ece.go):
// per element as model-supplied, labelled hypothetical, never as inflow.
//
// Kron's own precedents (Tensors for Circuits, Dover 1959):
//   - ch. 32, eq. 32.15: the field tensor in an accelerated frame keeps
//     only the antisymmetric part (torsion) of the connection, which is the
//     pattern of ECE's "field = torsion";
//   - eq. 32.12: e = L·δi/dt = L·di/dt + Γ·i·ẋ. With the time coordinate
//     ẋ⁰ = c this is the torsion term s·cω₀L·i, power booked because
//     nothing pays it;
//   - eqs. 31.9–31.11: the standard counterpart is the commutator machine,
//     G independent of L, whose power goes through the shaft
//     (analysis.TestECECommutatorTwin);
//   - eq. 31.5: a skew F does no work. K's w0skew places mutual torsion
//     that way (analysis.TestECETorsionPlacement);
//   - eq. 31.14: Γ is not a tensor. A frame change adds Cᵀ·L·dC/dt, whose
//     symmetric part is ½·d(CᵀLC)/dt. So with a constant transformed
//     inductance a frame generates only skew terms, never the symmetric r_t
//     (analysis.TestFrameTermSymmetricPart);
//   - p. 245: an equivalent circuit needs a tensor equation. Of the
//     scalar terms, only a κ that differs between an element's leads
//     depends on the datum; a uniform κ scales every impedance by
//     1/(1 + s·κ) (analysis.TestECEDatumShift, TestECEUniformKappaIsScaling).

// hasECE reports whether branch k carries spin-connection terms.
func (b branch) hasECE() bool { return b.ca >= 0 || b.cb >= 0 || b.rt != 0 || b.ge != 0 }

// bodyVoltage is v_b of branch k in a solution.
func (n *Network) bodyVoltage(sol *Solution, k int) complex128 {
	b := n.branches[k]
	v := sol.E.At(k)
	if b.ca >= 0 {
		v += complex(b.sa, 0) * sol.E.At(b.ca)
	}
	if b.cb >= 0 {
		v -= complex(b.sb, 0) * sol.E.At(b.cb)
	}
	return v
}

// eceOps completes the operating points of the elements that carry
// spin-connection terms, as device.Operating / OperatingAC would: body
// voltage, storage rate and the power the terms deliver (ECEPower). At DC
// the values are instantaneous; in AC, ECEPower is the average ½·Re(·).
func (n *Network) eceOps(ops []device.OperatingPoint, sol *Solution, dc bool) {
	for ei, k := range n.elemBranch {
		if k < 0 || !n.branches[k].hasECE() {
			continue
		}
		b := n.branches[k]
		E, I := sol.E.At(k), sol.I.At(k)
		vb := n.bodyVoltage(sol, k)
		vl := vb - complex(b.rt, 0)*I // the voltage that stores energy (vb for R, C)
		op := &ops[ei]
		// UFT 456's shunt g_e on a capacitor carries g_e·v_b of the terminal
		// current; its power g_e·v_b² is absorbed by the term (booked as
		// negative model-supplied, as device/reactive.go does).
		ge := complex(b.ge, 0)
		if dc {
			op.BodyVoltage = real(vb)
			op.ECEPower = real(vl-E)*real(I) - b.ge*real(vb)*real(vb)
			if b.kind != bR {
				op.StoragePower = real(vl) * real(I-ge*vb)
			}
			continue
		}
		op.ECEPower = 0.5*real((vl-E)*cmplx.Conj(I)) - 0.5*b.ge*real(vb*cmplx.Conj(vb))
	}
	for _, m := range n.eceK {
		i1, i2 := sol.I.At(m.a), sol.I.At(m.b)
		if dc {
			// The coils' storage rates include r_m·i_other·i; the K element
			// takes it back out and books it as model-supplied.
			p := -2 * m.rm * real(i1) * real(i2)
			ops[m.elem].StoragePower, ops[m.elem].ECEPower = p, p
			continue
		}
		ops[m.elem].ECEPower = -m.rm * real(i1*cmplx.Conj(i2))
	}
}
