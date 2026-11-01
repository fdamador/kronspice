// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package device

import "math"

// TruncContext carries the tolerances and step data for local truncation
// error (LTE) estimates.
type TruncContext struct {
	Time   float64 // time of the candidate solution [s]
	H      float64 // candidate step [s]
	Method Method  // integration order used for the candidate
	RelTol float64 // relative tolerance (SPICE reltol, default 1e-3)
	AbsTol float64 // absolute current tolerance (SPICE abstol, default 1e-12 A)
	VnTol  float64 // absolute voltage tolerance (SPICE vntol, default 1e-6 V)
	ChgTol float64 // charge/flux tolerance (SPICE chgtol, default 1e-14)
	TrTol  float64 // truncation-error overestimate factor (SPICE trtol, default 7)
}

// history keeps the last accepted values of a state variable (capacitor
// charge or inductor flux) for divided-difference LTE estimates.
type history struct {
	t [3]float64
	q [3]float64
	n int // number of valid points (≤ 3), most recent at index 0
}

func (h *history) push(t, q float64) {
	h.t[2], h.t[1], h.t[0] = h.t[1], h.t[0], t
	h.q[2], h.q[1], h.q[0] = h.q[1], h.q[0], q
	if h.n < 3 {
		h.n++
	}
}

func (h *history) reset() { h.n = 0 }

// lteStep proposes a new step size for one state variable.
//
// The (k+1)-th derivative of the state q is estimated from a divided
// difference over the candidate point and the last k+1 accepted points
// (k = integration order). The local truncation error in q is
// C_k·h^(k+1)·|q^(k+1)| with C_1 = 1/2 (backward Euler) and C_2 = 1/12
// (trapezoidal); dividing by h gives the error in the element current. The
// allowed current error is trtol·max(reltol·max|i| + abs,
// reltol·max(|q|, chgtol)/h), the same tolerance structure as SPICE's
// truncation-error control (Nagel 1975, ch. 6; SPICE3 CKTterr). The proposal
// is h·(allowed/error)^(1/k), since the current error scales as h^k.
func lteStep(hist *history, qn, in, inPrev, absTol float64, tc *TruncContext) float64 {
	k := int(tc.Method)
	if k < 1 || k > 2 || hist.n < k+1 {
		return inf // not enough history yet; the engine limits growth instead
	}
	// Points: (tc.Time, qn), then history[0..k].
	ts := [4]float64{tc.Time, hist.t[0], hist.t[1], hist.t[2]}
	qs := [4]float64{qn, hist.q[0], hist.q[1], hist.q[2]}
	// Divided differences up to order k+1.
	var dd [4]float64
	copy(dd[:], qs[:k+2])
	for order := 1; order <= k+1; order++ {
		for i := 0; i+order <= k+1; i++ {
			den := ts[i] - ts[i+order]
			if den == 0 {
				return inf
			}
			dd[i] = (dd[i] - dd[i+1]) / den
		}
	}
	deriv := dd[0] // (k+1)-th divided difference ≈ q^(k+1)/(k+1)!
	fact := 2.0
	c := 0.5
	if k == 2 {
		fact, c = 6, 1.0/12
	}
	h := tc.H
	errI := c * math.Pow(h, float64(k)) * math.Abs(deriv*fact)
	if errI == 0 || math.IsNaN(errI) {
		return inf
	}
	qmax := math.Max(math.Abs(qn), math.Abs(hist.q[0]))
	allowed := tc.TrTol * math.Max(tc.RelTol*math.Max(math.Abs(in), math.Abs(inPrev))+absTol,
		tc.RelTol*math.Max(qmax, tc.ChgTol)/h)
	return h * math.Pow(allowed/errI, 1/float64(k))
}
