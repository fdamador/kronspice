// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"fmt"
	"math"
	"sort"
)

// The momentum ledger (CLAUDE.md §3f "Momentum ledger and push–pull
// designs"; doc/em-push-pull.md). Beside the energy ledger, every force
// the engine computes is booked with the body it acts on and the source
// it comes from: another body of the system, or an exterior channel (the
// source of a declared external field, the plasma, radiation). Internal
// forces must cancel in pairs (Newton's third law, which holds exactly for
// the magnetostatic forces between closed circuits; Jackson §5.7); the
// net force on the system then equals the sum over exterior sources, and
// each such row names where the reaction momentum goes.

// ForceRow is one booked force: On is the body it acts on, From the body
// or exterior source that exerts it (and receives the reaction).
type ForceRow struct {
	On, From     string
	F            [3]float64 // [N]
	Exterior     bool       // From is an exterior source
	Hypothetical bool       // from a hypothetical model (kept apart)
}

// MomentumLedger collects ForceRows.
type MomentumLedger struct {
	Rows []ForceRow
}

// AddInternal books the force of body `from` on body `on`.
func (l *MomentumLedger) AddInternal(on, from string, f [3]float64) {
	l.Rows = append(l.Rows, ForceRow{On: on, From: from, F: f})
}

// AddExterior books a force on a body from an exterior source; the source
// receives −f.
func (l *MomentumLedger) AddExterior(on, source string, f [3]float64, hypothetical bool) {
	l.Rows = append(l.Rows, ForceRow{On: on, From: source, F: f, Exterior: true, Hypothetical: hypothetical})
}

// Net is the net force on the system from real (non-hypothetical) rows:
// internal rows should cancel, so it equals the exterior sum.
func (l *MomentumLedger) Net() [3]float64 {
	var s [3]float64
	for _, r := range l.Rows {
		if !r.Hypothetical {
			s = add3(s, r.F)
		}
	}
	return s
}

// ByExterior is the force from each exterior source (real rows only),
// with the reaction each source receives being its negative.
func (l *MomentumLedger) ByExterior() map[string][3]float64 {
	out := map[string][3]float64{}
	for _, r := range l.Rows {
		if r.Exterior && !r.Hypothetical {
			out[r.From] = add3(out[r.From], r.F)
		}
	}
	return out
}

// Hypothetical is the net force from hypothetical rows, never part of Net.
func (l *MomentumLedger) Hypothetical() [3]float64 {
	var s [3]float64
	for _, r := range l.Rows {
		if r.Hypothetical {
			s = add3(s, r.F)
		}
	}
	return s
}

// Check verifies Newton's third law for every pair of internal bodies:
// |F(a←b) + F(b←a)| ≤ tol·max(|F|). It returns an error naming the worst
// pair otherwise — an internal force that does not cancel would be thrust
// from nowhere, a numerical or modelling error.
func (l *MomentumLedger) Check(tol float64) error {
	pair := map[[2]string][3]float64{}
	scale := 0.0
	for _, r := range l.Rows {
		if r.Exterior || r.Hypothetical {
			continue
		}
		a, b := r.On, r.From
		if a > b {
			a, b = b, a
		}
		pair[[2]string{a, b}] = add3(pair[[2]string{a, b}], r.F)
		scale = math.Max(scale, math.Sqrt(dot(r.F, r.F)))
	}
	keys := make([][2]string, 0, len(pair))
	for k := range pair {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		return keys[i][0] < keys[j][0] || keys[i][0] == keys[j][0] && keys[i][1] < keys[j][1]
	})
	for _, k := range keys {
		f := pair[k]
		if m := math.Sqrt(dot(f, f)); m > tol*scale {
			return fmt.Errorf("em: momentum ledger: forces between %s and %s do not cancel (%.3g N of %.3g N)", k[0], k[1], m, scale)
		}
	}
	return nil
}
