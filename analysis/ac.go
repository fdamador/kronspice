// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/output"
)

// ErrBadACSpec is returned for an invalid AC sweep.
var ErrBadACSpec = errors.New("analysis: invalid AC specification")

// ACSpec is an .ac request.
type ACSpec struct {
	Variation     string // "dec", "oct" or "lin"
	Points        int    // per decade / per octave, or total for lin
	FStart, FStop float64
}

// MaxACPoints bounds the number of frequencies in one sweep.
const MaxACPoints = 1_000_000

// Frequencies returns the sweep points, as SPICE3 generates them.
func (s ACSpec) Frequencies() ([]float64, error) {
	if s.Points < 1 || s.FStop < s.FStart || s.FStart < 0 || math.IsNaN(s.FStart+s.FStop) || math.IsInf(s.FStop, 0) {
		return nil, fmt.Errorf("%w: %+v", ErrBadACSpec, s)
	}
	var out []float64
	switch s.Variation {
	case "lin":
		if s.Points == 1 {
			return []float64{s.FStart}, nil
		}
		step := (s.FStop - s.FStart) / float64(s.Points-1)
		for k := 0; k < s.Points; k++ {
			out = append(out, s.FStart+float64(k)*step)
		}
	case "dec", "oct":
		if s.FStart <= 0 {
			return nil, fmt.Errorf("%w: logarithmic sweep needs fstart > 0", ErrBadACSpec)
		}
		base := 10.0
		if s.Variation == "oct" {
			base = 2
		}
		ratio := math.Pow(base, 1/float64(s.Points))
		n := int(math.Floor(math.Log(s.FStop/s.FStart)/math.Log(ratio)+1e-9)) + 1
		if n > MaxACPoints {
			return nil, fmt.Errorf("%w: %d points", ErrBadACSpec, n)
		}
		for k := 0; k < n; k++ {
			out = append(out, s.FStart*math.Pow(ratio, float64(k)))
		}
	default:
		return nil, fmt.Errorf("%w: sweep type %q", ErrBadACSpec, s.Variation)
	}
	if len(out) > MaxACPoints {
		return nil, fmt.Errorf("%w: %d points", ErrBadACSpec, len(out))
	}
	return out, nil
}

// AC runs a small-signal AC analysis: the DC operating point first (devices
// are linearized there; all M3 devices are linear), then one complex MNA
// solve per frequency (CLAUDE.md §3, AC). Each frequency gets an
// average-power ledger: power delivered by AC sources equals the power
// dissipated, and reactive power sums to zero (complex Tellegen theorem).
func AC(c *circuit.Circuit, spec ACSpec, opt Options) (*output.ACResult, error) {
	freqs, err := spec.Frequencies()
	if err != nil {
		return nil, err
	}
	s, err := New(c)
	if err != nil {
		return nil, err
	}
	if _, err := s.OP(opt); err != nil {
		return nil, fmt.Errorf("AC operating point: %w", err)
	}
	if err := s.polish(); err != nil {
		return nil, fmt.Errorf("AC operating point: %w", err)
	}
	n := s.m.Size()
	m, err := opt.cmatrix(n)
	if err != nil {
		return nil, err
	}
	rhs, x := make([]complex128, n), make([]complex128, n)
	absBuf := make([]float64, n)
	actx := &device.ACContext{M: m, RHS: rhs}
	ops := make([]device.OperatingPoint, len(s.devs))

	r := &output.ACResult{
		Title:     c.Title,
		Freq:      freqs,
		NodeNames: append([]string(nil), c.NodeNames...),
		V:         make([][]complex128, c.NumNodes()),
		I:         make([][]complex128, len(s.devs)),
	}
	for _, d := range s.devs {
		r.ElementNames = append(r.ElementNames, d.Name())
		r.ElementKinds = append(r.ElementKinds, d.Kind().String())
	}
	for i := 1; i < c.NumNodes(); i++ {
		r.V[i] = make([]complex128, len(freqs))
	}
	for i := range s.devs {
		r.I[i] = make([]complex128, len(freqs))
	}
	l := &r.Ledger
	l.Inflow = make([]float64, len(freqs))
	l.Outflow = make([]float64, len(freqs))
	l.Dissipation = make([]float64, len(freqs))
	l.Reactive = make([]float64, len(freqs))
	l.Residual = make([]float64, len(freqs))
	l.ModelSupplied = make([]float64, len(freqs))

	for k, f := range freqs {
		if err := opt.interrupted(); err != nil {
			return r, err
		}
		m.Clear()
		for i := range rhs {
			rhs[i] = 0
		}
		actx.Omega = 2 * math.Pi * f
		for _, d := range s.devs {
			if err := d.LoadAC(actx); err != nil {
				return r, fmt.Errorf("load AC %s: %w", d.Name(), err)
			}
		}
		if err := m.Solve(rhs, x); err != nil {
			return r, fmt.Errorf("f = %g Hz: %w", f, err)
		}
		for i := 1; i < c.NumNodes(); i++ {
			r.V[i][k] = x[i-1]
		}
		var acc em.Accumulator
		var sumS complex128
		scaleS := 0.0
		for i, d := range s.devs {
			p := d.OperatingAC(x)
			r.I[i][k] = p.Current
			sAbs := 0.5 * p.Voltage * cmplx.Conj(p.Current) // complex power absorbed
			if mp, ok := d.(device.ACPowerer); ok {
				sAbs = mp.ACPower(x)
			}
			sumS += sAbs
			scaleS += cmplx.Abs(sAbs)
			// Average power drives the exterior ports; the element models'
			// ExteriorPower functions are linear in absorbed power.
			ops[i] = device.OperatingPoint{Voltage: cmplx.Abs(p.Voltage), Current: cmplx.Abs(p.Current), PowerAbsorbed: real(sAbs)}
			if _, ok := d.(device.Reactive); ok {
				l.Reactive[k] += imag(sAbs)
			}
		}
		for _, pr := range s.ports {
			acc.AddPort(pr.ch, pr.status, s.devs[pr.dev].ExteriorPower(pr.ch, ops[pr.dev]))
		}
		for _, i := range s.supply {
			acc.AddModelSupplied(s.devs[i].(device.ModelSupplier).ModelSuppliedPower(ops[i]))
		}
		l.ModelSupplied[k] = acc.ModelSupplied
		l.Inflow[k], l.Outflow[k], l.Dissipation[k] = acc.Inflow, acc.Outflow, acc.Dissipation
		l.Residual[k] = acc.Residual()
		// Round-off floor, as for the DC ledger (output.Ledger.Roundoff).
		m.AbsProduct(x, absBuf)
		ro := 0.0
		for i, v := range absBuf {
			ro += cmplx.Abs(x[i]) * (v + cmplx.Abs(rhs[i]))
		}
		ro *= 64 * 0x1p-52
		rr := 0.0
		if scaleS > 0 && cmplx.Abs(sumS) > ro {
			rr = cmplx.Abs(sumS) / scaleS
		}
		if rr > l.MaxRelativeResidual {
			l.MaxRelativeResidual = rr
		}
		if !opt.SkipEnergyCheck && rr > opt.ledgerTol() {
			return r, fmt.Errorf("%w at f = %g Hz: relative %.3g", em.ErrLedgerResidual, f, rr)
		}
	}
	return r, nil
}
