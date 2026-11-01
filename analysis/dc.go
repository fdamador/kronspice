// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"strings"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/output"
)

// Errors from DC sweeps.
var (
	ErrBadDCSpec     = errors.New("analysis: invalid DC sweep")
	ErrUnknownSource = errors.New("analysis: swept source not found")
)

// DCSweep is one source swept by .dc: start to stop (inclusive) by step.
type DCSweep struct {
	Source            string
	Start, Stop, Step float64
}

// MaxDCPoints bounds the total number of operating points in one sweep.
const MaxDCPoints = 1_000_000

func (s DCSweep) values() ([]float64, error) {
	if s.Step == 0 || math.IsNaN(s.Start+s.Stop+s.Step) || math.IsInf(s.Start+s.Stop+s.Step, 0) {
		return nil, fmt.Errorf("%w: %+v", ErrBadDCSpec, s)
	}
	n := (s.Stop - s.Start) / s.Step
	if n < 0 || n > MaxDCPoints {
		return nil, fmt.Errorf("%w: %+v", ErrBadDCSpec, s)
	}
	count := int(math.Floor(n+1e-9)) + 1
	out := make([]float64, count)
	for k := range out {
		out[k] = s.Start + float64(k)*s.Step
	}
	return out, nil
}

// DC sweeps one or two sources (sweeps[0] is the inner loop, as in SPICE) and
// solves the operating point at every value. The ledger is checked at every
// point.
func DC(c *circuit.Circuit, sweeps []DCSweep, opt Options) (*output.DCSweepResult, error) {
	if len(sweeps) < 1 || len(sweeps) > 2 {
		return nil, fmt.Errorf("%w: need 1 or 2 sweeps", ErrBadDCSpec)
	}
	s, err := New(c)
	if err != nil {
		return nil, err
	}
	vals := make([][]float64, len(sweeps))
	srcs := make([]device.DCOverrider, len(sweeps))
	total := 1
	for i, sw := range sweeps {
		if vals[i], err = sw.values(); err != nil {
			return nil, err
		}
		total *= len(vals[i])
		for _, d := range s.devs {
			if strings.EqualFold(d.Name(), sw.Source) {
				if o, ok := d.(device.DCOverrider); ok {
					srcs[i] = o
				}
			}
		}
		if srcs[i] == nil {
			return nil, fmt.Errorf("%w: %s", ErrUnknownSource, sw.Source)
		}
	}
	if total > MaxDCPoints {
		return nil, fmt.Errorf("%w: %d points", ErrBadDCSpec, total)
	}
	defer func() {
		for _, o := range srcs {
			o.ClearDCOverride()
		}
	}()

	r := &output.DCSweepResult{
		Title:     c.Title,
		NodeNames: append([]string(nil), c.NodeNames...),
		V:         make([][]float64, c.NumNodes()),
		I:         make([][]float64, len(s.devs)),
		Values:    make([][]float64, len(sweeps)),
	}
	for _, sw := range sweeps {
		r.Sources = append(r.Sources, sw.Source)
	}
	for _, d := range s.devs {
		r.ElementNames = append(r.ElementNames, d.Name())
	}
	outer := []float64{0}
	if len(sweeps) == 2 {
		outer = vals[1]
	}
	for _, vo := range outer {
		if len(sweeps) == 2 {
			srcs[1].SetDCOverride(vo)
		}
		for _, vi := range vals[0] {
			srcs[0].SetDCOverride(vi)
			res, err := s.OP(opt)
			if err != nil {
				return r, fmt.Errorf("%s = %g: %w", sweeps[0].Source, vi, err)
			}
			r.Values[0] = append(r.Values[0], vi)
			if len(sweeps) == 2 {
				r.Values[1] = append(r.Values[1], vo)
			}
			for n := 1; n < c.NumNodes(); n++ {
				r.V[n] = append(r.V[n], res.NodeVoltages[n])
			}
			for i, e := range res.Elements {
				r.I[i] = append(r.I[i], e.Current)
			}
			if rr := res.Ledger.RelativeResidual(); rr > r.MaxRelativeResidual {
				r.MaxRelativeResidual = rr
			}
		}
	}
	return r, nil
}
