// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
	"sort"
	"strings"

	"github.com/fdamador/kronspice/circuit"
)

// The lumped (quasi-static) limit (CLAUDE.md §3c, "What em/ does" item 3).
//
// Kirchhoff's laws follow from Maxwell's equations when the changing
// fields stay inside the components: charge is stored only in capacitors
// and flux only in inductors, so the scalar part of ∇ applied to the total
// current vanishes at a node and the vector part of ∇𝔈 integrates to the
// inductor voltages around a loop. That holds when every part of the
// circuit is small against the shortest wavelength present: a conductor of
// length d carries a signal of frequency f with a phase delay
// 2π·d/λ, λ = c/(f·√εr), which a lumped model sets to zero. The usual
// engineering limit is d < λ/10 (D. M. Pozar, Microwave Engineering, 4th
// ed., Wiley 2012, §2.1: lumped versus distributed; H. W. Johnson,
// M. Graham, High-Speed Digital Design, Prentice Hall 1993, §1.1 and §4.1
// for the edge-rate form). Beyond it the element should be a transmission
// line (em/tline, M9i) or a field model.
//
// For a transient the shortest wavelength is set by the fastest edge: a
// rise time t_r has most of its spectrum below the knee frequency
// f_knee = 0.5/t_r (Johnson & Graham §1.1).

// ErrNotQuasiStatic means part of a circuit is electrically large at the
// highest frequency of the analysis, so the lumped model is not valid.
var ErrNotQuasiStatic = errors.New("em: circuit is not quasi-static (electrically large)")

// SpeedOfLight is c in vacuum [m/s] (exact, SI 2019).
const SpeedOfLight = 299792458.0

// SizeOptions configures CheckQuasiStatic.
type SizeOptions struct {
	// Threshold is the largest allowed d/λ (default 0.1).
	Threshold float64
	// EpsR is the effective relative permittivity of the medium around the
	// conductors (default 1, vacuum or air; ≈ 3.3 for microstrip on FR-4).
	EpsR float64
	// WarnOnly reports the electrically large parts without an error.
	WarnOnly bool
}

// SizeEntry is one checked part.
type SizeEntry struct {
	Name    string  // element name, or "circuit" for the overall size
	Length  float64 // [m]
	Ratio   float64 // d/λ at the highest frequency
	TooLong bool
}

// SizeReport is the result of CheckQuasiStatic.
type SizeReport struct {
	FMax       float64 // highest frequency checked [Hz]
	Wavelength float64 // λ at FMax in the medium [m]
	Entries    []SizeEntry
}

// CheckQuasiStatic compares the declared physical lengths (element name →
// conductor length in metres, and the overall circuit size, 0 if not
// declared) with the wavelength at fmax. Entries are sorted by d/λ,
// largest first. With SizeOptions.WarnOnly it never fails; otherwise it
// returns ErrNotQuasiStatic naming the parts above the threshold.
func CheckQuasiStatic(c *circuit.Circuit, lengths map[string]float64, overall, fmax float64, opt SizeOptions) (*SizeReport, error) {
	if opt.Threshold <= 0 {
		opt.Threshold = 0.1
	}
	if opt.EpsR <= 0 {
		opt.EpsR = 1
	}
	if !(fmax > 0) || math.IsInf(fmax, 0) {
		return nil, fmt.Errorf("em: CheckQuasiStatic: highest frequency %g Hz", fmax)
	}
	lam := SpeedOfLight / (fmax * math.Sqrt(opt.EpsR))
	rep := &SizeReport{FMax: fmax, Wavelength: lam}
	add := func(name string, d float64) error {
		if !(d >= 0) || math.IsInf(d, 0) {
			return fmt.Errorf("%w: %s: length %g m", ErrOutsideModelValidity, name, d)
		}
		r := d / lam
		rep.Entries = append(rep.Entries, SizeEntry{Name: name, Length: d, Ratio: r, TooLong: r > opt.Threshold})
		return nil
	}
	names := make([]string, 0, len(lengths))
	for n := range lengths {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, n := range names {
		if c != nil && c.Index(n) < 0 {
			return nil, fmt.Errorf("em: CheckQuasiStatic: no element %q", n)
		}
		if err := add(n, lengths[n]); err != nil {
			return nil, err
		}
	}
	if overall > 0 {
		if err := add("circuit", overall); err != nil {
			return nil, err
		}
	}
	sort.SliceStable(rep.Entries, func(i, j int) bool { return rep.Entries[i].Ratio > rep.Entries[j].Ratio })
	var long []string
	for _, e := range rep.Entries {
		if e.TooLong {
			long = append(long, fmt.Sprintf("%s (%.3g m = %.3g λ)", e.Name, e.Length, e.Ratio))
		}
	}
	if len(long) > 0 && !opt.WarnOnly {
		return rep, fmt.Errorf("%w at %.4g Hz (λ = %.4g m, limit %.3g λ): %s — model it as a transmission line or a field region",
			ErrNotQuasiStatic, fmax, lam, opt.Threshold, strings.Join(long, ", "))
	}
	return rep, nil
}

// SourceBandwidth is the highest significant frequency of the circuit's
// independent sources in a transient: SIN at its frequency; PULSE and PWL
// edges at the knee 0.5/t_r of their shortest rise, fall or segment; EXP at
// the knee of its 10–90 % rise 2.2·τ. defaultEdge stands in for a PULSE rise
// or fall left at its SPICE default (the time step); 0 if none.
func SourceBandwidth(c *circuit.Circuit, defaultEdge float64) float64 {
	f := 0.0
	knee := func(tr float64) {
		if tr > 0 {
			f = math.Max(f, 0.5/tr)
		}
	}
	for _, e := range c.Elements {
		if e.Kind != circuit.VoltageSource && e.Kind != circuit.CurrentSource {
			continue
		}
		w := e.Wave
		arg := func(i int) (float64, bool) {
			if i < len(w.Args) {
				return w.Args[i], true
			}
			return 0, false
		}
		switch w.Kind {
		case circuit.WaveSin:
			if fr, ok := arg(2); ok && fr > 0 {
				f = math.Max(f, fr)
			}
		case circuit.WavePulse:
			for _, k := range []int{3, 4} {
				if tr, ok := arg(k); ok && tr > 0 {
					knee(tr)
				} else {
					knee(defaultEdge)
				}
			}
		case circuit.WaveExp:
			for _, k := range []int{3, 5} {
				if tau, ok := arg(k); ok && tau > 0 {
					knee(2.2 * tau)
				}
			}
		case circuit.WavePWL:
			for i := 2; i+1 < len(w.Args); i += 2 {
				if w.Args[i+1] != w.Args[i-1] {
					knee(w.Args[i] - w.Args[i-2])
				}
			}
		}
	}
	return f
}
