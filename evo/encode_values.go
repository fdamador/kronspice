// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"fmt"
	"math"
	"strings"
)

// Preferred-number series, IEC 60063:2015 (mantissas of one decade).
var (
	e6  = []float64{1.0, 1.5, 2.2, 3.3, 4.7, 6.8}
	e12 = []float64{1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2}
	e24 = []float64{1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0,
		3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1}
	e96 = []float64{1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30,
		1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74,
		1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32,
		2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09,
		3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12,
		4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49,
		5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32,
		7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76}
	e48 = everyOther(e96)
)

func everyOther(s []float64) []float64 {
	out := make([]float64, 0, len(s)/2)
	for i := 0; i < len(s); i += 2 {
		out = append(out, s[i])
	}
	return out
}

// series returns the mantissas of a named series (nil for none).
func series(name string) ([]float64, error) {
	switch strings.ToUpper(strings.TrimSpace(name)) {
	case "", "NONE":
		return nil, nil
	case "E6":
		return e6, nil
	case "E12":
		return e12, nil
	case "E24":
		return e24, nil
	case "E48":
		return e48, nil
	case "E96":
		return e96, nil
	}
	return nil, fmt.Errorf("unknown series %q (want E6, E12, E24, E48 or E96)", name)
}

// snap rounds a positive magnitude to the nearest series value in log
// distance, preferring values inside [min, max]. With no series (or when no
// series value lies in the range) the clamped value is returned unchanged.
func snap(v float64, ser []float64, min, max float64) float64 {
	v = math.Min(math.Max(v, min), max)
	if ser == nil {
		return v
	}
	lv := math.Log10(v)
	dec := math.Floor(lv)
	best, bestD := v, math.Inf(1)
	for d := dec - 1; d <= dec+1; d++ {
		p := math.Pow(10, d)
		for _, m := range ser {
			// Round to the series' own precision so 4.7e3 is exactly 4700.
			c := roundSig(m*p, 3)
			if c < min*(1-1e-12) || c > max*(1+1e-12) {
				continue
			}
			if dd := math.Abs(math.Log10(c) - lv); dd < bestD {
				best, bestD = c, dd
			}
		}
	}
	return best
}

// roundSig rounds x to n significant digits.
func roundSig(x float64, n int) float64 {
	if x == 0 {
		return 0
	}
	e := math.Floor(math.Log10(math.Abs(x))) - float64(n-1)
	p := math.Pow(10, e)
	r := math.Round(x/p) * p
	// Remove representation noise: format and reparse.
	var out float64
	fmt.Sscanf(fmt.Sprintf("%.*e", n-1, r), "%g", &out)
	return out
}
