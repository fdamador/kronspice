// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package netlist

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// ErrBadValue is returned when a numeric field cannot be parsed.
var ErrBadValue = errors.New("netlist: bad numeric value")

// ParseValue parses a SPICE number with an optional engineering suffix.
//
// Suffixes are case-insensitive, as in SPICE3: T=1e12, G=1e9, MEG=1e6,
// K=1e3, M=1e-3, MIL=25.4e-6, U=1e-6, N=1e-9, P=1e-12, F=1e-15. Any letters
// after the suffix (for example the "ohm" in "10kohm") are ignored. This is
// the only place in the code base where engineering suffixes are handled
// (CLAUDE.md §8 rule 3). Results must be finite.
func ParseValue(s string) (float64, error) {
	t := strings.TrimSpace(s)
	if t == "" {
		return 0, fmt.Errorf("%w: empty", ErrBadValue)
	}
	// Longest numeric prefix (digits, sign, point, exponent).
	end := numericPrefix(t)
	if end == 0 {
		return 0, fmt.Errorf("%w: %q", ErrBadValue, s)
	}
	lit := t[:end]
	rest := strings.ToLower(t[end:])
	// Power-of-ten suffixes are applied as a decimal exponent before the
	// conversion, so "10u" is the correctly rounded double nearest 1e-5
	// (multiplying by 1e-6 afterwards would round twice).
	exp10, scale := 0, 1.0
	switch {
	case strings.HasPrefix(rest, "meg"):
		exp10 = 6
	case strings.HasPrefix(rest, "mil"):
		scale = 25.4e-6
	case rest == "":
	default:
		switch rest[0] {
		case 't':
			exp10 = 12
		case 'g':
			exp10 = 9
		case 'k':
			exp10 = 3
		case 'm':
			exp10 = -3
		case 'u':
			exp10 = -6
		case 'n':
			exp10 = -9
		case 'p':
			exp10 = -12
		case 'f':
			exp10 = -15
		default:
			if rest[0] < 'a' || rest[0] > 'z' {
				return 0, fmt.Errorf("%w: %q", ErrBadValue, s)
			}
			// Other letters are a unit name only (e.g. "5V"): scale 1.
		}
	}
	if exp10 != 0 {
		mant, e := lit, 0
		if k := strings.IndexAny(lit, "eE"); k >= 0 {
			n, err := strconv.Atoi(lit[k+1:])
			if err != nil || n > 100000 || n < -100000 {
				return 0, fmt.Errorf("%w: %q", ErrBadValue, s)
			}
			mant, e = lit[:k], n
		}
		lit = mant + "e" + strconv.Itoa(e+exp10)
	}
	v, err := strconv.ParseFloat(lit, 64)
	if err != nil && !errors.Is(err, strconv.ErrRange) {
		return 0, fmt.Errorf("%w: %q", ErrBadValue, s)
	}
	v *= scale
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, fmt.Errorf("%w: %q is not finite", ErrBadValue, s)
	}
	return v, nil
}

// numericPrefix returns the length of the longest prefix of s that forms a
// floating-point literal (without suffix).
func numericPrefix(s string) int {
	i := 0
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		i++
	}
	digits := 0
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
		digits++
	}
	if i < len(s) && s[i] == '.' {
		i++
		for i < len(s) && s[i] >= '0' && s[i] <= '9' {
			i++
			digits++
		}
	}
	if digits == 0 {
		return 0
	}
	// Exponent: e or E followed by optional sign and at least one digit.
	if i < len(s) && (s[i] == 'e' || s[i] == 'E') {
		j := i + 1
		if j < len(s) && (s[j] == '+' || s[j] == '-') {
			j++
		}
		k := j
		for k < len(s) && s[k] >= '0' && s[k] <= '9' {
			k++
		}
		if k > j {
			i = k
		}
	}
	return i
}
