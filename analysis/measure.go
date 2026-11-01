// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/output"
)

// ErrNoMeasurement means a measurement is undefined for the given data (for
// example a level that is never crossed). The GA maps it to a penalty.
var ErrNoMeasurement = errors.New("analysis: measurement not found")

// Transfer returns the AC transfer function H(f) = V(out)/V(ref); ref may be
// "0" to use V(out) itself (for a 1 V AC source).
func Transfer(r *output.ACResult, out, ref string) ([]complex128, error) {
	vo, ok := r.NodeVoltage(out)
	if !ok {
		return nil, ErrNoMeasurement
	}
	if ref == "0" || ref == "gnd" || ref == "" {
		return vo, nil
	}
	vr, ok := r.NodeVoltage(ref)
	if !ok {
		return nil, ErrNoMeasurement
	}
	h := make([]complex128, len(vo))
	for k := range vo {
		if vr[k] == 0 {
			return nil, ErrNoMeasurement
		}
		h[k] = vo[k] / vr[k]
	}
	return h, nil
}

// GainAt returns |H| at frequency f, interpolated linearly in log-frequency
// and in dB between sweep points.
func GainAt(freq []float64, h []complex128, f float64) (float64, error) {
	if len(freq) == 0 || f < freq[0] || f > freq[len(freq)-1] || f <= 0 {
		return 0, ErrNoMeasurement
	}
	for k := 1; k < len(freq); k++ {
		if f <= freq[k] {
			d0, d1 := output.DB(h[k-1]), output.DB(h[k])
			x := (math.Log(f) - math.Log(freq[k-1])) / (math.Log(freq[k]) - math.Log(freq[k-1]))
			return math.Pow(10, (d0+x*(d1-d0))/20), nil
		}
	}
	return cmplx.Abs(h[0]), nil
}

// Cutoff3dB finds the frequencies where |H| falls 3.0103 dB (a factor √2)
// below its maximum, on each side of the maximum. fLow is 0 when the
// response does not fall below the threshold at low frequencies (low-pass);
// fHigh is +Inf when it does not fall at high frequencies (high-pass).
// Crossings are interpolated in log-frequency and dB.
func Cutoff3dB(freq []float64, h []complex128) (fLow, fHigh float64, err error) {
	if len(freq) < 2 || len(h) != len(freq) {
		return 0, 0, ErrNoMeasurement
	}
	peak, kp := 0.0, 0
	for k, z := range h {
		if a := cmplx.Abs(z); a > peak {
			peak, kp = a, k
		}
	}
	if peak == 0 {
		return 0, 0, ErrNoMeasurement
	}
	thr := output.DB(complex(peak, 0)) - 10*math.Log10(2)
	cross := func(k0, k1 int) float64 {
		d0, d1 := output.DB(h[k0]), output.DB(h[k1])
		x := (thr - d0) / (d1 - d0)
		l0, l1 := math.Log(freq[k0]), math.Log(freq[k1])
		return math.Exp(l0 + x*(l1-l0))
	}
	fLow, fHigh = 0, math.Inf(1)
	for k := kp; k > 0; k-- {
		if output.DB(h[k-1]) < thr {
			fLow = cross(k, k-1)
			break
		}
	}
	for k := kp; k < len(h)-1; k++ {
		if output.DB(h[k+1]) < thr {
			fHigh = cross(k, k+1)
			break
		}
	}
	return fLow, fHigh, nil
}

// Crossing returns the first time after tFrom at which y crosses level in
// the given direction (+1 rising, −1 falling, 0 either), by linear
// interpolation between time points.
func Crossing(t, y []float64, level float64, dir int, tFrom float64) (float64, error) {
	for k := 1; k < len(t); k++ {
		if t[k] < tFrom {
			continue
		}
		a, b := y[k-1]-level, y[k]-level
		rising := a < 0 && b >= 0
		falling := a > 0 && b <= 0
		if (dir >= 0 && rising) || (dir <= 0 && falling) {
			tc := t[k-1] + (t[k]-t[k-1])*a/(a-b)
			if tc >= tFrom {
				return tc, nil
			}
		}
	}
	return 0, ErrNoMeasurement
}

// RiseTime returns the time for y to go from lo to hi fraction (e.g. 0.1 and
// 0.9) of the step from its first value to its final value.
func RiseTime(t, y []float64, lo, hi float64) (float64, error) {
	if len(y) < 2 || lo >= hi {
		return 0, ErrNoMeasurement
	}
	y0, y1 := y[0], y[len(y)-1]
	if y1 == y0 {
		return 0, ErrNoMeasurement
	}
	dir := 1
	if y1 < y0 {
		dir = -1
	}
	t1, err := Crossing(t, y, y0+lo*(y1-y0), dir, t[0])
	if err != nil {
		return 0, err
	}
	t2, err := Crossing(t, y, y0+hi*(y1-y0), dir, t1)
	if err != nil {
		return 0, err
	}
	return t2 - t1, nil
}

// Overshoot returns the peak excursion beyond the final value as a fraction
// of the step from the first to the final value (0 if none).
func Overshoot(y []float64) (float64, error) {
	if len(y) < 2 {
		return 0, ErrNoMeasurement
	}
	y0, y1 := y[0], y[len(y)-1]
	step := y1 - y0
	if step == 0 {
		return 0, ErrNoMeasurement
	}
	worst := 0.0
	for _, v := range y {
		if over := (v - y1) / step; over > worst {
			worst = over
		}
	}
	return worst, nil
}

// SettlingTime returns the last time y is outside a band of ±tol (fraction
// of the step) around its final value, measured from t[0].
func SettlingTime(t, y []float64, tol float64) (float64, error) {
	if len(y) < 2 || tol <= 0 {
		return 0, ErrNoMeasurement
	}
	y1 := y[len(y)-1]
	band := tol * math.Abs(y1-y[0])
	if band == 0 {
		return 0, ErrNoMeasurement
	}
	for k := len(y) - 1; k >= 0; k-- {
		if math.Abs(y[k]-y1) > band {
			if k == len(y)-1 {
				return 0, ErrNoMeasurement
			}
			return t[k+1] - t[0], nil
		}
	}
	return 0, nil
}
