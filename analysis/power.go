// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"errors"
	"fmt"
	"math"
	"sort"

	"github.com/fdamador/kronspice/hyper/ga"
)

// Power in non-sinusoidal conditions, two ways (CLAUDE.md §3b item 3,
// M9f): geometric algebra, and the conventional quantities of IEEE Std
// 1459-2010 ("Definitions for the Measurement of Electric Power Quantities
// Under Sinusoidal, Nonsinusoidal, Balanced, or Unbalanced Conditions").
//
// Geometric algebra (M. Castro-Núñez, R. Castro-Puche, "The IEEE Standard
// 1459, the CPC power theory, and geometric algebra in circuits with
// nonsinusoidal sources and linear loads," IEEE Trans. Circuits Syst. I
// 59(12), 2012; F. G. Montoya et al., Mathematics 9(11):1259, 2021). Each
// harmonic k contributes two orthonormal functions over a period,
// σ_ck = √2·cos kωt and σ_sk = √2·sin kωt, the basis vectors e_{2k−1},
// e_{2k} of the Euclidean algebra Cl(2H). A periodic voltage is the vector
// V = Σ (V_ck·σ_ck + V_sk·σ_sk) with RMS components, so the mean of v·i
// is the dot product V·I. The geometric product
//
//	M = V·I = V⌋I + V∧I
//
// has a scalar part, the active power P, and a bivector part. Its
// σ_ck∧σ_sk components (the same harmonic) are the reactive powers
// Q_k = V_ck·I_sk − V_sk·I_ck — with the phasor convention
// x(t) = Re{√2·X·e^(jkωt)}, Q_k = Im(V_k·I_k*), positive for an inductive
// load. The components between different harmonics are distortion power,
// each tied to its pair of harmonics. Since |M| = |V|·|I| for vectors,
//
//	S² = (V_rms·I_rms)² = P² + Σ_k Q_k² + Σ_{h≠k} D_hk²,
//
// and grouping V and I into fundamental and harmonic parts reproduces the
// IEEE 1459 decomposition S² = S₁² + D_I² + D_V² + S_H²: |V₁·I₁| = S₁,
// |V₁·I_H| = D_I, |V_H·I₁| = D_V, |V_H·I_H| = S_H.

// ErrPowerWindow is a waveform too short for the requested periods.
var ErrPowerWindow = errors.New("analysis: waveform shorter than the power-analysis window")

// Harmonics are the RMS cosine and sine components of a voltage and a
// current at k·F0, k = 1 … len(V) (index 0 is k = 1), and their DC parts.
type Harmonics struct {
	F0         float64
	V, I       [][2]float64 // [k−1] = (C_k, S_k), x = Σ C_k·√2cos + S_k·√2sin
	VDC, IDC   float64
	VRMS, IRMS float64 // total RMS over the window from the samples themselves
	Samples    int     // uniform samples per window used
}

// HarmonicsOf resamples v(t) and i(t) uniformly over the last `periods`
// periods of 1/f0 (samples per window: at least 64·hmax and 4096) and
// takes their components up to harmonic hmax (of f0) by the discrete
// Fourier transform; the RMS values are from the same samples, so for
// band-limited waveforms Parseval makes them agree with the components.
// A DC component is reported separately (VDC, IDC) and is not a GA basis
// vector.
func HarmonicsOf(t, v, i []float64, f0 float64, periods, hmax int) (Harmonics, error) {
	if len(t) < 2 || len(v) != len(t) || len(i) != len(t) || !(f0 > 0) || periods < 1 || hmax < 1 {
		return Harmonics{}, fmt.Errorf("analysis: HarmonicsOf: bad arguments")
	}
	T := float64(periods) / f0
	t1 := t[len(t)-1]
	t0 := t1 - T
	if t0 < t[0]-1e-12*T {
		return Harmonics{}, fmt.Errorf("%w: need %g s, have %g s", ErrPowerWindow, T, t1-t[0])
	}
	n := 64 * hmax * periods
	if n < 4096 {
		n = 4096
	}
	vs, is := resample(t, v, t0, T, n), resample(t, i, t0, T, n)
	h := Harmonics{F0: f0, V: make([][2]float64, hmax), I: make([][2]float64, hmax), Samples: n}
	var sv, si, sv2, si2 float64
	for s := 0; s < n; s++ {
		sv += vs[s]
		si += is[s]
		sv2 += vs[s] * vs[s]
		si2 += is[s] * is[s]
	}
	h.VDC, h.IDC = sv/float64(n), si/float64(n)
	h.VRMS, h.IRMS = math.Sqrt(sv2/float64(n)), math.Sqrt(si2/float64(n))
	for k := 1; k <= hmax; k++ {
		m := k * periods // harmonic k of f0 is bin m of the window
		var vc, vsn, ic, isn float64
		for s := 0; s < n; s++ {
			sn, cs := math.Sincos(2 * math.Pi * float64(m) * float64(s) / float64(n))
			vc += vs[s] * cs
			vsn += vs[s] * sn
			ic += is[s] * cs
			isn += is[s] * sn
		}
		// RMS components on √2cos, √2sin: (√2/n)·Σ x·cos.
		f := math.Sqrt2 / float64(n)
		h.V[k-1] = [2]float64{f * vc, f * vsn}
		h.I[k-1] = [2]float64{f * ic, f * isn}
	}
	return h, nil
}

// resample samples x(t) linearly at t0 + s·T/n, s = 0 … n−1 (one window,
// the endpoint excluded so the window is exactly periodic).
func resample(t, x []float64, t0, T float64, n int) []float64 {
	out := make([]float64, n)
	j := 0
	for s := 0; s < n; s++ {
		ts := t0 + T*float64(s)/float64(n)
		for j < len(t)-2 && t[j+1] < ts {
			j++
		}
		dt := t[j+1] - t[j]
		a := 0.0
		if dt > 0 {
			a = (ts - t[j]) / dt
		}
		out[s] = x[j] + a*(x[j+1]-x[j])
	}
	return out
}

// HarmonicPair is the distortion power between two harmonics: the norm
// of the bivector components that couple harmonic H with harmonic K (the
// four products of their cosine and sine functions), in VA.
type HarmonicPair struct {
	H, K  int // harmonic orders, H < K
	Value float64
}

// GAPower is the geometric power M = V·I of one voltage–current pair.
type GAPower struct {
	P float64 // scalar part: active power (AC part; DC power is PDC)
	Q float64 // √Σ Q_k²: the same-harmonic bivector part
	D float64 // the cross-harmonic bivector part
	S float64 // |M| = |V|·|I|
	// PerHarmonic are P_k and Q_k for k = 1 … H.
	PerHarmonic []struct{ P, Q float64 }
	// Distortion lists the distortion power of each pair of harmonics,
	// largest first: which pairs cause it. Σ Value² = D².
	Distortion []HarmonicPair
	// IEEE-1459 groups from the same product: S1 = |V₁·I₁|,
	// DI = |V₁·I_H|, DV = |V_H·I₁|, SH = |V_H·I_H|.
	S1, DI, DV, SH float64
	// PDC is the DC power VDC·IDC (outside the AC algebra).
	PDC float64
	// M is the multivector itself, in Alg.
	M   ga.Multivector
	Alg *ga.Algebra
}

// GeometricPower forms M = V·I in Cl(2H) (H ≤ 32 harmonics).
func GeometricPower(h Harmonics) (*GAPower, error) {
	H := len(h.V)
	if H < 1 || H > 32 || len(h.I) != H {
		return nil, fmt.Errorf("analysis: GeometricPower: 1 … 32 harmonics, got %d", H)
	}
	alg, err := ga.New(2*H, 0, 0)
	if err != nil {
		return nil, err
	}
	vec := func(x [][2]float64, from, to int) ga.Multivector {
		c := make([]float64, 2*H)
		for k := from; k < to; k++ {
			c[2*k], c[2*k+1] = x[k][0], x[k][1]
		}
		m, _ := alg.Vector(c...)
		return m
	}
	V, I := vec(h.V, 0, H), vec(h.I, 0, H)
	M := alg.Mul(V, I)
	out := &GAPower{M: M, Alg: alg, P: M.Coef(0), PDC: h.VDC * h.IDC}
	out.S = alg.Norm(M)
	out.PerHarmonic = make([]struct{ P, Q float64 }, H)
	var q2, d2 float64
	pairs := map[[2]int]float64{}
	for _, b := range M.Blades() {
		if b.Grade() != 2 {
			continue
		}
		c := M.Coef(b)
		lo := trailing(uint64(b))
		hi := trailing(uint64(b) &^ (1 << lo))
		kl, kh := lo/2+1, hi/2+1
		if kl == kh {
			// σ_ck∧σ_sk: Q_k.
			out.PerHarmonic[kl-1].Q = c
			q2 += c * c
			continue
		}
		d2 += c * c
		pairs[[2]int{kl, kh}] += c * c
	}
	for hk, v := range pairs {
		out.Distortion = append(out.Distortion, HarmonicPair{H: hk[0], K: hk[1], Value: math.Sqrt(v)})
	}
	for k := 0; k < H; k++ {
		out.PerHarmonic[k].P = h.V[k][0]*h.I[k][0] + h.V[k][1]*h.I[k][1]
	}
	out.Q, out.D = math.Sqrt(q2), math.Sqrt(d2)
	sort.Slice(out.Distortion, func(i, j int) bool {
		a, b := out.Distortion[i], out.Distortion[j]
		if a.Value != b.Value {
			return a.Value > b.Value
		}
		return a.H < b.H || a.H == b.H && a.K < b.K
	})
	V1, VH := vec(h.V, 0, 1), vec(h.V, 1, H)
	I1, IH := vec(h.I, 0, 1), vec(h.I, 1, H)
	out.S1 = alg.Norm(alg.Mul(V1, I1))
	out.DI = alg.Norm(alg.Mul(V1, IH))
	out.DV = alg.Norm(alg.Mul(VH, I1))
	out.SH = alg.Norm(alg.Mul(VH, IH))
	return out, nil
}

func trailing(x uint64) int {
	n := 0
	for x&1 == 0 && n < 64 {
		x >>= 1
		n++
	}
	return n
}

// IEEE1459 are the single-phase quantities of IEEE Std 1459-2010 §3.2,
// computed conventionally: RMS values from the time samples, the
// fundamental phasors from the DFT, and the harmonic parts by subtraction.
type IEEE1459 struct {
	V, I, V1, I1, VH, IH float64 // RMS (AC parts: DC removed)
	P, P1, PH            float64 // active: total, fundamental, harmonic
	Q1                   float64 // fundamental reactive power V1·I1·sin θ1
	S, S1, SN            float64 // apparent: total, fundamental, non-fundamental
	DI, DV, SH           float64 // current and voltage distortion powers, harmonic apparent power
	THDV, THDI           float64 // VH/V1, IH/I1
	PF, PF1              float64 // P/S and P1/S1
}

// IEEE1459Of computes the IEEE 1459 quantities from the same uniform
// samples as HarmonicsOf (it calls it): V and I are RMS values of the
// samples with the DC part removed, P the mean of (v − VDC)(i − IDC).
func IEEE1459Of(t, v, i []float64, f0 float64, periods int) (IEEE1459, error) {
	h, err := HarmonicsOf(t, v, i, f0, periods, 1)
	if err != nil {
		return IEEE1459{}, err
	}
	T := float64(periods) / f0
	n := h.Samples
	vs, is := resample(t, v, t[len(t)-1]-T, T, n), resample(t, i, t[len(t)-1]-T, T, n)
	var q IEEE1459
	var sv, si, sp float64
	for s := 0; s < n; s++ {
		dv, di := vs[s]-h.VDC, is[s]-h.IDC
		sv += dv * dv
		si += di * di
		sp += dv * di
	}
	q.V, q.I, q.P = math.Sqrt(sv/float64(n)), math.Sqrt(si/float64(n)), sp/float64(n)
	v1, i1 := complex(h.V[0][0], -h.V[0][1]), complex(h.I[0][0], -h.I[0][1]) // RMS phasors
	q.V1, q.I1 = cabs(v1), cabs(i1)
	s1 := v1 * complex(real(i1), -imag(i1))
	q.P1, q.Q1 = real(s1), imag(s1)
	q.VH = math.Sqrt(math.Max(q.V*q.V-q.V1*q.V1, 0))
	q.IH = math.Sqrt(math.Max(q.I*q.I-q.I1*q.I1, 0))
	q.PH = q.P - q.P1
	q.S, q.S1 = q.V*q.I, q.V1*q.I1
	q.SN = math.Sqrt(math.Max(q.S*q.S-q.S1*q.S1, 0))
	q.DI, q.DV, q.SH = q.V1*q.IH, q.VH*q.I1, q.VH*q.IH
	if q.V1 > 0 {
		q.THDV = q.VH / q.V1
	}
	if q.I1 > 0 {
		q.THDI = q.IH / q.I1
	}
	if q.S > 0 {
		q.PF = q.P / q.S
	}
	if q.S1 > 0 {
		q.PF1 = q.P1 / q.S1
	}
	return q, nil
}

func cabs(z complex128) float64 { return math.Hypot(real(z), imag(z)) }
