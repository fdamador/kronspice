// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package em

import (
	"errors"
	"fmt"
	"math"
)

// Multiconductor lossless lines: modes and transition points.
//
// L. V. Bewley, Tensor Analysis of Electric Circuits and Machines, Ronald
// Press 1961, Ch. 6 "Multivelocity Traveling Waves" ("Bewley 1961"); C. R.
// Paul, Analysis of Multiconductor Transmission Lines, 2nd ed., Wiley-IEEE
// 2008, Ch. 7 ("Paul 2008"). For n conductors over a reference with
// per-unit-length inductance L and capacitance K (Bewley's K = C = p⁻¹,
// Maxwell's potential coefficients inverted), R = G = 0, the waves travel
// in n modes whose velocities solve |L·K·v² − 1| = 0 (Bewley Eq. 6‑11):
// v_i = 1/√λ_i for the eigenvalues λ_i of L·K. In a homogeneous medium
// L·K = 1/c² and every mode travels at c (Eq. 6‑15). The surge
// (characteristic) impedance matrix is
//
//	Z_c = K^(−½)·(K^(½)·L·K^(½))^(½)·K^(−½)       (Paul 2008 §7.2)
//
// so Z_c·K·Z_c = L and Y_c = Z_c⁻¹ is Bewley's surge admittance (Eq. 6‑14).
//
// Where a line meets a network whose impedance seen from the line's end is
// Z′ (a Thévenin matrix, found by Kron's elimination, Bewley Eq. 3‑18), an
// incident wave e is reflected as e_r = R·e with (Bewley Eqs. 6‑22 to
// 6‑25)
//
//	R = (1 + Z′·Y_c)⁻¹·(Z′·Y_c − 1),   total voltage (1 + R)·e.
//
// These are the pieces a SPICE-style multiconductor T element needs
// (modal decoupling, then one delay per mode); the element itself is not
// in the circuit engine yet (doc/bewley.md).

// MTLModes is the modal solution of a lossless multiconductor line.
type MTLModes struct {
	// Velocity of each mode [m/s], ascending.
	Velocity []float64
	// T holds the modal voltage vectors as columns: L·K·T = T·diag(1/v²).
	T [][]float64
	// Zc and Yc are the characteristic impedance and admittance matrices.
	Zc, Yc [][]float64
}

// LosslessModes solves the modes of a lossless line with per-unit-length
// inductance l [H/m] and capacitance k [F/m] (both symmetric positive
// definite, n×n).
func LosslessModes(l, k [][]float64) (*MTLModes, error) {
	n := len(l)
	if n == 0 || len(k) != n {
		return nil, fmt.Errorf("em: LosslessModes: L is %d×?, K is %d×?", n, len(k))
	}
	for i := 0; i < n; i++ {
		if len(l[i]) != n || len(k[i]) != n {
			return nil, errors.New("em: LosslessModes: L and K must be square and the same size")
		}
		for j := 0; j < n; j++ {
			if math.Abs(l[i][j]-l[j][i]) > 1e-12*math.Abs(l[i][i]) || math.Abs(k[i][j]-k[j][i]) > 1e-12*math.Abs(k[i][i]) {
				return nil, errors.New("em: LosslessModes: L and K must be symmetric")
			}
		}
	}
	kh, khi, err := sqrtSPD(k)
	if err != nil {
		return nil, fmt.Errorf("em: LosslessModes: K: %w", err)
	}
	s := matMul(matMul(kh, l), kh)
	lam, u, err := symEig(s)
	if err != nil {
		return nil, err
	}
	for _, v := range lam {
		if v <= 0 {
			return nil, errors.New("em: LosslessModes: L is not positive definite")
		}
	}
	// Order modes by ascending velocity (descending λ).
	idx := make([]int, n)
	for i := range idx {
		idx[i] = i
	}
	for i := 1; i < n; i++ {
		for j := i; j > 0 && lam[idx[j]] > lam[idx[j-1]]; j-- {
			idx[j], idx[j-1] = idx[j-1], idx[j]
		}
	}
	m := &MTLModes{Velocity: make([]float64, n)}
	us := make([][]float64, n)
	sroot := make([][]float64, n)
	for i := range us {
		us[i] = make([]float64, n)
		sroot[i] = make([]float64, n)
	}
	for c, j := range idx {
		m.Velocity[c] = 1 / math.Sqrt(lam[j])
		for i := 0; i < n; i++ {
			us[i][c] = u[i][j]
		}
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			var v float64
			for q := 0; q < n; q++ {
				v += u[i][q] * math.Sqrt(lam[q]) * u[j][q]
			}
			sroot[i][j] = v
		}
	}
	m.T = matMul(khi, us)
	m.Zc = matMul(matMul(khi, sroot), khi)
	m.Yc, err = invert(m.Zc)
	if err != nil {
		return nil, errors.New("em: LosslessModes: singular surge impedance")
	}
	return m, nil
}

// Reflection is Bewley's transition-point matrix R (Eqs. 6‑22 to 6‑25):
// a wave e arriving on a line of surge impedance zc at a network of
// impedance zterm (seen from the line's end) is reflected as R·e.
func Reflection(zc, zterm [][]float64) ([][]float64, error) {
	n := len(zc)
	if n == 0 || len(zterm) != n {
		return nil, errors.New("em: Reflection: matrices must be the same size")
	}
	yc, err := invert(zc)
	if err != nil {
		return nil, errors.New("em: Reflection: singular surge impedance")
	}
	zy := matMul(zterm, yc)
	a := make([][]float64, n)
	b := make([][]float64, n)
	for i := 0; i < n; i++ {
		a[i] = make([]float64, n)
		b[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			a[i][j], b[i][j] = zy[i][j], zy[i][j]
		}
		a[i][i]++
		b[i][i]--
	}
	ai, err := invert(a)
	if err != nil {
		return nil, errors.New("em: Reflection: singular transition (1 + Z′·Y)")
	}
	return matMul(ai, b), nil
}

func matMul(a, b [][]float64) [][]float64 {
	n, m, p := len(a), len(b), len(b[0])
	out := make([][]float64, n)
	for i := 0; i < n; i++ {
		out[i] = make([]float64, p)
		for k := 0; k < m; k++ {
			if a[i][k] == 0 {
				continue
			}
			for j := 0; j < p; j++ {
				out[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return out
}

// sqrtSPD returns A^(½) and A^(−½) of a symmetric positive definite A.
func sqrtSPD(a [][]float64) (h, hi [][]float64, err error) {
	lam, u, err := symEig(a)
	if err != nil {
		return nil, nil, err
	}
	n := len(a)
	h = make([][]float64, n)
	hi = make([][]float64, n)
	for i := 0; i < n; i++ {
		if lam[i] <= 0 {
			return nil, nil, errors.New("not positive definite")
		}
		h[i] = make([]float64, n)
		hi[i] = make([]float64, n)
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			var s, si float64
			for q := 0; q < n; q++ {
				r := math.Sqrt(lam[q])
				s += u[i][q] * r * u[j][q]
				si += u[i][q] / r * u[j][q]
			}
			h[i][j], hi[i][j] = s, si
		}
	}
	return h, hi, nil
}

// symEig diagonalizes a symmetric matrix by cyclic Jacobi rotations:
// a = u·diag(lam)·uᵀ, eigenvectors in the columns of u.
func symEig(a0 [][]float64) (lam []float64, u [][]float64, err error) {
	n := len(a0)
	a := make([][]float64, n)
	u = make([][]float64, n)
	for i := range a {
		a[i] = append([]float64(nil), a0[i]...)
		u[i] = make([]float64, n)
		u[i][i] = 1
	}
	for sweep := 0; sweep < 100; sweep++ {
		off, scale := 0.0, 0.0
		for i := 0; i < n; i++ {
			scale += a[i][i] * a[i][i]
			for j := i + 1; j < n; j++ {
				off += a[i][j] * a[i][j]
			}
		}
		if off <= 1e-30*scale || off == 0 {
			lam = make([]float64, n)
			for i := range lam {
				lam[i] = a[i][i]
			}
			return lam, u, nil
		}
		for p := 0; p < n; p++ {
			for q := p + 1; q < n; q++ {
				if a[p][q] == 0 {
					continue
				}
				th := (a[q][q] - a[p][p]) / (2 * a[p][q])
				t := math.Copysign(1, th) / (math.Abs(th) + math.Sqrt(th*th+1))
				c := 1 / math.Sqrt(t*t+1)
				s := t * c
				for k := 0; k < n; k++ {
					akp, akq := a[k][p], a[k][q]
					a[k][p], a[k][q] = c*akp-s*akq, s*akp+c*akq
				}
				for k := 0; k < n; k++ {
					apk, aqk := a[p][k], a[q][k]
					a[p][k], a[q][k] = c*apk-s*aqk, s*apk+c*aqk
				}
				for k := 0; k < n; k++ {
					ukp, ukq := u[k][p], u[k][q]
					u[k][p], u[k][q] = c*ukp-s*ukq, s*ukp+c*ukq
				}
			}
		}
	}
	return nil, nil, errors.New("em: eigenvalues did not converge")
}
