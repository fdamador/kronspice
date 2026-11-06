// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import "math"

// The light: the standard account of single-bubble sonoluminescence is
// thermal emission from the hot, weakly ionized gas at the collapse,
// optically thin (S. Hilgenfeldt, S. Grossmann, D. Lohse, "A simple
// explanation of light emission in sonoluminescence," Nature 398, 402,
// 1999; Brenner, Hilgenfeldt & Lohse 2002, §IV). The estimate here keeps
// its simplest part:
//
//   - ionization of argon from the Saha equation at the adiabatic core
//     temperature and the mean atom density, χ = 15.7596 eV (NIST ASD),
//     partition functions 1 (Ar) and 4 + 2·exp(−0.1775 eV/kT) (Ar⁺ ²P₃/₂,
//     ²P₁/₂, NIST ASD); no lowering of the ionization potential;
//   - electron–ion free–free (bremsstrahlung) emission,
//     ε = 1.4×10⁻⁴⁰·√T·n_e·n_i·ḡ W/m³ with ḡ = 1.2 (G. B. Rybicki,
//     A. P. Lightman, Radiative Processes in Astrophysics, Wiley, 1979,
//     eq. 5.15b, converted to SI);
//   - a uniform core, every photon escaping.
//
// It leaves out electron–neutral bremsstrahlung, recombination radiation,
// heat conduction out of the gas and the lowered ionization potential,
// which Hilgenfeldt et al. include; so it is an order-of-magnitude
// estimate, compared below with the measured photon numbers and flash
// widths. The light is not taken out of the dynamics: it is far below the
// other ledger terms, and its size is reported beside them.

const (
	eV        = 1.602176634e-19
	chiAr     = 15.7596 * eV
	me        = 9.1093837015e-31
	hPlanck   = 6.62607015e-34
	gauntFF   = 1.2
	ffCoeffSI = 1.4e-40 // W·m³·K^(−1/2)
)

// Ionization returns the electron density [m⁻³] of argon at atom density n
// and temperature T from the Saha equation (charge neutrality n_e = n_i).
func Ionization(n, T float64) float64 {
	kt := kB * T
	ui := 4 + 2*math.Exp(-0.1775*eV/kt)
	s := 2 * ui * math.Pow(2*math.Pi*me*kt/(hPlanck*hPlanck), 1.5) * math.Exp(-chiAr/kt) // n_e·n_i/n_0
	// n_e² = s·(n − n_e)
	return (-s + math.Sqrt(s*s+4*s*n)) / 2
}

// Light is the flash estimate.
type Light struct {
	Energy  float64   // emitted energy per flash [J]
	Peak    float64   // peak power [W]
	FWHM    float64   // full width at half maximum of the power [s]
	PeakIon float64   // ionization fraction at the peak
	Photons float64   // Energy / mean photon energy of the free–free spectrum (kT at the peak)
	Power   []float64 // emitted power at each sample [W]
}

// FlashEstimate integrates the free–free emission over one cycle's samples.
func (b Bubble) FlashEstimate(s []Sample) Light {
	N := b.Atoms()
	var L Light
	pw := make([]float64, len(s))
	ipk := 0
	for i, x := range s {
		V := 4 * math.Pi / 3 * x.R * x.R * x.R
		n := N / V
		ne := Ionization(n, x.Temp)
		pw[i] = ffCoeffSI * math.Sqrt(x.Temp) * ne * ne * gauntFF * V
		if pw[i] > pw[ipk] {
			ipk = i
			L.PeakIon = ne / n
		}
		if i > 0 {
			L.Energy += 0.5 * (pw[i] + pw[i-1]) * (x.T - s[i-1].T)
		}
	}
	L.Peak = pw[ipk]
	// FWHM by linear interpolation on either side of the peak
	half := L.Peak / 2
	cross := func(dir int) float64 {
		for i := ipk; i+dir >= 0 && i+dir < len(pw); i += dir {
			if pw[i+dir] < half {
				a, c := pw[i], pw[i+dir]
				return s[i].T + (s[i+dir].T-s[i].T)*(a-half)/(a-c)
			}
		}
		return s[ipk].T
	}
	L.FWHM = cross(1) - cross(-1)
	L.Photons = L.Energy / (kB * s[ipk].Temp)
	L.Power = pw
	return L
}

// Spectral form of the same emission (Rybicki & Lightman eq. 5.14b in SI):
//
//	dP/dV dν = 6.8×10⁻⁵¹ · n_e·n_i · T^(−1/2) · e^(−hν/kT) · ḡ  W/(m³·Hz)
//
// whose integral over ν is the 1.4×10⁻⁴⁰ √T law used above. The spectrum
// of the whole flash is its time integral; per unit wavelength,
// dE/dλ = (dE/dν)·c/λ².
const ffNuCoeffSI = 6.8e-51

// Spectrum is the flash's time-integrated spectral energy [J/m] at the
// given wavelengths [m], and the power in a band [W] at each sample.
func (b Bubble) Spectrum(s []Sample, lambdas []float64) []float64 {
	N := b.Atoms()
	out := make([]float64, len(lambdas))
	prev := make([]float64, len(lambdas))
	for i, x := range s {
		V := 4 * math.Pi / 3 * x.R * x.R * x.R
		ne := Ionization(N/V, x.Temp)
		pre := ffNuCoeffSI * ne * ne / math.Sqrt(x.Temp) * gauntFF * V
		for k, l := range lambdas {
			nu := cLight / l
			p := pre * math.Exp(-hPlanck*nu/(kB*x.Temp)) * cLight / (l * l)
			if i > 0 {
				out[k] += 0.5 * (p + prev[k]) * (x.T - s[i-1].T)
			}
			prev[k] = p
		}
	}
	return out
}

const cLight = 299792458.0

// BandWidth is the full width at half maximum [s] of the flash seen
// through a band filter from lo to hi [m] (flat transmission).
func (b Bubble) BandWidth(s []Sample, lo, hi float64) float64 {
	const nb = 16
	ls := make([]float64, nb)
	for k := range ls {
		ls[k] = lo + (hi-lo)*(float64(k)+0.5)/nb
	}
	N := b.Atoms()
	pw := make([]float64, len(s))
	ipk := 0
	for i, x := range s {
		V := 4 * math.Pi / 3 * x.R * x.R * x.R
		ne := Ionization(N/V, x.Temp)
		pre := ffNuCoeffSI * ne * ne / math.Sqrt(x.Temp) * gauntFF * V
		for _, l := range ls {
			pw[i] += pre * math.Exp(-hPlanck*cLight/(l*kB*x.Temp)) * cLight / (l * l)
		}
		if pw[i] > pw[ipk] {
			ipk = i
		}
	}
	half := pw[ipk] / 2
	cross := func(dir int) float64 {
		for i := ipk; i+dir >= 0 && i+dir < len(pw); i += dir {
			if pw[i+dir] < half {
				a, c := pw[i], pw[i+dir]
				return s[i].T + (s[i+dir].T-s[i].T)*(a-half)/(a-c)
			}
		}
		return s[ipk].T
	}
	return cross(1) - cross(-1)
}

// PlanckFit returns the temperature [K] whose Planck curve B_λ(T) best
// matches the shape of a spectrum (least squares in log, free scale),
// searched from 2,000 to 100,000 K.
func PlanckFit(lambdas, spec []float64) float64 {
	cost := func(T float64) float64 {
		var lr []float64
		mean := 0.0
		for k, l := range lambdas {
			bb := 1 / (math.Pow(l, 5) * (math.Exp(hPlanck*cLight/(l*kB*T)) - 1))
			r := math.Log(spec[k] / bb)
			lr = append(lr, r)
			mean += r
		}
		mean /= float64(len(lr))
		c := 0.0
		for _, r := range lr {
			c += (r - mean) * (r - mean)
		}
		return c
	}
	lo, hi := math.Log(2000.0), math.Log(1e5)
	for it := 0; it < 200; it++ { // golden section on log T
		a := hi - (hi-lo)*0.618
		bb := lo + (hi-lo)*0.618
		if cost(math.Exp(a)) < cost(math.Exp(bb)) {
			hi = bb
		} else {
			lo = a
		}
	}
	return math.Exp((lo + hi) / 2)
}
