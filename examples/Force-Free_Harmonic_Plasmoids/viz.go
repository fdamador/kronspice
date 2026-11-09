// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// -viz dir: the numbers behind the pictures (viz/render.py draws them).

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

func writeViz(dir string) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	w := func(name string, body func(f *os.File)) error {
		f, err := os.Create(filepath.Join(dir, name))
		if err != nil {
			return err
		}
		defer f.Close()
		body(f)
		return nil
	}
	// lifetimes vs radius (1–50 cm), both recipes, two energies
	if err := w("lifetimes.csv", func(f *os.File) {
		fmt.Fprintln(f, "R_m,W_J,tauTE_recipe,tauTM_recipe,beta_recipe,f_recipe,tauTE_nonrel,tauTM_nonrel,f_nonrel")
		for i := 0; i <= 40; i++ {
			r := 0.01 * math.Pow(50, float64(i)/40)
			for _, wj := range []float64{1e4, 1e6} {
				a := DesignFromC(r, wj, 0.01)
				b := DesignFromBeta(r, wj, 0.1)
				fmt.Fprintf(f, "%g,%g,%g,%g,%g,%g,%g,%g,%g\n", r, wj, a.TauTE, a.TauTM, a.Beta, a.F, b.TauTE, b.TauTM, b.F)
			}
		}
	}); err != nil {
		return err
	}
	// electron density needed to last tau
	if err := w("needed.csv", func(f *os.File) {
		fmt.Fprintln(f, "R_m,tau_s,kR,f_Hz,ne_over_loschmidt")
		for i := 0; i <= 40; i++ {
			r := 0.01 * math.Pow(50, float64(i)/40)
			for _, tau := range []float64{1e-6, 1e-3, 1, 3, 10} {
				y, fr, ne := ForLifetime(r, tau)
				fmt.Fprintf(f, "%g,%g,%g,%g,%g\n", r, tau, y, fr, ne/Loschmidt)
			}
		}
	}); err != nil {
		return err
	}
	// harmonic comb at 12 cm
	if err := w("comb.csv", func(f *os.File) {
		fmt.Fprintln(f, "n,kR,f_Hz,ne,tauTE,tauTM")
		z := J1Zeros(60)
		for i := 1; i < len(z); i++ {
			d := DesignAt(0.12, 1e5, z[i])
			fmt.Fprintf(f, "%d,%g,%g,%g,%g,%g\n", i+1, z[i], d.F, d.Ne, d.TauTE, d.TauTM)
		}
	}); err != nil {
		return err
	}
	// kronspice spectrum of the Appendix B plasmoid: one sweep over all 8 modes
	al := J1Zero / 0.1
	e := AppendixB("printed", 35.87/(al*al), 2686)
	modes := Spectrum(e.R, e.Wp, 4)
	ref := real(modes[0].Omega)
	b := circuit.NewBuilder("spectrum")
	b.Add(circuit.Element{Name: "IGEN", Kind: circuit.CurrentSource, ACMag: 1}, "0", "n0")
	for i, m := range modes {
		wn := real(m.Omega) / ref
		a, c := fmt.Sprintf("n%d", i), fmt.Sprintf("n%d", i+1)
		if i == len(modes)-1 {
			c = "0"
		}
		b.Capacitor(fmt.Sprintf("C%d", i), a, c, 1)
		b.Inductor(fmt.Sprintf("L%d", i), a, c, 1/(wn*wn))
		b.Resistor(fmt.Sprintf("R%d", i), a, c, m.Q/wn)
	}
	ckt, err := b.Build()
	if err != nil {
		return err
	}
	lo, hi := math.Inf(1), 0.0
	for _, m := range modes {
		lo, hi = math.Min(lo, real(m.Omega)), math.Max(hi, real(m.Omega))
	}
	span := hi - lo
	r, err := analysis.AC(ckt, analysis.ACSpec{Variation: "lin", Points: 40001,
		FStart: (lo - 0.05*span) / ref / (2 * math.Pi), FStop: (hi + 0.05*span) / ref / (2 * math.Pi)}, analysis.Options{})
	if err != nil {
		return err
	}
	v0, _ := r.NodeVoltage("n0")
	if err := w("spectrum.csv", func(f *os.File) {
		fmt.Fprintln(f, "f_Hz,absZ_ohm_normalized")
		for k := range r.Freq {
			fmt.Fprintf(f, "%.12g,%g\n", r.Freq[k]*ref, cmplx.Abs(v0[k]))
		}
	}); err != nil {
		return err
	}
	if err := w("modes.csv", func(f *os.File) {
		fmt.Fprintln(f, "kind,f_Hz,Q,lifetime_s")
		for _, m := range modes {
			k := "TE"
			if !m.TE {
				k = "TM"
			}
			fmt.Fprintf(f, "%s,%.12g,%g,%g\n", k, real(m.Omega)/(2*math.Pi), m.Q, m.Lifetime)
		}
	}); err != nil {
		return err
	}
	// kronspice ring-down, shown with Q lowered to 40 so the decay is visible
	ck, err := circuit.NewBuilder("ring-down").CapacitorIC("C1", "a", "0", 1, 1).Inductor("L1", "a", "0", 1).
		Resistor("R1", "a", "0", 40).BuildTransient()
	if err != nil {
		return err
	}
	T := 2 * math.Pi
	wf, err := analysis.Tran(ck, analysis.TranSpec{Step: T / 64, Stop: 30 * T, MaxStep: T / 64, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return err
	}
	va, _ := wf.NodeVoltage("a")
	return w("ringdown.csv", func(f *os.File) {
		fmt.Fprintln(f, "t_periods,v,stored")
		for k := range wf.Time {
			fmt.Fprintf(f, "%g,%g,%g\n", wf.Time[k]/T, va[k], wf.Ledger.StoredEnergy[k])
		}
	})
}
