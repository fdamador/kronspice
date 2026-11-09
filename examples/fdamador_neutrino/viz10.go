// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render_part10.py (weak.png), as the run computed it.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
)

func writeViz10(dir string, p Part10) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	write := func(name string, b *strings.Builder) error {
		return os.WriteFile(filepath.Join(dir, name), []byte(b.String()), 0o644)
	}
	ring0, ring1 := p.Nu[0].MuEff2, p.Nu[1].MuEff2
	var b strings.Builder
	// recoil spectra
	b.WriteString("e_nu,t,weak_e,weak_mu,weak_ebar,dip_xenon,ring0_point,ring0,ring1\n")
	for _, ee := range []float64{0.862, 10} {
		tm := TMax(ee)
		for t := 1e-3; t < tm; t *= 1.04 {
			fmt.Fprintf(&b, "%g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", ee, t, WeakDSigmaDT(NuE, false, ee, t, 0), WeakDSigmaDT(NuMu, false, ee, t, 0),
				WeakDSigmaDT(NuE, true, ee, t, 0), DipoleDSigmaDT(boundXe*boundXe, ee, t, 0), DipoleDSigmaDT(ring0, ee, t, 0),
				DipoleDSigmaDT(ring0, ee, t, p.ANu), DipoleDSigmaDT(ring1, ee, t, p.ANu))
		}
	}
	if err := write("p10_recoil.csv", &b); err != nil {
		return err
	}
	// total cross-section against energy
	b.Reset()
	b.WriteString("e_nu,weak_e,weak_mu,weak_ebar,ring0,ring1\n")
	for ee := 0.2; ee <= 40; ee *= 1.15 {
		fmt.Fprintf(&b, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", ee, SigmaWeak(NuE, false, ee, 1e-3), SigmaWeak(NuMu, false, ee, 1e-3), SigmaWeak(NuE, true, ee, 1e-3),
			SigmaDipole(ring0, ee, 1e-3, p.ANu), SigmaDipole(ring1, ee, 1e-3, p.ANu))
	}
	if err := write("p10_sigma.csv", &b); err != nil {
		return err
	}
	// form factor and deuteron E1
	b.Reset()
	b.WriteString("x,ff2\n")
	for x := 0.01; x < 300; x *= 1.08 {
		fmt.Fprintf(&b, "%.6g,%.6g\n", x, FormFactor2(x))
	}
	if err := write("p10_formfactor.csv", &b); err != nil {
		return err
	}
	b.Reset()
	b.WriteString("w,sigma_e1\n")
	for w := bDeut + 0.01; w < 20; w += 0.05 {
		fmt.Fprintf(&b, "%.6g,%.6g\n", w, SigmaE1(w))
	}
	if err := write("p10_e1.csv", &b); err != nil {
		return err
	}
	// kronspice twin (every 4th point)
	b.Reset()
	b.WriteString("t_periods,i_ring,p_electron,work\n")
	T := 1 / p.Twin.F0
	var work float64
	for k := range p.Twin.Time {
		if k > 0 {
			work += 0.5 * (p.Twin.PElectron[k] + p.Twin.PElectron[k-1]) * (p.Twin.Time[k] - p.Twin.Time[k-1])
		}
		if k%4 == 0 {
			fmt.Fprintf(&b, "%.6g,%.6g,%.6g,%.6g\n", p.Twin.Time[k]/T, p.Twin.IRing[k], p.Twin.PElectron[k], work)
		}
	}
	if err := write("p10_twin.csv", &b); err != nil {
		return err
	}
	b.Reset()
	b.WriteString("key,value\n")
	kv := func(k string, v float64) { fmt.Fprintf(&b, "%s,%.8g\n", k, v) }
	kv("a_nu", p.ANu)
	kv("mu0_amp", p.Nu[0].MAmp/MuB)
	kv("mu0_rms", math.Sqrt(ring0))
	kv("pc1_amp", p.Nu[1].PAmp*C0/MuB)
	kv("eff_xe0", p.EffXe[0])
	kv("eff_xe1", p.EffXe[1])
	kv("eff_bx0", p.EffBx[0])
	kv("bound_xe", boundXe)
	kv("bound_bx", boundBx)
	kv("sno_ring0", p.SNORing[0])
	kv("sno_ring1", p.SNORing[1])
	kv("sno_point1", p.SNOPoint1)
	kv("sig_nc_ref", sigNCref)
	kv("anti_ratio", p.AntiRatioWeak)
	kv("r2_ring", p.Nu[0].R2Max/E*1e4)
	kv("q_ring", p.Nu[1].QMax/E)
	return write("p10_meta.csv", &b)
}
