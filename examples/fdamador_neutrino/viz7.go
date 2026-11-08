// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Data for viz/render_part7.py (Part 7, spin.png), as the run computed it.

import (
	"fmt"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/spacetime"
)

func writeViz7(dir string, p7a Part7a, p7b []Part7b) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	// 7a: δ_CP, J/J_max and the angles against ν (Part 6 defects), both signs.
	var nu, dcp, jr, t12, t23, t13, ratio []float64
	for e := -6.0; e <= -1.5+1e-9; e += 0.05 {
		for _, sgn := range []float64{-1, 1} {
			x := sgn * math.Pow(10, e)
			c := RotClusterMix(0.2064, 0.1628, 2, 23, x, [3]int{0, 1, 2})
			nu, dcp, jr = append(nu, x), append(dcp, c.DeltaCP), append(jr, c.J/c.JMax)
			t12, t23, t13, ratio = append(t12, c.Th12), append(t23, c.Th23), append(t13, c.Th13), append(ratio, c.Ratio)
		}
	}
	if err := writeCSV(dir, "p7_scan.csv", "nu,delta_cp,j_over_jmax,th12,th23,th13,ratio", nu, dcp, jr, t12, t23, t13, ratio); err != nil {
		return err
	}
	m := NuFIT6()
	meta := fmt.Sprintf("key,value\nfit_nu,%.9g\nfit_delta,%.6g\npinned_nu,%.9g\npinned_delta,%.6g\nburst_nu,%.9g\n"+
		"th12,%.6g\ns12,%.6g\nth23,%.6g\ns23,%.6g\nth13,%.6g\ns13,%.6g\nratio,%.6g\nsratio,%.6g\ndelta,%.6g\n"+
		"tmeas,%.6g\ntpred,%.6g\ntzero,%.6g\n",
		p7a.Fit.Nu, p7a.Fit.DeltaCP, p7a.Pinned.Nu, p7a.Pinned.DeltaCP, 2*spacetime.BurstSpeed(spacetime.T1000G)/C0,
		m.Th12, m.S12, m.Th23, m.S23, m.Th13, m.S13, m.Ratio, m.SRatio, m.DeltaCP, p7a.TMeas, p7a.TPred, p7a.TZeroMeas)
	if err := os.WriteFile(filepath.Join(dir, "p7_meta.csv"), []byte(meta), 0o644); err != nil {
		return err
	}
	// T violation: P(e→μ) and P(μ→e), with and without rotation, on the e run's times.
	e, mu := p7a.Pinned.Assign[0], p7a.Pinned.Assign[1]
	thin := func(a, b RotThree) (t, pem, pme, pemP, pmeP []float64) {
		r := DefaultFRing(false)
		w := r.OmegaK(r.Kappa(r.N/2 - 1))
		for i := 0; i < len(a.TAll); i += 10 {
			tt := a.TAll[i]
			t = append(t, tt*w/(2*math.Pi))
			pem, pemP = append(pem, a.PAll[mu][i]), append(pemP, a.PredAll[mu][i])
			pme, pmeP = append(pme, interp(b.TAll, b.PAll[e], tt)), append(pmeP, interp(b.TAll, b.PredAll[e], tt))
		}
		return
	}
	t, a1, a2, a3, a4 := thin(p7a.ThreeE, p7a.ThreeMu)
	if err := writeCSV(dir, "p7_tviol.csv", "periods,p_e_mu,p_mu_e,p_e_mu_pred,p_mu_e_pred", t, a1, a2, a3, a4); err != nil {
		return err
	}
	t, a1, a2, a3, a4 = thin(p7a.ZeroE, p7a.ZeroMu)
	if err := writeCSV(dir, "p7_tviol_zero.csv", "periods,p_e_mu,p_mu_e,p_e_mu_pred,p_mu_e_pred", t, a1, a2, a3, a4); err != nil {
		return err
	}
	// 7b: per run of the finest grid.
	if len(p7b) > 0 {
		fine := p7b[len(p7b)-1]
		if err := os.WriteFile(filepath.Join(dir, "p7b_grid.csv"), []byte(fmt.Sprintf("cells\n%g\n", fine.Cells)), 0o644); err != nil {
			return err
		}
		names := "index,name,q,jzw0,hw0\n"
		for i, r := range fine.Runs {
			u := make([]float64, len(r.U))
			for k := range r.U {
				u[k] = r.U[k] / r.U[0]
			}
			if err := writeCSV(dir, fmt.Sprintf("p7b_run%d.csv", i), "periods,u_over_u0,jz_w_over_u,h_w_over_u,eb,jz_mink_w_over_u", r.T, u, r.JzW, r.HW, r.EB, r.JzMW); err != nil {
				return err
			}
			names += fmt.Sprintf("%d,%s,%.6g,%.6g,%.6g\n", i, r.Name, r.Q, r.JzW0, r.HW0)
		}
		if err := os.WriteFile(filepath.Join(dir, "p7b_runs.csv"), []byte(names), 0o644); err != nil {
			return err
		}
	}
	return nil
}
