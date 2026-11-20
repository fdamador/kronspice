// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"fmt"
	"math"
	"path/filepath"
	"sort"
)

// RegenResult is one regenerative-load run with its plain-load twin.
type RegenResult struct {
	N, C3Factor float64
	R           Result
}

// Regen runs design A with its load replaced by the regenerative load:
// turns ratio n from 3 to 6 and C3 retuned by 0.8–1.2×. It returns every
// run, the best (largest power through the load port) and design A with
// its plain load for comparison.
func Regen(ctx context.Context) (all []RegenResult, best RegenResult, plain Result) {
	var ds []Design
	var keys [][2]float64
	for n := 3.0; n <= 6.0+1e-9; n += 0.25 {
		for _, f := range []float64{0.8, 0.9, 1.0, 1.1, 1.2} {
			d := EvolvedSync()
			d.Regen, d.RegenN, d.C3 = true, n, d.C3*f
			ds = append(ds, d)
			keys = append(keys, [2]float64{n, f})
		}
	}
	ds = append(ds, EvolvedSync())
	rs := parallel(ctx, ds, false)
	plain = rs[len(rs)-1]
	for i, r := range rs[:len(rs)-1] {
		x := RegenResult{N: keys[i][0], C3Factor: keys[i][1], R: r}
		all = append(all, x)
		if r.Err == nil && r.PPort > best.R.PPort {
			best = x
		}
	}
	return
}

func runRegen(ctx context.Context, dir string) error {
	all, best, plain := Regen(ctx)
	fmt.Println("== Regenerative load on design A (synchronized) ==")
	fmt.Println("  turns ratio n:1 and C3 retune; P_port = power into the load port, P_return = power back into the 24 V rail")
	var rows []string
	for _, x := range all {
		r := x.R
		if r.Err != nil {
			fmt.Printf("  n %.2f C3×%.1f: %v\n", x.N, x.C3Factor, r.Err)
			continue
		}
		net := r.PInHV + r.PInRoy
		fmt.Printf("  n %.2f C3×%.1f  P_port %7.1f W  P_return %7.1f W (%.3f)  supplies net %6.1f W (HV %5.1f, Royer %6.1f; Royer before return %6.1f)  mismatch %.1e\n",
			x.N, x.C3Factor, r.PPort, r.PReturn, safeDiv(r.PReturn, r.PPort), net, r.PInHV, r.PInRoy, r.PRoyGross, r.Mismatch)
		rows = append(rows, fmt.Sprintf("%g,%g,%g,%g,%g,%g,%g,%g", x.N, x.C3Factor, r.PPort, r.PReturn, r.PInHV, r.PInRoy, r.PRoyGross, r.Mismatch))
	}
	r := best.R
	net := r.PInHV + r.PInRoy
	fmt.Printf("\n  best: n %.2f, C3 %.3g F (×%.1f)\n", best.N, best.R.Design.C3, best.C3Factor)
	fmt.Printf("  power through the load port %.1f W, returned to the 24 V rail %.1f W (recovery %.3f)\n", r.PPort, r.PReturn, r.PReturn/r.PPort)
	fmt.Printf("  the supplies now deliver %.1f W net (HV %.1f W + Royer %.1f W); the Royer supply delivers %.1f W and takes %.1f W back\n",
		net, r.PInHV, r.PInRoy, r.PRoyGross, r.PReturn)
	fmt.Printf("  with the plain %g Ω load: %.1f W to the load from %.1f W in (η %.3f), losses %.1f W\n",
		plain.Design.RL, plain.POut, plain.PInHV+plain.PInRoy, plain.Eta, plain.PInHV+plain.PInRoy-plain.POut)
	fmt.Printf("  ledger: integration mismatch %.1e of input; nothing is kept at the load, so net input = losses\n", r.Mismatch)
	printLosses(r)
	if dir == "" {
		return nil
	}
	if err := writeCSV(filepath.Join(dir, "regen.csv"), "n,c3factor,p_port,p_return,p_hv,p_roy_net,p_roy_gross,mismatch", rows); err != nil {
		return err
	}
	var lr []string
	keys := make([]string, 0, len(r.PDiss))
	for k := range r.PDiss {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		lr = append(lr, fmt.Sprintf("%s,%g", k, r.PDiss[k]))
	}
	lr = append(lr, fmt.Sprintf("PPort,%g", r.PPort), fmt.Sprintf("PReturn,%g", r.PReturn), fmt.Sprintf("PInHV,%g", r.PInHV),
		fmt.Sprintf("PInRoy,%g", r.PInRoy), fmt.Sprintf("PRoyGross,%g", r.PRoyGross), fmt.Sprintf("N,%g", best.N),
		fmt.Sprintf("C3,%g", r.Design.C3), fmt.Sprintf("Mismatch,%g", r.Mismatch),
		fmt.Sprintf("PlainPOut,%g", plain.POut), fmt.Sprintf("PlainPIn,%g", plain.PInHV+plain.PInRoy))
	return writeCSV(filepath.Join(dir, "regen_best.csv"), "name,value", lr)
}

func safeDiv(a, b float64) float64 {
	if b == 0 || math.IsNaN(b) {
		return 0
	}
	return a / b
}
