// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
)

// WriteViz writes the data behind the report's figures:
//
//	calib.csv        Xue's rotor: σ and the collapse rate rpm·g(σ) against speed
//	gcurve.csv       the vapour law g(σ)
//	front.csv        the dimple search's feasible Pareto designs
//	sensitivity.csv  best η with each declared constant changed
//	ring.csv         |p| of the gap ring (kronspice AC), clear water
//	collapse.csv     R(t), Ṙ(t), radiated pressure of the best design's cavity
//	drives.csv       energy books, motor and pump + turbine
//	pern.csv         best Mandelbulb η for each n
//	bulb_n*.csv      equator profiles of the n = 3, 8, 12 bulbs
//	best.csv         the best dimple and bulb designs
func WriteViz(dir string, p Params, seed uint64) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	w := func(name string, head []string, rows [][]float64) error {
		f, err := os.Create(filepath.Join(dir, name))
		if err != nil {
			return err
		}
		defer f.Close()
		c := csv.NewWriter(f)
		_ = c.Write(head)
		for _, r := range rows {
			s := make([]string, len(r))
			for i, v := range r {
				s[i] = strconv.FormatFloat(v, 'g', 8, 64)
			}
			_ = c.Write(s)
		}
		c.Flush()
		return c.Error()
	}
	fmt.Println("writing", dir)

	// 1–2: calibration
	var cb [][]float64
	for rpm := 2500.0; rpm <= 6000; rpm += 25 {
		sg := Sigma(Xue2024.PRef, omega(rpm)*Xue2024.R)
		cb = append(cb, []float64{rpm, sg, rpm * p.G(sg)})
	}
	if err := w("calib.csv", []string{"rpm", "sigma", "rate"}, cb); err != nil {
		return err
	}
	var gc [][]float64
	for sg := 0.0; sg <= 1.4; sg += 0.01 {
		gc = append(gc, []float64{sg, p.G(sg)})
	}
	if err := w("gcurve.csv", []string{"sigma", "g"}, gc); err != nil {
		return err
	}

	// 4: front and sensitivity
	front, err := Evolve(p, dimpleSpace, DimpleFrom, dimpleLand, 80, 80, seed, nil)
	if err != nil || len(front) == 0 {
		return fmt.Errorf("GA: %v", err)
	}
	var fr [][]float64
	for _, s := range front {
		fr = append(fr, []float64{s.E.Eta, s.E.PCollapse, s.E.PShaft, s.E.Sigma, s.X[0], s.X[2], s.X[7], s.X[8]})
	}
	if err := w("front.csv", []string{"eta", "Pcollapse", "Pshaft", "sigma", "R", "W", "gap", "rpm"}, fr); err != nil {
		return err
	}
	best := Heater(front)
	var se [][]float64
	names := []string{}
	for _, r := range Sensitivity(p, seed) {
		if r.Err == nil {
			se = append(se, []float64{float64(len(names)), r.Best.E.Eta})
			names = append(names, r.Name)
		}
	}
	if err := w("sensitivity.csv", []string{"i", "eta"}, se); err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(dir, "sensitivity_names.txt"), []byte(fmt.Sprint(names)), 0o644); err != nil {
		return err
	}

	// 5: ring spectrum
	m, e := best.M, best.E
	bandLen := float64(m.Rows) * m.RowLen
	f0 := cLiq / (2 * math.Pi * m.R)
	fq, mag, err := RingSpectrum(m.R, m.Gap, bandLen, cLiq, 64, 0.3*f0, 4.3*f0, 4000)
	if err != nil {
		return err
	}
	var rg [][]float64
	for i := range fq {
		rg = append(rg, []float64{fq[i], mag[i]})
	}
	if err := w("ring.csv", []string{"f", "abs_p"}, rg); err != nil {
		return err
	}

	// 6: collapse trace
	c0 := m.Rotor[0][0]
	pulse := 0.0
	if len(m.Housing) > 0 {
		pulse = p.Cpulse * 0.5 * rho * e.U * e.U * math.Pow(p.G0/m.Gap, p.GapExp)
	}
	dp := e.P - pVap + pulse
	r0 := math.Cbrt(3 * p.KappaV * e.G * c0.V / (4 * math.Pi))
	var tr [][]float64
	k := 0
	RunCollapseTrace(r0, dp, 1e-3, 0.0717, p.Mu, func(t, r, v, prad float64) {
		if k%20 == 0 {
			tr = append(tr, []float64{t, r, v, prad / m.R})
		}
		k++
	})
	if err := w("collapse.csv", []string{"t", "R", "Rdot", "p_at_axis"}, tr); err != nil {
		return err
	}

	// 7: drives
	var dr [][]float64
	for i, d := range Drives(e) {
		dr = append(dr, []float64{float64(i), d.PElec, d.PShaft, d.PWaterHeat, d.PRoomHeat, d.PCollapse, d.CavPerJoule, d.SteamKgH, d.HeatPumpW})
	}
	if err := w("drives.csv", []string{"drive(0 motor,1 turbine)", "Pelec", "Pshaft", "water_heat", "room_heat", "Pcollapse", "cav_per_J", "steam_kg_h", "heat_pump_W"}, dr); err != nil {
		return err
	}

	// 8: per n and profiles
	var pn [][]float64
	for _, r := range PerN(p, seed) {
		pn = append(pn, []float64{float64(r.N), r.Eta, float64(r.Cavities)})
	}
	if err := w("pern.csv", []string{"n", "eta", "cavities"}, pn); err != nil {
		return err
	}
	for _, n := range []float64{3, 8, 12} {
		prof := Bulb{N: n, MaxIt: bulbMaxIt}.Profile(0, 720)
		var pr [][]float64
		for i, r := range prof {
			pr = append(pr, []float64{2 * math.Pi * float64(i) / 720, r})
		}
		if err := w(fmt.Sprintf("bulb_n%d.csv", int(n)), []string{"phi", "r"}, pr); err != nil {
			return err
		}
	}
	bulbFront, err := Evolve(p, bulbSpace, BulbFrom, bulbLand, 60, 60, seed, nil)
	if err != nil || len(bulbFront) == 0 {
		return fmt.Errorf("bulb GA: %v", err)
	}
	bb := bulbFront[0]

	// 10: closed steam loop
	var le [][]float64
	for _, pd := range []float64{1.18e5, 1.5e5, 3e5} {
		_, rpm := pathEta(p, m, pd, 60, 60)
		mm := m
		mm.POut, mm.RPM = pd, rpm
		for t := 20.0; t <= 135; t += 2.5 {
			q := p
			q.PVap = SatPressure(t + 273.15)
			eta := math.NaN()
			if q.PVap < pd {
				if ev, err := q.Evaluate(mm); err == nil {
					eta = ev.Eta
				}
			}
			le = append(le, []float64{pd, t, eta})
		}
	}
	if err := w("loop_eta.csv", []string{"p_drum", "T_C", "eta"}, le); err != nil {
		return err
	}
	var lp [][]float64
	for i, r := range Loops(p, best) {
		lp = append(lp, []float64{float64(i), r.SteamKgH, r.PTurbine, r.PMotorShaft, r.PElec, r.PHome, r.PRoom, r.PCollapse, r.EtaLoop, r.CavPerElecJ, r.FracSensible})
	}
	if err := w("loop.csv", []string{"cycle(0 A,1 B)", "steam_kg_h", "Pturbine", "Pmotor_shaft", "Pelec", "Phome", "Proom", "Pcollapse", "eta_loop", "cav_per_elecJ", "frac_sensible"}, lp); err != nil {
		return err
	}
	if err := w("best.csv", []string{"kind(0 dimple,1 bulb)", "eta", "Pcollapse", "Pshaft", "Pfriction", "Pdrag", "sigma", "U", "cavities"}, [][]float64{
		{0, e.Eta, e.PCollapse, e.PShaft, e.PFriction, e.PDrag, e.Sigma, e.U, float64(e.Cavities)},
		{1, bb.E.Eta, bb.E.PCollapse, bb.E.PShaft, bb.E.PFriction, bb.E.PDrag, bb.E.Sigma, bb.E.U, float64(bb.E.Cavities)},
	}); err != nil {
		return err
	}
	return nil
}
