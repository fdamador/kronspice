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
	"strings"

	"github.com/fdamador/kronspice/analysis"
)

// spiceNetlist returns a plain SPICE3 netlist of a measured design, with
// the EMF switched on at t = 0 (PWL), an .ac sweep two decades either side
// of the resonance (the 1 A probe makes v(a) the tank impedance) and a
// .tran of the switch-on. It runs as is in kronspice, ngspice, LTspice and
// QSpice; the comment block maps every element to the physical part.
func spiceNetlist(r Result, label string) (string, error) {
	d := r.Design
	c, err := d.Build(r.Rcoil, true)
	if err != nil {
		return "", err
	}
	c.Title = fmt.Sprintf("Bi-Al %s element with %s coil, %s (kronspice examples/bial_bifilar)",
		strings.Fields(d.Variant.String())[0], coilKind(d.Coil.Bifilar), label)
	periods := math.Min(2000, math.Max(20, 3*r.Q))
	ac := fmt.Sprintf(".ac dec 200 %s %s", num(r.F0/100), num(r.F0*100))
	tran := fmt.Sprintf(".tran %s %s 0 %s", num(1/(40*r.F0)), num(periods/r.F0), num(1/(40*r.F0)))
	var body strings.Builder
	if err := c.WriteNetlist(&body, ac, tran); err != nil {
		return "", err
	}
	lines := strings.SplitN(body.String(), "\n", 2)
	el := r.Element
	var h strings.Builder
	p := func(f string, x ...any) { fmt.Fprintf(&h, "* "+f+"\n", x...) }
	p("")
	p("Coil (%s): %d turns, diameter %s, wire %s, pitch %s, insulation eps %.3g, copper",
		coilKind(d.Coil.Bifilar), d.Coil.Turns, eng(d.Coil.Diameter, "m"), eng(d.Coil.WireDiameter, "m"),
		eng(d.Coil.Pitch, "m"), d.Coil.InsulationEps)
	p("  Lcoil  = %s   Wheeler (Proc. IRE 1928)", eng(r.L, "H"))
	p("  Cself  = %s   turn-to-turn capacitance, %s", eng(r.Cself, "F"),
		map[bool]string{true: "V/2 per gap (bifilar, Tesla US 512,340)", false: "V/N per gap (ordinary)"}[d.Coil.Bifilar])
	p("  Rcoil  = %s   copper, skin effect at f0", eng(r.Rcoil, "Ohm"))
	p("Bi-Al element: plates %s square, gap %s", eng(d.PlateSide, "m"), eng(d.Gap, "m"))
	switch d.Variant {
	case Dry, Thermopile:
		pairs := 1
		if d.Variant == Thermopile {
			pairs = d.Pairs
		}
		p("  Vbial  = %s   Seebeck EMF %d x (S_Al - S_Bi)*dT = %d x 75.5 uV/K * %g K; energy from HEAT (thermal port)",
			eng(el.EMF, "V"), pairs, pairs, d.A.DeltaT)
		if d.Variant == Thermopile {
			p("  Rint   = %s   %d Bi and Al legs (%s dia. x %s) plus plates; legs conduct %s of heat at dT",
				eng(el.Rint, "Ohm"), 2*d.Pairs, eng(d.LegDiameter, "m"), eng(d.LegLength, "m"), eng(el.HeatLeak, "W"))
		} else {
			p("  Rint   = %s   Bi and Al plate resistance", eng(el.Rint, "Ohm"))
		}
		p("  Cplate = %s   plates with dielectric eps %.3g, in series with the EMF (DC is blocked)",
			eng(el.Cplate, "F"), d.DielectricEps)
	case Wet:
		p("  Vbial  = %s   Al|electrolyte|Bi cell EMF (%g x standard 1.970 V); energy from CHEMISTRY (chemical port)",
			eng(el.EMF, "V"), d.A.EMFFactor)
		p("  Relec  = %s   electrolyte, rho %g Ohm*m", eng(el.Rint, "Ohm"), d.A.ElectrolyteRho)
		p("  Rct    = %s   charge transfer, two electrodes (assumed %g Ohm*cm2 each)", eng(el.Rct, "Ohm"), d.A.ChargeTransfer*1e4)
		p("  Cdl    = %s   double layers, two in series (%g uF/cm2 each)", eng(el.Cdl, "F"), d.A.DoubleLayer*100)
		p("  Cgeo   = %s   geometric plate capacitance, water eps %g", eng(el.Cgeo, "F"), d.A.WaterEps)
	}
	p("Iprobe: 1 A AC current into node a (DC 0), so vm(a) = |Z| of the tank in the .ac run.")
	p("Vbial switches on at t = 0 (PWL; an idealized step - a real thermal or chemical source is slower).")
	p("Simulated: f0 = %s, Q = %.4g, stored energy at steady state %s (%s per volume).",
		eng(r.F0, "Hz"), r.Q, eng(r.StoredJ, "J"), eng(r.EnergyDensity, "J/m^3"))
	p("Note: SPICE itself books Vbial as an ordinary source; kronspice's ledger attributes its energy to the %s port.",
		d.channelName())
	p("ngspice: run, then  plot vm(a)  (ac)  or  plot v(a) i(Vbial)  (tran)")
	p("")
	// Plain ASCII, so every simulator and editor reads the comments.
	header := strings.NewReplacer("µ", "u", "–", "-", "³", "^3", "²", "^2").Replace(h.String())
	return lines[0] + "\n" + header + lines[1], nil
}

func num(v float64) string { return strconv.FormatFloat(v, 'g', 6, 64) }

// readFront reads a front CSV written by writeCSV and returns the gene
// vectors.
func readFront(path string, v Variant) ([][]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		return nil, err
	}
	var out [][]float64
	for _, row := range rows[1:] {
		x := make([]float64, len(geneNames(v)))
		for i := range x {
			if x[i], err = strconv.ParseFloat(row[i], 64); err != nil {
				return nil, fmt.Errorf("%s: %v", path, err)
			}
		}
		out = append(out, x)
	}
	return out, nil
}

// writeNetlists writes SPICE netlists (.cir), equivalent-circuit drawings
// (.svg) and editable LTspice schematics (.asc) of the highest-Q and the highest
// energy-per-volume designs of each front CSV in dir (front_dry*.csv,
// front_wet*.csv), plus the highest-Q dry design with an ordinary coil for
// comparison. It returns the files written.
func writeNetlists(dir string, a Assumptions) ([]string, error) {
	var written []string
	for _, v := range variants {
		name := v.short()
		matches, _ := filepath.Glob(filepath.Join(dir, "front_"+name+"*.csv"))
		if len(matches) == 0 {
			return written, fmt.Errorf("no front_%s*.csv in %s", name, dir)
		}
		xs, err := readFront(matches[0], v)
		if err != nil {
			return written, err
		}
		var hiQ, hiE Result
		for _, x := range xs {
			r, err := Measure(designFrom(v, a, true, x), analysis.Options{})
			if err != nil {
				continue
			}
			if r.Q > hiQ.Q {
				hiQ = r
			}
			if r.EnergyDensity > hiE.EnergyDensity {
				hiE = r
			}
		}
		type job struct {
			r     Result
			label string
			file  string
		}
		jobs := []job{{hiQ, "highest-Q design", name + "_highest_q.cir"}, {hiE, "highest energy-per-volume design", name + "_highest_energy.cir"}}
		if v == Dry {
			if o, err := Measure(withBifilar(hiQ.Design, false), analysis.Options{}); err == nil {
				jobs = append(jobs, job{o, "highest-Q design wound ordinarily (comparison)", name + "_highest_q_ordinary.cir"})
			}
		}
		for _, j := range jobs {
			s, err := spiceNetlist(j.r, j.label)
			if err != nil {
				return written, err
			}
			path := filepath.Join(dir, j.file)
			if err := os.WriteFile(path, []byte(s), 0o644); err != nil {
				return written, err
			}
			svgPath := strings.TrimSuffix(path, ".cir") + ".svg"
			if err := os.WriteFile(svgPath, []byte(schematicSVG(j.r, true)), 0o644); err != nil {
				return written, err
			}
			ascText, err := schematicASC(j.r, j.label)
			if err != nil {
				return written, err
			}
			ascPath := strings.TrimSuffix(path, ".cir") + ".asc"
			if err := os.WriteFile(ascPath, []byte(ascText), 0o644); err != nil {
				return written, err
			}
			written = append(written, path, svgPath, ascPath)
		}
	}
	return written, nil
}
