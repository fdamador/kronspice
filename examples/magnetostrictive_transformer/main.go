// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command magnetostrictive_transformer: a transformer whose core is a
// magnetostrictive material carrying ultrasonic standing waves from two
// transducers, as a mechanical relative of Barrett's oscillator-shuttle
// mixer (register entry H37). Terfenol-D and Galfenol cores, as a rod, a
// ring and a helical spring, with a DC bias coil; the primary is driven at
// 5 kHz, the two transducers at the core's modes, and the secondary shows
// what reaches it: the transformer line, the ultrasound read through the
// Villari effect, and the products of the two (the nonlinear part). The
// swap test drives the secondary instead and watches the primary.
package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
)

// Case is one core's runs.
type Case struct {
	Core           Core
	Fwd, Rev       Result
	LinFwd, LinRev Result // small drive, transducers off: reciprocity
	NoSound        Result // forward drive with the transducers off
}

// RunCase runs one core.
func RunCase(c Core) (Case, error) {
	var k Case
	k.Core = c
	r := DefaultRig()
	var err error
	if k.Fwd, err = c.Run(r); err != nil {
		return k, err
	}
	rr := r
	rr.Reverse = true
	if k.Rev, err = c.Run(rr); err != nil {
		return k, err
	}
	q := r
	q.Acoustic = false
	if k.NoSound, err = c.Run(q); err != nil {
		return k, err
	}
	q.Vdrive = 1e-3
	if k.LinFwd, err = c.Run(q); err != nil {
		return k, err
	}
	q.Reverse = true
	if k.LinRev, err = c.Run(q); err != nil {
		return k, err
	}
	return k, nil
}

// Reciprocity is |I_b(fwd) − I_a(rev)| / |I_b(fwd)| at fe for the small drive.
func (k Case) Reciprocity() float64 {
	a, b := k.LinFwd.Lines["fe"], k.LinRev.Lines["fe"]
	d := a - b
	return math.Hypot(real(d), imag(d)) / math.Hypot(real(a), imag(a))
}

func main() {
	data := flag.String("data", "", "write figure data (CSV) to this directory")
	flag.Parse()
	if err := run(*data); err != nil {
		fmt.Fprintln(os.Stderr, "magnetostrictive_transformer:", err)
		os.Exit(1)
	}
}

func run(data string) error {
	r := DefaultRig()
	fmt.Println("H37 — magnetostrictive transformer with two ultrasonic transducers (Terfenol-D, Galfenol; rod, ring, helix)")
	fmt.Println("  materials (declared, at the bias Hb = 0.66 a):")
	for _, m := range []Material{Terfenol, TerfenolComposite, Galfenol} {
		d, mur, kk := m.Small()
		fmt.Printf("    %-22s d = %4.1f nm/A, μr = %5.1f, k = %.2f, E = %.0f GPa, c = %.0f m/s, Q = %.0f, bias %.0f kA/m\n",
			m.Name, d*1e9, mur, kk, m.E/1e9, math.Sqrt(m.E/m.Rho), m.Q, m.Bias()/1e3)
	}
	fmt.Printf("  coils %.0f and %.0f turns wound along the core; primary port %.0f Ω, secondary (load) %.0f Ω; drive %.0f V at %.0f kHz;\n",
		r.N1, r.N2, r.Ra, r.Rb, r.Vdrive, r.Fe/1e3)
	fmt.Printf("  two transducers of %.0f N each at the core's coupled modes\n", r.Force)
	var cs []Case
	for _, c := range Cores() {
		k, err := RunCase(c)
		if err != nil {
			return err
		}
		cs = append(cs, k)
	}
	fmt.Println()
	fmt.Println("  forward (drive the primary, watch the secondary): load current lines [mA]")
	fmt.Println("  core                         tones fa, fb [kHz]   H peak/a   fe       fa       fb      | fa±fe          fb±fe          fa+fb    fb−fa   | ledger")
	for _, k := range cs {
		f := k.Fwd
		fmt.Printf("  %-28s %5.1f, %5.1f        %5.2f    %7.3f  %7.4f  %7.4f | %.1e %.1e  %.1e %.1e  %.1e %.1e | %.0e\n",
			k.Core.Mat.Name+" "+k.Core.Shape.String(), f.Fa/1e3, f.Fb/1e3, f.Hac/k.Core.Mat.A,
			1e3*f.Amp("fe"), 1e3*f.Amp("fa"), 1e3*f.Amp("fb"),
			1e3*f.Amp("fa+fe"), 1e3*f.Amp("fa-fe"), 1e3*f.Amp("fb+fe"), 1e3*f.Amp("fb-fe"), 1e3*f.Amp("fa+fb"), 1e3*f.Amp("fb-fa"),
			f.Residual/f.In)
	}
	fmt.Println()
	fmt.Println("  swap test: drive the secondary (same source in series with the 20 Ω side), watch the primary")
	fmt.Println("  core                         reciprocity (small drive)   fe fwd → rev [mA]     products fa±fe, fb±fe summed: fwd → rev [mA]")
	for _, k := range cs {
		sum := func(r Result) float64 {
			return r.Amp("fa+fe") + r.Amp("fa-fe") + r.Amp("fb+fe") + r.Amp("fb-fe")
		}
		fmt.Printf("  %-28s %.1e                    %7.3f → %7.3f     %.2e → %.2e\n",
			k.Core.Mat.Name+" "+k.Core.Shape.String(), k.Reciprocity(), 1e3*k.Fwd.Amp("fe"), 1e3*k.Rev.Amp("fe"), 1e3*sum(k.Fwd), 1e3*sum(k.Rev))
	}
	fmt.Println()
	if data != "" {
		if err := writeData(data, cs); err != nil {
			return err
		}
	}
	return nil
}

func writeData(dir string, cs []Case) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, "lines.csv"))
	if err != nil {
		return err
	}
	defer f.Close()
	names := []string{"fe", "fa", "fb", "fa+fe", "fa-fe", "fb+fe", "fb-fe", "fa+fb", "fb-fa", "2fe"}
	fmt.Fprint(f, "core,direction")
	for _, n := range names {
		fmt.Fprint(f, ","+n)
	}
	fmt.Fprintln(f)
	for _, k := range cs {
		for _, d := range []struct {
			n string
			r Result
		}{{"forward", k.Fwd}, {"reverse", k.Rev}, {"no sound", k.NoSound}} {
			fmt.Fprintf(f, "%s %s,%s", k.Core.Mat.Name, k.Core.Shape, d.n)
			for _, n := range names {
				fmt.Fprint(f, ","+strconv.FormatFloat(d.r.Amp(n), 'g', 6, 64))
			}
			fmt.Fprintln(f)
		}
		w, err := os.Create(filepath.Join(dir, fmt.Sprintf("wave_%s_%s.csv", short(k.Core.Mat.Name), k.Core.Shape)))
		if err != nil {
			return err
		}
		fmt.Fprintln(w, "t,with_sound,no_sound")
		for j := range k.Fwd.T {
			fmt.Fprintf(w, "%g,%g,%g\n", k.Fwd.T[j], k.Fwd.Iout[j], k.NoSound.Iout[j])
		}
		w.Close()
	}
	return nil
}

func short(s string) string {
	switch s {
	case Terfenol.Name:
		return "terfenol"
	case Galfenol.Name:
		return "galfenol"
	default:
		return "composite"
	}
}
