// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command three_coil_sync: three coils, one output. Coil 1 is an air-core
// high-voltage step-up transformer whose HV terminal is wired to coil 3;
// coil 2 is the bifilar, centre-tapped tank coil of a resonant (ZVS,
// Mazzilli) Royer oscillator, wound over coil 3; coil 3 is a series
// resonant RLC loop feeding a load. The HV coil drives coil 3 with high
// voltage at low current; the Royer drives it through the magnetic
// coupling with high current at low voltage. The genetic algorithm
// evolves the coil turns, the spacing, the capacitors, the load and the
// HV generator's frequency and phase for the largest power in the load,
// which needs the two drives synchronized (same frequency, phase that
// adds in coil 3).
//
// Every inductance is computed from the winding geometry (Maxwell's
// coaxial-circle formula, turn by turn), so the evolved design is a set
// of turn counts and dimensions someone can wind. kronspice runs every
// candidate as a full transient, MOSFETs and diodes included, with the
// energy ledger on: the load's power is traced to the two supplies.
package main

import (
	"context"
	"flag"
	"fmt"
	"os"
)

// Baseline is a hand-picked starting design (not evolved).
func Baseline() Design {
	d := Design{N3: 40, N2: 6, Gap23: 5e-3, N1s: 400, N1p: 8, C2: 680e-9, C3: 100e-9, RL: 50, PhiHV: 0}
	d.FHV = d.RoyerFreq0()
	return d
}

func main() {
	dataDir := flag.String("data", "", "write CSV data for the report into this directory")
	gens := flag.Int("gens", 30, "GA generations")
	pop := flag.Int("pop", 48, "GA population")
	seed := flag.Uint64("seed", 1, "GA seed")
	quick := flag.Bool("quick", false, "baseline only")
	flag.BoolVar(&syncMode, "sync", false, "clock the HV amplifier from the Royer's tank (synchronized search)")
	waves := flag.Bool("waves", false, "only rewrite the waveform files and netlists of the two evolved designs into -data and -data-free")
	dataFree := flag.String("data-free", "", "with -waves: data directory of the free-running design")
	flag.Parse()
	ctx := context.Background()

	if *waves {
		for _, x := range []struct {
			dir string
			d   Design
		}{{*dataFree, EvolvedFree()}, {*dataDir, EvolvedSync()}} {
			if x.dir == "" {
				continue
			}
			if err := writeWaves(ctx, x.dir, x.d); err != nil {
				fmt.Fprintln(os.Stderr, err)
				os.Exit(1)
			}
		}
		return
	}
	b := Baseline()
	fmt.Println("== Baseline (hand-picked) ==")
	report(Simulate(ctx, b, false))
	if *quick {
		return
	}
	if err := run(ctx, *dataDir, *gens, *pop, *seed); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func report(r Result) {
	d := r.Design
	if r.Err != nil {
		fmt.Println("  error:", r.Err)
		return
	}
	e := r.E
	fmt.Printf("  coils: N3 %d, N2 %d per strand, gap %.1f mm, N1s %d, N1p %d; C2 %.3g F, C3 %.3g F, RL %.3g Ω\n",
		d.N3, d.N2, d.Gap23*1e3, d.N1s, d.N1p, d.C2, d.C3, d.RL)
	fmt.Printf("  L3 %.3g H, L2 strand %.3g H, k(2 strands) %.3f, k(2–3) %.3f, L1s %.3g H, L1p %.3g H, k1 %.3f, C1self %.3g F\n",
		e.L3, e.L2a, e.K2ab, e.K2a3, e.L1s, e.L1p, e.K1, e.C1self)
	if d.Sync {
		fmt.Printf("  HV amplifier slaved to the tank: phase %.0f°, gain %.3g, output %.0f V peak\n", d.SyncAngle(), d.SyncGain, r.VgPk)
	}
	fmt.Printf("  HV gen %.4g Hz, %.0f°; Royer tank (unloaded) %.4g Hz, measured %.4g Hz; output %.4g Hz; locked %v, beat %.3f\n",
		d.FHV, d.PhiHV, r.FR0, r.FRoyer, r.FOut, r.Locked, r.Beat)
	fmt.Printf("  P_out %.4g W; in: HV %.4g W + Royer %.4g W; η %.3f\n", r.POut, r.PInHV, r.PInRoy, r.Eta)
	fmt.Printf("  peaks: Vds %.1f V, Id %.1f A, V(hv) %.0f V, V(C3) %.0f V, V(load) %.0f V; I3 %.2f A rms, I2 %.2f A rms, supply %.2f A\n",
		r.VdsPk, r.IdPk, r.VHVPk, r.VC3Pk, r.VOutPk, r.I3rms, r.I2rms, r.IChoke)
	fmt.Printf("  ledger: max residual %.2g, integration mismatch %.2g of input; %d steps\n", r.Residual, r.Mismatch, r.Steps)
}

// EvolvedFree is the design the free-running search found (seed 1,
// population 48, 25 generations; run_free.txt, viz/data_free).
func EvolvedFree() Design {
	return Design{N3: 101, N2: 7, Gap23: 0.004238877084804275, N1s: 55, N1p: 12, C2: 2.2e-6, C3: 100e-9, RL: 24,
		FHV: 24016.988656209804, PhiHV: 157.25148166666554}
}

// EvolvedSync is the design the synchronized search found (flag -sync;
// seed 1, population 48, 25 generations; run_output.txt, viz/data).
func EvolvedSync() Design {
	d := Design{N3: 184, N2: 8, Gap23: 0.007253941848441254, N1s: 50, N1p: 4, C2: 4.7e-6, C3: 150e-9, RL: 24,
		Sync: true, SyncGain: 0.1670924419549027}
	d.SetSyncAngle(287.3050662877631)
	return d
}
