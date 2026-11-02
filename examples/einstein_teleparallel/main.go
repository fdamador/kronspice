// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command einstein_teleparallel puts A. Einstein, "New possibility for a
// unified field theory of gravitation and electricity," Sitzungsber.
// Preuss. Akad. Wiss. 1928, 224–227 (A. Unzicker's translation), beside
// Kron's machines: does Einstein's distant parallelism answer Kron's
// multi-dimensional, non-Riemannian work?
//
//  1. Einstein's 4-bein and torsion are Kron's non-holonomic frame and
//     Hoffmann's Ω: computed both ways on Kron's motor.
//
//  2. Einstein's electromagnetic potential φ = trace of the torsion, on the
//     same motor in two frames.
//
//  3. Curvature: Einstein's connection has none; the machine's space does.
//
//  4. Einstein's first approximation (eqs. 5–9) checked numerically.
//
//  5. The motor built in kronspice from Einstein's torsion alone, against
//     Kron's frames and against the torsion removed.
//
//     go run ./examples/einstein_teleparallel [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/kron"
)

func maxDiff(a, b [4][4][4]float64, sign float64) float64 {
	d := 0.0
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				d = math.Max(d, math.Abs(a[i][j][k]-sign*b[i][j][k]))
			}
		}
	}
	return d
}

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()
	f := Frame{}
	th := 0.7

	fmt.Println("1. Einstein's 4-bein is Kron's frame (Hoffmann's single-phase motor: q1, q2, q3, θ — four coordinates)")
	A := f.Anholonomy(th)
	E := f.ToFrame(f.EinsteinTorsion(th), th)
	fmt.Printf("   Einstein's torsion Λ (coordinate form, his eq. 9–10 of June 7) in the frame's components equals minus\n")
	fmt.Printf("   the frame's anholonomy (sign convention): difference %.0e.\n", maxDiff(A, E, -1))
	omS, err := kron.Torsion(Motor(), J, f.Moving(), th)
	if err != nil {
		log.Fatal(err)
	}
	var om [4][4][4]float64
	for i := range om {
		for j := range om[i] {
			copy(om[i][j][:], omS[i][j])
		}
	}
	lb := f.NewMetric(th)
	var low [4][4][4]float64
	for l := 0; l < 4; l++ {
		for m := 0; m < 4; m++ {
			for n := 0; n < 4; n++ {
				for r := 0; r < 4; r++ {
					low[l][m][n] += lb[n][r] * A[r][l][m]
				}
			}
		}
	}
	fmt.Printf("   Hoffmann's Ω (kron.Torsion, eq. 19) is that torsion with one index lowered by the inductances: Ω_λμν = L̄_νρ·Λ^ρ_λμ,\n")
	fmt.Printf("   difference %.0e. Same object: Kron's non-Riemannian frame is a teleparallel 4-bein.\n", maxDiff(om, low, 1))
	fmt.Printf("   Einstein's own connection has torsion 2Λ (check %.0e).\n", TorsionOfW(f, th))

	fmt.Println("\n2. Einstein's electromagnetic potential φ_μ = Λ^α_μα (his eq. 2) on the motor")
	p0 := Phi(f.EinsteinTorsion(th))
	fs := Frame{S: 0.3}
	p1 := Phi(fs.EinsteinTorsion(th))
	fmt.Printf("   Kron's frame (stationary rotor coils, a rotation): φ = (%.0e, %.0e, %.0e, %.0e) — a \"pure gravitational\" field to Einstein\n", p0[0], p0[1], p0[2], p0[3])
	fmt.Printf("   the same motor in a frame whose rotor legs also grow as e^{0.3θ}: φ_θ = %.3f\n", p1[3])
	fmt.Println("   In Kron's machine the frame is a choice of current paths; the motor is the same in both. So φ here is not")
	fmt.Println("   a field the machine has but a property of the frame chosen — Einstein's own worry that \"the separation of the")
	fmt.Println("   gravitational and electromagnetic field seems arbitrary\".")

	fmt.Println("\n3. Curvature")
	ce, ck := Curvatures(f, th)
	ce2, _ := Curvatures(fs, th)
	fmt.Printf("   Einstein's distant-parallelism connection: largest curvature component %.0e (and %.0e in the scaled frame) —\n", ce, ce2)
	fmt.Println("   flat, all its geometry is torsion.")
	fmt.Printf("   The machine's own connection (Christoffel symbols of its inductances and inertia): largest component %.3g\n", ck)
	fmt.Println("   (mixed units of henry and kg·m²; what matters is that it is not zero — the machine's space is curved).")
	fmt.Println("   Kron's torsion is that curved space seen from a moving frame (Hoffmann eq. 17); it cannot remove the curvature.")
	fmt.Println("   Einstein's torsion replaces curvature. Same mathematics of frames, different geometry.")

	fmt.Println("\n4. Einstein's first approximation (§2), plane waves e^{p·x} with x⁴ imaginary")
	for _, c := range []struct {
		name string
		p    [4]complex128
	}{
		{"p² ≠ 0", [4]complex128{1, 0.3, -0.2, 0.5}},
		{"p² ≠ 0, complex", [4]complex128{0.4, 0.2 + 0.3i, -0.5, 0.1 - 0.7i}},
		{"p² = 0 (light-like)", [4]complex128{1, 1i, 0, 0}},
	} {
		lc := CheckWave(c.p)
		fmt.Printf("   %-22s %2d solutions of eq. 5; largest φ %.2g; eq. 8 %.0e, eq. 9 %.0e, eq. 7 %.0e; R_αβ when φ = 0: %.0e\n",
			c.name, lc.Solutions, lc.MaxPhi, lc.Eq8, lc.Eq9, lc.Eq7, lc.Ricci)
	}
	fmt.Println("   Away from the light cone the only solutions are coordinate changes (φ = 0, R = 0). On it, φ ≠ 0 obeys")
	fmt.Println("   ∂φ = 0 and □φ = 0 — Maxwell's vacuum equations — and with φ = 0 the metric obeys R_αβ = 0. Einstein's")
	fmt.Println("   first-approximation claims hold.")

	fmt.Println("\n5. The motor in kronspice, built from Einstein's torsion alone (120 V, 60 Hz, pushed to 40 rad/s, 4 s)")
	s1, s2, s3, s4, err := Specs()
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("   speed matrix G = 2(Ω_θλν + Ω_νλθ) from the torsion: %v (kron.RotatingFrame: %v)\n", rnd(s3.G), rnd(s2.G))
	var runs []Run
	for _, c := range []struct {
		name string
		spec *circuit.MachineSpec
	}{{"frame I (coils turn with the rotor)", s1}, {"Kron frame II(b) (kron.RotatingFrame)", s2}, {"built from Einstein's torsion", s3}, {"same inductances, torsion removed", s4}} {
		r, err := Start(c.name, c.spec, 40, 4)
		if err != nil {
			log.Fatal(err)
		}
		runs = append(runs, r)
		fmt.Printf("   %-40s final speed %9.4f rad/s; ledger mismatch %.0e\n", c.name, r.Final, r.Mismatch)
	}
	fmt.Printf("   The torsion motor and Kron's frame II(b) agree to %.0e rad/s; frame I to %.2f rad/s (time step).\n",
		math.Abs(runs[2].Final-runs[1].Final), math.Abs(runs[2].Final-runs[0].Final))
	fmt.Printf("   Without the torsion it only coasts: 40·e^{−4 s·friction/J} = %.2f rad/s.\n", 40*math.Exp(-4*1e-4/J))

	fmt.Println("\nAnswer: Einstein's paper gives Kron's non-holonomic frames their geometric language — a 4-bein whose")
	fmt.Println("torsion carries the motion — and kronspice runs a motor from that torsion alone. It does not make Kron's")
	fmt.Println("machines electromagnetism in Einstein's sense: their torsion is trace-free (φ = 0) and sits on a curved space.")
	fmt.Println("Einstein's theory itself was dropped by 1931; its gravity survives as the teleparallel equivalent of GR.")

	if *viz == "" {
		return
	}
	g, err := os.Create(filepath.Join(*viz, "speed.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(g, "run,t,speed")
	for _, r := range runs {
		step := len(r.Time)/1500 + 1
		for k := 0; k < len(r.Time); k += step {
			fmt.Fprintf(g, "%q,%g,%g\n", r.Name, r.Time[k], r.Speed[k])
		}
	}
	g.Close()
	h, err := os.Create(filepath.Join(*viz, "torsion.csv"))
	if err != nil {
		log.Fatal(err)
	}
	fmt.Fprintln(h, "lambda,mu,nu,einstein_frame,kron_anholonomy,hoffmann_omega")
	for l := 0; l < 4; l++ {
		for m := 0; m < 4; m++ {
			for n := 0; n < 4; n++ {
				if A[l][m][n] != 0 || om[l][m][n] != 0 {
					fmt.Fprintf(h, "%d,%d,%d,%g,%g,%g\n", l, m, n, E[l][m][n], A[l][m][n], om[l][m][n])
				}
			}
		}
	}
	h.Close()
}

func rnd(m [][]float64) [][]float64 {
	o := make([][]float64, len(m))
	for i := range m {
		o[i] = make([]float64, len(m[i]))
		for j := range m[i] {
			o[i][j] = math.Round(m[i][j]*1e6) / 1e6
		}
	}
	return o
}
