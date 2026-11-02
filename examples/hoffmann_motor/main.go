// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command hoffmann_motor works through B. Hoffmann, "Kron's Non-Riemannian
// Electrodynamics," Rev. Mod. Phys. 21, 535–540 (1949): the single-phase
// squirrel-cage induction motor written two ways.
//
//  1. Frame I (holonomic, Riemannian): the cage as two rotor coils that turn
//     with the rotor; inductances depend on θ (eqs. 24, 42), Lagrange–
//     Maxwell equations with Christoffel symbols (eq. 5).
//  2. Stationary coils 2̄, 3̄ (non-holonomic, eq. 35): constant inductances
//     (eq. 29), but the space now has torsion Ω (eq. 19), and the equations
//     gain the term 2Ω·ẋ·ẋ (eq. 5′). Without it — ordinary Lagrange
//     equations on the constant inductances — the results are wrong (p. 540).
//
// Each is solved three ways — kronspice's machine element in frame I, the
// same element in Kron's frame II(b) (kron.RotatingFrame), and eq. 5′
// integrated directly (eq5.go) — from rest, pushed forward and backward,
// held at fixed speeds for the torque–speed curve, and with a salient pole.
//
//	go run ./examples/hoffmann_motor [-viz dir]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/kron"
)

func main() {
	viz := flag.String("viz", "", "directory for picture data")
	flag.Parse()
	must := func(err error) {
		if err != nil {
			log.Fatal(err)
		}
	}
	p := Default()
	s1 := p.FrameI()
	s2, err := p.FrameIIb()
	must(err)

	fmt.Println("Hoffmann's single-phase squirrel-cage motor (Rev. Mod. Phys. 21, 535, 1949); values declared:")
	fmt.Printf("   stator L11 %.2f H, rotor coils L_S %.2f H, mutual M %.2f H, salient pole L_D %.2f H; Rs %.0f Ω, rotor coils %.0f Ω;\n", p.L11, p.LS, p.M, p.LD, p.Rs, p.Rr)
	fmt.Printf("   J %.0e kg·m², friction %.0e N·m·s; 120 V, 60 Hz; two poles, synchronous speed %.1f rad/s\n", p.J, p.Friction, p.Sync())

	fmt.Println("\n1. The two frames")
	lb, _, err := kron.NewMetric(s1, p.J, Frame(), 0)
	must(err)
	fmt.Println("   Stationary rotor coils (eq. 35): the inductances no longer depend on θ (eq. 29):")
	names := []string{"1 ", "2̄ ", "3̄ ", "θ "}
	for i, row := range lb {
		fmt.Printf("      %s", names[i])
		for _, v := range row {
			fmt.Printf(" %8.4f", v)
		}
		fmt.Println()
	}
	om, err := kron.Torsion(s1, p.J, Frame(), 0)
	must(err)
	fmt.Print("   but the frame has torsion Ω_λμν (eq. 19), non-zero components:")
	for a := range om {
		for b := range om[a] {
			for c := range om[a][b] {
				if a < b && math.Abs(om[a][b][c]) > 1e-12 {
					fmt.Printf(" Ω(%s%s%s) = %+.3f", trim(names[a]), trim(names[b]), trim(names[c]), om[a][b][c])
				}
			}
		}
	}
	fmt.Printf("\n   (and Ω(μλν) = −Ω(λμν)); it does not change with θ (spread %.0e).\n", torsionSpread(p))
	worst, wrong := twin(p, 200)
	fmt.Printf("   At 200 random states, eq. 5′ with the torsion gives the frame-I forces to %.0e;\n", worst)
	fmt.Printf("   ordinary Lagrange–Maxwell equations on the constant inductances miss by up to %.0f %% of the largest force.\n", 100*wrong)

	fmt.Println("\n2. Starting (5 s runs)")
	var runs []Run
	for _, c := range []struct {
		name string
		mode int // 0 frame I, 1 frame II(b), 2 eq. 5′, 3 naive
		w0   float64
	}{
		{"from rest, frame I", 0, 0},
		{"pushed to +40 rad/s, frame I", 0, 40},
		{"pushed to +40 rad/s, frame II(b)", 1, 40},
		{"pushed to +40 rad/s, eq. 5′ (Runge–Kutta)", 2, 40},
		{"pushed to −40 rad/s, frame I", 0, -40},
		{"pushed to +40 rad/s, no torsion (Hoffmann's warning)", 3, 40},
	} {
		var r Run
		switch c.mode {
		case 0:
			r, err = p.Start(c.name, s1, c.w0, 5)
		case 1:
			r, err = p.Start(c.name, s2, c.w0, 5)
		case 2:
			r, err = p.Eq5Run(c.w0, 5)
		case 3:
			naive, e := p.Naive()
			must(e)
			r, err = p.Start(c.name, naive, c.w0, 5)
		}
		must(err)
		r.Name = c.name
		runs = append(runs, r)
		ledger := ""
		if c.mode != 2 {
			ledger = fmt.Sprintf("; mains %.0f J = heat %.0f J + stored %.1f J (mismatch %.0e)", r.Supplied, r.Heat, r.Stored, r.Mismatch)
		}
		fmt.Printf("   %-52s final %8.3f rad/s%s\n", c.name, r.Final, ledger)
	}
	fmt.Println("   A single-phase motor has no starting torque: its pulsating field is two equal fields turning opposite ways.")
	fmt.Println("   Pushed either way it runs up that way to just under synchronous speed. The three solutions agree")
	fmt.Printf("   (frame I vs eq. 5′: %.1e rad/s). Without the torsion term there is no speed voltage and no torque:\n", math.Abs(runs[1].Final-runs[3].Final))
	fmt.Printf("   the pushed rotor only coasts down on friction (40 · e^(−5 s · friction/J) = %.3f rad/s).\n", 40*math.Exp(-5*p.Friction/p.J))

	fmt.Println("\n3. Torque–speed curve (shaft held at each speed; mean torque over the last 10 cycles)")
	var speeds []float64
	for k := -12; k <= 12; k++ {
		speeds = append(speeds, float64(k)/12*p.Sync()*0.995)
	}
	sal := p
	sal.LD = 0.02
	sals := sal.FrameI()
	var tq, tqs, tq2 []float64
	for _, w := range speeds {
		a, err := p.Torque(s1, w)
		must(err)
		b, err := sal.Torque(sals, w)
		must(err)
		c, err := p.Torque(s2, w)
		must(err)
		tq, tqs, tq2 = append(tq, a), append(tqs, b), append(tq2, c)
	}
	dev := 0.0
	top := 0.0
	for k := range tq {
		dev = math.Max(dev, math.Abs(tq[k]-tq2[k]))
		top = math.Max(top, math.Abs(tq[k]))
	}
	fmt.Printf("   %10s %14s %18s\n", "rad/s", "torque N·m", "salient L_D 0.02")
	for k := 0; k < len(speeds); k += 2 {
		fmt.Printf("   %10.1f %14.4f %18.4f\n", speeds[k], tq[k], tqs[k])
	}
	fmt.Printf("   Zero at standstill, odd in speed, largest near %0.f %% of synchronous speed. Frames I and II(b) agree to %.1e of the peak.\n",
		100*math.Abs(speeds[argmax(tq)])/p.Sync(), dev/top)
	fmt.Println("   The salient pole (L_D = 0.02 H) makes the stationary coil 2̄ couple more tightly to the stator (L̄22 = L_S − L_D")
	fmt.Println("   = M), which roughly doubles the torque at most speeds. With no load it runs a little closer to synchronous speed:")
	for _, ld := range []float64{0, 0.01, 0.02} {
		q := p
		q.LD = ld
		r, err := q.Start("salient", q.FrameI(), 40, 5)
		must(err)
		fmt.Printf("      L_D = %.2f H: final %.3f rad/s (slip %.3f %%)\n", ld, r.Final, 100*(1-r.Final/p.Sync()))
	}
	fmt.Println("   It does not pull into synchronism at these values (a reluctance motor needs a much larger L_D).")
	fmt.Println("\nThe motor adds no energy: every joule of shaft motion and heat comes from the 120 V supply (or the push).")

	if *viz == "" {
		return
	}
	f, err := os.Create(filepath.Join(*viz, "speed.csv"))
	must(err)
	fmt.Fprintln(f, "run,t,speed")
	for _, r := range runs {
		step := len(r.Time)/1500 + 1
		for k := 0; k < len(r.Time); k += step {
			fmt.Fprintf(f, "%q,%g,%g\n", r.Name, r.Time[k], r.Speed[k])
		}
	}
	f.Close()
	g, err := os.Create(filepath.Join(*viz, "torque.csv"))
	must(err)
	fmt.Fprintln(g, "speed,torque,torque_frame2,torque_salient")
	for k := range speeds {
		fmt.Fprintf(g, "%g,%g,%g,%g\n", speeds[k], tq[k], tq2[k], tqs[k])
	}
	g.Close()
}

func trim(s string) string { return s[:len(s)-1] }

func argmax(x []float64) int {
	k := 0
	for i := range x {
		if x[i] > x[k] {
			k = i
		}
	}
	return k
}

// twin compares eq. 5′ with the frame-I equations at random states and
// returns the largest relative difference with the torsion term and
// without it.
func twin(p Params, n int) (worst, wrong float64) {
	s, f := p.FrameI(), Frame()
	rng := rand.New(rand.NewSource(1949))
	for k := 0; k < n; k++ {
		th := 2 * math.Pi * rng.Float64()
		ip := []float64{rng.NormFloat64(), rng.NormFloat64(), rng.NormFloat64()}
		dip := []float64{100 * rng.NormFloat64(), 100 * rng.NormFloat64(), 100 * rng.NormFloat64()}
		om, dom := 300*rng.NormFloat64(), 1000*rng.NormFloat64()
		c, d := f.C(th), f.D(th)
		i, di := make([]float64, 3), make([]float64, 3)
		for a := 0; a < 3; a++ {
			for b := 0; b < 3; b++ {
				i[a] += c[a][b] * ip[b]
				di[a] += c[a][b]*dip[b] + d[a][b]*ip[b]*om
			}
		}
		old := kron.Geodesic(s, p.J, th, i, di, om, dom)
		want := make([]float64, 4)
		for la := 0; la < 3; la++ {
			for a := 0; a < 3; a++ {
				want[la] += c[a][la] * old[a]
			}
		}
		want[3] = old[3]
		got, _ := kron.NonRiemannian(s, p.J, f, th, ip, dip, om, dom, true)
		bad, _ := kron.NonRiemannian(s, p.J, f, th, ip, dip, om, dom, false)
		scale := 0.0
		for _, v := range want {
			scale = math.Max(scale, math.Abs(v))
		}
		for la := range want {
			worst = math.Max(worst, math.Abs(got[la]-want[la])/scale)
			wrong = math.Max(wrong, math.Abs(bad[la]-want[la])/scale)
		}
	}
	return worst, wrong
}
