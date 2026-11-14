// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command void_shells tests proposed hypothesis H50: a particle as a void
// in a nonlinear ether, its light pushing out and the ether pushing in,
// built up in shells (nested void and ether layers). It is the simplest
// structure left standing after examples/ether_horizon Part 13, which found
// that every donut-shaped void either travels or settles into a round bubble
// of light.
//
// The fields (field.go): a Gross–Pitaevskii ether ψ and the double-loop
// light (m = ½, spin ½), one photon per occupied mode, coupled both ways on
// an axisymmetric grid. No vortex, no added structure: the shapes come out
// of the balance.
//
// Parts:
//  1. the round bubble: one photon in the lowest mode, across the light's
//     strength ε;
//  2. one photon in a higher mode: does its node let ether in (a shell)?
//  3. a prepared shell: an ether core inside a void shell; relaxed, and in
//     real time;
//  4. several photons, one per mode (they are spin ½), filling up;
//  5. a prepared onion: two photons in two voids with an ether layer between;
//  6. g for every state, and why none reaches 2.
//
// Standard prediction beside it: the electron is pointlike below 10⁻¹⁸ m
// with g = 2.00231930436 from Dirac's equation and QED; there is no ether.
package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// The run's grid and constants (healing-length units; sound speed 1/√2).
const (
	GridH  = 0.4
	GridLR = 24.0
	GridLZ = 24.0
	V0     = 20.0 // the ether's opacity to light
	EpsRef = 300.0
	Relax  = 8000
	// AeMeasured is (g − 2)/2 (Fan et al. 2023).
	AeMeasured = 1.15965218059e-3
)

func grid() *Grid { return NewGrid(GridH, GridLR, GridLZ) }

// Result is one relaxed state.
type Result struct {
	Name    string
	Ee, EL  float64
	G0      float64
	L       Layout
	Lam     []float64
	RLbar   float64 // void reach at z = 0 over the reduced Compton length ƛ = ε/M
	EtherSh float64 // share of the mass in the ether
}

func relaxed(name string, eps float64, occ []int, amp func(r, z float64) float64, steps int) (Result, *State) {
	s := NewState(grid(), eps, V0, occ, amp)
	s.Relax(steps)
	ee, el := s.Energies()
	l := s.Measure()
	m := ee + el
	return Result{Name: name, Ee: ee, EL: el, G0: m / el, L: l, Lam: s.Lam, RLbar: l.REq * m / eps, EtherSh: ee / m}, s
}

func main() {
	viz := flag.String("viz", "", "write CSVs and animation frames to DIR")
	flag.Parse()
	var csvRows []string
	add := func(part string, r Result) {
		csvRows = append(csvRows, fmt.Sprintf("%s,%q,%.6f,%.6f,%.8f,%d,%.4f,%.3f,%.3f,%.6f,%.4f", part, r.Name, r.Ee, r.EL, r.G0, r.L.Voids, r.L.CoreDens, r.L.REq, r.L.RAxis, r.L.LightIn, r.RLbar))
	}

	fmt.Println("void_shells (proposed H50): light pushing out, ether pushing in, built up in shells")
	fmt.Printf("Grid %.1f ξ, box ρ < %.0f, |z| < %.0f; ether opacity V0 = %.0f; relaxed %d steps (energy only falls).\n", GridH, GridLR, GridLZ, V0, Relax)
	fmt.Println("g = Mc²/E_L (rest energy over the light's, the charge riding with the light; ether_horizon Part 11).")
	fmt.Println("Standard prediction: a pointlike electron, g = 2.00231930436 from Dirac + QED, no ether.")

	fmt.Println()
	fmt.Println("Part 1. The round bubble: one double-loop photon in the lowest mode, across its strength ε (ε plays ħc)")
	fmt.Println("  ε       void ρ / z reach   light in void   ether share   g = M/E_L   size: ρ reach / ƛ")
	for _, eps := range []float64{30, 100, 300, 1000, 3000} {
		r, _ := relaxed(fmt.Sprintf("ball ε=%g", eps), eps, []int{0}, Ball(2+math.Cbrt(eps)/2), Relax)
		add("1", r)
		fmt.Printf("  %-6.0f  %4.1f / %-11.1f  %-14.4f  %-12.4f  %-10.4f  %.2f\n", eps, r.L.REq, r.L.RAxis, r.L.LightIn, r.EtherSh, r.G0, r.RLbar)
	}
	fmt.Println("  The light carries about three quarters of the mass and the ether the rest, so g sits between 4/3 and 3/2")
	fmt.Println("  for every strength (Part 6). The bubble is a few reduced Compton lengths across.")

	fmt.Println()
	fmt.Println("Part 2. One photon in a higher mode (ε = 300): does a node let the ether in and make a shell?")
	fmt.Println("  mode   voids   void ρ / z reach   centre density   g")
	for _, k := range []int{1, 2, 3, 4, 5, 6} {
		r, _ := relaxed(fmt.Sprintf("mode %d", k), EpsRef, []int{k - 1}, Ball(4+1.5*float64(k)), 2*Relax)
		add("2", r)
		fmt.Printf("  %-5d  %-6d  %4.1f / %-11.1f  %-15.2e  %.4f\n", k, r.L.Voids, r.L.REq, r.L.RAxis, r.L.CoreDens, r.G0)
	}
	fmt.Println("  No. Every excited mode keeps one void: the ether cannot enter a node thinner than its healing length,")
	fmt.Println("  so the void just stretches along the axis into a chain of lobes. g goes down, not up.")

	fmt.Println()
	fmt.Println("Part 3. A prepared shell: an ether core (r < 3) inside a void shell (3 < r < 8), one photon (ε = 300)")
	sh, _ := relaxed("prepared shell", EpsRef, []int{0}, ShellAmp(3, 8), Relax)
	ball, _ := relaxed("ball", EpsRef, []int{0}, Ball(5), Relax)
	add("3", sh)
	fmt.Printf("  relaxed:  voids %d, centre density %.1e, M = %.6f, g = %.6f\n", sh.L.Voids, sh.L.CoreDens, sh.Ee+sh.EL, sh.G0)
	fmt.Printf("  the plain round bubble: M = %.6f, g = %.6f — the same state: the core dissolved.\n", ball.Ee+ball.EL, ball.G0)
	fmt.Println("  In real time (no friction), the same start:")
	var frameDir string
	if *viz != "" {
		frameDir = filepath.Join(*viz, "frames")
		if err := os.MkdirAll(frameDir, 0o755); err != nil {
			log.Fatal(err)
		}
	}
	rt := realTime("shell", ShellAmp(3, 8), []int{0}, 40, frameDir)
	for _, row := range rt {
		fmt.Println(row)
	}
	fmt.Println("  The light squeezes the core out as sound: by t ≈ 10 the centre is empty and a breathing round bubble is left.")

	fmt.Println()
	fmt.Println("Part 4. Several photons, one per mode (spin ½, so no two share a mode), ε = 300 each")
	fmt.Println("  N   voids   void ρ / z reach   g        mass per photon")
	for n := 1; n <= 6; n++ {
		occ := make([]int, n)
		for i := range occ {
			occ[i] = i
		}
		r, _ := relaxed(fmt.Sprintf("N=%d", n), EpsRef, occ, Ball(4+1.2*float64(n)), 2*Relax)
		add("4", r)
		fmt.Printf("  %-2d  %-6d  %4.1f / %-11.1f  %-7.4f  %.2f\n", n, r.L.Voids, r.L.REq, r.L.RAxis, r.G0, (r.Ee+r.EL)/float64(n))
	}
	fmt.Println("  They share one growing void; no ether layer forms between them. g edges toward 4/3.")

	fmt.Println()
	fmt.Println("Part 5. A prepared onion: void core r < 3, ether layer 3 … 5.5, void shell 5.5 … 9, two photons")
	on, _ := relaxed("onion", EpsRef, []int{0, 1}, OnionAmp(3, 5.5, 9), 2*Relax)
	two, _ := relaxed("N=2 from a ball", EpsRef, []int{0, 1}, Ball(6.4), 2*Relax)
	add("5", on)
	fmt.Printf("  relaxed:  voids %d, M = %.4f, g = %.6f;  two photons from a ball: M = %.4f, g = %.6f\n", on.L.Voids, on.Ee+on.EL, on.G0, two.Ee+two.EL, two.G0)
	fmt.Println("  The ether layer is squeezed out from both sides and the onion becomes the two-photon bubble.")
	fmt.Println("  In real time the round bubble breathes but holds (energy conserved below; ε = 300):")
	for _, row := range realTime("ball", Ball(5), []int{0}, 40, frameDir) {
		fmt.Println(row)
	}

	fmt.Println()
	fmt.Println("Part 6. Why no shell reaches g = 2")
	fmt.Println("  g = Mc²/E_L = 1 + E_ether/E_L. Scaling the void by λ, the light's energy goes as 1/λ, the ether's")
	fmt.Println("  pressure work as λ³ and its wall as λ². At balance E_L = 3E_P + 2E_W, so")
	fmt.Println("  g = 1 + (E_P + E_W)/(3E_P + 2E_W): between 4/3 (pressure) and 3/2 (wall). The round bubbles lie there; the")
	fmt.Println("  excited modes dip just below 4/3 (1.31–1.32) because the ether's opacity is finite and the light leaks into")
	fmt.Println("  the wall, so its energy is not purely 1/λ. Every state computed here has g = 1.31–1.37, about 32–35 % below 2.")
	fmt.Println("  g = 2 needs an ether energy that grows as λ¹ — a constant-force string. In this ether only a vortex")
	fmt.Println("  line does that, and ether_horizon Part 13 found vortex rings travel or unwind. Shells only add walls")
	fmt.Println("  (λ²) and layers (λ³), so they cannot move g toward 2; the light squeezes them out anyway.")

	if *viz != "" {
		f, err := os.Create(filepath.Join(*viz, "states.csv"))
		if err != nil {
			log.Fatal(err)
		}
		fmt.Fprintln(f, "part,name,Ee,EL,g,voids,core_density,rho_reach,z_reach,light_in,rho_reach_over_lbar")
		for _, r := range csvRows {
			fmt.Fprintln(f, r)
		}
		f.Close()
		for _, c := range []struct {
			name string
			occ  []int
			amp  func(r, z float64) float64
		}{
			{"ball", []int{0}, Ball(5)}, {"mode3", []int{2}, Ball(8.5)}, {"mode6", []int{5}, Ball(13)},
			{"N4", []int{0, 1, 2, 3}, Ball(8.8)}, {"shell_start", []int{0}, ShellAmp(3, 8)}, {"onion_start", []int{0, 1}, OnionAmp(3, 5.5, 9)},
		} {
			steps := 2 * Relax
			if c.name == "shell_start" || c.name == "onion_start" {
				steps = 0
			}
			s := NewState(grid(), EpsRef, V0, c.occ, c.amp)
			if steps > 0 {
				s.Relax(steps)
			}
			writeSnapshot(filepath.Join(frameDir, "state_"+c.name+".bin"), s)
		}
	}
}

// realTime evolves a start without friction and returns printed rows; with
// frameDir set it writes frames every 0.5 time units.
func realTime(tag string, amp func(r, z float64) float64, occ []int, T float64, frameDir string) []string {
	s := NewState(grid(), EpsRef, V0, occ, amp)
	e0, l0 := s.Energies()
	dt := 0.01
	n := int(math.Round(T / dt))
	var rows []string
	var fr *os.File
	if frameDir != "" {
		var err error
		fr, err = os.Create(filepath.Join(frameDir, "rt_"+tag+".bin"))
		if err != nil {
			log.Fatal(err)
		}
		defer fr.Close()
		binary.Write(fr, binary.LittleEndian, []int32{int32(s.G.NR), int32(s.G.NZ)})
	}
	for k := 0; k <= n; k++ {
		t := float64(k) * dt
		if k%50 == 0 && fr != nil {
			writeFrame(fr, t, s)
		}
		if k%500 == 0 {
			ee, el := s.Energies()
			l := s.Measure()
			rows = append(rows, fmt.Sprintf("    t = %-3.0f centre density %.3f  voids %d  void ρ / z reach %4.1f / %-4.1f  light in void %.4f  energy drift %+.5f",
				t, l.CoreDens, l.Voids, l.REq, l.RAxis, l.LightIn, (ee+el-e0-l0)/(e0+l0)))
		}
		if k < n {
			s.Step(dt)
		}
	}
	return rows
}

func writeFrame(f *os.File, t float64, s *State) {
	binary.Write(f, binary.LittleEndian, float32(t))
	d := make([]float32, len(s.dens))
	x := make([]float32, len(s.dens))
	for i := range d {
		d[i] = float32(s.dens[i])
		for _, q := range s.Occ {
			x[i] += float32(s.Chi[q][i] * s.Chi[q][i])
		}
	}
	binary.Write(f, binary.LittleEndian, d)
	binary.Write(f, binary.LittleEndian, x)
}

func writeSnapshot(path string, s *State) {
	f, err := os.Create(path)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	s.updateDens()
	binary.Write(f, binary.LittleEndian, []int32{int32(s.G.NR), int32(s.G.NZ)})
	writeFrame(f, 0, s)
}
