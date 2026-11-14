// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
)

// Part13Grid is the grid of the Part 13 runs (healing-length units).
func Part13Grid() *EtherGrid { return NewEtherGrid(0.4, 20, 26) }

// Part13V0 is the ether's opacity to light (V0 > λ everywhere below, so
// the light is bound in the void and cannot travel through dense ether).
const Part13V0 = 20

// RunCase is one real-time run.
type RunCase struct {
	Tag    string
	Eps    float64
	R0, A0 float64
	T, DT  float64
	Frames bool
}

// Part13Cases: the thin donut at three light strengths, a half-fat donut
// and the horn torus of Part 12.
var Part13Cases = []RunCase{
	{"thin, ε = 100", 100, 8, 2, 60, 0.01, false},
	{"thin, ε = 300", 300, 8, 2, 60, 0.01, true},
	{"thin, ε = 1000", 1000, 8, 2, 60, 0.01, true},
	{"half, ε = 300", 300, 8, 4, 60, 0.01, false},
	{"horn, ε = 300", 300, 8, 8, 60, 0.01, true},
}

// RunSummary is what a real-time run reports.
type RunSummary struct {
	Speed       float64 // dz_c/dt of the void, fitted over t ∈ [5, T]
	Drift       float64 // relative change of E_e + E_L by the end
	MinLightIn  float64
	EndShape    VoidShape
	AxisShare   float64 // share of samples (t ≥ 10) with the void on the axis
	VortexShare float64 // share of samples (t ≥ 10) with winding ±1 round the void
	AnomalyWall float64 // mean α⟨1/d⟩/(2√λ) over t ≥ 10 (charge on the computed wall)
	MeanAoverR  float64 // mean (outer − inner)/(outer + inner) over samples off the axis
	Rows        []string
	CSV         []string // t, inner, outer, height, on_axis, winding, light_in, z_c, drift, anomaly
}

// RunRealTime evolves one case and samples it every unit of time.
func RunRealTime(c RunCase, frameDir string) RunSummary {
	g := Part13Grid()
	s := NewRingState(g, c.Eps, Part13V0, c.R0, c.A0)
	ee0, el0, _ := s.Energies()
	var sum RunSummary
	sum.MinLightIn = 1
	var ts, zs []float64
	var an, ar []float64
	nAxis, nVort, nLate := 0, 0, 0
	every := int(math.Round(1 / c.DT))
	n := int(math.Round(c.T / c.DT))
	var fr *os.File
	if c.Frames && frameDir != "" {
		var err error
		fr, err = os.Create(filepath.Join(frameDir, frameName(c)))
		if err != nil {
			log.Fatal(err)
		}
		defer fr.Close()
		binary.Write(fr, binary.LittleEndian, []int32{int32(g.NR), int32(g.NZ)})
	}
	for k := 0; k <= n; k++ {
		if k%every == 0 {
			t := float64(k) * c.DT
			sh := s.Shape()
			ee, el, _ := s.Energies()
			if sh.Any {
				ts = append(ts, t)
				zs = append(zs, sh.ZC)
				sum.MinLightIn = math.Min(sum.MinLightIn, sh.LightIn)
				if t >= 10 {
					nLate++
					if sh.OnAxis {
						nAxis++
					} else {
						ar = append(ar, (sh.Outer-sh.Inner)/(sh.Outer+sh.Inner))
					}
					if math.Abs(math.Abs(sh.Winding)-1) < 0.01 {
						nVort++
					}
					an = append(an, Alpha*s.WallInvDGrid()/(2*math.Sqrt(s.Lam))) // the light's void, set by Shape above
				}
			}
			if k%(10*every) == 0 {
				sum.Rows = append(sum.Rows, fmt.Sprintf("    t = %-4.0f void ρ %5.1f–%-5.1f height %5.1f  on axis %-5v winding %4.1f  light in void %.4f  z_c %7.2f  energy drift %+.4f",
					t, sh.Inner, sh.Outer, sh.Height, sh.OnAxis, sh.Winding, sh.LightIn, sh.ZC, (ee+el-ee0-el0)/(ee0+el0)))
			}
			if fr != nil {
				binary.Write(fr, binary.LittleEndian, float32(t))
				d := make([]float32, len(s.dens))
				x := make([]float32, len(s.Chi))
				for i := range d {
					d[i] = float32(s.dens[i])
					x[i] = float32(s.Chi[i] * s.Chi[i])
				}
				binary.Write(fr, binary.LittleEndian, d)
				binary.Write(fr, binary.LittleEndian, x)
			}
			anom := math.NaN()
			if sh.Any {
				anom = Alpha * s.WallInvDGrid() / (2 * math.Sqrt(s.Lam))
			}
			sum.CSV = append(sum.CSV, fmt.Sprintf("%.2f,%.3f,%.3f,%.3f,%v,%.2f,%.6f,%.4f,%.6f,%.6e", t, sh.Inner, sh.Outer, sh.Height, sh.OnAxis, sh.Winding, sh.LightIn, sh.ZC, (ee+el-ee0-el0)/(ee0+el0), anom))
			sum.EndShape = sh
			sum.Drift = (ee + el - ee0 - el0) / (ee0 + el0)
		}
		if k < n {
			s.Step(c.DT)
		}
	}
	sum.Speed = slope(ts, zs, 5)
	sum.AxisShare = float64(nAxis) / math.Max(1, float64(nLate))
	sum.VortexShare = float64(nVort) / math.Max(1, float64(nLate))
	sum.AnomalyWall = mean(an)
	sum.MeanAoverR = mean(ar)
	return sum
}

func frameName(c RunCase) string {
	switch {
	case c.A0 >= c.R0:
		return "frames_horn.bin"
	default:
		return fmt.Sprintf("frames_thin_e%.0f.bin", c.Eps)
	}
}

func slope(t, y []float64, from float64) float64 {
	var n, st, sy, stt, sty float64
	for i := range t {
		if t[i] < from {
			continue
		}
		n++
		st += t[i]
		sy += y[i]
		stt += t[i] * t[i]
		sty += t[i] * y[i]
	}
	return (n*sty - st*sy) / (n*stt - st*st)
}

func mean(x []float64) float64 {
	if len(x) == 0 {
		return math.NaN()
	}
	s := 0.0
	for _, v := range x {
		s += v
	}
	return s / float64(len(x))
}

// RelaxResult is the end of a dissipative run.
type RelaxResult struct {
	Eps                 float64
	Shape               VoidShape
	G0, Lam, EnergyFall float64
}

// RunRelax lets the energy drain (imaginary time) from the thin donut.
func RunRelax(eps float64, steps int) RelaxResult {
	g := Part13Grid()
	s := NewRingState(g, eps, Part13V0, 8, 2)
	e0, l0, _ := s.Energies()
	s.Relax(steps, 20)
	ee, el, _ := s.Energies()
	return RelaxResult{Eps: eps, Shape: s.Shape(), G0: (ee + el) / el, Lam: s.Lam, EnergyFall: (ee + el) / (e0 + l0)}
}

func part13(vizDir string) {
	fmt.Println()
	fmt.Println("Part 13. The full ether simulation: Gross–Pitaevskii ether + the m = ½ light mode on an axisymmetric grid")
	fmt.Println("  (ρ, z) grid 0.4 ξ, box 20 × 52 ξ; ether opaque to light where dense (V0 = 20); no thin-ring formulas.")
	fmt.Println("  Light re-solved every step (inverse iteration, conjugate gradients); ether by RK4 in real time.")

	fmt.Println()
	fmt.Println("Part 13a. With friction (energy drains, as Part 12 assumed): start from the thin donut R = 8, a = 2")
	fmt.Println("  ε       end: void ρ        on axis   winding   light in void   g0 = M/E_L   energy kept")
	for _, eps := range []float64{30, 100, 300} {
		r := RunRelax(eps, 30000)
		fmt.Printf("  %-6.0f  %4.1f–%-10.1f  %-8v  %-8.1f  %-14.4f  %-11.4f  %.3f\n", eps, r.Shape.Inner, r.Shape.Outer, r.Shape.OnAxis, r.Shape.Winding, r.Shape.LightIn, r.G0, r.EnergyFall)
	}
	fmt.Println("  Every run loses the vortex: the ring slides in through the void and unwinds on the axis. What is left")
	fmt.Println("  is a round bubble of light with no swirl, held by the ether's pressure (g0 ≈ 1.4, near Part 9's 4/3).")

	fmt.Println()
	fmt.Println("Part 13b. Without friction (real time, energy conserved): five starts, 60 time units each")
	var frameDir string
	if vizDir != "" {
		frameDir = filepath.Join(vizDir, "part13_frames")
		if err := os.MkdirAll(frameDir, 0o755); err != nil {
			log.Fatal(err)
		}
	}
	var csv []string
	csv = append(csv, "case,t,inner,outer,height,on_axis,winding,light_in,z_c,drift,anomaly")
	for _, c := range Part13Cases {
		r := RunRealTime(c, frameDir)
		for _, row := range r.CSV {
			csv = append(csv, fmt.Sprintf("%q,%s", c.Tag, row))
		}
		fmt.Printf("  %s (start R = %.0f, a = %.0f):\n", c.Tag, c.R0, c.A0)
		for _, row := range r.Rows {
			fmt.Println(row)
		}
		fmt.Printf("    speed along the axis %.3f (sound 0.707); light kept ≥ %.4f; void on the axis %.0f%% of t ≥ 10;\n",
			r.Speed, r.MinLightIn, 100*r.AxisShare)
		fmt.Printf("    vortex round the void %.0f%% of t ≥ 10; mean a/R off axis %.3f; (g − 2)/2 from the wall charge %.3e (%.2f of measured)\n",
			100*r.VortexShare, r.MeanAoverR, r.AnomalyWall, r.AnomalyWall/AeMeasured)
	}
	if vizDir != "" {
		f, err := os.Create(filepath.Join(vizDir, "part13_runs.csv"))
		if err != nil {
			log.Fatal(err)
		}
		for _, row := range csv {
			fmt.Fprintln(f, row)
		}
		f.Close()
	}
	fmt.Println()
	fmt.Println("  What it shows:")
	fmt.Println("  - The horn torus is not a resting state. Started there, it breaks up into sound and ends as a small round")
	fmt.Println("    bubble of light on the axis; the vortex leaves it.")
	fmt.Println("  - With weak or moderate light (ε = 100, 300) the thin donut keeps its shape (a/R ≈ 0.24–0.32), its vortex and")
	fmt.Println("    its light, but it never swells toward the edge.")
	fmt.Println("  - With strong light (ε = 1000) the core does swell to the edge, as Part 12 said, but does not stop: the hole")
	fmt.Println("    closes, the void becomes a round bubble on the axis, and the vortex trails off outside the light.")
	fmt.Println("  - Part 12's stopping argument fails: the m = ½ light vanishes on the axis by itself (χ ∝ √ρ), so it needs")
	fmt.Println("    no hole in the ether to keep winding.")
	fmt.Println("  - The anomaly from the computed walls is 0.28–0.67 of measured, not Part 12's 0.959.")
	fmt.Println("  - A vortex ring always travels (here ≈ a third of the ether's sound speed); an electron at rest cannot")
	fmt.Println("    be a single vortex ring.")
}
