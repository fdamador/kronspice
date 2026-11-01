// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command kron_nmr solves magnetic resonance (NMR for the proton, ESR for
// the electron) as one of Kron's machines. Bloch's equations
//
//	dM/dt = γ·M × B − (Mx, My, 0)/T2 − (0, 0, Mz − M0)/T1
//
// are a rotation of the magnetization vector about the field, and a
// rotation of currents at a fixed rate is Kron's frame II(b) (Equivalent
// Circuits of Electric Machinery, 1951, Epilogue eq. 9, p. 257): windings
// with constant inductance and a skew-symmetric speed matrix G, so that
// L·di/dt = −ω·G·i. The three components Mx, My, Mz are the currents of
// three windings of unit inductance, G = L·γ·[B]ₓ (the matrix of M × B),
// and the mechanical port is turned at ω = 1. The torque iᵀ·G·i of a skew
// G is zero, so the precession takes no power from the shaft, as it must.
//
// The NMR rotating frame is Kron's rotating reference frame (§3b, Park's
// transform): in the frame turning at the RF frequency ω_rf the RF field
// B₁ is constant along x and the static field becomes B₀ − ω_rf/γ, so the
// whole problem has a constant G and resonance is where the effective
// z-field vanishes. Relaxation is resistance on the windings: T2 = L/R₂
// on x and y, T1 = L/R₁ on z with a source R₁·M₀ in series that restores
// the equilibrium magnetization; the energy the lattice receives is on
// the ledger as heat. A fourth, open winding coupled to x is the pickup
// coil, whose voltage M_p·dMx/dt is the signal a spectrometer sees.
//
//  1. Free precession in the laboratory frame: frequency against γB₀, T2
//     decay of the transverse magnetization, T1 recovery of Mz.
//  2. Rabi nutation at resonance in the rotating frame: Mz = M₀ cos(γB₁t).
//  3. A π/2 pulse then the free-induction decay, with the pickup voltage.
//  4. Continuous-wave line: steady-state absorption and dispersion against
//     the detuning, compared with Bloch's steady-state solution.
//  5. The resonance frequencies for the proton and the electron, measured,
//     and what the H21 ring (g = 1) would give instead.
//
// Animation (-gif): the magnetization vector in the rotating frame through
// the π/2 pulse and the decay, with the pickup signal drawn as it arrives.
//
//	go run ./examples/kron_nmr [-gif viz/nmr.gif]
package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/gif"
	"log"
	"math"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// Gyromagnetic ratios [rad/s/T] (CODATA 2018): proton γ_p = 2.675 221 87e8
// (42.577 MHz/T), electron |γ_e| = 1.760 859 63e11 (28.025 GHz/T).
const (
	GammaProton   = 2.6752218744e8
	GammaElectron = 1.76085963023e11
	QE            = 1.602176634e-19
	ME            = 9.1093837015e-31
	MP            = 1.67262192369e-27
)

// Spin is one magnetic-resonance run: gyromagnetic ratio Gamma, static
// field B0 along z, RF field B1 along x of the rotating frame turning at
// Wrf (0 = laboratory frame), relaxation times, equilibrium M0 and the
// initial magnetization. L is the winding inductance (1 H: currents are
// the magnetization in its own units); Mp the pickup mutual inductance.
type Spin struct {
	Gamma, B0, B1, Wrf float64
	T1, T2             float64
	M0                 float64
	Minit              [3]float64
	L, Mp              float64
}

// Spec is the machine: windings x, y, z, pickup; G = L·γ·[B_eff]ₓ with
// B_eff = (B1, 0, B0 − ω_rf/γ); unit shaft speed.
func (s Spin) Spec() *circuit.MachineSpec {
	m := circuit.NewMachineSpec(4)
	for k := 0; k < 4; k++ {
		m.L0[k][k] = s.L
	}
	m.L0[0][3], m.L0[3][0] = s.Mp, s.Mp
	bz := s.B0 - s.Wrf/s.Gamma
	bx := s.B1
	// L·di/dt = −G·i must equal L·γ·(M × B): (M × B)_x = My·Bz − Mz·By,
	// (M × B)_y = Mz·Bx − Mx·Bz, (M × B)_z = Mx·By − My·Bx, with By = 0.
	g := s.L * s.Gamma
	m.G[0][1] = -g * bz
	m.G[1][0] = g * bz
	m.G[1][2] = -g * bx
	m.G[2][1] = g * bx
	copy(m.IC, []float64{s.Minit[0], s.Minit[1], s.Minit[2], 0})
	return m
}

// Circuit: the machine with its windings closed through the relaxation
// resistances, the z winding's series source restoring M0, the pickup
// open through 1 MΩ, the shaft turned at 1 rad/s.
func (s Spin) Circuit() (*circuit.Circuit, error) {
	r2, r1 := s.L/s.T2, s.L/s.T1
	return circuit.NewBuilder("magnetic resonance as Kron's frame II(b)").
		Domain(exchange.Mechanical, "w").
		MachineModel("bloch", s.Spec()).
		Machine("N1", "bloch", [][2]string{{"x", "0"}, {"y", "0"}, {"z", "0"}, {"p", "0"}}, "w", "0").
		Resistor("RX", "x", "0", r2).
		Resistor("RY", "y", "0", r2).
		Resistor("RZ", "zs", "z", r1).
		VoltageSource("VZ", "zs", "0", r1*s.M0).
		Resistor("RP", "p", "0", 1e6).
		VoltageSource("VW", "w", "0", 1).
		BuildTransient()
}

// Run solves the circuit to stop with the given step and returns the
// magnetization series and the pickup voltage.
type Run struct {
	Time       []float64
	Mx, My, Mz []float64
	Pickup     []float64 // pickup winding voltage [V]
	Ledger     output.TranLedger
}

func (s Spin) Run(stop, step float64) (Run, error) {
	c, err := s.Circuit()
	if err != nil {
		return Run{}, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: step, Stop: stop, UIC: true},
		analysis.TranOptions{Options: analysis.Options{RelTol: 1e-6}})
	if err != nil {
		return Run{}, err
	}
	// The winding currents return through the relaxation resistors:
	// i_x = −I(RX)·… simpler: i_x = v(x)/(−R2) since v(x) = −R2·i_x.
	vx, _ := wf.NodeVoltage("x")
	vy, _ := wf.NodeVoltage("y")
	vz, _ := wf.NodeVoltage("z")
	vp, _ := wf.NodeVoltage("p")
	r := Run{Time: wf.Time, Ledger: wf.Ledger}
	r2, r1 := s.L/s.T2, s.L/s.T1
	for k := range wf.Time {
		r.Mx = append(r.Mx, -vx[k]/r2)
		r.My = append(r.My, -vy[k]/r2)
		r.Mz = append(r.Mz, (r1*s.M0-vz[k])/r1)
		r.Pickup = append(r.Pickup, vp[k])
	}
	return r, nil
}

func last(x []float64) float64 { return x[len(x)-1] }

// Final returns the magnetization at the end of the run.
func (r Run) Final() [3]float64 { return [3]float64{last(r.Mx), last(r.My), last(r.Mz)} }

// Frequency estimates the precession frequency [rad/s] of Mx from its
// zero crossings.
func Frequency(t, x []float64) float64 {
	var cross []float64
	for k := 1; k < len(x); k++ {
		if (x[k-1] < 0) != (x[k] < 0) && x[k-1] != 0 {
			f := x[k-1] / (x[k-1] - x[k])
			cross = append(cross, t[k-1]+f*(t[k]-t[k-1]))
		}
	}
	if len(cross) < 2 {
		return 0
	}
	return math.Pi * float64(len(cross)-1) / (last(cross) - cross[0])
}

// BlochSteady is Bloch's steady-state solution in the rotating frame for
// detuning dw = ω_rf − γB₀ (sign by this program's convention): the
// transverse components (dispersion Mx, absorption My) and Mz.
func BlochSteady(s Spin, dw float64) (mx, my, mz float64) {
	w1 := s.Gamma * s.B1
	d := 1 + dw*dw*s.T2*s.T2 + w1*w1*s.T1*s.T2
	mx = -s.M0 * w1 * s.T2 * s.T2 * dw / d
	my = s.M0 * w1 * s.T2 / d
	mz = s.M0 * (1 + dw*dw*s.T2*s.T2) / d
	return
}

// ---------- animation ----------

type canvas struct {
	img  *image.RGBA
	w, h int
}

func newCanvas(w, h int) *canvas {
	c := &canvas{img: image.NewRGBA(image.Rect(0, 0, w, h)), w: w, h: h}
	draw.Draw(c.img, c.img.Bounds(), image.White, image.Point{}, draw.Src)
	return c
}

func (c *canvas) line(x0, y0, x1, y1 float64, col color.Color) {
	n := int(math.Max(math.Abs(x1-x0), math.Abs(y1-y0))) + 1
	for k := 0; k <= n; k++ {
		f := float64(k) / float64(n)
		x, y := x0+f*(x1-x0), y0+f*(y1-y0)
		if x >= 0 && y >= 0 && int(x) < c.w && int(y) < c.h {
			c.img.Set(int(x), int(y), col)
		}
	}
}

func (c *canvas) thick(x0, y0, x1, y1 float64, col color.Color) {
	for _, d := range []float64{-1, 0, 1} {
		c.line(x0+d, y0, x1+d, y1, col)
		c.line(x0, y0+d, x1, y1+d, col)
	}
}

func (c *canvas) dot(x, y, r float64, col color.Color) {
	for dy := -r; dy <= r; dy++ {
		for dx := -r; dx <= r; dx++ {
			if dx*dx+dy*dy <= r*r {
				c.img.Set(int(x+dx), int(y+dy), col)
			}
		}
	}
}

// DrawGIF animates the rotating-frame magnetization (left: vector on the
// unit sphere, oblique view, z up, x towards the viewer's lower left) over
// the pulse run then the decay run, and the pickup voltage of the decay
// run drawn progressively (right).
func DrawGIF(path string, pulse, fid Run, wrf float64, frames int) error {
	const W, H = 1000, 500
	blue := color.RGBA{30, 60, 200, 255}
	red := color.RGBA{200, 40, 40, 255}
	grey := color.RGBA{160, 160, 160, 255}
	black := color.RGBA{0, 0, 0, 255}
	pal := color.Palette{color.White, blue, red, grey, black}
	// Rotating-frame components of the decay run (solved in the lab).
	type pt struct{ x, y, z float64 }
	var seq []pt
	var tseq []float64
	for k := range pulse.Time {
		seq = append(seq, pt{pulse.Mx[k], pulse.My[k], pulse.Mz[k]})
		tseq = append(tseq, pulse.Time[k])
	}
	tp := last(pulse.Time)
	for k := range fid.Time {
		// In the lab M⊥ turns as e^(−iγB₀t); multiplying by e^(+iω_rf t)
		// brings it to the rotating frame.
		c, s := math.Cos(wrf*fid.Time[k]), math.Sin(wrf*fid.Time[k])
		seq = append(seq, pt{fid.Mx[k]*c - fid.My[k]*s, fid.Mx[k]*s + fid.My[k]*c, fid.Mz[k]})
		tseq = append(tseq, tp+fid.Time[k])
	}
	total := last(tseq)
	proj := func(p pt) (float64, float64) {
		const cx, cy, r = 250.0, 260.0, 170.0
		return cx + (p.y-0.5*p.x)*r, cy - (p.z-0.3*p.x)*r
	}
	vmax := 0.0
	for _, v := range fid.Pickup {
		vmax = math.Max(vmax, math.Abs(v))
	}
	g := &gif.GIF{}
	for f := 0; f < frames; f++ {
		tt := total * float64(f) / float64(frames-1)
		cv := newCanvas(W, H)
		// Sphere: equator and the x–z and y–z meridians.
		for k := 0; k < 360; k++ {
			a0, a1 := 2*math.Pi*float64(k)/360, 2*math.Pi*float64(k+1)/360
			for _, ring := range []func(a float64) pt{
				func(a float64) pt { return pt{math.Cos(a), math.Sin(a), 0} },
				func(a float64) pt { return pt{math.Cos(a), 0, math.Sin(a)} },
				func(a float64) pt { return pt{0, math.Cos(a), math.Sin(a)} },
			} {
				x0, y0 := proj(ring(a0))
				x1, y1 := proj(ring(a1))
				cv.line(x0, y0, x1, y1, grey)
			}
		}
		ox, oy := proj(pt{})
		zx, zy := proj(pt{0, 0, 1.25})
		cv.line(ox, oy, zx, zy, black) // B₀ axis
		// Trail and vector.
		var idx int
		for idx = 0; idx < len(tseq)-1 && tseq[idx+1] <= tt; idx++ {
		}
		for k := 1; k <= idx; k++ {
			x0, y0 := proj(seq[k-1])
			x1, y1 := proj(seq[k])
			cv.line(x0, y0, x1, y1, blue)
		}
		mx, my := proj(seq[idx])
		cv.thick(ox, oy, mx, my, red)
		cv.dot(mx, my, 5, red)
		// Pickup trace of the decay run, time from the end of the pulse.
		const x0, x1, ymid, amp = 540.0, 980.0, 260.0, 180.0
		cv.line(x0, ymid, x1, ymid, grey)
		cv.line(x0, ymid-amp, x0, ymid+amp, grey)
		for k := 1; k < len(fid.Time) && tp+fid.Time[k] <= tt; k++ {
			ax := x0 + (x1-x0)*fid.Time[k-1]/last(fid.Time)
			bx := x0 + (x1-x0)*fid.Time[k]/last(fid.Time)
			cv.line(ax, ymid-amp*fid.Pickup[k-1]/vmax, bx, ymid-amp*fid.Pickup[k]/vmax, blue)
		}
		// Progress bar.
		cv.line(0, H-3, W*tt/total, H-3, black)
		cv.line(0, H-2, W*tt/total, H-2, black)
		fr := image.NewPaletted(image.Rect(0, 0, W, H), pal)
		draw.Draw(fr, fr.Bounds(), cv.img, image.Point{}, draw.Src)
		g.Image = append(g.Image, fr)
		g.Delay = append(g.Delay, 6)
	}
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	fh, err := os.Create(path)
	if err != nil {
		return err
	}
	defer fh.Close()
	return gif.EncodeAll(fh, g)
}

// Proton is the main program's spin: a proton at 1 T (42.58 MHz) with
// short relaxation times so the runs stay small (the equations do not
// care), M0 = 1, L = 1 nH so that the winding voltages ωL·M stay at
// volts, pickup mutual 0.1 L.
func Proton() Spin {
	return Spin{Gamma: GammaProton, B0: 1, T1: 5e-6, T2: 2e-6, M0: 1, Minit: [3]float64{0, 0, 1}, L: 1e-9, Mp: 1e-10}
}

func main() {
	gifPath := flag.String("gif", "", "write the animation to this GIF")
	flag.Parse()
	fmt.Println("Magnetic resonance as Kron's frame II(b): Bloch's equations are L·dM/dt = −G·M with G = Lγ[B]ₓ, shaft at ω = 1")
	fmt.Println()

	s := Proton()
	w0 := s.Gamma * s.B0
	period := 2 * math.Pi / w0
	fmt.Printf("1. Free precession, laboratory frame: proton, B₀ = 1 T, T1 = %.0f µs, T2 = %.0f µs, M tipped to x\n", s.T1*1e6, s.T2*1e6)
	lab := s
	lab.Minit = [3]float64{1, 0, 0}
	r, err := lab.Run(4*s.T2, period/400)
	if err != nil {
		log.Fatal(err)
	}
	f := Frequency(r.Time, r.Mx)
	k := len(r.Time) - 1
	mxy := math.Hypot(r.Mx[k], r.My[k])
	fmt.Printf("   precession %.6f MHz (γB₀/2π = %.6f MHz, %.1e); |Mxy| at 4T2 %.5f (e⁻⁴ = %.5f); Mz %.5f (1 − e⁻⁴ᵀ²ᐟᵀ¹ = %.5f); %d points\n",
		f/2/math.Pi/1e6, w0/2/math.Pi/1e6, f/w0-1, mxy, math.Exp(-4), r.Mz[k], 1-math.Exp(-4*s.T2/s.T1), len(r.Time))
	fmt.Printf("   ledger: source %.4g J, heat %.4g J, stored %.4g → %.4g J, mismatch %.1e J\n",
		r.Ledger.EnergyIn-r.Ledger.EnergyOut, r.Ledger.EnergyDissipated, r.Ledger.StoredStart, r.Ledger.StoredEnd, r.Ledger.IntegrationMismatch)
	fmt.Println()

	fmt.Println("2. Rabi nutation at resonance, rotating frame: B₁ = 1 mT along x, relaxation off (T1 = T2 = 1 s)")
	rabi := s
	rabi.B1, rabi.Wrf, rabi.T1, rabi.T2 = 1e-3, w0, 1, 1
	wr := rabi.Gamma * rabi.B1
	trabi := 2 * math.Pi / wr
	r, err = rabi.Run(2*trabi, trabi/400)
	if err != nil {
		log.Fatal(err)
	}
	worst := 0.0
	for k := range r.Time {
		worst = math.Max(worst, math.Abs(r.Mz[k]-math.Cos(wr*r.Time[k])))
	}
	fmt.Printf("   Rabi period %.4f µs (2π/γB₁); Mz − cos(γB₁t) at worst %.1e over two periods; shaft power = 0 (iᵀGi = 0)\n", trabi*1e6, worst)
	fmt.Println()

	fmt.Println("3. π/2 pulse, then free-induction decay with the pickup coil")
	pulse := rabi
	tp := math.Pi / (2 * wr)
	pr, err := pulse.Run(tp, tp/400)
	if err != nil {
		log.Fatal(err)
	}
	m := pr.Final()
	fid := s
	fid.Minit = m
	fr, err := fid.Run(4*s.T2, period/400)
	if err != nil {
		log.Fatal(err)
	}
	vmax := 0.0
	for _, v := range fr.Pickup {
		vmax = math.Max(vmax, math.Abs(v))
	}
	fmt.Printf("   after the pulse (%.4f µs): M = (%.4f, %.4f, %.4f); FID at %.4f MHz, pickup peak %.4g V (M_p·γB₀·M₀ = %.4g V), decaying with T2\n",
		tp*1e6, m[0], m[1], m[2], Frequency(fr.Time, fr.Pickup)/2/math.Pi/1e6, vmax, s.Mp*w0*s.M0)
	fmt.Println()

	fmt.Println("4. Continuous-wave line: steady state in the rotating frame against the detuning (B₁ = 10 µT, T1 = 5 µs, T2 = 2 µs)")
	fmt.Printf("   %12s %12s %12s %12s %12s %12s\n", "Δω·T2", "My (abs.)", "Bloch", "Mx (disp.)", "Bloch", "Mz")
	cw := s
	cw.B1 = 1e-5
	for _, dT2 := range []float64{-3, -1, -0.5, 0, 0.5, 1, 3} {
		dw := dT2 / s.T2
		cw.Wrf = w0 + dw
		r, err := cw.Run(12*s.T1, s.T2/50)
		if err != nil {
			log.Fatal(err)
		}
		bx, by, _ := BlochSteady(cw, dw)
		fmt.Printf("   %12.1f %12.5f %12.5f %12.5f %12.5f %12.5f\n", dT2, last(r.My), by, last(r.Mx), bx, last(r.Mz))
	}
	fmt.Println()

	fmt.Println("5. Resonance frequencies at 1 T: measured γ/2π, and the H21 ring with the particle's charge and mass (g = 1, γ = e/2m)")
	fmt.Printf("   proton   NMR %10.4f MHz measured, ring %10.4f MHz (g_p = 5.5857 measured)\n", GammaProton/2/math.Pi/1e6, QE/(2*MP)/2/math.Pi/1e6)
	fmt.Printf("   electron ESR %10.4f GHz measured, ring %10.4f GHz (g_e = 2.0023 measured)\n", GammaElectron/2/math.Pi/1e9, QE/(2*ME)/2/math.Pi/1e9)
	fmt.Println()

	if *gifPath != "" {
		if err := DrawGIF(*gifPath, pr, fr, w0, 100); err != nil {
			log.Fatal(err)
		}
		fmt.Printf("Animation written to %s\n", *gifPath)
	} else {
		fmt.Println("Animation: run with -gif viz/nmr.gif")
	}
}
