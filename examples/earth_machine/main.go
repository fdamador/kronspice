// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command earth_machine models the Earth's core dynamo as a Kron machine:
// the Rikitake two-disk dynamo (T. Rikitake, "Oscillations of a system of
// disk dynamos," Proc. Camb. Phil. Soc. 54, 89, 1958), two homopolar disks
// each excited by the other's field coil. Part I of the "Earth as an
// electromagnetic machine" study (CLAUDE.md in this folder).
//
// The equations (dimensional; L coil inductance, R loop resistance, M the
// disk-coil coupling, C shaft inertia, G drive torque):
//
//	L·İ1 = −R·I1 + M·ω1·I2        C·ω̇1 = G1 − M·I1·I2
//	L·İ2 = −R·I2 + M·ω2·I1        C·ω̇2 = G2 − M·I1·I2
//
// With L = R = M = C = G = 1 these are Rikitake's dimensionless equations
// (x = I1, y = I2, z = ω1, A = ω1 − ω2, μ = 1).
//
// The circuit: each disk is kronspice's machine element (N) with two
// windings — its own radial path and the OTHER loop's field coil — and a
// speed matrix entry G[0][1] = ±M, which gives the speed voltage ω·M·i_coil
// in the disk and the back torque M·i_disk·i_coil on the shaft. Each shaft
// is a mechanical node with inertia C (a capacitor) and drive torque G (a
// current source on the mechanical channel), as in eckardt_homopolar.
//
// The program runs the circuit and an independent RK4 integration of the
// same equations, compares them, and runs three cases:
//
//  1. co-rotating disks (Rikitake's machine);
//  2. disk 2 driven the other way, same wiring;
//  3. disk 2 driven the other way, coil 2 rewound.
//
// The model is dimensionless: it shows that a self-excited two-disk
// dynamo reverses its field on its own, not when the Earth's field does.
//
//	go run ./examples/earth_machine [-t 1000] [-csv out.csv]
package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Params are one two-disk machine.
type Params struct {
	L, R, M, C float64 // coil inductance, loop resistance, coupling, shaft inertia
	G1, G2     float64 // drive torques on the two shafts
	Wind2      float64 // +1: coil 2 wound as Rikitake; −1: rewound
	I1, I2     float64 // initial currents
	W1, W2     float64 // initial shaft speeds
	// K adds a feedback coil around each disk carrying the disk's own
	// current (paper107 §14.5.2.3; Bullard's single-disk dynamo inside
	// Rikitake's): speed EMF K·M·ωᵢ·Iᵢ and back torque K·M·Iᵢ². K > 0 aids
	// the field, K < 0 opposes it.
	K float64
	// GMod and Wd modulate both drive torques, Gᵢ(t) = Gᵢ + GMod·sin(Wd·t)
	// (paper107's time-varying rotation, §14.2–14.3).
	GMod, Wd float64
	// Speed, when set, drives each shaft at a prescribed speed instead of
	// a torque (paper107 §14.5.2.2: each disk on its own variable-frequency
	// drive): a voltage source on the mechanical node, whose current is the
	// torque the drive supplies.
	Speed *[2]circuit.Waveform
	// Stacked puts the two disks on one axis, each inside the other loop's
	// coil (stacked.go): every coil's stray field also threads its own disk
	// (KS, a fraction of M, a built-in feedback coil) and the two coils
	// couple (KM, a fraction of L). StackedCouplings derives both from the
	// geometry.
	Stacked bool
	KS, KM  float64
}

// Rikitake returns μ = 1, A = 1 from (x, y, z) = (1, 0.5, 1).
func Rikitake() Params {
	return Params{L: 1, R: 1, M: 1, C: 1, G1: 1, G2: 1, Wind2: 1, I1: 1, I2: 0.5, W1: 1, W2: 0}
}

// Circuit builds the two-loop machine.
//
//	loop 1: a1 →[N1 disk 1]→ b1 →R1→ c1 →[N2 coil 1]→ a1
//	loop 2: a2 →[N2 disk 2]→ b2 →R2→ c2 →[N1 coil 2]→ a2
func (p Params) Circuit() (*circuit.Circuit, error) {
	if p.Stacked {
		return p.stackedCircuit()
	}
	d1 := circuit.NewMachineSpec(2) // disk 1; winding 1 is coil 2 (field on disk 1)
	d1.L0[0][0], d1.L0[1][1] = 1e-9, p.L
	d1.G[0][1] = p.Wind2 * p.M
	// The element books ω·G·i as a voltage drop along each winding. The
	// cross term reads as a source because winding 1 (coil 2) runs against
	// loop 2; winding 0 runs with loop 1, so an aiding self-feedback is
	// −K·M, and its torque iᵀ·G·i = −K·M·I1² is the back torque.
	d1.G[0][0] = -p.K * p.M
	d1.IC[0], d1.IC[1] = p.I1, p.I2
	d2 := circuit.NewMachineSpec(2) // disk 2; winding 1 is coil 1 (field on disk 2)
	d2.L0[0][0], d2.L0[1][1] = 1e-9, p.L
	d2.G[0][1] = p.M
	d2.G[0][0] = -p.K * p.M
	d2.IC[0], d2.IC[1] = p.I2, p.I1

	b := circuit.NewBuilder("Rikitake two-disk dynamo").Domain(exchange.Mechanical, "w1", "w2")
	b.MachineModel("disk1", d1).MachineModel("disk2", d2)
	b.Machine("N1", "disk1", [][2]string{{"a1", "b1"}, {"a2", "c2"}}, "w1", "0")
	b.Machine("N2", "disk2", [][2]string{{"a2", "b2"}, {"a1", "c1"}}, "w2", "0")
	b.Resistor("R1", "b1", "c1", p.R).Resistor("R2", "b2", "c2", p.R)
	// The loops float; reference every loop node (UIC fixes winding currents).
	for _, n := range []string{"a1", "b1", "c1", "a2", "b2", "c2"} {
		b.Resistor("RREF"+n, n, "0", 1e7)
	}
	p.shafts(b)
	return b.Build()
}

// shafts adds the two shafts: inertia C, torque G (optionally modulated),
// or a prescribed speed.
func (p Params) shafts(b *circuit.Builder) {
	for k, w := range []struct {
		node   string
		w0, tq float64
	}{{"w1", p.W1, p.G1}, {"w2", p.W2, p.G2}} {
		n := fmt.Sprint(k + 1)
		if p.Speed != nil {
			e := circuit.Element{Name: "VFD" + n, Kind: circuit.VoltageSource, Value: p.Speed[k].Args[0], HasDC: true,
				Channel: exchange.Mechanical}
			if p.Speed[k].Kind != circuit.WaveNone {
				e.Wave = p.Speed[k]
			}
			b.Add(e, w.node, "0")
			continue
		}
		b.CapacitorIC("C"+n, w.node, "0", p.C, w.w0)
		b.Resistor("RW"+n, w.node, "0", 1e9) // DC path; windage 1e-9, negligible
		td := circuit.Element{Name: "TD" + n, Kind: circuit.CurrentSource, Value: w.tq, HasDC: true, Channel: exchange.Mechanical}
		if p.GMod != 0 {
			td.Wave = circuit.Waveform{Kind: circuit.WaveSin, Args: []float64{w.tq, p.GMod, p.Wd / (2 * math.Pi), 0, 0, 0}}
		}
		b.Add(td, "0", w.node)
	}
}

// Trace is a run: time, the two loop currents and the two shaft speeds.
type Trace struct{ T, I1, I2, W1, W2 []float64 }

// RunCircuit runs the kronspice transient.
func (p Params) RunCircuit(tEnd, dt float64) (Trace, float64, error) {
	c, err := p.Circuit()
	if err != nil {
		return Trace{}, 0, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: tEnd, MaxStep: dt, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return Trace{}, 0, err
	}
	var tr Trace
	tr.T = wf.Time
	tr.W1, _ = wf.NodeVoltage("w1")
	tr.W2, _ = wf.NodeVoltage("w2")
	i1, _ := wf.ElementCurrent("R1")
	i2, _ := wf.ElementCurrent("R2")
	// Side by side the resistor's current is booked against the loop (the
	// winding ICs set the loop running b → c through R as −I); stacked, the
	// sense windings run along the loops and R carries +I. Both are checked
	// against RunRK4.
	sign := -1.0
	if p.Stacked {
		sign = 1
	}
	for k := range i1 {
		tr.I1, tr.I2 = append(tr.I1, sign*i1[k]), append(tr.I2, sign*i2[k])
	}
	return tr, wf.Ledger.IntegrationMismatch, nil
}

// Books is an energy balance of a trace from its samples after t = 0
// (the UIC start point, where the winding currents are not yet in the
// resistors, is left out): drive work in, ohmic heat out, change in
// stored (shaft kinetic + coil magnetic) energy, and the residual.
func (p Params) Books(tr Trace) (in, heat, dStored, residual float64) {
	stored := func(k int) float64 {
		return 0.5*p.C*(tr.W1[k]*tr.W1[k]+tr.W2[k]*tr.W2[k]) + 0.5*p.L*(tr.I1[k]*tr.I1[k]+tr.I2[k]*tr.I2[k])
	}
	for k := 2; k < len(tr.T); k++ {
		dt := tr.T[k] - tr.T[k-1]
		in += 0.5 * dt * (p.G1*(tr.W1[k]+tr.W1[k-1]) + p.G2*(tr.W2[k]+tr.W2[k-1]))
		heat += 0.5 * dt * p.R * (tr.I1[k]*tr.I1[k] + tr.I1[k-1]*tr.I1[k-1] + tr.I2[k]*tr.I2[k] + tr.I2[k-1]*tr.I2[k-1])
	}
	dStored = stored(len(tr.T)-1) - stored(1)
	residual = in - heat - dStored
	return
}

// RunRK4 integrates the same equations directly, as the reference.
func (p Params) RunRK4(tEnd, dt float64) Trace {
	f := p.rhs
	s := [4]float64{p.I1, p.I2, p.W1, p.W2}
	var tr Trace
	n := int(math.Round(tEnd / dt))
	const sub = 10 // substeps per sample: the counter-rotating case spins up to ω·dt ≫ 1
	h := dt / sub
	for k := 0; k <= n; k++ {
		tr.T = append(tr.T, float64(k)*dt)
		tr.I1, tr.I2, tr.W1, tr.W2 = append(tr.I1, s[0]), append(tr.I2, s[1]), append(tr.W1, s[2]), append(tr.W2, s[3])
		for q := 0; q < sub; q++ {
			s = rk4(f, float64(k)*dt+float64(q)*h, s, h)
		}
	}
	return tr
}

// rhs is the right-hand side of the machine equations at time t.
func (p Params) rhs(t float64, s [4]float64) [4]float64 {
	i1, i2, w1, w2 := s[0], s[1], s[2], s[3]
	g := 0.0
	if p.GMod != 0 {
		g = p.GMod * math.Sin(p.Wd*t)
	}
	ks, km := p.K, 0.0
	if p.Stacked {
		ks, km = p.K+p.KS, p.KM
	}
	e1 := -p.R*i1 + p.Wind2*p.M*w1*i2 + ks*p.M*w1*i1
	e2 := -p.R*i2 + p.M*w2*i1 + ks*p.M*w2*i2
	// L·[[1, km], [km, 1]]·İ = e
	d := p.L * (1 - km*km)
	return [4]float64{
		(e1 - km*e2) / d,
		(e2 - km*e1) / d,
		(p.G1 + g - p.Wind2*p.M*i1*i2 - ks*p.M*i1*i1) / p.C,
		(p.G2 + g - p.M*i1*i2 - ks*p.M*i2*i2) / p.C,
	}
}

func rk4(f func(float64, [4]float64) [4]float64, t float64, s [4]float64, dt float64) [4]float64 {
	k1 := f(t, s)
	var s2, s3, s4 [4]float64
	for j := range s {
		s2[j] = s[j] + dt/2*k1[j]
	}
	k2 := f(t+dt/2, s2)
	for j := range s {
		s3[j] = s[j] + dt/2*k2[j]
	}
	k3 := f(t+dt/2, s3)
	for j := range s {
		s4[j] = s[j] + dt*k3[j]
	}
	k4 := f(t+dt, s4)
	for j := range s {
		s[j] += dt / 6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
	}
	return s
}

// Reversals returns the times I1 changes sign after tSkip, ignoring
// chatter about zero (|I1| below floor between crossings).
func Reversals(tr Trace, tSkip, floor float64) []float64 {
	var out []float64
	last := 0.0
	for k := range tr.T {
		if tr.T[k] < tSkip || math.Abs(tr.I1[k]) < floor {
			continue
		}
		s := math.Copysign(1, tr.I1[k])
		if last != 0 && s != last {
			out = append(out, tr.T[k])
		}
		last = s
	}
	return out
}

// Stats are the reversal count and interval min, mean, max and
// coefficient of variation (1 for a Poisson process).
func Stats(rev []float64) (n int, mn, mean, mx, cv float64) {
	n = len(rev)
	if n < 2 {
		return
	}
	mn, mx = math.Inf(1), 0
	var s, s2 float64
	for k := 1; k < n; k++ {
		d := rev[k] - rev[k-1]
		mn, mx = math.Min(mn, d), math.Max(mx, d)
		s += d
		s2 += d * d
	}
	m := float64(n - 1)
	mean = s / m
	cv = math.Sqrt(s2/m-mean*mean) / mean
	return
}

// Equilibrium is Rikitake's fixed point for μ = R = 1 (L = M = C = G = 1):
// I1 = k, I2 = 1/k, ω1 = k², with k² − k⁻² = A; and the oscillation
// period 2π/|Im λ| of the linearization about it.
func Equilibrium(A float64) (k2, period float64) {
	k2 = (A + math.Sqrt(A*A+4)) / 2
	// Characteristic polynomial of the Jacobian at the fixed point (μ = 1):
	// λ³ + 2λ² + (k² + k⁻²)λ + 2(k² + k⁻²) = (λ + 2)(λ² + k² + k⁻²).
	period = 2 * math.Pi / math.Sqrt(k2+1/k2)
	return
}

func absMax(v []float64, from int) float64 {
	m := 0.0
	for _, x := range v[from:] {
		m = math.Max(m, math.Abs(x))
	}
	return m
}

func main() {
	tEnd := flag.Float64("t", 1000, "run length (dimensionless time)")
	dt := flag.Float64("dt", 0.01, "time step")
	csv := flag.String("csv", "", "write the co-rotating circuit trace to this CSV")
	layer := flag.Int("layer", 0, "run one part only (1 geodynamo, 2 atmospheric circuit, 3 cavity, 4 layer-1 studies, 5 the layers linked, 6 stacked design); 0 runs all")
	data := flag.String("data", "", "write the figure data (CSV) to this directory and exit")
	flag.Parse()
	if *data != "" {
		if err := writeFigData(*data); err != nil {
			log.Fatal(err)
		}
		return
	}
	if *layer == 0 || *layer == 2 {
		if err := runAtmosphere(); err != nil {
			log.Fatalf("layer 2: %v", err)
		}
		if err := runAtmosphereRefined(); err != nil {
			log.Fatalf("layer 2 refinements: %v", err)
		}
	}
	if *layer == 0 || *layer == 3 {
		if err := runCavity(); err != nil {
			log.Fatalf("layer 3: %v", err)
		}
		if err := runCavityTwoHeight(); err != nil {
			log.Fatalf("layer 3 refinement: %v", err)
		}
	}
	if *layer == 0 || *layer == 5 {
		if err := runCoupled(); err != nil {
			log.Fatalf("links: %v", err)
		}
	}
	if *layer == 0 || *layer == 6 {
		if err := runStacked(); err != nil {
			log.Fatalf("stacked: %v", err)
		}
	}
	if *layer == 0 || *layer == 4 {
		if err := runStudies(); err != nil {
			log.Fatalf("studies: %v", err)
		}
	}
	if *layer != 0 && *layer != 1 {
		return
	}

	k2, per := Equilibrium(1)
	fmt.Println("Rikitake two-disk dynamo (μ = 1, A = 1) — the geodynamo as a Kron machine")
	fmt.Printf("  equilibrium: I1² = k² = %.6f (the golden ratio), natural period 2π/√(k²+k⁻²) = %.4f\n\n", k2, per)

	cases := []struct {
		name string
		p    Params
	}{
		{"1  co-rotating (Rikitake)", Rikitake()},
		{"2  disk 2 reversed, same wiring", func() Params { p := Rikitake(); p.G2, p.W2 = -1, -1; return p }()},
		{"3  disk 2 reversed, coil 2 rewound", func() Params { p := Rikitake(); p.G2, p.W2, p.Wind2, p.I2 = -1, 0, -1, -0.5; return p }()},
		{"4  co-rotating, coil 2 rewound (Tesla's opposed fields)", func() Params { p := Rikitake(); p.Wind2 = -1; return p }()},
	}
	for ci, c := range cases {
		tr, mis, err := c.p.RunCircuit(*tEnd, *dt)
		if err != nil {
			log.Fatalf("%s: %v", c.name, err)
		}
		ref := c.p.RunRK4(*tEnd, *dt)
		skip := 50.0
		i0 := 0
		for i0 < len(tr.T) && tr.T[i0] < skip {
			i0++
		}
		n, mn, mean, mx, cv := Stats(Reversals(tr, skip, 1e-2))
		nr, _, meanr, _, _ := Stats(Reversals(ref, skip, 1e-2))
		// Agreement over the first 10 time units (chaos separates the two
		// integrations after that; the statistics are what must agree).
		var dmax float64
		for k := 1; k < len(tr.T); k++ {
			if tr.T[k] > 10 {
				break
			}
			j := int(math.Round(tr.T[k] / *dt))
			if j < len(ref.T) {
				dmax = math.Max(dmax, math.Abs(tr.I1[k]-ref.I1[j]))
			}
		}
		fmt.Printf("Case %s\n", c.name)
		fmt.Printf("  kronspice: %d reversals, interval min %.2f mean %.2f max %.2f, CV %.2f (Poisson: 1)\n", n, mn, mean, mx, cv)
		fmt.Printf("  RK4 ref:   %d reversals, mean interval %.2f\n", nr, meanr)
		in, heat, ds, res := c.p.Books(tr)
		fmt.Printf("  max |I1| after t=50: %.3g; |ΔI1| vs RK4 over t ≤ 10: %.2g\n", absMax(tr.I1, i0), dmax)
		fmt.Printf("  energy from the trace (t > 0): drive in %.4g = heat %.4g + Δstored %.4g; residual %.2g (%.2g%%)\n", in, heat, ds, res, 100*res/in)
		fmt.Printf("  engine ledger (whole run, UIC start included): integration mismatch %.2g J (%.1g of the drive work)\n", mis, mis/in)
		fmt.Printf("  final shaft speeds ω1 = %.3g, ω2 = %.3g\n\n", tr.W1[len(tr.W1)-1], tr.W2[len(tr.W2)-1])
		if ci == 0 && *csv != "" {
			f, err := os.Create(*csv)
			if err != nil {
				log.Fatal(err)
			}
			fmt.Fprintln(f, "t,I1,I2,w1,w2")
			for k := range tr.T {
				fmt.Fprintf(f, "%g,%g,%g,%g,%g\n", tr.T[k], tr.I1[k], tr.I2[k], tr.W1[k], tr.W2[k])
			}
			f.Close()
		}
	}
}
