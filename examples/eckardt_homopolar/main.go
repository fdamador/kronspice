// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command eckardt_homopolar analyses H. Eckardt, "Simulation of resonances
// in the homopolar generator" (A.I.A.S./UPITEC, 13 February 2012), with
// kronspice.
//
// The machine (paper Fig. 1): two disks of radius r on one shaft, in series
// with a resistor R and a coil of N turns, length l, core permeability μr.
// The "static" disk sits in a permanent-magnet field Bm; the "dynamic" disk
// sits in the coil's field B_I = μ0·μr·N·I/l, so the machine's own current
// excites it. The paper's equations:
//
//	L·İ + R·I = ½·ω·r²·(Bm + B_I)            (6)
//	R_eff = R − ½·ω·r²·μ0·μr·N/l             (8)
//	ω_res = 2·R·l/(r²·μ0·μr·N)               (9)
//	J·ω̇ = α (constant drive torque)          (12)
//
// Eq. (12) leaves out the torque the disk currents feel in the fields they
// cross: the Lorentz force that gives the EMF ½ωr²B on the moving charges
// gives the torque ½r²B·I on the disk (per disk). With it,
//
//	J·ω̇ = α − ½·r²·(Bm + B_I)·I
//
// and ω·(torque) = EMF·I exactly: every joule into R and into the coil's
// field is taken from the shaft. The program runs:
//
//  1. the paper's equations as printed (RK4), with their energy books;
//  2. the same machine as a kronspice circuit — one machine element with the
//     two disk paths, the coil and the magnet's field as windings, the
//     shaft a mechanical node with inertia J and drive α — where the back
//     torque is part of the element;
//  3. the paper's second run (disk spun above ω_res, then connected).
//
// The physics of R_eff < 0 is standard: it is the self-excitation threshold
// of a series-wound dynamo (W. Siemens 1867; E. C. Bullard, "The stability
// of a homopolar dynamo," Proc. Camb. Phil. Soc. 51, 744, 1955; H. K.
// Moffatt, Magnetic Field Generation in Electrically Conducting Fluids,
// 1978, ch. 6). Above it the current grows by taking power from the drive.
//
//	go run ./examples/eckardt_homopolar
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Params are the paper's simulation values (p. 4), SI.
type Params struct {
	R, L, N, Len, Rad, Mu0, MuR, Bm, J, Alpha float64
}

// Paper returns the values printed on p. 4.
func Paper() Params {
	return Params{R: 0.3, L: 0.001, N: 800, Len: 0.05, Rad: 0.1, Mu0: 4 * math.Pi * 1e-7, MuR: 100, Bm: 0.1, J: 0.01, Alpha: 1}
}

// KI is the coil field per ampere, μ0·μr·N/l [T/A].
func (p Params) KI() float64 { return p.Mu0 * p.MuR * p.N / p.Len }

// K1 is the static disk's EMF per unit speed, ½r²Bm [V·s].
func (p Params) K1() float64 { return 0.5 * p.Rad * p.Rad * p.Bm }

// K2 is the dynamic disk's EMF per unit speed and ampere, ½r²·μ0μrN/l [V·s/A].
func (p Params) K2() float64 { return 0.5 * p.Rad * p.Rad * p.KI() }

// OmegaRes is eq. (9).
func (p Params) OmegaRes() float64 { return p.R / p.K2() }

// CoilInductance is the inductance of the stated coil if its field B_I
// fills the disk's area πr² (as eq. 4 assumes): μ0·μr·N²·πr²/l [H].
func (p Params) CoilInductance() float64 {
	return p.Mu0 * p.MuR * p.N * p.N * math.Pi * p.Rad * p.Rad / p.Len
}

// Trace is a time history.
type Trace struct {
	T, I, W                []float64
	EOhm, EDrive, ERot, EL float64 // energies at the end [J]
	EBackTorque            float64 // work done against the back torque [J]
	Blowup                 bool
}

// RunPaper integrates eqs. (6) and (12) as printed (no back torque), from
// current i0 and speed w0, with RK4, to tEnd or until |I| > iMax.
func RunPaper(p Params, i0, w0, alpha, tEnd, dt, iMax float64) Trace {
	f := func(i, w float64) (float64, float64) {
		return (0.5*w*p.Rad*p.Rad*(p.Bm+p.KI()*i) - p.R*i) / p.L, alpha / p.J
	}
	var tr Trace
	i, w := i0, w0
	for t := 0.0; t <= tEnd; t += dt {
		tr.T, tr.I, tr.W = append(tr.T, t), append(tr.I, i), append(tr.W, w)
		if math.Abs(i) > iMax {
			tr.Blowup = true
			break
		}
		a1, b1 := f(i, w)
		a2, b2 := f(i+dt/2*a1, w+dt/2*b1)
		a3, b3 := f(i+dt/2*a2, w+dt/2*b2)
		a4, b4 := f(i+dt*a3, w+dt*b3)
		in := i + dt/6*(a1+2*a2+2*a3+a4)
		wn := w + dt/6*(b1+2*b2+2*b3+b4)
		tr.EOhm += 0.5 * dt * p.R * (i*i + in*in)
		tr.EDrive += 0.5 * dt * alpha * (w + wn)
		// What the back torque would have absorbed: ω·½r²(Bm+B_I)·I.
		tr.EBackTorque += 0.5 * dt * (w*0.5*p.Rad*p.Rad*(p.Bm+p.KI()*i)*i + wn*0.5*p.Rad*p.Rad*(p.Bm+p.KI()*in)*in)
		i, w = in, wn
	}
	tr.ERot = 0.5 * p.J * (w*w - w0*w0)
	tr.EL = 0.5 * p.L * (i*i - i0*i0)
	return tr
}

// Circuit builds the machine for kronspice. Windings of the one machine
// element (one shaft, node "w"):
//
//	1  static disk path   a → b    EMF ω·K1·i_f
//	2  dynamic disk path  b → c    EMF ω·K2·I (from the coil current)
//	3  coil               c → d    self-inductance L
//	4  magnet equivalent  f → 0    1 A (a permanent magnet as a fixed current)
//
// R closes the loop d → a. The coil's polarity is the one that excites
// (R_eff = R − ωK2, as in the paper). The shaft carries inertia J (a
// capacitor in the mechanical domain) and the drive torque α (a current
// source booked on the mechanical channel), started at ω0.
func (p Params) Circuit(l, alpha, w0 float64) (*circuit.Circuit, error) {
	s := circuit.NewMachineSpec(4)
	s.L0[0][0], s.L0[1][1], s.L0[2][2], s.L0[3][3] = 1e-9, 1e-9, l, 1e-3
	s.G[0][3] = p.K1()
	s.G[1][2] = -p.K2()
	s.IC[3] = 1 // the magnet's field is present from t = 0
	b := circuit.NewBuilder("Eckardt resonant homopolar").Domain(exchange.Mechanical, "w")
	b.MachineModel("hg", s)
	b.CurrentSource("IF", "0", "f", 1).Resistor("RF", "f", "0", 1e3)
	b.Machine("NH", "hg", [][2]string{{"a", "b"}, {"b", "c"}, {"c", "d"}, {"f", "0"}}, "w", "0")
	b.Resistor("R", "d", "a", p.R)
	// Reference the floating loop. At the start of a UIC run the windings
	// are fixed currents, so each loop node needs its own reference; 100 kΩ
	// to ground (100 kΩ) leaks little; what it takes is reported.
	b.Resistor("RREF", "a", "0", 1e5).Resistor("RREFB", "b", "0", 1e5).Resistor("RREFC", "c", "0", 1e5).Resistor("RREFD", "d", "0", 1e5)
	b.CapacitorIC("J", "w", "0", p.J, w0)
	// Shaft node needs a DC path: windage of 1e-9 N·m per rad/s, negligible
	// (the paper has no friction).
	b.Resistor("RW", "w", "0", 1e9)
	if alpha != 0 {
		b.Add(circuit.Element{Name: "TD", Kind: circuit.CurrentSource, Value: alpha, HasDC: true, Channel: exchange.Mechanical}, "0", "w")
	}
	return b.Build()
}

// RunCircuit runs the kronspice transient and returns its trace and ledger.
func (p Params) RunCircuit(l, alpha, w0, tEnd float64) (Trace, float64, error) {
	c, err := p.Circuit(l, alpha, w0)
	if err != nil {
		return Trace{}, 0, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: tEnd / 2000, Stop: tEnd, MaxStep: tEnd / 2000, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return Trace{}, 0, err
	}
	var tr Trace
	w, _ := wf.NodeVoltage("w")
	ir, _ := wf.ElementCurrent("R")
	tr.T, tr.W = wf.Time, w
	for k := range ir {
		tr.I = append(tr.I, -ir[k]) // sign so I > 0 when R_eff < 0 drives it, as in the paper's figures
	}
	n := len(wf.Time) - 1
	for k := 1; k <= n; k++ {
		dt := wf.Time[k] - wf.Time[k-1]
		tr.EOhm += 0.5 * dt * p.R * (ir[k]*ir[k] + ir[k-1]*ir[k-1])
		tr.EDrive += 0.5 * dt * alpha * (w[k] + w[k-1])
	}
	tr.ERot = 0.5 * p.J * (w[n]*w[n] - w0*w0)
	tr.EL = 0.5 * l * ir[n] * ir[n]
	return tr, wf.Ledger.IntegrationMismatch, nil
}

// SteadyCurrent is the current at which the back torque balances α at the
// self-excitation speed: α = K1·I + K2·I².
func (p Params) SteadyCurrent(alpha float64) float64 {
	a, b := p.K2(), p.K1()
	return (-b + math.Sqrt(b*b+4*a*alpha)) / (2 * a)
}

func at(tr Trace, t float64) (float64, float64) {
	for k := range tr.T {
		if tr.T[k] >= t {
			return tr.I[k], tr.W[k]
		}
	}
	n := len(tr.T) - 1
	return tr.I[n], tr.W[n]
}

func main() {
	p := Paper()
	fmt.Println("H. Eckardt, \"Simulation of resonances in the homopolar generator\" (2012) — kronspice analysis")
	fmt.Printf("parameters (p. 4): R %.2f Ω, L %.3g H, N %.0f, l %.2f m, r %.2f m, μr %.0f, Bm %.2f T, J %.2f kg m², α %.1f N m\n",
		p.R, p.L, p.N, p.Len, p.Rad, p.MuR, p.Bm, p.J, p.Alpha)
	fmt.Printf("coil field %.3f T per ampere; K1 = ½r²Bm = %.3g V s; K2 = ½r²μ0μrN/l = %.4g V s/A\n", p.KI(), p.K1(), p.K2())
	fmt.Printf("ω_res = R/K2 = %.2f rad/s (paper: about 30; eq. 9)\n\n", p.OmegaRes())

	fmt.Println("1. The paper's equations as printed (J ω̇ = α, no back torque)")
	pp := RunPaper(p, 0, 0, p.Alpha, 0.6, 1e-6, 1e6)
	for _, t := range []float64{0.2, 0.3, 0.4, 0.45, 0.5} {
		i, w := at(pp, t)
		fmt.Printf("   t %.2f s: ω %6.2f rad/s, I %12.4g A\n", t, w, i)
	}
	fmt.Printf("   stopped at t = %.3f s (|I| > 1e6 A: %v)\n", pp.T[len(pp.T)-1], pp.Blowup)
	fmt.Printf("   energy into R %.4g J and the coil %.4g J; drive work α∫ω dt only %.3g J, of which %.3g J went to spin-up\n",
		pp.EOhm, pp.EL, pp.EDrive, pp.ERot)
	fmt.Printf("   the missing %.4g J is exactly the work the omitted back torque ½r²(Bm+B_I)·I would have taken from the shaft (%.4g J):\n",
		pp.EOhm+pp.EL-(pp.EDrive-pp.ERot), pp.EBackTorque)
	fmt.Println("   the \"energy from spacetime\" in these equations is the shaft work they leave out.")
	fmt.Println()

	fmt.Println("2. The same machine in kronspice (back torque inside the machine element)")
	for _, l := range []float64{p.L, p.CoilInductance()} {
		tEnd := 3.0
		if l > 1 {
			tEnd = 600
		}
		tr, mm, err := p.RunCircuit(l, p.Alpha, 0, tEnd)
		if err != nil {
			log.Fatal(err)
		}
		fmt.Printf("   L = %.3g H%s, %0.f s:\n", l, map[bool]string{true: " (the stated coil filling the disk area: μ0μrN²πr²/l)", false: " (as printed)"}[l > 1], tEnd)
		for _, f := range []float64{0.1, 0.2, 0.4, 0.7, 1} {
			i, w := at(tr, f*tEnd)
			fmt.Printf("     t %7.2f s: ω %7.3f rad/s, I %8.3f A, B_I %7.3f T\n", f*tEnd, w, i, p.KI()*math.Abs(i))
		}
		fmt.Printf("     drive work %.4g J = R %.4g J + coil %.4g J + spin-up %.4g J + reference resistors/windage %.3g J (ledger closes to %.1e J)\n",
			tr.EDrive, tr.EOhm, tr.EL, tr.ERot, tr.EDrive-tr.EOhm-tr.EL-tr.ERot, mm)
	}
	is := p.SteadyCurrent(p.Alpha)
	fmt.Printf("   steady state: ω → ω_res = %.2f rad/s, α = K1·I + K2·I² → I = %.3f A, P = α·ω_res = %.2f W = R·I² = %.2f W\n",
		p.OmegaRes(), is, p.Alpha*p.OmegaRes(), p.R*is*is)
	fmt.Printf("   that needs B_I = %.1f T in the iron; iron saturates near 1.5–2 T, so a real coil limits the excitation well before\n\n", p.KI()*is)

	fmt.Println("3. Paper's second run: disk at ω0 = 40 rad/s, then connected (paper Figs. 6–9)")
	for _, alpha := range []float64{0, p.Alpha} {
		tr, mm, err := p.RunCircuit(p.L, alpha, 40, 0.5)
		if err != nil {
			log.Fatal(err)
		}
		imax, tmax := 0.0, 0.0
		for k := range tr.I {
			if math.Abs(tr.I[k]) > imax {
				imax, tmax = math.Abs(tr.I[k]), tr.T[k]
			}
		}
		_, wEnd := at(tr, 0.5)
		fmt.Printf("   α = %.0f N m: peak I %.2f A at %.3f s; ω 40 → %.2f rad/s; R took %.3f J = flywheel %.3f J + drive %.3f J − coil %.3f J (books %.1e J)\n",
			alpha, imax, tmax, wEnd, tr.EOhm, -tr.ERot, tr.EDrive, tr.EL, mm)
	}
	fmt.Println("   the pulse is the flywheel's energy discharged through the self-excited dynamo: what R receives is what the shaft gave.")
}
