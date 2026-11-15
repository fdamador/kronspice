// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// Part 4: the same machine with the hypothetical ECE and ECE2
// spin-connection terms switched on (register entry H48, CLAUDE.md §3c;
// circuit/ece.go). The circuit engine attaches those terms to R, L, C and K
// only, not to the N machine element, so the coil's self-inductance is
// moved out of winding 3 into a series L element (the loop current is the
// same everywhere, so the circuit is unchanged; ECEPlainTwin checks it).
// Torsion then acts as a series resistance s·c·ω₀·L in the coil, and the
// scalar term as the factor 1 + s·κ on the resistor's terminal potentials
// (s = +1 ECE, −2 ECE2). No measured ω₀ or κ exists: the values are declared,
// sized like examples/ece_rlc (torsion = a quarter of R, κ = 0.05).

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// ECECase is one spin-connection setting.
type ECECase struct {
	Name  string
	Form  circuit.ECEForm
	W0    float64 // torsion ω₀ [1/m] on the coil
	Kappa float64 // scalar κ on both terminals of R
}

// ECETorsionOhms is the declared size of the torsion term, c·ω₀·L [Ω].
const ECETorsionOhms = 0.075

// ECEKappa is the declared scalar term on R.
const ECEKappa = 0.05

// ECECases lists the runs: plain, then torsion and scalar in each form.
func (p Params) ECECases() []ECECase {
	w0 := ECETorsionOhms / (circuit.SpeedOfLight * p.L)
	return []ECECase{
		{Name: "plain (Maxwell)"},
		{Name: "ECE torsion", Form: circuit.ECE, W0: w0},
		{Name: "ECE2 torsion", Form: circuit.ECE2, W0: w0},
		{Name: "ECE scalar", Form: circuit.ECE, Kappa: ECEKappa},
		{Name: "ECE2 scalar", Form: circuit.ECE2, Kappa: ECEKappa},
	}
}

// ECECircuit is Circuit with the coil's inductance as a series L element
// and the spin-connection terms of k.
func (p Params) ECECircuit(k ECECase, alpha, w0 float64) (*circuit.Circuit, error) {
	s := circuit.NewMachineSpec(4)
	s.L0[0][0], s.L0[1][1], s.L0[2][2], s.L0[3][3] = 1e-9, 1e-9, 1e-9, 1e-3
	s.G[0][3] = p.K1()
	s.G[1][2] = -p.K2()
	s.IC[3] = 1
	b := circuit.NewBuilder("Eckardt resonant homopolar, "+k.Name).Domain(exchange.Mechanical, "w")
	if k.Form != circuit.ECENone {
		b.ECE(circuit.ECESpec{Form: k.Form})
	}
	b.MachineModel("hg", s)
	b.CurrentSource("IF", "0", "f", 1).Resistor("RF", "f", "0", 1e3)
	b.Machine("NH", "hg", [][2]string{{"a", "b"}, {"b", "c"}, {"c", "e"}, {"f", "0"}}, "w", "0")
	lc := circuit.Element{Name: "LC", Kind: circuit.Inductor, Value: p.L, HasIC: true}
	if k.W0 != 0 {
		lc.Params = map[string]float64{"w0": k.W0}
	}
	b.Add(lc, "e", "d")
	r := circuit.Element{Name: "R", Kind: circuit.Resistor, Value: p.R}
	if k.Kappa != 0 {
		r.Params = map[string]float64{"kappa": k.Kappa}
	}
	b.Add(r, "d", "a")
	for _, n := range []string{"a", "b", "c", "d", "e"} {
		b.Resistor("RREF"+n, n, "0", 1e5)
	}
	b.CapacitorIC("J", "w", "0", p.J, w0)
	b.Resistor("RW", "w", "0", 1e9)
	if alpha != 0 {
		b.Add(circuit.Element{Name: "TD", Kind: circuit.CurrentSource, Value: alpha, HasDC: true, Channel: exchange.Mechanical}, "0", "w")
	}
	return b.Build()
}

// ECEResult is one run's trace and books.
type ECEResult struct {
	Case                     ECECase
	Tr                       Trace
	EHeatR, EModel, Mismatch float64 // heat in R, model-supplied (ECE terms), ledger mismatch [J]
	Label                    string
}

// RunECE runs case k as RunCircuit does.
func (p Params) RunECE(k ECECase, alpha, w0, tEnd float64) (ECEResult, error) {
	c, err := p.ECECircuit(k, alpha, w0)
	if err != nil {
		return ECEResult{}, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: tEnd / 2000, Stop: tEnd, MaxStep: tEnd / 2000, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return ECEResult{}, err
	}
	res := ECEResult{Case: k}
	for _, e := range c.Elements {
		if t := c.ECEFor(e); t.On {
			res.Label += e.Name + ": " + t.Label() + "\n"
		}
	}
	w, _ := wf.NodeVoltage("w")
	ir, _ := wf.ElementCurrent("R")
	var pr []float64
	for j, n := range wf.ElementNames {
		if n == "R" {
			pr = wf.P[j]
		}
	}
	res.Tr.T, res.Tr.W = wf.Time, w
	for i := range ir {
		res.Tr.I = append(res.Tr.I, -ir[i])
	}
	n := len(wf.Time) - 1
	for i := 1; i <= n; i++ {
		dt := wf.Time[i] - wf.Time[i-1]
		res.Tr.EOhm += 0.5 * dt * p.R * (ir[i]*ir[i] + ir[i-1]*ir[i-1])
		res.Tr.EDrive += 0.5 * dt * alpha * (w[i] + w[i-1])
		if pr != nil {
			res.EHeatR += 0.5 * dt * (pr[i] + pr[i-1])
		}
	}
	res.Tr.ERot = 0.5 * p.J * (w[n]*w[n] - w0*w0)
	res.Tr.EL = 0.5 * p.L * ir[n] * ir[n]
	res.EModel = wf.Ledger.EnergyModelSupplied
	res.Mismatch = wf.Ledger.IntegrationMismatch
	return res, nil
}

// ECEOmegaRes is the self-excitation speed with the terms: the loop
// resistance becomes R/(1 + s·κ) + s·c·ω₀·L.
func (p Params) ECEOmegaRes(k ECECase) float64 {
	s := k.Form.Factor()
	return (p.R/(1+s*k.Kappa) + s*circuit.SpeedOfLight*k.W0*p.L) / p.K2()
}

func runECE(p Params) {
	fmt.Println("4. ECE and ECE2 spin-connection terms switched on (hypothetical, H48; circuit/ece.go)")
	fmt.Println("   the terms attach to R, L, C, K, not to the N element, so the coil's L is a series L element here")
	fmt.Printf("   declared (no measured value exists): torsion c·ω₀·L = %.3f Ω (ω₀ = %.4g 1/m) on the coil; κ = %.2f on both terminals of R\n",
		ECETorsionOhms, ECETorsionOhms/(circuit.SpeedOfLight*p.L), ECEKappa)
	fmt.Println("   the loop resistance becomes R/(1 + s·κ) + s·c·ω₀·L, s = +1 ECE, −2 ECE2; the steady current is still set by α = K1·I + K2·I²")
	for _, run := range []struct {
		alpha, w0, tEnd float64
		title           string
	}{{p.Alpha, 0, 3, "driven, α = 1 N m, from rest, 3 s"}, {0, 40, 0.5, "flywheel, α = 0, ω0 = 40 rad/s, 0.5 s"}} {
		fmt.Printf("   %s:\n", run.title)
		fmt.Printf("     %-16s %8s %9s %9s %9s %10s %10s %10s %10s %9s\n", "case", "ω_res", "ω end", "I end", "peak I", "drive J", "flywheel J", "heat R J", "model J", "books J")
		for _, k := range p.ECECases() {
			r, err := p.RunECE(k, run.alpha, run.w0, run.tEnd)
			if err != nil {
				fmt.Printf("     %-16s error: %v\n", k.Name, err)
				continue
			}
			n := len(r.Tr.T) - 1
			pk := 0.0
			for _, i := range r.Tr.I {
				pk = math.Max(pk, math.Abs(i))
			}
			fmt.Printf("     %-16s %8.2f %9.3f %9.3f %9.3f %10.4g %10.4g %10.4g %+10.4g %9.1e\n",
				k.Name, p.ECEOmegaRes(k), r.Tr.W[n], r.Tr.I[n], pk, r.Tr.EDrive, -r.Tr.ERot, r.EHeatR, r.EModel, r.Mismatch)
		}
	}
	is := p.SteadyCurrent(p.Alpha)
	fmt.Printf("   steady state with drive: I = %.3f A in every case; the drive gives α·ω_res, R takes R_body·I², the ECE terms book the difference as model-supplied\n", is)
	fmt.Println("   ECE2 torsion is a negative series resistance: it lowers ω_res and R heats more than the drive supplies, the surplus booked to the spin connection, never as inflow")
}
