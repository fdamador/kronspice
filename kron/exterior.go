// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package kron

import (
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/device"
	"github.com/fdamador/kronspice/em"
	"github.com/fdamador/kronspice/exchange"
	"github.com/fdamador/kronspice/output"
)

// OpenPath is one element's open path to the exterior junction Ω of one
// channel (CLAUDE.md §3a.1 items 1–2): a junction-pair axis from the
// element to everything outside the circuit on that channel. Power is the
// power through it, + into the circuit from the exterior, − out. An
// inactive path (no coupling model) is present in the frame and carries
// exactly zero (Kron 1939 p. 485: no current I flows through it).
type OpenPath struct {
	Element string
	Channel exchange.ChannelID
	Status  exchange.PortStatus
	Active  bool
	PowerW  float64
}

// Exterior is the exchange of a solved network with the outside: one open
// path per element and channel, the energy ledger they make up, and Kron's
// invariants of the same solution.
type Exterior struct {
	Paths  []OpenPath
	Ledger *output.Ledger
	// Result holds the element rows the ledger was built from (voltages,
	// currents, absorbed and stored power, ports), as from analysis.OP.
	Result *output.Result
	// Invariants and the frame they were computed in (DC only).
	Invariants Invariants
}

// ExteriorDC solves the operating point in the given form and builds its
// exterior: each element's port powers come from the element's own port
// model (device.ExteriorPower, the same model the MNA ledger uses) applied
// to the voltage and current Kron's frame gives for its branch; a
// radiation edge's power leaves through its element's light/radiation
// path. The ledger is the sum over the open paths — Kron's power through
// the junction pairs to the exterior — plus storage and model-supplied
// power.
func ExteriorDC(n *Network, form Form) (*Exterior, error) {
	sol, err := n.Solve(form, 0, true)
	if err != nil {
		return nil, err
	}
	f, err := n.Frame()
	if err != nil {
		return nil, err
	}
	iv, err := n.Invariants(sol, f, 0, true)
	if err != nil {
		return nil, err
	}
	ops := make([]device.OperatingPoint, len(n.devs))
	for ei, k := range n.elemBranch {
		if k < 0 {
			continue
		}
		E, I := real(sol.E.At(k)), real(sol.I.At(k))
		op := device.OperatingPoint{Voltage: E, Current: I, PowerAbsorbed: real(n.absorbed(sol, k))}
		if _, ok := n.devs[ei].(device.Reactive); ok {
			op.StoragePower = E * I
		}
		ops[ei] = op
	}
	rad := make([]float64, len(n.devs)) // radiated power per element, ≥ 0
	for _, r := range n.radiation {
		rad[r.elem] += real(sol.E.At(r.branch) * sol.I.At(r.branch))
	}
	ext := n.exterior(ops, rad)
	ext.Invariants = iv
	v, err := n.voltages(sol)
	if err != nil {
		return nil, err
	}
	for i := range ext.Result.NodeVoltages {
		ext.Result.NodeVoltages[i] = real(v[i])
	}
	return ext, nil
}

// ExteriorAC is the average-power exterior at one frequency [Hz]: phasors
// are peak amplitudes, so each element's average absorbed power is
// ½·Re(E·i*), and the port models are applied to it as the MNA AC ledger
// does (they are linear in absorbed power).
func ExteriorAC(n *Network, form Form, freq float64) (*Exterior, error) {
	sol, err := n.Solve(form, 2*math.Pi*freq, false)
	if err != nil {
		return nil, err
	}
	ops := make([]device.OperatingPoint, len(n.devs))
	for ei, k := range n.elemBranch {
		if k < 0 {
			continue
		}
		E, I := sol.E.At(k), sol.I.At(k)
		ops[ei] = device.OperatingPoint{Voltage: cmplx.Abs(E), Current: cmplx.Abs(I),
			PowerAbsorbed: n.absorbedAC(sol, k)}
	}
	rad := make([]float64, len(n.devs))
	for _, r := range n.radiation {
		rad[r.elem] += 0.5 * real(sol.E.At(r.branch)*cmplx.Conj(sol.I.At(r.branch)))
	}
	return n.exterior(ops, rad), nil
}

// exterior assembles the open paths and the ledger from per-element
// operating points.
func (n *Network) exterior(ops []device.OperatingPoint, rad []float64) *Exterior {
	c := n.c
	res := &output.Result{Title: c.Title, Analysis: "op", NodeNames: append([]string(nil), c.NodeNames...),
		NodeVoltages: make([]float64, c.NumNodes())}
	ext := &Exterior{Result: res}
	for ei, d := range n.devs {
		e := c.Elements[ei]
		op := ops[ei]
		ports := d.Ports()
		if rad[ei] != 0 {
			// The radiation edge makes the element's light/radiation path
			// active.
			ports = append([]exchange.ExternalPort(nil), ports...)
			for i := range ports {
				if ports[i].Channel == exchange.Light {
					ports[i].Status = exchange.Modeled
					ports[i].Direction = exchange.Out
					ports[i].Note = "radiation-resistance edge (kron)"
				}
			}
		}
		er := output.ElementResult{Name: e.Name, Kind: e.Kind.String(),
			Nodes:   [2]string{c.NodeNames[e.Nodes[0]], c.NodeNames[e.Nodes[1]]},
			Voltage: op.Voltage, Current: op.Current, PowerAbsorbed: op.PowerAbsorbed, Ports: ports}
		if rd, ok := d.(device.Reactive); ok {
			er.Stores = true
			er.StoragePower = op.StoragePower
			if _, mutual := d.(*device.Mutual); !mutual {
				er.StoredEnergy = rd.StoredEnergy(op)
			}
		}
		if ms, ok := d.(device.ModelSupplier); ok {
			er.ModelSupplied = ms.ModelSuppliedPower(op)
		}
		for _, p := range ports {
			pw := 0.0
			if p.Status != exchange.NotModeled {
				pw = d.ExteriorPower(p.Channel, op)
				if p.Channel == exchange.Light {
					pw -= rad[ei]
				}
			}
			er.ExteriorPower[p.Channel] = pw
			ext.Paths = append(ext.Paths, OpenPath{Element: e.Name, Channel: p.Channel, Status: p.Status,
				Active: p.Status != exchange.NotModeled, PowerW: pw})
		}
		res.Elements = append(res.Elements, er)
	}
	res.Ledger = em.BuildLedger(res)
	ext.Ledger = res.Ledger
	return ext
}

// Junction returns the open paths that meet at one channel's exterior
// junction Ω_channel, and the net power through them (+ into the circuit).
func (x *Exterior) Junction(c exchange.ChannelID) (paths []OpenPath, net float64) {
	for _, p := range x.Paths {
		if p.Channel == c {
			paths = append(paths, p)
			net += p.PowerW
		}
	}
	return paths, net
}
