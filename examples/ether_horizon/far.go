// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"math"
	"strconv"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Parts 6 and 7: where the ether goes, and what a sink electron does to
// the rest of the world.
//
// Far from its core a sink draws the ether in radially: through every
// sphere the same volume flux Q, so v = Q/(4π r²). A body sitting in that
// inflow and swallowing ether itself takes the inflow's momentum: the
// Lagally force on sink 2 in the flow of sink 1 is ρe·Q2·v1, towards sink 1
// (L. M. Milne-Thomson, Theoretical Hydrodynamics, §§ 9.xx; Lagally 1922).
// Two sinks attract with ρe·Q1·Q2/(4π d²) — an inverse square, with the
// sign of gravity (like attracts like), the opposite of Coulomb's.

// Constants for the bounds.
const (
	EarthHeat      = 47e12                 // W, Earth's surface heat flow (Davies & Davies, Solid Earth 1, 5, 2010)
	EarthElectrons = 5.972e24 * 0.49 / Amu // about 0.49 electrons per nucleon mass
	Amu            = 1.66053906660e-27     // kg
	GMSunRate      = 1e-13 / Year          // 1/s: planetary ephemerides fit d(GM_sun)/dt/GM_sun to ~1e-13 per year (Pitjeva & Pitjev 2012, 2021)
	RhoCrit        = 8.5e-27               // kg/m³, critical density of the universe (H0 ≈ 67.7 km/s/Mpc)
)

// HeatBoundPerElectron is the most steady power an electron at rest can
// put out without heating the Earth past its measured heat flow (all of
// which already has known sources — radiogenic and primordial heat).
func HeatBoundPerElectron() float64 { return EarthHeat / EarthElectrons }

// HawkingPowerPaidByFlow is the one-channel power ħc²/(24π r_h²) of a
// horizon of radius rh, c_ether = c. In a sink flow the flow pays for it
// (as the pump does in laboratory analogues), so the electron need not
// lose rest energy — but the power still leaves as radiation.
func HawkingPowerPaidByFlow(rh float64) float64 {
	return OneChannelPower(HawkingT(rh, C0))
}

// HeatHorizonBound is the smallest horizon whose Hawking power meets the heat bound.
func HeatHorizonBound() float64 {
	return math.Sqrt(Hbar * C0 * C0 / (24 * math.Pi * HeatBoundPerElectron()))
}

// SinkElectron is an electron as a point sink of ether of density RhoE
// [kg/m³] whose mutual attraction is fitted to gravity: ρe·Q²/(4π) = G m_e².
type SinkElectron struct{ RhoE float64 }

// Q is the volume flux [m³/s].
func (s SinkElectron) Q() float64 { return Me * math.Sqrt(4*math.Pi*G/s.RhoE) }

// MassRate is the ether mass swallowed per second, ρe·Q = m_e·√(4πGρe) [kg/s].
func (s SinkElectron) MassRate() float64 { return s.RhoE * s.Q() }

// GrowthRate is the fractional rate √(4πGρe) at which the electron would
// grow if it kept the ether it swallows [1/s]. The same rate holds for any
// body whose gravity comes from the same sinks (Q ∝ mass), the Sun too.
func (s SinkElectron) GrowthRate() float64 { return math.Sqrt(4 * math.Pi * G * s.RhoE) }

// EnergyColumn is the energy taken per second if each kilogram of ether
// carries c² (a superfluid whose waves move at c has μ/m = c²) [W].
func (s SinkElectron) EnergyColumn() float64 { return s.MassRate() * C0 * C0 }

// Horizon is where the inflow Q/(4π r²) reaches c, √(Q/4πc) [m]: a
// horizon exists only if the sink's core is smaller than this.
func (s SinkElectron) Horizon() float64 { return math.Sqrt(s.Q() / (4 * math.Pi * C0)) }

// RhoForGrowth is the ether density whose growth rate equals rate.
func RhoForGrowth(rate float64) float64 { return rate * rate / (4 * math.Pi * G) }

// ----------------------------------------------------------- kronspice twin

// SinkGrid is the far field as a circuit: incompressible potential flow
// obeys Laplace's equation, as a resistor lattice does — node voltage =
// velocity potential, branch current = volume flow, KCL = continuity, and
// a sink = a current source drawing Q out of its node. An N³ lattice of
// unit resistors with the outer faces at 0 (ether at rest far away).
type SinkGrid struct {
	N      int
	Sinks  [][3]int // node coordinates
	Q      []float64
	Volt   []float64 // solved potentials, index (i*N + j)*N + k
	PIn    float64   // power the sink sources deliver [W]
	PRes   float64   // power in the resistors [W]
	Bound  float64   // current into the faces (the inflow from "infinity")
	Ledger float64   // kronspice's own relative residual
}

func gridName(i, j, k int) string {
	return "n" + strconv.Itoa(i) + "_" + strconv.Itoa(j) + "_" + strconv.Itoa(k)
}

// Solve builds the lattice in kronspice and runs a DC operating point.
func (g *SinkGrid) Solve() error {
	n := g.N
	b := circuit.NewBuilder("ether sink lattice")
	name := func(i, j, k int) string {
		if i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n {
			return "0" // the faces: potential 0
		}
		return gridName(i, j, k)
	}
	id := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				for _, d := range [][3]int{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}} {
					ii, jj, kk := i+d[0], j+d[1], k+d[2]
					inside := ii >= 0 && jj >= 0 && kk >= 0 && ii < n && jj < n && kk < n
					if inside && (d[0] < 0 || d[1] < 0 || d[2] < 0) {
						continue // each inner branch once
					}
					b.Resistor("R"+strconv.Itoa(id), name(i, j, k), name(ii, jj, kk), 1)
					id++
				}
			}
		}
	}
	for s, p := range g.Sinks {
		// The sink draws Q from its node to ground (the ether leaves the lattice there).
		b.CurrentSource("IQ"+strconv.Itoa(s), name(p[0], p[1], p[2]), "0", g.Q[s])
	}
	c, err := b.Build()
	if err != nil {
		return err
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return err
	}
	idx := make(map[string]float64, len(r.NodeNames))
	for i, nm := range r.NodeNames {
		idx[nm] = r.NodeVoltages[i]
	}
	g.Volt = make([]float64, n*n*n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				g.Volt[(i*n+j)*n+k] = idx[gridName(i, j, k)]
			}
		}
	}
	at := func(i, j, k int) float64 {
		if i < 0 || j < 0 || k < 0 || i >= n || j >= n || k >= n {
			return 0
		}
		return g.Volt[(i*n+j)*n+k]
	}
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				v := at(i, j, k)
				for _, d := range [][3]int{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} {
					w := at(i+d[0], j+d[1], k+d[2])
					g.PRes += (v - w) * (v - w)
				}
				for _, d := range [][3]int{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}} {
					ii, jj, kk := i+d[0], j+d[1], k+d[2]
					if ii < 0 || jj < 0 || kk < 0 {
						g.PRes += v * v
						g.Bound += 0 - v
					}
				}
				for _, d := range [][3]int{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}} {
					if i+d[0] >= n || j+d[1] >= n || k+d[2] >= n {
						g.Bound += 0 - v
					}
				}
			}
		}
	}
	for s, p := range g.Sinks {
		// A source drawing Q from a node at potential V delivers −V·Q.
		g.PIn += -at(p[0], p[1], p[2]) * g.Q[s]
	}
	g.Ledger = r.Ledger.RelativeResidual()
	return nil
}

// Speed is the inflow speed at node (i, j, k) along x from the central
// difference of the potential (flow = −∇φ in lattice units).
func (g *SinkGrid) Speed(i, j, k int) float64 {
	n := g.N
	at := func(i, j, k int) float64 { return g.Volt[(i*n+j)*n+k] }
	return -(at(i+1, j, k) - at(i-1, j, k)) / 2
}

// RhoForHeat is the ether density whose energy column, if it all became
// heat, would just meet the heat bound per electron.
func RhoForHeat() float64 {
	rate := HeatBoundPerElectron() / (Me * C0 * C0)
	return RhoForGrowth(rate)
}

// CoulombOverGravity is e²/(4πε0·G·m_e²) for two electrons.
func CoulombOverGravity() float64 { return E * E / (4 * math.Pi * Eps0 * G * Me * Me) }
