// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"errors"
	"fmt"
	"math"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
)

// Stage 2: a reduced cavitation model, calibrated on the published
// rotational reactors (LITERATURE.md §1), for a rotor band of cavities
// turning inside a housing. It is not CFD: every constant below is either
// fitted to a published number or declared, and the declared ones are run
// at both ends of a range in Sensitivity.

// Params are the model constants.
type Params struct {
	Cal Calibration
	// SigmaI is the inception cavitation number, the middle of the
	// calibration bracket; Delta sets how fast the vapour volume grows
	// below it and then chokes. Delta is fitted so that the
	// collapse-energy rate of Xue's rotor peaks at its best speed
	// (FitDelta).
	SigmaI, Delta float64
	// KappaV is the vapour volume at the peak of g(σ) as a fraction of the
	// cavity volume. Declared 0.2: Sun et al. 2021 report 5.44e-8 m³ of
	// vapour beside Ø 10 mm units (≈ 2.6e-7 m³ as hemispheres).
	KappaV float64
	// Cd is the drag coefficient of one cavity on its projected area
	// (½ρU²·Cd·A). Declared 0.25: the order Sun et al. 2021's 1.63 N·m
	// baseline torque gives if it is one 1/32 sector at 46 m/s.
	Cd float64
	// Cpulse is the pressure pulse of a passing housing cavity, as a
	// fraction of ½ρU², at gap G0. Declared 0.3 from Petkovšek et al.
	// 2013's 0.7 bar at 3.5 mm and 1.7 bar at 0.8 mm (as cited by Sun et
	// al. 2020) at an assumed 20–30 m/s; GapExp = 0.60 is fitted to those
	// two points (ln(1.7/0.7)/ln(3.5/0.8)).
	Cpulse, G0, GapExp float64
	// St are the cloud Strouhal numbers fW/U the design must work for:
	// 0.28 (cloud shedding, Stutz & Legoupil) and Xue's 0.446. A design is
	// scored at the worse of the two.
	St []float64
	// Mu is water's viscosity at 300 K (IAPWS R12-08: 853.8 µPa·s).
	Mu float64
	// PVap is the water's vapour pressure [Pa]; 0 means water at 300 K
	// (pVap). The closed steam loop (loop.go) sets it from the water's
	// temperature. Density and viscosity stay at their 300 K values there
	// (declared simplification).
	PVap float64
}

// pv returns the vapour pressure in use.
func (p Params) pv() float64 {
	if p.PVap > 0 {
		return p.PVap
	}
	return pVap
}

// DefaultParams calibrates on Xue 2024 and fills the declared constants.
func DefaultParams() Params {
	cal := Calibrate(Xue2024)
	si := 0.5 * (cal.SigmaOn + cal.SigmaNone)
	return Params{
		Cal: cal, SigmaI: si, Delta: FitDelta(si, cal.SigmaBest),
		KappaV: 0.2, Cd: 0.25, Cpulse: 0.3, G0: 1.5e-3, GapExp: math.Log(1.7/0.7) / math.Log(3.5/0.8),
		St: []float64{0.28, cal.StLock}, Mu: 853.8e-6,
	}
}

// G is the vapour-volume law: 0 above inception, x·e^(1−x) below it with
// x = (σ_i − σ)/Δ, so it rises to 1 at σ_i − Δ and falls again as the
// cavity chokes (supercavitation collapses outside the cavity, weakly).
func (p Params) G(sigma float64) float64 {
	if sigma >= p.SigmaI {
		return 0
	}
	x := (p.SigmaI - sigma) / p.Delta
	return x * math.Exp(1-x)
}

// FitDelta returns Δ such that rpm·g(σ(rpm)), the collapse-energy rate of
// a plain-housing rotor at fixed pressure (rate ∝ U ∝ rpm, energy per
// event ∝ (p − p_v)·V·g), peaks at σ_best. With σ ∝ rpm⁻²,
// d/d rpm [rpm·g] = 0 gives Δ = 2σ(σ_i − σ)/(σ_i + σ).
func FitDelta(sigmaI, sigmaBest float64) float64 {
	return 2 * sigmaBest * (sigmaI - sigmaBest) / (sigmaI + sigmaBest)
}

// Cavity is one kind of cavity in a row.
type Cavity struct {
	W     float64 // mouth width along the motion [m]
	V     float64 // volume [m³]
	A     float64 // projected area facing the flow [m²]
	Count int     // per row
}

// Hemicap returns a spherical-cap dimple of mouth Ø w and depth h ≤ w/2.
func Hemicap(w, h float64, count int) Cavity {
	a := w / 2
	rs := (a*a + h*h) / (2 * h) // sphere radius
	v := math.Pi * h * (3*a*a + h*h) / 6
	area := rs*rs*math.Acos((rs-h)/rs) - (rs-h)*math.Sqrt(2*rs*h-h*h) // circular segment
	return Cavity{W: w, V: v, A: area, Count: count}
}

// Machine is the rotor: two dimple drums, one at each outlet (the user's
// choice, 2026-11-08), each a cylinder of radius R carrying the cavity rows,
// with two annular end faces from the hub to R, turning in a fixed housing.
// The water leaves at the drum's end, where the outlet pressure is set.
type Machine struct {
	Name    string
	R       float64    // band radius [m]
	Rows    int        // rows on each band
	RowLen  float64    // axial length per row [m]
	Rotor   [][]Cavity // per row: the rotor's cavity kinds
	Housing []Cavity   // per row, the housing's cavities (nil: plain)
	Gap     float64    // rotor-housing gap [m]
	RPM     float64
	POut    float64 // pressure at the outlet, the drum's end [Pa]
	Hub     float64 // drum end faces run from Hub to R [m] (0: no end faces)
	Drums   int     // 2
}

// Eval is the machine's power split and line set.
type Eval struct {
	U, P, Sigma, G   float64
	PFriction, PDrag float64 // shaft power to wall friction and cavity drag [W]
	PCollapse        float64 // collapse-energy rate, a part of PDrag [W]
	PShaft           float64 // = PFriction + PDrag
	Eta              float64 // PCollapse / PShaft
	Cavities         int
	FPass            float64 // rotor-cavity passing frequency at the housing [Hz]
	FCloud           []float64
	TCollapse        float64 // Rayleigh collapse time of the largest cavity [s]
	Consistent       bool    // PCollapse ≤ PDrag
}

// tau is the wall shear stress at relative speed U across a gap: the
// larger of laminar Couette μU/gap and a turbulent flat-plate estimate at
// the mean speed U/2 (c_f = 0.0592 Re⁻⁰·², Schlichting; x = one
// circumference). An estimate, labelled.
func (p Params) tau(U, gap, x float64) float64 {
	lam := p.Mu * U / gap
	re := rho * (U / 2) * x / p.Mu
	if re < 1 {
		return lam
	}
	turb := 0.5 * rho * 0.0592 * math.Pow(re, -0.2) * (U / 2) * (U / 2)
	return math.Max(lam, turb)
}

// Evaluate scores a machine.
func (p Params) Evaluate(m Machine) (Eval, error) {
	if m.R <= 0 || m.Rows < 1 || m.Gap <= 0 || m.RPM <= 0 || m.Drums < 1 || len(m.Rotor) != m.Rows {
		return Eval{}, fmt.Errorf("bad machine %+v", m)
	}
	w := omega(m.RPM)
	U := w * m.R
	pr := m.POut
	s := (pr - p.pv()) / (0.5 * rho * U * U)
	e := Eval{U: U, P: pr, Sigma: s, G: p.G(s), FCloud: make([]float64, len(p.St)), Consistent: true}
	pulse := 0.0
	if len(m.Housing) > 0 {
		pulse = p.Cpulse * 0.5 * rho * U * U * math.Pow(p.G0/m.Gap, p.GapExp)
	}
	pDrive := pr - p.pv() + pulse
	if pDrive <= 0 {
		return e, errors.New("no collapse pressure")
	}
	nHousing := 0
	for _, c := range m.Housing {
		nHousing += c.Count
	}
	nRotor := 0
	if len(m.Rotor) > 0 {
		for _, c := range m.Rotor[0] {
			nRotor += c.Count
		}
	}
	e.FPass = float64(nRotor) * m.RPM / 60
	for i, st := range p.St {
		e.FCloud[i] = st * U / 0.01 // reference 10 mm, for the report
	}

	// collapse-energy rate of one cavity, at the worse Strouhal number;
	// the event rate is the cloud rate on a plain housing, and the passing
	// rate of the other side's cavities (capped by the cloud rate) on a
	// dimpled one; never faster than the Rayleigh collapse allows.
	one := func(c Cavity, others int) float64 {
		vc := p.KappaV * e.G * c.V
		if vc <= 0 {
			return 0
		}
		req := math.Cbrt(3 * vc / (4 * math.Pi))
		tc := 0.915 * req * math.Sqrt(rho/pDrive) // Rayleigh 1917
		if tc > e.TCollapse {
			e.TCollapse = tc
		}
		worst := math.Inf(1)
		for _, st := range p.St {
			f := st * U / c.W
			if others > 0 {
				f = math.Min(f, float64(others)*m.RPM/60)
			}
			f = math.Min(f, 1/(2*tc))
			worst = math.Min(worst, f)
		}
		return pDrive * vc * worst
	}

	// friction on the drum surface (rotor and housing walls, less the
	// cavity mouths) and on its two end faces, per drum
	bandArea := 2 * math.Pi * m.R * float64(m.Rows) * m.RowLen
	mouths := 0.0
	for _, row := range m.Rotor {
		for _, c := range row {
			mouths += float64(c.Count) * math.Pi * c.W * c.W / 4
			e.Cavities += c.Count
			e.PDrag += float64(c.Count) * 0.5 * rho * U * U * p.Cd * c.A * U
			e.PCollapse += float64(c.Count) * one(c, nHousing)
		}
		for _, c := range m.Housing {
			mouths += float64(c.Count) * math.Pi * c.W * c.W / 4
			e.Cavities += c.Count
			e.PDrag += float64(c.Count) * 0.5 * rho * U * U * p.Cd * c.A * U
			e.PCollapse += float64(c.Count) * one(c, nRotor)
		}
	}
	smooth := math.Max(0, 2*bandArea-mouths) / 2 // rotor-side equivalent of both walls' shear
	e.PFriction = p.tau(U, m.Gap, 2*math.Pi*m.R) * smooth * U
	const slices = 40
	if m.Hub > 0 && m.Hub < m.R {
		for k := 0; k < slices; k++ { // two end faces, hub to R
			r := m.Hub + (m.R-m.Hub)*(float64(k)+0.5)/slices
			dr := (m.R - m.Hub) / slices
			u := w * r
			e.PFriction += 2 * p.tau(u, m.Gap, 2*math.Pi*r) * 2 * math.Pi * r * dr * u
		}
	}
	n := float64(m.Drums)
	e.PFriction *= n
	e.PDrag *= n
	e.PCollapse *= n
	e.Cavities *= m.Drums
	e.PShaft = e.PFriction + e.PDrag
	e.Eta = e.PCollapse / e.PShaft
	e.Consistent = e.PCollapse <= e.PDrag*(1+1e-12)
	return e, nil
}

// ShaftCircuit is the machine's power split as a kronspice circuit, in the
// mobility analogy (CLAUDE.md §3a.1 item 3: angular speed across, torque
// through): a speed source ω drives one conductance per sink, G = T/ω, so
// each element absorbs ω·T. The operating point's ledger must return the
// shaft power with every watt booked.
func ShaftCircuit(e Eval, rpm float64) (*circuit.Circuit, error) {
	w := omega(rpm)
	b := circuit.NewBuilder("harmonic_cavitation shaft")
	b.VoltageSource("VSHAFT", "w", "0", w)
	add := func(name string, pw float64) {
		if pw > 0 {
			b.Resistor(name, "w", "0", w*w/pw)
		}
	}
	add("RFRICTION", e.PFriction)
	add("RDRAG_HEAT", e.PDrag-e.PCollapse)
	add("RCOLLAPSE", e.PCollapse)
	return b.Build()
}

// LedgerTwin solves the shaft circuit and returns the source's power and
// the ledger's relative residual.
func LedgerTwin(e Eval, rpm float64) (pIn, residual float64, err error) {
	c, err := ShaftCircuit(e, rpm)
	if err != nil {
		return 0, 0, err
	}
	r, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return 0, 0, err
	}
	src, _ := r.Element("VSHAFT")
	return -src.PowerAbsorbed, r.Ledger.RelativeResidual(), nil
}

// DimpleMachine builds two dimple drums of hemispherical (or shallower)
// dimples: rotor pitch pr·W, housing pitch ph·W (ph = 0: plain housing).
// The outlet is at the drum's end, at radius R; hub 10 mm (declared).
func DimpleMachine(r float64, rows int, w, hw, pr, ph, gap, rpm, pOut float64) Machine {
	nr := int(2 * math.Pi * r / (pr * w))
	m := Machine{Name: "dimple drums", R: r, Rows: rows, RowLen: w * 1.6, Gap: gap, RPM: rpm, POut: pOut,
		Hub: 0.010, Drums: 2}
	m.Rotor = make([][]Cavity, rows)
	for i := range m.Rotor {
		m.Rotor[i] = []Cavity{Hemicap(w, hw*w, nr)}
	}
	if ph > 0 {
		m.Housing = []Cavity{Hemicap(w, hw*w, int(2*math.Pi*r/(ph*w)))}
	}
	return m
}
