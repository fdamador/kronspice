// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// The three layers linked.
//
// (A) Layers 2 and 3 are one network. The global circuit is the DC response
// of the Earth–ionosphere system and the Schumann resonances its ELF
// response, so one Kron network carries both: every colatitude band of the
// cavity ladder gets, as its shunt branch, the atmospheric column of that
// band (the leaky shells of layer 2, two-scale profile), and the tops of the
// columns are joined band to band by the ionosphere's series path (μ0·h_m,
// with its loss). At DC the series inductors are shorts, the tops are one
// equipotential ionosphere and the columns are Wilson's circuit; at ELF the
// lower air is an insulator (ωε0 ≫ σ) and the column becomes a capacitor
// whose effective height is where σ(z) ≈ ωε0: the cavity's electric height
// is no longer fitted, it comes out of the conductivity profile.
//
// (B) Layer 1 drives layer 2. The air's conductivity is set by ions that
// galactic cosmic rays make; the dynamo's dipole shields the atmosphere
// from them, and the shielding (the cutoff rigidity) scales with the dipole
// moment (Störmer). With ion–ion recombination σ ∝ √q, so when the dipole
// weakens during a reversal the air conducts better and V_I = I·R falls.
// How much the global ion production rises without shielding is declared
// (QUnshielded), not measured: the chain is standard, its size is a
// parameter.
//
// (C) One ledger: the dynamo is driven through its mechanical channel, the
// storms through the thermal channel; the cavity is fed by lightning, a
// part of the storms' work.

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/exchange"
)

// Unified is the combined layer-2/3 network.
type Unified struct {
	Bands, Layers int
	Top           float64               // column height [m]
	Sigma         func(float64) float64 // fair-weather profile
	Hm            float64               // magnetic height for the series path [m]
	RLoss         float64               // series loss per unit of μ0·h_m·ω (Im h_m / Re h_m at the reference frequency)
	FRef          float64
}

func (u Unified) shellZ(k int) float64 {
	const stretch = 3.0
	return u.Top * (math.Exp(stretch*float64(k)/float64(u.Layers)) - 1) / (math.Exp(stretch) - 1)
}

// Circuit builds the network. stormDC is the storm current (DC, ground to
// the pole band's top); strokeAC is a unit AC stroke at the same place.
func (u Unified) Circuit(stormDC, strokeAC float64) (*circuit.Circuit, error) {
	d := math.Pi / float64(u.Bands)
	b := circuit.NewBuilder("Earth: global circuit and cavity in one network")
	top := func(k int) string { return fmt.Sprintf("t%d", k) }
	for k := 0; k < u.Bands; k++ {
		th0, th1 := float64(k)*d, float64(k+1)*d
		frac := (math.Cos(th0) - math.Cos(th1)) / 2 // fraction of the sphere
		for j := 0; j < u.Layers; j++ {
			z0, z1 := u.shellZ(j), u.shellZ(j+1)
			const m = 32
			h := (z1 - z0) / m
			var integ float64
			for i := 0; i <= m; i++ {
				w := 2.0
				if i == 0 || i == m {
					w = 1
				} else if i%2 == 1 {
					w = 4
				}
				integ += w / u.Sigma(z0+float64(i)*h)
			}
			integ *= h / 3
			r0, r1 := rEarth+z0, rEarth+z1
			rm := 0.5 * (r0 + r1)
			area := frac * 4 * math.Pi * rm * rm
			lo, hi := fmt.Sprintf("n%d_%d", k, j), fmt.Sprintf("n%d_%d", k, j+1)
			if j == 0 {
				lo = "0"
			}
			if j == u.Layers-1 {
				hi = top(k)
			}
			b.Resistor(fmt.Sprintf("R%d_%d", k, j), hi, lo, integ/area)
			b.Capacitor(fmt.Sprintf("C%d_%d", k, j), hi, lo, frac*4*math.Pi*eps0*r0*r1/(r1-r0))
		}
		if k < u.Bands-1 {
			per := d / (2 * math.Pi * math.Sin(th1))
			l := mu0 * u.Hm * per
			mid := fmt.Sprintf("s%d", k)
			b.Resistor(fmt.Sprintf("RS%d", k), top(k), mid, 2*math.Pi*u.FRef*l*u.RLoss)
			b.Inductor(fmt.Sprintf("LS%d", k), mid, top(k+1), l)
		}
	}
	b.Add(circuit.Element{Name: "ISTORM", Kind: circuit.CurrentSource, Value: stormDC, HasDC: true, ACMag: strokeAC, Channel: exchange.Thermal}, "0", top(0))
	return b.Build()
}

// ElectricHeight is where σ(z) = 2πf·ε0 in the column profile.
func (u Unified) ElectricHeight(f float64) float64 {
	target := 2 * math.Pi * f * eps0
	lo, hi := 0.0, 200e3
	for it := 0; it < 100; it++ {
		mid := 0.5 * (lo + hi)
		if u.Sigma(mid) < target {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

func runCoupled() error {
	fmt.Println("Linking the layers")
	a := StandardAtmosphere()
	a.Top = 100e3 // the column runs up into the lower ionosphere
	prof, _ := StandardAtmosphere().TwoScale(1500, 5000)
	t2, _ := fitTwoHeight(5)
	_, hm := t2.Heights(complex(8, 0))
	u := Unified{Bands: 61, Layers: 40, Top: 100e3, Sigma: prof, Hm: real(hm), RLoss: imag(hm) / real(hm), FRef: 8}
	fmt.Printf("  (A) one network: %d bands × %d shells to %.0f km (two-scale profile of layer 2), series path at h_m = %.0f km\n",
		u.Bands, u.Layers, u.Top/1e3, u.Hm/1e3)
	I := a.GlobalCurrent()
	c, err := u.Circuit(I, 1)
	if err != nil {
		return err
	}
	op, err := analysis.OP(c, analysis.Options{})
	if err != nil {
		return err
	}
	var vmin, vmax = math.Inf(1), 0.0
	for k := 0; k < u.Bands; k++ {
		v, _ := op.Voltage(fmt.Sprintf("t%d", k))
		vmin, vmax = math.Min(vmin, v), math.Max(vmax, v)
	}
	res := 0.0
	if op.Ledger != nil && op.Ledger.Scale > 0 {
		res = op.Ledger.Residual / op.Ledger.Scale
	}
	fmt.Printf("    DC: ionosphere %.2f–%.2f kV across all bands (one equipotential; layer 2 gives 249.7 kV), %.0f MW, ledger %.0e\n",
		vmin/1e3, vmax/1e3, vmin*I/1e6, res)
	ac, err := analysis.AC(c, analysis.ACSpec{Variation: "lin", Points: 391, FStart: 3, FStop: 42}, analysis.Options{})
	if err != nil {
		return err
	}
	vp, _ := ac.NodeVoltage("t0")
	ve, _ := ac.NodeVoltage(fmt.Sprintf("t%d", u.Bands/2))
	var mp, me []float64
	for k := range ac.Freq {
		mp, me = append(mp, cmplx.Abs(vp[k])), append(me, cmplx.Abs(ve[k]))
	}
	pp, pe := Peaks(ac.Freq, mp), Peaks(ac.Freq, me)
	fmt.Printf("    ELF, same network (stroke at the pole): spectrum peaks")
	for _, p := range pp {
		fmt.Printf(" %.1f", p)
	}
	fmt.Printf(" Hz; equator")
	for _, p := range pe {
		fmt.Printf(" %.1f", p)
	}
	fmt.Printf(" Hz; measured %v; AC ledger ≤ %.0e\n", MeasuredSchumann, ac.Ledger.MaxRelativeResidual)
	fmt.Printf("    electric height from the profile, σ(z) = ωε0: %.1f km at 8 Hz, %.1f km at 33 Hz (layer-3 fit: %.1f km at 8 Hz)\n",
		u.ElectricHeight(8)/1e3, u.ElectricHeight(33)/1e3, t2.He0/1e3)

	// (B) Dipole strength from the Rikitake run → conductivity → V_I.
	const qUnshielded = 2.0 // declared: global ion production without geomagnetic shielding / today's
	tr := Rikitake().RunRK4(400, 0.05)
	var rms float64
	k0 := len(tr.T) / 4
	for _, x := range tr.I1[k0:] {
		rms += x * x
	}
	rms = math.Sqrt(rms / float64(len(tr.I1)-k0))
	vOf := func(mRel float64) float64 { // quasi-static: V_I ∝ 1/σ ∝ q^(−1/2)
		q := 1 + (qUnshielded-1)*(1-math.Min(1, mRel))
		return a.VI / math.Sqrt(q)
	}
	// Twin: the layer-2 circuit with σ scaled gives the same V_I.
	sc := StandardAtmosphere()
	s0, H := sc.Calibrate()
	scale := math.Sqrt(qUnshielded)
	g, err := sc.SolveGlobal([]Column{{"g", 1, func(z float64) float64 { return scale * s0 * math.Exp(z/H) }}}, Storms{Is: sc.GlobalCurrent()})
	if err != nil {
		return err
	}
	fmt.Printf("  (B) dynamo → atmosphere: unshielded ion production ×%.1f (declared), σ ∝ √q:\n", qUnshielded)
	fmt.Printf("    dipole gone: V_I %.1f kV (circuit with σ × √%.0f: %.1f kV), relaxation time %.0f s; dipole at its mean: V_I %.0f kV\n",
		vOf(0)/1e3, qUnshielded, g.VI/1e3, eps0/(scale*s0), vOf(1)/1e3)
	var vmn, sum float64 = math.Inf(1), 0
	n := 0
	for _, x := range tr.I1[k0:] {
		v := vOf(math.Abs(x) / rms)
		vmn = math.Min(vmn, v)
		sum += v
		n++
	}
	fmt.Printf("    over a Rikitake run (|dipole| ∝ |I1|, normalized to its RMS): V_I mean %.0f kV, minimum %.0f kV at each reversal\n", sum/float64(n)/1e3, vmn/1e3)
	fmt.Println("    the two time scales are far apart (minutes for the air, kyr for a reversal): the link is quasi-static.")
	fmt.Println("  (C) one ledger: the dynamo's drive enters on the mechanical channel (core convection; geodynamo ohmic")
	fmt.Println("    losses are estimated at 0.2–0.5 TW, Christensen & Tilgner, Nature 429, 169, 2004), the storms on the")
	fmt.Println("    thermal channel (layer 2: 0.26–3.9 GW with the storm paths), and the cavity draws on lightning, part of")
	fmt.Println("    the storms' work; no layer supplies another with energy, the dynamo only modulates the air's conductivity.")
	fmt.Println()
	return nil
}
