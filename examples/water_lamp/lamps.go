// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"math/rand/v2"
)

// Constants (CODATA 2018 / SI exact).
const (
	Faraday = 96485.33212 // C/mol
	Erev    = 1.229       // V, reversible water-splitting voltage at 25 °C
)

// Lamp is the switching element across the capacitor.
type Lamp struct {
	Vb   float64 // static breakdown (strike) voltage [V]
	Vm   float64 // maintaining (glow or arc) voltage while on [V]
	Rg   float64 // resistance while on [Ω]
	Iext float64 // goes out when its current falls below this [A]
	// Seeds is the rate of free starting electrons [1/s]; 0 means a
	// deterministic strike at Vb (no statistical wait).
	Seeds float64
	// Water gap (liquid): breakdown from E·τ^(1/3) = K (MV/cm, µs) with a
	// per-flash scatter; Gap > 0 switches this on.
	Gap, K, Scatter float64
	Rleak           float64 // leakage resistance through the medium [Ω]; 0 = none
	Electrolysis    bool    // leakage current splits water
	Recovery        float64 // hold-off recovery time constant after a spark [s]; 0 = instant
	// Spark gaps conduct both ways while the arc lasts (Part 3, with the
	// coil). The arc goes out when its current has stayed below Iext for
	// Quench seconds (de-ionization, declared) or at a current zero where
	// the capacitor voltage is below the arc voltage, or — with QuenchAt > 0 —
	// at the first current zero after QuenchAt seconds (a quenched gap).
	Bipolar  bool
	Quench   float64
	QuenchAt float64
}

// Case is one lamp in its flashing-lamp circuit.
type Case struct {
	Name, Tag, Note string
	Vs, R, C        float64
	Lamp            Lamp
	Seconds         float64 // how long to simulate (Part 1)
	DriveCoil       bool
	CoilSeconds     float64
}

// PaschenAir is the breakdown voltage of air at pressure p [Torr] and gap
// d [cm]: V = B·pd / ln(A·pd / ln(1 + 1/γ)), A = 15 /(cm·Torr),
// B = 365 V/(cm·Torr), γ = 0.01 (Yu. P. Raizer, Gas Discharge Physics,
// Springer 1991, Table 4.1 and eq. 7.3; γ declared).
func PaschenAir(p, d float64) float64 {
	pd := p * d
	return 365 * pd / math.Log(15*pd/math.Log(1+1/0.01))
}

// PaschenWater is the breakdown voltage of water vapour, V = B·pd /
// ln(pd/x0), with B and x0 fixed by the measured Paschen minimum of
// Škoro et al. (Phys. Rev. E 84, 055401, 2011): 480 V at 0.6 Torr·cm. At
// the minimum V = B·e·x0 and pd = e·x0, so B = 480/0.6 V/(Torr·cm) and
// x0 = 0.6/e Torr·cm. Valid near the minimum (the measured range).
func PaschenWater(pd float64) float64 {
	B, x0 := 480/0.6, 0.6/math.E
	return B * pd / math.Log(pd/x0)
}

// cosmicSeeds is the free-electron rate in a gas volume vol [cm³] at
// pressure p [Torr] from cosmic rays at sea level: about 2 ion pairs per
// cm³ per second in air at 1 atm (declared; cosmic-ray part of the
// ground-level ionization), scaled with gas density.
func cosmicSeeds(p, vol float64) float64 { return 2 * p / 760 * vol }

// Cases are the lamps compared, each in its own flashing-lamp circuit.
func Cases() []Case {
	ne := Lamp{Vb: 76, Vm: 58, Rg: 5e3, Iext: 20e-6} // NE-2: TekWiki "Neon Lamps" (76 V, 58 V); Rg, Iext declared
	nePrimed, neDark := ne, ne
	nePrimed.Seeds = 1e4                 // krypton-85 starter: declared, fast enough to make the wait negligible
	neDark.Seeds = cosmicSeeds(60, 0.03) // fill ~60 Torr, ~0.03 cm³ (declared)
	vap := Lamp{Vb: PaschenWater(0.6), Vm: 0.76 * PaschenWater(0.6), Rg: 20e3, Iext: 20e-6}
	vapDark, vapUV := vap, vap
	vapDark.Seeds = cosmicSeeds(1.2, 2) // 1.2 Torr, 5 mm gap, 2 cm³ tube
	vapUV.Seeds = 1e3                   // ultraviolet LED on the cathode (declared photoelectron rate)
	sigma := 5.5e-6                     // S/m: pure water at 25 °C, 0.055 µS/cm
	gap, area := 1e-3, 0.2e-4
	liq := Lamp{Gap: gap, K: 0.135, Scatter: 0.12, Vm: 50, Rg: 0.5, Iext: 10e-3,
		Rleak: gap / (sigma * area), Electrolysis: true}
	dry := Lamp{Vb: PaschenAir(760, 0.1), Vm: 30, Rg: 0.5, Iext: 10e-3, Seeds: 1e3, Recovery: 1e-3}
	mist := dry
	mist.Vb = dry.Vb * 29.4 / 21.7
	return []Case{
		{Name: "neon NE-2, radioactive starter", Tag: "neon_primed", Vs: 120, R: 10e6, C: 10e-9, Lamp: nePrimed, Seconds: 2,
			Note: "reference: strike 76 V, maintain 58 V (TekWiki NE-2); krypton-85 starter rate declared", DriveCoil: true, CoilSeconds: 0.5},
		{Name: "neon NE-2, no starter, in the dark", Tag: "neon_dark", Vs: 120, R: 10e6, C: 10e-9, Lamp: neDark, Seconds: 36000,
			Note: "cosmic rays only (≈2 ion pairs/cm³/s at 1 atm, scaled to ~60 Torr, ~0.03 cm³; declared)"},
		{Name: "water-vapour tube, in the dark", Tag: "vapour_dark", Vs: 600, R: 22e6, C: 10e-9, Lamp: vapDark, Seconds: 36000,
			Note: "1.2 Torr water vapour, 5 mm gap, at the measured Paschen minimum 480 V (Škoro 2011); maintain 0.76 × strike (NE-2 ratio, declared)"},
		{Name: "water-vapour tube, ultraviolet LED starter", Tag: "vapour_uv", Vs: 600, R: 22e6, C: 10e-9, Lamp: vapUV, Seconds: 2,
			Note: "same tube; a UV LED frees electrons from the cathode (rate declared) instead of a radioactive gas", DriveCoil: true, CoilSeconds: 1.5},
		{Name: "spark gap in liquid water (pure)", Tag: "liquid", Vs: 3000, R: 1e6, C: 10e-9, Lamp: liq, Seconds: 2,
			Note: "1 mm, 0.2 cm² electrodes; E·τ^(1/3) = 0.135 ± 12 % (Stygar 2006) EXTRAPOLATED from ≤ 0.6 µs to ms; leakage 0.055 µS/cm", DriveCoil: true, CoilSeconds: 0.3},
		{Name: "air spark gap, dry", Tag: "air_dry", Vs: 10000, R: 1e6, C: 10e-9, Lamp: dry, Seconds: 2,
			Note: "1 mm at 1 atm, Paschen's law (Raizer); room-air small ions as starters and 1 ms recovery (declared)"},
		{Name: "air spark gap with water mist", Tag: "mist", Vs: 10000, R: 1e6, C: 10e-9, Lamp: mist, Seconds: 2,
			Note: "fog raises breakdown ×1.355 (Liu 2018, 5 cm rod–plane, EXTRAPOLATED to 1 mm); recovery as dry (faster recovery = H5, untested)", DriveCoil: true, CoilSeconds: 0.3},
	}
}

// BreakdownText describes the strike condition.
func (c Case) BreakdownText() string {
	if c.Lamp.Gap > 0 {
		return fmt.Sprintf("E·τ^(1/3) = %.3f MV/cm·µs^(1/3) ± %.0f %% over %.1f mm", c.Lamp.K, 100*c.Lamp.Scatter, c.Lamp.Gap*1e3)
	}
	return fmt.Sprintf("%.0f V", c.Lamp.Vb)
}

// SeedText describes the starting electrons.
func (c Case) SeedText() string {
	switch {
	case c.Lamp.Gap > 0:
		return "not needed (breakdown starts at the electrode in the liquid)"
	case c.Lamp.Seeds == 0:
		return "always present"
	default:
		return fmt.Sprintf("%.3g per second", c.Lamp.Seeds)
	}
}

// Run is a Part 1 (or Part 2) result.
type Run struct {
	Case                                 Case
	T                                    float64
	Flashes                              int
	MeanStrike, SDStrike                 float64
	MeanInterval, JitterInterval         float64
	EnergyPerFlash                       float64
	Supply, RHeat, LampEnergy, Leak      float64
	Stored, Residual, Chemical, Hydrogen float64
	TraceT, TraceV                       []float64
	StrikeTimes                          []float64 // [s]
}

// Run simulates the flashing-lamp circuit for `seconds` with seed.
func (c Case) Run(seconds float64, seed uint64) (Run, error) {
	r, err := simulate(c, nil, seconds, seed)
	return r.Run, err
}

// mean and standard deviation
func meanSD(x []float64) (float64, float64) {
	if len(x) == 0 {
		return 0, 0
	}
	m := 0.0
	for _, v := range x {
		m += v
	}
	m /= float64(len(x))
	s := 0.0
	for _, v := range x {
		s += (v - m) * (v - m)
	}
	if len(x) > 1 {
		s /= float64(len(x) - 1)
	}
	return m, math.Sqrt(s)
}

// simulate is the one engine for Parts 1 and 2: RK4 on
// y = [V, Ip, Is, Vc2, supply, R heat, leak, lamp, coil heat], with the
// energies integrated as states so the ledger closes to RK4 accuracy.
func simulate(c Case, coil *Coil, seconds float64, seed uint64) (CoilRun, error) {
	l := c.Lamp
	rng := rand.New(rand.NewPCG(seed, 0x6c616d70))
	inv := func(x float64) float64 {
		if x == 0 {
			return 0
		}
		return 1 / x
	}
	gl := inv(l.Rleak)
	const (
		iV = iota
		iIp
		iIs
		iVc
		iSup
		iR
		iLeak
		iLamp
		iCoil
		nY
	)
	var y [nY]float64
	on := false
	der := func(y [nY]float64) [nY]float64 {
		var d [nY]float64
		V := y[iV]
		iR_ := (c.Vs - V) / c.R
		ip := 0.0
		if on {
			if coil == nil {
				ip = (V - l.Vm) / l.Rg
			} else {
				ip = y[iIp]
			}
		}
		d[iV] = (iR_ - V*gl - ip) / c.C
		d[iSup] = c.Vs * iR_
		d[iR] = iR_ * iR_ * c.R
		d[iLeak] = V * V * gl
		if on {
			d[iLamp] = l.Vm*math.Abs(ip) + l.Rg*ip*ip
		}
		if coil != nil {
			// [Lp M; M Ls]·[dIp; dIs] = [V − Vm − Rg·Ip (on) ; Vc2]
			a := 0.0
			if on {
				sg := 1.0
				if y[iIp] < 0 || (y[iIp] == 0 && V < 0) {
					sg = -1
				}
				a = V - l.Vm*sg - l.Rg*y[iIp]
			}
			det := coil.Lp*coil.Ls - coil.M*coil.M
			dIp := (coil.Ls*a - coil.M*y[iVc]) / det
			dIs := (coil.Lp*y[iVc] - coil.M*a) / det
			if !on {
				dIp, dIs = 0, y[iVc]/coil.Ls
			}
			d[iIp], d[iIs] = dIp, dIs
			d[iVc] = (-y[iIs] - y[iVc]/coil.Rs) / coil.Cs
			d[iCoil] = y[iVc] * y[iVc] / coil.Rs
		}
		return d
	}
	stored := func(y [nY]float64) float64 {
		w := 0.5 * c.C * y[iV] * y[iV]
		if coil != nil {
			w += 0.5*coil.Lp*y[iIp]*y[iIp] + coil.M*y[iIp]*y[iIs] + 0.5*coil.Ls*y[iIs]*y[iIs] + 0.5*coil.Cs*y[iVc]*y[iVc]
		}
		return w
	}
	coilStored := func(y [nY]float64) float64 {
		if coil == nil {
			return 0
		}
		return 0.5*coil.Lp*y[iIp]*y[iIp] + coil.M*y[iIp]*y[iIs] + 0.5*coil.Ls*y[iIs]*y[iIs] + 0.5*coil.Cs*y[iVc]*y[iVc]
	}
	w0 := stored(y)
	tauOff := c.R * c.C
	if l.Rleak > 0 {
		tauOff = c.C / (1/c.R + gl)
	}
	dtOff := tauOff / 400
	dtOn := c.C * l.Rg / 100
	if coil != nil {
		per := 1 / coil.F1
		dtOff = math.Min(dtOff, per/60)
		dtOn = math.Min(math.Max(dtOn, 2*math.Pi*math.Sqrt(coil.Lp*c.C)/400), per/400)
		// the primary's leakage inductance with the lamp's resistance is the
		// fastest time constant while the lamp is on
		dtOn = math.Min(dtOn, 0.2*(coil.Lp-coil.M*coil.M/coil.Ls)/l.Rg)
	}
	var res CoilRun
	res.Case = c
	t := 0.0
	hazard, need := 0.0, rng.ExpFloat64()
	scatter := 1 + l.Scatter*rng.NormFloat64()
	tExt := -1.0 // time the last spark went out
	var strikes, times []float64
	// liquid gap: effective time = time spent above 63 % of the present voltage
	var histT, histV []float64
	var lampOnStart, flashE0, prevIp, lastHigh float64
	maxIp := 0.0
	nextTrace := 0.0
	traceDt := seconds / 4000
	if coil != nil {
		traceDt = math.Min(traceDt, 1/coil.F1/20)
	}
	for t < seconds {
		dt := dtOff
		if on {
			dt = dtOn
		}
		// Settled and waiting for a starting electron: the voltage no longer
		// changes, so the wait is drawn directly and the energies advance at
		// constant power.
		if !on && coil == nil && l.Seeds > 0 && l.Recovery == 0 && l.Gap == 0 {
			vInf := c.Vs / c.R / (1/c.R + gl)
			if V := y[iV]; math.Abs(vInf-V) < 1e-9*vInf && V > l.Vb {
				lam := l.Seeds * (1 - l.Vb/V)
				jump := math.Min((need-hazard)/lam, seconds-t)
				iR_ := (c.Vs - V) / c.R
				y[iSup] += c.Vs * iR_ * jump
				y[iR] += iR_ * iR_ * c.R * jump
				y[iLeak] += V * V * gl * jump
				hazard += lam * jump
				t += jump
				if t >= seconds {
					break
				}
				dt = math.Min(dtOff, 1e-12*t+1e-15) // tiny step, then the strike test below fires
			}
		}
		// Between flashes, once the coil has rung down to nothing, stop
		// following it (its last bit of energy is booked as coil heat) and
		// take the capacitor's own step.
		if !on && coil != nil {
			if w := coilStored(y); w < 1e-15 {
				y[iCoil] += w
				y[iIp], y[iIs], y[iVc] = 0, 0, 0
				dt = tauOff / 400
			}
		}
		if t+dt > seconds {
			dt = seconds - t
		}
		k1 := der(y)
		var y2, y3, y4 [nY]float64
		for i := range y {
			y2[i] = y[i] + dt/2*k1[i]
		}
		k2 := der(y2)
		for i := range y {
			y3[i] = y[i] + dt/2*k2[i]
		}
		k3 := der(y3)
		for i := range y {
			y4[i] = y[i] + dt*k3[i]
		}
		k4 := der(y4)
		for i := range y {
			y[i] += dt / 6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
		}
		t += dt
		V := y[iV]
		if t >= nextTrace {
			res.TraceT = append(res.TraceT, t)
			if coil != nil {
				res.TraceV = append(res.TraceV, y[iVc])
			} else {
				res.TraceV = append(res.TraceV, V)
			}
			nextTrace = t + traceDt
		}
		if coil != nil {
			res.CoilPeak = math.Max(res.CoilPeak, math.Abs(y[iVc]))
		}
		if on {
			ip := (V - l.Vm) / l.Rg
			if coil != nil {
				ip = y[iIp]
			}
			maxIp = math.Max(maxIp, ip)
			res.PeakIp = math.Max(res.PeakIp, math.Abs(ip))
			out := false
			switch {
			case coil != nil && l.Bipolar && l.QuenchAt > 0:
				out = t-lampOnStart >= l.QuenchAt && ip*prevIp <= 0 && prevIp != 0
			case coil != nil && l.Bipolar:
				if math.Abs(ip) >= l.Iext {
					lastHigh = t
				}
				// out once the current has stayed small for the de-ionization
				// time, or at a current zero when the capacitor can no longer
				// drive the arc's own voltage
				out = t-lastHigh >= l.Quench || (ip*prevIp <= 0 && prevIp != 0 && math.Abs(V) < l.Vm)
			default:
				out = t-lampOnStart > dt && ip < l.Iext && ip < maxIp
			}
			prevIp = ip
			if out {
				on = false
				tExt = t
				res.FlashEnergy += flashE0 - 0.5*c.C*V*V
				if coil != nil {
					// the opening gap stops the primary's small remaining
					// current; its magnetic energy goes into the spark
					w := coilStored(y)
					y[iIp] = 0
					y[iLamp] += w - coilStored(y)
				}
				hazard, need = 0, rng.ExpFloat64()
				scatter = 1 + l.Scatter*rng.NormFloat64()
				histT, histV = histT[:0], histV[:0]
			}
			continue
		}
		// strike?
		vb := l.Vb
		if l.Recovery > 0 && tExt >= 0 {
			vb *= -math.Expm1(-(t - tExt) / l.Recovery)
		}
		strike := false
		switch {
		case l.Gap > 0:
			histT, histV = append(histT, t), append(histV, V)
			// first time the voltage was ≥ 0.63·V (the rise is monotonic)
			lo := 0
			for lo < len(histV) && histV[lo] < 0.63*V {
				lo++
			}
			tau := (t - histT[min(lo, len(histT)-1)]) * 1e6 // µs
			if tau > 0 && V > 0 {
				eb := l.K * scatter / math.Cbrt(tau) * 1e8 // V/m
				strike = V/l.Gap >= eb
			}
		case l.Seeds == 0:
			strike = V >= vb
		default:
			if V > vb {
				hazard += l.Seeds * (1 - vb/V) * dt
				strike = hazard >= need
			}
		}
		if strike {
			on = true
			lampOnStart, maxIp, prevIp, lastHigh = t, 0, 0, t
			flashE0 = 0.5 * c.C * V * V
			strikes = append(strikes, V)
			times = append(times, t)
		}
	}
	r := &res.Run
	r.Case, r.T = c, seconds
	r.Flashes = len(strikes)
	r.StrikeTimes = times
	r.MeanStrike, r.SDStrike = meanSD(strikes)
	var iv []float64
	for i := 1; i < len(times); i++ {
		iv = append(iv, times[i]-times[i-1])
	}
	r.MeanInterval, r.JitterInterval = meanSD(iv)
	if r.Flashes > 0 {
		r.EnergyPerFlash = res.FlashEnergy / float64(r.Flashes)
	}
	r.Supply, r.RHeat, r.Leak, r.LampEnergy = y[iSup], y[iR], y[iLeak], y[iLamp]
	r.Stored = stored(y) - w0
	res.CoilEnergy = y[iCoil] + coilStored(y)
	r.Residual = (r.Supply - r.RHeat - r.Leak - r.LampEnergy - y[iCoil] - r.Stored) / math.Max(r.Supply, 1e-30)
	if l.Electrolysis && l.Rleak > 0 {
		// charge through the water × 1.23 V is chemical energy (only while V > 1.23 V,
		// which holds for nearly all the time here); H₂: one mole per 2 F
		q := 0.0
		for i := 1; i < len(res.TraceT); i++ {
			v := 0.5 * (res.TraceV[i] + res.TraceV[i-1])
			if v > Erev {
				q += v / l.Rleak * (res.TraceT[i] - res.TraceT[i-1])
			}
		}
		r.Chemical = q * Erev
		r.Hydrogen = q / (2 * Faraday)
	}
	r.TraceT, r.TraceV = res.TraceT, res.TraceV
	res.Supply, res.RHeat, res.LampEnergy, res.Leak, res.Stored, res.Residual = r.Supply, r.RHeat, r.LampEnergy, r.Leak, r.Stored, r.Residual
	res.Flashes = r.Flashes
	return res, nil
}
