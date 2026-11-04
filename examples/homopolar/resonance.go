// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Resonance and harmonics in the stacked dual-motor machine (dual.go), and
// what a load needs to lock a resonance in.
//
// A homopolar loop on its own is first order: resistance and inductance,
// one kind of stored energy (magnetic). Oscillating its speed only makes the
// speed EMF M·ω(t)·I a time-varying resistance: harmonics, no resonance. A
// resonance needs a second store. The machine already has one: the shafts'
// kinetic energy. In Kron's analogy an inertia is a capacitor on the
// mechanical channel, and the speed coupling joins it to the coils' L, so
// the self-excited machine is an L–J oscillator. About its steady state
// (symmetric mode, I1 = I2 = I, ω1 = ω2 = ω* = R/M_s) the loop's net
// resistance −R + M_s·ω* is zero, and
//
//	L_s·δİ = M_s·I·δω,   J·δω̇ = −2·M_s·I·δI − b·δω
//	⇒ f_em = (1/2π)·√(2·M_s²·I²/(L_s·J)) = (1/2π)·√(2·M_s·G/(L_s·J)),  Q = J·2πf_em/b,
//
// with M_s = M·(1 + KS) and L_s = L·(1 + KM) the stacked couplings. The load
// sets the steady speed (ω* = R/M_s) and takes its share of the power, but to
// first order it does not damp this resonance: only friction does. (The
// side-by-side machine has the same oscillation: Rikitake's neutral period.)
//
// The other road is the RLC one: a capacitor across each load. The load and
// the capacitor make a parallel tank fed by the coils; DC still reaches the
// load, and the capacitor's ESR stays out of the DC path. With the speed held,
// the loop's characteristic equation is
//
//	L·C·s² + (L/R_L + r·C)·s + (1 + r/R_L) = 0,   r = R_int − M_s·ω,
//
// which goes to DC at r = −R_L but oscillates first once C > L/R_L² (the
// tank's Q = R_L/√(L/C) > 1). With free shafts the inertia joins the tank and
// the machine locks into AC at a much smaller C (lockInC).

package main

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"

	"github.com/fdamador/kronspice/analysis"
)

// stackedLoad is the stacked machine with a load in each loop, started at
// its symmetric steady state.
func stackedLoad(rLoad, windage float64) Dual {
	m := DefaultDual()
	m.Stacked, m.RLoad, m.Windage = true, rLoad, windage
	_, _, ks, r := m.loop()
	ms := (1 + ks) * m.Coupling()
	w := r / ms                                    // ω* (friction neglected in the start point)
	i := math.Sqrt((m.Torque[0] - windage*w) / ms) // G = b·ω + M_s·I²
	m.I0, m.W0 = [2]float64{i, i}, [2]float64{w, w}
	return m
}

// EMResonance returns the electromechanical resonance f_em [Hz] and Q of
// the symmetric mode about the steady state.
func (m Dual) EMResonance() (f, q float64) {
	l11, l12, ks, _ := m.loop()
	ms, ls := (1+ks)*m.Coupling(), l11+l12
	i := m.I0[0]
	w := math.Sqrt(2 * ms * ms * i * i / (ls * m.Inertia()))
	return w / (2 * math.Pi), m.Inertia() * w / m.Windage
}

// amplitude is the complex amplitude of x at frequency f over the samples
// with t ≥ t0 (a whole number of periods is taken).
func amplitude(t, x []float64, t0, f float64) complex128 {
	period := 1 / f
	k0 := 0
	for k0 < len(t) && t[k0] < t0 {
		k0++
	}
	span := t[len(t)-1] - t[k0]
	tEnd := t[k0] + math.Floor(span/period)*period
	var s complex128
	var n int
	for k := k0; k < len(t) && t[k] <= tEnd; k++ {
		s += complex(x[k], 0) * cmplx.Exp(complex(0, -2*math.Pi*f*t[k]))
		n++
	}
	return 2 * s / complex(float64(n), 0)
}

// mean is the average of x over t ≥ t0.
func mean(t, x []float64, t0 float64) float64 {
	var s float64
	var n int
	for k := range t {
		if t[k] >= t0 {
			s += x[k]
			n++
		}
	}
	return s / float64(n)
}

// rmsGrowth is the growth rate of the RMS of x between the two halves of
// t ≥ t0.
func rmsGrowth(t, x []float64, t0 float64) float64 {
	k0 := 0
	for k0 < len(t) && t[k0] < t0 {
		k0++
	}
	mid := (k0 + len(t)) / 2
	rms := func(a, b int) float64 {
		var s float64
		for k := a; k < b; k++ {
			s += x[k] * x[k]
		}
		return math.Sqrt(s / float64(b-a))
	}
	return math.Log(rms(mid, len(t))/rms(k0, mid)) / (t[(mid+len(t))/2] - t[(k0+mid)/2])
}

func runResonance() error {
	const rLoad, windage = 50e-6, 5e-3 // declared: a 50 µΩ busbar load; bearing + air friction
	const esrCell = 0.29e-3            // declared: typical ESR of one 3000 F, 2.7 V supercapacitor cell
	m := stackedLoad(rLoad, windage)
	_, _, ks, r := m.loop()
	ms := (1 + ks) * m.Coupling()
	l11, l12, _, _ := m.loop()
	ls := l11 + l12
	fmt.Println("Resonance and harmonics, stacked dual motor with a load (50 µΩ in each loop, friction 5e-3 N·m·s, declared)")
	w, i := m.W0[0], m.I0[0]
	fmt.Printf("  steady state at 20 N·m: ω* = R/M_s = %.0f rad/s (rim %.0f m/s, limit %.0f), I = %.2f kA,\n", w, w*m.Disk.Ro, burst(), i/1e3)
	fmt.Printf("  load %.2f kW per loop of %.2f kW from each motor (%.0f %%)\n", rLoad*i*i/1e3, m.Torque[0]*w/1e3, 100*rLoad*i*i/(m.Torque[0]*w))

	// A. Electromechanical resonance (no added capacitor).
	fem, q := m.EMResonance()
	fmt.Printf("  A. no capacitor: the shafts' inertia is the second store (L–J). f_em = √(2·M_s·G/(L_s·J))/2π = %.3f Hz, Q = %.0f\n", fem, q)
	fmt.Println("     torque oscillating ±0.5 N·m (±2.5 %) at f, after 40 s: speed and current oscillation amplitudes")
	fmt.Println("       f / f_em   δω [rad/s]   δI [A]")
	var fPeak, aPeak float64
	for _, x := range []float64{0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.5, 2} {
		q := m
		q.TorqueAmp, q.TorqueFreq = 0.5, x*fem
		tr := q.RunRK4(60, 2e-3)
		aw := cmplx.Abs(amplitude(tr.T, tr.W1, 40, x*fem))
		ai := cmplx.Abs(amplitude(tr.T, tr.I1, 40, x*fem))
		if aw > aPeak {
			aPeak, fPeak = aw, x*fem
		}
		fmt.Printf("       %-6.2f     %8.2f     %8.1f\n", x, aw, ai)
	}
	// static response for comparison: δω at f → 0 ~ the torque swing over the speed stiffness
	fmt.Printf("     peak at %.3f Hz (f_em %.3f): the machine resonates on its own L and J\n", fPeak, fem)
	// Twin: the circuit at the peak.
	{
		q := m
		q.TorqueAmp, q.TorqueFreq = 0.5, fPeak
		tr, err := q.Run(60, 2e-3)
		if err != nil {
			return err
		}
		ac := cmplx.Abs(amplitude(tr.T, tr.W1, 40, fPeak))
		fmt.Printf("     circuit (sine torque sources) at the peak: δω %.2f rad/s against RK4 %.2f; ledger mismatch %.0e of the drive work\n",
			ac, aPeak, tr.Mismatch/tr.In)
	}
	// Harmonics at resonance, larger drive.
	fmt.Println("     harmonics of the loop current at resonance (fraction of the fundamental):")
	fmt.Println("       torque swing   n = 2     n = 3     n = 4     n = 5     reversals")
	for _, amp := range []float64{0.5, 2, 5} {
		q := m
		q.TorqueAmp, q.TorqueFreq = amp, fem
		tr := q.RunRK4(80, 1e-3)
		a1 := cmplx.Abs(amplitude(tr.T, tr.I1, 50, fem))
		fmt.Printf("       ±%-4.1f N·m    ", amp)
		for n := 2; n <= 5; n++ {
			fmt.Printf("%.4f    ", cmplx.Abs(amplitude(tr.T, tr.I1, 50, float64(n)*fem))/a1)
		}
		fmt.Printf("%d\n", len(reversals(tr.T, tr.I1, 0, 50)))
	}
	fmt.Println("     the speed EMF (M·ω·I) and the torque (M·I²) are products: driven hard, the oscillation makes")
	fmt.Println("     harmonics of f_em; the load does not damp it (the steady speed sits where R = M_s·ω).")

	// B. A capacitor across the load: a parallel tank fed by the coils.
	fmt.Println("  B. a capacitor across each load (parallel): the load and C make a tank fed by the coils' L.")
	fmt.Println("     DC still reaches the load, and the capacitor's ESR is outside the DC path. Holding the speed,")
	fmt.Println("     the loop obeys L·C·s² + (L/R_L + r·C)·s + (1 + r/R_L) = 0 with r = R_int − M_s·ω: it settles")
	fmt.Println("     to DC when r reaches −R_L, but oscillates first, at the tank, once L/R_L + r·C < 0 can come")
	fmt.Println("     earlier: when C > L/R_L², i.e. tank Q = R_L/√(L/C) > 1.")
	rint := r - rLoad
	cCrit := ls / (rLoad * rLoad)
	fmt.Printf("     here L_s = %.1f µH, R_L = %.0f µΩ: C_crit = L/R_L² = %.0f F (%.1f cells of 3000 F)\n", ls*1e6, rLoad*1e6, cCrit, cCrit/3000)
	fmt.Println("       cells   C [F]    ESR [mΩ]  tank f0 [Hz]  Q = R_L/√(L/C)  AC onset ω [rad/s] at f [Hz]  DC onset ω  first")
	for _, n := range []float64{1, 2, 4, 8} {
		c := 3000 * n
		zc := math.Sqrt(ls / c)
		f0 := 1 / (2 * math.Pi * math.Sqrt(ls*c))
		wDC := (rint + rLoad) / ms
		wAC := (rint + ls/(rLoad*c)) / ms
		first, fH := "DC", "  —  "
		if wAC < wDC {
			first = "AC (tank)"
			fH = fmt.Sprintf("%.3f", f0*math.Sqrt(1-cCrit/c)) // f0·√(1 − L/(R_L²·C)) at the onset
		}
		fmt.Printf("       %-5.0f   %-6.0f   %.4f    %.3f         %.2f            %6.0f at %s             %6.0f      %s\n",
			n, c, esrCell*1e3/n, f0, rLoad/zc, wAC, fH, wDC, first)
	}
	fmt.Println("     (onsets without the ESR; the runs below include it.)")

	// Torque-driven with the tank: the shafts are free, so the electromechanical
	// mode of A joins the tank. Every run starts 1 % off the DC state, so an
	// unstable DC state is left and a stable one returned to.
	fmt.Println("     torque-driven at 20 N·m (free shafts), from 1 % off the DC state, 80 s (ESR of whole 3000 F cells):")
	fmt.Println("       C [F]    tank f0 [Hz]  Q      outcome  loop DC [kA]  loop AC [kA] at f [Hz]  harmonics 2, 3 [%]  load [kW]")
	for _, c := range []float64{30, 300, 1000, 3000, 12000, 24000} {
		q := tankCase(m, c, esrCell)
		tr := q.RunRK4(80, 2e-3)
		dc := mean(tr.T, tr.I1, 50)
		fBest, aBest := 0.0, 0.0
		for f := 0.05; f <= 3.0; f += 0.01 {
			if a := cmplx.Abs(amplitude(tr.T, tr.I1, 50, f)); a > aBest {
				aBest, fBest = a, f
			}
		}
		out, h := "DC", "    —     "
		if aBest > 0.1*math.Abs(dc) || aBest > 1000 {
			out = "AC"
			h = fmt.Sprintf("%4.1f, %4.1f", 100*cmplx.Abs(amplitude(tr.T, tr.I1, 50, 2*fBest))/aBest, 100*cmplx.Abs(amplitude(tr.T, tr.I1, 50, 3*fBest))/aBest)
		}
		pl, err := q.loadPower(80, 2e-3, 50)
		if err != nil {
			return err
		}
		fmt.Printf("       %-7.0f  %6.3f        %5.2f  %-7s  %8.2f      %8.2f at %.2f        %s          %5.2f\n",
			c, 1/(2*math.Pi*math.Sqrt(ls*c)), rLoad/math.Sqrt(ls/c), out, dc/1e3, aBest/1e3, fBest, h, pl/1e3)
	}
	cLock := lockInC(m, esrCell)
	fmt.Printf("     with free shafts the DC state gives way to AC at C ≈ %.0f F (bisection), far below the held-speed\n", cLock)
	fmt.Printf("     rule C > L/R_L² = %.0f F: the shafts' inertia joins the tank and supplies the gain.\n", cCrit)
	{
		q := tankCase(m, 3000, esrCell)
		tr := q.RunRK4(40, 2e-3)
		ct, err := q.Run(40, 2e-3)
		if err != nil {
			return err
		}
		var d float64
		for k := 1; k < len(ct.T); k++ {
			j := int(math.Round(ct.T[k] / 2e-3))
			d = math.Max(d, math.Abs(ct.I1[k]-tr.I1[j]))
		}
		fmt.Printf("     circuit vs RK4 with the tank (3000 F, 40 s): |ΔI1| ≤ %.1f A of ~17 kA; ledger mismatch %.0e of the drive work\n", d, ct.Mismatch/ct.In)
	}

	// Oscillating speed with the tank.
	fmt.Println("     speed drives oscillating ±30 % at f_p, mean speed 0.9 of the first onset (8 cells, ESR included):")
	q8 := m
	q8.CParallel, q8.ESR = 24000, esrCell/8
	f8 := 1 / (2 * math.Pi * math.Sqrt(ls*q8.CParallel))
	wOn := math.Min((rint+rLoad)/ms, (rint+ls/(rLoad*q8.CParallel))/ms)
	fmt.Println("       f_p / f0   no capacitor   capacitor across the load")
	for _, x := range []float64{0.5, 1, 1.5, 2, 2.5, 3} {
		var g [2]float64
		for j, c := range []float64{0, q8.CParallel} {
			d := DefaultDual()
			d.Stacked, d.RLoad, d.CParallel, d.ESR = true, rLoad, c, q8.ESR
			d.Speed = &[2]float64{0.9 * wOn, 0.9 * wOn}
			d.SpeedAmp, d.SpeedFreq = 0.3, x*f8
			d.I0 = [2]float64{100, 100}
			tr := d.RunRK4(60, 1e-3)
			g[j] = rmsGrowth(tr.T, tr.I1, 15)
		}
		fmt.Printf("       %-5.2f      %+8.3f       %+8.3f\n", x, g[0], g[1])
	}
	fmt.Println("     with the tank, a slow speed oscillation (f_p ≲ f0) lifts the field, because the speed spends part")
	fmt.Println("     of each cycle above the AC onset; faster ones average out. The oscillation that locks in is the")
	fmt.Println("     machine's own: speed swings are not needed to start it once C is past the lock-in value.")
	fmt.Println()
	fmt.Println("  C. locking a resonance in, by the RLC rule (two stores, one loss), with the load:")
	fmt.Println("     L: already there (the coils, 13.5 µH per loop with the stacked mutual).")
	fmt.Println("     C: across the load. It keeps DC in the load and makes a parallel tank. With the shafts held at a")
	fmt.Printf("        speed the tank locks in when Q = R_L/√(L/C) > 1, i.e. C > L/R_L² = %.0f F; with free shafts the\n", cCrit)
	fmt.Printf("        inertia joins in and it locks in at C ≈ %.0f F (Q ≈ %.2f). A larger load resistance needs less C\n",
		cLock, rLoad/math.Sqrt(ls/cLock))
	fmt.Println("        (C ∝ 1/R_L²) but raises the DC threshold speed.")
	fmt.Println("        Locked in, the machine is an AC generator at 0.14–0.27 Hz, below both the tank and f_em: the")
	fmt.Println("        loop current reverses every half cycle, with a 3rd harmonic of 37–78 %. It gives the load less")
	fmt.Println("        (2.3 kW on DC; 2.0, 1.8, 1.0, 0.67 kW at 1000, 3000, 12 000, 24 000 F from the same 20 N·m):")
	fmt.Println("        much of the AC circulates between the coils and the capacitor, heating the ESR and the coils.")
	fmt.Println("     both: an inductor added in series lowers the frequency and raises √(L/C), so it needs more C;")
	fmt.Println("        added L is only worth it to place f0, never to lock the resonance in.")
	fmt.Println("     Without a capacitor, the shafts' inertia already gives the resonance of A.")
	fmt.Println("  Every watt in the load, the coils and the oscillation is paid by the motors or the speed drives (ledger).")
	fmt.Println()
	return nil
}

// writeResonanceData writes the CSVs for viz/plot_resonance.py.
func writeResonanceData(dir string) error {
	m := stackedLoad(50e-6, 5e-3)
	fem, _ := m.EMResonance()
	var rows [][]float64
	for x := 0.6; x <= 1.4001; x += 0.02 {
		q := m
		q.TorqueAmp, q.TorqueFreq = 0.5, x*fem
		tr := q.RunRK4(60, 2e-3)
		rows = append(rows, []float64{x * fem, cmplx.Abs(amplitude(tr.T, tr.W1, 40, x*fem)), cmplx.Abs(amplitude(tr.T, tr.I1, 40, x*fem))})
	}
	if err := writeResCSV(dir, "res_curve.csv", "f,dw,di", rows); err != nil {
		return err
	}
	q := m
	q.TorqueAmp, q.TorqueFreq = 5, fem
	tr := q.RunRK4(80, 1e-3)
	rows = nil
	for k := range tr.T {
		if tr.T[k] >= 70 {
			rows = append(rows, []float64{tr.T[k], tr.I1[k], tr.W1[k]})
		}
	}
	if err := writeResCSV(dir, "res_wave.csv", "t,i1,w1", rows); err != nil {
		return err
	}
	// The lock-in: the loop current with 300 F (DC), 1000 F and 12 000 F (AC)
	// across the load, torque-driven at 20 N·m, from 1 % off the DC state.
	var tr3 [3]DualTrace
	for j, c := range []float64{300, 1000, 12000} {
		tr3[j] = tankCase(m, c, 0.29e-3).RunRK4(60, 2e-3)
	}
	rows = nil
	for k := range tr3[0].T {
		if k%5 == 0 {
			rows = append(rows, []float64{tr3[0].T[k], tr3[0].I1[k], tr3[1].I1[k], tr3[2].I1[k]})
		}
	}
	return writeResCSV(dir, "res_tank.csv", "t,i1_300F,i1_1000F,i1_12000F", rows)
}

func writeResCSV(dir, name, header string, rows [][]float64) error {
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	f, err := os.Create(filepath.Join(dir, name))
	if err != nil {
		return err
	}
	defer f.Close()
	fmt.Fprintln(f, header)
	for _, r := range rows {
		for i, v := range r {
			if i > 0 {
				fmt.Fprint(f, ",")
			}
			fmt.Fprintf(f, "%.8g", v)
		}
		fmt.Fprintln(f)
	}
	return nil
}

// loadPower runs the circuit and returns the mean power in loop 1's load
// resistor over t ≥ t0 (with a capacitor across the load, the capacitor
// carries part of the AC, so the load's own current is used).
func (m Dual) loadPower(tEnd, dt, t0 float64) (float64, error) {
	c, err := m.Circuit()
	if err != nil {
		return 0, err
	}
	wf, err := analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: tEnd, MaxStep: dt, UIC: true}, analysis.TranOptions{})
	if err != nil {
		return 0, err
	}
	name := "R1"
	if m.CParallel > 0 && m.RLoad > 0 {
		name = "RL1"
	}
	i, _ := wf.ElementCurrent(name) // without a capacitor the load is lumped into R1, which carries the loop current
	var s float64
	var n int
	for k := range wf.Time {
		if wf.Time[k] >= t0 {
			s += m.RLoad * i[k] * i[k]
			n++
		}
	}
	return s / float64(n), nil
}

// tankCase is m with a capacitor c across each load, made of 3000 F cells
// (at least one) of the given ESR each, started 1 % off the DC state.
func tankCase(m Dual, c, esrCell float64) Dual {
	q := m
	q.CParallel, q.ESR = c, esrCell*3000/math.Max(c, 3000)
	q.I0[0], q.I0[1] = 1.01*q.I0[0], 1.01*q.I0[1]
	return q
}

// lockIn reports whether the torque-driven machine with c across each load
// leaves its DC state for an AC oscillation.
func lockIn(m Dual, c, esrCell float64) bool {
	tr := tankCase(m, c, esrCell).RunRK4(80, 2e-3)
	dc := mean(tr.T, tr.I1, 50)
	var hi float64
	for k := range tr.T {
		if tr.T[k] >= 50 {
			hi = math.Max(hi, math.Abs(tr.I1[k]-dc))
		}
	}
	return hi > 0.1*math.Abs(m.I0[0])
}

// lockInC is the capacitance at which the DC state gives way to AC, by
// bisection on a log scale between 30 F and 30 000 F.
func lockInC(m Dual, esrCell float64) float64 {
	lo, hi := 30.0, 30000.0
	for it := 0; it < 12; it++ {
		mid := math.Sqrt(lo * hi)
		if lockIn(m, mid, esrCell) {
			hi = mid
		} else {
			lo = mid
		}
	}
	return math.Sqrt(lo * hi)
}

// pumpBest scans the speed-oscillation frequency f_p (0.25–4 f0) at swing h
// for a capacitor c across each load (ESR = τ/C, τ = 0.87 s as the 3000 F
// cell; τ = 0: ideal) and returns the largest growth rate and its f_p/f0.
func pumpBest(c, h, tau float64) (best, at float64) {
	m := stackedLoad(50e-6, 5e-3)
	l11, l12, _, _ := m.loop()
	f0 := 1 / (2 * math.Pi * math.Sqrt((l11+l12)*c))
	best = math.Inf(-1)
	for x := 0.25; x <= 4.01; x += 0.25 {
		d := DefaultDual()
		d.Stacked, d.RLoad, d.CParallel, d.ESR = true, 50e-6, c, tau/c
		wb := 0.95 * m.W0[0]
		d.Speed = &[2]float64{wb, wb}
		d.SpeedAmp, d.SpeedFreq = h, x*f0
		d.I0 = [2]float64{100, 100}
		step := math.Min(math.Min(1e-3, 0.02/(x*f0)), 50e-6*c/5)
		tr := d.RunRK4(math.Max(20, 40/f0), step)
		if g := rmsGrowth(tr.T, tr.I1, tr.T[len(tr.T)-1]/4); g > best {
			best, at = g, x
		}
	}
	return
}

// runFastPump answers: can a faster speed oscillation let the capacitor
// shrink, and how high can the resonance go?
func runFastPump() error {
	m := stackedLoad(50e-6, 5e-3)
	l11, l12, ks, r := m.loop()
	ls, ms := l11+l12, (1+ks)*m.Coupling()
	rint := r - m.RLoad
	wBurst := burst() / m.Disk.Ro
	wb := 0.95 * m.W0[0]
	hMax := wBurst/wb - 1
	fmt.Println("  D. a faster speed oscillation and a smaller capacitor across the load")
	fmt.Printf("     speed held at %.0f rad/s (0.95 ω*) and swung ±h at f_p; the rim limit allows h ≤ %.2f\n", wb, hMax)
	fmt.Println("       C [F]   tank f0 [Hz]  Q      best growth [1/s] at f_p/f0: h = 0.16 (rim-safe)   h = 0.5 (past the rim)")
	for _, c := range []float64{3000, 1000, 300, 100, 30} {
		g1, x1 := pumpBest(c, 0.16, 0.87)
		g2, x2 := pumpBest(c, 0.5, 0.87)
		fmt.Printf("       %-6.0f  %6.2f        %5.2f  %+.2f at %.2f                     %+.2f at %.2f\n",
			c, 1/(2*math.Pi*math.Sqrt(ls*c)), m.RLoad/math.Sqrt(ls/c), g1, x1, g2, x2)
	}
	fmt.Println("     (ESR = 0.87 s / C, the time constant of a 3000 F cell.) Below ~1000 F nothing grows, at any f_p:")
	fmt.Println("     a smaller C raises f0 but the load damps the tank at 1/(R_L·C), faster than f0 rises (Q = R_L·√(C/L)")
	fmt.Println("     falls as √C). Where it grows, the best pump is slow (f_p ≈ f0/4): the speed spends part of each")
	fmt.Println("     cycle above the onset, a quasi-static gain, not a resonance. Faster rpm swings do not buy a smaller C.")

	fmt.Println("     how high the resonance can go:")
	// Electrical: tank Q ≥ q needs R_L ≥ q·2πf0·L; R_L is capped by the burst-limited onset speed.
	rlMax := ms*wBurst - rint
	for _, q := range []float64{1, 0.29} {
		f := rlMax / (q * ls) / (2 * math.Pi)
		fmt.Printf("       tank, Q ≥ %.2f: R_L ≤ M_s·ω_burst − R_int = %.0f µΩ ⇒ f0 ≤ R_L/(2π·Q·L) = %.2f Hz (C ≥ %.0f F)\n",
			q, rlMax*1e6, f, 1/(math.Pow(2*math.Pi*f, 2)*ls))
	}
	fmt.Println("         (Q ≥ 1: speed held; 0.29: the free-shaft lock-in found in B)")
	// Mechanical: swinging Δω at f_p takes J·Δω·2πf_p of drive torque.
	dw := hMax * wb
	fmt.Printf("       drives: a rim-safe swing ±%.0f rad/s at f_p takes J·Δω·2πf_p of torque per shaft:", dw)
	for _, f := range []float64{1, 10, 100} {
		fmt.Printf("  %g Hz %.0f N·m", f, m.Inertia()*dw*2*math.Pi*f)
	}
	fmt.Println()
	// Electromechanical: f_em = √(2·M_s·G/(L_s·J))/2π.
	fem, _ := m.EMResonance()
	fmt.Printf("       no capacitor: f_em = √(2·M_s·G/(L_s·J))/2π = %.3f Hz at 20 N·m; ∝ √G: %.2f Hz at 200 N·m, %.2f Hz at 2000 N·m\n",
		fem, fem*math.Sqrt(10), fem*math.Sqrt(100))
	fmt.Println("       (the drive power G·ω* grows with G: 19 kW, 190 kW, 1.9 MW per shaft), and ∝ 1/√J: the disk alone")
	jd := m.Inertia() - m.JMotor
	fmt.Printf("       (no motor rotor) gives %.3f Hz.\n", fem*math.Sqrt(m.Inertia()/jd))
	fmt.Println("     So this copper machine resonates up to about 1–3 Hz electrically and a few Hz mechanically; higher")
	fmt.Println("     needs a lower-L, higher-R_L design whose onset speed stays under the rim limit (stronger coupling M).")
	fmt.Println()
	return nil
}

// runSchumann: the shafts oscillated at 7.8 Hz, the first Schumann frequency.
func runSchumann() error {
	const fs = 7.8
	m := stackedLoad(50e-6, 5e-3)
	l11, l12, ks, _ := m.loop()
	ls, ms := l11+l12, (1+ks)*m.Coupling()
	wb := m.W0[0]
	dwMax := burst()/m.Disk.Ro - wb
	fmt.Println("  E. the shafts oscillated at 7.8 Hz (the first Schumann frequency)")
	tq := m.Inertia() * dwMax * 2 * math.Pi * fs
	fmt.Printf("     a rim-safe swing (±%.0f rad/s about ω* = %.0f) at 7.8 Hz needs ±%.0f N·m per shaft (J·Δω·2πf)\n", dwMax, wb, tq)
	q := m
	q.TorqueAmp, q.TorqueFreq = tq, fs
	tr := q.RunRK4(10, 2e-4)
	aw := cmplx.Abs(amplitude(tr.T, tr.W1, 5, fs))
	ai := cmplx.Abs(amplitude(tr.T, tr.I1, 5, fs))
	fmt.Printf("     torque-driven ±%.0f N·m at 7.8 Hz: speed swing ±%.0f rad/s, loop current ±%.0f A at 7.8 Hz on %.2f kA DC\n",
		tq, aw, ai, mean(tr.T, tr.I1, 5)/1e3)
	fmt.Printf("     (forced, not resonant: δI ≈ M_s·I·δω/(2πf·L_s) = %.0f A); harmonics 2f %.1f %%, 3f %.1f %%\n",
		ms*m.I0[0]*aw/(2*math.Pi*fs*ls),
		100*cmplx.Abs(amplitude(tr.T, tr.I1, 5, 2*fs))/ai, 100*cmplx.Abs(amplitude(tr.T, tr.I1, 5, 3*fs))/ai)
	fmt.Printf("     7.8 Hz power into the load: %.0f W (½·R_L·δI²)\n", 0.5*m.RLoad*ai*ai)
	ct, err := q.Run(10, 2e-4)
	if err != nil {
		return err
	}
	fmt.Printf("     circuit: ±%.0f A at 7.8 Hz (RK4 ±%.0f); ledger mismatch %.0e of the drive work\n",
		cmplx.Abs(amplitude(ct.T, ct.I1, 5, fs)), ai, ct.Mismatch/ct.In)
	fmt.Println("     resonant at 7.8 Hz instead? a tank needs C = 1/((2π·7.8)²·L) = 31 F, Q = 0.07 (overdamped, D);")
	fmt.Printf("     the machine's own f_em = 7.8 Hz needs G = %.0f N·m (%.1f MW per shaft): neither fits this machine.\n",
		20*math.Pow(fs/0.587, 2), 20*math.Pow(fs/0.587, 2)*wb/1e6)

	// Coupling to the Earth–ionosphere cavity: the machine is a small magnetic dipole.
	area := math.Pi * m.Rc * m.Rc
	mdip := 2 * m.N * ai * area // two coils, N turns, the 7.8 Hz current
	w := 2 * math.Pi * fs
	const mu0, c0 = 1.25663706212e-6, 299792458.0
	prad := mu0 * mdip * mdip * math.Pow(w, 4) / (12 * math.Pi * math.Pow(c0, 3))
	fmt.Printf("     as a source for the Schumann cavity: the two coils make a magnetic dipole of %.0f A·m² at 7.8 Hz;\n", mdip)
	fmt.Printf("     its field falls as μ0·m/(4π·r³): %.2g T at 10 m, %.2g T at 1 km, against a natural Schumann field of\n",
		1e-7*mdip/1e3, 1e-7*mdip/1e9)
	fmt.Printf("     about 1 pT (Sentman 1995): above it only within %.0f m of the machine. It radiates μ0·m²·ω⁴/(12π·c³) =\n",
		math.Cbrt(1e-7*mdip/1e-12))
	fmt.Printf("     %.1g W at a wavelength of %.0f km. The cavity is fed by ~50 lightning strokes a second, each\n", prad, c0/fs/1e3)
	fmt.Println("     a vertical current of kiloamperes over kilometres; this machine is a local 7.8 Hz magnet, not a")
	fmt.Println("     Schumann source. (examples/earth_machine, layer 3, has the cavity's response per ampere of stroke.)")
	fmt.Println()
	return nil
}
