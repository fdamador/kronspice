// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// A parallel RLC output: every load element across the loop's terminals,
// R_L ∥ C ∥ L_p, with L_p the primary of a transformer whose secondary feeds
// an AC load R_ac. The rpm follows 1768·(1 + h·sin 2π·7.8t) (rpm_schumann.go).
//
// At DC the primary is a short: it carries the DC current, R_L sees none,
// and the loop's DC resistance becomes R_int + R_p (the primary's winding).
// At 7.8 Hz the machine is a source of EMF M_s·I·ω0·h behind
// Z_s = R_s + j·ω·L_s (its coils). The most AC power reaches the outputs when
// they present the conjugate, R_s − j·ω·L_s: a parallel conductance
// G = R_s/|Z_s|² (split equally here between R_L and the secondary's load as
// seen through the transformer, R_ac/n²) and a capacitive susceptance
// B = ω·L_s/|Z_s|², made by C against L_p: ω·C − 1/(ω·L_p) = B. The tank's own
// Q = R_eff/(ω·L_p) is chosen (3); the coils and the output then form a
// series RLC at 7.8 Hz with Q ≈ ω·L_s/(2·R_s), underdamped. Matched, each
// loop can give EMF²/(8·R_s) at 7.8 Hz.

package main

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/output"
)

// OutputDesign is one loop's matched R ∥ C ∥ L_p output. The machine uses
// ONE tank between the two machines' coils (Dual.Output = OutSeries): these
// values in series twice, so 2·R_L, C/2, 2·L_p and 2·R_ac, with the same turns
// ratio. Shaft 2 counter-rotates; the shafts are the common, and the
// secondary's common end is tied to it.
type OutputDesign struct {
	RL, LP, C, LS, RAC, N, Reff float64
	Rs, Xs                      float64 // the source seen at 7.8 Hz
}

// outputMachine returns the iron machine whose onset, with the primary
// carrying the DC, is 1768 rpm, and its matched output design.
func outputMachine(rPrim, qTank, rAC, kT float64) (Dual, OutputDesign) {
	w0 := rpmBase * 2 * math.Pi / 60
	m := DefaultDual()
	m.Stacked, m.Windage, m.RLoad = true, 5e-3, 1 // placeholder load, replaced below
	_, _, ks, _ := m.loop()
	rint := m.Disk.DiskResistance() + 2*m.Disk.RShaft
	_, rc := m.Coil()
	rint += rc
	rdc := rint + rPrim
	m.IronGain = rdc / ((1 + ks) * m.Coupling() * w0)
	l11, l12, _, _ := m.loop()
	om := 2 * math.Pi * fSchu
	rs, xs := rdc, om*(l11+l12)
	z2 := rs*rs + xs*xs
	g, bsus := rs/z2, xs/z2
	reff := 1 / g
	d := OutputDesign{Reff: reff, Rs: rs, Xs: xs}
	d.RL = 2 * reff // half the conductance in R_L, half through the transformer
	d.LP = reff / (qTank * om)
	d.C = (bsus + 1/(om*d.LP)) / om
	rRefl := 2 * reff
	d.RAC = rAC
	d.N = math.Sqrt(rAC / rRefl)
	d.LS = d.N * d.N * d.LP
	m.RLoad, m.CParallel, m.LPrim, m.RPrim, m.LSec, m.KT, m.RAC = d.RL, d.C, d.LP, rPrim, d.LS, kT, rAC
	m.Output = OutSeries // one tank between the coils: the per-loop values above, two in series
	return m, d
}

func runOutputTank() error {
	const rPrim, qTank, rAC, kT = 10e-6, 3.0, 10.0, 0.98 // declared: primary winding, tank Q, AC load, coupling
	m, d := outputMachine(rPrim, qTank, rAC, kT)
	w0 := rpmBase * 2 * math.Pi / 60
	zt, it := m.TankScale()
	fmt.Println("  G. parallel RLC output: ONE R_L ∥ C ∥ L_p between the two machines' coils, L_p the primary of a transformer")
	fmt.Println("     to an AC load; one current runs shaft → disk 1 → coil 1 → tank → coil 2 → disk 2 → shaft (shaft 2 counter-rotates);")
	fmt.Println("     the two shafts are the common (ground), and the secondary's common end is tied to it")
	fmt.Printf("     iron machine (gain %.2f) with the DC through the primary (R_p = %.0f µΩ, declared): onset %.0f rpm\n",
		m.IronGain, rPrim*1e6, rpmBase)
	fmt.Printf("     source at 7.8 Hz: R_s = %.3f mΩ, X_s = ω·L_s = %.3f mΩ; matched output: G = R_s/|Z_s|² → R_eff = %.2f mΩ\n",
		d.Rs*1e3, d.Xs*1e3, d.Reff*1e3)
	fmt.Printf("     per loop: R_L = %.2f mΩ ∥ C = %.2f F ∥ L_p = %.1f µH (tank f = %.2f Hz, Q = %.1f); transformer 1:%.1f (k = %.2f) to R_ac = %.0f Ω\n",
		d.RL*1e3, d.C, d.LP*1e6, 1/(2*math.Pi*math.Sqrt(d.LP*d.C)), qTank, d.N, kT, rAC)
	fmt.Printf("     the one tank (two of those in series): R_L = %.2f mΩ ∥ C = %.2f F ∥ L_p = %.1f µH, 1:%.1f to R_ac = %.0f Ω\n",
		d.RL*zt*1e3, d.C/zt, d.LP*zt*1e6, d.N, rAC*zt)
	fmt.Printf("     coils + output: a series RLC at 7.8 Hz with Q ≈ X_s/(2·R_s) = %.1f (underdamped)\n", d.Xs/(2*d.Rs))
	ms := m.Coupling()
	_, _, ks, _ := m.loop()
	ms *= 1 + ks
	const idc = 10e3
	fmt.Println("     rpm = 1768·x·(1 + h·sin 2π·7.8t), speed drives. With the matched output the 7.8 Hz sine pumps the")
	fmt.Println("     field (it grows at x = 1); the steady point is the mean speed x where it neither grows nor decays.")
	fmt.Println("     Powers at that point, scaled to 10 kA DC (the equations are linear: P ∝ I²):")
	fmt.Println("     (loop AC per machine; primary, R_L, secondary and R_ac are the one tank's, i.e. both machines together)")
	fmt.Println("       h     growth at x = 1   steady x (rpm)    loop AC [A]   primary AC [A]   R_L AC [W]   secondary [V / A peak]   R_ac [W]   primary DC heat [W]")
	run := func(h, x, tEnd float64) DualTrace {
		q := m
		q.Speed = q.held(x * w0)
		q.SpeedAmp, q.SpeedFreq = h, fSchu
		q.I0 = [2]float64{idc, idc}
		return q.RunRK4(tEnd, math.Min(2e-4, q.RLoad*q.CParallel/5))
	}
	for _, h := range []float64{0.1, 0.5, 1} {
		g1 := rmsGrowth(run(h, 1, 8).T, run(h, 1, 8).I1, 2)
		lo, hi := 0.3, 1.0 // growth rises with x: bisect for zero
		for it := 0; it < 30; it++ {
			mid := 0.5 * (lo + hi)
			tr := run(h, mid, 6)
			if rmsGrowth(tr.T, tr.I1, 2) > 0 {
				hi = mid
			} else {
				lo = mid
			}
		}
		x := 0.5 * (lo + hi)
		tr := run(h, x, 8)
		dc := mean(tr.T, tr.I1, 3)
		sc := idc / dc
		ia := cmplx.Abs(amplitude(tr.T, tr.I1, 3, fSchu)) * sc
		ipa := cmplx.Abs(amplitude(tr.T, tr.IP1, 3, fSchu)) * sc
		var prl, pac, isPk float64
		var cnt int
		for k := range tr.T {
			if tr.T[k] >= 3 {
				v := tr.V1[k] * sc            // R_L sees the whole output voltage (its DC part is only R_p·I)
				prl += v * v / (m.RLoad * zt) // the one tank
				is := tr.IS1[k] * sc
				pac += m.RAC * zt * is * is
				isPk = math.Max(isPk, math.Abs(is))
				cnt++
			}
		}
		prl /= float64(cnt)
		pac /= float64(cnt)
		fmt.Printf("       %-4.1f  %+.4f 1/s        %.4f (%4.0f)    %7.0f       %7.0f          %7.1f      %6.1f / %5.1f           %7.1f    %5.0f\n",
			h, g1, x, x*rpmBase, ia, ipa, prl, isPk*m.RAC*zt, isPk, pac, rPrim*zt*(it*idc)*(it*idc))
	}
	fmt.Println("     the sine's work at 7.8 Hz enters through the drives: it pays the AC outputs and part of the DC field,")
	fmt.Println("     so the field holds at a lower mean rpm than the DC onset (ledger below).")
	_ = ms
	{
		q := m
		q.Speed = q.held(w0)
		q.SpeedAmp, q.SpeedFreq = 1, fSchu
		q.I0 = [2]float64{idc, idc}
		tr := q.RunRK4(4, math.Min(2e-4, q.RLoad*q.CParallel/5))
		ct, err := q.Run(4, 2e-4)
		if err != nil {
			return err
		}
		fmt.Printf("     circuit at h = 1: loop AC %.0f A (RK4 %.0f); ledger mismatch %.0e of the drive work\n",
			cmplx.Abs(amplitude(ct.T, ct.I1, 1, fSchu)), cmplx.Abs(amplitude(tr.T, tr.I1, 1, fSchu)), ct.Mismatch/ct.In)
	}
	fmt.Printf("     the primary carries the loop's DC, %.0f kA: %.0f J stored at %.1f µH; its core needs an air gap (or the DC flux cancelled)\n",
		it*idc/1e3, 0.5*(d.LP*zt)*(it*idc)*(it*idc), d.LP*zt*1e6)
	fmt.Printf("     with L_p across the terminals the DC goes to the primary, not R_L: a DC output needs a DC load in series\n")
	fmt.Printf("     with the primary, or on the DC side of the machine.\n")
	fmt.Println()
	return nil
}

// runRectifier: the circulating current across the transformer against the
// power it carries, and the secondary rectified into a capacitor bank.
func runRectifier() error {
	const rPrim, rAC, kT, h, idc = 10e-6, 10.0, 0.98, 1.0, 10e3
	w0 := rpmBase * 2 * math.Pi / 60
	fmt.Println("  H. the circulating current across the transformer, and rectifying the secondary into a capacitor bank")
	fmt.Println("     tank Q swept at the same match (h = 1, steady x of the Q = 3 design, scaled to 10 kA DC; RK4):")
	fmt.Println("       tank Q   L_p [µH]   C [F]    primary AC [A peak]   secondary [A peak]   R_L + R_ac AC [W]   (the one tank)")
	var xs float64
	for _, qt := range []float64{1, 3, 10, 30} {
		m, d := outputMachine(rPrim, qt, rAC, kT)
		run := func(x, tEnd float64) DualTrace {
			q := m
			q.Speed = q.held(x * w0)
			q.SpeedAmp, q.SpeedFreq = h, fSchu
			q.I0 = [2]float64{idc, idc}
			leak := (1 - q.KT*q.KT) * q.LSec / q.RAC // the transformer's leakage time constant
			return q.RunRK4(tEnd, math.Min(math.Min(2e-4, q.RLoad*q.CParallel/5), leak/5))
		}
		lo, hi := 0.3, 1.0
		for it := 0; it < 25; it++ {
			mid := 0.5 * (lo + hi)
			tr := run(mid, 6)
			if rmsGrowth(tr.T, tr.I1, 2) > 0 {
				hi = mid
			} else {
				lo = mid
			}
		}
		x := 0.5 * (lo + hi)
		zq, _ := m.TankScale()
		if qt == 3 {
			xs = x
		}
		tr := run(x, 8)
		sc := idc / mean(tr.T, tr.I1, 3)
		var p float64
		var cnt int
		for k := range tr.T {
			if tr.T[k] >= 3 {
				v, is := tr.V1[k]*sc, tr.IS1[k]*sc
				zt, _ := m.TankScale()
				p += v*v/(m.RLoad*zt) + m.RAC*zt*is*is
				cnt++
			}
		}
		fmt.Printf("       %-6.0f   %7.1f    %6.2f    %8.0f              %6.1f               %7.0f   (steady x %.3f)\n",
			qt, d.LP*zq*1e6, d.C/zq, cmplx.Abs(amplitude(tr.T, tr.IP1, 3, fSchu))*sc, cmplx.Abs(amplitude(tr.T, tr.IS1, 3, fSchu))*sc, p/float64(cnt), x)
	}
	fmt.Println("     a higher tank Q raises the current circulating between C and the primary, but that current is the")
	fmt.Println("     primary's magnetizing current: it stays in the tank. Only the in-phase (load) part crosses the")
	fmt.Println("     transformer, and that is set by the match, at most EMF²/(8·R_s) per machine. The circulating current")
	fmt.Println("     sizes the copper and the capacitor's ripple rating; the secondary current sizes the rectifier.")

	return rectified(xs)
}

func runRectifierOnly() error { return rectified(0.756) }

// rectified runs the Q = 3 design with a bridge, a bank and a DC load on
// the secondary, at the steady mean speed x.
func rectified(xs float64) error {
	const rPrim, rAC, kT, h, idc = 10e-6, 10.0, 0.98, 1.0, 10e3
	w0 := rpmBase * 2 * math.Pi / 60
	m, d := outputMachine(rPrim, 3, rAC, kT)
	zt, _ := m.TankScale()
	// declared, per loop: 50 mF bank, pre-charged to 200 V; R_dc = 2·R_ac loads a capacitor-input bridge like R_ac
	m.CBank, m.RDC, m.BankV0 = 0.05, 2*rAC, 200*zt
	m.Speed = m.held(xs * w0)
	m.SpeedAmp, m.SpeedFreq = h, fSchu
	m.I0 = [2]float64{idc, idc}
	c, err := m.Circuit()
	if err != nil {
		return err
	}
	wf, err := analysisTran(c, 5, 5e-5)
	if err != nil {
		return err
	}
	vp, _ := wf.NodeVoltage("dp") // the one bank
	vn, _ := wf.NodeVoltage("dn")
	rdc := m.RDC * zt
	il, _ := wf.ElementCurrent("R1")
	var vb, vmin, vmax float64 = 0, math.Inf(1), math.Inf(-1)
	var cnt int
	for k := range wf.Time {
		if wf.Time[k] >= 4 {
			v := vp[k] - vn[k]
			vb += v
			vmin, vmax = math.Min(vmin, v), math.Max(vmax, v)
			cnt++
		}
	}
	vb /= float64(cnt)
	dcI := mean(wf.Time, il, 4)
	sc := idc / dcI
	fmt.Printf("     rectified (Q = 3 design, one transformer 1:%.1f, one bridge, %.0f mF bank, %.0f Ω DC load, speed held at x = %.3f, h = 1):\n", d.N, m.CBank/zt*1e3, rdc, xs)
	fmt.Printf("       at the run's %.1f kA DC per loop: bank %.1f V (ripple %.1f V p-p at 15.6 Hz), DC load %.0f W\n",
		dcI/1e3, vb, vmax-vmin, vb*vb/rdc)
	fmt.Printf("       scaled to 10 kA DC (P ∝ I², V ∝ I): bank %.0f V, DC load %.0f W\n", vb*sc, vb*vb*sc*sc/rdc)
	fmt.Printf("       engine ledger mismatch %.0e of the drive work (diode and bank losses included)\n",
		wf.Ledger.IntegrationMismatch/wf.Ledger.EnergyIn)
	fmt.Println("       the bank charges to about the secondary's peak; the load sees DC; the energy came from the speed")
	fmt.Println("       drives' 7.8 Hz work, through the field, the tank and the transformer.")
	fmt.Println()
	return nil
}

func analysisTran(c *circuit.Circuit, tEnd, dt float64) (*output.Waveforms, error) {
	return analysis.Tran(c, analysis.TranSpec{Step: dt, Stop: tEnd, MaxStep: dt, UIC: true}, analysis.TranOptions{})
}
