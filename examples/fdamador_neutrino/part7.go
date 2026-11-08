// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The report for Part 7 (H40): a spinning ring.

import (
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"time"

	"github.com/fdamador/kronspice/spacetime"
)

// Part7a holds the rotating-ring-line results.
type Part7a struct {
	Sagnac     []SagnacRun
	Two        TwoState
	Scan       []RotCluster
	Fit        RotCluster // δ_CP free within the NuFIT band
	Pinned     RotCluster // δ_CP pinned at 212°
	ThreeE     RotThree   // pinned point, from e
	ThreeMu    RotThree   // pinned point, from μ
	ZeroE      RotThree   // same defects, no rotation, from e
	ZeroMu     RotThree
	TMeas      float64 // max |P(e→μ) − P(μ→e)|, kronspice, with rotation
	TPred      float64 // the same from the exact expansion
	TZeroMeas  float64 // without rotation
	TZeroPred  float64
	Rims       []RimRow
	FitSeconds float64
}

// RimRow is one rim-speed estimate.
type RimRow struct {
	Case  string
	Guide string
	Rim   SagnacRim
}

// Part7b holds the field results for one grid.
type Part7b struct {
	Cells   float64
	Modes   []Mode
	Runs    []SpinRun
	Seconds float64
}

const bench7a = 10e-3 // ring radius for the rim-speed estimates [m] (the bench ring of Part 1)

func reportPart7a(cpAll bool) (Part7a, error) {
	var p Part7a
	t0 := time.Now()
	fmt.Println("\n7. A spinning ring (H40, user hypothesis 2026-11-07): \"If the torus spins about its axis and spins around the E")
	fmt.Println("   axis (the ring), would this also add to the missing information of neutrino data?\" Missing after Parts 1–6: mass,")
	fmt.Println("   weak interaction, a CP phase (the reciprocal ring gives only δ_CP = 0 or 180°), left-handed helicity, three")
	fmt.Println("   generations, confinement. Two readings, kept apart: 7a a material ring that turns (the moving medium makes the")
	fmt.Println("   ring non-reciprocal: Sagnac); 7b a field that spins (circulating E×B, in the field solver).")
	fmt.Println("7a. Rotation on the ring line. A ring turning at Ω delays the co-rotating wave by Δt = 4AΩ/c² per turn for any")
	fmt.Println("    co-moving medium (Sagnac 1913; Post, Rev. Mod. Phys. 39, 475, 1967), so a resonance splits by Δω/ω = 2aΩ·v_g/c²")
	fmt.Println("    (spacetime.SagnacSplitting × v_g/c). In the ring's frame this is a gyroscopic (Coriolis) term, antisymmetric and")
	fmt.Println("    workless (Landau & Lifshitz, Mechanics §39; Kron's speed tensor G): a kronspice gyrator Gy across every section,")
	fmt.Println("    Cmat·du/dt = −g + 2·Gy·Q·u. Uniform ring: c_κω² ± 4·Gy·sin κ·ω − 1/Lg = 0, so CW and CCW split by 4·Gy·sin κ/c_κ")
	fmt.Println("    exactly. Gy is set from the non-reciprocity ν = Δω/ω on a chosen pair; ν is turned into a rim speed below.")
	fmt.Println("    Test 1, harmonic 1 split by rotation, started standing (S at c = ½), kronspice transient, 60 steps a period:")
	fmt.Printf("    %-8s %7s %10s %12s %12s %12s %12s %12s %9s %9s %9s %8s\n", "ring", "ν", "Gy [S]", "Δω exact", "Δω modes", "Δω 1st ord.", "exact, trap.", "kronspice", "S→D amp", "gyr pass", "gyr abs", "ledger")
	for _, crossed := range []bool{false, true} {
		s, err := Sagnac(crossed, 0.01, 120)
		if err != nil {
			return p, err
		}
		p.Sagnac = append(p.Sagnac, s)
		fmt.Printf("    %-8s %7.3g %10.4g %12.6g %12.6g %12.6g %12.6g %12.6g %9.5f %9.2e %9.1e %8.1e\n", s.Ring, s.Nu, s.Gyro, s.DwExact, s.DwModes, s.DwFirst,
			s.DwWarped, s.DwKron, s.ConvAmp, s.GyroPass, s.GyroPowerMax, s.LedgerRel)
	}
	s0 := p.Sagnac[0]
	fmt.Printf("    (Δω in rad/s, ω_CW − ω_CCW; f = %.5g MHz. 'exact, trap.' = the exact frequencies after the trapezoidal rule's warping\n", s0.F0/1e6)
	fmt.Printf("    ω → (2/h)·atan(ωh/2); kronspice/that = %.7f. gyr pass = largest power one gyrator passes from section to section,\n", s0.DwKron/s0.DwWarped)
	fmt.Printf("    gyr abs = largest power any gyrator absorbs, both / (ω·W): the rotation terms move energy round the ring and do no\n")
	fmt.Printf("    net work, ∫ΣP dt = %.1e of the stored energy; booked as the rotation's (mechanical) column, it is zero.)\n", s0.GyroWork)
	fmt.Println("    Rotation alone turns S into D and back (the pattern turns): amplitude 1, the two-state beat with CW/CCW as modes.")
	two, err := TwoStateRot(0.1, 1, 1.25)
	if err != nil {
		return p, err
	}
	p.Two = two
	fmt.Println("    Test 2, defect and rotation together (δ = 0.1 on section 0, rotation splitting the pair by the same amount), started CW:")
	fmt.Printf("      splittings: defect alone %.6g Hz, rotation alone %.6g Hz, both (exact modes) %.6g Hz, √(Δd² + Δr²) %.6g Hz\n", two.DfDefect, two.DfRot, two.DfExact, two.DfQuad)
	fmt.Printf("      CW → CCW: sin²2θ predicted Δd²/(Δd² + Δr²) = %.4f (θ %.2f°); kronspice %.4f (θ %.2f°), Δf %.6g Hz; max |res| %.1e, ledger %.1e\n",
		two.AmpPred, two.ThetaPred, two.AmpFit, two.ThetaFit, two.DfFit, two.MaxResid, two.LedgerRel)
	fmt.Println("      Part 6 found θ = 45° for every reciprocal defect in the CW/CCW basis; rotation gives any θ: the modes become")
	fmt.Println("      elliptical (complex mixtures of CW and CCW).")

	fmt.Println("    Test 3, the three-state cluster of Part 6 with rotation (ν on the uniform pair at κ = 11π/12; Part 6 defects 0.2064 on")
	fmt.Println("    section 2, 0.1628 on section 23; e, μ, τ = f1, f2, f3). The modes are complex, so U can carry a CP phase:")
	fmt.Printf("    %10s %9s %9s %8s %8s %8s %8s %8s %9s %6s\n", "ν", "δ_CP", "J", "J/J_max", "θ12", "θ23", "θ13", "ratio", "cost", "leak")
	for _, nu := range []float64{0, -1e-5, -3e-5, -1e-4, -3e-4, -1e-3, -3e-3, -1e-2, 1e-4, 1e-3} {
		c := RotClusterMix(0.2064, 0.1628, 2, 23, nu, [3]int{0, 1, 2})
		p.Scan = append(p.Scan, c)
		fmt.Printf("    %10.3g %8.1f° %9.5f %8.4f %7.2f° %7.2f° %7.2f° %8.3f %9.4g %6.3f\n", nu, c.DeltaCP, c.J, c.J/c.JMax, c.Th12, c.Th23, c.Th13, c.Ratio, c.Cost, c.Leak)
	}
	fmt.Println("    (J = Im(U_e1 U_μ2 U*_e2 U*_μ1) of the unitary closest to the projected U (Löwdin); J_max = c12 s12 c23 s23 c13² s13;")
	fmt.Println("    reversing the rotation gives 360° − δ_CP. Measured (NuFIT 6.0, NO with SK-atm): δ_CP = 212° (+26/−41).)")
	t1 := time.Now()
	nus := []float64{-1e-4, -3e-4, -1e-3}
	if cpAll {
		p.Fit = FitCPAll(nus, 0)
		p.Pinned = FitCPAll(nus, 1)
	} else {
		p.Fit = FitCP(2, 23, nus, 0)
		p.Pinned = FitCP(2, 23, nus, 1)
	}
	p.FitSeconds = time.Since(t1).Seconds()
	m := NuFIT6()
	for _, x := range []struct {
		name string
		c    RotCluster
	}{{"δ_CP free within the NuFIT band", p.Fit}, {"δ_CP pinned at 212° (±1°)", p.Pinned}} {
		c := x.c
		fmt.Printf("    fit, %s (sections %d, %d): δ = %.4f, %.4f, ν = %.4g (Gy = %.4g S)\n", x.name, c.P1, c.P2, c.D1, c.D2, c.Nu, c.Gyro)
		fmt.Printf("      δ_CP %.1f°, J %.5f (J_max %.5f); θ12 %.2f° (%.2f ± %.2f), θ23 %.2f° (%.1f ± %.1f), θ13 %.2f° (%.2f ± %.2f), ratio %.2f (%.2f); cost %.3g\n",
			c.DeltaCP, c.J, c.JMax, c.Th12, m.Th12, m.S12, c.Th23, m.Th23, m.S23, c.Th13, m.Th13, m.S13, c.Ratio, m.Ratio, c.Cost)
	}
	if cpAll {
		fmt.Println("    (every section pair searched)")
	} else {
		fmt.Println("    (Part 6's sections; a search over every section pair, run once with -cpall in 5.4 min, found only the mirror pair")
		fmt.Println("    (1, 22) with the same cost)")
	}
	fmt.Printf("    Without the pin the fit stays near Part 6 and takes δ_CP inside the measured 1σ band; pinned at 212° the defects cannot\n")
	fmt.Printf("    hold θ13 (%.1fσ) and the splitting ratio (%.1f %%) at the same time. (fits %.0f s)\n",
		math.Abs(p.Pinned.Th13-m.Th13)/m.S13, 100*(p.Pinned.Ratio/m.Ratio-1), p.FitSeconds)
	U := p.Pinned.UU
	fmt.Println("    pinned point, U (rows e, μ, τ; columns ν1, ν2, ν3), |U| and arg U:")
	for a := 0; a < 3; a++ {
		fmt.Printf("      %.4f∠%+7.1f°  %.4f∠%+7.1f°  %.4f∠%+7.1f°\n", cmplx.Abs(U[a][0]), cmplx.Phase(U[a][0])*180/math.Pi,
			cmplx.Abs(U[a][1]), cmplx.Phase(U[a][1])*180/math.Pi, cmplx.Abs(U[a][2]), cmplx.Phase(U[a][2])*180/math.Pi)
	}

	// T violation in kronspice at the pinned point, and with the rotation off.
	cl := p.Pinned
	e, mu := cl.Assign[0], cl.Assign[1]
	var errs [4]error
	p.ThreeE, errs[0] = RotThreeFlavour(cl, e, 1)
	p.ThreeMu, errs[1] = RotThreeFlavour(cl, mu, 1)
	z := RotClusterMix(cl.D1, cl.D2, cl.P1, cl.P2, 0, cl.Assign)
	p.ZeroE, errs[2] = RotThreeFlavour(z, e, 1)
	p.ZeroMu, errs[3] = RotThreeFlavour(z, mu, 1)
	for _, err := range errs {
		if err != nil {
			return p, err
		}
	}
	p.TMeas, p.TPred = TAsym(p.ThreeE, p.ThreeMu, e, mu)
	p.TZeroMeas, p.TZeroPred = TAsym(p.ZeroE, p.ZeroMu, e, mu)
	fmt.Printf("    kronspice at the pinned point, started in e and in μ, %.0f periods each (one slow beat):\n", p.ThreeE.Periods)
	fmt.Printf("      against the exact complex mode expansion (all %d modes, same readout): largest deviation %.2e (e), %.2e (μ)\n",
		len(cl.Modes.Omega), p.ThreeE.MaxDevAll, p.ThreeMu.MaxDevAll)
	fmt.Printf("      against the 3×3 formula with the unitarized U: %.3f (the %.1f %% of each mode outside the three states)\n", p.ThreeE.MaxDev3, 100*cl.Leak)
	fmt.Printf("      max |P(e→μ)(t) − P(μ→e)(t)|: kronspice %.4f, exact expansion %.4f; rotation off (same defects): %.4f, %.4f\n",
		p.TMeas, p.TPred, p.TZeroMeas, p.TZeroPred)
	fmt.Printf("      ledgers %.1e, %.1e; gyrators' net work %.1e, %.1e of the ring's energy\n", p.ThreeE.LedgerRel, p.ThreeMu.LedgerRel, p.ThreeE.GyroRel, p.ThreeMu.GyroRel)
	fmt.Println("      Without rotation the difference is the readout's (the flavour projections are not exact energy eigen-projections) and")
	fmt.Println("      the exact expansion has it too; with rotation P(e→μ) and P(μ→e) separate: T (and so CP) violation in the ring.")

	// rim speeds
	r := DefaultFRing(false)
	k11 := r.Kappa(r.N/2 - 1)
	f11 := r.OmegaK(k11) / (2 * math.Pi)
	burst := spacetime.BurstSpeed(spacetime.T1000G)
	fmt.Printf("    What turning the ring takes, on the bench ring a = %.0f mm: aΩ = ν·c²/(2v_g); burst = T1000G carbon-fibre hoop, %.0f m/s\n", bench7a*1e3, burst)
	fmt.Println("    (spacetime.BurstSpeed; the strongest material in the package, a ceramic ring bursts far lower). v_g = c is the most")
	fmt.Println("    favourable passive guide; c/10 the εr = 100 ring; the circuit line's own wave laid on the ring is a slow wave.")
	fmt.Printf("    %-36s %-26s %11s %11s %10s %10s %11s %8s\n", "case", "guide (v_g)", "ν", "aΩ [m/s]", "(v/c)²", "burst ×", "½v² [J/kg]", "possible")
	type cs struct {
		name string
		nu   float64
	}
	for _, c := range []cs{{"δ_CP in the 1σ band (fit)", p.Fit.Nu}, {"δ_CP = 212° (pinned)", p.Pinned.Nu}, {"Test 1 (harmonic 1)", 0.01}} {
		for _, gd := range []struct {
			name string
			vg   float64
		}{{"c", C0}, {"c/10 (εr = 100)", C0 / 10}, {fmt.Sprintf("circuit line (%.3g m/s)", r.GroupVelocityRing(k11, bench7a)), r.GroupVelocityRing(k11, bench7a)}} {
			rim := Rim(c.nu, gd.vg, bench7a, f11)
			p.Rims = append(p.Rims, RimRow{c.name, gd.name, rim})
			b2 := fmt.Sprintf("%10.2e", rim.Beta2)
			bm := fmt.Sprintf("%10.3g", rim.BurstMargin)
			ke := fmt.Sprintf("%11.3g", rim.KinPerKg)
			if !rim.BelowC {
				b2, bm, ke = "       ≥ c", "         —", "          —"
			}
			fmt.Printf("    %-36s %-26s %11.3g %11.3g %s %s %s %8v\n", c.name, gd.name, math.Abs(c.nu), rim.RimSpeed, b2, bm, ke, rim.Possible)
		}
	}
	fmt.Printf("    (twin: the same ν from spacetime.SagnacSplitting(πa², 2πa, c/f, Ω)·v_g/c, relative difference ≤ %.1e)\n", maxTwin(p.Rims))
	fmt.Printf("    At the burst speed (v_g = c) the reachable ν is 2·%.0f/c = %.2e, %.0f× below the 212° point; δ_CP there is %.2f° from 180°.\n",
		burst, 2*burst/C0, math.Abs(p.Pinned.Nu)/(2*burst/C0), math.Abs(RotClusterMix(p.Pinned.D1, p.Pinned.D2, p.Pinned.P1, p.Pinned.P2,
			math.Copysign(2*burst/C0, p.Pinned.Nu), p.Pinned.Assign).DeltaCP-180))
	fmt.Println("    Energy: the gyrators do no net work at constant Ω (above); the ring's own kinetic energy, ½v² per kilogram of ring,")
	fmt.Println("    is paid by whatever turns it, and in the lab a co-rotating wave carries ΩJ_z more energy than in the ring's frame,")
	fmt.Println("    also paid by the drive. Nothing enters from the rotation itself.")
	fmt.Println("    Standard prediction: rotation-induced non-reciprocity is ordinary Sagnac physics (ring-laser gyroscopes); in the")
	fmt.Println("    Standard Model CP violation in neutrino mixing is a complex phase of the weak-interaction (PMNS) matrix, not rotation.")
	fmt.Printf("   (Part 7a: %s)\n", time.Since(t0).Round(time.Second))
	return p, nil
}

func maxTwin(rs []RimRow) float64 {
	var m float64
	for _, r := range rs {
		if r.Rim.BelowC {
			m = math.Max(m, r.Rim.TwinRel)
		}
	}
	return m
}

func reportPart7b(cells []float64) ([]Part7b, error) {
	var out []Part7b
	fmt.Println("7b. The spinning field (em/fdtd.go), the torus of Part 5b (εr = 100, lossless, absorber from 3a, diagnostics inside the")
	fmt.Println("    ±2a box). J_z = ε0∫ρ(E×B)_φ dV (Abraham; Minkowski ∫ρ(D×B)_φ beside it), helicity h = (1/2c)∫(A·B/μ0 − C·D/ε0) dV with")
	fmt.Println("    B = curl A, D = −curl C (temporal gauges, carried beside the solver; h·ω/U = ±1 for circular light, tested), ∫E·B.")
	fmt.Println("    For modes of azimuthal order m_k, J_z = Σ m_k U_k/ω_k: J_z·ω/U is an energy-weighted mean of whole numbers.")
	for _, c := range cells {
		t0 := time.Now()
		sf := DefaultSpinField()
		sf.Ring.Cells = c
		st := time.Now()
		sf.Log = func(f string, a ...any) {
			fmt.Fprintf(os.Stderr, "[7b r/%g %s] "+f+"\n", append([]any{c, time.Since(st).Round(time.Second)}, a...)...)
		}
		pb := Part7b{Cells: c}
		g := sf.grid()
		mE, err := sf.FindMode("E-ring (TE01δ, m = 0)", "E", 1.0e9, 4e9)
		if err != nil {
			return out, err
		}
		mB, err := sf.FindMode("B-ring (B round the ring, m = 0)", "B", 1.0e9, 6e9)
		if err != nil {
			return out, err
		}
		m1, err := sf.FindMode("m = 1 (cos φ)", "m1", 1.0e9, 6e9)
		if err != nil {
			return out, err
		}
		m1s := g.Rot90(m1)
		pb.Modes = []Mode{mE, mB, m1}
		fmt.Printf("    grid r/%g (dx = %.3g mm, %d nodes): modes found (probe spectrum peak, then filtered): E-ring %.5g GHz, B-ring %.5g GHz,\n",
			c, g.dx*1e3, g.n*g.n*g.n, mE.F/1e9, mB.F/1e9)
		fmt.Printf("      m = 1 %.5g GHz (lumped E-ring of Part 1: 2.0606 GHz)\n", m1.F/1e9)
		s2 := complex(math.Sqrt2/2, 0)
		i := complex(0, 1)
		type rc struct {
			name string
			ms   []Mode
			a    []complex128
		}
		cases := []rc{
			{"E-ring alone", []Mode{mE}, []complex128{1}},
			{"B-ring alone", []Mode{mB}, []complex128{1}},
			{"E + B, quadrature (+90°)", []Mode{mE, mB}, []complex128{s2, i * s2}},
			{"E + B, quadrature (−90°)", []Mode{mE, mB}, []complex128{s2, -i * s2}},
			{"E + B, in phase", []Mode{mE, mB}, []complex128{s2, s2}},
			{"m = 1 standing (cos φ)", []Mode{m1}, []complex128{1}},
			{"m = 1 travelling (cos φ + i sin φ)", []Mode{m1, m1s}, []complex128{s2, i * s2}},
			{"m = 1 travelling (cos φ − i sin φ)", []Mode{m1, m1s}, []complex128{s2, -i * s2}},
			{"½ E-ring + ½ m = 1 travelling", []Mode{mE, m1, m1s}, []complex128{s2, 0.5, i * 0.5}},
		}
		for _, x := range cases {
			r, err := sf.RingDown(x.name, x.ms, x.a)
			if err != nil {
				return out, err
			}
			pb.Runs = append(pb.Runs, r)
			sf.Log("%s done", x.name)
		}
		pb.Seconds = time.Since(t0).Seconds()
		fmt.Printf("      %d ring-downs of %.0f lumped periods each (%.0f s for the grid). J_z·ω/U and h·ω/U at the start, the end, and their range:\n",
			len(cases), sf.RunT, pb.Seconds)
		fmt.Printf("      %-36s %8s %18s %8s %18s %8s %18s %8s %7s %8s\n", "field", "Jzω/U A", "range", "Jzω/U M", "range", "hω/U 0", "range", "circ.", "Q", "ledger")
		for _, r := range pb.Runs {
			circ := r.JzAbs[0] * r.Omega / r.U[0]
			fmt.Printf("      %-36s %+8.4f [%+7.4f, %+7.4f] %+8.4f [%+7.4f, %+7.4f] %+8.4f [%+7.4f, %+7.4f] %8.4f %7.3g %8.1e\n", r.Name, r.JzW0, r.JzWMin, r.JzWMax,
				r.JzMW0, r.JzMWMin, r.JzMWMax, r.HW0, r.HWMin, r.HWMax, circ, r.Q, r.Residual)
		}
		fmt.Println("      (A = Abraham ε0∫ρ(E×B)_φ, M = Minkowski ∫ρ(D×B)_φ; start, and range while U ≥ 1e-3 U0; circ. = ε0∫ρ|(E×B)_φ| dV·ω/U,")
		fmt.Println("      the azimuthal circulation whatever its sign; ω = U/Σ(U_k/ω_k); Q from the decay of the energy inside the box;")
		fmt.Println("      ledger = field-ledger residual / start energy. h·ω/U with this vacuum-normalized h is ±1 for circular light in")
		fmt.Printf("      vacuum and (1 + n²)/2n = %.2f inside εr = 100: values above 1 are the medium's)\n", (1+100.0)/(2*10))
		r := pb.Runs
		fmt.Printf("      potentials at the start, max |curl A − B| / max |B| inside: E + B %.2f, m = 1 travelling %.2f. They are built from the\n", r[2].CurlARes, r[6].CurlARes)
		fmt.Println("      filtered phasors (A = iÊ/ω), which are not exact single-frequency solutions, so h is read only for its sign and its")
		fmt.Println("      swings. The potential-free pseudoscalar χ = c·ε0·∫E·B dV / U carries the handedness without that error:")
		for _, k := range []int{0, 2, 3, 4, 6} {
			lo, hi := math.Inf(1), math.Inf(-1)
			for q := 0; q <= r[k].Last; q++ {
				x := C0 * Eps0 * r[k].EB[q] / r[k].U[q]
				lo, hi = math.Min(lo, x), math.Max(hi, x)
			}
			fmt.Printf("        %-36s χ at the start %+.4f, range [%+.4f, %+.4f]\n", r[k].Name, C0*Eps0*r[k].EB[0]/r[k].U[0], lo, hi)
		}
		az := g.AzimuthalContent([]Mode{m1, m1s}, []complex128{s2, i * s2})
		fmt.Printf("      azimuthal content of the travelling m = 1 (share of Σ|E_φ|² on the ring's centre circle): m = 1 %.3f, −3 %.3f, 5 %.3f, −7 %.3f;\n",
			az[1], az[-3], az[5], az[-7])
		fmt.Println("      the cubic grid has only fourfold symmetry, which couples m to m ± 4: on the staircased torus J_z is conserved only")
		fmt.Println("      modulo 4 and the 'm = 1' lattice mode is not a pure order, so its J_z·ω/U is not a clean 1 in either form.")
		fmt.Printf("      the m = 1 mode: %.1f %% of its energy inside the torus; overlap with itself turned 90° %+.4f, 180° %+.4f (odd order)\n",
			100*m1.Inside, m1.Turn90, m1.Turn180)
		out = append(out, pb)
	}
	fmt.Println("    Reading 1 of 'spin about the axis' (E-ring + B-ring, m = 0): the azimuthal E×B circulation is there, but it integrates")
	fmt.Println("    to J_z = 0 (the field is unchanged by a rotation about the axis); the handedness (χ, h) starts with the sign of the")
	fmt.Println("    relative phase and then swings through both signs as the two modes (different frequencies) beat: it is not conserved,")
	fmt.Println("    as expected where the dielectric breaks duality symmetry (Fernandez-Corbaton et al., PRL 111, 060401, 2013).")
	if len(out) > 0 {
		r := out[0].Runs
		fmt.Printf("    Reading 2 (a pattern that turns, m = 1): net J_z, its sign set by the phase (the ± runs are mirror images: Minkowski\n")
		fmt.Printf("    %+.4f / %+.4f, Abraham %+.4f / %+.4f at the start); its size is not a clean 1 on this grid (above), the Abraham form\n",
			r[6].JzMW0, r[7].JzMW0, r[6].JzW0, r[7].JzW0)
		fmt.Printf("    is small because the field's momentum sits in the dielectric, and Q is the standing mode's (%.4g against %.4g).\n", r[6].Q, r[5].Q)
	}
	fmt.Println("    ½ per quantum: J_z·ω/U is an energy-weighted mean of whole m, so a value of ½ can only be a mixture of states")
	fmt.Println("    (the ½ E-ring + ½ m = 1 run), never one coherent classical mode.")
	fmt.Println("    Measured: neutrinos are left-handed (helicity −½), antineutrinos right-handed. Maxwell's equations are parity-even, so the")
	fmt.Println("    sign of J_z or h is set by how the field is started; nothing in the model fixes it. Standard prediction: a field's")
	fmt.Println("    helicity and angular momentum are ordinary electromagnetism; neutrino handedness comes from the V−A weak interaction.")
	return out, nil
}

// reportMass7 states what spinning adds to the energy and what it does to
// confinement, from 7a and 7b.
func reportMass7(p7a Part7a, p7b []Part7b) {
	fmt.Println("7c. Mass and confinement. Rotation adds energy, so mass E/c², but no binding:")
	burst := spacetime.BurstSpeed(spacetime.T1000G)
	fmt.Printf("    a material ring at the burst speed carries ½v² = %.3g J/kg, ΔM/M = v²/2c² = %.2e; a wave of order m in a ring turning\n",
		0.5*burst*burst, 0.5*burst*burst/(C0*C0))
	e := DefaultERing()
	om := burst / bench7a
	fmt.Printf("    at Ω has energy U(1 + mΩ/ω) in the lab: at the burst Ω = %.3g rad/s on the bench ring (ω = %.4g rad/s), mΩ/ω = %.2e for m = 1.\n",
		om, e.Omega0(), om/e.Omega0())
	if len(p7b) > 0 {
		r := p7b[0].Runs
		fmt.Printf("    Field solver (r/%g): Q of the m = 1 mode standing %.4g, travelling %.4g and %.4g; E-ring %.4g, E + B %.4g (decay times\n",
			p7b[0].Cells, r[5].Q, r[6].Q, r[7].Q, r[0].Q, r[2].Q)
		w := 2 * math.Pi * p7b[0].Modes[2].F
		fmt.Printf("    Q/ω: m = 1 standing %.3g ns, travelling %.3g ns). Spinning neither lengthens nor shortens the decay: the turning\n",
			r[5].Q/w*1e9, r[6].Q/w*1e9)
		fmt.Println("    pattern radiates as the standing one (to the grid's resolution), and in vacuum nothing holds either.")
	} else {
		fmt.Println("    (field-solver decay times: run with -field)")
	}
	fmt.Println("    A centrifugal pull outward needs a wall to push against; the dielectric is that wall here, and a neutrino has none.")
	_ = p7a
}
