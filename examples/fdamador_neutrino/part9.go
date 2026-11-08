// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

// The report for Part 9 (H42): how two of these rings interact.

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

// Part9 holds what the report computed (also for viz9.go).
type Part9 struct {
	Single       AMode
	SingleRun    ArrayRun
	Pairs        []PairRow
	Beats        []Beats
	Forces       []ForceAt
	ForceGeoms   []Geom
	Best         PairRow // largest antiphase Q (with loss), coaxial
	BestRad      PairRow // largest antiphase radiation-only Q, coaxial
	Chain        []ChainRow
	ChainD       float64
	ChainExp     float64 // fitted exponent of Q_rad(dark) ∝ N^p over the largest N
	ChainCoef    float64
	NNeeded      float64 // rings for Q = Part 2's Q needed at that scaling
	Field        []PairFieldResult
	FieldLabel   []string
	FieldPred    []float64 // circuit radiation-only Q for the same configuration
	FieldPredF   []float64
	Nu           NuPair
	Seconds      float64
	FieldSeconds float64
}

// PairRow is one distance of the pair table.
type PairRow struct {
	G                 Geom
	D, KD             float64
	ML, RmR, RmDip    float64 // (M + ΔM)/(L + ΔL), R_m/R11 exact and dipole form, at f0
	In, Anti          AMode
	QInK, QAntiK      float64 // kronspice mode runs
	LedgerK, MismK    float64 // worst kronspice ledger residual and energy mismatch of the two runs
	MinRadK           float64 // smallest Iᵀ·R_rad·I / its largest over both runs
	SplitRel          float64 // (f_anti − f_in)/f_single
	QAntiOverSingle   float64
	QRadAntiOverSingl float64
}

// ChainRow is one chain length.
type ChainRow struct {
	N            int
	Bright, Dark AMode
	DarkRad      AMode // the same chain without dielectric loss (exact radiation-only mode)
	QDarkK       float64
	QBrightK     float64
	LedgerK      float64
	MismK        float64
	Kron         bool
}

// NuPair is test 9e: two model rings at the Part 2 neutrino scale.
type NuPair struct {
	A, KA, Omega0, EpsR float64
	Kappa2a             float64 // (M + ΔM)/(L + ΔL) coaxial at d = 2a
	Shift2a             float64 // ½ħω0·κ: the energy of each mode moves by this [eV]
	KappaSide           float64 // side by side at centre spacing 2.5a (gap 0.5a)
	ShiftSide           float64
	Rm2a                float64 // R_m/R11 coaxial at 2a
	QPairGain           float64 // antiphase/single radiation Q at 2a, Part 2 ring
	EpsPair             float64 // effective εr the antiphase pair (d = 2r) would need for Q needed
	YukawaExp           float64 // r/λ_Z at r = 2a (the weak potential falls as e^{−r/λ_Z})
	ContactV            float64 // √2·G_F(ħc)³ averaged over one ring's volume [eV]
	SigmaWeak           float64 // ~G_F²s/π for two 10 MeV neutrinos head on [cm²]
	SigmaModel          float64 // perturbative energy-exchange cross-section of two model rings passing at c [cm²]
	BHalf               float64 // impact parameter (in a) at which the transfer phase is π/2
	PhaseContact        float64 // transfer phase at b = a
}

const (
	gF     = 1.1663787e-5               // Fermi constant G_F/(ħc)³ [GeV⁻²] (PDG 2024)
	hbarC  = 0.1973269804               // ħc [GeV·fm]
	gev2mb = 0.3893793721               // (ħc)² [GeV²·mb]
	lamZ   = 0.1973269804e-15 / 91.1876 // ħ/(M_Z c) [m] (M_Z = 91.1876 GeV)
)

func runPart9(field bool, cells []float64, viz, fieldCSV string) error {
	t0 := time.Now()
	var p Part9
	e := DefaultERing()
	fmt.Println("\n9. Two neutrinos (H42, user question 2026-11-07): \"how would two neutrinos interact with each other?\" In the model:")
	fmt.Println("   two identical E-rings (Part 1's bench ring, TE01δ) a distance d apart. They carry no charge, so their whole")
	fmt.Println("   coupling is the vector potential of one ring's displacement current acting round the other: the retarded")
	fmt.Println("   Neumann integral Z12 = jω(μ0/4π)∮∮dl1·dl2 e^{−jkR}/R = R_m + jω(M + ΔM) (static M by em.MutualCircles for")
	fmt.Println("   coaxial rings, em.Neumann for side-by-side rings; ΔM and R_m from the cos and sin parts of the kernel). The")
	fmt.Println("   same kernel gives the self terms R11, ΔL, so subradiance (R11 − R_m) is computed consistently. Each mode is")
	fmt.Println("   solved at its own frequency; kronspice: loops C–L–R11–R_d, K elements for M, H elements (V = R_m·I_j) for the")
	fmt.Println("   mutual radiation resistance, whose power is the interference term of the radiated power (model-supplied column).")

	// single ring in the Part 9 model
	one := NewArray(e, Coaxial, e.A, 1)
	p.Single = one.SolveMode([]float64{1})
	sr, err := one.RunMode(p.Single, 40, 200)
	if err != nil {
		return err
	}
	p.SingleRun = sr
	fmt.Printf("   One ring with the kernel's self terms: R11 = %.4g Ω (Part 1's small-loop formula %.4g Ω at f0), ΔL/L = %.4f;\n",
		one.R11(e.Omega0()), e.RRad(), one.LSelf(e.Omega0())/e.L()-1)
	fmt.Printf("   f0 = %.5g GHz, Q = %.4g (radiation alone %.4g); kronspice ring-down Q %.4g, ledger %.1e. Part 1: 2.0606 GHz, Q 40.52;\n",
		p.Single.F/1e9, p.Single.Q, p.Single.QRad, sr.Q, sr.LedgerRel)
	fmt.Println("   the first retardation terms move the lumped ring by a few per cent; every pair below is compared with this ring.")

	// 9a
	fmt.Println("9a. Coupled loops: normal modes against d (matrix: self-consistent modes; kronspice: each mode started alone).")
	fmt.Println("    In phase = currents the same way round (I2 = I1), antiphase = opposite (I2 = −I1); k = (M + ΔM)/(L + ΔL) at f0.")
	fmt.Println("    geometry      d/a   kd     k        R_m/R11 (dipole)  f_in GHz  Q_in   kron   f_anti GHz  Q_anti  kron    Q_rad,anti  split")
	type cfg struct {
		g  Geom
		dd []float64
	}
	for _, c := range []cfg{{Coaxial, []float64{0.5, 0.75, 1, 1.5, 2, 3, 5, 10}}, {Coplanar, []float64{2.5, 3, 4, 5, 7.5, 10}}} {
		for _, dd := range c.dd {
			d := dd * e.A
			a := NewArray(e, c.g, d, 2)
			in, an := a.PairModes()
			w0 := e.Omega0()
			m, rm := a.Coupling(0, 1, w0)
			row := PairRow{G: c.g, D: d, KD: w0 * d / C0, ML: m / a.LSelf(w0), RmR: rm / a.R11(w0), RmDip: real(DipoleZ12(c.g, w0*d/C0)), In: in, Anti: an}
			ri, err := a.RunMode(in, 20, 200)
			if err != nil {
				return err
			}
			ra, err := a.RunMode(an, 40, 200)
			if err != nil {
				return err
			}
			row.QInK, row.QAntiK = ri.Q, ra.Q
			row.LedgerK = math.Max(ri.LedgerRel, ra.LedgerRel)
			row.MismK = math.Max(math.Abs(ri.Mismatch), math.Abs(ra.Mismatch))
			row.MinRadK = math.Min(ri.MinRadP, ra.MinRadP)
			row.SplitRel = (an.F - in.F) / p.Single.F
			row.QAntiOverSingle = an.Q / p.Single.Q
			row.QRadAntiOverSingl = an.QRad / p.Single.QRad
			p.Pairs = append(p.Pairs, row)
			if c.g == Coaxial && an.Q > p.Best.Anti.Q {
				p.Best = row
			}
			if c.g == Coaxial && an.QRad > p.BestRad.Anti.QRad {
				p.BestRad = row
			}
			fmt.Printf("    %-12s %5.2f %5.3f %+8.4f %+8.4f (%+7.4f) %8.4f %6.4g %6.4g %10.4f %7.4g %7.4g %10.4g %+7.4f\n",
				c.g, dd, row.KD, row.ML, row.RmR, row.RmDip, in.F/1e9, in.Q, row.QInK, an.F/1e9, an.Q, row.QAntiK, an.QRad, row.SplitRel)
		}
	}
	var led, mis, minr float64 = 0, 0, 1
	for _, r := range p.Pairs {
		led, mis, minr = math.Max(led, r.LedgerK), math.Max(mis, r.MismK), math.Min(minr, r.MinRadK)
	}
	fmt.Printf("    kronspice: every mode run's ledger residual ≤ %.1e, energy books (start − end − radiated − heat − interference) ≤ %.1e\n", led, mis)
	fmt.Printf("    of the start (trapezoidal step T/200), total radiated power Iᵀ·R_rad·I ≥ %.1e of its peak (never negative).\n", minr)
	fmt.Println("    Coaxial rings: M > 0, so the in-phase mode (L + M) lies below f0 and the antiphase mode (L − M) above; side by side")
	fmt.Println("    M < 0 and the order turns over. The dipole form of R_m/R11 agrees with the exact kernel from d ≈ 3a outward.")
	for _, dd := range []float64{1, 2} {
		bt, err := RunBeats(e, Coaxial, dd*e.A, 60)
		if err != nil {
			return err
		}
		p.Beats = append(p.Beats, bt)
		fmt.Printf("    Beats, coaxial d = %ga: ring 0 started alone (capacitor charged), ring 1 at rest, 60 periods: E0 − E1 (carrier-averaged)\n", dd)
		fmt.Printf("      crosses zero %d times, spacing → Δf = %.5g MHz; the circuit as integrated (matrices frozen at the mean frequency)\n",
			bt.Crossings, bt.DfMeasured/1e6)
		fmt.Printf("      has f_anti − f_in = %.5g MHz (%+.1e); the self-consistent modes %.5g MHz (dispersion of ΔL, ΔM). Largest share\n",
			bt.DfFrozen/1e6, bt.DfMeasured/bt.DfFrozen-1, bt.DfMatrix/1e6)
		fmt.Printf("      in ring 1 %.3f at %.2f periods (the beat period 1/Δf = %.2f periods); at the end %.3g of the energy is left, %.4f\n",
			bt.MaxTransfer, bt.TMaxTransfer, (bt.FIn+bt.FAnti)/2/bt.DfMatrix, bt.EnergyLeft, bt.DarkShare)
		fmt.Printf("      of it in the antiphase pattern (Q_in %.3g, Q_anti %.4g: the bright half leaves, the dark half stays); ledger %.1e,\n",
			bt.QIn, bt.QAnti, bt.Run.LedgerRel)
		fmt.Printf("      books %.1e of the start.\n", bt.Run.Mismatch)
	}
	fmt.Println("    Standard prediction: this is the classical physics of two coupled resonators (normal-mode splitting, beats,")
	fmt.Println("    collective radiation: Dicke super- and subradiance; R. H. Dicke, Phys. Rev. 93, 99, 1954). Neutrinos: no measured")
	fmt.Println("    neutrino–neutrino interaction exists; the Standard Model gives only Z exchange (9e).")

	// 9b
	fmt.Println("9b. Force, time-averaged: F = ⟨I1·I2⟩·∂(M + ΔM)/∂d (energy method at constant currents; the induction field is the")
	fmt.Println("    whole electric near field, so there is no separate electric term), in units of U/a (U the pair's stored energy);")
	fmt.Println("    + = repulsion. Analytic for a mode: ±∂M/∂d/(2(L ± M)); kronspice: ⟨I1I2⟩/⟨U⟩ over 3 periods of each mode run.")
	fmt.Println("    geometry      d/a   ∂M/∂d [H/m]   in phase (analytic, kron)   antiphase (analytic, kron)   random phases: first period, last period")
	for _, x := range []struct {
		g  Geom
		dd float64
	}{{Coaxial, 0.75}, {Coaxial, 1}, {Coaxial, 1.5}, {Coaxial, 2}, {Coaxial, 3}, {Coplanar, 2.5}, {Coplanar, 4}} {
		f, err := Force(e, x.g, x.dd*e.A, 30, 9001)
		if err != nil {
			return err
		}
		p.Forces = append(p.Forces, f)
		p.ForceGeoms = append(p.ForceGeoms, x.g)
		fmt.Printf("    %-12s %5.2f %+12.4e   %+9.4f %+9.4f         %+9.4f %+9.4f         %+.4f ± %.4f, %+.4f ± %.4f (energy left %.3f)\n",
			x.g, x.dd, f.DMdd, f.FInA, f.FInK, f.FAntiA, f.FAntiK, f.FRandEarly, f.FRandEarlySigma, f.FRandLate, f.FRandLateSigma, f.LateEnergyLeft)
	}
	f1 := p.Forces[1]
	fmt.Println("    Coaxial: in phase attracts (parallel currents), antiphase repels; side by side the signs swap (M < 0). Random")
	fmt.Printf("    relative phases (32 pairs, seeded): no mean force at the start (%+.4f ± %.4f at d = a), but the pairs do not stay\n", f1.FRandEarly, f1.FRandEarlySigma)
	fmt.Printf("    incoherent: the bright (in-phase) half radiates away and what is left after 30 periods is mainly the dark antiphase\n")
	fmt.Printf("    pattern, which repels (%+.4f, %.2f of the coherent antiphase value). Ledger of the random runs ≤ %.1e.\n", f1.FRandLate, f1.CoherentInRatio, f1.RandomMaxMismach)
	fmt.Println("    Standard prediction: two neutral classical resonators exchange forces through their fields only when their phases")
	fmt.Println("    are correlated (optical binding, Burns, Fournier & Golovchenko, Phys. Rev. Lett. 63, 1233, 1989). Neutrinos: the")
	fmt.Println("    weak force is a contact (Z-exchange, range ħ/M_Z c = 2.2e-18 m) interaction, unmeasured between neutrinos (9e).")

	// 9c
	fmt.Println("9c. Subradiance: Q of both modes against d (9a table; coaxial).")
	fmt.Printf("    One ring: Q %.4g (radiation alone %.4g). Antiphase pair best with dielectric loss at d = %.2ga: Q %.4g (×%.3g);\n",
		p.Single.Q, p.Single.QRad, p.Best.D/e.A, p.Best.Anti.Q, p.Best.QAntiOverSingle)
	fmt.Printf("    radiation alone best at d = %.2ga: Q_rad %.4g (×%.3g). In phase at d = a: Q %.4g (superradiant: R11 + R_m ≈ 2R11, partly\n",
		p.BestRad.D/e.A, p.BestRad.Anti.QRad, p.BestRad.QRadAntiOverSingl, p.Pairs[2].In.Q)
	fmt.Println("    offset by the lower frequency, R ∝ ω⁴).")
	fmt.Println("    Scaling: for a far coaxial pair R11 − R_m = R11·(kd)²/10 + …, so Q_rad(anti)/Q_rad(single) ≈ 10/(kd)² at fixed frequency:")
	for _, r := range p.Pairs {
		if r.G != Coaxial || r.D > 3*e.A {
			continue
		}
		kdA := r.Anti.W * r.D / C0
		fmt.Printf("      d = %4.2fa: kd at f_anti %.3f, 10/(kd)² = %7.4g, measured Q_rad ratio %7.4g, f_anti/f0 = %.3f\n",
			r.D/e.A, kdA, 10/(kdA*kdA), r.QRadAntiOverSingl, r.Anti.F/p.Single.F)
	}
	c0 := p.Pairs[0]
	fmt.Printf("    The gain stops growing below d ≈ %.2ga: close rings couple strongly (k = %.2f at contact, d = 2r), the antiphase mode\n", p.BestRad.D/e.A, c0.ML)
	fmt.Printf("    moves up in frequency (f_anti/f0 = %.2f there) and R11 ∝ ω⁴ grows faster than the cancellation improves; with\n", c0.Anti.F/p.Single.F)
	fmt.Println("    tan δ = 1e-4 the dielectric (Q_d = 1/tan δ) caps it as well.")
	if field {
		tf := time.Now()
		if fieldCSV != "" {
			var err error
			p.Field, p.FieldLabel, p.FieldPred, p.FieldPredF, p.FieldSeconds, err = readPairFields(fieldCSV)
			if err != nil {
				return err
			}
			fmt.Println("    (field rows from an earlier -only 9 -field run, read back from its p9_field.csv with -pairfieldcsv)")
		} else {
			p.Field, p.FieldLabel, p.FieldPred, p.FieldPredF = runPairFields(e, cells)
			p.FieldSeconds = time.Since(tf).Seconds()
		}
		reportPairFields(p)
	} else {
		fmt.Println("    Field solver: run with -only 9 -field (about an hour); see run_output.txt for the recorded run.")
	}

	// 9d
	p.ChainD = e.A
	fmt.Printf("9d. Chains: N coaxial rings %ga apart. Darkest mode (largest Q) against the brightest; matrix, and kronspice for N ≤ 16\n", p.ChainD/e.A)
	fmt.Println("    (darkest started alone, without dielectric loss, so its Q is the radiation Q; 200 periods at T/100).")
	fmt.Println("    N    bright f GHz  Q      dark f GHz  Q (tan δ 1e-4)  Q_rad (exact, no loss)  kronspice Q_rad  ledger   net moment")
	eL := e
	eL.TanD = 0
	for _, n := range []int{1, 2, 4, 8, 16, 32, 64} {
		a := NewArray(e, Coaxial, p.ChainD, n)
		var row ChainRow
		row.N = n
		if n == 1 {
			row.Bright, row.Dark = p.Single, p.Single
		} else {
			row.Bright, row.Dark, _ = a.ChainModes()
		}
		al := NewArray(eL, Coaxial, p.ChainD, n)
		row.DarkRad = al.SolveMode(realVec(row.Dark.V))
		if n <= 16 {
			r, err := al.RunMode(row.DarkRad, 200, 100)
			if err != nil {
				return err
			}
			row.QDarkK, row.LedgerK, row.MismK, row.Kron = r.Q, r.LedgerRel, r.Mismatch, true
		}
		p.Chain = append(p.Chain, row)
		kq := "      —"
		if row.Kron {
			kq = fmt.Sprintf("%11.5g", row.QDarkK)
		}
		led := "   —"
		if row.Kron {
			led = fmt.Sprintf("%.0e", row.LedgerK)
		}
		fmt.Printf("    %-4d %9.4f %8.4g %10.4f %12.5g %20.5g %s       %s    %.1e\n", n, row.Bright.F/1e9, row.Bright.Q, row.Dark.F/1e9, row.Dark.Q,
			row.DarkRad.Q, kq, led, row.Dark.Dipole)
	}
	// exponent over the largest three N
	var xs, ys []float64
	for _, r := range p.Chain[len(p.Chain)-3:] {
		xs, ys = append(xs, math.Log(float64(r.N))), append(ys, math.Log(r.DarkRad.Q))
	}
	p.ChainExp = fitSlope(xs, ys)
	last := p.Chain[len(p.Chain)-1]
	p.ChainCoef = last.DarkRad.Q / math.Pow(float64(last.N), p.ChainExp)
	b := DefaultBench().Run()
	p.NNeeded = math.Pow(b.QNeeded/p.ChainCoef, 1/p.ChainExp)
	fmt.Printf("    The darkest mode is the band-edge mode (neighbours opposite, smooth envelope, net moment 0). Q_rad ∝ N^%.2f over N = 16–64:\n", p.ChainExp)
	fmt.Println("    it compounds as a power of N, not exponentially (the darkest finite-chain modes of sub-wavelength arrays scale as")
	fmt.Println("    N⁻³: A. Asenjo-Garcia et al., Phys. Rev. X 7, 031024, 2017). The brightest mode's Q falls toward the in-phase limit.")
	fmt.Printf("    With dielectric loss the chain saturates at Q ≈ %.4g (the dielectric, Q_d·f/f0).\n", last.Dark.Q)
	fmt.Printf("    Needed (Part 2): Q = %.3g for SN1987A's proper time, i.e. a suppression %.1e of the single ring's radiation (Part 4).\n", b.QNeeded, p.Single.QRad/b.QNeeded)
	fmt.Printf("    At the bench ring's scaling Q_rad ≈ %.3g·N^%.2f that is N ≈ %.2g rings in one coherent band-edge mode, a chain %.2g m\n",
		p.ChainCoef, p.ChainExp, p.NNeeded, p.NNeeded*p.ChainD)
	fmt.Println("    long, and lossless to 1e-19. For rings in vacuum (εr = 1) there is no subradiance to use: each ring is electrically")
	fmt.Println("    large at its own frequency (ka ≈ 4, Part 2), so neighbours sit beyond kd ≈ 1 and the cancellation is gone.")
	fmt.Println("    Standard prediction: subradiance is limited by size and disorder (a random frequency spread mixes the dark mode with")
	fmt.Println("    bright ones); one neutrino is one particle, not a phase-locked chain of 10⁵–10⁶ rings.")

	// 9e
	p.Nu = nuPair(b)
	nu := p.Nu
	fmt.Println("9e. Against measurement: two model rings at the Part 2 neutrino scale (0.1 eV, a/r = 4, the effective εr that holds")
	fmt.Printf("    the ring for SN1987A's proper time): a = %.3g m, ka = %.3g, εr = %.3g. The coupling k(d/a) depends only on shape\n", nu.A, nu.KA, nu.EpsR)
	fmt.Printf("    (retardation negligible at this ka): coaxial at d = 2a, k = %.4f, so each mode moves by ½ħω0·k = %.3g eV; side by\n", nu.Kappa2a, nu.Shift2a)
	fmt.Printf("    side with a 0.5a gap k = %.4f (%.3g eV). Mutual radiation R_m/R11 = %.8f at 2a: the antiphase pair's radiation Q\n", nu.KappaSide, nu.ShiftSide, nu.Rm2a)
	fmt.Printf("    is %.3g× the single ring's, so a touching antiphase pair (d = 2r) would need only εr ≈ %.3g instead of %.3g.\n", nu.QPairGain, nu.EpsPair, nu.EpsR)
	fmt.Println("    Standard: neutrinos interact with each other only by Z exchange. Its potential is a Yukawa of range")
	fmt.Printf("    ħ/M_Z c = %.3g m: at r = 2a it carries e^{−%.3g} — zero for every purpose. Even the contact interaction √2·G_F(ħc)³\n", lamZ, nu.YukawaExp)
	fmt.Printf("    spread over one model ring's volume gives %.2g eV, %.2g× below the model's coaxial shift at 2a.\n", nu.ContactV, nu.Shift2a/nu.ContactV)
	fmt.Printf("    Cross-sections: weak νν scattering ~ G_F²s/π = %.2g cm² for two 10 MeV neutrinos head on (s = 4E², order of\n", nu.SigmaWeak)
	fmt.Println("    magnitude; never measured in a laboratory). Two model rings passing at c exchange energy through the near-field")
	fmt.Println("    coupling with a transfer phase φ(b) = ∫½ω0 k dt; for point dipoles side by side φ = π ka a²/(4Λ b²) (b the impact")
	fmt.Printf("    parameter, Λ = ln(8a/r) − 7/4): φ(a) = %.2g, π/2 only at b = %.2g a (inside the ring), and σ = ∫φ² 2πb db (b ≥ a) =\n", nu.PhaseContact, nu.BHalf)
	fmt.Printf("    π a²(π ka/4Λ)² = %.2g cm², %.2g× the weak value. Two model neutrinos would interact far more strongly than\n", nu.SigmaModel, nu.SigmaModel/nu.SigmaWeak)
	fmt.Println("    measured ones are allowed to: they would stay in thermal contact in the early universe long after the observed")
	fmt.Println("    decoupling at ~1 s (T ~ 1 MeV; the effective number of species N_eff = 2.99 ± 0.17, Planck 2018, fits the standard")
	fmt.Println("    decoupling), and a coherent pair would attract or repel with forces ~0.1·U/a. Collective flavour conversion in")
	fmt.Println("    supernovae comes from forward νν scattering with potential √2·G_F·n_ν (Pantaleone, Phys. Lett. B 287, 128, 1992),")
	fmt.Println("    not from any EM coupling. The rings are classical: two in-phase rings simply add (no exclusion), while neutrinos")
	fmt.Println("    are fermions (Pauli exclusion; the Tremaine–Gunn bound, Phys. Rev. Lett. 42, 407, 1979, uses it). The model")
	fmt.Println("    gives none of this; its neutral-current analogue is missing as in Parts 2–8.")

	fmt.Println("    Status H42: bounded. Two rings interact as two coupled resonators: splitting k·f0, beats at Δf (kronspice to the")
	fmt.Println("    matrix), forces only for correlated phases (coaxial in phase attracts, antiphase repels); an antiphase pair is")
	fmt.Printf("    subradiant (Q ×%.3g at best, scaling 10/(kd)² far out), a chain's band-edge mode as N^%.1f. The confinement Part 2 needs\n",
		p.BestRad.QRadAntiOverSingl, p.ChainExp)
	fmt.Printf("    is still out of reach (~%.1g phase-locked rings), and the interaction is %.0e× the weak one in cross-section.\n", p.NNeeded, nu.SigmaModel/nu.SigmaWeak)
	p.Seconds = time.Since(t0).Seconds()
	fmt.Printf("   (Part 9: %s)\n", time.Duration(p.Seconds*1e9).Round(time.Second))
	if viz != "" {
		if err := writeViz9(viz, p); err != nil {
			return err
		}
		fmt.Println("\nviz data for Part 9 written to", viz, "— run python3 viz/render_part9.py")
	}
	return nil
}

func realVec(v []complex128) []float64 {
	out := make([]float64, len(v))
	for i := range v {
		out[i] = real(v[i])
	}
	return out
}

// nuPair evaluates test 9e at the Part 2 bench ring.
func nuPair(b BenchResult) NuPair {
	bench := DefaultBench()
	var n NuPair
	n.A, n.KA, n.Omega0, n.EpsR = b.RadiusNeeded, b.KANeeded, b.Omega0, b.EpsRNeeded
	ring := ERing{A: b.RadiusNeeded, R: b.RadiusNeeded / bench.AoverR, EpsR: b.EpsRNeeded, N: 1}
	co := NewArray(ring, Coaxial, 2*ring.A, 2)
	m, rm := co.Coupling(0, 1, b.Omega0)
	n.Kappa2a = m / co.LSelf(b.Omega0)
	evOmega := Hbar * b.Omega0 / E
	n.Shift2a = 0.5 * evOmega * n.Kappa2a
	n.Rm2a = rm / co.R11(b.Omega0)
	sd := NewArray(ring, Coplanar, 2.5*ring.A, 2)
	ms, _ := sd.Coupling(0, 1, b.Omega0)
	n.KappaSide = ms / sd.LSelf(b.Omega0)
	n.ShiftSide = 0.5 * evOmega * math.Abs(n.KappaSide)
	_, anti := co.PairModes()
	one := NewArray(ring, Coaxial, ring.A, 1).SolveMode([]float64{1})
	n.QPairGain = anti.QRad / one.QRad
	// a touching antiphase pair (d = 2r): Q_pair ≈ Q_rad·10/(kd)², Q_rad = 6Λ/(π(ka)³), solve for ka, then εr
	lam := ring.Lambda()
	dOverA := 2 / bench.AoverR
	ka5 := 60 * lam / (math.Pi * dOverA * dOverA * b.QNeeded)
	kaP := math.Pow(ka5, 0.2)
	n.EpsPair = 2 * bench.AoverR * bench.AoverR / (lam * kaP * kaP)
	n.YukawaExp = 2 * ring.A / lamZ
	vol := 2 * math.Pi * math.Pi * ring.A * ring.R * ring.R / 1e-45 // fm³
	n.ContactV = math.Sqrt2 * gF * math.Pow(hbarC, 3) / vol * 1e9   // eV
	s := 4 * math.Pow(bench.EnergyEV/1e9, 2)
	n.SigmaWeak = gF * gF * s / math.Pi * gev2mb * 1e-27
	c := math.Pi * n.KA / (4 * lam)
	n.PhaseContact = c
	n.BHalf = math.Sqrt(2 * c / math.Pi)
	n.SigmaModel = math.Pi * ring.A * ring.A * c * c * 1e4
	return n
}

// runPairFields runs the field-solver configurations of 9c and returns
// them with labels and the circuit's radiation-only predictions.
func runPairFields(e ERing, cells []float64) (res []PairFieldResult, labels []string, pred, predF []float64) {
	logf := func(format string, args ...any) { fmt.Printf("    [field] "+format+"\n", args...) }
	type job struct {
		label   string
		d       float64
		w       []float64
		cells   float64
		halfA   float64
		periods float64
	}
	var jobs []job
	for ci, c := range cells {
		jobs = append(jobs, job{"one ring", e.A, []float64{1}, c, 5, 20})
		ds := []float64{0.75, 1, 1.5, 2}
		if ci > 0 {
			ds = []float64{1}
		}
		for _, dd := range ds {
			jobs = append(jobs, job{fmt.Sprintf("in phase d=%ga", dd), dd * e.A, []float64{1, 1}, c, 5, 20})
			jobs = append(jobs, job{fmt.Sprintf("antiphase d=%ga", dd), dd * e.A, []float64{1, -1}, c, 5, 60})
		}
		if ci == 0 {
			// absorber distance check, and a four-ring band-edge chain
			jobs = append(jobs, job{"one ring (absorber at 8a)", e.A, []float64{1}, c, 7, 20})
			jobs = append(jobs, job{"antiphase d=1a (absorber at 8a)", e.A, []float64{1, -1}, c, 7, 60})
			jobs = append(jobs, job{"chain N=4 band edge d=1a", e.A, nil, c, 5, 80})
		}
	}
	for _, j := range jobs {
		w := j.w
		var ar *Array
		if w == nil { // band-edge weights from the lossless circuit's darkest mode
			el := e
			el.TanD = 0
			ar = NewArray(el, Coaxial, j.d, 4)
			_, dk, _ := ar.ChainModes()
			w = realVec(dk.V)
		} else {
			el := e
			el.TanD = 0
			ar = NewArray(el, Coaxial, j.d, len(w))
		}
		pf := DefaultPairField(j.d, w)
		pf.Cells, pf.HalfA, pf.Periods = j.cells, j.halfA, j.periods
		pf.Log = logf
		r, err := pf.Run()
		if err != nil {
			fmt.Println("    field error:", err)
			continue
		}
		m := ar.SolveMode(w)
		res = append(res, r)
		labels = append(labels, j.label)
		pred = append(pred, m.Q)
		predF = append(predF, m.F)
	}
	return
}

func reportPairFields(p Part9) {
	fmt.Println("    Field solver (em/fdtd.go): tori of Part 5b (εr = 100, lossless) on a common axis, each started with Part 5b's D_φ")
	fmt.Println("    bump times its weight, B = 0; the mode's f and Q from the demodulated E_φ probe at (a, 0, z0) after 4 periods")
	fmt.Println("    (the strongest peak in 1–4.5 GHz, so a slower higher mode is filtered out); circuit = the same weights, no loss.")
	fmt.Println("    configuration                      grid  nodes      f GHz   Q field   circuit f, Q_rad   field/one ring  circuit/one ring  ledger")
	// each row against the single ring on the same grid and in the same box
	box := func(l string) string {
		if i := strings.Index(l, "(absorber"); i >= 0 {
			return l[i:]
		}
		return ""
	}
	type key struct {
		cells float64
		box   string
	}
	oneF, oneC := map[key]float64{}, map[key]float64{}
	for i, r := range p.Field {
		if strings.HasPrefix(p.FieldLabel[i], "one ring") {
			k := key{r.Cells, box(p.FieldLabel[i])}
			oneF[k], oneC[k] = r.Q, p.FieldPred[i]
		}
	}
	gain := map[string]float64{}
	for i, r := range p.Field {
		k := key{r.Cells, box(p.FieldLabel[i])}
		gf, gc := r.Q/oneF[k], p.FieldPred[i]/oneC[k]
		gain[fmt.Sprintf("%s r/%g", p.FieldLabel[i], r.Cells)] = gf
		fmt.Printf("    %-34s r/%g  %8d  %7.4f  %8.4g   %6.4f %8.4g    %8.4g        %8.4g      %.0e\n",
			p.FieldLabel[i], r.Cells, r.Nodes, r.F/1e9, r.Q, p.FieldPredF[i]/1e9, p.FieldPred[i], gf, gc, r.Residual)
	}
	fmt.Printf("    (%.0f s.) Ratios are to the single ring on the same grid and in the same box. The single ring's own Q depends on\n", p.FieldSeconds)
	fmt.Printf("    the box: absorber at 3a from the axis (Part 5b) 41.7, at 6a %.4g, at 8a %.4g (r/2), at 6a on r/3 %.4g. At ka ≈ 0.4\n",
		oneF[key{2, ""}], oneF[key{2, "(absorber at 8a)"}], oneF[key{3, ""}])
	fmt.Println("    the absorber sits in the ring's reactive near field and returns a little of what it takes, so the absolute field")
	fmt.Println("    Q is uncertain by about ±10 %, Part 5b's 41.7 included. The subradiant gain does not move with it: antiphase d = a")
	fmt.Printf("    gains ×%.3g (r/2, 6a), ×%.3g (r/2, 8a), ×%.3g (r/3), against ×21.1 in the circuit; d = 0.75a ×%.3g (circuit ×22.3);\n",
		gain["antiphase d=1a r/2"], gain["antiphase d=1a (absorber at 8a) r/2"], gain["antiphase d=1a r/3"], gain["antiphase d=0.75a r/2"])
	fmt.Printf("    the four-ring band-edge chain ×%.3g (circuit ×54.0). Field frequencies lie 2.5–5 %% above the lumped ones (r/3 closer), as for one ring.\n",
		gain["chain N=4 band edge d=1a r/2"])
}

// readPairFields reads the field rows that writeViz9 saved (p9_field.csv),
// so the report can be printed again without the hour of field runs.
func readPairFields(path string) (res []PairFieldResult, labels []string, pred, predF []float64, secs float64, err error) {
	f, err := os.Open(path)
	if err != nil {
		return
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		return
	}
	col := map[string]int{}
	for i, h := range rows[0] {
		col[h] = i
	}
	num := func(r []string, k string) float64 { v, _ := strconv.ParseFloat(r[col[k]], 64); return v }
	e := DefaultERing()
	for _, r := range rows[1:] {
		res = append(res, PairFieldResult{Cells: num(r, "cells"), D: num(r, "d_over_a") * e.A, Nodes: int(num(r, "nodes")),
			F: num(r, "f"), Q: num(r, "q"), QEnergy: num(r, "q_energy"), Residual: num(r, "residual"), Seconds: num(r, "seconds")})
		labels = append(labels, r[col["label"]])
		pred = append(pred, num(r, "circuit_q"))
		predF = append(predF, num(r, "circuit_f"))
		secs += num(r, "seconds")
	}
	return
}
