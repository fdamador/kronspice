# CLAUDE.md — examples/earth_machine (Earth as an electromagnetic machine)

Notes for Claude sessions working in this example. The repo's root CLAUDE.md
rules apply here too (standard prediction alongside any hypothesis, energy
ledger closes, hypotheses register in NOMENCLATURE/register).

## Goal
Model the Earth as an electromagnetic machine with Kron's machine theory in
kronspice, and write a paper. Three layers, then the links between them:
1. **Geodynamo** — Rikitake two-disk dynamo (two cross-coupled homopolar
   machines). <- DONE, with convergence, (μ, A) and paper107 studies
2. **Global atmospheric electric circuit** — thunderstorms charging the
   ionosphere to about 250 kV through leaky fair-weather air. <- DONE, with
   two-scale profile, land/ocean columns, storm generators
3. **Earth–ionosphere cavity** — Schumann resonances on a transmission line in
   colatitude. <- DONE, with frequency-dependent (two-height) ionosphere
4. **Links** — one network for layers 2 and 3; dynamo → conductivity through
   cosmic-ray shielding; one ledger. <- DONE (first pass)
5. **Stacked design** — the two disks on one axis, each in the other's coil
   (`stacked.go`). <- DONE

Paper (Claude Docs): "The Earth as an Electromagnetic Machine: A Kron-Network
Approach" — https://claude.ai/artifact/KVPScKu6m7MaRCTeu7VEk3. A PDF export,
the figures and a README are kept in `reference/earth_machine/` (local:
`reference/` is git-ignored); re-export there whenever the draft changes.
Sections 1–9,
Figures 0–10 placed (§9 = stacked design) and a coupling diagram, references. Figure numbers in the
paper follow reading order; the repo files keep their build order:
paper Fig. 4 = `fig6_reversal_statistics.png`, Fig. 5 = `fig7_paper107_studies.png`,
Fig. 6 = `fig4_global_circuit.png`, Fig. 7 = `fig5_schumann.png`,
Fig. 8 = `fig8_layers_linked.png` (Figs. 0–3 match).

## Files
- `main.go` — layer 1: two machine elements N1, N2 (shafts as mechanical nodes,
  inertia C, torque G), RK4 twin, reversal statistics, trace books, four cases
  (co-rotating; disk 2 reversed; disk 2 reversed + coil 2 rewound; co-rotating
  + coil 2 rewound = Tesla's opposed fields). Params also carry K (feedback
  coil), GMod/Wd (torque modulation) and Speed (prescribed shaft speeds).
  Flags: `-layer 0|1|2|3|4|5`, `-t`, `-dt`, `-csv`, `-data DIR`.
- `studies.go` — task 2 (convergence, two-sample KS), task 3 ((μ, A), bootstrap
  KS vs exponential, precursor ratio, polarity record), task 4a–e.
- `atmosphere.go` — layer 2: 60 leaky shells (R ∥ C), storm source on the
  thermal channel, calibration σ0 = J/E0 and H, OP, switch-on vs exact ladder,
  Carnegie curve by AC at 1/day.
- `atmosphere2.go` — columns with any profile (numerical ∫dz/σ), two-scale
  profile, land/ocean, storm generators with R_up, R_down, R_cloud.
- `cavity.go` — layer 3: 181-band line in colatitude, Sturm-bisection modes vs
  Legendre, AC sweep (pole, equator), constant-height fit to mode 1.
- `cavity2.go` — two-height ionosphere: complex heights, complex mode solver,
  Nelder–Mead fit, ladder rebuilt per frequency (series R+L, shunt G+C).
- `coupled.go` — the links (Unified network, ElectricHeight, shielding chain).
- `stacked.go` — stacked design: StackedCouplings(gap) from em.MutualCircles /
  em.Nagaoka (homopolar geometry in disk radii), stackedCircuit (coils as
  inductors + Mutual, 1 nH sense windings along the loops, G01 = −M,
  G02 = −(K+KS)·M; loop current = +I(R), unlike side by side), gap and
  unequal-torque studies (`-layer 6`). `stacked_test.go`.
- `viz/plot_fluid.py` (Fig. 9, the machine as two eddies), `viz/plot_stacked.py`
  (Fig. 10, reads `viz/data/stacked_traces.csv`).
- `figdata.go` — CSVs for the figures (`-data viz/data`, about 7 minutes).
- `viz/plot_layers.py` — figures 4–8 from the CSVs (matplotlib; 300 dpi).
  `viz/fig0–3` are the earlier figures.
- `main_test.go`, `layers_test.go`, `studies_test.go` — tests (about 7 s).
- `run_output.txt` — output of `go run ./examples/earth_machine -t 1000`.
- `rikitake_reference.py` — the original SciPy reference.

## Model (layer 1)
`L İ1 = −R I1 + M ω1 I2 + K M ω1 I1`, `L İ2 = −R I2 + M ω2 I1 + K M ω2 I2`,
`C ω̇i = Gi(t) − M I1 I2 − K M Ii²`. All constants 1 ⇒ Rikitake, μ = 1, A = 1.
Disk i = machine element, winding 0 its radial path (L 1e-9), winding 1 the
OTHER loop's coil. The element books ω·G·i as a voltage DROP along a
winding: the cross term G[0][1] = +M reads as a source only because winding 1
runs against its loop, so the aiding self term is G[0][0] = −K·M (checked
against RK4 to 1.4e-3). Loop current = −ElementCurrent("Ri"). At the UIC
point (k = 0) the resistor currents are not yet the loop currents: skip k = 0
when comparing traces.

## Results
Layer 1
- Equilibrium I1² = golden ratio; neutral oscillation, period 4.20 = the floor
  of the intervals. Circuit, t 50–1000: 71 reversals, mean 12.17, CV 0.92.
- Convergence: RK4 t = 20 000, h 0.004 → 0.0005: 1385 → 1473 reversals; KS vs
  finest p = 0.01, 0.65, 0.88; circuit dt 0.02/0.01/0.005: p = 0.13/0.38/0.39.
- (μ, A): exponential rejected everywhere (bootstrap p ≤ 0.0025);
  longest/mean ≤ 6 (record ≳ 100); last swing before reversal 8–23 % above
  median (record: weakening); μ = 0.5, A ≥ 2 periodic.
- 4a reversal rate 0.03–0.14 (0.075 unmodulated), max near 0.8 ωn.
  4b growth (−R + M√(ω1ω2))/L exact. 4c K ≥ 0.125 locks polarity; K < 0 rare
  reversals, shafts driven by I² torque. 4d opposed fields: λ = −1 ± i, decays.
  4e machine = exact under dω/dt ≠ 0 (7e-6); inertial terms ≈ Ω/ω_c ≈ 1e-9.
Layer 2
- σ0 = 1.67e-14 S/m, H = 2.08 km; I = 1020 A, R = 244.9 Ω; Q = −5.42e5 C;
  half the balancing charge below 1.42 km; C_eff 2.17 F vs 0.076 F;
  R·C_eff = ε0/σ0 = 531 s for ANY profile; 255 MW; Carnegie lag 4.4 min.
- Two-scale: half below 0.72 km. Land/ocean: same E0 (108 V/m), J ∝ σ0.
  Storms: 23–72 % of the generator current reaches the ionosphere, 0.57–3.9 GW.
Layer 3
- Ladder = Legendre to 4e-4, second order. Constant heights (fit n = 1):
  modes 2–5 up to 10 % low, Q ∝ f. Two-height fit: all five within 2 %,
  h_e 52.3 km, h_m 90.2 km at 8 Hz, Q2–5 = 4.7–4.2 (reported 4–6).
Links
- One network (61 bands × 40 shells): DC one ionosphere 249.6–250.3 kV;
  ELF peaks 7.2, 13.0, 18.9, 24.9, 30.9 Hz with no fitted h_e (8 % low);
  emergent h_e 45.6 km at 8 Hz.
- Shielding, q ×2 declared: V_I 250 → 177 kV, τ 531 → 376 s with no dipole.

Stacked
- KS = 0.519, KM = 0.387 at gap 0.53 r_o (= homopolar); circuit = RK4 1.2e-3.
- No gap reverses: self-coupling breaks the conservation of A = ω1 − ω2, A
  drains (gap 6: 0.005 by t = 1000) and the field locks. Unequal torques bring
  back only clock-like reversals (gap 0.53, δ 0.5: interval 3.87, CV 0.13).
  Same mechanism explains task 4c (aiding feedback coil locks).

## Engine issue (fixed)
At a UIC start with inconsistent winding ICs, the imposed point's powers
(currents routed through reference resistors) were integrated: −3.75e4 J
mismatch. `analysis/tran.go` now integrates the first interval at its right
end when the first point is imposed (backward-Euler start), with
`analysis/tran_uic_test.go` (fails without the fix: −3e4 J). Mismatch now
1.1e-5 J (4e-9 of the drive work). `Books` kept as an independent check.

## Declared (not measured) inputs — keep labelled
Boundary layer 1.5 km and free-air scale height 5 km; land/ocean σ0 1.0e-14 /
2.2e-14 S/m; storm resistances; unshielded ion production ×2; Q1 = 5 and the
measured Schumann frequencies 7.8, 14.3, 20.8, 27.3, 33.8 Hz (Sentman 1995).

## Next tasks
1. Layer 3: day–night asymmetry (non-zonal modes, 2-D sphere grid); lightning
   source spectrum and the H_φ spectrum a station records; fit to spectral
   peaks rather than mode frequencies; break the ζe/h_e0 trade-off of the fit.
2. Links: a measured upper-stratosphere conductivity profile in the unified
   network (does the emergent h_e reach 52 km?); latitude columns with their
   own cosmic-ray cutoff, so the shielding link acts band by band; a measured
   ion-production-vs-cutoff relation to replace the declared ×2.
3. Layer 1: larger Kron networks (more coupled loops) toward superchrons and
   pre-reversal weakening; Glatzmaier & Roberts (1995) as the benchmark.
4. Paper: review pass by the author; the circuit diagram of the unified
   network (one band's column + series path) as a figure.

## Notes
- Layer 1 is dimensionless: no claim about reversal times in years.
- Figures at paper quality (300 dpi); matplotlib was run from a scratch venv,
  not a repo dependency.
- The full run takes several minutes (studies and the two-height sweep).
