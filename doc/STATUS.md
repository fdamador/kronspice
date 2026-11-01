# kronspice — implementation status and usage details

Moved from the README (2026-10-27) when the README became the project's front page.


A SPICE-compatible analog circuit simulator in pure Go, built as the
fitness-evaluation engine for genetically evolved circuits. The full design,
rules and milestones are in [CLAUDE.md](../CLAUDE.md).

## Status: milestones M1–M9, M9a, M9a.1, M9b, M9c and M9c.1 complete

- `exchange/` — every element carries one external energy port per channel
  (electrical, light, thermal, magnetic, …, neutrino, reserved). Unmodeled
  ports report exactly 0 W and are listed as "not modeled".
- `netlist/` — SPICE netlists: title, `*` and `;` comments, `+` continuation,
  engineering suffixes, R / C / L (with `IC=`), V / I sources with DC, AC,
  PULSE, SIN, EXP and PWL, K (mutual inductance), D (diode), Q (BJT) and M (MOSFET, with instance
  parameters `L= W= AD= …`) with `.model`, E / F / G / H
  controlled sources, B (nonlinear PWL
  elements, ngspice syntax), negative resistance, and `.op`, `.dc`, `.ac`,
  `.tran`, `.title`, `.end`.
- `circuit/` — circuit model, validation (ground, floating nodes, duplicate
  names, bad values and waveforms) and the `Builder` API used by the GA.
- `device/` — resistor, capacitor, inductor, voltage and current sources,
  SPICE3 level-1 diode (series resistance, junction and diffusion charge,
  breakdown, `pnjlim` limiting), voltage- and current-controlled nonlinear
  elements (N-type and S-type), linear controlled sources E, F, G, H,
  Gummel–Poon bipolar transistors (NPN and PNP, charges, area), MOSFETs
  levels 1–3 (Shichman–Hodges, Grove–Frohman, semi-empirical; NMOS and
  PMOS, process parameters, Meyer and junction capacitances, weak
  inversion, velocity saturation; exact conductances by dual numbers), SPICE3
  waveforms with breakpoints, truncation-error estimates.
- `solver/` — real and complex LU, no external dependencies, allocation-free
  solve: dense with partial pivoting, and sparse Gilbert–Peierls with
  minimum-degree ordering and diagonal-preference pivoting; chosen by size
  (dense below 30 unknowns) or forced with `Options.Backend` /
  `kronspice -solver`.
- `analysis/` — Newton–Raphson with gmin stepping, source stepping and
  pseudo-transient continuation, node sets, open-circuit / short-circuit
  port measurement;
  operating point (`.op`), DC sweeps (`.dc`), AC small-signal
  (`.ac`) and transient (`.tran`: trapezoidal / backward Euler, adaptive step
  from local truncation error, exact breakpoints, UIC); measurements (gain,
  −3 dB, rise time, overshoot, settling, crossings).
- `kron/` — Gabriel Kron's tensor analysis of networks: the primitive
  network, the mesh and node-pair connection tensors, and the transformed
  equations `Z′ = Cₜ·Z·C`, `Y′ = Aᵀ·Y·A` with typed reference frames, for
  linear DC and AC; agrees with MNA to 1e−15 on the linear golden circuits.
  Kron's orthogonal frame (closed and open paths), one open path per element
  and channel to each exterior junction with the energy ledger built from
  them (equal to the MNA ledger row by row), TAN chords (mutual and delayed
  coupling) and radiation-resistance edges. Diakoptics: named branches
  are torn out, each piece is factored on its own and the pieces are
  joined through the torn-branch currents (Aitchison's form), with
  factored pieces cached by content and shared between goroutines; equal
  to MNA to 1e−10, and cuts that leave a floating piece are reported.
  Nonlinear circuits: every Newton Jacobian (operating point, DC sweep,
  AC, transient) can be solved torn (`kron.JacobianTearing` as
  `analysis.Options.Solver`), cut at resistors, capacitors or the
  nonlinear devices themselves; all golden circuits match ngspice that way.
  Nonlinear circuits in Kron's own frame: every device a primitive network
  with a state-dependent tensor, joined by a constant connection tensor
  (`kron.NonlinearOP`, `kron.NonlinearAC`; `kron.OP`/`AC` node form);
  equal to MNA on every golden circuit. Kron's machine dynamics: the
  Christoffel-symbol form of the machine equations and rotating
  (non-holonomic) frames (`kron.Christoffel`, `kron.RotatingFrame`).
- `em/` — open-system energy ledger, element × channel, including stored
  energy, with a residual check (Tellegen's theorem).
- `output/` — text reports and CSV (transients resampled on the print step).
- `evo/` — genetic algorithm, value genome (E-series, signed genes, model
  and instance parameters) and topology genome (parts grown into an
  embryo), weighted-sum or NSGA-II multi-objective selection (Pareto
  fronts, constraints), objectives including frequency masks, gain,
  bandwidth, supply power and parsimony, worker pool,
  deterministic seeds, checkpoint and resume, run directories.
- `material/` — cited material constants (Seebeck coefficients,
  resistivities, electrode potentials, emission lines); `em/coil.go` —
  analytic coil and plate models (Wheeler, turn-to-turn and bifilar
  capacitance, skin effect).
- `cmd/kronspice` — command-line simulator; `cmd/evolve` — GA experiments.
- `examples/amplifier_gain` — a 2N3904 amplifier evolved to 30 dB with
  100 Hz–200 kHz bandwidth from standard part values.
- `examples/filter_topology` — a low-pass filter grown from an empty test
  fixture; evolution finds the four-part elliptic (Cauer) low-pass.
- Energy domains: mechanical (velocity, force) and acoustic (pressure,
  volume velocity) networks in the same circuit, ideal transducers (piezo,
  piston), radiation into water, all on the ledger (`.domain`, `Y`, `EXT=`;
  [doc/domains.md](doc/domains.md)); `fluid/` — H₂O and D₂O properties,
  piston radiation impedance, the spherical flask's breathing mode.
- Nonlinear and time-varying capacitors and inductors, defined by charge
  q = s(t)·Q(v) and flux λ = s(t)·Λ(i), so i = dq/dt = C dv/dt + v dC/dt
  and v = dλ/dt = L di/dt + i dL/dt; the work of the variation,
  ṡ·∫₀ˣ Q, is booked on a mechanical port (or `EXT=`):
  `C1 a 0 Q=pwl(v1 q1 …)`, `L1 a 0 FLUX=pwl(i1 λ1 …)`,
  `C2 a 0 1u SCALE=SIN(1 0.1 10k)`; try
  `go run ./cmd/kronspice examples/parametric_pump/pump.cir`.
- `examples/acoustic_drive` — generator → piezo → flask mode, the drive
  chain of a single-bubble acoustic resonator: pressure per volt and the
  power budget for 1.3 atm, in H₂O and D₂O (the bubble itself is out of
  scope for now).
- Moving and time-varying inductors (`N` element, `.model … machine`):
  windings whose inductances depend on a rotor angle, a speed matrix, and
  a mechanical port (torque, angular velocity); a coordinate driven at a
  set speed makes it a parametric inductance, its pump power on the
  ledger. `examples/kron_machine` — an induction motor started in Kron's
  frames I and II, and a pumped LC tank against Mathieu's growth rate
  ([doc/M9c1.md](doc/M9c1.md)).
- `examples/torn_ga` — a pre-filter evolved in front of a fixed plant,
  with AC solved by MNA, by Kron untorn and by Kron torn and cached: same
  best design, the torn run 11× faster than the untorn one
  ([doc/M9b.md](doc/M9b.md)).
- `examples/mos_amplifier_pareto` — an NMOS common-source stage evolved
  with NSGA-II into a 99-design gain / bandwidth / supply-power Pareto
  front (up to 34.5 dB, up to 34 MHz, down to 22 µW).
- `examples/bial_bifilar` — the first experiment; its first run is in
  [doc/experiments/bial_bifilar_run1.md](doc/experiments/bial_bifilar_run1.md)
  (dry and wet) and [run 2](doc/experiments/bial_bifilar_run2.md) (adds a
  Bi–Al thermopile).

Milestone notes: [doc/M1.md](doc/M1.md), [doc/M2.md](doc/M2.md), [doc/M3.md](doc/M3.md), [doc/M4.md](doc/M4.md), [doc/M5.md](doc/M5.md), [doc/M6.md](doc/M6.md), [doc/M7.md](doc/M7.md), [doc/M8.md](doc/M8.md), [doc/M9.md](doc/M9.md), [doc/M9a.md](doc/M9a.md), [doc/M9a1.md](doc/M9a1.md), [doc/M9b.md](doc/M9b.md), [doc/M9c.md](doc/M9c.md), [doc/M9c1.md](doc/M9c1.md); requested future work: [doc/future-work.md](doc/future-work.md).

## Use

```bash
go run ./cmd/kronspice testdata/netlists/divider.cir
go run ./cmd/kronspice -all-channels -o out.csv testdata/netlists/bridge.cir
go run ./cmd/kronspice -o rlc.csv testdata/netlists/rlc_step.cir
go run ./cmd/kronspice -o bode.csv testdata/netlists/rlc_bandpass_ac.cir
go run ./cmd/kronspice -o rect.csv testdata/netlists/rectifier_rc.cir
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json
go run ./cmd/evolve -config examples/negative_resistance_gain/config.json
go run ./examples/bial_bifilar          # first experiment: Bi–Al element ∥ bifilar coil
```

From Go:

```go
c, err := circuit.NewBuilder("divider").
    VoltageSource("V1", "in", "0", 10).
    Resistor("R1", "in", "out", 1e3).
    Resistor("R2", "out", "0", 1e3).
    Build()
r, err := analysis.OP(c, analysis.Options{})
v, _ := r.Voltage("out")       // 5
fmt.Println(r.Ledger.Inflow)   // 0.05 W from V1
```

## Tests

```bash
go test ./...
go test -race ./...
go test -bench=. ./solver ./analysis
go test ./netlist -run=^$ -fuzz=FuzzParse -fuzztime=30s
go test ./analysis -run=^$ -fuzz=FuzzSimulate -fuzztime=30s
```

Golden references in `testdata/golden` are exact: `.op` values from an
independent rational-arithmetic nodal analysis, `.tran` tables from analytic
step responses, `.ac` tables from exact transfer functions. To replace the `.op` references with ngspice results,
install ngspice and run:

```bash
go run ./internal/golden/cmd/regen -netlists testdata/netlists -out testdata/golden
```

The nonlinear and time-varying C and L tables (`nlcap_tran`, `nlind_tran`,
`varcap_tank_tran`, `varind_tran`) come from ngspice: each netlist has an
ngspice twin in `internal/golden/ref/ngspice/` (ngspice's charge-defined
`C … Q='…'`, and for inductors its gyrator dual). A container with Go, gcc
and ngspice runs those references and the race detector on any host:

```bash
docker build -t kronspice-tools internal/golden/docker
docker run --rm -v "$PWD":/src kronspice-tools go test -race ./...
docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_ref.sh
```

## License

kronspice is dual-licensed: **GPL-3.0-only** ([LICENSE](../LICENSE)) or a
**commercial license** for proprietary use ([COMMERCIAL-LICENSE.md](../COMMERCIAL-LICENSE.md)).
See [LICENSING.md](../LICENSING.md) for which applies, and
[CONTRIBUTING.md](../CONTRIBUTING.md) for the contributor agreement.

Copyright (C) 2026 Franklin Amador.
