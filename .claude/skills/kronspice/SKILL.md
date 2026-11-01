---
name: kronspice
description: Simulate, check or evolve an analog circuit with kronspice (SPICE netlist in, JSON or text report out, with the energy ledger), run one of its examples, or interpret its output and errors. Use when asked to simulate a circuit, compute an operating point, sweep, Bode plot or transient, evolve component values or a topology, or read a kronspice ledger.
---

# kronspice

kronspice is this repository's SPICE-compatible simulator (pure Go). Use
the CLI for one-off simulations and the Go API inside programs. Everything
below runs from the repository root.

## Simulate a netlist

1. Write a SPICE netlist: title line, elements, `.model` cards, analyses,
   `.end`. Ground is node `0`. Engineering suffixes work (`1k`, `4.7u`,
   `10meg`). Supported elements: R, C, L, K, V, I (DC, AC, PULSE, SIN,
   EXP, PWL), E/F/G/H, D, Q, M (levels 1–3), B (`I=pwl(V(a,b)) …` or
   `V=pwl(I(Vx)) …`), and the extensions Y (gyrator), N (Kron machine),
   P (bearing motor). Analyses: `.op`, `.dc`, `.ac`, `.tran`. Not
   supported: `.subckt`, `.param`, `.include`, switches, JFETs, BSIM,
   transmission-line devices, noise, pole-zero.
2. Run it and read JSON, which is the reliable interface:

   ```
   go run ./cmd/kronspice -json circuit.cir
   printf '...netlist...' | go run ./cmd/kronspice -json -
   ```

   The text report (no `-json`) is for humans; `-o out.csv` writes the
   last analysis as CSV; `-solver dense|sparse` forces a backend.
3. Check `ok`. If false, branch on `error.class`:
   - `syntax`: fix the netlist line named in `error.message`.
   - `invalid`: no ground, a floating node (every node needs a path to
     ground), a bad value or model, a bad analysis spec, or an
     electrically large circuit.
   - `singular`: a voltage-source loop, a capacitor-only cut set with no
     DC path, or an inductor loop. Add a path or a small resistance.
   - `no-convergence`: Newton or the time step failed. Try a smaller
     `.tran` step or a gentler source edge; from the Go API, a starting
     guess through `analysis.Options.NodeSet`.
   - `ledger`: the energy balance failed. This is a bug to find, not a
     result to report.
   - `timeout`, `error`: see the message.
4. Read results by name: `analyses[k].op.nodes["out"]`,
   `analyses[k].tran.v["out"]` aligned with `tran.time`,
   `analyses[k].ac.v["out"].mag` and `.phase_deg` aligned with `ac.freq`,
   `analyses[k].dc.v["out"]` aligned with `dc.values[0]`.

## Read the ledger honestly

Every run has an energy ledger. `residual` (op) or `max_relative_residual`
(dc, ac, tran) must be at round-off. `model_supplied` is power from
idealized elements (negative resistance, controlled sources) and is never
inflow from outside. `hypothetical` is from hypothesis-mode models and is
never counted as real. Never report a surplus as energy gain: an
imbalance means a numerical or modelling error, and CLAUDE.md §3c and §1
say so. When a result tests a register hypothesis (H-numbers in
CLAUDE.md §1), state the standard prediction beside it.

## Use the Go API

```go
c, err := circuit.NewBuilder("name").
    VoltageSource("V1", "in", "0", 10).
    Resistor("R1", "in", "out", 1e3).
    Resistor("R2", "out", "0", 1e3).
    Build()
r, err := analysis.OP(c, analysis.Options{Context: ctx})
wf, err := analysis.Tran(c, analysis.TranSpec{Step: 1e-5, Stop: 5e-3}, analysis.TranOptions{})
ac, err := analysis.AC(c, analysis.ACSpec{Variation: "dec", Points: 10, FStart: 10, FStop: 1e5}, analysis.Options{})
```

Errors are sentinels for `errors.Is`; `evo.Penalty(err)` turns any of them
into a penalty fitness and class. Always pass a context with a deadline
when evaluating generated circuits. Never share a circuit between
goroutines; build one per worker.

## Evolve values or a topology

```
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json
go run ./cmd/evolve -config examples/filter_topology/config.json -workers 4
```

Copy an existing config (examples/*/config.json) and change the netlist,
genes and objectives. Outputs land in `runs/<timestamp>-<name>/`: a copy
of the config, per-generation statistics as CSV, the best netlist per
generation (opens in ngspice, LTspice and QSpice) and `checkpoint.json`
for `-resume`. A seed gives a bit-identical run; quote the seed with any
result.

## Run or add an example

`examples/README.md` indexes every example by group with its scope and
run command; each folder's README has the construction, results and
sources. To add an example: a `main.go` with a package doc comment, tests
that check known values and a closed ledger, an `export_test.go` writing
`export/` (netlist, KiCad schematic, results sheet), a README in the same
shape as its neighbours, and a row in `examples/README.md`. Every new
file starts with the SPDX dual-licence header. Run `go test -race ./...`
before calling it done.

## Where to look

- CLAUDE.md: design rules, Kron's methods, the ledger, the hypothesis
  register, references, milestones.
- llms.txt: the one-page summary of the above for any tool.
- doc/STATUS.md and doc/M*.md: what each package does and the measured
  results behind each milestone.
- testdata/netlists/*.cir: small netlists that are known to run, with
  ngspice twins in testdata/golden.
