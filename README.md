# kronspice

**A SPICE-compatible circuit simulator you embed in Go programs.** Run
thousands of circuits in parallel, in memory and safely, with results
checked against ngspice.

## Welcome, universal engineer!

kronspice is named after Gabriel Kron (1901–1968), the General Electric engineer
who treated a circuit as a geometric object. In two 1956 articles he gave a name to
where he thought this led: **universal engineering**.

Kron looked at an electric machine, a bridge or an airplane and saw the same thing:

> a collection of numerous multi-dimensional spaces connected together into one
> unit

So he looked for one method that carries across them. Describe the parts on their
own, connect them by a change of reference frame, and tear large systems into
pieces that can be solved separately and rejoined exactly. He took those ideas from
motors to power grids, electromagnetic fields, elastic structures, reactors and the
Schrödinger equation.

If you think a problem in one field might be a network problem in disguise, this
project is for you, whatever your field and whether or not you have ever drawn a
circuit. Electrical engineers, mechanical engineers, physicists, chemists, students
and curious tinkerers are all welcome here as universal engineers.

**Try it across domains in kronspice:**

| Domain | As a network | Example |
|---|---|---|
| Circuits | Torn, solved piece by piece and rejoined inside a running evolution | [examples/torn_ga](examples/torn_ga) |
| Rotating machines | Kron's own reference frames, with an induction motor starting up | [examples/kron_machine](examples/kron_machine) |
| Mechanics and acoustics | Generator → piezo → flask, with force and pressure as network variables | [examples/acoustic_drive](examples/acoustic_drive) |
| Molecules | Bond vibrations of water as AC resonances | [examples/water_vibration](examples/water_vibration) |
| Quantum mechanics | Kron's circuit model of the Schrödinger equation | [examples/kron_schrodinger](examples/kron_schrodinger) |

Kron's program was never finished, and neither is this one: his late work on
multi-dimensional networks is not implemented here. To go deeper, read
[the history of universal engineering](HISTORICAL.md), and use
[the nomenclature guide](NOMENCLATURE.md) to translate Kron's terms into modern
ones.

## Why kronspice

At its core, kronspice is a circuit simulator. Tools that design
circuits automatically (genetic algorithms, optimizers, AI-assisted design, cloud
design services) need a simulator inside their loop. Established simulators are
separate programs built for one engineer and one netlist at a time. Driving them
means writing files, starting processes and parsing output for every candidate, and
generated candidates are often broken in ways that stall them. kronspice is a
library built for that loop.

| | |
|---|---|
| **SPICE-compatible** | 36 test circuits agree with ngspice-39 within the stated tolerances: DC, AC and transient; diode, BJT, MOSFET levels 1–3, coupled inductors. Four charge- and flux-defined C/L circuits have ngspice twins too. |
| **Safe on broken circuits** | Floating nodes, shorted sources, singular matrices and non-convergence return typed errors, never a panic or a hang. Every analysis takes a `context` deadline. |
| **Parallel and deterministic** | No shared mutable state, and the full suite passes Go's race detector. The same netlist, options and seed give bit-identical results. |
| **Pure Go** | No cgo, no external program, no third-party dependencies (`go.mod` has none). Cross-compiles everywhere. |
| **Built-in evolution** | Value and topology genomes, weighted-sum and NSGA-II multi-objective fitness, checkpoints that resume exactly. Best circuits are written as plain SPICE netlists (checked in ngspice; LTspice and QSpice not yet tested). |
| **Kron diakoptics** | Tears a circuit into pieces, solves them separately and joins them. This agrees with the full solve to 1e-9, and power is the same torn and untorn. A mutation re-solves only its piece. |
| **More than circuits** | Mechanical and acoustic networks join the circuit through transducers, and machines run as moving inductors with a mechanical port, all on the same energy ledger. |
| **Every joule accounted for** | A per-element, per-channel energy ledger reports energy appearing from nowhere as an error, so evolution cannot exploit a simulator flaw. |
| **KiCad export** | Any netlist becomes a wired KiCad 8 schematic plus a results sheet. |

## Install

```bash
go get github.com/fdamador/kronspice
```

Go 1.24 or later. Command-line tools and container images are attached to
each tagged release.

## Quick start

```go
c, err := circuit.NewBuilder("divider").
    VoltageSource("V1", "in", "0", 10).
    Resistor("R1", "in", "out", 1e3).
    Resistor("R2", "out", "0", 1e3).
    Build()
if err != nil {
    log.Fatal(err)
}
r, err := analysis.OP(c, analysis.Options{})
if err != nil {
    log.Fatal(err)
}
v, _ := r.Voltage("out")
fmt.Println(v, r.Ledger.Inflow) // 5 0.05 (volts at out, watts from V1)
```

Evaluate candidates in parallel. Each goroutine owns its circuit, and a bad
candidate gets a penalty rather than a crash:

```go
fitness := make([]float64, len(candidates))
var wg sync.WaitGroup
for i, r2 := range candidates {
    wg.Add(1)
    go func() {
        defer wg.Done()
        ctx, cancel := context.WithTimeout(context.Background(), time.Second)
        defer cancel()
        c, err := divider(r2).Build()
        if err != nil {
            fitness[i] = penalty // invalid circuit
            return
        }
        r, err := analysis.OP(c, analysis.Options{Context: ctx})
        if err != nil {
            fitness[i] = penalty // singular, no convergence or timeout
            return
        }
        v, _ := r.Voltage("out")
        fitness[i] = math.Abs(v - 2.5)
    }()
}
wg.Wait()
```

From the command line:

```bash
go run ./cmd/kronspice -o bode.csv testdata/netlists/rlc_bandpass_ac.cir   # simulate a netlist
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json       # evolve R and C for a 1 kHz cutoff
go run ./cmd/cir2kicad examples/water_vibration/netlists/h2o.cir         # KiCad schematic + results sheet
```

## What it simulates

- **Elements:** R, C, L, coupled inductors (K), V and I sources (DC, AC,
  PULSE, SIN, EXP, PWL), controlled sources E/F/G/H, piecewise-linear
  nonlinear elements (B), negative resistance, SPICE3 diode, Gummel–Poon
  BJT, MOSFET levels 1–3, and charge- and flux-defined and time-varying C
  and L.
- **Analyses:** operating point, DC sweep, AC small-signal, transient
  (trapezoidal, with backward Euler at start and breakpoints; adaptive
  step from truncation error), measurements (gain, −3 dB bandwidth, rise
  time, overshoot, settling time), IEEE 1459 power quantities, and
  geometric-algebra and quaternion power theory.
- **Solvers:** dense LU and sparse Gilbert–Peierls LU, both native. Newton
  with gmin stepping, source stepping and pseudo-transient continuation.
- **Netlists:** SPICE syntax with `.op`, `.dc`, `.ac`, `.tran`, `.model`
  and engineering suffixes.

**Not yet:** `.subckt`, `.param` and `.include`; temperature other than
27 °C; switches, transmission-line devices, JFETs and BSIM models; noise and
pole–zero analysis. For those, use [ngspice](https://ngspice.sourceforge.io/),
which kronspice uses as its reference. No speed comparison with ngspice has
been measured yet.

## Verified against ngspice

The tools container (Go, gcc, ngspice) generates ngspice references for
every test netlist, and `go test` compares kronspice with them. That check is
alongside exact and independent references (rational-arithmetic nodal
analysis, analytic step responses, Python re-implementations of the device
equations).

```bash
docker build -t kronspice-tools internal/golden/docker
docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_all.sh
docker run --rm -v "$PWD":/src kronspice-tools go test -race ./...
```

On Windows, `internal/golden/docker/run.ps1` and the VS Code tasks in
`.vscode/tasks.json` run the same steps.

## Learn more

- [HISTORICAL.md](HISTORICAL.md): Gabriel Kron, his method and universal engineering.
- [NOMENCLATURE.md](NOMENCLATURE.md): Kron's terms in modern electrical engineering language.
- [doc/STATUS.md](doc/STATUS.md): every package, milestone and command in detail.
- [CLAUDE.md](CLAUDE.md): design, numerical rules, references and milestones.
- [examples/](examples/): runnable examples, each with an `export/` folder
  holding its netlist, KiCad schematic and Kron results sheet.
- **Research test bench.** kronspice also has field solvers in Maxwell's
  quaternion form, and test benches for proposals beyond standard physics.
  Each is labeled hypothetical and reports the standard-physics prediction
  beside its own. See CLAUDE.md, "Hypotheses under test".

## License

Dual-licensed: **GPL-3.0-only** ([LICENSE](LICENSE)) or a **commercial
license** for proprietary use ([COMMERCIAL-LICENSE.md](COMMERCIAL-LICENSE.md)).
See [LICENSING.md](LICENSING.md) for which applies, and
[CONTRIBUTING.md](CONTRIBUTING.md) for the contributor agreement.

Copyright (C) 2026 Franklin Amador.
