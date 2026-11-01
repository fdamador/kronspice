# A divider with gain: evolution finds negative resistance

An M5 example (CLAUDE.md §4, value genome) built to show the `sign: either`
gene and the ledger's model-supplied column (CLAUDE.md §3c, "Negative
resistance, duality and sign inversion"). A two-resistor divider is asked
for a voltage gain of 2 (6.02 dB) at 1 kHz. No passive divider can do it:
the only solution is a lower leg of −2 kΩ, so the search has to cross zero
and settle on a negative value.

```bash
go run ./cmd/evolve -config examples/negative_resistance_gain/config.json   # ≈ 2 s
go run ./cmd/kronspice runs/<run>/best.cir                                   # re-simulate the result
```

## Circuit, gene and objective

```
V1 in  0   DC 1 AC 1
R1 in  out 1k
R2 out 0   1k          gene: 10 Ω … 100 kΩ in magnitude, sign either
.ac lin 1 1k 1k
```

| Objective | Target | Error |
|---|---|---|
| gain of V(out)/V(in) at 1 kHz | 6.0206 dB (×2) | \|gain − 6.0206\| in dB / 20 |

With `sign: either` the sign is its own gene, flipped by mutation, so the
log-scale magnitude and the sign evolve separately. GA settings: population
30, 30 generations, seed 5, mutation step 0.3 decades, `target_fitness`
1e-4 (not reached, so the run goes the full 30 generations).

## Result

Seed 5 after 30 generations (842 evaluations): R2 = −2000.86 Ω, a gain of
6.017 dB (error 1.9e-4). The exact answer is R2 = −2 kΩ, since
H = R2/(R1 + R2) = −2000/(1000 − 2000) = 2.

The run also prints

```
model-supplied power 0.000999 W (idealized active elements: a real circuit needs a supply for it)
```

A negative resistor delivers v²/|R| = 2²/2000 = 1 mW. The ledger books it
as **model-supplied**, per element, never as inflow from the exterior: an
ideal negative resistance stands for a tunnel diode, a Gunn diode or an
op-amp negative-impedance converter whose bias supply the model leaves
out (doc/M5.md). A fitness that measures harvested energy must count this
column as an input; `evo.ModelSuppliedLimit` is the objective that bounds
it.

`best.cir` opens in ngspice, LTspice and QSpice as it is (ngspice accepts
negative resistors too).

## Files

| File | What it holds |
|---|---|
| `circuit.cir` | the template netlist, both legs 1 kΩ (gain ½) |
| `config.json` | the gene, the objective and the GA settings |
| `export/` | the template's netlist, KiCad schematic and AC results sheet (MNA against Kron's torn solve), written by `KRONSPICE_EXPORT=1 go test ./internal/export` |

## See also

- [doc/M5.md](../../doc/M5.md), "Negative resistance and nonlinear elements".
- [rc_lowpass_values](../rc_lowpass_values) and
  [rectifier_filter](../rectifier_filter) — the other config-driven M5 examples.
