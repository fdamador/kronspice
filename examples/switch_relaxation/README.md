# switch_relaxation

Two relaxation oscillators from one parameterized subcircuit. The netlist
([relax.cir](relax.cir)) uses the netlist preprocessing and the switch
device together:

- `.param vcc=5 rc=10k cc=1u` at the top level;
- `.subckt relax top out params: r={rc} c={cc} vt=2.5 vh=0.5`, with defaults
  that refer to the top-level parameters, and a local
  `.model swh SW(VT={vt} VH={vh} …)` (renamed `xa.swh`, `xb.swh`);
- two instances, `XA` with the defaults and `XB` with `r={2*rc} vh=1`;
- an `S` switch across each capacitor, controlled by its own voltage.

Each cell charges C through R until the switch turns on at VT + VH, then
discharges through RON (1 kΩ) until it turns off at VT − VH. The exact
period is

    T = R·C·ln((Vcc − lo)/(Vcc − hi)) + τ·ln((hi − v∞)/(lo − v∞)),
    v∞ = Vcc·RON/(R + RON),  τ = C·R·RON/(R + RON).

| cell | R | window | exact | simulated |
|---|---|---|---|---|
| XA | 10 kΩ | 2.0–3.0 V | 4.50828 ms | 4.50869 ms (+9e-5, 15 cycles) |
| XB | 20 kΩ | 1.5–3.5 V | 17.8504 ms | 17.8551 ms (+2.6e-4, 3 cycles) |

The switching instants are located by the switch itself
(`device.Switch.LimitStep`: a step in which the control crosses a threshold
is repeated to end just after the crossing), so the period does not drift
by one step per cycle. The remaining difference is the resolution of the
turn-on instants read from the recorded points. The ledger's integration
mismatch is 3.8e-4 of the energy supplied by V1; all of it comes from V1.

Run:

    go run ./examples/switch_relaxation
    go test ./examples/switch_relaxation
