# Parametric pump: an LC tank whose capacitor is modulated at twice resonance

A one-netlist example of the time-varying capacitor, the kronspice
`SCALE=` extension (doc/STATUS.md, "Nonlinear and time-varying capacitors
and inductors"). A 1 mH, 1 µF tank (f0 = 5.033 kHz) starts with 1 V on the
capacitor, and the capacitance is driven as C(t) = 1 µF·(1 − 0.1·sin 2ω0t):
it falls while |v| peaks, so the tank's energy grows. It is the capacitive
twin of the pumped inductance in [kron_machine](../kron_machine) (M9c.1,
"temporal rate of change"), and the simplest form of the parametric
circuits of H16 ([parametric_lcr](../parametric_lcr)).

```bash
go run ./cmd/kronspice examples/parametric_pump/pump.cir
go run ./cmd/kronspice -o pump.csv examples/parametric_pump/pump.cir     # waveforms and ledger columns
go test ./analysis -run 'TestParametricPumping|TestParametricSelfExcitation'
```

## The netlist

```
C1 a 0 1u SCALE=SIN(1 0.1 10.066k 0 0 180) IC=1
L1 a b 1m
R1 b 0 0.1
.tran 1u 2m UIC
```

`SCALE=wave(...)` multiplies the element's charge curve by s(t): the
charge is q = s(t)·C·v, so the current is i = C·s·dv/dt + C·v·ds/dt. The
work of changing s, ṡ·∫Q, is the work of whatever agent squeezes the
capacitor, and it is booked on the element's **mechanical port** (an
`EXT=` channel can redirect it). The 180° phase makes s fall as |v| peaks;
with 0° the tank decays and the port takes energy out instead
(`analysis.TestParametricPumping`).

## Result

The envelope of v(a) grows from 1.02 V at the first peak to 4.44 V at
1.992 ms, an exponential rate of 740.8 /s. Mathieu's first-order rate for
a pump depth m at 2ω0, less the series damping, is
m·ω0/4 − R/(2L) = 790.6 − 50 = 740.6 /s (Kouřil & Vrba §10.2.2, eq. 10.101;
`analysis.TestParametricSelfExcitation` holds the engine to it).

Energy ledger over the 2 ms (from the run above):

| | |
|---|---|
| in through the mechanical port | 11.70 µJ |
| returned to the port | 1.33 µJ |
| dissipated in R1 | 0.61 µJ |
| stored, start → end | 0.50 µJ → 10.25 µJ |
| integration mismatch | 1.2 nJ |

The net work of the pump, 10.36 µJ, equals the stored-energy gain plus the
heat: nothing comes from nowhere, and the ledger's worst point residual is
0. This is the standard prediction beside which H16's "energy from the
time domain" is read (CLAUDE.md §1, H16): the agent that varies C pays.

## Files

| File | What it holds |
|---|---|
| `pump.cir` | the netlist |
| `export/` | the same netlist, a KiCad schematic and the operating-point sheet, written by `KRONSPICE_EXPORT=1 go test ./internal/export` |

## Sources

- F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits*, Ellis Horwood
  1988, §10.2.2 (parametric self-excitation), and `doc/kouril-vrba-verification.md`.
- L. Mandelstam, N. Papalexi, *Z. Phys.* 73, 223, 1932 — parametric excitation.
- [doc/M9c1.md](../../doc/M9c1.md) — the pumped inductance against Mathieu's rate.
