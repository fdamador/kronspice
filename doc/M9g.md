# M9g — the open-system energy ledger and the lumped limit

CLAUDE.md §3c ("Open systems and the energy ledger", "What em/ does"
items 3–4). The ledger itself (inflow, outflow, dissipation, storage,
model-supplied power, per element × channel) was built with M1–M2 and is
checked inside every analysis (`em.CheckLedger`, the transient
recorder, the AC balance). M9g back-fills the checks over all golden
circuits and adds what was missing.

## New

- **`em/limits.go`** — `CheckQuasiStatic(c, lengths, overall, fmax, opt)`
  compares declared conductor lengths (element → metres, and the overall
  size) with λ = c/(f·√εr) at the highest frequency; above `Threshold`
  (default λ/10) it returns `ErrNotQuasiStatic` naming each part and its
  d/λ, or only reports with `WarnOnly`. `SourceBandwidth(c, defaultEdge)`
  gives the highest frequency of a transient's sources: SIN at its
  frequency, PULSE and PWL edges and EXP rises at the knee 0.5/t_r
  (Johnson & Graham, High-Speed Digital Design, §1.1). em does not import
  analysis, so the caller passes fmax (an AC stop frequency, or
  `SourceBandwidth`). `evo.Penalty` maps `ErrNotQuasiStatic` to "invalid",
  so a `Problem.Build` that checks the limit penalizes electrically large
  candidates.
- **`em/efficiency.go`** — `OPEfficiency(result, loads…)` and
  `TranEfficiency(waveforms, from, to, loads…)`: energy absorbed by the
  named elements over the net energy the exterior delivers
  (∫(Inflow − Outflow): energy a reactive load hands back to its source
  within a cycle is not counted as delivered; model-supplied power is
  never inflow). `Tellegen(c, nodeV, elemI)` is Σ v_k·i_k for circuits of
  two-terminal elements, for the theorem's general form.
- **`evo.OPEfficiency`, `evo.TranEfficiency`** — 1 − η as objectives.

## Checks

| Check | Result |
|---|---|
| Ledger on every golden circuit, every analysis in its netlist plus the operating point | 40 OP, 12 DC sweeps, 6 AC sweeps, 13 transients: no residual above the round-off estimate (the largest absolute residual is 1e−15 W, at points whose powers are themselves ~1e−11 W) |
| Tellegen at the operating point of every two-terminal golden circuit | 20 circuits, Σ v·i = 0 to 1e−12 of Σ\|v·i\| |
| Tellegen's general form: voltages of one circuit × currents of another with the same graph; voltages and currents of one transient at two different instants | 50 random pairs, 1e−12; 1e−9 |
| DC efficiency R_L/(R_L + 50) from the ledger | exact to 1e−12 |
| Transient efficiency over whole periods = AC phasor efficiency of the same circuit (R_S, R_L ∥ C) | 1e−4 |
| λ/10 limit: 1 m at 0.99·c/10 passes, at 1.01·c/10 fails naming the part; εr = 4 halves λ; WarnOnly | as stated |
| Source bandwidth for SIN, PULSE (given and defaulted edges), EXP, PWL | exact |
| A 100-turn, 5 cm coil (15.7 m of wire) is lumped only to 1.9 MHz | as stated |
| GA: efficiency objective drives R_L to its 1 kΩ bound (η = 0.9524); a lead longer than λ/10 at 100 MHz is penalized "invalid" and never wins | as stated |

## Not done here

- Automatic sizes: elements carry no physical length yet, so lengths are
  declared by the caller (or come from em.Coil and, later, em/extract).
- The rule of at least 20 points per rise time for arc and pulse runs
  (§3g) belongs with M9o.
