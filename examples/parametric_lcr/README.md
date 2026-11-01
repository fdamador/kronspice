# Parametric LCR circuit with a time-varying capacitor (H16)

F. D. Amador and H. Eckardt, "Simulation of a Parametric Resonance
Circuit" (A.I.A.S./UPITEC, 30 March 2012) drive a series R, L, C(t) loop
and report a coefficient of performance of about 200 against the source's
input, read as energy "from the time domain" (spin connection resonance).
The example reproduces every case in the paper with kronspice's
time-varying capacitor and puts the one agent the paper's energy sum leaves
out on the ledger: the pump that varies C (hypothesis H16, CLAUDE.md §1;
Kouřil & Vrba §10.2 in `doc/kouril-vrba-verification.md`). The same
parametric pumping in Kron's machine frame is in
[kron_machine](../kron_machine) and [parametric_pump](../parametric_pump).

```
go run ./examples/parametric_lcr
go test ./examples/parametric_lcr
```

## The circuit

Paper Fig. 2, eqs. 19–27: a sine source U = U0·sin(ωt), inductance L,
resistance R and a capacitance C_eff(t) = C·(1 + a0·sin(ω1·t + φ)) with
q = C_eff·u_C (eq. 25). Parameters (eqs. 28–36): U0 = 5 V, C = 230 nF,
L = 20 mH, R = 20 Ω, so f0 = 2346.6 Hz, Q = 14.74, α = R/2L = 500 1/s;
ω1 = 2ω0, a0 = 0.2, φ = 0.

kronspice's `VaryingCapacitor`, q = s(t)·C·v, is eq. (25) exactly. The
work the agent that moves the plates (or a varactor's bias supply) does on
the circuit is P_pump = −ṡ·½·C·v², the pump of Mandelstam and Papalexi
(1932, the paper's ref. [13]) and of Kouřil & Vrba §10.2.2. Every case is
checked as source + pump = heat in R + change of stored energy.

Standard theory beside the paper: a0 = 0.2 gives an elastance depth
m = 0.2020 against the tuned self-excitation threshold m_c = 2/Q = 0.1356
(K&V 10.101), so pumping at 2ω0 self-excites. C(1 + a sin) raises the
mean elastance by 1/√(1 − a²), so the circuit's own frequency is
1.0103 f0; pumped at 2 × 0.97 f0 the detuning is 4.2 % and the threshold
rises to 0.2144 > m, below self-excitation: a parametric amplifier of the
drive (K&V 10.103).

## Results

| Paper case | kronspice | Books (source + pump vs R + stored) |
|---|---|---|
| Fig. 4, constant C at f0 | peak 0.25 A = U0/R; source 0.0175 J → R 0.0169 J + stored 0.0006 J | 0.001 % |
| §3.1 Fig. 5, C(1 + 0.2 sin 2ω0t) at f0, 30 ms | 87.3 A and growing without limit, σ = 245 1/s as K&V predict; source −0.085 J, pump 256 J → R 165 J + stored 91 J | 0.07 % |
| §3.1 Figs. 6–7, drive and pump at 0.97 f0, φ = 0 | bounded 7.34 A, 504 W; pump 241.5 J, source −3.8 J | 0.04 % |
| same, φ = −90° | bounded 16.1 A, 2.43 kW; source 6.75 J, pump 1139 J → R 1143 J | 0.04 % |
| same, φ = +90° | bounded 5.72 A, 307 W; pump 138 J | 0.04 % |
| §3.2 Fig. 8, C(1 + 0.2 sin⁷ ω1t), ω1 = drive = 1.07 ω0 | 0.123 A, no growth | 0.4 % |
| §3.2 Fig. 8, a0 = 0.8 | 207 A; pump 607 J → R 262 J + stored 345 J | 0.3 % |
| §3.2 Fig. 9, rectangular C 0.8C ↔ 1.2C at 0.95 f0, drive off after 2 periods | decays, 1.1 mA | 0.4 % |
| §3.2 Fig. 9, rectangular C 0.2C ↔ 1.8C | 7.6 kA from a 2-period drive, or 443 A from a 1 mA start with no drive; pump 2360 J → R 614 J + stored 1854 J | few % (integration error, falls with the step: 0.97 % → 0.32 % from 2000 to 8000 steps per period) |

The paper's COP divides the heat in R by the source energy alone: 0.964
for Fig. 4, 169 for the φ = −90° case the paper highlights, and negative
where the source absorbs. Divided by source plus pump it is 0.998 in the
bounded cases and below 1 everywhere. Pumping at the drive frequency
(§3.2) grows only for a0 ≳ 0.7, as the second Mathieu tongue requires.

In every case the heat in R and the stored energy are paid by the source
plus the pump that varies C. The pump's share is what the paper's "energy
from spacetime" counts. **Bench test:** measure the power that drives
C(t), the motor shaft or the varactor bias supply. Status in the register:
*untested* (reproduced; the COP omits the pump).

## Files

| File | What it holds |
|---|---|
| `main.go` | the cases, `Run` with the pump's power on the ledger, K&V thresholds (`ElastanceIndex`, `Threshold`, `MeanShift`) |
| `main_test.go` | classical amplitude, tuned growth paid by the pump, detuned bounded case, second tongue (a0 = 0.2 vs 0.8) |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |

## Sources

- F. D. Amador, H. Eckardt, "Simulation of a Parametric Resonance Circuit," A.I.A.S./UPITEC, 30 March 2012 (`reference/amador_eckardt_2012_parametric_lcr.pdf`). Hypothesis source.
- L. Mandelstam, N. Papalexi, *Z. Phys.* 73, 223 (1932).
- F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits*, Ellis Horwood 1988, §10.2 (`reference/kouril-vrba/`, `doc/kouril-vrba-verification.md`).
- J. M. Manley, H. E. Rowe, *Proc. IRE* 44, 904 (1956).
