# A craft's double hull as a ring waveguide

The ring-tower experiment of [ring_towers](../ring_towers) moved inside a
craft's double wall: the gap between an inner and an outer conducting hull
is a closed parallel-plate waveguide around the craft, driven at one of its
ring resonances by radiators in the gap, with a DC bias between the hulls
added to the RF. The questions are what gap voltage a given RF power
reaches, whether the vacuum gap multipactors on the way up (CLAUDE.md §3i),
and what the least RF power is for a target peak voltage. The same
multipactor model drives the electron pumps of
[photon_electron](../photon_electron).

```
go run ./examples/double_hull
go test ./examples/double_hull
```

## The construction

The gap is the parallel-plate TEM line of the telegrapher's equations
(Paul §4.2): per unit length around the hull L′ = μ0·d/w, C′ = ε0·w/d and
Z0 = η0·d/w, with d the gap and w the hull height; wall loss
R′ = 2·Rs/w with Rs = √(π·f·μ0·ρ) the skin-effect surface resistance of
each copper wall, so α = Rs/(η0·d). N sections (24 per wavelength) close
the ring. The AC analysis gives the steady state at the drive frequency,
the operating point gives the bias, and the gap voltage is their sum.
Expected traveling-wave amplitude at mode n: √(Z0·P_in·Q/(π·n)); a
standing wave doubles it at the antinodes.

The DC bias charges the outer hull through a high resistance, stores
½·C_gap·V_dc² once and draws no steady power. The gap then swings between
V_dc ± V_rf. A DC bias across a coaxial gap is the measured cure for
multipacting in the LEP2 352 MHz couplers (2.5 kV) and the KEKB 509 MHz
couplers (2 kV).

**Multipactor check** (`em/multipactor.go`): two-surface resonance between
parallel walls (Hatch & Williams 1954; Vaughan 1988). Electrons launched
from either wall at every RF phase with 2 eV are tracked through the
combined RF + DC field impact after impact, multiplying Furman's yield
δ = δmax·s·x/(s − 1 + x^s) at each one, with the measured copper values of
Baglin et al. 2002 (as received δmax 2.03 at 262 eV, s 1.39; conditioned
1.13 at 318 eV, s 1.35). A phase whose mean yield per impact stays above 1
is a growing discharge. Space charge and the spread of emission energies
are left out, so the check marks where a discharge can start, not how
large it grows. Ring-up passes through every amplitude below the final
one, so each is checked. Kilpatrick's RF sparking limit
(f[MHz] = 1.64·E²·exp(−8.5/E), E in MV/m) bounds the peak gap field.

## Results

Example hull: radius 5 m, height 10 m, gap 0.10 m, copper, mode 4,
96 sections, 1 kW from the radiators.

| Run | f | Q | RF gap amplitude | Formula | Peak gap (RF + bias) |
|---|---|---|---|---|---|
| traveling, 0 V | 37.685 MHz | 9410 | 1681–1683 V | 1680 V | 1683 V |
| traveling, 2500 V | 37.685 MHz | 9410 | 1681–1683 V | 1680 V | 4183 V |
| standing, 0 V | 37.685 MHz | 9410 | 0–2379 V | 1680·√2 V | 2379 V |
| standing, 2500 V | 37.685 MHz | 9410 | 0–2379 V | 1680·√2 V | 4879 V |

Ledger in every run: radiators 1000.0 W = wall loss 999.9 W (999.7 W
standing) + bias feed; 39.7 mJ of RF energy in the gap (Q·P/ω); the 2.5 kV
bias stores 87.8 mJ once. The tests hold the resonance to n·c/ℓ times the
ladder's lattice factor to 1e-9, the amplitude to the formula within 1 %,
and the exported 193-node circuit torn into 6 pieces equals MNA to 5.6e-9
with the same power.

**Multipactor at f·d = 3.77 MHz·m**, worst-phase mean yield per impact
(as-received copper; > 1 grows):

| V_rf | bias 0 V | 500 V | 2500 V | 5000 V |
|---|---|---|---|---|
| 300 V | 1.10 at 53 eV, susceptible | 0.07 | 0.06 | 0.06 |
| 1000 V | 1.78 at 630 eV, susceptible | 0.11 | 0.07 | 0.06 |
| 1700 V | 0.29 | 0.16 | 0.07 | 0.06 |
| 10000 V | 0.06 | 0.06 | 0.07 | 0.07 |

Ring-up to 1680 V (traveling) or 2375 V (standing) without bias meets a
worst growth of 2.02 at about 760 V for as-received copper (1.13 for
conditioned); the least bias on the grid that clears the whole ring-up is
500 V. Kilpatrick's limit at 37.7 MHz is 8.1 MV/m against 0.049 MV/m for
the biased standing wave.

**Evolved** (gap, mode, RF power and bias for a 20 kV peak with at least
5 kV of it RF, least RF power, multipactor-free ring-up, under
Kilpatrick): best of 332 evaluations is a 0.470 m gap at mode 1
(9.09 MHz, Q 21729), 451 W of RF and 12496 V of bias, giving
7467 V + 12496 V = 19963 V peak, worst ring-up growth 0.06 and
0.042 MV/m against Kilpatrick's 5.3. Ledger: 451.0 W = wall loss
448.3 W + feed 2.720 W; the bias supply delivers C·V² = 0.97 J once, half
stored and half lost in the feed resistor; no steady DC current (no gap
leakage is modeled).

## Files

| File | What it holds |
|---|---|
| `main.go` | the hull as a ring ladder, the AC and OP solve, the ring-up and bias scans, Kilpatrick's limit, the GA problem |
| `main_test.go` | circuit against the waveguide formulas, Furman's yield, multipactor zones and bias suppression, Kilpatrick's inverse |
| `export_test.go` | traveling and standing netlists, KiCad schematics and results sheets in `export/` |
| `run_output.txt` | the run above |

## Sources

- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed., §4.2 (`reference/`).
- A. J. Hatch, H. B. Williams (1954); J. R. M. Vaughan, *IEEE Trans. Electron Devices* 35, 1172 (1988) — two-surface multipactor.
- V. Baglin, I. Collins, B. Henrist, N. Hilleret, G. Vorlaufer, CERN LHC Project Report 472 (2002) — copper secondary-electron yield, in Furman's fit.
- H. P. Kindermann, E. Haebel, M. Stirbet, V. Veshcherevich, EPAC 1996 (LEP2 couplers); Y. Kijima et al., EPAC 2000 (KEKB couplers) — DC bias against multipacting.
- W. D. Kilpatrick, *Rev. Sci. Instrum.* 28, 824 (1957); T. J. Boyd, LANL LA-UR-82-1166 (1982) — the sparking criterion.
