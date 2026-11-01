# US 10,144,532 B2 (Pais) against the ring towers — Kron analysis

S. C. Pais, "Craft using an inertial mass reduction device", U.S. Patent
10,144,532 B2, assignee the U.S. Navy, filed 28 Apr 2016, issued 4 Dec 2018
(Justia: patents.justia.com/patent/10144532). Code:
`examples/pais_cavity` (run: `run_output.txt`). Register entry: H18.

## What the patent claims

1. An inner and an outer resonant-cavity wall form an annular cavity.
   Microwave emitters drive it, which makes the outer wall vibrate "in an
   accelerated mode" and "create a local polarized vacuum" outside it.
2. The cavity is filled with a noble gas (xenon).
3. The outer wall is electrically charged.
4. The cavity is rotated "in an accelerated mode".

Figures from the description:

- emitters at 300 MHz–300 GHz;
- cavity Q of 10⁴–10⁹;
- wall vibration at 10⁹–10¹⁸ Hz;
- example: σ = 50,000 C/m² on a 2 m disc at 30,000 rpm;
- energy flux S = f_G·(σ²/ε0)·[R_r·ω + R_v·ν + v·R], giving about 10²⁴ W/m² for rotation and 10³³ W/m² at σ = 1 C/m² and 10¹¹ Hz;
- steering by driving sections of the shroud separately.

The patent gives no wall gap, height, voltage or drive power.

## Same circuit as the ring towers and the double hull

| | `ring_towers` | `double_hull` | Pais cavity (`pais_cavity`) |
|---|---|---|---|
| Ring | Earth–ionosphere (normalized LC ladder) | gap between a craft's two hulls | gap between inner and outer cavity walls |
| Line | 96 sections, Z0 31.6 Ω | parallel-plate TEM | parallel-plate TEM, 640 sections, Z0 37.7 Ω |
| Mode | 4 wavelengths, 1.31 MHz | chosen ring mode | 16 wavelengths, 381 MHz (declared R 2 m, gap 5 cm, height 0.5 m) |
| Loss per lap | 0.2 Np (Q ≈ 63) | copper walls | 3.4×10⁻³ Np, copper walls (Q 1.48×10⁴) |
| Drive | launcher pairs at 0/90/180/270°, phase-synchronized | radiators in the gap, plus DC bias | emitters at 0/90/180/270° (patent: "emitters arranged in the cavity", sections driven separately) |
| Charge | DC bias per pair | DC bias between hulls | outer wall "electrically charged" |
| Claim | none (energy thought experiment) | none (voltage, multipactor) | inertial mass reduction |

So the Pais hull is the `double_hull` circuit: the closed ring
transmission line of `ring_towers`, at craft scale and microwave frequency.

## Kron analysis

The ring is torn into 16 sectors. Each sector is a subdivision hinged at
the common bus, which is the inner wall. The sectors are joined by 16 link
branches. This uses Happ's two-level single-common-bus diakoptics
(`kron.TwoLevel`; Happ 1971 Ch. XII §4, Fig. 12.12). The interconnection
system has order 16, against 1,280 unknowns for MNA.

**Twin check against MNA AC.**

- Off resonance, Kron and MNA agree to 3.3×10⁻¹⁰.
- On resonance, the system is ill-conditioned (by about Q·n). Both solutions have backward error ≤ 10⁻¹³, and they differ by 5.8×10⁻⁸.

**Energy ledger.** Input power equals wall loss to 3×10⁻¹⁰.

Results at 1 MW total drive:

| Case | Q | Stored | V_max | E_max | Energy in emitter's quadrant |
|---|---|---|---|---|---|
| A. 4 emitters, in phase (ring-tower sync) | 1.48×10⁴ | 6.17 J | 149 kV | 2.98 MV/m | 0.25 |
| B. top emitter only, same current, on resonance | 1.48×10⁴ | 0.385 J (62.5 kW in) | 37 kV | 0.74 MV/m | 0.25 |
| C. top emitter only, between modes | — | 1.1 µJ (0.18 W in) | 63 V | — | 0.25 |
| D. top only, walls 2,000× lossier | 7.4 | — | — | — | 0.82 |
| E. niobium walls, Rs 20 nΩ | 3.8×10⁹ | 6.17 J for 3.9 W | — | — | — |

**Sectional steering does not work on a high-Q ring.** One wave decays
over 1/α = 2·Z0/R′ ≈ 3.7 km, and the ring is only 12.6 m around. So
energy fed into one section fills the whole ring evenly, whether the drive
is on or off resonance. Only a ring with Q in the single digits keeps the
energy near its emitter, as in case D. This matches the ring-tower result:
synchronized pairs change *where the antinodes are*, not *which sector
holds the energy*.

## The patent's figures against physical limits

| Figure | Value | Limit |
|---|---|---|
| σ = 1 C/m² | E = σ/ε0 = 1.13×10¹¹ V/m | Kilpatrick 1.9×10⁷ V/m at 381 MHz; Fowler–Nordheim emission 2.5×10¹⁵ A/m² drains the charge in 4×10⁻¹⁶ s; wall pressure σ²/2ε0 = 5.6×10¹⁰ Pa vs copper's 2.1×10⁸ Pa; RF drive to reach it in this cavity 1.4×10¹⁵ W |
| σ = 50,000 C/m² | 5.6×10¹⁵ V/m | 4×10⁻³ of the Schwinger field (1.32×10¹⁸ V/m); pressure 1.4×10²⁰ Pa |
| 10²⁴ W/m² (rotation) | reproduced: σ²/ε0·Rω = 1.77×10²⁴ | This is the magnitude of E×H next to a moving charged surface. That flux circulates (∇·S = 0), so no energy leaves. |
| 10³³ W/m² (vibration) | needs R_v·ν = 8.9×10²¹ m/s | 3×10¹³ c |
| 30,000 rpm, R = 2 m | rim 6,283 m/s | 3.3× the burst speed of a T1000G carbon-fibre hoop (`spacetime.AnalyzeRotor`); (v/c)² = 4.4×10⁻¹⁰ |
| wall vibration 10⁹–10¹⁸ Hz | — | No lattice vibration exists above the Debye frequency: Cu 7.1×10¹² Hz, Al 8.9×10¹² Hz |
| thrust | — | A closed cavity exerts no net time-averaged force on itself. Beaming all 1 MW out gives P/c = 3.3 mN. |

The patent's experimental support is the 1989 gyroscope weight reduction
(Hayasaka & Takeuchi, PRL 63, 2701). Three groups did not reproduce it:

- Faller et al., PRL 64, 825 (1990);
- Quinn & Picard, Nature 343, 732 (1990);
- Nitschke & Wilmarth, PRL 64, 2115 (1990).

## What a bench test would measure (H18)

The standard prediction is computed beside the patent's, as the register
requires:

- weigh a sealed, microwave-driven, charged cavity on a balance, against a dummy load dissipating the same power (to separate heating and convection);
- the standard prediction is no weight change beyond P/c of any escaping radiation, and the null results above bound the gyroscope version;
- the ledger books every watt into the emitters and out as heat.
