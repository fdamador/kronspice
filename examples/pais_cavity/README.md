# The Pais annular cavity (US 10,144,532 B2) as a Kron ring line

S. C. Pais's patent "Craft using an inertial mass reduction device" (U.S.
Navy, filed 2016, issued 2018) claims an annular resonant cavity between
two hull walls, filled with xenon, driven by microwave emitters, with the
outer wall charged, vibrating and rotated, to "create a local polarized
vacuum" and reduce the craft's inertial mass. This example puts that
cavity through a Kron circuit analysis and the energy ledger, beside the
same ring circuit at other scales: [ring_towers](../ring_towers) and
[double_hull](../double_hull). Register entry H18 (CLAUDE.md §1); the full
comparison is in `doc/pais-us10144532-comparison.md`.

```
go run ./examples/pais_cavity
go test ./examples/pais_cavity
```

**Standard prediction beside the claim:** a closed cavity's fields exert
no net time-averaged force on the craft (momentum conservation), mass is
unchanged, and the only way to push is to emit momentum, at most P/c. The
gyroscope weight reduction the patent cites (Hayasaka & Takeuchi 1989) was
not reproduced by three groups.

## The construction

The annulus is a closed parallel-plate TEM line around the craft (Paul
§4.2): L′ = μ₀d/h, C′ = ε₀h/d, R′ = 2R_s/h with R_s the copper surface
resistance at the drive frequency. The patent gives no gap, height,
voltage or drive power, so the geometry is declared: R = 2 m (its example
radius), gap 5 cm, walls 0.5 m, 640 sections (40 per wavelength), ring
mode n = 16 at 381.315 MHz (inside its 300 MHz–300 GHz band), Z₀ = 37.7 Ω,
1 MW total drive. Emitters are current injections.

**Kron analysis.** Happ's two-level diakoptics with a single common bus
(`kron.TwoLevel`; Happ 1971 Ch. XII §4, Fig. 12.12): the ring is torn into
16 sectors, each hinged at the inner wall, joined by 16 link branches, so
the interconnection system has order 16 against 1,280 MNA unknowns. The
twin is the same netlist through kronspice's MNA AC analysis.

## Results

| Case (1 MW scaling) | Q | Stored | V_max | E_max | Energy in the emitter's quadrant |
|---|---|---|---|---|---|
| A. 4 emitters at 0/90/180/270°, in phase | 1.477×10⁴ | 6.167 J | 148.9 kV | 2.98 MV/m | 0.250 |
| B. top emitter only, same current, on resonance | 1.477×10⁴ | 0.385 J (62.5 kW in) | 37.2 kV | 0.74 MV/m | 0.250 |
| C. top emitter only, between modes n and n + 1 | 1.5×10⁴ | 1.1 µJ (0.18 W in) | 63 V | — | 0.250 |
| D. top only, walls 2000× copper's R_s | 7.4 | 3.3×10⁻⁷ J (107 W in) | 64 V | — | 0.817 |
| E. niobium walls, R_s = 20 nΩ (declared) | 3.76×10⁹ | 6.17 J for 3.93 W | — | — | — |

- **Twin and ledger.** Off resonance Kron and MNA agree to 3.3×10⁻¹⁰. On
  the Q ≈ 1.5×10⁴ resonance the system is ill-conditioned by about Q·n,
  so both are checked by backward error (8.4×10⁻¹⁴ and 2.4×10⁻¹⁴) and
  differ by 5.8×10⁻⁸. Input power equals wall loss to 2.5×10⁻¹⁰; Q equals
  ωL′/R′ to 2×10⁻³; the resonance follows the lumped-ladder dispersion
  sin x/x exactly.
- **Sectional steering does not work on a high-Q ring.** One wave decays
  over 1/α = 2Z₀/R′ ≈ 3.7 km and the ring is 12.6 m around, so energy fed
  into one section fills the ring evenly on or off resonance. Only a ring
  with single-digit Q keeps the energy near its emitter (case D), as the
  ring-tower runs found: synchronized drives move the antinodes, not
  which sector holds the energy.

**The patent's figures against physical ceilings:**

| Figure | Value | Limit |
|---|---|---|
| σ = 1 C/m² | E = σ/ε₀ = 1.13×10¹¹ V/m | Kilpatrick 1.91×10⁷ V/m at 381 MHz; Fowler–Nordheim emission 2.45×10¹⁵ A/m² drains the charge in 4.1×10⁻¹⁶ s; wall pressure 5.65×10¹⁰ Pa against copper's 2.1×10⁸ Pa; RF drive to reach it here 1.44×10¹⁵ W |
| σ = 50,000 C/m² | 5.65×10¹⁵ V/m | 4×10⁻³ of the Schwinger field 1.32×10¹⁸ V/m; pressure 1.41×10²⁰ Pa |
| 10²⁴ W/m² (rotation) | reproduced: σ²/ε₀·Rω = 1.77×10²⁴ W/m² | the magnitude of E×H beside a moving charged surface; it circulates (∇·S = 0), so no energy leaves |
| 10³³ W/m² (vibration) | needs R_v·ν = 8.85×10²¹ m/s | 3×10¹³ c |
| 30,000 rpm at R = 2 m | rim 6,283 m/s, (v/c)² = 4.4×10⁻¹⁰ | burst margin 3.34 against a T1000G carbon-fibre hoop (`spacetime.AnalyzeRotor`): flies apart |
| wall vibration 10⁹–10¹⁸ Hz | — | lattice vibrations end at the Debye frequency: Cu 7.1×10¹² Hz, Al 8.9×10¹² Hz |
| thrust | — | zero net time-averaged force on a closed cavity; beaming all 1 MW out gives P/c = 3.34 mN |

**Bench test (H18).** Weigh a sealed, microwave-driven, charged cavity on
a balance against a dummy load dissipating the same power; the standard
prediction is no weight change beyond P/c of any escaping radiation, with
every watt booked in and out as heat. Status: *untested*.

## Files

| File | What it holds |
|---|---|
| `main.go` | the ring line, Happ two-level solve, MNA twin netlist, backward-error residual, the limits |
| `main_test.go` | Kron vs MNA on and off resonance; resonance, Q, ledger and localization; the patent figures |
| `export_test.go` | netlist, KiCad schematic and results sheet for case A in `export/` |
| `run_output.txt` | the run above |

## Sources

- S. C. Pais, "Craft using an inertial mass reduction device," U.S. Patent
  10,144,532 B2 (2018).
- H. H. Happ, *Diakoptics and Networks*, Academic Press, 1971, Ch. XII §4
  (`reference/`).
- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed.,
  §4.2 — parallel-plate line parameters.
- R. G. Forbes, *J. Vac. Sci. Technol. B* 26, 788 (2008), eq. 1 —
  Fowler–Nordheim emission; W. D. Kilpatrick, *Rev. Sci. Instrum.* 28,
  824 (1957), with T. J. Boyd's 1982 fit; C. Kittel, *Introduction to Solid
  State Physics*, 8th ed., Table 5.1 — Debye temperatures.
- H. Hayasaka, S. Takeuchi, *Phys. Rev. Lett.* 63, 2701 (1989); not
  reproduced by J. E. Faller et al., *PRL* 64, 825 (1990), T. J. Quinn,
  A. Picard, *Nature* 343, 732 (1990), and J. M. Nitschke, P. A. Wilmarth,
  *PRL* 64, 2115 (1990).
