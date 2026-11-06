# Acoustic drive of a sonoluminescence flask

Electrical → mechanical → acoustic, one circuit: a 50 Ω generator, an
optional series inductor that tunes out the transducer's capacitance, a
piezoelectric ring (clamped capacitance, force factor, moving mass, loss),
its coupling into the breathing mode of a 3 cm liquid sphere, and the mode
itself, in H₂O and D₂O. The output is the center pressure per volt, the
drive and power for a 1.3 atm amplitude, and where every watt goes. It is
the first circuit that joins three of the energy domains of CLAUDE.md §3.1
(electrical, mechanical, acoustic) through Kron's network analogies
(§3a.1 items 3 and 7), built after M9a.1 and described, with the solver
changes it needed, in [doc/domains.md](../../doc/domains.md). The liquids
come from [fluid/](../../fluid/) (§3h).

```
go run ./examples/acoustic_drive
go test ./examples/acoustic_drive ./fluid
```

## The construction

Every node belongs to one domain, and ground is the common reference of all
of them (zero velocity, ambient pressure). The chain is one netlist with
two non-electrical domains declared on it:

| Stage | Elements | Domain |
|---|---|---|
| Generator | `VGEN` 1 V AC, `RGEN` 50 Ω | electrical |
| Matching | `LMATCH` = 1/(ω0²·C0) with 0.5 Ω, or `RWIRE` 0.5 Ω alone | electrical |
| Piezo ring, electrical side | `C0` 4 nF clamped capacitance | electrical |
| Transducer | gyrator `YPZT`, force factor φ = 0.6 N/V | electrical ↔ mechanical |
| Piezo ring, mechanical side | `CMASS` 4 g moving mass, `LSPRING` 1/k with k = ω0²·m (resonance placed at the mode), `RLOSS` 1/b, b = 3 N·s/m | mechanical |
| Coupling | gyrator `YCOUPLE`, effective area S_eff = 2 cm² | mechanical ↔ acoustic |
| Flask breathing mode | `CMODE` Ca, `LMODE` Ma, `RMODE` Rpar from `fluid.SphericalResonator` | acoustic |

The mechanical domain uses the mobility analogy (velocity across, force
through), so a mass is a capacitor to ground and Kirchhoff's current law is
the force balance at a point; the acoustic domain is pressure across, volume
velocity through. A transducer is an ideal gyrator, `i1 = g·V2`,
`i2 = −g·V1`, which passes power from one port to the other without loss:
a piezo element of force factor φ has `g = φ`, a piston of area S has
`g = S`. In Kron's primitive network it is two axes coupled by an
antisymmetric mutual admittance, a non-reciprocal element the mesh and
node-pair forms carry like any other.

The breathing mode of a liquid sphere of radius R with a pressure-release
wall sits at f0 = c/(2R) and is modeled as a parallel acoustic resonator
whose across variable is the center pressure: Ca = 2R³/(πρc²) from the
mode's potential energy, Ma = 1/(ω0²·Ca), Rpar = Q/(ω0·Ca). The program
sweeps 2001 points over ±3 % of f0 with `analysis.AC`, finds the pressure
peak, scales the drive to 1.3 atm at the center and reads the ledger.

## Results

For a 1.3 atm center amplitude, the range where single-bubble
sonoluminescence is observed (about 1.2–1.5 atm at 20–40 kHz, Brenner,
Hilgenfeldt & Lohse 2002):

| Run | f0 (peak) | Pressure per volt | Generator EMF | Source power | Into the liquid | Transducer loss | Generator's own 50 Ω |
|---|---|---|---|---|---|---|---|
| H₂O, no matching | 25.03 kHz (24.27) | 3.22 kPa/V | 40.9 V | 10.9 mW | 6.96 mW | 3.50 mW | 0.42 mW |
| H₂O, matched | 25.03 kHz (25.02) | 69.4 kPa/V | 1.90 V | 26.2 mW | 6.96 mW | 0.42 µW | 19.0 mW |
| D₂O, no matching | 23.40 kHz (22.69) | 3.20 kPa/V | 41.2 V | 10.5 mW | 6.72 mW | 3.26 mW | 0.51 mW |
| D₂O, matched | 23.40 kHz (23.40) | 72.2 kPa/V | 1.82 V | 23.5 mW | 6.72 mW | 0.39 µW | 16.6 mW |

The ledger closes in every run: source power in equals heat out, residual 0
of the scale. The power the liquid takes at a given pressure is fixed by the
mode, P²/(2·Rpar); what matching changes is the drive voltage (×21 more
pressure per volt) and where the rest of the power goes. Without the
inductor the transducer's own loss takes a third of the input; with it,
most of the input is burned in the generator's 50 Ω, because the matched
chain is a low impedance at resonance. D₂O resonates 6.5 % lower than H₂O
(sound speed 1404 against 1502 m/s) and takes slightly less power.

`main_test.go` checks, for both liquids: liquid heat equals P²/(2·Rpar) to
1e-9; ledger residual below 1e-9; matching puts the peak back at f0 within
0.1 % and raises the pressure per volt by more than ten times; D₂O's f0 is
below H₂O's. The export sheets show the chain torn automatically at
`LMATCH` into two pieces, the torn and MNA solutions agreeing to 5e-12.

**What is assumed.** The flask Q (1500), the coupling area S_eff and the
transducer values (C0, φ, m, b) are declared assumptions in `main.go`, not
measurements. A hydrophone calibration of the real flask sets them; after
that the chain is a design tool, and an evolvable one (matching network,
transducer and drive frequency for the least power at a target pressure).
The bubble itself (Rayleigh–Plesset or Keller–Miksis dynamics, collapse,
light) was milestone M9p and stays out of the engine (CLAUDE.md §1,
non-goals); it enters as a cited equivalent circuit in
[bubble_collapse](../bubble_collapse), which puts it on this chain's node p.
[doc/future-work.md](../../doc/future-work.md) notes that these acoustic
modes, about 25 kHz, sit ten orders of magnitude below molecular
vibration frequencies.

## Files

| File | What it holds |
|---|---|
| `main.go` | the chain builder (`chain`), the four runs and the ledger report |
| `main_test.go` | `TestChain`: P²/2R heat, ledger, matching, D₂O below H₂O |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` (`KRONSPICE_EXPORT=1`) |
| `export/` | `acoustic_drive_water`, `_water_matched`, `_d2o_matched`: `.cir`, `.kicad_sch`, `_results.md/.csv` |

## Sources

- M. P. Brenner, S. Hilgenfeldt, D. Lohse, "Single-bubble sonoluminescence,"
  *Rev. Mod. Phys.* 74, 425 (2002), §II.A — the driving-pressure range and
  the flask resonator.
- L. E. Kinsler, A. R. Frey, A. B. Coppens, J. V. Sanders, *Fundamentals of
  Acoustics*, 4th ed., Wiley (2000), §7.5 (piston radiation impedance,
  eqs. 7.5.11–7.5.13) and §9.3 (spherical resonator), cited in
  `fluid/acoustic.go`.
- F. A. Firestone, "A new analogy between mechanical and electrical
  systems," *J. Acoust. Soc. Am.* 4, 249 (1933) — the mobility analogy.
- W. Wagner, A. Pruß, "The IAPWS formulation 1995 for the thermodynamic
  properties of ordinary water substance for general and scientific use,"
  *J. Phys. Chem. Ref. Data* 31, 387 (2002), Table 7 — H₂O at 300 K.
- S. Herrig, M. Thol, A. H. Harvey, E. W. Lemmon, "A reference equation of
  state for heavy water," *J. Phys. Chem. Ref. Data* 47, 043102 (2018),
  Table 7 — D₂O at 300 K.
- G. Kron, *Tensor Analysis of Networks* (1939), Ch. V p. 133 — separate
  networks of any kind (electrical, mechanical, acoustical) joined by a
  connection tensor (CLAUDE.md §3a.1).
