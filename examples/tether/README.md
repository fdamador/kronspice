# Electrodynamic tether in low Earth orbit

A conductor pushing and pulling against a planetary magnetic field
(`doc/em-push-pull.md`; CLAUDE.md §3f; milestone M9n). A vertical aluminium
wire moves east at orbital speed through the northward geomagnetic field;
the Lorentz force on its charges, by grade selection (`em.MotionalEMF`),
gives an EMF (v×B)·L along it. The circuit closes through the ionosphere
by plasma contactors at the two ends. The same tether coiled up is
[coiled_tether](../coiled_tether); the force's reaction is booked in the
momentum ledger as every push–pull design must (§3f). NASA flew this
circuit three times (PMG 1993, TSS-1 1992, TSS-1R 1996) and designed it
twice more (ProSEDS, the ISS reboost study); their numbers are set beside
the example's below.

```
go run ./examples/tether
go test ./examples/tether
```

Declared, not looked up: v = 7.7 km/s, B = 30 µT (order of magnitude in
LEO; a real orbit takes B from a cited field model such as IGRF), contactor
resistances 100 Ω each (a linear stand-in for the contactors' nonlinear
current–voltage characteristics), 5 km of 1 mm aluminium wire (169 Ω).

## The construction

Two modes of the same loop, solved as a kronspice operating point with the
ledger on:

```
EMF — tether resistance — contactors — load          (generator)
EMF — tether resistance — contactors — supply > EMF  (thruster)
```

The current feels ∫I·dl×B. As a generator the force is drag: the power
the load receives comes from the orbit's kinetic energy, F·v = −EMF·I.
With a supply driving the current against the EMF the force turns into
thrust. `em.MomentumLedger` books the force on the tether with its
reaction on the geomagnetic field and ionosphere, the exterior source.

This is the circuit NASA's *Tethers in Space Handbook* writes for the
electrodynamic tether (Cosmo & Lorenzini 1997, §4): the EMF
`V = ∫(v×B)·dl`, which for a straight tether perpendicular to the field is
`V = L·v·B·sin θ`; the force `F = ∫I·dl×B`; and the loop equation
`V_ind = I·R + ΔV_low + ΔV_up + ΔV_ion + ΔV_load`, where the two ΔV terms at
the ends are the drops across the space-charge sheaths or plasma clouds
of the contactors and ΔV_ion is the drop across the ionosphere, of order
1 to 20 Ω. The handbook notes that the sheath regions are the largest
impedances in the loop; the example's two 100 Ω resistors stand in for
them, and the ionospheric drop is folded into the same resistors. The
handbook also states the reversal the example shows: fed from an on-board
supply against the induced EMF the same tether is a thruster, at the
expense of that electrical power.

## Results

EMF from the Lorentz force 1155.0 V, equal to vBL.

| | Generator, 500 Ω load | Thruster, 3000 V supply |
|---|---|---|
| Current up the tether | 1.330 A | −5.004 A |
| Electrical power | 884 W to the load | 15,012 W from the supply |
| Force along x | −0.199 N (drag) | +0.751 N (thrust) |
| F·v | −1536 W: the orbit pays | +5780 W into the orbit |
| Thrust per power | — | 0.05 N/kW |
| Momentum ledger | reaction +0.199 N on the field and ionosphere | reaction −0.751 N on the field and ionosphere |

Checks in the test (all to 10⁻⁹): F·v = −EMF·I in both modes; generator
load power plus I²R losses equals EMF·I; thruster supply power equals
losses plus F·v; the momentum ledger's net force equals the force on the
tether; the energy ledger closes to 10⁻¹².

## Against NASA's flights and studies

| | This example | PMG (1993) | TSS-1R (1996) | ProSEDS (designed) | ISS reboost study (1998) |
|---|---|---|---|---|---|
| Tether | 5 km, 1 mm Al, 169 Ω | 500 m, #18 AWG Teflon-insulated Cu | 19.7 km deployed of 20.7 km; 10 strands of 34 AWG Cu on a Nomex core, 0.305 mm FEP insulation, 2.54 mm overall | 5 km of 7-strand 28 AWG 1350-O Al, mostly bare, plus 10 km Dyneema | 10 km, bare lower section, < 200 kg |
| Orbit | v = 7.7 km/s, B = 30 µT, v ⟂ B | 193 × 869 km, 25.7° | 297 km at the break, near the equator | 400 km circular, 36° | ISS orbit |
| EMF | 1155 V (231 V/km) | bias ±65 to −130 V applied in series with the I×B EMF | 3482 V (177 V/km) | — | about 100 V/km |
| Current | 1.33 A generator, 5.0 A thruster | 100–300 mA by day, 10–50 mA at night; 0.3 A at 130 V bias | 0.97 A just before the break; 1.1 A for 75 s after it | — | — |
| Force | 0.199 N drag, 0.751 N thrust | — | — | about 0.4 N drag; 14 km/day descent | 0.5 N from 5 kW, 0.8 N from 10 kW |
| Power | 884 W load; 15 kW supply | — | — | about 100 W to the spacecraft | 5–10 kW from the arrays |
| Ends | two 100 Ω linear contactors | hollow-cathode xenon plasma contactors at both ends, current fully reversible | 1.6 m conducting sphere as anode, electron guns and the Orbiter ground as cathode | bare wire collects, hollow cathode emits | bare wire collects, plasma contactor on the ISS emits |

What the comparison says about the example's declared inputs:

- **EMF per kilometre.** The example's 231 V/km is the upper end: it takes
  v perpendicular to B at 30 µT. TSS-1R measured 177 V/km near the
  equator at 297 km, and the ISS study quotes about 100 V/km for the
  station's 51.6° orbit. The handbook's `sin θ` is the difference, which
  is why the code says a real orbit must take B from a field model.
- **Thrust per watt.** The example's thruster gives 0.05 N/kW against the
  study's 0.1 N/kW at 5 kW and 0.08 N/kW at 10 kW. Of the 15 kW drawn,
  9.2 kW is I²R heat in the 369 Ω of tether and contactor resistance and
  5.8 kW goes into the orbit. NASA's design cuts both: a bare tether
  collects electrons along its length with no end sheath to cross
  (Sanmartín, Martínez-Sánchez & Ahedo 1993; the "breakthrough" of the
  study's §6.2), and a hollow cathode emits them at a drop of tens of
  volts.
- **The contactors.** The linear 100 Ω stand-ins are the weakest part of
  the model. TSS-1R's sphere collected currents up to three times what
  the Parker–Murphy model predicted (handbook §1.1, Fig. 1.4), at
  satellite potentials lower than predicted, and Stone & Raitt conclude
  that the system's motion through the plasma affects collection and
  must be modelled. A nonlinear contactor characteristic, fitted to these
  flight curves, is the next step for this example.
- **Reversibility.** PMG demonstrated in flight what the two columns of
  the results table assume: the same tether and contactors run as a
  generator with electron current down the tether and as a motor with
  electron current driven up it.

**The TSS-1R break.** On 25 February 1996 (057/01:29:16 GMT), during
STS-75, with 19.7 km of the planned 20.7 km deployed and the EMF at
3482 V, the tether parted at the Orbiter end. The failure board
(Szalai 1996) traced it to an arc through a flaw in the insulation that
burned the tether through, and Gilchrist et al. (1998) report that the
current, which had reached 0.97 A, rose to 1.1 A as the bare broken end
entered the plasma and stayed there for about 75 s: a few short strands
of copper emitted more current than the Orbiter's conducting surfaces,
by a gas-enhanced discharge later reproduced in a plasma chamber. For
this example the lesson is the one §3e applies to coils: every design is
checked against insulation breakdown at its peak voltage, and here the
peak is the full 1155 V between the tether and the surrounding plasma
when the far end is open.

## Files

| File | What it holds |
|---|---|
| `main.go` | the tether, its motional EMF, the two circuits and both ledgers |
| `main_test.go` | the energy and momentum identities above |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` (`KRONSPICE_EXPORT=1`) |

The results sheets record that automatic tearing finds no cut of this
single loop whose pieces each keep a reference (Aitchison 1983, Lemma 2),
so only the MNA column is filled.

## Sources

NASA flights and studies (all on the NASA Technical Reports Server,
ntrs.nasa.gov, by document ID):

- M. L. Cosmo, E. C. Lorenzini (eds.), *Tethers in Space Handbook*, 3rd
  ed., Smithsonian Astrophysical Observatory for NASA Marshall Space
  Flight Center, grant NAG8-1160, December 1997, NTRS 19980018321 —
  §1.1 TSS-1 and TSS-1R flights and findings; §1.3 the PMG flight (date,
  orbit, tether, hollow cathodes, bias and load sweeps, day and night
  currents, reversibility); §4 the electrodynamic tether equations (EMF,
  force, loop equation with sheath and ionospheric drops, 1–20 Ω
  ionosphere), hollow cathode versus passive and electron-gun contactors.
- L. Johnson, M. Herrmann, *International Space Station Electrodynamic
  Tether Reboost Study*, NASA/TM-1998-208538, Marshall Space Flight
  Center, July 1998, NTRS 19980223479 — about 100 V/km at the ISS orbit;
  10 km bare tether under 200 kg; 0.5 N from 5 kW and 0.8 N from 10 kW;
  0.43/0.7 N for a 7 km tether at 5–10 kW; the bare-tether anode (§6.2);
  Appendix A on ProSEDS (5 km bare wire, about 0.4 N, Delta II deorbit in
  about 17 days from 400 km, about 100 W extracted).
- N. H. Stone, W. J. Raitt, "The TSS-1R Electrodynamic Tether Experiment:
  Scientific and Technological Results," NASA Marshall Space Flight
  Center, *Adv. Space Res.* (1998), NTRS 19990084046 — EMF above 3500 V
  and current above 1 A; orbital motion affects current collection.
- B. E. Gilchrist, C. Bonifazi, S. G. Bilén, W. J. Raitt, W. J. Burke,
  N. H. Stone, J. P. Lebreton, "Enhanced electrodynamic tether currents
  due to electron emission from a neutral gas discharge: Results from the
  TSS-1R mission," *Geophys. Res. Lett.* 25(4), 437–440, 1998,
  doi:10.1029/97GL03023 — 19.7 km, 3482 V, 0.97 A before and 1.1 A for
  75 s after the break, 297 km altitude, the tether's construction, and
  the break time.
- K. Szalai, *TSS-1R Mission Failure Investigation Board: Final Report*,
  NASA, 1996 (as cited by Gilchrist et al.) — the insulation flaw and arc
  that cut the tether.
- L. Johnson, R. D. Estes, E. Lorenzini, M. Martínez-Sánchez,
  J. Sanmartín, I. Vas, "Electrodynamic Tethers for Spacecraft
  Propulsion," AIAA 36th Aerospace Sciences Meeting, Reno, January 1998,
  NTRS 19990107381 — propellantless thrust against a planetary field;
  SEDS-1, PMG and TSS-1 milestones; ProSEDS; Jupiter.
- J. Ballance, L. Johnson, "Propulsive Small Expendable Deployer System
  (ProSEDS)," NASA Marshall Space Flight Center, NTRS 20010020385 — the
  tether sections (10 km Dyneema, 5 km of seven-strand 28 AWG 1350-O
  aluminium, 0.8–1.6 mm, 250 N breaking strength), 400 km at 36°, about
  0.4 N, 14 km/day descent, the hollow cathode plasma contactor and the
  high-voltage switching between open, short, resistor box and battery.
  L. Curtis, L. Johnson, N. S. Brown, "Propulsive Small Expendable
  Deployer System (ProSEDS)," STAIF 2002, NTRS 20020048608 — the same
  mission as planned for 2002. ProSEDS was built but never flown.
- I. Katz et al., "Plasma turbulence enhanced current collection: Results
  from the plasma motor generator electrodynamic tether flight,"
  *J. Geophys. Res.* 100, 1995, doi:10.1029/94JA03142 — PMG's 500 m
  tether carried up to 0.3 A through the ionosphere; the contactor plasma
  scatters electrons across field lines.

Contactor and sheath theory the NASA documents rely on:

- J. R. Sanmartín, M. Martínez-Sánchez, E. Ahedo, "Bare wire anodes for
  electrodynamic tethers," *J. Propulsion and Power* 9(3), 353–360, 1993
  (reference 1 of NASA/TM-1998-208538) — current collection by an
  uninsulated tether in the orbital-motion-limited regime.
- L. W. Parker, B. L. Murphy, "Potential buildup on an electron-emitting
  ionospheric satellite," *J. Geophys. Res.* 72, 1631–1636, 1967 — the
  magnetically limited collection model that TSS-1R exceeded.

This project:

- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16 —
  Lorentz force by grade selection (reference only, §3f).
- P. W. Aitchison, "Diakoptics or tearing — a mathematical approach,"
  *Quart. Appl. Math.* 41, 265 (1983) — Lemma 2, quoted by the export
  sheet.
- `doc/em-push-pull.md` — design notes for pushing against fields in
  space; `doc/M9n.md` — the M9n checks this example belongs to.
