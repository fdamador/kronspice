# Plasmoids as toroidal current rings (H26)

J. Nachamkin, "Force-Free Time-Harmonic Plasmoids," Phillips Laboratory
report PL-TR-92-3044 (University of Dayton Research Institute, 1992).

Nachamkin proposes a plasmoid built from three parts:
- a sphere of plasma holding a standing electromagnetic wave whose
  electric and magnetic fields are parallel everywhere ("force-free", or
  Beltrami);
- a swirling plasma flow that balances the field's pressure;
- a "resonant size" at which it should stop radiating.

He offers it as a model of ball lightning, with a worked example: 1 MJ in
a 10 cm sphere.

Seen as rings:
- **Toroidal structure.** The field lines of his l = 1 field lie on nested
  tori around a ring of current: a toroidal structure with an empty core.
- **Spheromak.** At zero frequency the same field is the **spheromak**, a
  laboratory plasma ring.
- **Collapsing ring.** A ring of electrons squeezed by a rising magnetic
  field is the electron-ring compressor of the electron-ring accelerators.

This example tests the paper and models both rings with kronspice. It sits
beside the other electron-ring models: [kron_electron](../kron_electron)
(H21), [space_vortex](../space_vortex) (H13) and
[photon_electron](../photon_electron) (H23).

```
go run ./examples/plasmoid_rings -viz examples/plasmoid_rings/viz   # about 8 seconds
python3 examples/plasmoid_rings/viz/render.py
go test ./examples/plasmoid_rings
```

## 1. The field

His l = 1 field (eq. 8, identical to the spheromak field):

- **Force-free.** It obeys ∇×F = αF to 2×10⁻¹⁰ and has zero divergence.
- **Nested tori.** A traced field line keeps its poloidal flux to 10⁻¹⁴,
  so it stays on its torus.
- **The ring.** The tori nest around a current ring at 0.611 of the
  sphere's radius (61 mm for a 10 cm plasmoid).

Picture: left panel of `viz/plasmoid_rings.png`.

## 2. His critical frequency is ordinary plasma physics

His Appendix B, recomputed (10 cm radius, k = 2686 /m):

| Quantity | Here | Nachamkin |
|---|---|---|
| Free charge density | −32.63 C/m³ | −32.63 C/m³ |
| Free electrons | 8.532×10¹⁷ | 8.532×10¹⁷ |
| Field for his C = 0.01 at 126 GHz | 3.566×10⁹ V/m | 3.566×10⁹ V/m |

- **Nothing exotic.** His condition on the electron density is exactly
  ω² = ωp² + α²c² (to 10⁻¹⁶), the textbook dispersion relation of an
  electromagnetic wave in a cold plasma. His "critical frequency" αc is
  that wave's cutoff, and the plasmoid is a standing wave of an ordinary
  plasma.
- **His own example breaks his assumptions.** At his chosen field the
  electrons would oscillate at ω·C·R = 7.9×10⁸ m/s, **2.6 times the speed
  of light**. That breaks the non-relativistic electron motion his whole
  analysis assumes. The field is also about 1,000 times the strength that
  breaks down air.

## 3. Does the "resonant size" stop radiation?

The exact l = 1 resonances of a plasma sphere with his radius and density,
with nothing feeding it from outside:

| Part of his field | Frequency | Size it needs (αR) | Q | Energy lifetime |
|---|---|---|---|---|
| TE (his M part) | 128.158 GHz | 4.493 — his choice | 4.8×10⁵ | **0.6 µs** |
| TM (his N part) | 128.147 GHz | **2.744** | 9.4×10⁵ | **1.2 µs** |

- **It radiates.** The plasmoid's energy leaks away in about a microsecond:
  about 2×10¹² W at 1 MJ.
- **The two halves need different sizes.** An equal mix of them — the mix
  that makes E parallel to H — is not one free mode of the sphere. Holding
  it needs outside driving or surface currents, the weak point Nachamkin
  names himself.

## 4. The spheromak ring in kronspice

At zero frequency his field is a force-free current ring. There is no
Lorentz force on the current, so the ring neither pinches nor bursts.
Electrical resistance ends it.

- **Decay time.** The field decays as one mode with lifetime μ0/(ηα²), with
  Spitzer resistivity η.
- **As a circuit.** The ring is an inductor L = 2W/I² carrying its ring
  current I, in series with a resistor R = P/I².

For 1 MJ in a 10 cm sphere: field scale 80 T, ring current 7.4 MA, and a
magnetic pressure of 25,000 atmospheres that only a wall or an outside
field can hold.

| Electron temperature | Ring R | Energy lifetime (theory) | kronspice | Heat |
|---|---|---|---|---|
| 1 eV | 30 mΩ | 0.60 µs | 0.60 µs | all 1 MJ |
| 10 eV | 0.96 mΩ | 18.9 µs | 18.9 µs | all 1 MJ |
| 100 eV | 30 µΩ | 0.60 ms | 0.60 ms | all 1 MJ |

To last one second the ring would need 14 keV — 163 million K, fusion
conditions. A ball-lightning plasmoid lasting seconds in air cannot be a
resistive spheromak. Picture: middle panel.

## 5. A collapsing electron ring, driven by kronspice

The setup (declared values):
- **Bank and coil.** A 1 mF, 10 kV capacitor bank is discharged into a
  10 µH solenoid with 16 turns per metre; kronspice solves this circuit.
- **The ring.** 2 MeV electrons circle the axis in a 0.05 T bias field.
- **Compression.** As the coil field rises the induced electric field
  speeds the electrons up and the stronger field bends them tighter: the
  ring collapses.

| | Start | End | Law |
|---|---|---|---|
| Field | 0.050 T | 1.984 T | |
| Ring radius | 16.4 cm | 2.60 cm | R ∝ B^−½: 2.60 cm |
| Momentum | 1 | ×6.3007 | ×6.2989 |
| Electron energy | 2.0 MeV | 15.0 MeV | |
| Ring current (10¹³ electrons) | 456 A | 2,933 A | |

- **Why it follows the law.** The electrons' angle around the axis is a
  cyclic coordinate (the same idea as in [kaluza_kron](../kaluza_kron)),
  so their canonical angular momentum is conserved, to 6×10⁻⁴ here.
- **Ledger.** The induced field's work equals the energy gained exactly,
  and kronspice's ledger closes to 10⁻⁹. The ring takes 21 J of the bank's
  50 kJ.

The "void" inside the ring is field, not a source: every joule the ring
gains comes from the bank. Picture: right panel.

## Where it stands

- **Sound pieces.** Nachamkin's field mathematics is sound, and his
  electron density follows from standard plasma physics.
- **Two failures.** His worked example puts the electrons above the speed
  of light, and the "resonant size" does not stop radiation: a free
  plasmoid of this kind lasts about a microsecond.
- **His suggested ingredients.** Long-lived ball-lightning plasmoids, if
  they are electromagnetic, need something the paper lacks; Nachamkin
  himself suggests quantum effects and outside fields.
- **Real rings.** The real, measured force-free rings are spheromaks
  (lifetimes of µs to ms in the lab) and compressed electron rings — both
  modelled here with kronspice.

## Files

| File | What it holds |
|---|---|
| `beltrami.go` | the l = 1 force-free field, curl check, field-line tracing, ring radius |
| `plasma.go` | his Appendix B and the cold-plasma dispersion relation |
| `radiation.go` | the plasma sphere's leaky TE and TM resonances |
| `spheromak.go` | the spheromak ring and its resistive decay in kronspice |
| `ring.go` | the electron-ring compressor: kronspice bank and coil, relativistic electron |
| `main.go` | the report |
| `main_test.go` | field, Appendix B, modes, spheromak decay, ring compression |
| `viz/render.py` | `plasmoid_rings.png` |

## Sources

- J. Nachamkin, *Force-Free Time-Harmonic Plasmoids*, PL-TR-92-3044 (Phillips Laboratory, 1992), DTIC AD-A257 765.
- S. Chandrasekhar, P. C. Kendall, *Astrophys. J.* 126, 457 (1957) — force-free fields.
- M. N. Rosenbluth, M. N. Bussac, *Nucl. Fusion* 19, 489 (1979) — the spheromak.
- F. F. Chen, *Introduction to Plasma Physics and Controlled Fusion* — electromagnetic waves in a cold plasma.
- J. A. Stratton, *Electromagnetic Theory* (1941), §9.22 — spherical boundary-value problems.
- NRL Plasma Formulary — Spitzer resistivity.
- V. I. Veksler, *Proc. CERN Symposium on High Energy Accelerators* (Geneva, 1956) — coherent (collective) acceleration, the idea behind electron-ring compressors.
