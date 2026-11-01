# One rotor in three gravitational fields

Rotation (special relativity) and curvature (general relativity) together,
for one machine: a 0.5 m T1000G carbon-fibre disk, 2 cm thick, spinning at
half its burst speed and carrying a 1 m² fibre loop (CLAUDE.md §3l,
"Rotation and curvature together"; package `spacetime/`). Then the reverse
question: how much does the machine itself curve and twist spacetime? The
rotor is the kind of machine [homopolar](../homopolar) and
[kron_machine](../kron_machine) solve electrically; this example asks only
what its clocks, its spin axis and a light loop riding it would show, and
what its material allows. The rotor energy-shift anomaly (H17) is the one
hypothesis on this bench and is computed in `spacetime.RotorShift` with
the standard k = 1/2 beside a declared k.

```
go run ./examples/machine_spacetime
go test ./examples/machine_spacetime ./spacetime
```

## The construction

CLAUDE.md §3l keeps three "curvatures" apart, and this example uses the
first two:

1. **Gravitational curvature.** The hub sits in the Schwarzschild field of
   the body at its site, moving with the site: clock rate
   √(1 − r_s/r − v²/c²), exact form in strong fields, with tidal
   stretching and the horizon and photon-sphere limits as typed errors.
2. **Rotation in flat spacetime.** In the hub's local inertial frame the
   rim is ordinary special relativity (the Langevin metric): rim clocks run
   at √(1 − ω²r²/c²), a loop on the rotor shows the Sagnac delay 4Aω/c²,
   and the rim speed is bounded by the material's burst speed √(σ/ρ).

By the equivalence principle the two compose: the exact rim rate equals the
first-order weak-field metric's rate to second order in the small
quantities (`TestRimRateTwin`). The spin axis precesses by Thomas +
geodetic + Lense–Thirring (Schiff 1960); the loop's Sagnac delay is taken
relative to the local gyroscopes, and the local frame's precession changes
it by the "local-frame correction" (`spacetime/machine.go`,
`AnalyzeMachine`).

**The rotor.** Toray T1000G fibre, 6,370 MPa and 1.80 g/cm³, bursts as a
thin hoop at √(σ/ρ) = 1,881 m/s (a solid disk bursts lower; the fibre
value is an upper bound for a composite rotor). The example runs at half
that, 941 m/s at 0.5 m, so ω = 1,881 rad/s; the disk's mass is
1800·π·0.5²·0.02 = 28.3 kg and its angular momentum ½MR²ω = 6.65×10³
kg·m²/s.

**The three sites.**

| Site | Body | Position and motion |
|---|---|---|
| Earth laboratory | Earth (GM, J2-free, rotating) | surface at 45° latitude, moving with the ground |
| Gravity Probe B orbit | Earth | circular polar orbit at 642 km above the equatorial radius, inertially pointed |
| Neutron star (illustrative) | 1.4 solar masses, 12 km radius, 10 ms spin, I = 0.4·M·R² (assumed values, not a measured object) | hovering 1,000 km from the centre, held up by a rocket |

Solar GM is the IAU 2015 Resolution B3 nominal value.

## Results

| | Earth lab, 45° | GP-B orbit, 642 km | Neutron star, 1,000 km (illustrative) |
|---|---|---|---|
| Hub clock rate | 0.999999999304050 | 0.999999999052363 | 0.997930583704828 |
| Rim clock rate | 0.999999999299128 | 0.999999999047441 | 0.997930583699917 |
| Rim loses against hub | 0.425 µs/day | 0.425 µs/day | 0.425 µs/day |
| GM/rc² | 6.95×10⁻¹⁰ | 6.32×10⁻¹⁰ | 2.07×10⁻³ |
| v²/c² at the rim | 9.84×10⁻¹² | 9.84×10⁻¹² | 9.84×10⁻¹² |
| Burst margin | 2.00 | 2.00 | 2.00 |
| Tidal stress on the disk (margin) | 0.0028 Pa (2.3×10¹²) | 0.0021 Pa (3.1×10¹²) | 3.3×10⁵ Pa (1.9×10⁴) |
| Thomas precession, here | 116 mas/yr | 0 | 0 |
| Geodetic precession, here | 350 mas/yr | 6,620 mas/yr | 0 (strong field, first-order terms indicative only) |
| Frame dragging, here | 172 mas/yr | 81.7 mas/yr | 4.9×10¹¹ mas/yr |
| Sagnac on the 1 m² loop | 8.372×10⁻¹⁴ s per pass | 8.372×10⁻¹⁴ s | 8.372×10⁻¹⁴ s |
| Local-frame correction to Sagnac | −1.0×10⁻³⁰ s | 4.5×10⁻²⁹ s | 0 |

Averaged over the GP-B orbit the model gives 6,621.0 mas/yr geodetic and
40.85 mas/yr frame dragging along Earth's axis, of which 39.1 mas/yr is
seen against the mission's guide star (IM Pegasi, declination +16°50′;
frame dragging lies along Earth's axis and GP-B measured its component
perpendicular to the star line). The rim lag is the same at every site
because it depends only on the rim speed; the hub rate carries the
gravity. At the neutron star the program warns that GM/rc² = 0.002 is a
strong field: the clock rates use the exact Schwarzschild form, and the
precession and Sagnac terms are first-order and indicative.

**The machine's own spacetime**, 1 m from the 28.3 kg rotor on its axis:

| | Rotor | Earth (surface, or at 45°) | Ratio |
|---|---|---|---|
| Clock shift GM/dc² | 2.1×10⁻²⁶ | 6.95×10⁻¹⁰ | 3×10⁻¹⁷ |
| Tidal curvature GM/c²d³ | 2.1×10⁻²⁶ /m² | 1.71×10⁻²³ /m² | 0.0012 |
| Frame dragging on the axis | 9.9×10⁻²⁴ rad/s = 6.4×10⁻⁸ mas/yr | 172 mas/yr | 3.7×10⁻¹⁰ |

A machine changes curvature through its mass and energy, and its spin
drags frames, but at laboratory scale the tidal curvature is the only part
within reach (torsion balances weigh that pull), and the frame dragging is
ten orders below Earth's own. The rotor's rim speed is limited by its
material long before any relativistic effect becomes large: at the burst
speed (v/c)² is 4×10⁻¹¹. No energy enters a circuit from spacetime
curvature; the gravitational channel of §3.1 carries clock rates and
frames, not power.

## Validation

The example's own test checks the self-field arithmetic (a 1 kg point mass
at 1 m shifts clocks by G/c² = 7.4261×10⁻²⁸, and a disk's J is ½MR²ω) and
that the program runs. The physics is held in `spacetime/`:

| Test | Check | Held to |
|---|---|---|
| `TestGravityProbeB` | geodetic 6,606.1 and frame dragging 39.2 mas/yr (GR, 642 km polar orbit) | 0.5 % and 3 % (spherical Earth, no J2); also within 2σ of the measured 6,601.8 ± 18.3 and 37.2 ± 7.2 |
| `TestRimRateTwin` | exact rim rate = weak-field first-order rate | second order in the small quantities |
| `TestSagnacEarthGyro` | 1 m² loop on the ground at 45°: 4Aω sin φ/c²; local-frame correction 10⁻¹¹ to 10⁻⁸ of it | 1e-12 |
| `TestGPSClock` | GPS satellite clock net +38.6 µs/day against the rotating Earth (model 38.54 without J2); gravitational part 45.7 | Ashby 2003 |
| `TestPoundRebka` | gh/c² over 22.5 m = 2.45×10⁻¹⁵ | Pound & Snider 1965 |
| `TestShapiroDelay`, `TestSolarDeflection` | closed form against the Plebanski index path integral; solar-limb deflection 1.751″ | Shapiro 1964; Will 2014 |
| `TestHorizonAndOrbitErrors`, `TestTidal` | r = 3r_s static rate √(2/3); typed errors at the horizon and inside the photon sphere; Earth-surface tidal gradient 3.08×10⁻⁶ s⁻² | exact |
| `TestRotorShift`, `TestSagnac`, `TestBurstAndRotor` | rotor shift with k = 1/2 equals exact rim dilation (H17's k beside it); exact Sagnac equals 4Aω/c²; T1000G burst speed | exact |

## Files

| File | What it holds |
|---|---|
| `main.go` | the rotor, the three sites, the per-site report and the self-field report |
| `main_test.go` | self-field known values; the program runs |
| `run_output.txt` | the run quoted above |
| [spacetime/machine.go](../../spacetime/machine.go) | `AnalyzeMachine`, `OrbitAveragedPrecession`, `MachineSelfField` |
| [spacetime/rotation.go](../../spacetime/rotation.go) | Langevin rim clocks, Sagnac, materials and `BurstSpeed` |

## Sources

- C. W. F. Everitt et al., "Gravity Probe B: final results of a space
  experiment to test general relativity," *Phys. Rev. Lett.* 106, 221101
  (2011) — the measured geodetic and frame-dragging rates.
- N. Bartel et al., arXiv:1509.07529, Table 1 — the guide star IM Pegasi's
  position (Hipparcos), used to project the frame dragging as GP-B saw it.
- L. I. Schiff, "Possible new experimental test of general relativity
  theory," *Phys. Rev. Lett.* 4, 215 (1960) — gyroscope precession.
- N. Ashby, "Relativity in the Global Positioning System," *Living Rev.
  Relativ.* 6, 1 (2003) — clocks near Earth.
- C. W. Misner, K. S. Thorne, J. A. Wheeler, *Gravitation*, Freeman, 1973
  — Schwarzschild clocks, orbits and tidal forces.
- J. Plebanski, *Phys. Rev.* 118, 1396 (1960); L. D. Landau, E. M.
  Lifshitz, *The Classical Theory of Fields*, §90 — gravity as an
  effective index, used by the Shapiro-delay twin.
- C. M. Will, *Living Rev. Relativ.* 17, 4 (2014) — PPN γ and the
  deflection and delay bounds.
- G. Sagnac, *C. R. Acad. Sci.* 157, 708 (1913); E. J. Post, *Rev. Mod.
  Phys.* 39, 475 (1967) — the Sagnac delay.
- J. Bailey et al., *Nature* 268, 301 (1977) — the clock hypothesis at
  10¹⁸ g in the CERN muon ring; B. Mashhoon, *Phys. Lett. A* 145, 147
  (1990) — the locality hypothesis behind the rim's co-moving frame.
- IAU 2015 Resolution B3 — nominal solar GM.
- Toray Composite Materials America, T1000G data sheet, 2024-02-16 —
  6,370 MPa, 1.80 g/cm³; S. Timoshenko, *Strength of Materials II*, 1956,
  §8 — the thin-hoop burst speed.
