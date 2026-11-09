# Literature: rotational cavitation, cavity flows, and the Mandelbulb

Collected 2026-11-08 for sizing the Mandelbulb cavity (see CLAUDE.md). Values
marked *derived* are computed here from a paper's stated geometry, not stated
by its authors. Check every number against the full paper before it goes
into code; several were read from abstracts or from PMC full text only.

## 1. Rotational hydrodynamic cavitation reactors (the closest machines)

A rotor with cavities (dimples, grooves, "cavitation generation units",
CGUs) turns past a stator, often with matching cavities, across a gap of a
few mm. Cavitation forms at the cavity walls where the water is
accelerated, and collapses as cavities move apart and come together again.
These are the machines whose numbers size our bulb.

| Paper | Geometry | Speed | Result |
|---|---|---|---|
| Xue, Hao, Ren et al., *Ultrason. Sonochem.* 109, 106999 (2024), [PMC11295464](https://pmc.ncbi.nlm.nih.gov/articles/PMC11295464/) | Cylindrical rotor Ø 80 mm, 20 rectangular grooves 6 mm wide × 3 mm deep, 18° oblique teeth, stator Ø 90 mm (gap 5 mm) | 3000–4200 rpm | No bubbles at 3000 rpm, cavitation from 3400. Optimum pulsation amplitude at 3800 rpm; more speed does not help. Attached cavity on the groove's **front wall**, cloud moves rearward and collapses. Dominant frequency 1122 Hz at 3600 rpm, rising linearly with speed; St = fW/u ≈ 0.47–0.48. *Derived:* surface speed 14.1 m/s at 3400 rpm, 15.9 m/s at 3800 rpm; groove-passing frequency 20 × 60 Hz = 1200 Hz at 3600 rpm |
| Sun, Xuan, Song et al., *Ultrason. Sonochem.* 70, 105311 (2020), [PMC7786598](https://pmc.ncbi.nlm.nih.gov/articles/PMC7786598/) | Cone-cylinder CGUs (Ø ≈ 10 mm, *derived*), 32 on rotor and on each cover, gap 2 mm | 3600 rpm, 1.4 m³/h | Local rotor speed 46 m/s. Sheet cavitation grows as moving and static units part ("leaving", "approaching") and breaks up when they line up. Period 0.5 ms (≈ 32 × 60 Hz). Schnerr–Sauer cavitation model, SST k-ω |
| Sun, Yang, Wei et al., *Ultrason. Sonochem.* (2021), doi:10.1016/j.ultsonch.2021.105771, [PMC8551246](https://pmc.ncbi.nlm.nih.gov/articles/PMC8551246/) | **GA + CFD** on the CGU: diameter D 9–11 mm, interaction distance s 1.5–2.5 mm, height h 1.5–2.5 mm, inclination θ 6.25–18.75° | 3600 rpm | Objectives: maximize vapour volume, minimize rotor torque (NSGA-II type). Pareto: s → 1.5 mm and θ → 18.75° always; D and h trade vapour against torque. Vapour 5.44 → 20.8 × 10⁻⁸ m³ (torque 1.63 → 1.94 N·m); a balanced point gives +156 % vapour at −14 % torque. Vapour volume is a proxy, not intensity: needs experiment |
| Kwon & Yoon, *Proc. IMechE E* 227, 67–73 (2012), doi:10.1177/0954408912451535 | Rotational cavitation **heat generator** | variable | Thermal efficiency up to ~90 % (heat into water / electrical input). Cavitation needs inlet > 1.25 bar and > 2700 rpm (as cited by Zhang et al. 2025) |
| Yoon group pilot reactors (cited in Zhang et al. 2025) | Hybrid rotor | 3600 rpm | Heat 48.15 MJ/h at 82.18 % thermal efficiency; pilot scale 57.21 MJ/h at 79.99 %. Efficiency **falls** as water temperature rises |
| Badve, Gogate, Pandit, Csoka, *Sep. Purif. Technol.* (2013–14) | Rotor with indentations; KI dosimetry | 2200–2700 rpm | Peak cavitation activity near 2200–2500 rpm; higher speed stabilizes the flow and lowers activity |
| Petkovšek, Zupanc, Dular et al., *Sep. Purif. Technol.* 118, 415–423 (2013) | Rotation generator, two rotors with opposite radial grooves | — | Pharmaceutical removal (with H₂O₂). An 8° inclined rotor gives a larger cavity zone, stronger pressure pulses, more erosion; pressure amplitude 0.7 bar at 3.5 mm gap against 1.7 bar at 0.8 mm (as cited by Sun 2020) |
| Zhang, Zheng, Ezekoye et al., *Ultrason. Sonochem.* 124, 107727 (2025), [PMC12827028](https://pmc.ncbi.nlm.nih.gov/articles/PMC12827028/) | **Review** of rotational cavitators | — | Dimple diameter drives intensity (11–19 mm widened low-pressure zones in CFD); six cavitation-number and three rotational-Reynolds definitions in use (Tables 1–2); no single optimum speed |
| Review, *Processes* 11(2), 514 (2023), [MDPI](https://www.mdpi.com/2227-9717/11/2/514) | Rotary generators for wastewater | — | Speed range reported ideal for one design 4320–5760 rpm; a critical speed exists below which no cavitation forms |
| *Fluids* 9(10), 244 (2024), [MDPI](https://www.mdpi.com/2311-5521/9/10/244) | Optimization of rotational reactor geometry (CFD) | — | At low speed cavitation starts at points in the narrow channel; higher speed lengthens the cavity zone. Not yet read in full (403) |

**What these give us (design rules, to be fitted, not assumed):**
1. **Inception speed.** Cavitation starts when the relative speed at the
   cavity edge is about **14–16 m/s** with the outlet near 1 atm (Xue 2024),
   or above 2700 rpm with inlet > 1.25 bar (Kwon & Yoon). Working designs run
   at 15–46 m/s.
2. **Optimum, not maximum.** Each reactor has a best speed (2200–2500,
   3800 rpm). Beyond it the flow stabilizes or the cavity chokes (supercavity)
   and the collapse intensity falls.
3. **Cavity size.** 5–20 mm across, 1.5–17 mm deep. Larger diameter gives more
   vapour and more torque.
4. **Gap.** Smaller is stronger: 1.5 mm (Sun 2021 optimum); pulses
   1.7 bar at 0.8 mm against 0.7 bar at 3.5 mm.
5. **Inclination.** A leading edge inclined about 8–19° helps (Petkovšek;
   Sun 2021).
6. **Frequency.** The pulsation follows the cavity-passing frequency N·Ω/2π
   (Xue: 1122 Hz near 20 × 60 Hz; Sun: 0.5 ms with 32 units at 60 Hz).
   St = fW/U ≈ 0.47.
7. **Heat.** Measured thermal efficiencies are **80–90 %** of the
   electrical input, below 100 %. Efficiency falls as the water gets hotter.

## 2. Cavity flows and cloud shedding (the physics of one cavity)

- **Open cavity, inception at the trailing corner.** X. Liu & J. Katz,
  *Phys. Fluids* 20, 041702 (2008); cavity L = 38.1 mm, D = 30.0 mm. The
  lowest pressure, and so inception, occurs periodically above the trailing
  corner as shear-layer vortices hit it. LES: Shams & Apte (Cav2009,
  [PDF](https://web.engr.oregonstate.edu/~sva/archive/Shams_Cav09_symp.pdf));
  Liu & Katz flow physics [arXiv:1110.3978](https://arxiv.org/pdf/1110.3978).
- **Cavity oscillation (Rossiter).** Incompressible limit St = fL/U ≈
  κ(m − γ) with κ ≈ 0.57, γ ≈ 0.25 (Rossiter 1964; see the generalized form,
  [arXiv:2601.09099](https://arxiv.org/html/2601.09099)). Axisymmetric cavity
  in water: Gharib & Roshko, *J. Fluid Mech.* 177, 501 (1987). In-plane
  modes of each lobe of the bulb follow this.
- **Cloud shedding.** St ≈ 0.2–0.3 for orifices and nozzles, 0.28 (Stutz &
  Legoupil); re-entrant-jet shedding at higher σ, shock-driven below about
  σ ≈ 0.8. Hatzissawidis, Kuhr & Pelz, *J. Fluid Mech.* (2025),
  "the transition from re-entrant flow driven to shockwave driven cloud
  cavitation". Shedding frequency falls as σ falls.
- **Vortex cavitation.** Arndt, *Annu. Rev. Fluid Mech.* 34, 143 (2002).
- **Textbook.** C. E. Brennen, *Cavitation and Bubble Dynamics*, OUP 1995
  (ch. 3 inception, ch. 6 Wood's equation, ch. 8 cavitating flows).

## 3. The nearest commercial heater

- J. L. Griggs, Hydrosonic Pump, later "ShockWave Power Reactor" (Hydro
  Dynamics Inc.; US 5,188,090, 1993: check the patent text). NASA Spinoff
  [Industrial Productivity](https://spinoff.nasa.gov/node/10190) and
  [NTRS 20020066442](https://ntrs.nasa.gov/citations/20020066442): Marshall
  engineers redesigned its bearings, which overheated. The company described
  it as converting mechanical energy to heat with high efficiency.
  No independent calorimetry found; claims of > 100 % are unverified
  (forum sources only).

## 4. The Mandelbulb

- White–Nylander power-n map: v ↦ vⁿ + c with
  vⁿ = rⁿ(sin nθ cos nφ, sin nθ sin nφ, cos nθ), θ = arccos(z/r),
  φ = atan2(y, x). Conventions differ (some use arcsin), so state the one
  used. [Wikipedia](https://en.wikipedia.org/wiki/Mandelbulb),
  [MathWorld](https://mathworld.wolfram.com/Mandelbulb.html),
  [P. Bourke](https://paulbourke.net/fractals/bulb),
  [arXiv:2206.06332](https://arxiv.org/pdf/2206.06332) (algebraic
  foundation).
- Distance estimate for rendering and meshing: D. Hvidtfeldt,
  ["Distance Estimated 3D Fractals (V)"](https://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/),
  DE = ½·r·ln r/|dr|.
- History: White and Nylander, 2007–2009; an earlier version attributed to
  J. Ruis (1997). Sources disagree.
- Properties tested in code (`model_test.go`): both conventions have
  (n − 1)-fold symmetry about the z-axis (turning c by 2π/(n − 1) turns
  every iterate by the same angle, since n·α ≡ α). Only the latitude
  convention (θ = arcsin z/r) has the degree-n Multibrot set as its
  equator. In the polar convention (θ = arccos z/r) an equatorial point
  maps toward the pole. So n sets the number of main lobes round the spin
  axis, 7 at n = 8, in both.

## 5. Still to find

- Liu & Katz 2008 full text (inception index values).
- Kwon & Yoon 2012 full text (heater geometry and calorimetry).
- Sun 2021 inception measurements to calibrate the vapour-volume proxy.
- Any cavitation study of fractal or multiscale cavity shapes (none found
  yet).
- Griggs patent text and any third-party calorimetry.
