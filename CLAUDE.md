# CLAUDE.md — kronspice

A SPICE-compatible analog circuit simulator written in pure Go, built to be the
fitness-evaluation engine for genetically evolved circuits.

This file tells Claude how the project is organized, which libraries and
references to use, and the rules to follow when changing code.

---

## 1. Project goals

1. **Fast, embeddable simulation.** Simulate a circuit in memory from a Go
   struct or a netlist string. No separate program to launch and no result
   files to read back for each candidate.
2. **Safe for evolved, often broken circuits.** A genetic algorithm produces
   floating nodes, shorted sources, singular matrices and circuits that never
   converge. The simulator must *return an error* for every one of these and
   never panic, hang or exhaust memory.
3. **Massively parallel.** Thousands of independent simulations run across
   goroutines with no shared mutable state.
4. **SPICE3-compatible where it matters.** Netlist syntax, device equations
   and analysis results match ngspice within stated tolerances, verified by
   automated golden tests.
5. **Deterministic.** The same netlist, options and GA seed always give
   bit-identical results.
6. **Consistent with Maxwell's electromagnetism in quaternion form.** Every
   circuit law and device model must be derivable from Maxwell's equations as
   he wrote them with Hamilton's quaternions, in the lumped (quasi-static)
   limit, and the engine must check that it is operating inside that limit
   (§3c).

### Long-term vision: neutrinos as an energy source

The end goal of this project is to genetically evolve circuits — for craft
travelling between planets, stars and galaxies — that use neutrinos as an
energy source. This is a research aspiration, not something any known circuit
can do today: neutrinos interact only through the weak force and gravity, and
the power a circuit can absorb from them under established physics is
vanishingly small (§3.1, neutrino channel).

The engine serves this vision by measuring it honestly:

1. **Neutrino inflow is always computed and reported** from the established
   weak-interaction model in `exchange/neutrino.go` (flux × cross-section ×
   energy deposited, with cited sources), a coupling model on the neutrino
   port that every element carries (§3.1), for every circuit and run. The
   flux at the circuit is a declared, cited input (default: the solar flux
   at 1 AU).
2. **The GA may optimize for neutrino-derived inflow** as an objective, and
   results are reported exactly as computed — including when they are
   effectively zero.
3. **No unvalidated coupling.** Claude must not add a neutrino coupling
   mechanism that isn't supported by published experimental evidence. If such
   evidence appears, the new interaction is added as a separate, cited model
   in `exchange/` with its own known-value tests.
4. **Surplus is never credited to neutrinos.** Any energy-ledger residual is
   investigated as an undeclared input, model error or numerical error
   (§3c) — evolution exploits simulator flaws, and a flaw must never look
   like progress toward this goal.
5. **Near-term work toward the vision:** resonant energy-storage circuits such
   as Tesla's bifilar coil (§3e), ultra-low-power operation, and harvesting
   of every real, declared energy source.

### Hypotheses under test

Ideas that go beyond established physics are kept in the project, stated
plainly and tested, not dropped. Each one is listed here with where it is
modeled, what would count as evidence, and its current status. Rules:

- **Nothing is removed for being unconventional, and nothing is promoted
  without a test.** A hypothesis leaves this list only when it has been
  tested: either a deviation from the standard prediction is found and
  reproduced, or a null result sets a bound — and even then it stays, with
  that bound recorded.
- Every hypothesis has a **standard prediction computed alongside it**, so a
  bench result can be read against both.
- Hypothetical models are off by default, labeled in all output, and keep
  their own ledger columns (§3c, §3f).
- **Status values:** *untested* · *bounded* (null result, limit recorded) ·
  *deviation found* (with data and conditions) · *reproduced*.

| # | Hypothesis | Modeled / tested in | Status |
|---|---|---|---|
| H1 | Ether-stream premise: pulses send outward HV electric fields that create toroidal magnetic fields along them | §3g pulse-duration experiment, checks 1–7 | untested |
| H2 | Magnetic charge, dual-symmetric Maxwell | `em/duality.go`, §3f | untested |
| H3 | Field-generated magnetic monopoles (Giardino 2020) | `em/hed.go`, §3f; pulse experiment | untested |
| H4 | Temporal scalar field linking EM pulses and heat (Jack 2003) | `em/temporal.go`, §3f; §3e Bi–Al test; pulse check 2 | untested |
| H5 | Water mist speeds spark-gap recovery and sharpens the pulse train | §3j misted gaps | untested |
| H6 | Solar energy and formation through high-voltage pulse fields | see below | untested |
| H7 | Neutrinos as a usable energy source | §1 long-term vision, §3.1 neutrino port | computed ≈ 0 under established physics |
| H8 | Extended electrodynamics with a physical scalar field (van Vlaenderen 2003; Hively 2012) | `em/scalar/`, §3k tier 2 | untested |
| H9 | Scalar longitudinal waves observed (Monstein & Wesley 2002) — replication | §3k tier 3, simulated setup and bench repeat | untested |
| H10 | Bearden: scalar EM; Motionless Electromagnetic Generator, COP > 1 | §3k tier 3, MEG replica with full ledger | untested |
| H11 | Meyl: longitudinal vortex scalar waves | §3k tier 3, delay, distance and resonance tests | untested |
| H12 | Electrogravitic thrust (Biefeld–Brown) beyond ion wind | §3k tier 3, air versus vacuum thrust | untested |
| H13 | Tewari's space-vortex electron: the electron as a vortex of an incompressible space fluid around a void (Physics Essays 31, 108, 2018); his reactionless generator | `examples/space_vortex` electron test bench; `examples/space_power_generator` ledger of his published SPG test; `doc/tewari-references.md` | untested — Zeeman prediction needed; model gives g = 15/8 (15/16 with the standard ring moment) against 2.0023, proton charge 12e, α = 15π²/16 against 1/137, a density limit white dwarfs exceed; SPG test (`examples/space_power_generator`) rests on an unmeasured loaded EMF: ≤ 1.37 V predicted, 3.41 V assumed; bearing (Huber) torque could close its 24.8 N·m torque gap only in EM mode, whose back-EMF leaves the same 1.37 V, and is excluded in thermal mode at the measured η ≈ 0.10 (`examples/bearing_motor`) |
| H14 | Marinov: a running ball-bearing motor yields heat beyond its electrical input (≈ 8 % more with rotation, *Electronics & Wireless World*, April 1989) | `examples/bearing_motor` registry; calorimetric repeat with the P element's ledger as the standard prediction (heat = V·I − ω·T) | untested — no data published |
| H15 | ECE theory (Evans; A.I.A.S. UFT paper 107): homopolar resonances drawing energy from spacetime; Eckardt 2012 "giant resonance" of a self-excited twin-disk machine at R_eff = 0 | `examples/eckardt_homopolar`: paper's equations and the kronspice machine element side by side, full ledger | untested — the classical resonance (R_eff = 0) is the series-dynamo self-excitation threshold and is reproduced; the paper's energy surplus equals, to 1e-3, the work of the back torque its eq. (12) omits; with it the current stays bounded (≈ 9.9 A at ω_res, 30 W from the drive) and the ledger closes; the two ECE-only resonances of UFT 107 are not modeled (no equations supplied) |
| H16 | Spin connection resonance (AIAS): a parametric LCR circuit with C(t) draws energy "from the time domain" (Amador & Eckardt 2012, COP ≈ 200 by source input) | `examples/parametric_lcr`: all of the paper's cases, with the power of the agent that varies C on the ledger (Mandelstam–Papalexi; K&V §10.2) | untested — reproduced: U0/R = 0.25 A classical case; unbounded growth at 2ω0 pumping (σ = 245 1/s); bounded ≈ 16 A, 2.4 kW at 0.97 f0, φ = −90° (source-only ratio 169); §3.2 pumping at ω0 grows only for a0 ≳ 0.7. In every case the source plus the pump pays for R and the stored energy (books close to ≤ 0.4 % where resolved); the COP omits the pump. Bench test: measure the power that drives C(t) (motor shaft or varactor bias supply) |
| H17 | Rotor Mössbauer anomaly: the energy shift between centre and rim of a fast rotor is ΔE/E = k·u²/c² with k ≈ 0.68, not the relativistic 1/2 | `spacetime.RotorShift` with k as a parameter; standard k = 1/2 beside it (§3l) | deviation reported, not independently reproduced — k = 0.68 ± 0.03 (Kholmetskii et al., Phys. Scr. 79, 2009), 0.69 ± 0.02 (Yarman et al., Can. J. Phys., doi:10.1139/cjp-2015-0063), 0.596 ± 0.006 from a reanalysis of Kündig 1963 (Phys. Scr. 77, 2008) against Kündig's own agreement with 1/2 within 1.1 %; interpretation disputed: clock desynchronization giving k = 2/3 (Corda, arXiv:2203.02282) against slit-broadening artifacts (Friedman & Nowik, Phys. Scr. 85, 065702, 2012). Needs an independent group with full line-shape measurement |
| H18 | Inertial mass reduction of a craft by a "polarized vacuum" from a charged, vibrating, microwave-driven resonant-cavity hull (Pais, US 10,144,532 B2, 2018) | `examples/pais_cavity`: the annular cavity as a ring line solved by Happ two-level diakoptics (twin MNA), ledger, and the patent's σ, flux, rotation and vibration figures against physical limits | untested — standard prediction: no net force on a closed cavity, mass unchanged, thrust ≤ P/c (3.3 mN per MW); the patent's σ = 1 C/m² needs 1.1×10¹¹ V/m (wall pressure 5.6×10¹⁰ Pa, field emission drains it in ~10⁻¹⁶ s); its 10³³ W/m² needs wall speed R_v·ν ≈ 3×10¹³ c; 10⁹–10¹⁸ Hz wall vibration lies above copper's Debye frequency (7.1×10¹² Hz); 30,000 rpm at 2 m exceeds carbon-fibre burst speed 3.3×; the cited gyroscope weight loss (Hayasaka & Takeuchi 1989) was not reproduced (Faller et al. 1990; Quinn & Picard 1990; Nitschke & Wilmarth 1990) |
| H19 | Bopp–Podolsky generalized electrodynamics: higher derivatives in the field equations themselves, (1 + a²□)□A = μ0J (fourth order), one new length a (F. Bopp 1940; B. Podolsky, Phys. Rev. 62, 68, 1942) | `em/podolsky.go` (point-charge potential and field, finite self-energy q²/(8πε0a), Maxwell and massive modes, Coulomb-test bound a < r/ln(1/ε)); `examples/rate_of_change` | untested — standard prediction a = 0 (Maxwell) beside every result; potential falls below Coulomb by e^(−r/a), so only very short-range experiments see it; electron self-energy equals m_e c² at a = 1.409×10⁻¹⁵ m; a lower bound on the Podolsky mass from measured Bhabha scattering is given by Bufalo, Pimentel & Soto, Phys. Rev. D 90, 085012 (2014) — value to be entered from the paper |
| H20 | Kron's polyhedral waves describe a physical field region: a region of space, with or without plasma, "crystallizes" into Kron's sequence of magnetic (even-dimensional) and dielectric (odd-dimensional) p-networks plus their dual, carrying complete e, h, b, d waves across dimensions, and can take a self-adaptive oscillatory state; crystals realize it (Kron, *J. Math. Anal. Appl.* 5, 46–69, 1962, pp. 46, 50, 53–54, 66–68; stated as physics in his 1958 and 1960 texts, where the waves are "accompanied by … magnetic poles" and the plasma field is "the stage for all cognitive activities", as quoted by Balasubramanian 1968 pp. 8, 15–17) | `kron/polyhedron.go`, `kron/curvefit.go`, `examples/kron_polyhedron`; standard twins `em/fdtd.go`, `em/tlm.go`; `reference/PII_ 0022-247X(62)90005-7.pdf`; `doc/kron-figure-catalogue.md` | untested — physical claims have no test yet; Kron himself called the physical interpretation "only tentative" (p. 66); the 2-D lattice and first planar frame reproduce his numbers, and his oscillatory column equals its input (fits any values, random numbers included), so its fit is no evidence for the physical reading; standard prediction: in three dimensions the polyhedron is a discretization of Maxwell's equations (Kron 1944; Tonti's cell method) and agrees with FDTD and TLM to discretization error, a passive linear structure has fixed modes, and X-ray-excited crystal waves follow Ewald–Laue dynamical diffraction |
| H21 | Kron's machine equations model the ring electron: a closed ring of circulating charge is a one-winding machine (self-inductance L, circulating mass as rotor inertia, speed voltages G), with voltage and current the pressure and flow of energy at that scale; Kron's eq. 11 (1951 Epilogue p. 259) is the formal link, the machine equation with e = R = 0 being the motion of a charged particle in an electromagnetic field (user hypothesis, 2026-10-31) | `examples/kron_electron` on the electron test bench of `doc/kelvin-references.md`, beside `examples/kelvin_ring` (Lorentz-force form) and `examples/space_vortex` (H13); `device/machine.go`, `em/inductance.go` | bounded — step 1 passes: Kelvin's ring as a three-winding machine with the inertia in the metric gives ∓eΦ/2π, Larmor's rate and the bench's kinetic energy to 1e-5, ledger residual 1e-13, and frame I's torque ½iᵀ(∂L/∂θ)i is −μB sin θ (period to 2e-5); step 2 computed: a ring whose mass circulates with its charge has g = 1/(1 + Le²/4π²ma²) ≤ 1 and precesses at Larmor's eB/2m, half the electron's; a massless ring carries the measured moment at S = ħ/2 only at a = 8.9 fm with v = 43.6 c (v/c = π/(2α(ln(8a/r) − 7/4)) for every massless thin ring, set by the fine-structure constant), magnetic energy 475 m_e c²; standard prediction: the Dirac electron with QED (g = 2.002 319 304 36, μ = 1.001 16 μ_B); scattering bounds the electron below 10⁻¹⁸ m, and Kelvin's full ring already gives the Zeeman lines with the wrong polarisation (`examples/kelvin_ring`) |
| H22 | "Quantumtricity": the quark model and the principles of circuit electricity (V, I, R, P) are one structure, the latter being electricity at the quantum level (user hypothesis, 2026-11-01) | `kron/gauge.go` and `examples/quantumtricity` (see below): the U(1) gauge constraint on Kron's 1944 field network (KCL as Gauss's law), then SU(3) link variables on the same graph | bounded — step 1 passes: on a two-mesh RLC in Kron's mesh frame KVL is zero holonomy of the branch voltage (2e-16), the voltage is pure gauge, KCL is the rate of the Gauss generator Div Q − ρ (1e-17), the mesh EMF is −jω × the flux holonomy, and the Hamiltonian's terms are the ledger's ½Cv² and ½Li²; the U(1) field on a Yee grid equals `em/fdtd.go` edge by edge to 6.6e-13 (TM110, 4000 steps; ω = Yee's to 1e-7, energy to 9e-15); step 2 null: SU(3) loops on the same graph compose with defect ‖[W_p, W_q]‖/3 = 0.75 (min 0.37 over 200 fields), 0 for the abelian (Cartan) part, and every quantity the simulator carries on a branch — voltage, current, charge, flux, power, in the linear circuit and in a diode pair — composes with defect ≤ 3e-16, the round-off of addition; standard prediction: circuit theory is the classical, abelian (U(1)) limit of a lattice gauge theory and the quark model is the non-abelian (SU(3)) one on the same kind of graph; the quantum level of V, I and R is set by electrons and photons through e and h (K_J = 2e/h = 483 597.848 4 GHz/V, R_K = h/e² = 25 812.807 459 Ω, I = e·f), with no quark content; conduction current is electrons only, quarks fix only the lattice charge Z; the e/3 charges seen in the fractional quantum Hall effect are collective electron excitations |
| H23 | The electron as a photon with toroidal topology: one Compton-wavelength photon circulating twice round a torus (Williamson & van der Mark, *Ann. Fond. Louis de Broglie* 22, 133, 1997) | `examples/photon_electron`: the paper's charge, coupling and g against measurement; plain and crossed two-wire ring lines as the circuit form of its 720° double loop (AC modes twin the ladder dispersion, pulse round trips, ledger) | bounded — model charge 0.910 e (eq. 5), α′ = (q/e)²α, anomaly α′/2π = 0.000962 against measured 0.00115965218059 (Fan et al., *PRL* 130, 071801, 2023), 17 % low; spin ½ from the double loop, reproduced as the crossed ring's half-integer modes and its sign flip per trip; the self-confinement is assumed, not derived: tested with `em/fdtd.go` — one wavelength set on the (2, 1) double loop in empty space keeps 5.8 % of its energy within 2R after T/4 and 1e-6 after T/2, all of it leaving through the walls (ledger residual 2e-13); so under linear Maxwell the configuration is not self-confined and the model needs a confining mechanism it does not give |
| H24 | Hydrogen's ground state as a zero-point-fluctuation-determined state: in stochastic electrodynamics the orbiting electron absorbs from the random zero-point field exactly what it radiates when m·ω·r² = ħ (Puthoff, *Phys. Rev. D* 35, 3266, 1987) | `examples/zpf_hydrogen`: the paper's balance with CODATA values (gives a0 and −13.6057 eV, i.e. Bohr's condition); the averaged balance is restoring about a0 (τ = 9.3×10⁻¹¹ s); the paper's 1-D oscillator in a simulated zero-point field holds ½ħω0 with absorbed = radiated (0.03 %); kronspice LC tank in its resistor's zero-point (Callen–Welton) noise holds the expected energy, ledger closes | reproduced at Bohr-theory level, as the paper claims; beyond it untested here — the QM ground state has L = 0 and a spread of radii; full SED simulations: Cole & Zou 2003 report the Schrödinger density, Nieuwenhuizen & Liska 2015 find self-ionization at long times; zero net energy exchange (not an energy source) |
| H25 | Electrical energy from the vacuum by the Casimir cohesion of charged foliated conductors, the "vacuum-fluctuation battery" (Forward, *Phys. Rev. B* 30, 1700, 1984) | `examples/forward_casimir`: Lifshitz energy for real Al (λp = 107 nm: 0.23 of the ideal at 20 nm); one leaf pair as a Kron (q, a) system under Forward's servo, ledger closed, kronspice twin of the electrical side; reopening the leaves | one-shot release of surface (Casimir) energy, consistent with energy conservation (Cole & Puthoff 1993); 12 µJ/m² per Al stroke, 0.02 J/kg; a close–reopen cycle is a net loss (computed) — not a continuous source |
| H26 | Long-lived plasmoids (ball lightning) as force-free time-harmonic electromagnetic fields in a plasma, stabilized by vortical flow and non-radiating at a resonant size (Nachamkin, PL-TR-92-3044, 1992) | `examples/plasmoid_rings`: his l = 1 Beltrami field (field lines on nested tori around a current ring); his electron density = cold-plasma dispersion ω² = ωp² + α²c²; Appendix B recomputed; exact leaky l = 1 resonances of the plasma sphere; zero-frequency limit as a spheromak ring decaying in kronspice; collapsing electron ring driven by a kronspice bank | field mathematics sound; his example needs electron quiver 2.6c (non-relativistic model broken); resonant size does not stop radiation (energy lifetime ~1 µs); TE and TM halves need different sizes; a 1 MJ, 10 cm spheromak lives 19 µs at 10 eV and needs 14 keV for 1 s — not ball lightning in air |
| H27 | Ken Shoulders' "EVs": micrometre beads of about 10¹¹ electrons at near-solid density, few ions, moving at c/10, held together by zero-point fluctuations (Shoulders, *EV: A Tale of Discovery*, 1987; patents 5,018,180 et seq.) | `examples/shoulders_ev`: bead electrostatics and relativistic Coulomb explosion; Casimir, Boyer sphere, magnetic pinch and Brillouin trap as containers; Gauss's law against his kilovolt probe readings; charge (C·V) and energy budgets of his 1 pF and 40 pF sources; kronspice source discharge and EV strike on a 1 pF anode into 50 Ω with ledgers | standard prediction: a 1 µm, 10¹¹-electron bead needs 2.77 J (3×10⁶ × his 1 µJ source), pushes out at 1.5×10¹⁸ Pa and doubles in radius in 1.7×10⁻¹⁵ s; Casimir pressure is 10²¹ too weak and pushes outward on a sphere; kV escape energies need a ≥14 cm "container" or ≤3×10⁵ net electrons; charge collected per shot can never exceed C·V of the source. Test: Faraday cup + ledger on one shot. Mainstream counterpart: Mesyats ectons |

**H6 — the Sun and high-voltage pulse fields (theoretical).** The working
idea is that the Sun's energy, and its formation, involve high-voltage
pulse fields, with the temporal field of H4 (and the field structures of
H1, H3) as the link between electromagnetic pulses and heat and light. The
established account is thermonuclear fusion in the core. To stand alongside
it, H6 has to reproduce what that account already matches:

- the solar luminosity, about `3.83 × 10²⁶ W`;
- the measured solar-neutrino fluxes and energy spectra, chain by chain
  (pp, ⁷Be, ⁸B and CNO; Homestake, Super-Kamiokande, SNO, Borexino);
- the Sun's interior sound-speed profile from helioseismology.

How KronSpice contributes: H6 depends on mechanisms that can be tested at
bench scale first — electric fields without magnetic fields and energy
beyond the source (pulse checks 2 and 7), and temporal-field heating (H4).
A reproduced deviation at bench scale is the prerequisite; if one is
found, a stellar-scale model is added at that point (as a new package),
with the three observational constraints above as its validation cases. Until then H6 is
recorded as *untested*, with its standard counterpart stated beside it.

**H20 — Kron's polyhedral waves as a physical field region.** In his 1962
curve-fitting paper Kron built a k-dimensional polyhedron whose q-networks are
copies of his 1944 Maxwell network, even-dimensional ones magnetic and
odd-dimensional ones dielectric, joined to a dual polyhedron with the roles
exchanged, so that one complete electromagnetic wave spans two neighbouring
dimensions and a k-dimensional polyhedron carries k/2 full waves (pp. 53–54).
He wrote that a field region filled with plasma "crystallizes" into such a
structure energized by k-dimensional generators (p. 46), that the structure
can oscillate self-adaptively, adjusting its modes and frequencies to its
environment (pp. 50, 67), and that a crystal's X-ray-excited dipole waves are
"surprisingly complete analogues" of these polyhedral waves (pp. 67–68). The
paper itself uses the polyhedron only as an estimator, in the space of a data
matrix's k variables, and calls the physics "still puzzling" and the
interpretation "only tentative" (p. 66); the hypothesis is the physical reading
the user has placed in the register. It is tested in three steps, each with
its standard counterpart beside it:

1. **Build the polyhedron** (`kron/polyhedron.go`): p-networks and the
   incidence matrices between them, the even/odd material rule, the dual
   polyhedron and its conductive and inductive ties, and the tearing of the
   polyhedron into its q-networks (pp. 54–55), on the `kron/` primitive and
   connection tensors. Validation: for k = 2 and 3 the structure is a
   discretization of Maxwell's equations and must reproduce the §3g
   validation fields (straight wire, coaxial TEM pulse, rectangular cavity
   modes) to the same tolerance as `em/fdtd.go` and `em/tlm.go`; the 1962
   curve-fitting examples (six functions at four points, p. 46 and the
   tables) are reproduced as the paper's own check. This is the standard
   prediction: nothing beyond Maxwell at k = 3.
2. **Self-adaptive oscillation.** The claim that the modes and frequencies
   adjust themselves to the excitation. Standard prediction: a passive linear
   network has fixed modes, and Kron's own "self-adaptation" in the paper is
   the re-choice of reference frame (tearing and reinterconnection) made by
   the analyst or by a nonlinear program (p. 67), not a field effect. The
   test is whether a built polyhedron, driven as in the paper, shifts its
   resonances with the drive beyond what its nonlinear elements explain; the
   ledger must show where any such adjustment's energy comes from.
3. **Crystals.** Standard counterpart: Ewald–Laue dynamical theory of X-ray
   diffraction, which already describes a crystal's dipole and diffracted
   waves as one self-consistent system with fixed modes set by the lattice.
   A deviation would be a measured mode or frequency adaptation of the
   diffracted field beyond that theory. Bench data, when available, are
   imported and compared as in §3g.

Higher dimensions (k > 3) have no standard physical counterpart; polyhedra
with k > 3 remain legitimate as estimators for the GA (the paper's own use)
and are reported as such, never as fields. Status stays *untested* until step
1 is built; a null result at steps 2 and 3 records a bound, as for every entry.

*What the references in hand supply, and what they do not* (every drawing
and table read from the page images, 2026-11-01; transcriptions, numbers
and the build plan in `doc/kron-figure-catalogue.md`). The 1962 paper gives
everything step 1 needs at k = 3: Fig. 1 (p. 49) is the 1944 Maxwell cell
with every element's impedance written out (dielectric lines
`Z = (h₁/h₂h₃)(Δu¹/Δu²Δu³)(c/εp)`, conductive `c/σ`, magnetic
`(h₂h₃/h₁)(Δu²Δu³/Δu¹)(μp/4c)`, ideal transformers `L = M = ∞`; SI values
and the working of the transformers in Balasubramanian 1968 eqs. 5-3/5-4
and Figs. 10–11), Fig. 2 the primal and dual polyhedra, Fig. 3 the
four-point example with its data matrix and six trial functions, pp. 54–55
the tearing and regression algorithm in eight steps (`β̂ = (XₜX)⁻¹Xₜy`,
`ŷ = Xβ̂`, errors `E′ − I′ = y − ŷ`), and Tables I–VI the six functions and
their divided differences to eight figures, the reproduction target
(Table I transcribed in the catalogue). The q-network rule for q ≥ 2,
first thought missing, is held at second hand: Balasubramanian's Fig. 9
(p. 152) reproduces Kron's 1959 algebraic diagram for space filters (an
orthogonal network per dimension, ε on odd and μ on even dimensions,
incidence matrices M^q_{q+1} as div and curl⁻¹ between them) and
Nicholson 1974 pp. 4–6 give its equations (the orthogonal solution N, the
four incidence relations, the propagation recursion); Kron 1959 (*AIEE
Trans.* 78, 554–561) is still wanted to confirm signs and the 1/4 at
q = 2, 3. Numerical targets beyond the tables: Kron 1945 Fig. 14 (a 4 × 4
TE cavity with its analyzer potentials and residuals) and Balasubramanian
Figs. 15–16 (waveguide matching-impedance error, 1 % at Δz/λ = 5 %, 4 % at
10 %). No source gives the oscillatory reference frame,
which Kron says "has not been reduced as yet to a routine computer procedure"
(p. 67). `examples/kron_polyhedron` infers it from the tables: a frame is a
residual direction of whole numbers (the GA recovers frame 1 from the printed
column), and the oscillatory column is the identity, met exactly by any frame
orthogonal to the 4-dimensional range of the characteristics. The GA does not
reach that on Kron's objective; a closed form does. The second planar frame,
the LS + linear column and the full 3-D cell (step 1's k = 3 validation) are not
yet reproduced, so step 1 is partly done. The physical claims (plasma, crystals) have
no equations in any source and stay hypothesis.

**H21 — Kron's machine equations as a model of the ring electron.** The
hypothesis reads a closed ring of circulating charge as one of Kron's
machines: one winding with the ring's self-inductance L, the circulating mass
as the rotor inertia J, an external field entering as the speed-voltage tensor
G (Epilogue eq. 9, G·i = B), and the ring's voltage and current as the
across and through variables of energy at that scale (the pressure–flow pair
of §3a.1 item 3). The established part is Kron's own: eq. 11 (p. 259) writes
the machine equation with the field tensor F_αβ and notes that with e and R
zero it is the equation of motion of a charged particle in an electromagnetic
field (doc/M9c1.md, "Charge and force"). The hypothesis is the step beyond:
that the electron *is* such a ring, so that its stored energy, angular
momentum, magnetic moment and Zeeman behaviour follow from L, I and J. Three
steps, each with its standard counterpart:

1. **Twin test (standard).** Kelvin's 1899 charged ring, already solved in
   `examples/kelvin_ring` charge by charge with the Lorentz force, is written
   as an N machine element: L from the core-averaged `em.MutualCircles` (as
   in `examples/vortex_ring`), J = w·a², the coil's field as G. Kron's frame
   must reproduce the bench: moments of momentum ∓eΦ/2π on switch-on, the
   fluid at Larmor's rate eB/2w, both ledgers closed, and the Zeeman spectra
   of the full ring. This is eq. 11 made concrete, and must agree to the
   bench's own tolerance before anything else is read from the model.
2. **The electron as the ring (hypothesis).** Fix the ring by two measured
   values (charge e and spin angular momentum ħ/2, or e and the magnetic
   moment) and compute the rest from Kron's equations: stored energy ½LI²
   against m_e c², magnetic moment against μ_B and the measured 1.001 16 μ_B,
   g-factor against 2.002 319 304 36, the Zeeman splitting and its
   polarisation, and the ring radius and speed implied. Standard
   prediction beside each: Dirac + QED values; for the classical ring, g = 1
   when mass circulates with charge, a Compton-scale radius for μ_B at
   v ≤ c, Larmor radiation and no stable state.
3. **Bounds.** The ring size against scattering limits (< 10⁻¹⁸ m) and the
   radiated power against the electron's stability record a bound, as for
   H13, which sits on the same bench (`examples/space_vortex`: g = 15/8
   against 2.0023).

Status: *bounded* (2026-11-01). Step 1 passes and step 2 is computed in
`examples/kron_electron`: no ring with the electron's charge-to-mass ratio
reaches g = 2 (g ≤ 1, Larmor precession eB/2m), and the massless ring that
carries the measured moment at S = ħ/2 needs v = 43.6 c at a = 8.9 fm. The
bench measurement that discriminates it is magnetic resonance
(`examples/kron_nmr`): a ring electron (γ = e/2m) would resonate at
14.0 GHz/T against the measured ESR 28.0 GHz/T, and a ring proton at
7.6 MHz/T against the measured NMR 42.6 MHz/T. Any later result is
reported next to the Dirac value and the classical-ring value, never
alone. Step 3's Zeeman polarisation stays with the 3-D bench of
`examples/kelvin_ring`; the 1-DOF N element cannot precess.

**H22 — Quantumtricity: the quark model as circuit electricity at the quantum
level.** The user's hypothesis, from setting the Ohm's-law wheel (P, V, I, R)
beside the Standard Model wheel: that the two describe one structure, the
quark model being the principles of electricity at the quantum level. As
stated it has no prediction, so it is held in the testable form below, with
the established links and the standard prediction beside it.

Established links, each cited:

- **Circuit theory is a gauge theory on a graph.** Kirchhoff's current law
  is charge conservation, the Noether consequence of U(1) gauge invariance;
  Kron's 1944 field network and Tonti's cell method are Maxwell's U(1)
  theory on a lattice (§3a.1). Lattice QCD is the same construction with
  SU(3) link variables (Wilson 1974; Kogut & Susskind 1975), where the
  "currents" are matrices that do not add linearly. Superconducting
  circuits have been designed to simulate such dynamical gauge fields
  (Marcos et al. 2013): a circuit can emulate the structure, which is not
  the claim that it contains it.
- **The quantum level of V, I and R is electrons and photons.** Since the
  2019 SI the volt follows from the Josephson constant 2e/h (Josephson
  1962), the ohm from the von Klitzing constant h/e² = 25 812.807 Ω (von
  Klitzing, Dorda & Pepper 1980), and the ampere can be realized by
  single-electron pumps, I = e·f (Pekola et al. 2013). One lepton and one
  gauge boson of the Standard Model wheel; no quark.
- **Quarks do not flow in a circuit.** Conduction current in a metal is
  electrons only; confined quarks fix the lattice charge Z and nothing
  else a circuit measures. The strong force and the gluon play no part in
  circuit theory; the weak bosons enter kronspice only as the neutrino
  coupling of `exchange/neutrino.go` (H7).
- **Fractional charge in circuits has a standard explanation.** Charges
  e/3 and e/5 are measured by shot noise in the fractional quantum Hall
  effect (de-Picciotto et al. 1997; Saminadayar et al. 1997); they are
  collective excitations of electrons in a two-dimensional gas, not quarks.

Test, in two steps on `kron/gauge.go` and `examples/quantumtricity`:

1. **Abelian step (standard).** Write the U(1) gauge constraint on Kron's
   1944 network and show it is KCL; show the Wilson loop (holonomy) of a
   mesh is its flux linkage, and that the lattice Hamiltonian's electric
   and magnetic terms are the `½Cv²` and `½Li²` of the ledger. This must
   agree with `em/fdtd.go` and `em/tlm.go` to discretization error before
   anything else is read from the model.
2. **Non-abelian step (the hypothesis).** Put SU(3) link variables on the
   same graph (Kogut–Susskind) and look for a circuit quantity that
   composes as those currents do. Prediction of the hypothesis: one
   exists. Standard prediction: none does, because every Kirchhoff sum
   commutes, so the resemblance ends at "both are gauge theories on a
   network", circuit electricity being the abelian classical limit.

A null result at step 2 records that bound; the entry stays, as for every
hypothesis. The name "quantumtricity" is the user's and is used for the
example directory only, never in output as an established term.

Status: *bounded*, work in progress (2026-11-01). The gauge sector is
built and tested; the matter sector is not: the dictionary matches photon
to gluon, and a quark's circuit counterpart is the electron (a charge field
at the nodes, entering the Gauss law as ρ), not V, I, R or P. V is the
electric field on a link, I the matter current along it, P = V·I the work
term J·E, and R has no gauge-theory counterpart (a gauge theory is
lossless). The node charge field with its link hopping current is the next
construction. `kron/gauge.go` builds the cell complex
of a Kron network (Grad = A, Curl = Cᵀ, Curl∘Grad = 0 is Aᵀ·C = 0), the
U(1) field with its Gauss generator, Hamiltonian and leapfrog, a Yee grid
as Kron's 1944 network, and SU(3) link variables with Wilson loops and the
non-abelian Gauss law. `examples/quantumtricity` runs both steps: step 1
holds to round-off on the circuit and the field (the U(1) field on the
Yee grid equals `em/fdtd.go` edge by edge to 6.6e-13, because Hamilton's
equations on the lattice are Yee's scheme; `em/tlm.go` is not yet
written, so that twin is pending); step 2 finds SU(3) loops composing with
defect 0.75 where every circuit quantity, linear or nonlinear, composes
with defect ≤ 3e-16. What would move the entry: a circuit quantity that
carries an index transported along branches and whose node sums do not
commute. None of the simulator's quantities does.

### Non-goals (for now)

- A full ngspice clone. BSIM3/4, XSPICE, CIDER and digital simulation are out
  of scope until the core is proven.
- A schematic GUI.
- A command-line compatible replacement for QSpice, LTspice or ngspice.
- Environment profiles (stars, planets, radiation belts), discharge-plasma
  physics, bubble dynamics and sonoluminescence. Removed from scope to keep the focus on Kron's
  linear and nonlinear circuit analysis, the genetic algorithm and the EM
  field solver. Where such effects are needed they enter as declared input
  waveforms on element ports (§3.1) or as cited equivalent-circuit models
  built from the simulator's linear and nonlinear elements (Kron's and
  Branin's network analogies, §3a.1). Gravity and fast rotation are back
  in scope (user decision, 2026-10-27) as the `spacetime/` package (§3l).

---

## 2. Repository layout

```
kronspice/
├── CLAUDE.md
├── go.mod                      # module github.com/<owner>/kronspice
├── Dockerfile                  # runtime image: static kronspice + evolve (§7)
├── .github/workflows/release.yml  # v* tag: tests, binaries, ghcr.io images (§7)
├── cmd/
│   ├── kronspice/              # CLI: simulate a netlist, write CSV / raw output
│   ├── evolve/                 # CLI: run a GA experiment from a YAML/JSON config
│   └── cir2kicad/              # CLI: netlist → KiCad 8 schematic + results sheet (MNA and Kron torn solutions)
├── kicad/                      # SPICE netlist → .kicad_sch: drawn wires to one rail per net (≤ 1000 elements) or net labels; Read checks connectivity from geometry
├── report/                     # results sheet: MNA vs Kron torn solution, links (torn branches) with current and voltage, power invariance
├── netlist/                    # SPICE netlist lexer, parser, .subckt expansion, .param
│   └── units.go                # engineering suffixes: f p n u m k meg g t, mil
├── exchange/                   # external energy channels + ports carried by EVERY element (§3.1)
│   ├── channel.go              # channel registry: light, thermal, magnetic, …, neutrino, reserved slots
│   ├── port.go                 # ExternalPort: channel, direction, coupling model, ledger column
│   └── base.go                 # device.Base embeds this: all channels present, default "not modeled"
├── circuit/                    # in-memory circuit model (nodes, elements, models)
│   └── builder.go              # programmatic API used by the GA (no text parsing)
├── device/                     # one file per device family
│   ├── device.go               # Device interface, stamping context
│   ├── passive.go              # R, C, L, mutual inductance K
│   ├── source.go               # V, I; DC, PULSE, SIN, PWL, EXP, AC
│   ├── environment.go          # environment-driven sources: photovoltaic cell, thermoelectric (later)
│   ├── controlled.go           # E, F, G, H
│   ├── switch.go               # S, W
│   ├── arc.go                  # arc / spark-gap models: Cassie, Mayr, combined (§3g)
│   ├── core.go                 # saturable magnetic core with hysteresis (Jiles–Atherton) (§3j)
│   ├── royer.go                # Royer self-oscillating converter stage (§3j)
│   ├── magamp.go               # self-saturating magnetic amplifier on the core model (§3j)
│   ├── diode.go                # D (SPICE3 level-1 diode)
│   ├── behavioral.go           # B: PWL nonlinear elements — N-type i=f(v), S-type v=g(i), controlled sources
│   ├── bjt.go                  # Q (Gummel–Poon)
│   ├── machine.go              # N: windings with L(θ), speed matrix G, mechanical port — Kron's moving inductors
│   ├── bearing.go              # P: ball-bearing (Huber) motor — measured contact drop VF, I² torque aiding rotation, back-EMF or thermal mode; κ, ω0 from measurement fits only
│   └── mosfet.go               # M levels 1–3 (later)
├── solver/
│   ├── matrix.go               # MNA matrix abstraction (dense + sparse backends)
│   ├── sparse.go               # sparse LU wrapper, symbolic reuse between iterations
│   ├── newton.go               # Newton–Raphson with damping and limiting
│   └── convergence.go          # gmin stepping, source stepping, pseudo-transient
├── analysis/
│   ├── op.go                   # .OP
│   ├── dc.go                   # .DC sweep
│   ├── ac.go                   # .AC (complex small-signal)
│   ├── tran.go                 # .TRAN, adaptive timestep, LTE control
│   ├── integrate.go            # backward Euler, trapezoidal, Gear-2
│   ├── measure.go              # .MEAS-style measurements in Go (rise time, gain, BW …)
│   └── power.go                # quaternion + GA power theory (active / reactive / distortion parts)
├── output/
│   ├── result.go               # in-memory vectors returned to callers
│   ├── csv.go
│   └── raw.go                  # SPICE3 raw file writer (binary + ASCII) for viewers
├── kron/                       # Gabriel Kron's tensor network analysis + diakoptics (§3a)
│   ├── primitive.go            # primitive (unconnected) network: branch Z / Y tensors
│   ├── connection.go           # connection tensor C built from the circuit graph
│   ├── transform.go            # covariant/contravariant transformations Z' = Cᵀ Z C
│   ├── mesh.go                 # mesh (loop) formulation
│   ├── node.go                 # node-pair formulation (cross-checks MNA)
│   ├── orthogonal.go           # Kron's orthogonal networks: closed paths + open paths in one frame
│   ├── exterior.go             # open paths to the exterior: each element's external ports (§3a.1)
│   ├── chords.go               # TAN interaction terms ("chords"), affine source terms, radiation edges
│   ├── diakoptics.go           # tearing into subnetworks, solving, interconnecting
│   ├── autotear.go             # automatic cut sets: BFS level partition, ground as common bus, retry with fewer pieces (§3a rule 5)
│   ├── polyhedron.go           # Kron's Maxwell cell in SI, TE field lattice, his 1945 cavity, thesis waveguide (H20)
│   ├── curvefit.go             # Kron 1962 polyhedron: characteristics, frames as residual directions, one-mesh circuit, closed-form oscillatory frame (H20)
│   ├── happ.go                 # Happ 1971: Table 12.1 piecewise algorithm, two-level single-common-bus diakoptics, Ch. XI mesh/nodal with paths ignored
│   ├── jacobian.go             # tearing of every Newton Jacobian, as an analysis solver backend (M9c)
│   ├── nonlinear.go            # nonlinear circuits in Kron's frame: state-dependent primitive device tensors (M9c.1)
│   ├── machine.go              # Kron's machine dynamics: Christoffel symbols, rotating (non-holonomic) frames (M9c.1)
│   ├── gaframe.go              # Kron's network in geometric algebra: branch space, loop/cut-set blades (§3b)
│   ├── gauge.go                # gauge theory on Kron's network (H22): cell complex (Grad, Curl, Div; Aᵀ·C = 0), U(1) field with Gauss generator, Hamiltonian ½Cv² + ½Li², Yee grid = em/fdtd.go; SU(3) link variables, Wilson loops, non-abelian Gauss law, composition defect
│   └── polyphase.go            # quaternion polyphase networks, Clarke/Park as rotors, Kron's rotating frames
├── em/                         # Maxwell's equations in quaternion form (§3c)
│   ├── nabla.go                # Hamilton's ∇ operator on discrete quaternion fields
│   ├── maxwell.go              # Maxwell's quaternion field equations (Treatise form)
│   ├── biquat.go               # single-equation biquaternion / spacetime-algebra form ∇F = J
│   ├── vcalc.go                # plain vector calculus (div, curl, grad, Maxwell residuals): the third twin, M9h
│   ├── frames.go               # parity / orthogonal frames and Lorentz boosts, vector rules vs GA (§3c rule 3)
│   ├── extract.go              # quasi-static field solves → lumped R, L, C, K parameters
│   ├── tline.go                # transmission lines derived from the field equations
│   ├── limits.go               # electrical-size check: is the lumped model valid here?
│   ├── lorentz.go              # multivector Lorentz force + power by grade selection ⟨F Q†⟩ (§3f)
│   ├── duality.go              # dual-symmetric Maxwell with optional magnetic sources (hypothetical mode)
│   ├── biotsavart.go           # Biot–Savart segments and loops, force on a loop (M9n)
│   ├── momentum.go             # momentum ledger: forces with their reactions on exterior sources (M9n)
│   ├── hed.go                  # quaternionic electrodynamics (Giardino 2020), hypothetical mode, free parameter β (§3f)
│   ├── temporal.go             # temporal scalar field T (Jack 2003), hypothetical mode, coupling from Thomson coefficients (§3f)
│   ├── scalar/                 # common interface for tier-2 scalar-field theories; comparisons on shared runs (§3k)
│   ├── fdtd.go                 # local Yee FDTD field solver (§3g): conducting walls, graded matched absorbing layer, materials (εr, σ), edge current sources, exact discrete energy ledger (source work = field + heat + absorbed, residual ~1e-13), div B / Gauss checks; validated: TM110 cavity (Yee ω to 1e-7), Hertzian dipole vs Larmor (+0.5 %), flying doughnut vs exact (3.8 % at 8 cells per q1)
│   ├── doughnut.go             # Hellwarth–Nouchi exact toroidal "flying doughnut" pulse (fields and potentials) and its FDTD validation run
│   ├── tlm.go                  # TLM field solver: Kron-style field network, lumped elements at nodes (§3a.1, §3g)
│   ├── topology.go             # field diagnostics: ∇·B, field-line tracing, toroidal moment, detachment
│   ├── probes.go               # virtual B-dot / D-dot / current probes; import of measured probe data
│   ├── dispersion.go           # rate-dependent materials: Debye and Lorentz polarization equations in time and ε(ω), plate-cell circuit twin; water at 25 °C (ε 78.6 → 5.4, τ 8.3 ps)
│   ├── podolsky.go             # Bopp–Podolsky electrodynamics, fourth-order field equations, hypothetical (H19)
│   ├── multipactor.go          # two-surface vacuum multipactor: Furman SEY fit, measured copper (Baglin 2002), growth per impact (from double_hull)
│   ├── coil.go                 # analytic coil (Wheeler L, turn-to-turn C, bifilar V/2, skin effect) and plate C
│   ├── inductance.go           # Neumann mutual inductance, coaxial circles (elliptic, AGM), Nagaoka solenoid (M9i)
│   ├── helix.go                # Knight helical-line model: free-coil SRF, velocity factor, self-capacitance (§3e)
│   ├── dipole.go               # oscillating point dipole: Jackson 9.18 fields, Larmor power, Poynting by GA vs vector (molecules ↔ light)
│   ├── efficiency.go           # ledger efficiency (DC, transient) and Tellegen's sum, M9g
│   ├── hoffmann.go             # Kron's non-Riemannian form (Hoffmann, Rev. Mod. Phys. 21, 535, 1949): torsion Ω of a non-holonomic frame (eq. 19), new metric (eq. 15), eq. 5′ with the 2Ω term; twin vs frame I to 1e-12, and without Ω the "wrong results" of p. 540
│   └── energy.go               # open-system energy ledger (Poynting / Tellegen balance)
├── fluid/                      # liquid media (§3h); leaf package, gonum only; to be expanded later
│   ├── media.go                # property tables: water (H₂O), heavy water (D₂O), others; vs. temperature
│   ├── mist.go                 # two-phase gap medium: gas plus water droplets, for misted spark gaps (§3j)
│   └── acoustic.go             # acoustic drive fields and pressure at distance
├── spacetime/                  # gravity and fast rotation (§3l), stdlib only: Schwarzschild clocks and redshift, Plebanski index, Shapiro delay, light bending, rotor clocks, Sagnac, burst limit
├── material/                   # cited material constants: Seebeck, resistivity, electrode potentials, spectral lines (leaf); to be expanded later
├── hyper/                      # hypercomplex number systems (§3b); standard library only
│   ├── quat/                   # Hamilton quaternions: products, rotors, matrices, exp/log, slerp (M9d)
│   └── ga/                     # geometric (Clifford) algebra Cl(p,q,r): multivectors, products, rotors
├── evo/                        # genetic algorithm layer
│   ├── genome.go               # genome interfaces + encodings
│   ├── encode_values.go        # fixed topology, evolve component values (E12/E24/E96 snapping)
│   ├── encode_topology.go      # variable topology (graph genome / Koza-style developmental)
│   ├── fitness.go              # fitness functions built from analysis + measure
│   ├── operators.go            # mutation, crossover, selection
│   ├── engine.go               # generational loop, worker pool, elitism, logging
│   └── checkpoint.go           # save/resume populations as JSON
├── internal/
│   └── golden/                 # helpers to run ngspice and compare waveforms
│       ├── ref/                    # independent references (Python), ngspice twins, ngspice_ref.sh, ngspice_all.sh
│       └── docker/                 # tools image (Go, gcc, ngspice) + run.ps1 for VS Code tasks
│   └── export/                 # every example's circuits → export/NAME.cir, .kicad_sch, _results.md/.csv (TestExport; KRONSPICE_EXPORT=1 writes)
├── testdata/
│   ├── netlists/               # *.cir test circuits
│   └── golden/                 # reference outputs (*.csv), committed; golden/ngspice/ = ngspice-computed twins
└── examples/
    ├── rc_lowpass_values/      # evolve R/C values to hit a cutoff
    ├── filter_topology/        # evolve a filter topology from scratch
    ├── amplifier_gain/         # evolve a BJT amplifier for gain + bandwidth
    ├── power_filter/           # evolve a three-phase harmonic filter scored with GA power theory
    ├── pcb_layout/             # evolve trace geometry; L/C/K extracted by em/ from Maxwell's equations
    ├── solar_harvester/        # evolve a solar-cell harvesting circuit from a declared irradiance waveform; fitness = ledger efficiency
    ├── bifilar_coil/           # reproduce Tesla US 512,340 energy comparison; evolve coil geometry (§3e)
    ├── bial_bifilar/           # first experiment: Bi–Al capacitor element (Seebeck / galvanic) ∥ bifilar coil (§3e)
    ├── homopolar/              # Tesla US 406,968 homopolar dynamo from the Lorentz force; evolve geometry (§3f)
    ├── kron_machine/           # induction motor in Kron's frames I and II; parametric inductance vs Mathieu (M9c.1)
    ├── rotating_efield/        # the dual case: B produced by a rotating E field, real vs. magnetic-charge prediction (§3f)
    ├── kelvin_ring/            # Kelvin 1899 charged-ring molecule: induction (eΦ/2π, Larmor) and his Zeeman test, gyrostat vs full Lorentz force
    ├── vortex_ring/            # Kelvin 1867 vortex atoms via his ring-electromagnet analogy: ring speed from inductance, leapfrogging rings
    ├── space_vortex/           # Tewari's space-vortex electron (H13) on the electron test bench: charge, moment, g = 15/8, hydrogen orbit
    ├── kron_electron/          # H21: a ring of circulating charge as a Kron machine (inertia in the metric); Kelvin's ring twin test, frame-I torque, electron rings vs Dirac (g ≤ 1, massless ring at 43.6 c); viz/: field-line picture and an animated GIF of the switch-on, the tilt and its field (stdlib image/gif)
    ├── kron_polyhedron/        # H20 step 1: Kron 1945 cavity and Fig. 15, thesis waveguide, Kron 1962 Tables I–VI; planar frames as integer vectors recovered by the GA through a one-mesh circuit (dissipation = error-square); oscillatory frame in closed form
    ├── kron_nmr/               # magnetic resonance (NMR, ESR) as Kron's frame II(b): Bloch's equations as a 3-winding machine with skew G = Lγ[B]ₓ, rotating frame = Park's frame, T1/T2 as resistors, pickup winding; precession, Rabi, π/2 pulse + FID, CW line vs Bloch; ring (g = 1) vs measured frequencies; GIF
    ├── quantumtricity/         # H22: circuit theory as U(1) gauge theory on Kron's network (two-mesh RLC in the mesh frame: KVL = zero holonomy, KCL = Gauss generator rate, mesh flux in Φ₀, ½Cv² + ½Li²; Yee grid beside em.FDTD to 7e-13), then SU(3) on the same graph (composition defect 0.75 vs 0 abelian; every circuit quantity ≤ 3e-16, linear and diode); K_J, R_K, I = e·f
    ├── space_power_generator/  # Tewari's SPG test (UPSM pp. 159–160) through the energy ledger: the loaded-EMF assumption
    ├── bearing_motor/          # Huber ball-bearing motor: measurement registry, κ/ω0 fits, can bearing torque explain the SPG? (H13, H14)
    ├── homopolar_bearing/      # Tesla two-disk homopolar with load current through its ball bearings: bore vs magnetic area, measured bearing drop and damage limit, least drive torque per load watt (circuit only)
    ├── eckardt_homopolar/      # Eckardt 2012 self-excited twin-disk homopolar (H15): paper's equations vs the machine element with its back torque, ledger
    ├── parametric_lcr/         # Amador & Eckardt 2012 parametric LCR with C(t) (H16): every case reproduced with the pump's work on the ledger
    ├── ring_towers/            # thought experiment: towers on the Earth–ionosphere ring as an LC ladder; one pair (traveling / standing / cancel) and two pairs at 90° synchronized so their standing waves add at 0°, 90°, 180°, 270°; eight towers at 45° (pairs A–B, C–D, E–F, G–H); DC bias per pair (continuous or sectored ring); viz/ animations; ledger
    ├── double_hull/            # craft double wall as a ring waveguide: gap voltage per section, DC bias between hulls, vacuum multipactor check (Furman copper SEY, Baglin 2002), Kilpatrick limit, GA for a target voltage
    ├── pais_cavity/            # US 10,144,532 B2 (Pais) annular resonant-cavity hull (H18): ring line by Happ two-level diakoptics vs MNA, sectional drive vs energy spread, patent figures against limits
    ├── water_vibration/        # chemistry (doc/chemistry-plan.md): H₂O bonds as a circuit — atoms = mass capacitors, bonds = spring inductors, Wilson's B as E/F transformers; AC peaks = IR fundamentals; fitted to H₂O, predicts D₂O/HDO within 1.3 %; .cir netlists
    ├── machine_spacetime/      # one rotor in three gravitational fields (Earth lab, GP-B orbit, illustrative neutron star): clocks, precession, Sagnac, burst and tidal margins; and the rotor's own curvature and frame dragging
    ├── kron_schrodinger/       # chemistry: Kron 1945 Schrödinger network on the (ρ, z) plane; H atom, H₂⁺ bond (R_e 1.998 a₀, D_e 2.786 eV vs exact 1.997, 2.79); LC-circuit twin in kronspice AC
    ├── dna_vibration/          # chemistry: the four DNA bases (methyl-capped) and the A·T, G·C Watson–Crick pairs as circuits (AMBER OL15); AC absorption, Einstein/Larmor release, pairing modes; model vs measured; stacks/ of 2–4 stacked pairs (stacking as its own coupling group G_ST, poly(dA) runs, ApT); sound/ (×3 Hz per cm⁻¹)
    ├── pulsed_bifilar/         # §3e step 3 with pulses: coils as 20 coupled sections (Maxwell coaxial-circle mutuals, K elements) driven through a 5-turn primary by one-sign pulses; ordinary and Tesla bifilar solenoid, conical and flat-spiral bifilar (flat-spiral L within 1.7 % of Wheeler); build-up only when a pulse-train harmonic lands on f1 (ordinary 1.17 MHz: 100 V → 7.4 kV; bifilar 35–42 kHz, V/2 between neighbouring turns); secondary.go: tuned low-voltage pickup at the coil's middle (C = 1/(ω1²L), 3.8 A at 0.77 V, 66–91 % of the power past the source); ledger; viz/ pictures and GIFs
    ├── rate_of_change/         # where higher time derivatives enter: water's Debye response (slow field ε 78.6, 1 ps ramp ε 9.7), time = frequency = kronspice AC circuit twin (7.9e-13); Bopp–Podolsky (H19) beside Maxwell
    ├── coiled_tether/          # the 5-mile tether coiled (ordinary and Tesla bifilar) vs straight: motional EMF along the real wire path (coil = its axial length; out-and-back = 0), solar-driven resonance (×112 voltage, bifilar turn voltage V/2 holding N²/4 the turn energy), ring-up multipactor check between turns, ledger: every joule from the Sun
    ├── photon_electron/        # Williamson & van der Mark toroidal-photon electron (H23) on the electron test bench; light on the double loop in the field solver (leaves in < T/2); pumps.go: three working electron pumps with ledgers — single-electron pump (Pothier 1992, orthodox master equation: 1.000000000 electrons per cycle at 1 MHz, uphill with the gates paying), multipactor pump in a 1 GHz resonator (copper SEY, Kishek–Lau beam-loading bursts, 9 µA collected), resonant-driven negative Cockcroft–Walton ground pump (10 V → 200 V → −1,195 V); magnetron.go: the multipactor pump driven by an LG 2M214 oven magnetron (2.45 GHz, 900 W) through a circulator into a copper re-entrant cavity, evolved by evo (gap, area ≤ λ/10, coupling β, grid, surface; Kilpatrick and 10 kV limits; every transit resolved by ≥ 40 steps so coarse steps cannot invent impacts): 7.5 µA → 0.66 mA collected for 1.32 kW mains, about 1.5 W into the cavity, the rest to the water load; Part 6: arrays of gaps sharing one magnetron through an ideal divider (N as a gene): 1.95 A from the oven magnetron (2,686 gaps, 0.40 m²) and 216 A from a Stellant CWM-100kW (915 MHz, 100 kW; 786,127 gaps, 311 m²), with current set by total copper gap area (m²) and microwave power; plain vs crossed two-wire ring lines: whole-wave vs half-wave modes, pulse upright every trip vs inverted on odd trips (720°), ledger
    ├── water_lamp/             # the neon flashing-lamp oscillator with water instead of neon + radioactive starter: NE-2 reference (primed / dark), sealed water-vapour glow tube at the measured Paschen minimum (480 V at 0.6 Torr·cm, Škoro 2011; dark = cosmic-ray wait ~11 min, UV-LED starter = 0.17 s ± 5 %), spark gap in pure liquid water (Stygar 2006 E·τ^(1/3) = 0.135, extrapolated and labelled; leakage electrolysis on the chemical channel), air gap dry vs water mist (Liu 2018 fog factor, extrapolated; recovery = H5 parameter); seeded statistical strike delay; each drives the bifilar coil's first mode (twin vs kronspice 8e-4): glow lamps ≈ mV, sparks 44–313 V, < 1 % of flash energy until the primary is tuned to f1; tuned.go: flash circuit tuned to the coil (Lp·Cp = Ls·Cs) with ordinary vs quenched gap, and evo-tuned (0.25 µF, 20 turns, quench at 0.65 of the transfer time: 75 % of each flash into the coil vs 0.6 % untuned; coil voltage then far above enamel ratings, flagged); parallel.go: parallel gaps on one capacitor (only one arc; more re-strikes) vs parallel capacitor plates (more energy and current per flash, retuned primary); ledgers to 1e-5
    ├── sparkgap_diode/         # a water gap as a diode from its electrode materials: valve metal (aluminium, anodic oxide) vs ordinary metal (lead, iron) = the Pollak/Nodon electrolytic rectifier (EB 1911: ~140 V per cell, 75–80 %); I–V of pairs (Al|Pb rectifies 13,700:1, Fe|Cu symmetric, Al|Al blocks), oxide sparking past 140 V; half-wave and bridge in kronspice (bridge 85 %, film capacitance fails at 20 kHz); polarity point–plane spark gaps in water (Stygar 2006, extrapolated; negative/positive ratio 2 declared): half-wave works, a bridge needs a ratio above 2; ledgers; pulse test (10 J: water gap 4.2 kA in 0.11 µs, cells 470 A in 9 µs) and bridge frequency sweep as amps per volt vs frequency (cells steady to 3 kHz then the oxide film leaks; water gaps only to 100 Hz: breakdown rises with f^(1/3) and hot gaps, recovering above the 480 V vapour Paschen minimum, re-strike backwards)
    ├── sparkgap_diodebridge/   # water cells + spark gaps combined for six jobs: mains DC supply (twin vs kronspice 0.01 %), surge protection (0.2 mm air-gap arrester + Al|Al two-way clamp: bridge cells 117 V vs 221 V unprotected under a 4 kV 1.2/50 µs surge), spark pulser (water-gap bridge + air gap fired at mains zero crossing: 120 pulses/s, 505 A, 50 % — self-breaking gaps fail by back-strike or power follow), Greinacher doubler (2 cells per position: 174 V, 85 %), oxide-film varactor with memory (C ∝ 1/V(form), 643–2,405 Hz in a 1 mH tank), frequency doubler; own nodal solver for gaps (sim.go) with exact Tellegen ledger
    ├── tether/                 # electrodynamic tether in LEO: push/pull against a planetary field, momentum ledger (§3f)
    ├── arc_surge/              # nanosecond arc pulses, toroidal field surges, field topology, probe comparison (§3g)
    ├── pulse_duration/         # unidirectional impulse sweep; standard vs premise checks at probes (§3g)
    ├── radiant_receiver/       # Tesla US 685,957: plate, condenser to ground, threshold discharge; input-channel ledger (§3d)
    └── charge_pump/            # ground charge pump, HV storage, staged Royer step-down (§3j)
```

Package dependency direction (never import upward):

```
cmd → evo → analysis → solver → device → circuit → netlist
          ↘  kron  ↗        ↘ output
              ↓
            hyper  (quat, ga — leaf packages, import nothing from kronspice)
              ↑
             em    (imports hyper + gonum; produces device parameters and checks)

  fluid, material  (leaf property packages, gonum / stdlib only)

  spacetime        (leaf package, stdlib only: gravity and rotating frames, §3l;
                    may be imported by em, kron, analysis, evo and examples)

  exchange         (leaf package, stdlib only: external energy channels and ports, §3.1)
   ↑ imported by device, em, analysis, evo — never imports anything from kronspice
```

External inputs (irradiance, ambient temperature, external field, neutrino
flux) are declared waveforms or parameter values with cited sources. `evo`
and `cmd` wire them into a circuit through `circuit.Builder`, so devices
never know where their inputs come from.

`em` never imports `analysis` or `solver`. It hands lumped parameters to
    ├── hoffmann_motor/         # Hoffmann 1949 single-phase squirrel-cage motor: frame I (rotor coils turn, L(θ), eq. 42) vs stationary coils 2̄, 3̄ (constant L̄, eq. 29, torsion Ω eq. 19); kronspice frame I, frame II(b) and eq. 5′ integrated directly agree (376.587 rad/s); no start from rest, runs either way when pushed; without Ω only friction coast-down; torque–speed curve, salient pole L_D; ledgers
    ├── zpf_hydrogen/           # Puthoff 1987 (H24): hydrogen ground state as a zero-point-fluctuation balance — eqs. 18–20 with CODATA (a0, −13.6057 eV), restoring averaged balance, 1-D oscillator in a simulated zero-point field (½ħω0, absorbed = radiated), kronspice LC tank in Callen–Welton resistor noise with ledger; later SED simulations and the standard answer stated
    ├── kaluza_kron/            # Kaluza 1921 five-dimensional unification: eq. 3 checked from his metric; exact 5-D geodesic = Lorentz force with q/m = 2p_w/√(1+p_w²) (charge = fifth-dimension momentum, capped at 2√G per mass; electron needs 2×10²¹); Klein circle R = 2l_P/√α = 3.8×10⁻³⁴ m; the same construction in a Kron machine (cyclic winding charges, conserved flux linkages, Routh-reduced rotor) — geodesic, reduced and kronspice agree
    ├── forward_casimir/        # Forward 1984 vacuum-fluctuation battery (H25): Casimir energy ideal and Lifshitz for real Al; one charged leaf pair as a Kron (q, a) electromechanical system under Forward's servo with ledger; kronspice varying-capacitor twin; close–reopen cycle is a net loss; as an electron pump: a one-shot Casimir-driven variable-capacitor pump
    ├── einstein_teleparallel/  # Einstein 1928 distant parallelism vs Kron: the motor's stationary-coil frame is a 4-bein; Einstein's torsion = frame anholonomy = Hoffmann's Ω lowered by L̄; φ (his EM potential) = 0 for Kron's frame, frame-dependent; Einstein's connection flat, the machine's curved; his eqs. 5–9 checked on plane waves; kronspice motor built from the torsion alone equals frame II(b)
    ├── plasmoid_rings/         # Nachamkin 1992 force-free plasmoids (H26) as toroidal current rings: Beltrami field on nested tori, cold-plasma dispersion, Appendix B (quiver 2.6c), exact leaky resonances (µs lifetimes), spheromak ring decay in kronspice (Spitzer), collapsing electron ring (R ∝ B^-½) driven by a kronspice bank, all with ledgers
    ├── shoulders_ev/           # Shoulders 1987 EVs (H27): 10¹¹-electron 1 µm bead tested — self-energy 2.77 J, fs Coulomb explosion, Casimir/Boyer/magnetic/Brillouin containers fail, Gauss's law vs his kV readings, C·V charge and energy budgets of his sources; kronspice source discharge and EV strike into a 50 Ω scope with ledgers; Faraday-cup test
`circuit`/`device` (via the builder) and receives simulation results for its
energy and validity checks through `output.Result`.

`analysis/power.go` (quaternion and GA power measurements) may also import
`hyper`. `solver` and `device` never do: the inner loop stays real and complex.

`kron` sits beside `analysis`: it reads `circuit`, uses `solver` for linear
algebra, and returns the same `output.Result` types so its answers can be
compared directly with MNA.

`device` must not import `analysis`. `evo` talks to the simulator only through
`circuit.Builder` and the `analysis` entry points.

---

## 3. Core numerical design

### Modified Nodal Analysis (MNA)

- Unknowns are node voltages (ground = node 0, eliminated) and branch currents
  for voltage sources, inductors and current-controlled elements.
- Every device **stamps** its contributions into the matrix `A` and the
  right-hand side `b` through a `StampContext`. Devices never touch matrix
  internals directly.
- Matrix positions are reserved once during setup (`Setup` phase) so the
  sparsity pattern is fixed and sparse LU symbolic analysis is reused for every
  Newton iteration and timestep.

### Device interface (sketch)

```go
type Device interface {
    Name() string
    Setup(ctx *SetupContext) error          // allocate nodes, branch eqs, matrix slots
    Load(ctx *LoadContext) error            // stamp for current iterate (DC / transient)
    LoadAC(ctx *ACContext) error            // complex small-signal stamp
    Converged(ctx *LoadContext) bool        // device-level convergence check
    Truncate(ctx *TranContext) float64      // suggested max timestep from LTE, +Inf if none
}
```

Nonlinear devices implement SPICE3-style voltage limiting (`pnjlim`, `fetlim`)
in `Load`. Port the *equations*, not the C structure.

As implemented (M2): every device has `Setup`, `Load`, `Operating`,
`Accept` (commit an accepted solution) and `ExteriorPower`; energy-storing
devices also implement `device.Reactive` (`StoredEnergy`, `Truncate`) and
sources implement `device.Breakpointer`. `LoadAC` and `Converged` arrive
with M3 and M4.

Every device also satisfies `exchange.Exchanger` (by embedding `device.Base`),
so every element carries external energy ports — see §3.1.

---

## 3.1 External energy ports on every element (`exchange/`)

**Vision:** every circuit element — resistor, capacitor, coil, transistor,
arc, solar cell, wire — is an open system that can take energy in from, or
give energy out to, the world outside the circuit. The interface therefore
gives *every* element a full set of **external energy ports**, one per
channel, as potential input and output variables — even for channels whose
coupling has not been defined yet. When a coupling model is added later, it
plugs into the existing port; no interface change, no device rewrite.

### Channels

`exchange/channel.go` defines a registry of channels. Each element has one
port per channel:

| Channel | Examples of in / out | Status at start |
|---|---|---|
| Electrical (terminals) | The normal circuit connection | Modeled (MNA) |
| Light / radiation | Photocurrent in; light emitted by LEDs and arcs | Modeled where defined (photovoltaic cell, photoemission §3d), from a declared irradiance |
| Thermal | Heat absorbed from surroundings; self-heating dissipated out | Modeled (TEMP, ledger) |
| Magnetic field | Induced EMF from external B(t); field radiated out | Modeled where geometry exists (`em/`), from a declared B(t) |
| Electric field | Coupling from external E(t); field emitted | Modeled where geometry exists (`em/`), from a declared E(t) |
| Mechanical | Shaft power (homopolar, N machine element), piezo strain, magnetostriction | Modeled where defined (§3f, M9c.1) |
| Acoustic | Sound/shock in and out (piezo) | Modeled where defined (gyrator ports, `fluid/acoustic.go`) |
| Ionizing radiation | Dose and particle strikes in | Not modeled (a declared pulse source can stand in) |
| Chemical | Battery, fuel cell, electrolysis | Not modeled yet |
| Neutrino | Energy deposited by neutrinos in; none out | Established model (`exchange/neutrino.go`), computed ≈0 (§1) |
| Gravitational / frame | Time dilation, redshift of the element's signals; rotor (rim) clocks, Sagnac | Frame transform (`spacetime/`, §3l), no energy |
| Hypothetical magnetic charge | Dual-pair force/power (§3f) | Hypothetical mode only, separate ledger column |
| Reserved (`Future1…FutureN`) | Channels not yet named | Not modeled |

New channels are added to the registry; every element automatically gains a
port for them.

### Port contract (`exchange/port.go`)

```go
type Direction uint8 // In, Out, Both

type ExternalPort struct {
    Channel   ChannelID
    Direction Direction
    Model     CouplingModel // nil ⇒ "not modeled"
    Status    PortStatus    // Modeled | NotModeled | Hypothetical
}

type CouplingModel interface {
    // Power entering (+) or leaving (−) the element through this port at time t,
    // plus any stamp contributions (e.g. an induced EMF source) for the solver.
    Exchange(ctx *ExchangeContext) (powerW float64, err error)
    Citation() string // published source for the model
}

type Exchanger interface {
    Ports() []ExternalPort // one per registered channel, always complete
}
```

- **Every element exposes every channel.** `device.Base` creates the full
  port set; devices override only the ports they model.
- **Not modeled means zero and visible.** A port with no model contributes
  exactly 0 W and is listed as `NotModeled` in reports, so it is always
  clear which channels were considered and which were not.
- **Models must be cited.** A `CouplingModel` is accepted only with a
  `Citation()`; hypothetical models carry `Status = Hypothetical` and their
  power goes to the hypothetical ledger column only (§3f).
- **Ledger per element, per channel.** The energy ledger (§3c) is organized
  as element × channel, so a report can show, for example, "R3: thermal out
  2.1 mW; neutrino in 0 (computed); chemical: not modeled".
- **Inputs and outputs are variables.** Each port's power (and any
  channel-specific quantity, e.g. irradiance or temperature) is available in
  `output.Result` as a time series, and can be used as a GA objective or
  constraint.
- **Zero cost when unused.** Ports with no model are skipped in the inner
  loop (a precomputed list of active ports per device), so carrying the full
  channel set does not slow simulation.

### Kron's view of the ports

In Kron's terms (§3a.1), each port is an **open path** from the element to the
exterior junction of its channel. A modeled port is an active open path; an
unmodeled port is an inactive junction pair. The ledger is the power flowing
through these open paths.

### How external inputs reach the ports

An experiment declares its external inputs — an irradiance waveform, an
ambient temperature, an external B(t) or E(t), the neutrino flux — each
with a cited source, and the builder attaches the matching cited coupling
model to each element's port: for example, a photocurrent model on the
light port of a solar cell, a thermal model on every element. Elements
remain unaware of where the energy comes from; they only see their ports.

### Rules

1. Every new device embeds `device.Base` and passes `TestAllChannelsPresent`
   (every registered channel has a port).
2. Every modeled port has unit tests with known values and is included in the
   ledger residual check.
3. Unmodeled ports must report exactly zero; a test enforces it.
4. Adding a channel never requires changing an existing device.

### Nonlinear solve

- Newton–Raphson with `reltol=1e-3`, `abstol=1e-12`, `vntol=1e-6` (SPICE
  defaults), configurable through `Options`.
- Iteration limits: `itl1=100` (DC), `itl4=10` (per timestep). Hitting a limit
  returns `ErrNoConvergence`; it is never an infinite loop.
- DC fallback chain: plain Newton → gmin stepping → source stepping →
  pseudo-transient. Each step is logged at debug level.
- As implemented (M4–M5): all four steps exist (`analysis/newton.go`).
  Pseudo-transient continuation adds a unit capacitance from every node to
  ground and a unit inductance in series with every branch, and integrates
  to steady state; it is the step that can follow an S-type element past
  its turning point (a glow lamp striking), where Newton and stepping stall.
  `Options.NodeSet` (SPICE `.nodeset`) gives the starting guess: a circuit
  with negative-slope elements can have several operating points, and the
  guess selects one. `Options.Context` bounds every analysis.

### Transient

- Integration: trapezoidal by default, Gear-2 as an option, backward Euler for
  the first step and after breakpoints.
- Adaptive timestep from local truncation error, with source breakpoints
  (PULSE edges, PWL corners) honored exactly.
- Hard guard: `minStep` and `maxSteps`. Exceeding either returns
  `ErrTimestepTooSmall` or `ErrTooManySteps`.

### AC

- Linearize at the DC operating point, then solve the complex MNA system per
  frequency point. Use `complex128` throughout.

### Linear solver backends

- `dense`: LU with partial pivoting. Default for circuits under about 30
  unknowns, which covers most evolved circuits and is often faster than
  sparse at that size. Since M1 this is a native, dependency-free
  implementation (`solver/dense.go`); a gonum-backed variant may be added
  behind the same interface and must pass the same tests.
- `sparse`: Gilbert–Peierls sparse LU with partial pivoting (see §5). Since
  M6 this is also native and dependency-free (`solver/sparse.go`):
  minimum-degree ordering, KLU-style diagonal preference, pattern learned
  from the first stamping and replayed in stamp order (`doc/M6.md`).
- The backend is chosen automatically by size and can be forced in `Options`.
  Both backends must pass the same tests.

---

## 3a. Kron tensor analysis and diakoptics (`kron/`)

Gabriel Kron treated a circuit as a geometric object. Start from the
**primitive network** — every branch on its own, unconnected, with its
impedance tensor `Z` (or admittance `Y`) and its sources. Interconnecting the
branches is a **change of reference frame**, described by a **connection
tensor** `C`. The connected network's equations follow from tensor
transformation laws rather than from writing Kirchhoff's laws by hand:

```
Primitive (mesh):   e = Z · i            (impressed e, response i)
Primitive (general, E + e = Z · (i + I)  (impressed e, I; responses i, E;
 orthogonal form)                        Kron 1939 p. 417, eq 19.3 p. 480)
Connection:         i = C · i'           (branch currents from mesh currents)
Transformed:        Z' = Cᵀ · Z · C      (covariant transform of impedance)
                    e' = Cᵀ · e
Solve:              i' = Z'⁻¹ · e'       then  i = C · i'
```

The mesh form is the special case with no impressed currents `I` and no
junction-pair voltages `E` on the axes kept. The general form is the one
§3a.1 needs, because open paths carry `E` and `I`. `C` is usually rectangular
(singular); the transformation laws that do not need `C⁻¹` still hold.
Sources (Kron 1939 = *Tensor Analysis of Networks*): Kron 1939 Ch. IV pp.
101–104 (eqs 4.23–4.37, square `C`) and Ch. V
pp. 119–121 (singular `C`); *Tensors for Circuits* Ch. 6 p. 43 (eqs 6.6–6.7).
`C` can be built from the Kirchhoff junction constraints `B · i = 0` by
solving for dependent currents (*Tensors for Circuits* Ch. 7 pp. 46–50); its
columns span the null space of `B`.

The dual node-pair form uses `Y' = Aᵀ · Y · A` with an incidence-based
connection tensor `A` (in Kron's notation `A = Cₜ⁻¹`, so this is his
`Y' = C⁻¹ · Y · Cₜ⁻¹`, Kron 1939 Ch. XIV p. 360; *Tensors for Circuits* eq 6.8).
For linear circuits it must agree with MNA, which makes it a built-in
cross-check.

### Diakoptics ("tearing")

Kron's diakoptics splits a large network into subnetworks, solves each one
independently, then joins the partial solutions through a small
interconnection system. The cost is dominated by the subnetworks, which can be
solved in parallel and **cached**.

**Algorithm (Aitchison 1983, Sec. 2).** For a system `A z = c` whose
matrix is block-diagonal apart from a few coupling entries:

1. Reorder `A` into diagonal blocks plus sparse off-block entries
   (partitioning algorithms: Steward 1965).
2. Build an enlarged, purely block-diagonal matrix `B` (one block per torn
   piece, plus blocks holding the moved coupling entries) with
   `A = P′ B Q′`, where `P′`, `Q′` are 0/1 selection matrices.
3. From them form `P`, `Q` (`P′P = I`, `QQ′ = I`) and `K`, `L` (rank
   `n − m`, built directly from `P′`, `Q′` by the paper's rules).
4. Solve `z = Q (I − B⁻¹K (L B⁻¹K)⁻¹ L) B⁻¹ P c`.

`B⁻¹` is applied block by block (the torn pieces, in parallel); the only
coupled solve is the interconnection system `(L B⁻¹K) y = L B⁻¹ P c`, of
order `n − m` = the number of torn coupling entries. `L B⁻¹K` is
non-singular if and only if `A` is (Aitchison Lemma 1), so tearing never
loses solvability; but the algorithm also needs every block of `B` to be
non-singular (Lemma 2). A torn piece that floats (no reference node in a
nodal formulation) or contains a voltage-source loop makes its block
singular, so the tearing step must give each piece its own reference or
choose cuts that avoid this, and report the reason when it cannot. Cost
estimate (Aitchison): about `N·t + r·t′` against `m·t` untorn, with `N`
the largest block size, `t` the non-zeros, `r` and `t′` the size and
non-zeros of the interconnection matrix.

**Conditions and invariants (Roth 1959).** In the algebraic-topology
formulation, "ohmicness" of the branch operator is necessary and
sufficient for existence and uniqueness of the generalized network problem
met in tearing (Sec. 2); power is invariant under tearing (Sec. 5); and any
linear system can be written as a network problem, though not always in a
form suited to efficient tearing (Sec. 3). Tests check the power
invariance directly: total power of the torn-and-interconnected solution
equals that of the untorn one.

**Caching for the GA.** A mutation inside one torn piece invalidates only
that block's factorization and the rows and columns of `L B⁻¹K` that
involve it; the other blocks' factors are reused from the parent, and the
interconnection matrix is updated rather than rebuilt.

### Why this matters for the genetic algorithm

- A mutation usually changes one part of a circuit. With tearing, only the
  torn piece containing that change is re-solved; the other pieces' factored
  solutions are reused from the parent. This is the main reason `kron` exists
  in this project, and the `evo` engine should use it when it pays off.
- Torn subnetworks map naturally onto goroutines.
- The connection tensor gives a clean, algebraic genome for topology
  evolution: mutating `C` rewires the circuit while the primitive branch
  library stays fixed.

### Rules for `kron/`

1. **Represent tensors as matrices with explicit index roles.** Kron's network
   tensors are rank-2 objects with covariant/contravariant indices; implement
   them on `gonum/mat` (`*mat.Dense`, `*mat.CDense` for AC), and record each
   tensor's frame (primitive, mesh, node, torn piece) in its Go type so a
   transform between the wrong frames fails to compile or returns an error.
2. **Never form `Z'⁻¹` explicitly.** Factor and solve, as in `solver`.
3. **Start linear.** Implement DC and AC for linear circuits first. For
   nonlinear circuits, apply tearing to the linearized Jacobian at each Newton
   iteration (and each transient step), reusing torn factorizations for pieces
   whose devices did not change operating region beyond tolerance.
4. **Every `kron` result has a twin MNA test.** The same circuit solved by
   `kron` (mesh, node and torn) and by `analysis` must agree to `1e-9`
   relative for linear DC and AC.
5. **Tearing strategy is pluggable**: manual cut sets from the netlist (a
   `.tear` directive or builder API), plus an automatic graph partitioner
   added later.
6. Cite Kron's text (book, chapter, equation) in comments when implementing
   one of his constructions.

### 3a.1 Open paths: Kron's model of external inputs and outputs

**What the literature says.**

- **Closed and open paths.** In Kron's orthogonal networks the branch space
  is described in one reference frame made of **closed paths** (meshes,
  carrying mesh currents) and **open paths** (junction pairs, carrying the
  voltages between junctions). For a network with `b` branches, `N` junctions
  and `S` separate sub-networks there are `p = N − S` open paths and
  `m = b − N + S` closed paths, and together they span the whole branch
  space (Kron 1939 p. 75, eqs 3.1–3.2; *Tensors for Circuits* p. 45,
  eqs 6.12–6.13). A connected network has `S = 1`; a torn network has one
  sub-network per piece. Kron 1939 develops this across the primitive mesh and
  junction networks
  (Ch. III), impressed and response quantities (Ch. III p. 76),
  potentials across open-circuited coils and interconnection of several
  networks (Ch. V pp. 128–129, 133), opening of circuits as transformations
  (Ch. VI p. 162), impressed voltages in eliminated axes (Ch. X pp. 248,
  254), junction networks (Ch. XIV), orthogonal networks with "active"
  junction pairs (Ch. XVI), inactive junction pairs and compound networks
  (Ch. XIX p. 485) and input, output and inactive junction pairs (Ch. XXII
  pp. 566–567).
- **Open paths carry the exchange with the outside.** Impressed voltages and
  currents — energy supplied from outside the network — act across open
  paths. Power is invariant under every Kron transformation; in the orthogonal
  frame Kron keeps two invariants, the power input `e·i` and the power
  output `E·I` (Kron 1939 p. 409, pp. 412–413). The power through the open
  paths —
  the `E·I` term — is exactly the power exchanged with the outside. An
  **inactive** junction pair is an open path through which no current `I`
  flows (Kron 1939 p. 485): present in the frame, carrying zero power.
- **Tearing creates open paths.** In diakoptics, removing branches opens
  paths between subnetworks; interconnection reconnects them through those
  open paths (Kron, *Diakoptics*; Roth's algebraic-topology account). The
  idea is already in Kron 1939: a system is divided by opening branches, the
  parts are analysed separately and recombined by a `C` (p. 133), and
  networks are subdivided physically into component networks (Ch. XIX
  p. 492).
- **Modern extension: Tensorial Analysis of Networks (TAN).** Maurice,
  Reineix, Durand, Dubois and co-workers extend Kron's method for
  electromagnetic compatibility: currents as 1-chains and potentials as
  0-cochains, with the electromotive force a cochain of degree 2 on a cell
  complex (arXiv:1412.0187, §5 and after); **chords** — the off-diagonal
  terms of the impedance tensor "replaced by functions of any complexity",
  including radiated coupling from free fields (Durand et al., WCE 2013);
  **affine source terms** — a delayed source emf `E·e^(−x p / v)` added on
  the receiving mesh (WCE 2013); and **radiation-resistance edges** for
  energy radiated away. The radiation resistance is intrinsic to the
  radiating edge but is modified by the environment: inside a shielded room
  the back-induction from the walls changes it, up to cancelling it
  (arXiv:1412.0187). Validated against FDTD (cavity with an aperture,
  shielding effectiveness; arXiv:1412.0187 §13), against measurement for a
  10 GHz horn 50 cm from a metal wall or absorbers (1.2% difference; §14),
  and against measured S21 between monopoles in a cavity, 100 MHz–2 GHz
  (Leman et al., C. R. Physique 2009).
- **Kron's field networks.** Kron also built equivalent-circuit networks of
  Maxwell's field equations (*Proc. IRE*, 1944) and of the Schrödinger
  equation (*Phys. Rev.*, 1945): a region of space becomes a network whose
  boundary is where the outside acts on it. Two worked examples are in
  `reference/` (website transcriptions; text and figures complete, displayed
  equations missing, so take formulas from the journal versions):
  - *J. Appl. Phys.* 16, 172–186 (1945), section "The Field Equations of
    Maxwell": a square waveguide or cavity resonator divided into a grid of
    cells, each with two inductors and a capacitor to ground. Inductor
    currents represent Ex and Ez, capacitor voltages represent Hy. Boundary
    conditions are known currents in some coils and known potentials at
    some junctions; at the open end a sinusoidal current distribution is
    impressed. The same paper solves such networks by iterating on the
    unbalanced junction currents and voltages, and refines the mesh only
    where needed (around corners).
  - *Phys. Rev.* 67, 39–43 (1945): the Schrödinger equation as a ladder of
    series inductors (kinetic energy), shunt inductors (potential energy) and
    shunt capacitors (total energy); ψ is the junction voltage to ground.
    Eigenvalues are found by driving the network from a generator and
    varying E until the generator current is zero. Boundaries: terminate by
    an equivalent impedance where ψ attenuates, otherwise by a short circuit
    at a point of zero voltage or an open circuit at a point of zero current.
  In Kron 1939 terms these boundary quantities are impressed currents and
  voltages on junction pairs, i.e. open paths, but the 1945 papers do not use
  that phrase. *Tensors for Circuits* Ch. 14 gives only Minkowski's tensor
  form of Maxwell's equations, not an equivalent circuit, and the full 3-D
  network is in the 1944 paper (not yet in `reference/`).

**Verified in Kron's text** (*Tensor Analysis of Networks* (Kron 1939), Internet
Archive scan `dli.ernet.8037`, and *Tensors for Circuits*, Dover 1959, both
in `reference/`; printed page numbers; full check in
`doc/kron-3a-verification.md`):

- **p. 75, eqs. 3.1–3.2:** junction-pairs = junctions − sub-networks;
  coils = meshes + junction-pairs.
- **p. 76:** an **open mesh** is "any circuit through the coils that starts
  at one junction of a junction-pair and ends at the other"; it may replace
  the junction-pair in all analysis. Superimposed quantities are of two
  kinds: **impressed quantities that originate outside the network**, and
  **response quantities** (the network's reaction). This is the exterior
  exchange that every port represents.
- **p. 133 (Ch. V):** separate networks are interconnected by a `C`, and the
  apparatus may be electrical, mechanical, acoustical, optical and so on.
  This is the basis for one exterior per channel and for tearing at the
  exterior.
- **pp. 407–408 (Ch. XVI):** impressed currents are also described as
  currents withdrawn by outside loads. A branch with no impedance carrying a
  known impressed voltage (or a coil with no admittance carrying a known
  impressed current) is an **apparent coil**. Kron then removes apparent
  coils from the analysis (p. 408, p. 412): their effect is known components
  of E′ and i′ (p. 417). So an element's external port is a junction-pair
  axis, not an extra coil in the primitive `Z`.
- **p. 408:** any mesh network is an orthogonal network in which currents of
  **zero value** are impressed across its junction-pairs, "for the sake of
  completeness". This is exactly an unmodeled port: present in the frame,
  zero impressed quantity, zero power.
- **p. 411:** a closed mesh can be opened, and an open mesh closed, by an
  apparent coil; an axis is a junction-pair when its E is unknown and a mesh
  when its i is unknown.
- **p. 429:** in setting up z, **inactive junction-pairs may be ignored**
  (and inactive meshes in setting up Y), and the singular transformation of
  a mesh or junction network is part of a non-singular square one obtained
  by supplying the missing orthogonal axes. This justifies skipping
  unmodeled ports in the inner loop while keeping them in the frame. The
  same page says **no mesh may be ignored** in setting up z, so a port
  modelled as a closed mesh cannot be skipped this way.
- **p. 485 (Ch. XIX):** junction pairs through which no current `I` flows are
  **inactive**; in setting up z they are ignored, while inactive meshes need
  the Ch. X reduction formulas.
- **pp. 566–567 (Ch. XXII):** junction pairs are classed as **input, output
  and inactive**, and the inactive ones are eliminated so the analysis runs
  on active axes only. This is Kron's own version of the port model.

Code implementing these constructions cites the page.

**How the engine applies it (`kron/orthogonal.go`, `kron/exterior.go`,
`kron/chords.go`).**

1. **One exterior junction per channel.** For every channel in the
   `exchange/` registry (§3.1) the network gains an **exterior junction**
   `Ω_channel` — Kron's datum for "everything outside the circuit".
2. **Every element port is an open path.** Each element's external port on a
   channel is an open path from that element to `Ω_channel`. It is a
   junction-pair axis, not an extra coil: a driven port enters as known
   components of E′ and i′ (Kron 1939 pp. 408, 412, 417), so no zero-impedance
   rows appear in `Z`. So every element, for every channel, has an open
   path to the outside, exactly as §3.1 requires.
3. **Generalized across/through variables.** Non-electrical channels use
   Kron's network analogies (formalized by Branin): each channel has an
   across variable and a through variable whose product is power — thermal
   (temperature, heat flow), mechanical (velocity, force or angular velocity,
   torque), acoustic (pressure, volume velocity). Channels with only a power
   balance (light, ionizing radiation, neutrino) use a power-only open path.
4. **Unmodeled = inactive junction pair.** An open path with no coupling
   model is an inactive junction pair: present in the frame, open-circuited,
   no current `I`, zero power (Kron 1939 p. 485). Inactive junction pairs may be
   skipped when setting up `Z`; meshes may not (Kron 1939 p. 429). This is the
   Kron-frame meaning of
   "not modeled".
5. **Couplings as impressed sources and chords.**
   - Environment inputs (§3d) are **impressed sources on open paths**
     (affine source terms when delayed).
   - Element-to-element field coupling and coupling from external fields are
     **chords** in the impedance tensor.
   - Energy radiated away is a **radiation-resistance edge** to
     `Ω_electromagnetic`. Its value depends on the environment: reflecting
     surfaces near the radiator return energy as chords (delayed emf) and
     change the effective radiation resistance, as in Maurice's antenna–
     wall case.
6. **The ledger is Kron's power invariant.** The per-element × per-channel
   ledger (§3c, §3.1) is computed in the orthogonal frame as the power
   through each open path to the exterior (Kron's output invariant `E·I`).
   Because power is invariant under the connection transformation, the
   open-path total must equal the
   MNA-computed total: this is an extra twin test of the ledger.
7. **Tearing at the exterior.** Diakoptics can tear along the open paths, so
   the circuit and each channel network are solved separately and joined
   only through the exterior junctions — the natural way to add thermal,
   mechanical or acoustic networks without enlarging the electrical solve.
8. **Field regions as Kron networks.** For `arc_surge` (§3g), the local field
   region can also be expressed as a Kron field network (1944) and compared
   with `em/fdtd.go`, with the region boundary as open paths: impressed
   boundary currents and potentials as in Kron's 1945 cavity example, and
   terminations (matching impedance, short at zero voltage, open at zero
   current) as in his 1945 Schrödinger networks. The 1945 cavity example
   (two inductors and a ground capacitor per cell, Ex/Ez as currents, Hy as
   voltages) is a good first validation case.

   **Implementation route: TLM (`em/tlm.go`).** The transmission-line matrix
   method is the modern numerical descendant of Kron's 1944 field network
   (Russer, *The Transmission Line Matrix Method*, cites Kron 1944 as the
   origin of the field-as-mesh analogy; copy in `reference/`). Space is
   divided into cells; wave pulses travel between cells and scatter at the
   nodes. It complements `em/fdtd.go` and maps onto §3a.1 as follows:
   - **Boundaries:** electric wall ρ = −1, magnetic wall ρ = +1, absorbing
     wall no reflection (Russer §VI). A boundary radiating into free space
     is non-local and needs a discrete Green's function, not a local rule.
   - **Boundary as exchange with the outside:** the discrete Green's
     function relates the waves entering a region's boundary to those
     leaving it (Russer §VII). This is the field-region form of an open
     path to the exterior, and of tearing: a subdomain can be replaced by
     its boundary response and joined to the rest only there.
   - **Circuit elements in the field:** lumped and nonlinear elements attach
     to a node through a half-cell line stub (Russer §X, distributed diode
     example). This is how the arc device of §3g couples into the field.
   - **Skin effect:** lumped R–L ladders at metal surfaces (Russer §IX).
   - **Variables:** TLM works in incident and reflected wave amplitudes, not
     Kron's `e = Z·i`. Where a TLM region meets the `kron/` network, convert
     explicitly between wave amplitudes and voltage/current, and check that
     the power through the interface is the same in both forms.

**Tests.**

- Closed-path and open-path counts satisfy `m = b − N + S` and
  `p = N − S` (S = number of sub-networks), both before and after tearing,
  and the orthogonal frame is non-singular.
- Power through all open paths to all exterior junctions equals the ledger
  computed from MNA, to `1e-9` relative.
- Inactive open paths carry exactly zero power.
- TAN chord coupling for two loops with mutual inductance matches the MNA
  `K` element; a radiation-resistance edge matches the analytic small-loop
  radiated power.
- Torn solutions: for every linear golden circuit, the Aitchison-form torn
  solve equals the untorn solve to `1e-9` relative, total power is the
  same (Roth Sec. 5), and a cut that leaves a singular block is detected
  and reported rather than solved.
- Radiation edge in an environment: an antenna edge facing a metal wall
  shows the delayed reflected emf and the changed input behavior of
  Maurice's 10 GHz horn case, and the same edge facing absorbers shows the
  free-space radiation resistance (arXiv:1412.0187 §14).

---

## 3b. Quaternions and geometric algebra (`hyper/`, `kron/gaframe.go`, `kron/polyphase.go`)

This layer extends Kron's work. Kron already described networks with change of
reference frame (connection tensors) and, for machines, with rotating frames.
Geometric algebra (GA) and quaternions express both ideas directly: a frame
change is an outer-morphism of a vector space, and a rotation is a rotor
`R · v · R̃`. The goal is a richer *formulation and analysis* layer, not a
replacement for the numerical core.

### Ground rule: the solver core stays real/complex

`solver`, `device` and the hot paths of `analysis` keep `float64` and
`complex128`. Complex numbers are exactly the even subalgebra of the 2-D GA,
so nothing is lost. GA and quaternion code runs **before** the solve (building
or checking the formulation) and **after** it (interpreting results, power
measures, fitness). Every GA/quaternion result must reduce to the same
answer as the MNA path — see the twin-test rule below.

### Where each tool is used

**1. Geometric algebra view of Kron's network (`kron/gaframe.go`)**
- Branch currents and voltages live in a branch vector space `B` with one
  basis vector per branch, giving the algebra `Cl(n_branches)`.
- Kron's connection tensor `C` and the incidence matrix `A` become linear maps
  on that space. The Kirchhoff laws are the Roth/Branin
  boundary–coboundary pair: currents are 1-chains, node potentials are
  0-cochains, meshes are cycles.
- The **outer product** of the mesh (or cut-set) vectors is a blade whose
  magnitude is the determinant of Kron's transformed tensor `Z' = Cᵀ Z C`.
  A near-zero blade magnitude flags dependent loops or a singular circuit
  **before** factorization. The GA uses this to reject degenerate evolved
  topologies cheaply and with a precise error.
- Orthogonal decomposition of `B` into the loop space and cut-set space by
  projection/rejection is the GA form of Kron's orthogonal networks.

**2. Quaternions for polyphase networks and Kron's rotating frames (`kron/polyphase.go`)**
- Three-phase quantities `(v_a, v_b, v_c)` map to a pure quaternion
  `v = v_a·i + v_b·j + v_c·k`.
- The Clarke (abc → αβ0) and Park (αβ0 → dq0) transforms are rotations, so
  they are implemented as unit-quaternion rotors `q · v · q*`, with Park's
  angle `θ(t)` making the rotor time-varying. This is the quaternion form of
  Kron's generalized rotating-machine frames.
- **Quaternion power**: `s = v · i*` gives a scalar part (instantaneous active
  power) and a vector part (instantaneous non-active power with a direction),
  following the three-phase quaternion power literature.
- `hyper/quat` does the arithmetic (standard library only, M9d); gonum's
  `num/dualquat` is a reference if rigid motions (rotation + translation)
  are ever needed for coupled electromechanical models.

**3. Geometric algebra power theory (`analysis/power.go`)**
- For non-sinusoidal waveforms, each harmonic (cos and sin component) is a
  basis vector of a GA. Voltage `V` and current `I` become multivectors, and
  the geometric product `M = V · I` splits into a scalar (active power) plus
  bivector parts (reactive and distortion power), each tied to a specific
  pair of harmonics.
- Results are computed from `.TRAN` output via FFT, or from multi-tone AC.
- These measures feed GA fitness functions: power factor, distortion power and
  which harmonic pairs cause it, three-phase imbalance. This makes it possible
  to evolve power-electronics filters and compensators for power quality.

**4. Hypercomplex genome operators (`evo/`, optional)**
- Parameters that are rotations or phase relationships (e.g. phase shifts in
  a polyphase network) are encoded as unit quaternions or rotors. Mutation is
  a small random rotor applied by composition, so mutated genes stay valid
  by construction.

### `hyper/ga` design

- **Basis blades as bitmasks**: blade `e₁∧e₃` is mask `0b101`. Products use
  bit operations plus a sign from reordering and the metric signature
  `(p, q, r)`.
- **Two tiers**:
  - Fixed, allocation-free types for small algebras used in hot-ish code:
    `G2` (complex check), `G3` (3-D vectors, rotors, quaternion isomorphism).
  - A general sparse `Multivector` for `Cl(p,q,r)` of any size, used by
    `gaframe` and power theory (harmonic spaces can have dozens of
    dimensions; store only non-zero blades).
- Operations: geometric, outer, left/right contraction, scalar product,
  reverse, grade involution, dual (where the pseudoscalar is invertible),
  grade projection, norm, inverse, rotor exponential/log.
- `hyper/quat` and `hyper/ga.G3` even subalgebra must interconvert exactly;
  test it.

### Rules for this layer

1. **Twin tests.** Each GA or quaternion formulation is checked against the
   plain matrix/complex result: blade magnitude vs. `det(Z')`; quaternion
   Clarke/Park vs. the standard matrix transforms; GA power vs. IEEE-1459-style
   P, Q and distortion computed conventionally. Tolerance `1e-9` relative.
2. **Property tests** (`testing/quick` or fuzzing) for algebra laws:
   associativity, `R R̃ = 1` for rotors, `|ab| = |a||b|` for quaternions,
   grade rules of the outer product.
3. **No hypercomplex types in `solver` or `device`.** If a GA idea needs a new
   numerical method, implement it as a separate solver backend behind
   `solver.Matrix` and benchmark it against LU before adopting it.
4. **Cite the source** (book, paper, equation) for each construction, as with
   Kron.
5. Keep sign and ordering conventions written down in `hyper/ga/doc.go`
   (basis order, metric signature, reverse/conjugate notation) — most GA bugs
   are convention bugs.

---

## 3c. Maxwell compliance: quaternionic electromagnetism (`em/`)

James Clerk Maxwell's *Treatise on Electricity and Magnetism* (1873) states
the general equations of the electromagnetic field both in components and in
Hamilton's quaternion notation. The whole engine is built to be consistent
with that formulation.

### What "compliance" means here — read this first

- Maxwell's quaternion equations, the modern Gibbs–Heaviside vector equations
  and the spacetime-algebra equation `∇F = J` are **mathematically equivalent
  descriptions of the same physics**. Using the quaternion form does not
  change any prediction; it is the formulation and the language of this code.
  Do not implement anything claiming extra physics beyond the standard
  equations (e.g. "scalar wave" effects); such claims are not supported by
  experiment.
- A circuit simulator is the **lumped, quasi-static limit** of Maxwell's
  equations: components much smaller than a wavelength, fields confined to
  components, no radiation. Compliance therefore has three concrete parts:
  1. **Derivation** — every circuit law and device model is traceable to
     Maxwell's quaternion equations plus a stated approximation.
  2. **Validity** — the engine checks each circuit against the assumptions
     of that approximation and warns or errors when they fail.
  3. **Energy balance for open systems** — energy may flow into and out of
     the circuit, but every joule is accounted for (see "Open systems and
     the energy ledger" below).

### Open systems and the energy ledger

Circuits are **open systems**. Like the Earth receiving sunlight and radiating
heat to space, a circuit takes energy in through sources (batteries, mains,
solar cells, signal generators, antennas) and sends it out through loads,
heat and radiation. Total energy inside the circuit is *not* constant and the
engine must never assume it is.

What Maxwell's equations do require is **local balance**, Poynting's theorem:
inside any boundary, the change in stored energy equals the energy that
crossed the boundary minus the energy dissipated. In circuit form:

```
dW_stored/dt  =  P_in(sources, ports, environment)  −  P_out(ports, radiation)  −  P_dissipated
```

so energy is neither created nor destroyed *inside* the boundary, while the
boundary itself is open. Tellegen's theorem is the instantaneous version: when
sources are counted as branches, the sum of all branch powers is zero, which
just says "power delivered by sources = power absorbed by everything else".

The engine therefore keeps an **energy ledger** for every analysis instead of
a "conservation" flag:

| Ledger entry | Elements | Sign |
|---|---|---|
| Inflow | independent sources, environment-driven sources (photovoltaic irradiance, thermal, RF pickup), external ports driven from outside | + |
| Outflow | external ports delivering to the outside world, radiation estimate when `limits.go` flags large electrical size | − |
| Dissipation | resistors, diode/transistor losses, switch losses | − |
| Storage | capacitors `½Cv²`, inductors `½Li²`, coupled inductors | Δ over time |
| Controlled sources | E/F/G/H are idealized and can supply power; their power is reported separately as "model-supplied" so it is visible, not hidden | ± |
| Negative resistance, nonlinear elements in the 2nd/4th quadrant | `R < 0`, and B elements whose curve makes `v·i < 0`, deliver power; booked as model-supplied, per element, never as inflow from the exterior | + |

The ledger is checked, not enforced: a residual (`inflow − outflow −
dissipation − Δstorage`) larger than the integration tolerance means a
numerical or modeling bug. A circuit that shows net energy appearing from
nowhere inside the boundary is reported as an error, not a discovery.

The ledger is also a **fitness source** for the GA: efficiency
(`useful outflow / inflow`), standby loss, energy harvested from an
environmental input, energy stored at the end of a run.

**Negative resistance, duality and sign inversion.** Two different
transformations relate the passive elements, and the engine keeps them apart:

- *Duality* (exchanging the roles of voltage and current, as a gyrator does)
  maps a capacitor to an inductor and back, and a resistance `R` to a
  conductance `1/R` — still a positive resistor.
- *Sign inversion* (reversing the current's sign, i.e. time reversal, or a
  negative-impedance converter) maps `R` to `−R`: dissipation becomes
  generation, while `C` and `L` map to themselves.

A negative resistance is accepted in netlists and genomes (as in ngspice) and
delivers `v²/|R|`. It is an idealized active element: a real one (tunnel or
Gunn diode, op-amp converter) takes that energy from a bias supply the ideal
element leaves out. So its power is **model-supplied** — listed per element,
balancing the ledger, never counted as inflow — and its thermal port is not
modeled. Fitness that measures harvested energy must count model-supplied
power as an input (`evo.ModelSuppliedLimit` constrains it).

**Nonlinear elements in Kron's frame.** Node-pair voltages are open-path
quantities (§3a.1). A voltage-controlled nonlinearity (`i = f(v)`, N-type:
diode, tunnel diode) is written directly in that frame; a current-controlled
one (`v = g(i)`, S-type: arc, glow lamp) enters by adding its current as an
unknown. The open-path frame therefore covers both, which a pure closed-path
(mesh) formulation cannot. Both forms are available (`device/behavioral.go`;
`circuit.Builder.NonlinearVC` / `NonlinearCC`; ngspice `B … I=pwl(V(…))` and
`B … V=pwl(I(V…))`), with negative-slope regions allowed. Every terminal pair
can be measured at its two extremes, **open-circuit voltage and short-circuit
current** (`analysis.OpenShort`), with the full nonlinear solve; both are GA
objectives.

### Maxwell's quaternion form

Hamilton's operator `∇ = i ∂/∂x + j ∂/∂y + k ∂/∂z` applied to a pure
quaternion (vector) field `v` gives

```
∇v = S.∇v + V.∇v = −div v + curl v      (Hamilton sign convention: i² = j² = k² = −1)
```

so a single quaternion product carries both divergence (scalar part) and curl
(vector part). In this notation the Treatise's main relations include, with
`𝔄` the vector potential ("electromagnetic momentum"), `ψ` the electric
potential, `𝔅` magnetic induction, `𝔈` electromotive intensity, `𝔇`
displacement, `ℌ` magnetic force and `ℭ` total current:

```
𝔅 = V.∇𝔄   (= ∇𝔄, since S.∇𝔄 = 0)   (A) induction from vector potential
𝔈 = V.𝔊𝔅 − ∂𝔄/∂t − ∇ψ               (B) electromotive intensity (𝔊 = velocity)
𝔅 = ℌ + 4π𝔍                          (D) magnetization (𝔍 = intensity of magnetization)
4π ℭ = V.∇ℌ                          (E) currents from magnetic force
𝔎 = C 𝔈                              Ohm's law
𝔇 = (1/4π) K 𝔈                       displacement
ℭ = 𝔎 + ∂𝔇/∂t                         total current: conduction + displacement
𝔅 = μ ℌ                              only when magnetization arises from induction
ℌ = −∇Ω                              where a magnetic potential exists
```

Source: Treatise Vol. II, Art. 618 (symbols) and Art. 619 (equations),
pp. 236–238 of the 1873 first edition in `reference/`; letters (A)–(E) are
Maxwell's. Equation (A) carries Maxwell's own condition `S.∇𝔄 = 0`, under
which `∇𝔄` is a pure vector and the `V.` is unnecessary — the Coulomb
gauge; the code states the gauge it uses. `𝔅 = μℌ` is Maxwell's relation
for induced magnetization only; the general relation is (D), so permanent
magnets and hysteretic cores (`device/core.go`) use (D). Maxwell calls the
scalar part `S.∇σ` the *convergence* of σ (Vol. I, Art. 25), i.e. `−div σ`.
Art. 619 also gives the electric and magnetic volume densities as scalar
parts of `∇𝔇` and `∇𝔍`; their signs relative to `div` must be read from the
page images before they are coded, because the scan's text layer drops
signs and the Fraktur symbols.

(Gaussian-style units as in the Treatise; the code works in SI and documents
the conversion in `em/doc.go`.)

The modern compact equivalent, used for implementation and cross-checks, is
the biquaternion / spacetime-algebra form: with `F = E + I c B` and the
four-gradient `∇`, all four equations become `∇F = J/(ε₀c)`. This is the
same algebra as `hyper/ga` in signature `Cl(1,3)`.

### How each circuit law traces back to Maxwell

| Circuit element of the engine | Derived from (quaternion form) | Approximation |
|---|---|---|
| KCL at every node | scalar part of `∇` applied to `ℭ` vanishes (total current is solenoidal); equivalently charge continuity | Charge stored only in capacitors |
| KVL around every loop | vector part of `∇𝔈` = `−∂𝔅/∂t`; integrated around a loop | Changing flux confined to inductors |
| Resistor `v = R i` | `𝔎 = C 𝔈` (Ohm's law, conductivity) | Uniform conductor, quasi-static |
| Capacitor `i = C dv/dt` | displacement current `∂𝔇/∂t`, `𝔇 = K𝔈/4π` | Field confined between plates |
| Inductor `v = L di/dt`, coupling `K` | `𝔅 = V.∇𝔄`, flux linkage of `𝔄` around the loop | Field confined to the coil |
| Transmission line | full field equations on a 2-D cross-section (TEM) | Uniform line, TEM mode |
| Tellegen's theorem / power balance | Poynting's theorem (energy conservation) | Lumped limit |

Each device file must carry a comment citing the row it implements.

### What `em/` does

1. **Parameter extraction (`extract.go`).** Small quasi-static field solves in
   quaternion form (finite differences on a grid, using `nabla.go`) compute
   R, L, C and mutual inductance K for physical geometries — PCB traces,
   coils, plates. These become device parameters, so evolved circuits can
   include physical layout and parasitic coupling.
2. **Transmission lines (`tline.go`).** Per-unit-length L and C come from a
   2-D field solve; the line is then a device in the circuit.
3. **Validity check (`limits.go`).** Each circuit can declare physical sizes.
   Using the highest frequency in the analysis (AC stop frequency, or the
   bandwidth implied by the fastest transient edge), the engine computes the
   electrical size `d/λ`. Above a configurable threshold (default `λ/10`) it
   returns `ErrNotQuasiStatic` or a warning, telling the user to use a
   transmission-line or field model. The GA penalizes such circuits when the
   experiment requests physical realism.
4. **Energy ledger and charge checks (`energy.go`).** For every analysis,
   build the open-system ledger above: inflow, outflow, dissipation, storage
   and model-supplied power, per element and in total, at each operating
   point and timestep. Check that the residual stays within the integration
   tolerance and that charge is conserved at each node (KCL). Expose the
   ledger in `output.Result` for fitness functions and reports.
5. **Field equations (`maxwell.go`, `biquat.go`).** Both forms are
   implemented on discrete quaternion fields and must agree with each other
   and with a standard vector-calculus implementation.

### Rules for `em/`

1. **Quaternion conventions are fixed and documented** in `em/doc.go`:
   Hamilton's sign convention, scalar part `= −div`, Treatise symbols mapped to
   Go identifiers, and the unit conversion Gaussian ↔ SI.
2. **Triple-twin tests.** Every field computation is checked three ways:
   Treatise quaternion form, biquaternion/spacetime-algebra form, and plain
   vector calculus. Tolerance `1e-9` relative on smooth analytic fields.
3. **Analytic validation cases** must pass before use: parallel-plate and
   coaxial capacitance, coaxial and solenoid inductance, a plane wave
   (checks `∇F = 0` in free space), and the telegrapher's equations for a
   lossless line. Two frame tests are also required:
   - **Parity.** Inverting the coordinate frame (`e₁,e₂,e₃ → −e₁,−e₂,−e₃`)
     flips E and leaves B unchanged. B is an axial vector, carried as the
     bivector `I·B` in `hyper/ga`; in plain vector calculus the curl and
     cross product carry the determinant of the transformation. Both forms
     must give identical fields after the inversion and after any proper
     or improper rotation.
   - **Moving observer.** The field seen by a moving probe is obtained by
     the Lorentz transformation of E and B (in GA, by a rotor applied to
     `F`). The GA result must match the textbook transformation; for a
     probe moving past a static charge it must show the expected magnetic
     field, and for `v ≪ c` it must reduce to `B' ≈ B − v×E/c²`.
   These guard against mistaking a bookkeeping error in frames for a
   difference between the quaternion/GA form and Maxwell: the forms are
   equivalent, so any difference is a bug.
4. **The solver core is unchanged.** `em` supplies parameters and checks; it
   never runs inside Newton iterations or the LU solve. A full-wave field
   simulation of an entire evolved circuit is out of scope. The one
   exception is `em/fdtd.go`: a *local* full-wave solve of a bounded region
   around a pulse source, driven by currents from the circuit simulation
   (one-way coupling first, §3g).
5. **Ledger checks are always on in tests** and can be enabled in GA runs
   (`Options.CheckEnergy`) to catch numerically suspicious candidates. Never
   write a check that assumes total energy inside the circuit is constant;
   always account for inflow and outflow through the boundary.
6. **Environment-driven sources** (e.g. a photovoltaic cell driven by an
   irradiance waveform, a thermoelectric source driven by a temperature
   difference) are first-class inflow devices. Their input from the outside
   world is a declared waveform, recorded in the ledger as inflow.
7. Cite the Treatise (volume, article) or the paper for every equation.

---

## 3d. Radiant-energy receiver (`examples/radiant_receiver/`)

Tesla's U.S. Patent 685,957 (1901, "Apparatus for the Utilization of
Radiant Energy") connects an insulated, preferably polished plate exposed to
sunlight, ultraviolet light or X-rays to one terminal of a condenser, the
other terminal to ground (or to an elevated conductor). Charge accumulates
in the condenser while the plate is exposed, and a threshold switch (two
light plates that attract and touch, or a small gap that breaks down)
discharges it periodically through a receiver (US 685,957, pp. 1–3). Tesla
explains the charging as electrified particles thrown off by the source;
he notes the discharge effect is strongest on a negatively electrified
conductor and that the energy collected per unit time is nearly
proportional to the exposed area (p. 1–2).

**Standard explanation modeled by the engine.** All energy enters through
declared input channels (§3.1), each with its own ledger column:

- **Light / radiation:** photoemission from the plate. Photons above the
  plate's work function (mostly ultraviolet for common metals) eject
  electrons, so the plate charges *positive*, which matches Tesla's
  observation that the rays convey a positive charge to the plate.
  Emission current ∝ area × spectral irradiance × quantum yield, falling
  as the plate's own positive potential pulls electrons back.
- **Ionizing radiation:** X-rays or particles ionize the surrounding air
  and the plate surface, adding conduction between plate and air.
- **Electric field:** an elevated plate sits in the atmosphere's
  fair-weather field (order 100 V/m near the ground, larger in disturbed
  weather), so it collects charge from the air–earth current through the
  ionized air; this is the "elevated conductor" option in the patent.
- **Charging:** for space or high-altitude settings, a declared charging
  current onto the floating plate.

**Circuit.** Plate (capacitance to surroundings, leakage resistance),
condenser (high-quality dielectric; Tesla specifies mica), ground return,
threshold switch (`device/arc.go` for a gap, or a voltage-controlled switch
for the attracting plates), and a receiver load. Optional: a primary and
secondary transformer in the discharge circuit, as in the patent's Fig. 4.

**Outputs.** Charging current from each channel, condenser voltage over
time, discharge rate and energy per discharge, and the ledger showing every
joule delivered to the load traced to light, ionizing radiation,
atmospheric field or charging current.

**Validation (must pass before use).**

- Charging current scales with exposed area (Tesla's observation, p. 2)
  and with ultraviolet irradiance; a plate in visible light only, with a
  work function above the photon energy, shows no photoemission.
- Polarity: under ultraviolet, a negatively biased plate discharges and a
  positively biased one barely does (photoelectric behavior; Tesla, p. 1).
- With the light channel off, an elevated plate in a 100 V/m field charges
  at the rate set by the air–earth current and its collecting area.
- Energy per discharge equals ½CV² at the switch threshold, and the total
  delivered energy equals the integrated inflow minus leakage.

**Genetic evolution.** Genome: plate area, material (work function),
height, condenser value, threshold voltage, receiver load. Fitness: power
delivered to the load for a given set of declared inputs, or discharge
rate for signaling (the patent's interrupted-radiation signaling idea).

**Note on terms.** In this patent "radiant energy" means external
radiation and charged particles collected by the plate. It is a different
use from the pulse-driven "radiant energy" of the §3g pulse-duration
experiment; the two are kept separate.

---

## 3e. Tesla bifilar coil study (`examples/bifilar_coil/`)

Nikola Tesla's U.S. Patent 512,340 (1894, "Coil for Electro-Magnets") winds
two insulated conductors side by side and connects the end of the first to
the start of the second. In his example, a 1,000-turn coil with 100 V across
it has about 0.1 V between adjacent turns when wound normally, but about 50 V
when wound bifilar (p. 1, col. 2 to p. 2, l. 3). Because the capacity effect
goes with the square of that difference, he states that the energy stored in
the coil as a whole becomes 250,000 times as great, `(50 / 0.1)² = 250,000`
(p. 2, ll. 3–7). His purpose was to build the capacitance needed to cancel
self-inductance into the coil itself, instead of adding condensers (p. 1,
col. 1); he also notes that capacity obtained this way is evenly distributed
along the coil (p. 2, ll. 16–19). Claims: p. 2, ll. 39–51.

Tesla's U.S. Patent 593,138 (1897, "Electrical Transformer") gives the
companion condition for wave operation: the secondary's length is about one
quarter of the wavelength, based on the velocity of propagation through the
coil and its circuit, so that the potential is zero at the grounded end and
maximum at the far terminal (p. 1, col. 2, top; repeated with a worked
example on p. 2: 185,000 mi/s and 925 Hz give a 200-mile wavelength and a
50-mile secondary). It also takes the opposite approach to turn-to-turn
voltage: the flat spiral keeps the potential difference between adjacent
turns small so very high terminal voltages can be held.

Page and line numbers for both patents come from the printed margin
numbering and are approximate (±2 lines); check against the scans in
`reference/` when quoting.

### What the example does

1. **Reproduce the patent's idealized comparison.** Model both windings
   (ordinary and bifilar, same turn count and terminal voltage) and compute
   the turn-to-turn electrostatic energy. The idealized model must give
   Tesla's 250,000 ratio.
2. **Compute the full picture with realistic geometry.** Use `em/extract` to
   get each coil's inductance, distributed turn-to-turn capacitance and
   inter-winding coupling. Report:
   - electric (capacitive) and magnetic (`½LI²`) stored energy separately,
     and the total;
   - self-resonant frequency and Q;
   - impedance versus frequency, showing where self-inductance is cancelled.
   The ledger shows that all stored energy comes from the driving source.
3. **Use the reflected wave for maximum amplitude.** Treat the winding as a
   distributed line, not a lumped L and C: a wave launched at the drive end
   travels along the wire, reflects at the far end and returns. When the
   reflected wave arrives back in phase with the drive, the two add and a
   standing wave builds up to an amplitude limited only by losses (Q).
   - **Model.** The coil is a transmission line with distributed
     inductance and capacitance per unit length (`em/tline.go`), or a field
     region solved by TLM (`em/tlm.go`, §3a.1 item 8) when the geometry
     matters. A bifilar winding is a coupled two-conductor line: adjacent
     turns of the two strands are strongly coupled, so it has two
     propagation modes with different velocities, and both are modeled.
     The wave velocity along the wire is set by the geometry and insulation,
     not assumed to be the speed of light (Tesla 593,138 bases the quarter
     wavelength on the velocity through the coil itself). It is not
     necessarily slower than c: measured along the wire, a long free coil
     shows a phase velocity above c (up to about 2c in Drude's 1902 data),
     falling with frequency, while long coils in circuit tend to the speed
     set by the surrounding medium (Knight 2016, §1.2 p. 13; §13 points
     2–4, pp. 101–102). `em/helix.go` computes it with Knight's helical
     transmission-line model, which treats the coil as a helical line
     terminated by its own axial-field and fringe-field capacitances;
     inputs are length/diameter, former permittivity and surrounding-medium
     permittivity, so it also covers coils in liquids (`fluid/media.go`) (Knight §11,
     procedure pp. 95–96). The model is for single-layer coils on uniform
     cylindrical formers; bifilar windings use it for each strand's helical
     propagation and Paul's coupled-line model for the coupling between
     strands.
     Use Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed.
     (copy in `reference/`) for the line model: per-unit-length L and C of
     wires (§4.2, §5.2; classical cross-check in Steinmetz Lecture X,
     pp. 119 ff.), decoupling the coupled-line equations into modes
     by similarity transformation (§7.2.1, p. 284 ff.), and inhomogeneous
     media (§7.2.2.3) — insulated wires in air are inhomogeneous, which is
     why the two bifilar modes travel at different speeds. Paul's lines are
     straight TEM lines (§1.5); a wound coil is not, so the straight-line
     model is a first approximation and TLM is the check.
   - **Condition for an amplitude maximum.** The reflected wave adds to the
     drive when the total phase around one round trip is a whole number of
     cycles: `2βℓ + φ_far + φ_near = 2πn`, with `β = 2π/λ` and `φ` the phase
     of the reflection coefficient at each end. This gives both families of
     resonance:
     - **Odd quarter-waves**, `ℓ = λ/4, 3λ/4, 5λ/4, …`, when the two ends
       reflect with opposite sign (one end open, ρ = +1, and the other
       shorted or driven from a low impedance, ρ = −1). The voltage maximum
       is at the open end.
     - **Half-waves**, `ℓ = λ/2, λ, 3λ/2, …` (even multiples of λ/4), when
       both ends reflect with the same sign (open–open or short–short).
       Voltage maxima are at open ends, current maxima at shorted ends.
     - **Other fractions such as `λ/6` or `λ/8`** occur when an end is
       terminated reactively (a capacitance or inductance, e.g. a top load
       or the plate capacitor of §3e's first experiment), because the
       reflection then has a phase between 0 and π. The engine computes the
       resonant lengths from the actual terminations instead of assuming a
       fixed fraction, so the GA can use all of these.
     For pulse drive, the equivalent condition is repeating the pulse at the
     round-trip time `2ℓ/v` (or a whole multiple of it), so each new pulse
     adds to the returning reflection. A bifilar coil has a separate set of
     resonances for each of its two modes. For a single-layer coil Knight
     confirms both basic cases: a free coil (no connections) first resonates
     with wire length λ/2, and one end grounded gives λ/4 at about half the
     frequency, because grounding removes the reflection at that end
     (Knight §13 point 1, p. 101; §1.4). He also reports resonances near
     sub-multiples of the self-resonant frequency from multiple internal
     reflections, each with its own field pattern (§1.4, p. 19); pulse-
     driven runs look for these. Steinmetz gives the same families for
     line oscillations in 1911: ends of different character (one open, one
     shorted) give a quarter wave and its odd harmonics `(2k−1)f`, ends of
     the same character the half-wave family (Steinmetz, Lecture VII §31,
     pp. 81–82). Tesla's 593,138 transformer is the odd quarter-wave case: one end of the secondary grounded, the far
     terminal at maximum potential. Terminal conditions and reflection
     coefficients for coupled lines: Paul §7.3, §7.6 (p. 327), wave tracing
     §8.1.1 (p. 346), capacitive and inductive loads §8.1.5.1 (p. 370).
      - **Step-up at a junction.** A second mechanism, separate from resonant
     build-up: where a wave passes from a section of natural impedance `Z₁`
     into one of `Z₂`, voltage and current must match at the junction, so
     the wave's maximum voltage scales by `√(Z₂/Z₁)` and its current by
     `√(Z₁/Z₂)`, with partial reflection at the junction (Steinmetz,
     Lecture IX §41–42; his example: 200 Ω line into 800 Ω transformer,
     voltage doubled). Feeding a high-impedance coil from a low-impedance
     line therefore raises the coil voltage before any resonance. The
     natural impedance of each section is `√(L/C)` (Lecture VI §25), which
     also gives peak voltage from peak current as a quick check.
   - **Outputs.** Standing-wave voltage and current along the winding,
     peak voltage and its position, build-up time to steady state,
     resonance frequencies of both bifilar modes, and the ratio of peak
     amplitude to drive amplitude (the resonant step-up).
   - **Energy.** The step-up is resonant storage: amplitude grows until the
     loss per cycle equals the input per cycle. The ledger must show the
     stored and dissipated energy coming from the driving source, and the
     step-up must fall as Q falls (wire resistance, skin effect, dielectric
     and radiation loss).
4. **Evolve coil geometry.** Genome: turn count, pitch, wire diameter,
   insulation thickness and permittivity, number of strands and winding
   pattern (ordinary, bifilar, multi-filar, sectioned), plus far-end
   termination (open, shorted, loaded), feed-line natural impedance, and
   drive frequency or pulse repetition. Fitness options: target self-resonant frequency, maximum Q,
   maximum stored energy per volume at a given drive, **maximum standing-wave
   amplitude from the reflected wave** (peak voltage or current, or resonant
   step-up, at a given drive power), and wireless power-transfer efficiency
   between two evolved coils.
5. **Report against the patent.** Each run outputs the evolved coil next to
   the ordinary and Tesla bifilar baselines with the same turn count.

### Rules

- Cite the patent (page, lines) for every claim reproduced from it.
- The idealized 250,000 ratio test and the realistic-geometry results are
  separate tests; never present the idealized number as the coil's total
  energy gain.
- The neutrino channel is computed for every coil (as for every circuit) and
  reported; it is not assumed to couple to the coil beyond the established
  model (§1, long-term vision).
- Reflected-wave designs need the distributed model. The lumped coil
  formulas (`em/coil.go`) and their `λ/10` wire-length limit (Paul §6.4,
  pp. 265–266, which also suggests `λ/20` for accuracy) apply only to
  lumped runs; a run that evolves for reflected-wave amplitude must use the
  line or TLM model, where wire lengths from a fraction of `λ` (with
  reactive end loading) up to several `λ/2` are expected.
- Validation for the helical model (`em/helix.go`): free air-cored coils
  reproduce Knight's and Drude's measured velocity factors versus
  length/diameter (Knight §1.2); grounding one end roughly halves the
  first self-resonant frequency; pulse-driven runs show the sub-multiple
  resonances Knight reports.
- Validation for the distributed model: a uniform line shorted at one end
  and open at the other resonates at `ℓ = λ/4, 3λ/4, …`; open–open and
  short–short lines resonate at `ℓ = λ/2, λ, …`; a line ended in a known
  capacitor resonates at the length given by the round-trip phase
  condition; in each case the step-up equals the analytic value for its Q; pulses
  repeated at `2ℓ/v` build up and pulses at other intervals do not; the
  two bifilar mode velocities match the coupled-line analytic values for a
  simple geometry (Paul §7.2).
- Turn-to-turn voltage is a design trade-off, not a free gain: 512,340
  raises it to build capacitance, 593,138 lowers it to hold high terminal
  voltage. Every evolved coil is checked against insulation breakdown at
  the peak standing-wave voltage between adjacent turns, and designs that
  exceed it are rejected.

### As implemented: first experiment (`examples/bial_bifilar/`)

A bismuth–aluminium capacitor element in parallel with a bifilar coil,
self-driven, evolved for Q at resonance against stored energy per volume
(weighted-sum sweep, Pareto front), in two variants:

- **dry** — Bi and Al plates with a solid dielectric; the Bi–Al leads are a
  thermocouple, EMF `(S_Al − S_Bi)·ΔT` (75.5 µV/K) in series with the plate
  capacitance; energy enters through the **thermal** port;
- **wet** — Al|electrolyte|Bi galvanic cell (≤ 1.970 V standard EMF),
  electrolyte resistance, Randles electrodes; energy enters through the
  **chemical** port, capacity by Faraday's law;
- **thermopile** — the dry element charged by `n` Bi–Al pairs in series
  (EMF × n, stored energy × n²); leg resistance, leg volume and the legs'
  conductive heat leak `n(k_Bi + k_Al)·A/ℓ·ΔT` are modeled (legs ≥ 5
  diameters long).

Coil formulas are analytic (`em/coil.go`: Wheeler inductance, turn-to-turn
capacitance with `V/N` ordinary and `V/2` bifilar per gap, skin effect), a
first approximation until `em/extract`; wire longer than `λ/10` at
resonance is rejected.

**Temporal-field test (hypothetical, §3f).** The Bi–Al element is also the
bench test of Jack's temporal-field model: under rapid heating, with an
imposed anisotropy (temperature gradient or isothermal current), the model
predicts an EMF from the time variation of temperature, with magnitude set
by the Thomson coefficients and conductivities of bismuth and aluminium.
The standard prediction (Seebeck EMF from the instantaneous temperature
difference only) is computed alongside, and the difference between the two
is the measurement target.

**Superseded for prediction.** Knight (2016, §6 pp. 53–54; §13 point 8)
shows that turn-to-turn capacitance theories have no predictive power for
single-layer solenoids: the coil's end-to-end (axial-field) capacitance
does not change with the number of turns. The `V/N` ordinary-coil model is
therefore kept only to reproduce Tesla's idealized 250,000× comparison
(§3e step 1). For predictions, ordinary-coil self-capacitance and
self-resonance come from `em/helix.go` (Knight §5 corrected formula, §11
SRF model), and the coupling between bifilar strands — a real capacitance
between two conductors at different potentials — comes from the coupled-
line model (Paul). The first run's Q and resonance results used the
superseded model and must be re-run before they are relied on. Each front design is re-simulated with an ordinary
coil, and its resonance is reported as the nearest whole-number
subharmonic of the strongest Cu, Al or Bi emission line (NIST; adjacent
subharmonics are only ~1/n apart, n ~ 10⁹–10¹¹, so this is reported, not a
real constraint). A source can declare the channel it draws on
(`circuit.Element.Channel`), so the ledger attributes the energy to heat or
chemistry. Runs: `doc/experiments/bial_bifilar_run1.md` (dry, wet),
`bial_bifilar_run2.md` (adds thermopile: 31 J/m³ at 5000 pairs, at a
~48 W heat leak).

---

## 3f. Lorentz force, the homopolar dynamo and its electric dual (`em/lorentz.go`, `em/duality.go`)

Reference: Peeter Joot, *Geometric Algebra for Electrical Engineers*
(GAelectrodynamics, v0.3.6.4), Theorem 3.16 "Lorentz force and power", source
repository `github.com/peeterjoot/GAelectrodynamics`. **No license file was
found in that repository, so treat it as reference only:** cite theorem and
equation numbers, implement the equations ourselves, and do not copy its
LaTeX or text into this project (or add it as a vendored dependency) unless
the author grants a license.

### The theorem

With field multivector `F = E + I η H` (equivalently `E + I c B`) and a
charge/velocity multivector `Q`, force and power on a charged particle come
from one product, split by grade:

```
(1/c) dT/dt = ⟨F Q†⟩₀,₁ = ½ (F† Q + F Q†)
grade 0 (scalar) → power dℰ/dt,   grade 1 (vector) → force c dp/dt
```

**Electric charge `q`, velocity `v`:**

```
dp/dt = q (E + v × B)          dℰ/dt = q E · v
```

**Magnetic charge `q_m`, velocity `v_m` (the dual):**

```
dp/dt = q_m (c B − (1/c) v_m × E)     dℰ/dt = q_m c B · v_m
(q_m here in coulombs. Erratum in the source, checked against
 reference/GAelectrodynamics.color.V0.3.6.4.pdf, pp. 239–240: Joot defines
 q_m = ∫ρ_m dV with ρ_m in Wb/m³ (§3.1), so Definition 3.10 should read
 Q = I(q_m/η)(1 + v_m/c). Theorem 3.16 as printed omits that 1/η in the
 force and gives the power as (1/η) q_m B·v_m/c; with q_m in Wb the force
 is q_m(H − ϵ v_m×E), agreeing with his Theorem 3.14, and the power is
 q_m H·v_m = (1/η) q_m c B·v_m — his own proof, eq. 3.183, gives grade 0 =
 q_m B·v_m before the scaling by c. The code uses q_m [A·m] = q_m[Wb]/μ₀
 = c × the coulomb value: q_m(B − v×E/c²) and q_m B·v.)
```

`em/lorentz.go` implements the grade-selection form with `hyper/ga` and must
agree with the conventional vector forms above (twin tests, `1e-12`
relative).

### Homopolar dynamo — Tesla U.S. Patent 406,968 (electric pair)

Tesla's 1889 "Dynamo-Electric Machine" patent is a homopolar (unipolar)
generator: a conducting disk spinning in an axial magnetic field. The
electric pair explains it directly. Charges in the disk move with
`v = ω × r`, feel `q v × B`, and are pushed radially, giving an EMF between
axle and rim:

```
EMF = ½ B ω (r_outer² − r_inner²)
```

The mechanical torque needed to keep the disk turning equals the electrical
power drawn, so the ledger shows inflow from the shaft (a declared mechanical
port) and outflow to the load. `examples/homopolar/` reproduces the patent's
arrangement, checks the EMF formula, then evolves disk geometry, field
arrangement and brush placement for output voltage, current and efficiency.

### The inverse analogy: a rotating electric field and magnetic fields

The dual pair suggests swapping roles: a rotating *electric* field acting on
*magnetic* charges. The engine treats two different physical questions
separately and never mixes them:

1. **A changing or rotating E field produces a B field — real, standard
   physics.** By the Ampère–Maxwell law (Treatise form `4π ℭ = V.∇ℌ` with
   displacement current `∂𝔇/∂t`), any time-varying electric field produces a
   magnetic field. A rotating charged disk also produces a magnetic field
   because its charge moves (Rowland's experiment, 1876).
   `examples/rotating_efield/` computes these real fields with `em/` and
   reports their size.
2. **Force on magnetic charges — the dual pair — requires magnetic
   monopoles.** No magnetic monopole has ever been observed, despite
   dedicated searches (e.g. the MoEDAL experiment at CERN, and searches for
   relativistic monopoles in the IceCube Neutrino Observatory, whose
   principal investigator Francis Halzen received the 2026 Nobel Prize in
   Physics for IceCube's high-energy neutrino work). Monopole results will be added here
   as cited evidence when the user supplies them. In normal mode,
   `q_m = 0` everywhere, so the dual pair gives zero force and zero power.

**Hypothetical mode.** `em/duality.go` implements dual-symmetric Maxwell
equations with magnetic charge and current, `∇F = J_e + I J_m`, so the dual
pair and the inverse analogy can be explored as theory. Rules:

- Off by default; enabled only by an explicit `Options.HypotheticalMagneticCharge`.
- Every result produced with it is labeled **hypothetical** in output, reports
  and GA logs, and its energy is kept in a separate ledger column that is
  never counted as real inflow.
- Neutrinos carry no electric or magnetic charge, so neither pair couples to
  them; the neutrino channel stays as defined in §1 and §3.1.

**Second hypothetical model: quaternionic electrodynamics (`em/hed.go`).**
Giardino (arXiv:2010.07748, 2020) extends Maxwell by letting the vector
potential carry extra quaternion components, `A = αi + βj` (α real, β
complex). The magnetic field gains a self-interaction term,
`B = ∇×A + A×A` (eq. 3), and magnetic charge and current arise from it
dynamically: `∇·(A×A) = −4πρ_M` (eq. 14) and `J_M = (1/4π) ∂(A×A)/∂t`
(eq. 16), with `E = −∇Φ − (1/c)∂A/∂t` unchanged (eq. 18); constant
factors are to be checked against the PDF's typeset equations, as the
text layer garbles them. Unlike the
point-source mode above, the monopole density here is produced by the
field's own structure. Rules:

- Off by default; enabled by `Options.HypotheticalHED`; labeled
  **hypothetical** everywhere, with its own ledger column, as above.
- The paper gives no mechanism by which a laboratory source creates β, so
  β enters as a declared free parameter (amplitude and spatial profile, or
  a coupling coefficient to a chosen source quantity, stated per run). An
  experiment yields either a fitted β or an upper bound on it.
- Tests: with β = 0 the model equals Maxwell to round-off; the magnetic
  continuity equation `∇·J_M + ∂ρ_M/∂t = 0` (eq. 17) holds; Giardino's
  Poynting theorem `∂U/∂t + ∇·S = ⟨E,J⟩ + ⟨B,J_M⟩` (eq. 28), with
  `U = (|E|² + |B|²)/8π`, closes in the ledger; the energy flux is
  orthogonal to E and B, `⟨S,E⟩ = ⟨S,B⟩ = 0` (eq. 30); the radial monopole
  solution of the paper's Sec. III is reproduced.

**Third hypothetical model: the temporal field (`em/temporal.go`).**
P. M. Jack (arXiv:math-ph/0307038, 2003) keeps both the left and the right
quaternion derivative of the potential `A = U + A₁i + A₂j + A₃k`, defining
`E = −{d/dr, A}` (symmetric) and `B = [d/dr, A]` (antisymmetric) (eqs.
9–10). This adds a seventh field component, the scalar **temporal field**
`T = −(1/c)∂U/∂t + div(A)` — the quantity that the Lorenz gauge sets to
zero in standard theory — and modifies two equations (eqs. 25–31):

```
curl(B) = (1/c)∂E/∂t + grad(T) + 4πJ/c      div(E) = (1/c)∂T/∂t + 4πρ
curl(E) = −(1/c)∂B/∂t                       div(B) = 0
```

Jack links T to reversible heat: a charge in T exchanges energy without
moving, which he identifies with Peltier, Thomson and Seebeck heat, and he
relates the medium's `dT/dK` (K = temperature) to the Thomson coefficient,
`h_T = −(c/4πσ)·dT/dK` (eq. 37). Rules:

- Off by default; enabled by `Options.HypotheticalTemporal`; labeled
  **hypothetical** everywhere, with its own ledger column for T-energy.
- The coupling is not free: in each material `dT/dK = −4πσ·h_T/c` is set
  from tabulated Thomson coefficient `h_T` and conductivity `σ` (sources
  cited per material), so the model makes numerical predictions.
- The paper gives no energy theorem for T and does not show charge
  conservation with T physical; both must be derived and checked before
  results are used.
- Tests: with `T = 0` the model equals Maxwell to round-off; charge
  conservation `∂ρ/∂t + div(J) = 0` is checked explicitly in every run and
  any violation reported; the derived energy theorem closes in the ledger;
  the Thomson-heat equation `dQ/dt = J²/σ − h_T J·grad(K)` (eq. 36) is
  reproduced for a wire in a temperature gradient.
- Predicted signatures: a longitudinal E pulse with no accompanying B
  (pulse-duration check 2, §3g); and an EMF induced by a temperature that
  varies in time in an anisotropic medium — the effect Bridgman noted as
  expected but unobserved (Jack Sec. III.B) — predicted in magnitude for
  the §3e Bi–Al element under rapid heating.

### Tests

- Grade-selection force/power equals the vector forms, for both pairs.
- Electric–magnetic duality transform (`E → cB`, `cB → −E`, `q → q_m`; `q → q_m/c` with q_m in A·m)
  maps the electric pair onto the magnetic pair.
- Homopolar EMF matches `½Bω(r_o² − r_i²)` and the shaft-power balance.
- Displacement-current B field of a known time-varying E field matches the
  analytic result (e.g. a charging parallel-plate capacitor).
- Röntgen–Eichenwald, the electric mirror of the homopolar disk: a dielectric
  disk spinning in a fixed axial E field is magnetized, M = P × v, giving
  B = −(εr − 1)·E·v/c² in its mid-plane (ring sum of the moving bound face
  charges ±P vs the thin-disk value, `examples/rotating_efield` case 1e). The
  1/c² is μ₀ε₀, not a fitted factor; it is why the mirror effect is nT where
  the homopolar one is volts: 1 T at 100 m/s gives 100 V/m, while 1 MV/m at
  100 m/s gives about 1 nT.

### Momentum ledger and push–pull designs

Beside the energy ledger, every force the engine computes is booked in a
**momentum ledger**: force on each element or body, and the reaction it
hands to each exterior channel (the source of an external field, the
plasma, radiation). An isolated system's internal forces sum to zero; a
predicted net thrust must show where its momentum goes, just as surplus
energy must show its source. Designs that push or pull against fields in
space — electrodynamic tethers, magnetic and electric sails, charged
bodies used for buoyancy — follow the design notes in
`doc/em-push-pull.md`, recorded at the user's request for future circuit
designs.

---

## 3g. Arc surges and toroidal field structures (`examples/arc_surge/`)

When a circuit arcs or is driven by a very fast pulse (nanoseconds and
shorter), a sharp current front travels down the conductors, surrounded by a
ring-shaped (toroidal) magnetic field that moves with it. Some pulse
geometries also launch field structures that detach from the circuit and
travel on their own. This example simulates these events end to end and
measures the field's shape.

### Pipeline

1. **Arc model (`device/arc.go`).** Arc conductance evolves by the Cassie
   model (high-current arcs), the Mayr model (near current zero), or a
   combination. An arc's falling voltage–current characteristic acts as a
   negative resistance that can sustain or grow oscillations; Steinmetz
   treats this as a traveling wave with negative power-transfer constant,
   with the arc as its source (Lecture VIII §35). The ledger must show the
   energy of such oscillations drawn from the supply feeding the arc. Parameters: time constant, power-loss coefficient, breakdown
   voltage. The arc is an ordinary nonlinear device in the circuit; for a
   gap in a liquid or a mist its parameters come from published
   measurements in that medium.
2. **Circuit transient.** `.TRAN` with timesteps fine enough to resolve the
   fastest edge (rule: at least 20 points per rise time; `limits.go` warns
   otherwise). Conductors that are long compared with the edge length are
   modeled as transmission lines (`em/tline.go`), so the pulse travels,
   reflects and decays.
3. **Local field solve (`em/fdtd.go`).** Branch currents and node voltages
   from the circuit drive a Yee finite-difference time-domain solve of the
   space around the pulse source (wires, gap, board and surrounding air or
   vacuum), with absorbing boundaries. This gives **E(t)** and **B(t)**
   everywhere in the region, including fields that leave the circuit.
   A TLM solve (`em/tlm.go`, §3a.1 item 8) can be used instead or
   alongside: it treats the region as a Kron-style network, attaches the
   arc and other lumped elements directly to field nodes, and exposes the
   region boundary as an exchange with the outside.
4. **Field diagnostics (`em/topology.go`).** At every saved time:
   - **Divergence of B** (the scalar part of the quaternion `∇𝔅`), mapped over
     the whole region. The Yee scheme keeps it at round-off level by
     construction, so any value above that is reported with its location as
     a numerical-quality warning.
   - **Field-line tracing:** closed loops, ring (toroidal) structures, their
     radius, thickness and position.
   - **Toroidal moment** of the current and field distribution.
   - **Detachment and propagation:** whether a ring structure separates from
     the conductors, how fast it moves, how long it persists, and how its
     energy decays (Poynting flux through surfaces around it).
   - **Energy accounting:** field energy in each structure, tracked in the
     ledger against the circuit energy that produced it.
5. **Virtual probes (`em/probes.go`).** B-dot loops, D-dot sensors and
   current probes placed anywhere in the region output exactly what a real
   probe would record, including probe bandwidth and loop area.
6. **Measurement comparison.** Real oscilloscope and probe data (CSV) can be
   imported, time-aligned with the simulation and compared probe by probe,
   with the difference reported per channel.
7. **Optional hypothetical comparison.** With
   `Options.HypotheticalMagneticCharge` (§3f), the same run can also compute
   the field pattern that a given magnetic charge would add, so its probe
   signature can be compared with the standard result and with measurement.
   Labeled hypothetical, as in §3f.

### Validation cases (must pass before use)

- Cassie and Mayr arc models reproduce published current-zero behavior.
- A step pulse on a lossless line travels at the expected speed and reflects
  correctly from open, shorted and matched ends.
- The B field around a straight wire matches Ampère's law; a pulse on a
  coaxial line matches the analytic TEM fields.
- The FDTD solver reproduces the analytic Hellwarth–Nouchi "flying doughnut"
  toroidal pulse propagating in free space, including its ring structure
  and speed.
- `∇·B` stays at round-off level in all of the above.
- When `em/tlm.go` is used: TLM and FDTD agree on the resonant modes of a
  closed rectangular cavity and on the step-pulse line case above; a lumped
  resistor attached through a half-cell stub absorbs the same power in TLM
  as in the circuit solve.

### Genetic evolution

The genome covers the pulse circuit (source, storage capacitor, switch or
gap, inductance, transmission-line lengths) and the electrode/antenna
geometry around it. Fitness options, all computed by the diagnostics above:

- energy carried by detached ring structures, as a fraction of input energy;
- persistence time and travel distance of a detached structure;
- toroidal moment at a target position;
- sharpness of the pulse edge, peak dB/dt at a probe;
- for experiments: the probe positions and circuit design that make the
  standard and hypothetical predictions most different, giving the clearest
  bench test.

### Pulse-duration experiment (`examples/pulse_duration/`)

Tesla's 1891 patent 462,418 charges a condenser from a high-tension source
and discharges it disruptively through an air gap into a working circuit,
with the stated aim of making the duration of each impulse extremely small;
the discharges "may be of the same direction or alternating" (US 462,418,
p. 1). His 1892 lecture separates the "fundamental sudden discharges" of a
disruptive-discharge coil from the oscillations superimposed on them, and
reports that instruments respond to the former and that electrodes wear far
faster under them than under harmonic currents (lecture, pp. 12, 33). His
1901 radiant-energy patent adds that a Roentgen or Lenard tube is effective
only when the exciting impulses are "wholly or at least preponderatingly of
one sign" (US 685,957, p. 3). Specific duration thresholds are reported only
secondhand (Vassilatos), from Tesla's later work; their primary source is
still to be found. This experiment reproduces that setup, sweeps the
impulse duration and repetition, and records what the field does, so the
standard prediction and the project's working premise can be compared at
real probes.

**Working premise (labeled hypothesis).** The pulse launches outward
"ether streams" that appear as high-voltage electric fields moving away from
the source, and those electric fields in turn create toroidal magnetic
fields wrapped around the direction of the electric field, travelling
outward with it. The premise is recorded as a hypothesis, not built into
the field solver.

**What standard Maxwell already predicts.** A changing electric field is a
displacement current (`4π ℭ = V.∇ℌ` with `ℭ = 𝔎 + ∂𝔇/∂t`, §3c), so a
pulsed E field along a direction produces ring-shaped B around that
direction, moving outward with the pulse. The engine computes this
standard field (FDTD or TLM). The experiment therefore asks where the
premise differs from it, and looks for those differences:

1. **Speed.** Standard fields propagate at the speed set by the medium.
   A front measured faster or slower at spaced probes is a deviation.
2. **Unaccompanied electric field.** Standard radiated fields carry E and B
   together (in vacuum `|E| = c|B|`). An outward E pulse with a weaker or
   absent ring B, or with the ring B arriving later than E, is a deviation.
3. **Ring geometry.** Standard theory puts the B rings around the local
   ∂E/∂t. Rings whose axis follows E where ∂E/∂t predicts otherwise, or
   rings that persist after ∂E/∂t has ended, are deviations.
4. **Duration thresholds.** Standard response follows the pulse's
   frequency spectrum smoothly. Abrupt changes in effect at particular
   durations that the spectrum does not explain are deviations. The ranges
   to watch (around 100 µs, 1 µs and shorter) are reported secondhand by
   Vassilatos, not in the Tesla sources in `reference/`; treat them as
   search targets, not as Tesla's stated values.
5. **Polarity.** Standard theory treats unidirectional and alternating
   pulses of the same spectrum the same way (for the standard distinction
   between unidirectional and oscillating discharges of a double-energy
   circuit, see Steinmetz, Lecture VI §24–27). A difference that the
   spectrum does not explain is a deviation. Two specific tests from
   Tesla: (a) a brush discharge rotating in a fixed direction in a
   magnetic field, reversed by reversing a permanent magnet, which he took
   as a sign that positive and negative impulses are unequal (1892
   lecture, p. 18); the run records rotation direction against pulse
   polarity and magnet orientation; (b) an X-ray-type tube driven by
   one-sign versus symmetrical impulses (US 685,957, p. 3).
6. **Shielding.** Standard fields are blocked by a closed conducting
   enclosure at these frequencies. A probe signal inside one is a deviation.
7. **Energy.** The ledger tracks field energy against source energy (§3c).
   Energy at probes or in detached structures beyond what the source
   supplied is a deviation. The ledger also separates energy lost locally
   in the surrounding gas (molecular impacts near charged surfaces) from
   energy radiated away: Tesla argued that most of the energy leaving a
   rapidly alternating charged body goes the first way, not into "long
   ether waves" (1892 lecture, pp. 36–37). Each configuration is run in air
   and in vacuum (or a gas of different pressure) to separate the two.

Each deviation is defined as a measured value outside the simulation's
predicted value by more than the combined probe and model uncertainty.

**Setup and sweep.**

- **Source:** a storage capacitor charged by a supply, discharged through a
  spark gap (`device/arc.go`) or a fast
  switch, giving unidirectional pulses. Optional series diode or magnetic
  quench to suppress reversal. A matching alternating source with the same
  spectrum is provided for the polarity comparison.
- **Load:** a straight discharge path, a ball terminal, or the coil of §3e
  (including its reflected-wave tuning).
- **Sweep:** impulse duration from 100 µs down to the shortest edge the
  circuit allows (decades, finer near the reported thresholds), repetition
  rate, charge voltage and polarity.
- **Probes (`em/probes.go`):** E (D-dot) and B (B-dot) sensors in pairs at
  several distances and angles, including inside a closed conducting box
  for check 6; probe timing is recorded for check 1.

**Outputs per run.** E(t) and B(t) at each probe, the E/B ratio and arrival
times, ring structures and their axes (`em/topology.go`), field energy
against source energy, and, for each of checks 1–7, the predicted value,
the measured value when bench data are imported, and whether it is a
deviation. When a discharge lamp or glow is the load, its light output is
also compared with Tesla's stated expectation that heat and light go with
frequency times the square of the potential (1892 lecture, p. 9).

**Optional comparisons.** The same runs can be repeated in the hypothetical
magnetic-charge mode (§3f), so its probe signature can also be compared,
and in the quaternionic-electrodynamics mode (`em/hed.go`, §3f), whose
signature for a declared β is a measurable `∇·B ≠ 0` and a mismatch between
`∇×E` and `−∂B/∂t` at paired probes. With no deviation measured, the run
reports the upper bound on β the probes allow. The temporal-field mode
(`em/temporal.go`, §3f) predicts a longitudinal E component with no
accompanying B, with its magnitude fixed by the materials' Thomson
coefficients.
Any further hypothesis the user specifies is added as its own predicted
signature, labeled hypothetical, alongside the standard one.

**Validation (must pass before use).**

- Capacitor discharge through a resistor gives the analytic RC pulse;
  through a gap, the Cassie/Mayr validation of §3g.
- A short current element radiates the analytic Hertzian-dipole field,
  including ring B around the element and `|E| = c|B|` far away.
- A closed conducting box attenuates an external pulse by the analytic
  shielding factor.
- Steinmetz's impulse on an artificial 130-mile line (Lecture VIII §37,
  Figs. 51–53): r = 93.6 Ω, L = 0.3944 H, C = 1.135 µF (natural impedance
  about 590 Ω), impulse duration about 3.6 ms, far end terminated in a
  resistance equal to the natural impedance (Fig. 51) and in the other
  terminations of Figs. 52–53. The simulated wave front, attenuation and
  reflections match the published oscillograms in shape and timing.
- Unidirectional and alternating sources with identical spectra give
  identical standard predictions (so checks 4 and 5 measure only
  deviations).

**Genetic evolution.** Genome: pulse circuit, gap or switch, load geometry
and probe placement. Fitness: the design and probe placement that maximize
the predicted separation for any of checks 1–7 at a given source energy,
so bench time goes to the clearest test.

**Sources.** Tesla's own descriptions are the primary sources, all in
`reference/`: the 1892 lecture "Experiments with Alternate Currents of High
Potential and High Frequency" (transcription from tesla.hu; page numbers are
that PDF's pages), US 462,418 (1891) and US 685,957 (1901). Patent text was
read by OCR; check quotations against the scans.
G. Vassilatos, *Secrets of Cold War Technology* (copy in `reference/`), is
background history and the source of the premise's "ether stream" wording;
it gives no equations or measurements and is not cited for physics.

---

## 3h. Liquid media (`fluid/`)

Properties are tables versus temperature (and pressure where it matters),
each with a cited source: density, viscosity, surface tension, sound speed,
vapor pressure, permittivity, electrical conductivity (set by purity and
dissolved ions), heat capacity and thermal conductivity. **Heavy water (D₂O)
and light water (H₂O) are separate media** — D₂O is denser, more viscous and
has a lower sound speed and vapor pressure, all of which change breakdown
and acoustic behavior. Any other liquid can be added as a table. The
liquid's permittivity and conductivity are also used by `em/fdtd.go` for
the field solve, and by `em/coil.go` for coils in liquids. `fluid/` and
`material/` are kept small for now and will be expanded later.

### 3i. Multipactor discharge (planned)

In vacuum, an RF field drives electrons across a gap; when the impact energy
gives a secondary-emission yield above 1 and the transit time is an odd
number of RF half-periods, the electron population grows into a resonant
discharge. Farnsworth built active devices on it; the patents to review
(details checked on Google Patents):

- US 1,969,399, "Electron multiplier" (filed 3 Mar 1930, issued 7 Aug
  1934; Television Laboratories) — secondary emission from repeated impacts
  along a hollow resistor with a graduated electrostatic field.
- US 2,091,439, "Multipactor oscillator and amplifier" (filed 24 Feb 1936,
  issued 31 Aug 1937) — an electron cloud oscillating between two opposed
  surfaces with secondary yield above 1; part of the cloud is diverted to
  deliver power, giving an oscillator or amplifier.
- US 2,135,615, "Multipactor" (filed 11 Feb 1936, issued 8 Nov 1938) —
  photocathode and a series of perforated secondary cathodes, each stage
  releasing more electrons than it receives.
- US 2,174,487, "Self-energized alternating current multiplier" (filed
  22 Mar 1937, issued 26 Sep 1939) — one electrode structure self-oscillates
  and its AC energizes a second multiplier (e.g. a photocell) in the same
  envelope. "Self-energized" means the tube generates its own RF drive from
  its DC supply; the ledger must show that supply as the energy source.
- US 2,179,996, "Electron multiplier" (Farnsworth and R. L. Snyder; filed
  9 Nov 1936, issued 14 Nov 1939) — perforate anode and potential gradients
  drift the electron cloud along secondary-emissive surfaces for many
  impacts; gain set by the gradient voltage.
- US 2,141,837, "Multistage multipactor", and US 2,180,279, "Method of
  operating electron multipliers" (already listed; still to be read).

**What the patents establish** (texts in `reference/`: US1969399, US2091439,
US2135615, US2174487, US2179996; read from their text layers):

- **Mechanism (US 2,091,439, p. 1).** Operation depends on energy transfer
  between an oscillating electron cloud and an oscillating EMF between two
  opposed cathodes whose secondary-emission ratio exceeds 1. When energy
  flows from the EMF to the cloud, the electrons strike the cathodes and
  multiply; when energy is taken from the cloud, it reinforces the
  oscillating field and the tube self-oscillates. The energy in the
  self-oscillating case is "supplied from a D. C. source applied to an
  anode" between the cathodes — the ledger books it there.
- **Transit-time tuning (US 2,091,439).** The anode potential sets the
  electrons' time of flight between cathodes so that it matches the RF
  period; electrons in the decelerating phase stop striking the cathodes
  and deliver their energy to the circuit, and an equilibrium current
  results. Separate multiplying and collecting structures (a collector
  cathode with a non-emitting, carbonized-nickel surface and its own
  output resonant circuit) let power be taken out without upsetting the
  multiplication.
- **Surfaces.** Treated surfaces give up to about six secondaries per
  primary at the right impact velocity; even untreated nickel exceeds
  unity (US 2,091,439).
- **Staged and drift multipliers.** US 1,969,399: graduated field along a
  hollow resistor, linear amplification. US 2,135,615: photo-emissive
  primary cathode and a series of perforated secondary cathodes at a rising
  potential gradient; no filament, so it can run above ground potential,
  and DC or AC excitation. US 2,179,996 (with R. L. Snyder): electrons
  drift along continuous emissive surfaces under a potential gradient,
  reducing space-charge limits; gain set by the gradient voltage.
- **Space-charge limit (US 2,174,487).** Self-oscillating multipliers settle
  at an equilibrium where space charge destroys linearity, which is why
  this patent uses one self-oscillating structure to drive a separate
  linear multiplier (e.g. a photocell) in the same envelope.

**Model consequences.** The §3i gap model needs: secondary-emission yield
versus impact energy per surface (including a low-yield collector
surface), transit time set by gap voltage and spacing, phase sorting of
electrons, space-charge saturation of the equilibrium current, and an
anode DC supply as the energy source for self-oscillation. Validation adds:
equilibrium current saturates with space charge; RF output power never
exceeds DC anode input minus losses.

On spacecraft the
same effect is a failure mode of RF filters, waveguides and connectors
(ECSS-E-ST-20-01, multipaction design and test). Planned, after the patents
are read:

1. **Susceptibility check** in vacuum: for each RF gap
   of an evolved design, gap voltage, frequency and gap length against the
   resonance zones from the surface's secondary-emission-yield curve;
   designs inside a zone are rejected or penalized in vacuum runs.
2. **Device model**: a multipactor gap as a nonlinear element, discharge
   current from secondary-emission yield and transit-time resonance, all of
   its energy drawn from the RF field and booked in the ledger.

### 3j. Ground charge pump and staged Royer step-down with magnetic amplifiers (`examples/charge_pump/`, planned)

**Design intent.** An electron-pumping action draws electrons from the
ground (the negative side), stores them in high-voltage capacitors, and
then steps the stored energy down — from high voltage and high frequency to
low voltage and low frequency — through staged Royer oscillators, with a
magnetic amplifier between every pair of Royer stages.

**Energy versus charge.** The earth is a charge reservoir and the potential
reference: it can supply electrons, as in Tesla's 685,957, where the
grounded terminal "may be considered as a vast reservoir of negative
electricity". It does not by itself supply energy. Every joule stored in
the bank comes from a declared source — radiation and the atmospheric field
(§3d radiant receiver), a drive supply (pulse source §3g, multipactor anode
supply §3i), or the field of a resonator (§3e). The ledger tracks **charge
drawn from ground** and **energy by source** as separate quantities; a run
that stores energy with no source behind it is a bug, and is reported as
one.

**Stages.**

1. **Front end: spark gaps in water mist (the nonlinear stage).** A
   supply charges a tank capacitor that discharges disruptively through
   one or more spark gaps (US 462,418) with a fine mist of water between
   the electrodes. The arcs are the strongly nonlinear element: an arc's
   falling voltage–current characteristic acts as a negative resistance
   that sustains and sharpens oscillation (Steinmetz Lecture VIII §35),
   and the mist changes the gap in three ways the model must capture:
   - **Breakdown — depends on droplet size and conductivity.** Fine clean
     fog (mean droplet ≈ 30 µm) *raised* the AC discharge voltage of
     5–15 cm rod–plane gaps by 15–35% over dry air, by lowering the peak
     field near the rod; raising the fog conductivity from 18 to
     4000 µS/cm took back 4–9% of that (Liu et al. 2018). Water vapour
     alone raises breakdown by about 2% (humidity results summarized by
     Allen & Phillips 1959). Large drops and water streams act the other
     way: breakdown is significantly lower than dry, falling further with
     stream velocity, size and conductivity (Yuan et al. 2014). The model
     must therefore follow droplet size, number density and conductivity,
     not just "wet or dry".
   - **Channel:** water vapour and droplets in the channel change its
     conductivity and absorb energy (evaporation, dissociation).
   - **Recovery:** after a spark the channel de-ionizes in tens of
     microseconds, but the hot, low-density gas takes milliseconds to cool,
     which limits repetition to a few hundred hertz without gas flow;
     electron-attaching gases and flow shorten recovery (Bhoomireddy 2016).
     That mist cooling of the channel shortens thermal recovery, giving a
     faster, more repetitive pulse train, is a **design hypothesis to be
     tested** — no measurement of it was found yet; the experiment
     measures recovery time with and without mist.
   The resulting pulses drive the bifilar coil (step 2) and the first
   magnetic amplifier (step 6), saturating its core each pulse. Two
   topologies are built and compared on the same front end:
   - **Topology A (series):** gap → bifilar primary → first magnetic
     amplifier; the coil filters and stores the pulses before the
     amplifier sees them.
   - **Topology B (parallel, preferred):** the gap pulses drive the first
     magnetic amplifier directly, with the bifilar coil connected in
     parallel across the same pulse node, charging and storing alongside
     it. The amplifier sees the sharp gap pulses; the coil adds storage
     and resonant support on the same node.
   **Bench mist source.** An ultrasonic nebulizer (mist maker) gives
   droplets of a few micrometres: by Lang's relation (Lang 1962) the
   median diameter is about `0.34·(8πσ/(ρf²))^(1/3)`, ≈ 3 µm for water at
   1.7 MHz. Use reverse-osmosis or de-ionized water: with tap water an
   ultrasonic humidifier also emits mineral particles (Na, Mg, Si, S, Ca;
   peak ≈ 180 nm), which RO water largely removes (Umezawa et al. 2013).
   Water conductivity is then a controlled variable (add salt to raise
   it), and it is a GA gene. Theatrical fog machines make glycol or
   glycerin fog, not water mist, and are not a substitute. The bench
   measures droplet size and liquid water content in the gap so runs can
   be matched to the model.
   Both run in every experiment of this section with identical gaps,
   coil and amplifier, and the report compares them on energy delivered to
   the Royer chain per joule from the supply, peak and RMS pulse into the
   amplifier, core reset margin, coil stored energy, and losses by element.
   Topology B is the default until a comparison says otherwise. Other
   front ends remain available: the radiant-energy plate (§3d) and a
   multipactor (§3i).
2. **Bifilar coil primary as storage.** The gap discharges into the
   primary of the §3e bifilar coil (in series ahead of the amplifier in
   topology A, in parallel with it in topology B). The coil stores energy in its magnetic
   field (`½LI²`), its distributed and inter-strand capacitance (`½CV²`),
   and, when tuned, in the reflected-wave standing wave (§3e step 3). This
   adds storage capacity to the chain: energy delivered by the gap is held
   in the coil and passed on, smoothing the pulses. The ledger shows the
   coil's stored energy rising and falling with every pulse, with all of
   it traced to the supply; storage does not add to the total.
3. **Pump / rectifier.** High-voltage diode stage or diode ladder
   (Greinacher / Cockcroft–Walton) oriented so that electrons flow from
   earth into the storage bank's negative terminal. Models: HV diode with
   reverse recovery and junction capacitance at the front end's frequency.
4. **Storage bank.** HV capacitors with leakage, dielectric absorption and
   breakdown limits; bank voltage and charge are logged.
5. **Staged Royer step-down.** A chain of Royer converters (self-
   oscillating push-pull with a saturable transformer), each fed from the
   previous stage's capacitor and feeding the next through a rectifier and
   capacitor, or ending in an AC output. Each stage lowers the voltage by
   its turns ratio; its frequency follows the Royer relation
   `f = V_in / (4 N_p B_sat A_e)`, so frequency is set stage by stage by
   input voltage, primary turns, core saturation flux density and core
   area. Royer stages are normally built for low input voltage; high-
   voltage stages need switch devices rated for the stage voltage (or
   series-stacked stages), and the engine rejects designs that exceed a
   device's rating.
6. **Magnetic amplifier between every pair of Royer stages.** Each Royer
   stage's square-wave output passes through a self-saturating magnetic
   amplifier (saturable reactor in series with a rectifier, half- or
   full-wave; Mali pp. 34–35) before the rectifier and capacitor that feed
   the next stage. A small DC control current sets how much of each
   half-cycle the core blocks before it saturates, so the amplifier
   regulates the volt-seconds — and therefore the DC voltage — passed to
   the next stage. Because a Royer stage's frequency follows its input
   voltage (`f = V_in / (4 N_p B_sat A_e)`), each magnetic amplifier also
   sets the next stage's frequency: the chain of control currents tunes the
   step-down, stage by stage. A Royer square wave and a square-loop core
   suit this well.
   - **Control and feedback (Mali pp. 50–55).** Power gain is output power
     over control power, and cascaded gains multiply. Negative feedback
     improves linearity and shortens response time at the cost of gain;
     positive feedback raises gain and lengthens response, and with a
     square loop gives bistable (two-level) operation. Response time is
     counted in supply cycles, so the higher-frequency stages near the
     front of the chain respond fastest. Typical figure of merit (power
     gain per cycle): 300–400 for low-power and 100–200 for high-power
     amplifiers.
   - **Energy.** A magnetic amplifier controls power; it does not add any.
     Its output power comes from the AC supplied by the Royer stage before
     it, and the control winding supplies only the small control power.
     Mali gives overall efficiency usually above 50% and up to 80–90% in
     larger sizes (p. 63), so every amplifier in the chain is a loss as
     well as a regulator. The ledger books control power, core and winding
     loss, and rectifier loss for each one, and the GA weighs regulation
     against the added loss.
   - **Pulsed drive and core reset.** The first amplifier is driven by
     the gap's pulse train rather than a symmetric square wave. With
     unidirectional pulses the core does not reset itself each cycle, so
     the amplifier needs an explicit reset (control winding bias or a
     reset circuit); the model checks that the core returns below
     saturation between pulses and flags designs where it ratchets into
     permanent saturation.

**Models needed.** `device/core.go`: saturable magnetic core with
hysteresis (Jiles–Atherton), required for the Royer transformer to saturate
and set the frequency. `device/royer.go`: the Royer stage as a composite
(two switches, centre-tapped saturable transformer, feedback winding).
`device/magamp.go`: self-saturating magnetic amplifier (half- and
full-wave) built on `device/core.go`, with control winding, load windings,
series rectifiers and optional feedback winding. `fluid/mist.go`: a
two-phase gap medium (air or other gas plus water droplets) with droplet
size distribution, liquid water content, droplet conductivity (H₂O or D₂O,
from `fluid/media.go`) and evaporation, setting the arc model's (§3g)
breakdown and recovery parameters from published measurements.

**Validation (must pass before use).**

- Royer stage frequency matches `V_in/(4 N_p B_sat A_e)` within the core's
  tolerance, and its waveform and efficiency match a published design.
- Diode ladder output voltage, ripple and droop match the analytic
  Cockcroft–Walton values.
- Charge conservation: charge drawn from ground equals charge on the bank
  plus leakage.
- Spark gap in mist: with a clean fog of ≈ 30 µm droplets the model
  reproduces Liu et al. (2018) — AC discharge voltage of 5, 10 and 15 cm
  rod–plane gaps of 29.4, 43.3 and 61.2 kV against 21.7, 35.9 and 53.2 kV
  dry, and 26.7, 39.3 and 58.7 kV at 4000 µS/cm; with water streams it
  reproduces the downward trends of Yuan et al. (2014); recovery in dry
  air and nitrogen follows the deionization (tens of µs) and thermal (ms)
  phases of Bhoomireddy (2016). With the mist off, the gap reduces to the
  dry-air §3g case. Energy lost to evaporation and dissociation is
  booked separately in the ledger.
- Bifilar storage: the coil's stored magnetic and electric energy, summed
  over a pulse, matches the energy delivered by the gap minus losses.
- Magnetic amplifier: transfer curve (output versus control current) and
  power gain from the core model match a published characteristic; the
  effects of negative and positive feedback on gain, linearity and
  response time have the directions Mali describes (pp. 54–55); response
  time scales with the supply period; efficiency lies in the published
  range for the core size.
- Energy closure: energy at the final output plus all losses equals the
  energy delivered by the declared sources, per source.

**Genetic evolution.** Genome: front-end type and parameters (for the
spark gaps: number of gaps, gap length, electrode material, mist droplet
size and liquid water content, tank capacitor; topology A or B), bifilar coil geometry and
tuning, ladder stage
count, bank size, number of Royer stages and each stage's turns, core and
capacitor, and for each magnetic amplifier its core, turns, control
current and feedback. Fitness: power at the low-frequency output per joule from the
declared source (end-to-end efficiency), output frequency accuracy, and
bank voltage within component ratings.

**References.** P. Mali, *Magnetic Amplifiers: Principles and
Applications*, Rider, 1960 (two scans in `reference/`; see §6). To add to `reference/`: G. H. Royer, "A switching
transistor DC to AC converter having an output frequency proportional to
the DC input voltage," *AIEE Trans.*, 1955; D. C. Jiles and D. L.
Atherton, "Theory of ferromagnetic hysteresis," *J. Magn. Magn. Mater.* 61,
1986; J. D. Cockcroft and E. T. S. Walton, *Proc. R. Soc. A* 136, 1932.

### 3k. Scalar and longitudinal wave models (`em/scalar/`)

"Scalar wave" is used for three different things. The engine keeps them in
three tiers, so established results, extended theories and frontier claims
are each handled at the right level of evidence — and none is left out
until it has been tested (§1, Hypotheses under test).

**Tier 1 — longitudinal fields inside standard electromagnetism
(established; modeled as standard physics).**

- Near fields and longitudinal E components of small sources and antennas.
- TM waveguide and cavity modes (including the longitudinal fields of
  Kron's 1945 cavity network, §3a.1).
- Surface waves guided by wires and by the ground — Sommerfeld, Zenneck and
  Goubau waves — the standard physics closest to Tesla's earth-conduction
  ideas.

Each is implemented in `em/fdtd.go` / `em/tlm.go` with an analytic validation case (dipole near field, TM₁₁₀
cavity mode, Sommerfeld–Goubau line attenuation). Reference for longitudinal waves within standard
theory: arXiv:1606.01738.

**Tier 2 — extended electrodynamics with a physical scalar field
(hypotheses with published equations).** These theories keep the scalar
term that the Lorenz gauge removes and make it physical:

- **Jack (2003)** — the temporal field T (H4, `em/temporal.go`).
- **van Vlaenderen (2003)**, arXiv:physics/0305098 — the scalar field S
  with generalized Gauss and Ampère laws; cites Tesla's work.
- **Hively (2012)**, "Toward a more complete electrodynamic theory," *Int.
  J. Signals and Imaging Systems Eng.* 5(1), 3–10, and later comparisons of
  extended theories.

All tier-2 models implement one interface in `em/scalar/` so they are
swappable and compared on the same runs, and each must pass the same gates
before its results are used:

1. reduces exactly to Maxwell when its scalar field is zero;
2. charge conservation checked explicitly in every run;
3. a stated energy theorem that closes in its own ledger column;
4. a cited source for every equation, and declared couplings (free
   parameters stated per run, or fixed from measured material data).

**Tier 3 — frontier claims (in the register; modeled where a measurable
prediction can be stated).** Work that pushes past accepted science is kept
and tested, not dismissed. Each claim is entered with its own predicted
signature, the standard prediction beside it, and published critiques
listed alongside its sources:

- **T. E. Bearden — "scalar EM" and the Motionless Electromagnetic
  Generator** (US 6,362,718, 2002, Patrick, Bearden et al.). Bench target: a
  replica built and run with the full ledger on (§3c), measuring every input
  (supply, ambient field, thermal) and every output; the claim is a
  coefficient of performance above 1, the standard prediction below 1.
- **K. Meyl — longitudinal (vortex) scalar waves** (*Scalar Waves*, INDEL,
  2003; demonstration transmission kits). Bench targets: propagation delay
  between transmitter and receiver (Meyl predicts speeds differing from c),
  received power versus distance (standard near-field and radiative
  fall-off versus Meyl's prediction), and resonance shift with distance.
  Critiques to record with it (e.g. G. W. Bruhn's published refutations).
- **Monstein & Wesley (2002)** — "Observation of scalar longitudinal
  electrodynamic waves," *Europhys. Lett.* 59, 514–520: 433.59 MHz ball
  antennas, 4.0 m and 4.4 m above ground, longitudinality tested with a
  cubic array of nine half-wave wires (blocking when parallel to the
  propagation direction only). Later Comment and Reply in EPL. The
  simulation reproduces the full setup, including the ground image source,
  under standard theory and every tier-2 model, so a bench repeat
  discriminates between them.
- **Electrogravitic thrust (Biefeld–Brown asymmetric capacitors).** The
  standard explanation is ion wind. Bench target: thrust in air versus in
  vacuum (or at reduced pressure) for the same capacitor
  and voltage; the standard prediction is thrust falling with the ion
  population, the electrogravitic claim is thrust persisting in vacuum.

**Shared outputs.** Every tier-2 and tier-3 run reports, at each probe pair:
longitudinal and transverse E, B, the E/B ratio, propagation delay, and the
ledger by column, with the standard prediction alongside — the same
deviation checks as the §3g pulse experiment.

**Keeping up.** A monthly literature watch (scheduled task "KronSpice
scalar-wave literature watch", also run on request) searches for new work on
these models and writes a dated report to `doc/literature/` with proposed
register updates. It never edits this file; changes to the physics are made
here only after review.

### 3l. Gravity and fast rotation (`spacetime/`)

Kept in scope for craft that travel between stars and for machines that
spin very fast (user decision, 2026-10-27). Three different "curvatures"
are kept apart, because they are often confused:

1. **Gravitational curvature of spacetime.** Schwarzschild exterior of a
   non-rotating mass: clock rates `dτ/dt = √(1 − r_s/r − v²/c²)`, redshift
   between clocks, tidal stretching, horizon and photon-sphere limits
   (typed errors). Maxwell's equations in a static gravitational field are
   flat-space Maxwell in a medium with `ε_r = μ_r = n(r)` (Plebanski, Phys.
   Rev. 118, 1396, 1960; Landau & Lifshitz §90), so `em/` field solvers take
   gravity as an index; weak field `n = 1 + (1+γ)GM/(rc²)` with the PPN
   parameter γ (Cassini bound γ − 1 = (2.1 ± 2.3)×10⁻⁵), exact isotropic
   Schwarzschild index for strong fields. Frame dragging (Kerr,
   Lense–Thirring) is not modeled yet.
2. **Rotating frames in flat spacetime** (rotors, homopolar disks,
   flywheels, a spinning hull). Lorentz transformations hold only between
   inertial frames; a rotor is described by the Langevin metric. Clocks on
   it run at `√(1 − ω²r²/c²)` (the clock hypothesis, tested to 10⁻³ at
   10¹⁸ g in the CERN muon ring); light around it shows the Sagnac delay
   `4Aω/c²` and a ring resonator on it splits by `4Aω/(λP)`; a probe riding
   it sees the fields of the local co-moving inertial frame (`em/frames.go`
   boost at `ω×r`, Mashhoon's locality hypothesis). A rotor's rim speed is
   bounded by its material: a thin hoop bursts at `√(σ/ρ)` (≈ 1,880 m/s for
   T1000G carbon fibre), where `(v/c)² ≈ 4×10⁻¹¹`; every rotor report gives
   `(v/c)²`, the rim clock lag and the burst margin, and a design past its
   burst speed is rejected. The rotor energy-shift anomaly is register
   entry H17: the standard k = 1/2 and a declared k are computed together.
3. **Kron's curved configuration space** (`kron/machine.go`): windings and
   rotor angle as coordinates, inductances `L(θ)` as the metric, speed
   voltages and torque from Christoffel symbols, non-holonomic d-q frames
   (Kron, *J. Math. Phys.* 13, 103, 1934). Not spacetime curvature; the
   same mathematics, and the frame in which every machine element is
   solved.

**Validation (must pass, `spacetime/spacetime_test.go`).** GPS satellite
clock net +38.6 µs/day against the rotating Earth (Ashby, Living Rev.
Relativ. 6, 1, 2003; model 38.54 without J2); gravitational part 45.7
µs/day; Pound–Rebka–Snider `gh/c²` over 22.5 m (2.45×10⁻¹⁵); Shapiro delay
closed form against the index path integral (and the exact isotropic
index); solar-limb deflection 1.751″ (VLBI: (1+γ)/2 = 0.99992 ± 0.00023);
`r = 3r_s` static rate `√(2/3)`; Earth-surface tidal gradient
3.08×10⁻⁶ s⁻²; rotor shift with k = 1/2 equal to exact rim time dilation;
exact Sagnac equal to `4Aω/c²`; T1000G burst speed.

**Rotation and curvature together (`spacetime/machine.go`, test
`machine_test.go`).** By the equivalence principle a rotor's hub runs at
the Schwarzschild rate for its position and speed, and its rim is special
relativity in the hub's local inertial frame (exact form checked against
the first-order weak-field metric). Its spin axis precesses by Thomas +
geodetic + Lense–Thirring (Schiff, *PRL* 4, 215, 1960); a loop riding it
sees the Sagnac delay relative to the local gyroscopes. Validated against
Gravity Probe B (Everitt et al., *PRL* 106, 221101, 2011): model 6621
mas/yr geodetic (GR 6606.1, measured 6601.8 ± 18.3; spherical Earth, no
J2) and 39.1 mas/yr frame dragging seen against the guide star IM Pegasi
(GR 39.2, measured 37.2 ± 7.2). `MachineSelfField` gives the reverse — the
curvature and frame dragging a machine's own mass and spin produce
(`examples/machine_spacetime`: a 28 kg rotor at 1 m drags frames ~10⁻¹⁰ of
Earth's own frame dragging).

**Rules.** Gravity changes circuits only through what reaches them:
in its own rest frame a circuit obeys the same Maxwell equations
everywhere (equivalence principle). Results seen by a distant observer
are converted by the clock rates; tidal stress beyond a declared
mechanical limit, a radius inside the photon sphere for an orbit, or at
the horizon is an error, never an extrapolation. No energy enters a
circuit from spacetime curvature by itself; the gravitational / frame
channel carries no power (§3.1).

---

## 4. Genetic algorithm layer (`evo/`)

### Encodings

1. **Value genome (start here).** Fixed topology; genes are component values
   on a log scale, optionally snapped to E12/E24/E96 series. This is the most
   reliable starting point.
2. **Topology genome.** Variable-length list of elements with node indices, or
   a Koza-style developmental program that grows a circuit from an embryo.
   Every decoded circuit passes `circuit.Validate()` before simulation.

### Fitness evaluation rules

- Fitness is computed from in-memory `output.Result` plus `analysis/measure`
  helpers (gain at f, −3 dB bandwidth, rise time, overshoot, THD, power).
- **Invalid or failed circuits get a large penalty fitness, never a crash.**
  Map each error class (`ErrSingularMatrix`, `ErrNoConvergence`,
  `ErrFloatingNode`, timeout) to a penalty so the GA can learn from them.
- Every evaluation runs with a `context.Context` deadline so one bad circuit
  cannot stall a generation.
- Optional parsimony pressure: a small penalty per component to prevent bloat.
- Multi-objective fitness (e.g. gain vs. power) should support weighted sums
  first and NSGA-II later.

### Parallelism and determinism

- A worker pool sized to `runtime.NumCPU()` evaluates genomes. Each worker owns
  its own simulator instance; nothing mutable is shared.
- All randomness comes from a `*rand.Rand` seeded from the experiment config,
  with child seeds derived per generation and individual. Never use the global
  `math/rand` source.
- Results are collected and sorted by individual index before selection, so
  scheduling order never affects the outcome.

### As implemented (M5)

- Value genome: genes target an element value or a `.model` parameter; log
  scale; optional E6/E12/E24/E48/E96 snapping; sign mode `positive`,
  `negative` or `either` (the sign is then its own gene, so evolution can
  discover negative resistance).
- Operators: tournament selection, BLX-0.5 crossover on log magnitudes,
  Gaussian mutation in decades, sign flips, elitism.
- Objectives (weighted sum): cutoff, gain at f, operating-point voltage,
  transient mean and ripple, open-circuit voltage, short-circuit current,
  model-supplied power limit; each analysis runs once per candidate.
- Penalty classes: no-measurement < no-convergence < singular < invalid <
  ledger imbalance < timeout < other < panic (caught and reported as a bug).
- Determinism: one PCG stream per generation derived from the seed, results
  stored by index; a checkpointed run resumes bit-identically.

### As implemented (M8)

- Topology genome (`evo/encode_topology.go`): a direct netlist encoding —
  a variable-length list of R/L/C parts (kind, two nodes, log value,
  E-series per kind) grown into an embryo netlist that holds only the
  fixture. Parts may use ground, the embryo's declared terminals and a pool
  of internal nodes. Dangling parts and self-loops are pruned at decode
  (they carry no current); what remains invalid is penalized.
- Operators: change value, add, delete, rewire, change kind, series split,
  parallel add; one-point crossover on the part lists.
- Objectives: `mask` (frequency bands with gain limits; error 0 inside,
  ≥ 0.1 outside, so "fitness < 0.1" means spec met) and `parts`
  (parsimony, weight per part).

### Experiment outputs

Each run writes to `runs/<timestamp>-<name>/`: config copy, per-generation
stats CSV (best / mean / worst fitness, diversity), best netlist per
generation, and a final checkpoint. Best netlists must be valid SPICE that
ngspice, LTspice and QSpice can open, so results can be checked in other tools.

---

## 5. Libraries

Keep dependencies few and pure Go (no cgo) so the project cross-compiles
everywhere, including Windows next to QSpice.

| Purpose | Library | Notes |
|---|---|---|
| Dense linear algebra, LU, complex | `gonum.org/v1/gonum` (`mat`, `floats`) | BSD-3. Core dependency. |
| Sparse LU | native (`solver/sparse.go`, M6) | Gilbert–Peierls sparse LU written from the published algorithm; `github.com/rwl/lufact` (BSD-3) was the fallback plan and was not needed. |
| Tensors (Kron) | `gonum.org/v1/gonum/mat` | Kron's network tensors are rank-2, so gonum matrices with typed frames are enough. Consider `gorgonia.org/tensor` (Apache-2.0) only if higher-rank tensors are ever needed. |
| Quaternions | native (`hyper/quat`, M9d) | A few dozen lines; written here to keep the project free of third-party dependencies. gonum's `num/quat` has the same algebra. |
| Geometric algebra | native (`hyper/ga`, M9d) | No mature pure-Go Clifford algebra library. Checked against Pauli and Dirac matrix representations; ganja.js (MIT) and the Python `clifford` / `kingdon` packages can add test vectors (numbers only; check licenses). |
| FFT (GA power theory) | `gonum.org/v1/gonum/dsp/fourier` | For harmonic decomposition of transient waveforms. |
| Plotting (examples, reports) | `gonum.org/v1/plot` | Used only in `examples/` and `cmd/`, never in core packages. |
| GA operators (optional) | `github.com/MaxHalford/eaopt` | MIT. GA, differential evolution, PSO, parallel evaluation. Our `evo` engine is primary; eaopt can be used for quick experiments or comparison. Last tagged release is from 2019, so do not make it a hard dependency of core code. |
| Testing | stdlib `testing`, `github.com/google/go-cmp` | No assertion frameworks. |
| CLI | stdlib `flag` | Add a CLI framework only if commands multiply. |
| Config | stdlib `encoding/json` (+ `gopkg.in/yaml.v3` if YAML is wanted) | |

Before adding any new dependency, explain why in the PR/commit message and
prefer the standard library.

---

## 6. References (study, don't copy)

Read these for algorithms and expected behavior. **Do not copy code from GPL
projects into this repository.** Check each project's license before reusing
any code at all; when in doubt, reimplement from the published equations.

**Licensing of this project.** kronspice is dual-licensed, GPL-3.0-only or
commercial (`LICENSING.md`). Every new file starts with
`// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial`
and the copyright line. Only original code, or code under a permissive
license compatible with both options (MIT, BSD, Apache-2.0, notice kept),
may be added; copyleft code (GPL, LGPL, AGPL, MPL…) can never be included,
because it could not be offered under the commercial license.

**Simulators**
- **Berkeley SPICE3f5** — the original; BSD-licensed. Canonical device
  equations, limiting functions, convergence strategy.
- **ngspice** — maintained descendant of SPICE3f5 (plus CIDER and XSPICE).
  Our golden reference for results; check per-file licenses.
- **Xyce** (Sandia) — modern C++ design, good model of clean device and
  solver separation. GPLv3: reference only.
- **ahkab** — small, readable SPICE-like simulator in Python. Good for seeing
  the whole pipeline at once. GPL: reference only.
- **PySpice** — Python interface to ngspice and Xyce; useful for netlist
  parsing ideas. GPLv3: reference only.
- **spicelib / qspice (Python)** — for reading QSpice/LTspice raw files and
  netlists, and later for cross-checking results against QSpice on Windows.

**Books and papers**
- L. W. Nagel, *SPICE2: A Computer Program to Simulate Semiconductor
  Circuits*, UC Berkeley ERL Memo M520, 1975.
- C.-W. Ho, A. Ruehli, P. Brennan, "The modified nodal approach to network
  analysis," *IEEE Trans. Circuits and Systems*, 1975.
- J. Vlach, K. Singhal, *Computer Methods for Circuit Analysis and Design*.
- F. N. Najm, *Circuit Simulation*.
- K. Kundert, *The Designer's Guide to SPICE and Spectre* (convergence and
  timestep behavior).
- J. R. Gilbert, T. Peierls, "Sparse Partial Pivoting in Time Proportional to
  Arithmetic Operations," *SIAM J. Sci. Stat. Comput.*, 1988.
- G. Kron, *Tensor Analysis of Networks*, Wiley, 1939 — primitive networks,
  connection tensors, transformation laws.
- G. Kron, *Diakoptics: The Piecewise Solution of Large-Scale Systems*,
  Macdonald, 1963 — tearing and interconnection.
- G. Kron, "Tearing and interconnecting as a form of transformation,"
  *Quarterly of Applied Mathematics* 13, 147–159, 1955 — tearing stated as a
  transformation of reference frame, the mathematical core of diakoptics
  (§3a). Copy in `reference/` (`S0033-569X-1955-73431-7.pdf`).
- H. H. Happ, *Diakoptics and Networks*, Academic Press, 1971 — a clearer,
  matrix-oriented presentation of Kron's methods: orthogonal networks with
  open and closed paths (Ch. III–IV), ignoring or eliminating open or
  closed paths (Ch. XI), diakoptics (Ch. XII). Copy in `reference/`
  (`happ_1971_diakoptics_and_networks.pdf`); its Appendix I worked example
  is reproduced in `kron/happ_test.go`; chapter map in
  `doc/happ-1971-notes.md`.
- F. Kouřil, K. Vrba, *Non-linear and Parametric Circuits: Principles,
  Theory and Applications*, Ellis Horwood / John Wiley / SNTL, Chichester,
  New York and Prague, 1988, 413 pp. (Ellis Horwood Series in Electrical
  and Electronic Engineering; ISBN 0-85312-606-2, 0-470-20014-6) —
  reference for the nonlinear and parametric (time-varying) circuit work:
  M9c tearing of nonlinear circuits and the parametric elements to follow.
  Photographed pages in `reference/kouril-vrba/` (§4.4–4.5 pp. 128–143,
  §6.4.4–6.4.5 pp. 222–228, §10.2 pp. 370–393, references pp. 394–397);
  transcription and the tests built on it in
  `doc/kouril-vrba-verification.md` (Manley–Rowe and Pantell power
  relations; parametric self-excitation, amplitude and amplifiers; the
  non-linear tank, where eq. 10.37 overstates the frequency shift).
- J. M. Manley, H. E. Rowe, "Some General Properties of Nonlinear
  Elements — Part I," *Proc. IRE* 44(7), 904–913, 1956 — power relations
  of non-linear storage elements (K&V [43]).
- R. H. Pantell, "General Power Relationship for Positive and Negative
  Nonlinear Resistive Elements," *Proc. IRE* 46(12), 1910–1913, 1958
  (K&V [47]).
- W. H. Louisell, *Coupled Mode and Parametric Electronics*, Wiley, 1960;
  D. G. Tucker, *Circuits with Periodically-Varying Parameters*,
  Macdonald, 1964 — parametric amplifiers and periodically varying
  circuits (K&V [39], [58]).
- J. P. Roth, "An application of algebraic topology to numerical analysis: on
  the existence of a solution to the network problem," *PNAS* 41(7),
  518–521, 1955 — the algebraic-topology reading of Kron. Copy in
  `reference/` (`roth-1955-an-application-of-algebraic-topology-to-numerical-
  analysis-on-the-existence-of-a-solution-to-the-network.pdf`).
- F. H. Branin Jr., "The algebraic-topological basis for network analogies and
  the vector calculus," 1966.
- E. Tonti, "Finite formulation of the electromagnetic field," *Progress In
  Electromagnetics Research* (PIER) 32, 1–44, 2001 — Maxwell's equations
  written directly on a cell complex and its dual (global variables on
  points, lines, surfaces, volumes; inner and outer orientation; a time cell
  complex), "the natural extension of network theory to the electromagnetic
  field": the modern form of Kron's 1944 field network and of the Roth/Branin
  chains and cochains (§3a.1). Copy in `reference/`
  (`tonti_2001_pier32_finite_formulation.pdf`, from the Wayback Machine copy
  of the jpier.org PDF). `reference/Electrom ICS.pdf` is the 7-page short
  version dated 16 October 2000 (International Compumag Society Newsletter
  8(1), 2001), not the PIER article.
- G. Kron, "Equivalent circuit of the field equations of Maxwell — I,"
  *Proc. IRE* 32, 289–299, 1944 — field regions as networks (§3a.1).
- G. Kron, "Numerical solution of ordinary and partial differential
  equations by means of equivalent circuits," *J. Appl. Phys.* 16, 172–186,
  1945 — field networks for Maxwell's equations (waveguide/cavity example)
  and the elastic field; boundary-value and initial-value solution methods.
  Copy in `reference/`.
- G. Kron, "Electric circuit models of the Schrödinger equation," *Phys.
  Rev.* 67, 39–43, 1945 — eigenvalue networks and boundary terminations.
  Copy in `reference/`.
- G. Kron, "Multidimensional curve-fitting with self-organizing automata,"
  *J. Math. Anal. Appl.* 5, 46–69, 1962 — the excited polyhedron and its
  dual as a least-squares estimator of functions of k variables and their
  divided differences; Kron's own placement of the polyhedral program in the
  exterior-form calculus of Cartan, de Rham, Hodge and Whitney (§3a.1,
  HISTORICAL.md). Copy in `reference/` (`PII_ 0022-247X(62)90005-7.pdf`);
  Figs. 1–3, the algorithm and Table I transcribed from the page images in
  `doc/kron-figure-catalogue.md`.
- P. Russer, "The Transmission Line Matrix Method," TU München (chapter,
  21 pp.) — 2-D and 3-D TLM, boundaries, discrete Green's functions, skin
  effect, nonlinear lumped elements in the mesh; cites Kron 1944 as the
  origin of the method. Implementation route for field regions (§3a.1
  item 8, §3g). Copy in `reference/`.
- P. B. Johns, R. L. Beurle, "Numerical solution of 2-dimensional
  scattering problems using a transmission-line matrix," *Proc. IEE* 118,
  1203–1208, 1971 — the original TLM paper.
- G. Kron, *Tensors for Circuits*, Dover, 1959 (first published 1942 as *A
  Short Course in Tensor Analysis for Electrical Engineers*) — compact
  statement of the transformation laws and equations of constraint. Copy in
  `reference/`.
- J. P. Roth, "An application of algebraic topology: Kron's method of
  tearing," *Quarterly of Applied Mathematics* 17, 1959 — ohmicness as the
  existence condition, power invariance under tearing. Copy in
  `reference/`.
- P. W. Aitchison, "Diakoptics or tearing — a mathematical approach,"
  *Quarterly of Applied Mathematics* 41, 265–272, 1983 — explicit tearing
  algorithm and non-singularity proof (§3a). Copy in `reference/`.
- H. Nicholson, "System Concepts in Kron's Polyhedron Model and the
  Scattering Problem," ACSE Research Report 26, University of Sheffield,
  1974 (journal version *Int. J. Control* 20(4), 529–555) — Kron's
  polyhedron model, p-networks and the scattering problem; the one
  continuation of Kron's multidimensional networks we hold. Copy in
  `reference/` (`acse research report 26.pdf`, scanned, no text layer;
  read from the page images, its equations and figures transcribed in
  `doc/kron-figure-catalogue.md`).
  Its companion, ACSE Report 14 (1973; *Proc. IEE* 120, 1545–1557), is open
  at White Rose eprints 86081.
- N. V. Balasubramanian, *Differential structures in electromagnetic field
  and circuit theory*, PhD thesis, University of Liverpool, 1968 (supervisor
  J. W. Lynn) — Kron's orthogonal networks and his 1944 network model of
  Maxwell's equations in exterior differential forms (Ch. 2–5), applied to
  waveguides (Ch. 6) and to Park's and Kron's machine frames, hunting and the
  affine connection (Ch. 7–9). Its introduction (pp. 6–22) quotes at length
  from six of Kron's 1958–65 polyhedron papers that are otherwise
  unobtainable, so it is the project's text for them (HISTORICAL.md). Bridge
  between `kron/` and `em/` (§3a.1 item 8) and background for
  `kron/machine.go`. Copy in `reference/` (`642931.pdf`); its field-network
  figures (5–17, including Fig. 9, Kron's 1959 space-filter diagram) and
  the waveguide results are transcribed in `doc/kron-figure-catalogue.md`.
- D. V. Steward, "Partitioning and tearing systems of equations," *SIAM J.
  Numer. Anal. Ser. B* 2, 345–365, 1965 — partitioning algorithms (§3a).
- A. Brameller, M. N. John, M. R. Scott, *Practical Diakoptics for Electrical
  Networks*, 1969.
- O. Maurice, A. Reineix, P. Durand, F. Dubois, "Kron's method and cell
  complexes for magnetomotive and electromotive forces," *IAENG Int. J.
  Applied Mathematics* 44(4), 2014 — chains/cochains, radiated coupling,
  validation against FDTD and measurement.
- P. Durand, O. Maurice, A. Reineix, "Generalized interaction principle
  implemented in the Kron's method," *Proc. World Congress on Engineering*,
  2013 — chords and affine source terms.
- *Tensorial Analysis of Networks (TAN) Modelling for PCB Signal Integrity and
  EMC Analysis*, IET, 2020 — TAN applied to boards and EMC.
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley, 1951 — the
  primitive machine, reference frames of rotating machines, and the
  Epilogue (pp. 253–259): the dynamical equations in holonomic and
  non-holonomic, Riemannian and non-Riemannian frames, the
  Christoffel-symbol form (eqs. 3–4) and the speed-voltage tensor G
  (eq. 9). Basis of `device/machine.go` and `kron/machine.go`. Noted for
  later work (doc/M9c1.md, "Not done here"): eq. 11 (p. 259) rewrites the
  machine equation with the field tensor F_αβ and, with e and R zero, is
  the equation of motion of a charged particle in an electromagnetic field,
  the link between the machine's G·i·ω and the Lorentz force of
  `em/lorentz.go`; frames III and IV (eqs. 12–13, pp. 259–260) cover
  interconnected machines on a freely rotating frame; part 2 (pp.
  261–265) reduces the dynamical equations to equivalent circuits through
  the absolute time derivative and the absolute-frequency tensor n,
  `Z = n⁻¹R + jX` (eqs. 14–34), the route to an AC model of the N element
  at speed. Two copies
  in `reference/`: `Equivalent Circuits of Electric Machinery Kron.pdf`, the
  complete HathiTrust scan (mdp.39015000477441, Google-digitized, public
  domain; 305 pages with a text layer; printed page n is PDF page n + 23,
  so the Epilogue starts at PDF page 276), to be preferred; and
  `374375674-Gabriel-Kron-Equivalent-Circuits-of-Electric-Machinery-pdf.pdf`,
  a 146-page image-only scan without a text layer.
- G. Kron, "Non-Riemannian Dynamics of Rotating Electrical Machinery,"
  *Journal of Mathematics and Physics* 13, 103–194, 1934 — the
  non-holonomic form of Lagrange's equations for machines (cited by the
  1951 Epilogue).
- G. Kron, *The Application of Tensors to the Analysis of Rotating Electrical
  Machinery* — Kron's rotating reference frames, the basis for the quaternion
  polyphase work.
- W. R. Hamilton, *Elements of Quaternions*.
- W. Thomson (Lord Kelvin), *Mathematical and Physical Papers*, Cambridge
  University Press — Vol. I (1882; Art. XXVII, mechanical representation of
  electric, magnetic and galvanic forces, 1847) and Vol. V (1911; Art. 175,
  "Magnetism and Molecular Rotation", 1899, pp. 354–358: a charged ring
  filled with oppositely charged fluid as a gyrostatic molecule, discarded
  by Kelvin after Zeeman's splitting, p. 357; "Æpinus Atomized", 1902).
  Copies in `reference/`; notes in `doc/kelvin-references.md`; the 1899
  ring is simulated in `examples/kelvin_ring` (induction ledger; Zeeman
  test: gyrostat as written broadens, full ring gives sharp lines with the
  wrong polarisation). Vol. IV (1910, Hydrodynamics and General
  Dynamics; Internet Archive `in.ernet.dli.2015.212310`, not yet in
  `reference/`): Art. 1, "On Vortex Atoms" (1867), atoms as vortex rings
  in a perfect liquid, with flow ↔ magnetic field of a ring current and
  impulse ↔ magnetic moment (pp. 9–10); Art. 3, ring speed (p. 67).
  Simulated in `examples/vortex_ring` through that analogy with the M9i
  inductance code.
- P. Tewari, "Structural relation between the vacuum space and the
  electron," *Physics Essays* 31(1), 108–129, 2018 — the space-vortex
  electron. **Hypothesis source** for H13; notes and the test bench in
  `doc/tewari-references.md` and `examples/space_vortex`. Copy in
  `reference/tewari_2018_physics_essays.pdf`. P. Tewari, *Universal
  Principles of Space and Matter*, Grain and Chaff, Toronto, 2017 (1st ed.
  2009) — the electron (ch. 2, 5) and the Space Power Generator test
  (§6.7, pp. 149–161); copy in `reference/tewari_upsm_2017.pdf`. Still to
  obtain: *Discovering Universal Reality* (2007).
- F. D. Amador, H. Eckardt, "Simulation of a Parametric Resonance
  Circuit," A.I.A.S./UPITEC, 30 March 2012 — series R, L, C(t) with
  C(t) = C(1 + a0 sin(ω1t + φ)), eqs. (19)–(38). **Hypothesis source** for
  H16; analysed in `examples/parametric_lcr`; copy in
  `reference/amador_eckardt_2012_parametric_lcr.pdf`. Standard counterpart:
  L. Mandelstam, N. Papalexi, *Z. Phys.* 73, 223, 1932 (parametric
  excitation); Kouřil & Vrba §10.2.2.
- H. Eckardt, "Simulation of resonances in the homopolar generator,"
  A.I.A.S./UPITEC, 13 February 2012 — self-excited twin-disk homopolar,
- B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Rev. Mod. Phys.*
  21, 535–540 (1949) — Kron's machine theory for physicists: Lagrange–
  Maxwell equations as geodesics with the inductances as metric (eqs. 1–8),
  non-holonomic frames, the torsion Ω (eq. 19) and eq. 5′, and the
  single-phase squirrel-cage motor with stationary rotor coils (eqs. 22–43).
  Basis of `kron/hoffmann.go` and `examples/hoffmann_motor`. Copy in
  `reference/hoffmann-1949-krons-non-riemannian-electrodynamics.pdf`
  (6-page image scan).
- H. E. Puthoff, "Ground state of hydrogen as a zero-point-fluctuation-
  determined state," *Phys. Rev. D* 35, 3266–3269 (1987) — the SED balance
  of absorbed zero-point and radiated Larmor power, m·ω·r² = ħ (H24,
  `examples/zpf_hydrogen`). Later full SED simulations: D. C. Cole &
  Y. Zou, *Phys. Lett. A* (2003), doi:10.1016/j.physleta.2003.08.022;
  T. M. Nieuwenhuizen & M. T. P. Liska, *Phys. Scr.* T165, 014006 (2015).
  Copy in `reference/puthoff-1987-hydrogen-ground-state-zpf.pdf`.
- Th. Kaluza, "Zum Unitätsproblem der Physik," *Sitzungsber. Preuss. Akad.
  Wiss. Berlin (Math. Phys.)* 1921, 966–972 — the five-dimensional metric
  framed by the electromagnetic potential, the "cylinder condition" Hoffmann
  later applies to Kron's machines; O. Klein, *Z. Phys.* 37, 895 (1926), the
  compact fifth dimension. `examples/kaluza_kron`. Copy (V. T. Toth's
  translation) in `reference/kaluza-1921-unification-problem-toth.pdf`.
- R. L. Forward, "Extracting electrical energy from the vacuum by cohesion
  of charged foliated conductors," *Phys. Rev. B* 30, 1700 (1984) (H25,
  `examples/forward_casimir`); D. C. Cole & H. E. Puthoff, "Extracting
  energy and heat from the vacuum," *Phys. Rev. E* 48, 1562 (1993);
  A. Lambrecht & S. Reynaud, *Eur. Phys. J. D* 8, 309 (2000), plasma
  wavelengths for the Lifshitz force. Copy in
  `reference/forward-1984-vacuum-fluctuation-battery.pdf`.
- A. Einstein, "New possibility for a unified field theory of gravitation
  and electricity," *Sitzungsber. Preuss. Akad. Wiss.* 1928, 224–227
  (A. Unzicker's translation) — distant parallelism: 4-bein, torsion Λ,
  φ_μ = Λ^α_μα as electromagnetic potential (`examples/einstein_teleparallel`);
  T. Sauer, *Historia Mathematica* 33, 399 (2006), its history 1928–31.
  Copy in `reference/einstein-1928-distant-parallelism-unzicker.pdf`.
- J. Nachamkin, *Force-Free Time-Harmonic Plasmoids*, PL-TR-92-3044
  (Phillips Laboratory, 1992), DTIC AD-A257 765 (H26,
  `examples/plasmoid_rings`). Copy in
  `reference/nachamkin-1992-force-free-plasmoids.pdf`.
- K. R. Shoulders, *EV: A Tale of Discovery* (Jupiter Technologies, Austin,
  1987) — electron clusters ("EVs") from fast cathode discharges (H27,
  `examples/shoulders_ev`); H. E. Puthoff, M. A. Piestrup, "Charge
  confinement by Casimir forces," arXiv:physics/0408114; G. A. Mesyats,
  *Plasma Phys. Control. Fusion* 47, A109 (2005), ectons. Copy in
  `reference/shoulders-1987-ev-tale-of-discovery.pdf`.
  eqs. (1)–(12), parameters p. 4. **Hypothesis source** for H15; analysed
  in `examples/eckardt_homopolar`; copy in
  `reference/eckardt_2012_homopolar_resonances.pdf`. Standard counterpart: the self-exciting
  homopolar dynamo, E. C. Bullard, *Proc. Camb. Phil. Soc.* 51, 744, 1955.
- Ball-bearing (Huber-effect) motor — measured sources for the P element
  (`device/bearing.go`, `examples/bearing_motor`, `doc/bearing-motor.md`):
  J. Huber, *Elektrotechnik und Maschinenbau* 76(8), 169–174, 1959;
  R. A. Milroy, *J. Appl. Mech.* 34, 525, 1967; H. Gruenberg, *Am. J. Phys.*
  46, 1213, 1978; P. G. Moyssides, P. Hatzikonstantinou, *IEEE Trans.
  Magn.* 26, 1274, 1990 and 33, 4566, 1997; P. Hatzikonstantinou,
  P. G. Moyssides, *J. Phys. A* 23, 3183, 1990; D. B. Watson et al., *IEE
  Proc. A* 140, 1993 and *IEE Proc. Sci. Meas. Technol.* 146, 83, 1999;
  J. L. Choo, W. L. Soong, D. Abbott, *Proc. SPIE* 5649, 700, 2005;
  S. Marinov, *Electronics & Wireless World*, April 1989; K. T. McDonald,
  "Ball-Bearing Motor," Princeton, 2020. Coefficients κ, ω0 are fitted
  only to measured tables from these papers; none is assumed.
- J. C. Maxwell, *A Treatise on Electricity and Magnetism*, 1873 — Vol. II,
  Part IV, "General Equations of the Electromagnetic Field", stated in
  components and in quaternion notation. Primary source for `em/`.
- J. C. Maxwell, "A Dynamical Theory of the Electromagnetic Field,"
  *Phil. Trans. R. Soc.*, 1865 — the original 20 equations in components.
- P. G. Tait, *An Elementary Treatise on Quaternions* — the ∇ calculus as
  Maxwell used it.
- D. Hestenes, *Space-Time Algebra*, 1966 — `∇F = J`, the modern
  single-equation form.
- B. D. H. Tellegen, "A general network theorem, with applications,"
  *Philips Research Reports*, 1952 — circuit-level energy conservation.
- K. S. Yee, "Numerical solution of initial boundary value problems involving
  Maxwell's equations in isotropic media," *IEEE Trans. Antennas Propag.*,
  1966 — reference for any future time-domain field solver.
- D. Hestenes, G. Sobczyk, *Clifford Algebra to Geometric Calculus*, 1984.
- K. J. van Vlaenderen, "A generalisation of classical electrodynamics for
  the prediction of scalar field effects," arXiv:physics/0305098, 2003 —
  tier-2 scalar field S (§3k). Hypothesis source.
- L. M. Hively, "Toward a more complete electrodynamic theory," *Int. J.
  Signals and Imaging Systems Eng.* 5(1), 3–10, 2012 (§3k). Hypothesis
  source.
- C. Monstein, J. P. Wesley, "Observation of scalar longitudinal
  electrodynamic waves," *Europhys. Lett.* 59(4), 514–520, 2002,
  DOI 10.1209/epl/i2002-00136-9; with the later Comment and Reply in EPL
  (§3k replication target).
- "Longitudinal electromagnetic waves in the framework of standard
  classical electrodynamics," arXiv:1606.01738 (§3k tier 1).
- K. Meyl, *Scalar Waves*, INDEL, 2003 — longitudinal vortex waves (§3k
  tier 3). Critiques: G. W. Bruhn's published refutations (TU Darmstadt).
- S. L. Patrick, T. E. Bearden, J. C. Hayes, K. D. Moore, J. L. Kenny,
  "Motionless electromagnetic generator," U.S. Patent 6,362,718, 2002
  (§3k tier 3).
- P. M. Jack, "Physical Space as a Quaternion Structure, I: Maxwell
  Equations. A Brief Note," arXiv:math-ph/0307038, 2003 — left and right
  quaternion derivatives, the scalar temporal field T, links to
  thermoelectric heat. **Hypothesis source** for `em/temporal.go` (§3f); no
  experimental support. Copy in `reference/`.
- P. W. Bridgman, *The Thermodynamics of Electrical Phenomena in Metals*,
  Dover, 1961 — "working" and "driving" EMFs; the expected but unobserved
  EMF from time-varying temperature (cited by Jack).
- S. Giardino, "Quaternionic Electrodynamics," arXiv:2010.07748
  [physics.gen-ph], 2020 — extension of Maxwell with a quaternionic vector
  potential whose self-interaction generates magnetic charge and current;
  field equations, continuity, force law, Poynting theorem. **Hypothesis
  source** for `em/hed.go` (§3f); no experimental support. Copy in
  `reference/`.
- M. D. Marko, J. Schaff, "Quaternion Mathematics in Electromagnetic
  Modeling and Simulation," arXiv:2209.14765, 2022 — **cautionary example,
  not a physics source.** It correctly notes that B is an axial vector and
  that `F = E + jB` handles inversion naturally, but it treats this as a
  correction to Maxwell (Maxwell's equations are parity-covariant when
  applied correctly), compares moving observers without Lorentz
  transformation, adds ad hoc fractal terms, and has no experimental
  validation. Motivates the parity and moving-observer tests in §3c.
  Copy in `reference/`.
- C. Doran, A. Lasenby, *Geometric Algebra for Physicists*, 2003.
- L. Dorst, D. Fontijne, S. Mann, *Geometric Algebra for Computer Science*,
  2007 — the most practical guide to implementing GA (bitmask blades, products).
- M. Castro-Núñez, R. Castro-Puche, "The IEEE Standard 1459, the CPC power
  theory, and geometric algebra in circuits with nonsinusoidal sources and
  linear loads," *IEEE Trans. Circuits and Systems I*, 2012.
- A. S. F. Komeno et al., "Three-phase quaternion power in three-wire systems,"
  *Renewable Energy and Power Quality Journal* (RE&PQJ).
- Geometric algebra power-theory work by the University of Almería group,
  e.g. *Mathematics* 9(11):1259, 2021, and the GA framework for balanced
  three-phase systems with sinusoidal and non-sinusoidal supply.
- ASTM G173 — reference solar spectral irradiance (Earth surface and
  extraterrestrial).
- P. T. Farnsworth, multipactor and electron-multiplier patents (§3i):
  US 1,969,399 (1934), US 2,091,439 (1937), US 2,135,615 (1938),
  US 2,174,487 (1939), US 2,179,996 (1939, with R. L. Snyder),
  US 2,141,837, US 2,180,279. Copies in `reference/` for the first five;
  2,141,837 and 2,180,279 still to be added.
- N. Tesla, "Coil for Electro-Magnets," U.S. Patent 512,340, 1894 —
  bifilar winding and the 250,000× energy example (§3e). Copy in
  `reference/`.
- N. Tesla, "Experiments with Alternate Currents of High Potential and
  High Frequency," lecture, London, February 1892 — disruptive condenser
  discharges, sudden versus harmonic impulses, brush rotation, local gas
  losses versus radiation (§3g pulse-duration experiment). Copy in
  `reference/` (`experiments_with_alternate_currents.pdf`, tesla.hu
  transcription).
- N. Tesla, "Method of and Apparatus for Electrical Conversion and
  Distribution," U.S. Patent 462,418, 1891 — condenser charged from a
  high-tension source, disruptive discharge into a working circuit,
  impulses of extremely small duration, same-direction or alternating
  (§3g). Copy in `reference/`.
- N. Tesla, "Apparatus for the Utilization of Radiant Energy," U.S. Patent
  685,957, 1901 — insulated plate exposed to radiation, condenser to
  ground, threshold discharge; one-sign impulses for X-ray tubes (§3g).
  Copy in `reference/`.
- N. Tesla, "Electrical Transformer," U.S. Patent 593,138, 1897 —
  quarter-wave secondary based on the propagation velocity in the coil;
  grounded end and maximum-potential terminal (§3e). Copy in `reference/`.
- D. W. Knight, "The self-resonance and self-capacitance of solenoid
  coils: applicable theory, models and calculation methods," v1.00, 2016,
  DOI 10.13140/RG.2.1.1472.0887 — coils as helical transmission lines,
  λ/2 free and λ/4 grounded resonances, helical velocity factor, free-coil
  SRF model, why turn-to-turn capacitance theories fail (§3e). Copy in
  `reference/` (`self-res.pdf`).
- Y. Liu, S. Jiang, Z. Zhong, J. Geng, F. Lv, "Alternating current
  discharge characteristics and simulation analysis of rod-plane short air
  gaps under salt fog conditions," *Energies* 11(9), 2251, 2018 — fog
  raises AC discharge voltage; conductivity effect (§3j).
- Y. Yuan, X. Jiang, S. M. Rowland, X. Cheng, Q. Li, "Effect of water
  streams on the AC breakdown performance of short rod-plane air gaps,"
  *IEEE Trans. Dielectr. Electr. Insul.* 21(4), 1747–1756, 2014,
  DOI 10.1109/TDEI.2014.004239 (§3j).
- K. R. Allen, K. Phillips, "Effect of humidity on the spark breakdown
  voltage," *Nature* 183, 174–175, 1959 (§3j).
- C. S. R. Bhoomireddy, *Recovery times of a spark gap with high pressure
  gases as insulation suitable for repetitive applications*, PhD thesis,
  HBNI/BARC, 2016 — deionization and thermal recovery phases (§3j).
- R. J. Lang, "Ultrasonic atomization of liquids," *J. Acoust. Soc. Am.*
  34, 6–8, 1962 — droplet size from frequency (§3j).
- M. Umezawa et al., "Effect of aerosol particles generated by ultrasonic
  humidifiers on the lung in mouse," *Part. Fibre Toxicol.* 10:64, 2013 —
  mineral particles from tap water versus RO water (§3j).
- P. Mali, *Magnetic Amplifiers: Principles and Applications*, John F.
  Rider, 1960 — saturable reactors, self-saturating amplifiers, gain,
  feedback, response time, efficiency (§3j). Two scans of the same 1960
  edition in `reference/`: `Rider-Magnetic-Amplifiers-Principles-And-
  Applications-Paul-Mali.pdf` (cleaner text, preferred for quotes) and
  `MagneticAmplifiers_text.pdf` (Internet Archive).
- C. P. Steinmetz, *Elementary Lectures on Electric Discharges, Waves and
  Impulses, and Other Transients*, McGraw-Hill, 1911 — single- and
  double-energy transients, natural impedance, quarter- and half-wave line
  oscillations, traveling waves and impulses, the arc as a source of
  growing oscillations, compound circuits and transition-point step-up,
  L and C of parallel round conductors (§3e, §3g). Copy in `reference/`
  (`elementarylectur00steirich.pdf`).
- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed.,
  Wiley–IEEE Press, 2008 — per-unit-length parameters, modal decoupling,
  terminal conditions, reflection coefficients, time-domain solutions
  (§3e bifilar line model). Copy in `reference/`.
- N. Tesla, "Dynamo-Electric Machine," U.S. Patent 406,968, 1889 —
  homopolar (unipolar) dynamo (§3f).
- P. Joot, *Geometric Algebra for Electrical Engineers* (GAelectrodynamics,
  v0.3.6.4), https://github.com/peeterjoot/GAelectrodynamics — Theorem 3.16,
  Lorentz force and power for electric and magnetic charges by grade
  selection. Reference only: no license file found.
- H. A. Rowland, magnetic effect of electric convection (1876) — rotating
  charged disk produces a magnetic field.
- W. C. Röntgen (1888) and A. Eichenwald (1903) — a dielectric moving in an
  electric field is magnetized (Röntgen current, M = P × v); H. A. Wilson &
  M. Wilson (1913) and Hertzberg et al., arXiv:gr-qc/0104101 (2001) —
  rotating magnetic dielectric; J. D. Jackson, *Classical Electrodynamics*,
  3rd ed., §11.10 — fields and constitutive relations in moving media.
- A. M. Cassie, arc model (CIGRE, 1939) and O. Mayr, arc model (1943) —
  dynamic arc conductance (§3g).
- R. W. Hellwarth, P. Nouchi, "Focused one-cycle electromagnetic pulses,"
  *Phys. Rev. E*, 1996 — exact free-space toroidal ("flying doughnut") pulse
  solution; FDTD validation case.
- A. Zdagkas et al., "Observation of toroidal pulses of light," *Nature
  Photonics*, 2022 — experimental generation of toroidal pulses.
- Yu. P. Raizer, *Gas Discharge Physics*, Springer, 1991 — breakdown, arc and
  spark plasma (arc and gap parameters, §3g).
- J. C. Martin, *J. C. Martin on Pulsed Power*, 1996 — pulsed breakdown of
  liquids and gases.
- IAPWS releases for heavy water (thermodynamic and transport properties) and
  for ordinary water — property tables in `fluid/media.go`.
- A. Taflove, S. C. Hagness, *Computational Electrodynamics: The
  Finite-Difference Time-Domain Method* — Yee scheme, absorbing boundaries.
- J. Koza et al., *Genetic Programming III* — evolved analog circuits,
  developmental encodings.
- K. Deb et al., "A fast and elitist multiobjective genetic algorithm:
  NSGA-II," *IEEE Trans. Evol. Comput.*, 2002.

**Gauge theory and quantum electrical metrology (H22)**
- K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* 10, 2445, 1974;
  J. Kogut, L. Susskind, "Hamiltonian formulation of Wilson's lattice gauge
  theories," *Phys. Rev. D* 11, 395, 1975 — gauge fields on a lattice, the
  non-abelian counterpart of Kron's 1944 network.
- D. Marcos, P. Rabl, E. Rico, P. Zoller, "Superconducting circuits for
  quantum simulation of dynamical gauge fields," *Phys. Rev. Lett.* 111,
  110504, 2013 — circuits emulating gauge fields.
- B. D. Josephson, *Phys. Lett.* 1, 251, 1962; K. v. Klitzing, G. Dorda,
  M. Pepper, *Phys. Rev. Lett.* 45, 494, 1980; J. P. Pekola et al.,
  "Single-electron current sources: toward a refined definition of the
  ampere," *Rev. Mod. Phys.* 85, 1421, 2013; BIPM, *The International
  System of Units (SI)*, 9th ed., 2019 — the volt, ohm and ampere from e
  and h.
- R. de-Picciotto et al., "Direct observation of a fractional charge,"
  *Nature* 389, 162, 1997; L. Saminadayar, D. C. Glattli, Y. Jin,
  B. Etienne, "Observation of the e/3 fractionally charged Laughlin
  quasiparticle," *Phys. Rev. Lett.* 79, 2526, 1997 — fractional charge
  as a collective electron excitation.

**Gravity and rotation (§3l)**
- C. W. Misner, K. S. Thorne, J. A. Wheeler, *Gravitation*, Freeman, 1973 —
  Schwarzschild clocks and orbits (§25), tidal forces.
- J. Plebanski, "Electromagnetic waves in gravitational fields," *Phys.
  Rev.* 118, 1396, 1960; L. D. Landau, E. M. Lifshitz, *The Classical
  Theory of Fields*, §90 — gravity as an effective medium.
- C. M. Will, "The confrontation between general relativity and
  experiment," *Living Rev. Relativ.* 17, 4, 2014 — PPN parameters.
- N. Ashby, "Relativity in the Global Positioning System," *Living Rev.
  Relativ.* 6, 1, 2003.
- R. V. Pound, J. L. Snider, *Phys. Rev.* 140, B788, 1965; I. I. Shapiro,
  *Phys. Rev. Lett.* 13, 789, 1964; S. S. Shapiro et al., *Phys. Rev.
  Lett.* 92, 121101, 2004; B. Bertotti, L. Iess, P. Tortora, *Nature* 425,
  374, 2003.
- G. Sagnac, *C. R. Acad. Sci.* 157, 708, 1913; E. J. Post, "Sagnac
  effect," *Rev. Mod. Phys.* 39, 475, 1967; W. W. Chow et al., "The ring
  laser gyro," *Rev. Mod. Phys.* 57, 61, 1985.
- J. Bailey et al., *Nature* 268, 301, 1977 (muon ring, clock hypothesis);
  B. Mashhoon, *Phys. Lett. A* 145, 147, 1990 (locality hypothesis).
- Rotor Mössbauer experiments (H17): H. J. Hay et al., *Phys. Rev. Lett.*
  4, 165, 1960; W. Kündig, *Phys. Rev.* 129, 2371, 1963; A. L. Kholmetskii
  et al., *Phys. Scr.* 77, 2008 and 79, 2009; T. Yarman et
  al., *Can. J. Phys.*, doi:10.1139/cjp-2015-0063 (arXiv:1503.05853); Y. Friedman,
  I. Nowik, *Phys. Scr.* 85, 065702, 2012; C. Corda, arXiv:2203.02282.
- Toray Composite Materials America, T1000G data sheet, 2024-02-16 —
  rotor burst limit.

---

## 7. Commands

```bash
go build ./...
go test ./...                         # unit + golden tests (golden uses committed CSVs)
go test -race ./...                   # required before merging anything in evo/ or solver/
go test -bench=. ./solver ./analysis  # performance benchmarks
go vet ./... && staticcheck ./...

# Regenerate golden files (requires ngspice on PATH; Linux/macOS/WSL)
go run ./internal/golden/cmd/regen -netlists testdata/netlists -out testdata/golden

# Simulate a netlist
go run ./cmd/kronspice testdata/netlists/rc_lowpass.cir -o out.csv

# Run a GA experiment
go run ./cmd/evolve -config examples/rc_lowpass_values/config.json

# Container checks (tools image: Go, gcc for -race, ngspice). On Windows,
# VS Code: Terminal > Run Task > "kronspice: …" (internal/golden/docker/run.ps1,
# logs in internal/golden/docker/logs/). Elsewhere:
docker build -t kronspice-tools internal/golden/docker
docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_ref.sh   # ngspice twins → testdata/golden
docker run --rm -v "$PWD":/src kronspice-tools sh internal/golden/ref/ngspice_all.sh   # every netlist through ngspice → testdata/golden/ngspice (TestNgspiceTwin)
docker run --rm -v "$PWD":/src kronspice-tools go test -race -timeout 60m ./...

# Runtime image (static kronspice + evolve, ~10 MB) and releases: Dockerfile;
# a v* tag runs .github/workflows/release.yml (tests, binaries for
# linux/windows/darwin × amd64/arm64, ghcr.io images kronspice and kronspice-tools)
docker build -t kronspice .
```

`go test -race` is about 10–30× slower: GA tests that compare an evolved
result set a generous `evo.Options.Timeout` and check `GenStats.TimedOut`
is zero, because a wall-clock timeout makes the result machine-dependent
(§1 goal 5).

---

## 8. Coding rules for Claude

1. **No panics in library code.** `netlist`, `circuit`, `device`, `solver`,
   `analysis` and `evo` return errors. Use typed sentinel errors
   (`errors.Is`-compatible) so `evo` can map them to penalties.
2. **No global mutable state.** Every simulation owns its matrices, devices and
   options. This is what makes parallel evaluation safe.
3. **SI units in code, float64 everywhere** (`complex128` for AC). Engineering
   suffixes are parsed only in `netlist/units.go`.
4. **Every new device or analysis ships with**: unit tests for its stamps, at
   least one golden test against ngspice, and a benchmark.
5. **Golden tolerances**: node voltages within `max(1e-6 V, 1e-3 relative)` for
   DC/AC; transient compared after resampling to a common time grid with
   tolerances stated in the test. Never loosen a tolerance to make a test pass
   without explaining why in the test.
6. **Fuzz the parser and validator** (`go test -fuzz`) — the GA will feed them
   inputs no human would write.
7. **Allocation-free inner loops.** `Load` and the LU solve must not allocate
   per iteration; check with `-benchmem`.
8. **Keep public APIs small** and documented with Go doc comments and runnable
   `Example` functions.
9. When porting an algorithm from SPICE3f5 or a book, cite the source
   (file/function or page) in a comment above the Go function.
10. Run `go test -race ./...` before declaring a task finished.
11. **Every element carries the full external port set** (§3.1) from the
    first device written in M1 onward — never add a device without it.
12. **Every new file carries the SPDX dual-license header** (§6, licensing).

---

## 9. Milestones

| # | Milestone | Done when |
|---|---|---|
| M1 | `exchange/` channels and ports, `device.Base`; netlist parser, R/V/I, dense MNA, `.OP` | Resistor networks match ngspice; every device exposes all channels; per-element × channel ledger report (electrical + thermal modeled, others "not modeled") |
| M2 | C, L, sources (PULSE, SIN, PWL), `.TRAN` with trapezoidal + LTE | RC/RLC step responses match golden |
| M3 | `.AC`, `.DC` sweep, measurements | Filter Bode plots match golden |
| M4 | Diode, Newton + limiting, gmin/source stepping | Rectifier and clipper circuits converge and match |
| M5 | `evo` value-genome engine, worker pool, checkpointing | RC low-pass example hits target cutoff within 1% |
| M5a | First experiment: Bi–Al capacitor element ∥ bifilar coil (`examples/bial_bifilar`) | Done: dry, wet and thermopile Pareto fronts, ordinary-coil comparison, ledger closes (`doc/experiments/`) |
| M6 | Sparse LU backend, benchmarks | Sparse and dense backends agree; speedup on >50-node circuits. Done: Gilbert–Peierls + minimum degree, stdlib only; all golden tests on both backends; 1.5× at 64 nodes, 5.6× at 576 (`doc/M6.md`) |
| M7 | BJT (Gummel–Poon), E/F/G/H sources | Amplifier example evolves to target gain/bandwidth. Done: 2N3904 stage evolved to 29.92 dB, 100 Hz–199.8 kHz from E-series values; 6 golden circuits on both backends (`doc/M7.md`) |
| M8 | Topology genome, validation, penalties | Filter evolved from an empty embryo meets spec. Done: 4-part elliptic low-pass grown from an empty fixture, mask met on a 10× denser sweep; seeds 1–8 all meet it (`doc/M8.md`) |
| M9 | MOSFET levels 1–3, NSGA-II | Multi-objective amplifier experiment. Done: levels 1–3 with 13 golden circuits on both backends; NSGA-II with constrained dominance reaches the ZDT1/ZDT2 fronts; a Level-3 NMOS amplifier evolved into a 99-design gain / bandwidth / power front (`doc/M9.md`) |
| M9a | `kron`: primitive network, connection tensor, mesh + node forms (linear DC/AC) | Agrees with MNA to 1e-9 on all linear golden circuits. Done: both forms on all 11 linear golden circuits, worst 9.8e-16; 300 random circuits with every linear element, worst 2.1e-11 (`doc/M9a.md`) |
| M9a.1 | `kron/orthogonal`, `exterior`, `chords`: open paths to exterior junctions for every element port; TAN chords | §3a.1 tests pass; open-path power equals MNA ledger; inactive paths carry zero. Done: ledger rows equal MNA to 1.1e-12 over 133 circuits, inactive paths exactly 0, chord = MNA `K` (now a full element), radiation edge = Balanis small loop; torn and reflecting-wall cases are M9b / `em` (`doc/M9a1.md`) |
| M9b | `kron`: diakoptics with manual tearing, cached torn solutions in `evo` | Value-genome GA runs faster with tearing than without, same results. Done: Aitchison-form torn solve equals MNA to 5e-11 over 2 289 cuts, power invariant, `m = b − N + S` after tearing, singular and nearly floating pieces reported; `evo.Problem.Tear` with a shared content-keyed cache; `examples/torn_ga` gives the same best design with the torn run 11× (line) and 5.6× (grid) faster than untorn Kron, still 1.5–2× slower than sparse MNA (`doc/M9b.md`) |
| M9c | `kron`: tearing of the Newton Jacobian for nonlinear circuits | Diode/BJT golden circuits match via torn solve. Done: every Newton Jacobian (OP, DC sweep, AC, transient) solved by Aitchison-form tearing through `analysis.Options.Solver` (`kron.JacobianTearing`); all 36 golden circuits (diodes, BJTs, MOSFETs) under 78 cuts, torn at resistors, capacitors and the nonlinear devices themselves (35 circuits split into pieces), match ngspice within the golden tolerances and the untorn run to 2.8e-12, with no untorn fallback; unchanged pieces keep their factors across Newton iterations and time steps; not yet faster than MNA on these sizes (`doc/M9c.md`) |
| M9c.1 | Kron's frame for nonlinear circuits and moving inductors: `kron/nonlinear.go`, `kron/machine.go`, N element (`device/machine.go`) | Nonlinear golden circuits solved in Kron's node-pair frame equal MNA; machine element checked against analytic cases and across Kron's frames I and II. Done: all 36 golden circuits (23 nonlinear) by state-dependent primitive tensors equal MNA to 3.7e-14 (OP) and 1.4e-13 (AC); Christoffel form (Kron 1951 eqs. 3–4) equals the element to 1e-12; induction motor in frames I and II agrees to 6e-4 of peak; pendulum (spatial) and parametric pump vs Mathieu (temporal) match with closed ledgers (`doc/M9c1.md`) |
| M9d | `hyper/quat` and `hyper/ga` (G2, G3, sparse Cl(p,q,r)) | Algebra property tests pass; quaternion ↔ G3 rotor conversion exact. Done: stdlib only; laws in seven signatures including degenerate ones; products equal Pauli (Cl(3,0)) and Dirac (Cl(1,3)) matrix products to 1e−12; quaternion ↔ G3 conversion exact, products to 1e−15; allocation-free G3 at 18 ns per product (`doc/M9d.md`) |
| M9e | `kron/gaframe`: loop/cut-set blades, singularity pre-check | Blade magnitude matches det(Z'); degenerate evolved circuits rejected before LU. Done: blade² = det(Cᵀ·Z·C) = product of rejections on 300 random networks and the linear goldens (1e−9), frame duality det Z′ = det Y′·ΠR, K₄ spanning trees; `kron.Precheck` names zero-impedance loops in 5.8 µs (MNA OP 1.7 ms) with no false positives, wired in as `evo.Problem.Precheck`; the conditioning diagnostic costs about a dense LU and is off by default (`doc/M9e.md`) |
| M9f | `kron/polyphase` + `analysis/power`: quaternion Clarke/Park, quaternion and GA power | Matches conventional transforms and IEEE-1459 quantities; used as a GA fitness in a power-filter example. Done: Clarke/Park rotors equal the matrices (1e−14); quaternion power invariant across frames; Kron's three-phase `ParkFrame` turns a salient stator into Ld, Lq, L0 with the right speed voltages; GA power equals IEEE 1459 P, S, S₁, D_I, D_V, S_H to 1e−9 and names the harmonic pairs behind the distortion; `examples/power_filter` evolves a six-pulse rectifier's filters from PF 0.903 to 0.957 (`doc/M9f.md`) |
| M9g | `em/energy` open-system energy ledger + Tellegen check, `em/limits` electrical-size check | Ledger residual within tolerance on all golden circuits (back-filled from M2); efficiency usable as GA fitness. Done: all 40 golden circuits balance in every analysis (no residual above round-off); Tellegen's sum vanishes at every two-terminal golden OP and in its general cross-circuit form; `em.CheckQuasiStatic` (λ/10, εr) and `SourceBandwidth`; `em.OPEfficiency` / `TranEfficiency` equal the analytic and AC values and drive `evo.OPEfficiency` / `TranEfficiency`; electrically large candidates penalized as invalid (`doc/M9g.md`) |
| M9h | `em/nabla`, `maxwell`, `biquat`: field equations in quaternion forms | Triple-twin tests and analytic cases pass. Done: Treatise quaternion form, spacetime algebra ∇F = J/(ε₀c) and vector calculus agree to 1e−9 on arbitrary sourced fields; ∇ exact on polynomials; charged sphere and current cylinder exact; plane wave ∇F = 0 converging at fourth order (5e−6 of kE₀ at 40 points/λ); potentials and gauge invariance; parity and rotations (vector rule = GA = recomputed Biot–Savart, 1e−12); Lorentz boost by rotor = textbook (1e−12) (`doc/M9h.md`) |
| M9i | `em/extract`, `em/tline`: lumped R/L/C/K and transmission lines from geometry | Analytic capacitance/inductance within 1%; evolved layout example runs. Done: 2-D finite-volume Laplace solve with cut-cell boundaries (coax 0.07 %, Cohn stripline −0.006 % extrapolated, two-wire 0.1 %, microstrip within 2 % of Hammerstad–Jensen), L = μ₀ε₀C₀⁻¹; Neumann mutual inductance = Maxwell's elliptic form (0.007 %), current sheet = Nagaoka; telegrapher ladder (Z₀, velocity) and coupled-line ladder with K and mutual C; `pcb_layout` GA evolves a 50 Ω pair below −30 dB NEXT (`doc/M9i.md`) |
| M9j | `exchange/neutrino.go`: cited neutrino coupling model on every element's neutrino port, declared flux input | Known-value test: interaction rate and deposited power for 1 kg of silicon at 1 AU match the cited calculation (orders of magnitude below thermal noise); reported as computed ≈ 0 in every ledger; never credited with a residual |
| M9m | `bifilar_coil` example: idealized + realistic comparison, evolved geometry | Idealized model gives 250,000; realistic results match analytic coil formulas within 5%; evolved coils reported against both baselines |
| M9n | `em/lorentz`, `em/duality`, `homopolar` + `rotating_efield` examples | Grade-selection and duality tests pass; homopolar EMF and shaft balance match; displacement-current B matches analytic case; hypothetical results labeled and kept out of the real ledger. Done: grade selection = vector forms (1e−12), duality invariance; Tesla two-disk dynamo as machine elements (EMF, torque, shaft power and ledger exact; GA meets 2 kW at η > 0.95); Rowland disk, charging capacitor and rotating-E fields checked with the M9h operators; momentum ledger (Newton's third law, exterior reactions, hypothetical column) and an electrodynamic-tether example (`doc/M9n.md`, `doc/em-push-pull.md`) |
| M9o | `device/arc`, `em/fdtd`, `em/topology`, `em/probes`, `arc_surge` example | All §3g validation cases pass (including the flying-doughnut pulse); probe data import and comparison work; GA runs on at least one fitness option |
| M9o.1 | `pulse_duration` example: impulse source, duration sweep, paired E/B probes, checks 1–7 | §3g pulse-duration validation cases pass; each check reports predicted value and deviation status; bench CSV import works |
| M9o.2 | `radiant_receiver` example: photoemission, ionization and atmospheric-field charging, threshold discharge | §3d radiant-receiver validation cases pass; ledger traces load energy to each declared input channel |
| M9q | Multipactor susceptibility check and multipactor gap model (§3i) | Resonance zones match published susceptibility charts; ledger closes for a driven gap |
| M9r | `device/core`, `device/royer`, `device/magamp`, `charge_pump` example: ground charge pump, HV bank, staged Royer step-down with magnetic amplifiers | §3j validation cases pass; charge from ground and energy by source both close in the ledger |
| M10 | Cross-check best circuits in QSpice via `qspice` Python package (Windows) | Evolved circuits reproduce in QSpice |

Work milestones in order. Do not start device models beyond the current
milestone until its golden tests pass.
