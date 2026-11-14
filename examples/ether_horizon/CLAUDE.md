# CLAUDE.md — examples/ether_horizon (a void and a horizon in a nonlinear ether, H49)

Notes for Claude sessions working in this example. The repo's root
CLAUDE.md governs everything here: standard prediction beside every
hypothesis, the energy ledger closes (§3c), cite source and equation for
every formula, hypothetical results labelled. New example ⇒ a row in
`examples/README.md`, its own `README.md`, and a citation in a `book/`
chapter. Every Go file starts with the SPDX dual-license header.

## The question (user, 2026-11-13)

Can a nonlinear ether open an empty void and a horizon (a black region)?
Would such a horizon at the centre of the photon electron (H23) act on its
light ring as the membrane's 377 Ω? And if a horizon needs an inflow,
where does the ether go (the sink)?

## Files

- `main.go`: Parts 1–4 (vortex void, horizon, coupling, the ring circuit with horizon resistors).
- `sink.go`: Part 5, steady GP inflow into a localized loss −iγψ.
- `far.go`: Parts 6–7, the sink's columns, heat-flow bound, Lagally attraction, resistor-lattice twin.
- `viz/sink_cross_section.svg|png`: the setup in cross-section.
- `fountain.go`, `spin.go`, `part8.go`: Part 8, the fountain (steady swirling Bragg–Hawthorne flow with an axis vortex) and the spinning core as a battery in kronspice.
- `viz/plot_fountain.py` → `viz/fountain.png`; `viz/spin_core_smith.svg|png` from the Go run.
- `viz/plot_sink3d.py` → `viz/sink3d.gif`: the photon electron's light twice round its torus
  (from `examples/photon_electron/viz/plot_structure.py`) with the sink at its centre.
- `go run ./examples/ether_horizon -viz examples/ether_horizon/viz` writes `viz/sink_profile.csv`.

## Register entry (proposed H49; H47 and H48 are taken in the root register; add to root §1 when merged)

| H49 | A nonlinear ether can open an empty void and a horizon (a black region): Tewari's void (H13) as the empty core of a fast enough ether vortex, an acoustic "dumb hole" where the ether flows in faster than its own waves, and such a horizon at the centre of the photon electron (H23) acting on its light ring as a resistor, the 377 Ω of the black-hole membrane (user hypothesis, 2026-11-13; the tearing route to a cavity is `examples/ether_cavitation`) | `examples/ether_horizon`: Gross–Pitaevskii vortex profile by Newton FD (n = 1–3), ergoregion, draining-bathtub horizon and Hawking temperature with c_ether = c, electron-mass GR hole (naked Kerr–Newman, Hawking, Carter g = 2), membrane (Z0) and Page absorption of the ring's light against horizon radius, photon_electron's crossed ring line with a horizon resistor per section in kronspice AC (Q and frequency against the ladder's complex dispersion, horizon power booked, ledger) | bounded — void: the axis empties for every winding n ≥ 1 (f ≈ 0.58318ρ, literature 0.58278; line energy ln(1.464 R/ξ), Fetter 1.46), core f² = ½ at 1.56 ξ, but the void is a line one healing length thick (a thin torus for a ring), never a sphere, and it needs an order parameter that can vanish — the one that lets vortex knots untie (Kleckner 2016), unlike H29/H44's unit-vector ether; horizon: a vortex alone has an ergoregion (1.83 ξ) and no horizon, a horizon needs inflow into a sink (an unbooked ether column); with c_ether = c an electron-sized horizon (1.93e-13 m) has T_H = 1.9e9 K and radiates the rest energy in 2.4e-20 s (one channel, a lower bound on the rate), and lasting the age of the universe needs r_h ≥ 8.2e5 m; coupling: the ring keeps its light for the age of the universe only if r_h ≤ 2.9e-32 m (membrane) or 1.3e-22 m (Page), so no ether horizon size satisfies both, by 37 and 28 orders; GR at electron mass: r_s = 1.35e-57 m (8e-23 l_P), a/M = 2.9e44 and r_Q/M = 2.0e21 (naked, no horizon), Hawking 6e-107 s; kronspice: horizon resistors set Q_h = ∞, 1000, 100, 30 give AC Q 34.27, 33.14, 25.53, 16.02 against the ladder's 34.27, 33.14, 25.53, 16.00, frequency to 2e-7, power split closes to 5e-10, the crossed ring's sign flip (spin ½) kept, frequency (and g) moved only at second order in the losses; sink (`sink.go`, `far.go`): a localized loss −iγψ in the GP ether makes a true acoustic horizon outside the sink for γ ≳ μ/ħ (γ_c = 1.28 at r_s = 1 ξ, horizon at 1.52 ξ for γ = 10), flux 2γπr_s² when weak and falling past its maximum (Zeno), mass ledger 1e-7; the flow pays the Hawking radiation, so the bound is Earth's heat flow (47 TW, ≤ 2.7e-38 W per electron): a horizon must be ≥ 2.2e9 m, overlapping through all matter — a sink electron must be subsonic, with no horizon; two sinks attract (Lagally, inverse square, gravity's sign, not Coulomb's, 4.2e42 weaker); fitted to gravity, ether kept means every body grows at √(4πGρe), so ρe ≤ 1.2e-32 kg/m³ from d(GM_sun)/dt (1e-13/yr), destroyed means an energy column off the books (≤ 1.3e-40 kg/m³ if it became heat); kronspice 21³ resistor lattice twin: face inflow = Q to 2e-14, source = resistor power, v·4πd² → 1 (1.08 at 6 cells); fountain (Part 8): exact steady swirling flow (Bragg–Hawthorne, C = √(C0² + k²ψ²), Chandrasekhar–Kendall l = 1, 2) closes the books with no sink; without the axis vortex no streamline torus winds even once, with it a q = 2 torus of tube/ring ½ is fitted (β = 1.51, a = 316 fm, one equatorial ring for l = 1; β = 1.61, a = 352 fm, two rings at ±137 fm for the two-pole flow l = 2) — a fit, not a prediction; spinning core as a membrane battery in kronspice: +R_H with sources off, −R_L at its load (Γ_core·Γ_load = 1), maximum power at Ω_F = Ω_H/2; superradiance for the m = ½ ring needs Ω > 2ω_C (rim speed c·r_c/R), but a spin of ħ/2 feeds the Q = 0.55 ring for 1.8e-22 s (0.02 period); standard prediction: no horizon or void inside the electron, a point particle to < 1e-18 m with g − 2 from QED |
