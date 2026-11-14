# CLAUDE.md — examples/void_shells (light–ether void shells, proposed H50)

Notes for Claude sessions working in this example. The repo's root
CLAUDE.md governs everything here: standard prediction beside every
hypothesis, the energy ledger closes (§3c), cite source and equation for
every formula, hypothetical results labelled. New example ⇒ a row in
`examples/README.md`, its own `README.md`, and a citation in a `book/`
chapter (not yet done: the user has uncommitted edits in both). Every Go
file starts with the SPDX dual-license header.

## The question (user, 2026-11-13)

Simple structure: a particle as a void in a nonlinear ether, light pushing
out and ether pushing in, built up in void shells. Do shells form, and can
they give g = 2?

## Files

- `field.go`: the solver (GP ether + m = ½ light modes, one photon per mode; relax, RK4 real time, measurements).
- `main.go`: Parts 1–6 (round bubble vs ε, excited modes, prepared shell, N photons, prepared onion, the g bound).
- `main_test.go`: Jacobi eigen check; bubble against the Python reference; shell dissolves; onion merges; real-time energy.
- `ref/void_shells_ref.py`: independent Python reference (scipy eigsh).
- `viz/plot_shells.py` → `viz/shells_g.png`, `viz/shells_gallery.png`, `viz/shells.mp4|gif` (from `viz/states.csv`, `viz/frames/`, frames not in git).
- `viz/make_report.py` → `void_shells_report.pdf`.

## Register entry (proposed H50; add to root §1 when merged — H49 is proposed in examples/ether_horizon)

| H50 | Particles as voids in a nonlinear ether, light pushing out and ether pushing in, built up in void shells (user hypothesis, 2026-11-13; the simplest structure after ether_horizon Part 13) | `examples/void_shells`: GP ether + double-loop light modes (m = ½, one photon per mode, V0 = 20) on an axisymmetric 0.4 ξ grid, relaxed and in real time; Python reference agrees to 11 digits | bounded — one photon makes a stable round resting bubble that keeps its light (light share ~73 %, g = 1.356–1.373 for ε = 30–3000, about 6 ƛ across); higher modes do not let ether into their nodes (thinner than ξ), the void stretches into a chain of lobes, g 1.31–1.37; a prepared shell (ether core in a void shell) dissolves into the bubble (relaxed: same M to 9 digits; real time: core squeezed out as sound by t ≈ 10); N photons one per mode share one void (g 1.373 → 1.359); a prepared onion (void, ether layer, void; two photons) merges into the two-photon bubble; scaling gives g = 1 + (E_P + E_W)/(3E_P + 2E_W) ∈ [4/3, 3/2] for hard walls, g = 2 needs a λ¹ (constant-force) ether energy, which only a vortex line gives; standard prediction: pointlike electron, g = 2.00231930436 (Dirac + QED) |
