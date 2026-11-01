// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package spacetime is the engine's gravity and fast-rotation layer
// (CLAUDE.md §3l): how curved spacetime and rotating frames change clocks,
// signal frequencies, light paths and the fields a circuit sees.
//
// Three different "curvatures" are kept apart:
//
//  1. Gravitational curvature (this package, Schwarzschild exterior):
//     clock rates, redshift, light delay and bending. Maxwell's equations
//     in a static gravitational field are equivalent to flat-space
//     Maxwell in a medium with ε_r = μ_r = n(r) (J. Plebanski, "Electro-
//     magnetic waves in gravitational fields," Phys. Rev. 118, 1396,
//     1960; L. D. Landau, E. M. Lifshitz, The Classical Theory of Fields,
//     §90), so em/ field solvers can take gravity as a refractive index.
//  2. Rotating frames in flat spacetime (this package, rotation.go):
//     clocks on a rotor run slow by √(1 − ω²r²/c²); light around a
//     rotating loop shows the Sagnac delay 4Aω/c². A probe riding a rotor
//     sees the fields of the local co-moving inertial frame (Lorentz
//     transformation at the local velocity ω×r, em/frames.go — the
//     "locality hypothesis"; B. Mashhoon, Phys. Lett. A 145, 147, 1990).
//  3. Kron's curved configuration space of a machine (kron/machine.go):
//     windings and rotor angle as coordinates, inductances L(θ) as the
//     metric, speed voltages and torque from Christoffel symbols (Kron,
//     J. Math. Phys. 13, 103, 1934). Not spacetime curvature; same
//     mathematics.
//
// Rotation and curvature together (machine.go): by the equivalence
// principle a rotor's hub follows the field (Schwarzschild rate at its
// position and speed) while its rim is special relativity in the hub's
// local inertial frame; its spin axis precesses by Thomas + geodetic +
// Lense–Thirring (L. I. Schiff, Phys. Rev. Lett. 4, 215, 1960), checked
// against Gravity Probe B (C. W. F. Everitt et al., Phys. Rev. Lett. 106,
// 221101, 2011); a loop riding it sees the Sagnac delay relative to the
// local gyroscopes. MachineSelfField estimates the reverse: the curvature
// and frame dragging the machine's own mass and spin produce.
//
// Conventions: SI units; G and c from CODATA 2018; masses given as
// GM [m³/s²] where measured that way (planetary GM is known far better
// than G or M separately). Rates dτ/dt are proper time per coordinate
// (distant-observer) time. Fractional frequency shifts are
// (f_received − f_emitted)/f_emitted.
//
// Validity: exterior of a non-rotating spherical mass (Schwarzschild);
// a rotating body's frame dragging is modeled to first order
// (Lense–Thirring, weak field); the full Kerr metric is not.
// Weak-field results carry the PPN parameter γ (= 1 in general
// relativity) so measured bounds can be applied: Cassini,
// γ − 1 = (2.1 ± 2.3)×10⁻⁵ (B. Bertotti, L. Iess, P. Tortora, Nature
// 425, 374, 2003).
package spacetime
