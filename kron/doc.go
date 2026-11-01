// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package kron is Gabriel Kron's tensor analysis of networks (CLAUDE.md
// §3a): a circuit is a primitive network — every branch on its own, with
// its impedance or admittance and its impressed sources — and
// interconnecting the branches is a change of reference frame described by
// a connection tensor. The connected network's equations follow from the
// transformation laws:
//
//	mesh form       i = C·i′,  Z′ = Cₜ·Z·C,  e′ = Cₜ·e,  Z′·i′ = e′
//	node-pair form  E = A·v,   Y′ = Aᵀ·Y·A                (A = Cₜ⁻¹ in Kron)
//
// (G. Kron, Tensor Analysis of Networks, Wiley 1939, "Kron 1939": Ch. IV
// pp. 101–104, eqs 4.23–4.37; Ch. V pp. 119–121 for a rectangular C;
// Ch. XIV p. 360 for the node-pair form; Tensors for Circuits, Dover 1959,
// Ch. 6 p. 43 and Ch. 7 pp. 46–50.)
//
// # Frames
//
// Each tensor carries its frames in its Go type (Mat[R, C], Vec[F]):
// Primitive (one axis per branch), Mesh (one per fundamental loop,
// m = b − N + S of them, Kron 1939 p. 75), Node (one potential per node
// other than the datum), and two frames for the axes that a pure mesh or
// node formulation cannot carry: Junction, the voltages across branches
// whose current is impressed (current sources, controlled current
// sources, control coils), which Kron treats as junction pairs rather than
// meshes (p. 411), plus the capacitor voltages; and Apparent, the currents
// of branches with no admittance (voltage sources, controlled voltage
// sources), Kron's apparent coils (pp. 407–408), plus the inductor
// currents. The reactive state variables are axes because recovering them
// from the other frame divides by ωC or ωL and loses accuracy as ω → 0
// (at DC a capacitor is an open junction pair and an inductor an apparent
// coil exactly). A transformation between the wrong frames does not
// compile.
//
// # Conventions
//
// Branch k runs from node From to node To of its element, the SPICE
// n+ → n− direction: its voltage is E = V(n+) − V(n−) and its current
// flows n+ → n− through it (for a source, the SPICE branch current).
// Each branch obeys one linear relation in the branch voltages and
// currents — Kron's general primitive equation E + e = z·(i + I) (Kron
// 1939 p. 417; eq 19.3 p. 480) written without infinite entries, so
// voltage sources (z = 0) and current sources (y = 0) need no special
// cases. A voltage-controlled source gets a control coil, an open branch
// (no admittance) across its controlling node pair, so its gain is a
// mutual term between two primitive axes; a current-controlled source
// couples to the axis of the voltage source it senses.
//
// The package solves linear DC and AC (rule 3: start linear); every result
// has an MNA twin test to 1e−9 relative (rule 4). The connection tensor C
// is built from a breadth-first spanning tree rooted at ground: its
// columns are the fundamental loops, one per link, so Aᵀ·C = 0 exactly and
// KCL holds by construction. Systems are factored and solved, never
// inverted (rule 2).
package kron
