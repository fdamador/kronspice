# Kelvin's vortex atoms through his ring-electromagnet analogy

Kelvin's 1867 "On Vortex Atoms" proposes that atoms are closed vortex rings
in a perfect liquid, and notes that the liquid's velocity round such a ring
is the magnetic field of the "corresponding ring electromagnet", with the
ring's impulse equal to that electromagnet's magnetic moment (pp. 9–10).
This example takes the analogy literally: circulation Γ plays the current I,
velocity plays B/μ₀, and the ring's kinetic energy is ½ρΓ²·L/μ₀ with L a
self-inductance. So the engine's inductance and Biot–Savart code (M9i,
`em.MutualCircles`, `em.Loop`) computes the vortex ring directly. It sits
beside the charged ring of [kelvin_ring](../kelvin_ring) on the electron
test bench (`doc/kelvin-references.md`, Vol. IV section).

```
go run ./examples/vortex_ring
go test ./examples/vortex_ring
```

## The construction

1. **Ring speed.** The self-inductance of a ring conductor of radius R and
   round core a with uniform current is Maxwell's mutual inductance of
   coaxial circles (Treatise Art. 701) averaged over pairs of filaments in
   the core section (Gauss–Legendre in ρ², midpoint in φ). The ring's speed
   is U = dE/dP at fixed circulation and core volume (a²R constant), with
   E = ½ρΓ²L/μ₀ and impulse P = ρΓπR², and is compared with Kelvin's
   U = Γ/4πR·(ln 8R/a − ¼) (Vol. IV Art. 3, as quoted by Saffman §10.2).
2. **Leapfrogging** (Helmholtz 1858; Tait's smoke rings, Kelvin pp. 2,
   11–12). Two coaxial rings of equal Γ are integrated with RK4; each moves
   at its self speed plus the other ring's Biot–Savart flow at its centre
   line. The energy is the self terms plus the mutual-inductance
   interaction, so impulse and energy must stay constant.
3. **Where the analogy stops.** The flow and the field share their shape
   but not their dynamics.

## Results (Γ = 1, R = 1)

| Core a/R | L/μ₀R (numeric) | thin ring ln 8R/a − 7/4 | U (numeric) | Kelvin's U | Difference |
|---|---|---|---|---|---|
| 0.20 | 1.95837 | 1.93888 | 0.27100 | 0.27366 | −0.97 % |
| 0.10 | 2.63730 | 2.63203 | 0.32798 | 0.32882 | −0.26 % |
| 0.05 | 3.32624 | 3.32517 | 0.38369 | 0.38398 | −0.07 % |

The difference falls as (a/R)², the order of the terms the thin-ring
formula drops. The flow at the ring's centre is 0.50000 = Γ/2R, as the
field at the centre of a loop is μ₀I/2R, and the impulse ρΓπR² is the
counterpart of the moment IπR² = 3.14159 (Kelvin p. 10).

**Leapfrogging**, two rings of core 0.1, 0.6 apart, 60 time units:

| | Value |
|---|---|
| Passes | 21, period 2.92 |
| Radii swing | 0.730 … 1.211 |
| Pair advance | 0.4990 (single ring 0.3288) |
| Impulse drift | 2.3×10⁻⁵ |
| Energy drift | 2.1×10⁻⁵ |

The constant energy also checks that the Biot–Savart velocities and the
elliptic mutual inductance agree with each other.

**Same shape, different dynamics.** A vortex ring travels because its
vorticity is carried along by the flow it makes. A current loop with the
identical magnetic field is not moved by that field: in Maxwell's
equations the ring electromagnet stays put. A self-propelled toroidal
*electromagnetic* structure has to come from the full field equations,
not from this static analogy; the standard case is Hellwarth and Nouchi's
flying doughnut (`em/doughnut.go`), which travels at c.

A vortex ring in a neutral liquid has no charge and no coupling to an
external field, so the Zeeman bench of [kelvin_ring](../kelvin_ring)
cannot be applied to it until a charge–vortex coupling is specified.

## Files

| File | What it holds |
|---|---|
| `main.go` | core-averaged self-inductance, Kelvin's speed, Biot–Savart velocity, the leapfrog integrator |
| `main_test.go` | speed error shrinking as a², the centre-flow identity, leapfrog invariants |

## Sources

- W. Thomson (Lord Kelvin), "On Vortex Atoms," *Proc. R. Soc. Edinb.* 6
  (1867); *Mathematical and Physical Papers* Vol. IV (1910), Art. 1,
  pp. 1–12; Art. 3, "The Translatory Velocity of a Circular Vortex Ring,"
  p. 67 (Internet Archive `in.ernet.dli.2015.212310`, not yet in
  `reference/`).
- H. von Helmholtz, "Über Integrale der hydrodynamischen Gleichungen,
  welche den Wirbelbewegungen entsprechen," *J. reine angew. Math.* 55,
  25 (1858) — leapfrogging rings.
- P. G. Saffman, *Vortex Dynamics*, Cambridge, 1992, §10.2 — the
  thin-ring speed formula as quoted.
- J. C. Maxwell, *Treatise*, Vol. II, Art. 701 — mutual inductance of
  coaxial circles.
