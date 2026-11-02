# Einstein's distant parallelism and Kron's machines

A. Einstein, "New possibility for a unified field theory of gravitation
and electricity," *Sitzungsberichte der Preussischen Akademie der
Wissenschaften* 1928, 224–227 (A. Unzicker's translation).

Einstein builds space-time from a **4-bein** — four frame legs at every
point — and a rule that legs at different points are parallel ("distant
parallelism"). That geometry has no curvature but has **torsion** Λ. He
proposes:
- the trace of the torsion, φ_μ, is the electromagnetic potential;
- to first order the theory gives both Einstein's vacuum gravity and
  Maxwell's equations.

**Question:** does this answer Kron's multi-dimensional, non-Riemannian
work?

**Short answer:** it supplies the geometric language Kron's machines
already use, and kronspice can run a motor from Einstein's torsion alone.
But Kron's machine torsion is not electromagnetism in Einstein's sense,
and Einstein's own theory did not survive.

```
go run ./examples/einstein_teleparallel -viz examples/einstein_teleparallel/viz   # about 2 seconds
python3 examples/einstein_teleparallel/viz/render.py
go test ./examples/einstein_teleparallel
```

## 1. Einstein's 4-bein is Kron's frame

Hoffmann's single-phase motor ([hoffmann_motor](../hoffmann_motor)) has
exactly four coordinates: three winding charges and the rotor angle. The
change to stationary rotor coils (Hoffmann eq. 35) is therefore a 4-bein,
Einstein's case n = 4.

The same torsion, computed three ways:

| | Agreement |
|---|---|
| Einstein's torsion Λ (his definition), in the frame's components | equals minus the frame's anholonomy (sign convention), to 6×10⁻¹¹ |
| Hoffmann's Ω (`kron.Torsion`, his eq. 19) | equals Λ with one index lowered by the inductances, Ω = L̄·Λ, to 3×10⁻¹⁷ |
| Torsion of Einstein's own connection | 2Λ, to 10⁻¹⁰ |

Kron's non-Riemannian frame is a teleparallel 4-bein. Picture: left half
of `viz/einstein_teleparallel.png`.

## 2. Einstein's electromagnetic potential on the motor

| Frame of the same motor | φ = trace of the torsion |
|---|---|
| Kron's stationary rotor coils (a rotation) | **0** — a "pure gravitational field" to Einstein |
| rotor legs also growing as e^{0.3θ} | φ_θ = 0.300 |

In a Kron machine the frame is just a choice of current paths: both frames
describe the same motor. So φ here belongs to the frame, not to the
machine. This echoes Einstein's own remark that "the separation of the
gravitational and electromagnetic field seems arbitrary".

## 3. Curvature: the real difference

| Connection | Curvature |
|---|---|
| Einstein's distant parallelism, either frame | **0** — all its geometry is torsion |
| The machine's own (Christoffel symbols of its inductances and inertia) | **not zero** — the machine's space is curved |

Kron's torsion is that curved space seen from a moving frame (Hoffmann
eq. 17); a change of frame cannot remove curvature. Einstein's torsion
*replaces* curvature. The mathematics of frames is the same; the geometry
is different.

## 4. Einstein's first approximation, checked

Einstein's linearized field equation (his eq. 5) was solved for plane
waves, with x⁴ imaginary as in his paper:

| Wave | Solutions | φ | His eqs. 8, 9, 7 | Ricci tensor when φ = 0 |
|---|---|---|---|---|
| not light-like (two cases) | 4 (coordinate changes only) | 0 | hold to 10⁻¹⁶ | 0 |
| light-like | 10 | up to 0.5 | hold exactly | 0 |

On the light cone φ obeys ∂φ = 0 and □φ = 0, Maxwell's vacuum equations,
and with φ = 0 the metric obeys vacuum gravity. **His first-approximation
claims hold.**

## 5. kronspice: the motor built from Einstein's torsion

The motor's machine model was built using only the 4-bein and its
torsion: constant inductances L̄ and a speed matrix G = 2(Ω_θλν + Ω_νλθ).
No inductance derivatives were used. Run: 120 V, 60 Hz, pushed to
40 rad/s, 4 s.

| Model | Final speed | Ledger |
|---|---|---|
| frame I (coils turn with the rotor) | 376.7928 rad/s | 10⁻⁶ |
| Kron's frame II(b) (`kron.RotatingFrame`) | 376.8235 rad/s | 5×10⁻⁸ |
| **built from Einstein's torsion** | **376.8235 rad/s** | 5×10⁻⁸ |
| same inductances, torsion removed | 26.81 rad/s — coasts on friction only | 4×10⁻⁷ |

The torsion-built motor matches Kron's frame II(b) to 2×10⁻¹² rad/s.
Without the torsion there is no torque. Picture: right half of
`viz/einstein_teleparallel.png`.

## Where it stands

- **Einstein's theory was dropped.** He worked on distant parallelism
  from summer 1928 to spring 1931. He gave it up because no single set of
  field equations could be justified (Sauer 2006) — the "unsecurity" his
  own note at the end of this paper admits.
- **Its gravity part survives** as the teleparallel equivalent of general
  relativity. Electromagnetism from the torsion trace did not survive.
- **For Kron.** Kron's non-Riemannian machine dynamics (1934; Hoffmann
  1949) uses the same objects — n-beins, anholonomy and torsion — as a
  working engineering method. Hoffmann himself later worked with Einstein
  (Einstein, Infeld and Hoffmann, 1938).
- **What this example does not decide.** Whether Kron drew on Einstein's
  1928–31 papers is a history question; it would take Kron's 1934
  references to settle.

## Files

| File | What it holds |
|---|---|
| `frame.go` | the motor, the 4-bein (Hoffmann's frame), Einstein's torsion, anholonomy, φ, L̄ |
| `geometry.go` | Einstein's distant-parallelism connection and the machine's connection; curvature |
| `linear.go` | Einstein's eq. 5 for plane waves; checks of eqs. 7–9 and R_αβ |
| `machine.go` | the machine model built from the torsion; kronspice runs |
| `main.go` | the report |
| `main_test.go` | torsion three ways, φ and curvature, the linear check, the motor |
| `viz/render.py` | `einstein_teleparallel.png` |

## Sources

- A. Einstein, *Sitzungsber. Preuss. Akad. Wiss.* 1928, 217–221 (Riemann geometry with distant parallelism) and 224–227 (this paper); A. Unzicker's translation.
- B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Rev. Mod. Phys.* 21, 535 (1949).
- G. Kron, "Non-Riemannian Dynamics of Rotating Electrical Machinery," *J. Math. Phys.* 13, 103 (1934).
- T. Sauer, "Field equations in teleparallel space-time: Einstein's Fernparallelismus approach toward unified field theory," *Historia Mathematica* 33, 399 (2006), arXiv:physics/0405142.
- J. W. Maluf, "The teleparallel equivalent of general relativity," *Ann. Phys.* 525, 339 (2013), arXiv:1303.3897.
