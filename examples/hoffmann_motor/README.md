# Hoffmann's account of Kron's non-Riemannian electrodynamics: the single-phase motor

B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Reviews of Modern
Physics* 21, 535–540 (1949), explains for physicists what Kron did with
rotating machines, and works one example: a single-phase squirrel-cage
induction motor. This example builds that motor and solves it the ways the
paper describes.

```
go run ./examples/hoffmann_motor -viz examples/hoffmann_motor/viz   # about 10 seconds
python3 examples/hoffmann_motor/viz/render.py
go test ./examples/hoffmann_motor ./kron
```

## The idea in plain words

A machine's currents and its rotor angle together are the coordinates of
one moving point. The inductances and the inertia act as the "geometry"
(the metric) of the space it moves in, and Lagrange's equations are the
path of the point (Hoffmann eqs. 1–8).

- **Frame I: coils that turn with the rotor** (Fig. 2). The squirrel cage
  is replaced by two perpendicular coils fixed to the rotor. Their
  inductances to the stator change with the angle θ (eqs. 24, 42). The
  space is ordinary (Riemannian) and the equations are Lagrange–Maxwell
  with Christoffel symbols (eq. 5). kronspice already had this
  (`kron.Christoffel`, `kron.Geodesic`).
- **Stationary coils 2̄, 3̄** (Fig. 3, eq. 35). Engineers prefer current
  paths that stay still while the cage turns, because then every
  inductance is a constant (eq. 29). The price is that this change of
  frame is *non-holonomic*: no set of coordinates goes with it (eq. 38).
  The space gains **torsion** Ω (eq. 19). The equations gain a term
  2Ω·ẋ·ẋ (eq. 5′), and that term carries both the speed voltages and the
  torque.
- **Hoffmann's warning** (p. 540): putting the constant inductances into
  the ordinary Lagrange–Maxwell equations "would obtain wrong results".

New in kronspice: `kron.Torsion` (eq. 19), `kron.NewMetric` (eq. 15) and
`kron.NonRiemannian` (eq. 5′), in `kron/hoffmann.go`.

## The motor (all values declared)

Stator coil L11 = 0.25 H, 2 Ω; rotor coils L_S = 0.25 H, 3 Ω each; mutual
M = 0.23 H; salient pole L_D = 0 (and 0.02 H in one case); inertia
0.001 kg·m²; friction 10⁻⁴ N·m·s; 120 V, 60 Hz; two poles, so synchronous
speed is 377 rad/s.

## 1. The two frames

In the stationary-coil frame the inductances no longer depend on θ:

| | 1 | 2̄ | 3̄ | θ |
|---|---|---|---|---|
| 1 | 0.25 | 0.23 | 0 | 0 |
| 2̄ | 0.23 | 0.25 | 0 | 0 |
| 3̄ | 0 | 0 | 0.25 | 0 |
| θ | 0 | 0 | 0 | 0.001 |

But the frame has torsion: Ω(2̄θ3̄) = −0.125, Ω(3̄θ1) = +0.115,
Ω(3̄θ2̄) = +0.125 (and the same with the first two indices swapped and the
sign changed). Ω does not change with the angle.

At 200 random states (angles, currents, their rates, speed, acceleration),
eq. 5′ with the torsion gives the frame-I forces to 5×10⁻¹⁶. Without the
torsion term it is wrong by up to 293 % of the largest force.

## 2. Starting

Each run is 5 seconds long.

| Run | Final speed | Energy |
|---|---|---|
| from rest | **0** — it does not start | 1,409 J from the mains, all heat |
| pushed to +40 rad/s, frame I (kronspice) | **376.587 rad/s** | 867 J = 796 J heat + 70.5 J stored |
| same, frame II(b) (kronspice, Kron's speed matrix) | 376.618 rad/s | same |
| same, **eq. 5′ integrated directly** | **376.587 rad/s** | — |
| pushed to −40 rad/s | −376.587 rad/s | runs backwards just as well |
| pushed, without the torsion term | 24.261 rad/s | only coasts down on friction |

- **No start from rest.** A single-phase motor has no starting torque: its
  pulsating field is two equal fields turning in opposite directions. Real
  motors add a start winding or a capacitor.
- **Pushed, it runs either way.** It goes up to just under synchronous
  speed (slip 0.1 % with no load).
- **The three solutions agree.** Frame I and eq. 5′ agree to 0.0004 rad/s.
  Eq. 5′ is integrated by this example's own Runge–Kutta code, which
  shares nothing with kronspice's machine element.
- **Hoffmann's warning.** Without the torsion there is no speed voltage
  and no torque. The pushed rotor follows friction alone:
  40·e^(−5 s·friction/J) = 24.261 rad/s.

Every kronspice run's energy ledger closes to 2×10⁻⁶ or better. The motor
adds no energy: every joule comes from the 120 V supply (or the push).

## 3. Torque–speed curve

With the shaft held at each speed, mean torque in N·m:

| rad/s | 0 | 62.5 | 125 | 188 | 250 | 313 | 375 |
|---|---|---|---|---|---|---|---|
| L_D = 0 | 0 | 0.071 | 0.152 | 0.262 | 0.433 | 0.712 | 0.071 |
| L_D = 0.02 H (salient pole) | 0 | 0.231 | 0.485 | 0.785 | 1.144 | 1.383 | 0.085 |

- **Shape.** The torque is zero at standstill and odd in speed (the same
  curve the other way round). It is largest near 91 % of synchronous
  speed. Frames I and II(b) agree to 0.2 % of the peak.
- **Salient pole.** It makes the stationary coil 2̄ couple more tightly to
  the stator (L̄22 = L_S − L_D = M), which roughly doubles the torque.
- **No-load slip.** The salient pole lowers it from 0.107 % to 0.053 %.
  The motor still does not lock into synchronous speed at these values: a
  reluctance motor needs a much larger L_D.

Picture: `viz/hoffmann_motor.png`.

## Files

| File | What it holds |
|---|---|
| `motor.go` | the motor in frame I (eq. 42) and frame II(b); start-up and fixed-speed runs in kronspice |
| `eq5.go` | Hoffmann's eq. 5′ integrated directly (Runge–Kutta), torsion constancy check |
| `main.go` | the report; the random-state twin |
| `main_test.go` | twin, torsion constancy, start-up agreement, coast-down, torque curve |
| `../../kron/hoffmann.go` | `Torsion`, `NewMetric`, `NonRiemannian` (eqs. 15, 19, 5′), tested in `kron/hoffmann_test.go` |

## Sources

- B. Hoffmann, "Kron's Non-Riemannian Electrodynamics," *Rev. Mod. Phys.* 21, 535–540 (1949).
- G. Kron, "Non-Riemannian Dynamics of Rotating Electrical Machinery," *J. Math. Phys.* 13, 103–194 (1934).
- G. Kron, *Equivalent Circuits of Electric Machinery* (Wiley, 1951), Epilogue — frames I and II(b) as used by `kron.RotatingFrame`.
