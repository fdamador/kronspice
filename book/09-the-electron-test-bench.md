# 9. The electron test bench

Four models of the electron as a ring, a vortex or a loop of light have
been put through the engine, and this chapter is about the bench they
are measured on rather than about any one of them. The bench is a fixed
set of measured numbers that every model must reproduce before it is
taken further, and a fixed standard prediction, the Dirac electron with
QED corrections, that is printed beside every result. The five examples
are [kelvin_ring](../examples/kelvin_ring), which builds the bench from
Kelvin's 1899 charged-ring molecule with nothing but the Lorentz force,
[vortex_ring](../examples/vortex_ring), its hydrodynamic cousin,
[space_vortex](../examples/space_vortex) (register entry H13),
[kron_electron](../examples/kron_electron) (H21) and
[photon_electron](../examples/photon_electron) (H23). The reading notes
behind them are [doc/kelvin-references.md](../doc/kelvin-references.md)
and [doc/tewari-references.md](../doc/tewari-references.md).

## 9.1 The bench

None of these examples is a netlist. Each is a Go program that takes the
measured constants, builds the model's circuit or field with the engine's
own elements (the N machine element, Maxwell's mutual inductance of
coaxial circles, Biot–Savart loops, the FDTD solver), and prints the
model's value next to the measured one. The measured values are the same
for every model (CODATA 2018; Fan et al. 2023):

| Quantity | Measured | Held to |
|---|---|---|
| Charge | e | exact by definition since 2019 |
| Magnetic moment | 1.001 16 μ_B | spin g/2 = 1.001 159 652 180 59(13) |
| g-factor | 2.002 319 304 36 | about 10⁻¹³ |
| Spin angular momentum | ħ/2 | |
| Electron spin resonance | 28.025 GHz/T | from g |
| Proton magnetic resonance | 42.578 MHz/T | γ_p, for the ring proton of chapter 3 |
| Size | below 10⁻¹⁸ m | scattering (Bourilkov 2001) |
| Zeeman effect | splitting eB/4πm per unit of m_l, σ± circularly polarised | Zeeman 1897 |

The last row is the one Kelvin himself used, and it is the bench's
sharpest instrument: a model has to give sharp lines with the measured
polarisation before its g-factor is even asked for.

## 9.2 Run it: Kelvin's ring

Kelvin's 1899 molecule is a rigid non-conducting ring carrying +e,
filled with an incompressible fluid carrying −e that slides round it
without friction. The example integrates twenty-four charges per half
under the Lorentz force (RK4) with the field of a coaxial coil from
Biot–Savart, and runs the two things Kelvin did with it.

```
go run ./examples/kelvin_ring
```

```
1. Induction: ring a = 5 cm, ±1 µC, tube and fluid 10 g each; coil r = 10 cm, 5 cm below, 20 kA-turns, 1 s ramps
   flux through ring 0.000698736 Wb (Biot–Savart) vs 0.000698732 Wb (Maxwell, elliptic)
   moment of momentum: tube -1.11207e-10, fluid 1.11207e-10 kg m²/s; Kelvin's corrected eAM/2π = 1.11207e-10
   fluid turns at 4.44829e-06 rad/s = Larmor eB̄/2w 4.44829e-06 rad/s
   ledger: source work -4.94683e-16 J = kinetic 4.94683e-16 J + mutual-field energy -9.89365e-16 J (residual -8.0e-16)
```

Switching the coil on gives each half of the ring a moment of momentum
eΦ/2π in opposite senses, which is Kelvin's own result once the editor's
correction to his 1899 text is applied, and the fluid settles at Larmor's
rate eB̄/2w. The ring is a mechanical picture of Larmor's theorem. The
ledger closes at 10⁻¹⁶ of the energies involved, and after switch-off
everything is at rest again.

The second part is Kelvin's Zeeman test. The ring is held by isotropic
springs in a steady field, kicked, and the spectrum of its motion and of
its radiated moment is taken over twelve orientations, for three models:
Lorentz's bound charge, the standard classical account of the normal
Zeeman effect; Kelvin's gyrostat as he wrote it, with a fixed internal
spin; and the full ring, with nothing linearised.

| Model | Splitting against orientation | Ensemble line | Circular polarisation along B, radiated |
|---|---|---|---|
| Lorentz bound charge | 0.0513 at every angle | sharp, width 0.0038 | −1.00 / +1.00 |
| Kelvin gyrostat as written | 0 to 0.0243, as cos θ | broad, width 0.0188 | +0.72 / −0.71 |
| Kelvin ring, full Lorentz force | 0.0243 to 0.0260 | sharp, width 0.0037 | −0.02 / −0.03 |

The gyrostat reproduces Kelvin's reason for discarding the model in 1899:
with a fixed spin the splitting goes as cos θ, so randomly oriented
molecules broaden the line instead of splitting it. The full ring is
new. When it tilts, the flux through it changes and the fluid's spin
follows, because the canonical angular momentum w·a²·u − eΦ/2π is
conserved, so every ring splits by the same amount and the lines are
sharp. It still fails the test: the modes circle the reflection of B in
the ring's plane rather than B itself, so the light along the field is
nearly unpolarised where Zeeman measured ±1, and the splitting is set by
the total mass of tube and fluid rather than by e/m. A ring of two
opposite charges with different charge-to-mass ratios cannot give the
normal Zeeman effect. That is the bench's first rule for every toroidal
model that follows.

## 9.3 Vortex atoms

Kelvin's 1867 vortex atom is a closed vortex ring in a perfect liquid,
and he noted that the liquid's velocity round the ring is the magnetic
field of the "corresponding ring electromagnet". The example takes that
literally: circulation plays current, velocity plays B/μ₀, and the ring's
kinetic energy is ½ρΓ²L/μ₀ with L a self-inductance. So the engine's
inductance code computes the vortex ring directly.

```
go run ./examples/vortex_ring
```

| Core a/R | Ring speed from the inductance | Kelvin's thin-ring formula | Difference |
|---|---|---|---|
| 0.20 | 0.27100 | 0.27366 | −0.97 % |
| 0.10 | 0.32798 | 0.32882 | −0.26 % |
| 0.05 | 0.38369 | 0.38398 | −0.07 % |

The difference falls as (a/R)², the order of the terms the thin-ring
formula drops. Two rings leapfrog through each other 21 times in 60 time
units with the impulse and energy drifting by 2×10⁻⁵, which also checks
that the Biot–Savart velocities and the elliptic mutual inductance agree
with each other. The analogy stops at the dynamics: a vortex ring travels
because its vorticity is carried by the flow it makes, while a current
loop with the identical field is not moved by that field. A neutral
vortex ring has no charge, so the Zeeman bench cannot be applied to it
until a charge-vortex coupling is specified.

## 9.4 Tewari's space-vortex electron (H13)

The 2018 paper describes the electron as an irrotational vortex of a
massless "space fluid" around a spherical void whose surface moves at c,
and states the charge, mass, angular momentum and magnetic moment in
terms of the void radius. The example fixes that radius from the measured
charge, as the paper does, and computes the rest.

```
go run ./examples/space_vortex
```

| Quantity | Model | Measured | Model / measured |
|---|---|---|---|
| Charge (eq. 4) | 4.803×10⁻¹⁰ statC | 4.803×10⁻¹⁰ | 1 (input) |
| Angular momentum (eq. 15) | 8.802×10⁻²⁸ erg·s | ħ/2 = 5.273×10⁻²⁸ | 1.669 |
| Magnetic moment (eq. 16) | 1.451×10⁻²⁰ erg/G | 9.285×10⁻²¹ | 1.563 |
| g-factor (eqs. 15–16) | 15/8 = 1.875 | 2.00232 | 0.9364 |
| g with the ring moment ½qvρ | 15/16 | 2.00232 | 0.468 |
| Proton charge (eq. 30) | 12 e | e | 12 |
| Fine-structure constant (book eq. 5.2) | 15π²/16 = 9.253 | 7.297×10⁻³ | 1268 |
| Hydrogen orbit speed | 0.69 c | αc | 94.6 |
| Orbit frequency against Lyman-α | 3.405×10¹⁸ Hz | 2.466×10¹⁵ Hz | 1381 |

The g-factor is the model's one real prediction, because the void radius,
the charge and the mass all cancel in it: 15/8 for any radius, 6.4 % from
a value measured to about 10⁻¹³. The hydrogen orbit gives the wrong
light, X-rays from a charge circling at 0.69 c where hydrogen's lines are
ultraviolet. The paper does not treat the Zeeman effect, so the bench's
sharpest test cannot be applied, and H13 stays *untested* with these rows
as the values it must also meet. The standard counterpart beside every
row is the Dirac electron with QED, which gives the measured values.

## 9.5 A ring of charge as one of Kron's machines (H21)

The hypothesis is that Kron's machine equations describe the ring
electron: a ring of charge q and mass m turning at θ̇ carries the current
i = qθ̇/2π, so the mechanical and electrical coordinates are one axis and
the kinetic energy ½ma²θ̇² is a magnetic energy with the inertia entering
Kron's metric tensor as an inductance ma²(2π/q)². Kron defined the metric
that way himself, as the inductances together with the rotor inertia
(1951, Epilogue eqs. 3–4), and his eq. 11 says that with no applied
voltage and no resistance the machine equation is the motion of a charged
particle in a field. The example uses the N element of chapter 3
unchanged.

```
go run ./examples/kron_electron
```

Step one is a twin test on the bench of section 9.2: Kelvin's ring as a
three-winding machine, coil, tube and fluid, with Maxwell's mutual
inductance between coil and ring and the two inertias on the charge axes.

| Quantity | Kron machine | Lorentz-force bench |
|---|---|---|
| Flux through the ring | 6.98732×10⁻⁴ Wb | 6.98736×10⁻⁴ Wb |
| Moment of momentum, tube / fluid | ∓1.11207×10⁻¹⁰ kg m²/s | ∓1.11207×10⁻¹⁰ |
| Fluid angular velocity | 4.44827×10⁻⁶ rad/s | 4.44829×10⁻⁶ |
| Kinetic energy | 4.94677×10⁻¹⁶ J | 4.94683×10⁻¹⁶ J |

The frame reproduces the bench to 10⁻⁵, the difference being the bench's
Biot–Savart flux against Maxwell's elliptic value. The same ring,
spinning at 1000 rad/s and tilting in 1 T, has Kron's torque
½iᵀ(∂L/∂θ)i equal to −μB sin θ, and the pendulum it makes has a period of
19.873 s against the exact 4K(k)/ω₀ = 19.872 s, with the odd harmonics of
the Jacobi series to 10⁻⁴. That is the standard part, and it passes.

Step two is the hypothesis: charge e, total angular momentum of the
metric ħ/2, and the rest follows.

| Ring | a [m] | v/c | μ/μ_B | g | ½LI²/m_e c² |
|---|---|---|---|---|---|
| Compton ring, electron mass circulating | 3.86×10⁻¹³ | 0.494 | 0.494 | 0.989 | 0.0014 |
| Half-Compton ring, mass circulating | 1.93×10⁻¹³ | 0.978 | 0.489 | 0.978 | 0.011 |
| Compton ring, massless (field inertia) | 3.86×10⁻¹³ | 43.6 | 43.6 | 87.2 | 10.9 |
| Massless ring fitted to the measured moment | 8.86×10⁻¹⁵ | 43.6 | 1.001 | 2.002 | 475 |
| **Measured (Dirac + QED)** | < 10⁻¹⁸ | | 1.001 16 | 2.002 32 | |

A ring whose mass circulates with its charge has
g = 1/(1 + Le²/4π²ma²) ≤ 1, so it precesses at Larmor's eB/2m, half the
electron's rate. Only a massless ring reaches the measured moment at
S = ħ/2, and it does so at 8.9 fm with v = 43.6 c and a magnetic energy of
475 electron masses; every massless thin ring with S = ħ/2 moves at
v/c = π/(2α(ln(8a/r) − 7/4)), fixed by the fine-structure constant. The
measurement that discriminates is magnetic resonance, which chapter 3
solves in the same frame: a ring electron would resonate at 14.0 GHz/T
against the measured 28.025, and a ring proton at 7.6 MHz/T against the
measured 42.6. H21 is recorded as *bounded*: the ring-as-machine
reproduces the classical ring exactly and inherits its obstacles, g ≤ 1
or v > c. Its two genuine matches are classical, the orbital Larmor
precession of the normal Zeeman effect and the measured large-amplitude
pendulum.

## 9.6 The electron as a photon on a double loop (H23)

Williamson and van der Mark (1997) take one photon of the Compton
wavelength going twice round a torus and derive the charge, moment and
spin ½ from that. The paper's own numbers, computed by the example:

| Quantity | Model | Measured |
|---|---|---|
| Charge (eq. 5) | 0.910 e | e |
| Magnetic anomaly (g − 2)/2 (eq. 12) | 0.000 962 | 0.001 160 |
| Spin | ½ from the double loop | ½ |

The anomaly is 17 % low because the charge is 9 % low. The 720° double
loop is then built as a circuit, a two-wire ring line of 200 LC sections
with its conductors crossed once: the modes have half-integer numbers of
wavelengths round the ring, and a pulse comes back inverted after one
trip and upright after two, each within 10⁻⁴ of the ladder's own theory.

What the paper assumes and does not derive is the confinement, and that
is what the field solver tests. One wavelength of light is set on the
(2, 1) torus-knot path in empty space and left alone under linear
Maxwell equations (chapter 5's `em/fdtd.go`, with its ledger closing to
2×10⁻¹³). The run is about ten minutes, so the figures below are the
example's recorded output rather than a fresh run.

| Time (periods) | 0 | ¼ | ½ | 1 | 2 |
|---|---|---|---|---|---|
| Energy within twice the ring radius | 98.9 % | 5.8 % | 1.6×10⁻⁶ % | 1×10⁻⁴ % | 6×10⁻⁷ % |
| Left the region through the absorber | 0 | 4×10⁻¹¹ % | 9.6 % | 99.65 % | 99.99 % |

Nothing holds light in a loop; it leaves at the speed of light within
half a period. The model needs all of it to stay for ever, so it needs a
confining mechanism the paper does not give, and H23 is recorded as
*bounded*. The same example then builds three electron pumps that do
work, a single-electron pump delivering exactly one electron per cycle,
a multipactor pump in a resonator and a resonant ground charge pump, each
with its ledger naming the source of every joule; those belong to chapter
6's family of gap devices and are summarised in the example's README.

## 9.7 What to keep

- The bench is fixed before any model is run: the measured values and the
  Dirac prediction appear beside every result, and the Zeeman test comes
  before the g-factor.
- A standard twin test comes first. Kron's frame had to reproduce Kelvin's
  ring to 10⁻⁵ before anything was read from it about the electron.
- A bound is a result. H21 and H23 are *bounded*, with the number that
  bounds them recorded: g ≤ 1 for a mass ring, light gone in half a period.
- Nothing is removed for failing. Each entry stays in the register with
  its status, and the values it must still meet are written down.

## Run it yourself

```
go run ./examples/kelvin_ring            # induction and the Zeeman test (a few minutes)
go run ./examples/vortex_ring
go run ./examples/space_vortex
go run ./examples/kron_electron          # -png viz/ring_field.png, -gif viz/ring_motion.gif
go run ./examples/photon_electron        # about ten minutes
go test ./examples/kelvin_ring ./examples/vortex_ring ./examples/space_vortex ./examples/kron_electron ./examples/photon_electron
```

## Sources

- W. Thomson (Lord Kelvin), "Magnetism and Molecular Rotation," *Proc.
  R. Soc. Edinb.* 22, 631–635 (1899); *Mathematical and Physical Papers*
  Vol. V, Art. 175, pp. 354–358. "On Vortex Atoms," *Proc. R. Soc.
  Edinb.* 6 (1867); *Papers* Vol. IV, Art. 1 and Art. 3.
- P. Zeeman, *Phil. Mag.* 43, 226 (1897); J. Larmor, *Phil. Mag.* 44, 503
  (1897).
- P. Tewari, "Structural relation between the vacuum space and the
  electron," *Physics Essays* 31, 108–129 (2018); *Universal Principles of
  Space and Matter*, 2017.
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley 1951,
  Epilogue eqs. 3–4, 9 and 11.
- J. G. Williamson, M. B. van der Mark, "Is the electron a photon with
  toroidal topology?", *Ann. Fond. Louis de Broglie* 22, 133 (1997).
- CODATA 2018 (E. Tiesinga et al., *Rev. Mod. Phys.* 93, 025010, 2021);
  X. Fan et al., *Phys. Rev. Lett.* 130, 071801 (2023); D. Hanneke et al.,
  *Phys. Rev. Lett.* 100, 120801 (2008); D. Bourilkov, *Phys. Rev. D* 64,
  071701 (2001) for the size bound.
- J. C. Maxwell, *Treatise*, Vol. II, Arts. 701 and 706; P. G. Saffman,
  *Vortex Dynamics*, 1992, §10.2; R. A. Nelson, M. G. Olsson, *Am. J.
  Phys.* 54, 112 (1986).
- [doc/kelvin-references.md](../doc/kelvin-references.md),
  [doc/tewari-references.md](../doc/tewari-references.md).
