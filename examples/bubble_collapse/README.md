# Bubble collapse: the sonoluminescing bubble on the flask, beside the knot

The bubble that [acoustic_drive](../acoustic_drive) stopped short of,
put on that example's flask as an equivalent circuit (CLAUDE.md §1
non-goals: the bubble returns as a cited equivalent-circuit model), and
compared with the knot's implosion in [knotted_light](../knotted_light)
(H29, H30). The question, from the user (2026-11-05): is the energy at the
centre of a collapsing bubble a large-scale picture of the knot's
implosion, and does its light carry the knot's topology?

```
go run ./examples/bubble_collapse            # about 1 min
go run ./examples/bubble_collapse -nofield   # parts 1–4 only, a second
go test ./examples/bubble_collapse           # -short skips the field solver
go run ./examples/bubble_collapse -data report_data.json   # every series behind the report
go run ./examples/bubble_collapse -report examples/bubble_collapse/viz/report.html   # rebuild the interactive report
```

## The bubble as a circuit

Acoustic domain ([doc/domains.md](../../doc/domains.md)): pressure across,
volume velocity U = dV/dt through. From the far liquid to the gas:

```
p∞(t) ── [ M(V) ∥ R_rad(V) ] ── R_visc(V) ── wall ── gas p_g(V) − 2σ/R
```

| Element | Value | Source |
|---|---|---|
| liquid inertance | M = ρ/(4πR), so ½MU² = 2πρR³Ṙ² | Rayleigh 1917 |
| radiation resistance | R_rad = ρc/(4πR²) | pulsating sphere: M ∥ R_rad = (ρc/4πR²)·jkR/(1+jkR), Kinsler et al. §7.2 |
| wall viscosity | R_visc = μ/(πR³) (= 4μṘ/R per U) | Newtonian liquid |
| gas | argon, adiabatic γ = 5/3, hard core R0/h = 8.86 | Brenner, Hilgenfeldt & Lohse 2002 |
| surface | 2σ/R | IAPWS R1-76(2014) σ = 71.69 mN/m at 300 K |

M depends on the volume, the circuit's own "charge", so the inertance is a
Kron inductance on a coordinate. With momentum Π = M·U_M,
dΠ/dt = Δp + ½(dM/dV)·U·U_M, where the ½M′ term is the Christoffel
(speed-voltage) term of Kron 1951 (Epilogue eqs. 3–4). Without the two
losses this is the Rayleigh–Plesset equation exactly
(`TestRayleighPlessetLimit`). R_rad adds the liquid's compressibility to
first order in Ṙ/c, the job Keller–Miksis does; Keller–Miksis runs beside
it as the twin. At small amplitude the element's admittance equals
1/(Z_rad + R_visc + K/jω) of linear acoustics to 1e-4
(`TestLinearAdmittance`).

**Ledger.** It closes by construction:

```
∫ −p_ac·U dt = Δ[½MU_M² + p_g(V−Vh)/(γ−1) + 4πσR² + P0·V] + ∫R_visc·U² dt + ∫Δp²/R_rad dt
   in from the sound field   =   Δ stored                 +   viscous heat  +  sound radiated
```

The derivation: d(½Π²/M)/dt = U_M·Δp from the state equation; the gas,
surface and P0 terms give −p_gU, (2σ/R)U and P0U; and
U_M·Δp = U·Δp − Δp²/R_rad.

Water: ρ, c from IAPWS-95 (`fluid.Water`); μ = 853.8 µPa·s (IAPWS R12-08).
The drive is 1.3 atm at the matched flask's peak, 25.02 kHz.

## Results

**Part 1: one acoustic cycle** (the sixth from rest):

| | Circuit | Keller–Miksis |
|---|---|---|
| R_max | 32.60 µm | 32.73 µm |
| R_min | 0.730 µm | 0.769 µm |
| minimum at phase | 0.5170 | 0.5181 |
| peak wall Mach | 0.78 | 0.52 |

- Peak gas pressure is 2.4×10⁹ Pa and the adiabatic peak temperature is
  15,000 K, with R_max/R_min = 45.
- **Ledger:** 14.29 nJ in from the sound field per cycle =
  1.04 nJ viscous heat + 13.25 nJ radiated sound, with Δ stored = 0 over
  the periodic cycle. Residual 6×10⁻¹⁰.
- 91 % of the radiated sound leaves within ±10 ns of the minimum: the
  collapse turns a 25 kHz input into a nanosecond outgoing pulse.
- Energy density goes from 3.9 J/m³ in the sound field to 9×10⁹ J/m³ in
  the collapsed bubble, a concentration of 2×10⁹.

The two models agree on the slow quantities (R_max to 0.4 %, the
minimum's phase to 0.1 %). They differ at the collapse, where Ṙ/c is no
longer small. The circuit model is first order in Ṙ/c. Keller–Miksis keeps
(1 − Ṙ/c) factors that slow the wall more (Mach 0.52 against 0.78, R_min
5 % larger). Neither is exact there; Brenner et al. discuss the
alternatives (Gilmore, full gas dynamics).

**Part 2: the bubble as a load on the drive chain.** Its fundamental
admittance Y₁ = −U₁/P₁ is put at node p of the matched water chain (a
resistor and a capacitor). This is a harmonic balance at the fundamental.
It is exact for power, because only the fundamental of U does work against
the sinusoidal centre pressure: ½Re(Y₁)Pa² = 357.66 µW equals the cycle
work × f.

| | No bubble | With bubble |
|---|---|---|
| peak | 25.025 kHz | 25.024 kHz |
| pressure per volt | 69.40 kPa/V | 68.46 kPa/V |
| generator EMF for 1.3 atm | 1.898 V | 1.924 V |
| source power | 26.18 mW | 26.54 mW |
| into the liquid (mode) | 6.957 mW | 6.957 mW |
| into the bubble | — | 0.358 mW |

The bubble takes 5 % of what the flask's own losses take. kronspice's
ledger closes. `export/bubble_on_flask.cir` is the chain with the bubble
element (RBUBBLE, CBUBBLE), with its KiCad sheet and the MNA vs torn
results.

**Part 3: reversibility.** The knot's implosion is its explosion run
backwards: Maxwell's equations are time-symmetric, and H30 fed the inward
half through a surface. Here the bubble runs from R_max through the
collapse and back out. Then its momentum is flipped and the drive is run
backwards:

| Losses | \|ΔR\|/R on return |
|---|---|
| none | 2.5×10⁻¹² |
| viscosity only | 0.065 |
| radiation only | 0.90 |
| both | 0.80 |

Without losses the bubble is a conservative (Hamiltonian) system and
retraces its collapse, like the knot. The radiated sound makes it
irreversible, because the outgoing pulse does not come back. Running a
bubble backwards would need that pulse to converge on it again, which is
the job H30's Huygens surface did for the knot.

**Part 4: the light.** Free–free (bremsstrahlung) emission from argon
ionized per the Saha equation, assuming a uniform core where every photon
escapes (`light.go`; Rybicki & Lightman eq. 5.15b; the standard account is
Hilgenfeldt, Grossmann & Lohse, *Nature* 398, 402, 1999):

- 12.8 fJ per flash, about 6×10⁴ photons, with a full width at half
  maximum of 180 ps;
- 0.6 % ionization at the peak;
- 9×10⁻⁷ of the collapse energy.

The measured flashes are 10⁵–10⁷ photons wide 60–250 ps (B. Gompf et al.,
*PRL* 79, 1405, 1997). The estimate is low by one to two orders in
photons and right in width. It leaves out electron–neutral bremsstrahlung,
recombination, heat loss and the lowered ionization potential, all of
which the full account includes. The light is far below the other ledger
terms and is not taken out of the dynamics.

**The flash spectrum against measurement.** The same emission, resolved in
frequency (Rybicki & Lightman eq. 5.14b) and integrated over the flash:

| Quantity | Model | Measured | Agreement |
|---|---|---|---|
| Width, UV band (250–300 nm) / red band (600–700 nm) | 164 / 176 ps | equal (Gompf et al. 1997) | yes |
| Spectral shape | peaks near 475 nm, falls into the UV | rises into the UV to water's edge near 200 nm | no |
| Planck-fit temperature of the shape | 6,500 K | of order 10⁴ K and above | too cool |

The width test agrees. The colour does not: an optically thin gas at the
model's 15,000 K peak emits most in the blue-green. The processes left out
above are where the standard account finds the ultraviolet. None of the
shortfall points to a non-thermal source.

**Part 5: does the flash carry the knot's topology?** A knot of field
energy carries net optical helicity, while thermal light from random
currents carries none. The bench measurement is the whole flash's degree of
circular polarization,

```
P_c = ΣV/ΣI,  V = −2 Im(Ê_θ*·Ê_φ),  I = |Ê_θ|² + |Ê_φ|²
```

summed over every frequency and over 600 probes on a sphere (V > 0 is spin
along the direction of travel). The sources run in `em/fdtd.go`, scaled to
a central wavelength λ₀ = 450 nm at 8 cells: the collapsed core,
R_min ≈ 0.79 µm, is 1.75 λ₀, and the probes sit at 6 λ₀.

| Source | P_c | Field ledger residual |
|---|---|---|
| hopfion of knotted_light (exact field at the same probes: +0.9999) | +0.9996 | 1.6×10⁻¹³ |
| chiral dipole pair, p ∥ m = c·p in quadrature | +0.980 | 6×10⁻¹⁵ |
| rotating electric dipole (circular on its axis) | +0.0007 | 4×10⁻¹⁵ |
| thermal core: 300 random current elements, 6 seeds | mean +0.006, rms 0.018 | ≤ 4×10⁻¹⁴ |

- **The knot is one helicity.** Rañada's hopfion in the Bateman form
  satisfies i∂F/∂t = ∇×F, which is the condition for one helicity. The
  measurement sees it at +0.9996, and the exact field sees +0.9999.
- **An electric dipole alone cannot carry helicity**, even when it is
  circular on its axis: the two hemispheres cancel. Helicity needs
  electric and magnetic parts together, since it is the charge of
  electric–magnetic duality (Calkin, *Am. J. Phys.* 33, 958, 1965; Cameron,
  Barnett & Yao, *New J. Phys.* 14, 053050, 2012). Random hot-gas currents
  are electric currents.
- **The thermal core gives P_c consistent with zero.** The seed-to-seed
  spread falls as the number of independent radiators grows. A real flash
  of 10⁵–10⁷ photons has a counting floor of 1/√N ≈ 3×10⁻³ to 3×10⁻⁴.

**Verdict.**
- **Where the analogy holds:** the bubble is a driven implosion and
  explosion with its inward supply booked, one cycle per acoustic period.
  This is the structure H30 found for the pulsing knot, with the flask's
  standing wave as the "declared medium".
- **Where it fails:** the collapse is a nonlinear, irreversible flow of
  matter, stopped by gas pressure. Its flow is spherically symmetric and
  carries no topology. Standard theory predicts thermal light with P_c = 0
  within counting noise.
- **The bench test:** if the flash came from a knotted field, P_c would be
  of order 1, as for the hopfion. A polarimeter on single flashes
  separates the two by three orders of magnitude.

## Files

| File | Holds |
|---|---|
| `bubble.go` | the equivalent circuit, its ledger, Keller–Miksis twin, Dormand–Prince integrator, reversal |
| `chain.go` | acoustic_drive's matched chain with the bubble's admittance at node p |
| `light.go` | Saha ionization and free–free emission, flash energy and width |
| `helicity.go` | probe sphere, Stokes V/I, the four field-solver sources |
| `hopfion.go` | Rañada's knot, copied from knotted_light |
| `viz/report.html` | the interactive report: 3D animation of the bubble beside the knot, polarization globe, field slices, every chart; open it in a browser (needs internet for three.js and fonts) |
| `viz/report_template.html` | its template, embedded in the program; `-report` fills in the data |
| `data.go` | `-data`: the series for the report (cycle, ledger, collapse, spectrum, sweeps, reversal paths, per-probe polarization, field slices, the knot's field lines from t = −3 to 3) |
| `export/` | `bubble_on_flask` netlist, KiCad schematic, results sheet |
| `run_output.txt` | full output of `go run` |

## Tests

| Test | Checks | Held to |
|---|---|---|
| `TestRayleighPlessetLimit` | lossless circuit = Rayleigh–Plesset | 1e-9 |
| `TestLinearAdmittance` | small-amplitude Y₁ = linear acoustic admittance | 1e-4 |
| `TestCycleLedgerAndTwin` | ledger; ½Re(Y₁)Pa² = cycle work × f; Keller–Miksis R_max, R_min, phase; flash plausible | 1e-8; 1e-6; 1 %, 10 %, 0.5 % |
| `TestChainWithBubble` | bubble element takes ½Re(Y₁)Pa², mode heat unchanged, kronspice ledger | 1e-3; 1e-9 |
| `TestReversal` | lossless return, lossy failure | 1e-9; > 1e-2 |
| `TestStokesSign`, `TestIonization` | sign of V; Saha limits | exact |
| `TestHelicity` | hopfion = exact, chiral > 0.9, rotating and thermal ≈ 0, ledgers | 0.01; 1e-10 |
| `em.TestFDTDSharedEdgeSources` | two sources on one edge act as their sum, ledger closed (the solver fix this example needed) | 1e-10 |

## Sources

- M. P. Brenner, S. Hilgenfeldt, D. Lohse, "Single-bubble
  sonoluminescence," *Rev. Mod. Phys.* 74, 425 (2002) — the bubble,
  argon's hard core, parameters, the light.
- Lord Rayleigh, *Phil. Mag.* 34, 94 (1917); M. S. Plesset, *J. Appl.
  Mech.* 16, 277 (1949); J. B. Keller, M. Miksis, *J. Acoust. Soc. Am.* 68,
  628 (1980).
- L. E. Kinsler, A. R. Frey, A. B. Coppens, J. V. Sanders, *Fundamentals
  of Acoustics*, 4th ed., Wiley, 2000, §7.2 (pulsating sphere).
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley, 1951,
  Epilogue eqs. 3–4.
- S. Hilgenfeldt, S. Grossmann, D. Lohse, *Nature* 398, 402 (1999);
  G. B. Rybicki, A. P. Lightman, *Radiative Processes in Astrophysics*,
  Wiley, 1979, eq. 5.15b; NIST Atomic Spectra Database (argon).
- B. Gompf, R. Günther, G. Nick, R. Pecha, W. Eisenmenger, *Phys. Rev.
  Lett.* 79, 1405 (1997) — flash widths.
- IAPWS R1-76(2014) (surface tension); IAPWS R12-08, M. L. Huber et al.,
  *J. Phys. Chem. Ref. Data* 38, 101 (2009) (viscosity).
- M. G. Calkin, *Am. J. Phys.* 33, 958 (1965); R. P. Cameron,
  S. M. Barnett, A. M. Yao, *New J. Phys.* 14, 053050 (2012) — optical
  helicity and duality.
- A. F. Rañada, *Lett. Math. Phys.* 18, 97 (1989); H. Kedia et al., *PRL*
  111, 150404 (2013) — the knot.
