# Magnetic resonance as Kron's rotating frame

NMR for the proton and ESR for the electron, solved by the circuit engine
with one of Kron's machines. A companion to [kron_electron](../kron_electron)
(H21): the same N element, now with the magnetization as its currents.

```
go run ./examples/kron_nmr -gif viz/nmr.gif
go test ./examples/kron_nmr
```

## The construction

Bloch's equations

```
dM/dt = γ·M × B − (Mx, My, 0)/T2 − (0, 0, Mz − M0)/T1
```

are a rotation of the magnetization about the field. A rotation of currents
at a fixed rate is Kron's frame II(b) (*Equivalent Circuits of Electric
Machinery*, 1951, Epilogue eq. 9): constant inductances and a
skew-symmetric speed matrix G, so that L·di/dt = −ω·G·i. Here:

- the three windings x, y, z carry Mx, My, Mz as their currents (L = 1 nH,
  so the winding voltages ωL·M are volts);
- G = L·γ·[B]ₓ, the matrix of M × B, and the shaft turns at ω = 1. A skew
  G has zero torque iᵀ·G·i, so the precession takes no shaft power;
- the NMR rotating frame is Kron's rotating reference frame (CLAUDE.md
  §3b, Park's transform): at the RF frequency ω_rf the RF field B₁ is a
  constant x-field and the static field becomes B₀ − ω_rf/γ, so the whole
  problem has a constant G and resonance is where the effective z-field
  vanishes;
- relaxation is resistance on the windings, T2 = L/R₂ on x and y, T1 = L/R₁
  on z with a series source R₁·M₀ that restores the equilibrium; the heat
  is on the ledger;
- a fourth, open winding coupled to x by M_p is the pickup coil, whose
  voltage M_p·dMx/dt is the spectrometer's signal.

## Results (proton at 1 T, T1 = 5 µs, T2 = 2 µs chosen short so the runs stay small)

| Check | Machine | Expected |
|---|---|---|
| Free precession | 42.5766 MHz | γB₀/2π = 42.5775 MHz (2×10⁻⁵, the trapezoidal phase error at 400 points per period) |
| Transverse decay at 4T2 | 0.01832 | e⁻⁴ = 0.01832 |
| Mz recovery at 4T2 | 0.79810 | 1 − e^(−4T2/T1) = 0.79810 |
| Rabi nutation, B₁ = 1 mT | period 23.487 µs, Mz − cos(γB₁t) within 2×10⁻⁴ | 2π/γB₁ |
| π/2 pulse (5.87 µs) | M = (0, 1.0000, 0) | onto the transverse plane |
| Pickup peak | 0.02676 V | M_p·γB₀·M₀ = 0.02675 V |
| Continuous-wave line, B₁ = 10 µT | absorption, dispersion and Mz at seven detunings | Bloch's steady state, to 5 digits |
| Ledger | source − heat − Δstored = −3×10⁻¹⁴ J on 10⁻⁹ J | |

**Resonance frequencies at 1 T**, measured against the H21 ring with the
particle's charge and mass (g = 1, γ = e/2m):

| | Measured | Ring (g = 1) |
|---|---|---|
| Proton NMR | 42.578 MHz (g_p = 5.5857) | 7.623 MHz |
| Electron ESR | 28.025 GHz (g_e = 2.0023) | 13.996 GHz |

This is the bench measurement that discriminates H21: a ring electron
would resonate at half the ESR frequency.

## Animation

`viz/nmr.gif`: the magnetization vector in the rotating frame (left) through
the π/2 pulse and the decay that follows, with the pickup coil's signal
drawn as it arrives (right). The frames are the solver's time points,
the decay run converted to the rotating frame by Park's rotation.

## Against existing measurements

The resonance frequencies are measured inputs, so the match is exact by
construction: the proton's 42.577 MHz/T is the CODATA value of γ_p/2π
(from the shielded-proton measurements in water that underlie it), and
the electron's 28.025 GHz/T follows from the measured g-factor. What the
machine adds is the dynamics, and those are Bloch's 1946 equations, whose
predictions (precession, T1 and T2 relaxation, the Lorentzian absorption
and dispersion lines, saturation) were confirmed in the first NMR
experiments of Purcell, Torrey and Pound and of Bloch, Hansen and Packard
(both 1946) and in every spectrometer since; electron spin resonance was
first observed by Zavoisky in 1945. The scaled relaxation times and
pickup units are not measurements (see the results table).

## Sources

- F. Bloch, "Nuclear induction," *Phys. Rev.* 70, 460–474, 1946; F. Bloch,
  W. W. Hansen, M. Packard, "Nuclear induction," *Phys. Rev.* 69, 127,
  1946; E. M. Purcell, H. C. Torrey, R. V. Pound, "Resonance absorption
  by nuclear magnetic moments in a solid," *Phys. Rev.* 69, 37–38, 1946.
- E. Zavoisky, "Spin-magnetic resonance in paramagnetics," *J. Phys. USSR*
  9, 211, 1945 — the first ESR.
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley 1951,
  Epilogue (`reference/`).
- CODATA 2018 (E. Tiesinga, P. J. Mohr, D. B. Newell, B. N. Taylor, *Rev.
  Mod. Phys.* 93, 025010, 2021): γ_p, γ_e, g_p, g_e; the electron g from
  D. Hanneke, S. Fogwell, G. Gabrielse, *Phys. Rev. Lett.* 100, 120801,
  2008, and X. Fan et al., *Phys. Rev. Lett.* 130, 071801, 2023.
