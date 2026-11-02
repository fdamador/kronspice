# Hydrogen's ground state held up by the zero-point field (H24)

H. E. Puthoff, "Ground state of hydrogen as a zero-point-fluctuation-
determined state," *Physical Review D* 35, 3266 (1987), works in
**stochastic electrodynamics** (SED). SED treats empty space as filled
with a real, random classical electromagnetic field carrying ½ħω of energy
per mode, the "zero-point field".

A classical electron circling a proton radiates, and should spiral into the
proton within about 16 picoseconds. Puthoff shows that on the right orbit
it **absorbs from the zero-point field exactly what it radiates**, so it
stays. The right orbit is the one with m·ω·r² = ħ — Bohr's condition, which
gives the Bohr radius and −13.6 eV.

This example is beside [kron_schrodinger](../kron_schrodinger) (the
standard quantum-mechanical hydrogen atom as Kron's network) and the
other electron models on the test bench. It is hypothesis **H24** in the
CLAUDE.md register.

```
go run ./examples/zpf_hydrogen -viz examples/zpf_hydrogen/viz   # about 15 seconds
python3 examples/zpf_hydrogen/viz/render.py
go test ./examples/zpf_hydrogen
```

## 1. The paper's numbers

With CODATA 2018 values:

| Quantity | Result | Compare |
|---|---|---|
| Orbit radius from m·ω·r² = ħ and Coulomb (eq. 20) | 5.291772×10⁻¹¹ m | Bohr radius a₀, the same |
| Energy | −13.6057 eV | Bohr and Schrödinger for an infinitely heavy proton, the same (measured ionization energy 13.598 eV includes the proton's recoil) |
| Power absorbed from the zero-point field (eq. 18) | 4.669288×10⁻⁸ W | equal to the radiated power (eq. 19) |
| Orbit frequency | 6.58×10¹⁵ Hz | |
| Collapse time with no zero-point field | 1.56×10⁻¹¹ s | the classical problem Puthoff answers |

## 2. Does the balance hold the orbit?

Averaged as in the paper, on Coulomb orbits of other radii:

| Radius | Absorbed ÷ radiated |
|---|---|
| 0.5 a₀ | 1.414 — gains energy, grows |
| 0.9 a₀ | 1.054 |
| **1.0 a₀** | **1.000** |
| 1.1 a₀ | 0.953 |
| 2.0 a₀ | 0.707 — loses energy, shrinks |

The ratio is ħ/(mωr²), which falls as 1/√r. So the averaged balance pulls
the orbit back to a₀, with a time constant of 9.3×10⁻¹¹ s, about 600,000
orbits. Picture: left half of `viz/zpf_hydrogen.png`.

## 3. The paper's oscillator in a simulated zero-point field

Puthoff builds the orbit from two one-dimensional oscillators (eqs. 4–17).
Here that oscillator is simulated directly:

- a random field with the zero-point spectrum is made by random phases
  (eq. 2) and a Fourier transform;
- the radiating oscillator is integrated in it, 16 independent runs.

| | Simulated | Expected for this field | Formula |
|---|---|---|---|
| Mean energy | 0.5037 ± 0.0009 ħω₀ | 0.5036 | ½ħω₀ |
| Power absorbed | 0.005076 | 0.005076 | eq. 16: 0.005 |
| Power radiated | 0.005077 | | |

Absorbed and radiated agree to 0.03 %, and the oscillator holds ½ħω₀ —
the known SED result for a harmonic oscillator.

Declared choices:

- **Units:** scaled, m = ω₀ = ħ = 1.
- **Damping:** Γω₀ = 0.01 instead of the electron's 2.6×10⁻⁷ at the
  hydrogen orbit frequency, so a run covers 500 relaxation times. Every
  result is compared with the formula for that damping.
- **Field cut:** the field stops at 2ω₀. The full zero-point spectrum
  makes a classical particle's kinetic energy diverge.
- **Radiation reaction:** used in its reduced (Landau–Lifshitz) form,
  which is equal at resonance.

## 4. The same oscillator as a kronspice circuit

The circuit: a series tank of L = 1 H, C = 1 F and R = 0.01 Ω, driven by
the resistor's own **zero-point noise**. By Callen and Welton's
fluctuation–dissipation theorem (1951, cited in Puthoff's footnote 13), a
resistor at absolute zero has a voltage noise of ħωR/π per unit angular
frequency.

| | kronspice | Expected |
|---|---|---|
| Energy stored in the tank | 0.4916 ± 0.0035 ħω₀ | 0.4935 |
| Ledger | noise source 645.0 J = resistor 639.3 J + stored 4.6 J | closes to 0.2 % |

ħ is scaled to 1 J·s; the circuit is linear, so energies scale with ħ.
The expected value allows for kronspice joining the noise waveform's
points by straight lines.

One kronspice lesson from this run: after every corner of a stepwise
(PWL) waveform kronspice takes a short backward-Euler step (standard SPICE
practice). With 65,536 corners, that step damped the tank about as much as
its resistor did, and the energy came out 7 % low. A maximum step of 0.02 s
makes the extra damping negligible.

**Energy.** The noise source puts in what the resistor takes out. Net, the
zero-point field supplies nothing here, as in Puthoff's balance.

## 5. Where this stands

- **Bohr-theory level.** Puthoff's result is a single circular orbit with
  angular momentum ħ, as he says. The quantum-mechanical ground state has
  zero orbital angular momentum and a spread of radii: most probable a₀,
  mean 1.5 a₀. Extending SED to that state was "yet to be achieved" in
  1987 (his footnote 17).
- **Later full simulations of SED hydrogen:**
  - Cole & Zou (*Phys. Lett. A*, 2003) report the Schrödinger
    ground-state probability density.
  - Nieuwenhuizen & Liska (*Phys. Scr.* T165, 014006, 2015), with longer
    runs, find that the atom eventually ionizes itself.
- **Why the numbers match.** Eq. 20 is Bohr's condition, so the energy
  and size agree with Bohr and Schrödinger. The test of SED is everything
  beyond that.
- **Energy.** The zero-point exchange is zero net: absorbed equals
  radiated. Nothing here is an energy source.

## Files

| File | What it holds |
|---|---|
| `puthoff.go` | the paper's eqs. 18–20 with CODATA values, the balance around the orbit, collapse time |
| `oscillator.go` | the oscillator in a simulated zero-point field; the FFT |
| `circuit.go` | the kronspice tank in zero-point resistor noise |
| `main.go` | the report |
| `main_test.go` | Bohr numbers, stability ratio, oscillator energy and balance, tank energy and ledger, FFT |
| `viz/render.py` | `zpf_hydrogen.png` |

## Sources

- H. E. Puthoff, *Phys. Rev. D* 35, 3266 (1987).
- T. H. Boyer, *Phys. Rev. D* 11, 790 (1975) — SED and the harmonic oscillator in zero-point radiation.
- H. B. Callen, T. A. Welton, *Phys. Rev.* 83, 34 (1951) — fluctuation–dissipation theorem, zero-point noise of a resistor.
- D. C. Cole, Y. Zou, "Quantum mechanical ground state of hydrogen obtained from classical electrodynamics," *Phys. Lett. A* (2003), doi:10.1016/j.physleta.2003.08.022, arXiv:quant-ph/0307154.
- T. M. Nieuwenhuizen, M. T. P. Liska, "Simulation of the hydrogen ground state in stochastic electrodynamics," *Phys. Scr.* T165, 014006 (2015), arXiv:1502.06856.
- CODATA 2018 recommended values.
