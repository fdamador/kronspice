# Evolving a microstrip pair from Maxwell's equations

Two parallel microstrip traces on FR-4 are evolved for a 50 Ω
characteristic impedance and near-end crosstalk no worse than −30 dB at
100 MHz, on the least board width. Nothing about the traces is looked up:
for every candidate the per-unit-length L and C matrices are extracted from
the cross-section by the quasi-static field solver of milestone M9i
([doc/M9i.md](../../doc/M9i.md); CLAUDE.md §3c), the pair becomes a
coupled-line ladder in the circuit, and an AC analysis measures the
crosstalk onto the quiet trace.

```
go run ./examples/pcb_layout [-gens 12] [-pop 20]
go test ./examples/pcb_layout          # -short skips the GA run
```

## The construction

- **Extraction** (`em.Section`): the cross-section (strips of zero
  thickness at the top of a 1.6 mm substrate, ε_r = 4.4, air above, a box
  twelve substrate thicknesses wide) is solved as div(ε grad φ) = 0 with
  each conductor driven in turn, giving the 2×2 capacitance matrix C; the
  same solve with ε = 1 gives C₀ and L = μ₀ε₀C₀⁻¹.
- **Circuit** (`em.CoupledLadder`): the 10 cm pair as 20 sections of the
  coupled telegrapher's equations, mutual inductance as a K element and
  mutual capacitance as a capacitor between the traces. The aggressor is
  driven through Z₀ and both traces are terminated in Z₀ at both ends.
- **Measurement**: near-end crosstalk is |v(v_in)/v(a_in)| in dB from a
  one-point AC analysis at 100 MHz.
- **Genome**: trace width w (0.3–6 mm) and gap s (0.2–10 mm). Objectives,
  weighted sum: |Z₀ − 50|/50, crosstalk above −30 dB, and 0.01 × board
  width (2w + s) in units of 10 mm. Each candidate's extraction is
  memoized across the objectives.

Limits: zero-thickness strips; quasi-static (TEM) only, so dispersion and
radiation are left out; the ladder's own λ/10 rule sets the section count.

## Results

Default run (population 20, 12 generations, seed 1):

```
evolved (218 evaluations): w 2.766 mm, s 3.020 mm, board 8.55 mm
  Z₀ 49.82 Ω, ε_eff 3.257, L₁₁ 301 nH/m, C₁₁ 121.3 pF/m, K = L₁₂/L₁₁ 0.0857, NEXT -33.1 dB at 100 MHz
```

The shorter run recorded in doc/M9i.md (`-gens 8 -pop 16`, 114
evaluations) settled on w 2.769 mm with a wider gap, s 8.913 mm, at
Z₀ 49.86 Ω and NEXT −50.8 dB; the board-width weight is small, so either
satisfies the specification.

Coupling falls with the gap, as the test checks for a 3 mm trace:

| Gap s | K = L₁₂/L₁₁ | NEXT at 100 MHz |
|---|---|---|
| 0.5 mm | 0.287 | −20.8 dB |
| 2 mm | 0.128 | −29.1 dB |
| 6 mm | 0.030 | −43.3 dB |

The extraction itself is validated in `em/extract_test.go` against the
analytic coax, Cohn's stripline, the two-wire line and Hammerstad–Jensen's
microstrip (doc/M9i.md). The GA test runs without a per-evaluation timeout
and fails if any candidate timed out, so the evolved result is the same on
every machine (CLAUDE.md §1 goal 5).

## Files

| File | What it holds |
|---|---|
| `main.go` | the cross-section, the coupled ladder, the crosstalk measurement and the GA problem |
| `main_test.go` | coupling falls with the gap; the GA meets the specification |
| `export_test.go` | the w = 3 mm, s = 2 mm pair as a netlist, KiCad schematic and results sheet in `export/` (torn and MNA equal to 1.2e-15) |

## Sources

- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed. (2008) — per-unit-length parameters and the coupled-line equations (`reference/`).
- E. Hammerstad, Ø. Jensen, "Accurate models for microstrip computer-aided design," *IEEE MTT-S* (1980) — the microstrip check in `em/extract_test.go`.
- doc/M9i.md — the extraction method and its validation cases.
