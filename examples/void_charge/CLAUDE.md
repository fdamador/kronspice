# CLAUDE.md — examples/void_charge (the charged void, proposed H51)

Notes for Claude sessions working in this example. The repo's root CLAUDE.md
governs everything here:

- the standard prediction goes beside every hypothesis;
- the energy ledger closes (§3c);
- every formula cites its source and equation;
- hypothetical results are labelled.

A new example also needs a row in `examples/README.md`, its own `README.md`,
and a citation in a `book/` chapter. The row and the citation are not done
yet, because the user has uncommitted edits in both. Every Go file starts
with the SPDX dual-license header.

## The question (user, 2026-11-14)

QED's point electron (no size, no structure) resembles the void of a
nonlinear ether. Assign the void electrical properties (negative charge,
vacuum pressure) and test them experimentally. Re-read Tewari's book, which
covers this (UPSM 2017 §2.3, §3.1, §8; `reference/tewari_upsm_2017.pdf`).

## Rule used

The charge value is a *suggestion* (the standard model puts it in by hand),
so assigning it is allowed. What is tested is where the charge sits, whether
it is exact, and what force and light follow from it.

## Files

- `field.go`: the void_shells solver, copied (GP ether + m = ½ light modes).
- `charge.go`: measured constants, the hydrogen contact formula, LEP2 bound,
  Uehling running, and Tewari's numbers.
- `measure.go`:
  - charge placements (missing ether, wall, light);
  - ⟨r²⟩, form factor, enclosed charge;
  - pair and stretched states;
  - ringing readers.
- `main.go`: Parts 1–6.
- `main_test.go`: the unit and physics tests.
- `viz/plot_charge.py` → graph, collision image, mp4/gif. Frames live in
  `viz/frames/`, which is not in git.
- `viz/make_report.py` → `void_charge_report.pdf`.

## Register entry (proposed H51; add to root §1 when merged)

| H51 | The void as the electron's charge: give the void of a nonlinear ether electrical properties (charge −e, pressure) and test them electrically (user hypothesis, 2026-11-14), starting from Tewari's charged void (H13; UPSM §2.3, §3.1, §8) and the light bubble of H50 | `examples/void_charge`: Tewari's formulas against measurement; the void_shells bubble (ε = 300, V0 = 20) with charge placed where the ether is missing, on the wall, or with the light; hydrogen contact shift and LEP2 form factor; missing ether across ε, N photons and in real time; e⁺e⁻ and e⁻e⁻ voids, statically against separation and colliding in real time; ringing after a 3 % disturbance | refuted as built — Tewari: non-conservative field (2.68 keV per lap), 2p spread 4.08 eV by spin orientation, r_e = 1.04 ƛ shell gives 9.9 GHz extra 2S shift against ±9 kHz, LEP2 form factor ≤ 2.4e-6, annihilation 0.8 × 2mc² in one shell of wavelength r_e (measured: 2 × 511 keV back-to-back photons); computed bubble: charge rms 3–6 ƛ in every placement (246–1121 GHz extra 2S shift), seen charge falls inside (QED's rises); missing ether is not in units (2 photons 1.86×, ε-dependent) and flickers 0.961–1.035 after a 3 % jolt; only the light's winding m = ±½ is exact but nothing feels its sign; two voids attract alike (e^{−1.43D}, no 1/r, no sign), e⁺e⁻ gather both photons in one void while the other collapses into a sound shell, 17.2 % of 2M freed as sound, light kept (285 against escape 1342): no annihilation; soft internal modes below a keV; EDM 0 passes; standard prediction: pointlike charge, Darwin 0.866 ƛ, α rising at short distance, e⁺e⁻ → 2γ, no excited electron |
