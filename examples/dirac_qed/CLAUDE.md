# CLAUDE.md — examples/dirac_qed (the standard electron, benchmark)

Notes for Claude sessions working in this example. The repo's root
CLAUDE.md governs everything here. This example is **not a hypothesis**:
it is the standard prediction (Dirac + QED) computed in kronspice, the
benchmark row of the electron ledger. No register entry. New example ⇒
a row in `examples/README.md` and a citation in a `book/` chapter (not
yet done: the user has uncommitted edits in both). Every Go file starts
with the SPDX dual-license header.

## What is computed here and what is taken from the literature

- **Computed:** the Dirac spectrum and spin (Part 1); Landau levels and
  g = 2 (Part 2); hydrogen n = 1, 2 by shooting (Part 3); Schwinger's
  one-loop integral, ζ(3) and Petermann's A⁽⁴⁾ (Part 4); the trap
  circuits (Part 5).
- **Literature, cited in `qed.go`:** A⁽⁶⁾, A⁽⁸⁾, A⁽¹⁰⁾, the mass-dependent
  terms, the hadronic and weak terms, α (Rb, Cs) and the measured a_e.
- **Illustrative:** the trap's effective electrode distance d = 3.5 mm and
  the detector tank (L = 100 nH, Q = 1000).

## Files

`dirac.go`, `qed.go`, `trap.go`, `main.go`, `main_test.go`,
`viz/plot_dirac.py`, `viz/make_report.py`.
