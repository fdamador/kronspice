# Skyrme knots: the nucleon (B = 1) and helium-4 (B = 4)

The knot picture of H31, computed here rather than quoted. In Skyrme's
model the proton and the neutron are one knot of the pion field (winding
number B = 1) in two states of internal rotation; helium-4 is one knot of
winding 4, with cubic symmetry.

```
cd examples/nucleon_circuits/skyrme
python3 skyrme.py      # solves both knots, writes skyrme.json (a few seconds)
python3 render.py      # fig1–fig5 from skyrme.json
```

**Method.** Hedgehog profile F(r), F(0) = π, F(∞) = 0, on a stretched grid
r = 1.5·tan(πx/2) (400 points), energy minimized by L-BFGS-B; units
1/(eFπ) for length and Fπ/e for energy (Adkins, Nappi & Witten 1983).
The B = 4 cube uses the rational-map ansatz of Houghton, Manton &
Sutcliffe (1998), R(z) = (z⁴ + 2√3 i z² + 1)/(z⁴ − 2√3 i z² + 1), angular
integral I = 20.65. e and Fπ are fitted to the nucleon and Δ masses, as
ANW did.

## Results (`run_output.txt`)

| | This computation | ANW 1983 | Measured |
|---|---|---|---|
| B = 1 energy (Fπ/e) | 36.462 | 36.5 | — |
| e, Fπ | 5.437, 129.1 MeV | 5.45, 129 MeV | Fπ = 186 MeV (this convention) |
| Isoscalar radius | 0.595 fm | 0.59 fm | 0.72 fm |
| μ_p | +1.880 μN | +1.87 | +2.793 |
| μ_n | −1.324 μN | −1.31 | −1.913 |
| \|μ_p/μ_n\| | 1.419 | 1.43 | 1.460 |
| B = 4 energy per winding (12π² units) | 1.1366 | HMS 1.137 | — |
| Helium-4 binding | 267 MeV (7.7 %) | — | 28.3 MeV (0.75 %) |

The B = 1 knot reproduces ANW to 0.1–1 %. The cube is lighter per
winding than four separate knots, so helium-4 holds together, but about
ten times too tightly: the classical knot over-binds. This is the known
shortcoming of the original model; lightly bound Skyrme models (Gillard,
Harland & Speight 2015) and quantized vibrations bring it down.

## Figures

- `fig1_profiles.png`: F(r) and the radial winding density for B = 1 and B = 4.
- `fig2_charge.png`: proton and neutron charge densities, ρ = ½ρ_S ± ½ρ_V
  (ρ_S the winding density, ρ_V from the moment of inertia density, each
  normalized). The neutron has a positive core and a negative skin, the
  sign pattern electron scattering measures.
- `fig3_compare.png`: one light ring, quark light rings, Skyrme knot and
  measurement for μ_p, μ_n and size.
- `fig4_helium.png`: energy per winding and binding against measurement.
- `fig5_knots3d.png`: surfaces of high winding density in 3-D, a ball for
  B = 1 and a cube for B = 4; the integrated winding inside the plotted
  box is 0.95 and 3.99.

## Sources

- T. H. R. Skyrme, *Proc. R. Soc. A* 260, 127 (1961).
- G. S. Adkins, C. R. Nappi, E. Witten, *Nucl. Phys. B* 228, 552 (1983).
- C. J. Houghton, N. S. Manton, P. M. Sutcliffe, *Nucl. Phys. B* 510, 507 (1998).
- M. Gillard, D. Harland, M. Speight, *Nucl. Phys. B* 895, 272 (2015).
- Helium-4 binding: AME 2020 (28.296 MeV).
