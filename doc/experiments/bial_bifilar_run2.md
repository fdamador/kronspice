# Bi–Al capacitor element with a bifilar coil: first experiment

Generated 2026-10-25 07:58 UTC by `examples/bial_bifilar` (seed 1; 5 weighted-sum runs per variant, population 40, 30 generations).

## Set-up

A coil wound bifilar (Tesla, US 512,340) in parallel with a bismuth–aluminium element, self-driven: the element is the only energy source. Evolution trades the tank's quality factor Q at resonance against the energy stored per volume (coil + element) once the element's EMF is on. Each non-dominated design is also simulated with an ordinary coil of the same geometry.

| Assumption | Value | Source |
|---|---|---|
| Bi–Al Seebeck coefficient | 75.5 µV/K (S_Al − S_Bi, relative to Pt: +3.5 and −72 µV/K) | Electronics Cooling (2006), Seebeck coefficient table |
| Dry: plates above coil connections | ΔT = 100 K → EMF 7.55 mV | assumption |
| Wet: cell EMF | 1.97 V (1 × standard 1.970 V: Al³⁺/Al −1.662 V, Bi³⁺/Bi +0.308 V) | CRC electrochemical series; real cells are lower |
| Wet: electrolyte | ρ = 0.21 Ω·m (sea water), ε = 80 | UNESCO 1981 |
| Wet: double layer, charge transfer | 20 µF/cm², 10 Ω·cm² per electrode | Bard & Faulkner (typical); assumed |
| Plates | 1 mm Bi + 1 mm Al | assumption |
| Coil | copper; Wheeler inductance; turn-to-turn capacitance (two-wire formula), V/N per gap ordinary, V/2 bifilar; skin effect at f0 | Wheeler 1928; Tesla US 512,340; Terman 1943 |
| Spectral line | Cu I 324.754 nm (air), 922.9 THz | NIST Handbook of Basic Atomic Spectroscopic Data |

## Dry (Seebeck) variant

5710 evaluations, 271 penalized (171 × invalid: coil too short for Wheeler's formula; 66 × invalid: wire too thick for the coil; 34 × invalid: plate gap not small against the plates; 8 × no-measurement: no resonance; 3 × invalid: wire longer than λ/10 at resonance). 59 non-dominated designs (highest Q first):

| # | turns | coil Ø | wire Ø | pitch | plate | gap | f0 | Q | stored | per volume | ordinary coil: f0 / Q / per volume |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 39 | 150 mm | 2 mm | 1.05×d | 26.53 mm | 1.988 mm | 311.7 kHz | 1.05e+03 | 8.224 fJ | 5.392 pJ/m³ | 621.5 kHz / 1.48e+03 / 5.392 pJ/m³ |
| 2 | 36 | 150 mm | 1.851 mm | 1.26×d | 29.81 mm | 1.66 mm | 436.3 kHz | 1.05e+03 | 10.68 fJ | 6.872 pJ/m³ | 596.1 kHz / 1.22e+03 / 6.872 pJ/m³ |
| 3 | 26 | 150 mm | 2 mm | 1.29×d | 29.47 mm | 1.545 mm | 562.5 kHz | 1.04e+03 | 14.18 fJ | 11.26 pJ/m³ | 675.9 kHz / 1.14e+03 / 11.26 pJ/m³ |
| 4 | 30 | 87.06 mm | 1.889 mm | 1.29×d | 19.62 mm | 883.9 µm | 892.1 kHz | 995 | 10.54 fJ | 22.58 pJ/m³ | 1.069 MHz / 1.09e+03 / 22.58 pJ/m³ |
| 5 | 33 | 83.09 mm | 1.823 mm | 1.24×d | 26.57 mm | 1.306 mm | 778.5 kHz | 942 | 12.72 fJ | 29.15 pJ/m³ | 927.1 kHz / 1.03e+03 / 29.15 pJ/m³ |
| 6 | 33 | 82.26 mm | 1.866 mm | 1.35×d | 33.96 mm | 1.737 mm | 751.7 kHz | 876 | 16.47 fJ | 34.61 pJ/m³ | 851.1 kHz / 931 / 34.61 pJ/m³ |
| 7 | 34 | 63.2 mm | 1.838 mm | 1.29×d | 33.13 mm | 1.634 mm | 966.7 kHz | 850 | 15.11 fJ | 54.21 pJ/m³ | 1.071 MHz / 894 / 54.21 pJ/m³ |
| 8 | 35 | 63.43 mm | 1.823 mm | 1.32×d | 36.11 mm | 1.872 mm | 894.9 kHz | 811 | 16.57 fJ | 56.41 pJ/m³ | 1.007 MHz / 859 / 56.41 pJ/m³ |
| 9 | 33 | 77.85 mm | 1.936 mm | 1.35×d | 34.45 mm | 855.9 µm | 610.2 kHz | 767 | 31.72 fJ | 71.51 pJ/m³ | 651.3 kHz / 792 / 71.51 pJ/m³ |
| 10 | 34 | 52.88 mm | 1.877 mm | 1.23×d | 40.96 mm | 1.391 mm | 853.1 kHz | 729 | 28.29 fJ | 145.1 pJ/m³ | 905.4 kHz / 750 / 145.1 pJ/m³ |
| 11 | 33 | 70.65 mm | 1.861 mm | 1.18×d | 58.74 mm | 1.398 mm | 504.7 kHz | 707 | 49.23 fJ | 154.7 pJ/m³ | 534.2 kHz / 727 / 154.7 pJ/m³ |
| 12 | 33 | 69.68 mm | 1.864 mm | 1.19×d | 59.96 mm | 1.54 mm | 490.4 kHz | 686 | 54.65 fJ | 174 pJ/m³ | 514.9 kHz / 703 / 174 pJ/m³ |
| 13 | 29 | 86.5 mm | 1.862 mm | 1.67×d | 81.74 mm | 1.177 mm | 344.2 kHz | 506 | 124.5 fJ | 212.1 pJ/m³ | 347.9 kHz / 508 / 212.1 pJ/m³ |
| 14 | 10 | 6.482 mm | 846.6 µm | 1.05×d | 23.8 mm | 10 µm | 1.208 MHz | 116 | 1.412 pJ | 911.7 nJ/m³ | 1.208 MHz / 116 / 911.7 nJ/m³ |
| 15 | 10 | 6.589 mm | 758.4 µm | 1.07×d | 29.64 mm | 10 µm | 922.3 kHz | 98.2 | 2.196 pJ | 1.026 µJ/m³ | 922.4 kHz / 98.2 / 1.026 µJ/m³ |
| 16 | 10 | 5.558 mm | 1.023 mm | 1.05×d | 34.79 mm | 10 µm | 1.021 MHz | 95.3 | 3.014 pJ | 1.064 µJ/m³ | 1.021 MHz / 95.3 / 1.064 µJ/m³ |
| 17 | 10 | 5.025 mm | 719.5 µm | 1.05×d | 26.02 mm | 10 µm | 1.292 MHz | 94.4 | 1.708 pJ | 1.085 µJ/m³ | 1.292 MHz / 94.4 / 1.085 µJ/m³ |
| 18 | 10 | 5.068 mm | 905.1 µm | 1.05×d | 30.58 mm | 10 µm | 1.194 MHz | 94.2 | 2.36 pJ | 1.087 µJ/m³ | 1.194 MHz / 94.2 / 1.087 µJ/m³ |
| 19 | 10 | 5.231 mm | 520.6 µm | 1.05×d | 23.64 mm | 10 µm | 1.219 MHz | 88.2 | 1.411 pJ | 1.103 µJ/m³ | 1.219 MHz / 88.2 / 1.103 µJ/m³ |
| 20 | 10 | 6.563 mm | 792.7 µm | 1.06×d | 39.15 mm | 10 µm | 710 kHz | 86.5 | 3.828 pJ | 1.103 µJ/m³ | 710.1 kHz / 86.5 / 1.103 µJ/m³ |
| 21 | 10 | 5.005 mm | 822.7 µm | 1.05×d | 35.4 mm | 10 µm | 1.006 MHz | 83.7 | 3.163 pJ | 1.141 µJ/m³ | 1.006 MHz / 83.7 / 1.141 µJ/m³ |
| 22 | 10 | 5.01 mm | 618.5 µm | 1.06×d | 30.95 mm | 10 µm | 1.031 MHz | 80.9 | 2.417 pJ | 1.149 µJ/m³ | 1.031 MHz / 80.9 / 1.149 µJ/m³ |
| 23 | 10 | 5.02 mm | 593.3 µm | 1.05×d | 31.06 mm | 10 µm | 1.007 MHz | 79.7 | 2.435 pJ | 1.155 µJ/m³ | 1.007 MHz / 79.7 / 1.155 µJ/m³ |
| 24 | 10 | 5.127 mm | 657 µm | 1.05×d | 35.11 mm | 10 µm | 912.3 kHz | 78.3 | 3.094 pJ | 1.156 µJ/m³ | 912.4 kHz / 78.3 / 1.156 µJ/m³ |
| 25 | 10 | 5.084 mm | 624 µm | 1.06×d | 33.39 mm | 10 µm | 947.6 kHz | 78.2 | 2.813 pJ | 1.16 µJ/m³ | 947.7 kHz / 78.2 / 1.16 µJ/m³ |
| 26 | 10 | 5.041 mm | 632 µm | 1.05×d | 34.46 mm | 10 µm | 926.1 kHz | 77.5 | 2.996 pJ | 1.166 µJ/m³ | 926.2 kHz / 77.5 / 1.166 µJ/m³ |
| 27 | 10 | 5.067 mm | 581.4 µm | 1.05×d | 33.33 mm | 10 µm | 924.2 kHz | 76.4 | 2.804 pJ | 1.168 µJ/m³ | 924.3 kHz / 76.4 / 1.168 µJ/m³ |
| 28 | 10 | 5.077 mm | 596.8 µm | 1.05×d | 34.74 mm | 10 µm | 893.7 kHz | 75.6 | 3.046 pJ | 1.172 µJ/m³ | 893.7 kHz / 75.6 / 1.172 µJ/m³ |
| 29 | 10 | 5.037 mm | 584.8 µm | 1.05×d | 35.72 mm | 10 µm | 868.7 kHz | 73.7 | 3.22 pJ | 1.179 µJ/m³ | 868.7 kHz / 73.7 / 1.179 µJ/m³ |
| 30 | 10 | 5 mm | 611.4 µm | 1.05×d | 36.77 mm | 10 µm | 863.3 kHz | 73.7 | 3.412 pJ | 1.18 µJ/m³ | 863.4 kHz / 73.7 / 1.18 µJ/m³ |
| 31 | 10 | 5.63 mm | 557.8 µm | 1.05×d | 46.21 mm | 10 µm | 601.2 kHz | 65.4 | 5.37 pJ | 1.197 µJ/m³ | 601.2 kHz / 65.4 / 1.197 µJ/m³ |
| 32 | 10 | 5.002 mm | 606 µm | 1.05×d | 49.12 mm | 10 µm | 644 kHz | 62.9 | 6.088 pJ | 1.212 µJ/m³ | 644 kHz / 62.9 / 1.212 µJ/m³ |
| 33 | 10 | 5.016 mm | 720.4 µm | 1.05×d | 66.61 mm | 10 µm | 505.6 kHz | 56.9 | 11.2 pJ | 1.226 µJ/m³ | 505.6 kHz / 56.9 / 1.226 µJ/m³ |
| 34 | 10 | 5.1 mm | 796.2 µm | 1.07×d | 73.7 mm | 10 µm | 471.4 kHz | 55.5 | 13.71 pJ | 1.227 µJ/m³ | 471.4 kHz / 55.5 / 1.227 µJ/m³ |
| 35 | 10 | 5.475 mm | 756.1 µm | 1.05×d | 80.44 mm | 10 µm | 394.8 kHz | 54 | 16.33 pJ | 1.23 µJ/m³ | 394.8 kHz / 54 / 1.23 µJ/m³ |
| 36 | 10 | 5.191 mm | 721.7 µm | 1.05×d | 79.75 mm | 10 µm | 410 kHz | 52.3 | 16.05 pJ | 1.234 µJ/m³ | 410 kHz / 52.3 / 1.234 µJ/m³ |
| 37 | 10 | 5 mm | 866.8 µm | 1.05×d | 89.27 mm | 10 µm | 406.9 kHz | 51.3 | 20.11 pJ | 1.235 µJ/m³ | 406.9 kHz / 51.3 / 1.235 µJ/m³ |
| 38 | 10 | 5 mm | 733.6 µm | 1.05×d | 94.77 mm | 10 µm | 358.9 kHz | 47.1 | 22.67 pJ | 1.241 µJ/m³ | 358.9 kHz / 47.1 / 1.241 µJ/m³ |
| 39 | 10 | 5.161 mm | 657 µm | 1.05×d | 92.9 mm | 10 µm | 341.4 kHz | 46.4 | 21.78 pJ | 1.241 µJ/m³ | 341.4 kHz / 46.4 / 1.241 µJ/m³ |
| 40 | 10 | 5 mm | 542.9 µm | 1.08×d | 84.01 mm | 10 µm | 366.1 kHz | 44.4 | 17.81 pJ | 1.242 µJ/m³ | 366.1 kHz / 44.4 / 1.242 µJ/m³ |
| 41 | 10 | 5.06 mm | 513.9 µm | 1.05×d | 84.69 mm | 10 µm | 348.4 kHz | 44.1 | 18.1 pJ | 1.243 µJ/m³ | 348.4 kHz / 44.1 / 1.243 µJ/m³ |
| 42 | 10 | 5.074 mm | 612.6 µm | 1.10×d | 101.4 mm | 10 µm | 314.6 kHz | 41.8 | 25.96 pJ | 1.244 µJ/m³ | 314.6 kHz / 41.8 / 1.244 µJ/m³ |
| 43 | 10 | 5.052 mm | 498.8 µm | 1.06×d | 99.07 mm | 10 µm | 295.7 kHz | 39.9 | 24.77 pJ | 1.247 µJ/m³ | 295.7 kHz / 39.9 / 1.247 µJ/m³ |
| 44 | 10 | 5.109 mm | 510.6 µm | 1.05×d | 103.5 mm | 10 µm | 282 kHz | 39.6 | 27.05 pJ | 1.247 µJ/m³ | 282 kHz / 39.6 / 1.247 µJ/m³ |
| 45 | 11 | 5 mm | 768 µm | 1.07×d | 140.2 mm | 10 µm | 235.1 kHz | 37.9 | 49.61 pJ | 1.247 µJ/m³ | 235.1 kHz / 37.9 / 1.247 µJ/m³ |
| 46 | 10 | 5.136 mm | 734.1 µm | 1.05×d | 144.2 mm | 10 µm | 230.4 kHz | 37.6 | 52.47 pJ | 1.249 µJ/m³ | 230.4 kHz / 37.6 / 1.249 µJ/m³ |
| 47 | 10 | 5.103 mm | 613.4 µm | 1.09×d | 145.9 mm | 10 µm | 217.2 kHz | 34.3 | 53.75 pJ | 1.25 µJ/m³ | 217.2 kHz / 34.3 / 1.25 µJ/m³ |
| 48 | 10 | 5.091 mm | 499.8 µm | 1.10×d | 139.7 mm | 10 µm | 211.2 kHz | 29 | 49.27 pJ | 1.251 µJ/m³ | 211.2 kHz / 29 / 1.251 µJ/m³ |
| 49 | 10 | 5.136 mm | 504.1 µm | 1.05×d | 146.4 mm | 10 µm | 197.6 kHz | 28.4 | 54.07 pJ | 1.251 µJ/m³ | 197.6 kHz / 28.4 / 1.251 µJ/m³ |
| 50 | 11 | 5.048 mm | 493.5 µm | 1.05×d | 150 mm | 10 µm | 182.6 kHz | 26.2 | 56.78 pJ | 1.251 µJ/m³ | 182.6 kHz / 26.2 / 1.251 µJ/m³ |
| 51 | 10 | 5 mm | 458.3 µm | 1.05×d | 144.3 mm | 10 µm | 198.4 kHz | 25.1 | 52.58 pJ | 1.252 µJ/m³ | 198.4 kHz / 25.1 / 1.252 µJ/m³ |
| 52 | 12 | 5 mm | 231.6 µm | 1.05×d | 143.6 mm | 10 µm | 141.2 kHz | 7.83 | 52.05 pJ | 1.253 µJ/m³ | 141.2 kHz / 7.83 / 1.253 µJ/m³ |
| 53 | 15 | 5 mm | 155 µm | 1.43×d | 149.2 mm | 10 µm | 108.7 kHz | 2.93 | 56.2 pJ | 1.253 µJ/m³ | 108.7 kHz / 2.93 / 1.253 µJ/m³ |
| 54 | 14 | 5.048 mm | 149.4 µm | 1.23×d | 149.6 mm | 10 µm | 107.2 kHz | 2.91 | 56.47 pJ | 1.254 µJ/m³ | 107.2 kHz / 2.91 / 1.254 µJ/m³ |
| 55 | 15 | 5.012 mm | 131.1 µm | 1.32×d | 148.5 mm | 10 µm | 98.2 kHz | 2.05 | 55.66 pJ | 1.254 µJ/m³ | 98.2 kHz / 2.05 / 1.254 µJ/m³ |
| 56 | 12 | 5.015 mm | 131.1 µm | 1.64×d | 148.7 mm | 10 µm | 122.3 kHz | 2.04 | 55.82 pJ | 1.254 µJ/m³ | 122.3 kHz / 2.04 / 1.254 µJ/m³ |
| 57 | 15 | 5 mm | 116.2 µm | 1.43×d | 148.7 mm | 10 µm | 92.14 kHz | 1.38 | 55.83 pJ | 1.254 µJ/m³ | 92.14 kHz / 1.38 / 1.254 µJ/m³ |
| 58 | 16 | 5 mm | 111.3 µm | 1.34×d | 148.3 mm | 10 µm | 83.6 kHz | 1.2 | 55.53 pJ | 1.254 µJ/m³ | 83.6 kHz / 1.2 / 1.254 µJ/m³ |
| 59 | 15 | 5 mm | 100 µm | 1.42×d | 148.5 mm | 10 µm | 78.93 kHz | 0.754 | 55.69 pJ | 1.254 µJ/m³ | 78.93 kHz / 0.754 / 1.254 µJ/m³ |

**Highest Q** — 39 turns, L = 225.5 µH, coil self-capacitance 867.8 pF bifilar (ordinary would be 2.282 pF; ratio 380.2 = N²/4 at this N), coil resistance at f0 420.6 mΩ. Element EMF 7.55 mV. At the steady state no current flows (the plate capacitance blocks DC): the EMF only charges the plates, with time constant (R_int + R_coil)·C_plate = 121.7 ps; a switch-on transient shorter than that ends with the plates partly charged. Resonance 311.7 kHz is the n = 2.96075e+09 subharmonic of Cu I, offset -0.111 ppb (any frequency is within 1/(2n) = 0.169 ppb of some subharmonic).

  Switch-on transient over 6.416 ms (80023 steps): energy in from the thermal port 613.7 pJ; dissipated 613.7 pJ; stored at the end 8.216 fJ; integration mismatch -5.965 fJ; worst ledger residual 7.4e-16; model-supplied 0 J. It rings at 312.6 kHz (resonance 311.7 kHz, difference 0.3 %).

**Highest energy per volume** — 15 turns, L = 1.264 µH, coil self-capacitance 3.449 pF bifilar (ordinary would be 61.32 fF; ratio 56.25 = N²/4 at this N), coil resistance at f0 504 mΩ. Element EMF 7.55 mV. At the steady state no current flows (the plate capacitance blocks DC): the EMF only charges the plates, with time constant (R_int + R_coil)·C_plate = 987.3 ns; a switch-on transient shorter than that ends with the plates partly charged. Resonance 78.93 kHz is the n = 1.16929e+10 subharmonic of Cu I, offset 0.000756 ppb (any frequency is within 1/(2n) = 0.0428 ppb of some subharmonic).

  Switch-on transient over 253.4 µs (825 steps): energy in from the thermal port 1.645 nJ; dissipated 1.591 nJ; stored at the end 55.69 pJ; integration mismatch -1.298 pJ; worst ledger residual 1.6e-14; model-supplied 0 J. The ringing dies out within about two cycles (Q = 0.754).

## Wet (galvanic) variant

5710 evaluations, 599 penalized (262 × invalid: wire too thick for the coil; 216 × no-measurement: no resonance; 85 × invalid: coil too short for Wheeler's formula; 46 × invalid: plate gap not small against the plates). 38 non-dominated designs (highest Q first):

| # | turns | coil Ø | wire Ø | pitch | plate | gap | f0 | Q | stored | per volume | ordinary coil: f0 / Q / per volume |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 23 | 141.3 mm | 1.993 mm | 1.36×d | 5 mm | 50 µm | 11.16 kHz | 6.43 | 4.819 µJ | 4.658 mJ/m³ | 11.16 kHz / 6.43 / 4.658 mJ/m³ |
| 2 | 24 | 149.7 mm | 1.994 mm | 1.37×d | 5.469 mm | 50 µm | 9.481 kHz | 6.32 | 5.774 µJ | 4.719 mJ/m³ | 9.482 kHz / 6.32 / 4.719 mJ/m³ |
| 3 | 25 | 150 mm | 1.993 mm | 1.21×d | 5.291 mm | 50.09 µm | 9.209 kHz | 6.31 | 5.408 µJ | 4.771 mJ/m³ | 9.21 kHz / 6.31 / 4.771 mJ/m³ |
| 4 | 10 | 12.27 mm | 2 mm | 1.05×d | 54.01 mm | 50 µm | 12.57 kHz | 5.13 | 559.1 µJ | 57.85 J/m³ | 12.57 kHz / 5.13 / 57.85 J/m³ |
| 5 | 10 | 11.61 mm | 1.999 mm | 1.05×d | 53.8 mm | 50 µm | 13.25 kHz | 5.1 | 554.7 µJ | 59.78 J/m³ | 13.25 kHz / 5.1 / 59.78 J/m³ |
| 6 | 10 | 11.49 mm | 2 mm | 1.05×d | 54.32 mm | 50 µm | 13.25 kHz | 5.09 | 565.5 µJ | 60.56 J/m³ | 13.25 kHz / 5.09 / 60.56 J/m³ |
| 7 | 10 | 11.25 mm | 2 mm | 1.05×d | 53.81 mm | 50 µm | 13.64 kHz | 5.08 | 555 µJ | 60.93 J/m³ | 13.64 kHz / 5.08 / 60.93 J/m³ |
| 8 | 10 | 11.16 mm | 2 mm | 1.05×d | 54 mm | 50 µm | 13.69 kHz | 5.08 | 558.9 µJ | 61.38 J/m³ | 13.69 kHz / 5.08 / 61.38 J/m³ |
| 9 | 10 | 11.13 mm | 1.999 mm | 1.05×d | 56.14 mm | 50 µm | 13.2 kHz | 5.06 | 603.9 µJ | 63.08 J/m³ | 13.2 kHz / 5.06 / 63.08 J/m³ |
| 10 | 10 | 10.73 mm | 2 mm | 1.05×d | 55.88 mm | 50 µm | 13.71 kHz | 5.04 | 598.4 µJ | 64.14 J/m³ | 13.71 kHz / 5.04 / 64.14 J/m³ |
| 11 | 10 | 10.86 mm | 2 mm | 1.05×d | 57.83 mm | 50 µm | 13.1 kHz | 5.03 | 640.7 µJ | 65.09 J/m³ | 13.1 kHz / 5.03 / 65.09 J/m³ |
| 12 | 10 | 10.74 mm | 2 mm | 1.05×d | 58.52 mm | 50 µm | 13.07 kHz | 5.02 | 656 µJ | 65.92 J/m³ | 13.07 kHz / 5.02 / 65.92 J/m³ |
| 13 | 10 | 10.73 mm | 2 mm | 1.05×d | 60.86 mm | 50 µm | 12.57 kHz | 4.99 | 709.4 µJ | 67.44 J/m³ | 12.57 kHz / 4.99 / 67.44 J/m³ |
| 14 | 10 | 10.34 mm | 2 mm | 1.05×d | 59.31 mm | 50 µm | 13.36 kHz | 4.99 | 673.8 µJ | 67.65 J/m³ | 13.36 kHz / 4.99 / 67.65 J/m³ |
| 15 | 10 | 10.42 mm | 2 mm | 1.05×d | 60.47 mm | 50 µm | 13 kHz | 4.98 | 700.3 µJ | 68.1 J/m³ | 13 kHz / 4.98 / 68.1 J/m³ |
| 16 | 10 | 10.51 mm | 2 mm | 1.11×d | 66.54 mm | 50 µm | 11.97 kHz | 4.86 | 847.3 µJ | 70.36 J/m³ | 11.97 kHz / 4.86 / 70.36 J/m³ |
| 17 | 10 | 9.86 mm | 1.847 mm | 1.15×d | 68.88 mm | 56 µm | 12.07 kHz | 4.4 | 905.5 µJ | 73.96 J/m³ | 12.07 kHz / 4.4 / 73.96 J/m³ |
| 18 | 11 | 9.51 mm | 1.862 mm | 1.05×d | 79.87 mm | 60.87 µm | 9.811 kHz | 4.13 | 1.215 mJ | 78.26 J/m³ | 9.811 kHz / 4.13 / 78.26 J/m³ |
| 19 | 10 | 10.22 mm | 1.801 mm | 1.05×d | 94.56 mm | 50 µm | 8.066 kHz | 3.99 | 1.701 mJ | 82.28 J/m³ | 8.066 kHz / 3.99 / 82.28 J/m³ |
| 20 | 10 | 8.27 mm | 1.622 mm | 1.05×d | 87.67 mm | 50 µm | 10.1 kHz | 3.97 | 1.462 mJ | 85.01 J/m³ | 10.1 kHz / 3.97 / 85.01 J/m³ |
| 21 | 10 | 8.262 mm | 1.612 mm | 1.05×d | 96.21 mm | 50.01 µm | 9.183 kHz | 3.76 | 1.756 mJ | 86.09 J/m³ | 9.183 kHz / 3.76 / 86.09 J/m³ |
| 22 | 10 | 8.292 mm | 1.622 mm | 1.05×d | 110.3 mm | 50.27 µm | 7.978 kHz | 3.47 | 2.299 mJ | 87.17 J/m³ | 7.978 kHz / 3.47 / 87.17 J/m³ |
| 23 | 10 | 8.443 mm | 1.664 mm | 1.05×d | 117.5 mm | 52.9 µm | 7.442 kHz | 3.35 | 2.608 mJ | 87.22 J/m³ | 7.442 kHz / 3.35 / 87.22 J/m³ |
| 24 | 10 | 8.406 mm | 1.666 mm | 1.05×d | 120 mm | 50.01 µm | 7.316 kHz | 3.34 | 2.718 mJ | 87.54 J/m³ | 7.316 kHz / 3.34 / 87.54 J/m³ |
| 25 | 10 | 8.349 mm | 1.645 mm | 1.05×d | 121 mm | 50.02 µm | 7.257 kHz | 3.29 | 2.763 mJ | 87.66 J/m³ | 7.257 kHz / 3.29 / 87.66 J/m³ |
| 26 | 10 | 8.216 mm | 1.623 mm | 1.05×d | 119.6 mm | 50 µm | 7.414 kHz | 3.29 | 2.695 mJ | 87.73 J/m³ | 7.414 kHz / 3.29 / 87.73 J/m³ |
| 27 | 10 | 8.206 mm | 1.623 mm | 1.05×d | 121.3 mm | 50 µm | 7.307 kHz | 3.26 | 2.775 mJ | 87.81 J/m³ | 7.307 kHz / 3.26 / 87.81 J/m³ |
| 28 | 10 | 8.246 mm | 1.622 mm | 1.05×d | 123.3 mm | 50 µm | 7.155 kHz | 3.21 | 2.864 mJ | 87.85 J/m³ | 7.155 kHz / 3.21 / 87.85 J/m³ |
| 29 | 10 | 8.253 mm | 1.622 mm | 1.05×d | 125.3 mm | 50 µm | 7.032 kHz | 3.17 | 2.956 mJ | 87.91 J/m³ | 7.032 kHz / 3.17 / 87.91 J/m³ |
| 30 | 10 | 8.202 mm | 1.623 mm | 1.05×d | 126.4 mm | 50 µm | 7.012 kHz | 3.16 | 3.006 mJ | 87.99 J/m³ | 7.012 kHz / 3.16 / 87.99 J/m³ |
| 31 | 10 | 8.218 mm | 1.623 mm | 1.05×d | 131 mm | 50 µm | 6.75 kHz | 3.07 | 3.224 mJ | 88.09 J/m³ | 6.75 kHz / 3.07 / 88.09 J/m³ |
| 32 | 10 | 8.239 mm | 1.627 mm | 1.05×d | 133.5 mm | 50 µm | 6.614 kHz | 3.03 | 3.345 mJ | 88.13 J/m³ | 6.614 kHz / 3.03 / 88.13 J/m³ |
| 33 | 10 | 5 mm | 872.1 µm | 1.05×d | 55.54 mm | 50 µm | 19.5 kHz | 2.84 | 585.3 µJ | 88.75 J/m³ | 19.5 kHz / 2.84 / 88.75 J/m³ |
| 34 | 10 | 5 mm | 764.2 µm | 1.05×d | 56.61 mm | 50.2 µm | 17.99 kHz | 2.51 | 605.4 µJ | 89.03 J/m³ | 17.99 kHz / 2.51 / 89.03 J/m³ |
| 35 | 10 | 5.021 mm | 784.7 µm | 1.05×d | 61.16 mm | 50 µm | 16.71 kHz | 2.49 | 705.5 µJ | 89.21 J/m³ | 16.71 kHz / 2.49 / 89.21 J/m³ |
| 36 | 10 | 5.019 mm | 770.1 µm | 1.05×d | 61.39 mm | 50 µm | 16.5 kHz | 2.43 | 710.2 µJ | 89.22 J/m³ | 16.5 kHz / 2.43 / 89.22 J/m³ |
| 37 | 10 | 5.028 mm | 775.5 µm | 1.05×d | 63.04 mm | 50 µm | 16.07 kHz | 2.41 | 748.3 µJ | 89.26 J/m³ | 16.07 kHz / 2.41 / 89.26 J/m³ |
| 38 | 10 | 5.035 mm | 782.4 µm | 1.05×d | 64.89 mm | 50 µm | 15.63 kHz | 2.4 | 792.1 µJ | 89.29 J/m³ | 15.63 kHz / 2.4 / 89.29 J/m³ |

**Highest Q** — 23 turns, L = 82.66 µH, coil self-capacitance 169.1 pF bifilar (ordinary would be 1.279 pF; ratio 132.2 = N²/4 at this N), coil resistance at f0 54.96 mΩ. Element EMF 1.97 V. At the steady state the cell drives 24.48 mA through the coil; source power 48.22 mW, all of it dissipated (48.22 mW). Aluminium plate capacity 724.2 C (0.201 Ah): at the steady-state current it lasts 8h13m2s. Resonance 11.16 kHz is the n = 8.27112e+10 subharmonic of Cu I, offset 0.000863 ppb (any frequency is within 1/(2n) = 0.00605 ppb of some subharmonic).

  Switch-on transient over 1.792 ms (832 steps): energy in from the chemical port 118.2 µJ; dissipated 91.19 µJ; stored at the end 4.819 µJ; integration mismatch -28.78 nJ; worst ledger residual 8.2e-16; model-supplied 0 J. It rings at 11.05 kHz (resonance 11.16 kHz, difference -0.98 %).

**Highest energy per volume** — 10 turns, L = 238.1 nH, coil self-capacitance 14.56 pF bifilar (ordinary would be 582.4 fF; ratio 25 = N²/4 at this N), coil resistance at f0 5.527 mΩ. Element EMF 1.97 V. At the steady state the cell drives 4.079 A through the coil; source power 8.035 W, all of it dissipated (8.035 W). Aluminium plate capacity 122 kC (33.9 Ah): at the steady-state current it lasts 8h18m24s. Resonance 15.63 kHz is the n = 5.90365e+10 subharmonic of Cu I, offset 0.00301 ppb (any frequency is within 1/(2n) = 0.00847 ppb of some subharmonic).

  Switch-on transient over 1.279 ms (832 steps): energy in from the chemical port 13.1 mJ; dissipated 11.07 mJ; stored at the end 792.1 µJ; integration mismatch -4.686 µJ; worst ledger residual 3e-16; model-supplied 0 J. The ringing dies out within about two cycles (Q = 2.4).

## Thermopile (Seebeck, n pairs) variant

5710 evaluations, 423 penalized (168 × invalid: thermopile legs must be at least 5 diameters long (1-D conduction); 147 × invalid: coil too short for Wheeler's formula; 59 × invalid: wire longer than λ/10 at resonance; 44 × invalid: plate gap not small against the plates; 11 × invalid: wire too thick for the coil). 19 non-dominated designs (highest Q first):

| # | turns | coil Ø | wire Ø | pitch | plate | gap | pairs | leg Ø × length | f0 | Q | stored | per volume | ordinary coil: f0 / Q / per volume |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 11 | 20.41 mm | 1.703 mm | 1.31×d | 19.64 mm | 1.941 mm | 3 | 562.3 µm × 3.34 mm | 21.15 MHz | 1.25e+03 | 997.6 aJ | 85.84 pJ/m³ | em: outside model validity: wire length 0.705 m exceeds λ/10 = 0.512 m at the resonance 5.853e+07 Hz (not quasi-static) |
| 2 | 11 | 19.53 mm | 1.408 mm | 1.23×d | 20.74 mm | 1.908 mm | 3 | 516.3 µm × 3.276 mm | 19.19 MHz | 1.14e+03 | 1.039 fJ | 119.2 pJ/m³ | em: outside model validity: wire length 0.675 m exceeds λ/10 = 0.555 m at the resonance 5.402e+07 Hz (not quasi-static) |
| 3 | 10 | 20.3 mm | 1.439 mm | 1.31×d | 36.32 mm | 1.974 mm | 3 | 382.5 µm × 3.262 mm | 20.6 MHz | 1.07e+03 | 3.164 fJ | 247.4 pJ/m³ | 36.03 MHz / 1.05e+03 / 247.4 pJ/m³ |
| 4 | 10 | 18.41 mm | 1.321 mm | 1.32×d | 33.24 mm | 2 mm | 3 | 509.7 µm × 3.338 mm | 22.67 MHz | 1.05e+03 | 2.557 fJ | 252 pJ/m³ | 41.93 MHz / 1.2e+03 / 252 pJ/m³ |
| 5 | 11 | 30.69 mm | 1.271 mm | 1.31×d | 32.27 mm | 1.927 mm | 5 | 511.4 µm × 3.294 mm | 11.13 MHz | 1.05e+03 | 6.953 fJ | 350.4 pJ/m³ | 25.27 MHz / 1.31e+03 / 350.4 pJ/m³ |
| 6 | 11 | 19.59 mm | 1.16 mm | 1.31×d | 36.98 mm | 1.926 mm | 7 | 515.7 µm × 3.461 mm | 17.86 MHz | 937 | 17.62 fJ | 1.543 nJ/m³ | 32.16 MHz / 899 / 1.543 nJ/m³ |
| 7 | 16 | 7.092 mm | 465.8 µm | 1.07×d | 101.3 mm | 10 µm | 4945 | 51.42 µm × 50 mm | 30.98 MHz | 334 | 632.9 µJ | 28.73 J/m³ | em: outside model validity: wire length 0.356 m exceeds λ/10 = 0.121 m at the resonance 2.478e+08 Hz (not quasi-static) |
| 8 | 15 | 6.866 mm | 332.1 µm | 1.07×d | 104.4 mm | 10 µm | 4983 | 53.26 µm × 41.82 mm | 32.31 MHz | 260 | 677.3 µJ | 29.39 J/m³ | em: outside model validity: wire length 0.324 m exceeds λ/10 = 0.124 m at the resonance 2.424e+08 Hz (not quasi-static) |
| 9 | 15 | 6.99 mm | 296.8 µm | 1.07×d | 150 mm | 10 µm | 5000 | 53.21 µm × 44.9 mm | 31.67 MHz | 254 | 1.409 mJ | 30.34 J/m³ | em: outside model validity: wire length 0.329 m exceeds λ/10 = 0.126 m at the resonance 2.375e+08 Hz (not quasi-static) |
| 10 | 20 | 5 mm | 632.7 µm | 1.06×d | 130.6 mm | 10 µm | 5000 | 50 µm × 18.03 mm | 33.05 MHz | 253 | 1.075 mJ | 30.75 J/m³ | em: outside model validity: wire length 0.314 m exceeds λ/10 = 0.0907 m at the resonance 3.305e+08 Hz (not quasi-static) |
| 11 | 13 | 5 mm | 287.9 µm | 1.05×d | 147.5 mm | 10 µm | 5000 | 50.01 µm × 22.46 mm | 38.82 MHz | 228 | 1.372 mJ | 31.01 J/m³ | em: outside model validity: wire length 0.204 m exceeds λ/10 = 0.119 m at the resonance 2.524e+08 Hz (not quasi-static) |
| 12 | 15 | 5 mm | 378.5 µm | 1.07×d | 146.5 mm | 10 µm | 5000 | 50 µm × 17.45 mm | 39.42 MHz | 210 | 1.355 mJ | 31.03 J/m³ | em: outside model validity: wire length 0.236 m exceeds λ/10 = 0.101 m at the resonance 2.957e+08 Hz (not quasi-static) |
| 13 | 10 | 5.038 mm | 372.8 µm | 1.08×d | 146.8 mm | 10 µm | 5000 | 50 µm × 15.86 mm | 65.47 MHz | 209 | 1.36 mJ | 31.09 J/m³ | em: outside model validity: wire length 0.158 m exceeds λ/10 = 0.0916 m at the resonance 3.274e+08 Hz (not quasi-static) |
| 14 | 10 | 5.43 mm | 311.9 µm | 1.05×d | 150 mm | 10 µm | 5000 | 50 µm × 15.26 mm | 54.89 MHz | 188 | 1.419 mJ | 31.11 J/m³ | em: outside model validity: wire length 0.171 m exceeds λ/10 = 0.109 m at the resonance 2.745e+08 Hz (not quasi-static) |
| 15 | 10 | 5 mm | 230.3 µm | 1.05×d | 147 mm | 10 µm | 5000 | 50 µm × 15.85 mm | 50.65 MHz | 181 | 1.363 mJ | 31.12 J/m³ | em: outside model validity: wire length 0.157 m exceeds λ/10 = 0.118 m at the resonance 2.533e+08 Hz (not quasi-static) |
| 16 | 10 | 5 mm | 278.7 µm | 1.05×d | 150 mm | 10 µm | 5000 | 50 µm × 13.33 mm | 59.72 MHz | 167 | 1.419 mJ | 31.15 J/m³ | em: outside model validity: wire length 0.157 m exceeds λ/10 = 0.1 m at the resonance 2.986e+08 Hz (not quasi-static) |
| 17 | 24 | 5 mm | 190.5 µm | 1.05×d | 147 mm | 10 µm | 5000 | 55.34 µm × 8.966 mm | 16.22 MHz | 86.6 | 1.364 mJ | 31.16 J/m³ | em: outside model validity: wire length 0.377 m exceeds λ/10 = 0.154 m at the resonance 1.947e+08 Hz (not quasi-static) |
| 18 | 18 | 5.154 mm | 114.9 µm | 1.10×d | 150 mm | 10 µm | 5000 | 51.63 µm × 10 mm | 23.05 MHz | 81.8 | 1.419 mJ | 31.21 J/m³ | em: outside model validity: wire length 0.291 m exceeds λ/10 = 0.145 m at the resonance 2.075e+08 Hz (not quasi-static) |
| 19 | 19 | 5.13 mm | 122.4 µm | 1.10×d | 148.7 mm | 10 µm | 5000 | 53.94 µm × 5.851 mm | 27.75 MHz | 43.3 | 1.396 mJ | 31.25 J/m³ | em: outside model validity: wire length 0.306 m exceeds λ/10 = 0.114 m at the resonance 2.636e+08 Hz (not quasi-static) |

**Highest Q** — 11 turns, L = 1.469 µH, coil self-capacitance 34.66 pF bifilar (ordinary would be 1.146 pF; ratio 30.25 = N²/4 at this N), coil resistance at f0 156.2 mΩ. Element EMF 22.65 mV. At the steady state no current flows (the plate capacitance blocks DC): the EMF only charges the plates, with time constant (R_int + R_coil)·C_plate = 819 fs; a switch-on transient shorter than that ends with the plates partly charged. Thermopile: 3 pairs, legs 562.3 µm × 3.34 mm, internal resistance 54.44 mΩ; holding ΔT = 100 K costs 5.463 W of heat conducted through the legs (not part of the electrical ledger). Resonance 21.15 MHz is the n = 4.36304e+07 subharmonic of Cu I, offset 0.949 ppb (any frequency is within 1/(2n) = 11.5 ppb of some subharmonic).

  Switch-on transient over 94.55 µs (80017 steps): energy in from the thermal port 160.8 fJ; dissipated 159.8 fJ; stored at the end 996.5 aJ; integration mismatch 8.138 aJ; worst ledger residual 2.5e-16; model-supplied 0 J. It rings at 21.21 MHz (resonance 21.15 MHz, difference 0.29 %).

**Highest energy per volume** — 19 turns, L = 1.921 µH, coil self-capacitance 17.12 pF bifilar (ordinary would be 189.7 fF; ratio 90.25 = N²/4 at this N), coil resistance at f0 1.08 Ω. Element EMF 37.75 V. At the steady state no current flows (the plate capacitance blocks DC): the EMF only charges the plates, with time constant (R_int + R_coil)·C_plate = 33.02 ms; a switch-on transient shorter than that ends with the plates partly charged. Thermopile: 5000 pairs, legs 53.94 µm × 5.851 mm, internal resistance 16.86 kΩ; holding ΔT = 100 K costs 47.83 W of heat conducted through the legs (not part of the electrical ledger). Resonance 27.75 MHz is the n = 3.32603e+07 subharmonic of Cu I, offset -3.21 ppb (any frequency is within 1/(2n) = 15 ppb of some subharmonic).

  Switch-on transient over 4.681 µs (5207 steps): energy in from the thermal port 395.7 nJ; dissipated 395.7 nJ; stored at the end 32.86 pJ; integration mismatch -29.55 fJ; worst ledger residual 6.8e-16; model-supplied 0 J. It rings at 27.69 MHz (resonance 27.75 MHz, difference -0.21 %).

## Reading the results

- **dry (Seebeck), bifilar against ordinary.** On 13 of 59 front designs the bifilar winding lowers the resonance by more than 0.1 %; on 0 it raises Q by more than 0.1 %. The stored energy at the steady state is the same for both windings: at DC no voltage stands across the coil's turn-to-turn capacitance, so the winding style changes the resonance, not the stored energy.
- **dry (Seebeck):** the highest energy per volume sits at the edge of the allowed range for coil_diameter_m, wire_diameter_m, insulation_eps, gap_m, dielectric_eps; a wider range would move it further, so that optimum is set by the chosen limits.
- **dry (Seebeck), physical limit.** The only store is the field in the gap, ½·ε·(EMF/gap)², diluted by the plates: for the best design 1.255 µJ/m³, against 1.254 µJ/m³ simulated (coil volume adds to the denominator). With a millivolt-level Seebeck EMF the energy per volume is necessarily microscopic; a thermopile of many junctions in series, or a larger ΔT, raises it as EMF².
- **wet (galvanic), bifilar against ordinary.** On 0 of 38 front designs the bifilar winding lowers the resonance by more than 0.1 %; on 0 it raises Q by more than 0.1 %. The stored energy at the steady state is the same for both windings: at DC no voltage stands across the coil's turn-to-turn capacitance, so the winding style changes the resonance, not the stored energy.
- **wet (galvanic):** the highest energy per volume sits at the edge of the allowed range for turns, pitch_ratio, gap_m; a wider range would move it further, so that optimum is set by the chosen limits.
- **Wet, what limits it.** The stored energy is mostly ½·L·I² from the steady cell current; the current is set by the assumed charge-transfer resistance (10 Ω·cm² per electrode) and reaches 4.079 A (96.9 mA/cm²) in the best design. A real aluminium electrode passivates well below such current densities, so these numbers are upper bounds until the cell is measured. The aluminium plate's charge lasts the same few hours for every design because both capacity and current scale with plate area.
- **thermopile (Seebeck, n pairs), bifilar against ordinary.** On 4 of 19 front designs the bifilar winding lowers the resonance by more than 0.1 %; on 2 it raises Q by more than 0.1 %. The stored energy at the steady state is the same for both windings: at DC no voltage stands across the coil's turn-to-turn capacitance, so the winding style changes the resonance, not the stored energy.
- **thermopile (Seebeck, n pairs):** the highest energy per volume sits at the edge of the allowed range for gap_m, dielectric_eps, pairs; a wider range would move it further, so that optimum is set by the chosen limits.
- **thermopile (Seebeck, n pairs), physical limit.** The only store is the field in the gap, ½·ε·(EMF/gap)², diluted by the plates: for the best design 31.39 J/m³, against 31.25 J/m³ simulated (coil volume adds to the denominator). With 5000 pairs the EMF is 37.75 V, 5000 times one junction's, and the stored energy rises with its square; the legs' heat leak (47.83 W) is the price of holding ΔT.
- Energy accounting: 0 evaluated designs failed the operating-point ledger check. Over the 6 switch-on transients: worst per-point ledger residual 1.6e-14 (relative to the run's peak power), worst integration mismatch 0.00079 of the energy taken in, model-supplied energy 0 J. So all stored and dissipated energy is accounted for by the element's declared source — heat for the dry variant, the aluminium plate's chemistry for the wet one — and nothing else.
- The bifilar winding raises the coil's self-capacitance by about N²/4 compared with an ordinary winding of the same geometry (Tesla's 250,000 at N = 1,000), which lowers the self-resonance; the inductance is unchanged.
- The idealized turn-to-turn model is a first approximation (it underestimates an ordinary coil's measured self-capacitance; Medhurst 1947). Field-solved values arrive with `em/extract` (M9i).
- The neutrino channel is not yet computed (the environment model, M9k); it is listed as not modeled.
