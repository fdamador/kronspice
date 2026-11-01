# Ball-bearing (Huber) motor — P element

A current through a rotating ball bearing drives the shaft. This note
records what has been measured, how kronspice models it, and what it means for
the homopolar and SPG analyses. **Rule:** coefficients come from measured
tables only. The element has no default κ or ω0, and a netlist without them
is rejected.

## Measured facts (open-access record)

| Finding | Source | Status |
|---|---|---|
| Steady currents 43.5–70.15 A; self-rotation at high current with small efficiency; at low current a driving motor is needed | Moyssides & Hatzikonstantinou, *IEEE Trans. Magn.* 26, 1274 (1990), abstract | measured |
| Back-EMF that depends on shaft speed | same | measured |
| No generator action | same | measured |
| Self-rotation 70–140 A; aided 18–140 A (two discs) | Moyssides 1992, cited by Lauterbach, Soong & Abbott, *Proc. SPIE* 4236 (2001) | measured (secondary) |
| Torque ∝ I² | Choo, Soong & Abbott, *Proc. SPIE* 5649 (2005), per McDonald (2020) | measured |
| Mechanical output ≈ 10 % of electrical input | Marinov, *Electronics & Wireless World*, Apr 1989 (calorimeter) | measured, data unpublished |
| Heat ≈ 8 % higher with rotation | same | claim (H14) |
| Needs a hand spin to start; sense set by the spin | Harrison (electricstuff); Huber 1959 | measured |
| Thermal effects start it; Lorentz force cannot start it from rest | McDonald (2020) fn. 6 and his 2016 experiments | measured |
| 1000 r/min (Milroy 1959); ≥ 1500 r/min (Ogden 1989) | Milroy *J. Appl. Mech.* 34 (1967); Ogden 1989 | measured / estimate |

No open-access source gives a torque–current–speed table. The papers that
do are listed in `Needed` in `examples/bearing_motor/main.go`
(Gruenberg 1978, Moyssides 1990/1997, Watson 1993/1999, Choo 2005,
Shen 1999, Polivanov 1973, Huber 1959). Once one is attached, its table goes
into `Datasets`, and `FitTorque` or `FitBackEMF` gives κ and ω0 with an rms
residual. In EM mode the two fits must agree, since E·I = ω·T.

## Element

`Pname a b m+ m- R KAPPA=κ W0=ω0 [THERMAL=1 ETAMAX=η]` (`device/bearing.go`)

- **Torque:** T = κ I² ψ(ω), with ψ = tanh(ω/ω0). This form has zero torque
  at rest, torque in the sense of rotation, and I² scaling. The rest state is
  unstable above the threshold current I_c = √(ω0 / (κ R_d)), where R_d is
  the shaft's mechanical damping. Below I_c a spin decays, matching the
  measured threshold.
- **EM mode:** V = R I + κ ψ ω I (back-EMF), so V I = R I² + ω T exactly.
- **Thermal mode:** V = R I. Shaft work comes out of the heat, capped by the
  measured η: ω T < η R I².
- **Ledger:** the heat V I − ω T leaves through the thermal port. Every
  joule of shaft work is drawn from the electrical circuit.

Tests:
- `device/bearing_test.go`: Jacobian, and an allocation-free `Load`.
- `analysis/bearing_test.go`: no start from rest; run-up above I_c in either
  sense; decay below I_c; closed ledgers in both modes.
- `circuit/bearing_test.go`: netlist round trip; missing coefficients are
  rejected.

## Application: the SPG (H13)

All inputs are measured: the SPG test (UPSM pp. 159–160) and the bearing
efficiency above.

- The claim needs k·I = 41.5 N·m on the shaft, but the drive motor supplied
  16.7 N·m. A bearing torque aiding rotation would have to supply the
  24.8 N·m gap (7.0 kW). If the full 3438 A crosses the bearings, that means
  κψ ≥ 2.1 × 10⁻⁶ N·m/A².
- **EM mode** (back-EMF measured, so this is the measured mode): the
  bearings then drop 2.04 V, leaving 1.37 V and 4.7 kW for the rest of the
  loop. This is exactly the SPG ledger bound. A bearing motor moves shaft
  work around; it adds none. The engine run sweeps the unknown contact
  resistance, and the load never gets more.
- **Thermal mode** at η = 0.10: the bearings would need 70 kW of
  electrical input, which is 20.4 V across them, more than the 3.41 V the
  machine makes. Excluded.

**Bench tests that decide it:**
1. Insulate the bearings (ceramic balls) and repeat the loaded run.
2. Measure the voltage across each bearing under load.
3. Put a torque transducer between the drive motor and the rotor.

The same element carries over to a DePalma N-machine model, whose
brush/bearing current path it will be part of.

## Contact drop (VF)

A ball bearing carrying DC current does not behave like a resistor. It
drops an almost constant voltage: 0.73 V on average at 4 A, and
0.5–1.5 V across the reported tests. The cause is repeated breakdown of the
oxide and lubricant film. Above an apparent current density of
0.1–0.3 A/mm² of Hertzian contact area, the bearing is damaged. For the
51208 bearing that area is A_Hz = 1.0 mm² at 400 N and 2.5 mm² at 1600 N.
Source: Safdarzadeh, Farahi, Binder et al., *Lubricants* 12(5), 148 (2024).

The P element carries this drop as `VF`, so the contact heat is
H = VF·|I| + R·I². In thermal mode the measured efficiency caps the shaft
work at η·H.

## Tesla homopolar with the current through its bearings (`examples/homopolar_bearing`)

This is a circuit simulation only. Each disk is an N element; each shaft
bearing is a P element carrying the measured VF and driving the shaft;
bearing friction M = ½·f·F·d (f = 0.0015) is a mechanical resistance; the
load is the external R_L, or a capacitor C charged through R_L.

**The trade-off.** The bearing bore sets the axle radius r_i, and with it
the disk's magnetic area, EMF = ½Bω(r_o² − r_i²). Bigger bearings, or more
heavily loaded ones, carry more current before damage (I_max = J·A_Hz,
with A_Hz ∝ (F·d)^(2/3) by Hertz), but their friction grows as F·d.

**Results.** These are for 15 cm copper disks in 0.8 T, with the rim speed
held to the yield of annealed copper with a 1.5 safety factor.

- At the damage limit, a bearing bank passes only 0.2–1.3 A, so a
  through-bearing unit is a watts-to-tens-of-watts machine. The reference
  unit runs at 7.6× the limit.
- Evolved for at least 20 W with the least drive torque, the GA chose six
  108 mm bearings per shaft, exactly at the current limit. That gives
  0.0101 N·m/W at 13 % efficiency; friction dominates.
- With the bearings insulated and the current on a shaft collector (its
  resistance declared, swept 10 µΩ–1 mΩ with no effect), the evolved unit
  gives 689 W at 99.8 %, i.e. 0.00135 N·m/W = 1/(η·ω). That is 7.4× less
  torque per watt.
- **Huber torque, EM mode:** it lowers the drive torque only by lowering
  the output with its back-EMF, so it never helps the torque per watt.
- **Huber torque, thermal mode:** it recovers at most η·H, about 10 % of
  the bearing heat.
- In every case the drive torque per load watt is at least 1/ω.
- **Capacitor load:** the bearing drops cost 2·V_b·Q however slowly the
  field is ramped. Electrical conversion rises from 43 % with a step to
  73 % with a slow ramp, against a 75.8 % ceiling. Bearing friction keeps
  running throughout, so the fastest charge that keeps resistive loss low
  is best.

**Bench data still needed:**
- the Huber coefficient κψ (examples/bearing_motor);
- the belt's rolling-contact resistance;
- the real collector's contact drop. A carbon brush also drops about 1 V,
  which would narrow the insulated design's lead; a liquid-metal collector
  would not.
