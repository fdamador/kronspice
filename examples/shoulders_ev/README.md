# Ken Shoulders' EVs: charge clusters tested with standard physics (H27)

K. R. Shoulders, *EV: A Tale of Discovery* (Jupiter Technologies, Austin,
1987). Page numbers below are the book's chapter-page numbers.

Shoulders reports small beads of electric charge thrown off by fast
discharges from sharp cathodes. He calls them **EVs** (*electrum validum*,
"strong electron"). What he reports about them:
- **Size and charge.** About 1 µm across, often in chains and rings. One
  bead carries "around 10¹¹ charges" (8-49).
- **Content.** Electrons, at a density near that of a solid
  (6.6×10²³ per cm³), with fewer than 1 ion per 100,000 electrons (4-63).
- **Motion and life.** They move at about c/10, live about 10⁻¹¹ s, and
  bore pits, leave ring-shaped marks and give X-rays where they hit.
- **What he thinks holds them together.** The zero-point fluctuations of
  empty space (Preface; 9-2).
- **Source.** They can be made from a cathode of about 1 pF storing less
  than 10⁻⁶ J (2-7).

This example takes those numbers at face value and checks them against
Coulomb's and Gauss's laws, the known containers, and charge and energy
budgets. It then runs his source and an EV strike as kronspice circuits
with energy ledgers.

```
go run ./examples/shoulders_ev -viz examples/shoulders_ev/viz   # under a second
python3 examples/shoulders_ev/viz/render.py
go test ./examples/shoulders_ev
```

## 1. A bead of 10¹¹ electrons, 1 µm across

| Quantity | Value |
|---|---|
| Charge | 1.6×10⁻⁸ C |
| Density | 1.9×10²³ electrons/cm³ (his 3.5×10¹⁴ electrons in 10 µm gives his 6.6×10²³) |
| Energy needed to assemble it | **2.77 J** |
| Field at its surface | 5.8×10¹⁴ V/m |
| Outward pressure | 1.5×10¹⁸ Pa (10¹³ atmospheres) |
| Kinetic energy of all its electrons at c/10 | 4.1×10⁻⁵ J (the field holds 67,000 times more) |

**Left bare, it flies apart.** The outer electrons are pushed off at
nearly the speed of light:
- the radius doubles in 1.7×10⁻¹⁵ s and grows tenfold in 1.5×10⁻¹⁴ s;
- each outer electron leaves with about 260 MeV.

Shoulders' lifetime (1 mm of track at c/10) is 3.3×10⁻¹¹ s, 20,000 times
longer than the doubling time. Something would have to hold it.

## 2. What could hold it?

| Candidate | Pressure | Against the 1.5×10¹⁸ Pa push |
|---|---|---|
| Magnetic attraction of electrons moving together at c/10 | 1% of the repulsion | the push still wins |
| A magnetic trap at this density (Brillouin limit) | needs 2×10⁵ T | the strongest field made indoors is about 1,200 T |
| Casimir pull between ideal plates 1 µm apart | 1.3×10⁻³ Pa | 10²¹ times too weak |
| Zero-point stress on a 1 µm conducting sphere | 1.9×10⁻³ Pa, **outward** | pushes the wrong way (Boyer 1968) |
| Positive ions | would cancel the charge | his own measurement: under 1 ion per 100,000 electrons |

The zero-point field is the container Shoulders proposes. Its known,
measured effect — the Casimir force — is 21 orders of magnitude too weak.
For a sphere it even pushes outward. Picture: left panel of
`viz/shoulders_ev.png`.

## 3. Gauss's law and his kilovolt measurements

- **The puzzle he states (9-2).** An electron leaving a bare 1 µm bead of
  10¹¹ electrons gains 288 MeV. His probes saw "only kilovolts at best".
- **His answer.** A container hides the core's field.
- **Gauss's law rules out hiding it.** Outside any closed surface, the
  potential of a net charge q is k·q/r whatever is inside. For escaping
  electrons to gain only V volts:

| Escape energy seen | The "container" must be at least |
|---|---|
| 1 kV | 14.4 cm in radius |
| 10 kV | 1.4 cm in radius |

Those sizes are larger than the millimetre gaps of his apparatus. The
alternative is a net charge of at most about 3×10⁵ electrons for a 1 µm
bead seen at 1 kV. That is 300,000 times less than 10¹¹, unless positive
charge cancels it, which his ion measurement excludes.

## 4. Charge and energy budgets

| Source (2-7) | Electrons it holds (C·V) | Energy |
|---|---|---|
| smallest: 1 pF, < 1 µJ (at most 1.4 kV) | 8.8×10⁹ | 1.0 µJ |
| earlier: 40 pF at 10 kV | 2.5×10¹² | 2.0 mJ |

| What he reports | Electrons | ÷ small source | ÷ larger source |
|---|---|---|---|
| residual charge of a 3 µm EV on an electrode (4-63) | 2×10¹⁰ | 2.3 | 0.008 |
| one EV bead (8-49) | 10¹¹ | 11 | 0.04 |
| emitted by a 10 µm EV over 1 mm of track (4-9) | 3.5×10¹⁴ | 40,000 | **140** |

- **Charge.** Charge is conserved exactly. A 10¹¹-electron bead cannot
  come from the 1 pF source, and the 3.5×10¹⁴ electrons he infers from
  camera brightness cannot come from either source. That count is
  probably the inference to check first: it rests on a phosphor
  calibration and on isotropic emission.
- **Energy.** A 1 µm bead of 10¹¹ electrons needs 2.77 J to assemble:
  3 million times the small source's energy and 1,400 times the larger
  one's.

Picture: middle panel.

## 5. kronspice: the source and a strike, with ledgers

The bench (declared values): a 1 pF anode — the witness plate or a
collector — feeding a 50 Ω scope line.

**a. His smallest source.** A 1 pF cathode at 1.414 kV (1 µJ) discharges
across the gap (a 20 Ω arc with 5 nH):

| | Value |
|---|---|
| Charge through the scope | 8.83×10⁹ electrons (all of C·V) |
| Energy | 1.000 µJ stored → 0.638 µJ to the scope, 0.362 µJ heat in the arc |
| Ledger | closes to 3×10⁻⁵ |

**b. An EV strike.** An ideal current source brings the EV's charge to the
anode in 20 ps:

| Electrons | Peak at the scope | Energy the arriving charge must bring | ÷ source energy |
|---|---|---|---|
| 2×10¹⁰ | 2.7 kV | 4.7 µJ | 5 |
| 10¹¹ | 13.4 kV | 117 µJ | 117 |

The ledger closes to 10⁻⁵ in both. Picture: right panel.

**What would settle it.** A Faraday cup on the anode and this ledger on
one shot: total charge collected against C·V of the source. Standard
physics predicts it never exceeds C·V. Shoulders' 10¹¹-per-bead figure
predicts more, from his smallest source.

## 6. The rings, and the mainstream counterpart

- **Ring marks.** Ring-shaped witness marks (2-19) fit a toroidal current
  ring. A real electron ring needs an outside field to hold it: the
  collapsing ring in [plasmoid_rings](../plasmoid_rings) needs a 2 T coil.
  Its electrons are 10⁶ times sparser than an EV's.
- **Ectons.** Mainstream physics has a measured phenomenon with the same
  fingerprints: G. A. Mesyats' **ectons**. These are explosive-emission
  micro-explosions at cathode spots, each lasting nanoseconds and
  releasing about 10¹¹–10¹² electrons. They leave craters and fast current
  pulses and give X-rays.
- **The difference.** An ecton is a neutral plasma of cathode metal, with
  ions. Shoulders reports almost no ions, so his ion measurement is the
  point where the two pictures disagree. It is worth repeating.

## Where it stands

- **The observations are real and repeatable.** Pits, rings, chains,
  X-rays and very fast pulses from tiny discharges, carefully documented.
- **The EV interpretation fails standard checks.** At the reported size
  and charge an EV:
  - holds 1,400 to 3 million times its source's energy;
  - flies apart in femtoseconds;
  - shows a field at odds with Gauss's law and his own kilovolt
    measurements;
  - has no known container. The zero-point field's measured effect is
    21 orders too weak and pushes outward for a sphere.
- **What would change this.** A one-shot Faraday-cup charge measurement
  exceeding C·V of the source. That is the H27 test.

Shoulders' later patents (5,018,180; 5,054,046; 5,054,047; 5,123,039;
5,148,461), and the Puthoff & Piestrup review, carry the same claims.

## Files

| File | What it holds |
|---|---|
| `ev.go` | the bead: self-energy, field, pressure, relativistic Coulomb explosion; Casimir, Boyer, Brillouin, Gauss radius |
| `circuit.go` | the kronspice bench: the source discharge and an EV strike, with ledgers |
| `main.go` | the report |
| `main_test.go` | electrostatics, explosion against the exact solution, containment, kronspice charge and ledger |
| `viz/render.py` | `shoulders_ev.png` |

## Sources

- K. R. Shoulders, *EV: A Tale of Discovery* (Jupiter Technologies, Austin, 1987).
- H. E. Puthoff, M. A. Piestrup, "Charge confinement by Casimir forces," arXiv:physics/0408114 (2004).
- T. H. Boyer, "Quantum electromagnetic zero-point energy of a conducting spherical shell and the Casimir model for a charged particle," *Phys. Rev.* 174, 1764 (1968).
- H. B. G. Casimir, *Proc. K. Ned. Akad. Wet.* 51, 793 (1948).
- L. Brillouin, "A theorem of Larmor and its importance for electrons in magnetic fields," *Phys. Rev.* 67, 260 (1945) — the density limit for magnetically confined electrons.
- G. A. Mesyats, "Ecton mechanism of the cathode spot phenomena in a vacuum arc," *Plasma Phys. Control. Fusion* 47, A109 (2005).
- D. Nakamura et al., "Record indoor magnetic field of 1200 T generated by electromagnetic flux-compression," *Rev. Sci. Instrum.* 89, 095106 (2018).
