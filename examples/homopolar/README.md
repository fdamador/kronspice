# Tesla's homopolar dynamo from the Lorentz force

Tesla's "Dynamo-Electric Machine," U.S. Patent 406,968 (1889): two
conducting disks on one shaft in fields of opposite polarity, joined at
their rims by a flexible metal belt, with the current taken from the two
shafts so no brush bears on a fast rim. The example derives the disk's EMF
from the Lorentz force by grade selection (CLAUDE.md §3f, milestone M9n),
builds the machine from kronspice's N machine element, compares Tesla's
two-disk arrangement with a single disk and a rim brush, and evolves the
two-disk machine for efficiency at a required output. The versions with
the current through the shaft bearings are in
[homopolar_bearing](../homopolar_bearing); the self-excited twin-disk case
of H15 is in [eckardt_homopolar](../eckardt_homopolar).

```
go run ./examples/homopolar [-gens 20] [-pop 30] [-p 2000]
go test ./examples/homopolar
```

## The construction

A disk spinning at ω in an axial field B carries its charges at
v = ω×r. The force q(v×B), the grade-1 part of ⟨F·Q†⟩ (`em.LorentzGA`),
pushes them radially, so the disk is a source of EMF between axle and rim.
`em.MotionalEMF` integrates v×B along a radius and must equal

```
EMF = ½ B ω (r_o² − r_i²)
```

In the circuit each disk is one N element with two windings, the disk's
radial path and the field winding (driven at 1 A, the field strength in
the coupling), and one speed-matrix entry K = ½·B·(r_o² − r_i²). The
element gives the speed voltage ω·K·i_f in the disk and the torque
K·i_disk·i_f on the shaft. Each shaft is a mechanical node driven at ω, so
the ledger books the shaft power as inflow and shows it leaving as heat in
the load, the disks, the belt and the contacts.

Reference design: 15 cm copper disks, 5 mm thick, 1 cm axle, 0.8 T,
3000 r/min, 10 mΩ load. Contact and belt resistances (50 µΩ shaft, 500 µΩ
rim brush, 100 µΩ belt) are declared, illustrative inputs, not measured.

## Results

| | Two disks, belt (Tesla) | Single disk, rim brush |
|---|---|---|
| EMF per disk, Lorentz force | 2.8149 V | 2.8149 V |
| ½Bω(r_o² − r_i²) | 2.8149 V | 2.8149 V |
| Disk spreading resistance | 1.45 µΩ | 1.45 µΩ |
| Load | 5.5178 V, 3044.6 W | 2.6678 V, 711.7 W |
| Shaft (ledger inflow) | 3106.4 W, torques 4.94 + 4.94 N·m | 750.9 W, torque 2.39 N·m |
| Efficiency | 0.9801 | 0.9477 |
| Ledger residual | 0 | 0 |

The two disks' EMFs add (open circuit: 2 × EMF at the terminals, to 1e-6),
and the belt costs far less than a brush on a rim moving at 47 m/s, so the
two-disk machine puts more than four times the power into the same load.

The tests check that the Lorentz EMF equals the formula to 1e-12, that the
shaft power from the ledger equals Σ torque·ω to 1e-9, that the torque per
disk is K·I with I = 2·EMF/R_total, and that the ledger residual is below
1e-12.

**Evolved** (genes: outer radius, thickness, field up to 1.5 T, load) for
at least 2 kW with the least loss, 20 generations of 30, seed 1 (562
evaluations): r_o 0.361 m, 6.8 mm thick, 1.22 T, 1 Ω load, giving
24.93 V per disk, 2485 W into the load from 2486 W at the shaft,
efficiency 0.9998.

## Files

| File | What it holds |
|---|---|
| `main.go` | `Design`, the Lorentz EMF, the N-element circuit, the GA problem |
| `main_test.go` | EMF formula, open-circuit sum, torque, ledger, single-disk comparison, short GA run |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |

## Sources

- N. Tesla, "Dynamo-Electric Machine," U.S. Patent 406,968, 1889 (`reference/`).
- P. Joot, *Geometric Algebra for Electrical Engineers*, Theorem 3.16 (Lorentz force and power by grade selection; reference only, see CLAUDE.md §3f).
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley 1951, Epilogue eq. 9 (the speed-voltage tensor G).
- `doc/M9n.md` for the milestone record.
