# Eckardt's self-excited homopolar generator (H15)

H. Eckardt, "Simulation of resonances in the homopolar generator"
(A.I.A.S./UPITEC, 13 February 2012) describes a twin-disk homopolar machine
that excites its own field and reports an energy surplus at a "giant
resonance" where the effective resistance vanishes. The example runs the
paper's equations as printed, then the same machine as a kronspice circuit
built on the N machine element, and compares their energy books (CLAUDE.md
§1, hypothesis H15; §3f). Companion: [homopolar](../homopolar) for the
Lorentz-force disk itself.

```
go run ./examples/eckardt_homopolar
go test ./examples/eckardt_homopolar
```

## The machine and the paper's equations

Two disks of radius r on one shaft, in series with R and a coil of N
turns, length l, core permeability μr (paper Fig. 1). The "static" disk
sits in a permanent-magnet field Bm; the "dynamic" disk sits in the coil's
own field B_I = μ0·μr·N·I/l, so the loop current excites it.

```
L·İ + R·I = ½·ω·r²·(Bm + B_I)          (6)
R_eff = R − ½·ω·r²·μ0·μr·N/l           (8)
ω_res = 2·R·l/(r²·μ0·μr·N)             (9)
J·ω̇ = α   (constant drive torque)      (12)
```

Parameters (p. 4): R 0.3 Ω, L 1 mH, N 800, l 5 cm, r 10 cm, μr 100,
Bm 0.1 T, J 0.01 kg·m², α 1 N·m. Eq. (9) gives ω_res = 29.84 rad/s
(paper: about 30).

Eq. (12) leaves out the torque the disk currents feel in the fields they
cross. The Lorentz force that gives the EMF ½ωr²B also gives the back
torque ½r²B·I per disk, so that ω·(torque) = EMF·I exactly:

```
J·ω̇ = α − ½·r²·(Bm + B_I)·I
```

In kronspice the machine is one N element with four windings (static disk
path, dynamic disk path, coil, magnet equivalent at 1 A), speed-matrix
entries K1 = ½r²Bm = 0.0005 V·s and K2 = ½r²μ0μrN/l = 0.01005 V·s/A, the
shaft a mechanical node with inertia J and the drive torque α. The back
torque is part of the element, not an added term.

## Results

**1. The paper's equations as printed** (RK4, no back torque): the speed
rises linearly, and once ω passes ω_res the current runs away: 0.56 A at
0.30 s, 211 A at 0.40 s, 1.2×10⁵ A at 0.45 s, over 10⁶ A at 0.463 s. The
energy into R and the coil is 1.42×10⁹ J against 10.7 J of drive work,
all of it spin-up. The missing 1.423×10⁹ J equals, to 1e-3, the work the
omitted back torque ½r²(Bm + B_I)·I would have taken from the shaft. The
"energy from spacetime" in these equations is the shaft work they leave out.

**2. The same machine in kronspice** (back torque inside the element),
L = 1 mH, 3 s:

| t | ω | I | B_I |
|---|---|---|---|
| 0.30 s | 30.08 rad/s | 0.579 A | 1.16 T |
| 0.60 s | 25.77 rad/s | 5.05 A | 10.2 T |
| 1.20 s | 27.37 rad/s | 12.69 A | 25.5 T |
| 3.00 s | 30.33 rad/s | 9.68 A | 19.5 T |

The current stays bounded and the speed settles at ω_res. Drive work
84.75 J = R 80.1 J + coil 0.047 J + spin-up 4.6 J + reference resistors
0.0002 J; the ledger closes to 3e-5 J. Steady state: α = K1·I + K2·I²
gives I = 9.949 A and P = α·ω_res = 29.84 W = R·I² (29.69 W). That needs
B_I = 20 T in the iron, which saturates near 1.5–2 T, so a real coil
limits the excitation well before. With the coil inductance the stated
geometry implies (50.5 H if B_I fills the disk area) the run is much
slower and the shaft reverses as the field hunts; the books still close.

**3. The paper's second run** (disk spun to 40 rad/s, then connected,
Figs. 6–9):

| Drive | Peak current | Speed after 0.5 s | Energy into R |
|---|---|---|---|
| α = 0 | 32.34 A at 0.057 s | 19.53 rad/s | 6.092 J = flywheel 6.092 J |
| α = 1 N·m | 46.32 A at 0.049 s | 25.50 rad/s | 19.197 J = flywheel 4.749 J + drive 14.466 J − coil 0.017 J |

The pulse is the flywheel's energy discharged through the self-excited
dynamo: what R receives is what the shaft gave (books close to 1e-8 J).

**Standard reading.** R_eff = 0 is the self-excitation threshold of a
series-wound dynamo (Siemens 1867; Bullard 1955; Moffatt 1978, ch. 6).
Above it the current grows by taking power from the drive. The classical
resonance is reproduced; the paper's two ECE-only resonances are not
modeled because no equations are supplied for them. Status in the register:
*untested*.

## Files

| File | What it holds |
|---|---|
| `main.go` | paper parameters, `RunPaper` (eqs. 6 and 12 as printed), the N-element circuit, `RunCircuit`, steady state |
| `main_test.go` | ω_res, the missing energy = back-torque work, bounded current with closed books, flywheel discharge; `ECKARDT_DUMP=1` prints the comparison series |
| `export_test.go` | netlists, KiCad schematics and results sheets in `export/` |

## Sources

- H. Eckardt, "Simulation of resonances in the homopolar generator," A.I.A.S./UPITEC, 13 February 2012 (`reference/eckardt_2012_homopolar_resonances.pdf`). Hypothesis source.
- E. C. Bullard, "The stability of a homopolar dynamo," *Proc. Camb. Phil. Soc.* 51, 744 (1955).
- H. K. Moffatt, *Magnetic Field Generation in Electrically Conducting Fluids*, Cambridge 1978, ch. 6.
- G. Kron, *Equivalent Circuits of Electric Machinery*, Wiley 1951, Epilogue (the speed-voltage tensor).
