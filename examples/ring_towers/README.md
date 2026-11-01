# Ring towers: standing waves on a closed ring

A thought experiment as a circuit: towers at points around the
Earth–ionosphere ring drive the same frequency into it, each tower sending
one wave east and one west, and their phases are chosen so that the waves
add, cancel or circulate. The ring is scaled to MHz so each lap takes
microseconds; the patterns depend only on the number of wavelengths around
the ring and on the phases. The same experiment inside a craft's double
wall is [double_hull](../double_hull). Every joule in the ring is supplied
by the towers and the ledger says so (CLAUDE.md §3c).

```
go run ./examples/ring_towers [-json out.json] [-phasors viz/launchers.json]
go test ./examples/ring_towers
```

## The construction

The ring is a lumped transmission line (telegrapher's equations;
Steinmetz Lecture VI, Paul §2): 96 sections of series L = 1 µH and R with
shunt C = 1 nF, the last closing back on the first. Z0 = 31.6 Ω, loss
0.20 Np per lap (R = 0.132 Ω per section), run at its fourth resonance,
f = 1.314 MHz, four wavelengths around, so waves returning after each lap
add to the new ones and the build-up is set by the loss per lap.

Each tower is two 10 mA current sources a quarter wavelength (6 sections)
apart, a two-element end-fire array: the eastern element lagging 90°
launches east only, leading 90° launches west only. Tower B sits at the
antipode. Each state is run to steady state (8 amplitude-decay times) by
transient analysis.

**Two towers.** *Traveling*: both launch east, B phased to the arrival of
A's wave, so one wave circulates with the same amplitude everywhere.
*Standing*: A launches east, B west, so counter-running waves form fixed
nodes and antinodes with no net flow. *Cancel*: both launch east with B
180° from A's arrival, so the waves subtract to nothing.

**Four towers.** A second pair C (90°) and D (270°) is added. With four
wavelengths around the ring every quarter point is an antinode of A–B's
cos(kx) pattern, all swinging in step, so phasing C–D to the same reference
waves puts its pattern exactly on top of A–B's and the two add. Two
contrasts: C–D in antiphase cancels A–B, and C–D shifted by an eighth of a
wavelength moves the antinodes off the tower points.

**Eight towers** add the 45°-rotated pairs E–F and G–H.

**DC bias.** A supply per pair through 1 MΩ feeds. The ring is one
conductor at DC, so on a continuous ring the pairs cannot hold different
levels: the ring settles at the average and (V_AB − V_CD)/R_feed flows
between the supplies. A sectored ring puts 10 µF blocking capacitors at
45°, 135°, 225° and 315° (reactance 0.012 Ω against ωL = 8.3 Ω per
section), so each sector holds its own bias while the RF passes.

## Results

| State | Peak | Envelope min/max | Eastward power flow | Towers deliver = ring loss | Stored |
|---|---|---|---|---|---|
| traveling | 3.294 V | 0.925 | 0.155 W | 0.06255 W = 0.06254 W | 0.475 µJ |
| standing | 3.328 V | 0.018 | 0.000163 W | 0.03136 W = 0.03135 W | 0.238 µJ |
| cancel | 0.167 V | 0.904 | 0.0003 W | 0.000161 W = 0.0001609 W | 0.00122 µJ |

Ledger mismatch is at most 1.1e-9 of the energy in.

| Four towers | \|V\| at 0°, 90°, 180°, 270° | Peak | Towers | Stored |
|---|---|---|---|---|
| A–B standing | 3.012, 3.324, 3.328, 3.008 V | 3.328 V | 0.03136 W | 0.238 µJ |
| C–D standing | 3.008, 3.012, 3.324, 3.328 V | 3.328 V | 0.03136 W | 0.238 µJ |
| both, synchronized | 6.023, 6.339, 6.656, 6.340 V | 6.656 V | 0.1252 W | 0.951 µJ |
| both, C–D antiphase | 0.004, 0.315, 0.004, 0.323 V | 0.323 V | 0.0002405 W | 0.00183 µJ |
| both, C–D shifted 45° | 5.127, 5.433, 5.669, 5.338 V | 6.088 V | 0.1069 W | 0.812 µJ |

Synchronized pairs give twice the amplitude and four times the stored
energy and tower power; eight synchronized towers give 11.95–13.22 V
around the ring from 0.5006 W, sixteen times one pair's stored energy
(3.8 µJ), and the 45° set in antiphase cancels everything to 0.32 V.

DC bias with A–B at 100 V and C–D at 50 V: on the continuous ring every
quarter point sits at 75.0 V with 50 µA between the supplies (2.5 mW in
the feeds); on the sectored ring the points hold 100.0, 50.0, 100.0 and
50.0 V with no DC current. In both the RF rides unchanged on top
(5.97–6.60 V).

## Animation

`viz/` holds self-contained HTML pages (four towers, phase sliders, eight
towers, DC bias). Each embeds `launchers.json`, the steady-state phasors of
each launcher simulated alone; the ring is linear, so the browser sums
them with the chosen phases, matching full simulations to 0.25 %. See
`viz/README.md` for how to regenerate them.

## Files

| File | What it holds |
|---|---|
| `main.go` | the ring, the towers, the states, the bias, the steady-state runs and the launcher phasors |
| `main_test.go` | two, four and eight towers and the bias, each against the expected pattern with the ledger closed |
| `export_test.go` | four-tower netlists (RF and sectored bias), KiCad schematics and results sheets in `export/` |
| `viz/` | the animations and their builder |
| `run_output.txt` | the run above |

## Sources

- C. P. Steinmetz, *Elementary Lectures on Electric Discharges, Waves and Impulses*, Lecture VI (1911) (`reference/`).
- C. R. Paul, *Analysis of Multiconductor Transmission Lines*, 2nd ed., §2 (`reference/`).
