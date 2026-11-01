# Ring-tower animations

Self-contained HTML pages (open in any browser, no server needed):

- `four_towers.html` — the four-tower ring (A 0°, C 90°, B 180°, D 270°):
  each pair alone, both synchronized, C–D in antiphase, C–D shifted 1/8 wave.
- `phase_control.html` — sliders for the C–D timing and position phase and
  the A–B timing phase, with the quarter-point voltages and the towers' power.
- `eight_towers.html` — eight towers at 45° spacing (pairs A–B, C–D, E–F,
  G–H): layout with each tower's two radiators and launch direction, timing
  and position phase and on/off per pair, power per pair (`TestEightTowers`).
- `dc_bias.html` — a DC bias on each pair, on a continuous ring (one DC
  level, steady current between the supplies) or a sectored ring (blocking
  capacitors; each pair holds its own level).

Data: every page embeds `launchers.json`, the steady-state voltage phasors
of each launcher (A, C, E, G east; B, D, F, H west) simulated alone by kronspice:

    go run ./examples/ring_towers -phasors examples/ring_towers/viz/launchers.json

The ring is linear, so each page sums those four responses with the phases
chosen in the browser; this matches full simulations of the combined cases
to 0.25 % (`TestFourTowers`, `TestBias`). The DC levels in `dc_bias.html`
follow `DCLevels` and the current between supplies `(VAB − VCD)/RFeed`,
both checked against the kronspice operating point in `TestBias`.

`eight_towers.html` embeds the A-east and B-west responses and rotates
them to the other six positions (the ring is uniform, so each launcher's
response is A's or B's shifted around the ring, with sign −1 for the 45°
set's 180° base phase; checked equal to the simulated responses).

To rebuild the first three pages after regenerating the data:

    python3 examples/ring_towers/viz/build_viz.py examples/ring_towers/viz/launchers.json examples/ring_towers/viz
