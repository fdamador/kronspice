# Torn GA: diakoptics in a running evolution

A two-pole RC pre-filter (R1, C1, R2, C2 in E24/E12) is evolved in front of
a fixed plant, joined to it by the resistor RT, for a 5 kHz cutoff at the
plant's output. The same run (same seed) is made three times: AC by MNA, by
Kron's node-pair form untorn, and by Kron's diakoptics torn at RT with the
plant's factored pieces cached and shared by all workers (milestone M9b).

```bash
go run ./examples/torn_ga                      # 40-section LC line
go run ./examples/torn_ga -plant grid -n 12    # 12×12 RC mesh
go test ./examples/torn_ga -bench Generation   # one generation, fresh seeds
```

The program prints each run's time and best design and exits non-zero if
the three best designs differ. Results and the discussion of where the
torn path stands against sparse MNA are in [doc/M9b.md](../../doc/M9b.md).
