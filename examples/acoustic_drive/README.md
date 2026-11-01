# Acoustic drive of a sonoluminescence flask

Electrical → mechanical → acoustic, one circuit: a 50 Ω generator, an
optional series inductor that tunes out the transducer's capacitance, a
piezoelectric ring (clamped capacitance, force factor, moving mass, loss),
its coupling into the breathing mode of a 3 cm liquid sphere, and the mode
itself, in H₂O and D₂O. The output is the center pressure per volt, the
drive and power for a 1.3 atm amplitude, and where every watt goes.

```bash
go run ./examples/acoustic_drive
```

The domains, transducers and liquid models are described in
[doc/domains.md](../../doc/domains.md), with the results table. The flask Q,
the coupling area and the transducer values in `main.go` are assumptions to
be replaced by bench values; the bubble itself is milestone M9p.
