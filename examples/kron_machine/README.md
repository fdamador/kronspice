# Kron's moving inductors: frames I and II, and a parametric pump

```bash
go run ./examples/kron_machine
```

1. **The rotor frame.** Kron's two-phase primitive machine (stator α, β;
   rotor a, b turning with θ) is transformed by `kron.RotatingFrame` to
   axes fixed to the stator: the inductances become constant and the rotor
   rows gain the speed matrix G (Kron, *Equivalent Circuits of Electric
   Machinery*, 1951, Epilogue eq. 9).
2. **Start-up of an induction motor** from rest on a balanced two-phase
   100 V, 50 Hz supply, solved in the holonomic frame I (L(θ)
   trigonometric) and in the rotor frame II(b). The two runs agree on speed
   and currents; each energy ledger closes (supply = heat + kinetic and
   magnetic energy).
3. **A parametric (temporal) inductance** L(t) = L0·(1 + m·cos 2ω0t), the
   same element with its coordinate driven at constant speed, pumping an
   LC tank at ω0: the envelope grows at Mathieu's rate m·ω0/4 less the
   damping R/(2L0), and the energy it gains comes in through the pump's
   mechanical port.

Details and the test results are in [doc/M9c1.md](../../doc/M9c1.md).
