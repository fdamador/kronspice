// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package analysis

import (
	"math"
	"testing"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// FuzzSimulate feeds arbitrary netlists through parsing, validation and
// simulation. Any outcome is acceptable except a panic, a hang or a
// non-finite result (CLAUDE.md §1 goal 2, §8 rule 6).
func FuzzSimulate(f *testing.F) {
	f.Add("t\nV1 a 0 5\nR1 a b 1k\nR2 b 0 1k\n.op\n")
	f.Add("t\nV1 in 0 PULSE(0 1 0 1n 1n 1u 2u)\nR1 in out 1k\nC1 out 0 1n\n.tran 10n 5u\n")
	f.Add("t\nC1 a 0 1u IC=1\nL1 a 0 1m\n.tran 1u 100u uic\n")
	f.Add("t\nV1 a 0 1\nV2 a 0 2\n.op\n")
	f.Add("t\nI1 0 a pwl(0 0 1m 1)\nR1 a 0 1\nL1 a 0 1\n.tran 0.1m 2m\n")
	f.Add("t\nV1 in 0 AC 1\nR1 in out 1k\nC1 out 0 1u\n.ac dec 10 1 1meg\n")
	f.Add("t\nV1 in 0 1\nR1 in 0 1k\n.dc V1 0 5 1\n")
	f.Add("t\nV1 a 0 100\nD1 a 0 d\n.model d D(IS=1e-14 RS=0)\n.op\n")
	f.Add("t\nV1 in 0 0\nR1 in o 1k\nD1 o 0 z\n.model z D(BV=3 IBV=1m CJO=1p)\n.dc V1 -9 9 1\n")
	f.Add("t\nV1 in 0 SIN(0 5 1k)\nD1 in o d 2\nC1 o 0 1u\nR1 o 0 1k\n.model d D(TT=1n)\n.tran 10u 2m\n")
	f.Add("t\nV1 in 0 DC 0.7 AC 1\nD1 in 0 d\n.model d D(N=2 CJO=1p M=0.9 FC=0.99)\n.ac dec 3 1 1g\n")
	f.Add("t\nV1 in 0 0.4\nR1 in a 50\nB1 a 0 I=pwl(V(a), 0,0, 0.1,10m, 0.3,1m, 0.6,20m)\n.dc V1 0 1 0.05\n")
	f.Add("t\nV1 in 0 100\nR1 in a 10k\nVs a m 0\nB1 m 0 V=pwl(I(Vs), 0,0, 1u,90, 1m,60, 10m,70)\n.op\n")
	f.Add("t\nV1 in 0 SIN(0 1 1k)\nB1 o 0 V=pwl(V(in), -1,-2, 1,2)\nR1 o 0 1k\nC1 o 0 1u\n.tran 10u 2m\n")
	f.Add("t\nVCC c 0 5\nRB c b 100k\nQ1 c b 0 q\n.model q NPN(BF=100 VAF=50 IKF=10m RB=10 CJE=1p TF=1n)\n.op\n")
	f.Add("t\nV1 in 0 PULSE(0 5 0 1n 1n 50n 100n)\nRB in b 1k\nRC v c 1k\nV2 v 0 5\nQ1 c b 0 0 q\n.model q NPN(TF=1n TR=10n CJC=1p CJS=1p)\n.tran 1n 200n\n")
	f.Add("t\nV1 e 0 5\nQ1 c b e q\nRC c 0 1k\nVB b 0 4.3\n.model q PNP(BF=50 RE=2)\n.dc VB 0 5 0.1\n")
	f.Add("t\nV1 a 0 AC 1\nR1 a 0 1k\nE1 o 0 a 0 -100\nG1 0 p o 0 1m\nRp p 0 1k\n.ac dec 2 1 1k\n")
	f.Add("t\nV1 a 0 1\nVS a m 0\nR1 m 0 1k\nF1 0 f VS 2\nRF f 0 1k\nH1 h 0 VS 1k\nRH h 0 1k\n.op\n")
	f.Add("t\nVDD vdd 0 5\nVIN in 0 PULSE(0 5 1n 1n 1n 5n 12n)\nMN out in 0 0 n L=1u W=5u AD=10p\nMP out in vdd vdd p L=1u W=10u\nCL out 0 10f\n.model n NMOS(TOX=20n NSUB=1e16 CJ=1e-4 CGSO=1e-10)\n.model p PMOS(VTO=-0.7 KP=20u)\n.tran 0.1n 20n\n")
	f.Add("pz\n.domain mechanical v\n.domain acoustic p\nV1 e 0 DC 0 AC 1 SIN(0 1 25k)\nC0 e 0 4n\nY1 e 0 v 0 0.6\nCM v 0 4m\nLK v 0 1u\nRB v 0 0.3\nY2 v 0 p 0 2e-4\nCP p 0 1e-14\nLP p 0 4e6\nRP p 0 1e9 EXT=acoustic\n.ac lin 3 20k 30k\n.tran 1u 100u\n")
	f.Add("k\nV1 a 0 PULSE(0 1 0 1n 1n 1u 2u) AC 1\nR1 a p 10\nL1 p 0 1m\nL2 s 0 4m\nK1 L1 L2 0.95\nR2 s 0 100\nC2 s 0 1n\n.tran 10n 3u\n.ac dec 3 1k 1g\n")
	f.Add("t\nVDD d 0 3\nVG g 0 DC 1.5 AC 1\nM1 d g 0 0 n L=2u W=2u\n.model n NMOS(VTO=0.6 KP=100u LAMBDA=0.1 RD=10)\n.ac dec 3 1 1g\n")
	f.Add("t\nVDD vdd 0 5\nVIN in 0 PULSE(0 5 1n 1n 1n 5n 12n)\nMN out in 0 0 n L=2u W=5u\nMP out in vdd vdd p L=2u W=10u\nCL out 0 10f\n.model n NMOS(LEVEL=2 TOX=25n NSUB=2e16 XJ=0.3u NFS=1e11 VMAX=5e4 UCRIT=5e4 UEXP=0.1 DELTA=1)\n.model p PMOS(LEVEL=2 VTO=-0.7 KP=20u GAMMA=0.4)\n.tran 0.1n 20n\n")
	f.Add("t\nVDD d 0 3\nVG g 0 DC 1.2 AC 1\nVB b 0 -1\nM1 d g 0 b n L=1u W=4u\n.model n NMOS(LEVEL=3 TOX=20n NSUB=3e16 XJ=0.25u ETA=0.3 THETA=0.1 KAPPA=0.4 VMAX=1.5e5 NFS=4e11 DELTA=0.5 WD=0.05u)\n.ac dec 3 1 1g\n.dc VG 0 3 0.5\n")
	f.Fuzz(func(t *testing.T, src string) {
		d, err := netlist.ParseString(src)
		if err != nil {
			return
		}
		c, err := circuit.FromDeck(d)
		if err != nil {
			return
		}
		for _, a := range d.Analyses {
			switch a.Kind {
			case "op":
				r, err := OP(c, Options{})
				if err != nil {
					continue
				}
				for _, v := range r.NodeVoltages {
					if math.IsNaN(v) || math.IsInf(v, 0) {
						t.Fatalf("non-finite node voltage")
					}
				}
			case "ac":
				spec := ACSpec{Variation: a.Variation, Points: a.Points, FStart: a.FStart, FStop: a.FStop}
				if f, err := spec.Frequencies(); err != nil || len(f) > 5000 {
					continue
				}
				_, _ = AC(c, spec, Options{})
			case "dc":
				var sw []DCSweep
				for _, s := range a.Sweeps {
					if n := (s.Stop - s.Start) / s.Step; n > 1000 {
						return
					}
					sw = append(sw, DCSweep{Source: s.Source, Start: s.Start, Stop: s.Stop, Step: s.Step})
				}
				_, _ = DC(c, sw, Options{})
			case "tran":
				// Bound the work per input so the fuzzer cannot hang.
				if a.Stop/a.Step > 1e5 {
					continue
				}
				_, _ = Tran(c, TranSpec{Step: a.Step, Stop: a.Stop, Start: a.Start, MaxStep: a.MaxStep, UIC: a.UIC},
					TranOptions{MaxSteps: 20000})
			}
		}
	})
}
