// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"context"
	"math"
	"testing"

	"github.com/fdamador/kronspice/evo"
)

// M9n: the Lorentz-force EMF equals ½Bω(r_o² − r_i²); open-circuited,
// Tesla's two disks give twice one disk's EMF at the terminals; loaded,
// the shaft power from the ledger equals Σ torque·ω and the electrical
// power plus the losses (the ledger closes); and a short GA run meets the
// power target with an efficiency near 1.
func TestHomopolar(t *testing.T) {
	d := Defaults()
	want := 0.5 * d.B * d.Omega * (d.Ro*d.Ro - d.Ri*d.Ri)
	if math.Abs(d.EMF()-want) > 1e-12*want {
		t.Errorf("EMF %g, want %g", d.EMF(), want)
	}
	open := d
	open.Load = 1e6
	r, _, err := open.Solve()
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(math.Abs(r.VLoad)-2*want) > 1e-6*want {
		t.Errorf("two-disk open-circuit voltage %g, want 2·EMF = %g", r.VLoad, 2*want)
	}
	r, out, err := d.Solve()
	if err != nil {
		t.Fatal(err)
	}
	shaft := (r.Torque[0] + r.Torque[1]) * d.Omega
	if math.Abs(shaft-r.PShaft) > 1e-9*r.PShaft {
		t.Errorf("Σ torque·ω = %g W, ledger shaft inflow %g W", shaft, r.PShaft)
	}
	if rr := out.Ledger.RelativeResidual(); rr > 1e-12 {
		t.Errorf("ledger residual %g", rr)
	}
	// Current I = 2·EMF/R_total and torque per disk K·I with K = ½B(r_o²−r_i²).
	rtot := d.Load + 2*d.RShaft + 2*d.DiskResistance() + d.RBelt
	i := 2 * want / rtot
	if k := 0.5 * d.B * (d.Ro*d.Ro - d.Ri*d.Ri); math.Abs(r.Torque[0]-k*i) > 1e-9*k*i {
		t.Errorf("torque %g N·m, want K·I = %g", r.Torque[0], k*i)
	}
	single := d
	single.TwoDisk = false
	r1, _, err := single.Solve()
	if err != nil {
		t.Fatal(err)
	}
	if !(r.PLoad > 3*r1.PLoad) {
		t.Errorf("two-disk %g W vs single disk %g W into the same load", r.PLoad, r1.PLoad)
	}
	s, err := evo.Run(context.Background(), Problem(2000), evo.Options{Population: 16, Generations: 6, Seed: 1}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	if s.Best.Terms[1] > 0 || s.Best.Terms[0] > 0.05 {
		t.Errorf("evolved: 1 − η = %g, power shortfall %g", s.Best.Terms[0], s.Best.Terms[1])
	}
}
