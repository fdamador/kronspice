// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package circuit

import (
	"errors"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/netlist"
)

// A P line parses, validates, writes and parses back unchanged; missing
// measured coefficients are an error, never a default.
func TestBearingNetlist(t *testing.T) {
	src := "bb\nI1 0 a 50\nP1 a 0 w 0 2m KAPPA=4u W0=20 VF=0.73 THERMAL=1 ETAMAX=0.1\nRD w 0 100\n.domain mechanical w\n"
	d, err := netlist.ParseString(src)
	if err != nil {
		t.Fatal(err)
	}
	c, err := FromDeck(d)
	if err != nil {
		t.Fatal(err)
	}
	var sb strings.Builder
	if err := c.WriteNetlist(&sb); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(sb.String(), "P1 a 0 w 0 0.002 ETAMAX=0.1 KAPPA=4e-06 THERMAL=1 VF=0.73 W0=20") {
		t.Errorf("written:\n%s", sb.String())
	}
	d2, err := netlist.ParseString(sb.String())
	if err != nil {
		t.Fatal(err)
	}
	if _, err := FromDeck(d2); err != nil {
		t.Fatal(err)
	}
	for _, bad := range []string{
		"P1 a 0 w 0 2m W0=20",                     // no KAPPA
		"P1 a 0 w 0 2m KAPPA=4u",                  // no W0
		"P1 a 0 w 0 2m KAPPA=4u W0=20 THERMAL=1",  // thermal without ETAMAX
		"P1 a 0 w 0 2m KAPPA=4u W0=20 ETAMAX=0.1", // ETAMAX without thermal
		"P1 a 0 w 0 0 KAPPA=4u W0=20",             // R = 0
		"P1 a 0 w 0 2m KAPPA=4u W0=20 FOO=1",      // unknown
	} {
		d, err := netlist.ParseString("bb\nI1 0 a 50\n" + bad + "\nRD w 0 100\n.domain mechanical w\n")
		if err != nil {
			continue
		}
		if _, err := FromDeck(d); !errors.Is(err, ErrBadValue) {
			t.Errorf("%q: err = %v, want ErrBadValue", bad, err)
		}
	}
}
