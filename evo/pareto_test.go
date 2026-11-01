// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

import (
	"reflect"
	"testing"
)

func TestNonDominated(t *testing.T) {
	pts := [][]float64{{1, 5}, {2, 2}, {3, 3}, {5, 1}, {2, 2}, {4, 4}}
	if got := NonDominated(pts); !reflect.DeepEqual(got, []int{0, 1, 3}) {
		t.Fatalf("front = %v", got)
	}
}
