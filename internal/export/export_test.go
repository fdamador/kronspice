// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package export

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Examples that ship netlists (configs for cmd/evolve, and the generated
// netlists of water_vibration and kron_schrodinger) get the same export as
// the circuits built in Go: examples/NAME/export/, KRONSPICE_EXPORT=1.
func TestExport(t *testing.T) {
	var files []string
	for _, g := range []string{"../../examples/*/*.cir", "../../examples/*/netlists/*.cir"} {
		m, _ := filepath.Glob(g)
		files = append(files, m...)
	}
	if len(files) < 10 {
		t.Fatalf("only %d example netlists", len(files))
	}
	for _, f := range files {
		src, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		ex := filepath.Base(filepath.Dir(f))
		if ex == "netlists" {
			ex = filepath.Base(filepath.Dir(filepath.Dir(f)))
		}
		name := ex + "_" + strings.TrimSuffix(filepath.Base(f), ".cir")
		Test(t, filepath.Join("../../examples", ex, "export"), Item{Name: name, Src: string(src)})
	}
}
