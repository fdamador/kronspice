// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"os"
)

func main() {
	if err := run(); err != nil {
		fmt.Fprintln(os.Stderr, "barrett_osc:", err)
		os.Exit(1)
	}
}

func run() error {
	if err := printFig2(); err != nil {
		return err
	}
	if err := printFig3(); err != nil {
		return err
	}
	if err := printFig4(); err != nil {
		return err
	}
	return printFig5()
}
