// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"flag"
	"fmt"
	"os"
)

var dataDir = flag.String("data", "", "write figure data (CSV) to this directory")

func main() {
	flag.Parse()
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
	if err := printFig5(); err != nil {
		return err
	}
	cs, err := printBench()
	if err != nil {
		return err
	}
	tc, err := printTesla()
	if err != nil {
		return err
	}
	if *dataDir != "" {
		if err := writeBenchData(*dataDir, cs); err != nil {
			return err
		}
		if err := writeTeslaData(*dataDir, tc); err != nil {
			return err
		}
	}
	return nil
}
