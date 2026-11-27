// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

//go:build race

package evo

// raceOn is true under the race detector, where the long evolution
// acceptance runs (minutes at full speed) would take an hour; they are
// skipped there and run in the ordinary test pass.
const raceOn = true
