// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package evo

// NonDominated returns the indices of the points not dominated by any other
// (all objectives minimized): the Pareto front. A point dominates another
// when it is no worse in every objective and better in at least one.
// Duplicates keep the first occurrence. Order follows the input.
func NonDominated(points [][]float64) []int {
	var out []int
	for i, p := range points {
		dominated := false
		for j, q := range points {
			if i == j {
				continue
			}
			le, lt := true, false
			for k := range p {
				if q[k] > p[k] {
					le = false
					break
				}
				if q[k] < p[k] {
					lt = true
				}
			}
			if le && (lt || j < i) {
				dominated = true
				break
			}
		}
		if !dominated {
			out = append(out, i)
		}
	}
	return out
}
