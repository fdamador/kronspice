// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"encoding/csv"
	"io"
	"strconv"
	"strings"
)

// WriteCSV writes an operating-point result as "name,value" rows:
// v(node) for nodes, i(element) and p(element) for elements, and
// ledger(element,channel) for every modeled ledger row.
func WriteCSV(w io.Writer, r *Result) error {
	cw := csv.NewWriter(w)
	f := func(x float64) string { return strconv.FormatFloat(x, 'g', 12, 64) }
	rows := [][]string{{"name", "value"}}
	for i, n := range r.NodeNames {
		if i == 0 {
			continue
		}
		rows = append(rows, []string{"v(" + n + ")", f(r.NodeVoltages[i])})
	}
	for _, e := range r.Elements {
		name := strings.ToLower(e.Name)
		rows = append(rows, []string{"i(" + name + ")", f(e.Current)}, []string{"p(" + name + ")", f(e.PowerAbsorbed)})
	}
	if r.Ledger != nil {
		for _, row := range r.Ledger.Rows {
			if row.Status == 0 { // not modeled
				continue
			}
			rows = append(rows, []string{"ledger(" + strings.ToLower(row.Element) + "," + row.Channel.String() + ")", f(row.PowerW)})
		}
		rows = append(rows,
			[]string{"ledger.inflow", f(r.Ledger.Inflow)},
			[]string{"ledger.outflow", f(r.Ledger.Outflow)},
			[]string{"ledger.dissipation", f(r.Ledger.Dissipation)},
			[]string{"ledger.residual", f(r.Ledger.Residual)},
		)
	}
	if err := cw.WriteAll(rows); err != nil {
		return err
	}
	return cw.Error()
}
