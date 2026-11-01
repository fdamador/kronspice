// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Package golden loads reference results and compares simulation output
// against them using the tolerances in CLAUDE.md §8 rule 5.
package golden

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/fdamador/kronspice/output"
)

// Tolerances for DC operating points (CLAUDE.md §8 rule 5).
const (
	AbsVoltage = 1e-6 // V
	AbsCurrent = 1e-12
	Rel        = 1e-3
)

// Load reads a golden CSV: lines "name,value"; "#" lines and the header are skipped.
func Load(path string) (map[string]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return Read(f)
}

// Read parses golden CSV data.
func Read(r io.Reader) (map[string]float64, error) {
	out := map[string]float64{}
	sc := bufio.NewScanner(r)
	n := 0
	for sc.Scan() {
		n++
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") || line == "name,value" {
			continue
		}
		parts := strings.SplitN(line, ",", 2)
		if len(parts) != 2 {
			return nil, fmt.Errorf("golden: line %d: want name,value", n)
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(parts[1]), 64)
		if err != nil {
			return nil, fmt.Errorf("golden: line %d: %v", n, err)
		}
		out[strings.ToLower(strings.TrimSpace(parts[0]))] = v
	}
	return out, sc.Err()
}

// Values converts an operating-point result into the golden naming scheme:
// v(node) for every non-ground node and i(name) for every element: the
// current into its + terminal, as ngspice reports it for voltage sources
// (ngspice files carry only the voltage-source currents).
func Values(r *output.Result) map[string]float64 {
	out := map[string]float64{}
	for i, n := range r.NodeNames {
		if i == 0 {
			continue
		}
		out["v("+n+")"] = r.NodeVoltages[i]
	}
	for _, e := range r.Elements {
		out["i("+strings.ToLower(e.Name)+")"] = e.Current
	}
	return out
}

// Compare checks every golden value against got and returns one message per
// mismatch (empty when everything agrees).
func Compare(got, want map[string]float64) []string {
	var msgs []string
	keys := make([]string, 0, len(want))
	for k := range want {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		w := want[k]
		g, ok := got[k]
		if !ok {
			msgs = append(msgs, fmt.Sprintf("%s: missing (want %g)", k, w))
			continue
		}
		abs := AbsVoltage
		if strings.HasPrefix(k, "i(") {
			abs = AbsCurrent
		}
		if tol := math.Max(abs, Rel*math.Abs(w)); math.Abs(g-w) > tol {
			msgs = append(msgs, fmt.Sprintf("%s: got %.9g, want %.9g (tolerance %.3g)", k, g, w, tol))
		}
	}
	return msgs
}

// Table is a transient or AC golden reference: an independent-variable
// column (time or freq) and named series.
type Table struct {
	X       string    // "time", "freq" or "sweep" (DC sweep value)
	Time    []float64 // values of the independent variable
	Columns map[string][]float64
	Names   []string // column names in file order (excluding time)
}

// IsTable reports whether a golden file holds a time table (header starts
// with "time,") rather than name,value pairs.
func IsTable(path string) (bool, error) {
	f, err := os.Open(path)
	if err != nil {
		return false, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		low := strings.ToLower(line)
		return strings.HasPrefix(low, "time,") || strings.HasPrefix(low, "freq,") || strings.HasPrefix(low, "sweep,"), nil
	}
	return false, sc.Err()
}

// LoadTable reads a transient golden CSV: a "time,name1,name2…" header then rows.
func LoadTable(path string) (*Table, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	tb := &Table{Columns: map[string][]float64{}}
	n := 0
	for sc.Scan() {
		n++
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		parts := strings.Split(line, ",")
		if tb.Names == nil {
			x := strings.ToLower(strings.TrimSpace(parts[0]))
			if x != "time" && x != "freq" && x != "sweep" {
				return nil, fmt.Errorf("golden: line %d: table header must start with time, freq or sweep", n)
			}
			tb.X = x
			for _, p := range parts[1:] {
				tb.Names = append(tb.Names, strings.ToLower(strings.TrimSpace(p)))
			}
			continue
		}
		if len(parts) != len(tb.Names)+1 {
			return nil, fmt.Errorf("golden: line %d: want %d columns", n, len(tb.Names)+1)
		}
		vals := make([]float64, len(parts))
		for i, p := range parts {
			v, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
			if err != nil {
				return nil, fmt.Errorf("golden: line %d: %v", n, err)
			}
			vals[i] = v
		}
		tb.Time = append(tb.Time, vals[0])
		for i, name := range tb.Names {
			tb.Columns[name] = append(tb.Columns[name], vals[i+1])
		}
	}
	return tb, sc.Err()
}

// TranFullScale is the transient tolerance as a fraction of each series'
// full scale (largest |value|), per CLAUDE.md §8 rule 5.
const TranFullScale = 2e-3

// CompareTable resamples the simulated series onto the golden time points and
// returns one message per series whose worst error exceeds the tolerance.
// get returns the simulated series for a golden column name (v(node), i(elem)).
func CompareTable(simTime []float64, get func(name string) ([]float64, bool), want *Table) []string {
	var msgs []string
	for _, name := range want.Names {
		y, ok := get(name)
		if !ok {
			msgs = append(msgs, fmt.Sprintf("%s: not in simulation", name))
			continue
		}
		ref := want.Columns[name]
		full := 0.0
		for _, v := range ref {
			full = math.Max(full, math.Abs(v))
		}
		tol := math.Max(TranFullScale*full, AbsVoltage)
		worst, at := 0.0, 0.0
		for k, tq := range want.Time {
			if d := math.Abs(output.Interp(simTime, y, tq) - ref[k]); d > worst {
				worst, at = d, tq
			}
		}
		if worst > tol {
			msgs = append(msgs, fmt.Sprintf("%s: max error %.3g at t = %.4g s (tolerance %.3g)", name, worst, at, tol))
		}
	}
	return msgs
}

// CompareTableAtSim is CompareTable with the roles of the grids swapped:
// the reference is interpolated onto the simulation's own time points (only
// those inside the reference's time span). It is for dense references —
// ngspice twins run with reltol=1e-6 time-step control, several points per
// simulation step — where linear interpolation of the reference adds
// almost nothing, so the measured error is the simulator's own at the
// points it computed. Interpolating a coarse simulation onto a dense
// reference would instead measure the straight lines drawn between its
// steps (mutual_tran.cir: 3.7e-3 V between 2 µs steps, against at most
// 1.2e-3 V at the steps, from the exact solution). Same tolerance.
func CompareTableAtSim(simTime []float64, get func(name string) ([]float64, bool), want *Table) []string {
	var msgs []string
	if len(want.Time) == 0 {
		return []string{"empty reference"}
	}
	t0, t1 := want.Time[0], want.Time[len(want.Time)-1]
	for _, name := range want.Names {
		y, ok := get(name)
		if !ok {
			msgs = append(msgs, fmt.Sprintf("%s: not in simulation", name))
			continue
		}
		ref := want.Columns[name]
		full := 0.0
		for _, v := range ref {
			full = math.Max(full, math.Abs(v))
		}
		tol := math.Max(TranFullScale*full, AbsVoltage)
		worst, at := 0.0, 0.0
		for k, ts := range simTime {
			if ts < t0 || ts > t1 {
				continue
			}
			if d := math.Abs(y[k] - output.Interp(want.Time, ref, ts)); d > worst {
				worst, at = d, ts
			}
		}
		if worst > tol {
			msgs = append(msgs, fmt.Sprintf("%s: max error %.3g at t = %.4g s (tolerance %.3g)", name, worst, at, tol))
		}
	}
	return msgs
}

// AC tolerances: the AC solution is exact at each frequency (no time-stepping
// error), so magnitudes must agree to 1e-6 relative and phases to 1e-4°.
const (
	ACRel      = 1e-6
	ACPhaseDeg = 1e-4
)

// CompareAC checks vm(node)/vp(node) columns against simulated phasors.
// get returns the simulated phasor series for a node name.
func CompareAC(freq []float64, get func(node string) ([]complex128, bool), want *Table) []string {
	var msgs []string
	if len(freq) != len(want.Time) {
		return []string{fmt.Sprintf("frequency count %d, want %d", len(freq), len(want.Time))}
	}
	for k := range freq {
		if math.Abs(freq[k]-want.Time[k]) > 1e-9*want.Time[k] {
			return []string{fmt.Sprintf("frequency %d is %g, want %g", k, freq[k], want.Time[k])}
		}
	}
	for _, name := range want.Names {
		if len(name) < 4 || (name[:3] != "vm(" && name[:3] != "vp(") {
			msgs = append(msgs, fmt.Sprintf("%s: unsupported AC column", name))
			continue
		}
		node := strings.TrimSuffix(name[3:], ")")
		v, ok := get(node)
		if !ok {
			msgs = append(msgs, fmt.Sprintf("%s: node not in simulation", name))
			continue
		}
		for k, ref := range want.Columns[name] {
			if name[:3] == "vm(" {
				got := cmplx.Abs(v[k])
				if math.Abs(got-ref) > ACRel*math.Max(math.Abs(ref), 1e-12) {
					msgs = append(msgs, fmt.Sprintf("%s at %g Hz: got %.10g, want %.10g", name, freq[k], got, ref))
					break
				}
			} else {
				got := cmplx.Phase(v[k]) * 180 / math.Pi
				d := math.Mod(math.Abs(got-ref), 360)
				if d > 180 {
					d = 360 - d
				}
				if d > ACPhaseDeg {
					msgs = append(msgs, fmt.Sprintf("%s at %g Hz: got %.8g°, want %.8g°", name, freq[k], got, ref))
					break
				}
			}
		}
	}
	return msgs
}

// CompareSweep checks a DC sweep point by point with the DC operating-point
// tolerances: max(AbsVoltage or AbsCurrent, Rel·|want|). sweep holds the
// simulated sweep values, which must match the table's first column.
func CompareSweep(sweep []float64, get func(name string) ([]float64, bool), want *Table) []string {
	if len(sweep) != len(want.Time) {
		return []string{fmt.Sprintf("sweep has %d points, want %d", len(sweep), len(want.Time))}
	}
	var msgs []string
	for _, name := range want.Names {
		y, ok := get(name)
		if !ok {
			msgs = append(msgs, fmt.Sprintf("%s: not in simulation", name))
			continue
		}
		abs := AbsVoltage
		if strings.HasPrefix(name, "i(") {
			abs = AbsCurrent
		}
		for k, ref := range want.Columns[name] {
			if math.Abs(sweep[k]-want.Time[k]) > 1e-9*math.Max(1, math.Abs(want.Time[k])) {
				return append(msgs, fmt.Sprintf("sweep point %d is %g, want %g", k, sweep[k], want.Time[k]))
			}
			if tol := math.Max(abs, Rel*math.Abs(ref)); math.Abs(y[k]-ref) > tol {
				msgs = append(msgs, fmt.Sprintf("%s at %g: got %.9g, want %.9g (tolerance %.3g)", name, sweep[k], y[k], ref, tol))
				break
			}
		}
	}
	return msgs
}
