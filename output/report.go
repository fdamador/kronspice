// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"fmt"
	"io"
	"strings"
	"text/tabwriter"

	"github.com/fdamador/kronspice/exchange"
)

// ReportOptions controls WriteReport.
type ReportOptions struct {
	// AllChannels lists every channel of every element, including the ones
	// that are not modeled. Otherwise unmodeled channels are summarized.
	AllChannels bool
}

// WriteReport writes a human-readable operating-point report: node voltages,
// element operating points and the element × channel energy ledger.
func WriteReport(w io.Writer, r *Result, opt ReportOptions) error {
	ew := &errWriter{w: w}
	ew.printf("%s\n", r.Title)
	ew.printf("Analysis: %s\n\n", r.Analysis)

	ew.printf("Node voltages\n")
	tw := tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	for i, n := range r.NodeNames {
		if i == 0 {
			continue
		}
		fmt.Fprintf(tw, "  v(%s)\t%s\n", n, Eng(r.NodeVoltages[i], "V"))
	}
	tw.Flush()

	ew.printf("\nElements\n")
	tw = tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	fmt.Fprintf(tw, "  name\tkind\tnodes\tvoltage\tcurrent\tpower absorbed\n")
	for _, e := range r.Elements {
		fmt.Fprintf(tw, "  %s\t%s\t%s→%s\t%s\t%s\t%s\n", e.Name, e.Kind, e.Nodes[0], e.Nodes[1],
			Eng(e.Voltage, "V"), Eng(e.Current, "A"), Eng(e.PowerAbsorbed, "W"))
	}
	tw.Flush()

	if l := r.Ledger; l != nil {
		ew.printf("\nEnergy ledger (+ into the circuit from the exterior, − out)\n")
		tw = tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
		fmt.Fprintf(tw, "  element\tchannel\tstatus\tpower\tnote\n")
		unmodeled := map[string]int{}
		for _, row := range l.Rows {
			if row.Status == exchange.NotModeled && !opt.AllChannels {
				unmodeled[row.Element]++
				continue
			}
			fmt.Fprintf(tw, "  %s\t%s\t%s\t%s\t%s\n", row.Element, row.Channel, row.Status, Eng(row.PowerW, "W"), row.Note)
		}
		tw.Flush()
		if !opt.AllChannels && len(unmodeled) > 0 {
			ew.printf("  (every other channel of each element is not modeled: 0 W each; use -all-channels to list them)\n")
		}
		ew.printf("\n  inflow        %s\n", Eng(l.Inflow, "W"))
		ew.printf("  outflow       %s\n", Eng(l.Outflow, "W"))
		ew.printf("  dissipation   %s\n", Eng(l.Dissipation, "W"))
		ew.printf("  storage rate  %s\n", Eng(l.StorageRate, "W"))
		if l.ModelSupplied != 0 {
			ew.printf("  model-supplied %s (idealized active elements; not inflow from the exterior)\n", Eng(l.ModelSupplied, "W"))
			for _, row := range l.Supplied {
				ew.printf("    %s  %s\n", row.Element, Eng(row.PowerW, "W"))
			}
		}
		if l.Hypothetical != 0 {
			ew.printf("  hypothetical (not counted) %s\n", Eng(l.Hypothetical, "W"))
		}
		ew.printf("  residual      %s (relative %.2g)\n", Eng(l.Residual, "W"), l.RelativeResidual())
	}
	return ew.err
}

// Eng formats x with an engineering prefix, e.g. Eng(0.0047, "A") = "4.7 mA".
func Eng(x float64, unit string) string {
	if x == 0 {
		return "0 " + unit
	}
	prefixes := []struct {
		s string
		f float64
	}{{"T", 1e12}, {"G", 1e9}, {"M", 1e6}, {"k", 1e3}, {"", 1}, {"m", 1e-3}, {"µ", 1e-6}, {"n", 1e-9}, {"p", 1e-12}, {"f", 1e-15}}
	a := x
	if a < 0 {
		a = -a
	}
	for _, p := range prefixes {
		if a >= p.f*0.9999995 {
			s := fmt.Sprintf("%.6g", x/p.f)
			return strings.TrimSpace(s + " " + p.s + unit)
		}
	}
	return fmt.Sprintf("%.6g %s", x, unit)
}

type errWriter struct {
	w   io.Writer
	err error
}

func (e *errWriter) Write(p []byte) (int, error) {
	if e.err != nil {
		return 0, e.err
	}
	n, err := e.w.Write(p)
	e.err = err
	return n, err
}

func (e *errWriter) printf(format string, a ...any) { fmt.Fprintf(e, format, a...) }

// WriteTranReport writes a summary of a transient analysis: time-stepping
// statistics, final node voltages and the integrated energy ledger.
func WriteTranReport(w io.Writer, wf *Waveforms) error {
	ew := &errWriter{w: w}
	ew.printf("%s\nAnalysis: tran\n\n", wf.Title)
	n := len(wf.Time)
	if n == 0 {
		ew.printf("no time points recorded\n")
		return ew.err
	}
	s := wf.Stats
	ew.printf("Time points  %d recorded, %d steps accepted, %d rejected\n", n, s.Accepted, s.Rejected)
	ew.printf("Step size    %s to %s\n", Eng(s.MinStep, "s"), Eng(s.MaxStep, "s"))
	ew.printf("Interval     %s to %s\n\n", Eng(wf.Time[0], "s"), Eng(wf.Time[n-1], "s"))

	ew.printf("Node voltages at the end\n")
	tw := tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	for i, name := range wf.NodeNames {
		if i == 0 {
			continue
		}
		fmt.Fprintf(tw, "  v(%s)\t%s\n", name, Eng(wf.V[i][n-1], "V"))
	}
	tw.Flush()

	l := wf.Ledger
	ew.printf("\nEnergy ledger over the interval\n")
	ew.printf("  energy in            %s\n", Eng(l.EnergyIn, "J"))
	ew.printf("  energy out           %s (returned to sources and ports)\n", Eng(l.EnergyOut, "J"))
	ew.printf("  energy dissipated    %s\n", Eng(l.EnergyDissipated, "J"))
	if l.EnergyModelSupplied != 0 {
		ew.printf("  model-supplied       %s (idealized active elements; not inflow from the exterior)\n", Eng(l.EnergyModelSupplied, "J"))
	}
	ew.printf("  stored: start → end  %s → %s\n", Eng(l.StoredStart, "J"), Eng(l.StoredEnd, "J"))
	ew.printf("  integration mismatch %s (time-stepping error, not a physical channel)\n", Eng(l.IntegrationMismatch, "J"))
	ew.printf("  worst point residual %.2g (relative)\n", l.MaxRelativeResidual)
	return ew.err
}

// WriteACReport writes a summary of an AC analysis: the sweep, each node's
// magnitude and phase at the first and last frequency, and the ledger check.
func WriteACReport(w io.Writer, r *ACResult) error {
	ew := &errWriter{w: w}
	ew.printf("%s\nAnalysis: ac\n\n", r.Title)
	n := len(r.Freq)
	if n == 0 {
		ew.printf("no frequencies\n")
		return ew.err
	}
	ew.printf("Sweep  %d points, %s to %s\n\n", n, Eng(r.Freq[0], "Hz"), Eng(r.Freq[n-1], "Hz"))
	tw := tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	fmt.Fprintf(tw, "  node\tat %s\tat %s\n", Eng(r.Freq[0], "Hz"), Eng(r.Freq[n-1], "Hz"))
	for i, name := range r.NodeNames {
		if i == 0 {
			continue
		}
		a, b := r.V[i][0], r.V[i][n-1]
		fmt.Fprintf(tw, "  v(%s)\t%.4g dB ∠ %.4g°\t%.4g dB ∠ %.4g°\n", name, DB(a), PhaseDeg(a), DB(b), PhaseDeg(b))
	}
	tw.Flush()
	ew.printf("\nAverage-power ledger: worst relative residual %.2g over the sweep\n", r.Ledger.MaxRelativeResidual)
	return ew.err
}

// WriteDCReport writes a summary of a DC sweep.
func WriteDCReport(w io.Writer, r *DCSweepResult) error {
	ew := &errWriter{w: w}
	ew.printf("%s\nAnalysis: dc\n\n", r.Title)
	if len(r.Values) == 0 || len(r.Values[0]) == 0 {
		ew.printf("no points\n")
		return ew.err
	}
	n := len(r.Values[0])
	ew.printf("Swept %v, %d operating points\n\n", r.Sources, n)
	tw := tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	fmt.Fprintf(tw, "  node\tfirst point\tlast point\n")
	for i, name := range r.NodeNames {
		if i == 0 {
			continue
		}
		fmt.Fprintf(tw, "  v(%s)\t%s\t%s\n", name, Eng(r.V[i][0], "V"), Eng(r.V[i][n-1], "V"))
	}
	tw.Flush()
	ew.printf("\nLedger: worst relative residual %.2g over the sweep\n", r.MaxRelativeResidual)
	return ew.err
}
