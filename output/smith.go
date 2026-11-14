// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package output

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"html"
	"io"
	"math"
	"math/cmplx"
	"sort"
	"strconv"
	"strings"
	"text/tabwriter"
)

// SmithRegion classifies one point of a Smith sweep by the sign of the
// port resistance Re Z, which is the sign of the average power the port
// absorbs. On the chart it says where the point lies: inside the unit
// circle, on it, or outside it.
type SmithRegion int

const (
	// SmithPassive: Re Z > 0, |Γ| < 1. The port absorbs power.
	SmithPassive SmithRegion = iota
	// SmithLossless: Re Z = 0 (to round-off), |Γ| = 1. Purely reactive.
	SmithLossless
	// SmithActive: Re Z < 0, |Γ| > 1. Negative resistance: the port
	// returns more power than it is given, and the energy ledger names
	// the element or channel that pays for it (ModelSupplied, Inflow).
	SmithActive
)

func (r SmithRegion) String() string {
	switch r {
	case SmithPassive:
		return "passive"
	case SmithLossless:
		return "lossless"
	case SmithActive:
		return "active"
	}
	return "unknown"
}

// SmithPoint is the port seen at one frequency.
type SmithPoint struct {
	Freq float64    // [Hz]
	Z    complex128 // driving-point impedance [Ω]
	// Gamma = (Z − Z0)/(Z + Z0). Infinite (cmplx.Inf) when Z = −Z0, the
	// pole of the map where the r = −1 line leaves the chart.
	Gamma  complex128
	Region SmithRegion
	// PortPower is the average power the port absorbs from the 1 A (peak)
	// probe, ½·Re Z [W]. Negative in the active region.
	PortPower float64
	// ModelSupplied and Inflow are the AC ledger's columns at this
	// frequency [W]: power from idealized active elements (a negative
	// resistor) and from declared external inputs. In the active region
	// they are where the returned power comes from.
	ModelSupplied, Inflow float64
}

// SmithResult is a Smith-chart sweep of one port: the circuit's
// independent AC sources set to zero and a 1 A AC current probe driven
// into Plus and out of Minus.
type SmithResult struct {
	Title       string
	Plus, Minus string
	Z0          float64 // reference impedance [Ω]
	Points      []SmithPoint
	// MaxRelativeResidual is the AC ledger's worst residual over the sweep.
	MaxRelativeResidual float64
}

// SmithGamma is the reflection coefficient of Z against a real reference
// Z0: (Z − Z0)/(Z + Z0), infinite at Z = −Z0.
func SmithGamma(z complex128, z0 float64) complex128 {
	d := z + complex(z0, 0)
	if d == 0 {
		return cmplx.Inf()
	}
	return (z - complex(z0, 0)) / d
}

// ClassifySmith returns the region of an impedance; tol is relative to
// max(|Z|, Z0) (1e-9 when tol ≤ 0).
func ClassifySmith(z complex128, z0, tol float64) SmithRegion {
	if tol <= 0 {
		tol = 1e-9
	}
	scale := math.Max(cmplx.Abs(z), z0)
	switch {
	case math.Abs(real(z)) <= tol*scale:
		return SmithLossless
	case real(z) > 0:
		return SmithPassive
	}
	return SmithActive
}

// MaxGamma returns the largest finite |Γ| of the sweep and whether some
// point is at Γ = ∞.
func (r *SmithResult) MaxGamma() (float64, bool) {
	m, inf := 0.0, false
	for _, p := range r.Points {
		a := cmplx.Abs(p.Gamma)
		if math.IsInf(a, 0) || math.IsNaN(a) {
			inf = true
			continue
		}
		m = math.Max(m, a)
	}
	return m, inf
}

// ActiveBands returns the frequency ranges [f_lo, f_hi] over which the
// port is in the active region (|Γ| > 1), as measured at the sweep points.
func (r *SmithResult) ActiveBands() [][2]float64 {
	var out [][2]float64
	in := false
	for _, p := range r.Points {
		if p.Region == SmithActive {
			if !in {
				out = append(out, [2]float64{p.Freq, p.Freq})
				in = true
			}
			out[len(out)-1][1] = p.Freq
		} else {
			in = false
		}
	}
	return out
}

// WriteSmithReport writes a summary of a Smith sweep: the reference, the
// largest |Γ|, the negative-resistance bands with the power the ledger
// books for them, and the first, last and most extreme points.
func WriteSmithReport(w io.Writer, r *SmithResult) error {
	ew := &errWriter{w: w}
	ew.printf("Smith chart: port %s–%s, Z0 = %s\n", r.Plus, r.Minus, Eng(r.Z0, "Ω"))
	n := len(r.Points)
	if n == 0 {
		ew.printf("no frequencies\n")
		return ew.err
	}
	mg, inf := r.MaxGamma()
	ew.printf("Sweep  %d points, %s to %s; max |Γ| = %.4g", n, Eng(r.Points[0].Freq, "Hz"), Eng(r.Points[n-1].Freq, "Hz"), mg)
	if inf {
		ew.printf(" (and Γ = ∞ where Z = −Z0)")
	}
	ew.printf("\n")
	bands := r.ActiveBands()
	if len(bands) == 0 {
		ew.printf("No negative resistance: every point is inside or on the unit circle.\n")
	}
	for _, b := range bands {
		var worst SmithPoint
		for _, p := range r.Points {
			if p.Freq >= b[0] && p.Freq <= b[1] && p.PortPower < worst.PortPower {
				worst = p
			}
		}
		ew.printf("Negative resistance %s to %s: strongest at %s, Re Z = %s, |Γ| = %.4g; port returns %s, ledger: model-supplied %s, inflow %s\n",
			Eng(b[0], "Hz"), Eng(b[1], "Hz"), Eng(worst.Freq, "Hz"), Eng(real(worst.Z), "Ω"), cmplx.Abs(worst.Gamma),
			Eng(-worst.PortPower, "W"), Eng(worst.ModelSupplied, "W"), Eng(worst.Inflow, "W"))
	}
	idx := []int{0, n - 1}
	best := 0
	for i, p := range r.Points {
		if a := cmplx.Abs(p.Gamma); a > cmplx.Abs(r.Points[best].Gamma) || math.IsInf(a, 0) {
			best = i
		}
	}
	idx = append(idx, best)
	sort.Ints(idx)
	tw := tabwriter.NewWriter(ew, 2, 4, 2, ' ', 0)
	fmt.Fprintf(tw, "\n  f\tZ\tΓ\t|Γ|\tregion\n")
	last := -1
	for _, i := range idx {
		if i == last {
			continue
		}
		last = i
		p := r.Points[i]
		fmt.Fprintf(tw, "  %s\t%.5g %+.5gj Ω\t%.4g ∠ %.4g°\t%.4g\t%s\n", Eng(p.Freq, "Hz"), real(p.Z), imag(p.Z),
			cmplx.Abs(p.Gamma), PhaseDeg(p.Gamma), cmplx.Abs(p.Gamma), p.Region)
	}
	tw.Flush()
	ew.printf("\nAverage-power ledger: worst relative residual %.2g over the sweep\n", r.MaxRelativeResidual)
	return ew.err
}

// WriteSmithCSV writes every point of a Smith sweep.
func WriteSmithCSV(w io.Writer, r *SmithResult) error {
	cw := csv.NewWriter(w)
	f := func(v float64) string { return strconv.FormatFloat(v, 'g', 12, 64) }
	if err := cw.Write([]string{"freq", "re_z", "im_z", "re_gamma", "im_gamma", "mag_gamma", "phase_gamma_deg", "region", "port_power", "model_supplied", "inflow"}); err != nil {
		return err
	}
	for _, p := range r.Points {
		g := p.Gamma
		row := []string{f(p.Freq), f(real(p.Z)), f(imag(p.Z)), f(real(g)), f(imag(g)), f(cmplx.Abs(g)), f(PhaseDeg(g)),
			p.Region.String(), f(p.PortPower), f(p.ModelSupplied), f(p.Inflow)}
		if err := cw.Write(row); err != nil {
			return err
		}
	}
	cw.Flush()
	return cw.Error()
}

// JSONFloat is a float64 that is written as null when it is not finite
// (JSON has no infinity), as Γ is at Z = −Z0.
type JSONFloat float64

// MarshalJSON implements json.Marshaler.
func (f JSONFloat) MarshalJSON() ([]byte, error) {
	v := float64(f)
	if math.IsInf(v, 0) || math.IsNaN(v) {
		return []byte("null"), nil
	}
	return json.Marshal(v)
}

// SmithDoc is a Smith sweep in the JSON document.
type SmithDoc struct {
	Plus          string       `json:"plus"`
	Minus         string       `json:"minus"`
	Z0            float64      `json:"z0"`
	Freq          []float64    `json:"freq"`
	ReZ           []float64    `json:"re_z"`
	ImZ           []float64    `json:"im_z"`
	ReGamma       []JSONFloat  `json:"re_gamma"` // null where Γ = ∞
	ImGamma       []JSONFloat  `json:"im_gamma"`
	Region        []string     `json:"region"` // passive, lossless, active
	PortPower     []float64    `json:"port_power"`
	ModelSupplied []float64    `json:"model_supplied"`
	Inflow        []float64    `json:"inflow"`
	MaxGamma      float64      `json:"max_gamma"`
	ActiveBands   [][2]float64 `json:"active_bands"`
}

// SmithDocument converts a Smith sweep for the JSON document.
func SmithDocument(r *SmithResult) *SmithDoc {
	d := &SmithDoc{Plus: r.Plus, Minus: r.Minus, Z0: r.Z0, ActiveBands: r.ActiveBands()}
	if d.ActiveBands == nil {
		d.ActiveBands = [][2]float64{}
	}
	d.MaxGamma, _ = r.MaxGamma()
	for _, p := range r.Points {
		d.Freq = append(d.Freq, p.Freq)
		d.ReZ = append(d.ReZ, real(p.Z))
		d.ImZ = append(d.ImZ, imag(p.Z))
		d.ReGamma = append(d.ReGamma, JSONFloat(real(p.Gamma)))
		d.ImGamma = append(d.ImGamma, JSONFloat(imag(p.Gamma)))
		d.Region = append(d.Region, p.Region.String())
		d.PortPower = append(d.PortPower, p.PortPower)
		d.ModelSupplied = append(d.ModelSupplied, p.ModelSupplied)
		d.Inflow = append(d.Inflow, p.Inflow)
	}
	return d
}

// SmithTrace is one curve drawn on an extended Smith chart.
type SmithTrace struct {
	Name  string
	Gamma []complex128
	Color string // CSS colour; a default palette is used when empty
}

// SmithChartOptions controls WriteSmithSVG.
type SmithChartOptions struct {
	Title string
	// Extent is the half-width of the Γ plane drawn: the chart covers
	// −Extent ≤ Re Γ, Im Γ ≤ Extent. 0 chooses it from the traces: 1.1
	// for passive traces, up to 5 for active ones.
	Extent float64
	// Size is the plot square in pixels (default 640).
	Size float64
}

var smithPalette = []string{"#1f6fd1", "#d1361f", "#2a9d3c", "#8a3ffc", "#c97a00", "#007d8a"}

// WriteSmithSVG draws an extended Smith chart as a standalone SVG.
//
// The ordinary chart is the unit disk |Γ| ≤ 1, where Re z ≥ 0. A negative
// resistance has Re z < 0 and |Γ| > 1, so this chart draws the whole Γ
// plane out to ±Extent: left, right, top and bottom. Constant-r circles
// are drawn for r ≥ 0 (inside the unit circle, grey) and r < 0 (outside
// it, red), with r = −1 the straight line Re Γ = 1; constant-x circles are
// drawn whole, both inside and outside. The region |Γ| > 1 is shaded, and
// trace segments there are dashed, so energy leaving a port is never
// mistaken for energy entering it (the folded 1/Γ* chart would put
// r = −R on top of r = +R).
func WriteSmithSVG(w io.Writer, traces []SmithTrace, opt SmithChartOptions) error {
	ext := opt.Extent
	if ext <= 0 {
		ext = 1.1
		for _, t := range traces {
			for _, g := range t.Gamma {
				a := cmplx.Abs(g)
				if !math.IsInf(a, 0) && !math.IsNaN(a) && a*1.1 > ext {
					ext = a * 1.1
				}
			}
		}
		ext = math.Min(math.Ceil(ext*2)/2, 5)
	}
	size := opt.Size
	if size <= 0 {
		size = 640
	}
	const ml, mt, mr, mb = 64.0, 48.0, 180.0, 56.0
	W, H := ml+size+mr, mt+size+mb
	s := size / (2 * ext)
	cx, cy := ml+size/2, mt+size/2
	X := func(gx float64) float64 { return cx + gx*s }
	Y := func(gy float64) float64 { return cy - gy*s }
	ew := &errWriter{w: w}
	ew.printf(`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.0f %.0f" width="%.0f" height="%.0f" font-family="sans-serif" font-size="11">`+"\n", W, H, W, H)
	ew.printf(`<rect width="100%%" height="100%%" fill="#ffffff"/>` + "\n")
	ew.printf(`<defs><clipPath id="plot"><rect x="%.2f" y="%.2f" width="%.2f" height="%.2f"/></clipPath></defs>`+"\n", ml, mt, size, size)
	if opt.Title != "" {
		ew.printf(`<text x="%.1f" y="%.1f" font-size="15" font-weight="bold">%s</text>`+"\n", ml, mt-20, html.EscapeString(opt.Title))
	}
	ew.printf(`<g clip-path="url(#plot)">` + "\n")
	// Active region |Γ| > 1 shaded: the square with the unit disk cut out.
	ew.printf(`<path fill="#fdecea" fill-rule="evenodd" d="M%.2f %.2fH%.2fV%.2fH%.2fZ M%.2f %.2f a%.2f %.2f 0 1 0 %.2f 0 a%.2f %.2f 0 1 0 %.2f 0Z"/>`+"\n",
		ml, mt, ml+size, mt+size, ml, X(-1), Y(0), s, s, 2*s, s, s, -2*s)
	// Γ-plane grid.
	step := 0.5
	if ext > 2.5 {
		step = 1
	}
	for v := -math.Floor(ext/step) * step; v <= ext+1e-9; v += step {
		ew.printf(`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#e4e4e4" stroke-dasharray="3 3"/>`+"\n", X(v), mt, X(v), mt+size)
		ew.printf(`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#e4e4e4" stroke-dasharray="3 3"/>`+"\n", ml, Y(v), ml+size, Y(v))
	}
	circle := func(c0, r float64, stroke string, width float64) {
		ew.printf(`<circle cx="%.2f" cy="%.2f" r="%.2f" fill="none" stroke="%s" stroke-width="%.2f"/>`+"\n", X(c0), Y(0), r*s, stroke, width)
	}
	// Constant-x circles: centre (1, 1/x), radius 1/|x|.
	for _, x := range []float64{0.2, 0.5, 1, 2, 5} {
		for _, sg := range []float64{1, -1} {
			ew.printf(`<circle cx="%.2f" cy="%.2f" r="%.2f" fill="none" stroke="#9db8e8" stroke-width="0.8"/>`+"\n", X(1), Y(sg/x), s/x)
		}
	}
	// Im z = 0: the real axis.
	ew.printf(`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#9db8e8" stroke-width="0.8"/>`+"\n", ml, Y(0), ml+size, Y(0))
	// Constant-r circles: centre r/(1+r), radius 1/|1+r|; r = −1 is the line Re Γ = 1.
	for _, r := range []float64{0.2, 0.5, 1, 2, 5} {
		circle(r/(1+r), 1/(1+r), "#b8b8b8", 0.8)
	}
	for _, r := range []float64{-0.2, -0.5, -1.5, -2, -3, -5} {
		circle(r/(1+r), 1/math.Abs(1+r), "#e39a90", 0.9)
	}
	ew.printf(`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#c0392b" stroke-width="1.2"/>`+"\n", X(1), mt, X(1), mt+size)
	circle(0, 1, "#222222", 1.6) // r = 0, |Γ| = 1
	// Traces: solid inside the unit circle, dashed outside; broken where Γ
	// runs off to infinity (Z through −Z0).
	far := 3 * ext
	for ti, t := range traces {
		col := t.Color
		if col == "" {
			col = smithPalette[ti%len(smithPalette)]
		}
		var run []complex128
		runOut := false
		flush := func() {
			if len(run) > 1 {
				dash := ""
				if runOut {
					dash = ` stroke-dasharray="6 3"`
				}
				var b strings.Builder
				for i, g := range run {
					if i == 0 {
						fmt.Fprintf(&b, "M%.2f %.2f", X(real(g)), Y(imag(g)))
					} else {
						fmt.Fprintf(&b, "L%.2f %.2f", X(real(g)), Y(imag(g)))
					}
				}
				ew.printf(`<path d="%s" fill="none" stroke="%s" stroke-width="2"%s/>`+"\n", b.String(), col, dash)
			}
			run = run[:0]
		}
		for _, g := range t.Gamma {
			a := cmplx.Abs(g)
			if math.IsInf(a, 0) || math.IsNaN(a) || a > far {
				flush()
				continue
			}
			out := a > 1+1e-12
			if len(run) > 0 && out != runOut {
				prev := run[len(run)-1]
				flush()
				run = append(run, prev)
			}
			runOut = out
			run = append(run, g)
		}
		flush()
		if len(t.Gamma) > 0 {
			g := t.Gamma[0]
			if a := cmplx.Abs(g); !math.IsInf(a, 0) && a <= far {
				ew.printf(`<circle cx="%.2f" cy="%.2f" r="3.5" fill="%s"/>`+"\n", X(real(g)), Y(imag(g)), col)
			}
		}
	}
	ew.printf("</g>\n")
	// Labels of the r circles where they cross the real axis, Γ = (r−1)/(r+1).
	for _, r := range []float64{0, 0.2, 0.5, 1, 2, 5, -0.2, -0.5, -1.5, -2, -3, -5} {
		g := (r - 1) / (r + 1)
		if math.Abs(g) > ext || (r > 0 && s < 120) { // passive labels only where there is room
			continue
		}
		col := "#666666"
		if r < 0 {
			col = "#c0392b"
		}
		ew.printf(`<text x="%.2f" y="%.2f" fill="%s" font-size="9">r=%s</text>`+"\n", X(g)+2, Y(0)-3, col, strconv.FormatFloat(r, 'g', -1, 64))
	}
	if 1 <= ext {
		ew.printf(`<text x="%.2f" y="%.2f" fill="#c0392b" font-size="9">r=−1</text>`+"\n", X(1)+2, mt+12)
	}
	// Axes.
	ew.printf(`<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="none" stroke="#444"/>`+"\n", ml, mt, size, size)
	for v := -math.Floor(ext/step) * step; v <= ext+1e-9; v += step {
		lab := strconv.FormatFloat(v, 'g', 3, 64)
		ew.printf(`<text x="%.2f" y="%.2f" text-anchor="middle">%s</text>`+"\n", X(v), mt+size+16, lab)
		ew.printf(`<text x="%.2f" y="%.2f" text-anchor="end">%s</text>`+"\n", ml-6, Y(v)+4, lab)
	}
	ew.printf(`<text x="%.2f" y="%.2f" text-anchor="middle" font-size="12">Re Γ</text>`+"\n", cx, mt+size+38)
	ew.printf(`<text x="%.2f" y="%.2f" text-anchor="middle" font-size="12" transform="rotate(-90 %.2f %.2f)">Im Γ</text>`+"\n", ml-40, cy, ml-40, cy)
	// Legend.
	lx, ly := ml+size+16, mt+10
	key := func(col, dash, label string) {
		ew.printf(`<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="2"%s/>`, lx, ly, lx+22, ly, col, dash)
		ew.printf(`<text x="%.2f" y="%.2f">%s</text>`+"\n", lx+28, ly+4, html.EscapeString(label))
		ly += 18
	}
	for ti, t := range traces {
		col := t.Color
		if col == "" {
			col = smithPalette[ti%len(smithPalette)]
		}
		key(col, "", t.Name)
	}
	ly += 8
	key("#222222", "", "|Γ| = 1 (r = 0)")
	key("#b8b8b8", "", "r > 0: passive")
	key("#e39a90", "", "r < 0: negative R")
	key("#9db8e8", "", "constant x")
	key("#555555", ` stroke-dasharray="6 3"`, "trace with |Γ| > 1")
	ew.printf(`<text x="%.2f" y="%.2f" fill="#888" font-size="9">shaded: |Γ| &gt; 1, energy returned</text>`+"\n", lx, ly+6)
	ew.printf("</svg>\n")
	return ew.err
}
