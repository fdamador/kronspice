// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"bufio"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"testing"

	"github.com/fdamador/kronspice/analysis"
	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/netlist"
)

// The ledger closes on every twin case, the magnets supply no net power,
// and the load never receives more than the drive delivers.
func TestLedger(t *testing.T) {
	for _, s := range Twins() {
		r, _, err := s.Run()
		if err != nil {
			t.Fatalf("%s: %v", s.Name, err)
		}
		if math.Abs(r.RelResidual(s)) > 1e-3 {
			t.Errorf("%s: ledger residual %.2g of the drive energy", s.Name, r.RelResidual(s))
		}
		if math.Abs(r.Pmagnet) > 1e-3*r.Pdrive {
			t.Errorf("%s: magnets %.3g W against drive %.3g W", s.Name, r.Pmagnet, r.Pdrive)
		}
		if !(r.COP > 0 && r.COP < 1) {
			t.Errorf("%s: COP %.4f", s.Name, r.COP)
		}
	}
}

// Halving the time step four times leaves the energy split unchanged, and
// the magnets' small net power shrinks with the step (it is integration
// error, not a source).
func TestTimestep(t *testing.T) {
	for _, k := range []Kind{Knife, MEG} {
		var cops, mags []float64
		for _, h := range []float64{0.2e-6, 0.1e-6, 0.05e-6} {
			r, _, err := Spec{Kind: k, Step: h}.Run()
			if err != nil {
				t.Fatal(err)
			}
			cops = append(cops, r.COP)
			mags = append(mags, math.Abs(r.Pmagnet/r.Pdrive))
		}
		for i := 1; i < len(cops); i++ {
			if math.Abs(cops[i]-cops[0]) > 1e-3*cops[0] {
				t.Errorf("%s: COP %v moves with the step", k, cops)
			}
		}
		if mags[len(mags)-1] > 1e-4 {
			t.Errorf("%s: magnet share %v at the smallest step", k, mags)
		}
	}
}

// The ngspice twins model each gyrator as two VCCS. Run that form in
// kronspice too: it must match the native gyrator, so a difference against
// ngspice is the simulators', not the translation's.
func TestGyratorAsVCCS(t *testing.T) {
	re := regexp.MustCompile(`(?m)^(Y\S+) (\S+) (\S+) (\S+) (\S+) (\S+)$`)
	for _, s := range []Spec{{Kind: Knife}, {Kind: MEG}} {
		native, _, err := s.Run()
		if err != nil {
			t.Fatal(err)
		}
		text := re.ReplaceAllStringFunc(s.Kron(), func(l string) string {
			m := re.FindStringSubmatch(l)
			g, _ := strconv.ParseFloat(m[6], 64)
			return "G" + m[1] + "_1 " + m[2] + " " + m[3] + " " + m[4] + " " + m[5] + " " + m[6] + "\n" +
				"G" + m[1] + "_2 " + m[4] + " " + m[5] + " " + m[2] + " " + m[3] + " " + strconv.FormatFloat(-g, 'g', 10, 64)
		})
		deck, err := netlist.ParseString(text)
		if err != nil {
			t.Fatal(err)
		}
		c, err := circuit.FromDeck(deck)
		if err != nil {
			t.Fatal(err)
		}
		a := deck.Analyses[0]
		wf, err := analysis.Tran(c, analysis.TranSpec{Step: a.Step, Stop: a.Stop}, analysis.TranOptions{})
		if err != nil {
			t.Fatal(err)
		}
		_, r, _ := s.elements()
		col := map[string][]float64{}
		for _, name := range s.columns(r) {
			inner := name[2 : len(name)-1]
			if strings.HasPrefix(name, "v(") {
				col[name], _ = wf.NodeVoltage(inner)
			} else {
				col[name], _ = wf.ElementCurrent(strings.ToUpper(inner))
			}
		}
		vccs, err := s.FromColumns(wf.Time, col)
		if err != nil {
			t.Fatal(err)
		}
		if d := math.Abs(vccs.COP-native.COP) / native.COP; d > 1e-4 {
			t.Errorf("%s: COP %.6f with VCCS pairs, %.6f with Y", s.Kind, vccs.COP, native.COP)
		}
	}
}

// netlists/ and ngspice/ are what -write generates.
func TestNetlistsCurrent(t *testing.T) {
	for _, s := range Twins() {
		for dir, want := range map[string]string{"netlists": s.Kron(), "ngspice": s.Ngspice(s.Name)} {
			got, err := os.ReadFile(filepath.Join(dir, s.Name+".cir"))
			if err != nil {
				t.Fatalf("%v (run: go run ./examples/bearden_meg -write)", err)
			}
			if strings.ReplaceAll(string(got), "\r\n", "\n") != want {
				t.Errorf("%s/%s.cir is stale (run: go run ./examples/bearden_meg -write)", dir, s.Name)
			}
		}
	}
}

// Alignment with ngspice: each twin's output (ngspice/NAME.dat, written by
// ngspice_run.sh in the tools container) gives the same drive power, load
// power and COP within 2 %, and the same probe waveforms within 2 % RMS of
// their swing. Skipped for a case until its .dat exists.
func TestNgspiceTwin(t *testing.T) {
	for _, s := range Twins() {
		path := filepath.Join("ngspice", s.Name+".dat")
		if _, err := os.Stat(path); err != nil {
			t.Logf("%s: no %s yet (run ngspice_run.sh in the tools container)", s.Name, path)
			continue
		}
		_, r, _ := s.elements()
		names := s.columns(r)
		tm, col, err := readWrdata(path, names)
		if err != nil {
			t.Fatalf("%s: %v", s.Name, err)
		}
		ng, err := s.FromColumns(tm, col)
		if err != nil {
			t.Fatalf("%s: %v", s.Name, err)
		}
		kr, wf, err := s.Run()
		if err != nil {
			t.Fatal(err)
		}
		for _, q := range []struct {
			what   string
			ng, kr float64
		}{{"drive power", ng.Pdrive, kr.Pdrive}, {"load power", ng.Pload, kr.Pload}, {"COP", ng.COP, kr.COP}} {
			d := math.Abs(q.ng-q.kr) / math.Abs(q.kr)
			t.Logf("%s %s: ngspice %.6g, kronspice %.6g (%.2g)", s.Name, q.what, q.ng, q.kr, d)
			if d > 0.02 {
				t.Errorf("%s %s: ngspice %.6g, kronspice %.6g", s.Name, q.what, q.ng, q.kr)
			}
		}
		t0, t1 := s.Window()
		for _, p := range r.probes {
			kv, _ := wf.NodeVoltage(p)
			nv := col["v("+p+")"]
			var se, lo, hi float64
			n := 0
			lo, hi = math.Inf(1), math.Inf(-1)
			for k, x := range tm {
				if x < t0 || x > t1 {
					continue
				}
				d := nv[k] - at(wf.Time, kv, x)
				se += d * d
				n++
				lo, hi = math.Min(lo, nv[k]), math.Max(hi, nv[k])
			}
			swing := hi - lo
			if n == 0 || swing == 0 {
				continue
			}
			rms := math.Sqrt(se/float64(n)) / swing
			t.Logf("%s v(%s): RMS difference %.2g of the swing", s.Name, p, rms)
			if rms > 0.02 {
				t.Errorf("%s v(%s): RMS difference %.3g of the swing", s.Name, p, rms)
			}
		}
	}
}

// readWrdata reads ngspice wrdata output: per row, a time and value pair
// for each vector, in the order written.
func readWrdata(path string, names []string) ([]float64, map[string][]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, nil, err
	}
	defer f.Close()
	col := map[string][]float64{}
	var tm []float64
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 1<<20), 1<<24)
	for sc.Scan() {
		fs := strings.Fields(sc.Text())
		if len(fs) != 2*len(names) {
			continue
		}
		x, err := strconv.ParseFloat(fs[0], 64)
		if err != nil {
			continue
		}
		tm = append(tm, x)
		for i, n := range names {
			v, err := strconv.ParseFloat(fs[2*i+1], 64)
			if err != nil {
				return nil, nil, err
			}
			col[n] = append(col[n], v)
		}
	}
	return tm, col, sc.Err()
}
