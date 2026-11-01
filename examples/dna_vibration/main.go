// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

// Command dna_vibration models the four DNA bases and the two Watson–Crick
// base pairs as electric circuits and finds, from kronspice's AC analysis,
// the infrared light each one absorbs (incoming) and releases (outgoing).
// It is the chemistry branch's next step after examples/water_vibration
// (doc/chemistry-plan.md).
//
// Molecules: 9-methyladenine, 1-methylthymine, 9-methylguanine,
// 1-methylcytosine (the methyl stands in for the sugar, as in gas-phase
// spectroscopy) and the pairs A·T and G·C, cut from an all-atom B-DNA
// duplex and minimized in the AMBER DNA force field (OL15; parm99/bsc0
// base terms; prep/prep.py writes molecules.json from OpenMM's
// amber14/DNA.OL15.xml and its bdna.xyz test structure).
//
// Circuit (circuit.go): atoms are capacitors (mass), bonds and angles are
// inductor springs between their atoms' velocity nodes (Wilson's B matrix
// as VCVS/CCCS transformers, as for water), the rest of the force field
// (torsions, out-of-plane, electrostatics) is a coupling network, and in a
// pair every force between the two bases — the hydrogen bonds — is its own
// named group G_HB. Light enters as a current source q_a at every atom (the
// force of the light's electric field on the atom's charge).
//
// Outputs: absorption spectrum (power the field puts in, from the AC
// sweep), release (emission) per mode by two routes that must agree —
// the Einstein A coefficient from the IR intensity, and Maxwell's
// radiated power of the vibrating dipole (em.LarmorPower, itself checked
// against the Poynting flux of the full dipole fields in GA and vector
// form) — and the release spectra for one quantum per mode and for
// thermal equilibrium at 300 K (Kirchhoff). Pairing: the hydrogen bonds,
// the modes that live in them, and how each base's lines move when it
// pairs.
//
//	go run ./examples/dna_vibration [-out examples/dna_vibration]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

var order = []string{"A", "T", "G", "C", "AT", "GC"}

func main() {
	log.SetFlags(0)
	out := flag.String("out", "examples/dna_vibration", "directory for netlists/, spectra/, modes/")
	doStacks := flag.Bool("stacks", true, "also analyse 2-, 3- and 4-pair stacks (stacks/)")
	onlyStacks := flag.Bool("only-stacks", false, "analyse the stacks only")
	stackKeys := flag.String("stack-keys", "", "comma-separated stacks to analyse (default: all)")
	flag.Parse()
	ms, err := Molecules()
	if err != nil {
		log.Fatal(err)
	}
	for _, d := range []string{"netlists", "spectra", "modes", "summary"} {
		if err := os.MkdirAll(filepath.Join(*out, d), 0o755); err != nil {
			log.Fatal(err)
		}
	}
	if *onlyStacks {
		if err := RunStacks(*out, ms, splitKeys(*stackKeys)); err != nil {
			log.Fatal(err)
		}
		return
	}
	states := map[string]*State{}
	var summary strings.Builder
	var geoAtoms, geoBonds, geoModes []map[string]any
	summary.WriteString("molecule,name,atoms,vibrations,springs,couplings_inside_bases,pairing_couplings,size_A,d_over_lambda,thermal_release_W,lowest_cm-1,highest_cm-1\n")
	fmt.Println("DNA bases and Watson–Crick pairs as circuits — infrared absorption and release")
	fmt.Println("Force field: AMBER OL15 (parm99/bsc0 base terms), methyl caps at the glycosidic N; harmonic, gas phase.")
	fmt.Println()
	for _, key := range order {
		s, err := Prepare(ms[key])
		if err != nil {
			log.Fatal(err)
		}
		states[key] = s
		springs := s.Springs()
		res, rst := s.ResidualEntries(springs)
		conv := 2 * math.Pi * wavenumber
		src, st := s.Netlist(NetlistOptions{FWHM: FWHMcm, Field: "x",
			Analysis: fmt.Sprintf(".ac lin %d %.9g %.9g", SweepPts, math.Max(SweepFrom, 1)/conv, SweepTo/conv)})
		path := filepath.Join(*out, "netlists", strings.ToLower(key)+".cir")
		if err := os.WriteFile(path, []byte(src), 0o644); err != nil {
			log.Fatal(err)
		}
		nu, ab, err := s.Absorption(SweepFrom, SweepTo, SweepPts, FWHMcm)
		if err != nil {
			log.Fatal(err)
		}
		qe, te := Emission(nu, ab)
		absN := append([]float64(nil), ab...)
		normalize(absN)
		var csv strings.Builder
		csv.WriteString("wavenumber_cm-1,absorption,release_one_quantum,release_thermal_300K,absorbed_power\n")
		for k := range nu {
			fmt.Fprintf(&csv, "%.1f,%.6g,%.6g,%.6g,%.6g\n", nu[k], absN[k], qe[k], te[k], ab[k])
		}
		if err := os.WriteFile(filepath.Join(*out, "spectra", strings.ToLower(key)+".csv"), []byte(csv.String()), 0o644); err != nil {
			log.Fatal(err)
		}
		var mcsv strings.Builder
		mcsv.WriteString("mode,wavenumber_cm-1,ir_intensity_km_mol,release_rate_A_s-1,lifetime_ms,thermal_release_W,pairing_energy_share,assignment,pair_motion\n")
		lam := 1e8 / s.Modes[len(s.Modes)-1].Wavenumber // shortest wavelength, Å
		fmt.Printf("%s (%s): %d atoms, %d vibrations %.0f–%.0f cm⁻¹, RMS gradient %.1e kcal/mol/Å\n",
			key, s.Mol.Title, len(s.Mol.Atoms), len(s.Modes), s.Modes[0].Wavenumber, s.Modes[len(s.Modes)-1].Wavenumber, s.RMS)
		fmt.Printf("  circuit %s: %d bond/angle springs, %d coupling sources inside the bases, %d base-pairing sources, %d entries < 1e-9 dropped (largest %.1e mdyn/Å)\n",
			path, st.Springs, rst.Intra, rst.Inter, rst.Dropped, rst.DroppedMax)
		fmt.Printf("  lumped check: size %.1f Å against the shortest infrared wavelength %.0f Å: d/λ = %.1e (λ/10 rule: lumped)\n", s.Size(), lam, s.Size()/lam)
		var strong []string
		total := 0.0
		for k, md := range s.Modes {
			label, _ := s.Label(md, springs, res)
			r := ReleaseOf(md)
			total += r.Thermal
			motion := ""
			if md.Inter >= 0.25 {
				motion = s.PairMotion(md)
			}
			fmt.Fprintf(&mcsv, "%d,%.1f,%.2f,%.4g,%.4g,%.4g,%.3f,%q,%q\n", k+1, md.Wavenumber, md.Intensity, r.A, 1e3*r.Lifetime, r.Thermal, md.Inter, label, motion)
			if md.Intensity >= 100 {
				strong = append(strong, fmt.Sprintf("%.0f (%s; %.0f km/mol, A = %.0f s⁻¹)", md.Wavenumber, label, md.Intensity, r.A))
			}
		}
		if err := os.WriteFile(filepath.Join(*out, "modes", strings.ToLower(key)+".csv"), []byte(mcsv.String()), 0o644); err != nil {
			log.Fatal(err)
		}
		fmt.Printf("  strongest absorption = strongest release lines (cm⁻¹): %s\n", strings.Join(strong, "; "))
		fmt.Printf("  thermal infrared release at %.0f K: %.2e W per molecule\n\n", TempK, total)
		for i, a := range s.Mol.Atoms {
			geoAtoms = append(geoAtoms, map[string]any{"molecule": key, "index": i, "name": a.Name, "element": a.Element,
				"x": round4(s.X[3*i]), "y": round4(s.X[3*i+1]), "z": round4(s.X[3*i+2])})
		}
		for _, b := range s.Mol.Bonds {
			geoBonds = append(geoBonds, map[string]any{"molecule": key, "i": b.I, "j": b.J, "kind": "covalent"})
		}
		idx := map[string]int{}
		for i, a := range s.Mol.Atoms {
			idx[a.Name] = i
		}
		for _, h := range s.HBonds() {
			geoBonds = append(geoBonds, map[string]any{"molecule": key, "i": idx[h.H], "j": idx[h.Acceptor], "kind": "hbond"})
		}
		for k, md := range s.Modes {
			mx := 0.0
			for _, v := range md.Vector {
				mx = math.Max(mx, math.Abs(v))
			}
			vec := make([]float64, len(md.Vector))
			for i, v := range md.Vector {
				vec[i] = round4(v / mx)
			}
			geoModes = append(geoModes, map[string]any{"molecule": key, "mode": k + 1, "wavenumber_cm-1": round4(md.Wavenumber),
				"pairing_energy_share": round4(md.Inter), "ir_intensity_km_mol": round4(md.Intensity), "vector": vec})
		}
		fmt.Fprintf(&summary, "%s,%q,%d,%d,%d,%d,%d,%.2f,%.3g,%.4g,%.1f,%.1f\n", key, s.Mol.Title, len(s.Mol.Atoms), len(s.Modes),
			st.Springs, rst.Intra, rst.Inter, s.Size(), s.Size()/lam, total, s.Modes[0].Wavenumber, s.Modes[len(s.Modes)-1].Wavenumber)
	}
	if err := os.WriteFile(filepath.Join(*out, "summary", "circuits.csv"), []byte(summary.String()), 0o644); err != nil {
		log.Fatal(err)
	}
	for name, rows := range map[string][]map[string]any{"atoms.json": geoAtoms, "bonds.json": geoBonds, "mode_shapes.json": geoModes} {
		b, err := json.Marshal(rows)
		if err != nil {
			log.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(*out, "summary", name), b, 0o644); err != nil {
			log.Fatal(err)
		}
	}
	var bands strings.Builder
	bands.WriteString("molecule,band,model_cm-1,energy_share\n")
	for _, b := range Bands(states) {
		fmt.Fprintf(&bands, "%s,%s,%.1f,%.2f\n", b.Molecule, b.Name, b.Wavenumber, b.Share)
	}
	if err := os.WriteFile(filepath.Join(*out, "summary", "bands.csv"), []byte(bands.String()), 0o644); err != nil {
		log.Fatal(err)
	}
	hb := pairing(states)
	if err := os.WriteFile(filepath.Join(*out, "summary", "hbonds.csv"), []byte(hb), 0o644); err != nil {
		log.Fatal(err)
	}
	if *doStacks {
		if err := RunStacks(*out, ms, splitKeys(*stackKeys)); err != nil {
			log.Fatal(err)
		}
	}
}

// pairing reports the hydrogen bonds, the inter-base modes and the shifts
// of the hydrogen-bonding groups when the bases pair.
func pairing(states map[string]*State) string {
	var csv strings.Builder
	csv.WriteString("pair,donor,hydrogen,acceptor,h_acceptor_A,donor_acceptor_A,angle_deg\n")
	fmt.Println("Base pairing")
	for _, p := range []struct{ pair, a, b string }{{"AT", "A", "T"}, {"GC", "G", "C"}} {
		s := states[p.pair]
		fmt.Printf("%s:\n", s.Mol.Title)
		for _, h := range s.HBonds() {
			fmt.Printf("  H-bond %s–%s···%s: H···A %.2f Å, D···A %.2f Å, angle %.0f°\n", h.Donor, h.H, h.Acceptor, h.HA, h.DA, h.Angle)
			fmt.Fprintf(&csv, "%s,%s,%s,%s,%.3f,%.3f,%.1f\n", p.pair, h.Donor, h.H, h.Acceptor, h.HA, h.DA, h.Angle)
		}
		springs := s.Springs()
		res, _ := s.ResidualEntries(springs)
		var inter []string
		for _, md := range s.Modes {
			if md.Inter >= 0.25 {
				label, _ := s.Label(md, springs, res)
				inter = append(inter, fmt.Sprintf("%.0f (%.0f%% in the H-bonds; %s)", md.Wavenumber, 100*md.Inter, label))
			}
		}
		fmt.Printf("  modes living in the pairing (≥ 25%% of their energy in the inter-base forces): %s\n", strings.Join(inter, "; "))
		// Shifts of the donor and acceptor groups: each monomer mode dominated
		// by a stretch of an H-bonding bond, against the pair mode dominated
		// by the same bond.
		hb := map[string]bool{}
		for _, h := range s.HBonds() {
			hb[short(h.Donor)+"_"+short(h.H)] = true
			hb[short(h.H)+"_"+short(h.Donor)] = true
		}
		for _, sp := range springs {
			if strings.HasPrefix(sp.Name, "b_") {
				f := strings.Split(strings.TrimPrefix(sp.Name, "b_"), "_")
				if len(f) == 4 {
					for _, h := range s.HBonds() {
						acc := short(h.Acceptor)
						if f[0]+"_"+f[1] == acc || f[2]+"_"+f[3] == acc {
							hb[f[0]+"_"+f[1]+"_"+f[2]+"_"+f[3]] = true
						}
					}
				}
			}
		}
		fmt.Printf("  shifts on pairing (model, harmonic):\n")
		for _, mono := range []string{p.a, p.b} {
			m := states[mono]
			for _, sp := range m.Springs() {
				name := strings.TrimPrefix(sp.Name, "b_")
				if !strings.HasPrefix(sp.Name, "b_") || !(hb[name] || hb[revBond(name)]) {
					continue
				}
				fm := dominant(m, sp.Name)
				fp := dominant(s, sp.Name)
				if fm > 0 && fp > 0 {
					fmt.Printf("    %-18s %7.1f → %7.1f cm⁻¹ (%+.1f)\n", pretty(sp.Name), fm, fp, fp-fm)
				}
			}
		}
		fmt.Println()
	}
	return csv.String()
}

func revBond(n string) string {
	f := strings.Split(n, "_")
	if len(f) != 4 {
		return n
	}
	return f[2] + "_" + f[3] + "_" + f[0] + "_" + f[1]
}

// dominant returns the frequency of the mode with the largest potential
// energy share in the named spring (0 if none has 30 %).
func dominant(s *State, spring string) float64 {
	var sp *Spring
	for _, x := range s.Springs() {
		if x.Name == spring {
			x := x
			sp = &x
		}
	}
	if sp == nil {
		return 0
	}
	type c struct{ f, sh float64 }
	var cs []c
	for _, md := range s.Modes {
		q := 0.0
		for i, b := range sp.B {
			q += b * md.Vector[i]
		}
		tot := 0.0
		for i := range md.Vector {
			for j := range md.Vector {
				tot += md.Vector[i] * s.K[i][j] * md.Vector[j]
			}
		}
		cs = append(cs, c{md.Wavenumber, sp.K * q * q / tot})
	}
	sort.Slice(cs, func(i, j int) bool { return cs[i].sh > cs[j].sh })
	if cs[0].sh < 0.3 {
		return 0
	}
	return cs[0].f
}

func round4(v float64) float64 { return math.Round(v*1e4) / 1e4 }

func splitKeys(s string) []string {
	if s == "" {
		return nil
	}
	return strings.Split(s, ",")
}
