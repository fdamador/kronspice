// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
)

// Stacked base pairs. Two, three and four consecutive Watson–Crick pairs
// are cut from the B-DNA duplex with their duplex rise and twist, the
// backbone removed and each base methyl-capped (prep/prep.py). The
// forces between bases of different pairs — van der Waals and
// electrostatics, AMBER's whole description of stacking — form a third
// coupling group in the circuit, G_ST, beside the bases' own couplings
// (G_R) and the pairing (G_HB).
//
// What stacking does to the spectrum: identical pairs that push on their
// neighbours turn each single-pair line into a band of nearby lines (one
// per pair), and the pairs moving against each other as near-rigid bodies
// add new low-frequency modes that a single pair does not have.

var stackOrder = []string{"AT", "AT2", "AT3", "AT4", "GC", "GC2", "GC3", "GC4", "polyA2", "polyA3", "polyA4", "ApT"}

// pairOf returns the base-pair index of atom a (0 for a lone pair, n−1 for
// prefix "A<n>").
func pairOf(m *Molecule, a int) int {
	pre := strings.SplitN(m.Atoms[a].Name, ":", 2)[0]
	n, err := strconv.Atoi(pre[1:])
	if err != nil {
		return 0
	}
	return n - 1
}

func numPairs(m *Molecule) int {
	n := 0
	for a := range m.Atoms {
		if p := pairOf(m, a) + 1; p > n {
			n = p
		}
	}
	return n
}

// StackGeometry returns the rise [Å] and twist [°] of each step: rise is the
// distance between pair centroids along the stack axis, twist the angle
// about the axis between the C1'(methyl)–C1' vectors of neighbouring pairs
// (the standard step parameters, simplified: El Hassan & Calladine, J. Mol.
// Biol. 251, 648, 1995, define them from base-pair frames).
func StackGeometry(m *Molecule, x []float64) (rise, twist []float64) {
	np := numPairs(m)
	if np < 2 {
		return nil, nil
	}
	cen := make([][3]float64, np)
	cnt := make([]float64, np)
	cms := make([][][3]float64, np)
	for a, at := range m.Atoms {
		p := pairOf(m, a)
		if at.Element != "H" {
			for d := 0; d < 3; d++ {
				cen[p][d] += x[3*a+d]
			}
			cnt[p]++
		}
		if strings.HasSuffix(at.Name, ":CM") {
			cms[p] = append(cms[p], [3]float64{x[3*a], x[3*a+1], x[3*a+2]})
		}
	}
	for p := range cen {
		cen[p] = scale(cen[p], 1/cnt[p])
	}
	z := sub(cen[np-1], cen[0])
	z = scale(z, 1/norm(z))
	perp := func(v [3]float64) [3]float64 { return sub(v, scale(z, dot(v, z))) }
	for p := 0; p+1 < np; p++ {
		rise = append(rise, dot(sub(cen[p+1], cen[p]), z))
		u := perp(sub(cms[p][1], cms[p][0]))
		v := perp(sub(cms[p+1][1], cms[p+1][0]))
		twist = append(twist, math.Atan2(dot(cross(u, v), z), dot(u, v))*180/math.Pi)
	}
	return rise, twist
}

// PairBodies describes how a mode moves the base pairs as rigid bodies:
// the share of the mode's kinetic energy carried by rigid motions of the
// pairs (each pair's best-fit translation and rotation, mass-weighted least
// squares), split into motion along the stack axis (rise), sideways in the
// pair planes (shift and slide), rotation about the axis (twist) and
// rotation about in-plane axes (tilt and roll).
type PairBodies struct {
	Rigid, Rise, Slide, Twist, TiltRoll float64
}

func (s *State) PairBodies(md Mode) PairBodies {
	m := s.Mol
	np := numPairs(m)
	var out PairBodies
	if np < 2 {
		return out
	}
	// Stack axis from the minimized centroids.
	cen := make([][3]float64, np)
	mass := make([]float64, np)
	for a, at := range m.Atoms {
		p := pairOf(m, a)
		for d := 0; d < 3; d++ {
			cen[p][d] += at.Mass * s.X[3*a+d]
		}
		mass[p] += at.Mass
	}
	for p := range cen {
		cen[p] = scale(cen[p], 1/mass[p])
	}
	z := sub(cen[np-1], cen[0])
	z = scale(z, 1/norm(z))
	total := 0.0
	for a, at := range m.Atoms {
		v := [3]float64{md.Vector[3*a], md.Vector[3*a+1], md.Vector[3*a+2]}
		total += at.Mass * dot(v, v)
	}
	for p := 0; p < np; p++ {
		// Normal equations for t, w: d_a ≈ t + w × r_a (r_a from the centroid).
		var A [6][6]float64
		var b [6]float64
		for a, at := range m.Atoms {
			if pairOf(m, a) != p {
				continue
			}
			r := sub([3]float64{s.X[3*a], s.X[3*a+1], s.X[3*a+2]}, cen[p])
			d := [3]float64{md.Vector[3*a], md.Vector[3*a+1], md.Vector[3*a+2]}
			// Row block: J = [I | −[r]×] maps (t, w) to t + w×r.
			var J [3][6]float64
			for k := 0; k < 3; k++ {
				J[k][k] = 1
			}
			J[0][4], J[0][5] = r[2], -r[1]
			J[1][3], J[1][5] = -r[2], r[0]
			J[2][3], J[2][4] = r[1], -r[0]
			for i := 0; i < 6; i++ {
				for k := 0; k < 3; k++ {
					b[i] += at.Mass * J[k][i] * d[k]
					for j := 0; j < 6; j++ {
						A[i][j] += at.Mass * J[k][i] * J[k][j]
					}
				}
			}
		}
		sol := solve6(A, b)
		t := [3]float64{sol[0], sol[1], sol[2]}
		w := [3]float64{sol[3], sol[4], sol[5]}
		// Kinetic shares: translation m|t|², rotation wᵀIw (I from the same
		// normal equations' rotational block).
		tz := dot(t, z)
		out.Rise += mass[p] * tz * tz / total
		out.Slide += mass[p] * (dot(t, t) - tz*tz) / total
		// Rotation energy wᵀIw, split between the axial component (twist)
		// and the in-plane one (tilt/roll) in proportion to their own
		// energies; the cross term is shared the same way.
		wz := scale(z, dot(w, z))
		wp := sub(w, wz)
		quad := func(u, v [3]float64) float64 {
			q := 0.0
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					q += u[i] * A[3+i][3+j] * v[j]
				}
			}
			return q
		}
		iw, ez, ep := quad(w, w), quad(wz, wz), quad(wp, wp)
		if ez+ep > 0 {
			out.Twist += iw * ez / (ez + ep) / total
			out.TiltRoll += iw * ep / (ez + ep) / total
		}
	}
	out.Rigid = out.Rise + out.Slide + out.Twist + out.TiltRoll
	return out
}

// Describe names a stack mode: the rigid motion of the pairs when it
// carries at least half the kinetic energy, else "inside the pairs".
func (pb PairBodies) Describe() string {
	if pb.Rigid < 0.5 {
		return ""
	}
	parts := []struct {
		name string
		v    float64
	}{{"rise (pairs moving apart and together)", pb.Rise}, {"shift/slide (pairs sliding over each other)", pb.Slide},
		{"twist (pairs turning about the axis)", pb.Twist}, {"tilt/roll (pairs rocking)", pb.TiltRoll}}
	sort.Slice(parts, func(i, j int) bool { return parts[i].v > parts[j].v })
	out := fmt.Sprintf("%s %.0f%%", parts[0].name, 100*parts[0].v/pb.Rigid)
	if parts[1].v/pb.Rigid >= 0.2 {
		out += fmt.Sprintf(", %s %.0f%%", parts[1].name, 100*parts[1].v/pb.Rigid)
	}
	return out
}

func solve6(A [6][6]float64, b [6]float64) [6]float64 {
	for k := 0; k < 6; k++ {
		p := k
		for i := k + 1; i < 6; i++ {
			if math.Abs(A[i][k]) > math.Abs(A[p][k]) {
				p = i
			}
		}
		A[k], A[p] = A[p], A[k]
		b[k], b[p] = b[p], b[k]
		for i := k + 1; i < 6; i++ {
			f := A[i][k] / A[k][k]
			for j := k; j < 6; j++ {
				A[i][j] -= f * A[k][j]
			}
			b[i] -= f * b[k]
		}
	}
	var x [6]float64
	for i := 5; i >= 0; i-- {
		s := b[i]
		for j := i + 1; j < 6; j++ {
			s -= A[i][j] * x[j]
		}
		x[i] = s / A[i][i]
	}
	return x
}

// ModalAbsorption is the circuit's absorbed power from its normal modes:
// with mass-proportional damping γ the modes stay uncoupled, and a unit
// field along each axis gives ½·Σ_k |L_k·q|²·γ/(γ² + (ω − ω_k²/ω)²),
// averaged over x, y, z, in the circuit's units (the same numbers as
// Absorption and DirectResponse, without solving the circuit at every
// point). All 3N modes are summed, the six rigid-body ones included.
func (s *State) ModalAbsorption(nu []float64, fwhm float64) []float64 {
	m := s.Mol
	n := len(s.K)
	conv := 2 * math.Pi * wavenumber
	gamma := 2 * math.Pi * fwhm / conv
	h := make([][]float64, n)
	for i := range h {
		h[i] = make([]float64, n)
		for j := range h[i] {
			h[i][j] = s.K[i][j] / math.Sqrt(m.Atoms[i/3].Mass*m.Atoms[j/3].Mass)
		}
	}
	vals, vecs := Eigen(h)
	amp := make([]float64, n) // Σ over axes of |L_k·q_axis|²
	for k := 0; k < n; k++ {
		var a [3]float64
		for i := 0; i < n; i++ {
			a[i%3] += m.Atoms[i/3].Charge * vecs[i][k] / math.Sqrt(m.Atoms[i/3].Mass)
		}
		amp[k] = a[0]*a[0] + a[1]*a[1] + a[2]*a[2]
	}
	out := make([]float64, len(nu))
	for p, v := range nu {
		w := 2 * math.Pi * v / conv
		for k := 0; k < n; k++ {
			d := w - vals[k]/w
			out[p] += 0.5 * amp[k] * gamma / (gamma*gamma + d*d) / 3
		}
	}
	return out
}

// circuitPoints runs kronspice's AC analysis of the stack circuit at a few
// frequencies (a short linear sweep from nu0 to nu1) and returns the
// absorbed power averaged over the three field directions: the check that
// the circuit itself gives the modal spectrum.
func (s *State) circuitPoints(nu0, nu1 float64, points int, fwhm float64) ([]float64, []float64, error) {
	return s.Absorption(nu0, nu1, points, fwhm)
}

type stackResult struct {
	key   string
	s     *State
	err   error
	rise  []float64
	twist []float64
	rise0 []float64
	tw0   []float64
}

// RunStacks analyses the lone pairs and the 2-, 3- and 4-pair stacks and
// writes stacks/: per-stack modes and spectra, a summary, and netlists of
// the 2-pair stacks.
func RunStacks(out string, ms map[string]*Molecule, keys []string) error {
	if len(keys) == 0 {
		keys = stackOrder
	}
	dir := filepath.Join(out, "stacks")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		return err
	}
	res := make([]stackResult, len(keys))
	var wg sync.WaitGroup
	for i, key := range keys {
		wg.Add(1)
		go func(i int, key string) {
			defer wg.Done()
			m := ms[key]
			r := stackResult{key: key}
			if m == nil {
				r.err = fmt.Errorf("no molecule %q", key)
				res[i] = r
				return
			}
			r.rise0, r.tw0 = StackGeometry(m, m.Coords())
			r.s, r.err = Prepare(m)
			if r.err == nil {
				r.rise, r.twist = StackGeometry(m, r.s.X)
			}
			res[i] = r
		}(i, key)
	}
	wg.Wait()
	for _, r := range res {
		if r.err != nil {
			return r.err
		}
	}
	fmt.Println()
	fmt.Println("Stacked base pairs (bases cut from the B-DNA duplex, backbone removed, methyl caps; AMBER OL15, gas phase)")
	nu := make([]float64, SweepPts)
	for k := range nu {
		nu[k] = SweepFrom + (SweepTo-SweepFrom)*float64(k)/float64(SweepPts-1)
	}
	var sum strings.Builder
	sum.WriteString("stack,pairs,atoms,vibrations,springs,couplings_inside_bases,pairing_couplings,stacking_couplings,rise_start_A,rise_min_A,twist_start_deg,twist_min_deg,lowest_cm-1,modes_below_200,stacking_modes_below_200\n")
	var low strings.Builder
	low.WriteString("stack,mode,wavenumber_cm-1,ir_intensity_km_mol,rigid_pair_share,stacking_energy_share,pairing_energy_share,motion\n")
	for _, r := range res {
		s := r.s
		springs := s.Springs()
		entries, rst := s.ResidualEntries(springs)
		np := numPairs(s.Mol)
		var mcsv strings.Builder
		mcsv.WriteString("mode,wavenumber_cm-1,ir_intensity_km_mol,release_rate_A_s-1,lifetime_ms,pairing_energy_share,stacking_energy_share,rigid_pair_share,stack_motion,assignment\n")
		below, stackBelow := 0, 0
		for k, md := range s.Modes {
			label, _ := s.Label(md, springs, entries)
			pb := s.PairBodies(md)
			motion := pb.Describe()
			rel := ReleaseOf(md)
			fmt.Fprintf(&mcsv, "%d,%.1f,%.2f,%.4g,%.4g,%.3f,%.3f,%.3f,%q,%q\n", k+1, md.Wavenumber, md.Intensity, rel.A, 1e3*rel.Lifetime, md.Inter, md.Stack, pb.Rigid, motion, label)
			if md.Wavenumber < 200 {
				below++
				if motion != "" {
					stackBelow++
				}
				fmt.Fprintf(&low, "%s,%d,%.1f,%.2f,%.3f,%.3f,%.3f,%q\n", r.key, k+1, md.Wavenumber, md.Intensity, pb.Rigid, md.Stack, md.Inter, firstNonEmpty(motion, label))
			}
		}
		if err := os.WriteFile(filepath.Join(dir, "modes_"+strings.ToLower(r.key)+".csv"), []byte(mcsv.String()), 0o644); err != nil {
			return err
		}
		ab := s.ModalAbsorption(nu, FWHMcm)
		q, th := Emission(nu, ab)
		var sp strings.Builder
		sp.WriteString("wavenumber_cm-1,absorption,release_one_quantum,release_thermal_300K,absorbed_power\n")
		na := append([]float64(nil), ab...)
		normalize(na)
		for k := range nu {
			fmt.Fprintf(&sp, "%.1f,%.6g,%.6g,%.6g,%.6g\n", nu[k], na[k], q[k], th[k], ab[k])
		}
		if err := os.WriteFile(filepath.Join(dir, "spectra_"+strings.ToLower(r.key)+".csv"), []byte(sp.String()), 0o644); err != nil {
			return err
		}
		if np == 2 {
			src, _ := s.Netlist(NetlistOptions{FWHM: FWHMcm, Field: "x", Analysis: fmt.Sprintf(".ac lin %d %.9g %.9g", SweepPts, SweepFrom/(2*math.Pi*wavenumber), SweepTo/(2*math.Pi*wavenumber))})
			if err := os.WriteFile(filepath.Join(dir, strings.ToLower(r.key)+".cir"), []byte(src), 0o644); err != nil {
				return err
			}
		}
		f := func(v []float64) string {
			if v == nil {
				return ""
			}
			var p []string
			for _, x := range v {
				p = append(p, fmt.Sprintf("%.2f", x))
			}
			return strings.Join(p, " ")
		}
		fmt.Fprintf(&sum, "%s,%d,%d,%d,%d,%d,%d,%d,%q,%q,%q,%q,%.1f,%d,%d\n", r.key, np, len(s.Mol.Atoms), len(s.Modes), len(springs), rst.Intra, rst.Inter, rst.Stack,
			f(r.rise0), f(r.rise), f(r.tw0), f(r.twist), s.Modes[0].Wavenumber, below, stackBelow)
		fmt.Printf("  %-4s %d pair(s), %3d atoms, %3d vibrations; rise %s Å (duplex %s), twist %s° (duplex %s); %d modes below 200 cm⁻¹, %d of them rigid pair motions\n",
			r.key, np, len(s.Mol.Atoms), len(s.Modes), f(r.rise), f(r.rise0), f(r.twist), f(r.tw0), below, stackBelow)
	}
	tag := ""
	if len(keys) != len(stackOrder) {
		tag = "_" + strings.ToLower(strings.Join(keys, "_"))
	}
	if err := os.WriteFile(filepath.Join(dir, "summary"+tag+".csv"), []byte(sum.String()), 0o644); err != nil {
		return err
	}
	if err := os.WriteFile(filepath.Join(dir, "low_modes"+tag+".csv"), []byte(low.String()), 0o644); err != nil {
		return err
	}
	// Twin checks on the 2-pair stacks: modal spectrum = direct solve =
	// kronspice AC of the circuit, near the strongest low and high lines.
	fmt.Println("  checks (2-pair stacks): modal spectrum against the direct solve and against kronspice's AC analysis of the circuit")
	for _, r := range res {
		if numPairs(r.s.Mol) != 2 {
			continue
		}
		s := r.s
		strong := s.Modes[0]
		for _, md := range s.Modes {
			if md.Intensity > strong.Intensity {
				strong = md
			}
		}
		lowm := s.Modes[0]
		for _, md := range s.Modes {
			if md.Wavenumber < 200 && md.Intensity > lowm.Intensity {
				lowm = md
			}
		}
		worstD, worstC := 0.0, 0.0
		for _, c := range []float64{lowm.Wavenumber, strong.Wavenumber} {
			cnu, cab, err := s.circuitPoints(c-4, c+4, 5, FWHMcm)
			if err != nil {
				return err
			}
			modal := s.ModalAbsorption(cnu, FWHMcm)
			for k := range cnu {
				d := 0.0
				for ax := 0; ax < 3; ax++ {
					d += s.DirectResponse(cnu[k], FWHMcm, ax) / 3
				}
				worstD = math.Max(worstD, math.Abs(modal[k]-d)/math.Abs(d))
				worstC = math.Max(worstC, math.Abs(modal[k]-cab[k])/math.Abs(cab[k]))
			}
		}
		fmt.Printf("    %s at %.1f and %.1f cm⁻¹: modal vs direct %.1e, modal vs circuit %.1e (relative)\n", r.key, lowm.Wavenumber, strong.Wavenumber, worstD, worstC)
		if worstD > 1e-5 || worstC > 1e-5 {
			return fmt.Errorf("%s: modal spectrum disagrees with the circuit (%.2g, %.2g)", r.key, worstD, worstC)
		}
	}
	return nil
}

func firstNonEmpty(a, b string) string {
	if a != "" {
		return a
	}
	return b
}
