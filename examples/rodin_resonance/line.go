// SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
// Copyright (C) 2026 Franklin Amador

package main

import (
	"fmt"
	"math"

	"github.com/fdamador/kronspice/circuit"
	"github.com/fdamador/kronspice/em"
)

// The distributed (sectional) model of the pair, as pulsed_bifilar does
// for its coils: every pass of every circuit of each winding is a
// section, with the inductance matrix from the geometry (pass matrix
// repeated over the circuits of a bundle, plus the bundle leakage) and
// declared capacitances:
//
//   - wire to wire inside a bundle: the wire of circuit c and that of
//     circuit c+1 lie side by side along every pass (layer order), with the
//     two-wire capacitance c′ = πε0εr/acosh(D/2r_w) per metre (enamel
//     0.03 mm per wire, εr = 3), split half to each end of the pass;
//   - between the A and B bundles: they touch only where they cross, so
//     c′·ℓ·ABContact per pass (ABContact = 0.1, declared);
//   - to the surroundings: the capacitance of a sphere of the coil's outer
//     radius, 4πε0·R, shared equally by all the nodes.
//
// This is the same declared level of approximation as pulsed_bifilar
// (sections equipotential inside, no radiation).
const (
	EnamelThick = 0.03e-3
	EnamelEps   = 3.0
	ABContact   = 0.1
)

// WireC is the two-wire capacitance per metre of two touching enamelled
// wires.
func WireC() float64 {
	d := 2*WireRadius + 2*EnamelThick
	return math.Pi * em.Eps0 * EnamelEps / math.Acosh(d/(2*WireRadius))
}

// Topology of the sectional network.
type Topology int

const (
	PlainRing   Topology = iota // A closed on itself, B closed on itself
	CrossedRing                 // A's end to B's start, B's end to A's start: one loop through both
	OpenA                       // A alone: start grounded (driven by a stiff source), end open
	SeriesTor                   // A start grounded, A end → B start, B end open
	SeriesDip                   // A start grounded, A end → B end, B start open
	DrivenPair                  // A and B each from a source at the start (B: at the end for Dipole), the other end grounded
)

func (t Topology) String() string {
	return [...]string{"plain ring", "crossed ring", "winding A open", "series, toroidal sense", "series, dipole sense", "driven pair"}[t]
}

// Net is the sectional network: sections (branches) between nodes.
type Net struct {
	Top      Topology
	NS       int   // sections per winding (Nc·Sec)
	Nodes    int   // non-ground nodes after merging
	A, B     []int // node id of a_j, b_j (j = 0..NS); −1 = ground
	L        [][]float64
	R        []float64 // section resistance [Ω]
	Len      []float64 // section wire length [m]
	CNode    [][]float64
	Branches [][2]int // per section: node ids (from, to)
	UseB     bool
}

// NewNet builds the network of a model in a topology. sense selects the
// driven end of B for DrivenPair (Toroidal: b_0, Dipole: b_NS).
func (m *Model) NewNet(top Topology, sense int, f float64) *Net {
	p := m.P
	ns := p.Nc * p.Sec
	n := &Net{Top: top, NS: ns, UseB: top != OpenA}
	// raw node ids: a_j = j, b_j = ns+1+j
	raw := 2 * (ns + 1)
	id := make([]int, raw)
	for i := range id {
		id[i] = i
	}
	parent := func(i int) int {
		for id[i] != i {
			i = id[i]
		}
		return i
	}
	union := func(i, j int) { id[parent(i)] = parent(j) }
	const G = -2 // marker for ground
	ground := map[int]bool{}
	a := func(j int) int { return j }
	b := func(j int) int { return ns + 1 + j }
	switch top {
	case PlainRing:
		union(a(ns), a(0))
		union(b(ns), b(0))
	case CrossedRing:
		union(a(ns), b(0))
		union(b(ns), a(0))
	case OpenA:
		ground[a(0)] = true
	case SeriesTor:
		ground[a(0)] = true
		union(a(ns), b(0))
	case SeriesDip:
		ground[a(0)] = true
		union(a(ns), b(ns))
	case DrivenPair:
		ground[a(ns)] = true
		if sense == Toroidal {
			ground[b(ns)] = true
		} else {
			ground[b(0)] = true
		}
	}
	// number the surviving nodes
	num := map[int]int{}
	final := make([]int, raw)
	for i := 0; i < raw; i++ {
		if i >= ns+1 && !n.UseB {
			final[i] = G
			continue
		}
		r := parent(i)
		if ground[r] || ground[i] {
			final[i] = -1
			continue
		}
		if _, ok := num[r]; !ok {
			num[r] = len(num)
		}
		final[i] = num[r]
	}
	n.Nodes = len(num)
	n.A = final[:ns+1]
	n.B = final[ns+1:]
	// sections and inductance matrix
	nb := ns
	if n.UseB {
		nb = 2 * ns
	}
	n.L = newMat(nb)
	n.R = make([]float64, nb)
	n.Len = make([]float64, nb)
	n.Branches = make([][2]int, nb)
	pass := func(g int) (w, c, q int) { // winding, circuit, pass
		w = g / ns
		c = (g % ns) / p.Sec
		q = g % p.Sec
		return
	}
	for g := 0; g < nb; g++ {
		wg, cg, qg := pass(g)
		for h := 0; h < nb; h++ {
			wh, ch, qh := pass(h)
			n.L[g][h] = m.Pass[wg*p.Sec+qg][wh*p.Sec+qh]
			if g == h {
				n.L[g][h] += m.Delta[wg*p.Sec+qg]
			}
			_ = cg
			_ = ch
		}
		n.Len[g] = m.SecLen[wg*p.Sec+qg]
		n.R[g] = WireResistance(n.Len[g], f)
		j := g % ns
		if wg == 0 {
			n.Branches[g] = [2]int{n.A[j], n.A[j+1]}
		} else {
			n.Branches[g] = [2]int{n.B[j], n.B[j+1]}
		}
	}
	// capacitances
	n.CNode = newMat(n.Nodes)
	addC := func(i, j int, c float64) {
		if i == j {
			return
		}
		if i >= 0 {
			n.CNode[i][i] += c
		}
		if j >= 0 {
			n.CNode[j][j] += c
		}
		if i >= 0 && j >= 0 {
			n.CNode[i][j] -= c
			n.CNode[j][i] -= c
		}
	}
	cw := WireC()
	nodesRaw := ns + 1
	if n.UseB {
		nodesRaw *= 2
	}
	cg := 4 * math.Pi * em.Eps0 * OuterRadius / float64(nodesRaw)
	for j := 0; j <= ns; j++ {
		addC(n.A[j], -1, cg)
		if n.UseB {
			addC(n.B[j], -1, cg)
		}
	}
	for w := 0; w < 2; w++ {
		if w == 1 && !n.UseB {
			break
		}
		nd := n.A
		if w == 1 {
			nd = n.B
		}
		for c := 0; c+1 < p.Nc; c++ {
			for q := 0; q < p.Sec; q++ {
				ct := cw * m.SecLen[w*p.Sec+q]
				j, k := c*p.Sec+q, (c+1)*p.Sec+q
				addC(nd[j], nd[k], ct/2)
				addC(nd[j+1], nd[k+1], ct/2)
			}
		}
	}
	if n.UseB {
		for c := 0; c < p.Nc; c++ {
			for q := 0; q < p.Sec; q++ {
				ct := cw * ABContact * 0.5 * (m.SecLen[q] + m.SecLen[p.Sec+q])
				j := c*p.Sec + q
				addC(n.A[j], n.B[j], ct/2)
				addC(n.A[j+1], n.B[j+1], ct/2)
			}
		}
	}
	return n
}

// SetFreq sets the section resistances for frequency f (skin effect).
func (n *Net) SetFreq(f float64) {
	for g := range n.R {
		n.R[g] = WireResistance(n.Len[g], f)
	}
}

// Modes returns the lossless network's resonant frequencies (ascending)
// and node-voltage mode shapes (columns): ω²·C·v = Dᵀ·L⁻¹·D·v.
func (n *Net) Modes() ([]float64, [][]float64, error) {
	li, err := spdInverse(n.L)
	if err != nil {
		return nil, nil, err
	}
	nb := len(n.L)
	// D: branch voltage = v(from) − v(to)
	k := newMat(n.Nodes)
	for g := 0; g < nb; g++ {
		for h := 0; h < nb; h++ {
			if li[g][h] == 0 {
				continue
			}
			for _, x := range []struct {
				node int
				s    float64
			}{{n.Branches[g][0], 1}, {n.Branches[g][1], -1}} {
				if x.node < 0 {
					continue
				}
				for _, y := range []struct {
					node int
					s    float64
				}{{n.Branches[h][0], 1}, {n.Branches[h][1], -1}} {
					if y.node < 0 {
						continue
					}
					k[x.node][y.node] += x.s * y.s * li[g][h]
				}
			}
		}
	}
	vals, vecs, err := genEigen(k, n.CNode)
	if err != nil {
		return nil, nil, err
	}
	f := make([]float64, len(vals))
	for i, v := range vals {
		f[i] = math.Sqrt(math.Max(0, v)) / (2 * math.Pi)
	}
	return f, vecs, nil
}

// TripParity is, for a mode of the crossed ring, the correlation between
// the voltages along A and along B (node j of A against node j of B): +1
// when the wave comes back the same way up after one trip (one winding),
// −1 when it comes back inverted (the 720° modes).
func (n *Net) TripParity(v [][]float64, mode int) float64 {
	sab, saa, sbb := 0.0, 0.0, 0.0
	for j := 0; j < n.NS; j++ {
		x, y := n.A[j], n.B[j]
		if x < 0 || y < 0 {
			continue
		}
		a, b := v[x][mode], v[y][mode]
		sab += a * b
		saa += a * a
		sbb += b * b
	}
	return sab / math.Sqrt(saa*sbb)
}

// Builder writes the network into a kronspice circuit: sections as
// inductors with their resistances, every mutual as a K element, node
// capacitances as capacitors. Node names "n<id>".
func (n *Net) Builder(title string) *circuit.Builder {
	b := circuit.NewBuilder(title)
	nm := func(i int) string {
		if i < 0 {
			return "0"
		}
		return fmt.Sprintf("n%d", i)
	}
	nb := len(n.L)
	for g := 0; g < nb; g++ {
		b.Inductor(fmt.Sprintf("L%d", g), nm(n.Branches[g][0]), fmt.Sprintf("m%d", g), n.L[g][g])
		b.Resistor(fmt.Sprintf("R%d", g), fmt.Sprintf("m%d", g), nm(n.Branches[g][1]), n.R[g])
	}
	k := 0
	for g := 0; g < nb; g++ {
		for h := g + 1; h < nb; h++ {
			kk := n.L[g][h] / math.Sqrt(n.L[g][g]*n.L[h][h])
			if math.Abs(kk) < 1e-5 {
				continue
			}
			b.Mutual(fmt.Sprintf("K%d", k), fmt.Sprintf("L%d", g), fmt.Sprintf("L%d", h), kk)
			k++
		}
	}
	c := 0
	for i := 0; i < n.Nodes; i++ {
		row := 0.0
		for j := 0; j < n.Nodes; j++ {
			row += n.CNode[i][j]
			if j > i && n.CNode[i][j] < 0 {
				b.Capacitor(fmt.Sprintf("C%d", c), nm(i), nm(j), -n.CNode[i][j])
				c++
			}
		}
		if row > 1e-18 {
			b.Capacitor(fmt.Sprintf("C%d", c), nm(i), "0", row)
			c++
		}
	}
	return b
}

// NodeName is the kronspice name of node id i.
func NodeName(i int) string {
	if i < 0 {
		return "0"
	}
	return fmt.Sprintf("n%d", i)
}
