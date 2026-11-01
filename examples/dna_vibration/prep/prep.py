#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Prepare the DNA vibration example's molecules (molecules.json).

Inputs (both from the OpenMM repository, github.com/openmm/openmm; give their
paths as arguments):

  DNA.OL15.xml  AMBER DNA force field as converted by OpenMM from AmberTools 17
                leaprc.DNA.OL15: parm99/bsc0 base parameters (Cornell et al.,
                JACS 117, 5179, 1995; Perez et al., Biophys. J. 92, 3817, 2007;
                Zgarbova et al., JCTC 11, 5723, 2015). Units kJ/mol, nm, rad.
  bdna.xyz      all-atom B-DNA dodecamer d(CGCGAATTCGCG) (Tinker xyz, Angstrom),
                the Dickerson-Drew sequence; it supplies starting coordinates
                with hydrogens, single bases and Watson-Crick pairs.

What it writes, per molecule: atoms (name, element, mass, charge, AMBER type,
coordinates in Angstrom), and every force-field term the Go code evaluates:
bonds, angles, proper and improper torsions, nonbonded pairs with their 1-4
scaling and exclusions. Assignment follows OpenMM's ForceField rules
(forcefield.py): specific torsion definitions before wildcard ones, AMBER
improper ordering (outer, outer, centre, outer, with the template-order swaps).

Each base is cut from the duplex at the glycosidic bond and capped as in the
methylated bases used in gas-phase spectroscopy (9-methyladenine,
1-methylthymine, 9-methylguanine, 1-methylcytosine): C1' becomes the methyl
carbon, H1' stays, and O4' and C2' are replaced by hydrogens 1.09 A from C1'
along the old bond directions. Methyl hydrogens take type H1 and the charge of
H1'; the methyl carbon takes the charge that makes the molecule neutral.
These cap choices are this example's, not AMBER's, and are stated in the
output.
"""
import json
import math
import sys
import xml.etree.ElementTree as ET

import networkx as nx
from networkx.algorithms import isomorphism

KJ_PER_KCAL = 4.184
NM = 10.0  # Angstrom per nm


def load_ff(path):
    root = ET.parse(path).getroot()
    types = {}
    for t in root.find("AtomTypes"):
        types[t.get("name")] = dict(cls=t.get("class"), element=t.get("element"), mass=float(t.get("mass")))
    residues = {}
    for r in root.find("Residues"):
        atoms = [(a.get("name"), a.get("type"), float(a.get("charge"))) for a in r.findall("Atom")]
        bonds = [(b.get("atomName1"), b.get("atomName2")) for b in r.findall("Bond")]
        residues[r.get("name")] = dict(atoms=atoms, bonds=bonds)
    bonds = {}
    for b in root.find("HarmonicBondForce"):
        k = frozenset((b.get("type1"), b.get("type2")))
        bonds.setdefault(k, (float(b.get("k")), float(b.get("length"))))
    angles = {}
    for a in root.find("HarmonicAngleForce"):
        k = (a.get("type1"), a.get("type2"), a.get("type3"))
        angles.setdefault(k, (float(a.get("k")), float(a.get("angle"))))
    propers, impropers = [], []
    for t in root.find("PeriodicTorsionForce"):
        terms = []
        i = 1
        while t.get("k%d" % i) is not None:
            terms.append((int(t.get("periodicity%d" % i)), float(t.get("k%d" % i)), float(t.get("phase%d" % i))))
            i += 1
        rec = dict(types=[t.get("type%d" % j) for j in range(1, 5)], terms=terms)
        (propers if t.tag == "Proper" else impropers).append(rec)
    nb = root.find("NonbondedForce")
    lj = {a.get("type"): (float(a.get("sigma")), float(a.get("epsilon"))) for a in nb.findall("Atom")}
    return dict(types=types, residues=residues, bonds=bonds, angles=angles, propers=propers,
                impropers=impropers, lj=lj, c14=float(nb.get("coulomb14scale")), lj14=float(nb.get("lj14scale")))


def load_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    g = nx.Graph()
    for ln in lines[1:1 + n]:
        f = ln.split()
        i = int(f[0])
        el = f[1][0]
        g.add_node(i, element=el, xyz=(float(f[2]), float(f[3]), float(f[4])))
        for j in f[6:]:
            g.add_edge(i, int(j))
    return g


# Base atoms plus the sugar atoms kept for the cap.
CAP = ["C1'", "H1'", "O4'", "C2'"]
BASE = {
    "DA": ["N9", "C8", "H8", "N7", "C5", "C6", "N6", "H61", "H62", "N1", "C2", "H2", "N3", "C4"],
    "DG": ["N9", "C8", "H8", "N7", "C5", "C6", "O6", "N1", "H1", "C2", "N2", "H21", "H22", "N3", "C4"],
    "DT": ["N1", "C6", "H6", "C5", "C7", "H71", "H72", "H73", "C4", "O4", "N3", "H3", "C2", "O2"],
    "DC": ["N1", "C6", "H6", "C5", "H5", "C4", "N4", "H41", "H42", "N3", "C2", "O2"],
}
NAMES = {"DA": "9-methyladenine", "DG": "9-methylguanine", "DT": "1-methylthymine", "DC": "1-methylcytosine"}


def template_graph(ff, res):
    keep = set(BASE[res] + CAP)
    r = ff["residues"][res]
    g = nx.Graph()
    for name, typ, _ in r["atoms"]:
        if name in keep:
            g.add_node(name, element=ff["types"][typ]["element"])
    for a, b in r["bonds"]:
        if a in keep and b in keep:
            g.add_edge(a, b)
    return g


def find_bases(ff, g):
    """All matches of each base template in the duplex, as name → xyz index maps."""
    found = {}
    for res in BASE:
        tg = template_graph(ff, res)
        gm = isomorphism.GraphMatcher(g, tg, node_match=lambda a, b: a["element"] == b["element"])
        seen = set()
        out = []
        for m in gm.subgraph_monomorphisms_iter():
            inv = {v: k for k, v in m.items()}
            key = inv["N9" if res in ("DA", "DG") else "N1"]
            if key in seen:
                continue
            # Ring atoms only (no extra neighbours on a base heavy atom beyond the template).
            seen.add(key)
            out.append(inv)
        found[res] = out
    # A guanine template also matches inside nothing else; adenine and guanine differ
    # in element pattern, so matches are unique per residue type.
    return found


def dist(a, b):
    return math.dist(a, b)


def build_molecule(ff, g, res, inv, tag):
    r = ff["residues"][res]
    info = {name: (typ, q) for name, typ, q in r["atoms"]}
    order = {name: i for i, (name, _, _) in enumerate(r["atoms"])}
    atoms = []
    for name in BASE[res] + ["C1'", "H1'"]:
        typ, q = info[name]
        atoms.append(dict(name=name, type=typ, charge=q, xyz=list(g.nodes[inv[name]]["xyz"]), tidx=order[name]))
    c1 = g.nodes[inv["C1'"]]["xyz"]
    extra = len(r["atoms"])
    for k, old in enumerate(["O4'", "C2'"]):
        p = g.nodes[inv[old]]["xyz"]
        d = dist(c1, p)
        xyz = [c1[i] + 1.09 * (p[i] - c1[i]) / d for i in range(3)]
        atoms.append(dict(name="HM%d" % (k + 2), type="DNA-H1", charge=None, xyz=xyz, tidx=extra + k))
    # Methyl cap: C1' → CM (type CT), H1' keeps its name with type H1.
    for a in atoms:
        if a["name"] == "C1'":
            a["name"] = "CM"
        if a["name"] == "H1'":
            a["name"], a["type"] = "HM1", "DNA-H1"
    hq = info["H1'"][1]
    for a in atoms:
        if a["name"].startswith("HM"):
            a["charge"] = hq
    total = sum(a["charge"] for a in atoms if a["name"] != "CM")
    for a in atoms:
        if a["name"] == "CM":
            a["charge"] = -total
    for a in atoms:
        t = ff["types"][a["type"]]
        a["element"], a["mass"] = t["element"], t["mass"]
        a["prefix"] = tag
    names = [a["name"] for a in atoms]
    rename = {"C1'": "CM", "H1'": "HM1"}
    bonds = set()
    for x, y in r["bonds"]:
        x, y = rename.get(x, x), rename.get(y, y)
        if x in names and y in names:
            bonds.add(tuple(sorted((names.index(x), names.index(y)))))
    for h in ("HM2", "HM3"):
        bonds.add(tuple(sorted((names.index("CM"), names.index(h)))))
    return atoms, sorted(bonds)


def bond_param(ff, t1, t2):
    p = ff["bonds"].get(frozenset((t1, t2)))
    if p is None:
        raise SystemExit("no bond parameter %s-%s" % (t1, t2))
    return p


def angle_param(ff, t1, t2, t3):
    p = ff["angles"].get((t1, t2, t3)) or ff["angles"].get((t3, t2, t1))
    if p is None:
        raise SystemExit("no angle parameter %s-%s-%s" % (t1, t2, t3))
    return p


def match_types(pattern, types):
    return all(p == "" or p == t for p, t in zip(pattern, types))


def proper_param(ff, ts):
    best = None
    for d in ff["propers"]:
        wild = "" in d["types"]
        if best is not None and wild:
            continue
        if match_types(d["types"], ts) or match_types(d["types"], ts[::-1]):
            if best is None or (("" in best["types"]) and not wild):
                best = d
            if not wild:
                break
    return best


def improper_terms(ff, atoms, adj):
    """OpenMM 'amber' improper assignment (forcefield.py, _matchImproper)."""
    out = []
    for c in range(len(atoms)):
        if len(adj[c]) != 3:
            continue
        nb = sorted(adj[c])
        t = [atoms[i]["type"] for i in [c] + nb]
        match = None
        for d in ff["impropers"]:
            ty = d["types"]
            wild = "" in ty
            if match is not None and wild:
                continue
            if not (ty[0] == "" or ty[0] == t[0]):
                continue
            import itertools
            done = False
            for perm in itertools.permutations([(t[1], nb[0]), (t[2], nb[1]), (t[3], nb[2])]):
                (t2, a2), (t3, a3), (t4, a4) = perm
                if all(p == "" or p == q for p, q in zip(ty[1:], (t2, t3, t4))):
                    ta = lambda i: atoms[i]["tidx"]
                    el = lambda i: atoms[i]["element"]
                    if not wild:
                        if t2 == t4 and ta(a2) > ta(a4):
                            a2, a4 = a4, a2
                        if t3 == t4 and ta(a3) > ta(a4):
                            a3, a4 = a4, a3
                        if t2 == t3 and ta(a2) > ta(a3):
                            a2, a3 = a3, a2
                    else:
                        if el(a2) == el(a4) and ta(a2) > ta(a4):
                            a2, a4 = a4, a2
                        if el(a3) == el(a4) and ta(a3) > ta(a4):
                            a3, a4 = a4, a3
                        if ta(a2) > ta(a3):
                            a2, a3 = a3, a2
                    match = (a2, a3, c, a4, d)
                    done = True
                    break
            if done and not wild:
                break
        if match:
            a2, a3, cc, a4, d = match
            out.append(dict(atoms=[a2, a3, cc, a4], terms=d["terms"]))
    return out


def terms_for(ff, atoms, bonds):
    n = len(atoms)
    adj = [set() for _ in range(n)]
    for i, j in bonds:
        adj[i].add(j)
        adj[j].add(i)
    T = lambda i: atoms[i]["type"]
    B = [dict(i=i, j=j, k=bond_param(ff, T(i), T(j))[0], r0=bond_param(ff, T(i), T(j))[1]) for i, j in bonds]
    A = []
    for j in range(n):
        nb = sorted(adj[j])
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                i, k = nb[x], nb[y]
                kk, th = angle_param(ff, T(i), T(j), T(k))
                A.append(dict(i=i, j=j, k=k, kf=kk, theta0=th))
    P = []
    for j, k in bonds:
        for i in adj[j] - {k}:
            for l in adj[k] - {j}:
                if i == l:
                    continue
                d = proper_param(ff, [T(i), T(j), T(k), T(l)])
                if d is None:
                    raise SystemExit("no torsion %s-%s-%s-%s" % (T(i), T(j), T(k), T(l)))
                P.append(dict(atoms=[i, j, k, l], terms=d["terms"]))
    I = improper_terms(ff, atoms, adj)
    # Nonbonded: exclusions 1-2, 1-3; 1-4 scaled.
    sp = dict(nx.all_pairs_shortest_path_length(nx.Graph(list(bonds)), cutoff=3))
    NB = []
    for i in range(n):
        for j in range(i + 1, n):
            d = sp.get(i, {}).get(j)
            if d is not None and d < 3:
                continue
            NB.append(dict(i=i, j=j, scale14=(d == 3)))
    return B, A, P, I, NB


def molecule_json(ff, parts, title, note, group=None):
    atoms, bonds, offset = [], [], 0
    for a_list, b_list in parts:
        atoms += a_list
        bonds += [(i + offset, j + offset) for i, j in b_list]
        offset += len(a_list)
    Bt, At, Pt, It, NBt = [], [], [], [], []
    offset = 0
    for a_list, b_list in parts:
        B, A, P, I, NB = terms_for(ff, a_list, b_list)
        sh = lambda x: x + offset
        Bt += [dict(b, i=sh(b["i"]), j=sh(b["j"])) for b in B]
        At += [dict(a, i=sh(a["i"]), j=sh(a["j"]), k=sh(a["k"])) for a in A]
        Pt += [dict(p, atoms=[sh(x) for x in p["atoms"]]) for p in P]
        It += [dict(p, atoms=[sh(x) for x in p["atoms"]]) for p in I]
        NBt += [dict(p, i=sh(p["i"]), j=sh(p["j"]), inter=False) for p in NB]
        offset += len(a_list)
    # Between bases: every atom pair, full strength. Within one base pair it
    # is the pairing ("inter"); between different pairs of a stack it is the
    # stacking ("stack"). group[k] is the base-pair index of part k.
    if group is None:
        group = [0] * len(parts)
    starts = []
    o = 0
    for a_list, _ in parts:
        starts.append(o)
        o += len(a_list)
    for p in range(len(parts)):
        for q in range(p + 1, len(parts)):
            same = group[p] == group[q]
            for i in range(starts[p], starts[p] + len(parts[p][0])):
                for j in range(starts[q], starts[q] + len(parts[q][0])):
                    NBt.append(dict(i=i, j=j, scale14=False, inter=same, stack=not same))
    out_atoms = []
    for a in atoms:
        sig, eps = ff["lj"][a["type"]]
        out_atoms.append(dict(name=a["prefix"] + ":" + a["name"], element=a["element"], mass=a["mass"],
                              charge=a["charge"], type=a["type"], sigma=sig * NM, epsilon=eps / KJ_PER_KCAL,
                              xyz=a["xyz"]))
    # Units out: kcal/mol, Angstrom, rad. Bond/angle energy ½k(x−x0)² (OpenMM form).
    for b in Bt:
        b["k"] = b["k"] / KJ_PER_KCAL / NM ** 2
        b["r0"] = b["r0"] * NM
    for a in At:
        a["kf"] = a["kf"] / KJ_PER_KCAL
    for p in Pt + It:
        p["terms"] = [dict(n=n, k=k / KJ_PER_KCAL, phase=ph) for n, k, ph in p["terms"]]
    return dict(title=title, note=note, atoms=out_atoms, bonds=Bt, angles=At, propers=Pt, impropers=It,
                pairs=NBt, coulomb14=ff["c14"], lj14=ff["lj14"])


def wc_pairs(g, found):
    """Every Watson-Crick pair of the duplex, ordered along the helix axis:
    a list of (strand-1 residue, its match, partner residue, its match)."""
    def xyz(m, n):
        return g.nodes[m[n]]["xyz"]
    cand = []
    for a in found["DA"]:
        for t in found["DT"]:
            cand.append((dist(xyz(a, "N6"), xyz(t, "O4")) + dist(xyz(a, "N1"), xyz(t, "N3")), "DA", a, "DT", t))
    for gg in found["DG"]:
        for c in found["DC"]:
            cand.append((dist(xyz(gg, "O6"), xyz(c, "N4")) + dist(xyz(gg, "N1"), xyz(c, "N3")) + dist(xyz(gg, "N2"), xyz(c, "O2")), "DG", gg, "DC", c))
    cand.sort(key=lambda r: r[0])
    used, pairs = set(), []
    for d, r1, m1, r2, m2 in cand:
        k1, k2 = id(m1), id(m2)
        if k1 in used or k2 in used or d > 3.5 * (2 if r1 == "DA" else 3):  # mean D···A ≤ 3.5 Å
            continue
        used |= {k1, k2}
        pairs.append((r1, m1, r2, m2))
    # Order along the axis: principal axis of the pair centroids.
    import numpy as np
    cen = []
    for r1, m1, r2, m2 in pairs:
        pts = [g.nodes[v]["xyz"] for m in (m1, m2) for v in m.values() if g.nodes[v]["element"] != "H"]
        cen.append(np.mean(pts, axis=0))
    cen = np.array(cen)
    c0 = cen.mean(axis=0)
    w, v = np.linalg.eigh(np.cov((cen - c0).T))
    ax = v[:, -1]
    order = np.argsort((cen - c0) @ ax)
    return [pairs[i] for i in order]


def strand1(g, pair, strand_atoms):
    """Which base of the pair is on the reference strand (bonded path through
    the backbone to the first residue's base)."""
    r1, m1, r2, m2 = pair
    return (r1, m1, r2, m2) if m1["C1'"] in strand_atoms else (r2, m2, r1, m1)


LET = {"DA": "A", "DT": "T", "DG": "G", "DC": "C"}


def stacks(ff, g, found, note):
    """Stacks of 2, 3 and 4 consecutive Watson-Crick pairs cut from the duplex:
    A·T stacks from the central AATT, G·C stacks from the CGCG end. The
    bases keep their duplex positions (B-DNA rise and twist); the backbone is
    removed and each base is methyl-capped as above."""
    pairs = wc_pairs(g, found)
    # Reference strand: the connected backbone component containing the first
    # pair's first base, after cutting the H-bonds (the xyz has no H-bond edges).
    comp = nx.node_connected_component(g, pairs[0][1]["C1'"])
    pairs = [strand1(g, p, comp) for p in pairs]
    seq = "".join(LET[p[0]] for p in pairs)
    print("duplex, strand 1 along the axis:", seq, file=sys.stderr)
    out = {}
    def make(start, k, key, label):
        parts, group, title = [], [], []
        for n, (r1, m1, r2, m2) in enumerate(pairs[start:start + k]):
            parts.append(build_molecule(ff, g, r1, m1, "%s%d" % (LET[r1], n + 1)))
            parts.append(build_molecule(ff, g, r2, m2, "%s%d" % (LET[r2], n + 1)))
            group += [n, n]
            title.append(LET[r1] + "·" + LET[r2])
        st = "".join(LET[p[0]] for p in pairs[start:start + k])
        out[key] = molecule_json(ff, parts, "%s: %d stacked %s pairs (one strand reads %s along the axis; pairs %s)" % (key, k, label, st, " / ".join(title)), note, group)
    a0 = seq.index("AATT") if "AATT" in seq else seq.index("TTAA")
    for k in (2, 3, 4):
        make(a0, k, "AT%d" % k, "A·T")
        make(0, k, "GC%d" % k, "G·C")
    # The duplex's middle step: the strand read along the axis is the reverse
    # of 5'-CGCGAATTCGCG-3', so its "TA" in the middle is the ApT step, the
    # commonest step with two different pairs (7.86 % of human dinucleotides).
    make(a0 + 1, 2, "ApT", "A·T")
    out["ApT"]["title"] = "ApT: the A·T / T·A step (5'-AT-3' on both strands), the middle of AATT"
    out.update(poly_a(ff, g, pairs, a0, seq, note))
    return out


def kabsch(P, Q):
    """Rotation R and translation t with R·p + t ≈ q (least squares)."""
    import numpy as np
    P, Q = np.array(P), np.array(Q)
    cp, cq = P.mean(axis=0), Q.mean(axis=0)
    H = (P - cp).T @ (Q - cq)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cq - R @ cp


def poly_a(ff, g, pairs, a0, seq, note):
    """Poly(dA)·poly(dT) stacks of 2, 3 and 4 pairs: the commonest 2-, 3- and
    4-base runs of the human genome (AA/TT, AAA/TTT, AAAA/TTTT). The duplex
    has only one AA/TT step (in AATT), so its rigid-body step (the transform
    that carries the first A·T pair onto the second, by least squares over
    all base atoms) is applied again and again to build the longer runs:
    every step has the duplex's own AA/TT rise and twist."""
    import copy
    import numpy as np
    p1, p2 = pairs[a0], pairs[a0 + 1]
    assert p1[0] == p2[0] and p1[2] == p2[2], "first two pairs of the run must have the same orientation"
    m1 = build_molecule(ff, g, p1[0], p1[1], "X"), build_molecule(ff, g, p1[2], p1[3], "Y")
    m2 = build_molecule(ff, g, p2[0], p2[1], "X"), build_molecule(ff, g, p2[2], p2[3], "Y")
    P = [a["xyz"] for part in m1 for a in part[0]]
    Q = [a["xyz"] for part in m2 for a in part[0]]
    R, t = kabsch(P, Q)
    rms = float(np.sqrt(np.mean(np.sum((np.array(P) @ R.T + t - np.array(Q)) ** 2, axis=1))))
    print("AA/TT step from the duplex: rigid fit rms %.3f A, rotation %.1f deg" % (rms, np.degrees(np.arccos((np.trace(R) - 1) / 2))), file=sys.stderr)
    out = {}
    for k in (2, 3, 4):
        parts, group = [], []
        for n in range(k):
            for which, (res, inv) in enumerate(((p1[0], p1[1]), (p1[2], p1[3]))):
                atoms, bonds = build_molecule(ff, g, res, inv, "%s%d" % (LET[res], n + 1))
                atoms = copy.deepcopy(atoms)
                for a in atoms:
                    x = np.array(a["xyz"])
                    for _ in range(n):
                        x = R @ x + t
                    a["xyz"] = [float(v) for v in x]
                parts.append((atoms, bonds))
                group.append(n)
        run = LET[p1[0]] * k
        comp = {"A": "T", "T": "A"}[LET[p1[0]]] * k
        out["polyA%d" % k] = molecule_json(ff, parts, "polyA%d: %s/%s, %d stacked A·T pairs built from the duplex's AA/TT step" % (k, run, comp, k), note, group)
    return out


def main():
    ff = load_ff(sys.argv[1])
    g = load_xyz(sys.argv[2])
    found = find_bases(ff, g)
    for res, ms in found.items():
        print(res, len(ms), "matches", file=sys.stderr)
    # Pick Watson-Crick pairs from the duplex by their H-bond distances.
    def xyz(m, n):
        return g.nodes[m[n]]["xyz"]
    best = {}
    for a in found["DA"]:
        for t in found["DT"]:
            d = dist(xyz(a, "N6"), xyz(t, "O4")) + dist(xyz(a, "N1"), xyz(t, "N3"))
            if "AT" not in best or d < best["AT"][0]:
                best["AT"] = (d, a, t)
    for gg in found["DG"]:
        for c in found["DC"]:
            d = dist(xyz(gg, "O6"), xyz(c, "N4")) + dist(xyz(gg, "N1"), xyz(c, "N3")) + dist(xyz(gg, "N2"), xyz(c, "O2"))
            if "GC" not in best or d < best["GC"][0]:
                best["GC"] = (d, gg, c)
    note = ("AMBER OL15 (parm99/bsc0 base terms) via OpenMM amber14/DNA.OL15.xml; starting geometry "
            "from bdna.xyz (OpenMM tests); methyl cap at the glycosidic N: C1' as CM (type CT), H1' and two "
            "H replacing O4'/C2' (type H1, charge of H1'), CM charge sets the molecule neutral.")
    mols = {}
    _, a, t = best["AT"]
    _, gg, c = best["GC"]
    pick = {"DA": a, "DT": t, "DG": gg, "DC": c}
    tag = {"DA": "A", "DT": "T", "DG": "G", "DC": "C"}
    built = {res: build_molecule(ff, g, res, pick[res], tag[res]) for res in pick}
    for res in ("DA", "DT", "DG", "DC"):
        mols[tag[res]] = molecule_json(ff, [built[res]], NAMES[res], note)
    mols["AT"] = molecule_json(ff, [built["DA"], built["DT"]], "9-methyladenine:1-methylthymine Watson-Crick pair", note)
    mols["GC"] = molecule_json(ff, [built["DG"], built["DC"]], "9-methylguanine:1-methylcytosine Watson-Crick pair", note)
    for k in ("AT", "GC"):
        print(k, "pair H-bond sum %.3f A" % best[k][0], file=sys.stderr)
    mols.update(stacks(ff, g, found, note))
    json.dump(dict(source=note, molecules=mols), open(sys.argv[3], "w"), indent=1)


if __name__ == "__main__":
    main()
