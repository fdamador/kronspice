#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Prepare the RNA vibration example's molecules (molecules.json).

    python examples/rna_vibration/prep/prep.py RNA.OL3.xml 1eka.pdb examples/rna_vibration/molecules.json

Inputs (give their paths as arguments; neither is committed):

  RNA.OL3.xml  AMBER RNA force field OL3 as converted by OpenMM from AmberTools 17
               leaprc.RNA.OL3 (file's own DateGenerated 2017-11-25), from
               https://raw.githubusercontent.com/openmm/openmm/master/wrappers/python/openmm/app/data/amber14/RNA.OL3.xml
               downloaded 2026-11-07; last commit touching the file:
               e53bdc5eab600b4137abd80e544cb100d1acac68 (2021-02-01).
               Base terms are Cornell et al. (parm94/99: W. D. Cornell et al.,
               JACS 117, 5179, 1995), with the bsc0 alpha/gamma correction
               (A. Perez et al., Biophys. J. 92, 3817, 2007) and the OL3 chi
               correction (M. Zgarbova et al., JCTC 7, 2886, 2011).
               Units kJ/mol, nm, rad.
  1eka.pdb     PDB entry 1EKA, the solution NMR structure of the RNA duplex
               r(GAGUGCUC)2 with all hydrogens (X. Chen, J. A. McDowell,
               R. Kierzek, T. R. Krugh, D. H. Turner, Biochemistry 39, 8970,
               2000, doi:10.1021/bi992938e), from
               https://files.rcsb.org/download/1EKA.pdb. Self-complementary
               strands A and B pair residue i with 9 - i: G1·C8, A2·U7 and
               G3·C6 (Watson-Crick) and the tandem wobble pairs U4·G5 and
               G5·U4 (two hydrogen bonds each in this entry). Near A-form.

What it writes, per molecule: atoms (name, element, mass, charge, AMBER type,
coordinates in Angstrom), and every force-field term the Go code evaluates:
bonds, angles, proper and improper torsions, nonbonded pairs with their 1-4
scaling and exclusions. Assignment follows OpenMM's ForceField rules
(forcefield.py), as in examples/dna_vibration/prep/prep.py.

Molecules: A·U is A2(A)·U7(B), G·C is G1(A)·C8(B), the G·U wobble pair is
G5(A)·U4(B); the single bases are taken from the A·U and G·C pairs. Each
base is cut at the glycosidic bond and methyl-capped exactly as in the DNA
example (9-methyladenine, 1-methyluracil, 9-methylguanine,
1-methylcytosine): C1' becomes the methyl carbon (type CT), H1' stays, and
O4' and C2' are replaced by hydrogens 1.09 A from C1' along the old bond
directions; methyl hydrogens take type H1 and the charge of H1', and the
methyl carbon takes the charge that makes the molecule neutral. The cap is
this example's choice, not AMBER's.

Stacks. 1EKA has two steps between Watson-Crick pairs (G1·C8 -> A2·U7,
5'-GA-3'/5'-UC-3', and A2·U7 -> G3·C6) and no AA/UU, GG/CC or CG step, so
the stacks are built by repeating the first real A-form step: each pair gets a frame (origin at the
midpoint of its two C1' atoms, z the mean base-plane normal pointing along
strand A 5'->3', y from the strand-II C1' to the strand-I C1' made normal
to z, x = y x z; a simplified version of the standard reference frame of
Olson et al., J. Mol. Biol. 313, 229, 2001), the step is the rigid motion
from the G·C frame to the A·U frame, and a stack of n pairs places its
pairs at that step applied 0 ... n-1 times. A pair whose strand-I base is a
pyrimidine (U·A, C·G) is the A·U or G·C pair turned 180 degrees about its
own x axis (the duplex dyad). So every stack has the duplex's GA/UC rise,
twist and roll at every step, a stated approximation. The tandem G·U stack
GU2 is cut directly from 1EKA (U4·G5 / G5·U4, real coordinates).
"""
import copy
import itertools
import json
import math
import sys
import xml.etree.ElementTree as ET

import networkx as nx
import numpy as np

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


def load_pdb(path):
    """Residues of the first model: (chain, resnum) -> {"name": resname, "xyz": {atom: (x, y, z)}}."""
    res = {}
    for ln in open(path):
        if ln.startswith("ENDMDL"):
            break
        if not ln.startswith("ATOM"):
            continue
        name = ln[12:16].strip()
        rn = ln[17:20].strip()
        ch = ln[21]
        num = int(ln[22:26])
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        r = res.setdefault((ch, num), dict(name=rn, xyz={}))
        r["xyz"][name] = xyz
    return res


# Base atoms; the sugar atoms C1', H1', O4', C2' are kept for the cap.
BASE = {
    "A": ["N9", "C8", "H8", "N7", "C5", "C6", "N6", "H61", "H62", "N1", "C2", "H2", "N3", "C4"],
    "G": ["N9", "C8", "H8", "N7", "C5", "C6", "O6", "N1", "H1", "C2", "N2", "H21", "H22", "N3", "C4"],
    "U": ["N1", "C6", "H6", "C5", "H5", "C4", "O4", "N3", "H3", "C2", "O2"],
    "C": ["N1", "C6", "H6", "C5", "H5", "C4", "N4", "H41", "H42", "N3", "C2", "O2"],
}
RING = {
    "A": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
    "G": ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"],
    "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
    "C": ["N1", "C2", "N3", "C4", "C5", "C6"],
}
NAMES = {"A": "9-methyladenine", "G": "9-methylguanine", "U": "1-methyluracil", "C": "1-methylcytosine"}


def build_molecule(ff, res, coords, tag):
    """One methyl-capped base: res is the residue letter, coords its atom -> xyz map."""
    r = ff["residues"][res]
    info = {name: (typ, q) for name, typ, q in r["atoms"]}
    order = {name: i for i, (name, _, _) in enumerate(r["atoms"])}
    atoms = []
    for name in BASE[res] + ["C1'", "H1'"]:
        typ, q = info[name]
        atoms.append(dict(name=name, type=typ, charge=q, xyz=list(coords[name]), tidx=order[name]))
    c1 = coords["C1'"]
    extra = len(r["atoms"])
    for k, old in enumerate(["O4'", "C2'"]):
        p = coords[old]
        d = math.dist(c1, p)
        xyz = [c1[i] + 1.09 * (p[i] - c1[i]) / d for i in range(3)]
        atoms.append(dict(name="HM%d" % (k + 2), type="RNA-H1", charge=None, xyz=xyz, tidx=extra + k))
    for a in atoms:
        if a["name"] == "C1'":
            a["name"] = "CM"
        if a["name"] == "H1'":
            a["name"], a["type"] = "HM1", "RNA-H1"
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
    atoms, offset = [], 0
    for a_list, _ in parts:
        atoms += a_list
    Bt, At, Pt, It, NBt = [], [], [], [], []
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
    starts, o = [], 0
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
    for b in Bt:
        b["k"] = b["k"] / KJ_PER_KCAL / NM ** 2
        b["r0"] = b["r0"] * NM
    for a in At:
        a["kf"] = a["kf"] / KJ_PER_KCAL
    for p in Pt + It:
        p["terms"] = [dict(n=n, k=k / KJ_PER_KCAL, phase=ph) for n, k, ph in p["terms"]]
    return dict(title=title, note=note, atoms=out_atoms, bonds=Bt, angles=At, propers=Pt, impropers=It,
                pairs=NBt, coulomb14=ff["c14"], lj14=ff["lj14"])


# ---- base-pair frames and the A-form step ----

def plane_normal(pts):
    P = np.array(pts)
    _, _, vt = np.linalg.svd(P - P.mean(axis=0))
    return vt[-1]


def pair_frame(r1, c1, r2, c2, axis):
    """Frame (R, o) of a pair: columns of R are x, y, z (see the docstring)."""
    n1 = plane_normal([c1[a] for a in RING[r1]])
    n2 = plane_normal([c2[a] for a in RING[r2]])
    if np.dot(n1, n2) < 0:
        n2 = -n2
    z = n1 + n2
    z /= np.linalg.norm(z)
    if np.dot(z, axis) < 0:
        z = -z
    p1, p2 = np.array(c1["C1'"]), np.array(c2["C1'"])
    y = p1 - p2
    y -= np.dot(y, z) * z
    y /= np.linalg.norm(y)
    x = np.cross(y, z)
    return np.column_stack([x, y, z]), 0.5 * (p1 + p2)


def centroid(c, res):
    return np.mean([c[a] for a in RING[res]], axis=0)


def main():
    ff = load_ff(sys.argv[1])
    pdb = load_pdb(sys.argv[2])
    out_path = sys.argv[3]

    def res(ch, n):
        r = pdb[(ch, n)]
        return r["name"], r["xyz"]

    pairs = {"AU": (res("A", 2), res("B", 7)), "GC": (res("A", 1), res("B", 8)), "GU": (res("A", 5), res("B", 4))}
    for k, ((n1, _), (n2, _)) in pairs.items():
        assert n1 + n2 == k, (k, n1, n2)
    note = ("AMBER RNA OL3 (Cornell et al. base terms, bsc0, OL3 chi) via OpenMM amber14/RNA.OL3.xml; starting geometry "
            "from PDB 1EKA (NMR, r(GAGUGCUC)2 with hydrogens, Chen et al. 2000); methyl cap at the glycosidic N: C1' as CM "
            "(type CT), H1' and two H replacing O4'/C2' (type H1, charge of H1'), CM charge sets the molecule neutral.")
    built = {}
    for k, ((n1, c1), (n2, c2)) in pairs.items():
        built[k] = (build_molecule(ff, n1, c1, n1), build_molecule(ff, n2, c2, n2))
    mols = {}
    mols["A"] = molecule_json(ff, [built["AU"][0]], NAMES["A"], note)
    mols["U"] = molecule_json(ff, [built["AU"][1]], NAMES["U"], note)
    mols["G"] = molecule_json(ff, [built["GC"][0]], NAMES["G"], note)
    mols["C"] = molecule_json(ff, [built["GC"][1]], NAMES["C"], note)
    mols["AU"] = molecule_json(ff, list(built["AU"]), "9-methyladenine:1-methyluracil Watson-Crick pair", note)
    mols["GC"] = molecule_json(ff, list(built["GC"]), "9-methylguanine:1-methylcytosine Watson-Crick pair", note)
    mols["GU"] = molecule_json(ff, list(built["GU"]), "9-methylguanine:1-methyluracil wobble pair", note)

    # Frames: the duplex axis direction from G1·C8 to A2·U7 (strand A 5'->3').
    (_, cg), (_, cc) = pairs["GC"]
    (_, ca), (_, cu) = pairs["AU"]
    axis = (centroid(ca, "A") + centroid(cu, "U")) / 2 - (centroid(cg, "G") + centroid(cc, "C")) / 2
    Rg, og = pair_frame("G", cg, "C", cc, axis)
    Ra, oa = pair_frame("A", ca, "U", cu, axis)
    RS, oS = Rg.T @ Ra, Rg.T @ (oa - og)
    twist = math.degrees(math.atan2(RS[1, 0] - RS[0, 1], RS[0, 0] + RS[1, 1]))
    print("GA/UC step of 1EKA in the pair frame: displacement (%.2f, %.2f, %.2f) A, rotation %.1f deg (about z %.1f)"
          % (oS[0], oS[1], oS[2], math.degrees(math.acos((np.trace(RS) - 1) / 2)), twist), file=sys.stderr)
    local = {}
    for k, R, o in (("AU", Ra, oa), ("GC", Rg, og)):
        (n1, c1), (n2, c2) = pairs[k]
        local[k] = (built[k], R, o)
    flip = np.diag([1.0, -1.0, -1.0])

    def place(seq, key, title):
        """Stack read along strand I 5'->3' (e.g. "AAUU"): each letter's pair at the step applied n times."""
        parts, group = [], []
        R, o = Rg.copy(), og.copy()  # first pair sits where 1EKA's G1·C8 sits
        for n, s in enumerate(seq):
            ref = "AU" if s in "AU" else "GC"
            (m1, m2), R0, o0 = local[ref]
            flipped = s in "UC"
            order = (m2, m1) if flipped else (m1, m2)
            for atoms, bonds in order:
                atoms = copy.deepcopy(atoms)
                for a in atoms:
                    L = R0.T @ (np.array(a["xyz"]) - o0)
                    if flipped:
                        L = flip @ L
                    a["xyz"] = [float(v) for v in (o + R @ L)]
                    a["prefix"] = a["prefix"] + str(n + 1)
                parts.append((atoms, bonds))
                group.append(n)
            o = o + R @ oS
            R = R @ RS
        comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
        pairs_s = " / ".join(s + "·" + comp[s] for s in seq)
        mols[key] = molecule_json(ff, parts, title % dict(seq=seq, pairs=pairs_s, n=len(seq)), note, group)

    for n in (2, 3, 4):
        place("AAUU"[:n], "AU%d" % n, "AU%(n)d: %(n)d stacked A·U pairs (strand I reads %(seq)s; pairs %(pairs)s), 1EKA's GA/UC step repeated")
        place("CGCG"[:n], "GC%d" % n, "GC%(n)d: %(n)d stacked G·C pairs (strand I reads %(seq)s; pairs %(pairs)s), 1EKA's GA/UC step repeated")
        place("A" * n, "polyA%d" % n, "polyA%(n)d: %(seq)s/" + "U" * n + ", %(n)d stacked A·U pairs, 1EKA's GA/UC step repeated")
    place("AU", "ApU", "ApU: the A·U / U·A step (5'-AU-3' on both strands), 1EKA's GA/UC step")
    # Tandem wobble pairs, real coordinates: U4(A)·G5(B) and G5(A)·U4(B).
    p1 = (build_molecule(ff, "U", res("A", 4)[1], "U1"), build_molecule(ff, "G", res("B", 5)[1], "G1"))
    p2 = (build_molecule(ff, "G", res("A", 5)[1], "G2"), build_molecule(ff, "U", res("B", 4)[1], "U2"))
    mols["GU2"] = molecule_json(ff, [p1[0], p1[1], p2[0], p2[1]],
                                "GU2: the tandem wobble pairs U·G / G·U (5'-UG-3' on both strands) cut from 1EKA", note, [0, 0, 1, 1])
    json.dump(dict(source=note, molecules=mols), open(out_path, "w"), indent=1)
    print("wrote", ", ".join(mols), file=sys.stderr)


if __name__ == "__main__":
    main()
