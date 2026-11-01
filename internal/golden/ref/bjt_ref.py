# SPDX-License-Identifier: GPL-3.0-only OR LicenseRef-kronspice-Commercial
# Copyright (C) 2026 Franklin Amador
"""Independent references for the M7 golden files (BJT and controlled sources).

ngspice could not be installed where these files were made, so, as for the
M4 diode files, the references come from a second implementation that
shares no code with the Go simulator:

* the SPICE3 Gummel-Poon equations written out again from the published
  model (Gummel & Poon 1970; Massobrio & Antognetti, Semiconductor Device
  Modeling with SPICE, 2nd ed., ch. 2; SPICE3f3 User's Manual 4.3), as
  terminal currents and charges only -- no derivatives;
* a plain nodal formulation (KCL residuals plus voltage-source rows) whose
  Jacobian, small-signal conductances and capacitances are all taken by
  central finite differences, so the Go simulator's analytic derivatives
  and stamps are checked, not repeated;
* DC by damped Newton to machine precision; AC by solving (G + jwC) x = b
  with G and C from finite differences at the operating point; transient
  by the trapezoidal rule at a fixed step 50 to 200 times below the Go
  simulator's print step.

When ngspice is available, internal/golden/cmd/regen replaces these files
with ngspice's own output; the tolerances in the tests stay the same.

Run from the repository root:  python3 internal/golden/ref/bjt_ref.py
It needs Python 3 and numpy.
"""
import math
import os

import numpy as np

K, Q, T = 1.380649e-23, 1.602176634e-19, 300.15
VT = K * T / Q
GMIN = 1e-12
OUT = os.path.join("testdata", "golden")


# ---------------------------------------------------------------- circuit --

class Circuit:
    def __init__(self):
        self.nodes = {"0": -1}
        self.elems = []
        self.nbranch = 0

    def node(self, name):
        name = name.lower()
        if name not in self.nodes:
            self.nodes[name] = len([n for n in self.nodes.values() if n >= 0])
        return self.nodes[name]

    def internal(self, name):
        return self.node("#" + name)

    def size(self):
        return len([n for n in self.nodes.values() if n >= 0]) + self.nbranch

    def finish(self):
        nn = len([n for n in self.nodes.values() if n >= 0])
        for e in self.elems:
            if hasattr(e, "branch_local"):
                e.br = nn + e.branch_local
        return self

    def add(self, e):
        self.elems.append(e)
        if getattr(e, "needs_branch", False):
            e.branch_local = self.nbranch
            self.nbranch += 1
        return e

    def find(self, name):
        for e in self.elems:
            if e.name.lower() == name.lower():
                return e
        raise KeyError(name)

    def residual(self, x, t, dc_scale=1.0):
        """Static residual: currents leaving each node; source rows."""
        f = np.zeros(self.size())
        for e in self.elems:
            e.static(x, t, f, dc_scale)
        return f

    def charge(self, x):
        """Charge leaving each node into charge-storing branches."""
        qv = np.zeros(self.size())
        for e in self.elems:
            e.charge(x, qv)
        return qv

    def ac_rhs(self):
        b = np.zeros(self.size(), dtype=complex)
        for e in self.elems:
            e.ac(b)
        return b


def V(x, n):
    return 0.0 if n < 0 else x[n]


def addf(f, n, i):
    if n >= 0:
        f[n] += i


class Element:
    needs_branch = False

    def static(self, x, t, f, s):
        pass

    def charge(self, x, qv):
        pass

    def ac(self, b):
        pass


class R(Element):
    def __init__(self, c, name, a, b, r):
        self.name, self.a, self.b, self.r = name, c.node(a), c.node(b), r

    def static(self, x, t, f, s):
        i = (V(x, self.a) - V(x, self.b)) / self.r
        addf(f, self.a, i)
        addf(f, self.b, -i)

    def current(self, x):
        return (V(x, self.a) - V(x, self.b)) / self.r


class Cap(Element):
    def __init__(self, c, name, a, b, cv):
        self.name, self.a, self.b, self.c = name, c.node(a), c.node(b), cv

    def charge(self, x, qv):
        q = self.c * (V(x, self.a) - V(x, self.b))
        addf(qv, self.a, q)
        addf(qv, self.b, -q)


def pulse(args, t):
    v1, v2, td, tr, tf, pw, per = args
    if t <= td:
        return v1
    tt = (t - td) % per if per > 0 else t - td
    if tt < tr:
        return v1 + (v2 - v1) * tt / tr
    if tt < tr + pw:
        return v2
    if tt < tr + pw + tf:
        return v2 + (v1 - v2) * (tt - tr - pw) / tf
    return v1


class Vsrc(Element):
    needs_branch = True

    def __init__(self, c, name, p, n, dc=0.0, acmag=0.0, pulse_args=None):
        self.name, self.p, self.n = name, c.node(p), c.node(n)
        self.dc, self.acmag, self.pulse = dc, acmag, pulse_args

    def value(self, t):
        return self.dc if self.pulse is None or t is None else pulse(self.pulse, t)

    def static(self, x, t, f, s):
        i = x[self.br]  # current into the + terminal, through the source
        addf(f, self.p, i)
        addf(f, self.n, -i)
        f[self.br] += V(x, self.p) - V(x, self.n) - s * self.value(t)

    def ac(self, b):
        b[self.br] += self.acmag

    def current(self, x):
        return x[self.br]


class Isrc(Element):
    def __init__(self, c, name, p, n, dc):
        self.name, self.p, self.n, self.dc = name, c.node(p), c.node(n), dc

    def static(self, x, t, f, s):
        addf(f, self.p, s * self.dc)
        addf(f, self.n, -s * self.dc)


class Ctrl(Element):
    """E, G (voltage control) or F, H (current control through a V source)."""

    def __init__(self, c, kind, name, p, n, ctrl, gain):
        self.kind, self.name, self.p, self.n, self.gain = kind, name, c.node(p), c.node(n), gain
        self.needs_branch = kind in "EH"
        if kind in "EG":
            self.cp, self.cn = c.node(ctrl[0]), c.node(ctrl[1])
        else:
            self.src = ctrl
        self.c = c

    def control(self, x):
        if self.kind in "EG":
            return V(x, self.cp) - V(x, self.cn)
        return x[self.c.find(self.src).br]

    def static(self, x, t, f, s):
        u = self.control(x)
        if self.kind in "EH":
            i = x[self.br]
            addf(f, self.p, i)
            addf(f, self.n, -i)
            f[self.br] += V(x, self.p) - V(x, self.n) - self.gain * u
        else:
            i = self.gain * u
            addf(f, self.p, i)
            addf(f, self.n, -i)


# ------------------------------------------------------------ Gummel-Poon --

DEFAULTS = dict(IS=1e-16, BF=100.0, NF=1.0, VAF=0.0, IKF=0.0, ISE=0.0, NE=1.5, BR=1.0, NR=1.0, VAR=0.0,
                IKR=0.0, ISC=0.0, NC=2.0, RB=0.0, IRB=0.0, RBM=None, RE=0.0, RC=0.0, CJE=0.0, VJE=0.75,
                MJE=0.33, TF=0.0, XTF=0.0, VTF=0.0, ITF=0.0, CJC=0.0, VJC=0.75, MJC=0.33, XCJC=1.0, TR=0.0,
                CJS=0.0, VJS=0.75, MJS=0.0, FC=0.5)


def diode(i_s, v, nvt):
    """Ideal junction current with SPICE's cubic tail below -3 nvt."""
    if i_s == 0:
        return 0.0
    if v >= -3 * nvt:
        return i_s * (math.exp(min(v / nvt, 700)) - 1)
    a = (3 * nvt / (v * math.e)) ** 3
    return -i_s * (1 + a)


def depletion_q(v, cz, vj, m, fc):
    if cz == 0:
        return 0.0
    if v < fc * vj:
        return cz * vj * (1 - (1 - v / vj) ** (1 - m)) / (1 - m)
    # Linear capacitance continuation: integrate C(u) = cz/(1-fc)^(1+m) *
    # (1 - fc(1+m) + m u/vj) from fc*vj, plus the charge there.
    v0 = fc * vj
    q0 = cz * vj * (1 - (1 - fc) ** (1 - m)) / (1 - m)
    k = cz / (1 - fc) ** (1 + m)
    return q0 + k * ((1 - fc * (1 + m)) * (v - v0) + m / (2 * vj) * (v * v - v0 * v0))


class BJT(Element):
    def __init__(self, c, name, col, base, emit, model, typ, subst="0"):
        p = dict(DEFAULTS)
        p.update(model)
        if p["RBM"] is None:
            p["RBM"] = p["RB"]
        self.p, self.t, self.name = p, (1 if typ == "npn" else -1), name
        self.c, self.b, self.e, self.s = c.node(col), c.node(base), c.node(emit), c.node(subst)
        self.cp = c.internal(name + "c") if p["RC"] > 0 else self.c
        self.bp = c.internal(name + "b") if p["RB"] > 0 else self.b
        self.ep = c.internal(name + "e") if p["RE"] > 0 else self.e

    def core(self, vbe, vbc):
        p = self.p
        ibf = diode(p["IS"], vbe, p["NF"] * VT)
        ibr = diode(p["IS"], vbc, p["NR"] * VT)
        ile = diode(p["ISE"], vbe, p["NE"] * VT) + GMIN * vbe
        ilc = diode(p["ISC"], vbc, p["NC"] * VT) + GMIN * vbc
        inv = lambda a: 1 / a if a > 0 else 0.0
        q1 = 1 / (1 - vbc * inv(p["VAF"]) - vbe * inv(p["VAR"]))
        q2 = ibf * inv(p["IKF"]) + ibr * inv(p["IKR"])
        qb = q1 * (1 + math.sqrt(max(0.0, 1 + 4 * q2))) / 2
        return ibf, ibr, ile, ilc, qb

    def volts(self, x):
        t = self.t
        vb, vbp, vcp, vep, vs = V(x, self.b), V(x, self.bp), V(x, self.cp), V(x, self.ep), V(x, self.s)
        return t * (vbp - vep), t * (vbp - vcp), t * (vb - vcp), t * (vs - vcp)

    def static(self, x, t_, f, s):
        p, t = self.p, self.t
        vbe, vbc, _, _ = self.volts(x)
        ibf, ibr, ile, ilc, qb = self.core(vbe, vbc)
        it = (ibf - ibr) / qb
        ib1 = ibf / p["BF"] + ile
        ib2 = ibr / p["BR"] + ilc
        # Physical currents, leaving the internal nodes into the device.
        addf(f, self.cp, t * (it - ib2))
        addf(f, self.bp, t * (ib1 + ib2))
        addf(f, self.ep, t * (-it - ib1))
        # Substrate junction (ngspice bjtload.c): gmin (and ISS, default 0)
        # from the substrate to the collector (vertical, NPN default) or from
        # the base to the substrate (lateral, PNP default).
        sa, sb = (self.s, self.cp) if p.get("SUBS", t) > 0 else (self.bp, self.s)
        vsub = t * (V(x, sa) - V(x, sb))
        isub = GMIN * vsub
        addf(f, sa, t * isub)
        addf(f, sb, -t * isub)
        # Series resistances; the base one falls with current (RBM, IRB).
        for a, b, r in ((self.c, self.cp, p["RC"]), (self.e, self.ep, p["RE"])):
            if r > 0:
                i = (V(x, a) - V(x, b)) / r
                addf(f, a, i)
                addf(f, b, -i)
        if p["RB"] > 0:
            rbm, rb = p["RBM"], p["RB"]
            if getattr(self, "rbb_frozen", None) is not None:
                rbb = self.rbb_frozen
            elif p["IRB"] > 0 and rb != rbm:
                a1 = max((ib1 + ib2) / p["IRB"], 1e-9)
                z = (-1 + math.sqrt(1 + 14.59025 * a1)) / 2.4317 / math.sqrt(a1)
                rbb = rbm + 3 * (rb - rbm) * (math.tan(z) - z) / (z * math.tan(z) ** 2)
            else:
                rbb = rbm + (rb - rbm) / qb
            self.rbb_last = rbb
            i = (V(x, self.b) - V(x, self.bp)) / rbb
            addf(f, self.b, i)
            addf(f, self.bp, -i)

    def charges(self, x):
        p = self.p
        vbe, vbc, vbx, vcs = self.volts(x)
        ibf, ibr, _, _, qb = self.core(vbe, vbc)
        qbe = depletion_q(vbe, p["CJE"], p["VJE"], p["MJE"], p["FC"])
        if p["TF"] > 0:
            if vbe > 0:
                factor = 0.0
                if p["XTF"] > 0:
                    w = ibf / (ibf + p["ITF"]) if p["ITF"] > 0 else 1.0
                    e = math.exp(vbc / (1.44 * p["VTF"])) if p["VTF"] > 0 else 1.0
                    factor = p["XTF"] * w * w * e
                qbe += p["TF"] * ibf * (1 + factor) / qb
            else:
                qbe += p["TF"] * ibf
        qbc = p["TR"] * ibr + depletion_q(vbc, p["CJC"] * p["XCJC"], p["VJC"], p["MJC"], p["FC"])
        qbx = depletion_q(vbx, p["CJC"] * (1 - p["XCJC"]), p["VJC"], p["MJC"], p["FC"])
        cz, ps, ms = p["CJS"], p["VJS"], p["MJS"]
        if cz == 0:
            qcs = 0.0
        elif vcs < 0:
            qcs = ps * cz * (1 - (1 - vcs / ps) ** (1 - ms)) / (1 - ms)
        else:
            qcs = vcs * cz * (1 + ms * vcs / (2 * ps))
        return qbe, qbc, qbx, qcs

    def charge(self, x, qv):
        t = self.t
        qbe, qbc, qbx, qcs = self.charges(x)
        for a, b, qq in ((self.bp, self.ep, qbe), (self.bp, self.cp, qbc), (self.b, self.cp, qbx),
                         (self.s, self.cp, qcs)):
            addf(qv, a, t * qq)
            addf(qv, b, -t * qq)

    def collector_current(self, x):
        """External collector current into the device (DC)."""
        f = np.zeros(len(x))
        self.static(x, None, f, 1.0)
        if self.p["RC"] > 0:
            return (V(x, self.c) - V(x, self.cp)) / self.p["RC"]
        return f[self.c] if self.c >= 0 else 0.0


# ---------------------------------------------------------------- solvers --

def jac(fun, x, h=1e-7):
    n = len(x)
    J = np.zeros((n, n))
    for k in range(n):
        d = h * max(1.0, abs(x[k]))
        xp, xm = x.copy(), x.copy()
        xp[k] += d
        xm[k] -= d
        J[:, k] = (fun(xp) - fun(xm)) / (2 * d)
    return J


def newton(fun, x, maxit=400, vmax=0.1):
    for it in range(maxit):
        f = fun(x)
        J = jac(fun, x)
        dx = np.linalg.solve(J, -f)
        m = np.max(np.abs(dx))
        if m > vmax:
            dx *= vmax / m
        x = x + dx
        if m < 1e-14 * max(1.0, np.max(np.abs(x))):
            return x
    f = fun(x)
    if np.max(np.abs(f)) < 1e-15:
        return x
    raise RuntimeError("Newton did not converge: |f| = %g" % np.max(np.abs(f)))


def op(c, t=None, x0=None):
    x = np.zeros(c.size()) if x0 is None else x0
    if x0 is None:
        # Source stepping to a robust start, then the full solve.
        for s in np.linspace(0.1, 1.0, 10):
            x = newton(lambda y: c.residual(y, t, s), x)
    return newton(lambda y: c.residual(y, t), x)


def ac(c, x, freqs):
    # SPICE stamps the current-dependent base resistance in AC as the
    # conductance 1/rbb at the operating point (its DC secant value), with
    # no d(rbb)/dV term; hold rbb there while differentiating.
    c.residual(x, None)
    bjts = [e for e in c.elems if isinstance(e, BJT) and e.p["RB"] > 0]
    for e in bjts:
        e.rbb_frozen = e.rbb_last
    G = jac(lambda y: c.residual(y, None), x, 1e-6)
    for e in bjts:
        e.rbb_frozen = None
    C = jac(lambda y: c.charge(y), x, 1e-6)
    b = c.ac_rhs()
    return [np.linalg.solve(G + 2j * math.pi * f * C, b) for f in freqs]


def tran(c, stop, h, record_every):
    x = op(c, t=0.0)
    q = c.charge(x)
    iq = np.zeros_like(q)
    rows = [(0.0, x.copy())]
    nsteps = int(round(stop / h))
    for k in range(1, nsteps + 1):
        t = k * h
        qn, iqn = q, iq

        def fun(y):
            return c.residual(y, t) + 2 * (c.charge(y) - qn) / h - iqn

        x = newton(fun, x.copy(), vmax=0.5)
        qnew = c.charge(x)
        iq = 2 * (qnew - qn) / h - iqn
        q = qnew
        if k % record_every == 0:
            rows.append((t, x.copy()))
    return rows


def write(name, header, rows, note):
    path = os.path.join(OUT, name)
    with open(path, "w") as fh:
        fh.write("# reference: " + note + " (internal/golden/ref/bjt_ref.py)\n")
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(repr(float(v)) for v in r) + "\n")
    print("wrote", path)


def write_op(name, values, note):
    path = os.path.join(OUT, name)
    with open(path, "w") as fh:
        fh.write("# reference: " + note + " (internal/golden/ref/bjt_ref.py)\n")
        fh.write("name,value\n")
        for k, v in values:
            fh.write("%s,%r\n" % (k, float(v)))
    print("wrote", path)


def frange(a, b, s):
    n = int(round((b - a) / s))
    return [a + k * s for k in range(n + 1)]


# --------------------------------------------------------------- circuits --

QGP = dict(IS=2e-15, BF=150, NF=1.01, VAF=80, IKF=50e-3, ISE=5e-14, NE=1.6, BR=4, NR=1.005, VAR=25,
           IKR=20e-3, ISC=1e-13, NC=1.9, RB=120, RBM=40, IRB=0.5e-3, RE=0.8, RC=6)


def ce_op():
    c = Circuit()
    vcc = c.add(Vsrc(c, "VCC", "vcc", "0", 12))
    c.add(R(c, "R1", "vcc", "b", 47e3))
    c.add(R(c, "R2", "b", "0", 10e3))
    c.add(R(c, "RC", "vcc", "c", 3.3e3))
    c.add(R(c, "RE", "e", "0", 1e3))
    q = c.add(BJT(c, "Q1", "c", "b", "e", QGP, "npn"))
    c.finish()
    x = op(c)
    vals = [("v(%s)" % n, V(x, c.nodes[n])) for n in ("vcc", "b", "c", "e")]
    vals += [("i(vcc)", vcc.current(x)), ("i(q1)", q.collector_current(x))]
    write_op("bjt_ce_op.csv", vals, "damped Newton on the Gummel-Poon equations, finite-difference Jacobian")


def output_dc():
    model = dict(IS=1e-15, BF=120, VAF=60, IKF=30e-3, RC=10, RB=50, BR=2, ISE=1e-14, NE=1.5)
    rows = []
    x = None
    for vce in frange(0, 6, 0.1):
        c = Circuit()
        c.add(Isrc(c, "IB", "0", "b", 20e-6))
        vs = c.add(Vsrc(c, "VCE", "c", "0", vce))
        c.add(BJT(c, "Q1", "c", "b", "0", model, "npn"))
        c.finish()
        x = op(c, x0=x)
        rows.append((vce, V(x, c.nodes["b"]), vs.current(x)))
    write("bjt_output_dc.csv", ["sweep", "v(b)", "i(vce)"], rows,
          "damped Newton on the Gummel-Poon equations at each sweep point")


def pnp_dc():
    model = dict(IS=5e-15, BF=80, VAF=40, IKF=20e-3, BR=3, RB=100, RE=1, RC=5, ISE=2e-14, NE=1.7)
    rows = []
    x = None
    for vin in frange(-8, 10, 0.25):
        c = Circuit()
        c.add(Vsrc(c, "VCC", "vcc", "0", 10))
        c.add(Vsrc(c, "VN", "vn", "0", -5))
        vi = c.add(Vsrc(c, "VIN", "in", "0", vin))
        c.add(R(c, "RE", "vcc", "e", 2e3))
        c.add(R(c, "RL", "c", "vn", 1e3))
        c.add(BJT(c, "Q1", "c", "in", "e", model, "pnp"))
        c.finish()
        x = op(c, x0=x) if x is not None else op(c)
        rows.append((vin, V(x, c.nodes["e"]), V(x, c.nodes["c"]), vi.current(x)))
    write("bjt_pnp_dc.csv", ["sweep", "v(e)", "v(c)", "i(vin)"], rows,
          "damped Newton on the Gummel-Poon equations at each sweep point")


QAC = dict(IS=2e-15, BF=150, VAF=80, IKF=50e-3, RB=120, RBM=40, RE=0.8, RC=6, CJE=12e-12, VJE=0.75,
           MJE=0.35, TF=300e-12, XTF=3, VTF=4, ITF=40e-3, CJC=5e-12, VJC=0.6, MJC=0.4, XCJC=0.6, TR=10e-9,
           CJS=3e-12, VJS=0.7, MJS=0.3)


def ce_ac():
    c = Circuit()
    c.add(Vsrc(c, "VCC", "vcc", "0", 12))
    c.add(Vsrc(c, "VIN", "in", "0", 0, acmag=1))
    c.add(Cap(c, "CIN", "in", "b", 10e-6))
    c.add(R(c, "R1", "vcc", "b", 47e3))
    c.add(R(c, "R2", "b", "0", 10e3))
    c.add(R(c, "RC", "vcc", "c", 3.3e3))
    c.add(R(c, "RE", "e", "0", 1e3))
    c.add(Cap(c, "CE", "e", "0", 100e-6))
    c.add(BJT(c, "Q1", "c", "b", "e", QAC, "npn", subst="0"))
    c.finish()
    x = op(c)
    freqs = [10 ** (1 + k / 10) for k in range(0, 91)]  # .ac dec 10 10 10g
    xs = ac(c, x, freqs)
    n = c.nodes["c"]
    rows = [(f, abs(v[n]), math.degrees(math.atan2(v[n].imag, v[n].real))) for f, v in zip(freqs, xs)]
    write("bjt_ce_ac.csv", ["freq", "vm(c)", "vp(c)"], rows,
          "(G + jwC) solve with G, C by finite differences at the operating point")


def switch_tran():
    c = Circuit()
    c.add(Vsrc(c, "VCC", "vcc", "0", 5))
    c.add(Vsrc(c, "VIN", "in", "0", 0, pulse_args=(0, 5, 10e-9, 2e-9, 2e-9, 200e-9, 400e-9)))
    c.add(R(c, "RB", "in", "b", 4.7e3))
    c.add(R(c, "RC", "vcc", "c", 1e3))
    c.add(Cap(c, "CL", "c", "0", 2e-12))
    model = dict(IS=1e-15, BF=100, BR=1, TF=0.4e-9, TR=20e-9, CJE=4e-12, CJC=3e-12, VJC=0.65, MJC=0.4,
                 RB=30, RC=2)
    c.add(BJT(c, "Q1", "c", "b", "0", model, "npn"))
    c.finish()
    rows = tran(c, 600e-9, 5e-12, 200)  # record every 1 ns
    nc, nb = c.nodes["c"], c.nodes["b"]
    write("bjt_switch_tran.csv", ["time", "v(c)", "v(b)"], [(t, x[nc], x[nb]) for t, x in rows],
          "trapezoidal rule at a fixed 5 ps step, Newton with finite-difference Jacobian")


def controlled_op():
    # Exact by hand; checked against the nodal solve below.
    c = Circuit()
    c.add(Vsrc(c, "V1", "a", "0", 2))
    c.add(R(c, "R1", "a", "0", 1e3))
    c.add(R(c, "R2", "a", "m", 1e3))
    vs = c.add(Vsrc(c, "VS", "m", "0", 0))
    c.add(Ctrl(c, "E", "E1", "e", "0", ("a", "0"), 3))
    c.add(R(c, "RE", "e", "0", 1e3))
    c.add(Ctrl(c, "G", "G1", "0", "g", ("a", "0"), 2e-3))
    c.add(R(c, "RG", "g", "0", 500))
    c.add(Ctrl(c, "H", "H1", "h", "0", "VS", 500))
    c.add(R(c, "RH", "h", "0", 1e3))
    c.add(Ctrl(c, "F", "F1", "0", "f", "VS", 3))
    c.add(R(c, "RF", "f", "0", 100))
    c.finish()
    x = op(c)
    exact = {"e": 6.0, "g": 2.0, "h": 1.0, "f": 0.6, "a": 2.0, "m": 0.0}
    for k, v in exact.items():
        assert abs(V(x, c.nodes[k]) - v) < 1e-12, (k, V(x, c.nodes[k]))
    assert abs(vs.current(x) - 2e-3) < 1e-15
    vals = [("v(%s)" % k, exact[k]) for k in ("a", "e", "g", "h", "f")] + [("i(vs)", 2e-3)]
    write_op("controlled_op.csv", vals, "exact values (V1 = 2 V; E 3x, G 2 mS into 500 ohm, "
             "H 500 ohm and F 3x of the 2 mA in VS)")


if __name__ == "__main__":
    import sys
    only = set(sys.argv[1:])
    for fn in (controlled_op, ce_op, output_dc, pnp_dc, ce_ac, switch_tran):
        if not only or fn.__name__ in only:
            fn()
